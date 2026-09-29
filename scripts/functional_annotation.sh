#!/usr/bin/env bash
#
# functional_annotation.sh - functional annotation of the Clausena lansium proteome
#
# Layers, in order of how much they should be trusted:
#   1. eggNOG-mapper  -> orthology-based GO / KEGG / COG / Pfam  (primary)
#   2. InterProScan   -> protein domains, independent of homology (primary)
#   3. Swiss-Prot     -> curated gene names and descriptions      (naming)
#   4. NR (plants)    -> "best hit description" column for the paper (descriptive only)
#
# Run everything on ONE representative isoform per gene, then propagate to the
# other isoforms. Annotating all 117,486 mRNAs separately wastes hours and gives
# near-identical results for isoforms of the same gene.
#
set -euo pipefail

##====================================================================
## configuration
##
## Values are taken from, in order of precedence:
##   1. environment variables set for this run
##   2. the config file named by $CONFIG
##   3. the defaults in that config file
##
## No paths are hard-coded in this script; edit the config file instead.
##====================================================================
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
CONFIG="${CONFIG:-$HERE/../config/annotation.config.sh}"
if [ ! -f "$CONFIG" ]; then
  echo "ERROR: config file not found: $CONFIG" >&2
  echo "       copy config/annotation.config.sh and set CONFIG=<your copy>" >&2
  exit 1
fi
# shellcheck source=/dev/null
. "$CONFIG"
OUT="$FUNC_OUT"



mkdir -p "$OUT"

##--------------------------------------------------------------------
## 0. representative isoform per gene
##    The renaming scheme sorted isoforms by CDS length, so ".1" is the
##    longest one for every gene.
##--------------------------------------------------------------------
if [ ! -s "$OUT/CL.rep.pep.fa" ]; then
  awk '/^>/{keep = ($0 ~ /\.1( |$)/)} keep' "$PEP" > "$OUT/CL.rep.pep.fa"
fi
echo "-- all proteins : $(grep -c '>' $PEP)"
echo "-- representative: $(grep -c '>' $OUT/CL.rep.pep.fa)"

# gffread writes a trailing '.' or '*' for the stop codon; strip it, some
# tools (notably InterProScan) reject '*' inside a sequence.
sed '/^>/!s/[*.]$//' "$OUT/CL.rep.pep.fa" > "$OUT/CL.rep.clean.pep.fa"
REP="$OUT/CL.rep.clean.pep.fa"

##--------------------------------------------------------------------
## 1. eggNOG-mapper  (uses eggnog_proteins.dmnd + eggnog.db internally)
##    Do NOT run diamond against eggnog_proteins.dmnd yourself - the raw
##    hits are meaningless without eggnog.db to map them to orthogroups.
##--------------------------------------------------------------------
if [ ! -s "$OUT/CL.emapper.annotations" ]; then
  emapper.py -i "$REP" -o CL --output_dir "$OUT" \
      --data_dir "$EGGNOG_DATA" \
      -m diamond --dmnd_db "$EGGNOG_DATA/eggnog_proteins.dmnd" \
      --itype proteins \
      --cpu "$CPU" \
      --tax_scope 33090 \
      --go_evidence non-electronic \
      --target_orthologs all \
      --evalue 1e-5 --pident 30 --query_cover 30 --subject_cover 30 \
      --override
fi
echo "-- eggNOG: $(grep -vc '^#' $OUT/CL.emapper.annotations || true) rows"

##--------------------------------------------------------------------
## 2. InterProScan - domains, and a second independent GO source
##--------------------------------------------------------------------
if command -v interproscan.sh >/dev/null && [ ! -s "$OUT/CL.iprscan.tsv" ]; then
  interproscan.sh -i "$REP" -f TSV,GFF3 -goterms -pa -iprlookup \
      -cpu "$CPU" -d "$OUT" -T "$OUT/ipr_tmp"
  mv "$OUT/$(basename $REP).tsv" "$OUT/CL.iprscan.tsv" 2>/dev/null || true
else
  echo "-- interproscan.sh not found, skipping (strongly recommended though)"
fi

##--------------------------------------------------------------------
## 3. Swiss-Prot - curated names. --max-target-seqs 1 is deliberate:
##    for naming you want the single best curated match, not a list.
##--------------------------------------------------------------------
if [ ! -s "$OUT/CL.sprot.tsv" ]; then
  diamond blastp -q "$REP" -d "$SPROT" -o "$OUT/CL.sprot.tsv" \
      --evalue 1e-5 --max-target-seqs 1 --threads "$CPU" \
      --outfmt 6 qseqid sseqid pident length evalue bitscore stitle \
      --very-sensitive
fi
echo "-- Swiss-Prot: $(cut -f1 $OUT/CL.sprot.tsv | sort -u | wc -l) proteins with a hit"

##--------------------------------------------------------------------
## 4. NR (Magnoliophyta) - description only, lowest trust
##--------------------------------------------------------------------
if [ ! -s "$OUT/CL.nr.tsv" ]; then
  diamond blastp -q "$REP" -d "$NRDB" -o "$OUT/CL.nr.tsv" \
      --evalue 1e-5 --max-target-seqs 1 --threads "$CPU" \
      --outfmt 6 qseqid sseqid pident length evalue bitscore stitle \
      --sensitive
fi
echo "-- NR: $(cut -f1 $OUT/CL.nr.tsv | sort -u | wc -l) proteins with a hit"

##--------------------------------------------------------------------
## 5. merge into one table, keyed on the representative mRNA ID
##--------------------------------------------------------------------
python3 - "$REP" "$OUT" <<'PY'
import sys, os, csv, re
rep, out = sys.argv[1], sys.argv[2]

ids = [l[1:].split()[0] for l in open(rep) if l.startswith('>')]

def load_dmnd(path):
    d = {}
    if not os.path.exists(path): return d
    for line in open(path):
        f = line.rstrip('\n').split('\t')
        if len(f) < 7: continue
        if f[0] in d: continue          # first hit = best hit
        d[f[0]] = {'sseqid': f[1], 'pident': f[2], 'evalue': f[4], 'desc': f[6]}
    return d

sprot = load_dmnd(f'{out}/CL.sprot.tsv')
nr    = load_dmnd(f'{out}/CL.nr.tsv')

egg = {}
p = f'{out}/CL.emapper.annotations'
if os.path.exists(p):
    hdr = None
    for line in open(p):
        if line.startswith('#query'):
            hdr = line.lstrip('#').rstrip('\n').split('\t'); continue
        if line.startswith('#') or not line.strip(): continue
        f = line.rstrip('\n').split('\t')
        if hdr is None: continue
        r = dict(zip(hdr, f))
        egg[r['query']] = r

ipr = {}
p = f'{out}/CL.iprscan.tsv'
if os.path.exists(p):
    for line in open(p):
        f = line.rstrip('\n').split('\t')
        if len(f) < 12: continue
        e = ipr.setdefault(f[0], {'ipr': set(), 'desc': set(), 'go': set()})
        if len(f) > 11 and f[11].startswith('IPR'):
            e['ipr'].add(f[11])
            if len(f) > 12 and f[12] not in ('', '-'): e['desc'].add(f[12])
        if len(f) > 13 and f[13].startswith('GO:'):
            e['go'].update(g for g in re.split(r'[|,]', f[13]) if g.startswith('GO:'))

cols = ['mRNA_id','gene_id','eggNOG_OG','eggNOG_desc','preferred_name',
        'GO','KEGG_ko','KEGG_Pathway','COG_category','PFAMs',
        'InterPro','InterPro_desc','InterPro_GO',
        'SwissProt_hit','SwissProt_desc','SwissProt_pident',
        'NR_hit','NR_desc','annotated']
w = csv.writer(open(f'{out}/CL.functional_annotation.tsv','w',newline=''), delimiter='\t')
w.writerow(cols)

n_any = 0
for q in ids:
    g = egg.get(q, {})
    i = ipr.get(q, {})
    s = sprot.get(q, {})
    n = nr.get(q, {})
    row = [
        q, q.rsplit('.',1)[0],
        g.get('eggNOG_OGs','-').split(',')[0] if g else '-',
        g.get('Description','-'), g.get('Preferred_name','-'),
        g.get('GOs','-'), g.get('KEGG_ko','-'), g.get('KEGG_Pathway','-'),
        g.get('COG_category','-'), g.get('PFAMs','-'),
        ';'.join(sorted(i.get('ipr',[]))) or '-',
        ';'.join(sorted(i.get('desc',[]))) or '-',
        ';'.join(sorted(i.get('go',[]))) or '-',
        s.get('sseqid','-'), s.get('desc','-'), s.get('pident','-'),
        n.get('sseqid','-'), n.get('desc','-'),
        '', ]
    has = any(x not in ('-','', None) for x in
              [row[2], row[10], row[13], row[16]])
    row[-1] = 'yes' if has else 'no'
    n_any += has
    w.writerow(row)

print(f"merged table: {out}/CL.functional_annotation.tsv")
print(f"proteins with at least one annotation: {n_any} / {len(ids)} "
      f"({100.0*n_any/len(ids):.1f}%)")
PY

##--------------------------------------------------------------------
## 6. propagate to every isoform and write the final GFF3
##--------------------------------------------------------------------
python3 - "$GFF" "$OUT" <<'PY'
import sys, csv, urllib.parse
gff, out = sys.argv[1], sys.argv[2]

ann = {}
with open(f'{out}/CL.functional_annotation.tsv') as fh:
    r = csv.DictReader(fh, delimiter='\t')
    for row in r:
        ann[row['gene_id']] = row      # keyed on gene, so all isoforms inherit

def esc(s):
    return urllib.parse.quote(s, safe='') if s and s != '-' else None

n = 0
with open(f'{out}/CL.functional.gff3','w') as o:
    for line in open(gff):
        if line.startswith('#'):
            o.write(line); continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9: o.write(line); continue
        if f[2] in ('gene','mRNA'):
            gid = None
            for kv in f[8].split(';'):
                if kv.startswith('ID='):
                    gid = kv[3:].rsplit('.',1)[0] if f[2]=='mRNA' else kv[3:]
            a = ann.get(gid)
            if a:
                extra = []
                for key, col in (('product','eggNOG_desc'),
                                 ('gene_name','preferred_name'),
                                 ('Ontology_term','GO'),
                                 ('kegg_ko','KEGG_ko'),
                                 ('pfam','PFAMs'),
                                 ('interpro','InterPro')):
                    v = esc(a.get(col,'-'))
                    if v: extra.append(f'{key}={v}')
                if extra:
                    f[8] = f[8] + ';' + ';'.join(extra); n += 1
        o.write('\t'.join(f) + '\n')
print(f"annotated features written: {n}  ->  {out}/CL.functional.gff3")
PY

##--------------------------------------------------------------------
## 7. summary for the paper
##--------------------------------------------------------------------
T="$OUT/CL.functional_annotation.tsv"
echo
echo "== functional annotation summary =="
awk -F'\t' 'NR>1{
  t++;
  if($3!="-")  e++;
  if($6!="-")  go++;
  if($7!="-")  ko++;
  if($11!="-") ip++;
  if($14!="-") sp++;
  if($17!="-") nr++;
  if($19=="yes") any++;
} END{
  printf "  total genes        : %d\n", t;
  printf "  eggNOG orthogroup  : %d (%.1f%%)\n", e,  e*100/t;
  printf "  GO terms           : %d (%.1f%%)\n", go, go*100/t;
  printf "  KEGG KO            : %d (%.1f%%)\n", ko, ko*100/t;
  printf "  InterPro domain    : %d (%.1f%%)\n", ip, ip*100/t;
  printf "  Swiss-Prot hit     : %d (%.1f%%)\n", sp, sp*100/t;
  printf "  NR hit             : %d (%.1f%%)\n", nr, nr*100/t;
  printf "  ANY annotation     : %d (%.1f%%)\n", any, any*100/t;
}' "$T"
