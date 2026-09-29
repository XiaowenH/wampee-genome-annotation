#!/usr/bin/env bash
#
# ncRNA_annotation.sh - tRNA / rRNA / miRNA / snRNA / snoRNA annotation
#                       for the 27-chromosome triploid Clausena lansium genome
#
# Tools, and what each is actually for:
#   tRNAscan-SE 2.0   tRNA           (the standard; do NOT use Rfam for tRNA)
#   barrnap           rRNA           (fast HMM scan for 5S/5.8S/18S/28S)
#   RNAmmer / RfamScan   rRNA        (alternative / cross-check)
#   Infernal + Rfam   everything else (miRNA, snRNA, snoRNA, and rRNA too)
#
# Note on miRNA: Rfam/Infernal finds miRNA PRECURSOR families by covariance
# model. That is what genome papers report. Identifying which arm is the mature
# miRNA, and confirming expression, needs small-RNA sequencing - Rfam alone
# cannot do it. Say "miRNA precursors predicted by Rfam" in the paper, not
# "miRNAs identified".
#
set -euo pipefail

##====================================================================
## configuration  (see config/annotation.config.sh)
##====================================================================
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
CONFIG="${CONFIG:-$HERE/../config/annotation.config.sh}"
if [ ! -f "$CONFIG" ]; then
  echo "ERROR: config file not found: $CONFIG" >&2
  exit 1
fi
# shellcheck source=/dev/null
. "$CONFIG"
OUT="$NCRNA_OUT"
##====================================================================

mkdir -p "$OUT"
G=$(readlink -f "$GENOME")
[ -s "${G}.fai" ] || samtools faidx "$G"

# genome size in Mb, needed for cmscan's -Z (both strands => x2)
SIZE_MB=$(awk '{s+=$2} END{printf "%.2f", s*2/1000000}' "${G}.fai")
echo "-- genome: $(awk '{s+=$2} END{printf "%.1f Mb", s/1e6}' "${G}.fai"), -Z $SIZE_MB"

##--------------------------------------------------------------------
## 1. tRNA - tRNAscan-SE 2.0
##    -E = eukaryotic mode. Always keep the isotype/score output: the
##    pseudogene column matters, most raw hits in a plant genome are
##    tRNA-derived repeats rather than functional tRNAs.
##--------------------------------------------------------------------
if [ ! -s "$OUT/tRNA.out" ]; then
  tRNAscan-SE -E \
      -o "$OUT/tRNA.out" \
      -f "$OUT/tRNA.struct" \
      -s "$OUT/tRNA.isospecific" \
      -m "$OUT/tRNA.stats" \
      -b "$OUT/tRNA.bed" \
      -a "$OUT/tRNA.fa" \
      --thread "$CPU" \
      "$G"
fi

# split functional from pseudogenes
awk 'NR>3 && $0!~/^-/ {if ($10 ~ /pseudo/ || $9 < 20) print > "'"$OUT"'/tRNA.pseudo.txt";
                       else print > "'"$OUT"'/tRNA.high_conf.txt"}' "$OUT/tRNA.out" || true
echo "-- tRNA total     : $(awk 'NR>3' $OUT/tRNA.out | wc -l)"
echo "   high confidence: $(wc -l < $OUT/tRNA.high_conf.txt 2>/dev/null || echo 0)"
echo "   pseudo/low     : $(wc -l < $OUT/tRNA.pseudo.txt   2>/dev/null || echo 0)"

##--------------------------------------------------------------------
## 2. rRNA - barrnap
##--------------------------------------------------------------------
if [ ! -s "$OUT/rRNA.gff3" ]; then
  barrnap --kingdom euk --threads "$CPU" --reject 0.3 "$G" > "$OUT/rRNA.gff3"
fi
echo "-- rRNA:"
awk -F'\t' '!/^#/{split($9,a,";"); split(a[1],b,"="); print "   "b[2]}' "$OUT/rRNA.gff3" \
  | sort | uniq -c | sort -rn

##--------------------------------------------------------------------
## 3. everything else - Infernal against Rfam
##    Two-step: cmscan with --fmt 2 (clan competition) then filter.
##    This is the slow step; split the genome to parallelise properly.
##--------------------------------------------------------------------
if [ ! -s "$RFAM_DIR/Rfam.cm.i1f" ]; then
  echo "-- pressing Rfam.cm (once)"
  cmpress "$RFAM_DIR/Rfam.cm"
fi

if [ ! -s "$OUT/rfam.tblout" ]; then
  cmscan --cpu "$CPU" \
         --rfam --cut_ga --nohmmonly \
         --tblout "$OUT/rfam.tblout" \
         --fmt 2 \
         --clanin "$RFAM_DIR/Rfam.clanin" \
         -Z "$SIZE_MB" \
         -o "$OUT/rfam.cmscan.out" \
         "$RFAM_DIR/Rfam.cm" "$G"
fi

# drop hits marked as overlapping a better-scoring hit in the same clan
awk 'NR>2 && $0!~/^#/ && $20!="=" {print}' "$OUT/rfam.tblout" > "$OUT/rfam.dedup.tblout"
echo "-- Rfam hits: $(wc -l < $OUT/rfam.dedup.tblout) after clan competition"

##--------------------------------------------------------------------
## 4. summarise Rfam families by class
##--------------------------------------------------------------------
python3 - "$OUT" <<'PY'
import sys, collections, re
out = sys.argv[1]
cls = collections.Counter()
fam = collections.Counter()
rows = []
with open(f'{out}/rfam.dedup.tblout') as fh:
    for line in fh:
        if line.startswith('#') or not line.strip(): continue
        f = line.split()
        if len(f) < 18: continue
        target, acc, query = f[1], f[2], f[3]
        seq, start, end, strand = f[3], f[9], f[10], f[11]
        # --fmt 2 column order: idx target accession query accession clan ...
        rows.append(f)
        name = f[1]
        n = name.lower()
        if   n.startswith('mir') or 'mir-' in n: c = 'miRNA'
        elif 'rrna' in n or re.match(r'^(5s|5_8s|ssu|lsu)', n): c = 'rRNA'
        elif n.startswith('trna') or n == 'trna':  c = 'tRNA'
        elif n.startswith('u') and n[1:2].isdigit(): c = 'snRNA'
        elif 'sno' in n or n.startswith('sn'): c = 'snoRNA'
        elif ('ribozyme' in n or 'riboswitch' in n or 'hammerhead' in n
              or 'rnase_p' in n or n.startswith('rnasep')): c = 'ribozyme/riboswitch'
        else: c = 'other'
        cls[c] += 1
        fam[name] += 1
print("\n-- Rfam families by class:")
for k, v in cls.most_common():
    print(f"   {k:<22} {v}")
print("\n-- top 15 families:")
for k, v in fam.most_common(15):
    print(f"   {k:<22} {v}")
with open(f'{out}/rfam_class_counts.tsv','w') as o:
    o.write("class\tcount\n")
    for k, v in cls.most_common(): o.write(f"{k}\t{v}\n")
PY

##--------------------------------------------------------------------
## 5. one merged GFF3
##--------------------------------------------------------------------
python3 - "$OUT" <<'PY'
import sys, os
out = sys.argv[1]
recs = []

# tRNAscan-SE bed -> gff3
p = f'{out}/tRNA.bed'
if os.path.exists(p):
    for i, line in enumerate(open(p), 1):
        f = line.rstrip('\n').split('\t')
        if len(f) < 6: continue
        recs.append([f[0], 'tRNAscan-SE', 'tRNA', str(int(f[1])+1), f[2],
                     f[4] if len(f)>4 else '.', f[5], '.',
                     f'ID=tRNA{i};Name={f[3]}'])

# barrnap gff3
p = f'{out}/rRNA.gff3'
if os.path.exists(p):
    for line in open(p):
        if line.startswith('#'): continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9: continue
        recs.append(f)

# Rfam
p = f'{out}/rfam.dedup.tblout'
if os.path.exists(p):
    for i, line in enumerate(open(p), 1):
        if line.startswith('#') or not line.strip(): continue
        f = line.split()
        if len(f) < 18: continue
        name, acc = f[1], f[2]
        seq = f[3]
        s, e, strand = int(f[9]), int(f[10]), f[11]
        if strand == '-': s, e = e, s
        recs.append([seq, 'Infernal-Rfam', 'ncRNA', str(s), str(e),
                     f[16] if len(f)>16 else '.', strand, '.',
                     f'ID=rfam{i};Name={name};rfam_acc={acc}'])

def natkey(r):
    import re
    m = re.match(r'^\D*(\d+)(\D*)$', r[0])
    return (int(m.group(1)), m.group(2)) if m else (999, r[0])
recs.sort(key=lambda r: (natkey(r), int(r[3])))

with open(f'{out}/CL.ncRNA.gff3','w') as o:
    o.write('##gff-version 3\n')
    for r in recs: o.write('\t'.join(r) + '\n')
print(f"\nmerged: {out}/CL.ncRNA.gff3  ({len(recs)} features)")
PY

##--------------------------------------------------------------------
## 6. per-haplotype counts - the three should be similar
##--------------------------------------------------------------------
echo
echo "-- features per haplotype (A / B / C should be comparable):"
for t in tRNA rRNA ncRNA; do
  printf "   %-8s" "$t"
  for h in A B C; do
    n=$(awk -F'\t' -v t="$t" -v h="$h" '$3==t && $1 ~ h"$"' "$OUT/CL.ncRNA.gff3" | wc -l)
    printf " %s=%-6s" "$h" "$n"
  done
  echo
done
