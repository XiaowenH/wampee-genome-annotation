#!/usr/bin/env python3
"""
sync_annotation.py - make the GFF3 and CDS files consistent with a cleaned
                     protein FASTA.

After clean_pep.py removed transcripts with internal stop codons (and empty
sequences), the GFF3 and the CDS FASTA still contain those models. This script
removes them so all three files describe exactly the same transcript set.

What it does:
  1. reads the mRNA IDs that SURVIVED in the cleaned protein file
  2. drops every mRNA not in that set from the GFF3, together with its
     exon / CDS / UTR children
  3. drops a gene only when ALL of its transcripts were removed - a gene with
     one bad isoform and two good ones is kept, with the bad isoform gone
  4. filters the CDS FASTA to the same set
  5. reports what was removed and verifies the three files now agree

Usage:
  python3 sync_annotation.py -g CL.final.gff3 -p CL_clean.pep -c CL.cds \\
      -o CL.filtered

  # if you kept the discard table from clean_pep.py you can pass it instead
  python3 sync_annotation.py -g CL.final.gff3 -d CL.pep.discarded.tsv \\
      -c CL.cds -o CL.filtered

Outputs:
  <prefix>.gff3   <prefix>.cds.fa   <prefix>.removed.txt
"""
import argparse
import collections
import re
import sys

ap = argparse.ArgumentParser()
ap.add_argument('-g', '--gff', required=True, help='input GFF3')
ap.add_argument('-p', '--pep', help='CLEANED protein fasta (IDs to keep)')
ap.add_argument('-d', '--discarded', help='discard table from clean_pep.py (IDs to drop)')
ap.add_argument('-c', '--cds', help='CDS fasta to filter')
ap.add_argument('-o', '--out', default='filtered', help='output prefix')
ap.add_argument('--gene-suffix-strip', default=r'\.\d+$',
                help="regex turning an mRNA ID into its gene ID [default: '\\.\\d+$']")
args = ap.parse_args()

if not args.pep and not args.discarded:
    sys.exit("give either -p (cleaned protein fasta) or -d (discard table)")


def fasta_ids(path):
    ids = []
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                ids.append(line[1:].split()[0])
    return ids


# ------------------------------------------------------------------ keep set
keep = drop = None
if args.pep:
    keep = set(fasta_ids(args.pep))
    print(f"cleaned protein file : {len(keep)} transcripts to keep", file=sys.stderr)
else:
    drop = set()
    with open(args.discarded) as fh:
        for i, l in enumerate(fh):
            if i == 0 and l.lower().startswith('id\t'):
                continue
            f = l.rstrip('\n').split('\t')
            if f and f[0]:
                drop.add(f[0])
    print(f"discard table        : {len(drop)} transcripts to remove", file=sys.stderr)

# ------------------------------------------------------------------ pass 1
# Work out, for every gene, which of its transcripts survive.
mrna_of_gene = collections.defaultdict(list)
mrna_ids = []
with open(args.gff) as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9:
            continue
        if f[2] in ('mRNA', 'transcript'):
            m = re.search(r'\bID=([^;]+)', f[8])
            p = re.search(r'\bParent=([^;]+)', f[8])
            if not m:
                continue
            mid = m.group(1)
            gid = p.group(1) if p else re.sub(args.gene_suffix_strip, '', mid)
            mrna_of_gene[gid].append(mid)
            mrna_ids.append(mid)

print(f"GFF3                 : {len(mrna_ids)} mRNAs in {len(mrna_of_gene)} genes",
      file=sys.stderr)

if keep is None:
    keep = set(mrna_ids) - drop
keep_mrna = set(mrna_ids) & keep
drop_mrna = set(mrna_ids) - keep_mrna

# a gene goes only if every one of its transcripts goes
keep_gene, drop_gene = set(), set()
for gid, mids in mrna_of_gene.items():
    if any(m in keep_mrna for m in mids):
        keep_gene.add(gid)
    else:
        drop_gene.add(gid)

# genes that lost some but not all isoforms - worth reporting
partial = {g: [m for m in mids if m in drop_mrna]
           for g, mids in mrna_of_gene.items()
           if g in keep_gene and any(m in drop_mrna for m in mids)}

print(f"\nmRNAs kept           : {len(keep_mrna)}", file=sys.stderr)
print(f"mRNAs removed        : {len(drop_mrna)}", file=sys.stderr)
print(f"genes kept           : {len(keep_gene)}", file=sys.stderr)
print(f"genes removed        : {len(drop_gene)}  (all isoforms bad)", file=sys.stderr)
print(f"genes partly trimmed : {len(partial)}  (kept, one or more isoforms dropped)",
      file=sys.stderr)

# ------------------------------------------------------------------ pass 2
n_in = n_out = 0
kinds = collections.Counter()
with open(args.gff) as fh, open(f'{args.out}.gff3', 'w') as o:
    o.write('##gff-version 3\n')
    for line in fh:
        if line.startswith('#'):
            if not line.startswith('##gff-version'):
                o.write(line)
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9:
            continue
        n_in += 1
        mid = re.search(r'\bID=([^;]+)', f[8])
        pid = re.search(r'\bParent=([^;]+)', f[8])
        mid = mid.group(1) if mid else None
        pid = pid.group(1) if pid else None

        if f[2] == 'gene':
            ok = mid in keep_gene
        elif f[2] in ('mRNA', 'transcript'):
            ok = mid in keep_mrna
        else:
            # exon / CDS / UTR: keep only if its parent transcript survives
            ok = pid in keep_mrna
        if ok:
            o.write('\t'.join(f) + '\n')
            n_out += 1
            kinds[f[2]] += 1

print(f"\nGFF3 lines: {n_in} -> {n_out}", file=sys.stderr)
print(f"wrote {args.out}.gff3", file=sys.stderr)
for k, v in kinds.most_common():
    print(f"   {k:<18} {v}", file=sys.stderr)

# ------------------------------------------------------------------ CDS
if args.cds:
    n_seq = n_kept = 0
    with open(args.cds) as fh, open(f'{args.out}.cds.fa', 'w') as o:
        write = False
        for line in fh:
            if line.startswith('>'):
                sid = line[1:].split()[0]
                write = sid in keep_mrna
                n_seq += 1
                if write:
                    n_kept += 1
            if write:
                o.write(line)
    print(f"\nCDS: {n_seq} -> {n_kept}   wrote {args.out}.cds.fa", file=sys.stderr)

# ------------------------------------------------------------------ report
with open(f'{args.out}.removed.txt', 'w') as o:
    o.write("#type\tid\tnote\n")
    for m in sorted(drop_mrna):
        gid = re.sub(args.gene_suffix_strip, '', m)
        note = 'gene_removed' if gid in drop_gene else 'isoform_only'
        o.write(f"mRNA\t{m}\t{note}\n")
    for g in sorted(drop_gene):
        o.write(f"gene\t{g}\tall_isoforms_removed\n")
print(f"removed list: {args.out}.removed.txt", file=sys.stderr)

# ------------------------------------------------------------------ verify
print("\n== consistency check ==", file=sys.stderr)
gff_mrna = set()
with open(f'{args.out}.gff3') as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) >= 9 and f[2] in ('mRNA', 'transcript'):
            m = re.search(r'\bID=([^;]+)', f[8])
            if m:
                gff_mrna.add(m.group(1))
print(f"  mRNAs in filtered GFF3 : {len(gff_mrna)}", file=sys.stderr)
if args.pep:
    pep = set(fasta_ids(args.pep))
    print(f"  proteins               : {len(pep)}", file=sys.stderr)
    print(f"  GFF3 == protein set    : {gff_mrna == pep}", file=sys.stderr)
    if gff_mrna != pep:
        only_g = gff_mrna - pep
        only_p = pep - gff_mrna
        if only_g:
            print(f"    only in GFF3 ({len(only_g)}): {list(only_g)[:5]}", file=sys.stderr)
        if only_p:
            print(f"    only in pep  ({len(only_p)}): {list(only_p)[:5]}", file=sys.stderr)
if args.cds:
    cds = set(fasta_ids(f'{args.out}.cds.fa'))
    print(f"  CDS                    : {len(cds)}", file=sys.stderr)
    print(f"  GFF3 == CDS set        : {gff_mrna == cds}", file=sys.stderr)

# orphan check: every child must have a surviving parent
orph = 0
seen = set(gff_mrna) | keep_gene
with open(f'{args.out}.gff3') as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) >= 9:
            p = re.search(r'\bParent=([^;]+)', f[8])
            if p and p.group(1) not in seen:
                orph += 1
print(f"  orphan features        : {orph}  (must be 0)", file=sys.stderr)
