#!/usr/bin/env python3
"""
filter_repeat_gene.py - keep only whitelisted genes in a GFF3.

Reads a list of gene IDs to keep and writes a GFF3 containing those genes,
their transcripts and all transcript children (exon / CDS / UTR). Used to
remove recovered gene models that overlap transposable elements.

Usage:
  filter_repeat_gene.py -g missed.filtered.gff3 -k keep_gene.ids \\
                        -o missed.nonoverlap.gff3

  # keep IDs listed in the first column of a table, skipping a header line
  filter_repeat_gene.py -g in.gff3 -k keep.tsv --column 1 --skip-header \\
                        -o out.gff3
"""
import argparse
import collections
import re
import sys

ID_RE = re.compile(r"\bID=([^;]+)")
PARENT_RE = re.compile(r"\bParent=([^;]+)")


def read_ids(path, column=1, skip_header=False):
    """One ID per line, or a chosen column of a whitespace/TSV table."""
    ids = set()
    with open(path) as fh:
        for i, line in enumerate(fh):
            if skip_header and i == 0:
                continue
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t") if "\t" in line else line.split()
            if column - 1 < len(fields):
                ids.add(fields[column - 1])
    return ids


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-g", "--gff", required=True, help="input GFF3")
    ap.add_argument("-k", "--keep", required=True, help="file of gene IDs to keep")
    ap.add_argument("-o", "--out", required=True, help="output GFF3")
    ap.add_argument("--column", type=int, default=1,
                    help="column of --keep holding the ID [1]")
    ap.add_argument("--skip-header", action="store_true",
                    help="ignore the first line of --keep")
    ap.add_argument("--quiet", action="store_true", help="suppress the summary")
    args = ap.parse_args()

    keep_genes = read_ids(args.keep, args.column, args.skip_header)
    if not keep_genes:
        sys.exit(f"error: no IDs read from {args.keep}")

    kept_mrna = set()
    counts = collections.Counter()
    seen_genes = set()

    with open(args.gff) as fh, open(args.out, "w") as out:
        out.write("##gff-version 3\n")
        for line in fh:
            if line.startswith("#"):
                if not line.startswith("##gff-version"):
                    out.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            feature = fields[2]
            fid = ID_RE.search(fields[8])
            parent = PARENT_RE.search(fields[8])
            fid = fid.group(1) if fid else None
            parent = parent.group(1) if parent else None

            if feature == "gene":
                write = fid in keep_genes
                if write:
                    seen_genes.add(fid)
            elif feature in ("mRNA", "transcript"):
                write = parent in keep_genes
                if write and fid:
                    kept_mrna.add(fid)
            else:
                write = parent in kept_mrna

            if write:
                out.write("\t".join(fields) + "\n")
                counts[feature] += 1

    if not args.quiet:
        missing = keep_genes - seen_genes
        print(f"wrote {args.out}", file=sys.stderr)
        for feature, n in counts.most_common():
            print(f"   {feature:<18} {n}", file=sys.stderr)
        if missing:
            print(f"warning: {len(missing)} requested IDs were not found in the GFF3 "
                  f"(first 5: {sorted(missing)[:5]})", file=sys.stderr)


if __name__ == "__main__":
    main()
