#!/usr/bin/env python3
"""
lncRNA_filter.py - keep lncRNA transcripts whose mature length passes a cutoff.

Mature length is the sum of exon lengths of a transcript, which is the length
that matters for the >= 200 nt lncRNA definition; the genomic span is longer
because it includes introns. A gene is kept when at least one of its
transcripts passes, and only the passing transcripts are written.

Usage:
  lncRNA_filter.py -g CL.lncRNA.final.gff3 -o CL_lncRNA_clean.gff3
  lncRNA_filter.py -g in.gff3 -o out.gff3 --min-length 200 \\
                   --report dropped.tsv
"""
import argparse
import collections
import re
import sys

ID_RE = re.compile(r"\bID=([^;]+)")
PARENT_RE = re.compile(r"\bParent=([^;]+)")
TRANSCRIPT_TYPES = ("mRNA", "transcript", "lnc_RNA", "lncRNA", "ncRNA")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-g", "--gff", required=True, help="input lncRNA GFF3")
    ap.add_argument("-o", "--out", required=True, help="output GFF3")
    ap.add_argument("--min-length", type=int, default=200,
                    help="minimum mature (exonic) transcript length in nt [200]")
    ap.add_argument("--report", help="optional TSV listing the dropped transcripts")
    ap.add_argument("--gene-suffix-strip", default=r"\.\d+$",
                    help="regex turning a transcript ID into its gene ID, used only "
                         r"when a transcript has no Parent= [default: '\.\d+$']")
    args = ap.parse_args()

    exon_len = collections.defaultdict(int)
    gene_of = {}
    transcripts = []

    with open(args.gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            feature, attrs = fields[2], fields[8]
            if feature == "exon":
                parent = PARENT_RE.search(attrs)
                if parent:
                    try:
                        exon_len[parent.group(1)] += int(fields[4]) - int(fields[3]) + 1
                    except ValueError:
                        sys.exit(f"error: non-numeric coordinates: {line.rstrip()}")
            elif feature in TRANSCRIPT_TYPES:
                tid = ID_RE.search(attrs)
                if not tid:
                    continue
                tid = tid.group(1)
                parent = PARENT_RE.search(attrs)
                gene_of[tid] = parent.group(1) if parent else re.sub(
                    args.gene_suffix_strip, "", tid)
                transcripts.append(tid)

    if not transcripts:
        sys.exit(f"error: no transcript features found in {args.gff}")

    keep_tx = {t for t in transcripts if exon_len[t] >= args.min_length}
    drop_tx = set(transcripts) - keep_tx
    keep_genes = {gene_of[t] for t in keep_tx}

    counts = collections.Counter()
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
            feature, attrs = fields[2], fields[8]
            fid = ID_RE.search(attrs)
            parent = PARENT_RE.search(attrs)
            fid = fid.group(1) if fid else None
            parent = parent.group(1) if parent else None

            if feature == "gene":
                write = fid in keep_genes
            elif feature in TRANSCRIPT_TYPES:
                write = fid in keep_tx
            else:
                write = parent in keep_tx

            if write:
                out.write("\t".join(fields) + "\n")
                counts[feature] += 1

    if args.report:
        with open(args.report, "w") as rep:
            rep.write("transcript_id\tgene_id\tmature_length\treason\n")
            for t in sorted(drop_tx):
                rep.write(f"{t}\t{gene_of[t]}\t{exon_len[t]}\t"
                          f"shorter_than_{args.min_length}nt\n")

    print(f"kept {len(keep_genes)} genes / {len(keep_tx)} transcripts "
          f"(dropped {len(drop_tx)} shorter than {args.min_length} nt)", file=sys.stderr)
    for feature, n in counts.most_common():
        print(f"   {feature:<18} {n}", file=sys.stderr)


if __name__ == "__main__":
    main()
