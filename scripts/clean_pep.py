#!/usr/bin/env python3
"""
clean_pep.py - strip the terminal stop codon and drop proteins with an
               internal stop.

gffread writes the stop codon as a trailing '*' (or '.'), which several
downstream tools reject. A stop in the MIDDLE of a sequence means the gene
model is wrong, so those transcripts are removed and listed in a discard
table; feed that table, or the cleaned FASTA, to sync_annotation.py so the
GFF3 and CDS files are filtered to the same transcript set.

Usage:
  clean_pep.py -i CL.pep -o CL_clean.pep -d CL.pep.discarded.tsv

  # the original positional form still works:
  clean_pep.py CL.pep CL_clean.pep CL.pep.discarded.tsv
"""
import argparse
import sys


def read_fasta(path):
    """Yield (header, sequence) pairs. Header keeps everything after '>'."""
    name = None
    seq = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name = line[1:]
                seq = []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def parse_args():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--in", dest="inp", help="input protein FASTA")
    ap.add_argument("-o", "--out", help="cleaned protein FASTA")
    ap.add_argument("-d", "--discarded", help="TSV listing the removed sequences")
    ap.add_argument("--stop-chars", default="*.",
                    help="characters treated as a stop codon [*.]")
    ap.add_argument("--min-length", type=int, default=1,
                    help="drop proteins shorter than this after stripping [1]")
    ap.add_argument("positional", nargs="*",
                    help="legacy form: <in> <out> <discarded>")
    args = ap.parse_args()

    if args.positional:
        if len(args.positional) != 3:
            ap.error("positional form needs exactly three files: in out discarded")
        args.inp, args.out, args.discarded = args.positional
    if not (args.inp and args.out and args.discarded):
        ap.error("give -i, -o and -d (or the three positional files)")
    return args


def main():
    args = parse_args()
    stops = set(args.stop_chars)
    kept = dropped_internal = dropped_short = 0

    with open(args.out, "w") as out, open(args.discarded, "w") as bad:
        bad.write("id\treason\tlength\tstop_positions\n")
        for name, seq in read_fasta(args.inp):
            seq_id = name.split()[0]
            seq = seq.rstrip(args.stop_chars)          # terminal stop only
            positions = [i + 1 for i, c in enumerate(seq) if c in stops]

            if positions:
                dropped_internal += 1
                bad.write(f"{seq_id}\tinternal_stop\t{len(seq)}\t"
                          f"{','.join(map(str, positions))}\n")
                continue
            if len(seq) < args.min_length:
                dropped_short += 1
                reason = "empty" if not seq else f"shorter_than_{args.min_length}"
                bad.write(f"{seq_id}\t{reason}\t{len(seq)}\t-\n")
                continue

            out.write(f">{name}\n{seq}\n")
            kept += 1

    print(f"kept={kept}  dropped_internal_stop={dropped_internal}  "
          f"dropped_short_or_empty={dropped_short}", file=sys.stderr)
    if kept == 0:
        sys.exit("error: no sequences survived - check the input file")


if __name__ == "__main__":
    main()
