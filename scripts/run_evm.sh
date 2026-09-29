#!/usr/bin/env bash
#
# run_evm.sh - merge GeMoMa + ANNEVO + TransDecoder + StringTie (+ miniprot)
#              into one consensus gene set with EVidenceModeler v2
#
# Order of operations:
#   0. resolve every input to an ABSOLUTE path and check it exists
#   1-5. convert each source to EVM-compatible GFF3
#   6. validate each prediction file (EVM fails late and cryptically otherwise)
#   7. concatenate the gene predictions into ONE file
#   8. write weights.txt (TAB separated!)
#   9. run EVidenceModeler
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
GENOME_IN="$GENOME"



##====================================================================
## 0. resolve to absolute paths BEFORE any cd.
##    readlink -f fails under `set -e` when the file is missing, so this
##    doubles as an input check - far better than failing hours later.
##====================================================================
E=$EVM/EvmUtils
[ -x "$EVM/EVidenceModeler" ] || { echo "ERROR: EVM not found at $EVM"; exit 1; }

need() {   # need <label> <path>
  # NOTE: errors MUST go to stderr. This function's stdout is captured by
  # $( ), so an error echoed to stdout would end up inside the variable
  # instead of on your screen. Likewise `readlink -f` succeeds on paths
  # that do not exist, so the -s test is what actually catches typos.
  local p
  p=$(readlink -f "$2" 2>/dev/null) || {
      echo "ERROR: cannot resolve $1: $2" >&2; exit 1; }
  [ -s "$p" ] || {
      echo "ERROR: $1 is missing or empty: $2" >&2
      echo "       (resolved to $p)" >&2; exit 1; }
  echo "$p"
}

GENOME=$(need    "genome"       "$GENOME_IN")
GEMOMA=$(need    "GeMoMa"       "$GEMOMA_IN")
ANNEVO=$(need    "ANNEVO"       "$ANNEVO_IN")
TD=$(need        "TransDecoder" "$TD_IN")
STRINGTIE=$(need "StringTie"    "$STRINGTIE_IN")

# miniprot is optional: MINIPROT stays empty when the file is absent
MINIPROT=""
if [ -s "$MINIPROT_IN" ]; then
  MINIPROT=$(readlink -f "$MINIPROT_IN")
  echo "-- miniprot protein evidence: $MINIPROT"
else
  echo "-- no miniprot file at $MINIPROT_IN; running without protein evidence"
fi

WORK=$(pwd)
echo "-- work dir: $WORK"

mkdir -p evm_in
cd evm_in

##====================================================================
## helper: EVM's validator demands ID= on EVERY exon/CDS line,
##         not just Parent=. GeMoMa/ANNEVO/TransDecoder all omit it.
##====================================================================
cat > add_ids.pl <<'PERL'
#!/usr/bin/env perl
use strict; use warnings;
my %n;
while (<>) {
  if (/^#/ or !/\S/) { print; next }
  chomp; my @f = split /\t/;
  if (@f>=9 and ($f[2] eq 'exon' or $f[2] eq 'CDS') and $f[8] !~ /\bID=/) {
    my ($p) = $f[8] =~ /Parent=([^;]+)/;
    $f[8] = "ID=$p.$f[2]." . (++$n{"$p.$f[2]"}) . ";$f[8]" if defined $p;
  }
  print join("\t",@f), "\n";
}
PERL

##====================================================================
## 1. GeMoMa -> OTHER_PREDICTION
##====================================================================
echo "== converting GeMoMa =="
agat_convert_sp_gxf2gxf.pl -g "$GEMOMA" -o gemoma.agat.gff3
awk 'BEGIN{OFS="\t"} /^#/{print;next} NF>=9{$2="GeMoMa"; print}' gemoma.agat.gff3 \
  | perl add_ids.pl > gemoma.evm.gff3

##====================================================================
## 2. ANNEVO -> ABINITIO_PREDICTION
##====================================================================
echo "== converting ANNEVO =="
agat_convert_sp_gxf2gxf.pl -g "$ANNEVO" -o annevo.agat.gff3
awk 'BEGIN{OFS="\t"} /^#/{print;next} NF>=9{$2="ANNEVO"; print}' annevo.agat.gff3 \
  | perl add_ids.pl > annevo.evm.gff3

##====================================================================
## 3. TransDecoder -> OTHER_PREDICTION
##====================================================================
echo "== converting TransDecoder =="
agat_convert_sp_gxf2gxf.pl -g "$TD" -o td.agat.gff3
awk 'BEGIN{OFS="\t"} /^#/{print;next} NF>=9{$2="transdecoder"; print}' td.agat.gff3 \
  | perl add_ids.pl > transdecoder.evm.gff3

##====================================================================
## 4. StringTie GTF -> TRANSCRIPT alignments (a DIFFERENT input slot)
##====================================================================
echo "== converting StringTie =="
"$E/misc/align_GTF_to_align_GFF3.pl" "$STRINGTIE" assembler-stringtie \
  > transcripts.evm.gff3

##====================================================================
## 5. miniprot -> PROTEIN alignments (optional)
##    NOTE: the *_alignment_GFF3.py script, NOT miniprot_GFF_2_EVM_GFF3.py
##====================================================================
if [ -n "$MINIPROT" ]; then
  echo "== converting miniprot =="
  python3 "$E/misc/miniprot_GFF_2_EVM_alignment_GFF3.py" "$MINIPROT" \
    > proteins.evm.gff3
fi

##====================================================================
## 6. validate BEFORE running EVM
##====================================================================
for f in gemoma.evm.gff3 annevo.evm.gff3 transdecoder.evm.gff3; do
  echo "-- validating $f"
  "$E/gff3_gene_prediction_file_validator.pl" "$f" \
    || { echo "VALIDATION FAILED: $f"; exit 1; }
done

##====================================================================
## 7. one combined predictions file
##====================================================================
{ echo "##gff-version 3"
  grep -hv "^##gff-version" gemoma.evm.gff3 annevo.evm.gff3 transdecoder.evm.gff3
} > "$WORK/gene_predictions.gff3"

cp transcripts.evm.gff3 "$WORK/transcript_alignments.gff3"
[ -n "$MINIPROT" ] && cp proteins.evm.gff3 "$WORK/protein_alignments.gff3"

cd "$WORK"

##====================================================================
## 8. weights.txt - MUST be TAB separated; column 2 must match column 2
##    of the GFF3 exactly, or that evidence is silently ignored.
##====================================================================
{
  printf 'OTHER_PREDICTION\tGeMoMa\t%s\n' "$W_GEMOMA"
  printf 'OTHER_PREDICTION\ttransdecoder\t%s\n' "$W_TRANSDECODER"
  printf 'ABINITIO_PREDICTION\tANNEVO\t%s\n' "$W_ANNEVO"
  printf 'TRANSCRIPT\tassembler-stringtie\t%s\n' "$W_STRINGTIE"
  [ -n "$MINIPROT" ] && printf 'PROTEIN\tminiprot\t%s\n' "$W_MINIPROT"
} > weights.txt

echo "-- sources present in gene_predictions.gff3:"
awk -F'\t' '!/^#/ && NF>=9 {print "   " $2}' gene_predictions.gff3 | sort -u
echo "-- weights.txt (^I must appear between every column):"
cat -A weights.txt | sed 's/^/   /'

##====================================================================
## 9. run EVM
##====================================================================
EVM_ARGS=(
  --sample_id "$SAMPLE"
  --genome "$GENOME"
  --weights "$WORK/weights.txt"
  --gene_predictions "$WORK/gene_predictions.gff3"
  --transcript_alignments "$WORK/transcript_alignments.gff3"
  --segmentSize "$EVM_SEGMENT_SIZE"
  --overlapSize "$EVM_OVERLAP_SIZE"
  --CPU "$CPU"
)
[ -n "$MINIPROT" ] && EVM_ARGS+=( --protein_alignments "$WORK/protein_alignments.gff3" )

"$EVM/EVidenceModeler" "${EVM_ARGS[@]}"

echo "== done: ${SAMPLE}.EVM.gff3 / .pep / .cds / .bed =="
printf "genes: %s\n" "$(awk -F'\t' '$3=="gene"' "${SAMPLE}.EVM.gff3" | wc -l)"
