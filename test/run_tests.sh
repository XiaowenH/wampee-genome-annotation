#!/usr/bin/env bash
#
# run_tests.sh - end-to-end smoke test on the toy dataset in test/data.
#
# Every script is run on inputs small enough to check by eye, and the results
# are asserted. No external tools, databases or reference genome are needed,
# so this runs anywhere in a few seconds (and in CI).
#
#   bash test/run_tests.sh
#
set -euo pipefail

HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
ROOT=$(cd -- "$HERE/.." && pwd)
DATA="$HERE/data"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

pass=0
fail=0

check() {  # check <description> <expected> <actual>
  if [ "$2" = "$3" ]; then
    printf '  ok   %-58s %s\n' "$1" "$2"
    pass=$((pass + 1))
  else
    printf '  FAIL %-58s expected %s, got %s\n' "$1" "$2" "$3"
    fail=$((fail + 1))
  fi
}

echo "== rename_genes.pl =="
perl "$ROOT/scripts/rename_genes.pl" -i "$DATA/toy.gff3" -o "$TMP/renamed.gff3" \
     -p Cla -m "$TMP/idmap.tsv" 2>"$TMP/rename.log"

check "genes renamed"          3 "$(awk -F'\t' '$3=="gene"' "$TMP/renamed.gff3" | wc -l)"
check "mRNAs renamed"          4 "$(awk -F'\t' '$3=="mRNA"' "$TMP/renamed.gff3" | wc -l)"
check "first gene on Chr5A"    Cla5AG00010 \
      "$(awk -F'\t' '$1=="Chr5A" && $3=="gene"' "$TMP/renamed.gff3" | head -1 | sed 's/.*ID=\([^;]*\).*/\1/')"
check "serial step of 10"      Cla5AG00020 \
      "$(awk -F'\t' '$1=="Chr5A" && $3=="gene"' "$TMP/renamed.gff3" | sed -n 2p | sed 's/.*ID=\([^;]*\).*/\1/')"
check "prev_id kept"           4 "$(grep -c 'prev_id=' <(awk -F'\t' '$3=="mRNA"' "$TMP/renamed.gff3"))"
check "id map rows"            7 "$(($(wc -l < "$TMP/idmap.tsv") - 1))"
# the longest-CDS isoform must become .1
long=$(grep -P '\tmRNA\t' "$TMP/renamed.gff3" | grep 'prev_id=evm.model.Chr5A.2b' | sed 's/.*ID=\([^;]*\).*/\1/')
check "longest isoform is .1"  "${long%.*}.1" "$long"

echo
echo "== clean_pep.py =="
python3 "$ROOT/scripts/clean_pep.py" -i "$DATA/toy.pep" -o "$TMP/clean.pep" \
        -d "$TMP/discarded.tsv" 2> "$TMP/clean.log"
check "proteins kept"          3 "$(grep -c '^>' "$TMP/clean.pep")"
check "internal stop dropped"  1 "$(grep -c 'internal_stop' "$TMP/discarded.tsv")"
check "terminal stop stripped" 0 "$(grep -c '\*' "$TMP/clean.pep")"

echo
echo "== sync_annotation.py =="
python3 "$ROOT/scripts/sync_annotation.py" -g "$DATA/toy.gff3" -p "$TMP/clean.pep" \
        -c "$DATA/toy.cds" -o "$TMP/sync" 2>"$TMP/sync.log"
check "mRNAs after sync"       3 "$(awk -F'\t' '$3=="mRNA"' "$TMP/sync.gff3" | wc -l)"
check "genes kept (partial)"   3 "$(awk -F'\t' '$3=="gene"' "$TMP/sync.gff3" | wc -l)"
check "CDS records filtered"   3 "$(grep -c '^>' "$TMP/sync.cds.fa")"
check "no orphan features"     "orphan features        : 0  (must be 0)" \
      "$(grep 'orphan features' "$TMP/sync.log" | sed 's/^ *//')"

echo
echo "== filter_repeat_gene.py =="
python3 "$ROOT/scripts/filter_repeat_gene.py" -g "$DATA/toy.gff3" \
        -k "$DATA/keep_gene.ids" -o "$TMP/kept.gff3" 2>"$TMP/filter.log"
check "genes kept"             2 "$(awk -F'\t' '$3=="gene"' "$TMP/kept.gff3" | wc -l)"
check "children follow parent" 2 "$(awk -F'\t' '$3=="mRNA"' "$TMP/kept.gff3" | wc -l)"
check "dropped gene absent"    0 "$(grep -c 'evm.TU.Chr5A.2' "$TMP/kept.gff3" || true)"

echo
echo "== lncRNA_filter.py =="
python3 "$ROOT/scripts/lncRNA_filter.py" -g "$DATA/toy.lncRNA.gff3" -o "$TMP/lnc.gff3" \
        --min-length 200 --report "$TMP/lnc.dropped.tsv" 2>"$TMP/lnc.log"
check "lncRNA genes kept"      1 "$(awk -F'\t' '$3=="gene"' "$TMP/lnc.gff3" | wc -l)"
check "short transcript dropped" 1 "$(($(wc -l < "$TMP/lnc.dropped.tsv") - 1))"

echo
echo "== shell scripts: syntax =="
for f in "$ROOT"/scripts/*.sh; do
  bash -n "$f" && printf '  ok   %s\n' "$(basename "$f")" && pass=$((pass + 1))
done

echo
echo "== shell scripts: no hard-coded absolute paths =="
if grep -nE '(^|[^A-Za-z_])(/home/|/Users/|\$HOME/soft)' "$ROOT"/scripts/*.sh; then
  echo "  FAIL hard-coded path found"; fail=$((fail + 1))
else
  echo "  ok   none found"; pass=$((pass + 1))
fi

echo
echo "-------------------------------------------"
printf '%d passed, %d failed\n' "$pass" "$fail"
[ "$fail" -eq 0 ]
