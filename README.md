# clausena-genome-annotation

Annotation utilities used to build the gene, functional and non-coding RNA
annotation of the haplotype-resolved, telomere-to-telomere genome of seedless
wampee — *Clausena lansium* (Lour.) Skeels cv. 'Yunan'.

The repository holds the in-house scripts referenced in the Code Availability
statement of the accompanying Data Descriptor. Third-party tools (EVidenceModeler,
eggNOG-mapper, InterProScan, Infernal, …) are **not** bundled; the shell pipelines
drive them and are documented below.

[![tests](https://github.com/XiaowenH/wampee-genome-annotation/actions/workflows/ci.yml/badge.svg)](https://github.com/XiaowenH/wampee-genome-annotation/actions/workflows/ci.yml)

---

## What the scripts do

| Script | Language | Purpose |
|---|---|---|
| `scripts/rename_genes.pl` | Perl | Assign systematic IDs (`Cla05AG00160`, `.1` isoforms) to a PASA/EVM GFF3, keeping the old ID as `prev_id=` and writing an old→new mapping table |
| `scripts/clean_pep.py` | Python | Strip terminal stop codons; drop proteins containing an internal stop, with a discard table |
| `scripts/sync_annotation.py` | Python | Make GFF3, CDS and protein files describe exactly the same transcript set, then verify it |
| `scripts/filter_repeat_gene.py` | Python | Keep only whitelisted genes (used to drop recovered models overlapping transposable elements) |
| `scripts/lncRNA_filter.py` | Python | Keep lncRNA transcripts whose **mature** (exonic) length passes a cutoff |
| `scripts/run_evm.sh` | Bash | Convert GeMoMa / ANNEVO / TransDecoder / StringTie / miniprot evidence to EVM format, validate it, and run EVidenceModeler |
| `scripts/functional_annotation.sh` | Bash | eggNOG-mapper, InterProScan, Swiss-Prot and NR annotation of one representative isoform per gene, merged into one table |
| `scripts/ncRNA_annotation.sh` | Bash | tRNAscan-SE and Infernal/Rfam scan with clan competition, filtered to a high-confidence non-coding RNA set |

Order of use in the published pipeline:

```
run_evm.sh  →  PASA refinement  →  rename_genes.pl  →  clean_pep.py
     →  sync_annotation.py  →  functional_annotation.sh
     →  ncRNA_annotation.sh  →  lncRNA_filter.py
```

---

## Installation

The five stand-alone scripts need only Python ≥ 3.8, Perl ≥ 5.32 and coreutils.
For a reproducible environment, including the third-party tools the shell
pipelines call:

```bash
git clone https://github.com/XiaowenH/wampee-genome-annotation.git
cd clausena-genome-annotation

conda env create -f environment.yml     # or: mamba env create -f environment.yml
conda activate clausena-annotation

bash test/run_tests.sh                  # ~5 s, no external data needed
```

`test/run_tests.sh` runs every script on the toy dataset in `test/data/` and
asserts the results — 23 checks. Run it after any change; CI runs the same suite
plus `shellcheck`, `flake8` and `perl -c`.

---

## Configuration

**No paths are hard-coded in the scripts.** Everything the shell pipelines need
lives in `config/annotation.config.sh`. Copy it, edit the copy, and point the
scripts at it:

```bash
cp config/annotation.config.sh config/my_run.local.sh
$EDITOR config/my_run.local.sh
CONFIG=config/my_run.local.sh scripts/run_evm.sh
```

Precedence is **environment variable → config file → default**, so a single
setting can be changed for one run without touching the file:

```bash
CPU=32 CONFIG=config/my_run.local.sh scripts/functional_annotation.sh
```

`config/*.local.sh` is git-ignored, so cluster-specific paths never reach the
repository. The main settings:

| Variable | Meaning | Default |
|---|---|---|
| `SAMPLE` | prefix for output files | `CL` |
| `SPECIES_PREFIX` | gene ID prefix | `Cla` |
| `CPU` | threads | `8` |
| `GENOME` | soft-masked assembly FASTA | `data/CL_v1.fa` |
| `EVM` | EVidenceModeler installation directory | `/opt/EVidenceModeler-v2.1.0` |
| `GEMOMA_IN`, `ANNEVO_IN`, `TD_IN`, `STRINGTIE_IN`, `MINIPROT_IN` | evidence files for EVM (`MINIPROT_IN` optional) | `data/…` |
| `W_GEMOMA`, `W_TRANSDECODER`, `W_ANNEVO`, `W_STRINGTIE`, `W_MINIPROT` | EVM evidence weights | 8, 7, 4, 10, 5 |
| `EGGNOG_DATA`, `SPROT`, `NRDB`, `RFAM_DIR` | database locations | `db/…` |
| `LNCRNA_MIN_LENGTH` | mature-length cutoff for lncRNAs | `200` |

The Python and Perl scripts take all their settings as command-line options, so
they need no config file at all.

---

## Usage

### Rename gene models

```bash
perl scripts/rename_genes.pl -i CL_pasa.gff3 -o CL.renamed.gff3 -p Cla
```

Produces `Cla5AG00010`, `Cla5AG00020`, … numbered along each chromosome in steps
of ten so new models can be inserted later, with isoforms `.1`, `.2` ordered by
CDS length (`.1` is the representative model). The old identifier is kept as
`prev_id=` and an old→new table is written to `<out>.idmap.tsv`.

Options: `--step`, `--digits`, `--prefix`, `--source`, `--keep-name`, `-m`.

### Clean proteins and synchronise the three annotation files

```bash
python3 scripts/clean_pep.py -i CL.pep -o CL_clean.pep -d CL.pep.discarded.tsv

python3 scripts/sync_annotation.py -g CL.renamed.gff3 -p CL_clean.pep \
        -c CL.cds -o CL.filtered
```

`sync_annotation.py` removes a gene only when **all** of its transcripts were
discarded — a gene with one bad isoform and two good ones is kept with the bad
isoform gone. It ends with a consistency check that must report
`orphan features : 0` and `GFF3 == protein set : True`.

### Filter repeat-overlapping genes and short lncRNAs

```bash
python3 scripts/filter_repeat_gene.py -g missed.filtered.gff3 \
        -k keep_gene.ids -o missed.nonoverlap.gff3

python3 scripts/lncRNA_filter.py -g CL.lncRNA.final.gff3 \
        -o CL_lncRNA_clean.gff3 --min-length 200 --report dropped.tsv
```

`lncRNA_filter.py` sums **exon** lengths, not the genomic span, which is the
length the ≥ 200 nt definition refers to.

### Run the pipelines

```bash
CONFIG=config/my_run.local.sh scripts/run_evm.sh
CONFIG=config/my_run.local.sh scripts/functional_annotation.sh
CONFIG=config/my_run.local.sh scripts/ncRNA_annotation.sh
```

---

## Input and output formats

All GFF3 output follows the GFF3 specification: attributes escaped only where
required (`;` `=` `&` `,`), one `ID` per feature, and every child carrying a
`Parent`. Validate with:

```bash
gt gff3 -sort -tidy -retainids CL.filtered.gff3 > /dev/null
```

Gene identifiers use the format `Cla<chromosome><haplotype>G<serial>`, e.g.
`Cla1AG00010`, so the haplotype and chromosome of origin are readable from the
identifier itself.

---

## Development

```bash
bash test/run_tests.sh                                    # toy-data test suite
shellcheck --severity=warning --exclude=SC1090,SC1091 scripts/*.sh
flake8 --max-line-length=100 scripts/*.py
perl -c scripts/rename_genes.pl
```

Conventions kept throughout:

- Shell scripts run under `set -euo pipefail`; every variable expansion is quoted.
- Inputs are resolved and checked for existence **before** any long-running step.
- Errors go to stderr, results to stdout, so output can be piped safely.
- No credentials, absolute personal paths or data files in the repository —
  `.gitignore` blocks FASTA, BAM, database and `config/*.local.sh` files.
- Scripts never delete or overwrite their input; output paths are explicit
  arguments.

Contributions: open an issue or pull request. CI must pass before merge.

---

## Citation

If you use these scripts, please cite the Data Descriptor (see `CITATION.cff`):

> [Authors]. Haplotype-resolved T2T genome assembly and annotation of seedless
> wampee [*Clausena lansium* (Lour.) Skeels cv. 'Yunan']. *Scientific Data* (2026).
> doi:[DOI]

Please also cite the third-party tools the pipelines call — EVidenceModeler,
PASA, GeMoMa, ANNEVO, TransDecoder, StringTie, miniprot, eggNOG-mapper,
InterProScan, DIAMOND, tRNAscan-SE, Infernal/Rfam and CPC2 — as listed in the
Methods of the Data Descriptor.

## License

MIT — see [LICENSE](LICENSE).
