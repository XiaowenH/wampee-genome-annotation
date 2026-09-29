# ---------------------------------------------------------------------------
# annotation.config.sh - every path and parameter used by the shell pipelines.
#
# Copy this file, edit the copy, and point the scripts at it:
#
#     cp config/annotation.config.sh config/my_run.config.sh
#     $EDITOR config/my_run.config.sh
#     CONFIG=config/my_run.config.sh scripts/run_evm.sh
#
# Anything already set in the environment wins over the value here, so a single
# value can be overridden for one run without editing the file:
#
#     CPU=8 scripts/run_evm.sh
#
# Nothing in this file should be committed with real absolute paths from your
# machine; keep those in an untracked copy (config/*.local.sh is git-ignored).
# ---------------------------------------------------------------------------

# --- sample ---------------------------------------------------------------
SAMPLE="${SAMPLE:-CL}"                     # prefix used for output files
SPECIES_PREFIX="${SPECIES_PREFIX:-Cla}"    # gene ID prefix, e.g. Cla05AG00160

# --- compute --------------------------------------------------------------
CPU="${CPU:-8}"                            # threads; raise on a cluster node

# --- reference ------------------------------------------------------------
GENOME="${GENOME:-data/CL_v1.fa}"          # soft-masked assembly, FASTA

# --- gene prediction inputs (run_evm.sh) ----------------------------------
GEMOMA_IN="${GEMOMA_IN:-data/GeMoMa/final_annotation.gff}"
ANNEVO_IN="${ANNEVO_IN:-data/ANNEVO/ANNEVO.gff}"
TD_IN="${TD_IN:-data/TransDecoder/stringtie.cDNA.fasta.transdecoder.genome.gff3}"
STRINGTIE_IN="${STRINGTIE_IN:-data/stringtie/stringtie.gtf}"
MINIPROT_IN="${MINIPROT_IN:-data/miniprot/swiss_miniprot_filtered.gff}"   # optional

# EVidenceModeler installation directory (must contain ./EVidenceModeler)
EVM="${EVM:-/opt/EVidenceModeler-v2.1.0}"

# EVM evidence weights
W_GEMOMA="${W_GEMOMA:-8}"
W_TRANSDECODER="${W_TRANSDECODER:-7}"
W_ANNEVO="${W_ANNEVO:-4}"
W_STRINGTIE="${W_STRINGTIE:-10}"
W_MINIPROT="${W_MINIPROT:-5}"
EVM_SEGMENT_SIZE="${EVM_SEGMENT_SIZE:-1000000}"
EVM_OVERLAP_SIZE="${EVM_OVERLAP_SIZE:-200000}"

# --- functional annotation (functional_annotation.sh) ---------------------
PEP="${PEP:-results/${SAMPLE}.pep}"        # proteins from gffread on the renamed GFF3
GFF="${GFF:-results/${SAMPLE}.gff3}"
EGGNOG_DATA="${EGGNOG_DATA:-db/eggnog}"    # directory with eggnog.db + eggnog_proteins.dmnd
SPROT="${SPROT:-db/SwissProt.dmnd}"
NRDB="${NRDB:-db/nr_Magnoliophyta.dmnd}"
TAX_SCOPE="${TAX_SCOPE:-33090}"            # Viridiplantae
EVALUE="${EVALUE:-1e-5}"
FUNC_OUT="${FUNC_OUT:-results/func}"

# --- non-coding RNA (ncRNA_annotation.sh) ---------------------------------
RFAM_DIR="${RFAM_DIR:-db/rfam}"            # Rfam.cm + Rfam.clanin
NCRNA_OUT="${NCRNA_OUT:-results/ncRNA}"

# --- lncRNA ---------------------------------------------------------------
LNCRNA_MIN_LENGTH="${LNCRNA_MIN_LENGTH:-200}"
