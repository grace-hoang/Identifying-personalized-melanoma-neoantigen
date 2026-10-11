#!/usr/bin/env bash
# =============================================================================
#  config/config.sh  —  SINGLE SOURCE OF TRUTH for the melanoma neoantigen
#  pipeline. Every script sources this (via scripts/lib/common.sh).
#  RESEARCH USE ONLY. Not for clinical decision making.
# =============================================================================
# shellcheck disable=SC2034   # variables are consumed by the other scripts

PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"

# -----------------------------------------------------------------------------
# 1. PATHS
# -----------------------------------------------------------------------------
CONFIG_DIR="$PROJECT_ROOT/config"
SAMPLESHEET="${SAMPLESHEET:-$CONFIG_DIR/samplesheet.csv}"
REF_DIR="$PROJECT_ROOT/reference"
DATA_DIR="$PROJECT_ROOT/data"
RESULTS_DIR="$PROJECT_ROOT/results"
LOG_DIR="$PROJECT_ROOT/logs"
TMP_DIR="$DATA_DIR/tmp"                  # big scratch space (sort tmp, STAR tmp...)
MARK_DIR="$RESULTS_DIR/.markers"         # "stage finished" markers used for resume
DOCKER_CTX_DIR="$PROJECT_ROOT/docker"

# Extra host directories to bind-mount into every container (same path inside).
# FASTQ directories from the samplesheet are mounted read-only automatically.
EXTRA_RW_MOUNTS=()
EXTRA_RO_MOUNTS=()

# -----------------------------------------------------------------------------
# 2. HARDWARE BUDGET  (tuned for: 16 cores, 28 GB RAM, WSL2, ~200 GB disk)
#    Stages run ONE AT A TIME, so each stage may use most of the machine, but
#    the sum of concurrently running containers inside a stage must stay below
#    ~24 GB so WSL2 and Windows keep some headroom.
# -----------------------------------------------------------------------------
THREADS="${THREADS:-16}"

FASTQC_THREADS=6;   FASTQC_MEM=8g
FASTP_THREADS=12;   FASTP_MEM=4g                 # fastp scales poorly > 12
ALIGN_THREADS=12;   SORT_THREADS=4;  SORT_MEM_PER_THREAD=2G; ALIGN_DOCKER_MEM=20g
MARKDUP_XMX=10g;    MARKDUP_DOCKER_MEM=14g
MOSDEPTH_THREADS=4
MUTECT_SCATTER=4                                 # parallel Mutect2 shards
MUTECT_XMX=4g;      MUTECT_DOCKER_MEM=5500m;    PAIRHMM_THREADS=3   # 4 x 5.5 GB = 22 GB
GATK_POST_XMX=6g;   GATK_POST_DOCKER_MEM=8g
STAR_THREADS=16;    STAR_DOCKER_MEM=25g
STAR_SORT_RAM=6000000000                         # bytes for BAM sorting
STAR_SA_SPARSE=2                                 # 2 => ~18 GB index in RAM (1 => ~30 GB: too big)
STAR_TWOPASS="None"                              # None | Basic (Basic ~2x slower)
SALMON_THREADS=16;  SALMON_DOCKER_MEM=14g
VEP_FORK=8;         VEP_DOCKER_MEM=8g
PVAC_THREADS=8;     PVAC_DOCKER_MEM=16g
PY_DOCKER_MEM=6g

# Disk guard-rails (GB free required before heavy steps)
MIN_FREE_GB_START=60
KEEP_INTERMEDIATE="${KEEP_INTERMEDIATE:-0}"      # 1 = never delete trimmed FASTQ / merged BAMs

# -----------------------------------------------------------------------------
# 3. REFERENCE FILES  (downloaded/built by:  ./run.sh setup)
# -----------------------------------------------------------------------------
GENCODE_RELEASE=46          # GENCODE 46 == Ensembl 112. Keep these two in sync.
ENSEMBL_RELEASE=112

REF_FASTA="$REF_DIR/genome/Homo_sapiens_assembly38.fasta"     # Broad hg38 (chr-prefixed)
REF_FAI="${REF_FASTA}.fai"
REF_DICT="$REF_DIR/genome/Homo_sapiens_assembly38.dict"
PRIMARY_CONTIGS_LIST="$REF_DIR/genome/primary_contigs.list"
PRIMARY_CONTIG_REGEX='^chr([0-9]+|X|Y)$'                       # calling restricted to these
PRIMARY_FASTA="$REF_DIR/genome/GRCh38.primary_assembly.genome.fa"   # STAR + salmon decoys

GERMLINE_RESOURCE="$REF_DIR/gatk/af-only-gnomad.hg38.vcf.gz"
PON="$REF_DIR/gatk/1000g_pon.hg38.vcf.gz"
CONTAM_VARIANTS="$REF_DIR/gatk/small_exac_common_3.hg38.vcf.gz"

GENCODE_GTF="$REF_DIR/annotation/gencode.v${GENCODE_RELEASE}.annotation.gtf"
GENCODE_TX_FASTA="$REF_DIR/annotation/gencode.v${GENCODE_RELEASE}.transcripts.fa.gz"
TX2GENE="$REF_DIR/annotation/tx2gene.v${GENCODE_RELEASE}.tsv"

STAR_INDEX_DIR="$REF_DIR/index/star_gencode${GENCODE_RELEASE}"
SALMON_INDEX_DIR="$REF_DIR/index/salmon_gencode${GENCODE_RELEASE}"
STAR_SJDB_OVERHANG=100

VEP_CACHE_DIR="$REF_DIR/vep/cache"
VEP_PLUGIN_DIR="$REF_DIR/vep/plugins"

# Aligner. "bwa" (classic bwa-mem, ~6 GB RAM, identical alignments) is the default
# because `bwa-mem2 index` needs ~28 GB RAM per Gbp of reference (~87 GB for hg38)
# which this workstation cannot provide. Set "bwa-mem2" only if you already have a
# prebuilt bwa-mem2 index next to $REF_FASTA (e.g. built on another machine).
ALIGNER="bwa"

# Exome capture kit BED (hg38, chr-prefixed). Leave EMPTY to derive callable target
# regions automatically from normal-sample coverage (>= MIN_COVERED_DEPTH).
TARGET_BED=""
MIN_COVERED_DEPTH=8
TARGET_PAD=100

# Download sources
URL_BROAD_HG38="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0"
URL_GATK_SOMATIC="https://storage.googleapis.com/gatk-best-practices/somatic-hg38"
URL_GENCODE="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_RELEASE}"
URL_VEP_CACHE="https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/variation/indexed_vep_cache/homo_sapiens_vep_${ENSEMBL_RELEASE}_GRCh38.tar.gz"

# -----------------------------------------------------------------------------
# 4. DOCKER IMAGES
#    Images tagged melpipe/*:local are built from docker/Dockerfile.conda using the
#    bioconda package lists below (./run.sh setup builds them). The others are pulled.
#    If a pull fails because a tag moved, edit the tag here — nothing else.
# -----------------------------------------------------------------------------
MICROMAMBA_BASE="mambaorg/micromamba:1.5.8"

IMG_FASTQC="melpipe/fastqc:local"
IMG_FASTP="melpipe/fastp:local"
IMG_ALIGN="melpipe/align:local"
IMG_STAR="melpipe/star:local"
IMG_SALMON="melpipe/salmon:local"
IMG_BCFTOOLS="melpipe/bcftools:local"
IMG_OPTITYPE="melpipe/optitype:local"
IMG_PY="melpipe/py:local"

IMG_GATK="broadinstitute/gatk:4.6.1.0"
IMG_VEP="ensemblorg/ensembl-vep:release_${ENSEMBL_RELEASE}.0"
IMG_PVAC="griffithlab/pvactools:7.1.5"           # fallback: griffithlab/pvactools:latest
IMG_BAMRC="mgibio/bam_readcount_helper-cwl:1.1.1"

declare -A LOCAL_IMAGES=(
  ["$IMG_FASTQC"]="fastqc=0.12.1 multiqc"
  ["$IMG_FASTP"]="fastp=0.23.4"
  ["$IMG_ALIGN"]="bwa bwa-mem2 samtools mosdepth bedtools"
  ["$IMG_STAR"]="star=2.7.11b samtools"
  ["$IMG_SALMON"]="salmon=1.10.3"
  ["$IMG_BCFTOOLS"]="bcftools htslib"
  ["$IMG_OPTITYPE"]="optitype glpk razers3"
  ["$IMG_PY"]="python=3.11 pandas numpy biopython pysam viennarna"
)

# -----------------------------------------------------------------------------
# 5. QC / TRIMMING
# -----------------------------------------------------------------------------
FASTP_LENGTH_REQUIRED=36
FASTP_QUALITY=20
FASTP_EXTRA=(--detect_adapter_for_pe --correction --cut_tail --cut_tail_mean_quality 20)
CLEAN_TRIMMED_FASTQ=1        # delete trimmed FASTQ once downstream product exists
OPTICAL_DUP_DISTANCE=2500    # 100 for non-patterned flowcells, 2500 for patterned (NovaSeq/HiSeq4000)

# -----------------------------------------------------------------------------
# 6. HLA TYPING
# -----------------------------------------------------------------------------
# Optional manual override (patient<TAB>HLA-A*02:01,HLA-A*24:02,...) — if a patient is
# listed here, OptiType is skipped for them.
HLA_OVERRIDE_FILE="$CONFIG_DIR/hla_override.tsv"
HLA_OPTITYPE_THREADS=8
HLA_MHC_REGION="chr6:28510120-33480577"      # GRCh38 MHC extended region

# -----------------------------------------------------------------------------
# 7. NEOANTIGEN PREDICTION (pVACseq)
# -----------------------------------------------------------------------------
VEP_PICK=1                                   # 1 = --pick (one transcript per variant; far faster)
PVAC_CLASS_I_LENGTHS="8,9,10,11"
PVAC_ALGORITHMS="NetMHCpan NetMHCpanEL MHCflurry"   # space separated; see `pvacseq run -h`
# "auto": use /opt/iedb inside the pVACtools image if it exists, else IEDB REST API (needs internet)
# "none": force REST API.  Any other value: path *inside the container* of an IEDB install.
PVAC_IEDB_DIR="auto"
# Host directory with a licensed local IEDB / NetMHCpan install to mount at /iedb (optional).
PVAC_LOCAL_IEDB_HOST_DIR=""
PVAC_BINDING_THRESHOLD=500
PVAC_PERCENTILE_THRESHOLD=2.0
PVAC_MIN_FOLD_CHANGE=0
PVAC_NORMAL_VAF_MAX=0.02
PVAC_TUMOR_DNA_VAF_MIN=0.05
PVAC_EXPN_VAL=1.0
TUMOR_PURITY=""                              # 0-1, optional (improves clonality estimate)
PVAC_EXTRA_ARGS=()

# -----------------------------------------------------------------------------
# 8. RANKING (scripts/py/rank_neoantigens.py)
# -----------------------------------------------------------------------------
RANK_IC50_MAX=500
RANK_PERCENTILE_MAX=2.0
RANK_MIN_FOLD_CHANGE=1.0
RANK_MIN_TUMOR_VAF=0.05
RANK_MIN_TUMOR_DP=10
RANK_MAX_NORMAL_VAF=0.02
RANK_MIN_EXPR_TPM=1.0
RANK_TOP_N=30                                # rows written to the "top candidates" table

# -----------------------------------------------------------------------------
# 9. mRNA CONSTRUCT DESIGN (in silico only)
# -----------------------------------------------------------------------------
DESIGN_DIR="$CONFIG_DIR/design"              # UTR / signal / MITD / PADRE sequences live here
DESIGN_N_EPITOPES=20                         # max neoantigen units in the construct
DESIGN_FLANK=12                              # residues each side of a missense site (0 = bare epitope)
DESIGN_LINKER_SHORT="AAY"                    # between bare 8-11-mers (proteasomal cleavage)
DESIGN_LINKER_LONG="GGSGGGGSGG"              # between long (>=20 aa) units
DESIGN_INCLUDE_SIGNAL=1                      # N-terminal secretory signal peptide
DESIGN_INCLUDE_MITD=1                        # C-terminal MHC-I trafficking domain
DESIGN_INCLUDE_PADRE=0                       # universal CD4 helper epitope
DESIGN_POLYA_LEN=100
DESIGN_SEED=42
DESIGN_ANNEAL_ITERS=6000
DESIGN_GC_MIN=0.48
DESIGN_GC_MAX=0.62
DESIGN_AVOID_SITES="GGTCTC GAGACC CGTCTC GAGACG GAATTC AAGCTT"   # BsaI, BsmBI, EcoRI, HindIII
DESIGN_T7_PROMOTER="TAATACGACTCACTATAG"
