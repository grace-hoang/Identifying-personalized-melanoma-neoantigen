# Melanoma neoantigen → candidate mRNA construct pipeline (Bash + Docker)

**Research use only.** Computational predictions; nothing here is validated or intended for clinical use.

Tumor WES + matched normal WES + tumor RNA-seq FASTQ → somatic variants → HLA class I typing →
pVACseq (NetMHCpan / MHCflurry) binding predictions → ranked neoantigens → in-silico multi-epitope mRNA design → HTML report.

## 1. One-time setup

1. **WSL2 resources** (Windows, `%UserProfile%\.wslconfig`), then `wsl --shutdown`:
   ```
   [wsl2]
   memory=27GB
   processors=16
   swap=32GB          # safety net for STAR index build
   ```
2. Docker must work inside WSL (`docker run hello-world`). Keep this repo under `/home/...`, **not** `/mnt/c`.
3. Edit `config/samplesheet.csv` (already filled with Patient01) and, if you know it, set `TARGET_BED` in `config/config.sh`
   (exome kit BED, hg38, `chr`-prefixed). If empty, targets are derived from normal-sample coverage.
4. `./run.sh setup` — builds 8 local bioconda images, pulls GATK/VEP/pVACtools/bam-readcount images, downloads and indexes
   references (hg38, gnomAD/PoN, GENCODE 46, STAR, Salmon, VEP 112 cache). Resumable; expect 2–4 h and ~100 GB
   (pVACtools image ≈ 25 GB on disk, VEP cache ≈ 20 GB, STAR index ≈ 18 GB).

## 2. Running

```
./run.sh                               # preflight + everything, resumes automatically
./run.sh qc | align | variants | rna | hla | neoantigen | report
./run.sh variants rna                  # several stages
./run.sh --step 10,11                  # numbered steps;  --from 03 --to 06
./run.sh --patient Patient01 --force   # recompute (ignore markers)
./run.sh --dry-run align               # print Docker commands only
./run.sh status                        # finished units + disk usage
./run.sh cleanup --level intermediate  # tmp | intermediate | docker | all
```

| Step | Script | Tool(s) |
|---|---|---|
| 00 | `00_preflight.sh` (+ `00_build_images.sh`, `00_setup_reference.sh`) | checks, image build, references |
| 01 | `01_fastqc.sh` | FastQC, MultiQC |
| 02 | `02_fastp.sh` | fastp |
| 03 | `03_bwa_mem2.sh` | BWA-MEM (or BWA-MEM2) + samtools sort |
| 04 | `04_markdup.sh` | GATK MarkDuplicates, samtools, mosdepth |
| 05 | `05_mutect2.sh` | GATK Mutect2 (4 parallel shards), contamination, FilterMutectCalls |
| 06 | `06_annotation.sh` | bcftools norm, Ensembl VEP + pVACtools plugins |
| 07 | `07_star.sh` | STAR (sparse SA index) |
| 08 | `08_salmon.sh` | Salmon (decoy-aware) → transcript/gene TPM |
| 09 | `09_hla.sh` | OptiType on MHC reads from the normal exome |
| 10 | `10_netmhcpan.sh` | VAtools expression/readcount annotation, **pVACseq** |
| 11 | `11_ranking.sh` | composite scoring + selection (`scripts/py/rank_neoantigens.py`) |
| 12 | `12_mrna_design.sh` | polyepitope + codon optimisation + UTRs/poly(A) + RNAfold (`design_mrna.py`) |
| 13 | `13_report.sh` | MultiQC + per-patient HTML |

**Resume:** every unit writes `results/.markers/<unit>.done` *after* its output is complete (outputs are written as `*.partial`
then renamed). A unit is skipped only if marker and output files both exist. After a failure, just re-run the same command.
Intermediates (trimmed FASTQ, lane/merged BAMs, Mutect2 shards) are deleted once the next product exists; later stages detect
this and do not redo upstream work. `--keep-intermediate` disables deletion.

## 3. Hardware decisions

* **BWA-MEM, not BWA-MEM2, by default.** `bwa-mem2 index` needs ≈28 GB RAM per Gbp (≈87 GB for hg38), impossible on 28 GB.
  Alignments are equivalent. Set `ALIGNER="bwa-mem2"` only with a prebuilt index (`reference/genome/*.bwt.2bit.64` …).
* **STAR** uses `--genomeSAsparseD 2` (~18 GB in RAM). If index generation is OOM-killed, set `STAR_SA_SPARSE=3`.
* **Mutect2** runs 4 shards × 5.5 GB; **stages run strictly one at a time**, so peak RAM stays < 26 GB.
* Peak disk during alignment ≈ 60 GB per patient on top of references; `run.sh cleanup` frees scratch space.
  WSL's virtual disk does not shrink automatically after deleting files (`wsl --shutdown` then compact the VHDX if needed).

## 4. Things you must check

* **Image tags** live in `config/config.sh` (`IMG_*`). pVACtools is pinned to `7.1.5` (fallback `latest`). If a pull fails,
  change the tag; nothing else needs editing. Conda-built images are unpinned except the main tools.
* **NetMHCpan:** pVACseq uses a local IEDB install if the image has one (`/opt/iedb`), otherwise the IEDB **web API**
  (internet, slow, rate-limited). NetMHCpan itself is DTU-licensed; to use your own licensed install set
  `PVAC_LOCAL_IEDB_HOST_DIR`. Change algorithms with `PVAC_ALGORITHMS`.
* **HLA typing** via OptiType is class I only. If it fails, copy `config/hla_override.tsv.example` to `hla_override.tsv`.
* **Construct elements** in `config/design/` (signal peptide, MITD, UTRs) are literature-style examples: verify against
  primary sources and replace with your validated sequences before any synthesis planning.
* **Not covered:** junctional neoepitope screening (use pVACvector or NetMHCpan on junctions), class II, fusions/splice
  neoantigens, HLA LOH, mRNA chemistry (cap, N1-methylpseudouridine).
* Ranking matches pVACseq columns by regex; if a future pVACtools renames columns the script prints the columns it found.
  `tests/smoke_test.sh` exercises ranking → design → report offline on synthetic data.

## 5. Outputs (`results/`)

`qc/` FastQC/fastp/mosdepth · `alignment/<P>/*.md.bam` · `variants/<P>/*.somatic.pass.norm.vep.vcf.gz` ·
`rna/<P>/star|salmon|expression` · `hla/<P>/*.hla_alleles.txt` · `neoantigen/<P>/{pvacseq,*.ranked_epitopes.tsv,*.top_candidates.tsv,*.design_input.tsv}` ·
`mrna/<P>/*.construct.{protein,cds,mrna,template}.fa, *.construct.map.tsv, *.design_report.json` · `report/<P>/<P>.report.html`.
Logs: `logs/run_*.log` (whole run) and `logs/<step>.log`, plus per-tool logs.
