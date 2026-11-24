# K.-brevis-viromes
Re-analysis of Illumina sequences published by Booker et al. 2026. 


## Experimental background: 

Harmful algal blooms, or HABs pose a growing threat to coastal 
community well-being and economic stability globally. In Florida, the dinoflagellate *Karenia 
brevis* is exceptionally problematic. This species produces the potent neurotoxins called brevetoxins 
that can cause immediate damage and also chronic neurological issues. Despite their significant impact
on Florida residents, little is known about what drives the dynamics of *Karenia* blooms. Current theory 
predicts that blooms are influenced by controls on physical, chemical, and biological drivers, including virus-mediated 
mortality. More recent work from other systems has shown that viral infection can sometimes prolong blooms (either 
by providing hosts with products from viral auxiliary metabolic genes that support host survival, infecting competitor species, 
or providing hosts with superinfection fitness), terminate blooms by causing widespread and rapid lysis, or increase toxin 
release into the environment. Others have shown that viruses have the potential to be used as a biocontrol strategy and can even 
serve as early warning signals for growing algal blooms before they are detected by conventional methods (using virus-specific peptides).  Therefore, 
expanding our current understanding of the viral community associated with harmful algae is essential to comprehensively 
characterizing bloom drivers and developing more effective early warning systems. 

In this study, the authors collected 2-5L of seawater from a sampling location off the coast
of Cayo Costa, Fl, USA in December-September 2021 and again in March 2022. Upon collection, 
seawater was first filtered through a 100um mesh filter, followed by 10-3-0.22um polycarbonate 
filters to size-fractionate the sample. After flash-freezing and storage at -80˚C, only the 0.22um 
filters were used for subsequent metagenomic sequencing, as these would capture viral particles. It is important 
to note that only very large viruses will be retained by a 0.22uM filter. This is why their analysis focused 
on identifying giant viruses using a set of tools specifically designed for this purpose. 

In tandem, the authors inoculated a culture of K. brevis bloom water that was either filtered 
to only contain the viral fraction (VF, <0.22um filtrate) or viral and bacterial fraction (BVF, <1um filtrate). 
After 8 or 15 days, 50mL of the culture was passed through 0.22um filters and total DNA was extracted from 
these filters. The DNA samples were shipped to the University of Illinois at Urbana-Champaign Roy J. Carver 
Biotechnology Center DNA Services Lab. Shotgun DNA libraries were prepared with the xGenTM ssDNA-Low Input 
DNA library prep kit (IDT, IA, USA) with a target size of 150 bp paired-end reads. The libraries were pooled, 
quantitated by qPCR, and sequenced on one NovaSeq 6000 SP lane (Illumina, CA, USA) for 151 cycles from both fragment 
ends. 

Giant virus "species" captured from the incubation experiment were then mapped against the reads from the environmental 
samples to determine their prevalence through time. 

## Purpose: 
In this project, I sought to re-analyze the metagenomes published by Booker et al. 2026 using 
a slightly different pipeline. Most importantly, I was interested in knowing whether I would 
be able to identify more viruses than they did by including the results of Vasquez et al. 2026. In
this paper, the authors scanned public databases and used a novel gene prediction method that allowed 
for alternative genetic codes. This significantly expanded giant virus databases (4 fold or 6 fold relative 
to two previous studies). However, this paper was published after Booker et al. 2026, and so they were not able to 
include this data. 


----------------------------------------------------------------------------------------

# My re-analysis: how to reproduce it

Everything below this line was rewritten on 2026-11-12 by Claude Code (an AI coding agent,
prompt 8 in `PROMPT_LOG.md`) from what is in this repository, in
`~/scr10/data/kbrevis/`, in the conda environments, and in the Slurm logs at that time. The
sections above (Experimental background, Purpose, Pipeline followed by the manuscript) are
my own notes and were not changed. The previous version of this part of the README (my
hand-written steps 1–2 and Claude's section 3) is still in git history
(`git show 6f9af91:README.md`); everything in it is carried over below or listed under
"Earlier attempts".

Every script is reproduced in full below, copied from `scripts/` on 2026-11-12. If a script
and this README ever disagree, the script in `scripts/` is what actually runs.

## Pipeline overview

| step | what | script(s) | main tool (version) | conda env | output (in `~/scr10/data/kbrevis/`) | status (2026-11-12) |
|---|---|---|---|---|---|---|
| 00 | make folders | `00_directory_setup.sh` | bash | – | folder tree | done (originally by hand; script fixed 2026-11-12) |
| 01 | download reads from NCBI SRA | `01_download_files.sh/.slurm` | fasterq-dump (sra-tools 3.4.1) | `sra-toolkit-env` | `raw/` | done |
| 02 | trim adapters / low-quality bases | `02_fastp.sh/.slurm` | fastp 1.3.7 | `fastp-env` | `trim/` | done (job 577472) |
| 03 | download contaminant references | `03_download_refs.sh/.slurm`, `03b_download_phix.sh/.slurm` | NCBI datasets 18.39.0, seqtk 1.5 | `datasets-env`, `seqtk-env` | `reference/` | done |
| 04 | remove reads from phiX, human, mouse, cat, dog, *K. brevis* | `04_clean_refs.sh/.slurm` | bowtie2 2.5.5, samtools 1.24 | `host-screen-env` | `clean/` | done (job 607624) |
| 05a | viral assembly | `05_metaspades.sh/.slurm` | metaviralSPAdes (SPAdes 4.3.0) | `metaspades-env` | `contigs/` | done, all 20 samples |
| 05b | whole-metagenome assembly | `05_metaspades_nov.sh/.slurm` | metaSPAdes (SPAdes 4.3.0) | `metaspades-env` | `contigs_no_v/` | **running** on gulf `gz` (jobs 92235–92254) |
| 06 | compare + merge the two assemblies | `06_merge_assemblies.sh/.slurm`, `merge_assemblies.py` | BLAST+ 2.16.0, `anicalc` (a helper bundled with CheckV; **not** CheckV itself and **not** virus identification), seqkit 2.14.0 | `seqkit-env`, `checkv-env` | `contigs_merged/` | written + tested, waiting for step 05b |
| – | bin contigs into genomes (one-sample test; scripts removed, real binning is inside step 08) | none (removed 2026-11-12) | BBMap 40.05, samtools, MetaBAT 2.12.1 | `binning-env` | `binning_test/` | test done on 1 sample (job 611846) |
| – | CheckV database (scripts removed; the database is now made by step 07) | none (removed 2026-11-12) | CheckV 1.1.1 | `checkv-env` | `~/data10/databases/checkv_db/checkv-db-v1.5` | database done (one copy only) |
| 07 | one-time setup: environments, GVClass, databases | `07_gvmag_setup.slurm` | pixi, GVClass 1.1.1, CheckV 1.1.1, geNomad 1.7.4 | creates `gvmag-bin-env`, `genomad-env`, `aniclustermap-env` | `~/data10/databases/`, `~/data10/programs/` | done (jobs 73832, 73834, 73835) |
| 08 | giant-virus genomes (GVMAGs) following Booker et al. | `08_gvmags.sh/.slurm`, `select_gvmags.py`, `summarize_gvmags.py` | MetaBAT 2.15, GVClass 1.1.1, CheckV 1.1.1 (db v1.5), geNomad 1.7.4 (predicts genes with pyrodigal-gv) | `gvmag-bin-env`, `checkv-env`, `genomad-env` | `gvmags/` | written + tested on mock data; waiting for step 06 |
| 09 | group GVMAGs into viral OTUs ("viral species") | `09_gvmag_votus.slurm`, `votu_cluster.py`, `votu_nearest_ref.py` | fastANI 1.33, ANIclustermap 1.1.0 | `aniclustermap-env`, `checkv-env` | `votus/` | written + tested on mock data; waiting for step 08 |
| 10 | concatenated marker-gene alignment for a phylogeny | `10_gvmag_phylogeny.slurm`, `select_refs.py` | ncldv_markersearch, HMMER 3.4, Clustal Omega 1.2.4, trimAl 1.5.1 | `phylo-env` | `phylogeny/` | written + partly tested; waiting for step 09 |
| 11 | maximum-likelihood tree from that alignment | `11_gvmag_tree.slurm` | IQ-TREE 2.3.6, RAxML 8.2.12 | `tree-env` | `tree/` | written + tested end to end on small alignments; waiting for step 10 |

## A. Where everything lives

#### A1. Directory structure

The git repository holds only code, notes and logs. All sequence data are in scratch
(`~/scr10`), because home has a 50 GB quota and jobs must not read and write heavily there.
Long-term databases go to `~/data10`, which is backed up. **Scratch files that are not
touched for 90 days are deleted**, so copy anything important to `~/data10`.

```
~/K_brevis_viromes/                 git repo (this README)
├── README.md                       this lab notebook
├── PROMPT_LOG.md                   every prompt given to AI agents in this project
├── metadata/SraRunTable.csv        SRA metadata for the 20 runs (from SRA Run Selector)
├── envs/                           conda environment files (see section B)
├── host-screen-env.yml             re-exported by 04_clean_refs.slurm each time it runs
├── scripts/                        all scripts, numbered by pipeline step
├── logs/                           Slurm .out/.err files, <step>_<jobid>.out/.err
└── output/                         empty (early runs wrote here; that output was moved/deleted)

~/scr10/data/kbrevis/               all sequence data (scratch: purged after 90 days unused)
├── SraRunTable.csv                 copy of metadata/SraRunTable.csv (read by step 01)
├── raw/                            SRRxxx_1.fastq, SRRxxx_2.fastq            (40 files, 2.0 TB)
├── trim/                           SRRxxx_{1,2}_TRIMMED.fastq                (40 files, 1.6 TB)
├── reference/                      contaminant genomes + bowtie2 index host_index.*.bt2l
├── clean/                          SRRxxx_{1,2}_CLEAN.fastq + SRRxxx_CLEAN.log (60 files, 1.1 TB)
├── contigs/                        metaviralSPAdes: SRRxxx_contigs.fasta + SRRxxx/ run folder (90 GB)
├── contigs_no_v/                   metaSPAdes: same layout (being filled by step 05b)
├── contigs_merged/                 step 06 output: merged contigs per sample (6.6 GB)
├── gvmags/                         step 08 output: bins, GVMAGs, CheckV/geNomad/GVClass (115 GB)
├── votus/                          step 09 output: vOTUs, ANI, heatmap (16 MB)
├── phylogeny/                      step 10 output: the trimmed alignment + manifest (629 MB)
└── tree/                           step 11 output: IQ-TREE run (in progress)

(`binning_test/`, `host_filtered/` and four smoke-test directories were deleted on
2026-11-17; see "Step 13. Cleaning up scratch".)

~/data10/databases/checkv_db/checkv-db-v1.5     CheckV database v1.5 (the only copy; used by every script)
~/data10/databases/genomad_db                   geNomad database (step 07)
~/data10/databases/gvmagsV2/GVDB_stats.tsv      GVMAGs V2 metadata, 18,725 giant virus genomes (step 07)
~/data10/programs/gvclass-v1.1.1                GVClass + its database, installed with pixi (step 07)
~/data10/programs/ncldv_markersearch            giant virus marker-gene search + concatenation (step 07; step 10 uses it)
(no entry: fastp now comes from the fastp-env conda environment, not a manual binary)
```

#### A2. Running jobs from dev

All work is submitted from `dev.sciclone.wm.edu` (the shared front-end; no computation
there). Each cluster has its own Slurm module, loaded in a subshell so they don't clash:

```
cd ~/K_brevis_viromes        # the .slurm scripts use paths relative to the repo root
( module purge; module load slurm/bora; sbatch scripts/02_fastp.slurm )
```

Most `.slurm` scripts call `bash ./scripts/...`, so **submit them from the repo root**.
`sbatch` copies your current shell environment into the job; this matters for fastp
(step 02), which is found through the `PATH` set in `~/.bashrc`.

## B. Software and conda environments

All conda environments live in `~/.conda/envs/` and are loaded inside each job with:

```
module load miniforge3
eval "$(conda shell.bash hook)"
conda activate <env-name>
```

| env | key tools (version installed) | used in step | file |
|---|---|---|---|
| `sra-toolkit-env` | sra-tools 3.4.1 (`fasterq-dump`) | 01, 03 | `envs/sra-toolkit-env.yml` |
| `fastp-env` | fastp 1.3.7 | 02 | `envs/fastp-env.yml` |
| `datasets-env` | ncbi-datasets-cli 18.39.0 | 03 | `envs/datasets-env.yml` |
| `seqtk-env` | seqtk 1.5 | 03 | `envs/seqtk-env.yml` |
| `host-screen-env` | bowtie2 2.5.5, samtools 1.24 | 04 | `envs/host-screen-env.yml` |
| `metaspades-env` | SPAdes 4.3.0 (`metaspades.py`, `metaviralspades.py`) | 05a, 05b | `envs/metaspades-env.yml` |
| `seqkit-env` | seqkit 2.14.0 | 06 | `envs/seqkit-env.yml` |
| `checkv-env` | CheckV 1.1.1, BLAST+ 2.16.0, DIAMOND 2.1.8, HMMER 3.4, Python 3.12 | 06, 08, 09 (CheckV and Python) | `envs/checkv-env.yml` |
| `binning-env` | BBMap 40.05, samtools 1.25, MetaBAT 2.12.1 | – (removed binning-test step) | `envs/binning-env.yml` |
| `bbmpa-env` (sic) | BBMap 40.05, samtools 1.24 | earlier attempt only | `envs/bbmpa-env.yml` |
| `readcorr-env` | BBMap 40.05, bfc r181, seqtk 1.5 | earlier attempt only | `envs/readcorr-env.yml` |
| `gvmag-bin-env` | BBMap 35.85, samtools 1.15.1, MetaBAT 2.15, seqkit | 08 (created by step 07) | `envs/gvmag-bin-env.yml` |
| `genomad-env` | geNomad 1.7.4, pyrodigal 3.5.1, pyrodigal-gv (used via geNomad) | 08 (created by step 07) | `envs/genomad-env.yml` |
| `aniclustermap-env` | ANIclustermap 1.1.0, fastANI 1.33, matplotlib <3.9 | 09 | `envs/aniclustermap-env.yml` |
| `phylo-env` | HMMER 3.4, Clustal Omega 1.2.4, MUSCLE 5.3, trimAl 1.5.1, pyrodigal-gv, pandas/numpy/natsort/biopython | 10 | `envs/phylo-env.yml` |
| `tree-env` | IQ-TREE 2.3.6, RAxML 8.2.12, Python 3.12 | 11 | `envs/tree-env.yml` | `envs/aniclustermap-env.yml` |

How the environment files were made:
- The three `gvmag-bin`/`genomad`/`aniclustermap` files were **written by hand** (by Claude, prompt
  6) and pin the tool versions reported in Booker et al. 2026's Supplemental Methods.
  `07_gvmag_setup.slurm` built them (job 73832). A fourth, `gvmag-checkv-env` (CheckV 0.9.0,
  the paper's version), was built too but removed on 2026-11-12 (prompt 17): I chose to use
  the newer CheckV 1.1.1 in `checkv-env` everywhere instead.
- `envs/metaspades-env.yml` was exported on 2026-11-10.
- All other files in `envs/` were exported from the existing environments on 2026-11-12
  (prompt 8), with build strings and the machine-specific `prefix:` line removed:

```
cd ~/K_brevis_viromes
module load miniforge3
for e in sra-toolkit-env datasets-env seqtk-env host-screen-env binning-env seqkit-env checkv-env bbmpa-env readcorr-env; do
    nice -n 19 conda env export -n $e --no-builds | grep -v '^prefix:' > envs/$e.yml
done
```

The environments were originally created by hand with commands like these (from my earlier
notes):

```
module load miniforge3
eval "$(conda shell.bash hook)"
conda create -y -n sra-toolkit-env -c bioconda sra-tools
conda create -y -n host-screen-env -c conda-forge -c bioconda bowtie2 samtools
conda create -y -n metaspades-env  -c conda-forge -c bioconda spades
conda create -y -n metabat-env     -c conda-forge -c bioconda metabat2
conda create -y -n checkv-env      -c conda-forge -c bioconda checkv
# datasets-env (ncbi-datasets-cli), seqtk-env (seqtk), seqkit-env (seqkit),
# binning-env (bbmap samtools metabat2) were made the same way.
```

**To rebuild an environment on another account**, run this inside a batch job or an
`salloc` session, not on dev, because solving an environment uses several cores:

```
module load miniforge3
conda env create -f envs/host-screen-env.yml
```

`--no-builds` exports pin versions but not exact builds, so a rebuild months later may pick
slightly different builds of the same versions. Some pinned versions may also disappear from
the channels. If a rebuild fails, relax the pins on minor libraries (e.g. `libgcc`,
`openssl`) and keep the pins on the bioinformatics tools.

**fastp moved into conda on 2026-11-19** (prompt 51). It used to be the prebuilt binary
`~/programs/fastp` (1.3.7, from the fastp GitHub releases), put on `PATH` by `~/.bashrc`.
That was the one genuinely unreproducible dependency in the pipeline: `sbatch` does not
source `~/.bashrc`, so step 02 worked only because a job inherits the submitting shell's
environment - it would have failed on any other account, and depended on the state of the
shell it was launched from here.

It is now `envs/fastp-env.yml`, **pinned to fastp 1.3.7**, the identical version, which
bioconda carries. `02_fastp.slurm` activates it, and `00_create_envs.slurm` builds it.
Both binaries were checked side by side and report `fastp 1.3.7`, and the invocation in
`02_fastp.sh` is unchanged, so re-running step 02 produces the same trimmed reads. The pin
matters: trimming parameters and fastp version together determine `trim/`, which every
later step depends on, so this should not be floated without re-running step 02 onwards.
The old `~/programs/fastp` binary can stay where it is; nothing references it any more.

Other environments exist in `~/.conda/envs` but are not used by any current script:
`fastqc-env`, `flye-env`, `megahit-env`, `metabat-env` (MetaBAT 2.18), `samtools`, `sra`.

## C. Samples

All 20 runs in BioProject PRJNA1077797 (Illumina NovaSeq 6000, paired-end 2 × 150 bp,
DNA from 0.22 µm filters). The table comes from `metadata/SraRunTable.csv`. (My earlier
notes said 22 runs; the run table has 20.)

| run | library name | collected | what it is | Gbases | metaviralSPAdes contigs |
|---|---|---|---|---|---|
| SRR28004153 | 1220_metaG | 2020-12 | EH25 field time series, Gulf of Mexico | 12.1 | 574 |
| SRR28004152 | 0121_metaG | 2021-01 | EH25 field time series | 16.4 | 231 |
| SRR28004151 | 0221_metaG | 2021-02 | EH25 field time series | 13.0 | 415 |
| SRR28004150 | 0321_metaG | 2021-03 | EH25 field time series | 11.9 | 214 |
| SRR28004149 | 0421_metaG | 2021-04 | EH25 field time series | 8.6 | 207 |
| SRR28004148 | 0521_metaG | 2021-05 | EH25 field time series | 16.8 | 212 |
| SRR28004147 | 0621_metaG | 2021-06 | EH25 field time series | 11.4 | 131 |
| SRR28004146 | 0721_metaG | 2021-07 | EH25 field time series | 10.5 | 295 |
| SRR28004145 | 0821_metaG | 2021-08 | EH25 field time series | 14.6 | 197 |
| SRR34748664 | 0921_metaG | 2021-09-29 | EH25 field time series, coastal Florida | 12.7 | 103 |
| SRR34748663 | 0322_metaG | 2022-03-08 | EH25 field time series, coastal Florida | 9.2 | 132 |
| SRR28004142 | CH1_Control_metaG | 2021-05 | incubation 1, control | 156.5 | 364 |
| SRR28004144 | CH1_VF1_metaG | 2021-05 | incubation 1, viral fraction rep 1 | 181.5 | 333 |
| SRR28004143 | CH1_VF2_metaG | 2021-05 | incubation 1, viral fraction rep 2 | 156.0 | 828 |
| SRR28004140 | CH2_Control_T0_metaG | 2021-08 | incubation 2, control day 0 | 9.1 | 902 |
| SRR28004141 | CH2_Control_T15_metaG | 2021-08 | incubation 2, control day 15 | 15.9 | 30 |
| SRR28004137 | CH2_VF1_metaG | 2021-08 | incubation 2, viral fraction rep 1 | 16.5 | 161 |
| SRR28004154 | CH2_VF2_metaG | 2021-08 | incubation 2, viral fraction rep 2 | 16.4 | 66 |
| SRR28004139 | CH2_BF1_metaG | 2021-08 | incubation 2, bacteria+virus fraction rep 1 | 15.5 | 87 |
| SRR28004138 | CH2_BF2_metaG | 2021-08 | incubation 2, bacteria+virus fraction rep 2 | 16.4 | 143 |

The three CH1 libraries are about 10 times deeper than the rest (156–182 Gbases each).
That is why they ran out of memory in step 05a and needed a bigger node. "What it is" is my
reading of the library names together with the paper; check it against the paper's sample
table before using it in figures.

----------------------------------------------------------------------------------------

# Step 00. Directory setup

#### SCRIPT: `scripts/00_directory_setup.sh`

```bash
#!/bin/bash
set -euo pipefail
# set up folders in data directory
mkdir -p "$HOME/scr10/data/kbrevis"
mkdir -p "$HOME/scr10/data/kbrevis/raw"
mkdir -p "$HOME/scr10/data/kbrevis/clean"

# set up folders in main directory
mkdir -p "$HOME/K_brevis_viromes/logs"
mkdir -p "$HOME/K_brevis_viromes/output"
mkdir -p "$HOME/K_brevis_viromes/scripts"
mkdir -p "$HOME/K_brevis_viromes/metadata"
mkdir -p "$HOME/K_brevis_viromes/envs"
```

**What it does.** Creates the scratch data folders (`~/scr10/data/kbrevis/`, `raw/`,
`clean/`) and the repo folders (`logs/`, `output/`, `scripts/`, `metadata/`, `envs/`).
`mkdir -p` makes any missing parent folders and does nothing if a folder already exists, so
the script is safe to re-run. It takes no input and runs in a second, so it can be run
directly on dev:

```
cd ~/K_brevis_viromes
bash scripts/00_directory_setup.sh
```

**Fix (2026-11-12, Claude Code, prompt 9 in `PROMPT_LOG.md`).** Every path was written
`~/$HOME/...`. Both `~` and `$HOME` stand for the home folder
(`/sciclone/home/iadanforth`), so together they pointed at the non-existent
`/sciclone/home/iadanforth/sciclone/home/iadanforth/...`, and the script would have created
a second copy of the whole tree there. The paths now use `$HOME` alone, quoted so that a
space in a path can't split it. (An earlier missing `#` in the first line had already been
fixed.) Before/after for each line:

```
# before
mkdir -p ~/$HOME/scr10/data/kbrevis
mkdir -p ~/$HOME/K_brevis_viromes/logs        # same pattern for raw, clean, output, scripts, metadata, envs
# after
mkdir -p "$HOME/scr10/data/kbrevis"
mkdir -p "$HOME/K_brevis_viromes/logs"
```

Tested with `bash -n scripts/00_directory_setup.sh` (no syntax errors) and then
`bash scripts/00_directory_setup.sh` (exit 0; all folders already existed, so nothing
changed). I also checked that the old version had never left a stray
`~/sciclone/home/...` or `~/~` folder: neither exists. The script was also made executable
(`chmod +x scripts/00_directory_setup.sh`). Originally the folders were made by hand on bora
with the equivalent commands.

The other data folders (`trim/`, `reference/`, `contigs/`, ...) are created by the scripts
that write to them.

# Step 01. Download the reads from NCBI

#### 01a. Get the run table

On the NCBI SRA Run Selector, search `PRJNA1077797`, select all runs, and download the
"Metadata" table as `SraRunTable.csv`. I did this on my laptop and copied it over:

```
scp ~/Downloads/SraRunTable.csv iadanforth@bora.sciclone.wm.edu:~/scr10/data/kbrevis
cp ~/scr10/data/kbrevis/SraRunTable.csv ~/K_brevis_viromes/metadata/
```

#### 01b. Check disk space first

The raw FASTQs are about 2 TB. My first attempt ran out of space, which is why everything is
written to scr10 and fasterq-dump gets its own temporary folder there.

```
df -h ~/scr10/data/kbrevis    # free space on the filesystem
du -sh ~/scr10/data/kbrevis   # space used by this folder
quota -s                      # home quota
```

#### SCRIPT: `scripts/01_download_files.sh`

```bash
#!/bin/bash
set -euo pipefail

# 2. apply for loop to SraRunTable.csv
RUN_TABLE="${1:-$HOME/scr10/data/kbrevis/SraRunTable.csv}"
RAW_DIR="$HOME/scr10/data/kbrevis/raw"
TEMP_DIR="$HOME/scr10/data/kbrevis/temp"

mkdir -p "$RAW_DIR" "$TEMP_DIR"

# Skip the CSV header, then read the first column (Run).
{
    IFS=, read -r _ || true
    while IFS=, read -r accession _ || [[ -n "${accession:-}" ]]; do
        accession="${accession//\"/}"
        accession="${accession%$'\r'}"

        [[ -z "$accession" ]] && continue

        echo "Downloading ${accession}..."
        fasterq-dump "$accession" \
            --split-files \
            --outdir "$RAW_DIR" \
            --temp "$TEMP_DIR" \
            --progress
    done
} < "$RUN_TABLE"
```

#### SCRIPT: `scripts/01_download_files.slurm`

```bash
#!/bin/bash
#SBATCH --job-name=kbrevis_download
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1-00:00:00
#SBATCH --mem=120G
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/kbrevis_download_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/kbrevis_download_%j.err

set -euo pipefail

# make sure that you're working in the SRA toolkit environment
# Activate sra-toolkit environment
module load miniforge3 
eval "$(conda shell.bash hook)"
conda activate sra-toolkit-env
conda env export --no-builds > sra-toolkit-env.yml

bash ./scripts/01_download_files.sh ~/scr10/data/kbrevis/SraRunTable.csv

conda deactivate 
conda deactivate

rm -r ~/scr10/data/kbrevis/temp
```

**What it does.** `fasterq-dump` (from NCBI's SRA Toolkit) downloads a sequencing run from
the Sequence Read Archive and converts it from NCBI's compressed SRA format into plain FASTQ
text files. The `.sh` script reads the first column (`Run`) of `SraRunTable.csv`, skipping
the header and stripping quotes and Windows line endings, and calls `fasterq-dump` once per
accession. The `.slurm` script loads `sra-toolkit-env`, runs the `.sh`, and deletes the
temporary folder at the end.

**Parameters.** `--split-files` writes read 1 and read 2 of each pair to separate files
(`SRRxxx_1.fastq`, `SRRxxx_2.fastq`), which every later step expects. `--temp` puts
fasterq-dump's large scratch files in `~/scr10/data/kbrevis/temp` instead of the small
default location. `--progress` prints progress to the `.err` log. The job asked for 4 cores
and 120 GB for 1 day on bora `batch`; fasterq-dump uses little memory, so 8–16 GB would have
been enough. To use all 4 cores it would also need `--threads 4` (the default is 6 threads
regardless of the allocation).

**Run.**
```
chmod +x scripts/01_download_files.sh scripts/01_download_files.slurm
( module purge; module load slurm/bora; sbatch scripts/01_download_files.slurm )
tail -n 100 logs/kbrevis_download_<JOBID>.err   # watch progress
```

**Result.** 40 FASTQ files (2 per run × 20 runs) in `~/scr10/data/kbrevis/raw/`, 2.0 TB.
The job's log files were not kept. The `conda env export` line in the `.slurm` writes
`sra-toolkit-env.yml` into the repo root; that copy was later deleted, and the current
export is `envs/sra-toolkit-env.yml`.

**Caveat.** Uncompressed FASTQ is very large. Later steps could read gzipped files
(`fasterq-dump` output piped to `pigz`) to save about 75 % of the space.

# Step 02. Trim adapters and low-quality bases (fastp)

#### SCRIPT: `scripts/02_fastp.sh`

```bash
#!/bin/bash
set -euo pipefail

RAW_DIR="$HOME/scr10/data/kbrevis/raw"
CLEAN_DIR="$HOME/scr10/data/kbrevis/trim"

command -v fastp >/dev/null 2>&1 || {
    echo "fastp is not available; activate the environment that contains it." >&2
    exit 1
}

mkdir -p "$CLEAN_DIR"

shopt -s nullglob
R1_FILES=("$RAW_DIR"/*_1.fastq)

if ((${#R1_FILES[@]} == 0)); then
    echo "No *_1.fastq files found in $RAW_DIR" >&2
    exit 1
fi

for R1 in "${R1_FILES[@]}"; do
    R2="${R1%_1.fastq}_2.fastq"

    if [[ ! -s "$R2" ]]; then
        echo "Missing or empty mate file: $R2" >&2
        exit 1
    fi

    NAME="${R1##*/}"
    NAME="${NAME%_1.fastq}"
    OUT1="$CLEAN_DIR/${NAME}_1_TRIMMED.fastq"
    OUT2="$CLEAN_DIR/${NAME}_2_TRIMMED.fastq"

    echo "Trimming paired reads for $NAME"
    fastp \
        --in1 "$R1" --in2 "$R2" \
        --out1 "$OUT1" --out2 "$OUT2" \
        --trim_front1 8 --trim_front2 8 \
        --trim_tail1 20 --trim_tail2 20 \
        --n_base_limit 0 \
        --length_required 100 \
        --average_qual 20 \
        --json /dev/null \
        --html /dev/null
done
```

#### SCRIPT: `scripts/02_fastp.slurm`

```bash
#!/bin/bash
#SBATCH --job-name=kbrevis_fastp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1-00:00:00
#SBATCH --mem=120G
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/fastp_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/fastp_%j.err

bash ./scripts/02_fastp.sh 

```

**What fastp does.** fastp is a fast quality-control tool for Illumina reads. For paired
reads it detects and removes adapter sequence (by finding where the two reads of a pair
overlap), trims bases, and throws out pairs that are too short, too low quality, or contain
ambiguous bases (N). Removing these early stops sequencing artefacts from creating false
contigs in the assembly.

**Parameters and why.**
- `--trim_front1 8 --trim_front2 8`: cut the first 8 bases of every read. The libraries were
  made with the IDT xGen ssDNA & Low-Input kit, whose "Adaptase" step adds a short
  low-complexity tail (up to about 10 nt) to the DNA fragments. Removing the first bases
  also removes the biased first cycles of Illumina reads.
- `--trim_tail1 20 --trim_tail2 20`: cut the last 20 bases, where NovaSeq quality is lowest.
- `--n_base_limit 0`: discard any read containing an N.
- `--length_required 100`: discard reads shorter than 100 bp after trimming. Reads start at
  151 bp, so at most 28 bases are removed by the fixed trims; this mainly removes reads
  shortened by adapter trimming.
- `--average_qual 20`: discard reads whose mean Phred quality is below 20 (1 % error).
- `--json /dev/null --html /dev/null`: don't write per-sample reports. (Keeping them would
  let MultiQC summarise the trimming; worth turning on next time.)

**Deviation from the paper.** Booker et al. quality-filtered with BBDuk
(`qtrim=rl trimq=17`, phiX removal, GC 5–95 %) and then error-corrected with bfc, and the
reads in SRA are probably already filtered. These fastp settings are mine, are stricter than
the paper's, and were chosen by hand rather than taken from a publication. In particular the
fixed 8-base front trim and 20-base tail trim remove sequence from every read, even good
reads.

**Run and result.** `sbatch scripts/02_fastp.slurm` from the repo root on bora: job
**577472**, bora `batch`, 4 cores / 120 GB, ran 2 h 15 min, COMPLETED. Output: 40
`*_TRIMMED.fastq` files in `~/scr10/data/kbrevis/trim/` (1.6 TB). An earlier attempt (job
573521) failed; its log was deleted. fastp is single-sample-at-a-time here and uses its
default 3 worker threads; adding `--thread "$SLURM_CPUS_PER_TASK"` would use the 4 cores
requested.

# Step 03. Download reference genomes of likely contaminants

#### SCRIPT: `scripts/03_download_refs.sh`

```bash
#!/bin/bash
set -euo pipefail

REF_DIR="$HOME/scr10/data/kbrevis/reference"
mkdir -p "$REF_DIR"

download_and_extract() {
    local name="$1"
    local accession_type="$2"
    local accession="$3"
    local zip="$REF_DIR/${name}_ref_genome.zip"
    local stage

    stage=$(mktemp -d "$REF_DIR/.${name}.XXXXXX")

    if [[ "$accession_type" == "taxon" ]]; then
        datasets download genome taxon "$accession" \
            --reference \
            --include genome \
            --filename "$zip"
    else
        datasets download genome accession "$accession" \
            --include genome \
            --filename "$zip"
    fi

    unzip -q -o "$zip" -d "$stage"

    local -a fasta_files=()
    while IFS= read -r -d '' file; do
        fasta_files+=("$file")
    done < <(find "$stage" -type f -name '*.fna' -print0)

    if ((${#fasta_files[@]} != 1)); then
        echo "Expected one .fna for $name, found ${#fasta_files[@]}; check the package or use a specific assembly accession." >&2
        return 1
    fi

    mv "${fasta_files[0]}" "$REF_DIR/${name}_ref_genome.fna"
    rm -rf -- "$stage"
}

download_and_extract dog taxon "Canis familiaris"
download_and_extract cat accession "GCF_018350175.1"
download_and_extract mouse accession "GCF_000001635.27"
download_and_extract human accession "GCF_000001405.40"

# lastly, delete the zip file 
rm -rf dog_ref_genome.zip
rm -rf cat_ref_genome.zip
rm -rf mouse_ref_genome.zip
rm -rf human_ref_genome.zip


# get K. brevis seqs - unclear if this is whole genome but I doubt it 
fasterq-dump SRR14099984 \
    --split-files \
    --outdir "$REF_DIR" \
    --progress

mv "$REF_DIR/SRR14099984_1.fastq" "$REF_DIR/kbrevis_reads_1.fastq"
mv "$REF_DIR/SRR14099984_2.fastq" "$REF_DIR/kbrevis_reads_2.fastq"

# now need to turn the kbrevis reads into .fna files so that they can be used in the next 
# script where you want to demove reads that match them. 
# convert R1 into fasta

module load miniforge3
eval "$(conda shell.bash hook)"
conda activate seqtk-env

seqtk seq -a "$REF_DIR/kbrevis_reads_1.fastq" >> "$REF_DIR/kbrevis_reads_1.fa"
seqtk seq -a "$REF_DIR/kbrevis_reads_2.fastq" >> "$REF_DIR/kbrevis_reads_2.fa"
cat "$REF_DIR/kbrevis_reads_1.fa" "$REF_DIR/kbrevis_reads_2.fa" > "$REF_DIR/kbrevis_ref.fna"
# this results in a .fna file in your reference directory that you can now pull from. It's 
# not an assembled genome, just all of the reads that came from a Kbrevis culture. But it might 
# help to get rid of some of the complexity in your contig-making attempts when you run metaSPAdes.
```

#### SCRIPT: `scripts/03_download_refs.slurm`

```bash
#!/bin/bash
#SBATCH --job-name=download_refs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1-00:00:00
#SBATCH --mem=120G
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/download_refs_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/download_refs_%j.err

module load miniforge3 
echo "$(conda info --base)"
eval "$(conda shell.bash hook)" 
conda activate datasets-env
conda env export --no-builds > datasets-env.yml

bash ./scripts/03_download_refs.sh 
```

#### SCRIPT: `scripts/03b_download_phix.sh`

```bash
#!/bin/bash
set -euo pipefail

module load miniforge3

REF_DIR="$HOME/scr10/data/kbrevis/reference"
PHIX_FASTA="$REF_DIR/phix174.fna"
mkdir -p "$REF_DIR"

if [[ ! -s "$PHIX_FASTA" ]]; then
    echo "Downloading PhiX174 (NC_001422.1)..."
    TMP=$(mktemp -d)
    datasets download virus genome accession NC_001422.1 \
        --filename "$TMP/phix.zip"
    unzip -q "$TMP/phix.zip" -d "$TMP"
    cp "$TMP/ncbi_dataset/data/genomic.fna" "$PHIX_FASTA"
    rm -rf "$TMP"
    echo "Saved $PHIX_FASTA"
else
    echo "$PHIX_FASTA already exists"
fi
```

#### SCRIPT: `scripts/03b_download_phix.slurm`

```bash
#!/bin/bash
#SBATCH --job-name=download_phix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1-00:00:00
#SBATCH --mem=120G
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/download_phix_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/download_phix_%j.err

module load miniforge3 
echo "$(conda info --base)"
eval "$(conda shell.bash hook)" 
conda activate datasets-env
conda env export --no-builds > datasets-env.yml

bash ./scripts/03b_download_phix.sh 
```

**Why.** Booker et al. removed reads matching phiX and mammalian contaminants (human, dog,
mouse, cat) before assembly. **phiX174** is a small bacteriophage whose DNA Illumina adds to
every run as a calibration control; its reads show up in almost every dataset and would
assemble into a perfect "virus" that isn't from the sample. **Human, mouse, cat and dog**
DNA gets into samples from people handling them and from lab reagents. I also added
**reads from a *K. brevis* culture** (SRA run SRR14099984) so that algal host DNA would be
removed too. This step only downloads the references; step 04 does the removal.

**What the tools do.** NCBI's `datasets` command-line tool downloads genome assemblies from
NCBI as a zip file. `fasterq-dump` (step 01) downloads the *K. brevis* reads. `seqtk seq -a`
converts FASTQ to FASTA (drops the quality lines) so the reads can be used as a bowtie2
reference.

**References used.**

| organism | how requested | file in `~/scr10/data/kbrevis/reference/` | size |
|---|---|---|---|
| dog | `taxon "Canis familiaris" --reference` (NCBI's current reference assembly for dog) | `dog_ref_genome.fna` | 2.5 GB |
| cat | `GCF_018350175.1` (F.catus_Fca126_mat1.0) | `cat_ref_genome.fna` | 2.5 GB |
| mouse | `GCF_000001635.27` (GRCm39) | `mouse_ref_genome.fna` | 2.8 GB |
| human | `GCF_000001405.40` (GRCh38.p14) | `human_ref_genome.fna` | 3.3 GB |
| phiX174 | virus accession `NC_001422.1` | `phix174.fna` | 5.5 kb |
| *K. brevis* culture reads | SRA run `SRR14099984`, all reads 1 + 2 as FASTA | `kbrevis_ref.fna` (+ `kbrevis_reads_{1,2}.fastq/.fa`) | 12.4 GB |

**Run.** From the repo root on bora:
```
chmod +x scripts/03_download_refs.sh scripts/03b_download_phix.sh
( module purge; module load slurm/bora; sbatch scripts/03_download_refs.slurm )
( module purge; module load slurm/bora; sbatch scripts/03b_download_phix.slurm )   # job 573541, 4 s, COMPLETED
```
`03_download_refs.sh` needs `fasterq-dump` (step 01) as well as `datasets`. `fasterq-dump` was
on `PATH` from `~/programs/sratoolkit.3.4.1-alma_linux64/bin` (set in `~/.bashrc`), not from
`datasets-env`. The script activates `seqtk-env` itself for the last part. The log of the
`download_refs` job was not kept.

**Caveats.**
- The "lastly, delete the zip file" lines use bare file names (`rm -rf dog_ref_genome.zip`),
  so they look in the folder you submitted from, not in `reference/`. The four
  `*_ref_genome.zip` files (≈ 3 GB) are still in `reference/` and can be deleted.
- The dog reference is chosen by taxon, so a later run may get a newer assembly. Use the
  accession printed in the `datasets` output to pin it.
- **The *K. brevis* "reference" is raw reads from a culture, not an assembled genome.** Lab
  cultures of *K. brevis* are usually not axenic: they carry bacteria and possibly viruses.
  Anything in those reads, including any viruses that infect *K. brevis*, is removed from
  every sample in step 04. Given that this project is looking for *K. brevis* viruses, that
  is a real risk. In practice step 04 removed very few reads (≤ 1.4 %), but the
  *K. brevis*-associated reads removed are worth checking (or step 04 could be re-run without
  `kbrevis_ref.fna`).
- `seqtk seq -a ... >>` appends; re-running the script would double the FASTA files.

# Step 04. Remove contaminant reads (bowtie2 + samtools)

#### SCRIPT: `scripts/04_clean_refs.sh`

```bash
#!/bin/bash
set -euo pipefail

REF_DIR="$HOME/scr10/data/kbrevis/reference"
TRIM_DIR="$HOME/scr10/data/kbrevis/trim"
OUT_DIR="$HOME/scr10/data/kbrevis/clean"
INDEX_PREFIX="$REF_DIR/host_index"

REF_FASTAS=(
    "$REF_DIR/dog_ref_genome.fna"
    "$REF_DIR/cat_ref_genome.fna"
    "$REF_DIR/mouse_ref_genome.fna"
    "$REF_DIR/human_ref_genome.fna"
    "$REF_DIR/phix174.fna"
    "$REF_DIR/kbrevis_ref.fna"
)

for ref in "${REF_FASTAS[@]}"; do
    [[ -s "$ref" ]] || {
        echo "Missing or empty reference FASTA: $ref" >&2
        exit 1
    }
done

mkdir -p "$OUT_DIR"

# Build the combined index once; reuse it on later runs.
if ! bowtie2-inspect "$INDEX_PREFIX" >/dev/null 2>&1; then
    HOST_FASTAS=$(IFS=,; printf '%s' "${REF_FASTAS[*]}")
    bowtie2-build --threads "$SLURM_CPUS_PER_TASK" \
        "$HOST_FASTAS" "$INDEX_PREFIX"
fi

shopt -s nullglob
R1_FILES=("$TRIM_DIR"/*_1_TRIMMED.fastq)

if ((${#R1_FILES[@]} == 0)); then
    echo "No trimmed *_1_TRIMMED.fastq files found in $TRIM_DIR" >&2
    exit 1
fi

for R1 in "${R1_FILES[@]}"; do
    R2="${R1%_1_TRIMMED.fastq}_2_TRIMMED.fastq"

    [[ -s "$R2" ]] || {
        echo "Missing or empty mate file: $R2" >&2
        exit 1
    }

    NAME="${R1##*/}"
    NAME="${NAME%_1_TRIMMED.fastq}"

    OUT1="$OUT_DIR/${NAME}_1_CLEAN.fastq"
    OUT2="$OUT_DIR/${NAME}_2_CLEAN.fastq"

    if [[ -e "$OUT1" || -e "$OUT2" ]]; then
        echo "Output already exists for $NAME; choose another output directory or move the existing files." >&2
        exit 1
    fi

    echo "Screening $NAME against contaminating genomes..."
    bowtie2 -x "$INDEX_PREFIX" \
        -1 "$R1" -2 "$R2" \
        -p "$SLURM_CPUS_PER_TASK" \
        2> "$OUT_DIR/${NAME}_CLEAN.log" |
        samtools view -u - |
        samtools fastq -f 12 -F 2304 \
            -1 "$OUT1" \
            -2 "$OUT2" \
            -0 /dev/null -s /dev/null -n -
done
```

#### SCRIPT: `scripts/04_clean_refs.slurm`

```bash
#!/bin/bash
#SBATCH --job-name=clean_metagenomes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/clean_metagenomes_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/clean_metagenomes_%j.err

module load miniforge3
eval "$(conda shell.bash hook)"
conda activate host-screen-env
conda env export --no-builds > host-screen-env.yml

bash ./scripts/04_clean_refs.sh
```

**What the tools do.** **bowtie2** aligns short reads to a reference genome; here, to all
six contaminant references at once. **samtools** then reads the alignments and writes back
out only the read pairs that did **not** align to any contaminant, as FASTQ. Those are the
"clean" reads used for assembly.

**Parameters and why.**
- `bowtie2-build` on all six references combined into one index (`reference/host_index.*`).
  The references total about 23 GB, so bowtie2 builds a "large" index (`.bt2l`). It is built
  once and reused on re-runs.
- `bowtie2` default settings (end-to-end, `--sensitive`): a read counts as contaminant only
  if the whole read aligns. This is the usual way to screen for host reads.
- `samtools fastq -f 12`: keep only pairs where **both** the read and its mate are unmapped
  (SAM flags 4 + 8). If either mate hits a contaminant, the pair is removed.
- `-F 2304`: skip secondary (256) and supplementary (2048) alignments so each read is written
  once.
- `-n`: keep read names as they are (no `/1`, `/2` added).
- The job used 8 cores and 64 GB on bora's `hima` partition.

**Deviation from the paper.** Booker et al. used BBMap/BBDuk (k-mer matching) against
phiX and the BBTools masked mammalian references. I used bowtie2 alignment against the
full unmasked NCBI genomes and added *K. brevis* reads. Unmasked genomes contain regions
(e.g. repeats, rRNA-like and low-complexity sequence) that can catch some non-contaminant
reads; with the low removal rates below, this has little effect.

**Run and result.** `sbatch scripts/04_clean_refs.slurm` from the repo root: job **607624**,
bora `hima` (hi01), 8 cores / 64 GB, 8 h 28 min, COMPLETED. (An earlier attempt, job 573868,
failed; its log was deleted.) Output: `~/scr10/data/kbrevis/clean/SRRxxx_{1,2}_CLEAN.fastq`
(1.1 TB) and a bowtie2 log per sample. Percentage of reads that matched a contaminant
(`grep 'overall alignment rate' clean/*_CLEAN.log`):

| sample | removed | sample | removed | sample | removed | sample | removed |
|---|---|---|---|---|---|---|---|
| SRR28004137 | 0.04 % | SRR28004142 | 1.38 % | SRR28004147 | 0.03 % | SRR28004152 | 0.01 % |
| SRR28004138 | 0.01 % | SRR28004143 | 0.56 % | SRR28004148 | 0.01 % | SRR28004153 | 0.02 % |
| SRR28004139 | 0.00 % | SRR28004144 | 0.02 % | SRR28004149 | 0.16 % | SRR28004154 | 0.01 % |
| SRR28004140 | 0.25 % | SRR28004145 | 0.11 % | SRR28004150 | 0.08 % | SRR34748663 | 0.09 % |
| SRR28004141 | 0.01 % | SRR28004146 | 0.03 % | SRR28004151 | 0.03 % | SRR34748664 | 0.16 % |

So contamination was very low; nearly all reads went on to assembly. The log doesn't say
which reference the removed reads hit.

**Note.** `04_clean_refs.slurm` does not set `set -euo pipefail` itself, and
`04_clean_refs.sh` refuses to overwrite existing `_CLEAN.fastq` files, so delete or move them
before re-running.

# Step 05a. Viral assembly with metaviralSPAdes

#### SCRIPT: `scripts/05_metaspades.sh` (submits one job per sample)

```bash
#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(dirname "$SCRIPT_DIR")"
READ_DIR="$HOME/scr10/data/kbrevis/clean"
OUT_DIR="$HOME/scr10/data/kbrevis/contigs"
SLURM_SCRIPT="$SCRIPT_DIR/05_metaspades.slurm"


[[ -d "$READ_DIR" ]] || {
    echo "Corrected read directory not found: $READ_DIR" >&2
    exit 1
}
[[ -f "$SLURM_SCRIPT" ]] || {
    echo "Slurm script not found: $SLURM_SCRIPT" >&2
    exit 1
}

mkdir -p "$OUT_DIR" "$REPO_DIR/logs"

shopt -s nullglob
R1_FILES=(
    "$READ_DIR"/*_1_*.fastq
)

((${#R1_FILES[@]} > 0)) || {
    echo "No paired R1 FASTQ files found in $READ_DIR" >&2
    exit 1
}

for R1 in "${R1_FILES[@]}"; do
    case "$R1" in
        *_1_*.fastq)
            FILE_NAME="${R1##*/}"
            NAME="${FILE_NAME%_1_*}"
            SUFFIX="${FILE_NAME#"$NAME"_1_}"
            R2="$READ_DIR/${NAME}_2_${SUFFIX}"
    esac

    [[ -s "$R1" ]] || {
        echo "Missing or empty R1 file: $R1" >&2
        exit 1
    }
    [[ -s "$R2" ]] || {
        echo "Missing or empty R2 mate for $R1: $R2" >&2
        exit 1
    }

    if [[ -s "$OUT_DIR/${NAME}_contigs.fasta" ]]; then
        echo "[$NAME] already assembled, skipping"
        continue
    fi

    echo "[$NAME] submitting MetaSPAdes job"
    sbatch \
        --job-name="metaspades_${NAME}" \
        "$SLURM_SCRIPT" "$NAME" "$R1" "$R2"
done
```

#### SCRIPT: `scripts/05_metaspades.slurm` (assembles one sample)

```bash
#!/bin/bash
#SBATCH --partition=serial
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=2-00:00:00
#SBATCH --mem=360G
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/05_metaspades_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/05_metaspades_%j.err
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/05_metaspades_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/05_metaspades_%j.err

set -euo pipefail

THREADS="${SLURM_CPUS_PER_TASK:?Submit with sbatch}"
export OMP_NUM_THREADS="$THREADS"

if (($# != 3)); then
    echo "Usage: sbatch $0 SAMPLE_NAME R1.fastq R2.fastq" >&2
    exit 2
fi

NAME="$1"
R1="$2"
R2="$3"
OUT_DIR="$HOME/scr10/data/kbrevis/contigs"
SAMPLE_OUT="$OUT_DIR/$NAME"

mkdir -p "$OUT_DIR"

[[ -s "$R1" ]] || {
    echo "Missing or empty R1 file: $R1" >&2
    exit 1
}
[[ -s "$R2" ]] || {
    echo "Missing or empty R2 file: $R2" >&2
    exit 1
}

if [[ -s "$OUT_DIR/${NAME}_contigs.fasta" ]]; then
    echo "[$NAME] already assembled, skipping"
    exit 0
fi
if [[ -e "$SAMPLE_OUT" ]]; then
    echo "Output directory already exists for $NAME: $SAMPLE_OUT" >&2
    echo "Move or remove it before resubmitting this sample." >&2
    exit 1
fi

module load miniforge3
eval "$(conda shell.bash hook)"
conda activate metaspades-env

command -v metaviralspades.py >/dev/null 2>&1 || {
    echo "metaviralspades.py not found in metaspades-env" >&2
    exit 1
}

# Assemble paired metagenomic DNA reads with metaviralSPAdes.
echo "[$NAME] assembling paired reads"
metaviralspades.py \
    --only-assembler \
    -1 "$R1" \
    -2 "$R2" \
    -o "$SAMPLE_OUT" \
    -t "$THREADS" \
    -m 340

cp "$SAMPLE_OUT/contigs.fasta" "$OUT_DIR/${NAME}_contigs.fasta"
```

**What metaviralSPAdes does.** An assembler joins overlapping short reads into longer
sequences (contigs), ideally whole genomes. SPAdes builds a de Bruijn graph from the reads.
**metaviralSPAdes** is the SPAdes mode for finding viruses in metagenomes: it looks for
paths in the graph that look like complete viral genomes, especially **circular** ones
(many phage genomes are circular, or assemble as circles because they are terminally
redundant). It then checks each candidate's genes against a database of viral, plasmid and
chromosomal gene profiles (the viralVerify step) and **outputs only contigs it classifies as
viral**. Each contig name ends in `type_circular`, `type_linear` or `type_linearrepeat`.

**Parameters and why.**
- `--only-assembler`: skip SPAdes' built-in read error correction (BayesHammer). Booker et
  al. had already corrected the reads with bfc before uploading them to SRA (their
  Supplemental Methods), and error correction is the slowest, most memory-hungry part of
  SPAdes.
- `-t 32 -m 340` (32 threads, 340 GB memory cap) with `--mem=360G` in Slurm: `-m` is set a
  bit below the Slurm limit so SPAdes stops with an error message before Slurm kills it.
- k-mer sizes: SPAdes' defaults (the paper used metaSPAdes with `-k 21,33,55,77,99,127`;
  metaviralSPAdes chooses its own k-mers).
- One job per sample (20 separate jobs) because samples are assembled separately, as in the
  paper.

**Deviation from the paper.** Booker et al. assembled with **metaSPAdes 3.15.4**. I used
**metaviralSPAdes from SPAdes 4.3.0** first, because I'm looking for viruses. metaviralSPAdes
was designed and tested mostly on phages and other small viruses. Giant viruses (genomes of
hundreds of kb) are often fragmented into many contigs, and metaviralSPAdes may not report
them, which is why step 05b also runs plain metaSPAdes.

**Run.**
```
cd ~/K_brevis_viromes
bash scripts/05_metaspades.sh      # submits one sbatch per sample, skips finished samples
```

**What actually ran.**
- Jobs **610419–610438** (bora `batch`, 8 cores / 64 GB, `-m 60`, the script's settings at
  the time): 17 samples finished in 3 h 20 min – 14 h. (Each had a few failed attempts of a
  few seconds to minutes just before, while the script was being debugged.)
- The three deep CH1 samples (SRR28004142–144, ~160–180 Gbases) did not fit: reruns at 96 GB
  (611253, 611255) were cancelled after ~3 h, and reruns at ~50 GB (611872–611874) ran
  **out of memory**.
- The script was then changed to 32 cores / 360 GB / `-m 340` and these three were rerun on
  **kuro `serial`**: jobs **109630–109632**, COMPLETED in 25 h, 14 h and 10 h.
- Earlier versions of the script wrote into `~/K_brevis_viromes/output/contigs/` (in home);
  that output was removed and the scripts now write to scratch.

**Result.** `~/scr10/data/kbrevis/contigs/SRRxxx_contigs.fasta` for all 20 samples (5,625
contigs in total; per-sample counts are in the sample table in section C), plus each
sample's full SPAdes run folder `contigs/SRRxxx/` (90 GB in total). I expected every contig
to be circular, but counting contig types shows otherwise:

```
cat ~/scr10/data/kbrevis/contigs/*_contigs.fasta | grep '>' | grep -o 'type_[a-z]*' | sort | uniq -c
#    904 type_circular
#     12 type_linearrepeat
#   4709 type_linear
```

**Caveats.** metaviralSPAdes' "viral" call is a prediction; it can miss viruses and can
include plasmids or prophage-bearing bacterial fragments. "Circular" means the assembly
graph closes on itself, which usually but not always means a complete genome. CheckV
(the removed CheckV-database step) is needed to judge completeness.

# Step 05b. Whole-metagenome assembly with metaSPAdes (for comparison)

#### SCRIPT: `scripts/05_metaspades_nov.sh`

```bash
#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(dirname "$SCRIPT_DIR")"
READ_DIR="$HOME/scr10/data/kbrevis/clean"
OUT_DIR="$HOME/scr10/data/kbrevis/contigs_no_v"
SLURM_SCRIPT="$SCRIPT_DIR/05_metaspades_nov.slurm"


[[ -d "$READ_DIR" ]] || {
    echo "Corrected read directory not found: $READ_DIR" >&2
    exit 1
}
[[ -f "$SLURM_SCRIPT" ]] || {
    echo "Slurm script not found: $SLURM_SCRIPT" >&2
    exit 1
}

mkdir -p "$OUT_DIR" "$REPO_DIR/logs"

shopt -s nullglob
R1_FILES=(
    "$READ_DIR"/*_1_*.fastq
)

((${#R1_FILES[@]} > 0)) || {
    echo "No paired R1 FASTQ files found in $READ_DIR" >&2
    exit 1
}

for R1 in "${R1_FILES[@]}"; do
    case "$R1" in
        *_1_*.fastq)
            FILE_NAME="${R1##*/}"
            NAME="${FILE_NAME%_1_*}"
            SUFFIX="${FILE_NAME#"$NAME"_1_}"
            R2="$READ_DIR/${NAME}_2_${SUFFIX}"
    esac

    [[ -s "$R1" ]] || {
        echo "Missing or empty R1 file: $R1" >&2
        exit 1
    }
    [[ -s "$R2" ]] || {
        echo "Missing or empty R2 mate for $R1: $R2" >&2
        exit 1
    }

    if [[ -s "$OUT_DIR/${NAME}_contigs.fasta" ]]; then
        echo "[$NAME] already assembled, skipping"
        continue
    fi

    echo "[$NAME] submitting MetaSPAdes job"
    sbatch \
        --job-name="metaspades_no_v_${NAME}" \
        "$SLURM_SCRIPT" "$NAME" "$R1" "$R2"
done
```

#### SCRIPT: `scripts/05_metaspades_nov.slurm`

```bash
#!/bin/bash
#SBATCH --partition=serial
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=2-00:00:00
#SBATCH --mem=360G
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/05_metaspades_no_v_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/05_metaspades_no_v_%j.err

set -euo pipefail

THREADS="${SLURM_CPUS_PER_TASK:?Submit with sbatch}"
export OMP_NUM_THREADS="$THREADS"

# SPAdes' -m is a cap, in GB, that SPAdes applies to its OWN allocator. It used to be
# hardcoded at 340 while #SBATCH --mem said 360G, and the two have to be kept consistent:
# Slurm's cgroup refuses any allocation past --mem, and SPAdes' internal accounting runs
# roughly 10% below its true RSS (job 92248 reported "325G / 340G" at the moment the
# kernel refused it at 360 GiB RSS). So derive -m from the actual allocation and keep 10%
# headroom. That way raising --mem on the sbatch command line is enough, and SPAdes aborts
# with its own "memory limit exceeded" message instead of dying inside mimalloc.
if [[ -n "${SLURM_MEM_PER_NODE:-}" ]]; then
    SPADES_MEM=$(( SLURM_MEM_PER_NODE / 1024 * 9 / 10 ))
elif [[ -n "${SLURM_MEM_PER_CPU:-}" ]]; then
    SPADES_MEM=$(( SLURM_MEM_PER_CPU * THREADS / 1024 * 9 / 10 ))
else
    echo "Cannot read the Slurm memory allocation; submit with --mem." >&2
    exit 1
fi
(( SPADES_MEM >= 16 )) || { echo "Allocation too small: -m would be ${SPADES_MEM}G" >&2; exit 1; }

# Resume a run that died part-way instead of starting over. Set to a SPAdes check-point
# ('k21', 'k33', 'k55', 'last'); 'k55' redoes only the final k-mer iteration and keeps the
# finished K21/K33 output in $SAMPLE_OUT. -o, -t and -m may be changed on a restart, the
# input reads may not, so the same SAMPLE/R1/R2 arguments must be passed again. Example:
#   RESTART_FROM=k55 sbatch --export=ALL,RESTART_FROM --mem=950G ... 05_metaspades_nov.slurm ...
RESTART_FROM="${RESTART_FROM:-}"

if (($# != 3)); then
    echo "Usage: sbatch $0 SAMPLE_NAME R1.fastq R2.fastq" >&2
    exit 2
fi

NAME="$1"
R1="$2"
R2="$3"
OUT_DIR="$HOME/scr10/data/kbrevis/contigs_no_v"
SAMPLE_OUT="$OUT_DIR/$NAME"

mkdir -p "$OUT_DIR"

[[ -s "$R1" ]] || {
    echo "Missing or empty R1 file: $R1" >&2
    exit 1
}
[[ -s "$R2" ]] || {
    echo "Missing or empty R2 file: $R2" >&2
    exit 1
}

if [[ -s "$OUT_DIR/${NAME}_contigs.fasta" ]]; then
    echo "[$NAME] already assembled, skipping"
    exit 0
fi
if [[ -n "$RESTART_FROM" ]]; then
    [[ -d "$SAMPLE_OUT" ]] || {
        echo "RESTART_FROM=$RESTART_FROM needs the previous run's output directory," >&2
        echo "but $SAMPLE_OUT does not exist. Drop RESTART_FROM to start from scratch." >&2
        exit 1
    }
elif [[ -e "$SAMPLE_OUT" ]]; then
    echo "Output directory already exists for $NAME: $SAMPLE_OUT" >&2
    echo "Move or remove it before resubmitting this sample," >&2
    echo "or set RESTART_FROM=<check-point> to resume it." >&2
    exit 1
fi

module load miniforge3
eval "$(conda shell.bash hook)"
conda activate metaspades-env

command -v metaspades.py >/dev/null 2>&1 || {
    echo "metaspades.py not found in metaspades-env" >&2
    exit 1
}

# Assemble paired metagenomic DNA reads with metaSPAdes.
#
# A restart takes a DIFFERENT argument list, not an extra flag. SPAdes reads the original
# command line back out of $SAMPLE_OUT/params.txt and refuses to be handed it a second
# time: passing -1/-2 aborts with "you cannot specify input data (-1, -2, -12, --pe-1,
# --pe-2 ...) with --restart-from option!" and --only-assembler aborts the same way (see
# check_options_for_restart_from() in spades_pipeline/options_parser.py). Only -o, -t and
# -m may be given; SPAdes logs them as "with updated parameters", which is what makes the
# new -m take effect. The reads still come from params.txt, so the R1/R2 arguments to this
# script are only used for the sanity checks above on a restart.
if [[ -n "$RESTART_FROM" ]]; then
    [[ -s "$SAMPLE_OUT/params.txt" ]] || {
        echo "No $SAMPLE_OUT/params.txt, so there is no run to restart from." >&2
        echo "Remove $SAMPLE_OUT and resubmit without RESTART_FROM." >&2
        exit 1
    }
    echo "[$NAME] resuming from check-point $RESTART_FROM with ${THREADS} threads, -m ${SPADES_MEM}G"
    metaspades.py \
        --restart-from "$RESTART_FROM" \
        -o "$SAMPLE_OUT" \
        -t "$THREADS" \
        -m "$SPADES_MEM"
else
    echo "[$NAME] assembling paired reads with ${THREADS} threads, -m ${SPADES_MEM}G"
    metaspades.py \
        --only-assembler \
        -1 "$R1" \
        -2 "$R2" \
        -o "$SAMPLE_OUT" \
        -t "$THREADS" \
        -m "$SPADES_MEM"
fi

cp "$SAMPLE_OUT/contigs.fasta" "$OUT_DIR/${NAME}_contigs.fasta"
```

**Why.** metaviralSPAdes outputs only contigs it calls viral, so its output is a filtered
set. **metaSPAdes** assembles everything in the sample (bacteria, algae, viruses) and is
what Booker et al. used. Comparing the two (step 06) shows which viral contigs both
assemblers recover (higher confidence) and which come from only one, and gives the binning
step the full set of contigs it needs to rebuild large, fragmented giant-virus genomes.

**How the script differs from step 05a.** These are copies of `05_metaspades.sh/.slurm`
with three changes: `metaspades.py` instead of `metaviralspades.py`; output to
`~/scr10/data/kbrevis/contigs_no_v/`; job names `metaspades_no_v_<SAMPLE>`. Same reads
(`clean/*_1_CLEAN.fastq` + `*_2_CLEAN.fastq`), same `--only-assembler`, 32 cores, 360 GB,
`-m 340`. (The listing above is the current version, after the memory fix described in
"Step 05b-fix" below replaced the hardcoded `-m 340` with a value derived from `--mem`.)
The only edit Claude made (prompt 2) was a typo in the stdout log name:

```
# before
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/05_metaspades_no_v__%j.out
# after
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/05_metaspades_no_v_%j.out
```

**Run.** The 20 jobs were first submitted to kuro `serial` (jobs 109726–109745). Slurm
estimated they wouldn't start until about Nov 18, so they were moved to gulf's `gz`
partition (64 cores, 1 TB, fastest CPUs). Each gz job was submitted first, and its kuro job
was cancelled only once the gz submission succeeded (prompt 4):

```
module purge; module load slurm/gulf
sbatch -p gz -t 3-00:00:00 --job-name="metaspades_no_v_${NAME}" \
    scripts/05_metaspades_nov.slurm "$NAME" "$R1" "$R2"
# then: module purge; module load slurm/kuro; scancel <old kuro job id>
```

gz jobs **92235–92254**. The `-p gz` and `-t` on the command line override the
`#SBATCH --partition=serial` and 2-day limit inside the script. Status at 12:05 on
2026-11-12: 92235 and 92236 (SRR34748663, SRR34748664) running on gz01, the other 18
pending. Check with `( module purge; module load slurm/gulf; squeue -u $USER )`.

**Risk.** My account only has the `normal` QOS on gz, so any job submitted with the
`gz_preempt` QOS **cancels** mine (no requeue). If that happens, delete the partial
`contigs_no_v/<SAMPLE>/` folder (the script refuses to run if it exists) and resubmit that
sample.

**Expected output.** `~/scr10/data/kbrevis/contigs_no_v/SRRxxx_contigs.fasta` plus the run
folder per sample. The three CH1 samples will take the longest (the metaviralSPAdes runs
took 10–25 h on kuro).

# Step 05b-fix. Why job 92248 (SRR28004143) ran out of memory

Diagnosed by Claude Code on 2026-11-13 (prompt 30 in `PROMPT_LOG.md`).

**What happened.** 17 of the 20 metaSPAdes jobs finished. One failed:

```
module purge; module load slurm/gulf
sacct -u $USER -j 92248 -o jobid,jobname%26,partition,state,elapsed,reqmem,maxrss,exitcode,nodelist
# 92248        metaspades_no_v_SRR280041+   gz   FAILED   13:36:47   360G              68:0   gz01
# 92248.batch  batch                             FAILED   13:36:47          377484792K  68:0   gz01
```

MaxRSS was 377,484,792 KB = **360.0 GiB**, i.e. exactly the `--mem=360G` the job asked
for. The end of `logs/05_metaspades_no_v_92248.out` shows the allocator, not Slurm,
reporting the refusal:

```
mimalloc: error: thread 0x14f5fdbf2640: unable to allocate OS memory
          (67108864 bytes, error code: 12 [Cannot allocate memory], ...)
    rss: current: 456906416128, peak: 363720421376
== Error == system call for: [... spades-core ... K55/configs/config.info ...]
            finished abnormally, OS return value: 68
```

So this is a genuine out-of-memory, not a crash: SPAdes asked the kernel for another 64 MB
and the job's cgroup, capped at 360 GiB, said no. SPAdes' own last progress line read
`310G / 325G`, below its `-m 340` ceiling, which is why its internal limiter never fired —
SPAdes' accounting tracks the memory it has deliberately allocated and runs about 10% below
true RSS, so the cgroup was reached first.

**Which stage.** The job died in `Preliminary Distance Estimation` for the K=55 k-mer
iteration, at the line `Using SMOOTHING distance estimator`, having just logged
`Paired info stats: missing = 358921965`. Distance estimation holds the paired-read
distance histograms for the whole assembly graph in memory, and the smoothing estimator is
the most memory-hungry step in metaSPAdes. K=21 and K=33 had already completed:

```
ls ~/scr10/data/kbrevis/contigs_no_v/SRR28004143/
# K21 (541M, complete)  K33 (486M, complete)  K55 (104K, only configs/)
```

**Why metaviralSPAdes managed the same reads in 106 GB.** This was the puzzling part, and
the answer is that the two programs do not do the same amount of work, even though they
share the same `spades-core` binary. Compare the config files each one loads (first lines
of each `spades.log`):

```bash
grep -m1 "== Running" ~/scr10/data/kbrevis/contigs/SRR28004143/spades.log       # metaviralSPAdes
# ... K21/configs/config.info  mda_mode.info  meta_mode.info  plasmid_mode.info  metaviral_mode.info
grep -m1 "== Running" ~/scr10/data/kbrevis/contigs_no_v/SRR28004143/spades.log  # metaSPAdes
# ... K21/configs/config.info  mda_mode.info  meta_mode.info
```

metaviralSPAdes layers `plasmid_mode` and `metaviral_mode` on top of metagenomic mode, and
those add a **chromosome removal** stage that runs *before* repeat resolving. Listing the
stages of each run makes the consequence obvious:

```bash
grep -E "STAGE ==" ~/scr10/data/kbrevis/contigs/SRR28004143/spades.log | sed 's/.*STAGE ==/STAGE ==/'
```

The metaviralSPAdes run goes `... Repeat Resolving -> Contig Output -> Chromosome Removal
-> Repeat Resolving -> Contig Output -> Chromosome Removal -> ...`, whereas the metaSPAdes
run goes `... Preliminary Simplification -> Preliminary Paired Information Counting ->
Preliminary Distance Estimation` and then the same three stages again for real. In plain
terms: **metaviralSPAdes throws away the cellular (chromosomal) part of the assembly graph
early and does its expensive graph work on the small viral remainder, while metaSPAdes
keeps the entire metagenome graph and does distance estimation and repeat resolving twice
over.** Peak memory for the two runs on this one sample:

```bash
for f in contigs contigs_no_v; do
  grep -oE "^ *[0-9:.]+ +[0-9]+G / [0-9]+G" ~/scr10/data/kbrevis/$f/SRR28004143/spades.log \
    | awk '{print $4}' | tr -d G | sort -n | tail -1
done
# 106   <- metaviralSPAdes (completed)
# 325   <- metaSPAdes (killed, still rising)
```

So "metaviralSPAdes worked on this file" does not imply metaSPAdes can. Expect metaSPAdes
to need several times the memory of metaviralSPAdes on the same reads, and more on a
complex sample.

**Why only this sample.** Three samples have roughly ten times as much sequence as the
other 17:

```bash
ls -laS ~/scr10/data/kbrevis/clean/*_1_CLEAN.fastq | head -4
# SRR28004144  145.5 GB
# SRR28004142  125.9 GB
# SRR28004143  118.0 GB
# SRR28004148   13.7 GB   <- and everything below it
```

and those are exactly the three with high peak memory (178 GB, 120 GB-and-rising, 325 GB);
every other sample peaked between 27 and 66 GB. Note that **SRR28004144 is the largest file
and it finished** (job 92247, 14:29:50, 178 GB peak), so raw read count is not the whole
story: assembly-graph complexity is, and SRR28004143 happens to produce a graph with
~359 million unresolved paired-info records. That is why a per-sample memory request, not a
single global one, is the right approach here.

**Fix 1: derive SPAdes' `-m` from the Slurm allocation.** The old script hardcoded `-m 340`
while `#SBATCH --mem` said `360G`. Two numbers that have to agree, in two places, is a bug
waiting to happen, and raising `--mem` on the command line silently left `-m` at 340. The
script now computes it:

```bash
# before
    -m 340

# after
if [[ -n "${SLURM_MEM_PER_NODE:-}" ]]; then
    SPADES_MEM=$(( SLURM_MEM_PER_NODE / 1024 * 9 / 10 ))
elif [[ -n "${SLURM_MEM_PER_CPU:-}" ]]; then
    SPADES_MEM=$(( SLURM_MEM_PER_CPU * THREADS / 1024 * 9 / 10 ))
else
    echo "Cannot read the Slurm memory allocation; submit with --mem." >&2
    exit 1
fi
...
    -m "$SPADES_MEM"
```

The 10% headroom is what closes the gap between SPAdes' accounting and its true RSS, so
that SPAdes now hits its own limit and exits with a readable "memory limit exceeded"
message instead of being refused inside mimalloc. Checked:

```bash
SLURM_MEM_PER_NODE=972800 bash -c 'echo $(( SLURM_MEM_PER_NODE / 1024 * 9 / 10 ))'  # 855
SLURM_MEM_PER_NODE=368640 bash -c 'echo $(( SLURM_MEM_PER_NODE / 1024 * 9 / 10 ))'  # 324
```

The `#SBATCH --mem=360G` default is unchanged, because 17 of the 20 samples peak under
70 GB and raising the default would make every job queue for a 1 TB node.

**Fix 2: allow a restart instead of a 20-hour redo.** SPAdes can resume from a check-point,
and K21 and K33 are intact on disk, so there is no reason to repeat the ~6.7 hours they
took. `metaspades.py --continue` ignores every option except `-o` (so it would keep the old
340 GB cap and be useless here), but `--restart-from <check-point>` accepts new `-o`, `-t`
and `-m` values:

```bash
metaspades.py --help | grep -A2 restart-from
#   --restart-from <cp>   restart run with updated options and from the specified
#                         check-point ('ec', 'as', 'k<int>', 'mc', 'last')
```

A new `RESTART_FROM` environment variable exposes this, and relaxes the
"output directory already exists" guard when it is set:

```bash
RESTART_FROM="${RESTART_FROM:-}"
...
if [[ -n "$RESTART_FROM" ]]; then
    [[ -d "$SAMPLE_OUT" ]] || { ...exit 1; }
elif [[ -e "$SAMPLE_OUT" ]]; then
    ...exit 1
fi
```

**A restart takes a different argument list, not an extra flag.** The first version of this
change simply added `--restart-from` to the existing `metaspades.py --only-assembler -1
"$R1" -2 "$R2" ...` call. That cannot work, and reading SPAdes' own source is what caught
it before the job was submitted:

```bash
grep -n "cannot specify" ~/.conda/envs/metaspades-env/share/spades/spades_pipeline/options_parser.py
#  906: "you cannot specify input data (-1, -2, -12, --pe-1, --pe-2 ...) with --restart-from option!"
#  930: "you cannot specify --only-assembler with --restart-from option!"
```

`check_options_for_restart_from()` rejects the reads *and* `--only-assembler`, because
SPAdes recovers the original command line from `$SAMPLE_OUT/params.txt` and will not be
handed it twice. Only `-o`, `-t` and `-m` may be supplied. Both branches were then tested
against a throwaway directory holding nothing but a copy of `params.txt`, so the real
46 GB run directory could not be touched:

```bash
mkdir -p /tmp/restart_test
cp ~/scr10/data/kbrevis/contigs_no_v/SRR28004143/params.txt /tmp/restart_test/
conda activate metaspades-env

# A. the wrong way - reads and --only-assembler passed again:
metaspades.py --only-assembler --restart-from k55 -1 "$R1" -2 "$R2" -o /tmp/restart_test -t 1 -m 16
# == Error == you cannot specify input data (-1, -2, -12, --pe-1, --pe-2 ...) with --restart-from option!
# SPAdes finished with the following error code: 67.

# B. the right way - only -o/-t/-m; 'k99' is a deliberately invalid check-point so SPAdes
#    stops before doing any actual assembly:
metaspades.py --restart-from k99 -o /tmp/restart_test -t 1 -m 16
# Restored from Command line: ... metaspades.py --only-assembler -1 ... -2 ... -t 32 -m 340
# Restart-from=k99
# with updated parameters: -t 1 -m 16
#   Memory limit (in Gb): 16          <- the new -m really does take effect
# == Error == failed to restart from k99 because this stage was not specified!
```

Test B is the important one: it proves the restart reads the reads and `--only-assembler`
back out of `params.txt` by itself, and that the **new `-m` overrides the stored 340**,
which is the entire point of the exercise. The `k99` error also confirms the check-point
name is validated against `k: [21, 33, 55]`, so `k55` is accepted.

So the call is written as two separate branches rather than one call with an optional flag:

```bash
if [[ -n "$RESTART_FROM" ]]; then
    [[ -s "$SAMPLE_OUT/params.txt" ]] || { ...exit 1; }
    metaspades.py --restart-from "$RESTART_FROM" -o "$SAMPLE_OUT" -t "$THREADS" -m "$SPADES_MEM"
else
    metaspades.py --only-assembler -1 "$R1" -2 "$R2" -o "$SAMPLE_OUT" \
        -t "$THREADS" -m "$SPADES_MEM"
fi
```

On a restart the `R1`/`R2` arguments to the script are therefore used only for the
existence checks near the top; the reads SPAdes actually assembles come from `params.txt`.
`bash -n scripts/05_metaspades_nov.slurm` passes.

**Fix 3: move this one sample to a 1 TB non-preemptible node.** Per `~/AGENTS.md` section 3:

```bash
~/hpc/best_cluster.sh -c 32 -m 950 -H 30
# bora  hima  11-14 15:51  11-17 10:31  module load slurm/bora; sbatch -p hima -N 1 -c 32 \
#       -C hi --exclude=hi04,hi05,hi07 --mem=950G -t 72:00:00   [32c/~251GB (hi03 ~1TB)]
# gulf  gz    11-20 11:07  11-21 17:07  ... [64c/1TB; PREEMPTIBLE by gz_preempt QOS]
# kuro serial NO-FIT ...  [64c/~376GB]
```

`hima` wins on both counts. `hi03` really does have 1 TB (`RealMemory=1031460` MB, versus
257,317 MB on hi01–hi07's other nodes), `scontrol show partition hima` reports
`PreemptMode=OFF` and `MaxTime=3-00:00:00`, and gz is booked out for a week. Preemption
matters a great deal for a job of this length: a ~40 hour gz job spends 40 hours exposed to
cancellation with no requeue, and this one has already burned 13.6 hours once. kuro
`serial` (376 GB) and gulf `batch` highmem (`gu01-02`, 503 GB but only 16 cores) cannot fit
950 GB at all, and astral's 503 GB node `as01` is a GPU node that CPU jobs must stay off
(`~/AGENTS.md` rule 2).

**Submitted** on 2026-11-13 at 11:35 as job **625965** (prompt 31):

```bash
cd ~/K_brevis_viromes
R1=~/scr10/data/kbrevis/clean/SRR28004143_1_CLEAN.fastq
R2=~/scr10/data/kbrevis/clean/SRR28004143_2_CLEAN.fastq
( module purge; module load slurm/bora
  sbatch --job-name=metaspades_no_v_SRR28004143 \
         -p hima -N 1 -c 32 -C hi --exclude=hi04,hi05,hi07 \
         --mem=950G -t 72:00:00 \
         --export=ALL,RESTART_FROM=k55 \
         scripts/05_metaspades_nov.slurm SRR28004143 "$R1" "$R2" )
# Submitted batch job 625965
```

Checked immediately before submitting, `hi03` was completely free (`CPUAlloc=0`,
`AllocMem=0`, `State=IDLE`, `RealMemory=1031460`) and `scontrol show reservation` listed no
reservation covering `hi[01-07]` (the three standing bora reservations are on `bo11` and
`bo23` only). The job went to `PD (Priority)` behind one pending 16 GB job. Track it with:

```bash
( module purge; module load slurm/bora; squeue -j 625965 -o "%.9i %.28j %.2t %.9M %R" )
tail -f ~/K_brevis_viromes/logs/05_metaspades_no_v_625965.out
```

The first thing to look for in that log is the line
`[SRR28004143] resuming from check-point k55 with 32 threads, -m 855G`, followed by
`Restart-from=k55` and `with updated parameters:` in
`~/scr10/data/kbrevis/contigs_no_v/SRR28004143/spades.log`. Anything else on the first
screenful means the restart was rejected rather than started.

> **This job did not run to completion.** The restart started correctly and ran for 5h55m,
> then was cancelled on 2026-11-13 for time reasons along with SRR28004142. See
> "05b-fix-c. Both remaining metaSPAdes jobs cancelled for time" below for what that means
> for the rest of the pipeline.

`--export=ALL,RESTART_FROM=k55` passes the restart point in without editing the script.
`-c 32` matches hi03's 32 cores; `--mem=950G` leaves ~57 GB of the node's 1007 GB for the
OS and gives SPAdes `-m 855`; `-t 72:00:00` is hima's maximum, chosen because bora's
per-core speed is roughly 0.45 of gz's, so the ~15–20 hours of gz work still outstanding
becomes ~35–45 hours there.

**Caveats.**

- 950 GB may still not be enough. The job died *entering* the smoothing distance estimator
  at 360 GB, so 325 GB is a floor, not an estimate of the peak. If 950 GB also fails, the
  remaining options are to subsample the reads for this sample (which changes what is
  assembled and should be reported if done) or to run it with `--disable-rr`, which skips
  repeat resolution entirely. `--disable-rr` is **not** recommended: it would make this
  sample's contigs qualitatively different from the other 19 and invalidate the step 06
  comparison. `--cov-cutoff` is not available, as metaSPAdes rejects it in metagenomic mode.
- **SRR28004142 (job 92249) was heading for the same wall.** By 18 h in it had reached a
  186 GB peak at the `Preliminary Simplification` milestone, where SRR28004143 had been at
  105 GB and SRR28004144 (which finished) at 83 GB, so it needed ~400-640 GB against a
  360 GB cap. It was cancelled before reaching that point; see 05b-fix-c below.
- `seff` cannot be used to right-size a job that was killed for exceeding its limit; the
  MaxRSS it reports is the cap, not the demand.
- A restart re-runs the whole K=55 iteration, including the construction and simplification
  that previously succeeded. It saves the K=21 and K=33 work (~6.7 hours of gz time), not
  all 13.6 hours.

### 05b-fix-d. Resuming those two assemblies, and an OOM caused by over-trimming the memory request (2026-11-21/22)

Written by Claude Code on 2026-11-22 (prompt 57 in `PROMPT_LOG.md`).

The two cancelled assemblies were resumed so that all 20 samples get a real merge
(`MV_ONLY_SAMPLES` is now empty by default in `06_merge_assemblies.sh`). Both SPAdes run
directories had survived with `params.txt` intact and had reached K55, so
`RESTART_FROM=k55` picks up from the last check-point rather than repeating ~30 h of K21
and K33.

**The first attempt (jobs 122119, 122120) was killed by the memory limit after 7-8 h**, and
the mistake was mine. I sized the requests from the "peak so far" recorded in each
cancelled run's `spades.log` - 171 GB for SRR28004142 and 310 GB for SRR28004143 - and
asked for 220 GB and 360 GB. But those runs had been **cancelled part-way through K55 with
memory still climbing**, so their high-water marks were floors, not requirements:

```
SRR28004142:  Memory limit set to 198 Gb  ->  mimalloc: unable to allocate 15674114048 bytes
SRR28004143:  Memory limit set to 324 Gb  ->  mimalloc: unable to allocate 67108864 bytes
```

The original `--mem=360G` had worked only because those runs sat on 1 TB nodes where the
limit was never the binding constraint. **A peak from an incomplete run tells you a lower
bound and nothing more** - the general lesson for `seff`-style right-sizing, which is
sound only for jobs that ran to completion.

Resubmitted with real headroom:

```bash
# SRR28004143 exceeded a 324 GB SPAdes limit outright, so it needs a 1 TB-class node
( module purge; module load slurm/bora; \
  RESTART_FROM=k55 sbatch --export=ALL,RESTART_FROM \
    -p hima -C hi --exclude=hi04,hi05,hi07 -c 32 --mem=900G -t 72:00:00 \
    scripts/05_metaspades_nov.slurm SRR28004143 \
      ~/scr10/data/kbrevis/clean/SRR28004143_{1,2}_CLEAN.fastq )   # job 695880, hi03, -m 810G

# SRR28004142 needs >214 GB; 333 GB is the limit its earlier gz run reached without OOM
( module purge; module load slurm/kuro; \
  RESTART_FROM=k55 sbatch --export=ALL,RESTART_FROM \
    -p serial -c 32 --mem=370G -t 48:00:00 \
    scripts/05_metaspades_nov.slurm SRR28004142 \
      ~/scr10/data/kbrevis/clean/SRR28004142_{1,2}_CLEAN.fastq )   # job 122136, ku58, -m 333G
```

A warning about this trap is now in the `#SBATCH --mem` block of
`05_metaspades_nov.slurm`, since the default of 360 GB is fine for the other 18 samples
and badly wrong for these two.

----------------------------------------------------------------------------------------

### 05b-fix-c. Both remaining metaSPAdes jobs cancelled for time (2026-11-13)

**Decision.** On 2026-11-13 the user chose to **cancel both outstanding metaSPAdes jobs
rather than wait for them**, because the remaining runtime (an estimated ~23 h for
SRR28004143 and ~12.5 h for SRR28004142) was more delay than the project could absorb
before moving on to step 06 and beyond. This was a deliberate time-versus-completeness
trade, not a technical failure, and it is recorded here because it changes what the
downstream steps are computed from (prompt 32 in `PROMPT_LOG.md`).

```bash
( module purge; module load slurm/bora; scancel 625965 )   # SRR28004143, hi03
( module purge; module load slurm/gulf; scancel 92249 )    # SRR28004142, gz02
```

Final states:

```
625965  metaspades_no_v_SRR28004143  CANCELLED  05:55:42   (bora hima, hi03)
 92249  metaspades_no_v_SRR28004142  CANCELLED  18:22:03   (gulf gz, gz02)
 92252  metaspades_no_v_SRR28004137  COMPLETED  02:01:12
 92253  metaspades_no_v_SRR28004138  COMPLETED  02:06:04
 92254  metaspades_no_v_SRR28004139  COMPLETED  01:48:16
```

**What the time estimates were based on.** Both numbers came from stage timings rather than
the hardware table. The speed of `hi03` relative to `gz` was measured directly, because the
restart re-ran K=55 stages the failed gz run had already done:

| K=55 stage | gz01 | hi03 | ratio |
|---|---|---|---|
| de Bruijn graph construction | 1h35m | 2h24m | 1.51 |
| Gap Closer (early) | 30m52s | 33m51s | 1.10 |
| Raw Simplification | 4m12s | 7m21s | 1.75 |
| block total | 2h10m | 3h05m | **1.42** |

So `hi03` is about 1.4x slower than `gz` on this work, not the 2.2x implied by the per-core
speeds in `~/AGENTS.md` section 4 — `hi03` was unshared, whereas `gz02` was running 60 of
its 64 cores for other users. Worth remembering when planning future assemblies. The
restart also paid a one-off 2h31m re-doing `Binary Read Conversion`, which the original
K=55 had done in 0.1 s by reusing K21's `.bin` files.

**SRR28004142 was also heading for the same out-of-memory failure.** Peak memory at the
*same* K=55 milestone (`Preliminary Simplification`) for the three deep-sequenced samples:

```bash
# per-stage peak, extracted from each spades.log
#                                     SRR28004144   SRR28004143   SRR28004142
#  Gap Closer (early)                     71 GB         91 GB        156 GB
#  Raw Simplification                     83 GB        105 GB        186 GB
#  Preliminary Simplification             83 GB        105 GB        186 GB
#  Gap Closer (preliminary)              117 GB        115 GB          -
#  Preliminary Paired Information Count  117 GB        115 GB          -
#  Preliminary Distance Estimation       136 GB        284 GB          -
#  Preliminary Repeat Resolving          145 GB        360 GB (died)   -
#  final                                 178 GB (OK)       -           -
```

SRR28004144 grew 83 → 178 GB (2.1x from that milestone) and finished; SRR28004143 grew
105 → 360 GB (3.4x) and was killed. SRR28004142 begins that same climb from **186 GB**
against a 360 GB cap, so it needs roughly 400 GB on SRR28004144's growth factor and
~640 GB on SRR28004143's. It would have exceeded its limit on either trajectory, about four
hours later, during `Preliminary Distance Estimation` — the stage where SRR28004143 jumped
from 115 GB to 284 GB in one step. It was also still running with the old hardcoded
`-m 340`, having been submitted before the fix above.

**Consequence for the analysis: n = 18, not 20.**

```bash
ls ~/scr10/data/kbrevis/contigs/*_contigs.fasta       | wc -l   # 20  metaviralSPAdes
ls ~/scr10/data/kbrevis/contigs_no_v/*_contigs.fasta  | wc -l   # 18  metaSPAdes
```

**SRR28004142 and SRR28004143 have a metaviralSPAdes assembly but no metaSPAdes
assembly.** Every other sample has both. No script needs changing for this —
`06_merge_assemblies.sh` already skips a sample that is missing one of the two:

```bash
        echo "[$NAME] no metaSPAdes assembly yet ($MS), skipping"
        continue
```

**But note what that skip implies downstream.** Step 06 writes nothing for those two
samples, and `08_gvmags.sh` builds its array from `contigs_merged/*_merged_contigs.fasta`:

```bash
for f in "$MERGED_DIR"/*_merged_contigs.fasta; do
```

so **SRR28004142 and SRR28004143 drop out of the pipeline entirely from step 06 onward**,
and their existing metaviralSPAdes contigs are not used at all. That is more than a
degradation of two samples: together with SRR28004144 they are the three most deeply
sequenced libraries in the dataset (118-145 GB of cleaned reads each, versus 5-14 GB for
the rest), so the GVMAG, vOTU and phylogeny results will be built from the 17 shallow
libraries plus SRR28004144. **Any write-up of the results must state this**, and it should
not be described as a re-analysis of all 20 runs.

If those two samples should instead be carried forward on their metaviralSPAdes contigs
alone, step 06 would need a documented fallback that copies the metaviralSPAdes contigs
through as the "merged" set and labels them so they are never mistaken for two-assembler
consensus. That has not been done, because it changes what the merged set means and is the
user's call.

**The partial runs were kept, not deleted:**

```bash
du -sh ~/scr10/data/kbrevis/contigs_no_v/SRR28004142   # 54G
du -sh ~/scr10/data/kbrevis/contigs_no_v/SRR28004143   # 47G
```

Either can be resumed later without repeating K21 and K33, which is most of a day's work
each. For SRR28004143, K=55 would restart from the beginning; for SRR28004142, likewise:

```bash
cd ~/K_brevis_viromes
S=SRR28004143      # or SRR28004142
R1=~/scr10/data/kbrevis/clean/${S}_1_CLEAN.fastq
R2=~/scr10/data/kbrevis/clean/${S}_2_CLEAN.fastq
( module purge; module load slurm/gulf
  sbatch --job-name="metaspades_no_v_${S}" -p gz -N 1 -c 32 --mem=900G -t 3-00:00:00 \
         --export=ALL,RESTART_FROM=k55 \
         scripts/05_metaspades_nov.slurm "$S" "$R1" "$R2" )
```

`gz01`/`gz02` have ~1007 GB each, so `--mem=900G` fits there as well as on `hi03`; run
`~/hpc/best_cluster.sh -c 32 -m 900 -H 20` first to see which will finish soonest, and
remember `gz` is preemptible (`~/AGENTS.md` rule 1). **Two deadlines apply if these are to
be revived:** `/sciclone/scr10` purges files not accessed for 90 days, and the cleaned
reads in `clean/` are on the same filesystem, so both the partial run directories and their
input reads expire around **2026-02-11** unless touched.

# Step 06. Compare and merge the two assemblies

Written by Claude Code on 2026-11-12 (prompt 2). Waiting for step 05b to finish; tested on
stand-in data only.

#### SCRIPT: `scripts/06_merge_assemblies.sh` (submits one job per sample)

```bash
#!/bin/bash
# Submit one 06_merge_assemblies.slurm job per sample that has BOTH a
# metaviralSPAdes assembly (contigs/) and a metaSPAdes assembly (contigs_no_v/).
# Run this after the 05_metaspades_nov.sh jobs have finished.
#
# Exception: the samples listed in MV_ONLY_SAMPLES have no metaSPAdes assembly and never
# will (their metaSPAdes jobs were cancelled for time - README "05b-fix-c"), so they are
# submitted with MV_ONLY=1 and their metaviralSPAdes contigs are used as-is. Any OTHER
# sample missing a metaSPAdes assembly is still skipped, because that means step 05b has
# not finished for it yet. If one of these samples is ever reassembled with metaSPAdes, the
# slurm script refuses MV_ONLY=1 and the sample should be taken off this list for a normal
# merge. Override for a one-off run with, e.g., MV_ONLY_SAMPLES="" bash 06_merge_assemblies.sh
set -euo pipefail

MV_ONLY_SAMPLES="${MV_ONLY_SAMPLES-SRR28004142 SRR28004143}"

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(dirname "$SCRIPT_DIR")"
MV_DIR="$HOME/scr10/data/kbrevis/contigs"
MS_DIR="$HOME/scr10/data/kbrevis/contigs_no_v"
OUT_DIR="$HOME/scr10/data/kbrevis/contigs_merged"
SLURM_SCRIPT="$SCRIPT_DIR/06_merge_assemblies.slurm"

mkdir -p "$OUT_DIR" "$REPO_DIR/logs"

shopt -s nullglob
for MV in "$MV_DIR"/*_contigs.fasta; do
    FILE_NAME="${MV##*/}"
    NAME="${FILE_NAME%_contigs.fasta}"
    MS="$MS_DIR/${NAME}_contigs.fasta"

    MV_ONLY=0
    if [[ " $MV_ONLY_SAMPLES " == *" $NAME "* ]]; then
        MV_ONLY=1
    fi

    if [[ ! -s "$MS" ]] && ((!MV_ONLY)); then
        echo "[$NAME] no metaSPAdes assembly yet ($MS), skipping"
        continue
    fi
    if [[ -s "$OUT_DIR/${NAME}_merged_contigs.fasta" ]]; then
        echo "[$NAME] already merged, skipping"
        continue
    fi

    if ((MV_ONLY)); then
        echo "[$NAME] submitting metaviralSPAdes-only job (no metaSPAdes assembly)"
        sbatch --job-name="merge_${NAME}" --export=ALL,MV_ONLY=1 "$SLURM_SCRIPT" "$NAME"
    else
        echo "[$NAME] submitting merge job"
        sbatch --job-name="merge_${NAME}" "$SLURM_SCRIPT" "$NAME"
    fi
done
```

#### SCRIPT: `scripts/06_merge_assemblies.slurm` (compares one sample)

```bash
#!/bin/bash
#SBATCH --partition=serial
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/06_merge_assemblies_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/06_merge_assemblies_%j.err

# Compare metaviralSPAdes and metaSPAdes contigs for one sample and merge them.
# Same contig = >=95% ANI over >=85% of the shorter contig (blastn + CheckV anicalc).
#
# MV_ONLY=1 handles a sample with NO metaSPAdes assembly (SRR28004142, SRR28004143 - their
# metaSPAdes jobs were cancelled for time; see README "05b-fix-c"). Nothing can be compared,
# so the metaviralSPAdes contigs become that sample's final contig set unchanged. It has to
# be asked for explicitly, not inferred from a missing file, so that running this step
# before the metaSPAdes jobs finish can never silently throw an assembly away.

set -euo pipefail

THREADS="${SLURM_CPUS_PER_TASK:-4}"
MIN_LEN="${MIN_LEN:-1000}"   # minimum metaSPAdes contig length kept
MIN_ANI="${MIN_ANI:-95}" # Average nucleotide identity threshold for merging contigs
MIN_COV="${MIN_COV:-85}"
MV_ONLY="${MV_ONLY:-0}"

if (($# != 1)); then
    echo "Usage: sbatch $0 SAMPLE_NAME" >&2
    exit 2
fi

NAME="$1"
SCRIPT_DIR="/sciclone/home/iadanforth/K_brevis_viromes/scripts"
MV_FASTA="$HOME/scr10/data/kbrevis/contigs/${NAME}_contigs.fasta"
MS_FASTA="$HOME/scr10/data/kbrevis/contigs_no_v/${NAME}_contigs.fasta"
OUT_DIR="$HOME/scr10/data/kbrevis/contigs_merged"
WORK="$OUT_DIR/$NAME"

[[ -s "$MV_FASTA" ]] || { echo "Missing or empty: $MV_FASTA" >&2; exit 1; }

if ((MV_ONLY)); then
    # Refuse if a metaSPAdes assembly does exist: then a real merge is what is wanted, and
    # passing MV_ONLY=1 by mistake would quietly discard half the data.
    if [[ -s "$MS_FASTA" ]]; then
        echo "MV_ONLY=1 was given but $MS_FASTA exists." >&2
        echo "Drop MV_ONLY to merge both assemblies for $NAME." >&2
        exit 1
    fi
    echo "[$NAME] ==== metaviralSPAdes-only mode ===="
    echo "[$NAME] no metaSPAdes assembly; its contigs become the final set unchanged."
else
    [[ -s "$MS_FASTA" ]] || { echo "Missing or empty: $MS_FASTA" >&2; exit 1; }
fi

mkdir -p "$WORK"

module load miniforge3
eval "$(conda shell.bash hook)"
conda activate seqkit-env

# Tag contig names with their assembler so the two sets can be told apart. The MV_ prefix
# is applied in both modes so contig IDs mean the same thing across all samples.
seqkit replace -p '^' -r 'MV_' "$MV_FASTA" > "$WORK/mv.fasta"
if ((MV_ONLY)); then
    conda deactivate
    conda activate checkv-env
    # No --ms/--ani: merge_assemblies.py passes the contigs through as one cluster each and
    # labels every record assembly_source=metaviralspades_only. No length filter is applied,
    # which matches how metaviralSPAdes contigs are treated in the merged samples too
    # (MIN_LEN only ever applied to the metaSPAdes side).
    python -I "$SCRIPT_DIR/merge_assemblies.py" \
        --sample "$NAME" \
        --mv "$WORK/mv.fasta" \
        --out-dir "$OUT_DIR" | tee "$OUT_DIR/${NAME}_merge_summary.txt"
    exit 0
fi

seqkit seq -m "$MIN_LEN" "$MS_FASTA" | seqkit replace -p '^' -r 'MS_' > "$WORK/ms.fasta"

conda deactivate
conda activate checkv-env

echo "[$NAME] blastn: metaviralSPAdes contigs vs metaSPAdes contigs"
makeblastdb -in "$WORK/ms.fasta" -dbtype nucl -out "$WORK/ms_db" > /dev/null
blastn \
    -query "$WORK/mv.fasta" \
    -db "$WORK/ms_db" \
    -outfmt '6 std qlen slen' \
    -max_target_seqs 10000 \
    -num_threads "$THREADS" \
    -out "$WORK/mv_vs_ms.blast.tsv"

anicalc -i "$WORK/mv_vs_ms.blast.tsv" -o "$WORK/mv_vs_ms.ani.tsv"

python -I "$SCRIPT_DIR/merge_assemblies.py" \
    --sample "$NAME" \
    --mv "$WORK/mv.fasta" \
    --ms "$WORK/ms.fasta" \
    --ani "$WORK/mv_vs_ms.ani.tsv" \
    --out-dir "$OUT_DIR" \
    --min-ani "$MIN_ANI" \
    --min-cov "$MIN_COV" | tee "$OUT_DIR/${NAME}_merge_summary.txt"

rm -f "$WORK"/ms_db.*
```

#### SCRIPT: `scripts/merge_assemblies.py`

```python
#!/usr/bin/env python3
"""Merge metaviralSPAdes and metaSPAdes contigs from one sample.

Two contigs (one from each assembler) are treated as the same sequence when
they share >= MIN_ANI % average nucleotide identity over >= MIN_COV % of the
SHORTER contig (MIUViG-style 95% ANI / 85% AF by default). ANI and coverage
come from `anicalc` (CheckV) run on blastn output.

Clustering is greedy and centroid-based (like CheckV's `aniclust`): contigs are
ranked, and each unassigned contig in rank order becomes a representative and
absorbs every still-unassigned contig from the OTHER assembler that matches it.
Ranking:
    1. metaviralSPAdes circular contigs (likely complete genomes)
    2. longer contigs first
    3. on equal length, metaviralSPAdes before metaSPAdes
Contigs are only compared across assemblers; nothing is deduplicated within
one assembler's output.

metaviralSPAdes-only mode: omit --ms and --ani for a sample that has no metaSPAdes
assembly (SRR28004142 and SRR28004143, whose metaSPAdes jobs were cancelled; see
README "05b-fix-c"). The metaviralSPAdes contigs are then passed through unchanged, one
cluster each, and every record is labelled assembly_source=metaviralspades_only so it is
never mistaken for a two-assembler result. Nothing is compared or deduplicated, which is
correct: duplicates can only arise BETWEEN assemblers, and there is only one here.

Outputs (in OUT_DIR):
    SAMPLE_merged_contigs.fasta   one sequence per cluster (union of both)
    SAMPLE_shared_contigs.fasta   only clusters found by BOTH assemblers
                                  (empty in metaviralSPAdes-only mode)
    SAMPLE_clusters.tsv           cluster membership and match statistics
"""
import argparse
import sys


def read_fasta(path):
    seqs = {}
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def topology(name):
    # metaviralSPAdes headers end in _type_circular / _type_linear / _type_linearrepeat
    if "_type_" in name:
        return name.rsplit("_type_", 1)[1]
    return "NA"


def wrap(seq, width=80):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sample", required=True)
    p.add_argument("--mv", required=True, help="metaviralSPAdes contigs, headers prefixed MV_")
    p.add_argument("--ms", help="metaSPAdes contigs, headers prefixed MS_. "
                                "Omit (with --ani) for metaviralSPAdes-only mode.")
    p.add_argument("--ani", help="anicalc output (MV queries vs MS subjects)")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--min-ani", type=float, default=95.0)
    p.add_argument("--min-cov", type=float, default=85.0)
    a = p.parse_args()

    # --ms and --ani go together: with both, the two assemblies are merged; with neither,
    # the sample has no metaSPAdes assembly and its metaviralSPAdes contigs are passed
    # through unchanged. One without the other is always a mistake, so refuse it rather
    # than silently dropping a whole assembly.
    if bool(a.ms) != bool(a.ani):
        sys.exit("--ms and --ani must be given together (or both omitted for "
                 "metaviralSPAdes-only mode)")
    merged_mode = bool(a.ms)
    assembly_source = "merged" if merged_mode else "metaviralspades_only"

    mv = read_fasta(a.mv)
    ms = read_fasta(a.ms) if merged_mode else {}
    seqs = {**mv, **ms}
    if len(seqs) != len(mv) + len(ms):
        sys.exit("Duplicate contig names between assemblies; expected MV_/MS_ prefixes")
    if not seqs:
        sys.exit(f"No contigs read for {a.sample}")

    # Cross-assembler matches: edges[contig] -> {other: (ani, cov_of_shorter)}
    # In metaviralSPAdes-only mode there is nothing to compare against, so every contig
    # keeps an empty edge set and becomes its own cluster below.
    edges = {n: {} for n in seqs}
    if merged_mode:
        with open(a.ani) as fh:
            next(fh)  # header: qname tname num_alns pid qcov tcov
            for line in fh:
                q, t, _, pid, qcov, tcov = line.rstrip("\n").split("\t")
                if q not in mv or t not in ms:
                    continue
                pid, qcov, tcov = float(pid), float(qcov), float(tcov)
                cov = qcov if len(seqs[q]) <= len(seqs[t]) else tcov
                if pid >= a.min_ani and cov >= a.min_cov:
                    prev = edges[q].get(t)
                    if prev is None or (pid, cov) > prev:
                        edges[q][t] = edges[t][q] = (pid, cov)

    def rank(n):
        is_mv = n in mv
        return (not (is_mv and topology(n) == "circular"), -len(seqs[n]), not is_mv, n)

    assigned = set()
    clusters = []
    for rep in sorted(seqs, key=rank):
        if rep in assigned:
            continue
        assigned.add(rep)
        members = []
        for other, (pid, cov) in sorted(edges[rep].items(), key=lambda kv: rank(kv[0])):
            if other not in assigned:
                assigned.add(other)
                members.append((other, pid, cov))
        clusters.append((rep, members))

    s = a.sample
    out = f"{a.out_dir}/{s}"
    counts = {"both": 0, "metaviralspades": 0, "metaspades": 0}
    with open(f"{out}_merged_contigs.fasta", "w") as merged, \
         open(f"{out}_shared_contigs.fasta", "w") as shared, \
         open(f"{out}_clusters.tsv", "w") as tsv:
        tsv.write("\t".join([
            "contig_id", "representative", "rep_assembler", "rep_length",
            "mv_topology", "assembly_source", "found_by", "n_matched",
            "matched_contigs", "matched_ani", "matched_cov_shorter"]) + "\n")
        for rep, members in clusters:
            rep_is_mv = rep in mv
            if members:
                found_by = "both"
            else:
                found_by = "metaviralspades" if rep_is_mv else "metaspades"
            counts[found_by] += 1

            mv_names = [rep] if rep_is_mv else []
            mv_names += [m for m, _, _ in members if m in mv]
            topos = sorted({topology(n) for n in mv_names})
            mv_topo = ",".join(topos) if topos else "NA"

            cid = f"{s}_{rep}"
            # assembly_source distinguishes "metaSPAdes ran and did not find this contig"
            # (merged) from "metaSPAdes was never run for this sample"
            # (metaviralspades_only). Without it, found_by=metaviralspades would be
            # ambiguous between the two, which matters for SRR28004142 and SRR28004143.
            header = (f">{cid} sample={s} assembly_source={assembly_source} "
                      f"found_by={found_by} "
                      f"rep_assembler={'metaviralspades' if rep_is_mv else 'metaspades'} "
                      f"mv_topology={mv_topo} n_matched={len(members)}")
            rec = f"{header}\n{wrap(seqs[rep])}\n"
            merged.write(rec)
            if members:
                shared.write(rec)
            tsv.write("\t".join([
                cid, rep,
                "metaviralspades" if rep_is_mv else "metaspades",
                str(len(seqs[rep])), mv_topo, assembly_source, found_by, str(len(members)),
                ",".join(m for m, _, _ in members) or "NA",
                ",".join(f"{pid:.2f}" for _, pid, _ in members) or "NA",
                ",".join(f"{cov:.2f}" for _, _, cov in members) or "NA",
            ]) + "\n")

    if merged_mode:
        mv_shared = sum(1 for n in mv if edges[n])
        print(f"[{s}] metaviralSPAdes contigs: {len(mv)} "
              f"({mv_shared} matched a metaSPAdes contig)")
        print(f"[{s}] metaSPAdes contigs (after length filter): {len(ms)}")
        print(f"[{s}] merged clusters: {len(clusters)} "
              f"(both={counts['both']}, metaviralspades_only={counts['metaviralspades']}, "
              f"metaspades_only={counts['metaspades']})")
    else:
        print(f"[{s}] assembly_source=metaviralspades_only "
              f"(no metaSPAdes assembly for this sample)")
        print(f"[{s}] metaviralSPAdes contigs passed through: {len(mv)}")
        print(f"[{s}] clusters: {len(clusters)} (one per contig; nothing was compared "
              f"or deduplicated, and {s}_shared_contigs.fasta is empty by definition)")


if __name__ == "__main__":
    main()
```

**What it does, step by step (per sample).**
1. `seqkit-env`: **seqkit** is a toolkit for FASTA/FASTQ files. Here it prefixes contig
   names with `MV_` (metaviralSPAdes) or `MS_` (metaSPAdes) so the two sets can be told
   apart, and drops metaSPAdes contigs shorter than 1,000 bp (`MIN_LEN`). Short contigs are
   too short to classify or compare reliably, and metaSPAdes produces hundreds of thousands
   of them.
2. `checkv-env`: **BLAST** (`blastn`) finds every stretch of similar DNA sequence between
   each metaviralSPAdes contig and the metaSPAdes contigs (`makeblastdb` turns the metaSPAdes
   contigs into a searchable database). `-max_target_seqs 10000` keeps all hits rather than
   the default 500.
3. CheckV's **`anicalc`** helper turns the BLAST hits into one row per contig pair: average
   nucleotide identity (ANI) and the fraction of each contig covered by the alignment.
4. **`merge_assemblies.py`** decides which contigs are the same and writes the merged files.
5. Deletes the BLAST database; other working files stay in `contigs_merged/SAMPLE/`.

**When two contigs count as the same.** One from each assembler, sharing **≥ 95 % ANI over
≥ 85 % of the shorter contig**. That is the community standard for grouping viral sequences
into the same "species-level" viral OTU (vOTU): the MIUViG guidelines (Roux et al. 2019,
*Nature Biotechnology*), which are also CheckV `aniclust`'s defaults. Both can be changed
(`MIN_ANI`, `MIN_COV` environment variables, or `--min-ani`, `--min-cov`).

**How contigs are merged** (greedy clustering, like CheckV's `aniclust`):
1. Rank all contigs: metaviralSPAdes **circular** contigs first (likely complete genomes),
   then longest first; on a tie, metaviralSPAdes first.
2. Go down the list. Each contig not yet used becomes a cluster representative and absorbs
   every unused contig from the *other* assembler that matches it.
3. Write one sequence per cluster: the representative. Sequences are not stitched together;
   the best existing contig is kept and the others are recorded as matches. Contigs are
   only compared between assemblers, not within one.

**Output** in `~/scr10/data/kbrevis/contigs_merged/`:
- `SAMPLE_merged_contigs.fasta`: one contig per cluster, i.e. everything from both
  assemblers with cross-assembler duplicates removed. Header example:
  `>SRR28004141_MV_NODE_1_..._type_circular sample=SRR28004141 found_by=both rep_assembler=metaviralspades mv_topology=circular n_matched=2`
- `SAMPLE_shared_contigs.fasta`: only clusters found by **both** assemblers.
- `SAMPLE_clusters.tsv`: one row per cluster (representative, length, metaviralSPAdes
  topology, `found_by` = both / metaviralspades / metaspades, matched contigs, ANI, coverage).
- `SAMPLE_merge_summary.txt`: counts.

**Resources.** 16 cores, 64 GB, 12 h. The header says partition `serial` (kuro); pick the
cluster with `~/hpc/best_cluster.sh -c 16 -m 64 -H 2` and override `-p` on the command line
if needed.

**Run (after step 05b finishes).**
```
cd ~/K_brevis_viromes
( module purge; module load slurm/kuro; bash scripts/06_merge_assemblies.sh )
```
Samples without a metaSPAdes assembly yet, or already merged, are skipped, so it can be
re-run as more assemblies finish.

**Testing.** No metaSPAdes output existed yet, so the merge step was tested on dev (single
thread, small files) with a stand-in "metaSPAdes" file made of: SRR28004141's
metaviralSPAdes contigs with renamed headers (should all match), SRR28004137's contigs
(mostly shouldn't match), and a copy of SRR28004141's longest contig with 2 kb trimmed off
(should match). Result: all 30 SRR28004141 contigs were `found_by=both`; the trimmed copy
merged into its full-length circular contig (kept as representative); 154 unmatched contigs
were kept as `metaspades` only. SRR28004137's NODE_4 also merged: it's identical to
SRR28004141's NODE_4, so the same phage was assembled in both samples.

**Where these commands come from (asked 2026-11-13, prompt 29).** Three things in this
script look odd until you know their provenance, so to be explicit:

1. **`checkv-env` is activated for `anicalc`, not for CheckV.** No `checkv` command runs in
   this step. `anicalc` is a *separate executable installed by the CheckV conda package* -
   its file manifest puts exactly `bin/checkv`, `bin/anicalc` and `bin/aniclust` on the
   PATH - so it never shows up in `checkv --help`, which is why it is easy to conclude it
   does not exist. Activating `checkv-env` is simply how you get `anicalc`.
2. **The whole `makeblastdb` -> `blastn` -> `anicalc` chain is CheckV's own documented
   recipe**, from the section "Supporting code: Rapid genome clustering based on pairwise
   ANI" near the end of the CheckV README. Read it as the raw file
   (https://bitbucket.org/berkeleylab/checkv/raw/master/README.md); Bitbucket's rendered
   view needs JavaScript and returns an empty shell to command-line fetches. The flags here
   match that recipe exactly, including `-outfmt '6 std qlen slen'` and
   `-max_target_seqs 10000`, and the recipe is also the source of the 95% ANI / 85% AF
   thresholds, which it describes as "the MIUViG recommended-parameters".
3. **`blastn`/`makeblastdb` come from BLAST+, which is NOT a CheckV dependency.** CheckV
   depends on DIAMOND and HMMER, not BLAST+. BLAST+ is in `checkv-env` from earlier work, so
   this is a convenient coincidence rather than a designed dependency. `envs/checkv-env.yml`
   pins `blast=2.16.0`, so a rebuild from that file keeps it; a rebuild from CheckV's own
   dependencies alone would break this script with `makeblastdb: command not found`.

**A deliberate deviation from that recipe:** CheckV's version ends with `aniclust`, which
dereplicates *one* set of sequences into vOTUs. This step is doing something different -
matching contigs *between two assemblers* and keeping both sets, including contigs unique to
one - so `merge_assemblies.py` does the clustering itself using the same thresholds.
(`aniclust` is used nowhere in this pipeline; step 09 reimplements the same greedy centroid
algorithm for multi-contig genomes.)

**What this step does NOT do (common misreading, clarified 2026-11-12, prompt 21).** No
virus identification happens here. `anicalc` is a small helper script that ships alongside
CheckV and only converts BLAST output into identity/coverage numbers; it never uses CheckV's
reference genomes or HMM profiles, and CheckV itself (`checkv end_to_end`) is not run until
step 08, after GVClass has already picked out the giant virus bins, and only to estimate
their completeness. CheckV never decides whether something is viral in any case. The
`seqkit seq -m 1000` call drops metaSPAdes contigs under 1 kb to keep the all-vs-all BLAST
tractable; it is **not** applied to the metaviralSPAdes contigs, which are kept at any
length. So `SAMPLE_merged_contigs.fasta` is every contig in the sample — bacterial, algal
and viral — not a viral subset.

**Caveats.** Most `found_by=metaspades` contigs will be bacterial or algal, because
metaSPAdes assembles the whole metagenome. The merged file is "all contigs", not "all viral
contigs": run a virus identifier (geNomad or VirSorter2) on it before calling anything
viral.

# Removed step (was 8). Binning contigs into genomes (MetaBAT2) — one-sample test

**Status: all the removed binning-test step scripts were removed on 2026-11-12** (prompts 14 and 15 in
`PROMPT_LOG.md`). Binning for the real analysis is done per sample inside step 08
(`08_gvmags.slurm`). This section is kept as a record of the test, and the two test scripts
below are the only remaining copy of their code: they were never committed to git.
Removal command (run by Claude Code at my request):

```
rm scripts/07_binning_test.sh scripts/07_binning_test.slurm
```

The test's data are still in scratch (`~/scr10/data/kbrevis/binning_test/`, 5.2 GB:
copied contigs, reads and the 4 bins). They aren't used by any remaining script and can be
deleted, or left for scratch to purge after 90 days unused.

#### SCRIPT (removed): `scripts/07_binning_test.sh`

```bash
#!/bin/bash
set -euo pipefail

module load miniforge3
eval "$(conda shell.bash hook)"
conda activate binning-env

for tool in bbmap.sh samtools jgi_summarize_bam_contig_depths metabat2; do
    command -v "$tool" >/dev/null 2>&1 || { echo "$tool not found" >&2; exit 1; }
done

READ_DIR="$HOME/scr10/data/kbrevis/clean"
CONTIG_DIR="$HOME/scr10/data/kbrevis/binning_test/contigs"
OUT_DIR="$HOME/scr10/data/kbrevis/binning_test/bins"
THREADS="$SLURM_CPUS_PER_TASK"
mkdir -p "$OUT_DIR"

shopt -s nullglob
CONTIG_FILES=("$CONTIG_DIR"/*_contigs.fasta)
((${#CONTIG_FILES[@]} > 0)) || { echo "No *_contigs.fasta in $CONTIG_DIR" >&2; exit 1; }
 


for CONTIGS in "${CONTIG_FILES[@]}"; do
    NAME="${CONTIGS##*/}"
    NAME="${NAME%_contigs.fasta}"
    R1="$READ_DIR/${NAME}_1_CLEAN.fastq"
	R2="$READ_DIR/${NAME}_2_CLEAN.fastq"

	[[ -s "$R1" ]] || { echo "Missing or empty R1 reads: $R1" >&2; exit 1; }
	[[ -s "$R2" ]] || { echo "Missing or empty R2 reads: $R2" >&2; exit 1; }
	
    SAMPLE_OUT="$OUT_DIR/$NAME"
    BAM="$SAMPLE_OUT/${NAME}_SORTED.bam"
    DEPTH="$SAMPLE_OUT/${NAME}_depth.txt"
 
    # Skip samples with no contigs (common with tiny test files)
    if [[ ! -s "$CONTIGS" ]]; then
        echo "[$NAME] contig file is empty, skipping" >&2
        continue
    fi
 
 
    # Skip samples that already finished
    if [[ -f "$SAMPLE_OUT/.done" ]]; then
        echo "[$NAME] already binned, skipping"
        continue
    fi
 
    mkdir -p "$SAMPLE_OUT"
 
    # 5. Map reads back to this sample's contigs and sort the alignments
    #    nodisk=t builds the contig index in memory instead of writing a ref/ folder
    echo "[$NAME] mapping reads to contigs"
    bbmap.sh -Xmx50g \
    ref="$CONTIGS" nodisk=t \
    in1="$R1" in2="$R2" \
    out=stdout.sam \
    threads="$THREADS" \
    2> "$SAMPLE_OUT/${NAME}_bbmap.log" |
    samtools sort -@ "$THREADS" -o "$BAM" -
 
    samtools index "$BAM"
 
    # 6. Calculate per-contig read depth
    echo "[$NAME] calculating contig depths"
    jgi_summarize_bam_contig_depths --outputDepth "$DEPTH" "$BAM"
 
    # 7. Bin contigs (default minimum contig length is 2500 bp)
    echo "[$NAME] binning"
    metabat2 \
        -i "$CONTIGS" \
        -a "$DEPTH" \
        -o "$SAMPLE_OUT/${NAME}_bin" \
        -t "$THREADS" \
        --seed 1
 
    N_BINS=$(find "$SAMPLE_OUT" -maxdepth 1 -name "${NAME}_bin.*.fa" | wc -l)
    echo "[$NAME] produced $N_BINS bins"
 
    touch "$SAMPLE_OUT/.done"
done
```

#### SCRIPT (removed): `scripts/07_binning_test.slurm`

```bash
#!/bin/bash
#SBATCH --job-name=07_binning_test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/07_binning_test_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/07_binnning_test_%j.err

set -euo pipefail

bash ./scripts/07_binning_test.sh
```

#### Removed: `06_bin_contigs.sh/.slurm` and `07_bin_contigs.sh` (2026-11-12)

Removed by Claude Code on 2026-11-12 at my request (prompt 14 in `PROMPT_LOG.md`):

```
rm scripts/06_bin_contigs.sh scripts/06_bin_contigs.slurm scripts/07_bin_contigs.sh
```

**Why they were removed.** All three binned contigs into genomes, which step 08
(`08_gvmags.slurm`, steps 1–2) now does for every sample. Step 08 does it better for this
project in three ways:
1. It bins the **merged** contigs (step 06: metaviralSPAdes + metaSPAdes, duplicates
   removed). The old scripts binned only the metaviralSPAdes contigs, which are mostly
   already single small viral genomes.
2. It uses the tool versions from Booker et al. (BBMap 35.85, samtools 1.15.1,
   MetaBAT 2.15, in `gvmag-bin-env`), not the newer versions in `binning-env`.
3. It sets memory properly (`--mem=96G`; BBMap gets 80 % of it). `06_bin_contigs.slurm` had
   no `--mem`, so it would have got less memory than `bbmap.sh -Xmx50g` asks for.

**What they were.**
- `06_bin_contigs.sh`: the same commands as `07_binning_test.sh` (above), for every
  sample in `~/scr10/data/kbrevis/contigs/` (the metaviralSPAdes contigs), writing to
  `~/scr10/data/kbrevis/bins/`. `06_bin_contigs.slurm` ran it with 8 cores, 1 day and no
  `--mem`. Never run, and never committed to git, so the files are gone; the commands are
  identical to `07_binning_test.sh` apart from the two folder paths.
- `07_bin_contigs.sh`: an older version of the same script that read interleaved,
  bfc-corrected reads (`~/scr10/data/kbrevis/corrected/*_CORRECTED_interleaved.fastq`, a
  folder that no longer exists). It had a broken first line (`#/!/bash`). Its `.slurm` was
  already deleted. The committed version can be recovered with
  `git show HEAD:scripts/07_bin_contigs.sh`; later uncommitted edits were not kept.

`07_binning_test.sh/.slurm` were removed the same day (prompt 15); their code is printed
above.

**What binning is and why.** A large genome, such as a giant virus or a bacterium, usually
assembles into many separate contigs. **Binning** groups contigs that probably come from
the same genome, using two signals: (1) contigs from one genome have similar **coverage**
(read depth), because the genome is present at one abundance; (2) they have similar
**sequence composition** (tetranucleotide frequencies). A bin is a "metagenome-assembled
genome" (MAG). Booker et al. built giant virus MAGs (GVMAGs) this way.

**What each tool does.**
- **BBMap** (`bbmap.sh`) aligns the sample's clean reads back to that sample's contigs.
  `nodisk=t` builds the contig index in memory instead of writing a `ref/` folder.
  `-Xmx50g` gives Java 50 GB.
- **samtools** sorts and indexes the alignments (BAM file).
- **`jgi_summarize_bam_contig_depths`** (part of MetaBAT2) computes each contig's mean depth.
- **MetaBAT2** bins the contigs using depth + composition. Defaults: contigs shorter than
  2,500 bp are ignored, bins smaller than 200 kb are not written. `--seed 1` makes runs
  repeatable.
- Each finished sample gets a `.done` file so re-runs skip it.

**Test run.** I copied SRR28004137's metaviralSPAdes contigs to
`~/scr10/data/kbrevis/binning_test/contigs/SRR28004137_contigs.fasta` by hand and ran
`sbatch scripts/07_binning_test.slurm` from the repo root: job **611846**, bora `batch`, 8
cores / 96 GB, 10 min, COMPLETED. MetaBAT 2.12.1 formed **4 bins (3.56 Mb in total)** in
`binning_test/bins/SRR28004137/`.

**Caveats and differences from the paper.**
- Binning **metaviralSPAdes** contigs mixes up the logic: those contigs are already
  (mostly) individual small viral genomes, and the 200 kb minimum bin size means a bin must
  hold many contigs. Booker et al. binned **metaSPAdes** contigs. The planned approach is to
  bin the merged contigs from step 06 (that is what step 08 does).
- `binning-env` has **MetaBAT 2.12.1** and BBMap 40.05; the paper used MetaBAT 2.15 and
  BBMap 35.85. Step 08 uses a new environment pinned to the paper's versions.
- The all-sample binning scripts (`06_bin_contigs.*`, `07_bin_contigs.sh`) were removed
  on 2026-11-12; binning for the real analysis happens in step 08 (see "Removed" above).

# Removed step (was 9). CheckV database (and a failed test)

**Status: both the removed CheckV-database step scripts were removed on 2026-11-12** (prompt 17 in `PROMPT_LOG.md`).
The CheckV database is now downloaded (if missing) by step 07 script,
`07_gvmag_setup.slurm`, and every script uses **one** copy:
`~/data10/databases/checkv_db/checkv-db-v1.5`, with CheckV 1.1.1 from `checkv-env`. This
section is kept as a record of how that database was first made. Commands run by Claude
Code at my request:

```
rm scripts/07_download_checkv_db.slurm scripts/08_test_checkv.slurm
rm -rf ~/scr10/databases/checkv_db                 # the original download (scratch copy)
rm -rf ~/data10/databases/checkv_db_for_v0.9.0     # the copy made for CheckV 0.9.0 by setup job 73832
rmdir ~/scr10/databases
```

All three copies were the same CheckV database, v1.5, and the same size (6,674,501 KB),
so nothing was lost. The kept copy is the one in data10, which is backed up; scratch is not.

#### SCRIPT (removed): `scripts/07_download_checkv_db.slurm`

```bash
#!/bin/bash
#SBATCH --job-name=download_checkV_db
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/download_checkV_db__%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/download_checkV_db__%j.err


module load miniforge3
eval "$(conda shell.bash hook)"
conda activate checkv-env

mkdir -p "$HOME/scr10/databases"
mkdir -p "$HOME/scr10/databases/checkv_db"

checkv download_database "$HOME/scr10/databases/checkv_db"
```

**What CheckV does.** CheckV estimates how complete a viral genome is by comparing its
genes to a database of complete viral genomes (or, for novel viruses, from conserved
protein families' typical genome lengths). It also finds host (cellular) genes at the ends
of contigs (prophage boundaries) and flags likely complete circular genomes. It assigns a
quality tier (complete / high / medium / low / not determined). **It does not decide
whether a contig is viral**; use a virus identifier for that.

**Run and result.** `sbatch scripts/07_download_checkv_db.slurm`: job **621410**, bora
`batch`, 2 min 50 s, COMPLETED. CheckV 1.1.1 downloaded **database v1.5** to
`~/scr10/databases/checkv_db/checkv-db-v1.5` and built its DIAMOND database (peak
memory 1.5 GB).

Because scratch is purged after 90 days unused, Claude copied it to backed-up data10 on
2026-11-12 (prompt 1):

```
mkdir -p /sciclone/data10/iadanforth/databases
rsync -a /sciclone/scr10/iadanforth/databases/checkv_db /sciclone/data10/iadanforth/databases/
diff -rq /sciclone/scr10/iadanforth/databases/checkv_db /sciclone/data10/iadanforth/databases/checkv_db   # no output = identical
```

Result: `~/data10/databases/checkv_db/checkv-db-v1.5` (6.4 GB). This is now the only copy,
and every script uses it. The scr10 original was deleted on 2026-11-12 (see the status note
above).

#### SCRIPT (removed): `scripts/08_test_checkv.slurm`

```bash
#!/bin/bash
#SBATCH --job-name=test_checkV
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/test_checkV_db__%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/test_checkV_db__%j.err


module load miniforge3
eval "$(conda shell.bash hook)"
conda activate seqkit-env

OUTPUT_DIR="$HOME/scr10/data/kbrevis/contigs/SRR34748663"
REF_DIR="$HOME/scr10/databases/checkv_db/checkv-db-v1.5"
FILE="$HOME/scr10/data/kbrevis/contigs/SRR34748663/before_rr.fasta"
seqkit seq -m 5000 "$FILE" > "$OUTPUT_DIR/before_rr.5kb.fasta"

conda deactivate seqkit-env

conda activate checkv-env

checkv end_to_end "$OUTPUT_DIR/before_rr.5kb.fasta" checkv_out \
    -t 16 \
    -d checkv_db/"$REF_DIR"
```

**What it was meant to do.** Take SRR34748663's metaviralSPAdes `before_rr.fasta` (contigs
before SPAdes' repeat resolution), keep contigs ≥ 5 kb with `seqkit seq -m 5000` (CheckV is
unreliable on short fragments), and run `checkv end_to_end`.

**Result: failed** (job **621411**, 8 s). The `.err` log says:
`Directory 'checkv_db//sciclone/home/iadanforth/scr10/databases/checkv_db/checkv-db-v1.5' does not exist.`
The fixes needed:
- `-d checkv_db/"$REF_DIR"` → `-d "$REF_DIR"` (and point `REF_DIR` at the data10 copy).
- `conda deactivate seqkit-env` → `conda deactivate` (deactivate takes no environment name).
- `checkv_out` is relative to the submit folder (the repo, in home); use
  `"$OUTPUT_DIR/checkv_out"` so output goes to scratch.
- `-t 16` but the job only has 8 cores; use `-t "$SLURM_CPUS_PER_TASK"`.

Never re-run; the script was removed on 2026-11-12 (prompt 17).

#### 06d. metaviralSPAdes-only fallback for SRR28004142 and SRR28004143

Added 2026-11-13 (prompt 33 in `PROMPT_LOG.md`). Section 05b-fix-c above explains that
those two samples have no metaSPAdes assembly, and warned that `06_merge_assemblies.sh`
would therefore skip them and drop them out of the pipeline from step 06 onward. This
change keeps them in: for those two samples the **metaviralSPAdes contigs become the final
contig set unchanged**, while all 18 other samples still get the full two-assembler merge
(the union of both, with contigs found by both counted once).

**Why nothing needs deduplicating in that mode.** The merge exists to stop a contig that
both assemblers recovered from being counted twice. Duplicates can therefore only arise
*between* the two assemblers; neither `merge_assemblies.py` nor anything upstream
deduplicates within a single assembler's output, because SPAdes already emits one contig per
assembly-graph path. With only one assembler's output in hand there is nothing to compare,
so the contigs pass straight through, one cluster each. No length filter is applied either,
which matches the merged samples: `MIN_LEN` (1000 bp) has only ever been applied to the
metaSPAdes side, never to metaviralSPAdes contigs.

**Three scripts changed.**

1. **`scripts/merge_assemblies.py`** — `--ms` and `--ani` became optional. Giving neither
selects metaviralSPAdes-only mode; giving one without the other is refused rather than
silently dropping an assembly. Putting the fallback here rather than writing a separate
shell path means both modes share one code path, so the FASTA headers and the
`_clusters.tsv` columns cannot drift apart between samples:

```python
    if bool(a.ms) != bool(a.ani):
        sys.exit("--ms and --ani must be given together (or both omitted for "
                 "metaviralSPAdes-only mode)")
    merged_mode = bool(a.ms)
    assembly_source = "merged" if merged_mode else "metaviralspades_only"
...
    ms = read_fasta(a.ms) if merged_mode else {}
...
    edges = {n: {} for n in seqs}
    if merged_mode:
        with open(a.ani) as fh:
            ...
```

In metaviralSPAdes-only mode every contig keeps an empty edge set, so the existing greedy
clustering loop makes each one its own cluster with no special-casing.

2. **A new `assembly_source` field**, in both the FASTA headers and `_clusters.tsv`, with
values `merged` or `metaviralspades_only`. This is the important part for interpreting the
results later. Without it, `found_by=metaviralspades` would mean two completely different
things: "metaSPAdes ran on this sample and did not recover this contig" for the 18 merged
samples, versus "metaSPAdes never ran" for these two. Those must not be conflated in any
downstream comparison of assemblers. All 20 samples carry the field so the schema is
uniform; nothing had been merged yet when this was added, so no existing output needed
regenerating.

3. **`scripts/06_merge_assemblies.slurm`** — an `MV_ONLY` switch:

```bash
MV_ONLY="${MV_ONLY:-0}"
...
if ((MV_ONLY)); then
    if [[ -s "$MS_FASTA" ]]; then
        echo "MV_ONLY=1 was given but $MS_FASTA exists." >&2
        echo "Drop MV_ONLY to merge both assemblies for $NAME." >&2
        exit 1
    fi
    ...
```

and a short branch that prefixes the contigs with `MV_` (as in the merged path, so contig
IDs mean the same thing in every sample) and calls `merge_assemblies.py` with no
`--ms`/`--ani`. `blastn`, `makeblastdb` and `anicalc` are skipped, since there is nothing to
align against.

4. **`scripts/06_merge_assemblies.sh`** — an explicit list rather than an inference:

```bash
MV_ONLY_SAMPLES="${MV_ONLY_SAMPLES-SRR28004142 SRR28004143}"
...
    if [[ ! -s "$MS" ]] && ((!MV_ONLY)); then
        echo "[$NAME] no metaSPAdes assembly yet ($MS), skipping"
        continue
    fi
```

**Why a named list instead of just "use MV only when the metaSPAdes file is absent".** The
automatic rule is tempting and would be correct today, but it fails unsafely: run step 06
while the step 05b jobs are still queued and every unfinished sample would quietly get a
metaviralSPAdes-only result that looks like a finished one. The named list makes the two
exceptions deliberate and visible, keeps the old protective skip for every other sample, and
the guard in the `.slurm` script means that if either sample is ever reassembled with
metaSPAdes the job refuses to run until the sample is taken off the list. Override the list
for a one-off run with `MV_ONLY_SAMPLES="" bash scripts/06_merge_assemblies.sh`.

**Tested** before any job was submitted.

*Test 1 — metaviralSPAdes-only mode on the real SRR28004143 contigs:*

```
[SRR28004143] assembly_source=metaviralspades_only (no metaSPAdes assembly for this sample)
[SRR28004143] metaviralSPAdes contigs passed through: 828
[SRR28004143] clusters: 828 (one per contig; nothing was compared or deduplicated, and
              SRR28004143_shared_contigs.fasta is empty by definition)
in: 828  out: 828          # no contigs lost
SRR28004143_shared_contigs.fasta -> 0 bytes   # correct: nothing was found by both
>SRR28004143_MV_NODE_1_length_3352221_cov_184.007226_cutoff_0_type_circular sample=SRR28004143 \
   assembly_source=metaviralspades_only found_by=metaviralspades \
   rep_assembler=metaviralspades mv_topology=circular n_matched=0
```

*Test 2 — regression test of the merged path on a synthetic 2-vs-2 dataset* (two MV contigs,
two MS contigs, a hand-written `anicalc` table with one pair at 98% ANI / 90% coverage and
one at 80% / 45%), run with a fixed seed so it is reproducible and needs no `blastn`:

```
[TEST] metaviralSPAdes contigs: 2 (1 matched a metaSPAdes contig)
[TEST] metaSPAdes contigs (after length filter): 2
[TEST] merged clusters: 3 (both=1, metaviralspades_only=1, metaspades_only=1)
```

That is the behaviour the merge is supposed to have: the 98%/90% pair collapsed into one
cluster (`found_by=both`, counted once, not twice), the 80%/45% pair correctly stayed apart
as below the 95%/85% thresholds, and the unmatched contig from each assembler was kept. All
three rows carried `assembly_source=merged`.

*Test 3 — the misuse guard:* `MV_ONLY=1 bash scripts/06_merge_assemblies.slurm SRR28004144`
(a sample that does have a metaSPAdes assembly) exited 1 with
"MV_ONLY=1 was given but ... exists. Drop MV_ONLY to merge both assemblies".

*Test 4 — the submitter, with `sbatch` stubbed out so nothing was queued:* 18 normal merge
jobs, 2 metaviralSPAdes-only jobs (`--export=ALL,MV_ONLY=1`), 0 skipped, 20 total.
`bash -n` passes on both shell scripts.

**Caveats.**

- SRR28004142 and SRR28004143 are now in the pipeline, but **not on equal footing with the
  other 18**. Their contigs come from one assembler, so for them the metaviralSPAdes
  chromosome-removal filter is the only thing that shaped the contig set, and there is no
  second-assembler corroboration. They are also two of the three most deeply sequenced
  libraries. Any comparison of per-sample contig counts, GVMAG yields or vOTU membership has
  to account for this, and `assembly_source` in the headers and `_clusters.tsv` is how to do
  it. It should be stated in any write-up.
- metaviralSPAdes reports far fewer contigs than metaSPAdes by design (364 and 828 for these
  two samples, against thousands for a metaSPAdes assembly), because it discards the
  cellular fraction of the graph. Expect these two samples to contribute proportionally
  fewer bins at step 08 for that reason alone, independently of their biology.
- `_shared_contigs.fasta` is empty for these samples by construction. Any downstream script
  that treats "found by both assemblers" as a confidence filter will select nothing from
  them, which is correct but is a silent zero, not an error.

# Step 07. One-time setup for the GVMAG pipeline (environments, GVClass, databases)

`07_gvmag_setup.slurm`, `08_gvmags.slurm`, `select_gvmags.py` and the four environment
files were written by Claude Code on 2026-11-12 (prompts 6–7) in a session that was
interrupted before it finished. The two missing pieces, `08_gvmags.sh` (the array
submitter) and `summarize_gvmags.py` (the summary table), were added in prompt 11. The
On 2026-11-12 (prompt 17) the scripts were renumbered (`09_gvmag_setup.slurm` →
`07_gvmag_setup.slurm`, `10_gvmags.sh/.slurm` → `08_gvmags.sh/.slurm`) and switched from
CheckV 0.9.0 to CheckV 1.1.1 (see below). Setup has run (job 73832); the per-sample jobs
have **not** been submitted. Order: (1) wait for step 06; (2) run `08_gvmags.sh`.

Commands run for the renumbering and the CheckV switch (by Claude Code, prompt 17):

```
mv scripts/09_gvmag_setup.slurm scripts/07_gvmag_setup.slurm
mv scripts/10_gvmags.sh         scripts/08_gvmags.sh
mv scripts/10_gvmags.slurm      scripts/08_gvmags.slurm
# inside the scripts (and this README): 10_gvmags -> 08_gvmags, 09_gvmag_setup -> 07_gvmag_setup,
# 11_gvmag_votus -> 09_gvmag_votus (job names, log names, usage comments)
sed -i -e 's/10_gvmags/08_gvmags/g' -e 's/09_gvmag_setup/07_gvmag_setup/g' -e 's/11_gvmag_votus/09_gvmag_votus/g' \
    scripts/07_gvmag_setup.slurm scripts/08_gvmags.sh scripts/08_gvmags.slurm scripts/select_gvmags.py scripts/summarize_gvmags.py
rm envs/gvmag-checkv-env.yml envs/gvmag-checkv-env.lock.yml
module load miniforge3; eval "$(conda shell.bash hook)"; conda env remove -y -n gvmag-checkv-env
# then edited by hand: setup step 03 (download db v1.5 into ~/data10/databases/checkv_db),
# 08_gvmags.slurm (CHECKV_DB path, `conda activate checkv-env`), 08_gvmags.sh (setup check)
bash scripts/08_gvmags.sh --dry-run    # test: setup check passed; stops at "No merged contigs yet" (expected)
```

The script listings below are copied from the current files.

**Goal.** Rebuild giant-virus MAGs from my merged contigs (step 06) with the same tools and
versions as Booker et al. 2026, Supplemental Methods ("Binning contigs"): bbmap 35.85 +
samtools 1.15.1 + MetaBAT 2.15 for binning, GVClass to find giant-virus bins, GVClass + CheckV
for completeness, pyrodigal-gv for gene prediction, and geNomad 1.7.4 for annotation
(bitscore ≥ 50).

**Deliberate deviation: CheckV version.** The paper used CheckV 0.9.0. I chose to use the
newer **CheckV 1.1.1** (`checkv-env`) with database **v1.5** for every script in this
project (prompt 17). CheckV's job is the same in both versions: it estimates how complete a
viral genome is by comparing its genes to a database of complete viral genomes, assigns a
quality tier (complete / high / medium / low), and counts genes that look cellular rather
than viral. The 1.x releases fixed bugs and work with the current database. Database v1.5
removed likely non-viral sequences from earlier versions, to reduce false matches to
cellular DNA. Completeness numbers can therefore differ slightly from the paper's, so
compare GVMAGs by their GVClass completeness (the paper's main measure) as well as CheckV.

#### Environment files (hand-written, version-pinned)

`envs/gvmag-bin-env.yml`:
```yaml
# Read mapping + binning, pinned to the versions in Booker et al. 2026 (Supplemental Methods)
name: gvmag-bin-env
channels: [conda-forge, bioconda]
dependencies:
  - bbmap=35.85
  - samtools=1.15.1
  - metabat2=2.15
  - seqkit
```

(`envs/gvmag-checkv-env.yml`, which pinned CheckV 0.9.0, was deleted on 2026-11-12 along
with its environment; CheckV now comes from the existing `checkv-env`, see
`envs/checkv-env.yml`.)

`envs/genomad-env.yml`:
```yaml
# Gene prediction (pyrodigal-gv; paper lists "v3.5.1" = the pyrodigal version) and
# gene annotation (geNomad 1.7.4), pinned to Booker et al. 2026 (Supplemental Methods)
name: genomad-env
channels: [conda-forge, bioconda]
dependencies:
  - genomad=1.7.4
  - pyrodigal=3.5.1
  - pyrodigal-gv
```

`envs/aniclustermap-env.yml`:
```yaml
# Pairwise ANI between GVMAGs, pinned to Booker et al. 2026 (Supplemental Methods)
name: aniclustermap-env
channels: [conda-forge, bioconda]
dependencies:
  - aniclustermap=1.1.0
```

#### SCRIPT: `scripts/07_gvmag_setup.slurm` (one-time setup)

```bash
#!/bin/bash
#SBATCH --job-name=07_gvmag_setup
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=6:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/07_gvmag_setup_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/07_gvmag_setup_%j.err

# One-time setup for the GVMAG pipeline (08_gvmags.*, 09_gvmag_votus.slurm):
#   1. conda environments (envs/*.yml): paper versions, except CheckV 1.1.1 (newer, by choice)
#   2. GVClass v1.1.1 (installed with pixi) + its reference database
#   3. CheckV database v1.5 (for CheckV 1.1.1 in checkv-env)
#   4. geNomad database
# Runs on a compute node because environment solves and downloads are too heavy for dev.
# Safe to re-run: finished steps are skipped.
#
# Submit (any CPU partition works, e.g. gulf highmem nodes):
#   module load slurm/gulf; sbatch -p batch -C highmem scripts/07_gvmag_setup.slurm

set -euo pipefail

REPO="/sciclone/home/iadanforth/K_brevis_viromes"
DB_DIR="$HOME/data10/databases"
PROG_DIR="$HOME/data10/programs"
GVCLASS_VERSION="v1.1.1"
GVCLASS_DIR="$PROG_DIR/gvclass-$GVCLASS_VERSION"

mkdir -p "$DB_DIR" "$PROG_DIR"

# Home has a 50 GB quota; the new conda envs need roughly 5 GB.
FREE_MB=$(quota -s 2>/dev/null | awk '/sciclone\/home/{getline; gsub("M","",$1); gsub("M","",$2); print $2-$1}')
if [[ -n "${FREE_MB:-}" && "$FREE_MB" -lt 8000 ]]; then
    echo "Less than 8 GB free in home quota ($FREE_MB MB). Free space first." >&2
    exit 1
fi

module load miniforge3
eval "$(conda shell.bash hook)"

# ---- 1. conda environments ----------------------------------------------------
for ENV in gvmag-bin-env checkv-env genomad-env aniclustermap-env; do
    if conda env list | awk '{print $1}' | grep -qx "$ENV"; then
        echo "[setup] $ENV already exists"
    else
        echo "[setup] creating $ENV from envs/$ENV.yml"
        conda env create -y -f "$REPO/envs/$ENV.yml"
    fi
    # Record the exact versions that were solved, for the README / reproducibility.
    conda env export -n "$ENV" --no-builds > "$REPO/envs/${ENV}.lock.yml"
done

# ---- 2. GVClass via pixi ------------------------------------------------------
# pixi and its package cache go to data10, not home (quota).
export PIXI_HOME="$PROG_DIR/pixi"
export PIXI_CACHE_DIR="$PROG_DIR/pixi-cache"
export PIXI_NO_PATH_UPDATE=1
if [[ ! -x "$PIXI_HOME/bin/pixi" ]]; then
    echo "[setup] installing pixi"
    curl -fsSL https://pixi.sh/install.sh | bash
fi
export PATH="$PIXI_HOME/bin:$PATH"
pixi --version

if [[ ! -d "$GVCLASS_DIR" ]]; then
    echo "[setup] cloning GVClass $GVCLASS_VERSION"
    git clone --depth 1 --branch "$GVCLASS_VERSION" \
        https://github.com/NeLLi-team/gvclass.git "$GVCLASS_DIR"
fi
cd "$GVCLASS_DIR"
pixi install
# Runs the bundled example; downloads the GVClass database (~700 MB) on first use
# and confirms the install works.
if [[ ! -f "$GVCLASS_DIR/.example_ok" ]]; then
    pixi run run-example
    touch "$GVCLASS_DIR/.example_ok"
fi
cd "$REPO"

# ---- 3. CheckV database v1.5 (used with CheckV 1.1.1 from checkv-env) ----------
# Pinned to v1.5: downloaded from its fixed URL instead of `checkv download_database`,
# which always fetches the newest release. The DIAMOND index is then built the same
# way `checkv download_database` does it (CheckV README, "manual install").
CHECKV_DB="$DB_DIR/checkv_db/checkv-db-v1.5"
if [[ ! -s "$CHECKV_DB/genome_db/checkv_reps.dmnd" ]]; then
    echo "[setup] downloading CheckV database v1.5"
    mkdir -p "$DB_DIR/checkv_db"
    curl -fL -o "$DB_DIR/checkv_db/checkv-db-v1.5.tar.gz" \
        https://portal.nersc.gov/CheckV/checkv-db-v1.5.tar.gz
    tar -xzf "$DB_DIR/checkv_db/checkv-db-v1.5.tar.gz" -C "$DB_DIR/checkv_db"
    rm "$DB_DIR/checkv_db/checkv-db-v1.5.tar.gz"
    conda activate checkv-env
    diamond makedb --threads "${SLURM_CPUS_PER_TASK:-1}" \
        --in "$CHECKV_DB/genome_db/checkv_reps.faa" \
        --db "$CHECKV_DB/genome_db/checkv_reps"
    conda deactivate
fi
echo "[setup] CheckV db: $CHECKV_DB"

# ---- 4. geNomad database ------------------------------------------------------
if [[ ! -d "$DB_DIR/genomad_db" ]]; then
    conda activate genomad-env
    genomad download-database "$DB_DIR"
    conda deactivate
fi
cat "$DB_DIR/genomad_db/version.txt" 2>/dev/null || true

echo "[setup] done"
```

**What it sets up.**
1. The conda environments `gvmag-bin-env`, `checkv-env` (only created if missing; it already
   existed), `genomad-env` and `aniclustermap-env`, plus a `*.lock.yml` export of exactly what
   was installed. Separate environments are needed because the paper's old versions (e.g.
   bbmap 35.85) conflict with the newer tools in my existing environments.
2. **GVClass v1.1.1.** GVClass classifies giant-virus genomes or bins: it finds giant-virus
   marker genes (GVOGs) and places them in reference phylogenies to give a taxonomy, and
   estimates completeness and contamination for NCLDV (*Nucleocytoviricota*) bins. It is
   installed with **pixi** (a package manager GVClass uses instead of conda), from GitHub,
   into `~/data10/programs/` to stay out of the home quota. The bundled example run also
   downloads the GVClass database (~700 MB). The paper used GVClass v0.9.4; v1.1.1 is the
   current release and the old one is not easily installable, so this is a deviation.
3. **The CheckV database v1.5**, at `~/data10/databases/checkv_db/checkv-db-v1.5`, only if
   it is missing. It is downloaded from its fixed address
   (`https://portal.nersc.gov/CheckV/checkv-db-v1.5.tar.gz`), not with
   `checkv download_database`, because that command always fetches the newest release and a
   rebuild should get the same version. The DIAMOND index (`checkv_reps.dmnd`, used by
   CheckV to search proteins quickly) is then built with `diamond makedb`, as
   `checkv download_database` would do. (Changed in prompt 17. Job 73832 ran the earlier
   version of this step, which used CheckV 0.9.0 and wrote to `checkv_db_for_v0.9.0/`; that
   copy has been deleted.)
4. **The geNomad database** in `~/data10/databases/genomad_db`.
5. **The GVMAGs V2 metadata table** (added in prompt 20) at
   `~/data10/databases/gvmagsV2/GVDB_stats.tsv`, ~5 MB: taxonomy, GVClass results and
   ecosystem for all 18,725 giant virus genomes of Vasquez et al. 2025
   (doi:10.1101/2025.09.26.678796, https://yvasquez-lbl.github.io/GVMAGsV2/). Step 09 uses
   it to describe each vOTU's nearest sequenced relative. The full 1.18 GB genome+protein
   bundle is **optional** and only fetched when the job is submitted with
   `--export=ALL,GVMAGS_V2_FULL=1`; it is needed only for `09_gvmag_votus.slurm --with-refs`.
   Note GVClass v1.1.1 already ships this same collection as its reference database (8,272
   NCLDV + 835 *Mirusviricota* genomes, counted in
   `~/data10/programs/gvclass-v1.1.1/resources/gvclassSeptember25_labels.tsv`), so the
   taxonomy in step 08 is already based on it; this download is only for the ANI work.
6. **A repair check for `aniclustermap-env`** (added in prompt 20): ANIclustermap 1.1.0
   needs matplotlib <3.9, and the first build of this environment solved matplotlib 3.10,
   which made `ANIclustermap` crash on import with
   `module 'matplotlib.cm' has no attribute 'register_cmap'`. The script now runs
   `ANIclustermap --help` and rebuilds the environment from the (now pinned) `.yml` if it
   fails.

Resources: 4 cores, 16 GB, 6 h. Needs internet access from the compute node (downloads
worked from bora in steps 01 and 03).

**Submitted 2026-11-12 (prompt 12).** The cluster was chosen with the helper script, which
asks each partition when a 4-core / 16 GB job would start (`sbatch --test-only`) and ranks
them by estimated finish time:

```
~/hpc/best_cluster.sh -c 4 -m 16 -H 2
# astral batch  start 12:23  finish ~15:15   (as00, shared node)     <- chosen
# femto  batch  start 12:23  finish ~16:01   (whole node reserved for one small job)
# bora   batch  start 12:23  finish ~16:50   (whole node reserved for one small job)
# bora   hima   start 19:14 ...  gulf batch 21:48 ...  kuro serial next day
# gulf   gz     not used (my instruction: no gz; jobs there can be cancelled by gz_preempt)
quota -s     # 24.6 GB of 50 GB used in home; the new conda envs need ~5 GB (script requires 8 GB free)
```

astral's `as00` could start immediately and is shared, so a 4-core job doesn't block a
whole node (femto and bora `batch` give you the entire node). `as01`, astral's GPU node, is
excluded because this job doesn't use GPUs. Walltime is the script's own 6 h rather than the
helper's 4 h, because download speeds are unpredictable. The script skips finished steps, so
if it runs out of time it can simply be resubmitted.

```
cd ~/K_brevis_viromes
( module purge; module load slurm/astral
  sbatch -p batch -N 1 -c 4 --exclude=as01 --mem=16G -t 6:00:00 scripts/07_gvmag_setup.slurm )
# Submitted batch job 73832
```

Job **73832** (astral). Logs: `logs/09_gvmag_setup_73832.out/.err` (written under the script's old name, `09_gvmag_setup.slurm`). Check with
`( module purge; module load slurm/astral; sacct -j 73832 -o jobid,state,elapsed,maxrss,exitcode )`.

----------------------------------------------------------------------------------------

# Step 08. Giant virus genomes (GVMAGs), following Booker et al. — COMPLETE (2026-11-16)

> All 20 samples finished: **848 bins, 26 GVMAGs in 9 of 20 samples**, almost all
> Imitervirales. See section 08h for the results, 08f for the `JAVA_HOME` fix and node
> choice, and 08g for the two timeouts and the bbmap fix.

Runs once per sample on the merged contigs from step 06, as a Slurm array submitted by
`08_gvmags.sh`. Needs step 07 to have finished first. Scripts: `08_gvmags.slurm`,
`08_gvmags.sh`, `select_gvmags.py`, `summarize_gvmags.py`.

#### SCRIPT: `scripts/08_gvmags.slurm` (one sample; submitted by `08_gvmags.sh`)

```bash
#!/bin/bash
#SBATCH --job-name=08_gvmags
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/08_gvmags_%A_%a.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/08_gvmags_%A_%a.err

# Build giant virus metagenome-assembled genomes (GVMAGs) for ONE sample, following
# Booker et al. 2026 (Supplemental Methods):
#   1. map the sample's reads to its merged contigs        (bbmap 35.85, samtools 1.15.1)
#   2. per-contig read depth + bin contigs                 (MetaBAT 2.15)
#   3. find the giant virus bins                           (GVClass)
#   4. quality/completeness of giant virus bins            (GVClass + CheckV 1.1.1, db v1.5)
#   5. gene prediction + annotation, bitscore >= 50        (geNomad 1.7.4, which predicts
#                                                           genes with pyrodigal-gv itself)
# Run as an array by 08_gvmags.sh (one task per line of the sample list),
# or for one sample:  sbatch 08_gvmags.slurm SAMPLE_NAME

set -euo pipefail

THREADS="${SLURM_CPUS_PER_TASK:?Submit with sbatch}"
REPO="/sciclone/home/iadanforth/K_brevis_viromes"
DATA="$HOME/scr10/data/kbrevis"
READ_DIR="$DATA/clean"
MERGED_DIR="$DATA/contigs_merged"
OUT_ROOT="$DATA/gvmags"
GVCLASS_DIR="$HOME/data10/programs/gvclass-v1.1.1"
GENOMAD_DB="$HOME/data10/databases/genomad_db"
CHECKV_DB="$HOME/data10/databases/checkv_db/checkv-db-v1.5"
MIN_BITSCORE=50

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${1:?sample list file}")
else
    NAME="${1:?Usage: sbatch $0 SAMPLE_NAME}"
fi

CONTIGS_IN="$MERGED_DIR/${NAME}_merged_contigs.fasta"
R1="$READ_DIR/${NAME}_1_CLEAN.fastq"
R2="$READ_DIR/${NAME}_2_CLEAN.fastq"
OUT="$OUT_ROOT/$NAME"
for f in "$CONTIGS_IN" "$R1" "$R2"; do
    [[ -s "$f" ]] || { echo "[$NAME] missing or empty: $f" >&2; exit 1; }
done
[[ -d "$GVCLASS_DIR" && -d "$GENOMAD_DB" && -d "$CHECKV_DB" ]] || {
    echo "Databases/GVClass not found; run scripts/07_gvmag_setup.slurm first" >&2; exit 1; }
if [[ -f "$OUT/.done" ]]; then
    echo "[$NAME] already finished, skipping"; exit 0
fi
mkdir -p "$OUT"/{binning,bins,gvclass_in,gvmags,checkv,genomad}

module load miniforge3
eval "$(conda shell.bash hook)"

# conda's per-package activation hooks are not written to survive `set -u`. openjdk's
# etc/conda/activate.d/java_home.sh starts with
#     export JAVA_HOME_CONDA_BACKUP=$JAVA_HOME
# which aborts immediately with "JAVA_HOME: unbound variable" when JAVA_HOME is not
# already set in the environment - and `sbatch` jobs start from a bare environment, so it
# never is. Job 110103 died this way after 2 seconds, before doing any work. gvmag-bin-env
# is affected because MetaBAT2/bbmap pull in openjdk (checkv-env and genomad-env are fine:
# their hooks use "${VAR:-}"). Relax `set -u` for the duration of the activate/deactivate
# calls only, rather than editing files inside the conda packages, which would be undone by
# any `conda env create` and would make the environment non-reproducible from its .yml.
conda_on()  { set +u; conda activate "$1"; set -u; }
conda_off() { set +u; conda deactivate;    set -u; }

# ---- 1-2. Map reads, compute depth, bin (gvmag-bin-env) -----------------------
conda_on gvmag-bin-env
# Contig names only (drop the "found_by=..." descriptions) so names match the BAM.
CONTIGS="$OUT/binning/${NAME}_contigs.fna"
seqkit seq -i "$CONTIGS_IN" > "$CONTIGS"

BAM="$OUT/binning/${NAME}.sorted.bam"
if [[ ! -s "$BAM.bai" ]]; then
    echo "[$NAME] mapping reads to merged contigs (bbmap)"
    bbmap.sh -Xmx$(( SLURM_MEM_PER_NODE * 8 / 10 / 1024 ))g \
        ref="$CONTIGS" nodisk=t \
        in1="$R1" in2="$R2" \
        out=stdout.sam threads="$THREADS" \
        2> "$OUT/binning/${NAME}_bbmap.log" |
        samtools sort -@ "$THREADS" -m 2G -o "$BAM" -
    samtools index "$BAM"
fi

echo "[$NAME] contig depths + MetaBAT2 binning"
jgi_summarize_bam_contig_depths --outputDepth "$OUT/binning/${NAME}_depth.txt" "$BAM"
# MetaBAT2 defaults as in the paper (contigs >= 2,500 bp, bins >= 200 kb); fixed seed
# so re-runs give identical bins.
metabat2 -i "$CONTIGS" -a "$OUT/binning/${NAME}_depth.txt" \
    -o "$OUT/bins/${NAME}_bin" -t "$THREADS" --seed 1
conda_off

# GVClass wants .fna files with no "." in the base name: SAMPLE_bin.3.fa -> SAMPLE_bin_3.fna
shopt -s nullglob
BINS=("$OUT"/bins/"${NAME}"_bin.*.fa)
echo "[$NAME] MetaBAT2 produced ${#BINS[@]} bins"
if ((${#BINS[@]} == 0)); then
    echo "[$NAME] no bins, nothing else to do"; touch "$OUT/.done"; exit 0
fi
for b in "${BINS[@]}"; do
    base=$(basename "$b" .fa); cp "$b" "$OUT/gvclass_in/${base//./_}.fna"
done

# ---- 3-4a. GVClass: taxonomy, completeness, contamination of every bin --------
echo "[$NAME] GVClass"
export PIXI_HOME="$HOME/data10/programs/pixi" PIXI_CACHE_DIR="$HOME/data10/programs/pixi-cache"
export PATH="$PIXI_HOME/bin:$PATH"
rm -rf "$OUT/gvclass"
( cd "$GVCLASS_DIR" && pixi run gvclass "$OUT/gvclass_in" -o "$OUT/gvclass" -t "$THREADS" )

# Giant virus bins = GVClass domain-level call is NCLDV (phylum Nucleocytoviricota),
# the classification shared by all GVMAGs in the paper (Supplemental Data 4).
conda_on checkv-env   # any env with python3 works for the helper
python -I "$REPO/scripts/select_gvmags.py" \
    --gvclass "$OUT/gvclass/gvclass_summary.tsv" \
    --bin-dir "$OUT/gvclass_in" \
    --out-dir "$OUT/gvmags" \
    --sample "$NAME"
conda_off

GVMAGS=("$OUT"/gvmags/*.fna)
echo "[$NAME] ${#GVMAGS[@]} giant virus bins (GVMAGs)"
if ((${#GVMAGS[@]} == 0)); then
    touch "$OUT/.done"; exit 0
fi

# ---- 4b. CheckV 1.1.1 with database v1.5 -------------------------------------
# CheckV scores single sequences, so each GVMAG's contigs are joined into one record
# (separated by 100 Ns) to get one completeness/quality value per GVMAG.
conda_on checkv-env
JOINED="$OUT/checkv/gvmags_joined.fna"
: > "$JOINED"
for g in "${GVMAGS[@]}"; do
    awk -v name="$(basename "$g" .fna)" '
        BEGIN { sp = sprintf("%100s", ""); gsub(/ /, "N", sp); printf ">%s\n", name }
        /^>/  { if (n++) printf "%s", sp; next }
              { printf "%s", $0 }
        END   { print "" }' "$g" >> "$JOINED"
done
checkv end_to_end "$JOINED" "$OUT/checkv/out" -t "$THREADS" -d "$CHECKV_DB" --remove_tmp
conda_off

# ---- 5. Gene prediction + annotation (geNomad) --------------------------------
# `genomad annotate` predicts the genes with pyrodigal-gv itself (the paper's gene caller)
# and writes them to <gvmag>_annotate/<gvmag>_proteins.faa, so pyrodigal-gv is NOT called
# separately here: that would predict the same genes twice with the same tool. --cleanup
# only deletes geNomad's mmseqs2 scratch directory, so the proteins are kept.
conda_on genomad-env
for g in "${GVMAGS[@]}"; do
    gname=$(basename "$g" .fna)
    genomad annotate --cleanup --threads "$THREADS" "$g" "$OUT/genomad/$gname" "$GENOMAD_DB"
    # Keep annotations with bitscore >= 50, as in the paper.
    GENES_TSV="$OUT/genomad/$gname/${gname}_annotate/${gname}_genes.tsv"
    awk -F'\t' -v min="$MIN_BITSCORE" \
        'NR==1{for(i=1;i<=NF;i++) if($i=="bitscore") c=i; print; next} c && $c!="NA" && $c>=min' \
        "$GENES_TSV" > "$OUT/genomad/${gname}_genes_bitscore${MIN_BITSCORE}.tsv"
done
conda_off

# ---- Summary table (like the paper's Supplemental Data 4) ---------------------
conda_on checkv-env
python -I "$REPO/scripts/summarize_gvmags.py" \
    --sample "$NAME" \
    --gvclass "$OUT/gvclass/gvclass_summary.tsv" \
    --gvmag-dir "$OUT/gvmags" \
    --checkv "$OUT/checkv/out/quality_summary.tsv" \
    --genomad-dir "$OUT/genomad" \
    --min-bitscore "$MIN_BITSCORE" \
    --out "$OUT/${NAME}_gvmag_summary.tsv"
conda_off

touch "$OUT/.done"
echo "[$NAME] done: $OUT/${NAME}_gvmag_summary.tsv"
```

#### SCRIPT: `scripts/select_gvmags.py`

```python
#!/usr/bin/env python3
"""Pick the giant virus bins (GVMAGs) out of one sample's MetaBAT2 bins.

A bin is called a GVMAG when GVClass's majority-rule taxonomy places it in the
domain NCLDV (phylum Nucleocytoviricota), i.e. `taxonomy_majority` starts with
"d_NCLDV". Every GVMAG reported by Booker et al. 2026 (Supplemental Data 4) has
this classification. No completeness or contamination cut-off is applied here:
the paper kept GVMAGs from 66 % completeness up; quality columns are reported in
the summary table instead so they can be filtered later.

Selected bins are copied to OUT_DIR as <query>.fna and listed in
OUT_DIR/gvmag_list.tsv.
"""
import argparse
import csv
import shutil
import sys
from pathlib import Path


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gvclass", required=True, help="GVClass gvclass_summary.tsv")
    p.add_argument("--bin-dir", required=True, help="directory of bin .fna files given to GVClass")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--sample", required=True)
    a = p.parse_args()

    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    n_bins = n_gv = 0
    with open(a.gvclass) as fh, open(out / "gvmag_list.tsv", "w") as lst:
        lst.write("sample\tgvmag\ttaxonomy_majority\torder_completeness\tcellular_unique\n")
        for row in csv.DictReader(fh, delimiter="\t"):
            n_bins += 1
            tax = row.get("taxonomy_majority", "")
            if not tax.startswith("d_NCLDV"):
                continue
            query = row["query"]
            src = Path(a.bin_dir) / f"{query}.fna"
            if not src.exists():
                sys.exit(f"GVClass query {query} has no matching bin file {src}")
            shutil.copy(src, out / f"{query}.fna")
            lst.write(f"{a.sample}\t{query}\t{tax}\t{row.get('order_completeness', 'NA')}\t"
                      f"{row.get('cellular_unique', 'NA')}\n")
            n_gv += 1
    print(f"[{a.sample}] GVClass classified {n_bins} bins; {n_gv} are NCLDV (GVMAGs)")


if __name__ == "__main__":
    main()
```

**What `08_gvmags.slurm` does, per sample.**
1. **Map + bin** (`gvmag-bin-env`): strip the `found_by=...` descriptions from the merged
   contig names (so names match the BAM), map the clean reads with bbmap, sort/index, compute
   depths, and bin with MetaBAT2 at its defaults (≥ 2,500 bp contigs, ≥ 200 kb bins, as in
   the paper), `--seed 1`. bbmap gets 80 % of the job's memory.
2. **GVClass** on every bin. Bins are renamed `SAMPLE_bin.3.fa` → `SAMPLE_bin_3.fna`
   because GVClass wants `.fna` files without extra dots.
**Why identification happens AFTER binning here (added 2026-11-12, prompt 23).** This order
is dictated by genome size, and it differs from a typical phage pipeline. A phage genome
(30-200 kb) usually assembles into a single contig, so the contig *is* the genome and you
identify contigs first (geNomad `end-to-end`, VirSorter2). A giant virus genome
(400 kb-2.5 Mb) shatters into many contigs - the paper's GVMAGs span 2 to 79 contigs each -
and any single fragment may carry no hallmark genes, so a contig-level classifier would miss
most of them. The genome has to be reassembled by binning before there are enough marker
genes to classify against. Hence: bin everything, then classify whole bins with GVClass.

Note GVClass is **not** giant-virus-only, which is what makes this work: its reference set
holds 27,333 bacterial, 2,094 archaeal, 1,800 eukaryotic and 4,193 phage genomes alongside
8,272 NCLDV (counted in `resources/gvclassSeptember25_labels.tsv`). Those cellular
references are there so it can say "this bin is bacterial", which is exactly what
`select_gvmags.py` relies on. VirSorter2 is not an alternative here: it scores contigs for
viralness (NCLDV is one of its classes) rather than assigning taxonomy and quality to
genome-sized bins.

3. **`select_gvmags.py`** keeps bins whose GVClass majority taxonomy starts with `d_NCLDV`
   (giant viruses, phylum *Nucleocytoviricota*), the classification shared by every GVMAG in
   the paper's Supplemental Data 4. **No completeness cut-off is applied here**; the paper
   reported GVMAGs from 66 % completeness up, so quality columns are kept to filter later.
4. **CheckV 1.1.1** (database v1.5) on each GVMAG. CheckV scores single sequences, so each GVMAG's contigs
   are joined into one record separated by 100 Ns to get one value per GVMAG. (This is
   Claude's choice; the paper doesn't say how multi-contig bins were given to CheckV.)
5. **geNomad annotate** both predicts the genes and annotates them. Gene prediction uses
   **pyrodigal-gv**, which geNomad calls internally: Prodigal (a gene finder for prokaryotes
   and viruses) extended with models for giant viruses and phages that use alternative
   genetic codes. (The paper writes "pyrodigal-gv v3.5.1"; 3.5.1 is a pyrodigal version.)
   A separate `pyrodigal-gv` call was **removed on 2026-11-12 (prompt 24)**: it predicted
   the same genes with the same tool a second time, and nothing read its output
   (`summarize_gvmags.py` takes its gene counts from geNomad's `_genes.tsv`). geNomad writes
   the proteins to `<gvmag>_annotate/<gvmag>_proteins.faa`, and `--cleanup` only deletes
   geNomad's mmseqs2 scratch directory, so the protein sequences are still kept — I
   confirmed this in the installed geNomad source (`modules/annotate.py`, "Run pyrodigal-gv"
   at line 170; cleanup removes only `annotate_mmseqs2_dir`). The `genes/` output folder is
   therefore gone.
**Expect most giant virus genes to stay unannotated (checked 2026-11-12, prompt 22).** In
the paper's own Supplemental Data 4, geNomad annotated only **595 of 5,912 genes across
their 11 GVMAGs, about 10%** (range 6.1-14.4% per GVMAG). That is normal: giant virus
genomes are full of genes with no known homolog. geNomad is kept here because it matches
the paper and is good at flagging virus **hallmark** genes, which is useful evidence that a
bin really is viral — but it will not tell you much about what these viruses do. Note also
that geNomad's main purpose is virus *identification* (`genomad end-to-end`), which this
pipeline does not use; only `annotate` is run, and gene *prediction* is done by
pyrodigal-gv, not geNomad. For functional biology later, consider eggNOG-mapper, or the
GVOG marker HMMs already installed with GVClass at
`~/data10/programs/gvclass-v1.1.1/resources/models/combined.hmm`.

6. **geNomad annotate** assigns functions to the genes using its marker database; rows
   with bitscore ≥ 50 are kept, as in the paper.
7. **`summarize_gvmags.py`** writes `gvmags/SAMPLE/SAMPLE_gvmag_summary.tsv`, one row per
   GVMAG (see below).

#### SCRIPT: `scripts/08_gvmags.sh` (submits the array; new, prompt 11)

```bash
#!/bin/bash
# Submit 08_gvmags.slurm as ONE Slurm array job: one task per sample that has merged
# contigs (step 06) and hasn't finished GVMAG building yet.
#
# Usage (from anywhere; load the Slurm module of the cluster you want first):
#   ( module purge; module load slurm/kuro; bash scripts/08_gvmags.sh -p serial )
#   bash scripts/08_gvmags.sh --dry-run          # show what would be submitted
# Anything after the options is passed straight to sbatch (partition, -t, ...).
# Pick the cluster with:  ~/hpc/best_cluster.sh -c 16 -m 96 -H 12
#
# Env vars: MAX_PARALLEL (default 10) = array tasks allowed to run at once.
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(dirname "$SCRIPT_DIR")"
DATA="$HOME/scr10/data/kbrevis"
MERGED_DIR="$DATA/contigs_merged"
READ_DIR="$DATA/clean"
OUT_ROOT="$DATA/gvmags"
SLURM_SCRIPT="$SCRIPT_DIR/08_gvmags.slurm"
MAX_PARALLEL="${MAX_PARALLEL:-10}"

DRY_RUN=0
if [[ "${1:-}" == "--dry-run" ]]; then
    DRY_RUN=1
    shift
fi

# ---- Check the one-time setup (07_gvmag_setup.slurm) was done -------------------
missing=()
[[ -d "$HOME/data10/programs/gvclass-v1.1.1" ]]  || missing+=("GVClass ($HOME/data10/programs/gvclass-v1.1.1)")
[[ -d "$HOME/data10/databases/genomad_db" ]]     || missing+=("geNomad database ($HOME/data10/databases/genomad_db)")
[[ -s "$HOME/data10/databases/checkv_db/checkv-db-v1.5/genome_db/checkv_reps.dmnd" ]] \
                                                 || missing+=("CheckV database v1.5 ($HOME/data10/databases/checkv_db/checkv-db-v1.5)")
if ((${#missing[@]})); then
    printf 'Setup not finished; missing:\n' >&2
    printf '  - %s\n' "${missing[@]}" >&2
    echo "Run scripts/07_gvmag_setup.slurm first." >&2
    ((DRY_RUN)) || exit 1
fi

# ---- Build the sample list -----------------------------------------------------
[[ -d "$MERGED_DIR" ]] || { echo "No merged contigs yet ($MERGED_DIR); run step 06 first." >&2; exit 1; }
mkdir -p "$OUT_ROOT" "$REPO_DIR/logs"

samples=()
shopt -s nullglob
for f in "$MERGED_DIR"/*_merged_contigs.fasta; do
    name="${f##*/}"
    name="${name%_merged_contigs.fasta}"
    if [[ -f "$OUT_ROOT/$name/.done" ]]; then
        echo "[$name] already finished, skipping"
        continue
    fi
    if [[ ! -s "$f" || ! -s "$READ_DIR/${name}_1_CLEAN.fastq" || ! -s "$READ_DIR/${name}_2_CLEAN.fastq" ]]; then
        echo "[$name] merged contigs or clean reads missing/empty, skipping" >&2
        continue
    fi
    samples+=("$name")
done

if ((${#samples[@]} == 0)); then
    echo "Nothing to submit."
    exit 0
fi

# A new, dated list for every submission: the running array reads it line by line,
# so it must not change while the array runs.
LIST="$OUT_ROOT/sample_list_$(date +%Y%m%d_%H%M%S).txt"
if ((DRY_RUN)); then
    LIST="(dry run: list not written)"
else
    printf '%s\n' "${samples[@]}" > "$LIST"
fi

N=${#samples[@]}
CMD=(sbatch --array="1-${N}%${MAX_PARALLEL}" "$@" "$SLURM_SCRIPT" "$LIST")
echo "${N} samples: ${samples[*]}"
echo "Sample list: $LIST"
if ((DRY_RUN)); then
    echo "Would run: ${CMD[*]}"
else
    "${CMD[@]}"
fi
```

**What it does.** Instead of 20 separate `sbatch` calls, it submits **one Slurm array
job** with one task per sample, which is easier on the scheduler and lets you cap how many
run at once. Each task runs `08_gvmags.slurm` on one sample.
1. **Checks that `07_gvmag_setup.slurm` has finished** (GVClass in
   `~/data10/programs/gvclass-v1.1.1`, the geNomad database in
   `~/data10/databases/genomad_db`, and the CheckV database at
   `~/data10/databases/checkv_db/checkv-db-v1.5`), and refuses to submit if not.
2. **Builds the sample list** from `~/scr10/data/kbrevis/contigs_merged/*_merged_contigs.fasta`
   (step 06 output). A sample is skipped if it already has `gvmags/SAMPLE/.done` (finished)
   or if its merged contigs or clean reads are missing or empty. So it can be re-run safely
   as more merges finish, and it only resubmits what's left.
3. Writes the list to a new, **dated** file, `gvmags/sample_list_YYYYMMDD_HHMMSS.txt`. Array
   task *N* reads line *N* of that file, so the file must not change while the array runs;
   a new name for every submission guarantees that.
4. Submits `sbatch --array=1-N%10 [your sbatch options] 08_gvmags.slurm LIST`.

**Parameters and why.**
- `%10` (`MAX_PARALLEL=10`): at most 10 samples run at once, per the cluster rules for
  arrays, so the job doesn't take over a partition. Change it with
  `MAX_PARALLEL=5 bash scripts/08_gvmags.sh ...`.
- `08_gvmags.slurm` has **no partition** in its header, on purpose: you choose the
  cluster when submitting. Anything after the options is passed to `sbatch` unchanged
  (e.g. `-p serial -t 1-00:00:00`).
- Each task requests 16 cores, 96 GB, 2 days (set in `08_gvmags.slurm`). bbmap gets 80 % of
  the memory. 96 GB is a guess for the largest samples (the deep CH1 libraries); after the
  first tasks finish, check `seff <jobid>_<task>` and lower it if possible.
- `--dry-run` prints the sample list and the `sbatch` command without writing or submitting
  anything (it still warns if setup isn't finished).

**Run (after step 06 and `07_gvmag_setup.slurm` have finished).**
```
cd ~/K_brevis_viromes
~/hpc/best_cluster.sh -c 16 -m 96 -H 12          # pick the cluster
bash scripts/08_gvmags.sh --dry-run -p serial     # check what would be submitted
( module purge; module load slurm/kuro; bash scripts/08_gvmags.sh -p serial )
```
Logs: `logs/08_gvmags_<arrayjob>_<task>.out/.err`.

#### SCRIPT: `scripts/summarize_gvmags.py` (one summary table per sample; new, prompt 11)

```python
#!/usr/bin/env python3
"""Summarise one sample's giant virus MAGs (GVMAGs) in one table.

One row per GVMAG, similar to Supplemental Data 4 of Booker et al. 2026:

  assembly statistics  (from the GVMAG .fna)    n_contigs, length_bp, gc_percent,
                                                n50_bp, longest_contig_bp
  CheckV 1.1.1         (quality_summary.tsv)    checkv_quality, checkv_completeness,
                                                checkv_completeness_method,
                                                checkv_contamination, checkv_viral_genes,
                                                checkv_host_genes, checkv_warnings
  geNomad annotate     (<gvmag>_genes.tsv)      genomad_genes, genomad_annotated
                                                (bitscore >= MIN_BITSCORE),
                                                genomad_virus_hallmarks
  GVClass              (gvclass_summary.tsv)    every GVClass column, prefixed
                                                "gvclass_" (copied as-is, so the table
                                                works whatever columns this GVClass
                                                version writes)

No GVMAG is filtered out here; filter on the quality columns afterwards.
Missing values are written as NA. Only the Python standard library is used.

Inputs are the files written by 08_gvmags.slurm:
  --gvmag-dir    OUT/gvmags/           <gvmag>.fna (one per GVMAG)
  --checkv       OUT/checkv/out/quality_summary.tsv  (contig_id = GVMAG name)
  --genomad-dir  OUT/genomad/          <gvmag>/<gvmag>_annotate/<gvmag>_genes.tsv
  --gvclass      OUT/gvclass/gvclass_summary.tsv     (query = GVMAG name)
"""
import argparse
import csv
import sys
from pathlib import Path

NA = "NA"


def read_fasta_lengths(path):
    """Return (list of contig lengths, total G+C count) for a FASTA file."""
    lengths, gc, cur = [], 0, None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    lengths.append(cur)
                cur = 0
            else:
                s = line.strip().upper()
                cur = (cur or 0) + len(s)
                gc += s.count("G") + s.count("C")
    if cur is not None:
        lengths.append(cur)
    return lengths, gc


def n50(lengths):
    total, run = sum(lengths), 0
    for length in sorted(lengths, reverse=True):
        run += length
        if run * 2 >= total:
            return length
    return 0


def read_tsv_by_key(path, key):
    """Read a TSV into {row[key]: row}. Missing file -> ({}, [])."""
    if not path or not Path(path).is_file():
        return {}, []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = {row[key]: row for row in reader if row.get(key)}
        return rows, list(reader.fieldnames or [])


def genomad_counts(genes_tsv, min_bitscore):
    """Count genes, annotated genes (bitscore >= min) and virus hallmark genes."""
    if not genes_tsv.is_file():
        return NA, NA, NA
    n = n_ann = n_hall = 0
    with open(genes_tsv, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            n += 1
            try:
                if float(row.get("bitscore", "NA")) >= min_bitscore:
                    n_ann += 1
            except ValueError:
                pass  # "NA" = no marker hit
            if row.get("virus_hallmark", "0") not in ("0", "", "NA"):
                n_hall += 1
    return n, n_ann, n_hall


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sample", required=True)
    p.add_argument("--gvclass", required=True, help="GVClass gvclass_summary.tsv")
    p.add_argument("--gvmag-dir", required=True, help="folder of <gvmag>.fna files")
    p.add_argument("--checkv", required=True, help="CheckV quality_summary.tsv")
    p.add_argument("--genomad-dir", required=True, help="folder of geNomad annotate outputs")
    p.add_argument("--min-bitscore", type=float, default=50.0)
    p.add_argument("--out", required=True, help="output TSV")
    a = p.parse_args()

    fnas = sorted(Path(a.gvmag_dir).glob("*.fna"))
    if not fnas:
        sys.exit(f"[{a.sample}] no .fna files in {a.gvmag_dir}")

    checkv, _ = read_tsv_by_key(a.checkv, "contig_id")
    gvclass, gvclass_cols = read_tsv_by_key(a.gvclass, "query")
    gvclass_cols = [c for c in gvclass_cols if c != "query"]
    if not checkv:
        print(f"[{a.sample}] warning: no CheckV results read from {a.checkv}", file=sys.stderr)
    if not gvclass:
        print(f"[{a.sample}] warning: no GVClass results read from {a.gvclass}", file=sys.stderr)

    checkv_cols = [("checkv_quality", "checkv_quality"),
                   ("checkv_completeness", "completeness"),
                   ("checkv_completeness_method", "completeness_method"),
                   ("checkv_contamination", "contamination"),
                   ("checkv_viral_genes", "viral_genes"),
                   ("checkv_host_genes", "host_genes"),
                   ("checkv_warnings", "warnings")]
    bs = int(a.min_bitscore) if a.min_bitscore.is_integer() else a.min_bitscore
    header = (["sample", "gvmag", "n_contigs", "length_bp", "gc_percent", "n50_bp",
               "longest_contig_bp"]
              + [c for c, _ in checkv_cols]
              + ["genomad_genes", f"genomad_annotated_bitscore{bs}", "genomad_virus_hallmarks"]
              + [f"gvclass_{c}" for c in gvclass_cols])

    out = Path(a.out)
    tmp = out.with_name(out.name + ".tmp")
    with open(tmp, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        for fna in fnas:
            name = fna.stem
            lengths, gc = read_fasta_lengths(fna)
            total = sum(lengths)
            cv = checkv.get(name, {})
            gv = gvclass.get(name, {})
            genes_tsv = Path(a.genomad_dir) / name / f"{name}_annotate" / f"{name}_genes.tsv"
            row = [a.sample, name, len(lengths), total,
                   f"{100 * gc / total:.2f}" if total else NA,
                   n50(lengths), max(lengths) if lengths else 0]
            row += [cv.get(col) or NA for _, col in checkv_cols]
            row += list(genomad_counts(genes_tsv, a.min_bitscore))
            row += [gv.get(c) or NA for c in gvclass_cols]
            w.writerow(row)
    tmp.replace(out)
    print(f"[{a.sample}] wrote {len(fnas)} GVMAG rows to {out}")


if __name__ == "__main__":
    main()
```

**What it does.** Collects everything known about each GVMAG into one tab-separated table,
`~/scr10/data/kbrevis/gvmags/SAMPLE/SAMPLE_gvmag_summary.tsv`, so the GVMAGs can be
compared with Booker et al.'s Supplemental Data 4. It is called at the end of
`08_gvmags.slurm` (in `checkv-env`, used only for its Python 3.12). It uses only the
Python standard library. Columns:

| column(s) | from | meaning |
|---|---|---|
| `sample`, `gvmag` | – | sample and bin name (e.g. `SRR28004137_bin_3`) |
| `n_contigs`, `length_bp`, `gc_percent`, `n50_bp`, `longest_contig_bp` | the GVMAG `.fna` | size and fragmentation of the genome. **N50** = the contig length at which half of the total length is in contigs that long or longer (higher = less fragmented). Giant virus genomes are usually 100 kb – 2.5 Mb and often have low GC. |
| `checkv_quality`, `checkv_completeness`, `checkv_completeness_method`, `checkv_contamination`, `checkv_viral_genes`, `checkv_host_genes`, `checkv_warnings` | CheckV 1.1.1 `quality_summary.tsv` | CheckV's estimate of how complete the genome is (%) and how it got that number, its quality tier, and how many genes look viral vs. cellular (`host_genes`). |
| `genomad_genes`, `genomad_annotated_bitscore50`, `genomad_virus_hallmarks` | geNomad `<gvmag>_genes.tsv` | number of genes geNomad predicted, how many matched a geNomad marker with **bitscore ≥ 50** (the paper's cut-off), and how many are **virus hallmark genes** (genes like major capsid protein that are specific to viruses; strong evidence the bin is viral). |
| `gvclass_*` | GVClass `gvclass_summary.tsv` | **all** GVClass columns, copied as-is with a `gvclass_` prefix: taxonomy at each rank, GVClass's completeness/contamination estimates, etc. |

Design choices:
- **Nothing is filtered out.** Every GVMAG chosen by `select_gvmags.py` gets a row, so you
  can apply the paper's thresholds (or your own) afterwards, e.g. in R.
- **GVClass columns are copied, not named.** GVClass v1.1.1 isn't installed yet, so its
  exact column names couldn't be checked. Copying them all means the table works whatever
  this version writes. (`select_gvmags.py` does assume `query` and `taxonomy_majority`
  exist; check that after the first real run.)
- Missing results (e.g. CheckV or geNomad failed for a GVMAG) become `NA` with a warning,
  not a crash, so one bad bin doesn't lose the whole table.
- The table is written to a temporary name and renamed when complete.

#### Testing (prompt 11)

Neither script could be run for real yet (no merged contigs and no GVClass/geNomad
install). Both were tested on dev, single-threaded, on small made-up inputs in the session
scratchpad:

- **`summarize_gvmags.py`**: two fake GVMAGs (`S1_bin_1`: 3 contigs of 8, 10 and 4 bp;
  `S1_bin_4`: 1 contig), a CheckV `quality_summary.tsv` and a geNomad `genes.tsv` in the
  real column layouts of CheckV 0.9 and geNomad 1.7 for `S1_bin_1` only, and a GVClass
  table that also has a bacterial bin. Result: 2 rows. `S1_bin_1`: length 22, GC 36.36 %,
  N50 8, CheckV `Medium-quality` / 62.5 %, 3 genes / 1 with bitscore ≥ 50 / 1 hallmark,
  all checked by hand; `S1_bin_4`: CheckV and geNomad columns `NA`; the bacterial bin
  left out (it isn't in `gvmags/`). With a missing CheckV file: a warning and `NA`s, exit 0.
- **`08_gvmags.sh`**: `bash -n` (syntax OK), then run with `HOME` pointed at a fake folder
  tree and a stand-in `sbatch` that only prints its arguments, with 4 samples: A and B
  ready, C already `.done`, D without clean reads. Setup missing → refused, exit 1.
  Setup missing + `--dry-run` → warned, listed A and B, printed the command. Setup present
  with `MAX_PARALLEL=5 ... -p serial -t 1-00:00:00` → wrote a list with A and B and ran
  `sbatch --array=1-2%5 -p serial -t 1-00:00:00 .../08_gvmags.slurm .../sample_list_<date>.txt`.
  C and D were skipped with messages.

Not tested: the real GVClass, CheckV and geNomad outputs, so check the first finished
sample's summary table by eye.

**Not yet written** (remaining step from the paper): GVMAG prevalence across the field
samples with CoverM. The vOTU and ANI steps are step 09 and the marker-gene phylogeny is
step 10, both below.

----------------------------------------------------------------------------------------

#### 08f. Running step 08 for real: the `JAVA_HOME` crash, the pilot, and the array (2026-11-13)

Written by Claude Code on 2026-11-13 (prompt 34 in `PROMPT_LOG.md`). Everything above in
step 08 was written but never executed; this is the first real run.

**What was run, in order.** Step 06 had just finished for all 20 samples, so step 08 became
runnable. Rather than launch all 20 array tasks at once, one sample was run alone first as a
pilot, so that any mistake would cost one job instead of twenty:

```bash
# pilot attempt 1 (FAILED after 2 seconds)
( module purge; module load slurm/kuro; sbatch -p serial -N 1 -c 16 --mem=96G -t 48:00:00 \
    -o /sciclone/home/iadanforth/K_brevis_viromes/logs/08_gvmags_pilot_%j.out \
    -e /sciclone/home/iadanforth/K_brevis_viromes/logs/08_gvmags_pilot_%j.err \
    scripts/08_gvmags.slurm SRR28004149 )
# Submitted batch job 110103

# after the fix below
rm -rf ~/scr10/data/kbrevis/gvmags/SRR28004149      # remove the half-built output dir
( module purge; module load slurm/kuro; sbatch -p serial -N 1 -c 16 --mem=96G -t 48:00:00 \
    -o /sciclone/home/iadanforth/K_brevis_viromes/logs/08_gvmags_pilot_%j.out \
    -e /sciclone/home/iadanforth/K_brevis_viromes/logs/08_gvmags_pilot_%j.err \
    scripts/08_gvmags.slurm SRR28004149 )
# Submitted batch job 110104   -> COMPLETED, 1:12:54, MaxRSS 26.74 GB

# the remaining 19 samples, as one array
( module purge; module load slurm/kuro; MAX_PARALLEL=6 bash scripts/08_gvmags.sh -p serial -t 12:00:00 )
# Submitted batch job 110120   (array 1-19%6)

# smoke test for the part the pilot could not reach (see "The pilot's blind spot" below)
( module purge; module load slurm/kuro; sbatch -p serial \
    <scratch>/tailtest.slurm )
# Submitted batch job 110122
```

**Bug 1: every task would have died in 2 seconds (`JAVA_HOME: unbound variable`).**
Job 110103 failed immediately with

```
/sciclone/home/iadanforth/.conda/envs/gvmag-bin-env/etc/conda/activate.d/java_home.sh: line 1: JAVA_HOME: unbound variable
```

The cause is a collision between two things that are each individually reasonable.
`08_gvmags.slurm` runs under `set -euo pipefail`, where `-u` aborts the script if it ever
reads a variable that was never set — a good safety net that catches typos in variable
names instead of silently using an empty string. But conda's *per-package activation
hooks* are not written to survive `-u`. When `gvmag-bin-env` is activated, conda sources
every file in `$CONDA_PREFIX/etc/conda/activate.d/`, and openjdk ships one whose first line
is:

```bash
export JAVA_HOME_CONDA_BACKUP=$JAVA_HOME
```

Its intent is to save any pre-existing `JAVA_HOME` so it can be restored on deactivation.
If `JAVA_HOME` was never set, that read is exactly what `-u` is designed to abort on. Jobs
submitted with `sbatch` do not source `~/.bashrc` and start from a nearly bare environment,
so `JAVA_HOME` is never set, and the failure is deterministic. `gvmag-bin-env` is the only
one of the five environments affected, because MetaBAT2 and bbmap are Java tools and pull
in openjdk; `checkv-env` and `genomad-env` activate cleanly (geNomad's libxml2 hook writes
`${VAR:-}`, the `-u`-safe form). This was latent rather than new: nothing in this project
had ever activated that environment under `set -u` before.

The fix adds two one-line helpers to `scripts/08_gvmags.slurm` and routes all ten
activate/deactivate calls through them, so `-u` is lifted only for the duration of the
activation itself and is back in force for the actual pipeline commands:

```bash
conda_on()  { set +u; conda activate "$1"; set -u; }
conda_off() { set +u; conda deactivate;    set -u; }
```

The alternative — editing `java_home.sh` inside the conda package to use `${JAVA_HOME:-}` —
was rejected deliberately. It would work today, but any `conda env create`/`conda install`
that rewrites the environment silently undoes it, and the environment would no longer be
reproducible from `envs/gvmag-bin-env.yml`, so anyone rebuilding this project from the repo
would hit the same crash with no clue why. Keeping the workaround in our own script makes it
visible and version-controlled. Verified before resubmitting, in a `set -euo pipefail`
subshell: the environment activates, `JAVA_HOME` is set to
`/sciclone/home/iadanforth/.conda/envs/gvmag-bin-env`, all five binning tools
(`bbmap.sh`, `samtools`, `metabat2`, `jgi_summarize_bam_contig_depths`, `seqkit`) are on
`PATH`, and deactivation returns cleanly while still under `-u`.

**Pilot result (job 110104, SRR28004149): the pipeline works, and the sample has no giant
viruses.** 1 h 12 min 54 s wall, 26.74 GB peak RSS of the 96 GB requested. bbmap mapped
40.7 % of the 42.5 M clean reads to that sample's merged contigs (11.5 M reads "well
mapped" per `jgi_summarize_bam_contig_depths`), MetaBAT2 produced 21 bins, and GVClass
classified all 21. GVClass itself took 68 m 50 s of the 73 min, i.e. this step is almost
entirely GVClass, and its runtime scales with the number of bins rather than the number of
contigs.

**Zero of the 21 bins are NCLDV, so SRR28004149 yields no GVMAGs.** Because this is a
negative result that silently skips the rest of the script, it was checked three ways
rather than taken at face value:

1. `select_gvmags.py` keys off `taxonomy_majority` beginning with `d_NCLDV`. That is the
   right column and the right string: GVClass's `src/core/summarize_full.py:92` comments
   *"For domain level, use the full domain name (e.g. "NCLDV" -> "d_NCLDV")"*. Note the
   separate `domain` column is **not** a plain label but a vote string such as
   `BAC:7(100.00%)`, so a selector written against `domain` would need to parse votes;
   ours does not use it.
2. Not one line of `gvclass_summary.tsv` contains the substring `NCLDV` at all, and no bin
   has a single `MIRUS` vote. There is therefore no borderline bin that a stricter or
   looser rule would have treated differently.
3. `mcp_total` is 0 for all 21 bins. The major capsid protein is the defining structural
   marker of Nucleocytoviricota, so its complete absence is the strongest available
   evidence that no giant virus genome is present.

The bins do carry a few giant-virus orthologous groups (`gvog4_unique` 0–4 of 4,
`gvog8_unique` 1–6 of 8), which at a glance looks like giant-virus signal. It is not: every
such bin also carries 16–70 unique *cellular* markers (`cellular_unique`), and the taxonomy
votes are 73–100 % bacterial. GVOGs are orthologous groups that are *enriched* in giant
viruses, not exclusive to them, so a few hits in a 0.5–3.5 Mb bacterial bin are the
expected background. Reading `gvog4`/`gvog8` as a GVMAG criterion would have produced a
large number of false giant viruses; GVClass's majority-rule taxonomy is what the paper
used and what is used here.

This is a biologically plausible outcome, not a failure. Giant virus MAGs are rare, the
paper recovered relatively few across all its samples, and a single sample yielding none is
unremarkable. It does mean the per-sample GVMAG count is a real result worth tabulating
once the array finishes.

**The pilot's blind spot, and the extra smoke test.** Because GVMAGs = 0, `08_gvmags.slurm`
took its early exit (`touch "$OUT/.done"; exit 0`) and therefore **never executed CheckV,
geNomad, the bitscore filter, or `summarize_gvmags.py`** — roughly the last 40 % of the
script, including several exact CheckV column names (`contig_id`, `completeness`,
`checkv_quality`, …) and geNomad's nested output path
`<gvmag>/<gvmag>_annotate/<gvmag>_genes.tsv`. A mistake in any of those would only surface
on the first sample that *does* have a GVMAG, after hours of mapping and GVClass.

So a separate one-off smoke test (job 110122, script kept in the session scratchpad, not in
the repo) runs exactly those commands using one existing MetaBAT2 bin
(`SRR28004149_bin_7`, 65 contigs, 212 kb — the smallest, so it is fast) as a stand-in
GVMAG. CheckV and geNomad do not inspect taxonomy, so a bacterial bin exercises the same
code paths: the 100-N joining `awk`, `checkv end_to_end`, `genomad annotate`, the bitscore
≥ 50 filter, and the summary table's joins against both CheckV and GVClass. It validates
the *plumbing* — commands, databases, file layout, column names, parsing — and explicitly
**not** the biology, since the stand-in is not a giant virus. It is run in parallel with
the array rather than before it, because the array's long pole (bbmap + GVClass) is already
proven and does not depend on the tail.

**Array submission (job 110120).** `08_gvmags.sh` built a dated sample list
(`~/scr10/data/kbrevis/gvmags/sample_list_20261113_190844.txt`), skipped SRR28004149 on its
`.done` marker, and submitted the other 19 samples as `--array=1-19%6`.

| setting | value | why |
|---|---|---|
| cluster/partition | kuro `serial` | see below |
| cores | 16 | as the pilot; GVClass runs 4 workers × 4 threads |
| memory | 96 GB | **not** reduced — see below |
| walltime | 12 h | ~10× the pilot's 73 min, leaving room for samples with more bins plus CheckV/geNomad; also short enough to backfill well. Kuro's maximum is 48 h |
| concurrency | 6 (`MAX_PARALLEL=6`) | 6 × 16 = 96 cores, ~19 % of the 512-core partition |

`~/hpc/best_cluster.sh -c 16 -m 64 -H 1.5` ranked gulf `gz` first and kuro `serial` second,
about 10 minutes apart in estimated finish. Kuro was chosen anyway: `gz` is preemptible
under the `normal` QOS, so a `gz_preempt` job **cancels** our task outright with no requeue
(AGENTS.md §3 rule 1). Each task represents up to ~2 h of mapping and GVClass, and across
19 tasks the expected loss from preemption far outweighs a 10-minute scheduling gain. Kuro
`serial` is also shared rather than `ExclusiveUser`, so a 16-core task does not idle a
whole 64-core node, and `MAX_PARALLEL` was lowered from the script's default of 10 to 6 to
leave most of a shared partition for other users.

**Memory was deliberately *not* shrunk**, even though the pilot peaked at 26.74 GB of 96 GB
and AGENTS.md §5 says to shrink over-requested jobs. bbmap's heap is derived from the
allocation (`-Xmx` = 80 % of `SLURM_MEM_PER_NODE`) and its reference index scales with the
assembly size, and SRR28004149 is a mid-sized sample: SRR28004148 (110,615 merged contigs),
SRR34748664 (100,610) and SRR28004145 (100,369) are considerably larger. Sizing 19 jobs
from one mid-sized sample's peak is precisely the mistake that killed the metaSPAdes job
92248 (see "Step 05b-fix"), and 96 GB of a 376 GB node still allows three tasks per node, so
there is nothing to be gained by trimming it. Re-check with `seff` once samples of all sizes
have finished.

**Caveats.**
- The array is running at the time of writing; only task 1 (SRR28004137, on ku60) had
  started. Nothing below step 08 has been verified against real GVMAGs yet.
- CheckV, geNomad and `summarize_gvmags.py` have been exercised only by the stand-in smoke
  test, never on a genuine NCLDV bin. The first sample that produces a GVMAG should still
  have its `SAMPLE_gvmag_summary.tsv` read by eye.
- The smoke-test script lives in the session scratchpad, so it will not survive. If that
  tail ever needs re-testing, it has to be rewritten; the commands in it are copies of the
  ones in `08_gvmags.slurm` and would drift if that script changes.
- A sample with 0 GVMAGs still gets a `.done` marker and **no** `_gvmag_summary.tsv`. Step
  09 must therefore tolerate missing summary files rather than assume one per sample.
- `conda_on`/`conda_off` only exist in `08_gvmags.slurm`. The other scripts
  (`09_gvmag_votus.slurm`, `10_gvmag_phylogeny.slurm`, `11_gvmag_tree.slurm`) use
  `aniclustermap-env`, `checkv-env`, `phylo-env` and `tree-env`, none of which depend on
  openjdk, so none is exposed today. If any of them ever gains a Java dependency, the same
  two helpers will be needed there.


----------------------------------------------------------------------------------------

#### 08g. Why tasks 110120_6 and _7 timed out, and the bbmap fix (2026-11-14)

Written by Claude Code on 2026-11-14 (prompt 35 in `PROMPT_LOG.md`).

Array tasks `110120_6` (SRR28004142) and `110120_7` (SRR28004143) both ended in `TIMEOUT`
at 2:00:19, killed by the 2 h limit set on them in section 08f. Neither reached MetaBAT2 or
GVClass: the last line of both logs is `mapping reads to merged contigs (bbmap)`.

**The mistake was treating contig count as a proxy for sample size.** Those two samples were
given a short limit because their merged contig files hold only 364 and 828 contigs, which
looked like minutes of work. But those files are small only because their metaSPAdes jobs
were cancelled for time (section 05b-fix-c), leaving them metaviralSPAdes-only (section
06d) - the *samples* are among the largest in the study:

| sample | clean R1 | merged contigs | note |
|---|---|---|---|
| SRR28004144 | 145.5 GB | 38,954 | largest read set (array task 8) |
| **SRR28004142** | **125.9 GB** | **364** | metaviralSPAdes-only; timed out |
| **SRR28004143** | **118.0 GB** | **828** | metaviralSPAdes-only; timed out |
| typical sample | ~13.4 GB | 17k-110k | |

bbmap's runtime is driven by the number of **reads**, not the size of the reference, so
these are roughly 9x the mapping work of a normal sample. A 2 h limit was never plausible.

**An inefficiency that made it much worse.** The mapping command piped bbmap's `out=`
(which emits *every* read, mapped or not) straight into `samtools sort`:

```bash
# before
bbmap.sh -Xmx$(( SLURM_MEM_PER_NODE * 8 / 10 / 1024 ))g \
    ref="$CONTIGS" nodisk=t in1="$R1" in2="$R2" \
    out=stdout.sam threads="$THREADS" 2> ...log |
    samtools sort -@ "$THREADS" -m 2G -o "$BAM" -
```

Unmapped reads have no reference position, so `jgi_summarize_bam_contig_depths` ignores
them and they cannot influence binning - but `samtools sort` still has to sort them. For a
sample whose reference is a few hundred viral contigs, almost nothing maps, so nearly the
entire 126 GB read set was being sorted to no purpose. The two killed tasks had already
written 32 and 112 sort temp fragments (13 GB and 44 GB, 57 GB total, since deleted).

The fix is `outm=` instead of `out=`, which makes bbmap emit mapped reads only. Verified
with a tiny single-threaded smoke test on dev (a 3 kb reference, 2 mappable and 2
unmappable read pairs): `out=` produced 8 SAM records, `outm=` produced 4. The depth table
and the bins are unchanged, because the discarded records contributed nothing to either.

The same edit fixes a latent **memory over-commit** in that command. bbmap's heap was 8/10
of `--mem` and `samtools sort` took `-m 2G` per thread, so with `--mem=96G` and 16 threads
the pipeline could ask for 76 GB + 32 GB = 108 GB inside a 96 GB allocation. It had not
bitten yet only because bbmap never used its full heap (the pilot peaked at 26.74 GB), but
a large sample could have hit it. Now 6/10 plus `-m 1G` = 57 GB + 16 GB = 73 GB.

```bash
# after
bbmap.sh -Xmx$(( SLURM_MEM_PER_NODE * 6 / 10 / 1024 ))g \
    ref="$CONTIGS" nodisk=t in1="$R1" in2="$R2" \
    outm=stdout.sam threads="$THREADS" 2> ...log |
    samtools sort -@ "$THREADS" -m 1G -o "$BAM" -
```

**Binning is still worth doing for both samples**, so they were resubmitted rather than
dropped. Despite the small contig counts they hold substantial assembled sequence:
SRR28004142 has 227 contigs >= 2,500 bp totalling 16.08 Mb (longest 1.24 Mb) and
SRR28004143 has 561 totalling 14.28 Mb (longest 3.35 Mb), so MetaBAT2's 2,500 bp contig
minimum and 200 kb bin minimum are both reachable. (A 3.35 Mb contig is larger than any
known giant virus genome and is far more likely a bacterial chromosome or a mis-assembly;
GVClass will decide.)

```bash
# orphaned sort fragments from the killed tasks
rm -f ~/scr10/data/kbrevis/gvmags/SRR2800414{2,3}/binning/*.tmp.*.bam   # 57 GB

# resubmitted as standalone jobs, 24 h instead of 2 h
( module purge; module load slurm/kuro; for s in SRR28004142 SRR28004143; do
    sbatch -p serial -N 1 -c 16 --mem=96G -t 24:00:00 \
      -o /sciclone/home/iadanforth/K_brevis_viromes/logs/08_gvmags_bigreads_%j.out \
      -e /sciclone/home/iadanforth/K_brevis_viromes/logs/08_gvmags_bigreads_%j.err \
      -J "08_gv_$s" scripts/08_gvmags.slurm "$s"
  done )
# Submitted batch job 110921   (SRR28004142)
# Submitted batch job 110922   (SRR28004143)

# keep the total footprint near the 6 slots agreed in 08f
( module purge; module load slurm/kuro; scontrol update JobId=110120 ArrayTaskThrottle=4 )
```

**Caveats.**
- **Slurm snapshots a batch script at submit time**, so the `outm` fix applies only to jobs
  110921/110922 and to future submissions. Array tasks 8-19 of job 110120 still run the old
  `out=` command and will sort their unmapped reads. They are normal-sized samples where
  ~40 % of reads map, so the waste is tolerable and they were left alone rather than
  cancelled and resubmitted.
- The BAMs produced from here on contain **mapped reads only**. That is equivalent for depth
  and binning, but the BAMs are no longer a complete record of the read set, so they cannot
  be reused to compute an overall mapping rate. The rate is still in each
  `<sample>_bbmap.log`.
- 24 h is an estimate, not a measurement: ~2 h of mapping at 9x a normal sample plus GVClass
  on an unknown number of bins. If these two time out again, raise to kuro's 48 h maximum.

----------------------------------------------------------------------------------------

#### 08h. Step 08 results: 848 bins, 26 GVMAGs (2026-11-16)

Written by Claude Code on 2026-11-16 (prompt 36 in `PROMPT_LOG.md`). All 20 samples are
finished and every output was checked; no job failed and no log contains an error.

**Jobs.** Pilot 110104 (SRR28004149), array 110120 tasks 1-19, and standalone 110921/110922
replacing the two tasks that timed out (section 08g). Elapsed times ran from 1 h 01 m to
10 h 47 m per sample. Smoke test 110122 - the stand-in GVMAG that exercised CheckV,
geNomad, the bitscore filter and `summarize_gvmags.py` - **passed** in 1 m 53 s, so the
part of the script the pilot could not reach is now known to work.

**Yield.**

| | count |
|---|---|
| samples | 20 |
| MetaBAT2 bins | 848 |
| GVMAGs (GVClass `d_NCLDV`) | **26** |
| samples with >= 1 GVMAG | **9** of 20 |
| GVMAGs per positive sample | 1-7 |

Eleven samples yielded no GVMAG at all. That is a real result, not a failure: giant virus
MAGs are rare, and only 3.1 % of bins (26/848) were NCLDV. The samples with the most GVMAGs
are SRR28004154 (7) and SRR28004137 (6); SRR28004144, SRR28004152 and SRR34748664 have one
each.

**Taxonomy.** 24 of the 26 are dominated by **Imitervirales**, the order that contains
Mimiviridae and Mesomimiviridae - the expected group for algal-associated giant viruses and
the same group that dominates Booker et al. The other two are
`SRR28004151_bin_1` (Pimascovirales) and `SRR28004139_bin_7` (Pandoravirales).
`SRR28004152_bin_36` is majority Algavirales.

**Quality (CheckV 1.1.1, db v1.5).**

| CheckV quality | n |
|---|---|
| High-quality | 5 |
| Medium-quality | 14 |
| Low-quality | 6 |
| Not-determined | 1 |

Sizes run from 203 kb to 2.19 Mb. The five high-quality GVMAGs are `SRR28004148_bin_1`
(2.19 Mb), `SRR28004154_bin_22` (1.23 Mb), `SRR28004137_bin_6` (530 kb),
`SRR34748664_bin_2` (348 kb) and `SRR28004151_bin_29` (203 kb). CheckV contamination is
0.00 % for all but two GVMAGs (`SRR28004151_bin_1` 1.23 %, `SRR28004148_bin_29` 0.49 %) and
no GVMAG carries a CheckV warning.

**Memory: the decision not to "right-size" was correct.** Section 08f kept `--mem=96G`
rather than shrinking to the pilot's 26.74 GB peak, on the argument that the pilot was
mid-sized and bbmap's index scales with the assembly. The finished jobs confirm it: peak RSS
reached the full 96 GB on SRR28004144 (task 8, the largest read set) and on both standalone
jobs, with 73.0, 71.1, 63.5 and 62.4 GB on other tasks. Shrinking to 32 or 48 GB, as the
pilot alone would have suggested, would have killed a large part of the run. Do **not**
reduce this request. Note that the 96 GB figures mean those jobs ran with no headroom at
all; if a future sample is larger, raise the allocation rather than hoping.

**Caveats.**
- Peak RSS equalling the allocation exactly is the cgroup ceiling, so the true demand of
  those jobs is unknown - it is at least 96 GB. They completed correctly, but the margin
  is zero.
- Eleven samples have a `.done` marker and **no** `_gvmag_summary.tsv`, because the script
  exits early when a sample has no GVMAG. Step 09 must tolerate missing summary files
  rather than assume one per sample.
- `SRR28004139_bin_7` has CheckV quality `Not-determined` and no completeness value; it is
  counted as a GVMAG because GVClass called it Pandoravirales, but nothing is known about
  how complete it is.
- Several GVMAGs are many small contigs (`SRR28004148_bin_1` is 463 contigs, `SRR28004137_bin_21`
  is 139). CheckV scored them as one joined sequence, so a high completeness value does not
  mean a contiguous genome.
- Bins whose GVClass vote is only partly NCLDV (e.g. `SRR28004148_bin_1` at 55.6 % NCLDV
  with 22 % bacterial votes) still pass the majority-rule test. Their bacterial contigs have
  not been removed; CheckV's contamination estimate for that GVMAG is nevertheless 0.00 %,
  which is worth treating sceptically for a 463-contig bin.

----------------------------------------------------------------------------------------

# Step 09. Grouping GVMAGs into viral OTUs (written 2026-11-12, not run yet)

Written by Claude Code on 2026-11-12 (prompt 20 in `PROMPT_LOG.md`). Files:
`scripts/09_gvmag_votus.slurm`, `scripts/votu_cluster.py`, `scripts/votu_nearest_ref.py`.
Nothing has been run yet: it needs the GVMAGs from step 08, which needs the merged contigs
from step 06.

#### 09a. What a vOTU is, and how it differs from binning and identification

These three operations are easy to confuse because they all "group sequences", but they
work at different levels and answer different questions.

| step | what is grouped | what comes out | the question it answers |
|---|---|---|---|
| **binning** (step 08, parts 1-2) | contigs **within one sample** | a bin = one draft genome (a MAG) | "which of these contigs are pieces of the same genome?" |
| **identification** (step 08, part 3) | nothing; each bin is labelled | a bin is a GVMAG, or it isn't | "is this draft genome a giant virus?" |
| **vOTU clustering** (step 09) | whole genomes **across all samples** | a vOTU = one viral "species" | "which of these genomes are the same virus?" |

1. **Binning** is assembly cleanup. A giant virus genome of 400 kb - 2.5 Mb almost never
   assembles into one piece, so MetaBAT2 gathers the fragments back together using read
   depth and sequence composition. It works inside a single sample, because coverage is
   only comparable within one sample. Output: `bin_1`, `bin_2`, ... Most are bacterial.
2. **Identification** asks what each bin is. GVClass looks for giant virus marker genes
   and places them in reference trees, and the bins that come back as NCLDV
   (*Nucleocytoviricota*) are the GVMAGs. This labels bins; it does not merge them.
3. **vOTU clustering** is the cross-sample step, and the one you were asking about. After
   steps 1-2 the same virus appears once per sample it was found in: if a virus infected
   the cultures in samples A, B and C, you have three near-identical GVMAGs that are three
   observations of one virus, not three viruses. Clustering genomes at >=95% ANI collapses
   them into one vOTU with one representative genome, so "we found six giant viruses" is a
   defensible statement and the same virus is not counted three times.

So your understanding is right, with one refinement: a vOTU is a cluster of **whole
genomes (bins)**, not of individual contigs. Contigs were already grouped into genomes by
binning in the previous step. The paper shows this clearly: its 11 GVMAGs collapsed into 7
vOTUs, because e.g. `VF_B.bins8` and `VF_C.bins2` (two samples) were both vOTU5, and
vOTU1-6 are the six giant viruses they report.

"vOTU" rather than "species" because viral species boundaries are not settled; the 95% ANI
cut-off is an operational convention from MIUViG (Roux et al. 2019), the same standard used
in step 06 to decide when the two assemblers had found the same contig.

#### 09a-bis. GVMAGs V2 bundle installed, and step 09 checked against missing summary files (2026-11-16)

Written by Claude Code on 2026-11-16 (prompt 37 in `PROMPT_LOG.md`). Two things were asked
for: install the GVMAGs V2 reference bundle so `--with-refs` works, and prove that step 09
copes with the 11 samples that have no `_gvmag_summary.tsv` because they produced no GVMAG.

**1. GVMAGs V2 bundle (job 92893).**

```bash
( module purge; module load slurm/gulf; sbatch -p gz -N 1 -c 4 --mem=16G -t 2:00:00 \
    --export=ALL,GVMAGS_V2_FULL=1 scripts/07_gvmag_setup.slurm )
# Submitted batch job 92893   -> COMPLETED in 2 m 25 s
```

Run on gulf `gz` although that partition is preemptible (AGENTS.md §3 rule 1), because the
job is short and fully restartable - every step of `07_gvmag_setup.slurm` is guarded and
skipped if already done, so a preemption costs only a resubmission. It was not preempted.

The bundle unpacks to `~/data10/databases/gvmagsV2/000_GVMAGS_V2/` containing `fna_files/`
(**8,508** reference genomes - the ANI-clustered representatives of Vasquez et al. 2025),
`faa_files/` and `GVMAGS_V2_STATS.tsv`. Total on disk is now 2.7 GB, in data10 because it
is a reference database that must survive scratch purges.

Three things were verified rather than assumed, because a silent mismatch here would make
`--with-refs` produce a table of `NA`s instead of failing:

- `09_gvmag_votus.slurm` finds the genomes with
  `find "$GVMAGS_V2_DIR" -maxdepth 3 -type d -name 'fna*'`, which matches
  `000_GVMAGS_V2/fna_files`. **Works as written, no change needed.**
- `votu_nearest_ref.py` joins hits to metadata on `meta[row["genome"]]`, keyed by the
  fastANI filename stem. 8,507 of the 8,508 `.fna` stems are present in the `genome` column
  of `GVDB_stats.tsv`. The one exception (`IMGM3300006193_BIN97`) falls through
  `meta.get(ref, {})` and would simply report `NA` metadata, so it cannot crash the step.
- All nine metadata columns the script reads (`Phylum`, `Class`, `Order`, `gvclass_family`,
  `new_genus`, `LENbp`, `GCperc`, `order_completeness`, `Ecosystem`) exist in
  `GVDB_stats.tsv`.

**A bug was found and fixed in `07_gvmag_setup.slurm`.** The download was guarded by
`! -d "$GVMAGS_V2_DIR/database"`, but the tarball never creates a `database/` directory, so
the test was always true and **every re-run with `GVMAGS_V2_FULL=1` would re-download the
full 1.18 GB**. The guard now tests for the same `fna*` directory that step 09 looks for,
and the script reports how many reference genomes are installed:

```bash
# before
if [[ "${GVMAGS_V2_FULL:-0}" == "1" && ! -d "$GVMAGS_V2_DIR/database" ]]; then
# after
if [[ "${GVMAGS_V2_FULL:-0}" == "1" && -z "$(find "$GVMAGS_V2_DIR" -maxdepth 3 -type d -name 'fna*' 2>/dev/null)" ]]; then
```

**2. A real contamination bug, caused by the step 08 smoke test.** Step 09 pooled GVMAGs
with `for sample_dir in "$GVMAG_ROOT"/*/`, i.e. *every* directory under `gvmags/`. The
smoke test of section 08f had left `~/scr10/data/kbrevis/gvmags/_tailtest/`, whose
`gvmags/SRR28004149_bin_7.fna` is a **bacterial** bin used only to exercise CheckV and
geNomad. Step 09 would have pooled it as a 27th "GVMAG" named
`_tailtest__SRR28004149_bin_7`, with no summary row, and clustered it alongside the real
ones - silently corrupting the vOTU table. The directory was moved out of the GVMAG root:

```bash
mv ~/scr10/data/kbrevis/gvmags/_tailtest ~/scr10/data/kbrevis/_tailtest_smoketest
```

and the pooling loop now refuses any directory that step 08 did not finish:

```bash
    if [[ ! -f "${sample_dir}.done" ]]; then
        echo "[votus] skipping $sample: no .done marker from step 08"
        continue
    fi
```

`.done` is the right test because `08_gvmags.slurm` writes it on every successful sample,
including the 11 that produced no GVMAG, so all 20 real samples still pass.

**3. Missing summary files: tested, not assumed.** Eleven of the 20 samples have a `.done`
marker and no `*_gvmag_summary.tsv`. Three facts make this safe, and all three were checked
by replaying the real step 09 logic against the real data:

- A sample with no GVMAG contributes **neither** a `.fna` **nor** a summary row, so the two
  counts stay consistent. Replaying the pooling and combining loops gave *26 pooled GVMAGs
  from 9 samples, 26 combined summary rows*, with the header written exactly once and
  `_tailtest` correctly skipped.
- `votu_cluster.py` only uses the summary for **ranking**: `completeness.get(n, -1.0)`
  sorts a genome with no row last, and the output reports its completeness as `NA`. The
  summary file is optional (`if a.summary and Path(a.summary).exists()`), so even a
  completely absent one is fine.
- Running `votu_cluster.py` on the 26 real GVMAGs with a synthetic fastANI table, under
  three scenarios - full 26-row summary, a summary with 3 rows deliberately removed, and no
  summary file at all - produced **25 vOTUs with identical cluster membership in all
  three**. Only the vOTU numbering changed. Rows with no completeness correctly showed
  `NA` (1, 4 and 26 rows respectively; the single `NA` in the full-summary case is
  `SRR28004139_bin_7`, whose CheckV quality is `Not-determined`).

The step now also prints both counts and warns, without failing, if they disagree:

```
[votus] combined summary rows: 26 (pooled GVMAGs: 26)
```

**Caveats.**
- Cluster *membership* is unaffected by missing completeness, but the **representative
  chosen within a multi-member vOTU can change**, because ranking is
  (completeness, length, name). With completeness unknown a genome ranks last, so it is
  less likely to be picked as the representative. For this dataset only one vOTU has more
  than one member, and both members have completeness, so nothing is affected today.
- The dry run used a **synthetic** fastANI table (one forced pair at 97.5 % ANI / 90 % AF)
  purely to exercise the code paths. The real clustering will differ; the 25 vOTUs above
  are not a result.
- `_tailtest_smoketest` is kept outside the GVMAG root as a record of the 08 validation.
  Anything else written under `~/scr10/data/kbrevis/gvmags/` must either be a real sample
  directory or will now be skipped for lack of a `.done` marker.
- The bundle's own `GVMAGS_V2_STATS.tsv` (2.2 MB, inside `000_GVMAGS_V2/`) is a different
  file from the top-level `GVDB_stats.tsv` (4.9 MB) that `votu_nearest_ref.py` reads. The
  script uses the latter, which is correct; do not swap them.

----------------------------------------------------------------------------------------

#### 09b. SCRIPT: `scripts/09_gvmag_votus.slurm`

Runs **once for the whole project**, not per sample, because its whole purpose is to
compare samples with each other. 16 cores, 64 GB, 12 h.

```bash
#!/bin/bash
#SBATCH --job-name=09_gvmag_votus
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/09_gvmag_votus_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/09_gvmag_votus_%j.err

# Group the GVMAGs from ALL samples into viral OTUs ("viral species"), following
# Booker et al. 2026 (Supplemental Methods, "GVMAG comparisons"):
#   1. pool every sample's GVMAGs from the removed binning-test step
#   2. all-vs-all ANI                               (fastANI, via ANIclustermap's engine)
#   3. cluster into vOTUs at >=95% ANI / >=85% AF    (MIUViG standard; votu_cluster.py)
#   4. clustered ANI heatmap                         (ANIclustermap 1.1.0, the paper's tool)
#   5. optional: nearest relative in GVMAGs V2       (Vasquez et al. 2025, 8,508 genomes)
#
# This runs ONCE for the whole project, after 08_gvmags.sh has finished for every sample.
# Submit (pick the cluster with ~/hpc/best_cluster.sh -c 16 -m 64 -H 2):
#   ( module purge; module load slurm/astral; sbatch -p batch --exclude=as01 scripts/09_gvmag_votus.slurm )
#   ...add --with-refs to also compare against the GVMAGs V2 reference genomes:
#   ( module purge; module load slurm/astral; sbatch -p batch --exclude=as01 scripts/09_gvmag_votus.slurm --with-refs )

set -euo pipefail

THREADS="${SLURM_CPUS_PER_TASK:?Submit with sbatch}"
REPO="/sciclone/home/iadanforth/K_brevis_viromes"
DATA="$HOME/scr10/data/kbrevis"
GVMAG_ROOT="$DATA/gvmags"
OUT="$DATA/votus"
GVMAGS_V2_DIR="$HOME/data10/databases/gvmagsV2"
MIN_ANI="${MIN_ANI:-95}"
MIN_AF="${MIN_AF:-85}"

WITH_REFS=0
[[ "${1:-}" == "--with-refs" ]] && WITH_REFS=1

[[ -d "$GVMAG_ROOT" ]] || { echo "No GVMAGs yet ($GVMAG_ROOT); run the removed binning-test step first." >&2; exit 1; }
mkdir -p "$OUT"/{gvmags,ani,heatmap,refs}

module load miniforge3
eval "$(conda shell.bash hook)"

# ---- 1. Pool every sample's GVMAGs, renamed <sample>__<bin> --------------------
# Bin names repeat between samples (every sample has a bin_1), so the sample name is
# prefixed to keep each GVMAG traceable to the sample it came from.
rm -f "$OUT"/gvmags/*.fna
shopt -s nullglob
n_samples=0
for sample_dir in "$GVMAG_ROOT"/*/; do
    sample=$(basename "$sample_dir")
    found=0
    for g in "$sample_dir"gvmags/*.fna; do
        cp "$g" "$OUT/gvmags/${sample}__$(basename "$g")"
        found=1
    done
    ((found)) && n_samples=$((n_samples + 1))
done
GVMAGS=("$OUT"/gvmags/*.fna)
echo "[votus] pooled ${#GVMAGS[@]} GVMAGs from $n_samples samples"
if ((${#GVMAGS[@]} == 0)); then
    echo "[votus] no GVMAGs were found in any sample; nothing to cluster."
    exit 0
fi

# One combined table of every sample's per-GVMAG summary, with the same <sample>__<bin>
# names, so completeness can be used to pick each vOTU's representative.
COMBINED="$OUT/all_gvmag_summaries.tsv"
: > "$COMBINED"
for s in "$GVMAG_ROOT"/*/*_gvmag_summary.tsv; do
    sample=$(basename "$(dirname "$s")")
    awk -F'\t' -v s="$sample" -v first="$([[ -s $COMBINED ]] && echo 0 || echo 1)" '
        NR==1 { if (first) print; next }
        { $2 = s "__" $2; print }' OFS='\t' "$s" >> "$COMBINED"
done
echo "[votus] combined summary rows: $(( $(wc -l < "$COMBINED") - 1 ))"

# ---- 2. All-vs-all ANI (fastANI) ----------------------------------------------
# fastANI compares whole genomes without aligning them base by base, and handles
# multi-contig genomes, which matters because a GVMAG is usually several contigs.
conda activate aniclustermap-env
printf '%s\n' "${GVMAGS[@]}" > "$OUT/ani/gvmag_list.txt"
fastANI --ql "$OUT/ani/gvmag_list.txt" --rl "$OUT/ani/gvmag_list.txt" \
    -t "$THREADS" -o "$OUT/ani/gvmags_all_vs_all.tsv"
echo "[votus] fastANI pairs: $(wc -l < "$OUT/ani/gvmags_all_vs_all.tsv")"

# ---- 4. Clustered ANI heatmap (the paper's Supplemental Figure 3) -------------
# Run before deactivating the env; failure here must not lose the vOTU results.
if ((${#GVMAGS[@]} >= 2)); then
    ANIclustermap -i "$OUT/gvmags" -o "$OUT/heatmap" \
        --fig_width 15 --fig_height 12 \
        || echo "[votus] WARNING: ANIclustermap failed; vOTU results below are unaffected" >&2
else
    echo "[votus] only one GVMAG: skipping the heatmap"
fi
conda deactivate

# ---- 3. Cluster into vOTUs ----------------------------------------------------
conda activate checkv-env   # used only for its Python 3.12
python -I "$REPO/scripts/votu_cluster.py" \
    --fastani "$OUT/ani/gvmags_all_vs_all.tsv" \
    --gvmag-dir "$OUT/gvmags" \
    --summary "$COMBINED" \
    --out-prefix "$OUT/gvmags" \
    --min-ani "$MIN_ANI" --min-af "$MIN_AF" | tee "$OUT/votu_summary.txt"
conda deactivate

# ---- 5. Optional: nearest relative among the GVMAGs V2 references -------------
# Each vOTU representative is compared to the 8,508 ANI-clustered reference genomes
# of Vasquez et al. 2025. This says how novel each vOTU is: a high ANI to a reference
# means it is a known giant virus, while no hit at all means fastANI found nothing
# above ~80% ANI, i.e. no close relative has been sequenced.
if ((WITH_REFS)); then
    REF_FNA_DIR=$(find "$GVMAGS_V2_DIR" -maxdepth 3 -type d -name 'fna*' | head -1)
    if [[ -z "$REF_FNA_DIR" ]]; then
        echo "[votus] --with-refs needs the full GVMAGs V2 bundle, which is not downloaded." >&2
        echo "[votus] Get it with: sbatch --export=ALL,GVMAGS_V2_FULL=1 ... scripts/07_gvmag_setup.slurm" >&2
    else
        conda activate aniclustermap-env
        find "$REF_FNA_DIR" -name '*.fna' > "$OUT/refs/ref_list.txt"
        echo "[votus] comparing vOTU representatives to $(wc -l < "$OUT/refs/ref_list.txt") references"
        fastANI --ql "$OUT/gvmags_representatives.txt" --rl "$OUT/refs/ref_list.txt" \
            -t "$THREADS" -o "$OUT/refs/votus_vs_gvmagsV2.tsv"
        conda deactivate
        # Best hit per representative, annotated from the GVMAGs V2 metadata table.
        conda activate checkv-env
        python -I "$REPO/scripts/votu_nearest_ref.py" \
            --fastani "$OUT/refs/votus_vs_gvmagsV2.tsv" \
            --votus "$OUT/gvmags_votus.tsv" \
            --metadata "$GVMAGS_V2_DIR/GVDB_stats.tsv" \
            --out "$OUT/votu_nearest_reference.tsv"
        conda deactivate
    fi
fi

echo "[votus] done. Results in $OUT:"
echo "  gvmags_votus.tsv            vOTU of every GVMAG"
echo "  gvmags_representatives.txt  one representative genome per vOTU"
echo "  votu_summary.txt            how many vOTUs were formed"
echo "  heatmap/                    clustered ANI heatmap"
((WITH_REFS)) && echo "  votu_nearest_reference.tsv  closest GVMAGs V2 relative per vOTU"
exit 0
```

**What each tool does and why.**
- **fastANI** computes average nucleotide identity between two whole genomes without
  aligning them base by base: it chops each genome into 3 kb fragments and maps them by
  shared k-mers. It is used instead of the BLAST + `anicalc` route from step 06 because it
  handles **multi-contig** genomes natively, which matters here since a GVMAG is typically
  2-80 contigs. It also reports how many fragments mapped, which gives the aligned
  fraction. Its practical floor is around 80% ANI; below that it reports nothing, which is
  why "no hit" in the reference table means "no close relative", not "not a virus".
- **ANIclustermap 1.1.0** is the tool the paper used (Supplemental Methods, "GVMAG
  comparisons"). It calls fastANI internally and draws the clustered heatmap with a
  dendrogram, the equivalent of the paper's Supplemental Figure 3. It is run for the figure
  only; the vOTU assignments do not depend on it, and the script warns and carries on if it
  fails.
- **`votu_cluster.py`** (below) does the actual clustering. Written by hand rather than
  using CheckV's `aniclust` because `aniclust` takes single sequences and an `anicalc`
  table, whereas these are multi-contig genomes compared with fastANI. It applies the same
  greedy centroid algorithm and the same default thresholds as `aniclust`.

**Thresholds, and where they come from.**
- **>=95% ANI** and **>=85% of the shorter genome aligned**: the MIUViG standard (Roux et
  al. 2019, *Nat Biotechnol*) for calling two viral genomes the same species-level vOTU,
  and also the defaults of CheckV's `aniclust` (`--min_ani 95`, `--min_tcov 85`), which I
  confirmed by running `aniclust --help`. Booker et al. say only that vOTUs were defined
  "based on their average nucleotide identity" and give no number, so this is my choice,
  not theirs. Override with `MIN_ANI=` / `MIN_AF=` as environment variables.
- **Representative choice:** highest CheckV completeness, then longest genome. The paper
  does not state how it picked one; this keeps the best-quality version of each virus.
- **The aligned fraction is measured on the shorter genome.** This stops a small fragment
  being absorbed into a large genome just because the fragment aligns well: a 200 kb bin
  matching 95% of a 1.2 Mb genome covers only 17% of the longer one, and taking the
  shorter genome's fraction is what `aniclust --min_tcov` does.

**Why greedy centroid clustering, and what it means.** GVMAGs are ranked best-first; each
unassigned genome in turn becomes a representative and absorbs the unassigned genomes that
match *it*. Absorbed members do not themselves recruit. This is deliberate (it is what
`aniclust` does) and avoids chaining, where A-B and B-C matches drag A and C into one
cluster even though A and C are not similar. The cost is that cluster membership depends on
the ranking, so a genome can end up as its own vOTU despite matching a member of another.

**Optional `--with-refs`.** Compares each vOTU representative to the 8,508 reference
genomes of GVMAGs V2 (Vasquez et al. 2025, see step 08) to find each vOTU's nearest
sequenced relative and the environment it came from. This needs the 1.18 GB bundle, which
`07_gvmag_setup.slurm` only downloads when submitted with
`--export=ALL,GVMAGS_V2_FULL=1`; without it the script prints how to get it and skips the
step rather than failing.

**Outputs** (in `~/scr10/data/kbrevis/votus/`):

| file | contents |
|---|---|
| `gvmags_votus.tsv` | one row per GVMAG: its vOTU, its `sample_set` (incubation or environmental, added 2026-11-24 - see 09i), whether it is the representative, and its ANI / aligned fraction to that representative |
| `gvmags_representatives.txt` | path to one representative genome per vOTU, the list to use for downstream work |
| `votu_summary.txt` | how many GVMAGs collapsed into how many vOTUs |
| `all_gvmag_summaries.tsv` | every sample's step-10 summary rows in one table, GVMAG names prefixed with the sample, plus a `sample_set` column (added 2026-11-24, see 09h) |
| `heatmap/` | the clustered ANI heatmap |
| `votu_nearest_reference.tsv` | with `--with-refs`: closest GVMAGs V2 relative per vOTU, with its taxonomy and ecosystem |

#### 09c. SCRIPT: `scripts/votu_cluster.py`

```python
#!/usr/bin/env python3
"""Group GVMAGs from all samples into viral operational taxonomic units (vOTUs).

A vOTU is the metagenomics stand-in for a viral "species": a cluster of genomes
similar enough to be treated as the same virus. The community standard (MIUViG,
Roux et al. 2019 Nat Biotechnol) is >=95% average nucleotide identity (ANI) over
>=85% of the shorter genome, which is also the default in CheckV's `aniclust`.
Booker et al. 2026 say only that they grouped GVMAGs into vOTUs "based on their
average nucleotide identity" without giving a cut-off, so the MIUViG values are
used here and can be changed with --min-ani / --min-af.

Clustering is greedy and centroid-based, like `aniclust`: GVMAGs are ranked
best-first, and each unassigned GVMAG in turn becomes a vOTU representative and
absorbs every unassigned GVMAG that matches it. Ranking is by CheckV completeness
(if a summary table is given), then genome length, then name, so the most complete
and longest genome represents each vOTU. That mirrors the paper, where one vOTU
(e.g. vOTU5) was recovered from two samples and reported once.

Inputs:
  --fastani   all-vs-all fastANI output (query, ref, ANI, mapped_frags, total_frags)
  --gvmag-dir directory of <sample>__<bin>.fna, one file per GVMAG
  --summary   optional combined per-GVMAG summary table (needs columns `gvmag`
              and `checkv_completeness`) used only for ranking

Outputs (written next to --out-prefix):
  <prefix>_votus.tsv          one row per GVMAG: its vOTU, whether it represents it,
                              and its ANI / alignment fraction to the representative
  <prefix>_representatives.txt  the representative GVMAGs' .fna paths, one per line
"""
import argparse
import csv
import sys
from pathlib import Path


def genome_length(path):
    return sum(len(line.strip()) for line in open(path) if not line.startswith(">"))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fastani", required=True)
    p.add_argument("--gvmag-dir", required=True)
    p.add_argument("--summary")
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--min-ani", type=float, default=95.0)
    p.add_argument("--min-af", type=float, default=85.0,
                   help="minimum aligned fraction of the SHORTER genome, in percent")
    a = p.parse_args()

    lengths = {f.stem: genome_length(f) for f in sorted(Path(a.gvmag_dir).glob("*.fna"))}
    if not lengths:
        sys.exit(f"No .fna files in {a.gvmag_dir}")

    completeness = {}
    if a.summary and Path(a.summary).exists():
        with open(a.summary) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                try:
                    completeness[row["gvmag"]] = float(row.get("checkv_completeness", ""))
                except (ValueError, TypeError):
                    pass

    # fastANI rows: query_path  ref_path  ANI  mapped_fragments  total_query_fragments
    # The aligned fraction is mapped/total, i.e. the fraction of the QUERY that aligned.
    # Running all-vs-all gives both directions, so for each pair take the fraction
    # belonging to the shorter genome, as aniclust's --min_tcov does.
    edges = {n: {} for n in lengths}
    with open(a.fastani) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 5:
                continue
            q, r = Path(f[0]).stem, Path(f[1]).stem
            if q == r or q not in lengths or r not in lengths:
                continue
            # Judge each pair on the row whose QUERY is the shorter genome, so `af`
            # is the aligned fraction of the shorter one. The reciprocal row is
            # skipped; all-vs-all guarantees the row we want exists.
            if lengths[q] > lengths[r]:
                continue
            ani, af = float(f[2]), 100.0 * float(f[3]) / float(f[4])
            if ani < a.min_ani or af < a.min_af:
                continue
            prev = edges[q].get(r)
            if prev is None or ani > prev[0]:
                edges[q][r] = edges[r][q] = (ani, af)

    def rank(n):
        return (-completeness.get(n, -1.0), -lengths[n], n)

    assigned, clusters = set(), []
    for rep in sorted(lengths, key=rank):
        if rep in assigned:
            continue
        assigned.add(rep)
        members = []
        for other in sorted(edges[rep], key=rank):
            if other not in assigned:
                assigned.add(other)
                members.append(other)
        clusters.append((rep, members))

    width = len(str(len(clusters)))
    with open(f"{a.out_prefix}_votus.tsv", "w") as tsv, \
         open(f"{a.out_prefix}_representatives.txt", "w") as reps:
        tsv.write("votu\tgvmag\tsample\tis_representative\trepresentative\t"
                  "length_bp\tcheckv_completeness\tani_to_rep\taf_to_rep\tn_members\n")
        for i, (rep, members) in enumerate(clusters, start=1):
            votu = f"vOTU{str(i).zfill(width)}"
            reps.write(f"{Path(a.gvmag_dir) / (rep + '.fna')}\n")
            for name in [rep] + members:
                ani, af = ("100.00", "100.00") if name == rep else \
                    tuple(f"{v:.2f}" for v in edges[rep][name])
                comp = completeness.get(name)
                tsv.write("\t".join([
                    votu, name, name.split("__")[0],
                    "yes" if name == rep else "no", rep,
                    str(lengths[name]), f"{comp:.2f}" if comp is not None else "NA",
                    ani, af, str(len(members) + 1)]) + "\n")

    multi = sum(1 for _, m in clusters if m)
    print(f"{len(lengths)} GVMAGs -> {len(clusters)} vOTUs "
          f"({multi} found in more than one sample/bin, "
          f"{len(clusters) - multi} singletons) "
          f"at >={a.min_ani}% ANI over >={a.min_af}% of the shorter genome")


if __name__ == "__main__":
    main()
```

#### 09d. SCRIPT: `scripts/votu_nearest_ref.py`

```python
#!/usr/bin/env python3
"""For each vOTU representative, report its closest relative in GVMAGs V2.

GVMAGs V2 (Vasquez et al. 2025, doi:10.1101/2025.09.26.678796) is a collection of
18,725 giant virus genomes reduced to 8,508 clustered representatives. Comparing our
vOTU representatives to it answers "is this a virus someone has already sequenced,
and where was it found?".

Reading the output:
  * ANI >= 95% over most of the genome  -> the same viral species as a known genome
  * ANI 80-95%                          -> a relative of a known genome, new species
  * no row at all                       -> fastANI found nothing above ~80% ANI, so no
                                           close relative has been sequenced yet. This
                                           is NOT evidence the GVMAG is wrong; most
                                           environmental giant viruses have no close
                                           relative in any database.

Inputs:
  --fastani   fastANI output, query = vOTU representatives, ref = GVMAGs V2 genomes
  --votus     votu_cluster.py's <prefix>_votus.tsv
  --metadata  GVDB_stats.tsv from the GVMAGs V2 release
"""
import argparse
import csv
from pathlib import Path

# Columns copied from GVDB_stats.tsv, as (output name, column in the file).
META_COLS = [("ref_phylum", "Phylum"), ("ref_class", "Class"), ("ref_order", "Order"),
             ("ref_family", "gvclass_family"), ("ref_genus", "new_genus"),
             ("ref_length_bp", "LENbp"), ("ref_gc_percent", "GCperc"),
             ("ref_order_completeness", "order_completeness"),
             ("ref_ecosystem", "Ecosystem")]
NA = "NA"


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fastani", required=True)
    p.add_argument("--votus", required=True)
    p.add_argument("--metadata", required=True)
    p.add_argument("--out", required=True)
    a = p.parse_args()

    # Best hit per query: highest ANI, then largest aligned fraction.
    best = {}
    with open(a.fastani) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 5:
                continue
            q, r = Path(f[0]).stem, Path(f[1]).stem
            ani, af = float(f[2]), 100.0 * float(f[3]) / float(f[4])
            if q not in best or (ani, af) > (best[q][1], best[q][2]):
                best[q] = (r, ani, af)

    meta = {}
    with open(a.metadata) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            meta[row["genome"]] = row

    with open(a.votus) as fh:
        reps = [r for r in csv.DictReader(fh, delimiter="\t")
                if r["is_representative"] == "yes"]

    hits = 0
    with open(a.out, "w") as out:
        out.write("\t".join(["votu", "representative", "n_members", "length_bp",
                             "checkv_completeness", "nearest_reference",
                             "ani_percent", "aligned_fraction_percent"]
                            + [c for c, _ in META_COLS]) + "\n")
        for r in reps:
            name = r["gvmag"]
            row = [r["votu"], name, r["n_members"], r["length_bp"],
                   r["checkv_completeness"]]
            if name in best:
                ref, ani, af = best[name]
                hits += 1
                row += [ref, f"{ani:.2f}", f"{af:.2f}"]
                m = meta.get(ref, {})
                row += [(m.get(col) or NA) for _, col in META_COLS]
            else:
                row += ["none_above_fastani_threshold", NA, NA] + [NA] * len(META_COLS)
            out.write("\t".join(row) + "\n")

    print(f"{len(reps)} vOTU representatives: {hits} with a GVMAGs V2 relative, "
          f"{len(reps) - hits} with none (no close relative sequenced)")


if __name__ == "__main__":
    main()
```

#### 09e. How these were tested

Neither real GVMAGs nor real fastANI output exist yet, so both helpers were tested on dev
with hand-made input (single-threaded and tiny, so safe on the front end).

`votu_cluster.py`: four mock GVMAGs, with fastANI rows written by hand. `SRR01__bin_1`
(500 kb, 71% complete) and `SRR02__bin_1` (480 kb, 95.5% complete) at 99.1% ANI;
`SRR02__bin_7`, a 200 kb fragment; and `SRR03__bin_2` at 82% ANI. Result:

```
4 GVMAGs -> 3 vOTUs (1 found in more than one sample/bin, 2 singletons) at >=95.0% ANI over >=85.0% of the shorter genome
votu   gvmag         sample  is_representative  representative  length_bp  checkv_completeness  ani_to_rep  af_to_rep  n_members
vOTU1  SRR02__bin_1  SRR02   yes                SRR02__bin_1    480000     95.50                100.00      100.00     2
vOTU1  SRR01__bin_1  SRR01   no                 SRR02__bin_1    500000     71.00                99.10       98.75      2
vOTU2  SRR03__bin_2  SRR03   yes                SRR03__bin_2    450000     88.00                100.00      100.00     1
vOTU3  SRR02__bin_7  SRR02   yes                SRR02__bin_7    200000     40.00                100.00      100.00     1
```

Three things were checked here: the two near-identical genomes from different samples
merged into vOTU1; the **shorter but more complete** genome was chosen as representative
(95.5% vs 71% complete, even though the other is 20 kb longer), confirming the ranking; and
the 82% ANI genome stayed separate. `SRR02__bin_7` also stayed separate, which is the
greedy centroid behaviour described above: its only match was to `SRR01__bin_1`, which had
already been absorbed as a member and so could not recruit it.

`votu_nearest_ref.py`: run with the mock vOTU table above plus the **real**
`GVDB_stats.tsv` header and one real row, so the actual column names were exercised. It
reported the 96.4% ANI hit with its family (`f_IM_07`) and reference completeness, and
printed `none_above_fastani_threshold` for the vOTU with no reference hit:

```
3 vOTU representatives: 2 with a GVMAGs V2 relative, 1 with none (no close relative sequenced)
```

The shell parts of `09_gvmag_votus.slurm` were checked with `bash -n`, and the `awk` that
merges each sample's summary table into one (keeping a single header and prefixing GVMAG
names with the sample) was tested separately on two mock summary files.

**Caveats.**
- fastANI below ~80% ANI reports nothing, so distant relationships are invisible. For those
  use the marker-gene phylogeny, which is still unwritten.
- An all-vs-all fastANI run is quadratic in the number of genomes. For the tens of GVMAGs
  expected here that is minutes; it would need rethinking for thousands.
- vOTUs are a convention, not a biological truth. A vOTU is not proof of a distinct
  species, and the 95% cut-off will split or merge borderline cases.
- This step says nothing about **endogenous viral elements**: a GVMAG that is really a
  stretch of integrated virus inside the *K. brevis* genome would still cluster into a
  tidy vOTU. That check is still to do (see prompt 19 in `PROMPT_LOG.md`).

----------------------------------------------------------------------------------------

----------------------------------------------------------------------------------------

#### 09f. Running step 09: the production run and its results (2026-11-16)

Written by Claude Code on 2026-11-16 (prompt 38 in `PROMPT_LOG.md`).

**Pre-flight checks.** Three things were confirmed before submitting, because step 09 is
cheap to run but its inputs are easy to get silently wrong:

```bash
# 1. every sample finished step 08, and the smoke-test directory is gone
ls -d ~/scr10/data/kbrevis/gvmags/*/ | wc -l          # 20
ls ~/scr10/data/kbrevis/gvmags/*/.done | wc -l        # 20  (all finished)
ls ~/scr10/data/kbrevis/gvmags/*/gvmags/*.fna | wc -l # 26  (the GVMAGs)
ls -d ~/scr10/data/kbrevis/gvmags/_tailtest           # absent (moved out, see 09a-bis)

# 2. the GVMAGs V2 reference bundle is present
ls ~/data10/databases/gvmagsV2/000_GVMAGS_V2/fna_files | wc -l   # 8508

# 3. both conda environments activate cleanly under `set -euo pipefail`
#    (step 08 job 110103 died in 2 s on openjdk's JAVA_HOME hook; see 08f)
for e in aniclustermap-env checkv-env; do ls ~/.conda/envs/$e/etc/conda/activate.d; done
#   -> neither environment has an activate.d hook, so the bare `conda activate`
#      calls in 09_gvmag_votus.slurm are safe and the script needed no change
```

The third check is the one worth keeping. Step 08 needed `set +u` wrappers around every
`conda activate` because the openjdk package (pulled in for bbmap, which needs Java) ships
an activation hook starting `export JAVA_HOME_CONDA_BACKUP=$JAVA_HOME`, which aborts
instantly when `JAVA_HOME` is unset and `set -u` is on. Step 09's two environments
(`aniclustermap-env`, `checkv-env`) contain no Java and have no activation hooks at all, so
the same fix was not needed here. This was verified by actually running the activations
rather than by reading the package list.

**Choosing where to run it.** `~/hpc/best_cluster.sh -c 16 -m 64 -H 2` ranked gulf `gz`
first (estimated finish 13:51) and kuro `serial` second (14:04). **kuro was chosen anyway.**
The `gz` nodes belong to another group and are preemptible: this account holds only the
`normal` QOS, so a `gz_preempt` job cancels ours outright (`PreemptMode=CANCEL`, no
requeue). Thirteen minutes of estimated finish time is not worth having the all-vs-all
fastANI and the 8,508-genome reference comparison killed partway through (AGENTS.md rule 1).

**Submitted:**

```bash
( module purge; module load slurm/kuro; \
  sbatch -p serial -N 1 -c 16 --mem=64G -t 4:00:00 \
         scripts/09_gvmag_votus.slurm --with-refs )
# Submitted batch job 111457
```

The walltime was cut from the script's own `#SBATCH --time=12:00:00` to 4 h on the command
line (flags override `#SBATCH` lines). A shorter request schedules sooner and is better
behaviour on a shared queue; 4 h was still a large overestimate, as it turned out.

**It finished in 1 minute 7 seconds.**

```
JobID        JobName            State      Elapsed   MaxRSS   ExitCode
111457       09_gvmag_votus     COMPLETED  00:01:07           0:0
111457.batch batch              COMPLETED  00:01:07 15000164K 0:0
```

Peak memory was **14.3 GB** of the 64 GB requested, and the whole job used about 1/200th of
its walltime. `seff` is not installed on kuro (`seff: command not found`), so these figures
come from `sacct -o jobid,jobname,state,elapsed,maxrss,exitcode`. **For any re-run, `-c 8
--mem=24G -t 30:00` is ample.** The 8,508-reference comparison is cheap because fastANI
sketches the reference set once and streams the 19 queries past it.

**Results.** All of these are in `~/scr10/data/kbrevis/votus` (and copied into
`output/votus`, see Step 12):

```
[votus] pooled 26 GVMAGs from 9 samples
[votus] combined summary rows: 26 (pooled GVMAGs: 26)
[votus] fastANI pairs: 97
26 GVMAGs -> 19 vOTUs (6 found in more than one sample/bin, 13 singletons)
            at >=95.0% ANI over >=85.0% of the shorter genome
19 vOTU representatives: 12 with a GVMAGs V2 relative, 7 with none
```

The missing-summary handling predicted in 09a-bis behaved exactly as tested in production:
26 GVMAGs pooled, 26 summary rows found, counts agreed, no warning fired, and the eleven
samples with no `_gvmag_summary.tsv` contributed neither a genome nor a row. One GVMAG
(`SRR28004139_bin_7`, CheckV completeness "Not-determined") has no numeric completeness; it
was ranked last and reported as `NA`, which is the designed behaviour, not a failure.

**The six multi-member vOTUs.** These are the biologically interesting rows, because a vOTU
with more than one member is the same viral "species" recovered independently more than
once. Every one of the six spans **different samples**, not two bins of a single assembly:

| vOTU | members | samples | representative | rep. completeness |
|---|---|---|---|---|
| vOTU02 | 2 | SRR28004137, SRR28004154 | SRR28004154_bin_22 | 100.00 |
| vOTU03 | 2 | SRR28004137, SRR28004154 | SRR28004137_bin_6  | 100.00 |
| vOTU08 | 3 | SRR28004137, SRR28004139, SRR28004154 | SRR28004139_bin_22 | 73.43 |
| vOTU12 | 2 | SRR28004137, SRR28004154 | SRR28004154_bin_8  | 64.79 |
| vOTU17 | 2 | SRR28004137, SRR28004154 | SRR28004154_bin_6  | 46.62 |
| vOTU18 | 2 | SRR28004138, SRR28004139 | SRR28004138_bin_51 | 43.70 |

This reproduces the pattern Booker et al. report, where one vOTU (their vOTU5) was recovered
from two samples and reported once. SRR28004137 and SRR28004154 share four vOTUs, which is
the strongest cross-sample signal in the data.

**Novelty against GVMAGs V2.** `votu_nearest_reference.tsv` gives each representative's
closest relative among the 8,508 reference genomes:

| vOTU | n | completeness | nearest reference | ANI % | AF % | ref. family | ref. ecosystem |
|---|---|---|---|---|---|---|---|
| vOTU01 | 1 | 100.00 | SRX802077-79 | 90.71 | 35.52 | f_Mimiviridae | NA |
| vOTU02 | 2 | 100.00 | bin_3206 | 79.34 | 48.85 | f_Mimiviridae | Freshwater |
| vOTU03 | 2 | 100.00 | IMGM3300017967_BIN561 | 78.55 | 2.53 | f_Mesomimiviridae | Marine |
| vOTU04 | 1 | 100.00 | *none* | NA | NA | NA | NA |
| vOTU05 | 1 | 100.00 | *none* | NA | NA | NA | NA |
| vOTU06 | 1 | 82.18 | IMGM3300040751_BIN751 | 89.22 | 13.77 | f_Mesomimiviridae | Marine |
| vOTU07 | 1 | 77.75 | *none* | NA | NA | NA | NA |
| vOTU08 | 3 | 73.43 | IMGM3300024256_BIN200 | 78.82 | 29.18 | f_Mimiviridae | Marine |
| vOTU09 | 1 | 71.28 | IMGM3300004096_BIN14 | 75.67 | 2.00 | f_Mesomimiviridae | Freshwater |
| vOTU10 | 1 | 70.58 | IMGM3300020408_BIN215 | 94.25 | 45.56 | f_Mesomimiviridae | Marine |
| vOTU11 | 1 | 65.31 | *none* | NA | NA | NA | NA |
| vOTU12 | 2 | 64.79 | IMGM3300037329_BIN195 | 74.44 | 0.92 | f_Allomimiviridae | Non-marine Saline/Alkaline |
| vOTU13 | 1 | 60.77 | IMGM3300023179_BIN729 | 74.60 | 21.59 | f_Mimiviridae | Freshwater |
| vOTU14 | 1 | 58.58 | *none* | NA | NA | NA | NA |
| vOTU15 | 1 | 48.65 | IMGM3300020469_BIN799 | 76.10 | 18.18 | f_AG_04 (Algavirales) | Marine |
| vOTU16 | 1 | 47.83 | IMGM3300046734_BIN24 | 97.95 | 25.60 | f_Mimiviridae | Marine |
| vOTU17 | 2 | 46.62 | *none* | NA | NA | NA | NA |
| vOTU18 | 2 | 43.70 | IMGM3300020187_BIN340 | 77.01 | 5.13 | f_Mimiviridae | Marine |
| vOTU19 | 1 | NA | *none* | NA | NA | NA | NA |

**The headline result: not one of the 19 vOTUs matches a known genome at species level.**
Applying the same MIUViG criteria used to form the vOTUs (>=95% ANI over >=85% of the
shorter genome) to the reference comparison returns nothing:

```bash
awk -F'\t' 'NR>1 && $7!="NA" && $7+0>=95 && $8+0>=85' votu_nearest_reference.tsv
#   -> no rows
```

This needs care in interpretation, because ANI and aligned fraction must be read together.
vOTU16's 97.95% ANI looks like a species-level match until you see it covers only 25.60% of
the genome: the two genomes share a quarter of their sequence almost exactly and differ
elsewhere, which is a relative, not the same virus. vOTU10 (94.25% ANI over 45.56%) is the
closest thing to a known virus in the set. At the other extreme, vOTU12's nearest relative
aligns over 0.92% of the genome, which is barely evidence of anything.

The seven vOTUs with no hit at all are **not** suspect. fastANI reports nothing below about
80% ANI, so "none" means only that no sequenced giant virus is within ~80% ANI. Most
environmental giant viruses have no close sequenced relative, which is exactly why GVMAGs V2
exists. Combined with the step 08 taxonomy (23 Imitervirales, 1 Pimascovirales,
1 Pandoravirales), the picture is a set of novel *Imitervirales*, mostly nearest to marine
Mimiviridae and Mesomimiviridae, which is what a *Karenia brevis* bloom virome should look
like.

**The heatmap.** `output/votus/heatmap/` holds `ANIclustermap.png`, `ANIclustermap.svg`,
the UPGMA dendrogram (`ANIclustermap_dendrogram.nwk`) and the underlying ANI matrix
(`ANIclustermap_matrix.tsv`), the equivalent of the paper's Supplemental Figure 3.

**Caveats.**

- The 19 vOTUs come from 26 GVMAGs recovered from only 9 of 20 samples. This describes the
  giant viruses that **assembled and binned well enough to be recovered**, not the
  diversity that was present. A virus at low abundance leaves no bin.
- Six of 19 representatives have CheckV completeness below 70%, and one is undetermined.
  A fragmentary genome has less sequence available to align, which depresses both ANI and
  aligned fraction against the references, so novelty is partly confounded with
  incompleteness. The highest-completeness vOTUs (01-05, all 100%) are the ones to trust
  most, and three of those five have no reference relative at all.
- Clustering used the MIUViG defaults (`--min-ani 95 --min-af 85`) because Booker et al.
  state only that they grouped GVMAGs "based on their average nucleotide identity" without
  giving a cut-off. A different cut-off would give a different vOTU count; 95/85 is the
  community standard and CheckV `aniclust`'s default.
- Missing completeness changes which genome **represents** a multi-member vOTU (ranking is
  completeness, then length, then name), though not which genomes group together. That was
  demonstrated in 09a-bis and does not affect this run: the one GVMAG with no completeness
  value is a singleton.
- A vOTU is a convention, not a species. These numbers also say nothing about whether any
  GVMAG is actually an endogenous viral element inside the *K. brevis* genome; that check
  is still outstanding.

----------------------------------------------------------------------------------------

# Step 10. Concatenated marker-gene phylogeny of the vOTUs (COMPLETE 2026-11-16)

Written by Claude Code on 2026-11-12 (prompt 26 in `PROMPT_LOG.md`). Files:
`scripts/10_gvmag_phylogeny.slurm`, `scripts/select_refs.py`, `envs/phylo-env.yml`.
**Submitted 2026-11-16 as job 111459** (kuro `serial`, 32 cores, 128 GB, 24 h) on the 19
vOTU representatives from step 09, after fixing an order-parsing bug in `select_refs.py`
that would have silently built the tree with no in-order references at all - see 10i for
the bug, the fix and the run.

#### 10a. What this step is for

Steps 08 and 09 tell you *what* your giant viruses are (GVClass taxonomy) and *which are
the same virus* (vOTUs). Neither tells you how your viruses are **related** to each other
and to known giant viruses. That needs a phylogenetic tree, and for giant viruses the
standard way to build one is a **concatenated marker-gene tree**.

Why concatenated markers rather than whole genomes? Giant virus genomes share almost no
gene content across orders - there is no alignable "whole genome" to compare. But a small
set of genes is present in nearly all of them, mostly involved in replication and
transcription. Aligning those genes separately, then joining the alignments end to end into
one long alignment per genome, gives enough phylogenetic signal to place a new genome in the
*Nucleocytoviricota* tree. This is the approach used by the field's reference phylogenies
and by Booker et al.

The seven default markers (`ncldv_markersearch`'s defaults, which are exactly the seven the
paper lists):

| GVOG | gene | why it is useful |
|---|---|---|
| GVOGm0890 | Poxvirus late transcription factor VLTF3 | near-universal in NCLDV |
| GVOGm0760 | A32 packaging ATPase | packages DNA into the capsid |
| GVOGm0461 | DNA topoisomerase II | DNA replication |
| GVOGm0172 | Transcription initiation factor IIB | transcription |
| GVOGm0054 | DNA polymerase family B | the classic NCLDV phylogenetic marker |
| GVOGm0023 | DNA-directed RNA polymerase alpha subunit | transcription |
| GVOGm0013 | DEAD/SNF2-like helicase | DNA/RNA remodelling |

Two further markers exist in the HMM set but are **not** in the default seven, matching the
paper: GVOGm0003 (major capsid protein) and GVOGm0022 (RNA polymerase beta subunit). Add
them with `-m` if you want, but then the tree is no longer directly comparable to published
ones.

#### 10b. SCRIPT: `scripts/10_gvmag_phylogeny.slurm`

```bash
#!/bin/bash
#SBATCH --job-name=10_gvmag_phylogeny
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/10_gvmag_phylogeny_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/10_gvmag_phylogeny_%j.err

# Build a trimmed, concatenated giant virus marker-gene alignment from the step 09 vOTU
# representatives plus GVMAGs V2 reference genomes, following Booker et al. 2026
# (Supplemental Methods, "GVMAG comparisons"):
#   1. predict proteins for each vOTU representative       (pyrodigal-gv, Prodigal headers)
#   2. choose reference genomes from GVMAGs V2             (select_refs.py)
#   3. find, extract, concatenate and align the 7 default
#      giant virus marker genes                            (ncldv_markersearch.py -c)
#   4. trim the concatenated alignment                      (trimAl 1.5.1 -automated1)
# The output *.trimmed.aln is the input for tree inference (IQ-TREE; see
# README step 10).
#
# Needs: step 09 finished, and 07_gvmag_setup.slurm run with GVMAGS_V2_FULL=1 so the
# 1.18 GB reference bundle is present.
#
# Submit (16+ cores REQUIRED: ncldv_markersearch hardcodes `hmmsearch --cpu 16`):
#   ( module purge; module load slurm/kuro;   sbatch -p serial scripts/10_gvmag_phylogeny.slurm )
#   ( module purge; module load slurm/astral; sbatch -p batch --exclude=as01 scripts/10_gvmag_phylogeny.slurm )
#
# Options (all optional, after the script name):
#   --tier mcl|skani     reference tier: 2,606 (default) or 8,507 representatives
#   --orders LIST        'auto' (default: the orders of our own GVMAGs), 'all', or e.g.
#                        Imitervirales,Algavirales
#   --max-refs N         cap in-order references (useful to keep the alignment small;
#                        0 = no cap, the default)
#   --aligner clustalo|muscle   per-marker aligner (default clustalo, as in the paper)
#   --no-refs            our vOTUs only, no reference genomes

set -euo pipefail

THREADS="${SLURM_CPUS_PER_TASK:?Submit with sbatch}"
REPO="/sciclone/home/iadanforth/K_brevis_viromes"
DATA="$HOME/scr10/data/kbrevis"
VOTU_DIR="$DATA/votus"
GVMAG_ROOT="$DATA/gvmags"
OUT="$DATA/phylogeny"
NMS_DIR="$HOME/data10/programs/ncldv_markersearch"
GVMAGS_V2_DIR="$HOME/data10/databases/gvmagsV2"
PROJECT="kbrevis_gvmags"

TIER=mcl; ORDERS=auto; MAX_REFS=0; ALIGNER=clustalo; USE_REFS=1
while (($#)); do
    case "$1" in
        --tier) TIER="$2"; shift 2 ;;
        --orders) ORDERS="$2"; shift 2 ;;
        --max-refs) MAX_REFS="$2"; shift 2 ;;
        --aligner) ALIGNER="$2"; shift 2 ;;
        --no-refs) USE_REFS=0; shift ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done

if ((THREADS < 16)); then
    echo "WARNING: ncldv_markersearch hardcodes 'hmmsearch --cpu 16'. With only $THREADS" >&2
    echo "cores allocated, hmmsearch will oversubscribe the node. Submit with -c 16 or more." >&2
fi

REPS="$VOTU_DIR/gvmags_representatives.txt"
[[ -s "$REPS" ]] || { echo "No vOTU representatives ($REPS); run step 09 first." >&2; exit 1; }
[[ -d "$NMS_DIR/hmm" ]] || {
    echo "ncldv_markersearch not installed ($NMS_DIR); run 07_gvmag_setup.slurm." >&2; exit 1; }

mkdir -p "$OUT"/{proteins,work}

module load miniforge3
eval "$(conda shell.bash hook)"
conda activate phylo-env

# ---- 1. Proteins for each vOTU representative ---------------------------------
# ncldv_markersearch requires Prodigal-style protein IDs (contigname_1, contigname_2, ...)
# because it joins split marker genes that sit near each other on the same contig; giant
# virus genes are often fragmented across a MAG, so this is essential, not cosmetic.
# It takes the taxon name from the .faa file name, so names must contain no dots.
rm -f "$OUT"/proteins/*.faa
n_q=0
while read -r fna; do
    [[ -s "$fna" ]] || continue
    base=$(basename "$fna" .fna)
    # name the taxon by its vOTU, so tree tips are readable: vOTU003__SRR28004141_bin_2
    votu=$(awk -F'\t' -v g="$base" '$2==g && $4=="yes" {print $1}' \
             "$VOTU_DIR/gvmags_votus.tsv" | head -1)
    name=$(echo "${votu:-noVOTU}__${base}" | tr -c 'A-Za-z0-9_-' '_')
    pyrodigal-gv -i "$fna" -a "$OUT/proteins/${name}.faa" -p meta -j "$THREADS" -o /dev/null
    n_q=$((n_q + 1))
done < "$REPS"
echo "[phylo] predicted proteins for $n_q vOTU representatives"

# ---- 2. Reference genomes from GVMAGs V2 --------------------------------------
if ((USE_REFS)); then
    REF_FAA_DIR=$(find "$GVMAGS_V2_DIR" -maxdepth 3 -type d -name 'faa*' | head -1)
    if [[ -z "$REF_FAA_DIR" ]]; then
        echo "The GVMAGs V2 protein bundle is not downloaded." >&2
        echo "Get it with: ( module purge; module load slurm/astral; sbatch -p batch \\" >&2
        echo "  --exclude=as01 --export=ALL,GVMAGS_V2_FULL=1 scripts/07_gvmag_setup.slurm )" >&2
        exit 1
    fi
    python -I "$REPO/scripts/select_refs.py" \
        --metadata "$GVMAGS_V2_DIR/GVDB_stats.tsv" \
        --faa-dir "$REF_FAA_DIR" \
        --out-dir "$OUT/proteins" \
        --tier "$TIER" --orders "$ORDERS" --max-refs "$MAX_REFS" \
        --gvclass-summaries "$GVMAG_ROOT"/*/*_gvmag_summary.tsv
else
    echo "[phylo] --no-refs: building the alignment from our vOTUs alone"
fi

N_TAXA=$(ls "$OUT"/proteins/*.faa | wc -l)
echo "[phylo] $N_TAXA taxa in the alignment ($n_q ours + $((N_TAXA - n_q)) references)"

# Reference .faa files come from the GVMAGs V2 bundle and were not necessarily produced by
# Prodigal. If their headers are not contigname_N, markersearch cannot join split genes for
# those taxa (it still finds single-copy hits, so the run succeeds but fragmented reference
# genes may be missed). Report rather than fail.
BAD=$(for f in "$OUT"/proteins/*.faa; do
          head -1 "$f" | grep -qE '^>[^ ]+_[0-9]+' || basename "$f"; done | wc -l)
((BAD == 0)) && echo "[phylo] all taxa have Prodigal-style protein IDs" \
             || echo "[phylo] NOTE: $BAD of $N_TAXA .faa files lack Prodigal-style IDs; split marker genes in those genomes may not be joined" >&2

# ---- 3. ncldv_markersearch: find, extract, concatenate, align -----------------
# Two quirks drive the directory juggling below:
#  * the script hardcodes a RELATIVE HMM path ("hmm/gvogs.common.hmm"), so it must run from
#    a directory containing `hmm/` - we symlink the installed copy;
#  * it writes all outputs into the current directory, so that same directory is the output.
# Default marker set is exactly the paper's seven: VLTF3 (GVOGm0890), A32 packaging ATPase
# (0760), DNA topoisomerase II (0461), TFIIB (0172), DNA polymerase B (0054), RNA polymerase
# alpha (0023) and SNF2-like helicase (0013). The major capsid protein (0003) and RNAP beta
# (0022) are available but NOT in the default set, as in the paper.
cd "$OUT/work"
ln -sfn "$NMS_DIR/hmm" hmm
echo "[phylo] running ncldv_markersearch (aligner: $ALIGNER); this is the slow step"
python -I "$NMS_DIR/ncldv_markersearch.py" \
    -i "$OUT/proteins" \
    -n "$PROJECT" \
    -t "$THREADS" \
    -g "$ALIGNER" \
    -c

CONCAT="$OUT/work/${PROJECT}.concat.aln"
[[ -s "$CONCAT" ]] || { echo "No concatenated alignment was produced ($CONCAT)." >&2; exit 1; }
echo "[phylo] concatenated alignment: $(grep -c '^>' "$CONCAT") taxa"

# ---- 4. trimAl -automated1 ----------------------------------------------------
# -automated1 chooses between gap-based and similarity-based trimming using a heuristic
# tuned for maximum-likelihood phylogenetics; it is what the paper used. The HTML report
# shows exactly which columns were dropped, which is worth checking before tree building.
TRIMMED="$OUT/${PROJECT}.trimmed.aln"
trimal -in "$CONCAT" -out "$TRIMMED" -automated1 -fasta \
       -htmlout "$OUT/${PROJECT}.trimal_report.html"
trimal -in "$TRIMMED" -sident > "$OUT/${PROJECT}.trimmed.identity.txt" 2>&1 || true

# ---- summary ------------------------------------------------------------------
cp "$OUT/work/${PROJECT}.table.tsv" "$OUT/" 2>/dev/null || true
cp "$OUT/work/${PROJECT}.full_output.txt" "$OUT/" 2>/dev/null || true
aln_len () {   # length of the first sequence in a FASTA alignment
    python -I -c "
import sys
seq = []
for line in open(sys.argv[1]):
    if line.startswith('>'):
        if seq: break
    else:
        seq.append(line.strip())
print(len(''.join(seq)))" "$1"
}
{
  echo "project:        $PROJECT"
  echo "date:           $(date)"
  echo "our vOTUs:      $n_q"
  echo "total taxa:     $(grep -c '^>' "$TRIMMED")"
  echo "reference tier: $([[ $USE_REFS == 1 ]] && echo "$TIER, orders=$ORDERS, max=$MAX_REFS" || echo none)"
  echo "aligner:        $ALIGNER"
  echo "columns before trimming: $(aln_len "$CONCAT")"
  echo "columns after trimming:  $(aln_len "$TRIMMED")"
} | tee "$OUT/${PROJECT}.phylogeny_summary.txt"

echo
echo "[phylo] done. Trimmed alignment for tree building:"
echo "  $TRIMMED"
echo "Per-marker occupancy per taxon: $OUT/${PROJECT}.table.tsv"
echo "trimAl report (which columns were dropped): $OUT/${PROJECT}.trimal_report.html"
echo "Next: infer a tree with IQ-TREE (see README step 10)."
```

#### 10c. What each tool does, and the parameters chosen

**`pyrodigal-gv` (step 1 of the script).** Predicts proteins from each vOTU
representative's nucleotide FASTA. This is needed because `ncldv_markersearch` works on
proteins, not DNA. Critically, it must produce **Prodigal-style protein IDs**
(`contigname_1`, `contigname_2`, ...): `ncldv_markersearch` parses those IDs to work out
which proteins are neighbours on a contig, so that a marker gene **split across two
predicted ORFs** can be joined back together. Giant virus genes in MAGs are frequently
split that way, so this is essential rather than cosmetic. I verified pyrodigal-gv's output
format on a test sequence: `>my_contig_name_1 # 1 # 882 # 1 # ID=1_1;...`, which is correct.
`-p meta` is used because these are metagenomic sequences rather than a single organism.

**`select_refs.py` (step 2).** Chooses which reference genomes to include. See 10d.

**`ncldv_markersearch.py -c` (step 3).** In one pass it: runs `hmmsearch` with the nine
giant virus marker HMMs against every taxon's proteins; keeps hits above its built-in
bit-score thresholds; joins hits to the same marker that are non-overlapping and within
`-p` genes of each other on one contig (default 3); aligns each marker separately
(Clustal Omega by default, `-g muscle` for MUSCLE5); and concatenates the per-marker
alignments into `<project>.concat.aln`, one entry per taxon, padding missing markers with
runs of `X`. The `-c` flag is what produces the concatenated alignment, and the paper used
exactly this option.

Three quirks of this script that the job script works around, all found by reading its
source:
- It hardcodes a **relative** HMM path, `speci_db = "hmm/gvogs.common.hmm"`, so it only
  works when the working directory contains `hmm/`. The job script symlinks the installed
  `hmm/` into its own work directory.
- It writes **all outputs into the current working directory**, prefixed with `-n`, so that
  same work directory is also the output directory.
- It **hardcodes `hmmsearch --cpu 16`**, ignoring `-t` for that step. So the job must
  request at least 16 cores or it will oversubscribe the node, which is exactly the kind of
  thing the front-end watchdog and your node-sharing neighbours notice. The script warns if
  `--cpus-per-task` is below 16; it requests 32.

It also writes `.domout` files into the **input** directory, which is why the script
rebuilds `proteins/` from scratch on each run (`rm -f`): stale `.domout` files from a
previous run with different taxa would otherwise be silently reused.

**`trimal -automated1` (step 4).** Alignment columns that are mostly gaps, or too variable
to be homologous, add noise rather than signal and can actively mislead tree inference.
trimAl removes them. `-automated1` is a heuristic that inspects the alignment and chooses
between gap-based and similarity-based trimming automatically, optimised for maximum
likelihood phylogenetics - it is what the paper used. Version **1.5.1** here, the newest
release; the paper used 1.4.22. The script also writes `-htmlout` report showing exactly
which columns were dropped, worth looking at before you trust the tree, plus a `-sident`
pairwise identity summary.

#### 10d. Choosing reference genomes: the most extensive database available

The request was to use the most extensive giant virus database including Vasquez et al.
2025, which is **GVMAGs V2**: 18,725 genomes reduced to 8,508 ANI-cluster representatives,
whose protein FASTAs ship in `GVMAGS_V2_database.gz` (1.18 GB). For comparison, the
Booker et al. paper used the Aylward GVDB from September 2024, roughly 1,382
representatives - so this is about a six-fold larger reference set.

`select_refs.py` reads the release's `GVDB_stats.tsv` and offers two nested tiers, whose
sizes I measured directly from that file:

| tier | column | representatives | Imitervirales among them |
|---|---|---|---|
| `mcl` (default) | `mcl_rep` | 2,606 | 1,049 |
| `skani` | `skani_rep` | 8,507 | 5,461 |

`mcl` is the default because 2,606 is closest to the scale of the published tree (2,037
genomes after filtering) and because 8,507 taxa make the Clustal Omega step punishingly
slow.

On top of the tier, references are filtered **by viral order**, which matters more than the
tier does. `--orders auto` (the default) reads your own GVMAGs' orders from the step 08
GVClass summaries and keeps references from those orders. All 11 GVMAGs in the paper were
*Imitervirales*, so if yours are too, that is 1,049 references at the `mcl` tier. A tree of
Imitervirales is far more informative about Imitervirales vOTUs than one diluted with
distant orders.

Because an unrooted tree has no direction, the script also keeps a few genomes from **every
other order** (`--outgroup-per-order`, default 3) so the tree can be rooted - the paper
midpoint-rooted theirs. And `--max-refs N` caps the in-order references, sampling **evenly
across families** rather than at random so the subset stays taxonomically spread rather than
dominated by whichever family happens to be largest.

#### 10e. SCRIPT: `scripts/select_refs.py`

```python
#!/usr/bin/env python3
"""Choose giant virus reference genomes from GVMAGs V2 for the phylogenomic tree.

GVMAGs V2 (Vasquez et al. 2025, doi:10.1101/2025.09.26.678796) is the largest curated
giant virus genome collection: 18,725 genomes reduced to 8,508 ANI-cluster
representatives, whose protein FASTAs ship in GVMAGS_V2_database.gz. Putting all 8,508
into one concatenated alignment is possible but slow, and far beyond what a Bayesian
phylogenetic program can comfortably handle, so this script picks a defensible subset.

Two nested tiers are available, both read from the release's GVDB_stats.tsv:

  --tier mcl    the coarser MCL clustering (column `mcl_rep`): 2,606 representatives.
                Closest in scale to the published tree, which had 2,037 genomes after
                filtering. This is the default.
  --tier skani  every ANI-cluster representative (column `skani_rep`): 8,507. Finest
                resolution, slowest.

On top of the tier you can restrict by viral order, which is usually what you want: a
tree of Imitervirales genomes is far more informative about Imitervirales vOTUs than one
diluted with distant orders. `--orders auto` reads the orders of your own GVMAGs from the
step 08 GVClass summaries and keeps references from those orders, plus (unless
--no-outgroup) a small, even sample of other orders to root the tree.

--max-refs caps the number of in-order references, sampled evenly across families so the
subset stays diverse rather than dominated by one clade.

Output: the chosen references' .faa files copied into --out-dir with sanitised names
(no dots, because ncldv_markersearch derives the taxon name from the file name), and a
manifest TSV recording what was chosen and why.
"""
import argparse
import csv
import random
import re
import shutil
import sys
from collections import defaultdict
from pathlib import Path


def sanitise(name):
    """ncldv_markersearch uses the .faa basename as the taxon name and does regex
    substitutions on '.faa', so dots and whitespace in names must go."""
    return re.sub(r"[^A-Za-z0-9_-]", "_", name)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--metadata", required=True, help="GVDB_stats.tsv")
    p.add_argument("--faa-dir", required=True, help="directory of reference .faa files")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--tier", choices=["mcl", "skani"], default="mcl")
    p.add_argument("--orders", default="auto",
                   help="'auto' (from --gvclass-summaries), 'all', or a comma-separated list")
    p.add_argument("--gvclass-summaries", nargs="*", default=[],
                   help="step 08 *_gvmag_summary.tsv files, for --orders auto")
    p.add_argument("--max-refs", type=int, default=0, help="0 = no cap")
    p.add_argument("--outgroup-per-order", type=int, default=3,
                   help="references kept from each OTHER order, to root the tree (0 = none)")
    p.add_argument("--seed", type=int, default=1)
    a = p.parse_args()

    rows = list(csv.DictReader(open(a.metadata), delimiter="\t"))
    by_genome = {r["genome"]: r for r in rows}
    col = "mcl_rep" if a.tier == "mcl" else "skani_rep"
    reps = sorted({r[col] for r in rows if r[col] not in ("NA", "")})
    print(f"[refs] tier '{a.tier}' ({col}): {len(reps)} representatives in the metadata")

    # ---- which orders do we care about? ----
    if a.orders == "all":
        wanted = None
    elif a.orders == "auto":
        wanted = set()
        for f in a.gvclass_summaries:
            if not Path(f).exists():
                continue
            with open(f) as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    o = (row.get("gvclass_order") or "").strip()
                    if o and o not in ("o_", "NA"):
                        wanted.add(o.removeprefix("o_"))
        if not wanted:
            sys.exit("--orders auto found no orders in the GVClass summaries; pass --orders "
                     "explicitly (e.g. --orders Imitervirales) or 'all'")
        print(f"[refs] orders taken from our own GVMAGs: {', '.join(sorted(wanted))}")
    else:
        wanted = {o.strip() for o in a.orders.split(",") if o.strip()}
        print(f"[refs] orders requested: {', '.join(sorted(wanted))}")

    # ---- split representatives into in-order and other-order ----
    in_order, others_by_order = [], defaultdict(list)
    for g in reps:
        r = by_genome.get(g)
        order = (r or {}).get("Order", "NA")
        if wanted is None or order in wanted:
            in_order.append(g)
        else:
            others_by_order[order].append(g)
    print(f"[refs] {len(in_order)} representatives in those orders, "
          f"{sum(len(v) for v in others_by_order.values())} in others")

    # ---- cap, sampling evenly across families to keep the subset diverse ----
    rng = random.Random(a.seed)
    if a.max_refs and len(in_order) > a.max_refs:
        by_family = defaultdict(list)
        for g in in_order:
            by_family[(by_genome.get(g) or {}).get("gvclass_family", "NA")].append(g)
        for v in by_family.values():
            rng.shuffle(v)
        picked, families = [], sorted(by_family)
        while len(picked) < a.max_refs and any(by_family[f] for f in families):
            for f in families:
                if by_family[f] and len(picked) < a.max_refs:
                    picked.append(by_family[f].pop())
        print(f"[refs] capped to {len(picked)} references across {len(families)} families "
              f"(--max-refs {a.max_refs})")
        in_order = sorted(picked)

    # ---- outgroups: a few genomes from each other order, to root the tree ----
    outgroups = []
    if wanted is not None and a.outgroup_per_order > 0:
        for order, gs in sorted(others_by_order.items()):
            if order in ("NA", ""):
                continue
            rng.shuffle(gs)
            outgroups.extend(gs[:a.outgroup_per_order])
        print(f"[refs] plus {len(outgroups)} outgroup references from "
              f"{len(others_by_order)} other orders")

    # ---- copy the .faa files ----
    out = Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)
    faa_dir = Path(a.faa_dir)
    index = {f.stem: f for f in faa_dir.rglob("*.faa")}
    print(f"[refs] {len(index)} .faa files available in {faa_dir}")

    copied, missing = 0, []
    with open(out.parent / "reference_manifest.tsv", "w") as man:
        man.write("taxon_file\tgenome\trole\torder\tfamily\tgenus\tlength_bp\tecosystem\n")
        for g, role in [(x, "in_order") for x in in_order] + \
                       [(x, "outgroup") for x in outgroups]:
            src = index.get(g)
            if src is None:
                missing.append(g); continue
            name = sanitise(g)
            shutil.copy(src, out / f"{name}.faa")
            r = by_genome.get(g, {})
            man.write("\t".join([name, g, role, r.get("Order", "NA"),
                                 r.get("gvclass_family", "NA"), r.get("new_genus", "NA"),
                                 r.get("LENbp", "NA"), r.get("Ecosystem", "NA")]) + "\n")
            copied += 1

    print(f"[refs] copied {copied} reference .faa files into {out}")
    if missing:
        print(f"[refs] WARNING: {len(missing)} selected genomes had no .faa in the bundle "
              f"(first few: {', '.join(missing[:5])})", file=sys.stderr)
    if copied == 0:
        sys.exit("No reference proteins were copied; check --faa-dir points at the "
                 "unpacked GVMAGs V2 bundle.")


if __name__ == "__main__":
    main()
```

#### 10f. How to run it

First make sure the one-time setup has installed `phylo-env` and `ncldv_markersearch`, and
downloaded the 1.18 GB reference bundle (the metadata table alone is not enough for this
step):

```
# phylo-env + ncldv_markersearch (also idempotent; added in prompt 26)
( module purge; module load slurm/gulf; sbatch -p batch -N 1 -c 4 -C highmem --mem=16G \
    -t 2:00:00 scripts/07_gvmag_setup.slurm )                      # job 92396

# the 1.18 GB GVMAGs V2 genome+protein bundle, needed for reference genomes
( module purge; module load slurm/gulf; sbatch -p batch -N 1 -c 4 -C highmem --mem=16G \
    -t 2:00:00 --export=ALL,GVMAGS_V2_FULL=1 scripts/07_gvmag_setup.slurm )
```

Then, after step 09:

```
# default: mcl tier, orders auto, all in-order references, Clustal Omega
( module purge; module load slurm/kuro; sbatch -p serial scripts/10_gvmag_phylogeny.slurm )

# a smaller alignment: ~120 in-order references plus outgroups
( module purge; module load slurm/kuro; sbatch -p serial scripts/10_gvmag_phylogeny.slurm \
    --max-refs 120 )

# our vOTUs only, no references (quick sanity check of the marker search)
( module purge; module load slurm/kuro; sbatch -p serial scripts/10_gvmag_phylogeny.slurm \
    --no-refs )
```

Outputs in `~/scr10/data/kbrevis/phylogeny/`:

| file | contents |
|---|---|
| `kbrevis_gvmags.trimmed.aln` | **the trimmed concatenated alignment - the input for tree building** |
| `kbrevis_gvmags.trimal_report.html` | which alignment columns trimAl dropped |
| `kbrevis_gvmags.table.tsv` | how many of the 7 markers were found per taxon |
| `kbrevis_gvmags.full_output.txt` | every marker hit, with bit scores and whether genes were joined |
| `kbrevis_gvmags.phylogeny_summary.txt` | taxa counts, alignment length before and after trimming |
| `reference_manifest.tsv` | which references were chosen, their taxonomy and ecosystem |
| `work/` | `ncldv_markersearch` intermediates, including the untrimmed `.concat.aln` |

**Check before building a tree:** open `kbrevis_gvmags.table.tsv` and look at marker
occupancy. A taxon with only one or two of the seven markers contributes mostly `X` padding
and can distort the tree; such taxa are usually better dropped. The paper's GVMAGs had 7-9
markers each.

#### 10g. How this was tested

`select_refs.py` was tested on dev against the **real** `GVDB_stats.tsv`, with a mock
protein bundle built from every `mcl_rep` genome and a mock step 08 GVClass summary
declaring `o_Imitervirales`:

```
[refs] tier 'mcl' (mcl_rep): 2606 representatives in the metadata
[refs] orders taken from our own GVMAGs: Imitervirales
[refs] 1049 representatives in those orders, 1557 in others
[refs] capped to 60 references across 12 families (--max-refs 60)
[refs] plus 30 outgroup references from 12 other orders
[refs] copied 90 reference .faa files
```

Checked: order auto-detection worked, the `--max-refs` cap sampled evenly across all 12
families (6 per family rather than 60 from one), outgroups were drawn from every other
order, and no taxon name contained a dot.

The `ncldv_markersearch` + trimAl chain was smoke-tested as a **Slurm job on astral**
(job 73841, 16 cores, 14 s) using the three giant virus genomes bundled with the repo
(medusavirus, megavirus, tupanvirus), because `hmmsearch --cpu 16` must not run on dev:

```
=== concat alignment:
>tupanvirus concatenated alignment of GVOGm0890,GVOGm0760,GVOGm0461,GVOGm0172,GVOGm0054,GVOGm0023,GVOGm0013
>megavirus  concatenated alignment of ...
>medusavirus concatenated alignment of ...
columns before trimming: 6671
columns after trimming : 4106
=== marker occupancy table:
genome       GVOGm0013 GVOGm0023 GVOGm0054 GVOGm0172 GVOGm0461 GVOGm0760 GVOGm0890
medusavirus        0.0       0.0       1.0       0.0       0.0       1.0       1.0
megavirus          1.0       1.0       1.0       1.0       1.0       1.0       1.0
tupanvirus         1.0       1.0       1.0       1.0       1.0       1.0       1.0
```

Four things this confirmed: the default marker set is **exactly the paper's seven** (the
GVOG IDs are printed in the FASTA description); the relative-`hmm/` symlink and
write-to-CWD handling work; trimAl `-automated1` removed 38% of columns (6,671 to 4,106);
and the occupancy table does flag poorly covered taxa (medusavirus has only 3 of 7 markers,
so in a real run it would be a candidate for removal).

The test also revealed that the bundled test proteins use headers like
`medusavirus.fna|AP018495.1_1` - i.e. `<genome>|<contig>_<n>`, the format the repo's own
`prodigal_launcher.py` produces. Both that and our plain pyrodigal-gv `contig_1` format end
in `_<number>`, which is what the script's contig-name regex (`re.sub("_\d*$", "", ...)`)
needs, so both work. The taxon name in the alignment comes from the **.faa file name**, not
the headers, which is why the job script names the files after the vOTU.

`bash -n` passes on `10_gvmag_phylogeny.slurm`; the Python helpers parse. The full script
has **not** been run, because there are no vOTU representatives yet.

#### 10h. Caveats

- `ncldv_markersearch` is an **unversioned GitHub script** with no releases or tags, so the
  setup script records the installed commit in
  `~/data10/programs/ncldv_markersearch/INSTALLED_COMMIT` for reproducibility.
- Reference `.faa` files come from the GVMAGs V2 bundle and were not necessarily produced by
  Prodigal. If their headers are not `contigname_N`, split marker genes in those genomes
  cannot be joined - single-copy hits still work, so the run succeeds, but some reference
  genes may be missed. The script counts and reports how many taxa look non-Prodigal.
- The concatenated alignment pads missing markers with `X`. That is standard practice, but
  taxa with few markers are poorly placed.
- Clustal Omega on thousands of sequences is the slow step: with the full `mcl` tier expect
  many hours, which is why the job asks for 2 days.
- This step builds an **alignment**, not a tree. Tree inference is deliberately separate, so
  you can try several methods on the same alignment. The paper used IQ-TREE v2.2.2.6 with
  `-m MFP -bb 1000 -nt AUTO -wbt`, then RAxML v8.2.12 (`-f i -m GTRCAT`) for tree certainty
  scores; neither is installed yet.

----------------------------------------------------------------------------------------

#### 10i. A silent order-parsing bug in `select_refs.py`, found before the first run (2026-11-16)

Written by Claude Code on 2026-11-16 (prompt 39 in `PROMPT_LOG.md`).

Step 10 had never been run. Before submitting it, the reference-selection step was checked
against the real data, and it was broken in a way that would **not** have produced an error.

**What was wrong.** `select_refs.py --orders auto` is supposed to read the viral orders of
our own GVMAGs out of the step 08 GVClass summaries and keep reference genomes from those
same orders. It did this:

```python
o = (row.get("gvclass_order") or "").strip()
if o and o not in ("o_", "NA"):
    wanted.add(o.removeprefix("o_"))
```

But `gvclass_order` in a GVClass summary is **not an order name - it is a vote string.**
The real values in our summaries look like this:

```
NCLDV__Imitervirales:10(100.00%)
NCLDV__Imitervirales:20(55.56%),BAC__UM-FILTER-46-8:8(22.22%),BAC__JAJFKO01:3(8.33%),...
NCLDV__Pandoravirales:1(100.00%)
```

GVClass reports how many marker genes voted for each order and with what support, rather
than a single label. So `wanted` was filled with 24 distinct strings of that form. The next
stage compares them against GVDB_stats.tsv's `Order` column, whose values are bare names
(`Imitervirales`, `Algavirales`, `Pimascovirales`, ...). **Nothing ever matched.**

**Why that was dangerous rather than merely broken.** With `in_order` empty, the script did
not fail. It fell through to the outgroup logic, which keeps 3 genomes from each *other*
order, and happily copied ~24 outgroup proteomes. The job would then have run to completion
and produced a perfectly well-formed 43-taxon alignment and tree from our 19 vOTUs plus two
dozen distant outgroups - with **zero** *Imitervirales* references, the very genomes the
tree exists to place our viruses among. Nothing in the logs would have looked wrong. This
is the same class of error as the `domain` column in step 08 (section 08c): a GVClass
column that reads like a label but is a vote tally.

**The fix.** The resolved lineage lives in `gvclass_taxonomy_majority`, which is clean:

```
d_NCLDV;p_Nucleocytoviricota;c_Megaviricetes;o_Imitervirales;f_IM_01;g_;s_singleton
```

A new `order_of()` helper reads the `o_` field from that column, falls back to parsing the
vote string only if the lineage column is absent, and returns `None` for an unresolved
order (an empty `o_`), which correctly contributes nothing:

```python
def order_of(row):
    tax = (row.get("gvclass_taxonomy_majority")
           or row.get("taxonomy_majority") or "").strip()
    for field in tax.split(";"):
        field = field.strip()
        if field.startswith("o_"):
            name = field[2:].strip()
            return name or None
    vote = (row.get("gvclass_order") or "").strip()
    m = re.search(r"NCLDV__([A-Za-z_]+?):", vote)
    return m.group(1) if m else None
```

A second change makes this class of mistake impossible to repeat quietly. The script now
**refuses to continue** if an order filter matches no reference at all:

```python
if wanted is not None and not in_order:
    sys.exit("No reference genome matched the requested orders "
             f"({', '.join(sorted(wanted))}). Those names must match GVDB_stats.tsv's "
             "`Order` column exactly (e.g. 'Imitervirales', no 'o_' prefix). Refusing "
             "to continue: the alignment would contain only outgroups, which looks "
             "like a valid tree but is not the one you asked for.")
```

That guard matters more than the fix itself: the original bug was invisible precisely
because an outgroup-only alignment is structurally valid. The previous version of the file
was kept as `scripts/select_refs.py.bak`.

**How it was tested.** First a dry run that exercised the real selection logic against the
real metadata and the real summaries, but counted instead of copying 1,547 protein files
(copying is job work, not front-end work):

```bash
# imports scripts/select_refs.py and calls order_of() on every step 08 summary row,
# then replays the tier/order/outgroup split
nice -n 19 python3 -I <dry-run script>
```

```
orders detected from our GVMAGs: ['Imitervirales', 'Pandoravirales', 'Pimascovirales']
tier mcl  :  2606 reps ->  1523 in-order + 24 outgroups = 1566 taxa with our 19 vOTUs
   in-order breakdown: {'Imitervirales': 1049, 'Pimascovirales': 361, 'Pandoravirales': 113}
tier skani:  8507 reps ->  6223 in-order + 24 outgroups = 6266 taxa with our 19 vOTUs
```

Then the production job confirmed the same numbers, which is the real test:

```
[refs] tier 'mcl' (mcl_rep): 2606 representatives in the metadata
[refs] orders taken from our own GVMAGs: Imitervirales, Pandoravirales, Pimascovirales
[refs] 1523 representatives in those orders, 1083 in others
[refs] plus 24 outgroup references from 10 other orders
[refs] copied 1547 reference .faa files
[phylo] 1566 taxa in the alignment (19 ours + 1547 references)
[phylo] all taxa have Prodigal-style protein IDs
```

Compare the pre-fix behaviour: 0 in-order, ~24 outgroups, 43 taxa.

**Other pre-flight checks.** `phylo-env` has no `activate.d` hook, so the bare
`conda activate phylo-env` under `set -euo pipefail` is safe (the step 08 `JAVA_HOME`
problem does not apply). The `faa_files` half of the GVMAGs V2 bundle is present
(8,508 files), which `--tier mcl` needs.

**Submitted:**

```bash
( module purge; module load slurm/kuro; \
  sbatch -p serial -N 1 -c 32 --mem=128G -t 24:00:00 \
         scripts/10_gvmag_phylogeny.slurm )
# Submitted batch job 111459
```

Defaults were used (`--tier mcl --orders auto`), giving 1,566 taxa. That is deliberately
close to the published tree's scale (2,037 genomes after filtering), and section 10d
explains why the `mcl` tier rather than `skani` (6,266 taxa) is the right default. kuro
`serial` was chosen over gulf `gz` despite `gz` ranking first on estimated finish: a job of
this length must not sit on a preemptible partition where a `gz_preempt` job would cancel
it outright with no requeue. Walltime was cut from the script's `2-00:00:00` to 24 h, which
is still ample - `ncldv_markersearch` is the slow step, dominated by Clustal Omega on seven
marker alignments of up to ~1,500 sequences each.

Note that `ncldv_markersearch.py` emits several `SyntaxWarning: invalid escape sequence`
messages on Python 3.12 (`'\s+'`, `"_\d*$"`, `"\*"` written without raw-string prefixes).
These are harmless warnings from the third-party script, not errors, and appear in the
`.err` log of every step 10 run.

----------------------------------------------------------------------------------------

#### 10j. Running step 10: the production run and its results (2026-11-16)

Written by Claude Code on 2026-11-16 (prompt 42 in `PROMPT_LOG.md`).

Job **111459** (kuro `serial`, 32 cores, 128 GB, 24 h) **completed in 1 h 09 m 44 s**, exit 0.

```
JobID        State      Elapsed   MaxRSS
111459       COMPLETED  01:09:44
111459.batch COMPLETED  01:09:44  3634196K
```

Peak memory was **3.6 GB of the 128 GB requested** - a 35x over-request. `ncldv_markersearch`
runs hmmsearch one proteome at a time and Clustal Omega one marker at a time, so nothing ever
holds the whole dataset in memory. **A re-run needs `-c 32 --mem=16G -t 4:00:00`**, and the
`#SBATCH --mem=128G` in the script is far too generous; it was left as written because the
job had already been submitted, but it should be reduced.

**Results:**

```
[phylo] predicted proteins for 19 vOTU representatives
[refs] 1523 representatives in those orders, 1083 in others
[refs] plus 24 outgroup references from 10 other orders
[phylo] 1566 taxa in the alignment (19 ours + 1547 references)
[phylo] all taxa have Prodigal-style protein IDs
[phylo] concatenated alignment: 1564 taxa
columns before trimming: 46737
columns after trimming:  37739
```

Note **1566 taxa went in and 1564 came out.** `ncldv_markersearch` silently drops a taxon in
which it found none of the seven markers, and exactly two were dropped: **vOTU05**
(`SRR28004151_bin_29`) and **vOTU19** (`SRR28004139_bin_7`). They were the only 2 of all
1,566 taxa with zero markers. So 17 of our 19 vOTUs are in the tree, and no decision about
removing them was needed - the tool made it.

**vOTU05 deserves a second look on its own merits.** CheckV called it 100% complete at
203,021 bp (section 09f), it has no relative in GVMAGs V2, and it carries **none** of the
seven core *Nucleocytoviricota* markers - no VLTF3, no A32 packaging ATPase, no DNA
topoisomerase II, no TFIIB, no DNA polymerase B, no RNA polymerase alpha, no SNF2 helicase.

> **Caveat on those gene names, added 2026-11-24 (prompt 69).** The seven markers are
> identified in this run only by GVOG accession - `GVOGm0013`, `GVOGm0023`, `GVOGm0054`,
> `GVOGm0172`, `GVOGm0461`, `GVOGm0760`, `GVOGm0890`. The familiar names used in the
> sentence above are **not derived from anything in this project** and have never been
> checked against it. Verified on 2026-11-24 that they cannot be: the HMM models in
> `work/hmm/gvogs.common.hmm` carry only `NAME GVOGm####.trim`, with no `DESC`, `ACC` or
> other free-text field, and all 8,508 GVMAGs V2 proteomes have bare headers
> (`>Fadolivirus|NODE_1_length_...`) with no functional annotation. Resolving an accession
> to a gene name needs the Giant Virus Database's own GVOG annotation table, which this
> project does not hold. **Treat the names above as unconfirmed.** The count - zero of
> seven - is solid and is what the argument rests on. Note also that the same directory
> holds a *different*, properly named HMM set (`NCLDV.hmm`: A32, D5, SFII, mcp, mRNAc,
> PolB, RNAPL, RNAPS, RNR, VLTF3) which the phylogeny did **not** use, so those names do
> not map onto these seven columns either.
A genuinely complete giant virus genome should carry most of these. The likelier
explanations are that it is a mis-binned contig set, or something real but so divergent that
the GVOG HMMs do not detect it.

> **Correction, 2026-11-24 (prompt 67).** This paragraph originally continued: *"Its
> GVClass lineage stopped at `d_NCLDV;p_Nucleocytoviricota;c_Megaviricetes;o_;f_;g_;s_` -
> no order - which is consistent with very few marker hits."* **That is wrong, and the
> lineage quoted is vOTU15's, not vOTU05's.** Checked against
> `output/votus/all_gvmag_summaries.tsv`: vOTU05 is
> `d_NCLDV;p_Nucleocytoviricota;c_Megaviricetes;o_Imitervirales;f_;g_;s_` - GVClass **does**
> assign it an order. What is true is that the assignment rests on very little: 2 marker
> votes, 1 unique GVOG8, `order_completeness` 11.4%, and `avgdist` 3.15, the largest of any
> of our 26 GVMAGs. So the caution below still stands, but on weaker grounds than first
> written: the evidence is thin, not absent.

**It should not be reported as a giant virus genome without
further evidence.** vOTU19 is the CheckV "Not-determined" genome, so its weakness is less
surprising.

**Marker occupancy** (`kbrevis_gvmags.table.tsv`), across all 1,566 taxa before the drop:

| markers found | taxa |
|---|---|
| 7 | 583 |
| 6 | 326 |
| 5 | 258 |
| 4 | 204 |
| 3 | 137 |
| 2 | 56 |
| 0 | 2 |

Our own vOTUs do well: **11 have all 7 markers**, 5 have 6, and only vOTU13 is poor with 2.
That is a better occupancy profile than the reference set as a whole, which is reassuring -
our GVMAGs are not systematically more fragmentary than the published genomes they are being
placed among.

**Alignment quality, and why it needs stating plainly.** `trimal -sident` reports:

```
## MaxIdentity      0.9306
## AverageIdentity  0.1450
```

**Average pairwise identity is 14.5%.** That is expected when concatenating markers across
three viral orders spanning hundreds of millions of years of divergence, and it is the
regime these marker genes were chosen for - but it is near the floor of what maximum
likelihood can resolve, and it predicts a tree with poorly supported deep branches. Treat
deep relationships in the step 11 tree with caution and lean on the UFBoot and IC/TC scores
rather than the topology alone.

trimAl `-automated1` removed only 19% of columns (46,737 -> 37,739), which is a lenient
trim. Section 11h measures what stricter trimming would retain and explains why the
alignment was nonetheless used as the paper specifies.

**Outputs** in `~/scr10/data/kbrevis/phylogeny`:

| file | what it is |
|---|---|
| `kbrevis_gvmags.trimmed.aln` | 58 MB, the alignment step 11 consumes |
| `kbrevis_gvmags.table.tsv` | per-marker occupancy for every taxon |
| `kbrevis_gvmags.trimal_report.html` | 164 MB; which columns trimAl dropped. Large enough that it should be opened with care, not in a light editor |
| `kbrevis_gvmags.trimmed.identity.txt` | 27 MB pairwise identity matrix; one very long line per taxon, so do not `head` it casually |
| `reference_manifest.tsv` | which references were chosen, their role (in_order/outgroup), taxonomy and ecosystem |
| `kbrevis_gvmags.phylogeny_summary.txt` | the summary block quoted above |

**Caveats.**

- The 1,523 in-order references are the `mcl` tier, not every known giant virus. A vOTU's
  position is relative to the genomes included; adding the finer `skani` tier (6,266 taxa)
  could change local topology.
- Only 24 outgroup genomes root a 1,564-taxon tree. That is enough to orient it but thin,
  and the tree is produced **unrooted** regardless - it must be midpoint-rooted before
  interpretation.
- `ncldv_markersearch` pads absent markers with `X`. 44 taxa are >=70% X (section 11h), and
  those are placed on very little data.
- The dropped vOTU05 and vOTU19 are absent from the tree entirely. They are not "unplaced";
  they simply were not analysed, and nothing in the tree output will mention them.

----------------------------------------------------------------------------------------

# Step 11. Maximum-likelihood tree from the alignment (COMPLETE 2026-11-18)

Written by Claude Code on 2026-11-12 (prompt 28 in `PROMPT_LOG.md`). Files:
`scripts/11_gvmag_tree.slurm`, `envs/tree-env.yml`, plus the two diagnostics
`scripts/11_pre_tree_check.slurm` and `scripts/11_pre_tree_check2.slurm` added on
2026-11-16, and `scripts/votu_tree_context.py` added on 2026-11-18. **Job 111564 completed
2026-11-18** in 1 d 11 h 46 m (kuro `serial`, 64 cores, 96 GB; 42.5 GB peak). Sections 11h
and 11i cover the diagnostics that sized the job and the three deliberate deviations from
the paper's IQ-TREE settings; **11j-11l cover the finished tree, where the vOTUs landed, and
the two findings that need following up.**

#### 11a. What this step does and why

Step 10 produces an *alignment*; this step turns it into a *tree*. It follows the paper
(Supplemental Methods, "GVMAG comparisons"), which used IQ-TREE for the tree and RAxML only
for certainty scores.

**IQ-TREE** infers a maximum-likelihood phylogeny: it searches tree topologies for the one
that makes the observed alignment most probable under a model of amino acid substitution.
Three things happen in one invocation:

1. **ModelFinder (`-m MFP`)** tests substitution models and picks the best by BIC. For
   protein data that means empirical matrices (LG, WAG, JTT, rtREV, cpREV, FLU and others)
   combined with `+F` (empirical amino acid frequencies), `+I` (invariant sites) and
   `+G`/`+R` (rate heterogeneity across sites). This matters: in the smoke tests
   ModelFinder chose `LG+F` for one alignment and `rtREV+F+I+R5` for another, so hard-coding
   a model would have been wrong both times.
2. **Tree search** under the chosen model.
3. **Ultrafast bootstrap (`-B 1000`)**, which resamples alignment columns 1,000 times and
   reports how often each branch recurs.

**RAxML `-f i`** then computes **internode certainty (IC)** and **tree certainty (TC)**
scores (Salichos & Rokas 2013) from the bootstrap trees. These expose *conflicting* signal
that a single bootstrap percentage hides: a branch can have high bootstrap support while a
competing arrangement is almost as well supported, and IC reveals that.

#### 11b. SCRIPT: `scripts/11_gvmag_tree.slurm`

```bash
#!/bin/bash
#SBATCH --job-name=11_gvmag_tree
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=3-00:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/11_gvmag_tree_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/11_gvmag_tree_%j.err

# Infer a maximum-likelihood phylogenomic tree from the step 10 trimmed, concatenated
# marker-gene alignment, following Booker et al. 2026 (Supplemental Methods):
#   1. IQ-TREE: model selection + ML tree + ultrafast bootstrap   (paper: -m MFP -bb 1000 -nt AUTO -wbt)
#   2. RAxML -f i: internode certainty (IC) and tree certainty (TC) scores
#
# Flag names below are IQ-TREE 2.x, verified against the installed 2.3.6 help text. The
# paper's `-bb`/`-nt`/`-wbt` are IQ-TREE 1.x spellings; the 2.x equivalents are
# `-B`/`-T`/`--boot-trees` (we use `--wbtl`, which also writes branch lengths because
# RAxML's IC/TC step wants them).
#
# Submit (needs a 72 h partition; kuro caps at 48 h, so add -t 2-00:00:00 there):
#   ( module purge; module load slurm/femto;  sbatch -p batch scripts/11_gvmag_tree.slurm )
#   ( module purge; module load slurm/astral; sbatch -p batch --exclude=as01 scripts/11_gvmag_tree.slurm )
#   ( module purge; module load slurm/kuro;   sbatch -p serial -t 2-00:00:00 scripts/11_gvmag_tree.slurm )
#
# Options:
#   --alignment PATH   alignment to use (default: step 10's kbrevis_gvmags.trimmed.aln)
#   --model STRING     IQ-TREE model (default MFP = ModelFinder Plus, as in the paper).
#                      Use e.g. LG+F+G4 to skip model selection on a large alignment.
#   --bootstraps N     ultrafast bootstrap replicates (default 1000, as in the paper)
#   --threads N|AUTO   default: the Slurm allocation (the paper used AUTO)
#   --alrt N           also run SH-aLRT with N replicates (default 0 = off)
#   --skip-ic          skip the RAxML IC/TC step

set -euo pipefail

THREADS_ALLOC="${SLURM_CPUS_PER_TASK:?Submit with sbatch}"
REPO="/sciclone/home/iadanforth/K_brevis_viromes"
DATA="$HOME/scr10/data/kbrevis"
PHYLO_DIR="$DATA/phylogeny"
OUT="$DATA/tree"
PROJECT="kbrevis_gvmags"

ALN="$PHYLO_DIR/${PROJECT}.trimmed.aln"
MODEL=MFP; BOOTSTRAPS=1000; THREADS="$THREADS_ALLOC"; ALRT=0; DO_IC=1
while (($#)); do
    case "$1" in
        --alignment) ALN="$2"; shift 2 ;;
        --model) MODEL="$2"; shift 2 ;;
        --bootstraps) BOOTSTRAPS="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --alrt) ALRT="$2"; shift 2 ;;
        --skip-ic) DO_IC=0; shift ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done

[[ -s "$ALN" ]] || { echo "No trimmed alignment ($ALN); run step 10 first." >&2; exit 1; }
mkdir -p "$OUT"

module load miniforge3
eval "$(conda shell.bash hook)"
conda activate tree-env

# ---- Pre-flight: look at the alignment before spending days on it -------------
# IQ-TREE refuses to bootstrap fewer than 4 sequences ("It makes no sense to perform
# bootstrap with less than 4 sequences"), and it only says so AFTER model selection has
# run. Catch it up front. Found while smoke-testing on a 3-taxon alignment.
N_TAXA=$(grep -c '^>' "$ALN")
if ((N_TAXA < 4)) && ((BOOTSTRAPS > 0)); then
    echo "Only $N_TAXA taxa in $ALN. IQ-TREE cannot bootstrap fewer than 4 sequences." >&2
    echo "Either add reference genomes (step 10 without --no-refs) or pass --bootstraps 0." >&2
    exit 1
fi

# Taxa made up mostly of X padding (markers that were not found) are placed unreliably and
# slow everything down. This reports them so they can be removed and the step re-run; it
# does not drop them automatically, since that is a scientific decision.
python - "$ALN" <<'PY'
import sys
aln, name, seqs = sys.argv[1], None, {}
for line in open(aln):
    line = line.rstrip("\n")
    if line.startswith(">"):
        name = line[1:].split()[0]; seqs[name] = []
    elif name:
        seqs[name].append(line.strip())
seqs = {k: "".join(v) for k, v in seqs.items()}
if not seqs:
    sys.exit("alignment contains no sequences")
L = len(next(iter(seqs.values())))
bad = {k: v for k, v in seqs.items() if len(v) != L}
print(f"[tree] {len(seqs)} taxa, {L} alignment columns")
if bad:
    sys.exit(f"alignment is not rectangular: {len(bad)} sequences differ in length "
             f"(first: {next(iter(bad))})")
def pad(s):
    return 100.0 * (s.count("X") + s.count("x") + s.count("-")) / len(s)
poor = sorted(((pad(v), k) for k, v in seqs.items()), reverse=True)
over = [(p, k) for p, k in poor if p >= 70]
print(f"[tree] worst padding/gap fraction: {poor[0][0]:.1f}% ({poor[0][1]})")
if over:
    print(f"[tree] WARNING: {len(over)} taxa are >=70% X/gap and are placed unreliably:")
    for p, k in over[:10]:
        print(f"[tree]   {p:5.1f}%  {k}")
    print("[tree] Consider removing them (see phylogeny/*.table.tsv) and re-running step 10.")
PY

# ---- 1. IQ-TREE: ModelFinder + ML tree + ultrafast bootstrap ------------------
# -m MFP  : ModelFinder Plus - tests substitution models (for protein data: LG, WAG, JTT,
#           ... with +F/+G/+R variants), picks the best by BIC, then infers the tree.
# -B 1000 : ultrafast bootstrap (UFBoot), 1000 replicates. UFBoot support values are NOT
#           comparable to classical bootstrap: >=95 is roughly as strong as >=70 classical.
# --bnni  : optimise each UFBoot tree by NNI. Recommended by the IQ-TREE authors because it
#           reduces UFBoot's tendency to overestimate support when the model is misspecified
#           - a real risk with concatenated markers. (Not in the paper; a deliberate addition.)
# --wbtl  : write the bootstrap trees, with branch lengths, to <prefix>.ufboot. The paper's
#           -wbt did the same without branch lengths; RAxML's IC/TC step reads this file.
# -T      : explicit thread count rather than the paper's AUTO, because AUTO probes the
#           machine and can exceed the Slurm allocation.
PREFIX="$OUT/$PROJECT"
# --mem caps IQ-TREE's own memory use at 85% of the Slurm allocation, so it reports a
# shortfall rather than being OOM-killed. (On small alignments IQ-TREE prints harmless
# "Too low -mem, automatically increased to N MB" notes per model; ignore them.)
IQ_ARGS=(-s "$ALN" -m "$MODEL" --bnni
         -T "$THREADS" --prefix "$PREFIX" --seed 12345
         --mem "$(( ${SLURM_MEM_PER_NODE:-131072} * 85 / 100 / 1024 ))G")
((BOOTSTRAPS > 0)) && IQ_ARGS+=(-B "$BOOTSTRAPS" --wbtl)
((ALRT > 0)) && IQ_ARGS+=(--alrt "$ALRT")

echo "[tree] iqtree2 ${IQ_ARGS[*]}"
iqtree2 "${IQ_ARGS[@]}"

TREE="$PREFIX.treefile"
[[ -s "$TREE" ]] || { echo "IQ-TREE produced no tree ($TREE)" >&2; exit 1; }

# ---- 2. RAxML -f i: internode certainty / tree certainty ---------------------
# -f i computes IC and TC scores (Salichos & Rokas 2013): how consistently each internal
# branch is supported across the bootstrap trees, which exposes conflicting signal that a
# single bootstrap number can hide.
#
# DELIBERATE DEVIATION FROM THE PAPER: it reports `-f i -m GTRCAT`, but GTRCAT is RAxML's
# GTR *nucleotide* model and this alignment is amino acids. RAxML's protein models are
# PROTCAT<matrix>/PROTGAMMA<matrix>, so PROTCATLG is used here (LG is the modern general
# protein matrix). The IC/TC calculation compares tree topologies rather than optimising a
# likelihood, so the matrix choice has little effect - but the model must at least match the
# data type or RAxML will reject the alignment.
if ((DO_IC)); then
    UFBOOT="$PREFIX.ufboot"
    if [[ -s "$UFBOOT" ]]; then
        IC_DIR="$OUT/ic"; rm -rf "$IC_DIR"; mkdir -p "$IC_DIR"
        # RAxML writes RAxML_* files into -w, which must be an absolute path.
        raxmlHPC-PTHREADS-SSE3 -T "$THREADS" -f i \
            -t "$TREE" -z "$UFBOOT" -s "$ALN" \
            -m PROTCATLG -n ic -w "$IC_DIR" -C \
            || echo "[tree] WARNING: RAxML IC/TC step failed; the ML tree above is unaffected" >&2
        cp "$IC_DIR"/RAxML_info.ic "$OUT/${PROJECT}.IC_TC_info.txt" 2>/dev/null || true
    else
        echo "[tree] no .ufboot file, so the IC/TC step was skipped" >&2
    fi
fi

# ---- Summary -----------------------------------------------------------------
{
  echo "project:            $PROJECT"
  echo "date:               $(date)"
  echo "alignment:          $ALN"
  echo "model requested:    $MODEL"
  echo "model chosen:       $(grep -m1 'Best-fit model' "$PREFIX.log" 2>/dev/null | sed 's/.*: //' || echo 'see log')"
  echo "log-likelihood:     $(grep -m1 'BEST SCORE FOUND' "$PREFIX.log" 2>/dev/null | sed 's/.*: //' || echo 'see log')"
  echo "bootstrap reps:     $( ((BOOTSTRAPS>0)) && echo "$BOOTSTRAPS (ultrafast UFBoot, --bnni)" || echo "none (--bootstraps 0)" )"
  echo "threads:            $THREADS"
  echo "taxa in tree:       $(tr ',' '\n' < "$TREE" | grep -c ':' || true)"
  echo
  if ((BOOTSTRAPS > 0)); then
  echo "UFBoot support distribution (branch support values):"
  python - "$TREE" <<'PY'
import re, sys
vals = [float(x) for x in re.findall(r"\)(\d+(?:\.\d+)?):", open(sys.argv[1]).read())]
if not vals:
    print("  (no support values found in the tree file)")
else:
    vals.sort()
    n = len(vals)
    strong = sum(1 for v in vals if v >= 95)
    print(f"  {n} internal branches; median {vals[n//2]:.0f}")
    print(f"  >=95 (well supported): {strong} ({100*strong/n:.0f}%)")
    print(f"  <70  (poorly supported): {sum(1 for v in vals if v < 70)}")
PY
  fi
} | tee "$OUT/${PROJECT}.tree_summary.txt"

echo
echo "[tree] done. Key outputs in $OUT:"
echo "  ${PROJECT}.treefile          the ML tree with UFBoot support (Newick)"
echo "  ${PROJECT}.iqtree            human-readable report: model, likelihood, parameters"
echo "  ${PROJECT}.log               full IQ-TREE log"
[[ -s "$PREFIX.ufboot" ]] && echo "  ${PROJECT}.ufboot            the bootstrap trees"
[[ -s "$OUT/${PROJECT}.IC_TC_info.txt" ]] && echo "  ${PROJECT}.IC_TC_info.txt    RAxML internode/tree certainty scores"
echo "  ${PROJECT}.tree_summary.txt  this summary"
echo
echo "The tree is UNROOTED. Booker et al. midpoint-rooted theirs; do that in iTOL,"
echo "FigTree or PearTree, or with ape::midpoint in R, before interpreting it."
```

#### 11c. Parameters, and where they differ from the paper

The paper's command was `-m MFP -bb 1000 -nt AUTO -wbt`. Those are IQ-TREE **1.x** flag
spellings; the installed version is 2.3.6, whose equivalents I verified against its own
help text:

| paper (IQ-TREE 1.x) | used here (2.3.6) | meaning |
|---|---|---|
| `-m MFP` | `-m MFP` | unchanged: ModelFinder Plus |
| `-bb 1000` | `-B 1000` | ultrafast bootstrap, 1,000 replicates |
| `-nt AUTO` | `-T $SLURM_CPUS_PER_TASK` | **deliberate change**: `AUTO` probes the machine and can exceed the Slurm allocation, so the thread count is set explicitly |
| `-wbt` | `--wbtl` | write bootstrap trees; `--wbtl` also writes branch lengths, which the RAxML step wants |
| – | `--bnni` | **addition**: optimises each bootstrap tree by NNI. The IQ-TREE authors recommend it because UFBoot overestimates support when the model is misspecified, a real risk with concatenated markers |
| – | `--seed 12345` | reproducible runs |
| – | `--mem` | caps IQ-TREE at 85% of the Slurm allocation so it reports a shortfall rather than being OOM-killed |

**The RAxML model is a genuine correction to the paper.** The paper reports
`-f i -m GTRCAT`, but **GTRCAT is RAxML's GTR nucleotide model** and this alignment is amino
acids - RAxML's protein models are `PROTCAT<matrix>` / `PROTGAMMA<matrix>`. This script uses
**`PROTCATLG`** (LG being the modern general protein matrix). The IC/TC calculation compares
tree topologies rather than optimising a likelihood, so the matrix barely affects the
numbers, but the model must match the data type or RAxML rejects the alignment. The paper's
published combination cannot have been run as written on a protein alignment.

**Reading the support values.** UFBoot values are **not** on the same scale as classical
bootstrap: **>=95 is roughly as trustworthy as >=70 classical**. The script's summary
reports how many branches clear 95 and how many fall below 70.

#### 11d. Environment

`envs/tree-env.yml`:

```yaml
# Tree inference from the step 10 trimmed alignment (step 11).
#
# iqtree 2.3.6: the paper used IQ-TREE v2.2.2.6, which is not on bioconda (available 2.x
# releases are 2.2.3, 2.2.5, 2.2.6, 2.3.0-2.3.6). 2.3.6 is the newest 2.x, so it keeps the
# paper's command-line interface. IQ-TREE 3.x exists but is a major release, so it is
# avoided here for comparability with the published tree.
#
# raxml 8.2.12: exactly the version the paper used, needed only for the internode/tree
# certainty scores (-f i), which IQ-TREE does not compute.
name: tree-env
channels: [conda-forge, bioconda]
dependencies:
  - python=3.12   # the job script's inline alignment and support-value checks
  - iqtree=2.3.6
  - raxml=8.2.12
```

IQ-TREE **2.3.6**: the paper used 2.2.2.6, which is not on bioconda (available 2.x releases
are 2.2.3, 2.2.5, 2.2.6 and 2.3.0-2.3.6), so this is the newest 2.x and keeps the same
command-line interface. IQ-TREE 3.x exists but is a major release, avoided here for
comparability. RAxML **8.2.12** is exactly the paper's version. `python` is included because
the job script uses it for the pre-flight alignment check and the support-value summary -
the environment did not have Python at first, and the script failed until it was added.

#### 11e. How to run it

```
# default: the step 10 alignment, ModelFinder, 1000 UFBoot replicates
( module purge; module load slurm/femto;  sbatch -p batch scripts/11_gvmag_tree.slurm )

# kuro caps walltime at 48 h, so shorten it there
( module purge; module load slurm/kuro;   sbatch -p serial -t 2-00:00:00 scripts/11_gvmag_tree.slurm )

# skip model selection on a very large alignment (much faster)
( module purge; module load slurm/femto;  sbatch -p batch scripts/11_gvmag_tree.slurm --model LG+F+G4 )

# add SH-aLRT as a second support measure, and skip the RAxML step
( module purge; module load slurm/femto;  sbatch -p batch scripts/11_gvmag_tree.slurm --alrt 1000 --skip-ic )
```

Resources: 32 cores, 128 GB, 3 days. A ~1,000-taxon protein alignment with ModelFinder and
1,000 UFBoot replicates is genuinely expensive; if it times out, rerun with an explicit
`--model` to skip model selection, or cut taxa with `--max-refs` in step 10.

Outputs in `~/scr10/data/kbrevis/tree/`:

| file | contents |
|---|---|
| `kbrevis_gvmags.treefile` | **the ML tree with UFBoot support values (Newick)** |
| `kbrevis_gvmags.iqtree` | readable report: chosen model, likelihood, parameters |
| `kbrevis_gvmags.log` | full IQ-TREE log, including the whole ModelFinder table |
| `kbrevis_gvmags.contree` | consensus tree of the bootstrap replicates |
| `kbrevis_gvmags.ufboot` | the 1,000 bootstrap trees |
| `kbrevis_gvmags.IC_TC_info.txt` | RAxML IC/TC summary |
| `ic/RAxML_IC_Score_BranchLabels.ic` | the tree with IC scores as branch labels |
| `kbrevis_gvmags.tree_summary.txt` | model, likelihood, support distribution |

#### 11f. How this was tested

All tests ran as Slurm jobs on astral, never on dev.

**Test 1 (job 73846, 3 taxa, FAILED - a useful failure).** Run on the step 10 smoke-test
alignment. ModelFinder completed and chose `LG+F`, then IQ-TREE stopped with
`ERROR: It makes no sense to perform bootstrap with less than 4 sequences.` The pre-flight
check correctly flagged `medusavirus` as 84.8% X/gap - the taxon that had only 3 of the 7
markers. Fix: the script now counts taxa up front and refuses early with a clear message,
rather than discovering this after model selection.

**Test 2 (job 73847, 3 taxa, `--bootstraps 0`, COMPLETED in 2 s).** Confirmed the no-
bootstrap path: model `LG+F`, log-likelihood -17659.226, tree written, and the IC/TC step
correctly skipped with a message because no `.ufboot` file exists. This run also exposed
three misleading summary lines (reporting "1000 (ultrafast)" when bootstrapping was off,
"could not parse support values", and listing files that were never created), all fixed.

**Test 3 (job 73848, FAILED - environment lesson).** The test alignment was passed from the
agent scratchpad under `/tmp`, which is **node-local**, so the compute node could not see it:
`No trimmed alignment (...); run step 10 first.` Re-staged to `~/scr10` and resubmitted.

**Test 4 (job 73849, 6 taxa, COMPLETED in 17 s)** - the full path. A 6-taxon alignment was
built from the 3 real genomes plus mutated copies, so bootstrapping and RAxML could actually
run:

```
model chosen:       rtREV+F+I+R5 chosen according to BIC
log-likelihood:     -20020.086
bootstrap reps:     1000 (ultrafast UFBoot, --bnni)
taxa in tree:       6
UFBoot support distribution (branch support values):
  3 internal branches; median 100
  >=95 (well supported): 3 (100%)
  <70  (poorly supported): 0
```

and RAxML produced real certainty scores on protein data with `PROTCATLG`:

```
Tree certainty for this tree: 2.947085
Relative tree certainty for this tree: 0.982362
```

So the verified-working pieces are: ModelFinder on protein data, the `-B`/`--wbtl`/`--bnni`
flags, the thread and memory handling, the UFBoot support parser, the RAxML `-f i` call with
a protein model, and the summary. The 100% support is an artefact of the test (each mutated
copy pairs with its original), not a meaningful result.

#### 11g. Caveats

- **The tree is unrooted.** The paper midpoint-rooted theirs. Root it in iTOL, FigTree or
  PearTree, or with `ape::midpoint()` in R, before interpreting it. The script prints a
  reminder.
- **UFBoot is not classical bootstrap.** Treat >=95 as the threshold, not >=70.
- Taxa that are mostly `X` padding get placed unreliably. The pre-flight lists any taxon
  that is >=70% X/gap; consider removing them and rerunning step 10.
- IQ-TREE silently merges identical sequences and warns about it - check the log if the taxon
  count in the tree is lower than in the alignment.
- ModelFinder on a large alignment can take a large share of the total runtime. `--model
  LG+F+G4` skips it.
- A tree of concatenated markers is a tree of *those genes*, not necessarily of the genomes:
  giant viruses recombine and acquire genes horizontally, so individual marker trees can
  disagree. That is exactly what the IC/TC scores are for.


# Earlier attempts (not part of the current pipeline)

Kept so the history makes sense. Their outputs either no longer exist or aren't used.

- **BBDuk phiX cleaning (`clean_files.sh/.slurm`, env `bbmpa-env`).** My first plan copied
  the paper's BBDuk command
  (`ref=phix k=27 hdist=1 qtrim=rl trimq=17 cardinality=t mingc=0.05 maxgc=0.95 -Xmx30g`)
  into `clean_files.sh`, tested on a 2,500-read subset (`head -n 10000` of each raw FASTQ
  into `~/scr10/data/kbrevis/test_data/`). Replaced by fastp (step 02) + bowtie2 screening
  (step 04); the scripts are no longer in `scripts/`.
- **Assembling a *K. brevis* reference from SRR14099984** (`01_/02_/03_assemble_kbrevis.*`).
  The idea was to trim the culture reads with fastp, remove mammalian reads with bowtie2,
  and assemble them with metaSPAdes to get a *K. brevis* reference. Only the host-screening
  part ran: its output is in `~/scr10/data/kbrevis/host_filtered/` (5.03 % of reads matched
  the mammalian references). Dinoflagellates like *K. brevis* have enormous genomes (tens of
  gigabases or more), so assembling one from a single short-read run isn't feasible. Abandoned in favour of using the raw
  reads directly as a reference (step 03).
- **Read correction with bfc (`clean_reads`, env `readcorr-env`).** Job 578087 (bora, 8
  cores / 64 GB, 10.5 h, COMPLETED) interleaved read pairs and ran
  `bfc -1 -s 10g -k 21 -t 8`, the paper's error-correction step, into
  `~/scr10/data/kbrevis/corrected/`. That folder and the script are gone, and the current
  pipeline assembles the bowtie2-cleaned reads without bfc correction (a deviation from the
  paper; SRA reads were probably already corrected by the authors).
- **Early metaSPAdes runs** wrote to `~/K_brevis_viromes/output/contigs/` (in home, against
  the cluster rules). Output deleted (`output/contigs/SRR34748664/spades.log` was removed
  from git).
- **`.yml` exports in the repo root.** `sra-toolkit-env.yml` and `datasets-env.yml` were
  written there by the `.slurm` scripts and later deleted; current exports are in `envs/`.

# Known problems in the scripts (not yet fixed)

| script | problem | fix |
|---|---|---|
| `01_download_files.slurm` | fasterq-dump ignores the 4 cores; 120 GB is far more than needed | `--threads "$SLURM_CPUS_PER_TASK"`; `--mem=16G` |
| `02_fastp.slurm` | relies on `PATH` from `~/.bashrc`; no `--thread` | activate an env with fastp; add `--thread "$SLURM_CPUS_PER_TASK"` |
| `03_download_refs.sh` | zip clean-up uses bare file names; `>>` appends on re-run | use `"$REF_DIR"/...zip`; use `>` |
| `05_metaspades.slurm` | `--mail-type`, `--output`, `--error` lines duplicated; `BEGIN` mail for 20 jobs | delete the duplicates; `--mail-type=FAIL,END` |
| most `.slurm` | `--mail-type=...BEGIN...` | `FAIL,END` |

# Rebuilding everything from scratch (order of commands)

**Audited end to end on 2026-11-19 (prompt 50 in `PROMPT_LOG.md`)** against the actual
contents of `~/scr10/data/kbrevis` and `~/K_brevis_viromes`, to confirm that a new user
cloning the repository can reproduce both. The gaps found are listed in section R below;
all of them are fixed in the commands here.

Run everything **from the repository root on dev**. This is not a style preference: the
step 01-03 `.slurm` wrappers invoke their worker script by a relative path
(`bash ./scripts/02_fastp.sh`), and `sbatch` starts a job in the directory it was
submitted from, so submitting from anywhere else fails with "No such file or directory".

Wait for each step's jobs to finish before starting the next. Before each submission,
check which cluster will finish soonest with `~/hpc/best_cluster.sh -c CORES -m MEM_GB
-H HOURS` and override `-p` (and `-t`) on the `sbatch` line; the partitions below are
the ones that worked on the dates given, not permanent recommendations.

### 0. Prerequisites a clone does not provide

```bash
bash scripts/00_directory_setup.sh                     # creates the scratch and repo folders
cp metadata/SraRunTable.csv ~/scr10/data/kbrevis/      # step 01 reads it from scratch
( module purge; module load slurm/astral; sbatch -p batch --exclude=as01 scripts/00_create_envs.slurm )
```

`00_create_envs.slurm` builds the six environments steps 01-06 need from `envs/*.yml`
and then verifies that every tool those steps call is actually on `PATH`. Steps 07 and 14
build their own environments later, so they are not duplicated here.

`fastp` is included: it was moved out of a manual `~/programs/fastp` binary and into
`envs/fastp-env.yml` (pinned to the same version, 1.3.7) on 2026-11-19, which removed the
pipeline's last dependency on a shell-specific `PATH`. Nothing outside the repository and
these environment files is needed any more.

### 1-6. Reads, references, cleaning, assembly, merge

```bash
( module purge; module load slurm/bora;   sbatch scripts/01_download_files.slurm )      # ~2.0 TB into raw/
( module purge; module load slurm/bora;   sbatch scripts/02_fastp.slurm )               # -> trim/ (~1.6 TB)
( module purge; module load slurm/bora;   sbatch scripts/03_download_refs.slurm )       # -> reference/  (SEE BELOW)
( module purge; module load slurm/bora;   sbatch scripts/03b_download_phix.slurm )      # -> reference/phix174.fna
( module purge; module load slurm/bora;   sbatch -p hima scripts/04_clean_refs.slurm )  # -> clean/ (~1.1 TB)
( module purge; module load slurm/kuro;   bash scripts/05_metaspades.sh )               # 20 jobs -> contigs/
( module purge; module load slurm/kuro;   bash scripts/05_metaspades_nov.sh )           # 20 jobs -> contigs_no_v/
( module purge; module load slurm/kuro;   bash scripts/06_merge_assemblies.sh )         # 20 jobs -> contigs_merged/
```

**Skipping step 03 on a re-run.** Steps 03 and 03b download ~79 GB of host and contaminant
genomes and take a long time. They only need to run once. To keep them:

```bash
# keep ~/scr10/data/kbrevis/reference/ when clearing scratch, then skip 03 and 03b entirely
```

This is safe, and the reason is worth stating. Step 04 needs six FASTAs in `reference/`
(`dog_ref_genome.fna`, `cat_ref_genome.fna`, `mouse_ref_genome.fna`,
`human_ref_genome.fna`, `phix174.fna`, `kbrevis_ref.fna`) and exits with a clear error
naming any that are missing, so a partial `reference/` fails loudly rather than silently
producing badly cleaned reads. The 42 GB bowtie2 index (`host_index.*.bt2l`, 6 files) also
lives in `reference/`, and `04_clean_refs.sh` rebuilds it **only if absent**
(`if ! bowtie2-inspect "$INDEX_PREFIX"`), so keeping `reference/` also skips several hours
of `bowtie2-build`. Verified present and complete on 2026-11-19.

### 7-11. GVMAGs, vOTUs, phylogeny, tree

```bash
( module purge; module load slurm/astral; sbatch -p batch -c 4 --exclude=as01 --mem=16G -t 6:00:00 \
    --export=ALL,GVMAGS_V2_FULL=1 scripts/07_gvmag_setup.slurm )      # envs, GVClass, CheckV/geNomad/GVMAGs V2 dbs
( module purge; module load slurm/kuro;   bash scripts/08_gvmags.sh -p serial )          # array -> gvmags/
( module purge; module load slurm/kuro;   sbatch -p serial -c 16 --mem=64G -t 4:00:00 \
    scripts/09_gvmag_votus.slurm --with-refs )                                           # -> votus/
( module purge; module load slurm/kuro;   sbatch -p serial -c 32 --mem=16G -t 4:00:00 \
    scripts/10_gvmag_phylogeny.slurm )                                                   # -> phylogeny/
( module purge; module load slurm/kuro;   sbatch -p serial scripts/11_pre_tree_check.slurm )   # optional: sizes step 11
( module purge; module load slurm/kuro;   sbatch -p serial scripts/11_pre_tree_check2.slurm )  # optional: X vs gap
( module purge; module load slurm/kuro;   sbatch -p serial -c 64 --mem=96G -t 48:00:00 \
    scripts/11_gvmag_tree.slurm --mset LG,WAG,JTT --cmax 4 --safe )                      # -> tree/
```

Note `--export=ALL,GVMAGS_V2_FULL=1` on step 07: without it the 1.18 GB GVMAGs V2 genome
bundle is not downloaded and **step 10 exits with an error**, because it needs the
reference proteomes. The original run discovered this the slow way.

Step 10's `--mem=16G` is the measured requirement (3.6 GB peak); the `#SBATCH` line in the
script still says 128 GB and should be overridden. Step 11's flags are explained in
sections 11h-11i: unrestricted ModelFinder would consume the whole walltime before building
a tree. If step 11 hits the 48 h limit, **resubmit the identical command** - IQ-TREE resumes
from its checkpoint.

### 12. Collect results into the repository

```bash
( module purge; module load slurm/kuro; sbatch -p serial scripts/copy_results.slurm )
```

Copies `gvmags/`, `votus/` and the step 11 tree files into `output/`, excluding the large
regenerable intermediates (`binning/`, `gvclass_in/`, `.ufboot`, `.ckp.gz`, `.mldist`).

### 13. iTOL annotations and vOTU tree placement

These are small, pure-Python, and run on dev:

```bash
python3 -I scripts/make_itol_annotations.py \
    --tree     output/tree/kbrevis_gvmags.treefile \
    --manifest ~/scr10/data/kbrevis/phylogeny/reference_manifest.tsv \
    --votus    output/votus/gvmags_votus.tsv \
    --nearest  output/votus/votu_nearest_reference.tsv \
    --out-dir  output/tree

python3 -I scripts/votu_tree_context.py \
    --tree     output/tree/kbrevis_gvmags.treefile \
    --manifest ~/scr10/data/kbrevis/phylogeny/reference_manifest.tsv \
    --out      output/tree/votu_tree_context.tsv
```

`output/Kbrevis_virome_tree.png` **was** a hand-exported iTOL image with no generating
script. **Superseded 2026-11-24 (section 20):** it is now drawn by
`scripts/plot_tree_figure.py` directly from `output/tree/kbrevis_gvmags.treefile`, so it is
reproducible and its type size and line weights can be changed without a round trip through
a web tool. The iTOL annotation files are still written by `make_itol_annotations.py` and
still work if an interactive tree is wanted.

### 14. vOTU abundance in the environment through time

```bash
( module purge; module load slurm/astral; sbatch -p batch --exclude=as01 scripts/14_abundance_setup.slurm )
bash scripts/14_votu_reference.sh                       # 19 vOTU representatives -> abundance/votu_ref
( module purge; module load slurm/astral; sbatch -p batch --exclude=as01 scripts/14_abundance_smoketest.slurm )  # optional but recommended
( module purge; module load slurm/bora;   sbatch -p hima -C hi --exclude=hi04,hi05,hi07 \
    scripts/14_votu_abundance.slurm )                   # array of 11 -> abundance/per_sample

python3 -I scripts/merge_votu_abundance.py \
    --coverm-dir ~/scr10/data/kbrevis/abundance/per_sample \
    --manifest   ~/scr10/data/kbrevis/abundance/votu_reference_manifest.tsv \
    --sra-table  metadata/SraRunTable.csv \
    --out-long   output/abundance/votu_abundance_long.tsv \
    --out-matrix output/abundance/votu_abundance_matrix.tsv

module load miniforge3 && eval "$(conda shell.bash hook)" && conda activate r-plot-env
Rscript scripts/plot_votu_abundance.R output/abundance/votu_abundance_long.tsv output/abundance
conda deactivate
```

`14_votu_reference.sh` must run before the array; the array checks for its output and exits
with an instruction if it is missing.

### 15. ANI against the published GVMAGs

```bash
bash scripts/15_fetch_booker_gvmags.sh                  # downloads 74 GenBank contigs -> 10 genomes
( module purge; module load slurm/astral; sbatch -p batch --exclude=as01 scripts/15_booker_ani.slurm )
```

Writes the three CSVs in `output/` directly, so no separate collection step is needed.

### R. What the 2026-11-19 reproducibility audit found

| # | finding | status |
|---|---|---|
| 1 | The runbook stopped at step 11. Steps 12, 14 and 15 and every manual post-processing command were undocumented, so a new user would have reproduced the scratch directory but not most of `output/` | **fixed** - this section now covers 00-15 |
| 2 | `output/tree/` was populated by an ad-hoc `cp` loop typed into a terminal, not by any script, so a re-run would not recreate it | **fixed** - `copy_results.slurm` now copies the named tree files |
| 3 | Steps 01-06 assumed their conda environments already existed; nothing created them | **fixed** - new `scripts/00_create_envs.slurm` |
| 4 | `03_download_refs.slurm` wrote `datasets-env.yml` into the **repository root** on every run, leaving a stray duplicate of `envs/datasets-env.yml` | **fixed** - writes `envs/datasets-env.lock.yml`, with the machine-specific `prefix:` line stripped |
| 5 | Root-level `datasets-env.yml` and `sra-toolkit-env.yml` were tracked duplicates of the `envs/` copies | **fixed** - untracked |
| 6 | `fastp` was a standalone binary on `PATH` via `~/.bashrc`, not in any `envs/*.yml`. It worked here only because jobs inherit the submitting shell; a new account would have failed at step 02 | **fixed 2026-11-19** - now `envs/fastp-env.yml` pinned to the same 1.3.7, activated by `02_fastp.slurm` and built by `00_create_envs.slurm` |
| 7 | `output/Kbrevis_virome_tree.png` is a hand-exported iTOL image with no generating script | **fixed 2026-11-24** - `scripts/plot_tree_figure.py` now draws it from the treefile (section 20) |
| 8 | Step 07 needs `--export=ALL,GVMAGS_V2_FULL=1` or step 10 fails later | **documented** in the step 7-11 block |
| 9 | **23 scripts hardcode `/sciclone/home/iadanforth/`** (59 occurrences) and 17 hardcode `--mail-user=iadanforth@wm.edu` | **documented** below - fine for re-running on this account, a blocker for anyone else |

#### Porting to a different account

Everything above assumes the pipeline is being re-run **on this account**, which is what
the audit was asked to confirm. A genuinely different user cloning from GitHub hits one
more obstacle, and it is worth being precise about it rather than claiming the repository
is portable when it is not.

23 of the 44 scripts contain the literal path `/sciclone/home/iadanforth/`, 59 times.
**44 of those occurrences are unavoidable**: AGENTS.md requires every job's `.out`/`.err`
to land in `<project>/logs/`, and `#SBATCH` directives cannot expand `~` or `$HOME`, so the
absolute path has to be written out. The rest are `REPO=` assignments which could use
`$HOME` but do not. A further 17 scripts carry `--mail-user=iadanforth@wm.edu`, which would
send another user's job notifications to this account.

For another user on SciClone, two substitutions fix all of it:

```bash
cd <their clone>
grep -rl '/sciclone/home/iadanforth' scripts/ \
  | xargs sed -i 's|/sciclone/home/iadanforth|/sciclone/home/<their-username>|g'
grep -rl 'iadanforth@wm.edu' scripts/ \
  | xargs sed -i 's|iadanforth@wm.edu|<their-address>@wm.edu|g'
```

Off SciClone the scripts need more than substitution: the `module load slurm/<cluster>`
lines, the partition names, `/local/scr/$USER` and `~/scr10` / `~/data10` are all specific
to this system. The *analysis* is portable - the tools, parameters and thresholds are all
documented above - but the *job scripts* are not, and rewriting them for another cluster is
a real piece of work rather than a find-and-replace.

Checked and found **already correct**, so left alone:

- Every `.sh` worker script creates the scratch directory it writes to (`mkdir -p`), so the
  short `00_directory_setup.sh` is not a problem despite listing only `raw/` and `clean/`.
- All 44 scripts, all 25 `envs/*.yml` and `metadata/SraRunTable.csv` are committed.
- `scripts/__pycache__/`, `scripts/select_refs.py.bak` and `scripts/11_gvmag_tree.slurm.orig`
  are correctly ignored and absent from a clone.
- Every artifact in `output/` is produced by a committed script, with the single exception
  of the hand-exported iTOL PNG above.

----------------------------------------------------------------------------------------

#### 11h. Pre-flight diagnostics before committing to a multi-day tree (2026-11-16)

Written by Claude Code on 2026-11-16 (prompt 42 in `PROMPT_LOG.md`). Files:
`scripts/11_pre_tree_check.slurm`, `scripts/11_pre_tree_check2.slurm`.

Step 11 is the most expensive job in the pipeline, and the script as written asks for
32 cores, 128 GB and **3 days**. Guessing those numbers for a 1,564-taxon x 37,739-column
amino-acid alignment is a bad bet in both directions: too little memory and the job dies
after days of work, too much and it sits in the queue blocking other users for nothing.
IQ-TREE reports the RAM it needs within seconds of reading an alignment, so two short
diagnostic jobs were run first. Together they cost 26 minutes and changed every resource
number in the submission.

**Why a separate script rather than running it on dev:** reading a 58 MB alignment and
starting IQ-TREE is exactly the kind of work AGENTS.md section 1 forbids on the front end.
A first attempt put the probe script in the session scratchpad under `/tmp` and failed
instantly with `No such file or directory` - `/tmp` is **node-local**, so a compute node
cannot see a file dev wrote there. Diagnostic scripts belong in the repository (or scratch),
like any other script a job runs.

**Diagnostic 1 (`11_pre_tree_check.slurm`, job 111549, 26 m 12 s, 33.7 GB peak).** Runs
`iqtree2 -m LG+F+G4 -n 0`, which initialises and optimises parameters but does not search
for a tree. The decisive lines:

```
Alignment has 1564 sequences with 37739 columns, 32764 distinct patterns
21156 parsimony-informative, 8598 singleton sites, 7985 constant sites
WARNING: 566 sites contain only gaps or ambiguous characters.
WARNING: 1564 sequences contain more than 50% gaps/ambiguity
NOTE: 31686 MB RAM (30 GB) is required!
WARNING: Numerical underflow for lh-branch, use `-safe` option to avoid this warning
Parameters optimization took 13 rounds (1491.654 sec)
```

Four findings, each of which changed the submission:

1. **30 GB of RAM is required, not 128 GB.** The script's request was more than 4x too
   large. 96 GB was requested in the end (headroom for 64 threads, whose per-thread buffers
   the 4-thread probe did not exercise), which still leaves most of a kuro node free.
2. **21,156 of 37,739 columns are parsimony-informative.** Despite the 14.5% average
   identity, there is real signal here - this is not a hopeless alignment.
3. **Optimising ONE model's parameters took 1,492 seconds (25 minutes).** This is the
   finding that mattered most. The paper used `-m MFP`, which asks ModelFinder to test every
   protein matrix against every rate-heterogeneity variant - on the order of 80 models. At
   this cost per model, unrestricted ModelFinder would spend **days** before the tree search
   began, and the walltime would expire during model selection.
4. **IQ-TREE explicitly asked for `-safe`.** Numerical underflow at 14.5% identity across
   1,564 taxa is expected, and an underflowing likelihood can mislead both the tree search
   and the support values.

**Diagnostic 2 (`11_pre_tree_check2.slurm`, job 111550, 6 seconds).** The first diagnostic
reported that **all 1,564 taxa are >=70% "padded"**, which is true and useless - it is also
what step 11's own built-in pre-flight warning says, and it fires for every taxon including
perfectly good ones. The reason is that the metric added together two completely different
things:

- **`X`** - a marker gene **not found** in that genome. This is missing data, and a taxon
  made mostly of `X` cannot be placed.
- **`-`** - an alignment gap, which says nothing about the genome. It means Clustal Omega
  inserted a gap because some *other* sequence has residues there. At 14.5% identity across
  1,564 divergent genomes the alignment is mostly gaps **for everyone**, complete genomes
  included.

Separating them gives the picture that actually matters:

```
=== MISSING DATA (X only) ===
X% histogram (bucket: taxa): {0: 762, 10: 190, 20: 167, 30: 160, 40: 104,
                              50: 91, 60: 46, 70: 33, 80: 11}
  >=50% X: 181 taxa      >=70% X: 44 taxa      >=90% X: 0 taxa

=== ALIGNMENT GAPS (- only) ===
  median 76.4%   min 11.0%   max 90.8%

=== REAL RESIDUES per taxon ===
  median 11.7%  -> ~4413 amino acids of actual data
```

So the honest reading is: **762 taxa have no missing markers at all**, 44 are >=70% missing
data, and none is >=90%. The median taxon contributes ~4,400 real amino acids, which is a
reasonable amount of data for a phylogeny. The alarming "1,564 taxa >=70% padded" was an
artifact of a gappy alignment, not 1,564 bad genomes.

Our own vOTUs, by missing data rather than by gaps:

| vOTU | X% (missing) | residue % | real amino acids |
|---|---|---|---|
| vOTU01, 02, 03, 06, 07, 08, 11, 12, 16, 17, 18 | **0.0** | 12.2-17.2 | 4,590-6,482 |
| vOTU04, vOTU09, vOTU14 | 6.7 | 10.9-14.1 | 4,098-5,332 |
| vOTU10 | 8.5 | 9.8 | 3,713 |
| vOTU15 | 33.2 | 9.6 | 3,641 |
| **vOTU13** | **74.8** | **2.7** | **1,022** |

Eleven of our 17 vOTUs in the tree have **zero** missing markers. Only **vOTU13**
(`SRR28004137_bin_11`) is genuinely data-poor at 74.8% missing and ~1,022 real amino acids;
**its position in the tree should not be trusted**, and it is flagged rather than removed
because dropping a taxon is a scientific decision, not a technical one.

**Column occupancy, and the decision to keep trimAl `-automated1`:**

```
columns with data in >=  5% of taxa:  10638 ( 28.2% of 37739)
columns with data in >= 10% of taxa:   7867 ( 20.8% of 37739)
columns with data in >= 20% of taxa:   6074 ( 16.1% of 37739)
columns with data in >= 50% of taxa:   4165 ( 11.0% of 37739)
```

Stricter gap trimming would cut the alignment to a quarter or less of its length, which
would run far faster. **It was not done**, for two reasons: the paper specifies
`-automated1`, and this run exists to reproduce the paper; and 21,156 of the columns are
parsimony-informative, so the length is carrying signal rather than only gaps. This is
recorded as a deliberate choice with a measured basis, and a stricter trim
(`trimal -gt 0.1`) is the obvious first thing to try if the tree turns out poorly resolved.

**Changes made to `scripts/11_gvmag_tree.slurm`** as a result (the previous version was kept
as `scripts/11_gvmag_tree.slurm.orig`):

- `--mset LIST` - restrict ModelFinder to named matrices.
- `--cmax N` - cap the rate categories ModelFinder tries.
- `--safe` - pass IQ-TREE's `-safe`.
- A note, and a log line, documenting that the script passes **no** `--redo`, so IQ-TREE
  finds `<prefix>.ckp.gz` and **resumes from its checkpoint** if a previous run hit the
  walltime. This is what makes the job safe to run on a partition whose time limit is
  shorter than the whole analysis: resubmit the identical command and it continues.

#### 11i. Running step 11 (submitted 2026-11-16)

```bash
( module purge; module load slurm/kuro; \
  sbatch -p serial -N 1 -c 64 --mem=96G -t 48:00:00 \
         scripts/11_gvmag_tree.slurm --mset LG,WAG,JTT --cmax 4 --safe )
# Submitted batch job 111564
```

Which becomes:

```
iqtree2 -s .../kbrevis_gvmags.trimmed.aln -m MFP --bnni -T 64 \
        --prefix .../tree/kbrevis_gvmags --seed 12345 --mem 81G \
        -B 1000 --wbtl -mset LG,WAG,JTT -cmax 4 -safe
```

**Every parameter, and why:**

| choice | value | why |
|---|---|---|
| cluster | kuro `serial` | 64 cores and an immediate start. astral `batch` (the only other fitting non-preemptible option) would not have started for ~24 h, and gulf `gz` is preemptible with `PreemptMode=CANCEL` - unacceptable for a job of this length even though it ranked first on estimated finish |
| cores | 64 | ~512 distinct patterns per thread, above IQ-TREE's rule of thumb of a few hundred per core |
| memory | 96 GB | measured requirement is 30 GB at 4 threads; the margin covers 64 threads' buffers. The script's original 128 GB was unjustified |
| walltime | 48 h | kuro's maximum. The analysis may well need longer, which is fine: IQ-TREE checkpoints, so resubmitting the identical command resumes |
| `-m MFP` | kept | the paper's ModelFinder-by-BIC approach is preserved |
| `-mset LG,WAG,JTT` | **deviation** | bounds ModelFinder to the three standard protein matrices instead of ~80 model combinations. Justified by the measured 1,492 s per model (section 11h): unrestricted selection would consume the entire walltime before any tree was built. LG is the modern general protein matrix; WAG and JTT are the classic alternatives, so the best-fit model is still chosen by BIC from a defensible set |
| `-cmax 4` | **deviation** | caps rate categories at 4 (IQ-TREE's default is 10), for the same cost reason. +G4/+R4 is the standard choice in protein phylogenetics |
| `-safe` | **addition** | requested by IQ-TREE itself on this alignment |
| `-B 1000 --wbtl` | kept | 1,000 ultrafast bootstrap replicates, as the paper did; `--wbtl` writes the bootstrap trees with branch lengths for the RAxML IC/TC step |
| `--bnni` | kept | NNI-optimises each UFBoot replicate, which the IQ-TREE authors recommend because UFBoot otherwise overestimates support under model misspecification - a real risk with concatenated markers |
| `--seed 12345` | kept | reproducibility |

The job started immediately on ku59 and IQ-TREE confirmed the safe kernel and thread count:

```
Host:    ku59.sciclone.wm.edu (AVX512, FMA3, 376 GB RAM)
Kernel:  Safe AVX+FMA - 64 threads (64 CPU cores detected)
Alignment has 1564 sequences with 37739 columns, 32764 distinct patterns
```

**Note on step 11's own pre-flight warning.** The script prints
`WARNING: 1564 taxa are >=70% X/gap and are placed unreliably` and lists ten of them. As
section 11h establishes, that warning is **over-broad**: it combines `X` with alignment gaps,
so it fires for every taxon in a gappy alignment. The genuinely data-poor taxa are the 44
that are >=70% `X`, one of which is our vOTU13. The warning's underlying `pad()` function
should be split into separate X and gap figures; this is noted rather than changed, because
the script was already running when the problem was understood and editing a file does not
affect a submitted job - Slurm snapshots the script at submit time.

----------------------------------------------------------------------------------------

#### 11j. Step 11 results: the tree (2026-11-18)

Written by Claude Code on 2026-11-18 (prompt 47 in `PROMPT_LOG.md`).

Job **111564** **COMPLETED** in **1 d 11 h 46 m** (35.8 h), exit 0, peak memory **42.5 GB**
of the 96 GB requested. The 96 GB request was right: the section 11h probe measured 30 GB
at 4 threads, and 64 threads pushed it to 42.5 GB, so the measured figure plus headroom was
the correct call and the script's original 128 GB would have been wasteful. It finished
inside kuro's 48 h limit, so the checkpoint-resume path was never needed.

```
model chosen:       WAG+F+I+R4 chosen according to BIC
log-likelihood:     -9852356.848 (s.e. 99100.142)
AIC:                19711015.696
taxa in tree:       1564
bootstrap reps:     1000 (ultrafast UFBoot, --bnni)
constant sites:     7985 (21.2%)
parsimony-informative sites: 21156
distinct site patterns:      32764
```

**The model was `WAG+F+I+R4`, not an LG variant.** This is worth recording because
restricting ModelFinder to LG alone - which was considered in section 11h to save time -
would have forced the wrong matrix. Including WAG and JTT in `-mset` is what let BIC find
it. An interim reading of the log part-way through model selection suggested LG was winning;
that was wrong, and the final BIC table decided otherwise.

**Support across the tree:**

```
1561 internal branches; median 100
>=95 (well supported): 998 (64%)
<70  (poorly supported): 225 (14%)
```

RAxML internode/tree certainty (`kbrevis_gvmags.IC_TC_info.txt`):

```
Tree certainty for this tree:          1143.280485
Relative tree certainty:               0.732403
Relative TC including all conflicting bipartitions (TCA): 0.733119
```

A relative tree certainty of **0.73** across 1,564 taxa at 14.5% average pairwise identity
is a reasonable result: most bipartitions are consistent across the bootstrap trees, but a
real minority conflict. Remember UFBoot values are **not** comparable to classical
bootstrap - `>=95` UFBoot is roughly as strong as `>=70` classical - so the 64% figure is
the one to quote, not a classical-bootstrap equivalent. The 225 branches below 70 are
concentrated in deep relationships, exactly as the 14.5% identity predicted (section 10j).

#### 11k. Where the *K. brevis* vOTUs sit, and two findings worth following up

`scripts/votu_tree_context.py` (new) answers the question a 1,564-tip Newick file cannot
answer by eye: for each of our vOTU tips it walks up to the parent node, reports the UFBoot
support on the branch that places it, and identifies the sister group's order and family
from `reference_manifest.tsv`. Output: `output/tree/votu_tree_context.tsv`.

| vOTU | UFBoot | sister tips | sister order | sister family |
|---|---|---|---|---|
| vOTU01 | 97 | 5 | Imitervirales ×5 | f_Mesomimiviridae ×5 |
| vOTU02 | 100 | 1 | Imitervirales | f_Mimiviridae |
| vOTU03 | 100 | 2 | *our own vOTUs* | — |
| **vOTU04** | **100** | 1 | **Asfuvirales** | **f_AF_01** |
| vOTU06 | 100 | 1 | Imitervirales | f_Mesomimiviridae |
| vOTU07 | 100 | 1 | *our own vOTU* | — |
| vOTU08 | 100 | 1 | Imitervirales | f_Mimiviridae |
| vOTU09 | 100 | 1 | *our own vOTU* | — |
| vOTU10 | 100 | 1 | Imitervirales | f_Mesomimiviridae |
| vOTU11 | 100 | 1 | Imitervirales | f_IM_14 |
| vOTU12 | 100 | 1 | *our own vOTU* | — |
| vOTU13 | 100 | 1 | *our own vOTU* | — |
| vOTU14 | 100 | 1 | Pimascovirales | f_PM_01 |
| vOTU15 | 100 | 1 | Algavirales | f_AG_04 |
| vOTU16 | 100 | 1 | Imitervirales | f_Mimiviridae |
| vOTU17 | 100 | 16 | Imitervirales ×16 | f_IM_07 ×16 |
| vOTU18 | 100 | 1 | Imitervirales | f_Mimiviridae |

**Every placement is supported at 97 or above**, and all but vOTU01 at 100. The overall
picture matches steps 08 and 09: the *K. brevis* bloom virome is dominated by
*Imitervirales*, with nearest relatives in *Mimiviridae* and *Mesomimiviridae*, plus single
representatives in *Pimascovirales* (vOTU14) and *Algavirales* (vOTU15).

**Finding 1: the tree contradicts GVClass about vOTU04, and GVClass was already unsure.**
vOTU04 (`SRR34748664_bin_2`) is placed sister to an **Asfuvirales** reference (`f_AF_01`)
with **100 UFBoot support**. But step 08's GVClass assigned it
`o_Imitervirales;f_IM_07`. Looking at the underlying vote explains the conflict - the call
was never confident:

```
NCLDV__Imitervirales:4(66.67%),NCLDV__Asfuvirales:2(33.33%)
```

Four marker genes voted *Imitervirales*, two voted *Asfuvirales*, and GVClass reported the
majority. The concatenated seven-marker phylogeny - which uses all the markers jointly
rather than voting - puts it firmly with *Asfuvirales*. It also has **no GVMAGs V2 relative
above the fastANI threshold** (section 09f), so ANI cannot break the tie. **This is exactly
what a marker-gene tree is for, and the tree should be preferred over the vote here.** It is
also the one vOTU whose order assignment changes between steps 08 and 11, so any table
carrying per-vOTU taxonomy needs a footnote. Worth confirming by looking at which markers
voted which way before publishing either label.

**Finding 2: two small clades contain only *K. brevis* viruses, with no reference
interleaved.**

- **{vOTU03, vOTU09, vOTU13}** - vOTU03's sister group is vOTU13 + vOTU09
- **{vOTU07, vOTU12}**

Groups of our own genomes that are each other's closest relatives, with no GVMAGs V2
genome falling between them, mean lineages **not represented in the reference database** -
consistent with these vOTUs having no close sequenced relative (section 09f). If these hold
up they are the most interesting result in the analysis: candidate *K. brevis*-associated
giant virus lineages rather than novel members of known families.

**But treat the vOTU13 membership with suspicion.** vOTU13 has **74.8% missing data** and
only ~1,022 real amino acids (section 11h). Taxa with very little data can attach to each
other with spuriously high support, because there is too little signal to place them
anywhere else - so the 100 UFBoot on vOTU13's branch does not mean what 100 means for a
complete genome. The **{vOTU03, vOTU09}** pairing is credible (0.0% and 6.7% missing data
respectively); vOTU13 joining them is not established by this tree. The {vOTU07, vOTU12}
pair is credible, both having 0.0% missing data.

**Caveats on the whole tree.**

- **It is UNROOTED.** Booker et al. midpoint-rooted theirs. Nothing about the branching
  *order* means anything until it is rooted (iTOL: Advanced -> Midpoint root).
- vOTU05 and vOTU19 are **not in the tree at all** - they carried none of the seven markers
  and `ncldv_markersearch` dropped them (section 10j). vOTU05 remains the open question from
  step 10: 100% CheckV complete, 203 kb, no reference relative, yet no core marker genes.
- The sister group is only as informative as the references included. These placements are
  relative to the 1,523 `mcl`-tier in-order references; the finer `skani` tier (6,266 taxa)
  could insert a reference inside one of the two *K. brevis* clades and dissolve it. **That
  is the single most valuable follow-up**, and it is one flag: re-run step 10 with
  `--tier skani`.
- 225 of 1,561 branches are below 70 UFBoot. None of them is a vOTU placement, but they do
  affect the deep backbone the vOTUs hang from.

#### 11l. Files copied to `output/tree/`

Copied from `~/scr10/data/kbrevis/tree` (the scratch run directory is left intact):

| file | size | what it is |
|---|---|---|
| `kbrevis_gvmags.treefile` | 80 KB | **the ML tree with UFBoot support** - this is the file to upload to iTOL |
| `kbrevis_gvmags.contree` | 80 KB | the extended majority-rule consensus of the bootstrap trees |
| `kbrevis_gvmags.iqtree` | 732 KB | human-readable report: the full BIC model table, likelihood, parameters |
| `kbrevis_gvmags.IC_TC_info.txt` | 228 KB | RAxML internode and tree certainty scores |
| `kbrevis_gvmags.IC_branch_labels.nwk` | — | the tree with IC scores as branch labels, for viewing certainty instead of UFBoot |
| `kbrevis_gvmags.log` | 132 KB | the full IQ-TREE log |
| `kbrevis_gvmags.tree_summary.txt` | 558 B | the summary block quoted above |
| `votu_tree_context.tsv` | — | the placement table above |
| `itol_*.txt` | — | five iTOL annotation files, regenerated against the final tree |

**Deliberately not copied:** `kbrevis_gvmags.ufboot` (74 MB of bootstrap trees),
`kbrevis_gvmags.ckp.gz` (18 MB checkpoint) and `kbrevis_gvmags.mldist` (24 MB distance
matrix). None is needed to visualise or interpret the tree, and committing 116 MB to git
would be the mistake Step 12 exists to avoid. They remain in scratch, which is purged after
90 days without access - so if the bootstrap trees are ever needed again (for example to
recompute IC/TC differently), copy them to `~/data10` before then.

The earlier `kbrevis_gvmags.PRELIMINARY_no_support.treefile` (prompt 44) has been
**deleted**, since a support-free tree sitting beside the real one is a trap.

----------------------------------------------------------------------------------------

# Step 12. Copying results out of scratch into the repository (2026-11-16)

Written by Claude Code on 2026-11-16 (prompt 40 in `PROMPT_LOG.md`). File:
`scripts/copy_results.slurm`, plus a new `.gitignore`.

#### 12a. Why this step exists

Everything steps 08 and 09 produced lives in `~/scr10/data/kbrevis`. Scratch is **never
backed up and is purged after 90 days without access** (AGENTS.md section 6), so results
left there will eventually vanish. This step copies them into the repository at
`~/K_brevis_viromes/output`, where they sit beside the scripts that made them and are
covered by the weeknightly home backup.

The scratch originals are **copied, not moved** - `~/scr10/data/kbrevis` is untouched, so
steps 10 and 11 keep reading their inputs from scratch as they were written to do.

#### 12b. What could not be copied, and why

The request was to copy `gvmags` and `votus` whole. That is not possible, and the reason is
worth recording rather than quietly working around:

| directory | size |
|---|---|
| `votus` | 16 MB |
| `gvmags` | **115 GB** |

The home filesystem gives this account a **50 GB quota**, of which 21.2 GB was free, and the
underlying filesystem was 92% full. Copying 115 GB into home is physically impossible, not
merely inadvisable; attempting it would have filled the quota and died partway through,
while putting avoidable pressure on a nearly-full shared filesystem.

Measuring where the 115 GB actually sits shows the problem is confined to one directory:

```bash
for d in binning bins gvclass_in gvclass checkv genomad gvmags; do
  printf "%-12s %s\n" "$d" "$(du -shc gvmags/*/$d | tail -1 | cut -f1)"
done
```

```
binning      108G      <- excluded
bins         1.9G
gvclass_in   1.9G      <- excluded
gvclass      2.5G
checkv       27M
genomad      11M
gvmags       14M
```

**Two directories are excluded, both for good reasons:**

- **`binning/` (108 GB)** holds, per sample, the bbmap BAM of reads against the merged
  contigs, its index, and the `jgi_summarize_bam_contig_depths` table. These are the
  mapping intermediates MetaBAT2 consumed to decide which contigs co-vary in depth.
  Nothing downstream reads them again - step 09 works from the GVMAG FASTAs, and steps 10
  and 11 from the vOTU representatives. They are also fully regenerable from the step 06
  merged contigs plus the step 02 cleaned reads by re-running `08_gvmags.slurm`. This is
  the only directory that makes the copy impossible, and it contains the least irreplaceable
  data in the tree.
- **`gvclass_in/` (1.9 GB)** is a byte-for-byte duplicate of `bins/`, differing only in file
  naming because GVClass requires a `.fna` extension where MetaBAT2 writes `bin.10.fa`.
  Verified on 2026-11-16: both hold 41 files for SRR28004137, and
  `bins/SRR28004137_bin.10.fa` and `gvclass_in/SRR28004137_bin_10.fna` are both exactly
  244,266 bytes with distinct inodes (separate copies, not hardlinks). Excluding it loses
  no information whatsoever.

**`bins/` (1.9 GB) is deliberately kept**, even though it is an intermediate. It holds all
848 MetaBAT2 bins, which means the *classification* half of step 08 (GVClass, CheckV,
geNomad, GVMAG selection) can be re-run without redoing the 108 GB of read mapping. That is
a cheap insurance policy for 1.9 GB.

If the BAMs are ever worth keeping, the right destination is `~/data10` (251 TB free, weekly
backup), not home. That was **not** done here, because nothing has asked for it and 108 GB
of regenerable intermediates is poor use of archival space.

#### 12c. SCRIPT: `scripts/copy_results.slurm`

The copy runs as a batch job rather than on dev, because AGENTS.md section 6 requires
GB-scale copies between filesystems to go through Slurm. One core is correct: this is pure
I/O and nothing in it is parallel.

Two safeguards are built in. The script **checks the free quota and aborts before copying
anything** if the transfer plus a 2 GB margin would not fit, so a too-large copy fails
immediately and harmlessly instead of filling the filesystem and dying halfway. And it uses
`rsync -a` rather than `cp`, so an interrupted run is resumed by simply resubmitting the
same job - rsync skips what is already identical.

```bash
#!/bin/bash
#SBATCH --job-name=copy_results
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=4:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/copy_results_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/copy_results_%j.err

set -euo pipefail

SRC="$HOME/scr10/data/kbrevis"
DEST="/sciclone/home/iadanforth/K_brevis_viromes/output"

mkdir -p "$DEST"

# Refuse to start if the copy cannot fit, rather than filling a 92%-full shared
# filesystem and dying halfway through.
NEED_MB=$(du -sm --exclude=binning --exclude=gvclass_in "$SRC/gvmags" "$SRC/votus" \
          | awk '{t+=$1} END{print t}')
FREE_MB=$(quota -s 2>/dev/null \
          | awk '/sciclone\/home/{getline} /^ *[0-9]/{gsub("M","",$1); gsub("M","",$2); print $2-$1}')
echo "[copy] need ${NEED_MB} MB, home quota has ${FREE_MB:-unknown} MB free"
if [[ -n "${FREE_MB:-}" ]] && (( NEED_MB + 2048 > FREE_MB )); then
    echo "[copy] ABORT: not enough room in the home quota (need ${NEED_MB} MB + 2 GB margin)." >&2
    exit 1
fi

# -a preserves times/permissions; --exclude drops the two directories named above.
# rsync (not cp) so an interrupted copy can be resumed by simply resubmitting.
for d in gvmags votus; do
    echo "[copy] $SRC/$d -> $DEST/$d"
    rsync -a --info=stats2 \
        --exclude='binning/' --exclude='gvclass_in/' \
        "$SRC/$d/" "$DEST/$d/"
done

echo
echo "[copy] verifying: every file present in the source (minus exclusions) is in the copy"
MISSING=$(rsync -rn --itemize-changes \
            --exclude='binning/' --exclude='gvclass_in/' \
            "$SRC/gvmags/" "$DEST/gvmags/" | grep -c '^>f' || true)
echo "[copy] files still differing after the copy: $MISSING (0 = complete)"

echo
echo "[copy] result sizes:"
du -sh "$DEST"/gvmags "$DEST"/votus
# Count with shell globs, NOT `find -path`: in find, '*' matches '/' as well, so
# -path '*/gvmags/*.fna' also matches .fna files nested deeper (e.g. under gvclass/),
# which inflated this line to 53 on the first run when the real GVMAG count is 26.
shopt -s nullglob
gv=( "$DEST"/gvmags/*/gvmags/*.fna )
bn=( "$DEST"/gvmags/*/bins/*.fa )
echo "[copy] GVMAG .fna files copied: ${#gv[@]}"
echo "[copy] MetaBAT2 bins copied:    ${#bn[@]}"
echo "[copy] done. Scratch originals in $SRC are untouched."
```

#### 12d. How it was run, and the result

```bash
( module purge; module load slurm/kuro; sbatch -p serial scripts/copy_results.slurm )
# Submitted batch job 111460
```

Completed in **3 minutes 30 seconds**, exit 0, empty `.err`:

```
[copy] need 4512 MB, home quota has 21685 MB free
[copy] files still differing after the copy: 0 (0 = complete)
4.4G    /sciclone/home/iadanforth/K_brevis_viromes/output/gvmags
16M     /sciclone/home/iadanforth/K_brevis_viromes/output/votus
[copy] done. Scratch originals are untouched.
```

The verification step is the part that matters: a second `rsync -rn --itemize-changes` pass
re-walks the source and reports any file that would still need transferring. It found **0**,
so the copy is complete rather than merely "finished without an error".

Final contents of `output/`:

| path | what it is |
|---|---|
| `output/votus/` | all of step 09: `gvmags_votus.tsv`, `gvmags_representatives.txt`, `votu_summary.txt`, `votu_nearest_reference.tsv`, `all_gvmag_summaries.tsv`, `ani/`, `heatmap/`, and the 26 pooled GVMAG FASTAs |
| `output/gvmags/<sample>/gvmags/` | the **26 GVMAGs**, the headline output of step 08 |
| `output/gvmags/<sample>/bins/` | all **848 MetaBAT2 bins** |
| `output/gvmags/<sample>/gvclass/` | GVClass taxonomy output per sample |
| `output/gvmags/<sample>/checkv/` | CheckV completeness/contamination per sample |
| `output/gvmags/<sample>/genomad/` | geNomad output per sample |
| `output/gvmags/<sample>/*_gvmag_summary.tsv` | the per-sample summary tables (9 samples have one) |

Home quota went from 28.8 GB to **33.2 GB used of 50 GB** (16.8 GB free).

#### 12e. `.gitignore`

`output/` is deliberately **not tracked by git**. Several GB of FASTA and tabular data
committed into history would be stored forever and would make the repository effectively
unclonable. The files live in the working tree so they are backed up and sit next to the
code, but git ignores them:

```
# Analysis outputs copied out of scratch by scripts/copy_results.slurm (~5 GB).
# These are data, not code: git would store every byte forever and the repo would
# become unclonable. They live in the working tree so they survive the 90-day scratch
# purge and sit beside the scripts that made them, but they are deliberately untracked.
# To track one specific small file anyway: git add -f output/votus/gvmags_votus.tsv
output/

# Slurm logs are kept (AGENTS.md requires them in logs/), so logs/ is NOT ignored here.
```

Individual small result files can still be committed deliberately with `git add -f`, which
is the right way to put something like `gvmags_votus.tsv` under version control. **Five step
09 result tables were added this way on 2026-11-16** (prompt 41 in `PROMPT_LOG.md`), because
they are the actual findings of step 09, are small (125 B to 24 KB), and are what anyone
reading the repository will want without having to re-run the pipeline:

```bash
git add -f output/votus/all_gvmag_summaries.tsv \
            output/votus/gvmags_representatives.txt \
            output/votus/gvmags_votus.tsv \
            output/votus/votu_nearest_reference.tsv \
            output/votus/votu_summary.txt
git commit   # commit 9356303
git push origin main
```

| file | size | what it holds |
|---|---|---|
| `gvmags_votus.tsv` | 3.4 KB | one row per GVMAG: its vOTU, `sample_set`, whether it represents it, ANI/AF to the representative |
| `votu_nearest_reference.tsv` | 3.0 KB | each vOTU representative's closest GVMAGs V2 relative, with that reference's taxonomy and ecosystem |
| `all_gvmag_summaries.tsv` | 24 KB | every sample's per-GVMAG summary combined, with `<sample>__<bin>` names and a `sample_set` column |
| `gvmags_representatives.txt` | 1.8 KB | the 19 representative `.fna` paths, the input to step 10 |
| `votu_summary.txt` | 125 B | the one-line vOTU count |

Two artifacts were explicitly kept **out** of the repository at the same time and added to
`.gitignore`: `scripts/__pycache__/` (bytecode written whenever `select_refs.py` is
imported, e.g. by the step 10 dry run) and `*.bak` (git history is the real backup, so a
committed `select_refs.py.bak` beside the fixed file would only be noise - the file is still
on disk, untracked).

#### 12f. Caveats

- This is a **copy**, so `output/` and scratch will drift apart as soon as anything re-runs.
  `scripts/copy_results.slurm` is idempotent and safe to resubmit, which is the intended way
  to refresh it; note that `rsync` without `--delete` will not remove files from `output/`
  that have disappeared from scratch.
- The 108 GB of `binning/` BAMs are **not** preserved anywhere. If step 08's binning ever
  needs to be audited rather than re-run, that data is gone once scratch is purged.
- Home is backed up weeknightly with roughly one month of retention, and there is no
  archival storage on this system. For anything that must outlive the project, copy it off
  the cluster.
- `output/` counts against the 50 GB home quota. With 16.8 GB free, a second round of
  results this size would need either `~/data10` or a cleanup first.

----------------------------------------------------------------------------------------

----------------------------------------------------------------------------------------

# Step 13. Cleaning up scratch (2026-11-17)

Written by Claude Code on 2026-11-17 (prompt 45 in `PROMPT_LOG.md`).

`~/scr10/data/kbrevis` had accumulated eleven working directories plus four smoke-test
leftovers. **Scratch has no backup and no undo**, and step 11 (job 111564) was running at
the time against `phylogeny/` and `tree/`, so nothing was deleted until each candidate had
been checked against the scripts that use it.

**How the candidates were identified.** Rather than judging by name, every directory was
grepped for across all scripts:

```bash
cd scripts
for d in raw trim clean host_filtered reference contigs contigs_no_v \
         contigs_merged gvmags votus phylogeny tree; do
  printf "%-16s %s\n" "$d" "$(grep -lE "kbrevis/$d\b|/$d/" *.slurm *.sh *.py 2>/dev/null | tr '\n' ' ')"
done
```

Exactly one directory came back with **no script referencing it**: `host_filtered/`. Its
contents confirmed why - a single merged read pair
(`kbrevis_reads_{1,2}_hostfree.fastq`, 7.2 GB each) rather than the per-sample files the
pipeline uses, dated 25 October, from before the 20-sample workflow existed. The README's
own directory listing already described it as "leftover from an earlier attempt".

**Deleted (5 directories, 14 GB):**

| directory | size | what it was |
|---|---|---|
| `host_filtered/` | 14 GB | superseded by `clean/`; no script reads it; an early single-pair host-filtering attempt |
| `nms_smoketest/` | 1.9 MB | smoke test of `ncldv_markersearch` while writing step 10 |
| `_tailtest_smoketest/` | 860 KB | the step 08 CheckV/geNomad smoke test (the bacterial bin that had to be moved out of the GVMAG root, section 09a-bis) |
| `tree_probe/` | 612 KB | the IQ-TREE memory probe from section 11h |
| `tree_smoketest/` | 28 KB | 6-taxon alignment used to test step 11 |

```bash
cd ~/scr10/data/kbrevis
rm -rf tree_smoketest tree_probe _tailtest_smoketest nms_smoketest host_filtered
```

The four smoke-test directories were all created by this agent while developing steps 08,
10 and 11; none is an input to anything. Step 11 continued running throughout, unaffected.

**Kept, and why.** Everything else is referenced by at least one script:

| directory | size | why it stays |
|---|---|---|
| `raw/` | 2.0 TB | step 01's downloads; input to step 02. Re-downloadable from SRA, but slowly |
| `trim/` | 1.6 TB | step 02 output; input to step 04 |
| `clean/` | 1.1 TB | step 04 output. **Still actively needed** - step 08 maps these reads to contigs, so re-running binning requires them |
| `contigs_no_v/` | 322 GB | metaSPAdes assemblies; input to step 06 |
| `gvmags/` | 115 GB | step 08 output (results already copied to `output/`, Step 12) |
| `contigs/` | 90 GB | metaviralSPAdes assemblies; input to step 06 |
| `reference/` | 79 GB | contaminant genomes and the bowtie2 host index; needed to re-run step 04 |
| `contigs_merged/` | 6.6 GB | step 06 output; step 08's input |
| `phylogeny/` | 629 MB | step 10 output; **step 11 is reading it** |
| `tree/` | 42 MB | **step 11 is writing it** |
| `votus/` | 16 MB | step 09 results |

**There is no space pressure, so nothing further was deleted.** `/sciclone/scr10` is 219 TB,
72% used, with **62 TB free**. The obvious further candidates, if space ever matters, are:

- `gvmags/*/binning/` - **108 GB** of bbmap BAMs, indexes and depth tables that nothing
  downstream reads and that `08_gvmags.slurm` regenerates (already excluded from the
  `output/` copy for exactly this reason, Step 12b).
- `raw/` + `trim/` - **3.6 TB** of upstream reads, regenerable from SRA and fastp, needed
  only if the cleaning step is re-run from scratch.

Deleting either is a judgement call about re-computation cost versus disk, not a tidiness
question, so neither was done unprompted.

#### 13a. Log files: almost none of them are irrelevant

`logs/` held 234 files (18 MB), every one of them tracked in git. The question "which are
irrelevant?" was answered by checking, not by guessing from filenames, because these logs
are the evidence behind the README's "How this was tested" sections - deleting a log that a
section cites would leave a claim with nothing behind it.

Every job ID appearing in a log filename was grepped for in `README.md` and
`PROMPT_LOG.md`:

```bash
for id in 73841 73842 73843 73846 73847 73848 73849 110122 611846 621411 111548; do
  grep -c "\b$id\b" README.md PROMPT_LOG.md
done
```

The result was the opposite of what the directory listing suggests. Logs that *look*
disposable are mostly cited:

| log | cited? | what cites it |
|---|---|---|
| `nms_smoketest_73841` | yes | step 10's testing section |
| `beast_introspect_73842` | yes | the BEAST X introspection note |
| `tree_smoketest_73846` | yes | "73846 (3 taxa) failed usefully, revealing..." |
| `tree_smoketest2_73847` | yes | the `--bootstraps 0` test |
| `tree_smoketest3_73848` | yes | the alignment-too-short failure |
| `tree_smoketest4_73849` | yes | the 6-taxon test built from real data |
| `08_gvmags_tailtest_110122` | yes | 7 mentions; the step 08 tail smoke test |
| `07_binnning_test_611846` | yes | 4 mentions (the misspelt name is deliberate - it records a real job) |
| `test_checkV_db__621411` | yes | the CheckV database failure, quoted verbatim |
| `08_gvmags_pilot_110103` | yes | the `JAVA_HOME` crash |

**Only two jobs were referenced nowhere at all**, and only those were deleted (4 files):

| deleted | why |
|---|---|
| `tree_introspect_73843.{out,err}` | exploratory introspection of the tree tooling; superseded by the `tree_smoketest*` jobs, which are the ones the README actually cites |
| `11_pre_tree_check_111548.{out,err}` | my own failed first attempt at the step 11 probe. It died in 1 s with `No such file or directory` because the script had been written to the session scratchpad under `/tmp`, which is node-local and invisible to compute nodes. Resubmitted as job 111549 within minutes, and the lesson is recorded in prose in section 11h, so the log adds nothing |

```bash
cd logs
rm -f tree_introspect_73843.out tree_introspect_73843.err \
      11_pre_tree_check_111548.out 11_pre_tree_check_111548.err
```

**Nothing else was deleted, on purpose.** Checks that might have justified a larger cull all
came back negative:

- **No log is empty on both sides.** Not one job produced an empty `.out` *and* an empty
  `.err`. (50 individual files are empty, but an empty `.err` beside a full `.out` is
  itself information - it says the job ran without errors.)
- **The large logs are real output**, not noise: `05_metaspades_109630.out` (1.4 MB) is an
  assembly log, and `09_gvmag_votus_111457.err` (868 KB) is fastANI's verbose INFO stream
  from a successful run.
- **Size is not a reason.** The whole directory is 18 MB.

Since all 230 remaining files are tracked in git, any future deletion is recoverable from
history - which is the main reason deleting logs is far less risky than deleting scratch.

**A side effect worth knowing.** Scratch is purged after **90 days without access**, and the
`du` and `stat` commands used to inventory these directories updated their access times. The
purge clock for everything under `~/scr10/data/kbrevis` therefore now runs from 2026-11-17.
That is harmless - it delays purging rather than hastening it - but it means access times
in this tree currently reflect inspection, not pipeline use.

----------------------------------------------------------------------------------------

----------------------------------------------------------------------------------------

# Step 14. vOTU abundance in the environment through time (2026-11-18)

Written by Claude Code on 2026-11-18 (prompt 48 in `PROMPT_LOG.md`). Files:
`envs/coverm-env.yml`, `envs/r-plot-env.yml`, `scripts/14_abundance_setup.slurm`,
`scripts/14_votu_reference.sh`, `scripts/14_abundance_smoketest.slurm`,
`scripts/14_votu_abundance.slurm`, `scripts/merge_votu_abundance.py`,
`scripts/plot_votu_abundance.R`.

#### 14a. The question, and what the sample metadata actually says

The question: **the giant viruses recovered from the incubation experiment - how abundant
were they in the environment, month by month?** Steps 08-11 say what these viruses are and
how they are related; nothing so far says whether they are out there in the water.

Two things in `metadata/SraRunTable.csv` shaped the design and both need stating, because
neither is obvious from the sample list.

**First, the vOTUs do not all come from the incubation.** Of the 19 vOTUs, **12 were
discovered in incubation samples and 7 in environmental samples**:

| discovered in | vOTUs | source samples |
|---|---|---|
| incubation | 02, 03, 06, 07, 08, 09, 11, 12, 13, 17, 18, 19 | SRR28004137, 138, 139, 144, 154 |
| environmental | 01, 10, 16 · 05, 14 · 15 · 04 | SRR28004148 · SRR28004151 · SRR28004152 · SRR34748664 |

All 19 were mapped, not just the 12. The 7 environmental-derived vOTUs are a **positive
control that costs nothing**: each one must be detected in the sample it was assembled
from, and if it is not, the mapping is broken. That check is worth more than the compute it
saves to skip it.

**Second, the "environmental samples" are two different surveys.**

| samples | n | dates | source |
|---|---|---|---|
| SRR28004145-153 | 9 | 2020-12 to 2021-08, **one per month** | Gulf of Mexico |
| SRR34748663, SRR34748664 | 2 | 2021-09-29, 2022-03-08 | labelled **coastal Florida waters** in the SRA, with recorded depths and day-resolution dates. That label is a mislabelling - see section 18a; they are the same survey as the other nine |

The nine Gulf samples are a genuine monthly time series. The two Florida samples are not
part of it - different location, different study, different accession series. Drawing one
line through all eleven would imply a continuous 15-month survey that does not exist, so
the figures draw the connecting line through the Gulf series only and mark the Florida
samples as separate points beyond a dotted divider. For orientation, the incubation
experiment's two batches (CH1 2021-05, CH2 2021-08) sit inside the Gulf series, which is
what makes "were these viruses in the water before and after?" answerable at all.

#### 14b. Method, and the three parameter choices that matter

**CoverM 0.7.0** maps the reads and summarises coverage **per genome** rather than per
contig. That is the reason for using it rather than the bbmap/`jgi_summarize_bam_contig_depths`
combination already installed for step 08: a GVMAG is many contigs, and a per-contig table
would turn 19 vOTUs into hundreds of uninterpretable rows. CoverM bundles minimap2, so no
separate aligner was needed. Neither CoverM nor R was installed on this system; both
environments were created by `14_abundance_setup.slurm`.

Each vOTU is represented by **one genome - its step 09 representative - not all its
members**. Including every member would split reads from a shared virus across
near-identical genomes and make each look less abundant than it is.

| choice | value | why |
|---|---|---|
| read identity | **>=95%** | the MIUViG species boundary, and the same cut-off that defined the vOTUs in step 09. Without it, reads from a related-but-different virus inflate a vOTU's abundance and closely related vOTUs steal each other's reads |
| read aligned | **>=75%** | rejects reads that match only through a short conserved patch |
| presence | **>=75% genome breadth** | MIUViG (Roux et al. 2019). A real detection covers most of the genome; a few conserved genes attracting stray reads covers very little. **Depth alone cannot tell these apart**, which is the whole point - see 14d |
| CoverM's own breadth filter | **disabled** (`--min-covered-fraction 0`) | its default of 10 silently rewrites low-breadth abundances to zero, baking a threshold into the raw numbers. Keeping the raw values and applying the 75% rule in `merge_votu_abundance.py` makes the threshold explicit and changeable **without re-mapping ~200 GB of reads** |

#### 14c. Three bugs the smoke test caught (`14_abundance_smoketest.slurm`)

Before committing an 11-task array to ~200 GB of reads, the exact CoverM command was run on
200,000 read pairs. It cost 3 seconds and caught three faults, each of which would have been
silent or late:

1. **`--min-read-percent-identity` takes a PERCENTAGE, not a fraction.** The command had
   `0.95`. CoverM accepts that without complaint and reads it as **0.95% identity** - no
   filtering at all. Every abundance would have been inflated by cross-mapping between
   related vOTUs, with nothing in any log to show for it. This is the one that would have
   quietly corrupted the result. `coverm genome --full-help` is explicit: "e.g. 95 for 95%".
2. **There is no `covered_fraction` method in CoverM 0.7.0.** Its methods are
   `relative_abundance, mean, trimmed_mean, coverage_histogram, covered_bases, variance,
   length, count, reads_per_base, rpkm, tpm`. Breadth is now computed in the merge step as
   `covered_bases / length`, which is the same number with the arithmetic visible.
3. **CoverM writes LONG format** - one row per genome, with columns named
   `<sample> <Method>`. `merge_votu_abundance.py` had assumed the opposite (one row,
   `<genome> <Method>` columns) and would have matched nothing, producing an empty table.

The smoke test also failed once on a fault of its own - `grep -c` exits 1 when it counts
zero matches, which `set -e` treats as fatal, so the job reported FAILED even though CoverM
had succeeded. Fixed with `|| true`; the corrected run (job 73931) exits 0.

It also gave the first positive control: 200,000 reads from SRR28004148 returned vOTU01,
vOTU10 and vOTU16 - exactly the three vOTUs assembled from that sample - as the top three
by breadth and RPKM.

#### 14d. Results

Mapping: job **651382** (bora `hima`, array of 11, 8 cores, 24 GB), **3-6 minutes per
sample**, all COMPLETED, no errors. Resources came from the smoke test rather than a guess;
the script's original 16 cores / 32 GB / 8 h was cut to 8 / 24 GB / 2 h.

**Seven detections out of 209 vOTU-sample pairs. All seven are environmental-derived vOTUs,
each detected in exactly the one sample it was assembled from.**

| vOTU | sample | date | mean coverage | breadth | RPKM |
|---|---|---|---|---|---|
| vOTU05 | SRR28004151 | 2021-02 | 97.3x | 1.000 | 4072.9 |
| vOTU16 | SRR28004148 | 2021-05 | 27.9x | 1.000 | 650.4 |
| vOTU10 | SRR28004148 | 2021-05 | 17.9x | 0.999 | 414.9 |
| vOTU15 | SRR28004152 | 2021-01 | 17.2x | 1.000 | 3104.5 |
| vOTU14 | SRR28004151 | 2021-02 | 15.5x | 0.999 | 644.1 |
| vOTU04 | SRR34748664 | 2021-09 | 11.9x | 1.000 | 2553.6 |
| vOTU01 | SRR28004148 | 2021-05 | 9.7x | 0.999 | 222.6 |

**None of the 12 incubation-derived vOTUs is detected in any environmental sample.** That
is the answer to the question as asked, and it is a real result rather than a failed
analysis - the positive controls prove the mapping works.

**The threshold is not doing arbitrary work.** Breadth in these data is sharply bimodal:

| mean coverage | n pairs | median breadth | max breadth |
|---|---|---|---|
| 0 (no reads) | 12 | 0.000 | 0.000 |
| <0.1x | 157 | 0.001 | 0.051 |
| 0.1-0.5x | 28 | 0.045 | 0.152 |
| 0.5-1.4x | 1 | 0.088 | 0.088 |
| **>=1.4x** | 11 | **0.999** | 1.000 |

Only **1 of 209** pairs falls between 15% and 99% breadth, so any threshold in that range
gives the same seven detections. (The 1.4x boundary is not arbitrary either: breadth
saturates as `1 - e^-coverage`, so >=75% breadth is unreachable below ~1.4x coverage. Every
pair that clears 1.4x here has ~100% breadth, so nothing is being lost to depth limitation.)

**Why breadth, not RPKM: the vOTU13 artifact.** Four pairs reach >=1.4x mean coverage and
still fail. The clearest is **vOTU13 in SRR34748664: 32x mean coverage over 0.9% of its
genome** - roughly 2,800 bp of 309,356 bp. Reads are piling onto one small region, almost
certainly a conserved or repetitive element, not a virus that is present. Plotted as RPKM
alone, vOTU13 would have looked like one of the more abundant environmental viruses in the
dataset. This is the same genome that carried only 2 of 7 marker genes and whose tree
placement section 11k flags as untrustworthy; three independent lines now say it is a poor
genome.

#### 14e. What this means, and what it does not

The giant viruses that came out of the incubation experiment were **not detectable in the
monthly environmental survey**. Reasonable readings, which this analysis cannot separate:

- they are genuinely rare *in situ* and the incubation enriched them from below the
  detection limit - the most likely explanation, and the usual reason for incubating;
- the environmental sequencing is not deep enough to cover them at >=1.4x;
- they are incubation-specific, for example induced from a lysogenic state or favoured by
  bottle conditions.

**A limitation that constrains all of it:** each environmental vOTU is detected *only* in
its own source sample, never in an adjacent month. The nine monthly samples therefore share
essentially no giant viruses with one another. The most economical explanation is that this
analysis is **assembly-detection-limited** - a virus is recovered from a sample only when it
is abundant enough there to assemble, and at that abundance it is detectable, while
everywhere else it falls below ~1.4x. So these data support "vOTU X was abundant in month
Y" but are weak evidence for absence elsewhere. **Deeper sequencing, or pooling the monthly
samples, is the way to test persistence.** The figures should not be read as showing viral
turnover between months.

#### 14h. Why we do not see the incubation vOTUs in the environment, and Booker et al. do (2026-11-21)

> **Superseded in part by Step 17 (2026-11-23).** This section tested CoverM's defaults on
> 9 samples against the 19 vOTU representatives. Step 14 was then rebuilt to map against
> each of the 26 individual GVMAGs, as the paper's Methods describe, running both filter
> settings on all 11 environmental samples. The conclusion below is unchanged and the
> evidence is stronger; read section 17d for the current numbers.

Written by Claude Code on 2026-11-21 (prompt 56 in `PROMPT_LOG.md`). File:
`scripts/14b_abundance_relaxed.slurm`.

Section 14d reports that **none of the 12 incubation-derived vOTUs is detectable in any
environmental sample.** Booker et al. 2026 report the opposite for four of them: their
vOTU1, vOTU2, vOTU5 and vOTU6 are described as present at the coastal EH25 site and
correlated with *K. brevis* cell counts. By the step 15 correspondence those are **our
vOTU02, vOTU08, vOTU12 and vOTU03**. The disagreement is real and had to be run down.

**The first hypothesis was wrong.** The obvious suspect was the two incubation samples
whose metaSPAdes assemblies were cancelled (SRR28004142, SRR28004143; section 05b-fix-c).
It cannot be the cause:

- **Neither sample contributed a single GVMAG.** All 26 came from other samples, so none of
  the 12 incubation vOTUs contains any sequence from those two. Better assemblies there can
  *add* vOTUs; they cannot change whether the existing 12 are detected.
- **It is not a threshold effect.** Lowering the breadth rule from 75% to 10% yields **zero**
  extra detections.
- **The scale is wrong by orders of magnitude.** A real detection in our data recruits
  **34,000-173,000 reads**; the incubation vOTUs across the nine Gulf samples recruit a
  **median of 1 read**, maximum 406.

**What the paper actually specifies.** Their Methods give one sentence: *"The relative
abundance and coverage of environmental metagenomic reads recruited to each GVMAG was
calculated using coverm."* No identity cutoff, no breadth cutoff, no version, no settings.
Presence is judged on relative abundance, and they count values they themselves call
background (*"~0.1-0.8 ... and <0.1"*) as presence. **CoverM applies no read-identity
filter by default.** Ours uses `--min-read-percent-identity 95`, the MIUViG species
boundary and the same cut-off that defined the vOTUs in step 09.

**The test.** `14b_abundance_relaxed.slurm` repeats step 14 on the same nine EH25 monthly
samples with CoverM's defaults, writing to `abundance/per_sample_relaxed/` so the main
result is untouched.

Removing the identity filter recruits **100-500x more reads**:

| our vOTU | their vOTU | they report it? | reads at ≥95% | reads, no filter | best breadth |
|---|---|---|---|---|---|
| vOTU02 | vOTU1 | yes | 32 | 5,718 | 0.008 |
| vOTU08 | vOTU2 | yes | 615 | 8,643 | 0.021 |
| vOTU12 | vOTU5 | yes | 8 | 4,651 | 0.005 |
| vOTU03 | vOTU6 | yes | 467 | 5,957 | 0.039 |
| vOTU17 | vOTU3 | no | 121 | 1,848 | 0.016 |
| vOTU07 | vOTU4 | no | 32 | 6,618 | 0.006 |

**So a relative-abundance number computed this way would look like a detection** - thousands
of reads, a non-zero percentage, something to plot against cell counts.

**But look at where those reads land.** Summed over the nine samples:

| vOTU | genome | bases covered | breadth | reads | depth on the covered part |
|---|---|---|---|---|---|
| vOTU02 | 1,226,811 | 9,387 | 0.008 | 5,718 | **91x** |
| vOTU08 | 857,067 | 17,894 | 0.021 | 8,643 | **72x** |
| vOTU12 | 343,608 | 1,562 | 0.005 | 4,651 | **447x** |
| vOTU03 | 530,439 | 20,933 | 0.039 | 5,957 | **43x** |
| vOTU07 | 495,517 | 3,169 | 0.006 | 6,618 | **313x** |

Thousands of reads stacked 43-447 deep on **0.5-4% of each genome**. For comparison, the
seven genomes we know are present - the positive controls, assembled from the very sample
they were mapped back to - have breadth **0.999-1.000**. A virus that is actually there is
covered end to end.

**This is reads from a different, related virus latching onto a conserved gene.** It is the
same artifact as vOTU13 in section 14d, which reaches 32x mean coverage over 0.9% of its
genome. Depth alone cannot distinguish "this genome is present" from "something shares one
conserved gene with it" - **only breadth can**, which is exactly why the MIUViG standard
gates on breadth and why step 14 does too.

**Conclusion, stated carefully.** The most likely explanation for the disagreement is that
the published analysis computed relative abundance from read recruitment **without a
genome-breadth requirement**, so conserved-gene cross-recruitment was counted as presence.
Our stricter criterion says these four vOTUs are **not** present in the environmental
samples at any meaningful abundance.

Two honest caveats against over-claiming:

1. **We cannot show the paper is wrong**, only that its stated method would produce an
   apparent detection from data in which the genome is demonstrably absent. Their Methods
   omit the settings, so they may have applied filters not described; and their correlation
   with *K. brevis* cell counts is a separate line of evidence this test does not address.
2. **The relaxed run does not reproduce their pattern either.** vOTU07 - which they do *not*
   report - recruits more reads (6,618) than three of the four they do. So read recruitment
   alone does not pick out their four, and something else differs as well. Their
   environmental accessions are in Fei et al. and Supplemental Data 1, which we have not
   checked; confirming the sample sets are identical is the obvious next step.

**For a write-up**, the defensible claim is the methodological one: applying the MIUViG
breadth criterion overturns an apparent detection, and the figures above show precisely
why. That is a stronger and more interesting result than simply reproducing the paper.

#### 14i. Restricting to the environmental vOTUs, and the cell-count analysis (2026-11-22)

Written by Claude Code on 2026-11-22 (prompt 58 in `PROMPT_LOG.md`). Files:
`scripts/plot_votu_vs_cells.R`, and `--only-source-set` in `merge_votu_abundance.py`.

`output/abundance/` now covers only the **7 vOTUs found in the environmental samples**.
Both readings of that phrase pick out the same seven - the vOTUs assembled *from* an
environmental sample (vOTU01, 04, 05, 10, 14, 15, 16) are exactly the seven ever *detected*
in one - so there was no ambiguity to resolve. `merge_votu_abundance.py --only-source-set
environmental` applies the filter, giving 77 rows instead of 209.

| file | what it is |
|---|---|
| `votu_abundance_environmental.csv` | the data: one row per vOTU per environmental sample |
| `votu_abundance_timeseries.{pdf,png}` | abundance against time, rebuilt from that CSV alone |
| `karenia_cell_counts.csv` | **a template with the cell-count column left blank** - see below |

**The *K. brevis* cell counts are not in the paper.** Booker et al. plot them in Fig. 4A,
but the caption gives no numbers, no units and no axis ranges; the Results contain no EH25
concentrations; and nothing states where the counts came from. The only numeric counts in
the paper are from lab cultures (17,063-33,141 cells/ml) and one New Pass bloom sample
(1,330,000 cells/l) - **neither is station EH25**. Supplemental Data 14 is the only
supplementary file cited beside the EH25 analysis and is the most likely home for them,
but it is inside a ~26 MB zip that could not be retrieved.

**No values were estimated from the figure.** Reading points off a published plot is a
legitimate technique when done with the image in hand and a digitiser; inventing numbers
that then get a slope, an R² and a p-value is not, and those numbers would be
indistinguishable from real ones once written into a CSV. `karenia_cell_counts.csv` is
therefore shipped with the 11 samples and dates pre-filled and `karenia_cells_per_L`
**empty**, plus a `source` column to record provenance. Fill it from Supplemental Data 14
or the FWC/FWRI HAB Monitoring Database, which publishes counts for Florida coastal
stations and is properly citable.

`scripts/plot_votu_vs_cells.R` is written and tested, and produces
`votu_vs_karenia_cells.{pdf,png}` plus `votu_vs_karenia_lm.csv` (slope, intercept, R²,
p-value and n per vOTU) as soon as the counts are present:

```bash
module load miniforge3; eval "$(conda shell.bash hook)"; conda activate r-plot-env
Rscript scripts/plot_votu_vs_cells.R \
    output/abundance/votu_abundance_environmental.csv \
    output/abundance/karenia_cell_counts.csv \
    output/abundance
```

With the template still blank it stops with an explanation rather than failing obscurely
or fitting nothing. Both paths were tested: the empty case prints the instruction above,
and the working case was exercised with synthetic counts **in a scratch directory**, which
confirmed the figure and the per-vOTU model table are produced and caught a real bug - both
CSVs carry a `survey` column, so the join silently produced `survey.x`/`survey.y` and the
plot failed with `object 'survey' not found`. The synthetic data was deleted and never
entered the repository.

Caveats for when the counts arrive: with nine monthly points a per-vOTU regression has
very little power, and each vOTU is non-zero in exactly one month, so a fit is largely
driven by that single point. A correlation found this way is suggestive at best, and the
two coastal Florida samples are a different study and should probably be excluded from the
fit rather than pooled.

> **Corrected 2026-11-24 (section 18a).** The two samples are *not* a different study - the
> SRA's "coastal Florida" label on them is a mislabelling. They belong in the fit. The
> power caveat above still stands.

#### 14f. Outputs

In `output/abundance/`:

| file | what it is |
|---|---|
| `votu_abundance_long.tsv` | tidy, one row per vOTU per sample: dates, breadth, mean coverage, RPKM, TPM, relative abundance, read count, and the `present` flag |
| `votu_abundance_matrix.tsv` | vOTU x sample matrix of breadth-gated RPKM |
| `votu_abundance_timeseries.{pdf,png}` | **the requested figure** - abundance against time, one panel per vOTU, filled = detected, hollow = recruited but absent, with the 12 undetected panels labelled "never detected" |
| `votu_abundance_heatmap.{pdf,png}` | genome breadth for all 19 vOTUs x 11 samples; the clearest single view of the result |

Figure design notes worth keeping: the time series plots **raw** RPKM rather than the
breadth-gated value, so that read recruitment which fails the breadth test is still visible
as a hollow point - plotting the gated value would show twelve flat zero lines, true but
mute. The y-axes are free per panel, which makes a tall hollow curve (vOTU13) look like
abundance, so every never-detected panel carries a "never detected" label. An early version
of the script dropped never-detected vOTUs to avoid empty panels, which was exactly backwards:
those twelve panels **are** the answer. An early version also rendered with a transparent
background (black in most viewers) and invisible titles, fixed with `bg = "white"` and an
explicit `plot.background`.

#### 14g. Caveats

- **Absence of evidence**: see 14e. Non-detection here means "below ~1.4x coverage in this
  sample", not "absent from the Gulf of Mexico".
- ~~The two coastal Florida samples are a different study and location. They are shown, but
  they are not month 10 and month 16 of the Gulf series.~~ **Corrected 2026-11-24
  (section 18a): the SRA label is wrong and these two are part of the same survey.**
- Mapping used vOTU **representatives**. A vOTU member genome that differs from its
  representative by more than 5% in some region will lose those reads to the 95% identity
  filter.
- RPKM is normalised by library size and genome length, so it is comparable across samples
  and vOTUs, but it is still a *relative* measure - it says nothing about absolute viral
  particle counts.
- The incubation samples themselves were **not** mapped. Comparing abundance in the
  incubation against the environment would need that, and it is a one-line change to the
  sample list in `14_votu_abundance.slurm`.

----------------------------------------------------------------------------------------

----------------------------------------------------------------------------------------

# Step 15. Do our GVMAGs match the ones the original study published? (2026-11-18)

Written by Claude Code on 2026-11-18 (prompt 49 in `PROMPT_LOG.md`). Files:
`scripts/15_fetch_booker_gvmags.sh`, `scripts/group_booker_contigs.py`,
`scripts/15_booker_ani.slurm`, `scripts/summarize_booker_ani.py`.

This is the direct test of the whole reanalysis: Booker et al. 2026 deposited their GVMAGs
in GenBank, so we can ask whether our independent assembly and binning recovered the same
viruses.

#### 15a. The accessions are contigs, not genomes

The paper lists its GVMAGs as **GenBank PV483752-PV483756 and PV637767-PV637835** - 74
accessions. It is tempting to treat those as 74 genomes. They are not. Each record is a
**contig**, and the GenBank title carries an `isolate` field naming the GVMAG it belongs
to:

```
PV483752.1  128444 bp  Nucleocytoviricota sp. isolate GVMAG_022B_bin.8 sequence
PV483756.1   54326 bp  Nucleocytoviricota sp. isolate GVMAG_022B_bin.8 sequence
PV637767.1  351190 bp  MAG: ... isolate vOTU6_CE2_022B_bin14r genomic sequence
PV637800.1    2709 bp  MAG: ... isolate vOTU6_CE2_022B_bin14r genomic sequence
```

Grouping by isolate turns the 74 contigs into **10 published GVMAGs**:

| published GVMAG | their vOTU | contigs | total bp |
|---|---|---|---|
| GVMAG_022B_bin.8 | (none given) | 5 | 436,894 |
| vOTU1_GVMAG_022C_bin.20 | vOTU1 | 3 | 1,249,865 |
| vOTU2_CE2_022B_bin31 | vOTU2 | 2 | 623,527 |
| vOTU2_CE2_022C_bin.14 | vOTU2 | 2 | 623,463 |
| vOTU3_CE2_022B_bin19 | vOTU3 | 4 | 418,504 |
| vOTU3_GVMAG_022C_bin.21 | vOTU3 | 5 | 419,385 |
| vOTU4_GVMAG_022C_bin.27 | vOTU4 | 5 | 451,954 |
| vOTU5_CE2_022C_bin.2 | vOTU5 | 4 | 438,504 |
| vOTU6_CE2_022B_bin14r | vOTU6 | 38 | 687,815 |
| vOTU6_CE2_022C_bin16 | vOTU6 | 6 | 457,346 |

**Why this matters:** comparing our whole multi-contig GVMAGs against individual 2.7 kb
contigs would produce meaningless ANI and alignment fractions, and would have looked like
74 mostly-failed comparisons rather than 10 genome-to-genome ones.
`group_booker_contigs.py` writes one multi-FASTA per isolate so both sides have the same
shape. All 74 contigs were assigned; none was left unparsed.

The download sends **no email or API key** to NCBI. E-utilities does not require one, and
the user's contact details are not ours to hand to an external service.

#### 15b. Method

`fastANI` is run **in both directions**, because it reports the aligned fraction relative
to the *query*. A 200 kb genome wholly contained in a 2 Mb genome aligns over 100% of the
small one but only 10% of the large one, so quoting one direction can make a partial match
look complete or a complete match look partial. `summarize_booker_ani.py` reads both runs
and reports the fraction belonging to the **shorter** genome.

"Same species" uses the **MIUViG rule: >=95% ANI over >=85% of the shorter genome** - the
same rule used to build our own vOTUs in step 09, so the answer is consistent with the rest
of the project.

Job **73934** (astral, 4 cores, 8 GB) finished in **10 seconds**.

#### 15c. Result: all six published vOTUs were recovered

**12 of our 26 GVMAGs are the same species as a published GVMAG**, nearly all at
**99.4-100% ANI**:

| our GVMAG | our vOTU | published GVMAG | their vOTU | ANI | AF (shorter) |
|---|---|---|---|---|---|
| SRR28004154_bin_6 | vOTU17 | vOTU3_CE2_022B_bin19 | vOTU3 | 99.996 | 99.2 |
| SRR28004137_bin_10 | vOTU17 | vOTU3_CE2_022B_bin19 | vOTU3 | 99.999 | 98.7 |
| SRR28004137_bin_32 | vOTU12 | GVMAG_022B_bin.8 | - | 99.999 | 100.0 |
| SRR28004154_bin_8 | vOTU12 | vOTU5_CE2_022C_bin.2 | vOTU5 | 99.989 | 100.0 |
| SRR28004154_bin_25 | vOTU03 | vOTU6_CE2_022C_bin16 | vOTU6 | 99.994 | 97.0 |
| SRR28004137_bin_6 | vOTU03 | vOTU6_CE2_022B_bin14r | vOTU6 | 99.552 | 86.7 |
| SRR28004154_bin_22 | vOTU02 | vOTU1_GVMAG_022C_bin.20 | vOTU1 | 99.994 | 99.5 |
| SRR28004137_bin_21 | vOTU02 | vOTU1_GVMAG_022C_bin.20 | vOTU1 | 99.374 | 91.6 |
| SRR28004137_bin_2 | vOTU08 | vOTU2_CE2_022B_bin31 | vOTU2 | 99.989 | 100.0 |
| SRR28004154_bin_23 | vOTU08 | vOTU2_CE2_022B_bin31 | vOTU2 | 99.986 | 99.5 |
| SRR28004139_bin_22 | vOTU08 | vOTU2_CE2_022C_bin.14 | vOTU2 | 99.712 | 87.4 |
| SRR28004154_bin_7 | vOTU07 | vOTU4_GVMAG_022C_bin.27 | vOTU4 | 99.803 | 87.8 |

At vOTU level this is a clean one-to-one correspondence, and **every one of the paper's six
vOTUs is recovered**:

| our vOTU | their vOTU | ANI range |
|---|---|---|
| vOTU02 | vOTU1 | 99.37-99.99 |
| vOTU08 | vOTU2 | 99.71-99.99 |
| vOTU17 | vOTU3 | 100.00 |
| vOTU07 | vOTU4 | 99.80 |
| vOTU12 | vOTU5 + the unlabelled GVMAG_022B_bin.8 | 99.99-100.00 |
| vOTU03 | vOTU6 | 99.55-99.99 |

**ANI of 99.4-100% means these are the same viruses, reassembled independently.** That is
about as strong a validation of steps 01-09 as this project can get: a different assembler
invocation, a different binning run and a different vOTU clustering all converged on the
same genomes. Note also that our vOTU numbering is **not** theirs - our vOTU17 is their
vOTU3 - so the two numbering schemes must never be used interchangeably.

**Our vOTU12 matches two published genomes**, `GVMAG_022B_bin.8` and
`vOTU5_CE2_022C_bin.2`, both at ~100% ANI. The simplest reading is that those two published
genomes are themselves the same virus, and that the unlabelled `GVMAG_022B_bin.8` is the
other member of their vOTU5 - consistent with it being the only deposited genome without a
vOTU prefix.

#### 15d. What is new in our reanalysis

**13 of our 19 vOTUs have no counterpart in the published set**, and 7 of our 26 GVMAGs
return no fastANI hit at all:

| our GVMAG | our vOTU | length |
|---|---|---|
| SRR28004138_bin_37 | vOTU11 | 671,956 |
| SRR28004144_bin_70 | vOTU06 | 421,742 |
| SRR28004148_bin_93 | vOTU10 | 351,669 |
| SRR34748664_bin_2 | vOTU04 | 348,214 |
| SRR28004151_bin_1 | vOTU14 | 261,816 |
| SRR28004139_bin_7 | vOTU19 | 216,065 |
| SRR28004151_bin_29 | vOTU05 | 203,021 |

Two readings, and this comparison cannot separate them: these are genuinely additional
giant viruses our pipeline recovered, **or** the paper recovered them too and deposited
only a subset. The deposited set is 10 genomes for 6 vOTUs, which is small for a study of
this size, so **the second possibility should be checked against the paper's supplementary
tables before claiming novelty.** Note that `SRR34748664_bin_2` (vOTU04) comes from a
sample outside the original study's accession series altogether, so that one is expected to
be new.

Also worth noting: **vOTU05 and vOTU19 are in this no-match list**, and they are the same
two genomes that carried none of the seven marker genes and were dropped from the step 11
tree (section 10j). vOTU05 now has three strikes against it - no markers, no GVMAGs V2
relative, and no match to the published set - while still being called 100% complete by
CheckV at 203 kb. It should not be reported as a novel giant virus without further work.

Seven more of our GVMAGs hit a published genome **below** the species threshold. Most sit
at ~74.8% ANI over 8-30% alignment, which is background-level similarity between different
viruses. Two are more interesting: `SRR28004137_bin_11` (our vOTU13) at 97.3% ANI over
71.6%, and `SRR28004154_bin_24` (our vOTU09) at 93.9% over 64.0% - both just under the
cut-off against their vOTU6. These are close relatives rather than the same virus, and
vOTU13 is in any case the fragmentary genome flagged in sections 11k and 14d.

#### 15e. Output files

In `output/`, as requested:

| file | what it is |
|---|---|
| `booker_gvmag_ani_summary.csv` | **the main result** - one row per our 26 GVMAGs, with its best-matching published GVMAG, ANI, aligned fraction of the shorter genome, our vOTU and CheckV completeness, and a `same_species_miuvig` verdict. GVMAGs with no match appear with `none_above_fastani_threshold`, because that is a result |
| `booker_gvmag_ani_by_votu.csv` | the same rolled up to vOTU level - the clearest view of the correspondence |
| `booker_gvmag_ani_all_pairs.csv` | every fastANI pair in both directions, with both aligned fractions, for checking the summary or applying a different threshold |

#### 15f. Caveats

- **fastANI reports nothing below roughly 80% ANI.** An absent pair means "no relative
  within ~80% ANI", not "0% identity", and it is not evidence that two genomes are
  unrelated at a deeper taxonomic level. That question belongs to the step 11 tree.
- ANI is computed between **assemblies**, so a fragmentary genome can fail the 85%
  alignment-fraction test while being the same virus. `SRR28004137_bin_11` at 97.3% ANI
  over 71.6% is exactly that situation, and is better described as "probably the same
  virus, incompletely assembled" than as a different one.
- The comparison is only against what was **deposited**. Absence from GenBank is not
  absence from the paper.
- Our vOTU numbers and the paper's vOTU numbers are independent. The mapping in 15c is the
  only valid translation between them.

----------------------------------------------------------------------------------------

# Next steps

1. ~~Wait for the metaSPAdes jobs (step 05b), then run the merge (step 06).~~ **Done
   2026-11-13**: step 05b was cancelled for time for the last two samples (section 05b-fix-c)
   and step 06 merged all 20 samples (jobs 110082, 110084-110102).
2. ~~Run a virus identifier (geNomad `end-to-end` or VirSorter2) on the merged contigs.~~
   **Decided against on 2026-11-12 (prompt 24): not interested in phages.** Note the
   consequence: the only viral selection in this pipeline is GVClass's `d_NCLDV` call on
   bins (step 08), so any bacteriophage is binned and then discarded, and no contig-level
   virus identification happens at all. If that changes, `genomad end-to-end` on the merged
   contigs is the cheapest route, since geNomad and its database are already installed.
3. CheckV runs inside step 08 on the GVMAGs. There is no separate CheckV step on contigs.
4. ~~Step 08 (GVMAGs)~~ **Done 2026-11-16** (sections 08f-08h): 848 bins, 26 GVMAGs in 9
   of 20 samples, mostly Imitervirales; 5 high-quality. Keep `--mem=96G` - several jobs
   peaked at the full 96 GB. Eleven samples have no `_gvmag_summary.tsv` at all, which
   step 09 must tolerate.
5. ~~Step 09 (vOTUs)~~ **Done 2026-11-16** (section 09f, job 111457, 1 m 07 s): 26 GVMAGs
   -> **19 vOTUs**, 6 of them recovered from more than one sample. Against all 8,508 GVMAGs
   V2 references, **none reaches the MIUViG species threshold**, so all 19 are new species;
   12 have a relative above ~80% ANI, 7 have none. A re-run needs only `-c 8 --mem=24G
   -t 30:00`.
6. ~~Step 10 (phylogeny)~~ **Done 2026-11-16** (section 10j, job 111459, 1 h 10 m):
   **1,564-taxon x 37,739-column** alignment (17 of our vOTUs + 1,523 in-order references +
   24 outgroups). vOTU05 and vOTU19 were dropped by `ncldv_markersearch` for having **zero**
   markers - vOTU05 in particular should not be reported as a giant virus without more
   evidence (100% CheckV complete, yet none of the 7 core markers). A re-run needs only
   `--mem=16G -t 4:00:00`, not 128 GB. Average pairwise identity is **14.5%**, so deep
   branches in the tree will be weakly supported.
7. ~~Step 11 (tree)~~ **Done 2026-11-18** (sections 11j-11l, job 111564, 1 d 11 h 46 m):
   `WAG+F+I+R4`, logL -9,852,356.85, **64% of 1,561 branches at >=95 UFBoot**, relative tree
   certainty **0.73**. All 17 vOTU placements supported at >=97. Two things to follow up:
   **(a)** the tree places **vOTU04 with *Asfuvirales*** at 100 support, contradicting
   GVClass's *Imitervirales* call - which was only a 4-vs-2 marker vote, so prefer the tree;
   **(b)** two clades contain **only our vOTUs** with no reference interleaved
   ({vOTU03, vOTU09, vOTU13} and {vOTU07, vOTU12}), i.e. candidate *K. brevis* lineages
   absent from GVMAGs V2 - but vOTU13's membership is suspect at 74.8% missing data.
   **Highest-value next run: step 10 with `--tier skani`** (6,266 taxa), to test whether a
   reference falls inside those clades and dissolves them. The tree is **unrooted** -
   midpoint-root it before interpreting.
8. ~~CoverM prevalence across the field samples~~ **Done 2026-11-18** (Step 14, job
   651382): **none of the 12 incubation-derived vOTUs is detectable in any of the 11
   environmental samples**; the only 7 detections are the environmental-derived vOTUs in
   their own source samples, which are the positive controls. Most likely the incubation
   enriched viruses that are rare *in situ*, but the analysis is detection-limited (14e) -
   deeper sequencing or pooling the monthly samples is the way to test persistence.
   Still unwritten: a check for endogenous viral elements in the *K. brevis* genome
   (prompt 19), and mapping the incubation samples themselves, which is a one-line change
   to the sample list in `14_votu_abundance.slurm`.
9. ~~Copy final results out of scratch~~ **Partly done 2026-11-16** (Step 12, job 111460):
   steps 08 and 09 results are now in `output/` (4.4 GB), excluding the 108 GB of `binning/`
   BAMs, which cannot fit in the 50 GB home quota and are regenerable. Still to do: copy
   the step 10/11 phylogeny results once they exist, and decide whether anything needs to go
   to `~/data10` or off-cluster. `output/` is git-ignored; use `git add -f` for individual
   small result files.

10. **Validation achieved (Step 15, 2026-11-18):** all six of the paper's published vOTUs
    were independently recovered at **99.4-100% ANI**, confirming steps 01-09. Thirteen of
    our nineteen vOTUs have no deposited counterpart; **check the paper's supplementary
    tables before claiming these are novel**, since only 10 genomes for 6 vOTUs were
    deposited and the rest may simply be undeposited.
11. **vOTU05 needs resolving before any write-up.** CheckV calls it 100% complete at
    203 kb, yet it has no giant-virus marker genes (10j), no GVMAGs V2 relative (09f), and
    no match to the published set (15d). Three independent signals say it is not what
    CheckV's completeness implies. vOTU13 and vOTU19 need similar scrutiny.

----------------------------------------------------------------------------------------

# Step 16. Which vOTUs appeared in which samples (2026-11-23)

Written by Claude Code on 2026-11-23 (prompt 59 in `PROMPT_LOG.md`). Files:
`scripts/votu_sample_table.py`, `scripts/plot_votu_sample_matrix.R`.

A vOTU-by-sample presence table, 19 vOTUs x 20 samples, as a CSV and a figure.

#### 16a. "Appeared" means two different things, and the table keeps them apart

This is the whole difficulty of the table. There are two kinds of evidence, and they
cover different samples:

| | what it means | available for |
|---|---|---|
| **ASSEMBLED** | a GVMAG was binned out of that sample's assembly and clustered into this vOTU | **all 20 samples** |
| **DETECTED** | that sample's reads were mapped back and covered >=75% of the vOTU genome | **the 11 environmental samples only** |

Step 14 mapped the environmental samples only, so **the 9 incubation samples have no
mapping evidence at all**. Colouring those cells the same as a tested absence would assert
152 absences that were never tested. They are therefore a separate category - *not
assessed* - with its own grey in the figure and its own label in the CSV. This is the one
thing the figure has to get right, and it is why it is not a plain presence/absence
heatmap.

#### 16b. Output

```bash
python3 -I scripts/votu_sample_table.py \
    --votus output/votus/gvmags_votus.tsv \
    --coverm-dir ~/scr10/data/kbrevis/abundance/per_sample \
    --sra-table metadata/SraRunTable.csv \
    --incubation SRR28004137,SRR28004138,SRR28004139,SRR28004140,SRR28004141,SRR28004142,SRR28004143,SRR28004144,SRR28004154 \
    --out-long output/votu_sample_table.csv \
    --out-matrix output/votu_sample_matrix.csv

module load miniforge3; eval "$(conda shell.bash hook)"; conda activate r-plot-env
Rscript scripts/plot_votu_sample_matrix.R output/votu_sample_table.csv output
```

| file | what it is |
|---|---|
| `output/votu_sample_table.csv` | long form, one row per vOTU per sample, with `status`, `status_label`, `sample_set` and `read_mapped` |
| `output/votu_sample_matrix.csv` | wide form, vOTUs down, samples across, with `A`, `D`, `A+D`, `.` and blank |
| `output/votu_sample_matrix.{pdf,png}` | the figure |

Cell counts: **19 assembled**, **7 assembled and detected**, 202 mapped-but-absent, and
**152 not assessed**.

#### 16c. What it shows

The figure splits cleanly into two blocks with no overlap between them, which is the
result of steps 14 and 15 in one picture:

- **Left (incubation):** 12 vOTUs assembled from 5 of the 9 incubation samples. SRR28004154
  and SRR28004137 carry most of them. Nothing is detected here because nothing was mapped.
- **Right (environmental):** 7 vOTUs, each assembled from **and** detected in exactly one
  sample - vOTU15 in January, vOTU05 and vOTU14 in February, vOTU01/vOTU10/vOTU16 in May,
  vOTU04 in the September coastal-Florida sample. Every other environmental cell is a
  tested absence.

So no vOTU crosses the divide, and within the environmental block no vOTU is found in more
than one sample. Both facts are visible at a glance, and both are the detection limit
described in section 14e rather than a biological turnover signal.

**The obvious gap this figure makes plain:** the incubation half is grey. Mapping the 9
incubation samples against the same 19 vOTUs would fill it in and turn this into a real
presence/absence matrix - it is the same `14_votu_abundance.slurm` with a different sample
list, roughly 30 minutes of compute, and it would let the incubation and environmental
halves be compared directly rather than only by assembly origin.

----------------------------------------------------------------------------------------

# Reproducing this project with one command

Added 2026-11-21 (prompt 53 in `PROMPT_LOG.md`). Files: `run_pipeline.sh`,
`pipeline.conf`, `scripts/setup_new_user.sh`.

The `# Commands` section below lists every command individually, which is what you want
when something has gone wrong or you are running one stage by hand. This section is the
other thing: **three commands that queue the entire analysis.**

```bash
git clone https://github.com/iadanforth/K_brevis_viromes.git
cd K_brevis_viromes
nano pipeline.conf                      # set PROJECT_ROOT, SCRATCH, MAIL_USER
bash scripts/setup_new_user.sh          # adapt the paths baked into #SBATCH directives
bash run_pipeline.sh --dry-run all      # read what it would do
bash run_pipeline.sh all                # queue all 17 stages
```

## How it works

Each stage is submitted with `--dependency=afterok:<the previous stage's job IDs>`, so the
whole pipeline is queued in seconds and Slurm releases each stage only when the one before
it has **succeeded**. `run_pipeline.sh` then exits. It deliberately does not sit on the
front-end waiting: a process that polls for five days would trip the watchdog
(AGENTS.md section 1) and die at logout regardless.

Stages that submit many jobs - 05a, 05b, 06 (20 jobs each) and 08 (an array) - report
every job ID they created, and the next stage depends on **all** of them. Stages that are
not Slurm jobs at all, such as the iTOL annotations and the R plotting, are wrapped in a
small `sbatch --wrap` job so they take their place in the same chain and their logs land
in `logs/` like everything else.

| command | what it does |
|---|---|
| `bash run_pipeline.sh --list` | the 17 stages, which cluster each uses, what each produces |
| `bash run_pipeline.sh --dry-run all` | prints every `sbatch` it would run, submits nothing |
| `bash run_pipeline.sh all` | queues the whole pipeline |
| `bash run_pipeline.sh --from 08 --to 11` | queues part of it |
| `bash run_pipeline.sh --only 14` | queues one stage |
| `bash run_pipeline.sh --status` | what is queued or running, across all five clusters |
| `bash run_pipeline.sh --check` | which stages have produced their outputs |

## When a stage fails

Everything after it stays `PENDING` with reason `DependencyNeverSatisfied`, and Slurm may
cancel it. **That is the point** - it stops a failed assembly from being quietly binned and
phylogenised. Read the failed stage's `.err` in `logs/`, fix the cause, then resume:

```bash
bash run_pipeline.sh --from 08_gvmags
```

Resuming is cheap because the stages skip finished work: 05a/05b/06 skip samples whose
contigs exist, 08 skips samples with a `.done` marker, 14 skips samples with a CoverM
table, and 00 and 07 skip conda environments that already exist.

## What `pipeline.conf` controls

Everything site- and user-specific, in one file: `PROJECT_ROOT`, `SCRATCH`, `ARCHIVE`,
`NODE_TMP`, `MAIL_USER`, the conda module name, the cluster and `sbatch` arguments for
each of the 17 stages, and the analysis parameters that the scripts take on the command
line (`IQTREE_MSET`, `ABUNDANCE_MIN_BREADTH`, and so on). Set
`SKIP_REFERENCE_DOWNLOAD=1` to skip the ~79 GB step 03 when `reference/` is already
populated.

The stage table in `pipeline.conf` also carries the two deliberate resource overrides that
are easy to lose: step 10 is given `--mem=16G` rather than the 128 GB in its `#SBATCH`
(measured peak 3.6 GB), and step 11 is given `-t 48:00:00` because its own 3-day request
exceeds kuro's limit.

## What `scripts/setup_new_user.sh` is for

Two things cannot be read from a config file at run time. Slurm parses `#SBATCH`
directives before any shell runs, so `--output=`, `--error=` and `--mail-user=` cannot
expand `~`, `$HOME` or a variable - and AGENTS.md requires every job's logs to land in
`<project>/logs/`, which forces 44 absolute paths into the scripts. This script rewrites
those and the 18 `--mail-user` lines to match `pipeline.conf`, verifies every script still
parses afterwards, and is idempotent. `--dry-run` shows the edits first.

## Tested

- `--list`, `--check`, `--status` and `--dry-run all` all run clean.
- `--check` was run against the finished project and correctly reports all 17 stages
  present. It caught one of my own wrong assumptions in the process: `05b_assembly`
  produces **18** metaSPAdes assemblies, not 20, because SRR28004142 and SRR28004143 were
  cancelled for time and are carried by step 06's `MV_ONLY_SAMPLES` mechanism instead
  (README 05b-fix-c). The threshold and the comment now say so.
- `scripts/setup_new_user.sh` was tested in a throwaway clone posing as another user: it
  rewrote 23 files and 18 mail lines, left **0** references to the old values, and every
  script still parsed.
- **The live path was tested, not just the dry run.** `bash run_pipeline.sh --only 15`
  submitted the wrapped `fetch_booker` job and the dependent `15_booker_ani` job; both
  completed, and `booker_gvmag_ani_summary.csv` came back **byte-identical** to the
  committed copy.

## What this does and does not give you

It gives a **complete, one-command reproduction on a SciClone account**. Another user
clones, edits three values, runs two commands, and the entire analysis is queued.

It does **not** make the pipeline cluster-independent, and no driver script could. The job
plumbing assumes Slurm with SciClone's per-cluster `slurm/<name>` modules (14 scripts),
`module load miniforge3` (24 scripts), and the `~/scr10` and `~/data10` filesystems
(29 and 6 scripts). `pipeline.conf` centralises the paths and partitions so that porting
means editing one file plus the `module load` lines, rather than hunting through 44
scripts - but on a non-SciClone cluster that editing is still real work.

**The science is portable even where the plumbing is not.** Every tool, version,
parameter and threshold is pinned in `envs/*.yml` and documented in the step sections
above, so the analysis can be rebuilt on any system by someone willing to rewrite the
submission layer.


### Step 17. Re-doing step 14 the way the paper did it: map to each GVMAG, not to one representative per vOTU

Written by Claude Code on 2026-11-23, prompt 60 in `PROMPT_LOG.md`.

#### 17a. Why step 14 was rebuilt

Section 14h established that none of the 12 incubation-derived vOTUs was detectable in any
environmental sample, and that this disagreed with Booker et al. 2026. One possible source
of the disagreement had not been tested: **what we were mapping to.** Their Methods say

> "The relative abundance and coverage of environmental metagenomic reads recruited to
> **each GVMAG** was calculated using coverm."

*Each GVMAG* - not each vOTU. The original step 14 mapped against 19 vOTU
**representatives**: step 09 had clustered the 26 GVMAGs into 19 vOTUs at ≥95% ANI, and
only one genome per cluster was carried forward. The other 7 were discarded from the
mapping reference. That is a reasonable design - it stops two near-identical genomes
stealing reads from each other - but it is a different experiment from the published one.

Some vocabulary, because the two words are easy to run together:

| term | what it is |
|---|---|
| **GVMAG** | *giant-virus metagenome-assembled genome*. One genome, binned out of one sample's assembly. We have 26. |
| **vOTU** | *viral operational taxonomic unit*. A cluster of GVMAGs sharing ≥95% ANI over ≥85% of the shorter genome - roughly "one viral species". Our 26 GVMAGs form 19 vOTUs: 13 singletons, 6 clusters of 2-3. |

So the old reference had 19 entries and the new one has 26. The 7 extra are the non-
representative members of those 6 multi-member clusters.

#### 17b. What the new step 14 does

**`scripts/14_gvmag_reference.sh`** writes one FASTA per GVMAG into
`~/scr10/data/kbrevis/abundance/gvmag_ref/`, named exactly as the GVMAG is named in
`output/votus/gvmags_votus.tsv` (CoverM takes each genome's name from its filename), plus
`gvmag_reference_manifest.tsv` joining each name back to its vOTU, source sample and
source set. It runs on dev: 26 file copies, 14 MB.

By default it writes **all 26** GVMAGs, not only the 19 from the incubation experiments.
The incubation GVMAGs are the ones the question is about - are the viruses that bloomed in
the incubations present in the environment? - but the 7 environmental GVMAGs are kept as
**internal positive controls**. Each was assembled from an environmental sample, so each
*must* be recoverable by mapping that sample's own reads back onto it. Without them, a
result of "the incubation genomes recruit almost nothing" would be indistinguishable from
"the mapping command is broken". `GVMAG_SET=incubation bash scripts/14_gvmag_reference.sh`
builds the 19-genome reference instead.

**`scripts/14_gvmag_abundance.slurm`** is an 11-task array, one task per environmental
sample, and each task runs CoverM **twice**:

| arm | settings | why |
|---|---|---|
| `paper` | CoverM defaults - no read-identity filter, no aligned-length filter | The paper specifies no cutoffs, no version and no settings, so "what they did" can only mean the defaults. This is the primary result. |
| `strict` | `--min-read-percent-identity 95 --min-read-aligned-percent 75` | 95% ANI is the MIUViG species boundary and the same cutoff that defined the vOTUs in step 09, so a read must plausibly come from *this* virus rather than a relative. 75% aligned rejects reads that match only through a short conserved patch. |

Both arms pass `--min-covered-fraction 0`. CoverM's default for that flag is 10, which
silently rewrites a genome's abundance to zero when breadth falls below 10% - that would
bake a presence/absence threshold into the raw numbers. The threshold is applied later, in
`merge_gvmag_abundance.py`, so it stays explicit and can be changed without re-mapping
~100 GB of reads.

Running both arms is the point of the rebuild. The difference between the columns *is* a
measurement: recruitment that disappears once reads must match at 95% identity was never
that genome's to begin with.

Note the scale on those flags. CoverM 0.7.0 wants a **percentage**, not a fraction.
Passing `0.95` is accepted silently and means 0.95% identity - no filtering at all - with
no error to warn you. This was caught by a smoke test before the first run and is recorded
in the script's header.

#### 17c. Reading the output, and one trap that was found the hard way

`scripts/merge_gvmag_abundance.py` reads both arms and writes four files to
`output/abundance/`:

| file | contents |
|---|---|
| `gvmag_abundance_long.tsv` | one row per GVMAG per sample per filter arm, all CoverM metrics |
| `gvmag_abundance.csv` | the same, trimmed to what the figures read |
| `votu_rollup.tsv` | per vOTU per sample, members recombined |
| `gvmag_abundance_matrix.tsv` | genome × sample relative abundance, one block per arm |

The **presence rule** is `covered_fraction ≥ 0.75` - at least 75% of the genome covered by
at least one read (MIUViG, Roux et al. 2019). Breadth, not depth, is the gate. A genome can
show high average coverage while only a few percent of it is covered at all, which is what
happens when a handful of conserved genes attract reads from a relative. Breadth saturates
as `1 − e^−coverage`, so 75% breadth needs only about 1.4× coverage: it is a weak demand on
depth and a strong demand on evenness, which is exactly the right shape of test.

The **rollup** exists because of a cost of following the paper's design. Mapping to 26
individual GVMAGs means the 6 multi-member vOTUs split their recruited reads between
members that are ≥95% identical, so each member alone understates the cluster. The rollup
reports, per vOTU per sample, `read_count_sum` (the split undone), `best_covered_fraction`
(breadth of the best-covered member - breadths are *not* summed, since overlapping members
would double-count), and the representative's own values for a like-for-like comparison
with the old per-representative analysis.

**The trap.** Each array task writes its `paper` table before its `strict` table, so
part-way through the run the `paper` directory is ahead. Merging at that moment produces a
table that looks complete and is quietly wrong: `max breadth across samples` silently drops
the samples that have not landed yet. This happened during job 74013 - a mid-run merge
showed vOTU01 at breadth 1.000 under `paper` and 0.060 under `strict`, which looks like a
dramatic finding and was nothing but `SRR28004148` still being `.part` in the strict
directory, so that genome's own source sample was missing from its strict maximum.
`merge_gvmag_abundance.py` now refuses to run unless both arms cover the identical sample
set, and names the difference.


#### 17d. The result: mapping per GVMAG does not change the answer

Job 74013 on astral, 11 array tasks, 3-7 min each, finished cleanly. Summary over the
11 environmental samples, with the maximum across samples for each genome:

| | `paper` arm | `strict` arm |
|---|---|---|
| detections (breadth ≥75%) | **7** | **7** |
| of those, incubation GVMAGs | **0** | **0** |
| minimum breadth among the 7 | 0.9998 | 0.9991 |
| incubation GVMAGs: median reads per genome-sample pair | 347 | 1 |
| incubation GVMAGs: best breadth achieved by any of the 19 | 0.196 | 0.097 |
| total reads recruited by all 19 incubation GVMAGs, all samples | 424,857 | 5,914 |

**The positive controls pass.** All 7 environmental GVMAGs reach ≥99.9% breadth under
*both* filter settings. The mapping command works, and it works at 95% identity. That is
what makes the incubation result interpretable rather than merely negative.

**No incubation GVMAG comes close.** The best any of the 19 manages is 19.6% breadth, and
that is in the unfiltered arm. Lifting the identity filter multiplies their recruited reads
by 72× in total, and moves nothing above 20% breadth. The 7 detections are the same 7 as
before the rebuild, so **mapping to each GVMAG rather than to one representative per vOTU
does not change the conclusion of section 14d.**

**The clearest single illustration** is `SRR28004137__SRR28004137_bin_11` (our vOTU13,
from an incubation) in sample `SRR34748664`, in the unfiltered arm:

| reads recruited | genome breadth | mean depth | relative abundance |
|---|---|---|---|
| 264,155 | **0.012** | 76.6× | **0.358 %** |

A quarter of a million reads, piled 77 deep onto **1.2% of the genome**, with the other
98.8% untouched. The relative abundance that comes out, 0.36%, sits squarely inside the
range Booker et al. describe as presence (*"~0.1-0.8"*). Under the 95% identity filter the
same genome-sample pair keeps almost nothing.

This is the whole argument in one row. Relative abundance cannot distinguish a virus that
is present from a conserved gene family that is present in something else; only breadth
can, and breadth is unambiguous here. It does **not** prove Booker et al. wrong - we cannot
see their settings, their read set differs, and a different cutoff is a defensible choice,
not an error. What it does show is that the quantity their Methods describe can reach
reportable values without the genome being there, so the disagreement between our
section 14d and their Fig. 4 is explained by the measurement, not by a difference in the
underlying data.

**The vOTU-by-sample matrix is unchanged too.** `output/votu_sample_matrix.png`, rebuilt
from the per-GVMAG mapping: 19 ASSEMBLED, 7 ASSEMBLED+DETECTED, 202 tested absences, 152
not assessed, and **0 DETECTED-only**. Per-GVMAG mapping makes it possible for a vOTU to
turn up in a sample it was not assembled from - that is the cell the old per-representative
design could not really test - and no vOTU does.

#### 17e. What this does not settle

- **The 9 incubation samples are still not read-mapped.** The grey half of the matrix is
  unchanged. Mapping them against the same 26 GVMAGs is the same array with a different
  sample list, about 30 minutes, and would turn the matrix into a true presence/absence
  table.
- **Two assemblies are still missing.** `SRR28004142` and `SRR28004143` (section 05b-fix-c)
  have never produced contigs. Neither contributed a GVMAG, so they cannot change whether
  the existing 19 incubation GVMAGs are detected - but they could add new ones. 695880 is
  running on hi03; 700722 is queued behind it.
- **The cell-count regression is still blocked on data**, not on code.
  `scripts/plot_votu_vs_cells.R` now reads `votu_rollup.tsv` and fits all 19 vOTUs, but
  `output/abundance/karenia_cell_counts.csv` is an empty template and the script exits
  with a message saying so. The counts are not in the paper's text; they need
  Supplemental Data 14 or the FWC/FWRI HAB Monitoring Database.
- **Read-splitting within multi-member vOTUs is real but not load-bearing here.** The 6
  vOTUs with 2-3 members do divide their reads, which is why `votu_rollup.tsv` exists.
  Recombining them changes no detection call: the summed read counts for incubation vOTUs
  remain three orders of magnitude below the detected ones.


### Step 18. The figure that corresponds to the paper's Fig. 4 (2026-11-24)

Written by Claude Code on 2026-11-24, prompt 62 in `PROMPT_LOG.md`. File:
`scripts/plot_incubation_votu_timeseries.R`. Outputs
`output/abundance/incubation_votu_timeseries.{pdf,png}` and
`output/abundance/incubation_votu_abundance.csv`.

Step 17's figure answered a methodological question - how much does the identity filter
change the numbers, and does breadth agree? - across all 26 GVMAGs. That is not the figure
the analysis is for. **This one is**: collection date on the x axis, and the abundance of
the vOTUs recovered from the incubation experiments measured in the environmental
metagenomes, which is the comparison Booker et al. 2026 make in their Fig. 4.

| choice | what it is | why |
|---|---|---|
| series | the **12** vOTUs assembled from an incubation sample | The 7 environmental vOTUs are excluded. The question is whether the viruses that bloomed in the incubations are out in the Gulf; an environmental vOTU recruiting its own source sample's reads answers nothing. |
| y | relative abundance, % of the sample's reads, log scale | The quantity the paper reports. Log because the values run from 5×10⁻⁶% to 0.36%, five orders of magnitude. |
| mapping | CoverM **defaults** - no read-identity filter, no aligned-length filter | A lenient threshold is what allows an environmental read to recruit to a genome assembled from a different experiment, and the paper specifies no cutoffs. This is step 14's `paper` arm. |
| per vOTU | members summed via `votu_rollup.tsv` | Six vOTUs hold 2-3 near-identical GVMAGs that split the reads between them; summing recovers the cluster's true recruitment. |
| lines | join **all 11** samples in date order | All 11 are one survey (station EH25) - see 18a. Monthly from Dec 2020 to Aug 2021, then two later points; the gap between them is a gap in sampling, not a boundary between studies. |
| colour | Okabe-Ito recycled, paired with line type | Twelve series need twelve keys and the colour-blind-safe palette has eight; the last four vOTUs are dashed. Colour and line type are both mapped to `votu`, so they merge into one legend. |

**What it shows.** Every one of the 12 incubation vOTUs recruits environmental reads at
every time point - there are no zeros - across a band from about 10⁻⁵% to 10⁻²%, rising
roughly an order of magnitude between the February 2021 minimum and the April-August bloom
period. vOTU13 is the outlier throughout, peaking at **0.358%** in the September 2021
sample.

**The caveat, stated once and carried in the figure's own caption.** Genome breadth never
exceeds **19.6%** for any point on this plot, and is below 2% for most of them, so these
are read-recruitment values rather than confirmed genome detections - the same reads could
be coming from a relative's conserved genes. Section 17d works through what that means and
why it does not amount to showing the paper is wrong. The figure is the paper's
measurement, made faithfully; the breadth column in
`output/abundance/incubation_votu_abundance.csv` is there so the reader can check it.

To redraw it with the 95% identity filter instead, pass `strict` as the third argument.

#### 18a. A metadata correction: the two "coastal Florida" samples are not a separate study

Added 2026-11-24, prompt 63 in `PROMPT_LOG.md`.

`SRR34748664` (2021-09-29) and `SRR34748663` (2022-03-08) carry `isolation_source` =
**coastal Florida** in `metadata/SraRunTable.csv`. Every figure in this project used to
treat them as a separate study on that basis: drawn as open triangles, never joined to the
Gulf of Mexico series, on the reasoning that a line through all eleven would assert a
continuous 15-month survey that did not exist.

**That label is wrong.** The user confirmed on 2026-11-24 that these two runs come from the
same survey as the other nine and were mislabelled when submitted to the SRA. (The README's
own sample table at the top already describes both as "EH25 field time series", which is
consistent with this.) So all 11 environmental samples are **one series**: one marker
shape, one line through all of them in date order.

**Where the correction lives.** In the three plotting scripts, each as a named constant:

```r
MISLABELLED <- c("SRR34748664", "SRR34748663")
```

It is deliberately **not** applied upstream. `metadata/SraRunTable.csv`,
`gvmag_abundance_long.tsv` and `votu_rollup.tsv` still carry the SRA's own record, so the
raw metadata stays recoverable and a reader can see exactly what was changed and why. If
the submission is ever corrected at NCBI, deleting those constants restores agreement.
`incubation_votu_abundance.csv` carries both: `survey` (corrected) and
`sra_isolation_source` (what the SRA says).

**Figures rebuilt under the correction:** `incubation_votu_timeseries.{pdf,png}`, and at
the time also `gvmag_abundance_timeseries.{pdf,png}` (since deleted - see 18b).
`plot_votu_vs_cells.R` lost its `shape = survey` mapping, since there is no sample-origin
distinction left for shape to carry. No numbers changed - only which points are joined and
how they are drawn.

#### 18b. The step-17 comparison figure was deleted (2026-11-24)

Prompt 64 in `PROMPT_LOG.md`. `output/abundance/gvmag_abundance_timeseries.{pdf,png}` -
the 2x2 panel of abundance and breadth across both filter settings for all 26 GVMAGs - was
deleted at the user's request, and the `Rscript scripts/plot_gvmag_abundance.R` call was
removed from `run_pipeline.sh` and from the Commands section, so a re-run does not put it
back.

`scripts/plot_gvmag_abundance.R` itself is **kept**, with a note at the top saying it is no
longer part of the default run. It is the only thing that draws the lenient-vs-strict and
abundance-vs-breadth comparison, and the table it reads,
`output/abundance/gvmag_abundance.csv`, is still written by step 14, so the figure can be
rebuilt at any time with one command:

```bash
Rscript scripts/plot_gvmag_abundance.R output/abundance/gvmag_abundance.csv output/abundance
```

Nothing in the analysis depended on that figure. The numbers it displayed are in
`gvmag_abundance_long.tsv` and are reported in section 17d; the headline figure is
`incubation_votu_timeseries.{pdf,png}` from step 18.


#### 18c. Taxonomy columns in `incubation_votu_abundance.csv` (2026-11-24)

Prompt 67 in `PROMPT_LOG.md`. The CSV now carries two extra columns, so the figure can be
facetted by taxonomy:

| column | source | values among the 12 incubation vOTUs |
|---|---|---|
| `order` | `o_` from `gvclass_taxonomy_majority` | Imitervirales (11), Pandoravirales (1) |
| `family` | `f_` from the same string, taken from the vOTU's **representative** | IM_01 (4), IM_16 (2), IM_14 (1), PV_05 (1), unclassified (4) |

Both come from GVClass's assignment in `output/votus/all_gvmag_summaries.tsv`, parsed out of
the lineage string by `pull_rank()` in `scripts/plot_incubation_votu_timeseries.R`.

A rank GVClass could not assign is present but **empty** in that string (`f_;`), so an empty
value is written as `"unclassified"` rather than `NA` - it means "GVClass looked and could
not say", which is a finding, not missing data.

Why the two ranks are derived differently:

- **Order is unanimous** within every vOTU - all members of a cluster agree - so it is taken
  directly, and the script **stops with an error** if that ever stops being true rather than
  quietly reporting whichever member sorted first.
- **Family is not unanimous.** A less complete member often falls back to blank while the
  representative resolves to a family, so family is taken from the representative, the
  genome step 09 chose to stand for the cluster.

**Faceting on `order` is lopsided**: 11 of the 12 vOTUs are Imitervirales and only vOTU19 is
Pandoravirales, so one panel holds eleven lines and the other holds one. `family` splits
them more evenly and is probably the more informative facet.

### Step 19. A one-row-per-vOTU origin table (2026-11-24)

Written by Claude Code on 2026-11-24, prompt 65 in `PROMPT_LOG.md`. File:
`scripts/votu_origin_table.py`. Output: `output/votu_origin.csv`.

The plainest summary of what the pipeline recovered: one row per vOTU, the sample it came
from, and whether that sample was an incubation or an environmental one.

| column | meaning |
|---|---|
| `votu` | the vOTU ID assigned in step 09 |
| `sample` | the SRA run(s) whose assembly produced the GVMAG(s) in this vOTU |
| `sample_type` | `incubation` (one of the 9 CH1/CH2 experiment runs) or `environmental` |

**"Found in" means assembled from** - the sample whose assembly produced the GVMAG that was
binned, quality-checked by CheckV/GVClass and clustered into this vOTU. That is how a vOTU
enters the study at all. For the 7 environmental vOTUs it happens to be the same answer as
read mapping gives, because no vOTU is detected in a sample it was not assembled from
(section 17d); for the 12 incubation vOTUs the two questions have different answers, and
this table answers the first one.

**Six vOTUs were recovered from more than one sample** - the same virus assembled
independently in several incubations - so those rows list every sample rather than
splitting into extra rows, which keeps the table at one row per vOTU as asked. All six are
incubation-only, so `sample_type` is never mixed within a row. The script **checks** this
rather than assuming it: if a vOTU ever turns up assembled from both an incubation and an
environmental sample, it exits with an error naming the vOTU instead of silently reporting
one of the two kinds. That could genuinely happen if the two outstanding assemblies
(`SRR28004142`, `SRR28004143`) add GVMAGs to an existing environmental vOTU.

Result: **19 vOTUs, 12 incubation and 7 environmental.** `SRR28004154` and `SRR28004137`
between them account for 9 of the 12 incubation vOTUs, and 5 vOTUs were recovered from
that pair together - which is the same clustering that shows up as the multi-member vOTUs
in `votu_rollup.tsv`.


### Step 20. A publication-quality phylogeny figure, and why the tree holds 17 vOTUs and not 19 (2026-11-24)

Written by Claude Code on 2026-11-24, prompt 66 in `PROMPT_LOG.md`. File:
`scripts/plot_tree_figure.py`. Output: `output/Kbrevis_virome_tree.{png,pdf}`.

#### 20a. Why there are 17 vOTUs on the tree, not 19

**vOTU05 and vOTU19 are not in the tree at all** - not unplaced, not collapsed, absent.
`ncldv_markersearch` found **none of the seven** core *Nucleocytoviricota* marker genes in
either, and it silently drops a taxon with zero markers before the alignment is built. They
were the only 2 of all 1,566 input taxa to be dropped, which is why 1,566 went in and 1,564
came out (section 10j). Straight from `kbrevis_gvmags.table.tsv`:

| genome | GVOGm0013 | GVOGm0023 | GVOGm0054 | GVOGm0172 | GVOGm0461 | GVOGm0760 | GVOGm0890 |
|---|---|---|---|---|---|---|---|
| `vOTU05__SRR28004151_bin_29` | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| `vOTU19__SRR28004139_bin_7` | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

Every other vOTU scores 1 on at least one marker. The reason is mechanical: a phylogeny is
built from an alignment of shared genes, so a genome sharing **no** gene with the others
has nothing to align and no position to occupy. This is not a threshold that could be
loosened.

**It is also a result worth reporting, not just a technicality.** vOTU05 is the awkward
one: CheckV called it 100% complete at 203 kb, it has no relative in GVMAGs V2, and it
carries none of VLTF3, the A32 packaging ATPase, DNA topoisomerase II, TFIIB, DNA
polymerase B, RNA polymerase alpha or SNF2 helicase. A genuinely complete giant-virus
genome should carry most of those. Either it is a mis-binned set of contigs, or it is real
but too divergent for the GVOG HMMs to detect. GVClass does place it in Imitervirales, but
on very little: 2 marker votes, 1 unique GVOG8, `order_completeness` 11.4% and the largest
`avgdist` of any of our 26 GVMAGs (3.15). (An earlier version of this paragraph said
GVClass assigned it no order at all; that was wrong - see the correction in section 10j.)
**It should not be reported as a giant-virus genome without further evidence.** vOTU19 is the
CheckV "Not-determined" genome, so its failure is less surprising.

The figure states this in its own legend key - "K. brevis vOTU (17 of 19)" - and the script
re-derives the count from the tree rather than hard-coding it, printing which vOTUs are
missing on every run.

#### 20b. The figure

The previous `Kbrevis_virome_tree.png` was exported by hand from iTOL. That made it
unreproducible, and its type size and line weights were whatever iTOL's defaults gave: a
legend at roughly 4 pt and hairline red branches, both too weak for publication.
`scripts/plot_tree_figure.py` now draws the same tree from
`output/tree/kbrevis_gvmags.treefile`, so every visual parameter is under version control.

No tree-plotting library is installed in this project's environments (no `ete3`, no
`Bio.Phylo`, no `ggtree`/`ape`), so the script carries its own Newick parser and radial
layout in about 60 lines of matplotlib. That is a smaller cost than adding a dependency,
and the Newick the script has to read is IQ-TREE's own output, which is small and
well-defined.

What changed, and why:

| | before | now |
|---|---|---|
| legend | ~4 pt, translucent frame | 15 pt text, 17 pt bold title, opaque frame |
| vOTU branches | hairline red, traced only at the tip | saturated `#D7191C`, 1.9 pt, drawn over a white casing so they stay legible against dense grey neighbours |
| finding a vOTU | hunt for a red hairline among 1,547 grey tips | red leader line from tip to ring, tall red wedge in the ring, and a 13 pt bold label outside it |
| labels | none | all 17 vOTUs labelled, rotated to read outward |
| scale bar | none | 0.5 substitutions per site, with end ticks |
| reproducible | no | `python3 -I scripts/plot_tree_figure.py` |

**Red marks each vOTU's own terminal branch, not its lineage back to the root.** Tracing
ancestors was tried first and looked worse: the 17 vOTUs sit in many different clades, so
their paths merge low down and turn most of the backbone red, hiding exactly what the
colour is meant to pick out. The leader lines do the work of drawing the eye instead.

**Labels are nudged apart, branches are not.** Several vOTUs are near neighbours in the
tree - vOTU03, vOTU09 and vOTU13 are each other's sisters - so their labels would overprint
at the true tip angle. A relaxation pass spreads labels to a minimum 3.6 degree separation
while keeping tree order; the leader line still points at the real position, so nothing is
misplaced.

**Placements were checked, not assumed.** Every red tip in the figure lands in the order
band that `output/tree/votu_tree_context.tsv` predicts from its sister taxa: vOTU14 in
Pimascovirales, vOTU15 in Algavirales, vOTU04 beside Asfuvirales, the remaining twelve in
Imitervirales, and vOTU03/vOTU09/vOTU13 adjacent to one another.

#### 20c. GVClass order and tree placement are different questions, and can disagree

GVClass classifies a genome by comparing its marker genes to a reference database. The step
11 tree places it by maximum likelihood on a concatenated alignment. They answer related but
distinct questions, and for three vOTUs they do not agree:

| vOTU | GVClass order | nearest sister in the tree | note |
|---|---|---|---|
| vOTU04 | Imitervirales (4 votes, 2 for Asfuvirales) | Asfuvirales | the vote was already split; the tree sides with the minority |
| vOTU15 | *unclassified* - GVClass left `o_` blank after a 3/2/1 split across Algavirales, a eukaryote and split_Imitervirales | Algavirales | the tree resolves what the vote could not |
| vOTU19 | Pandoravirales | **not in the tree** | it carries none of the 7 phylogeny markers, but GVClass still classified it from a different gene set (2 unique GVOG8) |

vOTU19 is the instructive case: **a genome can be classifiable and still be unplaceable.**
GVClass and `ncldv_markersearch` use different gene sets, so "GVClass gave it an order" and
"it has enough markers to go in the tree" are not the same test. That is why the CSV's
`order` column has a value for vOTU19 even though section 20a explains it is absent from the
figure in `output/Kbrevis_virome_tree.png`.



#### 09h. A `sample_set` column in `all_gvmag_summaries.tsv` (2026-11-24)

Prompt 68 in `PROMPT_LOG.md`. `output/votus/all_gvmag_summaries.tsv` now carries a
`sample_set` column, inserted immediately after `sample`, saying whether that GVMAG came
from an **incubation** experiment or an **environmental** sample. Of the 26 GVMAGs,
**19 are incubation and 7 environmental**, which matches the split used everywhere else in
the project.

The whole analysis turns on that distinction - the question is whether the viruses that
bloomed in the incubations are present in the environment - so carrying it in the table
means each consumer no longer has to re-derive it from its own copy of the accession list.

**Written by the pipeline, not patched in afterwards.** The change is in the awk loop in
`scripts/09_gvmag_votus.slurm` that builds the combined table, so a re-run produces it. The
committed copy was regenerated on dev by running that same loop over the nine per-sample
`*_gvmag_summary.tsv` files still in scratch (26 rows, no compute), and checked: with the
new column removed it is **byte-identical** to the previous committed file, so nothing but
the added column changed.

The script also now refuses to run if its hard-coded incubation list does not hold exactly
9 accessions, which catches a typo rather than letting it silently reclassify a sample. That
list is duplicated in `pipeline.conf` (`INCUBATION_SAMPLES`), `scripts/14_gvmag_reference.sh`
and `scripts/votu_origin_table.py`; **they must be kept identical**, and the count check is
the only thing currently guarding that.

Both existing readers of this table - `scripts/votu_cluster.py` (`csv.DictReader`) and
`scripts/plot_incubation_votu_timeseries.R` (`read_tsv`) - address columns by name, so
inserting a column does not disturb them. Verified by re-running the R script: it produced
an unchanged `incubation_votu_abundance.csv`.


#### 10k. A per-vOTU marker-gene table (2026-11-24)

Prompt 69 in `PROMPT_LOG.md`. Files: `scripts/votu_marker_table.py`,
`output/votus/votu_marker_genes.csv`, `output/tree/marker_occupancy_raw.tsv`.

Step 10 already works out which of the seven marker genes each genome carries -
`ncldv_markersearch` writes an occupancy table as a by-product - but it landed in scratch,
which is purged after 90 days without access, and it is keyed by the long mangled alignment
taxon name. The raw table is now copied into `output/tree/marker_occupancy_raw.tsv` so the
repo does not depend on scratch, and `scripts/votu_marker_table.py` turns it into a readable
per-vOTU table:

| vOTU | set | 0013 | 0023 | 0054 | 0172 | 0461 | 0760 | 0890 | n | in tree |
|---|---|:-:|:-:|:-:|:-:|:-:|:-:|:-:|---|---|
| vOTU01 | environmental | Y | Y | Y | Y | Y | Y | Y | 7 | yes |
| vOTU02 | incubation | Y | Y | Y | Y | Y | Y | Y | 7 | yes |
| vOTU03 | incubation | Y | Y | Y | Y | Y | Y | Y | 7 | yes |
| vOTU04 | environmental | Y | Y | Y | . | Y | Y | Y | 6 | yes |
| **vOTU05** | environmental | . | . | . | . | . | . | . | **0** | **no** |
| vOTU06 | incubation | Y | Y | Y | Y | Y | Y | Y | 7 | yes |
| vOTU07 | incubation | Y | Y | Y | Y | Y | Y | Y | 7 | yes |
| vOTU08 | incubation | Y | Y | Y | Y | Y | Y | Y | 7 | yes |
| vOTU09 | incubation | Y | Y | Y | Y | Y | Y | . | 6 | yes |
| vOTU10 | environmental | Y | Y | . | Y | Y | Y | Y | 6 | yes |
| vOTU11 | incubation | Y | Y | Y | Y | Y | Y | Y | 7 | yes |
| vOTU12 | incubation | Y | Y | Y | Y | Y | Y | Y | 7 | yes |
| vOTU13 | incubation | Y | . | . | . | . | Y | . | 2 | yes |
| vOTU14 | environmental | Y | Y | Y | Y | Y | Y | . | 6 | yes |
| vOTU15 | environmental | Y | . | Y | Y | Y | Y | Y | 6 | yes |
| vOTU16 | environmental | Y | Y | Y | Y | Y | Y | Y | 7 | yes |
| vOTU17 | incubation | Y | Y | Y | Y | Y | Y | Y | 7 | yes |
| vOTU18 | incubation | Y | Y | Y | Y | Y | Y | Y | 7 | yes |
| **vOTU19** | incubation | . | . | . | . | . | . | . | **0** | **no** |

**11 of 19 carry all seven.** Two carry none and are therefore absent from the tree
(section 20a). **vOTU13 is the one to look at twice**: it has only 2 of 7, so although it is
on the tree, its position rests on two genes and should not be trusted the way a 7-marker
placement can be. It is also the biggest read recruiter among the incubation vOTUs in the
step 18 figure, which makes its weak placement worth knowing.

Two things the table does **not** say:

1. **The columns have no gene names**, only GVOG accessions. See the caveat added to
   section 10j: the names cannot be recovered from anything in this project, and the
   familiar list that appears elsewhere in this README is unverified.
2. **It describes one genome per vOTU** - the step 09 representative, which is the genome
   step 10 put into the alignment. A non-representative member of one of the 6 multi-member
   vOTUs could carry a marker its representative lacks; that is untested.


#### 09i. A `sample_set` column in `gvmags_votus.tsv` (2026-11-24)

Prompt 70 in `PROMPT_LOG.md`. The same change as 09h, applied to the other step 09 table:
`output/votus/gvmags_votus.tsv` now carries `sample_set` as **column 4**, straight after
`sample` - 19 incubation GVMAGs, 7 environmental. It is written by `scripts/votu_cluster.py`,
so a re-run of step 09 produces it.

**This one needed care, because two consumers read the file by column POSITION**, and a
positional reader fails silently when a column moves - it does not error, it just reads the
wrong field:

| consumer | what it did | what would have happened |
|---|---|---|
| `scripts/14_gvmag_reference.sh` | `while IFS=$'\t' read -r votu gvmag sample isrep rep ...` | `isrep` would have received `"incubation"`, so no genome would ever match `"yes"` |
| `scripts/10_gvmag_phylogeny.slurm` | `awk '$2==g && $4=="yes"'` | `$4` is no longer `is_representative`, so **every tree tip would have been named `noVOTU__...`** |

Both are fixed - the bash read now takes 11 variables and `awk` tests `$5` - and both were
exercised against the real file before committing. The other four consumers
(`votu_cluster.py`'s own output aside: `votu_origin_table.py`, `votu_marker_table.py`,
`make_itol_annotations.py`, `summarize_booker_ani.py`, and `plot_incubation_votu_timeseries.R`)
all address columns by name and were re-run to confirm unchanged output.

**One copy of the incubation list was removed.** `14_gvmag_reference.sh` used to carry its
own set of 9 accessions and re-derive the split; it now just reads the `sample_set` column.
`votu_cluster.py` takes the list as `--incubation`, which `09_gvmag_votus.slurm` passes from
the single list it already checks, and refuses to run unless it holds exactly 9 accessions.
That takes the number of independent copies from four down to three (`pipeline.conf`,
`09_gvmag_votus.slurm`, `votu_origin_table.py`), with `votu_cluster.py`'s literal a
documented fallback for standalone runs.

As in 09h, the committed table was regenerated rather than re-run on the cluster, and
checked: with the new column removed it is **byte-identical** to the previous file.


----------------------------------------------------------------------------------------

# Commands

Every command needed to reproduce this project from nothing, in order. Written
2026-11-20 (prompt 52 in `PROMPT_LOG.md`) and checked against the scripts as they stand on
that date: every script named here exists, every flag was verified against the script's own
argument parser, and every `#SBATCH` default quoted below was read from the file.

This section is the **copy-and-paste list**. The reasoning behind each choice -
why a tool, why a threshold, what the output means - is in the numbered step sections
above, and the shorter narrative version is under "Rebuilding everything from scratch".

**Read these five things first.**

1. **Run everything from `~/K_brevis_viromes` on dev.** The step 01-03 wrappers call their
   worker script by a relative path (`bash ./scripts/02_fastp.sh`) and `sbatch` starts a
   job in the directory it was submitted from. Submitting from elsewhere fails.
2. **Wait for each step to finish before starting the next.** Nothing here polls or blocks.
   Check with `squeue` (see "Monitoring" at the end).
3. **The `-p` partitions below are what worked on the dates in the step sections, not
   permanent advice.** Before each submission run
   `~/hpc/best_cluster.sh -c CORES -m MEM_GB -H HOURS` and override `-p` and `-t`. Only
   `05_metaspades*.slurm` and `06_merge_assemblies.slurm` carry a partition (`-p serial`)
   in the file; everything else takes it from the command line.
4. **Several commands deliberately override the script's own `#SBATCH` values**, because
   the originals were over-requested and later measured. Those are flagged inline. Do not
   "tidy" them away.
5. **Disk**: this produces about **5.3 TB** in scratch. `raw/` 2.0 TB, `trim/` 1.6 TB,
   `clean/` 1.1 TB, `contigs_no_v/` 322 GB, `gvmags/` 115 GB, `contigs/` 90 GB,
   `reference/` 79 GB. Scratch is never backed up and is purged after 90 days without
   access.

## 0. Directories, repository, environments, metadata

```bash
# --- the two directories, and the code -------------------------------------------------
mkdir -p ~/scr10/data/kbrevis
cd ~
git clone https://github.com/iadanforth/K_brevis_viromes.git
cd ~/K_brevis_viromes

# creates the remaining scratch and repo folders (each script also makes its own, so this
# is belt-and-braces rather than strictly required)
bash scripts/00_directory_setup.sh

# --- the sample table step 01 reads ----------------------------------------------------
cp metadata/SraRunTable.csv ~/scr10/data/kbrevis/

# --- conda environments ----------------------------------------------------------------
# Builds the seven environments steps 01-06 need (sra-toolkit, datasets, seqtk,
# host-screen, metaspades, seqkit, fastp) and then verifies every tool resolves.
# Steps 07 and 14 build their own environments later; do not duplicate them here.
( module purge; module load slurm/astral; sbatch -p batch --exclude=as01 scripts/00_create_envs.slurm )
```

**If you are not `iadanforth`**, also run these two substitutions before anything else. 23
scripts contain the literal path and 17 contain the e-mail; 44 of the path occurrences are
`#SBATCH` log directives, which cannot use `~` or `$HOME` and so must be absolute:

```bash
grep -rl '/sciclone/home/iadanforth' scripts/ \
  | xargs sed -i 's|/sciclone/home/iadanforth|/sciclone/home/YOUR_USERNAME|g'
grep -rl 'iadanforth@wm.edu' scripts/ \
  | xargs sed -i 's|iadanforth@wm.edu|YOUR_ADDRESS@wm.edu|g'
```

## 1-6. Reads, references, cleaning, assembly, merging

```bash
cd ~/K_brevis_viromes

# --- 01. download 20 SRA runs -> raw/  (~2.0 TB) ---------------------------------------
( module purge; module load slurm/bora; sbatch scripts/01_download_files.slurm )

# --- 02. adapter/quality trimming -> trim/  (~1.6 TB) ----------------------------------
# fastp 1.3.7 from fastp-env; the wrapper activates it.
( module purge; module load slurm/bora; sbatch scripts/02_fastp.slurm )

# --- 03 + 03b. host and contaminant genomes -> reference/  (~79 GB) --------------------
# SKIPPABLE ON A RE-RUN - see the note below.
( module purge; module load slurm/bora; sbatch scripts/03_download_refs.slurm )
( module purge; module load slurm/bora; sbatch scripts/03b_download_phix.slurm )

# --- 04. bowtie2 host removal -> clean/  (~1.1 TB) -------------------------------------
# Also builds the 42 GB bowtie2 index into reference/, but only if it is not already there.
( module purge; module load slurm/bora; sbatch -p hima -C hi --exclude=hi04,hi05,hi07 \
    scripts/04_clean_refs.slurm )

# --- 05a. metaviralSPAdes -> contigs/  (90 GB; 20 separate jobs) -----------------------
# The submitter skips any sample whose contigs already exist, so it is safe to re-run.
( module purge; module load slurm/kuro; bash scripts/05_metaspades.sh )

# --- 05b. metaSPAdes -> contigs_no_v/  (322 GB; 20 jobs) -------------------------------
( module purge; module load slurm/kuro; bash scripts/05_metaspades_nov.sh )

# --- 06. merge the two assemblies per sample -> contigs_merged/  (6.6 GB; 20 jobs) -----
( module purge; module load slurm/kuro; bash scripts/06_merge_assemblies.sh )
```

**To skip step 03 on a re-run**, keep `~/scr10/data/kbrevis/reference/` when you clear
scratch and simply do not submit `03_download_refs.slurm` or `03b_download_phix.slurm`.
This is safe and saves both the ~79 GB download and several hours of `bowtie2-build`:
step 04 needs six FASTAs (`dog_ref_genome.fna`, `cat_ref_genome.fna`,
`mouse_ref_genome.fna`, `human_ref_genome.fna`, `phix174.fna`, `kbrevis_ref.fna`) and exits
with an error naming any that are missing, and it rebuilds the index only when
`bowtie2-inspect` cannot find one. Check before relying on it:

```bash
ls -la ~/scr10/data/kbrevis/reference/*.fna
ls ~/scr10/data/kbrevis/reference/host_index.*.bt2l | wc -l    # expect 6
```

## 7. One-time setup for the GVMAG pipeline

```bash
# Builds 6 more conda environments, installs GVClass + pixi, and downloads the CheckV v1.5,
# geNomad and GVMAGs V2 databases into ~/data10/databases (about 5 GB of environments and
# several GB of databases).
#
# GVMAGS_V2_FULL=1 IS REQUIRED. Without it the 1.18 GB GVMAGs V2 genome/protein bundle is
# not downloaded and STEP 10 FAILS much later with a confusing error.
( module purge; module load slurm/astral; sbatch -p batch -c 4 --exclude=as01 --mem=16G \
    -t 6:00:00 --export=ALL,GVMAGS_V2_FULL=1 scripts/07_gvmag_setup.slurm )
```

## 8-11. Binning, GVMAGs, vOTUs, alignment, tree

```bash
# --- 08. map, bin, classify -> gvmags/  (115 GB; array, <=10 concurrent) ---------------
# Keep --mem=96G: several samples peaked at the full 96 GB.
( module purge; module load slurm/kuro; bash scripts/08_gvmags.sh -p serial )

# --- 09. cluster GVMAGs into vOTUs -> votus/  (16 MB; ~1 minute) -----------------------
# --with-refs also compares each vOTU to the 8,508 GVMAGs V2 reference genomes.
# 16c/64G is far more than needed (measured: 14.3 GB, 1 m 07 s); -c 8 --mem=24G -t 30:00 suffices.
( module purge; module load slurm/kuro; sbatch -p serial -c 16 --mem=64G -t 4:00:00 \
    scripts/09_gvmag_votus.slurm --with-refs )

# --- 10. marker genes + concatenated alignment -> phylogeny/  (629 MB; ~1 h 10 m) ------
# OVERRIDE THE MEMORY: the script asks for 128 GB, the measured peak is 3.6 GB.
( module purge; module load slurm/kuro; sbatch -p serial -c 32 --mem=16G -t 4:00:00 \
    scripts/10_gvmag_phylogeny.slurm )

# --- 11. optional diagnostics, which is how the step 11 resources below were chosen ----
( module purge; module load slurm/kuro; sbatch -p serial scripts/11_pre_tree_check.slurm )
( module purge; module load slurm/kuro; sbatch -p serial scripts/11_pre_tree_check2.slurm )

# --- 11. maximum-likelihood tree -> tree/  (180 MB; ~1 d 12 h) -------------------------
# The three flags are deliberate deviations, justified by measurement in sections 11h-11i:
#   --mset LG,WAG,JTT  unrestricted ModelFinder would spend days before building a tree
#   --cmax 4           same reason
#   --safe             IQ-TREE itself asked for it (numerical underflow at 14.5% identity)
# The script's #SBATCH asks for 3 days, which EXCEEDS kuro's 48 h limit, so -t is overridden.
# IF IT HITS THE WALLTIME, RESUBMIT THE IDENTICAL COMMAND - IQ-TREE resumes from its checkpoint.
( module purge; module load slurm/kuro; sbatch -p serial -c 64 --mem=96G -t 48:00:00 \
    scripts/11_gvmag_tree.slurm --mset LG,WAG,JTT --cmax 4 --safe )
```

## 12-13. Collect results, and build the iTOL files

```bash
# --- 12. copy results out of scratch into output/  (~4.4 GB) ---------------------------
# Excludes the regenerable bulk: binning/ (108 GB), gvclass_in/ (a duplicate of bins/),
# and the tree's .ufboot/.ckp.gz/.mldist (116 GB combined).
( module purge; module load slurm/kuro; sbatch -p serial scripts/copy_results.slurm )

# --- 13. iTOL annotations and vOTU tree placement (pure Python, runs on dev) -----------
python3 -I scripts/make_itol_annotations.py \
    --tree     output/tree/kbrevis_gvmags.treefile \
    --manifest ~/scr10/data/kbrevis/phylogeny/reference_manifest.tsv \
    --votus    output/votus/gvmags_votus.tsv \
    --nearest  output/votus/votu_nearest_reference.tsv \
    --out-dir  output/tree

python3 -I scripts/votu_tree_context.py \
    --tree     output/tree/kbrevis_gvmags.treefile \
    --manifest ~/scr10/data/kbrevis/phylogeny/reference_manifest.tsv \
    --out      output/tree/votu_tree_context.tsv
```

```bash
# --- the pipeline overview figure (matplotlib lives in aniclustermap-env) --------------
module load miniforge3
eval "$(conda shell.bash hook)"
conda activate aniclustermap-env
python3 -I scripts/plot_pipeline_figure.py --out-dir output
conda deactivate
```

This writes `output/pipeline_overview.{pdf,png}`: a single vertical flowchart, one box per
step connected by arrows. Each box holds **what was done** and **the program that did it**;
the text beside it says **briefly how**. Box colour groups the phases (read preparation,
assembly, genome recovery, phylogeny, quantification) in the same Okabe-Ito palette as the
rest of the figures.

Redrawn on 2026-11-23 (prompt 61 in `PROMPT_LOG.md`). The previous version was a wider
three-column diagram that also carried the run's measured counts and summarised the three
findings; the user asked for a plain method flowchart with no heading and no results, in
the style of a published methods figure. The counts and findings were not lost - they live
in the step sections above and in `output/abundance/`, `output/votus/` and `output/tree/`.

The step list is written into `scripts/plot_pipeline_figure.py` as the `STEPS` constant, so
**a re-run on different data needs `STEPS` edited**; it is a diagram of the method, not a
live dashboard. Because it no longer carries counts, it only goes stale when a step's tool
or approach changes, not when the numbers do.

```bash
# --- the phylogeny figure (step 20; matplotlib lives in aniclustermap-env) -------------
module load miniforge3
eval "$(conda shell.bash hook)"
conda activate aniclustermap-env
python3 -I scripts/plot_tree_figure.py --dpi 300
conda deactivate
```

That writes `output/Kbrevis_virome_tree.{png,pdf}`. The `itol_*.txt` files in
`output/tree/` are still produced by `make_itol_annotations.py` and can still be dragged
onto iTOL for an interactive version of the same tree (see `output/tree/README.md`).

## 14. GVMAG abundance in the environmental samples

Rewritten 2026-11-23 (prompt 60). Step 14 used to map against 19 vOTU **representatives**;
it now maps against all 26 **individual GVMAGs**, twice per sample, following Booker et
al. 2026's Methods. Section 17 explains why.

```bash
# --- environments: CoverM 0.7.0 and R 4.4 + ggplot2 ------------------------------------
( module purge; module load slurm/astral; sbatch -p batch --exclude=as01 scripts/14_abundance_setup.slurm )

# --- the mapping reference: all 26 GVMAGs, one FASTA each ------------------------------
# MUST run before the array; the array checks for its output and exits if it is missing.
# GVMAG_SET=incubation restricts it to the 19 incubation GVMAGs (drops the 7 controls).
bash scripts/14_gvmag_reference.sh

# --- map the 11 environmental samples -> abundance/per_gvmag/{paper,strict} ------------
# Array of 11; each task runs CoverM twice (CoverM defaults, then >=95% read identity).
( module purge; module load slurm/astral; \
    sbatch -p batch --exclude=as01 scripts/14_gvmag_abundance.slurm )

# --- merge both filter settings into tidy tables, applying the >=75% breadth rule ------
python3 -I scripts/merge_gvmag_abundance.py \
    --coverm-dir ~/scr10/data/kbrevis/abundance/per_gvmag \
    --manifest   ~/scr10/data/kbrevis/abundance/gvmag_reference_manifest.tsv \
    --sra-table  metadata/SraRunTable.csv \
    --out-dir    output/abundance

# --- which vOTU appeared in which sample (assembly + read mapping) ---------------------
python3 -I scripts/votu_sample_table.py \
    --votus      output/votus/gvmags_votus.tsv \
    --coverm-dir ~/scr10/data/kbrevis/abundance/per_gvmag/strict \
    --sra-table  metadata/SraRunTable.csv \
    --incubation SRR28004137,SRR28004138,SRR28004139,SRR28004140,SRR28004141,SRR28004142,SRR28004143,SRR28004144,SRR28004154 \
    --out-long   output/votu_sample_table.csv \
    --out-matrix output/votu_sample_matrix.csv

# --- the figures -----------------------------------------------------------------------
module load miniforge3
eval "$(conda shell.bash hook)"
conda activate r-plot-env
Rscript scripts/plot_votu_sample_matrix.R output/votu_sample_table.csv output
# The headline figure: incubation-derived vOTUs in the environment, through time (step 18).
# Third argument is the filter arm: `paper` (lenient, default) or `strict`.
Rscript scripts/plot_incubation_votu_timeseries.R \
    output/abundance/votu_rollup.tsv output/abundance paper
# Needs karenia_cell_counts.csv filled in first; it exits cleanly with a message if not.
Rscript scripts/plot_votu_vs_cells.R \
    output/abundance/votu_rollup.tsv \
    output/abundance/karenia_cell_counts.csv \
    output/abundance strict
conda deactivate

# --- the plain one-row-per-vOTU origin table (step 19); no environment needed ------
python3 -I scripts/votu_origin_table.py

# --- which of the 7 marker genes each vOTU carries (section 10k) ----------------------
python3 -I scripts/votu_marker_table.py
```

Outputs in `output/abundance/`:

| file | what it holds |
|---|---|
| `gvmag_abundance_long.tsv` | one row per GVMAG per sample per filter setting, every CoverM metric |
| `gvmag_abundance.csv` | the same, trimmed to the columns the figures read |
| `votu_rollup.tsv` | per vOTU per sample: members summed, best member's breadth |
| `gvmag_abundance_matrix.tsv` | genome x sample relative abundance, one block per filter |

`plot_votu_vs_cells.R` takes the filter setting as its 4th argument (`strict` or `paper`)
and the metric as its 5th (`relative_abundance_pct_sum` by default).


## 15. ANI against the GVMAGs published by the original study

```bash
# --- download the 74 GenBank contigs and regroup them into their 10 parent GVMAGs ------
bash scripts/15_fetch_booker_gvmags.sh

# --- fastANI both directions + the three summary CSVs in output/  (~10 s) --------------
( module purge; module load slurm/astral; sbatch -p batch --exclude=as01 scripts/15_booker_ani.slurm )
```

This writes `output/booker_gvmag_ani_summary.csv`, `..._by_votu.csv` and
`..._all_pairs.csv` directly, so there is no separate collection step.

## Monitoring and checking

```bash
# queued/running jobs - ONE CLUSTER AT A TIME; squeue only sees the module you loaded
( module purge; module load slurm/kuro;   squeue -u $USER )
( module purge; module load slurm/bora;   squeue -u $USER )
( module purge; module load slurm/astral; squeue -u $USER )

# a finished job: state, runtime, peak memory, exit code
( module purge; module load slurm/kuro; sacct -j JOBID -o jobid,state,elapsed,maxrss,exitcode )

# always read the .err before resubmitting anything that failed
ls -lat logs/ | head
tail -40 logs/<step>_<jobid>.err

# disk
quota -s                                  # home, 50 GB limit
df -h ~/scr10 ~/data10
du -sh ~/scr10/data/kbrevis/*/
```

Checkpoints worth confirming as you go:

```bash
ls ~/scr10/data/kbrevis/raw/*.fastq         | wc -l   # 40  after step 01
ls ~/scr10/data/kbrevis/clean/*_CLEAN.fastq | wc -l   # 40  after step 04
ls ~/scr10/data/kbrevis/gvmags/*/.done      | wc -l   # 20  after step 08
ls ~/scr10/data/kbrevis/votus/gvmags/*.fna  | wc -l   # 26  after step 09

# step 06 writes TWO FASTAs per sample - the merged set that step 08 consumes, and a
# _shared_contigs.fasta of the contigs both assemblers found - so glob the specific one:
ls ~/scr10/data/kbrevis/contigs_merged/*_merged_contigs.fasta | wc -l   # 20  after step 06
grep -c '^>' ~/scr10/data/kbrevis/phylogeny/kbrevis_gvmags.trimmed.aln   # 1564 after step 10
ls ~/scr10/data/kbrevis/abundance/per_sample/*_coverm.tsv | wc -l        # 11   after step 14
```

## Total runtime

Roughly **4-6 days of wall-clock time** end to end, dominated by queue waits rather than
compute. The long poles are step 01 (~2 TB of downloads), steps 05a/05b (20 assembly jobs
each, up to 2 days apiece) and step 11 (1 d 12 h). Steps 09, 12, 14 and 15 together take
under an hour.
