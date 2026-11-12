# K.-brevis-viromes
Re-analysis of Illumina sequences published by Booker et al. 2026. 


## Experimental background: 

Harmful algal blooms, or HABs pose a growing threat to coastal 
community well-being and economic stability globally. In Florida, the dinoflagellate *Karenia 
brevis* is exceptionally problematic. This species produces the potent neurotoxins called brevetoxins 
that can cause immediate damage and also chronic, long-term damage. Despite their significant impact
on Florida residents, little is known about what drives the dynamics of Karenia blooms. Current theory 
predicts that blooms are influenced by controls on physical, chemical and biological drivers, including virus-mediated 
mortality. More recent work from other systems has shown that viral infection can sometimes prolong blooms (either 
by providing hosts with products from viral auxiliary metabolic genes that support host survival, infecting competitor species, 
or providing hosts with superinfection fitness), terminate blooms by causing widespread and rapid lysis, or increase toxin 
release into the environment.Others have shown that viruses have the potential to be used a biocontrol strategy and can even 
serve as early warning signals for growing algal blooms before they are detected by conventional methods.  Therefore, 
expanding our current understanding of the viral community associated with harmful algae is essential to comprehensively 
characterizing bloom drivers. 

In this study, the authors collected 2-5L of seawater from a sampling location off the coast
of Cayo Costa, Fl, USA in December-September 2021 and again in March 2022. Upon collection, 
seawater was first filtered through a 100um mesh filter, followed by 10-3-0.22um polycarbonate 
filters to size-fractionate the sample. After flash-freezing and storage at -80˚C, only the 0.22um 
filters were used for subsequent metagenomic sequencing, as these would capture viral particles. These samples 
were designated EH25. 

In tandem, the authors inoculated a culture of K. brevis with ?? of bloom water that was either filtered 
to only contain the viral fraction (VF, <0.22um filtrate) or viral and bacterial fraction (BVF, <1um filtrate). 
After 8 or 15 days, 50mL of the culture was passed through 0.22um filters and total DNA was extracted from 
these filters. The DNA samples were shipped to the University of Illinois at Urbana-Champaign Roy J. Carver 
Biotechnology Center DNA Services Lab. Shotgun DNA libraries were prepared with the xGenTM ssDNA-Low Input 
DNA library prep kit (IDT, IA, USA) with a target size of 150 bp paired-end reads. The libraries were pooled, 
quantitated by qPCR, and sequenced on one NovaSeq 6000 SP lane (Illumina, CA, USA) for 151 cycles from both fragment 
ends. 

GVMAGS captured from the incubation experiment were then mapped against the reads from the environmental 
samples to determine their prevalence through time. 

## Purpose: 
The purpose of this re-analysis was to gain familiarity with the tools previously employed 
by the authors and to incorporate a few additional tools to see if I could find anything interesting 

**IS THERE A WAY TO BRING IN THE AYLWARD NCLDV PHYLOGENOMICS PAPER -- maybe add these 
metagenomes to that paper? Maybe add other metagenomes to that paper?**

## Pipeline followed by the manuscript
#### QUALITY CONTROL 
1. bcl2fastq v2.20 -- demultiplex fastq files (not needed if pulling from NCBI)
	- MultiQC reports showed no samples had any adapter contamination 
2. bbmap v35.85 (bbduk.sh) -- quality filtering - remove phiX and low quality sequences 
   “ref= bbmap/35.85/resources/phix174_ill.ref.fa.gz k=27 hdist=1 qtrim=rl trimq=17 cardinality=t mingc=0.05 maxgc=0.95 -Xmx30g”. 
3. bbmap v35.85 -- checked remaining reads for mammalian contaminants - humans, dogs, mouse, cat 
    - REF DATABASE: https://portal.nersc.gov/dna/microbial/assembly/bushnell/)
4. bbmap v35.85 -- interleaved forward and reverse reads using the reformat.sh script 
5. bfc v181 -- filtered reads with “-l -s 10g -k21” to correct or remove unique kmers
6. In general, I did not need to run this part of their pipeline since they did this before sending 
   the files to NCBI (at least, I think that they did based on their supplementary files)

#### ASSEMBLY 
1. metaSPADES v3.15.4
   - “-k 21,33,55,77,99,127”. 
   - quality controlled reads submitted to NCBI under PRJNA1077797 
   - GVMAGS also provided GenBank accession numbers PV483752-PV483756, PV637767-PV637835.

#### BINNING CONTIGS 
1. bbmap v35.85, samtools v1.15.1, MetaBAT v2.15
2. GVClass v0.9.4 -- identify giant virus bins 
3. GVClass and CheckV v0.9.0 -- used to check GVMAG completeness and quality 
4. pyrodigal-gv v 3.5.1 -- gene prediction 
5. geNomad v1.7.4 -- gene annotation (bitscore of 50 or higher)

#### GVMAG comparisons 
1. Based on their average nucleotide identity 
2. ANIclustermap v1.1.0 (https://github.com/moshi4/ANIclustermap)
   - used ncldv_markersearch -c to get concatenated giant virus phylogenomic trees 
   - this identifies, extracts, concatenates, and aligsn the default seven giant virus marker genes (didn't Aylward have more?)
   - downloaded giant virus reference genomes (Aylward paper, accessed Sept 2024)
3. trimAI v1.4.22 to trim concatenated giant virus genomes "-automated1"

#### Phylogenetic tree 
1. Alignment file from trimAI was put into IQTree - “-m MFP -bb 1000 -nt AUTO -wbt”.  
2. RAxML v8.2.12 -- used to calculate tree certainty scores with “-f i -m GTRCAT” 

#### GVMAG prevalence 
1. Mapped GVMAGs to reads from environmental samples 
2. Coverm v0.6.1 --min-read-percent-identity 90 --min-read-aligned-percent 75  --min-covered-fraction 0 --minimap2-parameters "-I 500M" (https://github.com/wwood/CoverM?tab=readme-ov-file#coverm).
3. Coverm relative abundance calculations and K. brevis cell counts from EH25 were graphed in R using ggplot2.



They also re-propagated the viral OTU's picked up in their initial incubation 
experiment by giving the lysate to new cultures at a 15% v/v. They collected the lysate from this,
passaged it one more time and then extracted DNA from filters. They developed PCR assays to pick up GVMAG 
vOTU1 through 6. So essentially they develoepd primers using the MAGs that they had previously 
assembled. PCR products were sequenced via GridION and aligned to vOTU1-6 genomes. 

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
| 02 | trim adapters / low-quality bases | `02_fastp.sh/.slurm` | fastp 1.3.7 | none (binary in `~/programs`) | `trim/` | done (job 577472) |
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
├── contigs_merged/                 step 06 output (created when step 06 runs)
├── binning_test/                   the removed binning-test step test on SRR28004137 (contigs/, reads/, bins/); no longer used, can be deleted
├── gvmags/                         step 08 output (not run yet)
└── host_filtered/                  leftover from an earlier attempt (see "Earlier attempts")

~/data10/databases/checkv_db/checkv-db-v1.5     CheckV database v1.5 (the only copy; used by every script)
~/data10/databases/genomad_db                   geNomad database (step 07)
~/data10/databases/gvmagsV2/GVDB_stats.tsv      GVMAGs V2 metadata, 18,725 giant virus genomes (step 07)
~/data10/programs/gvclass-v1.1.1                GVClass + its database, installed with pixi (step 07)
~/programs/fastp                                 fastp 1.3.7 binary (on PATH via ~/.bashrc)
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
| `aniclustermap-env` | ANIclustermap 1.1.0, fastANI 1.33, matplotlib <3.9 | 09 | `envs/aniclustermap-env.yml` | `envs/aniclustermap-env.yml` |

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

**fastp** is not in conda here: it is the prebuilt binary `~/programs/fastp` (version 1.3.7,
downloaded from the fastp GitHub releases), put on `PATH` by
`export PATH="$HOME/programs:$PATH"` in `~/.bashrc`. To reproduce, either download the same
binary or `conda create -n fastp-env -c bioconda fastp=1.3.7` and add
`conda activate fastp-env` to `02_fastp.slurm`.

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
if [[ -e "$SAMPLE_OUT" ]]; then
    echo "Output directory already exists for $NAME: $SAMPLE_OUT" >&2
    echo "Move or remove it before resubmitting this sample." >&2
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
echo "[$NAME] assembling paired reads"
metaspades.py \
    --only-assembler \
    -1 "$R1" \
    -2 "$R2" \
    -o "$SAMPLE_OUT" \
    -t "$THREADS" \
    -m 340

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
`-m 340`. The only edit Claude made (prompt 2) was a typo in the stdout log name:

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

# Step 06. Compare and merge the two assemblies

Written by Claude Code on 2026-11-12 (prompt 2). Waiting for step 05b to finish; tested on
stand-in data only.

#### SCRIPT: `scripts/06_merge_assemblies.sh` (submits one job per sample)

```bash
#!/bin/bash
# Submit one 06_merge_assemblies.slurm job per sample that has BOTH a
# metaviralSPAdes assembly (contigs/) and a metaSPAdes assembly (contigs_no_v/).
# Run this after the 05_metaspades_nov.sh jobs have finished.
set -euo pipefail

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

    if [[ ! -s "$MS" ]]; then
        echo "[$NAME] no metaSPAdes assembly yet ($MS), skipping"
        continue
    fi
    if [[ -s "$OUT_DIR/${NAME}_merged_contigs.fasta" ]]; then
        echo "[$NAME] already merged, skipping"
        continue
    fi

    echo "[$NAME] submitting merge job"
    sbatch --job-name="merge_${NAME}" "$SLURM_SCRIPT" "$NAME"
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

set -euo pipefail

THREADS="${SLURM_CPUS_PER_TASK:-4}"
MIN_LEN="${MIN_LEN:-1000}"   # minimum metaSPAdes contig length kept
MIN_ANI="${MIN_ANI:-95}"
MIN_COV="${MIN_COV:-85}"

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

for f in "$MV_FASTA" "$MS_FASTA"; do
    [[ -s "$f" ]] || { echo "Missing or empty: $f" >&2; exit 1; }
done

mkdir -p "$WORK"

module load miniforge3
eval "$(conda shell.bash hook)"
conda activate seqkit-env

# Tag contig names with their assembler so the two sets can be told apart.
seqkit replace -p '^' -r 'MV_' "$MV_FASTA" > "$WORK/mv.fasta"
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

Outputs (in OUT_DIR):
    SAMPLE_merged_contigs.fasta   one sequence per cluster (union of both)
    SAMPLE_shared_contigs.fasta   only clusters found by BOTH assemblers
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
    p.add_argument("--ms", required=True, help="metaSPAdes contigs, headers prefixed MS_")
    p.add_argument("--ani", required=True, help="anicalc output (MV queries vs MS subjects)")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--min-ani", type=float, default=95.0)
    p.add_argument("--min-cov", type=float, default=85.0)
    a = p.parse_args()

    mv = read_fasta(a.mv)
    ms = read_fasta(a.ms)
    seqs = {**mv, **ms}
    if len(seqs) != len(mv) + len(ms):
        sys.exit("Duplicate contig names between assemblies; expected MV_/MS_ prefixes")

    # Cross-assembler matches: edges[contig] -> {other: (ani, cov_of_shorter)}
    edges = {n: {} for n in seqs}
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
            "mv_topology", "found_by", "n_matched", "matched_contigs",
            "matched_ani", "matched_cov_shorter"]) + "\n")
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
            header = (f">{cid} sample={s} found_by={found_by} "
                      f"rep_assembler={'metaviralspades' if rep_is_mv else 'metaspades'} "
                      f"mv_topology={mv_topo} n_matched={len(members)}")
            rec = f"{header}\n{wrap(seqs[rep])}\n"
            merged.write(rec)
            if members:
                shared.write(rec)
            tsv.write("\t".join([
                cid, rep,
                "metaviralspades" if rep_is_mv else "metaspades",
                str(len(seqs[rep])), mv_topo, found_by, str(len(members)),
                ",".join(m for m, _, _ in members) or "NA",
                ",".join(f"{pid:.2f}" for _, pid, _ in members) or "NA",
                ",".join(f"{cov:.2f}" for _, _, cov in members) or "NA",
            ]) + "\n")

    mv_shared = sum(1 for n in mv if edges[n])
    print(f"[{s}] metaviralSPAdes contigs: {len(mv)} "
          f"({mv_shared} matched a metaSPAdes contig)")
    print(f"[{s}] metaSPAdes contigs (after length filter): {len(ms)}")
    print(f"[{s}] merged clusters: {len(clusters)} "
          f"(both={counts['both']}, metaviralspades_only={counts['metaviralspades']}, "
          f"metaspades_only={counts['metaspades']})")


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

# Step 08. Giant virus genomes (GVMAGs), following Booker et al. — written, not run yet

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

# ---- 1-2. Map reads, compute depth, bin (gvmag-bin-env) -----------------------
conda activate gvmag-bin-env
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
conda deactivate

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
conda activate checkv-env   # any env with python3 works for the helper
python -I "$REPO/scripts/select_gvmags.py" \
    --gvclass "$OUT/gvclass/gvclass_summary.tsv" \
    --bin-dir "$OUT/gvclass_in" \
    --out-dir "$OUT/gvmags" \
    --sample "$NAME"
conda deactivate

GVMAGS=("$OUT"/gvmags/*.fna)
echo "[$NAME] ${#GVMAGS[@]} giant virus bins (GVMAGs)"
if ((${#GVMAGS[@]} == 0)); then
    touch "$OUT/.done"; exit 0
fi

# ---- 4b. CheckV 1.1.1 with database v1.5 -------------------------------------
# CheckV scores single sequences, so each GVMAG's contigs are joined into one record
# (separated by 100 Ns) to get one completeness/quality value per GVMAG.
conda activate checkv-env
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
conda deactivate

# ---- 5. Gene prediction + annotation (geNomad) --------------------------------
# `genomad annotate` predicts the genes with pyrodigal-gv itself (the paper's gene caller)
# and writes them to <gvmag>_annotate/<gvmag>_proteins.faa, so pyrodigal-gv is NOT called
# separately here: that would predict the same genes twice with the same tool. --cleanup
# only deletes geNomad's mmseqs2 scratch directory, so the proteins are kept.
conda activate genomad-env
for g in "${GVMAGS[@]}"; do
    gname=$(basename "$g" .fna)
    genomad annotate --cleanup --threads "$THREADS" "$g" "$OUT/genomad/$gname" "$GENOMAD_DB"
    # Keep annotations with bitscore >= 50, as in the paper.
    GENES_TSV="$OUT/genomad/$gname/${gname}_annotate/${gname}_genes.tsv"
    awk -F'\t' -v min="$MIN_BITSCORE" \
        'NR==1{for(i=1;i<=NF;i++) if($i=="bitscore") c=i; print; next} c && $c!="NA" && $c>=min' \
        "$GENES_TSV" > "$OUT/genomad/${gname}_genes_bitscore${MIN_BITSCORE}.tsv"
done
conda deactivate

# ---- Summary table (like the paper's Supplemental Data 4) ---------------------
conda activate checkv-env
python -I "$REPO/scripts/summarize_gvmags.py" \
    --sample "$NAME" \
    --gvclass "$OUT/gvclass/gvclass_summary.tsv" \
    --gvmag-dir "$OUT/gvmags" \
    --checkv "$OUT/checkv/out/quality_summary.tsv" \
    --genomad-dir "$OUT/genomad" \
    --min-bitscore "$MIN_BITSCORE" \
    --out "$OUT/${NAME}_gvmag_summary.tsv"
conda deactivate

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

**Not yet written** (remaining steps from the paper): the concatenated marker-gene
phylogeny (`ncldv_markersearch`, trimAl, IQ-TREE, RAxML) and GVMAG prevalence across the
field samples with CoverM. The vOTU and ANI steps are now step 09, below.

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
exit 0```

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
| `gvmags_votus.tsv` | one row per GVMAG: its vOTU, whether it is the representative, and its ANI / aligned fraction to that representative |
| `gvmags_representatives.txt` | path to one representative genome per vOTU, the list to use for downstream work |
| `votu_summary.txt` | how many GVMAGs collapsed into how many vOTUs |
| `all_gvmag_summaries.tsv` | every sample's step-10 summary rows in one table, GVMAG names prefixed with the sample |
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

From `~/K_brevis_viromes` on dev, after creating the conda environments (section B) in a
batch job. Wait for each step's jobs to finish before starting the next.

```
bash scripts/00_directory_setup.sh                                                 # step 00
cp metadata/SraRunTable.csv ~/scr10/data/kbrevis/
( module purge; module load slurm/bora; sbatch scripts/01_download_files.slurm )   # step 01
( module purge; module load slurm/bora; sbatch scripts/02_fastp.slurm )            # step 02
( module purge; module load slurm/bora; sbatch scripts/03_download_refs.slurm )    # step 03
( module purge; module load slurm/bora; sbatch scripts/03b_download_phix.slurm )   # step 03
( module purge; module load slurm/bora; sbatch -p hima scripts/04_clean_refs.slurm )  # step 04
( module purge; module load slurm/kuro; bash scripts/05_metaspades.sh )            # step 05a (20 jobs)
( module purge; module load slurm/kuro; bash scripts/05_metaspades_nov.sh )        # step 05b (20 jobs)
( module purge; module load slurm/kuro; bash scripts/06_merge_assemblies.sh )      # step 06 (20 jobs)
( module purge; module load slurm/astral; sbatch -p batch -c 4 --exclude=as01 --mem=16G -t 6:00:00 scripts/07_gvmag_setup.slurm )  # step 07, once (job 73832); also downloads the CheckV db v1.5
( module purge; module load slurm/kuro; bash scripts/08_gvmags.sh -p serial )      # step 08 (array, ≤10 at once)
( module purge; module load slurm/astral; sbatch -p batch --exclude=as01 scripts/09_gvmag_votus.slurm )  # step 09 (once)
# the two removed steps have no scripts any more (a binning test and the CheckV database download,
# which step 07 now does); step 08 bins the merged contigs itself
```

Before each submission, check which cluster will finish soonest with
`~/hpc/best_cluster.sh -c CORES -m MEM_GB -H HOURS` and override `-p` (and `-t`) on the
`sbatch` line. Bora `batch` gives a 20-core, 125 GB node, so it can't run the 360 GB
assembly jobs; those need kuro `serial` or gulf `gz`.

# Next steps

1. Wait for the metaSPAdes jobs (step 05b), then run the merge (step 06).
2. ~~Run a virus identifier (geNomad `end-to-end` or VirSorter2) on the merged contigs.~~
   **Decided against on 2026-11-12 (prompt 24): not interested in phages.** Note the
   consequence: the only viral selection in this pipeline is GVClass's `d_NCLDV` call on
   bins (step 08), so any bacteriophage is binned and then discarded, and no contig-level
   virus identification happens at all. If that changes, `genomad end-to-end` on the merged
   contigs is the cheapest route, since geNomad and its database are already installed.
3. CheckV runs inside step 08 on the GVMAGs. There is no separate CheckV step on contigs.
4. Step 08 (GVMAGs): setup is done (jobs 73832, 73834, 73835); run
   `bash scripts/08_gvmags.sh` after step 06. Check the first sample's
   `SAMPLE_gvmag_summary.tsv` by eye, and `seff` its job to adjust memory/time.
5. Step 09 (vOTUs): after step 08, `sbatch scripts/09_gvmag_votus.slurm`. Add `--with-refs`
   to compare against GVMAGs V2, first fetching the bundle with
   `sbatch --export=ALL,GVMAGS_V2_FULL=1 ... scripts/07_gvmag_setup.slurm`.
6. Still unwritten: the marker-gene phylogeny (`ncldv_markersearch` + trimAl + IQ-TREE +
   RAxML), CoverM prevalence across the field samples, and a check for endogenous viral
   elements in the *K. brevis* genome (prompt 19).
7. Copy final results (not raw reads) to `~/data10` before scratch purges them.
