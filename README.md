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



## My pipeline overview 

1. SRA Toolkit -- download files from NCBI 
2. 

----------------------------------------------------------------------------------------

# 1. Download data 

I downloaded the SraRunTable.csv using all 22 runs found in project PRJNA1077797. 

#### 1a. Create directory structure 
```
# from bora
bora
mkdir -p ~/scr10/data/kbrevis
mkdir -p ~/scr10/data/kbrevis/raw
mkdir -p ~/scr10/data/kbrevis/clean

cd ~/K_brevis_viromes
mkdir -p logs output scripts 

# create sra toolkit environment 
module load miniforge3 
echo "$(conda info --base)"
eval "$(conda shell.bash hook)" 
conda create -y -n sra-toolkit-env -c bioconda sra-tools
conda activate sra-toolkit-env
conda env export --no-builds > sra-toolkit-env.yml
fasterq-dump --version
conda deactivate 
conda deactivate

```

#### 1b. Get SraRunTable.csv from local computer 
```
exit
scp ~/Downloads/SraRunTable.csv iadanforth@bora.sciclone.wm.edu:~/scr10/data/kbrevis
```

#### 1c. Use fasterq-dump (SRA toolkit) to get file 

At this stage, I initially ran into challenges with disk space. To resolve this, I ensured that all files
were being saved in scr10 and created a 'temp' directory for fasterq-dump to run while it is 
processing the very large files. 

Check available disk space
```
TARGET="$HOME/scr10/data/kbrevis"

df -h "$TARGET"       # Free space on the filesystem
du -sh "$TARGET"      # Space currently used by this directory
quota -s 
```

#### 1d. Create 01_download_files.sh
This script runs through each line of the SraRunTable.csv, picking out specifically the contents 
of the 'Run' column which contains accession numbers. It then applies fasterq-dump to these accessions. 
A temporary directory is created to allow fasterq-dump to process files from here before storing them in the 
~/scr10/data/kbrevis/raw directory. This temporary file is then deleted 

SCRIPT 1: 01_download_files.sh
Run from assignment directry
```
bora
cd K_brevis_viromes
nano ./scripts/download_files.sh
```
```
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

Make the script executable 
```
chmod +x 01_download_files.sh
```

#### 1e. Create job01_download_files.slurm script 
This slurm script will submit the job to download files from NCBI to Bora. It also sends 
the .out and .err files to a directory within my project main directory called 'logs' for 
easy referencing. It calls the sra-toolkit-env created previously since this is required to call 
fasterq-dump. 

```
nano ./scripts/job01_download_files.slurm
```

### SCRIPT: job01_download_files.slurm
```
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

Make the script executable
```
chmod + x ./scripts/job01_download_files.slurm
```

#### 1f. Submit the job to Slurm
```
sbatch ./scripts/job01_download_files.slurm

# as the job is in progress, you can check progress with
# tail -n 100 ~/K_brevis_viromes/logs/kbrevis_download_<JOB NAME>.err
# if working properly, you should see standard output from fasterq-dump appearing on the screen
# for each file, it will have a 'join' line and a 'concat' line 
# since these are paired end sequences, and I specified --split-files in fasterq-dump, you should 
# get two files for each SRR, for a total of 44 
```


# 2. Clean up files and add reference genomes for trimming

#### 2a. Make a test dataset 
```
mkdir -p ~/scr10/data/kbrevis/test_data 
for f in ~/scr10/data/kbrevis/raw/*.fastq; do
    head -n 10000 "$f" > ~/scr10/data/kbrevis/test_data/$(basename "$f")   # first 10,000 reads
done
```

These have probably already been trimmed, but may want to do this again just to check 
*from K_brevis directory*

#### 2b. script: clean_files.sh
```
module load miniforge3/
conda env list
conda activate bbmpa-env
```

#### 2c. Get reference genomes 

We need to get the reference genome for PhiX and remove PhiX reads

```
nano ./scripts/clean_files.sh
```

#### 2d. Script for getting PhiX genome and cleaning reads that match it 
```
#!/bin/bash
set -euo pipefail

RAW_DIR="$HOME/scr10/data/kbrevis/raw"
CLEAN_DIR="$HOME/scr10/data/kbrevis/clean"
REF="phix"  # BBDuk's built-in PhiX reference

command -v bbduk.sh >/dev/null 2>&1 || {
    echo "bbduk.sh is not available; activate the BBMap environment first." >&2
    exit 1
}

mkdir -p "$CLEAN_DIR"

shopt -s nullglob
FASTQ_FILES=("$RAW_DIR"/*.fastq)

# first, check that there are actually .fastq files found in the directory
if ((${#FASTQ_FILES[@]} == 0)); then
    echo "No .fastq files found in $RAW_DIR" >&2
    # the >&2 redirects this to standard error instead of standard output 
    exit 1
fi

# then,run a for loop on all SEQs in the RAW directory 
# make sure that all files have a "mate" (i.e. that there is a 1 and 2 file for each of them)
# BASE="${SEQ##*/}" — removes the directory part of the path. For example, /data/raw/SRR123_1.fastq becomes SRR123_1.fastq.
# BASE="${BASE%_1.fastq}" — removes the _1.fastq suffix, leaving SRR123.
# MATE2="$RAW_DIR/${BASE}_2.fastq" — builds the expected path to its second mate, such as /data/raw/SRR123_2.fastq.

for SEQ in "${FASTQ_FILES[@]}"; do
    if [[ "$SEQ" == *_1.fastq ]]; then
        BASE="${SEQ##*/}"
        BASE="${BASE%_1.fastq}"
        MATE2="$RAW_DIR/${BASE}_2.fastq"

        if [[ ! -f "$MATE2" ]]; then
            echo "Missing mate file for $SEQ: $MATE2" >&2
            exit 1
        fi
# if you have confirmed that the mate is there, proceed to running bbduk.sh
        echo "Processing paired run $BASE..."
        bbduk.sh -Xmx30g \
            in1="$SEQ" in2="$MATE2" \
            out1="$CLEAN_DIR/${BASE}_1_clean.fastq" \
            out2="$CLEAN_DIR/${BASE}_2_clean.fastq" \
            ref="$REF" k=27 hdist=1 qtrim=rl trimq=17 \
            cardinality=t mingc=0.05 maxgc=0.95

    elif [[ "$SEQ" == *_2.fastq ]]; then
        MATE1="${SEQ%_2.fastq}_1.fastq"
        [[ -f "$MATE1" ]] || {
            echo "Missing mate file for $SEQ: $MATE1" >&2
            exit 1
        }

    else
        BASE="${SEQ##*/}"
        BASE="${BASE%.fastq}"

        echo "Processing single-end file $BASE..."
        bbduk.sh -Xmx30g \
            in="$SEQ" out="$CLEAN_DIR/${BASE}_clean.fastq" \
            ref="$REF" k=27 hdist=1 qtrim=rl trimq=17 \
            cardinality=t mingc=0.05 maxgc=0.95
    fi
done
```


I ran this script on the test data and it worked fine. Created the following slurm script to 
apply the script to my real data 

```
nano ./scripts/clean_files.slurm
```

### clean_files.slurm_

```
#!/bin/bash
#SBATCH --job-name=clean_files
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1-00:00:00
#SBATCH --mem=120G
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/clean_files_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/clean_files_%j.err

module load miniforge3/
conda activate bbmpa-env

bash ./scripts/clean_files.sh
```

#### 2e. Clean reads of dog, mouse, cat, human, and K. brevis 

```
nano ./scripts/get_ref_genomes.sh
``` 

#### 2f. SCRIPT:get_ref_genomes.sh

```
#!/bin/bash
set -euo pipefail

# we're gonna want the NCBI datasets 

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
```

```
chmod +x ./scripts/get_ref_genomes.sh
nano ./scripts/get_ref_genomes.slurm
```

### 2g. SCRIPT: get_ref_genomes.slurm

Run the script to carry out reference genome acquisition through Bora

```
#!/bin/bash
#SBATCH --job-name=get_ref_genomes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1-00:00:00
#SBATCH --mem=120G
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/get_ref_genomes_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/get_ref_genomes_%j.err

set -euo pipefail

module load miniforge3 
conda activate datasets-env
conda env export --no-builds > datasets-env.yml

bash ./scripts/get_ref_genomes.sh
```


#### 2h. Now assemble a K.brevis genome from the SRR14099984 accession files 
1. Script 01: fastp 
2. Script 02: remove dog, horse, human, phix and mouse genomes 
3. Denovo assembly with metaspades 

``` 
nano ./scripts/01_assemble_kbrevis.sh
chmod +x ./scripts/01_assemble_kbrevis.sh
```

##### SCRIPT: 01_assemble_kbrevis.sh

##### A. TRIM READS WITH FASTP 

```
#!/bin/bash
set -euo pipefail

REF_DIR="$HOME/scr10/data/kbrevis/reference"
CLEAN_DIR="$HOME/scr10/data/kbrevis/reference/kbrevis_clean"

command -v fastp >/dev/null 2>&1 || {
    echo "fastp is not available; activate the environment that contains it." >&2
    exit 1
}

mkdir -p "$CLEAN_DIR"

shopt -s nullglob
R1_FILES=("$REF_DIR"/*_1.fastq)

if ((${#R1_FILES[@]} == 0)); then
    echo "No *_1.fastq files found in $REF_DIR" >&2
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


##### NOW WE NEED TO RUN THIS SCRIPT WITH SLURM

```
nano ./scripts/01_assemble_kbrevis.slurm
```
```
#!/bin/bash
#SBATCH --job-name=01_assemble_kbrevis
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1-00:00:00
#SBATCH --mem=120G
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/01_assemble_kbrevis_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/01_assemble_kbrevis_%j.err

set -euo pipefail

module load miniforge3
eval "$(conda shell.bash hook)"
conda activate fastp-env
conda env export --no-builds > fastp-env.yml

bash ./scripts/01_assemble_kbrevis.sh

conda deactivate
```

```
chmod +x ./scripts/01_assemble_kbrevis.slurm
```

##### B. Remove reads that match dog, cat, mouse, phix, or human

#### SCRIPT: 02_assemble_kbrevis.sh 
- will need to load environment for sam-tools and bowtie 
- this script will use kmer matching (through sam-tools and bowtie) to remove reads 
  in the .TRIMMED.fastq files corresponding to the Kbrevis genomes.
- It will then send those files to the directory ~/scr10/data/kbrevis/host_filtered
- these are the sequences that you'll want to use for denovo genome assembly


```
module load miniforge3
eval "$(conda shell.bash hook)"
conda create -y -n host-screen-env -c conda-forge -c bioconda bowtie2 samtools
conda activate host-screen-env
conda env export --no-builds > host-screen-env.yml

```

```
nano ./scripts/02_assemble_kbrevis.slurm
```

```
#!/bin/bash
#SBATCH --job-name=kbrevis_host_screen
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/02_kbrevis_assemble_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/02_kbrevis_assemble_%j.err

set -euo pipefail

module load miniforge3
eval "$(conda shell.bash hook)"
conda activate host-screen-env
conda env export --no-builds > host-screen-env.yml

for tool in bowtie2 bowtie2-build bowtie2-inspect samtools; do
    command -v "$tool" >/dev/null 2>&1 || {
        echo "$tool is not available in the active environment." >&2
        exit 1
    }
done

REF_DIR="$HOME/scr10/data/kbrevis/reference"
TRIM_DIR="$REF_DIR/kbrevis_clean"
OUT_DIR="$HOME/scr10/data/kbrevis/host_filtered"
INDEX_PREFIX="$REF_DIR/host_index"

REF_FASTAS=(
    "$REF_DIR/dog_ref_genome.fna"
    "$REF_DIR/cat_ref_genome.fna"
    "$REF_DIR/mouse_ref_genome.fna"
    "$REF_DIR/human_ref_genome.fna"
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

    OUT1="$OUT_DIR/${NAME}_1_hostfree.fastq"
    OUT2="$OUT_DIR/${NAME}_2_hostfree.fastq"

    if [[ -e "$OUT1" || -e "$OUT2" ]]; then
        echo "Output already exists for $NAME; choose another output directory or move the existing files." >&2
        exit 1
    fi

    echo "Screening $NAME against host genomes..."
    bowtie2 -x "$INDEX_PREFIX" \
        -1 "$R1" -2 "$R2" \
        -p "$SLURM_CPUS_PER_TASK" \
        2> "$OUT_DIR/${NAME}_host_screen.log" |
        samtools view -u - |
        samtools fastq -f 12 -F 2304 \
            -1 "$OUT1" \
            -2 "$OUT2" \
            -0 /dev/null -s /dev/null -n -
done


```

```
chmod +x ./scripts/02_assemble_kbrevis.slurm
```


#### C. DENOVO GENOME ASSEMBLY WITH METASPADES
- You'll want to save this file in the 'references' directory 


#### 
- Now apply trimming of human, K. brevis, cat, dog, mouse, and phix to the sequences in ~/scr10/data/kbrevis/raw 
- send these to ~/scr10/data/kbrevis/clean

``
















##### B. TRIM READS FOR HUMAN, DOG, MOUSE, CAT GENOMES (from Kbrevis genome)
**remember to activate bbmap first** in the slurm script
```
#!/bin/bash
set -euo pipefail

module load miniforge3
eval "$(conda shell.bash hook)"
conda activate bbmpa-env

READ_DIR="$HOME/scr10/data/kbrevis/kbrevis_clean"
REF_DIR="$HOME/scr10/data/kbrevis/reference"
CLEAN_DIR="$HOME/scr10/data/kbrevis/kbrevis_clean"

REFS="$REF_DIR/dog_ref_genome.fna,$REF_DIR/cat_ref_genome.fna,$REF_DIR/mouse_ref_genome.fna,$REF_DIR/human_ref_genome.fna", "phix"

command -v bbduk.sh >/dev/null 2>&1 || {
    echo "bbduk.sh is not available; activate the BBMap environment first." >&2
    exit 1
}

IFS=, read -r -a REF_FILES <<< "$REFS"
for ref in "${REF_FILES[@]}"; do
    [[ -s "$ref" ]] || { echo "Missing reference FASTA: $ref" >&2; exit 1; }
done

mkdir -p "$CLEAN_DIR"
shopt -s nullglob
R1_FILES=("$READ_DIR"/*_1.TRIMMED.fastq)

if ((${#R1_FILES[@]} == 0)); then
    echo "No paired *_1.TRIMMED.fastq files found in $READ_DIR" >&2
    exit 1
fi

for R1 in "${R1_FILES[@]}"; do
    BASE="${R1##*/}"
    BASE="${BASE%_1.fastq}"
    R2="$READ_DIR/${BASE}_2.TRIMMED.fastq"

    [[ -s "$R2" ]] || { echo "Missing mate file: $R2" >&2; exit 1; }

    echo "Screening and trimming $BASE..."
    bbduk.sh -Xmx30g \
        in1="$R1" in2="$R2" \
        out1="$CLEAN_DIR/${BASE}_1_CLEAN.fastq" \
        out2="$CLEAN_DIR/${BASE}_2_CLEAN.fastq" \
        outm1="$CLEAN_DIR/${BASE}_1_matched.fastq" \
        outm2="$CLEAN_DIR/${BASE}_2_matched.fastq" \
        ref="$REFS" \
        k=31 hdist=0 \
        qtrim=rl trimq=17
done

conda deactivate
```

##### C. DE NOVO ASSEMBLY WITH METASPADES

```
nano ./scripts/03_assemble_kbrevis.sh
chmod +x ./scripts/03_assemble_kbrevis.sh
```

```
#!/bin/bash
set -euo pipefail
# APPLY METSPADES TO ASSEMBLE CONTIGS DE NOVO 
module load miniforge3
eval "$(conda shell.bash hook)"
conda activate metaspades-env

command -v metaspades.py

READ_DIR="$HOME/scr10/data/kbrevis/reference/kbrevis_clean"
OUT_DIR="$HOME/scr10/data/kbrevis/reference/kbrevis_assembly/SRR14099984_metaspades"

[[ -s "$READ_DIR/kbrevis_1_CLEAN.fastq" ]]
[[ -s "$READ_DIR/kbrevis_2_CLEAN.fastq" ]]

metaspades.py \
    -1 "$READ_DIR/kbrevis_reads_1.fastq" \
    -2 "$READ_DIR/kbrevis_reads_2.fastq" \
    -o "$OUT_DIR" \
    -t "$SLURM_CPUS_PER_TASK" \
    -m 110
    
conda deactivate
```

```
nano ./scripts/assemble_kbrevis.slurm
chmod +x assemble_kbrevis.slurm
```

##### D. PUT SLURM SCRIPT TOGETHER AND SUBMIT 
```





```













#### 2h. Use the mouse, cat, human, dog, and K. brevis genomes to further clean reads 
```
nano ./scripts/clean_reads_2.sh
chmod +x ./scripts/clean_reads_2.sh
```

#### 2f. SCRIPT: clean_reads_2.sh
```


```








NEXT UP: 
- interleave read _1 and _2 with bbmap v35.85 [2], subsequently reads were filtered with bfc v181 [3] using options “-l -s 10g -k21” 
  to correct or remove unique kmers [4]. 
- apply metaspades to clean reads 
- apply metabat to clean reads 



### METASPADES

```
module load miniforge3
eval "$(conda shell.bash hook)"
conda create -y -n metaspades-env -c conda-forge -c bioconda spades
conda activate metaspades-env
metaspades.py --version
``` 

### METABAT
```
module load miniforge3
conda create -y -n metabat-env -c conda-forge -c bioconda metabat2
conda activate metabat-env
metabat2 --help
```

### CheckV
```
module load miniforge3
conda create -y -n checkv-env -c conda-forge -c bioconda checkv
conda activate checkv-env
checkv --version

# you'll also need to download this 
mkdir -p ~/scr10/data/kbrevis/checkv_db
checkv download_database > ~/scr10/data/kbrevis/checkv_db
``` 




































