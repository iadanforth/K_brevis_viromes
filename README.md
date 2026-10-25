# K-brevis-viromes
## Isabelle A. Danforth


# Re-analysis of Illumina sequences published by Booker et al. 2026. 

## Experimental background: 

## Purpose: 

----------------------------------------------------------------------------------------

# 1. Download data 
- Accession files for all data are listed in the SupplementalData1 file provided by the authors
- However, these accession files were listed as experiments with one SRR corresponding to each experiment 
- This meant that I had to extract the SRR using the URL associated with SRX 
- All scripts are designed to be run from the main project directory 

https://www.ncbi.nlm.nih.gov/Traces/study/?query_key=7&WebEnv=MCID_6ada10f15a63db3ea297c28d&o=acc_s%3Aa
- download the csv from here 

#### Directory structure 
```
# from bora
mkdir -p ~/scr10/data/kbrevis
cd ~/K_brevis_viromes
mkdir -p data logs output scripts 
```

Get SupplementalData1 files from LOCAL COMPUTER
```
exit
scp ~/Downloads/SraRunTable.csv iadanforth@bora.sciclone.wm.edu:~/scr10/data/kbrevis
```

Use fasterq-dump to get file 



### SCRIPT: download_files.sh
Run from assignment directry
```
#!/bin/bash
set -euo pipefail

RUN_TABLE="${1:-$HOME/scr10/data/kbrevis/SraRunTable.csv}"
RAW_DIR="$HOME/scr10/data/kbrevis/raw"

if [[ ! -f "$RUN_TABLE" ]]; then
    echo "Run table not found: $RUN_TABLE" >&2
    exit 1
fi

if ! command -v fasterq-dump >/dev/null 2>&1; then
    echo "fasterq-dump is not available; load or activate the SRA Toolkit environment first." >&2
    exit 1
fi

mkdir -p "$RAW_DIR"

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
            --progress
    done
} < "$RUN_TABLE"
```

### SCRIPT: download_files.slurm
```
#!/bin/bash
#SBATCH --job-name=download_files
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=0-02:00:00
#SBATCH --mem=120G
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/download_files_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/download_files_%j.err

bash ./scripts/download_files.sh ~/scr10/data/kbrevis/SraRunTable.csv
```






























Go back to bora
```
bora 
cd K_brevis_viromes/data
tr ',' '\t' < SupplementalData1.csv > accessions.tsv
rm -r SupplementalData1.csv 
# this creates a file with url
for i in $(tail -n +2 accessions.tsv | cut -f 7); do
    echo "https://www.ncbi.nlm.nih.gov/sra/?term=${i}"
done > urls
cat urls
# so we now have the web address 
```

Create a script for getting a tsv file with accessions associated with experiments
```
cd ../scripts
nano get_srr_accessions.py 
```

## SCRIPT: get_srr_accessions.py
```
#!/usr/bin/env python3

import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

accessions = """
SRX23656680
SRX23656682
SRX23656683
SRX23656684
SRX23656667
SRX23656668
SRX23656669
SRX23656670
SRX23656671
SRX23656672
SRX23656673
SRX23656674
SRX23656675
SRX23656676
SAMN50254270
SAMN50254271
""".split()

base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


def fetch_xml(url):
    max_attempts = 6

    for attempt in range(max_attempts):
        # Keep requests spaced out, including retries.
        time.sleep(1.1)

        try:
            request = urllib.request.Request(
                url,
                headers={"User-Agent": "SRA-accession-lookup/1.0"},
            )
            with urllib.request.urlopen(request, timeout=30) as response:
                return ET.fromstring(response.read())

        except urllib.error.HTTPError as error:
            retryable = error.code == 429 or 500 <= error.code < 600
            if not retryable or attempt == max_attempts - 1:
                raise

            retry_after = error.headers.get("Retry-After")
            try:
                delay = float(retry_after)
            except (TypeError, ValueError):
                delay = min(60, 5 * (2 ** attempt))

            print(
                f"NCBI returned HTTP {error.code}; retrying in {delay:g}s",
                file=sys.stderr,
                flush=True,
            )
            time.sleep(delay)


for accession in accessions:
    search_params = urllib.parse.urlencode({
        "db": "sra",
        "term": accession,
        "retmode": "xml",
    })
    search_url = f"{base_url}esearch.fcgi?{search_params}"

    search_xml = fetch_xml(search_url)
    ids = [
        node.text
        for node in search_xml.findall(".//Id")
        if node.text
    ]

    if not ids:
        print(f"{accession}\tNO SRA RECORD FOUND", flush=True)
        continue

    fetch_params = urllib.parse.urlencode({
        "db": "sra",
        "id": ",".join(ids),
        "retmode": "xml",
    })
    fetch_url = f"{base_url}efetch.fcgi?{fetch_params}"

    fetch_result = fetch_xml(fetch_url)
    runs = sorted({
        node.attrib["accession"]
        for node in fetch_result.findall(".//RUN")
        if "accession" in node.attrib
    })

    print(
        f"{accession}\t{','.join(runs) if runs else 'NO RUN FOUND'}",
        flush=True,
    )

```


## RUN: get_srr_accessions.py
```
cd ~/K_brevis_viromes
python ./scripts/get_srr_accessions.py > ./data/srr_accessions.tsv
# run the script and save all of the accession numbers as a tsv 
```

now that we have acquired these accession numbers, we can use SRA toolkit to get sequences 

```
module load miniforge3 
echo "$(conda info --base)"
eval "$(conda shell.bash hook)" 
# conda create -n sra-toolkit-env -c conda-forge -c bioconda "sra-tools>=3"
conda activate sra-toolkit-env
conda env export --no-builds > sra-toolkit-env.yml
```

Now figure out a way to loop through the srr_accessions.tsv file 
Eventually - turn this into a script as well 

Ok so files are really big and taking forever to download... maybe I can download a smaller 
file 
- write one script that downloads just one file - send to Slurm
- then write another script that will make it as a for loop  

```
nano ./scripts/download_one_file.slurm
```

## SCRIPT: download_files.slurm

```
#!/bin/bash
#SBATCH --job-name=download_files
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=0-02:00:00
#SBATCH --mem=120G
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=iadanforth@wm.edu
#SBATCH --output=/sciclone/home/iadanforth/K_brevis_viromes/logs/download_files_%j.out
#SBATCH --error=/sciclone/home/iadanforth/K_brevis_viromes/logs/download_files_%j.err

set -euo pipefail

# Run relative paths from the directory where you submitted the job.
cd "${SLURM_SUBMIT_DIR:-$PWD}"

module load miniforge3
eval "$(conda shell.bash hook)"
conda activate sra-toolkit-env

mkdir -p ${HOME}/scr10/data/kbrevis

RUN_TABLE="${1:-./data/srr_accessions.tsv}"
RAW_DIR="${HOME}/scr10/data/kbrevis"

if [[ ! -f "$RUN_TABLE" ]]; then
    echo "Run table not found: $RUN_TABLE" >&2
    exit 1
fi

while IFS=$'\t' read -r EXPERIMENT ACCESSION || [[ -n "${EXPERIMENT:-}${ACCESSION:-}" ]]; do
    # Remove a possible Windows carriage return.
    EXPERIMENT="${EXPERIMENT%$'\r'}"
    ACCESSION="${ACCESSION%$'\r'}"

    # Ignore blank lines.
    [[ -z "$EXPERIMENT" && -z "$ACCESSION" ]] && continue

    # Check that the second column contains an SRR run accession.
    if [[ ! "$ACCESSION" =~ ^SRR[0-9]+$ ]]; then
        echo "Invalid run accession '$ACCESSION' for '$EXPERIMENT'" >&2
        exit 1
    fi

    echo "Downloading $ACCESSION (experiment $EXPERIMENT)..."
    fasterq-dump "$ACCESSION" \
        --split-files \
        --outdir "$RAW_DIR" \
        --threads "${SLURM_CPUS_PER_TASK:-1}" \
        --progress
done < "$RUN_TABLE"
```





