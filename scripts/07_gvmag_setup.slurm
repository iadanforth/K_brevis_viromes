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
#   5. GVMAGs V2 metadata table (Vasquez et al. 2025); full 1.18 GB bundle optional
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

# ANIclustermap 1.1.0 needs matplotlib <3.9 (see envs/aniclustermap-env.yml). An
# environment solved before that pin was added imports but crashes, so check that it
# actually runs and rebuild it from the .yml if it does not.
if conda env list | awk '{print $1}' | grep -qx aniclustermap-env; then
    conda activate aniclustermap-env
    if ! ANIclustermap --help >/dev/null 2>&1; then
        echo "[setup] aniclustermap-env is broken; rebuilding from envs/aniclustermap-env.yml"
        conda deactivate
        conda env remove -y -n aniclustermap-env
        conda env create -y -f "$REPO/envs/aniclustermap-env.yml"
        conda env export -n aniclustermap-env --no-builds > "$REPO/envs/aniclustermap-env.lock.yml"
        conda activate aniclustermap-env
        ANIclustermap --help >/dev/null 2>&1 \
            && echo "[setup] aniclustermap-env now works" \
            || { echo "[setup] ANIclustermap still fails after rebuild" >&2; exit 1; }
    else
        echo "[setup] aniclustermap-env works"
    fi
    conda deactivate
fi

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

# ---- 5. GVMAGs V2 metadata (Vasquez et al. 2025) ------------------------------
# Metadata for all 18,727 giant virus genomes in the GVMAGs V2 release: taxonomy,
# GVClass results and the environment each genome came from. Used by
# 09_gvmag_votus.slurm to look up what our GVMAGs' nearest relatives are and where
# they were found. Small (~5 MB), so it is always downloaded.
#   https://yvasquez-lbl.github.io/GVMAGsV2/   doi:10.1101/2025.09.26.678796
GVMAGS_V2_DIR="$DB_DIR/gvmagsV2"
if [[ ! -s "$GVMAGS_V2_DIR/GVDB_stats.tsv" ]]; then
    echo "[setup] downloading GVMAGs V2 metadata (GVDB_stats.tsv)"
    mkdir -p "$GVMAGS_V2_DIR"
    curl -fL -o "$GVMAGS_V2_DIR/GVDB_stats.tsv.part" \
        https://github.com/yvasquez-lbl/GVMAGsV2/raw/main/additional_files/GVDB_stats.tsv
    mv "$GVMAGS_V2_DIR/GVDB_stats.tsv.part" "$GVMAGS_V2_DIR/GVDB_stats.tsv"
fi
echo "[setup] GVMAGs V2 metadata: $(wc -l < "$GVMAGS_V2_DIR/GVDB_stats.tsv") lines"

# Optional: the full GVMAGs V2 genome + protein bundle (1.18 GB), needed only to
# compare our GVMAGs against the 8,508 reference representatives by ANI
# (09_gvmag_votus.slurm --with-refs). Download it by resubmitting this script with
#   sbatch --export=ALL,GVMAGS_V2_FULL=1 ...
if [[ "${GVMAGS_V2_FULL:-0}" == "1" && ! -d "$GVMAGS_V2_DIR/database" ]]; then
    echo "[setup] downloading the full GVMAGs V2 database (1.18 GB)"
    curl -fL -o "$GVMAGS_V2_DIR/GVMAGS_V2_database.gz" \
        https://portal.nersc.gov/cfs/nelli/gvmagsV2/GVMAGS_V2_database.gz
    tar -xzf "$GVMAGS_V2_DIR/GVMAGS_V2_database.gz" -C "$GVMAGS_V2_DIR"
    rm "$GVMAGS_V2_DIR/GVMAGS_V2_database.gz"
fi

echo "[setup] done"
