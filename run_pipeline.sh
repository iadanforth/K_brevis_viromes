#!/bin/bash
# =======================================================================================
# run_pipeline.sh - reproduce the whole K. brevis virome analysis with one command.
#
#   bash run_pipeline.sh --list              what the stages are, and what each produces
#   bash run_pipeline.sh --dry-run all       print every submission without submitting
#   bash run_pipeline.sh all                 queue the entire pipeline
#   bash run_pipeline.sh --from 08 --to 11   queue part of it
#   bash run_pipeline.sh --only 14           queue one stage
#   bash run_pipeline.sh --status            what is queued or running, on every cluster
#   bash run_pipeline.sh --check             which stages have produced their outputs
#
# HOW IT WORKS. Each stage is submitted with `--dependency=afterok:<previous job ids>`,
# so the whole run is queued in seconds and Slurm releases each stage only when the one
# before it has succeeded. The script then exits - it does not sit on the front-end
# waiting, which would trip the watchdog and die at logout anyway.
#
# IF A STAGE FAILS, everything after it stays PENDING with reason DependencyNeverSatisfied
# (Slurm may also cancel them). That is the intended behaviour: it stops a bad assembly
# from propagating. Read the failed stage's .err in logs/, fix it, and resume with
# `--from <that stage>`.
#
# RESUMABLE. Stages skip work that is already done - 05a/05b/06 skip finished samples,
# 08 skips samples with a .done marker, 14 skips samples with a CoverM table, and 07 and
# 00 skip environments that exist. Re-running a completed stage is cheap and safe.
#
# BEFORE THE FIRST RUN: edit pipeline.conf, then `bash scripts/setup_new_user.sh`.
#
# SCOPE. This reproduces the pipeline on a SciClone account. It assumes Slurm with the
# per-cluster `slurm/<name>` modules, `module load miniforge3`, and the ~/scr10 and
# ~/data10 filesystems. On a different cluster the analysis is unchanged but the job
# plumbing is not - see "Porting to a different account" in README.md.
# =======================================================================================
set -uo pipefail

HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"
[[ -f pipeline.conf ]] || { echo "pipeline.conf not found in $HERE" >&2; exit 1; }
# shellcheck disable=SC1091
set -a; source pipeline.conf; set +a

DRY_RUN=0; FROM=""; TO=""; ONLY=""; ACTION=""
while (($#)); do
    case "$1" in
        --dry-run) DRY_RUN=1; shift ;;
        --from)    FROM="$2"; shift 2 ;;
        --to)      TO="$2";   shift 2 ;;
        --only)    ONLY="$2"; shift 2 ;;
        --list|--status|--check) ACTION="${1#--}"; shift ;;
        all)       shift ;;
        -h|--help) sed -n '2,33p' "$0"; exit 0 ;;
        *) echo "Unknown argument: $1 (try --help)" >&2; exit 2 ;;
    esac
done

# --- stage table: name | description | what proves it finished --------------------------
STAGES=(00_envs 01_download 02_fastp 03_refs 04_clean 05a_assembly 05b_assembly 06_merge
        07_setup 08_gvmags 09_votus 10_phylogeny 11_tree 12_collect 13_itol 14_abundance
        15_ani)

describe () {
    case "$1" in
        00_envs)      echo "conda environments for steps 01-06, + dirs + sample table" ;;
        01_download)  echo "download 20 SRA runs                     -> raw/ (~2.0 TB)" ;;
        02_fastp)     echo "adapter/quality trimming (fastp 1.3.7)   -> trim/ (~1.6 TB)" ;;
        03_refs)      echo "host + contaminant genomes, PhiX         -> reference/ (~79 GB)" ;;
        04_clean)     echo "bowtie2 host removal                     -> clean/ (~1.1 TB)" ;;
        05a_assembly) echo "metaviralSPAdes, 20 jobs                 -> contigs/ (90 GB)" ;;
        05b_assembly) echo "metaSPAdes, 20 jobs (18 succeed)         -> contigs_no_v/ (322 GB)" ;;
        06_merge)     echo "merge the two assemblies, 20 jobs        -> contigs_merged/ (6.6 GB)" ;;
        07_setup)     echo "GVClass, CheckV/geNomad/GVMAGs V2 dbs, 6 more envs" ;;
        08_gvmags)    echo "map, bin, classify; array                -> gvmags/ (115 GB)" ;;
        09_votus)     echo "cluster GVMAGs into vOTUs                -> votus/ (16 MB)" ;;
        10_phylogeny) echo "marker genes + concatenated alignment    -> phylogeny/ (629 MB)" ;;
        11_tree)      echo "IQ-TREE ML tree + RAxML IC/TC            -> tree/ (180 MB)" ;;
        12_collect)   echo "copy results out of scratch              -> output/" ;;
        13_itol)      echo "iTOL annotations + vOTU tree placement   -> output/tree/" ;;
        14_abundance) echo "CoverM abundance through time + figures  -> output/abundance/" ;;
        15_ani)       echo "ANI vs the published GVMAGs              -> output/*.csv" ;;
    esac
}

# A cheap test of whether a stage's output already exists. Used by --check and reported
# before submission so you can see what will be skipped.
verify () {
    case "$1" in
        00_envs)      [[ -s "$SCRATCH/SraRunTable.csv" ]] ;;
        01_download)  (( $(ls "$SCRATCH"/raw/*.fastq 2>/dev/null | wc -l) >= 40 )) ;;
        02_fastp)     (( $(ls "$SCRATCH"/trim/*_TRIMMED.fastq 2>/dev/null | wc -l) >= 40 )) ;;
        03_refs)      (( $(ls "$SCRATCH"/reference/*.fna 2>/dev/null | wc -l) >= 6 )) ;;
        04_clean)     (( $(ls "$SCRATCH"/clean/*_CLEAN.fastq 2>/dev/null | wc -l) >= 40 )) ;;
        05a_assembly) (( $(ls "$SCRATCH"/contigs/*_contigs.fasta 2>/dev/null | wc -l) >= 20 )) ;;
        # 18, not 20, and that is correct. SRR28004142 and SRR28004143 are the largest
        # read sets and their metaSPAdes jobs were cancelled for time (README 05b-fix-c).
        # Step 06 handles them through MV_ONLY_SAMPLES, using their metaviralSPAdes contigs
        # alone, so 18 metaSPAdes assemblies is the finished state of this stage. If your
        # run produces 20, so much the better - edit MV_ONLY_SAMPLES in 06_merge_assemblies.sh
        # to empty so those two samples get a real merge.
        05b_assembly) (( $(ls "$SCRATCH"/contigs_no_v/*_contigs.fasta 2>/dev/null | wc -l) >= 18 )) ;;
        06_merge)     (( $(ls "$SCRATCH"/contigs_merged/*_merged_contigs.fasta 2>/dev/null | wc -l) >= 20 )) ;;
        07_setup)     [[ -d "$ARCHIVE/databases/gvmagsV2" ]] ;;
        08_gvmags)    (( $(ls "$SCRATCH"/gvmags/*/.done 2>/dev/null | wc -l) >= 20 )) ;;
        09_votus)     [[ -s "$SCRATCH/votus/gvmags_votus.tsv" ]] ;;
        10_phylogeny) [[ -s "$SCRATCH/phylogeny/kbrevis_gvmags.trimmed.aln" ]] ;;
        11_tree)      [[ -s "$SCRATCH/tree/kbrevis_gvmags.treefile" ]] ;;
        12_collect)   [[ -s "$HERE/output/tree/kbrevis_gvmags.treefile" ]] ;;
        13_itol)      [[ -s "$HERE/output/tree/votu_tree_context.tsv" ]] ;;
        14_abundance) [[ -s "$HERE/output/abundance/votu_abundance_long.tsv" ]] ;;
        15_ani)       [[ -s "$HERE/output/booker_gvmag_ani_summary.csv" ]] ;;
        *) return 1 ;;
    esac
}

opts_for () { local v="STAGE_OPTS_$1"; echo "${!v:-${STAGE_OPTS_helper}}"; }
cluster_of () { opts_for "$1" | cut -d: -f1; }
sbargs_of  () { opts_for "$1" | cut -d: -f2- ; }

# --- submission helpers -----------------------------------------------------------------
# Every message goes to stderr; only job IDs go to stdout, so callers can capture them.
LOGS="$HERE/logs"; mkdir -p "$LOGS"

run_in_cluster () {   # run_in_cluster <cluster> <command...>
    local cluster="$1"; shift
    ( module purge >/dev/null 2>&1
      [[ -n "$cluster" && "$cluster" != "none" ]] && module load "slurm/$cluster" >/dev/null 2>&1
      "$@" )
}

submit_script () {    # submit_script <stage> <dep> <script> [script args...]
    local stage="$1" dep="$2"; shift 2
    local cluster; cluster=$(cluster_of "$stage")
    local -a args; read -r -a args <<< "$(sbargs_of "$stage")"
    local -a dflag=(); [[ -n "$dep" ]] && dflag=(--dependency="afterok:$dep")
    if (( DRY_RUN )); then
        echo "  [dry-run] ( module load slurm/$cluster; sbatch --parsable ${dflag[*]} ${args[*]} $* )" >&2
        echo "DRYRUN$RANDOM"; return 0
    fi
    run_in_cluster "$cluster" sbatch --parsable "${dflag[@]}" "${args[@]}" "$@"
}

submit_wrap () {      # submit_wrap <stage> <dep> <name> <shell command string>
    local stage="$1" dep="$2" name="$3" cmd="$4"
    local cluster; cluster=$(cluster_of "$stage")
    local -a args; read -r -a args <<< "$(sbargs_of "$stage")"
    local -a dflag=(); [[ -n "$dep" ]] && dflag=(--dependency="afterok:$dep")
    if (( DRY_RUN )); then
        echo "  [dry-run] ( module load slurm/$cluster; sbatch --parsable ${dflag[*]} ${args[*]} -J $name --wrap \"$cmd\" )" >&2
        echo "DRYRUN$RANDOM"; return 0
    fi
    run_in_cluster "$cluster" sbatch --parsable "${dflag[@]}" "${args[@]}" \
        -J "$name" -o "$LOGS/${name}_%j.out" -e "$LOGS/${name}_%j.err" \
        --wrap "set -euo pipefail; cd '$HERE'; $cmd"
}

submit_multi () {     # submit_multi <stage> <dep> <driver.sh> [args...]  -> colon-joined ids
    local stage="$1" dep="$2"; shift 2
    local cluster; cluster=$(cluster_of "$stage")
    local -a args; read -r -a args <<< "$(sbargs_of "$stage")"
    local -a dflag=(); [[ -n "$dep" ]] && dflag=(--dependency="afterok:$dep")
    if (( DRY_RUN )); then
        echo "  [dry-run] ( module load slurm/$cluster; bash $* ${args[*]} ${dflag[*]} )  # submits many" >&2
        echo "DRYRUN$RANDOM"; return 0
    fi
    local out
    out=$(run_in_cluster "$cluster" bash "$@" "${args[@]}" "${dflag[@]}" 2>&1)
    echo "$out" | sed 's/^/      /' >&2
    echo "$out" | sed -n 's/^SUBMITTED_JOB_IDS=//p' | tail -1
}

# --- the stages ---------------------------------------------------------------------------
submit_stage () {     # submit_stage <name> <dep>  -> job ids (colon-joined)
    local s="$1" dep="$2"
    case "$s" in
      00_envs)
        local a b
        a=$(submit_script "$s" "$dep" scripts/00_create_envs.slurm)
        b=$(submit_wrap "$s" "$a" copy_sra_table \
            "mkdir -p '$SCRATCH' && cp metadata/SraRunTable.csv '$SCRATCH/' && bash scripts/00_directory_setup.sh")
        echo "$a:$b" ;;
      01_download)  submit_script "$s" "$dep" scripts/01_download_files.slurm ;;
      02_fastp)     submit_script "$s" "$dep" scripts/02_fastp.slurm ;;
      03_refs)
        if [[ "${SKIP_REFERENCE_DOWNLOAD:-0}" == "1" ]]; then
            echo "  skipped: SKIP_REFERENCE_DOWNLOAD=1 in pipeline.conf" >&2
            echo ""; return 0
        fi
        local a b
        a=$(submit_script "$s" "$dep" scripts/03_download_refs.slurm)
        b=$(submit_script "$s" "$dep" scripts/03b_download_phix.slurm)
        echo "$a:$b" ;;
      04_clean)     submit_script "$s" "$dep" scripts/04_clean_refs.slurm ;;
      05a_assembly) submit_multi  "$s" "$dep" scripts/05_metaspades.sh ;;
      05b_assembly) submit_multi  "$s" "$dep" scripts/05_metaspades_nov.sh ;;
      06_merge)     submit_multi  "$s" "$dep" scripts/06_merge_assemblies.sh ;;
      07_setup)
        # GVMAGS_V2_FULL=1 is required: without it step 10 fails much later.
        local cluster; cluster=$(cluster_of "$s")
        local -a args; read -r -a args <<< "$(sbargs_of "$s")"
        local -a dflag=(); [[ -n "$dep" ]] && dflag=(--dependency="afterok:$dep")
        if (( DRY_RUN )); then
            echo "  [dry-run] sbatch ${args[*]} --export=ALL,GVMAGS_V2_FULL=1 scripts/07_gvmag_setup.slurm" >&2
            echo "DRYRUN$RANDOM"
        else
            run_in_cluster "$cluster" sbatch --parsable "${dflag[@]}" "${args[@]}" \
                --export=ALL,GVMAGS_V2_FULL=1 scripts/07_gvmag_setup.slurm
        fi ;;
      08_gvmags)    submit_multi  "$s" "$dep" scripts/08_gvmags.sh ;;
      09_votus)     submit_script "$s" "$dep" scripts/09_gvmag_votus.slurm --with-refs ;;
      10_phylogeny) submit_script "$s" "$dep" scripts/10_gvmag_phylogeny.slurm ;;
      11_tree)      submit_script "$s" "$dep" scripts/11_gvmag_tree.slurm \
                        --mset "$IQTREE_MSET" --cmax "$IQTREE_CMAX" --safe ;;
      12_collect)   submit_script "$s" "$dep" scripts/copy_results.slurm ;;
      13_itol)
        submit_wrap "$s" "$dep" itol_annotations "
python3 -I scripts/make_itol_annotations.py \
  --tree output/tree/kbrevis_gvmags.treefile \
  --manifest '$SCRATCH/phylogeny/reference_manifest.tsv' \
  --votus output/votus/gvmags_votus.tsv \
  --nearest output/votus/votu_nearest_reference.tsv \
  --out-dir output/tree
python3 -I scripts/votu_tree_context.py \
  --tree output/tree/kbrevis_gvmags.treefile \
  --manifest '$SCRATCH/phylogeny/reference_manifest.tsv' \
  --out output/tree/votu_tree_context.tsv" ;;
      14_abundance)
        local a b c d
        a=$(submit_script "14_abundance" "$dep" scripts/14_abundance_setup.slurm)
        b=$(submit_wrap   "14_abundance" "$a" gvmag_reference "bash scripts/14_gvmag_reference.sh")
        c=$(submit_script "14_abundance" "$b" scripts/14_gvmag_abundance.slurm)
        d=$(submit_wrap   "14_abundance" "$c" abundance_merge_plot "
mkdir -p output/abundance
python3 -I scripts/merge_gvmag_abundance.py \
  --coverm-dir '$SCRATCH/abundance/per_gvmag' \
  --manifest '$SCRATCH/abundance/gvmag_reference_manifest.tsv' \
  --sra-table metadata/SraRunTable.csv \
  --min-covered-fraction $ABUNDANCE_MIN_BREADTH \
  --out-dir output/abundance
python3 -I scripts/votu_sample_table.py \
  --votus output/votus/gvmags_votus.tsv \
  --coverm-dir '$SCRATCH/abundance/per_gvmag/strict' \
  --sra-table metadata/SraRunTable.csv \
  --incubation '$INCUBATION_SAMPLES' \
  --min-breadth $ABUNDANCE_MIN_BREADTH \
  --out-long output/votu_sample_table.csv \
  --out-matrix output/votu_sample_matrix.csv
module load $CONDA_MODULE
eval \"\\\$(conda shell.bash hook)\"
conda activate r-plot-env
Rscript scripts/plot_votu_sample_matrix.R output/votu_sample_table.csv output
Rscript scripts/plot_incubation_votu_timeseries.R output/abundance/votu_rollup.tsv output/abundance paper")
        echo "$d" ;;
      15_ani)
        local a b
        a=$(submit_wrap   "15_ani" "$dep" fetch_booker "bash scripts/15_fetch_booker_gvmags.sh")
        b=$(submit_script "15_ani" "$a" scripts/15_booker_ani.slurm)
        echo "$b" ;;
    esac
}

# --- actions --------------------------------------------------------------------------
if [[ "$ACTION" == "list" ]]; then
    printf "%-14s %-6s %s\n" STAGE CLUSTER DESCRIPTION
    for s in "${STAGES[@]}"; do
        printf "  %-12s %-6s %s\n" "$s" "$(cluster_of "$s")" "$(describe "$s")"
    done
    echo
    echo "Config: PROJECT_ROOT=$HERE"
    echo "        SCRATCH=$SCRATCH"
    echo "        MAIL_USER=$MAIL_USER"
    exit 0
fi

if [[ "$ACTION" == "check" ]]; then
    echo "Stage outputs present on disk:"
    for s in "${STAGES[@]}"; do
        if verify "$s"; then printf "  %-14s DONE\n" "$s"; else printf "  %-14s -\n" "$s"; fi
    done
    exit 0
fi

if [[ "$ACTION" == "status" ]]; then
    for c in kuro bora astral femto gulf; do
        out=$(run_in_cluster "$c" squeue -u "$USER" -h -o '%.12i %.22j %.10T %.10M %R' 2>/dev/null)
        [[ -n "$out" ]] && { echo "--- $c"; echo "$out"; }
    done
    exit 0
fi

# --- select the stages to run -----------------------------------------------------------
selected=()
if [[ -n "$ONLY" ]]; then
    for s in "${STAGES[@]}"; do [[ "$s" == "$ONLY"* ]] && selected+=("$s"); done
else
    started=0
    for s in "${STAGES[@]}"; do
        [[ -z "$FROM" || "$s" == "$FROM"* ]] && started=1
        (( started )) && selected+=("$s")
        [[ -n "$TO" && "$s" == "$TO"* ]] && break
    done
fi
if ((${#selected[@]} == 0)); then
    echo "No stages matched. Try --list." >&2; exit 2
fi

echo "K. brevis virome pipeline"
echo "  project : $HERE"
echo "  scratch : $SCRATCH"
echo "  stages  : ${selected[*]}"
(( DRY_RUN )) && echo "  MODE    : DRY RUN - nothing will be submitted"
echo

dep=""
for s in "${selected[@]}"; do
    status="will run"
    verify "$s" && status="output already present (the stage will skip finished work)"
    echo "[$s] $(describe "$s")"
    echo "  $status"
    ids=$(submit_stage "$s" "$dep")
    ids="${ids//[[:space:]]/}"
    if [[ -n "$ids" ]]; then
        echo "  job(s): $ids"
        dep="$ids"
    else
        echo "  nothing submitted; later stages keep the previous dependency"
    fi
    echo
done

if (( DRY_RUN )); then
    echo "Dry run complete - nothing was submitted."
else
    echo "All stages queued. Each starts only when the one before it succeeds."
    echo "  watch:  bash run_pipeline.sh --status"
    echo "  verify: bash run_pipeline.sh --check"
    echo "  logs:   $LOGS"
    echo
    echo "If a stage fails, later stages stay PENDING (DependencyNeverSatisfied)."
    echo "Read its .err, fix it, then resume with: bash run_pipeline.sh --from <stage>"
fi
