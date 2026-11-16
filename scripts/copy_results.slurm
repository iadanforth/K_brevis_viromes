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

# Copy the step 08 (gvmags) and step 09 (votus) results out of scratch and into the
# repository, so they survive the 90-day scratch purge and sit next to the scripts that
# made them.
#
#   ~/scr10/data/kbrevis/gvmags  ->  ~/K_brevis_viromes/output/gvmags
#   ~/scr10/data/kbrevis/votus   ->  ~/K_brevis_viromes/output/votus
#
# This COPIES; the scratch originals are left untouched.
#
# TWO DIRECTORIES ARE EXCLUDED, because the full gvmags tree is 115 GB and the home
# filesystem gives this account a 50 GB quota (21 GB free), so copying all of it is
# impossible, not merely unwise:
#   binning/     108 GB. The bbmap BAM, its index and the jgi depth table for each
#                sample. These are the mapping intermediates MetaBAT2 consumed; nothing
#                downstream reads them again, and they are regenerable from the step 06
#                contigs plus the step 02 reads by re-running 08_gvmags.slurm.
#   gvclass_in/  1.9 GB. A byte-for-byte duplicate of bins/, differing only in file
#                naming (bin.10.fa -> bin_10.fna) because GVClass requires .fna
#                extensions. Verified identical by size on 2026-11-16; excluding it
#                loses no information at all.
# Everything else is kept, including bins/ (all 848 MetaBAT2 bins, 1.9 GB), so the
# classification half of step 08 can be re-run without redoing the 108 GB of mapping.
#
# Submit (1 core is right: this is pure I/O, and AGENTS.md section 6 requires GB-scale
# copies between filesystems to run in a batch job rather than on dev):
#   ( module purge; module load slurm/kuro; sbatch -p serial scripts/copy_results.slurm )

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
