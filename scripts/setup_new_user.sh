#!/bin/bash
# Adapt this clone to the account running it.
#
# Most settings live in pipeline.conf and are read at run time. Two things cannot be:
#
#   1. #SBATCH directives. Slurm parses them before any shell runs, so they cannot expand
#      ~, $HOME or a variable. AGENTS.md requires every job's .out/.err to land in
#      <project>/logs/, which forces an absolute path into 44 lines across the scripts.
#   2. --mail-user, for the same reason.
#
# This script rewrites both to match the current account. It is idempotent - running it
# twice changes nothing the second time - and --dry-run shows the edits without making
# them.
#
# Usage:
#   bash scripts/setup_new_user.sh --dry-run
#   bash scripts/setup_new_user.sh
set -euo pipefail

DRY_RUN=0
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN=1

HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$HERE"
[[ -f pipeline.conf ]] || { echo "pipeline.conf not found in $HERE" >&2; exit 1; }
# shellcheck disable=SC1091
set -a; source pipeline.conf; set +a

# The values currently baked into the scripts, as committed.
OLD_ROOT="/sciclone/home/iadanforth/K_brevis_viromes"
OLD_MAIL="iadanforth@wm.edu"

NEW_ROOT="$PROJECT_ROOT"
NEW_MAIL="$MAIL_USER"

echo "Adapting this clone"
echo "  project root : $OLD_ROOT"
echo "               -> $NEW_ROOT"
echo "  mail-user    : $OLD_MAIL"
echo "               -> $NEW_MAIL"
echo

if [[ "$NEW_ROOT" != /* ]]; then
    echo "PROJECT_ROOT must be an absolute path (got '$NEW_ROOT')" >&2; exit 1
fi
if [[ "$NEW_ROOT" != "$HERE" ]]; then
    echo "WARNING: PROJECT_ROOT in pipeline.conf is '$NEW_ROOT' but this clone is at" >&2
    echo "         '$HERE'. #SBATCH log paths would point somewhere this clone is not," >&2
    echo "         and Slurm fails a job whose log directory does not exist." >&2
    echo "         Fix PROJECT_ROOT, or move the clone, before continuing." >&2
    exit 1
fi

# `grep ... | head` makes grep die of SIGPIPE once head has what it wants, which is fatal
# under `set -o pipefail`. Every such pipeline below is guarded with `|| true`.
n_root=$(grep -rl "$OLD_ROOT" scripts/ 2>/dev/null | wc -l || true)
n_mail=$(grep -rl "$OLD_MAIL" scripts/ 2>/dev/null | wc -l || true)
echo "  files containing the old root : $n_root"
echo "  files containing the old mail : $n_mail"

# The question is whether the values DIFFER, not whether the old strings are present. On
# the original account they are present in 24 files and there is still nothing to do.
if [[ "$NEW_ROOT" == "$OLD_ROOT" && "$NEW_MAIL" == "$OLD_MAIL" ]]; then
    echo
    echo "Nothing to change: pipeline.conf already matches what the scripts contain."
    echo "(This is the original account, or the clone has been adapted already.)"
    exit 0
fi

if (( DRY_RUN )); then
    echo
    echo "--dry-run: showing the first few lines that would change"
    { grep -rn "$OLD_ROOT" scripts/ 2>/dev/null || true; } | head -5 | sed 's/^/    /'
    { grep -rn "$OLD_MAIL" scripts/ 2>/dev/null || true; } | head -3 | sed 's/^/    /'
    echo
    echo "Re-run without --dry-run to apply."
    exit 0
fi

if [[ "$NEW_ROOT" != "$OLD_ROOT" ]]; then
    grep -rl "$OLD_ROOT" scripts/ | while read -r f; do
        sed -i "s|$OLD_ROOT|$NEW_ROOT|g" "$f"
    done
    echo "  rewrote project root in $n_root files"
fi
if [[ "$NEW_MAIL" != "$OLD_MAIL" ]]; then
    grep -rl "$OLD_MAIL" scripts/ | while read -r f; do
        sed -i "s|$OLD_MAIL|$NEW_MAIL|g" "$f"
    done
    echo "  rewrote mail-user in $n_mail files"
fi

# Prove the result still parses: a bad substitution that breaks a script is worse than
# leaving it alone, and sed cannot tell.
echo
echo "Verifying every script still parses"
bad=0
for f in scripts/*.sh scripts/*.slurm; do
    bash -n "$f" 2>/dev/null || { echo "  ** SYNTAX ERROR after rewrite: $f" >&2; bad=1; }
done
for f in scripts/*.py; do
    python3 -I -c "import ast,sys; ast.parse(open(sys.argv[1]).read())" "$f" 2>/dev/null \
        || { echo "  ** SYNTAX ERROR after rewrite: $f" >&2; bad=1; }
done
(( bad )) && { echo "Rewrite produced a broken script; use 'git checkout -- scripts/' to undo." >&2; exit 1; }
echo "  all scripts parse"

left=$(grep -rl "$OLD_ROOT\|$OLD_MAIL" scripts/ 2>/dev/null | wc -l || true)
echo "  files still mentioning the old values: $left (expect 0)"
echo
echo "Done. Next: bash run_pipeline.sh --list"
