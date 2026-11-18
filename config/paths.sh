# Central path configuration for shell scripts. Source it from the repo root:
#
#   source "${ADSEX_ROOT:-.}/config/paths.sh"

PROJECT_ROOT="$(cd "${ADSEX_ROOT:-.}" && pwd)"
if [ ! -f "$PROJECT_ROOT/config/paths.sh" ]; then
    echo "Project root not found. Run from the repository root or set ADSEX_ROOT." >&2
    exit 1
fi

DATA_DIR="${ADSEX_DATA:-$PROJECT_ROOT/data}"
RESULTS_DIR="${ADSEX_RESULTS:-$PROJECT_ROOT/results}"
