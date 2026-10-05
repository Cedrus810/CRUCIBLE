#!/usr/bin/env bash
# T0: run a command with a START line, per-line wall-clock timestamps and an
# END line (exit code + total wall time).  The command is received as an
# array ("$@"); shell functions and aliases are not supported (plan §1.7).
#
#   scripts/run_logged.sh python -m pytest --durations=0 > jobs/x/pytest.log 2>&1

if [ "$#" -lt 1 ]; then
    echo "usage: $0 <command> [args...]" >&2
    exit 64
fi

t0=$(date +%s)
stamp() {
    while IFS= read -r line || [ -n "$line" ]; do
        printf '[%s +%ss] %s\n' "$(date '+%F %T')" "$(( $(date +%s) - t0 ))" "$line"
    done
}

commit=$(git -c safe.directory="$PWD" rev-parse --short HEAD 2>/dev/null || echo unknown)
echo "START host=$(hostname) commit=${commit} argv=$*" | stamp
"$@" 2>&1 | stamp
status=${PIPESTATUS[0]}
echo "END exit=${status} wall=$(( $(date +%s) - t0 ))s" | stamp
exit "${status}"
