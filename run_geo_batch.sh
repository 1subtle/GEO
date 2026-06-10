#!/usr/bin/env bash
set -euo pipefail

# Batch compile + parallel run + result table for GEO.cpp.
#
# Usage:
#   chmod +x run_geo_batch.sh
#   JOBS=16 SEEDS="1 2 3 4 5" TIME_LIMIT=3600 ./run_geo_batch.sh
#
# Outputs:
#   batch_results_YYYYmmdd_HHMMSS/logs/*.log
#   batch_results_YYYYmmdd_HHMMSS/results.csv
#   batch_results_YYYYmmdd_HHMMSS/results.md

CXX="${CXX:-g++}"
JOBS="${JOBS:-}"
SEEDS="${SEEDS:-1}"
SRC="${SRC:-geo.cpp}"
EXE="${EXE:-geo}"
INSTANCE_DIR="${INSTANCE_DIR:-geo_instances}"
OUTDIR="${OUTDIR:-batch_results_$(date +%Y%m%d_%H%M%S)}"
TIME_LIMIT="${TIME_LIMIT:-}"

if [ -z "$JOBS" ]; then
    if command -v nproc >/dev/null 2>&1; then
        JOBS="$(nproc)"
    elif command -v sysctl >/dev/null 2>&1; then
        JOBS="$(sysctl -n hw.ncpu)"
    else
        JOBS=4
    fi
fi

mkdir -p "$OUTDIR/logs"

echo "Compiling $SRC -> $EXE"
"$CXX" -O3 "$SRC" -o "$EXE"

JOBS_FILE="$OUTDIR/jobs.tsv"
: > "$JOBS_FILE"

find "$INSTANCE_DIR" -type f -name "*.txt" | sort | while IFS= read -r inst; do
    base="$(basename "$inst" .txt)"
    for seed in $SEEDS; do
        log="$OUTDIR/logs/${base}_seed${seed}.log"
        printf '%s\t%s\t%s\n' "$inst" "$seed" "$log" >> "$JOBS_FILE"
    done
done

run_one() {
    inst="$1"
    seed="$2"
    log="$3"
    if [ -n "$TIME_LIMIT" ]; then
        "./$EXE" "$inst" "$seed" "$TIME_LIMIT" > "$log" 2>&1
    else
        "./$EXE" "$inst" "$seed" > "$log" 2>&1
    fi
}
export -f run_one
export EXE
export TIME_LIMIT

echo "Running $(wc -l < "$JOBS_FILE" | tr -d ' ') jobs with JOBS=$JOBS"
xargs -P "$JOBS" -n 3 bash -c 'run_one "$@"' _ < "$JOBS_FILE"

CSV="$OUTDIR/results.csv"
MD="$OUTDIR/results.md"

printf 'instance,seed,best_profit,best_time,elapsed_time,tasks_served,num_tasks,verified,log\n' > "$CSV"

while IFS=$'\t' read -r inst seed log; do
    instance="$(basename "$inst")"
    best_profit="$(awk '
        /^(ILS done|MSBTS-GEO done)/ {
            for(i=1;i<=NF;i++) if($i ~ /^bestProfit=/){split($i,a,"="); bp=a[2]}
        }
        /Total profit[[:space:]]*:/ {tp=$4}
        END { if(bp!="") print bp; else if(tp!="") print tp }
    ' "$log")"
    best_time="$(awk '
        /^(ILS done|MSBTS-GEO done)/ {
            for(i=1;i<=NF;i++) if($i ~ /^bestTime=/){split($i,a,"="); bt=a[2]}
        }
        END { if(bt!="") print bt }
    ' "$log")"
    elapsed="$(awk '/Elapsed time:/ {print $3}' "$log" | tail -1)"
    served="$(awk '/Tasks served/ {print $4}' "$log" | tail -1)"
    total_tasks="$(awk '/Tasks served/ {print $6}' "$log" | tail -1)"
    verified="NO"
    if grep -q "Solution verified OK" "$log"; then
        verified="YES"
    fi
    log_path="$(cd "$(dirname "$log")" && pwd)/$(basename "$log")"
    printf '%s,%s,%s,%s,%s,%s,%s,%s,%s\n' \
        "$instance" "$seed" "$best_profit" "$best_time" "$elapsed" \
        "$served" "$total_tasks" "$verified" "$log_path" >> "$CSV"
done < "$JOBS_FILE"

{
    echo "| instance | seed | best_profit | best_time_s | elapsed_s | tasks_served | num_tasks | verified | log |"
    echo "|---|---:|---:|---:|---:|---:|---:|---|---|"
    tail -n +2 "$CSV" | while IFS=',' read -r instance seed best_profit best_time elapsed served total_tasks verified log; do
        printf '| %s | %s | %s | %s | %s | %s | %s | %s | %s |\n' \
            "$instance" "$seed" "$best_profit" "$best_time" "$elapsed" \
            "$served" "$total_tasks" "$verified" "$log"
    done
} > "$MD"

echo "Done."
echo "CSV table: $CSV"
echo "Markdown table: $MD"
echo "Logs: $OUTDIR/logs"
