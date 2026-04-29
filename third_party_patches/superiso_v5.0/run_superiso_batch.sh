#!/bin/bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
NP_DIR="$ROOT_DIR/np_data"
LOG_DIR="$NP_DIR/superiso_batch_logs"
CHECK_SUPERISO="$ROOT_DIR/third_party_patches/superiso_v5.0/check_superiso.py"
JOBS="${1:-8}"
CHECKPOINT_EVERY="${2:-50}"

mkdir -p "$LOG_DIR"

cd "$ROOT_DIR"

find "$NP_DIR" -maxdepth 1 -type f -name 'scan_after*.dat' \
  ! -name '*_superiso*.dat' \
  ! -name '*_logcheck*.dat' \
  ! -name '*_test.dat' \
  ! -name '*_serial_*.dat' \
  ! -name '*_parallel_*.dat' | sort | while read -r input_file; do
    base_name="$(basename "$input_file" .dat)"
    output_file="$NP_DIR/${base_name}_superiso.dat"
    log_file="$LOG_DIR/${base_name}_superiso.log"

    if [ -f "$output_file" ] && [ -f "$log_file" ] && grep -q '\] DONE ' "$log_file"; then
      printf '[%s] SKIP completed output: %s\n' "$(date '+%F %T')" "$output_file" | tee -a "$log_file"
      continue
    fi

    if [ -f "$output_file" ]; then
      printf '[%s] RESUME partial output: %s\n' "$(date '+%F %T')" "$output_file" | tee -a "$log_file"
    else
      printf '[%s] START %s\n' "$(date '+%F %T')" "$input_file" | tee "$log_file"
      printf 'COMMAND: %s\n' "python3 -u third_party_patches/superiso_v5.0/check_superiso.py np_data/$(basename "$input_file") --mode chi2 --append-table --jobs $JOBS --checkpoint-every $CHECKPOINT_EVERY --output-file np_data/$(basename "$output_file")" | tee -a "$log_file"
    fi

    python3 -u "$CHECK_SUPERISO" "np_data/$(basename "$input_file")" \
      --mode chi2 \
      --append-table \
      --jobs "$JOBS" \
      --checkpoint-every "$CHECKPOINT_EVERY" \
      --output-file "np_data/$(basename "$output_file")" >> "$log_file" 2>&1

    printf '[%s] DONE %s\n' "$(date '+%F %T')" "$output_file" | tee -a "$log_file"
  done