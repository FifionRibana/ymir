#!/usr/bin/env bash
# Aggregate Phase 1-bis diagnostic metrics from per-scenario log files.
#
# Reads `logs/phase1bis_<S>_<NN>.log` for S in {A,B,C,D} and NN in {01,02,03},
# plus `logs/summary.txt` for wallclock.
#
# Emits tab-separated tables for each subsection of the report
# (wallclock, phase percentages, RHS spike, η distribution, residual
# localization) restricted to macro-step range [step_lo, step_hi] for the
# steady-state metrics. Range is 100..200 by default; override via
# PHASE1BIS_STEP_LO and PHASE1BIS_STEP_HI.

set -euo pipefail

LOG_DIR="${1:-logs}"
STEP_LO="${PHASE1BIS_STEP_LO:-100}"
STEP_HI="${PHASE1BIS_STEP_HI:-200}"

cd "$(dirname "$0")/.."

if [[ ! -f "$LOG_DIR/summary.txt" ]]; then
    echo "missing $LOG_DIR/summary.txt — run all scenarios first" >&2
    exit 1
fi

awk_prog='
# Helper: extract numeric value for a "key=value" token on a line.
function kv(line, key,   v, re) {
    re = "(^| )" key "=[^ ]+"
    if (match(line, re)) {
        v = substr(line, RSTART, RLENGTH)
        sub(".*=", "", v)
        return v + 0
    }
    return -1
}
# Helper: extract integer "step=N" field from the solver_step span.
function step_of(line,   v, re) {
    re = "step=[0-9]+"
    if (match(line, re)) {
        v = substr(line, RSTART, RLENGTH)
        sub("step=", "", v)
        return v + 0
    }
    return -1
}

BEGIN { lo = step_lo + 0; hi = step_hi + 0 }

# Accumulate per-metric per-scenario. sc comes from the filename.
function acc(metric, val,   key) {
    if (val < 0) return
    key = sc "|" metric
    sum[key] += val
    cnt[key] += 1
    if (!(key in maxv) || val > maxv[key]) maxv[key] = val
}

{
    step = step_of($0)
    if (step < lo || step > hi) next

    if ($0 ~ /rhs_breakdown/) {
        # Extract 5 GPE + 5 T_plates + 5 total + 3 spike ratios.
        gm = kv($0, "gpe_rhs_max_abs")
        gs = kv($0, "gpe_rhs_p95")
        tm = kv($0, "tplates_rhs_max_abs")
        ts = kv($0, "tplates_rhs_p95")
        tn = kv($0, "tplates_rhs_norm")
        gn = kv($0, "gpe_rhs_norm")
        if (gm >= 0 && gs > 0) acc("gpe_spike_p95", gm / gs)
        if (tm >= 0 && ts > 0) acc("tp_spike_p95", tm / ts)
        if (gm >= 0)            acc("gpe_max_abs",   gm)
        if (tm >= 0)            acc("tp_max_abs",    tm)
        if (gn >= 0)            acc("gpe_norm",      gn)
        if (tn >= 0)            acc("tp_norm",       tn)
    }
    if ($0 ~ /eta_breakdown/) {
        er = kv($0, "eta_ratio")
        yf = kv($0, "yielding_cells_fraction")
        sc2 = kv($0, "saturated_cells_count")
        if (er >= 0) acc("eta_ratio",  er)
        if (yf >= 0) acc("yield_frac", yf)
        if (sc2 >= 0) acc("saturated", sc2)
    }
    if ($0 ~ /residual_spatial/) {
        rl = kv($0, "residual_localization")
        if (rl >= 0) acc("resid_local", rl)
    }
    if ($0 ~ /phase_timings/) {
        tb = kv($0, "t_boundaries_us")
        ts = kv($0, "t_solve_us")
        ta = kv($0, "t_advection_us")
        tr = kv($0, "t_recycling_us")
        tp = kv($0, "t_plates_us")
        acc("t_boundaries", tb)
        acc("t_solve",       ts)
        acc("t_advection",   ta)
        acc("t_recycling",   tr)
        acc("t_plates",      tp)
    }
}

END {
    for (k in sum) {
        n = cnt[k]
        if (n > 0) printf "%s MEAN %.6g\n", k, sum[k]/n
        if (k in maxv) printf "%s MAX %.6g\n", k, maxv[k]
    }
}
'

echo "=== Wallclock (from SUMMARY lines) ==="
awk '/SUMMARY/ {
    sc = gensub(/.*scenario=([A-D]).*/, "\\1", 1)
    ms = gensub(/.*elapsed_ms=([0-9]+).*/, "\\1", 1) + 0
    s[sc] += ms; n[sc] += 1
    if (!(sc in lo) || ms < lo[sc]) lo[sc] = ms
    if (!(sc in hi) || ms > hi[sc]) hi[sc] = ms
}
END {
    printf "scenario\tmean_s\tmin_s\tmax_s\n"
    for (sc in s) printf "%s\t%.2f\t%.2f\t%.2f\n", sc, s[sc]/n[sc]/1000, lo[sc]/1000, hi[sc]/1000
}' "$LOG_DIR/summary.txt" | sort

echo
echo "=== Per-scenario aggregates (steps $STEP_LO..$STEP_HI, 3 reps merged) ==="
for sc in A B C D; do
    files=( $LOG_DIR/phase1bis_${sc}_*.log )
    if [[ ! -f "${files[0]}" ]]; then
        echo "scenario $sc: no log files found" >&2
        continue
    fi
    awk -v sc="$sc" -v step_lo="$STEP_LO" -v step_hi="$STEP_HI" "$awk_prog" "${files[@]}"
done
