#!/bin/bash
# Pell equation benchmark: Egypt vs PARI/GP
# Simulates auto-doubling precision for fair comparison
#
# Usage: ./scripts/bench_pell.sh [D] [runs]
#   D     - discriminant (default: 4729494 = cattle problem)
#   runs  - number of benchmark runs (default: 10)

D=${1:-4729494}
RUNS=${2:-10}
EGYPT=${EGYPT:-./target/release/egypt}

echo "=== Pell Equation Benchmark ==="
echo "D = $D"
echo "Runs = $RUNS"
echo ""

# Check dependencies
if ! command -v gp &> /dev/null; then
    echo "Warning: PARI/GP (gp) not found, skipping comparison"
    PARI_AVAILABLE=0
else
    PARI_AVAILABLE=1
fi

if [ ! -x "$EGYPT" ]; then
    echo "Error: $EGYPT not found. Run 'cargo build --release' first."
    exit 1
fi

# Egypt with auto-doubling
egypt_doubling() {
    local total=0
    for p in 64 128 256 512 1024 2048 4096 8192; do
        local start=$(date +%s%N)
        $EGYPT "$D sqrt" 1 --pell -p $p 2>&1 | grep -q "Fundamental"
        local found=$?
        local end=$(date +%s%N)
        local elapsed=$(( (end - start) / 1000000 ))
        total=$((total + elapsed))
        if [ $found -eq 0 ]; then
            echo $total
            return 0
        fi
    done
    echo "FAILED"
    return 1
}

# PARI/GP benchmark (quadunit)
pari_bench() {
    local start=$(date +%s%N)
    echo "quadunit(4*$D)" | gp -q > /dev/null 2>&1
    local end=$(date +%s%N)
    echo $(( (end - start) / 1000000 ))
}

# PARI/GP regulator only (same complexity)
pari_regulator() {
    local start=$(date +%s%N)
    echo "quadregulator(4*$D)" | gp -q > /dev/null 2>&1
    local end=$(date +%s%N)
    echo $(( (end - start) / 1000000 ))
}

# Run benchmarks
echo "--- Egypt (with auto-doubling) ---"
egypt_times=()
for i in $(seq 1 $RUNS); do
    t=$(egypt_doubling)
    egypt_times+=($t)
    echo "Run $i: ${t}ms"
done

# Calculate Egypt average
egypt_sum=0
for t in "${egypt_times[@]}"; do
    egypt_sum=$((egypt_sum + t))
done
egypt_avg=$((egypt_sum / RUNS))
echo "Egypt average: ${egypt_avg}ms"
echo ""

if [ $PARI_AVAILABLE -eq 1 ]; then
    echo "--- PARI/GP ---"
    pari_times=()
    for i in $(seq 1 $RUNS); do
        t=$(pari_bench)
        pari_times+=($t)
        echo "Run $i: ${t}ms"
    done

    # Calculate PARI average
    pari_sum=0
    for t in "${pari_times[@]}"; do
        pari_sum=$((pari_sum + t))
    done
    pari_avg=$((pari_sum / RUNS))
    echo "PARI/GP average: ${pari_avg}ms"
    echo ""

    # Summary
    echo "=== Summary ==="
    echo "Egypt (auto-doubling): ${egypt_avg}ms"
    echo "PARI/GP:              ${pari_avg}ms"
    speedup=$(echo "scale=1; $pari_avg / $egypt_avg" | bc)
    echo "Speedup: ${speedup}x"
fi
