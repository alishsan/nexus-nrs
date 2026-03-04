#!/bin/bash
# Simple script to run all benchmark calculations
# Usage: ./paper/run_benchmarks.sh

cd "$(dirname "$0")/.." || exit 1

echo "Running Nexus-NRS benchmark calculations..."
echo ""

clj -M -e "(load-file \"paper/run_benchmarks.clj\")" -e "(paper.run-benchmarks/run-all-benchmarks)"

echo ""
echo "Done! Check paper/plots/ for generated plots."
