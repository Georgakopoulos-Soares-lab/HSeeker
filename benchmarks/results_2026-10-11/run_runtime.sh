#!/bin/bash
# Runtime/memory re-measurement at the 2026-10-11 defaults (minrep 8), same protocol as
# benchmarks/benchmark_comparison.md: hseeker CLI, /usr/bin/time -v, scoring on, -purity-rmq.
set -u
export PYTHONPATH=/scratch/10899/kimopro/HSeeker/src
export OMP_NUM_THREADS=1
PY=/scratch/10899/kimopro/kouzine/env/bin/python
R=/scratch/10899/kimopro/kouzine/runtime
cd /scratch/10899/kimopro/HSeeker
lscpu > $R/lscpu.txt; free -m > $R/free.txt; uname -a > $R/uname.txt
$PY benchmarks/parallel_sweep.py --minrep 8 --workers 1,2,4,8,16 --runs 3 --json $R/sweep_greedy_minrep8.json > $R/sweep_greedy.log 2>&1
echo "sweep exit $?"
cli() {  # label fasta strategy run
  /usr/bin/time -v $PY -m hseeker -seq $2 -out $R/out_$1_$4 -minrep 8 -maxrep 1000 -maxspacer 10 \
     -purity 0.90 -mismatch 0.10 -at-threshold 0.80 -overlap-strategy $3 -workers 16 -purity-rmq \
     > $R/$1_$4.stdout 2> $R/$1_$4.time
  echo "$1 run $4 exit $? $(grep -E 'Elapsed|Maximum resident' $R/$1_$4.time | tr -s ' ' | tr '\n' ' ')"
}
for i in 1 2 3; do cli chr1_score benchmarks/data/chr1.fa score $i; done
for i in 1 2 3; do cli genome_greedy /scratch/10899/kimopro/kouzine/dl/hg38/hg38.primary.fa greedy $i; done
echo RUNTIME_DONE
