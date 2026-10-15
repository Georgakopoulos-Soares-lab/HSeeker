#!/bin/bash
# Comparator runtime/memory on hg38 chr1, same node and protocol as HSeeker (2026-10-11).
set -u
R=/scratch/10899/kimopro/kouzine/runtime/comparators
CHR1=/scratch/10899/kimopro/kouzine/dl/hg38/chr1.fa
RS=/home1/10899/kimopro/WORK/miniconda3/envs/r_triplex/bin/Rscript
T=$WORK/tools
cd /scratch/10899/kimopro/HSeeker
lscpu > $R/lscpu.txt
for i in 1 2 3; do
  /usr/bin/time -v $RS benchmarks/triplex_benchmark.R $CHR1 $R/triplex_chr1_$i.csv 16 > $R/triplex_$i.stdout 2> $R/triplex_$i.time
  echo "triplex run $i exit $? $(grep -E 'Elapsed|Maximum resident' $R/triplex_$i.time | tr -s ' ' | tr '\n' ' ')"
done
for i in 1 2 3; do
  mkdir -p $R/tx_$i
  /usr/bin/time -v $T/triplexator/bin/triplexator -ds $CHR1 -rm 1 -p 16 -of 0 -od $R/tx_$i -o tts > $R/triplexator_$i.stdout 2> $R/triplexator_$i.time
  echo "Triplexator run $i exit $? $(grep -E 'Elapsed|Maximum resident' $R/triplexator_$i.time | tr -s ' ' | tr '\n' ' ') hits $(grep -vc '^#' $R/tx_$i/tts)"
done
for i in 1 2 3; do
  mkdir -p $R/gfa_$i && cd $R/gfa_$i && ln -sf $CHR1 s.fa
  /usr/bin/time -v $T/non-B_gfa/gfa -skipWGET -skipAPR -skipSTR -skipDR -skipIR -skipGQ -skipZ -seq s.fa -out o > ../gfa_$i.stdout 2> ../gfa_$i.time
  echo "nBMST run $i exit $? $(grep -E 'Elapsed|Maximum resident' ../gfa_$i.time | tr -s ' ' | tr '\n' ' ') MR $(($(wc -l < o_MR.tsv)-1)) triplex $(awk -F'\t' 'NR>1 && $12==1' o_MR.tsv | wc -l)"
  cd /scratch/10899/kimopro/HSeeker
done
mkdir -p $R/nessie && cd $R/nessie
/usr/bin/time -v timeout 7200 bash -c "for k in 16 17 18 19 20; do $T/nessie/nessie -I $CHR1 -O out_k\$k -T -k \$k > /dev/null 2>&1 || exit \$?; done" > ../nessie.stdout 2> ../nessie.time
echo "NeSSie (k 16-20) exit $? $(grep -E 'Elapsed|Maximum resident' ../nessie.time | tr -s ' ' | tr '\n' ' ')"
echo COMPARATORS_DONE
