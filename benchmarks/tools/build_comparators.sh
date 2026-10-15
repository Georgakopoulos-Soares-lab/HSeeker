#!/bin/bash
# Build the comparator tools for the HSeeker tool comparison (2026-10-11).
# Installs into $WORK/tools (TACC /scratch is purged). Verified on Lonestar6, Rocky Linux 8.10.
set -euo pipefail
T=${TOOLS_DIR:-$WORK/tools}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$T" && cd "$T"

# Triplexator 1.3.2 (Buske et al. 2012)
[ -d triplexator ] || git clone https://github.com/Gurado/triplexator.git
(cd triplexator && git checkout -q 4505bba && git apply --check "$HERE/triplexator_seqan_guessFileFormat.patch" 2>/dev/null \
   && git apply "$HERE/triplexator_seqan_guessFileFormat.patch" || true)
module load gcc/9.4.0 cmake/4.1.1
(cd triplexator && rm -rf build_release && mkdir build_release && cd build_release \
   && cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DCMAKE_CXX_FLAGS="-std=gnu++98" \
   && make -j8)
# The patch changes only the return type of the never-instantiated seqan::guessFileFormat to a
# const reference; -std=gnu++98 restores the pre-C++11 make_pair overloads the code relies on.
(cd triplexator/demos && sh smoketest_triplexator.sh ../bin/triplexator)   # expect 14 passed, 0 failed

# NeSSie (Berselli et al. 2018)
[ -d nessie ] || git clone https://github.com/B3rse/nessie.git
(cd nessie && git checkout -q dbe6cb2 && make)

# non-B_gfa / nBMST (Cer et al. 2013)
[ -d non-B_gfa ] || git clone https://github.com/abcsFrederick/non-B_gfa.git
(cd non-B_gfa && git checkout -q a891b59 && make)
# Self-test: the 14 files of test_files.tar are reproduced exactly (run gfa with short relative paths).
