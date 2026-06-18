# chr1 HSeeker purity_rmq Runtime Comparison

Settings: 16 worker threads, scoring enabled, minrep=10, maxrep=1000, maxspacer=10, purity=0.90, mismatch=0.10.

- `purity_rmq=False`: 39.79s, 96,729 hits
- `purity_rmq=True`: 10.34s, 96,729 hits
- Speedup (`True` vs `False`): 3.85x
- Lost hits with RMQ: 0
- Extra hits with RMQ: 0
