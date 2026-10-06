# Anonymous review artifact

This self-contained artifact independently recomputes the seven paired Safe Effective Repair contrasts reported in the manuscript from frozen per-run records. It also reproduces the task-template-clustered sensitivity analysis (equal-weight template means, 100,000 cluster-bootstrap draws, exact sign-flip enumeration, and leave-one-template-out ranges). It does not import an experiment runner or the original analyzers.

## Run

```powershell
python -m pip install -r .\requirements.txt
python .\reproduce_all.py
```

The script writes episode-level and cluster-level results under `generated/`, then exits with code 2 if any planned count, discordance, risk difference, clustered interval, or exact sign-flip result differs from the frozen reference.

The clustered analysis uses the canonical task identifiers already present in each formal record; no result-dependent regrouping is performed. The frozen method note and task map are under `sensitivity/`.

ProtocolSuite implementation source is omitted because redistribution is not authorized. The included frozen records are sufficient for the paper's numerical result audit.
