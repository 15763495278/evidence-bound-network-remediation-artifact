# Anonymous review artifact

This self-contained artifact independently recomputes the seven paired Safe Effective Repair contrasts reported in the manuscript from frozen per-run records. It does not import an experiment runner or the original analyzers.

## Run

```powershell
python .\reproduce_all.py
```

The script writes `generated/recomputed_results.json`, `effect_estimates.csv`, and `verification_report.md`, and exits with code 2 if any planned count, discordance, or risk difference differs from the frozen reference.

ProtocolSuite implementation source is omitted because redistribution is not authorized. The included frozen records are sufficient for the paper's numerical result audit.
