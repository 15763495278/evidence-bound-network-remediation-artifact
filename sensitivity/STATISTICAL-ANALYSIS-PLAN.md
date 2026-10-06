# Frozen task-template sensitivity analysis

This secondary analysis treats the canonical task template, not the episode, as the resampling unit. It was run only after the F1--F4 records and the seven episode-level contrasts had been frozen.

For each contrast, paired episode differences are averaged within each canonical task template. The reported clustered risk difference is the equal-weight mean of those template means. Uncertainty is summarized with 100,000 nonparametric bootstrap resamples of task templates using NumPy PCG64 seed `20260930`. The two-sided sign-flip p-value is obtained by exact enumeration of all `2^G` sign assignments, where `G` is the number of task templates. Leave-one-template-out ranges are descriptive. These p-values do not replace the planned episode-level tests and are not entered into the original Holm family.

All seven pre-existing contrasts are reported. No contrast, task, or model request was added or removed on the basis of the sensitivity results.
