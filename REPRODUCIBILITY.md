# Reproducibility guide

A cbFBA result is determined by both the source code and the exact constrained metabolic model.

## Minimum experiment record

Archive:

- repository Git commit;
- MATLAB release;
- Optimization Toolbox version;
- LP solver used by `linprog`;
- input model filename;
- MAT variable name;
- model-file SHA-256 checksum;
- source/release of the metabolic reconstruction;
- all preprocessing/reconstruction steps;
- reaction-bound changes;
- irreversible-conversion method when used;
- tolerance factor;
- complex-decomposition tolerance;
- generated comparison CSV;
- metadata JSON;
- result MAT file.

## Study-specific reconstruction

The historical E. coli and yeast scripts encode publication-specific assumptions in `legacy/`.

Those transformations should not be hidden inside a general-purpose library function.

For a reproducible re-analysis:

1. reconstruct one exact model/submodel;
2. save that prepared model;
3. record its checksum;
4. run the maintained cbFBA/pFBA comparison on that saved model.

This separates biological-condition reconstruction from the optimization algorithm itself.

## Solver status

The maintained package requires successful `linprog` termination for:

- the primary cbFBA/pFBA objective;
- every reaction minimum;
- every reaction maximum.

Unexpected solver failure raises an error rather than creating a misleading zero range.

## Tolerance factor

The historical scripts used:

```text
1.00001
```

The maintained runner exposes this as an explicit input.

Report it because a larger factor permits a wider near-optimal feasible region and may change the resulting flux ranges.

## Complex decomposition

The maintained code constructs `Y` and `A` numerically and records the reconstruction error.

A run should preserve the complex-decomposition tolerance whenever exact reproducibility matters.

## Publication-era comparison

If maintained and historical results differ, compare:

- the same prepared model;
- the same reaction bounds;
- the same irreversible representation;
- the same near-optimal factor;
- the same solver;
- the historical scripts under `legacy/`;
- the maintained package.

Do not attribute a difference to the cbFBA method before ruling out model-reconstruction or solver differences.
