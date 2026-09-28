# Data and model provenance

This directory contains historical model and experimental-data snapshots used by the original study scripts.

The repository currently includes files such as E. coli and yeast MAT models, carbon-source tables, strain information, and experimental measurements.

## Reproducible use

Before reporting new results, record for every source file:

- biological/model source;
- database or reconstruction version;
- publication/reference;
- retrieval date when applicable;
- SHA-256 checksum;
- preprocessing performed before the maintained cbFBA analysis.

## Prepared models

The maintained `cbfba` package is designed to receive one **prepared model** with its final constraints.

For study-specific workflows, a strong reproducibility pattern is:

```text
historical/raw model + experimental condition
        |
        v
explicit reconstruction/preprocessing
        |
        v
saved prepared model.mat
        |
        v
cbfba.runAnalysis
```

This avoids mixing condition reconstruction with the cbFBA/pFBA optimization layer.

## External dependencies

Some historical reconstruction functions require COBRA Toolbox functions such as `convertToIrreversible`.

Record the exact COBRA Toolbox release/commit when reproducing those steps.
