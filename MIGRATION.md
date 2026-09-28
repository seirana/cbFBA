# Modernization notes

The original cbFBA MATLAB scripts are preserved under `legacy/`. The maintained implementation lives under `src/+cbfba/`.

## Historical structure

The original workflow mixed several responsibilities:

```text
study script
  -> reconstruct/load one submodel
  -> derive chemical-complex strings
  -> build Y and A
  -> classify complexes
  -> solve cbFBA
  -> solve pFBA
  -> write spreadsheets and MAT files
```

Several functions relied on current-working-directory filenames, mutable model fields, repeated string construction, or implicit assumptions about irreversible reaction ordering.

## Maintained structure

```text
validated model
  -> numerical complex decomposition
  -> verify S = Y*A
  -> build sparse optimization problem
  -> checked primary optimization
  -> near-optimal objective constraint
  -> checked reaction-wise range optimization
  -> matched cbFBA/pFBA comparison
  -> optional persistence
```

## Senior-level engineering changes

### Namespaced API

Maintained functions use the MATLAB package namespace `cbfba`.

This avoids polluting the global MATLAB path with generic function names such as `complexes`, `save_data`, or `A_matrix`.

### Numerical complex construction

The historical implementation encoded complexes as concatenated strings such as:

```text
1*A+2*B
```

and later parsed those strings back into the `Y` matrix.

The maintained implementation builds complexes numerically from stoichiometric columns and uses names only for reporting.

### Reconstruction invariant

After building `Y` and `A`, the maintained implementation checks:

```text
norm(S - Y*A, "fro")
```

against an explicit tolerance.

A decomposition that does not reconstruct the stoichiometric matrix fails immediately.

### Correct absolute-value formulations

The maintained LP introduces explicit non-negative auxiliary variables:

```text
-t <= A*v <= t       cbFBA
-t <= v   <= t       pFBA
```

This makes the objective definition explicit even when signed reaction variables are present.

### Solver-status handling

Every maintained `linprog` call checks the exit flag.

The workflow does not convert optimization failure into a zero flux bound.

### Separation of computation and output

Core optimization functions return structs/tables.

CSV, MAT, and JSON files are written only by `cbfba.writeResults` or the file-driven runner.

### Model loading

`cbfba.loadModelFromMat` either loads an explicitly named variable or automatically selects the only model-like struct in a MAT file.

Ambiguous files fail with a clear error instead of relying on workspace side effects from `load`.

## Scientific behavior

The near-optimal factor remains explicit and supports the historical value `1.00001`.

Historical dataset-specific reconstruction scripts are not silently replaced because they encode study-specific assumptions. They remain under `legacy/` and should be compared against a maintained/prepared model before claiming numerical equivalence.

## Compatibility boundary

The maintained core accepts a prepared model directly.

If a historical study requires:

- COBRA Toolbox `convertToIrreversible`;
- gene-knockout reconstruction;
- strain-specific experimental bounds;
- study-specific carbon-source handling;

perform those steps explicitly before calling `cbfba.runAnalysis` and record them as part of experiment provenance.
