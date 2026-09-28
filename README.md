# cbFBA — Complex-Balanced Flux Balance Analysis

cbFBA is a MATLAB research implementation for comparing a complex-based parsimonious flux objective with conventional parsimonious FBA (pFBA).

The repository now separates the publication-era exploratory scripts from a maintained package:

```text
src/+cbfba/       maintained, tested implementation
legacy/           original scripts and functions
```

The original code is preserved for provenance. New analyses should use the maintained package.

## Why the code was modernized

The historical implementation was useful research code, but it mixed model reconstruction, optimization, result writing, global file assumptions, and low-level matrix construction in one workflow.

The maintained implementation follows a senior-level MATLAB design:

- namespaced package API;
- explicit inputs and outputs;
- input and dimension validation;
- no hard-coded working-directory assumptions in core functions;
- numerical complex construction instead of string-concatenation logic;
- sparse optimization matrices;
- checked solver exit flags;
- reproducible near-optimal tolerance handling;
- computation separated from persistence;
- testable solver-independent components;
- machine-readable outputs and provenance metadata;
- publication-era code archived rather than silently rewritten.

## Method overview

For a stoichiometric model (S), cbFBA constructs a complex matrix (Y) and an incidence matrix (A) such that:

```text
S = Y * A
```

The maintained implementation compares two objectives.

### cbFBA

```text
minimize sum(abs(A * v))
subject to S * v = 0
           lb <= v <= ub
```

The auxiliary variables represent absolute complex-flow imbalance terms.

### pFBA baseline

```text
minimize sum(abs(v))
subject to S * v = 0
           lb <= v <= ub
```

After finding the optimum objective, the workflow constrains the solution space to:

```text
objective <= toleranceFactor * optimum
```

and computes reaction-wise feasible flux ranges.

The historical examples used:

```text
toleranceFactor = 1.00001
```

## Maintained architecture

```text
src/+cbfba/
├── validateModel.m
├── buildComplexMatrices.m
├── classifyComplexes.m
├── buildParsimoniousProblem.m
├── netFluxObjective.m
├── optimizeFluxRanges.m
├── compareMethods.m
├── runAnalysis.m
├── loadModelFromMat.m
└── writeResults.m
```

### Responsibilities

**`validateModel`**  
Validates stoichiometric dimensions, reaction/metabolite identifiers, bounds, and optional irreversible-pair metadata.

**`buildComplexMatrices`**  
Builds the numerical complex matrix `Y` and incidence matrix `A`, and verifies that `Y*A` reconstructs `S` within tolerance.

**`classifyComplexes`**  
Implements deterministic structural classification of trivially balanced, non-trivially balanced, and remaining non-balanced complexes.

**`buildParsimoniousProblem`**  
Builds sparse LP matrices for cbFBA or pFBA using explicit absolute-value auxiliary variables.

**`optimizeFluxRanges`**  
Solves the primary parsimonious objective, constrains the near-optimal space, and computes reaction-wise flux ranges with solver-status checks.

**`compareMethods`**  
Runs cbFBA and pFBA under matched settings and returns one comparison table.

**`writeResults`**  
Writes results only after computation completes.

## Requirements

Full optimization requires:

- MATLAB;
- Optimization Toolbox (`linprog`);
- a metabolic model struct with fields:
  - `S`
  - `rxns`
  - `mets`
  - `lb`
  - `ub`

The core method does not require COBRA Toolbox after the model has been prepared. Study-specific reconstruction scripts may still require COBRA Toolbox, for example when converting a reversible model to an irreversible representation.

## Installation

Clone the repository and add the maintained package:

```matlab
addpath("src")
```

## Programmatic use

```matlab
addpath("src")

model = cbfba.loadModelFromMat("data/ExampleModel.mat");

output = cbfba.runAnalysis( ...
    model, ...
    1.00001, ...
    struct());

head(output.comparison)
```

The comparison table reports signed minima/maxima and the absolute endpoint ranges used by the historical workflow.

## File-driven runner

A reproducible runner is provided:

```matlab
addpath("scripts")

output = runCbFBA( ...
    "data/ExampleModel.mat", ...
    "artifacts/example", ...
    1.00001);
```

If a MAT file contains more than one model-like struct, pass the variable name as the fourth argument.

The runner records:

- MATLAB release/version;
- Git commit when available;
- input MAT-file SHA-256;
- tolerance factor;
- objective optima;
- reaction and complex counts.

## Outputs

A file-driven run writes:

```text
<run>_comparison.csv
<run>_complexes.csv
<run>_metadata.json
<run>_results.mat
```

Generated outputs are ignored by Git.

## Historical scripts

The original code is retained under `legacy/`, including:

- E. coli and yeast study drivers;
- historical cbFBA/pFBA implementations;
- complex-construction scripts;
- strain/submodel reconstruction scripts;
- random flux sampling utilities.

These files are preserved to support comparison with publication-era behavior. They are not the maintained API.

See [MIGRATION.md](MIGRATION.md).

## Tests

Run:

```matlab
addpath("src")
results = runtests("tests");
assertSuccess(results)
```

The solver-independent test suite covers:

- model validation;
- complex decomposition and exact `Y*A = S` reconstruction on a toy network;
- deterministic complex classification;
- cbFBA LP construction;
- pFBA LP construction;
- forward/reverse net-flux objectives;
- safe MAT-model loading.

GitHub Actions runs these tests in a clean MATLAB environment and verifies that the maintained entry points resolve.

## Repository layout

```text
.
├── src/+cbfba/          maintained MATLAB package
├── tests/               MATLAB unit tests
├── scripts/             reproducible runners
├── legacy/              original implementation
├── data/                historical model/data snapshots
├── docs/
│   └── cbFBA with an example.docx
├── README.md
├── MIGRATION.md
├── REPRODUCIBILITY.md
└── CONTRIBUTING.md
```

## Scientific interpretation

cbFBA and pFBA are model-based optimization methods. Their outputs depend on:

- model reconstruction;
- reaction bounds;
- environmental/strain constraints;
- reversible-reaction representation;
- solver behavior and tolerances;
- the near-optimal objective factor.

A narrower or different flux range is a computational property of the specified model and optimization formulation. It is not experimental validation of a flux value by itself.

## Reproducibility

For any reported comparison, preserve the exact model file, model-variable name, SHA-256 checksum, MATLAB release, solver/toolbox version, tolerance factor, reconstruction steps, and generated metadata/results.

See [REPRODUCIBILITY.md](REPRODUCIBILITY.md).

## License

No explicit software license file is currently included. Public repository visibility alone does not grant reuse rights.
