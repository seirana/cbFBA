# Contributing

## Scope

Changes to the maintained implementation belong under:

```text
src/+cbfba/
```

Publication-era code under `legacy/` should remain unchanged except for archival/documentation maintenance.

## Development workflow

1. create a branch;
2. keep scientific and engineering changes reviewable;
3. add or update MATLAB unit tests;
4. run the test suite;
5. update README/MIGRATION documentation when behavior changes;
6. open a pull request.

## Tests

```matlab
addpath("src")
results = runtests("tests");
assertSuccess(results)
```

Tests should prefer small synthetic models and should not depend on large third-party toolboxes unless the behavior specifically requires them.

## Scientific changes

Any change that can alter numerical results should document:

- the previous behavior;
- the new behavior;
- the reason for the change;
- whether publication-era output is expected to change;
- how the change was tested.

Avoid silently changing model-specific biological assumptions.

## Style

Prefer:

- package-qualified public functions;
- descriptive names;
- explicit options structs;
- validated inputs;
- sparse matrices for LP construction;
- deterministic ordering;
- checked optimization exit flags;
- structured return values;
- file I/O at orchestration boundaries rather than deep inside algorithms.
