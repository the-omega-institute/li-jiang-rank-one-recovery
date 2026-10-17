# Earlier qubit verification

These scripts, CSV grid, and JSON outputs accompany `../previous-note.tex`
and `../previous-note.pdf`. The mathematical status in that earlier note is
historical: the current `../note.pdf` proves joint qubit optimality for all
`0 < p < 1`. The 100-point, 800-start search is an archived experiment.

Install the version-pinned `requirements.txt` here or the root requirements.
Run from the repository root:

```sh
.venv/bin/python previous-verification/symbolic.py
.venv/bin/python previous-verification/validate.py
(cd previous-verification && ../.venv/bin/python counterexample.py)
```

The symbolic checker verifies all-parameter identities with
`p = 3(1-t^2)/(3+t^2), 0 < t < 1`. The validator reconstructs all 24 physical
noise Kraus operators and checks the saved witnesses, both objective
vectorization conventions, covariance, and 4,000 random encoders/contractions.
The random tests provide consistency evidence. `counterexample.py` is an exact
fixed-parameter comparison at `p = 9/13` and writes `verification.json` in its
current working directory.

To rerun the costly grid and overwrite its saved outputs:

```sh
.venv/bin/python previous-verification/numerics.py --starts 8 --outer 80
```

The grid has 100 noise points with eight starts each. `joint_best` includes
the feasible family baseline; `searched_best` records the searched candidate,
which can be slightly below that baseline due to solver error. Approximate
Choi matrices are made positive and trace-preserving before recording their
fidelities. `optimal_inaccurate` solver statuses are permitted, and the
feasible-primal/adjusted-dual gap is recorded. These floating-point bounds
are not interval certificates.

For the decoder, `vec_row(D)` stacks the rows of its `2 x 4` Kraus matrix,
and its Choi matrix uses output tensor input. The decoder cost uses the
transpose and conjugation implemented in `model.py` and checked in `validate.py`.
The current qutrit scripts use their explicitly documented input-first convention.

`model.py` implements the literal source channel, exact-formula numerical
evaluation, objectives, and polar normalization. `decoder-grid.csv`,
`decoder-runs.json`, `decoder-witnesses.json`, and `numerics-summary.json`
record the original experiment, seeds, tolerances, and separate real/imaginary
matrix entries. The CSV `p` field is authoritative; witness dictionary keys are
rounded for lookup. `environment.json` records the executed package versions.

No new Lean result is included. For current proofs and evidence boundaries,
see the root README and `../verification/README.md`.
