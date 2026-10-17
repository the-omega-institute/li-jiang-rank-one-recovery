# Rank-one encoding and recovery for the Li--Jiang channel

Mathematical proofs, reproducible verification scripts, exact certificates, and
saved numerical witnesses for the noise family in Bikun Li and Liang Jiang,
[*High-Rank Encoding Can Improve Approximate Quantum Error Correction*,
arXiv:2609.00778v1](https://arxiv.org/abs/2609.00778v1).
The accompanying notes are by Haobo Ma and Wenlin Zhang, dated 17 October 2026.

## What is verified

- **Joint qubit optimality:** for logical dimension `d = 2` and every
  `0 < p < 1`, the decoder `tr_S` can be chosen in a globally optimal pair
  with a rank-one encoder, even when both encoder and decoder may vary over
  completely positive trace-nonincreasing maps. The proof reduces the recovery
  to a qubit channel and uses the unitary recovery theorem of Karimipour,
  Benatti, and Floreanini.
- **Fixed-decoder optimum in every dimension:** for every integer `d >= 2`
  and `0 < p < 1`, the full complex rank-one optimum with `tr_S` fixed is
  attained by an explicit two-amplitude family. A specified quartic root
  determines the optimum. Exact primal and dual SDP certificates and a
  positive quadratic gap identity certify the family optimum after the
  analytic reduction from arbitrary complex encoders.
- **Decoder certificates:** an explicit SDP dual certificate proves that
  `tr_S` is an optimal decoder for every encoder in the real attaining family.
  An exact qutrit instance at `d = 3, p = 8/11` is provided as algebraic data.

**Joint global optimality for `d >= 3` remains open.** Separate best-response
certificates do not settle that question. The notes contain the all-parameter
proofs; the scripts check symbolic identities and finite exact matrix instances.
Numerical searches and floating-point witnesses supply consistency evidence.
No Lean formalization or rigorous interval certificate is claimed.

Read [note.pdf](note.pdf) or its [LaTeX source](note.tex) for the complete
argument, and [verification/README.md](verification/README.md) for the scope of
each checker and certificate.

## Run the exact checks

Run from the repository root. The recorded environment uses Python 3.14.4;
use Python 3.14 to reproduce the pinned package environment.

```sh
python3.14 -m venv .venv
.venv/bin/python -m pip install -r verification/requirements-exact.txt
.venv/bin/python verification/check_exact.py
.venv/bin/python verification/independent/check_referee.py
.venv/bin/python previous-verification/symbolic.py
(cd previous-verification && ../.venv/bin/python counterexample.py)
```

The main checker prints nine `PASS` instance lines and a final `PASS` line.
The independent checker prints nine `OK` lines and imports no author checker.
The earlier symbolic and fixed-parameter scripts print JSON with their exact
checks. These commands regenerate the corresponding result JSON files.

## Validate saved numerical witnesses

All direct Python dependencies are version-pinned in [requirements.txt](requirements.txt)
and the requirement files it includes: SymPy 1.14.0, NumPy 2.5.2,
SciPy 1.18.1, CVXPY 1.9.3, and Clarabel 0.11.1.

```sh
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python verification/validate_qutrit.py
.venv/bin/python previous-verification/validate.py
```

These reconstruct the original physical noise Kraus operators and recompute
the saved feasible witnesses: twelve qutrit runs and a 100-point qubit grid
with eight starts per point. They check floating-point fidelity and feasibility
residuals; they do not prove global upper bounds.

To rerun the optional searches, which overwrite saved search outputs:

```sh
.venv/bin/python verification/search_qutrit.py
.venv/bin/python previous-verification/numerics.py --starts 8 --outer 80
```

Solver warnings and nominal tolerances are documented with the recorded runs.
The costly 800-start qubit search is preserved as an earlier experiment.

## Files and integrity

- `note.tex`, `note.pdf`: current all-parameter proofs and open questions.
- `verification/`: exact checker and certificates, independent implementation,
  qutrit numerical exploration, saved witnesses, and recorded outputs.
- `previous-note.tex`, `previous-note.pdf`: earlier fixed-decoder qubit note.
  Its statement that joint optimality is open records the earlier status;
  the current note proves the joint qubit result.
- `previous-verification/`: earlier scripts, 100-point numerical grid, and
  matrix witnesses; see its README for conventions and reproduction.
- `ORIGIN.txt`: provenance and preservation of verification source bytes.
- `SHA256SUMS`: SHA-256 hashes of all published content files except the
  manifest itself. Check with `shasum -a 256 -c SHA256SUMS`.

To rebuild either PDF, run `pdflatex` twice on its `.tex` source. The notes use
standard LaTeX packages listed in their preambles.

## References and license

The channel, unrestricted optimum, and high-rank obstruction are due to
Li and Jiang, [arXiv:2609.00778v1](https://arxiv.org/abs/2609.00778v1).
The qubit recovery theorem is from V. Karimipour, F. Benatti, and R. Floreanini,
*Quasi inversion of qubit channels*, Physical Review A **101**, 032109 (2020),
[doi:10.1103/PhysRevA.101.032109](https://doi.org/10.1103/PhysRevA.101.032109),
[arXiv:1909.06118v5](https://arxiv.org/abs/1909.06118v5).
See [verification/citation-checks.md](verification/citation-checks.md) for the
bounded source check. No exhaustive priority claim is made.

The Python source code is licensed under the [MIT License](LICENSE).
This code license does not relicense the cited papers or the mathematical notes.
