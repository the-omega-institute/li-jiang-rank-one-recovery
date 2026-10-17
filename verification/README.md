# Verification source

The readable all-parameter proofs are in ../note.tex and ../note.pdf.
The new mathematics is not Lean-verified. Scripts and certificates below make
the algebra and finite exact instances reproducible; they do not replace the
analytic proof by extrapolation.

## Exact checks

Requirements: Python 3.14 and SymPy 1.14.0, pinned in requirements-exact.txt.
Run from the repository root:

    python3 verification/check_exact.py

Expected output: nine PASS instance lines followed by a final PASS line.
The script regenerates exact-results.json and qutrit-certificate.json.
It checks symbolic secular/Schur/gap identities, singular-value invariants,
the full encoder quadratic matrix at nine rational channel instances,
trace preservation, quartic root isolation and normalization, decoder dual PSD
and complementarity at rational family witnesses, and the exact quadratic
relaxation obstruction. Exact matrices are finite checks; the all-d result
rests on the proofs in the note.

qutrit-certificate.json specifies the quartic root, W,H,h, optimizer and fidelity
as rational functions of that root. S=mu W-H and its zero-Schur-complement
3-by-3 slack are given in the note. The decoder dual for the same family is
given explicitly in the note, with a complete analytic PSD proof.

## Optional numerical exploration

The optional search requires the version-pinned NumPy, SciPy, CVXPY and
CLARABEL dependencies in requirements-numerical.txt:

    python3 verification/search_qutrit.py
    python3 verification/validate_qutrit.py

The first command performs twelve seeded starts at d=3,p=8/11; it regenerates
qutrit-search.json and qutrit-witnesses.json. It uses a reduced decoder SDP
and polar encoder updates, with the tolerances/caps recorded in the JSON and
source. Nominal solver tolerances are not rigorous objective-error bounds.
No positive improvement was found; the largest remaining fidelity deficit
was about 6.29e-10. Neither global optimality nor exact feasibility follows.

The second command independently reconstructs the original 99 physical noise
Kraus operators, uses the saved matrices, and recomputes every reported fidelity.
Its scope is floating-point consistency, not rigorous interval certification.
Recorded output is in qutrit-validation.json.

## Independent referee implementation

independent/check_referee.py was written by a separate reviewer without importing
the author checker. It independently constructs the rationalized channel and
checks nine exact full quadratic matrices and decoder certificates. See its
recorded output. Independence here means a separate derivation and
implementation, not external human peer review. Run it from the root with:

    python3 verification/independent/check_referee.py

## Environment and archives

Fresh exact run: Python 3.14, SymPy 1.14.0. Fresh optional numerical run:
NumPy 2.5.2, SciPy 1.18.1, CVXPY 1.9.3 and CLARABEL 0.11.1.
No Lean/lake invocation occurred.

The public package also preserves the earlier qubit verification scripts and
matrix witnesses in previous-verification/. Those archived numerical runs
are archived. Publication checks recompute their saved witnesses without
rerunning that costly search. Only selected recorded outputs are included;
Python bytecode and build files are omitted.
