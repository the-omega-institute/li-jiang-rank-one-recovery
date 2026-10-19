# Exact rank-one encoding and recovery for the Li--Jiang channel

For every logical dimension d >= 2 and every 0 < p < 1, partial trace can
be chosen in a jointly optimal rank-one encoder--decoder pair for the noise
family of Bikun Li and Liang Jiang,
[High-Rank Encoding Can Improve Approximate Quantum Error Correction,
arXiv:2609.00778v1](https://arxiv.org/abs/2609.00778v1).
The optimization allows all complex rank-one CPTNI encoders and arbitrary
CPTNI recoveries. The specified quartic-root formula in the note is therefore
the exact joint optimum, including the qutrit case d=3,p=8/11.

Read [note.pdf](note.pdf) or [note.tex](note.tex), by Haobo Ma and Wenlin Zhang,
dated 20 October 2026. The proof combines the rank-one replacement bound of
Li--Jiang with an explicit sufficient SDP dual criterion for unitary recovery
on square channels. A product encoder supplies the needed fidelity threshold.
The all-parameter result is an analytic proof, with exact symbolic checks;
it is not a new Lean formalization.

## Reproduce the exact checks

The new proof and exact qutrit certificate are documented in
[verification/JOINT-README.md](verification/JOINT-README.md).
Install SymPy 1.14.0 and run from the repository root:

    python3 -m pip install -r verification/requirements-exact.txt
    python3 verification/check_exact.py
    python3 verification/check_joint.py

The first checker verifies the all-dimensional fixed-decoder algebra and nine
finite full-matrix instances. The second verifies the new product-fidelity and
unitary-recovery threshold identities, positive polynomials for all d >= 4,
and exact algebraic certificates and a rational fidelity interval for the qutrit
instance. The finite instances support the analytic proof, without extrapolation.

The concise [qutrit joint certificate](verification/qutrit-joint-certificate.json)
uses the unique real root T in (1305,1306) of

    T^4 - 1054*T^3 - 440147*T^2 + 155349060*T - 12068784500.

The joint fidelity lies in (0.892249,0.892250).

## Independent implementation and archived experiments

The new independent implementation is
[verification/independent/check_joint_referee.py](verification/independent/check_joint_referee.py).
It reconstructs the channel, qutrit quartic, all-d inequalities and family decoder
actions without importing the author checkers. It also tests Choi conventions and
the recovery criterion numerically on nine seeded complex channels. Run it with:

    python3 -m pip install -r verification/requirements-referee.txt
    python3 verification/independent/check_joint_referee.py

These optional floating-point checks do not establish the universal theorem.
The older independent implementation, searches and numerical witnesses remain
in verification/ and previous-verification/ with their original documentation.
They record earlier experiments and are not represented as fresh numerical runs.

previous-note.tex and previous-note.pdf contain the earlier fixed-decoder qubit
note. Its open joint-optimality statement records the earlier scope; note.tex
contains the current complete result. Existing scripts and matrix witnesses for
the blocked original attachment remain available at their established paths.

## Remaining question and references

For maximally mixed qubits in physical dimension four, the sharp universal
small-error rank-gap coefficient still satisfies 1/9 <= c_(2,4) <= 1/4.
The current family has coefficient 1/9; the exact joint solution does not
establish the sharp coefficient over other noise channels.

Primary sources are Li--Jiang (arXiv:2609.00778v1),
Karimipour--Benatti--Floreanini (Phys. Rev. A 101, 032109 (2020),
DOI10.1103/PhysRevA.101.032109, arXiv:1909.06118v5), and
Shahbeigi et al. (J. Phys. A 54, 345301 (2021),
DOI10.1088/1751-8121/ac13db, arXiv:2104.06062v3).
Their specific contributions and the complete proofs appear in the note.
No exhaustive priority claim is made.

To compile the nine-page PDF, run pdflatex twice on note.tex.
SHA256SUMS hashes every published content file except the manifest itself;
check it with `shasum -a 256 -c SHA256SUMS`.
The Python source is licensed under LICENSE (MIT). This code license does not
relicense the cited papers or mathematical notes.
