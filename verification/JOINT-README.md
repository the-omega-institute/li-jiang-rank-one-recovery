# Exact rank-one encoding and recovery

The note proves F_joint(p,d)=f_d(p) for every integer d >= 2 and 0 < p < 1,
including arbitrary complex rank-one CPTNI encoders and arbitrary CPTNI decoders.
The fixed-decoder quartic formula is consequently the joint optimum.
The proof applies Li--Jiang's rank-one replacement theorem to square channels,
then certifies an optimal unitary recovery through an explicit recovery SDP dual.
The accompanying algebra uses exact rational arithmetic and an isolated real
quartic root. These results are mathematical proofs with exact symbolic checks.
No Lean/lake validation was run for them.

## Exact checks

Install Python 3 and SymPy 1.14.0, then run from the package root:

    python3 -m pip install -r verification/requirements-exact.txt
    python3 verification/check_exact.py
    python3 verification/check_joint.py

check_exact.py checks the inherited all-d fixed-decoder algebra and nine finite
full matrix instances, with d in {2,3,4}, q in {1/4,1/2,3/4}. It checks trace
preservation, the full complex quadratic form, root isolation, exact normalization,
decoder slack positivity and complementary slackness, and the quadratic moment
obstruction. These nine instances support the analytic proof and are not used to
extrapolate to all parameters.

check_joint.py checks the new qutrit product-fidelity gap, all-d product lower
bound, recovery threshold identities and positivity polynomial for d >= 4.
At d=3,p=8/11 it verifies exact root isolation, normalization modulo the quartic,
stationarity, Schur factorization, positive definiteness, decoder equality,
and the rational joint-fidelity enclosure (0.892249,0.892250).
The exact qutrit-joint-certificate.json specifies the algebraic value without
requiring numerical reconstruction.

The commands regenerate the following outputs, omitted from this archive.
Their expected SHA256 values, with the pinned SymPy version, are listed below.
Remove colons to recover each conventional hexadecimal digest:

    99b9d530:9da29143:4abe1525:9c494fc5:9864724b:66401c9b:a752157e:1f9e32cf  verification/exact-results.json
    b94de860:e02f7613:bc7c4129:d0ff6979:c0948f1c:2e967e37:92bd5a8e:40af9c64  verification/qutrit-certificate.json
    e3b919d3:8a9ecb73:21c0b9eb:7b79ea4b:f7179f7b:51d10ce0:0ba47978:ee55525c  verification/joint-results.json

## Independent implementation

The independent referee implemented independent/check_joint_referee.py without
importing either author checker. Its exact part reconstructs the noise action,
qutrit quartic, family decoder actions and all-dimensional identities separately.
Its numerical part checks the Choi convention, stationarity and SDP slack on
nine seeded complex channels in dimensions 2, 3 and 4. To run it:

    python3 -m pip install -r verification/requirements-referee.txt
    python3 verification/independent/check_joint_referee.py

The numerical checks use floating-point optimization. Two runs emitted solver
accuracy warnings in the recorded verification; the maximum objective difference
was below 5e-9. These checks validate conventions and examples, while the universal
result rests on the proof. No global numerical search is asserted.

## Mathematical source and scope

note.tex contains the complete argument, by Haobo Ma and Wenlin Zhang, dated
20 October 2026. Compile it twice with pdflatex. The delivered PDF has nine pages.
The note cites and distinguishes the channel construction and replacement bound
of Li--Jiang (arXiv:2609.00778v1), the qubit quasi-inverse theorem of Karimipour,
Benatti and Floreanini (arXiv:1909.06118v5), and the higher-dimensional context of
Shahbeigi et al. (arXiv:2104.06062v3).
The sharp dimension-four qubit small-error coefficient remains between 1/9
and 1/4. No sharp coefficient theorem is claimed.
