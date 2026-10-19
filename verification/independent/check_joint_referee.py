"""Independent referee checks, rederived without importing author code."""
import json
from pathlib import Path
import numpy as np
import sympy as s
from scipy.linalg import expm
from scipy.optimize import minimize
import cvxpy as cp

q=s.symbols('q', real=True)
d=3
p=8*(1-q*q)/(8+q*q)
l1=1-p/7;l4=(8+p)/9;l5=(1-p)/9
l2=s.cancel((1-l1)*(9*l1-1)/(9*l4));l3=s.cancel(2*(1-l1)**2/(9*l4))
a2=3*l4/8
O=s.zeros(9,1)
for j in range(3):O[3*j+j]=1
P=O*O.T/3
Qbar=s.eye(9)+(q-1)*P
Ds=[];Es=[]
for i in range(3):
    D=s.zeros(3,9)
    for j in range(3):D[j,3*j+i]=1
    ei=s.eye(3)[:,i]*O.T-q*a2*D*Qbar
    Ds.append(D*Qbar);Es.append(ei)

def ch(X):
    return (3*l3*s.trace(X)*s.eye(3)+3*(l1-l3)*a2*sum((D*X*D.T for D in Ds),s.zeros(3))+l2/l4*sum((E*X*E.T for E in Es),s.zeros(3))).applyfunc(s.cancel)

# Derive baseline fidelity directly from three noise groups.
V=s.zeros(9,3)
for j in range(3):V[3*j,j]=1
baseline=s.cancel(l3+(3*(l1-l3)*a2*sum(s.trace(D*V)**2 for D in Ds)+l2/l4*sum(s.trace(E*V)**2 for E in Es))/9)
claimed=(-128*q**4+256*q**3+711*q*q+752*q+3512)/(567*(q*q+8))
assert s.cancel(baseline-claimed)==0
numerator=440+3008*q+1143*q*q+512*q**3+512*q**3*(1-q)
assert s.cancel(baseline-s.Rational(3,4)-numerator/(2268*(q*q+8)))==0
# TP of reduced map, by coefficient matrix.
assert s.simplify(9*l3*s.eye(9)+3*(l1-l3)*a2*sum((D.T*D for D in Ds),s.zeros(9))+l2/l4*sum((E.T*E for E in Es),s.zeros(9))-s.eye(9))==s.zeros(9)
# Independent fixed-family fidelity at q=1/2 and Lagrange polynomial.
x,z,mu=s.symbols('x z mu',real=True)
V=s.zeros(9,3);V[0,0]=x;V[4,0]=z;V[8,0]=z;V[3,1]=1;V[6,2]=1
f=s.cancel(l3*(x*x+2*z*z+2)/3+(3*(l1-l3)*a2*sum(s.trace(D*V)**2 for D in Ds)+l2/l4*sum(s.trace(E*V)**2 for E in Es))/9)
f=s.factor(f.subs(q,s.Rational(1,2)))
# Independently confirm every matrix-unit action in the decoder-family certificate.
ll=s.Rational(4,1617); AA=s.Rational(527,1617); CC=s.Rational(578,4851)
alpha=3*ll; beta=3*AA; c=9*CC/2; tt=-z; ss=x-2*tt
av=alpha+beta*x*x-c*(2*ss*x/3-2*ss*ss/9)
bv=alpha+beta*tt*tt+c*(2*ss*tt/3+2*ss*ss/9)
gv=beta*x-c*ss/3
for ii in range(3):
    for jj in range(3):
        actual=ch(V[:,ii]*V[:,jj].T).subs(q,s.Rational(1,2))
        if ii==jj==0: expected_action=s.diag(av,bv,bv)
        elif ii==jj: expected_action=alpha*s.eye(3)+beta*s.eye(3)[:,ii]*s.eye(3)[:,jj].T
        else: expected_action=(gv if 0 in (ii,jj) else beta)*s.eye(3)[:,ii]*s.eye(3)[:,jj].T
        for diff in actual-expected_action:
            num=s.cancel(diff).as_numer_denom()[0]
            assert s.rem(num,x*x+2*z*z-1,x)==0
print('EXACT: all nine reduced-channel matrix units match the family decoder certificate.')

# Construct Hessian from actual channel independently.
qc=s.expand(3*f-s.Rational(4,1617)*(x*x+2*z*z+2))
H=s.hessian(qc,(x,z))/2
h=s.Matrix([qc.diff(v).subs({x:0,z:0})/2 for v in (x,z)])
W=s.diag(1,2)
w=(mu*W-H).inv()*h
poly=s.together((w.T*W*w)[0]-1).as_numer_denom()[0]
T=s.symbols('T')
poly=s.Poly(poly.subs(mu,T/s.Integer(1617)),T).monic().as_expr()
expected=T**4-1054*T**3-440147*T**2+155349060*T-12068784500
assert s.expand(poly-expected)==0
intervals=s.polys.polytools.intervals(expected,eps=s.Rational(1,10**12))
isolated=[(ab,m) for ab,m in intervals if ab[0]>1305 and ab[1]<1306]
assert len(isolated)==1 and isolated[0][1]==1
root=float(sum(isolated[0][0])/2)
wv=np.array(w.subs(mu,root/1617),float).ravel()
assert abs(wv@np.diag([1.,2.])@wv-1)<1e-12
fopt=float(f.subs({x:wv[0],z:wv[1]}))
# Exact source theorem constants, monotonic bootstrap endpoint.
assert s.Rational(3,4)-s.Rational(2,9)*(1-s.sqrt(1-s.Rational(3,4)))**2==s.Rational(25,36)
assert 4*s.Rational(25,36)-2*s.Rational(5,6)-1==s.Rational(1,9)

# All-dimensional product bound and bootstrap proposed during referee review.
dsym,k=s.symbols('d k',positive=True)
lam=s.symbols('lam',nonnegative=True)
Asym=1/dsym-dsym*lam
Flower=lam+dsym*Asym-(dsym**2-1)/dsym**2*(2*Asym/dsym)
md=(dsym**6-4*dsym**4+4*dsym**2-2)/(dsym**4*(dsym**2-2))
assert s.cancel(Flower.subs(lam,2/(dsym**2*(dsym**2-2)**2))-md)==0
N=dsym**6-2*dsym**5-8*dsym**4+4*dsym**3+12*dsym**2-2*dsym-6
assert s.cancel(md-(dsym+2)/(dsym+3)-N/(dsym**4*(dsym+3)*(dsym**2-2)))==0
assert s.expand(N.subs(dsym,k+4))==k**6+22*k**5+192*k**4+836*k**3+1852*k**2+1822*k+434
threshold=(dsym+2)/(dsym+3)
assert s.factor((dsym+1)*threshold**2-(dsym-1)*threshold-1)==(dsym+1)/(dsym+3)**2
print('EXACT: all-dimension product lower bound m_d, d>=4 polynomial positivity and unitary-optimum threshold checked.')

print('EXACT: baseline, positive 3/4 gap, TP, q=1/2 quartic and isolating interval all checked.')
print('Qutrit instance: x=%.12f z=%.12f fidelity=%.12f'%(wv[0],wv[1],fopt))

# Generic fixed-decoder algebra: solve the secular system without specializing d.
Dg,Ag,Cg,mg=s.symbols('D A C m', positive=True)
ng=Dg-1
Wg=s.diag(1,ng)
Hg=s.Matrix([[Ag-ng*Cg/Dg,-ng*(Dg-2)*Cg/(2*Dg)],[-ng*(Dg-2)*Cg/(2*Dg),ng*ng*Cg/Dg]])
hg=s.Matrix([ng*(Ag-Cg/2),-ng*ng*Cg/2])
delta=Cg*(Ag/Dg-Cg/4)
Rg=Ag*Ag-Ag*Cg+Dg*Cg*Cg/4
Dmu=mg*mg-Ag*mg+ng*delta
Sg=mg*Wg-Hg
assert s.cancel(Sg.det()-ng*Dmu)==0
wg=s.Matrix([ng*((Ag-Cg/2)*mg-ng*delta)/Dmu,-ng*Cg*(mg-(2*Ag/Dg-Cg/2))/(2*Dmu)])
assert (Sg*wg-hg).applyfunc(s.cancel)==s.zeros(2,1)
eta=ng*ng*(Rg*mg-ng*Ag*delta)/Dmu
assert s.cancel((hg.T*wg)[0]-eta)==0
assert s.cancel(((wg.T*Wg*wg)[0]-1)*Dmu**2+ Dmu**2-ng*ng*(Rg*mg*mg-2*ng*Ag*delta*mg+ng*Dg*delta*delta))==0
print('EXACT: generic-d secular equations, optimizer, quartic and objective correction checked.')

# Numerical complex-channel checks of Choi conventions and sufficient lemma.
rng=np.random.default_rng(78239)
def vec(a):return a.reshape(-1,order='F')
def fidelity(kraus):return sum(abs(np.trace(a))**2 for a in kraus)/len(kraus[0])**2
def random_kraus(dim,rank):
    Z=rng.normal(size=(dim*rank,dim))+1j*rng.normal(size=(dim*rank,dim))
    Q,_=np.linalg.qr(Z)
    return [Q[dim*i:dim*(i+1)] for i in range(rank)]
def choi_adjoint(k):return sum(np.outer(vec(a.conj().T),vec(a.conj().T).conj()) for a in k)/len(k[0])**2

def decoder_sdp(k):
    D=len(k[0]);K=choi_adjoint(k);X=cp.Variable((D*D,D*D),hermitian=True)
    constraints=[X>>0]
    for i in range(D):
        for j in range(D):constraints.append(sum(X[i*D+a,j*D+a] for a in range(D))==(i==j))
    prob=cp.Problem(cp.Maximize(cp.real(cp.trace(K@X))),constraints)
    prob.solve(solver='CLARABEL',tol_gap_abs=1e-9,tol_feas=1e-9,tol_gap_rel=1e-9)
    return float(prob.value)

def Hermitian(v,D):
    H=np.zeros((D,D),complex);H[np.diag_indices(D)]=v[:D]
    idx=D
    for i in range(D):
        for j in range(i+1,D):H[i,j]=v[idx]+1j*v[idx+1];H[j,i]=v[idx]-1j*v[idx+1];idx+=2
    return H

max_gap=0.;max_ann=0.;min_eig=1.; cases=0
for D in (2,3,4):
    for noise in (.01,.05,.15):
        raw=random_kraus(D,3)
        k=[np.sqrt(1-noise)*np.eye(D)]+[np.sqrt(noise)*a for a in raw]
        def objective(v):
            U=expm(1j*Hermitian(v,D));return -fidelity([U@a for a in k])
        sol=minimize(objective,np.zeros(D*D),method='BFGS',options={'gtol':1e-10,'maxiter':800})
        U=expm(1j*Hermitian(sol.x,D));k=[U@a for a in k]
        f0=fidelity(k)
        expression=(D+1)*f0-(D-1)*np.sqrt(f0)-1
        assert expression>0
        B=sum(np.trace(a).conj()*a for a in k)
        # Numerical stationarity checked separately from Hermitization.
        assert np.linalg.norm(B-B.conj().T)<2e-6
        Y=B.T/D**2
        K=choi_adjoint(k)
        slack=np.kron(Y,np.eye(D))-K
        ann=np.linalg.norm(slack@vec(np.eye(D)))
        eig=np.linalg.eigvalsh((slack+slack.conj().T)/2)
        assert ann<5e-7 and eig.min()>-5e-7
        # Independent arbitrary-recovery SDP.
        best=decoder_sdp(k)
        gap=abs(best-f0);assert gap<2e-7
        cases+=1;max_gap=max(max_gap,gap);max_ann=max(max_ann,ann);min_eig=min(min_eig,eig[0])
        # Choi objective equals direct Kraus objective for a complex random recovery.
        rk=random_kraus(D,2)
        X=sum(np.outer(vec(a),vec(a).conj()) for a in rk)
        direct=sum(abs(np.trace(r@a))**2 for r in rk for a in k)/D**2
        assert abs(np.trace(K@X).real-direct)<1e-12
print('NUMERICAL: %d complex channels (d=2,3,4), independent recovery SDP; max difference %.3g, slack-annihilation %.3g, min slack eigenvalue %.3g.'%(cases,max_gap,max_ann,min_eig))
print('All independent referee checks passed.')
