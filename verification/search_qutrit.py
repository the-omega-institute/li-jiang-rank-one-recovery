"""Floating-point search, NOT a global certificate. Optional NumPy/SciPy/CVXPY.

Reduced decoder SDP, complex isometry power/polar updates. All starts and
stopping rules are recorded. Exact claims use check_exact.py and note.tex.
"""
import json
from pathlib import Path
import numpy as np
import scipy
from scipy.optimize import brentq
import cvxpy as cp


def polar(V):
    U,_,Wh=np.linalg.svd(V,full_matrices=False)
    return U@Wh


def model(d,q):
    m=d*d-1;p=m*(1-q*q)/(m+q*q)
    l1=1-p/(d*d-2);l4=(m+p)/(d*d)
    l2=(1-l1)*(d*d*l1-1)/(d*d*l4);l3=2*(1-l1)**2/(d*d*l4)
    A=1/d-d*l3;C=2*A/d*(1-(d*d-2)*q/(m-q*q))
    v=np.eye(d).reshape(-1);H=np.eye(d*d)+(q-1)*np.outer(v,v)/d
    aa=d*l4/m
    D=[np.kron(np.eye(d),np.eye(d)[i:i+1]) for i in range(d)]
    M=[np.sqrt(d*l3)*np.outer(np.eye(d)[:,i],np.eye(d*d)[:,j])
       for i in range(d) for j in range(d*d)]
    M += [np.sqrt(d*(l1-l3)*aa)*Di@H for Di in D]
    M += [np.sqrt(l2/l4)*(np.outer(np.eye(d)[:,i],v)-q*aa*D[i]@H) for i in range(d)]
    M=np.array(M)
    assert np.linalg.norm(sum(Z.T@Z for Z in M)-np.eye(d*d))<1e-12
    n=d-1;W=np.diag([1,n]);h=np.array([n*(A-C/2),-n*n*C/2])
    H2=np.array([[A-n*C/d,-n*(d-2)*C/(2*d)],[-n*(d-2)*C/(2*d),n*n*C/d]])
    f=lambda mu: np.linalg.solve(mu*W-H2,h)@W@np.linalg.solve(mu*W-H2,h)-1
    mu=brentq(f,A,d*A,xtol=1e-14);x,z=np.linalg.solve(mu*W-H2,h)
    V=np.zeros((d*d,d),complex);V[0,0]=x
    for j in range(1,d):V[d*j+j,0]=z;V[d*j,j]=1
    optimum=l3+(A*n*n+mu+h@np.array([x,z]))/d
    return p,M,V,optimum


def decoder(M,V,d):
    # Input-first vec(D): D.T.flatten(), so trace(D L)=vec(D).T L.flatten().
    L=M@V;g=L.reshape(len(L),d*d)
    K=g.conj().T@g/d**2
    X=cp.Variable((d*d,d*d),hermitian=True)
    constraints=[X>>0]
    constraints += [sum(X[a*d+i,b*d+i] for i in range(d))==int(a==b)
                    for a in range(d) for b in range(d)]
    prob=cp.Problem(cp.Maximize(cp.real(cp.trace(K@X))),constraints)
    prob.solve(solver='CLARABEL',tol_gap_abs=1e-10,tol_feas=1e-10,tol_gap_rel=1e-10,max_iter=150)
    assert prob.status in ('optimal','optimal_inaccurate')
    J=(X.value+X.value.conj().T)/2
    ev,U=np.linalg.eigh(J)
    Ds=[np.sqrt(v)*U[:,k].reshape(d,d).T for k,v in enumerate(ev) if v>1e-13]
    B=sum(D.conj().T@D for D in Ds)
    b,Z=np.linalg.eigh(B);W=(Z/np.sqrt(b))@Z.conj().T
    return np.array([D@W for D in Ds])


def encoder(M,Ds,V,d):
    gs=np.array([(D@L).T.reshape(-1) for D in Ds for L in M])
    K=gs.conj().T@gs/d**2
    value=lambda Q: float(np.real(np.vdot(Q.reshape(-1),K@Q.reshape(-1))))
    old=value(V)
    for _ in range(150):
        U=polar((K@V.reshape(-1)).reshape(d*d,d));new=value(U)
        V=U
        if new-old<1e-13:break
        old=new
    return V


def fidelity(M,V,Ds,d):
    return float(sum(abs(np.trace(D@L@V))**2 for D in Ds for L in M)/d**2)


def pack(Z):
    return {'real':Z.real.tolist(),'imag':Z.imag.tolist()}


def main():
    d=3;q=.5;p,M,V0,benchmark=model(d,q)
    rng=np.random.default_rng(20261017)
    runs=[];witnesses=[]
    for j in range(12):
        V=V0.copy() if j==0 else polar(rng.normal(size=(d*d,d))+1j*rng.normal(size=(d*d,d)))
        old=-1
        for it in range(60):
            Ds=decoder(M,V,d)
            V=encoder(M,Ds,V,d)
            value=fidelity(M,V,Ds,d)
            if it>=3 and abs(value-old)<1e-10:break
            old=value
        Ds=decoder(M,V,d);value=fidelity(M,V,Ds,d)
        row={'start':j,'outer_steps':it+1,'fidelity':value,'excess_over_fixed':value-benchmark,
             'encoder_residual':float(np.linalg.norm(V.conj().T@V-np.eye(d))),
             'decoder_residual':float(np.linalg.norm(sum(D.conj().T@D for D in Ds)-np.eye(d)))}
        print(json.dumps(row),flush=True);runs.append(row)
        witnesses.append({'start':j,'V':pack(V),'decoder_kraus':[pack(D) for D in Ds]})
    result={'status':'floating point search, not a global or interval certificate',
            'd':d,'q':q,'p':p,'fixed_decoder_benchmark':benchmark,'seed':20261017,
            'outer_cap':60,'encoder_inner_cap':150,'outer_increment_threshold':1e-10,
            'encoder_increment_threshold':1e-13,'decoder_solver':'CLARABEL',
            'nominal_decoder_tolerances':1e-10,
            'numpy':np.__version__,'scipy':scipy.__version__,'cvxpy':cp.__version__,'runs':runs}
    Path(__file__).with_name('qutrit-search.json').write_text(json.dumps(result,indent=2)+'\n')
    Path(__file__).with_name('qutrit-witnesses.json').write_text(json.dumps(witnesses)+'\n')


if __name__=='__main__':main()
