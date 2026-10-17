"""Li--Jiang v1, Eqs. (18)--(23); explicit Kraus and fidelity conventions."""
import numpy as np
from scipy.optimize import brentq
I2=np.eye(2); I4=np.eye(4)
psi=np.array([1,0,0,1])/np.sqrt(2)
D=[np.kron(I2,I2[j:j+1,:]) for j in range(2)]

def channel(p):
    l1=1-p/2; l2=p*(3-2*p)/(2*(3+p)); l3=p*p/(2*(3+p))
    l4=(3+p)/4; l5=(1-p)/4
    a=np.sqrt((3+p)/6); b=np.sqrt((1-p)/2)
    Q=a*I4+(b-a)*np.outer(psi,psi)
    E=[(np.outer(I4[:,2*i+j],psi)-np.sqrt(l5)*D[j].T@D[i]@Q)/np.sqrt(l4) for i in range(2) for j in range(2)]
    K=[np.sqrt(l3)*np.outer(I4[:,i],I4[:,j]) for i in range(4) for j in range(4)]
    K += [np.sqrt(l1-l3)*np.kron(I2,np.outer(I2[:,i],I2[:,j]))@Q for i in range(2) for j in range(2)]
    K += [np.sqrt(l2)*e for e in E]
    assert np.linalg.norm(sum(k.T@k for k in K)-I4)<2e-14
    return np.array(K)

def family(p):
    k=(3-2*p)/(3+p); h=(1+p)/2; w=np.sqrt((3+p)*(1-p)/12)
    A=k*(h+w); B=k*(h-w)
    eta=brentq(lambda z:(A/(z-A))**2+(B/(z-B))**2-1,A+1e-14,4,xtol=5e-15)
    x=A/(eta-A); y=-B/(eta-B)
    V=np.array([[x,0],[0,0],[0,1],[y,0]],dtype=complex)
    F=(1+eta+A*A/(eta-A)+B*B/(eta-B))/4
    return V,F,eta,A,B

def fidelity(K,V,dec):
    return sum(abs(np.trace(d@k@V))**2 for d in dec for k in K)/4

def decoder_cost(K,V):
    # d=vec_row(D), trace(D L)=d^T vec_row(L^T).
    L=K@V
    g=L.transpose(0,2,1).reshape(-1,8)
    return g.conj().T@g/4

def encoder_cost(K,dec):
    g=np.array([(d@k).T.reshape(8) for d in dec for k in K])
    return g.conj().T@g/4

def polar(V):
    u,_,vh=np.linalg.svd(V,full_matrices=False)
    return u@vh

def kraus_from_choi(J):
    vals,vec=np.linalg.eigh((J+J.conj().T)/2)
    ds=[np.sqrt(max(v,0))*vec[:,j].reshape(2,4) for j,v in enumerate(vals) if v>1e-14]
    S=sum(d.conj().T@d for d in ds)
    ev,U=np.linalg.eigh(S)
    W=(U/np.sqrt(ev))@U.conj().T
    return np.array([d@W for d in ds])
