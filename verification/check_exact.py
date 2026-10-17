"""Exact certificate checks for the rank-one Li--Jiang problem.

Python 3 + SymPy. No numerical eigensolver or SDP solver is used here.
The analytic all-parameter proof is in note.tex; finite tests are labeled.
"""
import json
from pathlib import Path
import sympy as s


def zero(v):
    vals = list(v) if isinstance(v, s.MatrixBase) else [v]
    assert all(s.cancel(x) == 0 for x in vals)


def parameters(d, q):
    d, q = s.sympify(d), s.sympify(q)
    m = d*d-1
    p = m*(1-q*q)/(m+q*q)
    l1 = 1-p/(d*d-2)
    l4 = (m+p)/(d*d)
    l2 = (1-l1)*(d*d*l1-1)/(d*d*l4)
    l3 = 2*(1-l1)**2/(d*d*l4)
    A = 1/d-d*l3
    C = 2*A/d*(1-(d*d-2)*q/(m-q*q))
    return tuple(map(s.factor, (p,l1,l2,l3,l4,A,C)))


def rational_channel(d, q):
    p,l1,l2,l3,l4,A,C = parameters(d,q)
    v = s.zeros(d*d,1)
    for j in range(d):
        v[d*j+j] = 1
    H = s.eye(d*d)+(q-1)*v*v.T/d
    aa = d*l4/(d*d-1)
    Ds = [s.kronecker_product(s.eye(d),s.eye(d)[i,:]) for i in range(d)]
    L = [s.eye(d)[:,i]*v.T-q*aa*Ds[i]*H for i in range(d)]
    groups = [(d*(l1-l3)*aa, [D*H for D in Ds]), (l2/l4,L)]
    # The depolarizing term has coefficient d*l3 and all d x d^2 matrix units.
    complete = d*d*l3*s.eye(d*d)
    for weight, mats in groups:
        for M in mats:
            complete += weight*M.T*M
    zero(complete-s.eye(d*d))
    return groups, (p,l1,l2,l3,l4,A,C)


def scalar_identities():
    d,A,C,mu = s.symbols('d A C mu', positive=True)
    n=d-1
    W=s.diag(1,n)
    H=s.Matrix([[A-n*C/d,-n*(d-2)*C/(2*d)],
                [-n*(d-2)*C/(2*d),n*n*C/d]])
    h=s.Matrix([n*(A-C/2),-n*n*C/2])
    S=mu*W-H
    delta=C*(A/d-C/4)
    R=A*A-A*C+d*C*C/4
    D=mu*mu-A*mu+n*delta
    zero(S.det()-n*D)
    x=n*((A-C/2)*mu-n*delta)/D
    z=-n*C*(mu-(2*A/d-C/2))/(2*D)
    w=s.Matrix([x,z])
    zero(S*w-h)
    eta=n*n*(R*mu-n*A*delta)/D
    zero((h.T*w)[0]-eta)
    polynomial=s.expand(D*D-n*n*(R*mu*mu-2*n*A*delta*mu+n*d*delta*delta))
    zero((1-(w.T*W*w)[0])*D**2-polynomial)
    u,v=s.symbols('u v',real=True)
    ww=s.Matrix([u,v])
    gap=A*n*n+mu+eta-(A*n*n+(ww.T*H*ww)[0]+2*(h.T*ww)[0])
    zero(gap-((ww-w).T*S*(ww-w))[0]-mu*(1-(ww.T*W*ww)[0]))
    Z=s.Matrix.vstack(s.Matrix([[eta,-h[0],-h[1]]]),
                      s.Matrix.hstack(-h,S))
    zero(Z*s.Matrix([1,x,z]))
    # Schur factorization, valid wherever S is invertible.
    B=s.Matrix.hstack(-w,s.eye(2))
    zero(Z-B.T*S*B)
    # Rationalized positivity formula, checked symbolically in d and q.
    q=s.symbols('q',positive=True)
    p,l1,l2,l3,l4,AA,CC=parameters(d,q)
    C_original=2*AA/d-2*(l1-l2-l3)*q*(d*d-1+p)/(d*d*(d*d-1))
    zero(CC-C_original)
    # Exceptional singular-value invariants, verified in a two-dimensional span.
    a,b,c,k=s.symbols('a b c k',real=True)
    # Choose alpha=a e0, beta=(u+iv)e0+t e1, and compute Gram's 2x2 block.
    u,v,t=s.symbols('u v t',real=True)
    al=s.Matrix([a,0]); be=s.Matrix([u+s.I*v,t])
    Gram=a*a*s.eye(2)+s.conjugate(be)*be.T+(s.conjugate(al)*be.T+s.conjugate(be)*al.T)/s.sqrt(d)
    bb=u*u+v*v+t*t; cc=a*u; kk=a*a*(u*u+v*v)
    zero(s.trace(Gram)-(2*a*a+bb+2*cc/s.sqrt(d)))
    zero(s.expand(Gram.det())-(a**4+a*a*bb*(1-1/d)+2*a*a*cc/s.sqrt(d)+kk/d))
    return {'symbolic_scalar_and_schur_identities':True,'singular_value_invariants':True}


def finite_case(d,q):
    groups,params=rational_channel(d,q)
    p,l1,l2,l3,l4,A,C=params
    n=d-1
    # Complete Hermitian quadratic matrix on all d^3 complex encoder entries.
    K=l3/d*s.eye(d**3)
    for weight,mats in groups:
        for M in mats:
            g=s.Matrix(list(M.T))
            K += weight/d**2*g*g.T
    rrows=[];zrows=[]
    for i in range(d):
        rr=s.zeros(1,d**3);zz=s.zeros(1,d**3)
        for a in range(d):
            rr[0,(a*d+i)*d+a]=1
            zz[0,(a*d+a)*d+i]=1  # sqrt(d)*z_i
        rrows.append(rr);zrows.append(zz)
    expected=l3/d*s.eye(d**3)
    for r,z in zip(rrows,zrows):
        expected += (A*r.T*r-C/2*(r.T*z+z.T*r)+C/d*z.T*z)/d
    zero(K-expected)
    # Algebraic optimizer and exact quartic-field certificates.
    mu=s.symbols('mu');delta=C*(A/d-C/4);R=A*A-A*C+d*C*C/4
    D=mu*mu-A*mu+n*delta
    poly=s.Poly(D*D-n*n*(R*mu*mu-2*n*A*delta*mu+n*d*delta*delta),mu).clear_denoms()[1].primitive()[1]
    intervals=s.polys.polytools.intervals(poly,eps=s.Rational(1,10**8))
    relevant=[ab for ab,multiplicity in intervals if ab[0]>A and ab[1]<d*A]
    assert len(relevant)==1
    x=n*((A-C/2)*mu-n*delta)/D
    z=-n*C*(mu-(2*A/d-C/2))/(2*D)
    numerator=s.cancel(x*x+n*z*z-1).as_numer_denom()[0]
    zero(s.rem(numerator,poly.as_expr(),mu))
    # Direct decoder Choi certificate at an independent rational family witness.
    u=s.Rational(1,4)
    xx=(1-n*u*u)/(1+n*u*u);zz=-2*u/(1+n*u*u)
    V=s.zeros(d*d,d)
    V[0,0]=xx
    for j in range(1,d):V[d*j+j,0]=zz;V[d*j,j]=1
    zero(V.T*V-s.eye(d))
    KD=l3/d*s.eye(d*d)  # isometric depolarizing contribution
    for weight,mats in groups:
        for M in mats:
            L=M*V;g=s.Matrix(list(L))
            KD += weight/d**2*g*g.T
    alpha=d*l3;beta=d*A;cc=d*d*C/2;t=-zz;ss=xx-n*t
    a0=alpha+beta*xx*xx-cc*(2*ss*xx/d-2*ss*ss/d**2)
    b0=alpha+beta*t*t+cc*(2*ss*t/d+2*ss*ss/d**2)
    g=beta*xx-cc*ss/d
    y0=(a0+n*g)/d**2;y1=(g+n*beta+alpha)/d**2
    Y=s.diag(y0,*([y1]*n));slack=s.kronecker_product(Y,s.eye(d))-KD
    iv=s.Matrix([int(i==j) for i in range(d) for j in range(d)])
    zero(slack*iv)
    zero(s.trace(Y)-(iv.T*KD*iv)[0])
    assert g>=0 and a0>=alpha and b0<=alpha+beta/n
    vals=slack.eigenvals()
    assert all(v.is_nonnegative for v in vals)
    # Degree-two relaxation obstruction using the actual unrestricted encoder.
    # C_i=a H D_i^T. Its summed fidelity equals lambda_1 and summed C_i^T C_i=I.
    m=d*d-1;aa=d*l4/m
    vv=s.zeros(d*d,1)
    for j in range(d):vv[d*j+j]=1
    HH=s.eye(d*d)+(q-1)*vv*vv.T/d
    sum_f=0;sum_tp=s.zeros(d)
    for i in range(d):
        Ci=HH*s.kronecker_product(s.eye(d),s.eye(d)[:,i])
        vec=s.Matrix(list(Ci))
        sum_f += aa*(vec.T*K*vec)[0]
        sum_tp += aa*Ci.T*Ci
    zero(sum_tp-s.eye(d));zero(sum_f-l1)
    if d==3 and q==s.Rational(1,2):
        T=s.symbols('T')
        scaled=s.Poly(poly.as_expr().subs(mu,T/1617),T).clear_denoms()[1].primitive()[1]
        assert scaled.as_expr()==T**4-1054*T**3-440147*T**2+155349060*T-12068784500
        assert scaled.count_roots(1305,1306)==1
        HH=s.Matrix([[A-n*C/d,-n*(d-2)*C/(2*d)],
                     [-n*(d-2)*C/(2*d),n*n*C/d]])
        hh=s.Matrix([n*(A-C/2),-n*n*C/2])
        certificate={'d':3,'p':'8/11','q':'1/2','A':str(A),'C0':str(C),
                     'lambda3':str(l3),'mu_equals':'T/1617',
                     'T_polynomial':str(scaled.as_expr()),'isolating_interval':['1305','1306'],
                     'root_count_in_interval':1,
                     'W':[['1','0'],['0',str(n)]],
                     'H':[[str(v) for v in HH.row(j)] for j in range(2)],
                     'h':[str(v) for v in hh],
                     'x_of_mu':str(s.factor(x)),'z_of_mu':str(s.factor(z)),
                     'eta_of_mu':str(s.factor(n*n*(R*mu-n*A*delta)/D)),
                     'fidelity_of_mu':str(s.factor(l3+(A*n*n+mu+n*n*(R*mu-n*A*delta)/D)/d)),
                     'PSD_and_complementarity':'S=mu W-H; Z=[h^T S^-1 h,-h^T;-h,S]; see note.tex'}
        Path(__file__).with_name('qutrit-certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
    return {'d':d,'q':str(q),'p':str(p),'encoder_matrix_dimension':d**3,
            'trace_preservation':True,'full_complex_quadratic_matrix':True,
            'optimizer_quartic':str(poly.as_expr()),
            'root_interval':[str(z) for z in relevant[0]],
            'quartic_normalization':True,'decoder_dual_psd_complementarity':True,
            'quadratic_relaxation_obstruction':True}


def main():
    result={'sympy':s.__version__,'all_parameter_identities':scalar_identities(),
            'finite_exact_instances':[],'lean_run':False}
    for d in (2,3,4):
        for q in (s.Rational(1,4),s.Rational(1,2),s.Rational(3,4)):
            row=finite_case(d,q)
            result['finite_exact_instances'].append(row)
            print('PASS exact matrix/certificate instance',d,q,flush=True)
    Path(__file__).with_name('exact-results.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS symbolic identities and all 9 exact instances; no Lean run.')


if __name__=='__main__':main()
