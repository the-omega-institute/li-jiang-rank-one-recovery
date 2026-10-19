"""Exact algebra for joint rank-one optimality; Python 3 + SymPy.

The all-input proof is in note.tex. This script checks its algebraic
inequalities and the exact qutrit certificate, with rational arithmetic.
It does not infer a theorem from sampled channels or numerical optima.
"""
import json
from pathlib import Path
import sympy as s


def zero(expr):
    values = list(expr) if isinstance(expr, s.MatrixBase) else [expr]
    assert all(s.cancel(x) == 0 for x in values)


def interval(expr, symbol, lo, hi):
    """Rational enclosure using interval Horner arithmetic."""
    def mul(x, y):
        p = [a*b for a in x for b in y]
        return min(p), max(p)
    def polynomial(poly):
        bounds = (s.S.Zero, s.S.Zero)
        for c in s.Poly(poly, symbol).all_coeffs():
            bounds = mul(bounds, (lo, hi))
            bounds = bounds[0]+c, bounds[1]+c
        return bounds
    num, den = s.cancel(expr).as_numer_denom()
    nn, dd = polynomial(num), polynomial(den)
    assert dd[0]*dd[1] > 0
    return mul(nn, (1/dd[1], 1/dd[0]))


def main():
    d,q,k,r = s.symbols('d q k r', real=True)
    p=(d*d-1)*(1-q*q)/(d*d-1+q*q)
    l3=2*p*p/((d*d-2)**2*(d*d-1+p))
    A=1/d-d*l3
    C=2*A/d*(1-(d*d-2)*q/(d*d-1-q*q))
    F=s.factor(l3+d*A-(d*d-1)*C/d**2)
    fq=s.factor(F.subs(d,3))
    gap=(440+3008*q+1143*q*q+512*q**3+512*q**3*(1-q))/(2268*(q*q+8))
    zero(fq-s.Rational(3,4)-gap)
    # Every summand in the displayed numerator is nonnegative on [0,1],
    # and the constant 440 and denominator are strictly positive.
    assert gap.subs(q,0)>0 and gap.subs(q,1)>0
    assert 4*s.Rational(25,36)-2*s.Rational(5,6)-1==s.Rational(1,9)
    threshold=(d+2)/(d+3)
    zero(((d+1)*r*r-(d-1)*r-1).subs(r,threshold)-(d+1)/(d+3)**2)
    # The derivative in sqrt(f) is positive above the stated threshold.
    zero((2*(d+1)*threshold-(d-1))-(d*d+4*d+7)/(d+3))
    eps=s.symbols('eps',real=True)
    l3p=2*eps**2/((d*d-2)**2*(d*d-1+eps))
    zero(s.diff(l3p,eps)-2*eps*(2*d*d-2+eps)/((d*d-2)**2*(d*d-1+eps)**2))
    Ap=1/d-d*l3p
    lower=s.factor(l3p+d*Ap-2*(d*d-1)*Ap/d**3)
    m=(d**6-4*d**4+4*d*d-2)/(d**4*(d*d-2))
    zero(lower.subs(eps,1)-m)
    coeff=(d*d-1)*(d*d-2)/d**2
    zero(lower-(m+coeff*(l3p.subs(eps,1)-l3p)))
    N=d**6-2*d**5-8*d**4+4*d**3+12*d*d-2*d-6
    zero(m-threshold-N/(d**4*(d+3)*(d*d-2)))
    positive=s.Poly(s.expand(N.subs(d,k+4)),k)
    assert positive.all_coeffs()==[1,22,192,836,1852,1822,434]
    print('PASS qutrit positive product gap and recovery threshold')
    print('PASS all-d product bound and positive shifted polynomial')

    # Exact field certificate for the globally optimal qutrit instance.
    T=s.symbols('T')
    poly=T**4-1054*T**3-440147*T**2+155349060*T-12068784500
    assert s.Poly(poly,T).count_roots(1305,1306)==1
    ints=s.polys.polytools.intervals(poly,eps=s.Rational(1,10**9))
    lo,hi=next(ab for ab,multiplicity in ints if ab[0]>1305 and ab[1]<1306)
    aa=s.Rational(527,1617); cc=s.Rational(578,4851)
    ll=s.Rational(4,1617);mu=T/1617;n=2
    W=s.diag(1,n)
    H=s.Matrix([[aa-n*cc/3,-n*cc/6],[-n*cc/6,n*n*cc/3]])
    h=s.Matrix([n*(aa-cc/2),-n*n*cc/2])
    S=mu*W-H
    w=S.inv()*h
    norm=s.cancel((w.T*W*w)[0]-1).as_numer_denom()[0]
    zero(s.rem(norm,poly,T))
    zero(S*w-h)
    eta=(h.T*w)[0]
    Z=s.BlockMatrix([[s.Matrix([[eta]]),-h.T],[-h,S]]).as_explicit()
    B=s.Matrix.hstack(-w,s.eye(2))
    zero(Z-B.T*S*B)
    zero(Z*s.Matrix([1,*list(w)]))
    # Strict PSD of S at the isolated root by exact interval signs.
    assert interval(S[0,0],T,lo,hi)[0]>0
    assert interval(S.det(),T,lo,hi)[0]>0
    x,z=w
    assert interval(x,T,lo,hi)[0]>0
    assert interval(z,T,lo,hi)[1]<0
    f=ll+(4*aa+mu+eta)/3
    flo,fhi=interval(f,T,lo,hi)
    assert flo>s.Rational(892249,10**6)
    assert fhi<s.Rational(892250,10**6)
    # Matrix certificate for the optimizer's reduced recovery.
    alpha=3*ll;beta=3*aa;c=9*cc/2;t=-z;v=x-2*t
    a0=alpha+beta*x*x-c*(2*v*x/3-2*v*v/9)
    b0=alpha+beta*t*t+c*(2*v*t/3+2*v*v/9)
    g=beta*x-c*v/3
    blocks=[g,a0+2*g-alpha,g+2*beta+alpha-b0,g+2*beta]
    # Expressions may have large cancellations; independently bound their
    # reduced rational forms over the isolated root interval.
    for value in blocks:
        assert interval(s.factor(value),T,lo,hi)[0]>0
    y0=(a0+2*g)/9;y1=(g+2*beta+alpha)/9
    num=s.cancel(y0+2*y1-f).as_numer_denom()[0]
    zero(s.rem(num,poly,T))
    print('PASS qutrit root isolation, field normalization, PSD and dual equality')
    print('Certified joint fidelity interval: (0.892249, 0.892250)')
    result={'all_parameter_product_inequalities':True,
            'qutrit_root_count_1305_1306':1,
            'qutrit_fidelity_interval':['0.892249','0.892250'],
            'qutrit_exact_field_certificate':True,
            'sympy':s.__version__,'lean_run':False}
    Path(__file__).with_name('joint-results.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    main()
