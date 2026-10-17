"""Exact all-parameter identities for the fixed-decoder theorem.
Uses t=q_parallel/q_perp, p=3(1-t^2)/(3+t^2), 0<t<1.
No sampling is used in these identity checks.
"""
import sympy as s
import json
from pathlib import Path

t=s.symbols('t',positive=True)
p=3*(1-t*t)/(3+t*t)
a=s.sqrt(2/(t*t+3)); b=t*a
l1=1-p/2; l2=p*(3-2*p)/(2*(3+p)); l3=p*p/(2*(3+p)); l4=(3+p)/4; l5=(1-p)/4
simp=lambda X:s.simplify(s.expand(X))
zero=lambda X: all(simp(z)==0 for z in X) if isinstance(X,s.MatrixBase) else simp(X)==0
I=s.eye(2); J=s.eye(4); psi=s.Matrix([1,0,0,1])/s.sqrt(2)
D=[s.kronecker_product(I,I[j,:]) for j in range(2)]
Q=a*J+(b-a)*psi*psi.T
C=[Q*d.T for d in D]
G=[psi*I[j,:]-b/s.sqrt(2)*C[j] for j in range(2)]
B=[g/(s.sqrt(s.Rational(3,2))*a) for g in G]
E=[(J[:,2*i+j]*psi.T-b/s.sqrt(2)*D[j].T*D[i]*Q)/(s.sqrt(s.Rational(3,2))*a) for i in range(2) for j in range(2)]
groups=[(l3,[J[:,i]*J[j,:] for i in range(4) for j in range(4)]), (l1-l3,[s.kronecker_product(I,I[:,i]*I[j,:])*Q for i in range(2) for j in range(2)]),(l2,E)]
complete=sum((w*k.T*k for w,ks in groups for k in ks),s.zeros(4))
assert zero(complete-J)
for i in range(2):
 for j in range(2):
  assert zero(s.trace(C[i].T*C[j])-int(i==j))
  assert zero(s.trace(C[i].T*B[j]))
  assert zero(s.trace(B[i].T*B[j])-int(i==j))
# Full arbitrary-complex-encoder quadratic form in row-major order.
K=s.zeros(8)
for w,ks in groups:
 for k in ks:
  for d in D:
   row=s.Matrix(1,8,list((d*k).T))
   K+=w*row.T*row/4
cv=[s.Matrix(8,1,list(c)) for c in C]; bv=[s.Matrix(8,1,list(z)) for z in B]
beta=(l1-l3)/2; delta=l2/2
K2=l3/2*s.eye(8)+beta*sum((v*v.T for v in cv),s.zeros(8))+delta*sum((v*v.T for v in bv),s.zeros(8))
assert zero(K-K2)
# The complex two-column fidelity, using |x|^2+|y|^2=1.
xr,xi,yr,yi=s.symbols('xr xi yr yi',real=True)
x=xr+s.I*xi; y=yr+s.I*yi
V=s.Matrix([[x,0],[0,0],[0,1],[y,0]])
v=s.Matrix(8,1,list(V))
F=(v.conjugate().T*K2*v)[0]
k=(3-2*p)/(3+p); h=(1+p)/2; omega=(3+p)*t/6
A=k*(h+omega); BB=k*(h-omega)
target=(1+A*(xr*xr+xi*xi+2*xr)+BB*(yr*yr+yi*yi-2*yr))/4
# Residual is a multiple of the normalization constraint.
normal=xr*xr+xi*xi+yr*yr+yi*yi-1
res=s.factor(simp(F-target)); quotient=s.factor(res/normal)
assert zero(res-quotient*normal)
assert not (set(quotient.free_symbols)&{xr,xi,yr,yi})
# Correct determinant identity; both real and imaginary phase coefficients matter.
u,vv,X,Y,z=s.symbols('u v X Y z',real=True)
M=u*C[0]+vv*((X+s.I*Y)*B[0]+z*B[1]); gram=M.conjugate().T*M

def reduce_z(expr):
 num,den=s.together(s.expand(expr)).as_numer_denom()
 return s.factor(s.rem(num,z*z-(1-X*X-Y*Y),z)/den)
assert zero(reduce_z(s.trace(gram))-u*u-vv*vv)
D0=(3*t**4*u*u*vv*vv+2*t**4*vv**4+18*t*t*u**4-6*t*t*u*u*vv*vv+18*t*t*vv**4+18*u**4+27*u*u*vv*vv)/(9*(t*t+3)**2)
D1=4*s.sqrt(3)*t*u*vv*(3*(1-t*t)*u*u-(9-t*t)*vv*vv)/(9*(t*t+3)**2)
Dx=u*u*vv*vv*(1-t*t)*(9-t*t)/(3*(t*t+3)**2)
Dy=u*u*vv*vv/3
det=reduce_z(s.det(gram)); assert zero(det-D0-D1*X-Dx*X*X-Dy*Y*Y)
# Secular quartic and global upper-bound identity for the circle problem.
AA,BBB,eta,U,W=s.symbols('A B eta U W',real=True)
secular=s.expand((eta-AA)**2*(eta-BBB)**2-AA**2*(eta-BBB)**2-BBB**2*(eta-AA)**2)
quartic=eta**4-2*(AA+BBB)*eta**3+4*AA*BBB*eta**2-(AA*BBB)**2
assert zero(secular-quartic)
bound=eta+AA**2/(eta-AA)+BBB**2/(eta-BBB)
obj=AA*(U*U+2*U)+BBB*(W*W+2*W)
squares=(eta-AA)*(U-AA/(eta-AA))**2+(eta-BBB)*(W-BBB/(eta-BBB))**2
assert zero(bound-obj-squares-eta*(1-U*U-W*W))
mu=s.symbols('mu')
normalized=s.factor(quartic.subs({AA:A,BBB:BB,eta:k*mu})/k**4)
expect=mu**4-2*(1+p)*mu**3+s.Rational(4,3)*p*(p+2)*mu**2-p*p*(p+2)**2/9
assert zero(normalized-expect)
# Exact p=9/13 consistency with the earlier polynomial on the circle.
assert simp(A.subs(t,s.Rational(1,2)))==s.Rational(105,208)
assert simp(BB.subs(t,s.Rational(1,2)))==s.Rational(49,208)
xx,yy=s.symbols('x y',real=True)
f9=(1+s.Rational(105,208)*(xx*xx+2*xx)+s.Rational(49,208)*(yy*yy-2*yy))/4
earlier=(132*xx*xx+210*xx+76*yy*yy-98*yy+181)/832
assert zero(f9-earlier+s.Rational(27,832)*(xx*xx+yy*yy-1))
# Small-p coefficients are in the source parameter, not in normalized infidelity.
pp=s.symbols('p',real=True)
kk=(3-2*pp)/(3+pp); hh=(1+pp)/2; ww=s.sqrt((3+pp)*(1-pp)/12)
aser=s.series(kk*(hh+ww),pp,0,3).removeO(); bser=s.series(kk*(hh-ww),pp,0,3).removeO()
# v*=p/3+O(p^2), u*=1-p^2/18+O(p^3).
asym=s.series((1+aser*((1-pp*pp/18)**2+2*(1-pp*pp/18))+bser*((pp/3)**2+2*pp/3))/4,pp,0,3)
assert asym.removeO()==1-pp/2-pp*pp/36
# Stationarity fixes the optimal first-order W coefficient, with nonzero Hessian.
WW=s.symbols('W',real=True)
objective=kk*(hh+ww)*(1-WW*WW+2*s.sqrt(1-WW*WW))+kk*(hh-ww)*(WW*WW+2*WW)
stationary=s.diff(objective,WW)
assert zero(stationary.diff(WW).subs({pp:0,WW:0})+4)
assert zero(s.series(stationary.subs(WW,pp/3),pp,0,2).removeO())
# Endpoint formula, using the exact secular root at A=B=1/4.
end_eta=(1+s.sqrt(2))/4
endpoint=(1+end_eta+2*s.Rational(1,16)/(end_eta-s.Rational(1,4)))/4
assert zero(endpoint-(5+2*s.sqrt(2))/16)
out={'sympy':s.__version__,'parameterization':'p=3(1-t^2)/(3+t^2), 0<t<1','kraus_completeness':True,'orthonormality':True,'full_8_by_8_quadratic_form':True,'complex_family_formula_modulo_norm':True,'determinant_identity':True,'determinant_coefficients':{'D0':str(D0),'D1':str(D1),'Dx':str(Dx),'Dy':str(Dy)},'secular_quartic':True,'circle_upper_bound_sum_of_squares':True,'normalized_quartic':True,'p_9_over_13_consistency':True,'asymptotic_series':str(asym),'stationarity_small_p':True,'p_to_one_limit':True,'new_lean_run':False}
print(json.dumps(out,indent=2)); Path(__file__).with_name('symbolic-results.json').write_text(json.dumps(out,indent=2)+'\n')
