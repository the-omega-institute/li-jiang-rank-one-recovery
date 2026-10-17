import sympy as s,json
from pathlib import Path
x,y=s.symbols('x y', real=True)
e=s.eye(4);I=s.eye(2)
A=s.Matrix([[3,0,0,-1],[0,4,0,0],[0,0,4,0],[-1,0,0,3]])
D=[s.kronecker_product(I,s.eye(2)[j,:]) for j in range(2)]
delta=s.Matrix([[1,0,0,1]])
T=[(13*e[:,2*i+j]*delta-D[j].T*D[i]*A)/2 for i in range(2) for j in range(2)]
R=[s.kronecker_product(I,s.eye(2)[:,i]*s.eye(2)[j,:])*A for i in range(2) for j in range(2)]
M=[e[:,a]*e[b,:] for a in range(4) for b in range(4)]
groups=[(s.Rational(702,10816),M),(s.Rational(245,10816),R),(s.Rational(21,10816),T)]
V=s.Matrix([[x,0],[0,0],[0,1],[y,0]])
polys=[s.expand(w*sum(s.trace(d*m*V)**2 for d in D for m in mats)/4) for w,mats in groups]
f=s.factor(sum(polys));target=(132*x*x+210*x+76*y*y-98*y+181)/832
assert s.expand(f-target)==0
complete=s.zeros(4)
for w,mats in groups:
 for m in mats:complete+=w*m.T*m
assert complete==s.eye(4)
p=s.Rational(9,13);l1=1-p/2;l4=(3+p)/4;l5=(1-p)/4;l2=(1-l1)*(4*l1-1)/(4*l4);l3=2*(1-l1)**2/(4*l4)
assert l3==groups[0][0]
assert (l1-l3)/26==groups[1][0]
assert l2/78==groups[2][0]
# Check the channel on every matrix unit against literal Eq. (23).
Q=A/s.sqrt(26);E=[t/s.sqrt(78) for t in T]
for X in M:
 qxq=Q*X*Q
 trS=sum((d*qxq*d.T for d in D),s.zeros(2))
 literal=l3*s.trace(X)*s.eye(4)+(l1-l3)*s.kronecker_product(trS,I)+l2*sum((t*X*t.T for t in E),s.zeros(4))
 kraus=sum((w*m*X*m.T for w,mats in groups for m in mats),s.zeros(4))
 assert s.simplify(literal-kraus)==s.zeros(4)
fp=s.simplify(f.subs({x:3/s.sqrt(10),y:-1/s.sqrt(10)}))
fw=f.subs({x:s.Rational(24,25),y:-s.Rational(7,25)})
gain=s.simplify(fw-fp)
assert s.simplify(gain-7*(10279-3250*s.sqrt(10))/260000)==0
out={'method':'Exact SymPy rational/algebraic arithmetic; constructed from source Eqs. (18), (22), (23), independently of Lean execution.','noise_kraus_count':24,'composite_kraus_count':48,'lambda_values':[str(z) for z in [l1,l2,l3,l4,l5]],'kraus_completeness':True,'literal_eq23_matches_on_all_16_matrix_units':True,'group_fidelity_polynomials':[str(z) for z in polys],'fidelity_polynomial':str(f),'polar_fidelity':str(fp),'competitor_fidelity':str(fw),'gain':str(gain),'gain_decimal':str(s.N(gain,16)),'square_gap':10279**2-10*3250**2,'new_lean_run':False}
print(json.dumps(out,indent=2))
Path('verification.json').write_text(json.dumps(out,indent=2)+'\n')
