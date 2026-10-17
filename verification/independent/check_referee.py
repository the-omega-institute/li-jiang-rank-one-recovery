# Independently derived rational checks. Does not import the authors' checker.
import sympy as s
from itertools import product
R=s.Rational

def vec(m): return s.Matrix(list(m))
def chan(d,q):
 p=R(d*d-1)*(1-q*q)/(d*d-1+q*q)
 l1=1-p/(d*d-2); l4=(d*d-1+p)/(d*d);l5=(1-p)/(d*d)
 l2=(1-l1)*(d*d*l1-1)/(d*d*l4);l3=2*(1-l1)**2/(d*d*l4)
 a2=d*l4/(d*d-1)
 om=s.zeros(d*d,1)
 for i in range(d): om[i*d+i]=1
 P=s.eye(d*d)+(q-1)*om*om.T/d
 D=[]
 for i in range(d):
  di=s.zeros(d,d*d)
  for j in range(d):di[j,j*d+i]=1
  D.append(di)
 # reduced T as sum c M X M†, with rational c,M
 K=[]
 for i in range(d):
  for j in range(d*d):
   e=s.zeros(d,d*d);e[i,j]=1
   K.append((d*l3,e))
 for i in range(d):
  K.append((d*(l1-l3)*a2,D[i]*P))
  u=s.zeros(d,1);u[i]=1
  K.append((l2/l4,u*om.T-a2*q*D[i]*P))
 A=R(1,d)-d*l3
 C=2*A/d-2*(l1-l2-l3)*l4*q/(d*d-1)
 return p,l3,A,C,K,D,om

def check(d,q):
 p,l3,A,C,K,D,om=chan(d,q)
 # exact reduced channel TP
 assert sum((c*k.T*k for c,k in K),s.zeros(d*d))==s.eye(d*d)
 # Full encoder quadratic matrix direct versus invariant formula
 def erow(k): return vec(k.T).T
 Q=sum((c*erow(k).T*erow(k) for c,k in K),s.zeros(d**3))/d**2
 RR=s.zeros(d,d**3);ZZ=s.zeros(d,d**3)
 for i in range(d):
  RR[i,:]=erow(D[i])
  zz=s.zeros(d,d*d);zz[i,:]=om.T
  ZZ[i,:]=erow(zz)
 Q2=l3*s.eye(d**3)/d+(A*RR.T*RR-C*(RR.T*ZZ+ZZ.T*RR)/2+C*ZZ.T*ZZ/d)/d
 assert s.simplify(Q-Q2)==s.zeros(d**3)
 # Family decoder dual blocks at exact feasible point
 n=d-1
 if d==2:x,t=R(3,5),R(4,5)
 elif d==3:x,t=R(1,3),R(2,3)
 elif d==4:x,t=R(1,2),R(1,2)
 else:raise ValueError()
 V=s.zeros(d*d,d);V[0,0]=x
 for i in range(1,d):V[i*d+i,0]=-t;V[i*d,i]=1
 assert V.T*V==s.eye(d)
 KD=sum((c*vec(k*V)*vec(k*V).T for c,k in K),s.zeros(d*d))/d**2
 al=d*l3;be=d*A;cc=d*d*C/2;ss=x-n*t
 aa=al+be*x*x-cc*(2*ss*x/d-2*ss*ss/d**2)
 bb=al+be*t*t+cc*(2*ss*t/d+2*ss*ss/d**2)
 gg=be*x-cc*ss/d
 yy=[(aa+n*gg)/d**2]+[(gg+n*be+al)/d**2]*n
 slack=s.kronecker_product(s.diag(*yy),s.eye(d))-KD
 Id=vec(s.eye(d));assert slack*Id==s.zeros(d*d,1)
 assert sum(yy)==(Id.T*KD*Id)[0]
 # spectral multiset exact, all rational for this selected family point
 claimed=[0,(n+1)*gg/d**2]+[(aa+n*gg-al)/d**2]*n+[(gg+n*be+al-bb)/d**2]*n+[(gg+n*be)/d**2]*(n*n-1)
 actual=[]
 for ev,m in slack.eigenvals().items():actual += [ev]*m
 assert sorted(actual)==sorted(claimed),(actual,claimed)
 assert min(actual)>=0
 print('OK',d,q,'p',p,'matrix dimension',d**3,'decoder eigenvalues',actual)
for d,q in product([2,3,4],[R(1,4),R(1,2),R(3,4)]):check(d,q)
