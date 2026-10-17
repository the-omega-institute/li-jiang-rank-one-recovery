"""Validate saved numerical witnesses against the original physical noise.

This reconstructs Eq.23 directly, independently of search_qutrit.model.
It remains floating-point validation and is not an interval certificate.
"""
import json
from pathlib import Path
import numpy as np

d=3;p=8/11;N=d*d
l1=1-p/(N-2);l4=(N-1+p)/N;l5=(1-p)/N
l2=(1-l1)*(N*l1-1)/(N*l4);l3=2*(1-l1)**2/(N*l4)
eye=np.eye(d);Omega=eye.reshape(-1)/np.sqrt(d)
P=np.outer(Omega,Omega)
Q=np.sqrt(d*l4/(N-1))*(np.eye(N)-P)+np.sqrt(d*l5)*P
Ds=[np.kron(eye,eye[j:j+1]) for j in range(d)]
Ks=[np.sqrt(l3)*np.outer(np.eye(N)[:,a],np.eye(N)[:,b]) for a in range(N) for b in range(N)]
Ks += [np.sqrt(l1-l3)*np.kron(eye,np.outer(eye[:,i],eye[:,j]))@Q for i in range(d) for j in range(d)]
Ks += [np.sqrt(l2/l4)*(np.outer(np.eye(N)[:,d*i+j],Omega)-np.sqrt(l5)*Ds[j].T@Ds[i]@Q)
       for i in range(d) for j in range(d)]
unpack=lambda a: np.array(a['real'])+1j*np.array(a['imag'])
data=json.loads(Path(__file__).with_name('qutrit-witnesses.json').read_text())
runs=json.loads(Path(__file__).with_name('qutrit-search.json').read_text())['runs']
max_f_error=0;max_iso=0;max_tp=0
for witness,run in zip(data,runs):
    V=unpack(witness['V']);R=[unpack(z) for z in witness['decoder_kraus']]
    F=sum(abs(np.trace(r@Dj@K@V))**2 for r in R for Dj in Ds for K in Ks)/d**2
    max_f_error=max(max_f_error,abs(F-run['fidelity']))
    max_iso=max(max_iso,np.linalg.norm(V.conj().T@V-eye))
    max_tp=max(max_tp,np.linalg.norm(sum(r.conj().T@r for r in R)-eye))
out={'scope':'floating-point direct physical Kraus recomputation of all saved witnesses',
     'witness_count':len(data),'physical_kraus_count':len(Ks),
     'max_fidelity_difference':float(max_f_error),'max_encoder_residual':float(max_iso),
     'max_decoder_residual':float(max_tp),
     'physical_noise_tp_residual':float(np.linalg.norm(sum(K.T@K for K in Ks)-np.eye(N)))}
assert max_f_error<1e-12 and max_iso<1e-12 and max_tp<1e-12
assert out['physical_noise_tp_residual']<1e-12
Path(__file__).with_name('qutrit-validation.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
