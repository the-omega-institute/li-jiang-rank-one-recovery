"""Recompute saved numerical witnesses directly from the 24 noise Kraus operators."""
from pathlib import Path
import csv,json
import numpy as np
from model import *
root=Path(__file__).resolve().parent
rows=list(csv.DictReader((root/'decoder-grid.csv').open()))
witnesses=json.loads((root/'decoder-witnesses.json').read_text())
runs=json.loads((root/'decoder-runs.json').read_text())
rng=np.random.default_rng(817263)
max_iso=max_tp=max_ferr=max_cost=max_sym=0.; max_random_excess=-np.inf; worst_run_loss=0.; nrandom=0
for row in rows:
 p=float(row['p']); K=channel(p); Vf,Ff,*_=family(p); key=f'{p:.12g}'
 z=witnesses[key]
 V=np.array(z['V_real'])+1j*np.array(z['V_imag']); ds=np.array(z['D_real'])+1j*np.array(z['D_imag'])
 max_iso=max(max_iso,np.linalg.norm(V.conj().T@V-I2))
 max_tp=max(max_tp,np.linalg.norm(sum(d.conj().T@d for d in ds)-I4))
 F=fidelity(K,V,ds)
 max_ferr=max(max_ferr,abs(F-float(row['searched_best'])))
 j=sum(np.outer(d.reshape(8),d.reshape(8).conj()) for d in ds)
 max_cost=max(max_cost,abs(np.trace(j@decoder_cost(K,V)).real-F),abs(np.vdot(V.reshape(8),encoder_cost(K,ds)@V.reshape(8)).real-F))
 worst_run_loss=max(worst_run_loss,Ff-min(runs[key]['final_values']))
 for _ in range(20):
  W=polar(rng.normal(size=(4,2))+1j*rng.normal(size=(4,2)))
  U=polar(rng.normal(size=(2,2))+1j*rng.normal(size=(2,2)))
  transformed=np.kron(U,U.conj())@W@U.conj().T
  FW=fidelity(K,W,D)
  max_sym=max(max_sym,abs(FW-fidelity(K,transformed,D)))
  max_random_excess=max(max_random_excess,FW-Ff)
  contraction=W@U@np.diag(rng.uniform(size=2))@U.conj().T
  max_random_excess=max(max_random_excess,fidelity(K,contraction,D)-Ff)
  nrandom+=2
assert max(max_iso,max_tp,max_ferr,max_cost,max_sym)<2e-12
assert max_random_excess<1e-10
assert worst_run_loss<1e-8
out={'points':len(rows),'saved_witnesses':len(witnesses),'saved_runs':sum(len(z['final_values']) for z in runs.values()),'max_isometry_residual':float(max_iso),'max_decoder_tp_residual':float(max_tp),'max_fidelity_recompute_error':float(max_ferr),'max_choi_cost_error':float(max_cost),'max_covariance_error':float(max_sym),'random_encoders_checked':nrandom,'max_random_encoder_excess':float(max_random_excess),'worst_start_final_loss_vs_family':float(worst_run_loss)}
print(json.dumps(out,indent=2));(root/'validation-results.json').write_text(json.dumps(out,indent=2)+'\n')
