"""Multistart rank-one see-saw: exact decoder SDP, polar encoder updates.
Numerical lower bounds and decoder dual checks, never global certificates.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import argparse,csv,json,time
from pathlib import Path
import numpy as np
import cvxpy as cp
from model import *

class DecoderSDP:
    def __init__(self):
        self.H=cp.Parameter((8,8),hermitian=True)
        self.J=cp.Variable((8,8),hermitian=True)
        self.tp=self.J[:4,:4]+self.J[4:,4:]==np.eye(4)
        self.problem=cp.Problem(cp.Maximize(cp.real(cp.trace(self.J@self.H))),[self.J>>0,self.tp])
    def solve(self,H):
        self.H.value=(H+H.conj().T)/2
        value=self.problem.solve(solver='CLARABEL',tol_gap_abs=2e-10,tol_feas=2e-10,tol_gap_rel=2e-10,max_iter=150)
        if self.problem.status not in ('optimal','optimal_inaccurate'): raise RuntimeError(self.problem.status)
        J=self.J.value
        dec=kraus_from_choi(J)
        Y=self.tp.dual_value
        gap_eig=np.linalg.eigvalsh(np.kron(np.eye(2),Y)-H).min()
        upper=np.trace(Y).real+4*max(0,-gap_eig)
        lower=np.trace(sum(np.outer(d.reshape(8),d.reshape(8).conj()) for d in dec)@H).real
        return dec, dict(primal=float(value),lower=float(lower),dual_upper=float(upper),gap=float(upper-lower),tp_residual=float(np.linalg.norm(J[:4,:4]+J[4:,4:]-np.eye(4))),min_eigenvalue=float(np.linalg.eigvalsh(J).min()))

def run(grid,starts,outer,out):
    rng=np.random.default_rng(20261017); solver=DecoderSDP(); rows=[]; witnesses={}; runs={}; t0=time.time()
    for p in grid:
        K=channel(p); Vf,Ff,*_=family(p)
        # Independent direct Kraus check of the formula.
        assert abs(fidelity(K,Vf,D)-Ff)<2e-13
        dec,cert=solver.solve(decoder_cost(K,Vf))
        fixed_gain=fidelity(K,Vf,dec)-Ff
        bestF=Ff; bestV=Vf; bestD=np.array(D); searchF=-np.inf; searchV=Vf; searchD=np.array(D); max_cert_gap=cert['gap']; iterations=[]; tails=[]; final_values=[]
        for start in range(starts):
            V=Vf.copy() if start==0 else polar(rng.normal(size=(4,2))+1j*rng.normal(size=(4,2)))
            old=-1.; tail=np.inf
            for it in range(outer):
                dec,cc=solver.solve(decoder_cost(K,V)); max_cert_gap=max(max_cert_gap,cc['gap'])
                H=encoder_cost(K,dec)
                # Inner monotone RW / polar iteration for the encoder.
                for inner in range(100):
                    W=polar((H@V.reshape(8)).reshape(4,2))
                    delta=(np.vdot(W.reshape(8),H@W.reshape(8))-np.vdot(V.reshape(8),H@V.reshape(8))).real
                    V=W
                    if delta<2e-12: break
                F=fidelity(K,V,dec); tail=F-old
                if F>bestF: bestF=F; bestV=V.copy(); bestD=dec.copy()
                if F>searchF: searchF=F; searchV=V.copy(); searchD=dec.copy()
                if tail<2e-10 and it>=2: break
                old=F
            iterations.append(it+1); tails.append(tail); final_values.append(F)
        row=dict(p=float(p),family=Ff,family_decoder_gain=fixed_gain,joint_best=bestF,joint_gain=bestF-Ff,searched_best=float(searchF),searched_gain=float(searchF-Ff),max_decoder_gap=max_cert_gap,outer_max=max(iterations),max_final_increment=float(max(tails)),min_start_final=float(min(final_values)))
        rows.append(row)
        witnesses[f'{p:.12g}']={'V_real':searchV.real.tolist(),'V_imag':searchV.imag.tolist(),'D_real':searchD.real.tolist(),'D_imag':searchD.imag.tolist()}
        runs[f'{p:.12g}']={'final_values':[float(v) for v in final_values],'iterations':iterations,'last_increments':[float(v) for v in tails],'family_decoder_certificate':cert}
        print(json.dumps(row),flush=True)
        with open(out/'decoder-grid.csv','w') as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
        (out/'decoder-witnesses.json').write_text(json.dumps(witnesses))
        (out/'decoder-runs.json').write_text(json.dumps(runs,indent=2))
    meta={'grid':[float(p) for p in grid],'starts_per_p':starts,'outer_limit':outer,'inner_limit':100,'seed':20261017,'solver':'CLARABEL','solver_tolerances':2e-10,'outer_stop':2e-10,'inner_stop':2e-12,'elapsed_seconds':time.time()-t0,'max_joint_gain':max(r['joint_gain'] for r in rows),'max_searched_gain':max(r['searched_gain'] for r in rows),'min_searched_gain':min(r['searched_gain'] for r in rows),'max_family_decoder_gain':max(r['family_decoder_gain'] for r in rows),'max_decoder_gap':max(r['max_decoder_gap'] for r in rows),'unconverged_points':int(sum(r['max_final_increment']>=2e-10 for r in rows)),'versions':{'numpy':np.__version__,'cvxpy':cp.__version__}}
    (out/'numerics-summary.json').write_text(json.dumps(meta,indent=2)+'\n')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--smoke',action='store_true');ap.add_argument('--starts',type=int,default=4);ap.add_argument('--outer',type=int,default=80)
    args=ap.parse_args(); out=Path(__file__).resolve().parent
    grid=[.1,9/13,.9] if args.smoke else sorted(set([k/100 for k in range(1,100)]+[9/13]))
    run(grid,args.starts,args.outer,out)
