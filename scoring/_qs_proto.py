"""Robust ligand-inclusive interface-QS via CA superposition (no aln fragility).

Method (per model vs corrected reference, same protein):
 1. Superpose model protein onto reference via matched CA (Kabsch).
 2. Map each model residue -> nearest reference residue (post-superposition).
 3. Reference contact vector over ref residues: in contact with REF ligand
    (min heavy-atom distance <= contact_d).
 4. Model contact vector over the SAME ref residues: the mapped model residue
    is in contact with the MODEL ligand (min heavy-atom distance <= contact_d).
 5. QS-style score over union of contacting ref residues:
        QS = sum_shared w / (sum_shared w + n_nonshared)
    with w = 1 - |d_ref - d_mdl|/contact_d  (clamped to [0,1]).
Continuous, robust, reflects how well the protein-ligand interface is reproduced.
"""
import sys, tempfile, os
from pathlib import Path
import numpy as np
sys.path.insert(0, "/mnt/c/Users/Ryan/AI-CoFolding-Allostery-Benchmark/scoring/plb_bench")
from plb_bench.references import get_reference, read_reference_text
from plb_bench.scoring import _get_ref_lig_name
REFS = Path("/mnt/c/Users/Ryan/AI-CoFolding-Archive/references_ref_cifs")
DATA = Path("/mnt/c/Users/Ryan/AI-CoFolding-Archive/plb_bench_data/dynamicbind_pla")
from ost import io as ost_io

def load(text):
    with tempfile.NamedTemporaryFile("w", suffix=".cif", delete=False) as f:
        f.write(text); tmp=f.name
    try:
        r = ost_io.LoadMMCIF(tmp, seqres=True, info=True, fault_tolerant=True)
        return r[0] if isinstance(r, tuple) else r
    finally: os.unlink(tmp)

def poly_cas(ent):
    """list of (reskey, ca_xyz, heavy_atoms_xyz) for peptide residues."""
    out=[]
    for ch in ent.chains:
        for r in ch.residues:
            ca=r.FindAtom("CA")
            if not ca.IsValid(): continue
            heavy=np.array([[a.pos.x,a.pos.y,a.pos.z] for a in r.atoms if a.element!="H"])
            out.append(((ch.name,r.number.num),
                        np.array([ca.pos.x,ca.pos.y,ca.pos.z]), heavy))
    return out

def lig_heavy(ent, names):
    pts=[]
    for ch in ent.chains:
        for r in ch.residues:
            if r.name in names:
                for a in r.atoms:
                    if a.element!="H": pts.append([a.pos.x,a.pos.y,a.pos.z])
    return np.array(pts) if pts else None

def kabsch(P,Q):
    Pc,Qc=P.mean(0),Q.mean(0)
    H=(P-Pc).T@(Q-Qc); U,_,Vt=np.linalg.svd(H)
    d=np.sign(np.linalg.det(Vt.T@U.T))
    R=Vt.T@np.diag([1,1,d])@U.T
    return R, Qc-R@Pc

def min_dist(res_heavy, lig):
    if res_heavy.size==0 or lig is None: return np.inf
    from scipy.spatial import distance
    return distance.cdist(res_heavy, lig).min()

def qs_iface(mdl, ref, ref_lig, contact_d=12.0):
    rcas=poly_cas(ref); mcas=poly_cas(mdl)
    if len(rcas)<3 or len(mcas)<3: return None
    rl=lig_heavy(ref,{ref_lig}); ml=lig_heavy(mdl,{"LIG"})
    if rl is None or ml is None: return None
    # match CAs by residue number (same protein); fallback positional
    rmap={k:(ca,hv) for k,ca,hv in rcas}
    mmap={k:(ca,hv) for k,ca,hv in mcas}
    common=[k for k in rmap if k in mmap]
    if len(common)<3:
        # positional by order
        n=min(len(rcas),len(mcas))
        rsel=rcas[:n]; msel=mcas[:n]
    else:
        rsel=[(k,)+rmap[k] for k in common]; msel=[(k,)+mmap[k] for k in common]
        rsel=[(a,b,c) for a,b,c in rsel]; msel=[(a,b,c) for a,b,c in msel]
    Rca=np.array([x[1] for x in rsel]); Mca=np.array([x[1] for x in msel])
    R,t=kabsch(Mca,Rca)               # transform model -> ref frame
    # transform model ligand + model residue heavy atoms
    ml_t=(R@ml.T).T+t
    # reference contact set over ref residues
    ref_d={k:min_dist(hv,rl) for k,_,hv in rcas}
    # model contact set, mapped to ref residue by nearest CA post-superposition
    refkeys=[k for k,_,_ in rcas]; refca=np.array([ca for _,ca,_ in rcas])
    mdl_d_byref={}
    for k,ca,hv in mcas:
        ca_t=(R@ca)+t
        j=int(np.argmin(np.linalg.norm(refca-ca_t,axis=1)))
        hv_t=(R@hv.T).T+t
        d=min_dist(hv_t, ml_t)
        rk=refkeys[j]
        if rk not in mdl_d_byref or d<mdl_d_byref[rk]:
            mdl_d_byref[rk]=d
    # union of contacting residues
    keys=set(k for k,d in ref_d.items() if d<=contact_d) | \
         set(k for k,d in mdl_d_byref.items() if d<=contact_d)
    if not keys: return None
    shared_w=0.0; nonshared=0
    for k in keys:
        dr=ref_d.get(k,np.inf); dm=mdl_d_byref.get(k,np.inf)
        if dr<=contact_d and dm<=contact_d:
            shared_w+=max(0.0,1.0-abs(dr-dm)/contact_d)
        else:
            nonshared+=1
    denom=shared_w+nonshared
    return shared_w/denom if denom>0 else None

import requests, gzip
def reftext(pid):
    return read_reference_text(get_reference(pid,REFS,allow_download=False)[0])

# test: good poses (low bisy) should give higher QS than bad poses (high bisy)
import pandas as pd
bench=pd.read_parquet("/mnt/c/Users/Ryan/AI-CoFolding-Allostery-Benchmark/plb_bench_output/pla/benchmark.parquet")
for pid in ["b6ywb8_dtp","o14965_6f2","o15530_3q3","p00044_6vb","p00533_az1" ]:
    if pid not in set(bench.pdb_id):
        # pick a real one
        continue
    ref=load(reftext(pid)); rlig=_get_ref_lig_name(ref)
    rows=[]
    for i in range(10):
        mp=DATA/pid/f"model_{i:03d}.cif"
        if not mp.exists(): continue
        mdl=load(mp.read_text())
        bisy=bench[(bench.pdb_id==pid)&(bench.model_idx==i)]["bisy_rmsd"]
        b=float(bisy.iloc[0]) if len(bisy) and pd.notna(bisy.iloc[0]) else None
        try: v=qs_iface(mdl,ref,rlig)
        except Exception as e: v=f"ERR:{type(e).__name__}:{str(e)[:40]}"
        rows.append((i,round(v,3) if isinstance(v,float) else v, round(b,2) if b else None))
    print(pid,"reflig",rlig)
    for i,v,b in rows: print(f"    model {i}: QS={v}  bisy={b}")
