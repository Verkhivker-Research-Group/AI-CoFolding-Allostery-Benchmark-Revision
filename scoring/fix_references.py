"""Replace WRONG references for currently-failed DynamicBind targets.

For each failed-but-fixable target (reference_remap.csv category == uniprot+lig
AND not currently scored), move the existing wrong <UNIPROT>_<LIG>.cif.gz into a
backup dir and download the correct co-crystal PDB from RCSB under the same name.

Working targets are NOT touched. Read source of truth = reference_remap.csv.
"""
import sys, gzip, shutil, time
from pathlib import Path
import requests
sys.path.insert(0, "/mnt/c/Users/Ryan/AI-CoFolding-Allostery-Benchmark/scoring/plb_bench")
import pandas as pd

REFS = Path("/mnt/c/Users/Ryan/AI-CoFolding-Archive/references_ref_cifs")
BACKUP = Path("/mnt/c/Users/Ryan/AI-CoFolding-Archive/references_ref_cifs_backup_wrong")
BACKUP.mkdir(exist_ok=True)

remap = pd.read_csv("/mnt/c/Users/Ryan/AI-CoFolding-Allostery-Benchmark/scoring/reference_remap.csv")
df = pd.read_parquet("/mnt/c/Users/Ryan/AI-CoFolding-Allostery-Benchmark/plb_bench_output/pla/benchmark.parquet")
db = df[df.producer == "dynamicbind_pla"]
scored = db.groupby("pdb_id")["bisy_rmsd"].apply(lambda s: s.notna().any())
remap = remap.merge(scored.rename("scored"), left_on="pdb_id", right_index=True, how="left")

todo = remap[(~remap.scored) & (remap.category == "uniprot+lig")].copy()
print("targets to fix:", len(todo))

ok = fail = 0
for _, r in todo.iterrows():
    pid = r.pdb_id; pdb = r.chosen_pdb
    name = pid.upper() + ".cif.gz"          # get_reference key
    dest = REFS / name
    # 1) back up existing wrong file (MOVE)
    if dest.exists():
        bdest = BACKUP / name
        if not bdest.exists():
            shutil.move(str(dest), str(bdest))
    # 2) download correct PDB (deposited entry, contains the ligand)
    try:
        resp = requests.get(f"https://files.rcsb.org/download/{pdb}.cif.gz", timeout=60)
        if resp.status_code != 200:
            print("  HTTP", resp.status_code, pid, pdb); fail += 1; continue
        # sanity: ligand present?
        txt = gzip.decompress(resp.content).decode("utf-8", "replace")
        lig = r.ligand
        if not any(l.startswith("HETATM") and l.split()[5] == lig for l in txt.splitlines()):
            print("  WARN ligand", lig, "absent in", pdb, "for", pid);
        with open(dest, "wb") as fh:
            fh.write(resp.content)
        ok += 1
    except Exception as e:
        print("  ERR", pid, pdb, e); fail += 1
    time.sleep(0.08)
    if (ok + fail) % 25 == 0:
        print("  ...", ok + fail, "/", len(todo))

print("done. downloaded:", ok, " failed:", fail)
print("backups in:", BACKUP)
