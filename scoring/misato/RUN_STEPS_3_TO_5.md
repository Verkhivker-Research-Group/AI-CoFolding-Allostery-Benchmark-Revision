# Run MISATO steps 3–5 visibly on this Windows machine

Open PowerShell and select this checkout:

```powershell
Set-Location 'C:\Users\Ryan\AI-CoFolding-Allostery-Benchmark-Revision-W-Misato'
```

The static-receptor QS run is **already complete** at
`misato_output/static_qs_all_v1.csv` (15,679/15,679 formally scored poses).
The gnina CSV was deliberately stopped at 4,143 scored rows so the remaining
work can run in your visible terminal. There is no active gnina process.
The source poses and reference CIFs have not been changed.

## Resume gnina CNNscore for all delivered poses

The executable and its CUDA libraries were copied into a persistent WSL
cache at `/home/ryan/.cache/misato_gnina_runtime/`. This is the same official
gnina v1.3.3 binary as the ignored copy under `misato_output/tools/`.
Run this command in PowerShell and leave the terminal open until it prints
the final status summary:

```powershell
wsl -d Ubuntu -- bash -lic 'conda activate plb && cd /mnt/c/Users/Ryan/AI-CoFolding-Allostery-Benchmark-Revision-W-Misato && python -u scoring/misato/rescore_gnina_wsl.py --gnina /home/ryan/.cache/misato_gnina_runtime/gnina --cuda-libs /home/ryan/.cache/misato_gnina_runtime/libs --workers 10 --batch --timeout-seconds 180 --resume --out-csv misato_output/gnina_rescore_all_v1.csv'
```

`--resume` skips all rows already written, so it does not rescore the first
4,143. The CSV is flushed after each target. Ten workers were tested on
this RTX 5070 Ti machine; change `--workers` only if you want to tune CPU
and GPU load. The script uses `--score_only` for each delivered SDF, never
redocks. A failed row retains a reason and blank CNNscore.

To see live progress, open **a second PowerShell terminal**, set the same
location, and run this read-only check as often as you like:

```powershell
python -c "import csv,collections; p='misato_output/gnina_rescore_all_v1.csv'; rows=list(csv.DictReader(open(p, newline='', encoding='utf-8'))); print(len(rows), '/ 16984 delivered poses'); print(dict(collections.Counter(row['status'] for row in rows))); print([(row['target_id'], row['method'], row['error'][:120]) for row in rows if row['status'] != 'scored'][:10])"
```

If WSL reports that the cached gnina executable is missing, recreate only
that cache from the already downloaded ignored files, then retry the
resume command:

```powershell
wsl -d Ubuntu -- bash -lc 'mkdir -p /home/ryan/.cache/misato_gnina_runtime/libs && cp -a /mnt/c/Users/Ryan/AI-CoFolding-Allostery-Benchmark-Revision-W-Misato/misato_output/tools/gnina.1.3.3.cuda12.8.static /home/ryan/.cache/misato_gnina_runtime/gnina && cp -a /mnt/c/Users/Ryan/AI-CoFolding-Allostery-Benchmark-Revision-W-Misato/misato_output/tools/gnina_cuda_libs/. /home/ryan/.cache/misato_gnina_runtime/libs/'
```

## Join the three metric families

When the progress check reaches **16,984 rows** and the gnina command has
exited, run:

```powershell
python scoring/misato/merge_static_metrics.py --gnina misato_output/gnina_rescore_all_v1.csv --out-csv misato_output/misato_static_metrics_v1.csv
```

The joined table contains formal pose RMSD and lDDT-PLI where scored,
`static-receptor QS`, native DiffDock confidence (EquiBind blank), gnina
`rescore_confidence` for both methods where available, and a blank pocket
RMSD with `receptor_mode = static_crystal`. It rejects incomplete or
duplicated gnina/QS cohorts rather than silently creating a partial table.
Do not rerun it over an existing output path; choose a new versioned name.

The final `evalspreadsheets/misato_v3/` exports, provenance sidecars, and
figures are **step 6 and not yet generated**. `METHODS_RESULTS.md` now
documents the validated methodology and clearly marks that remaining work.
