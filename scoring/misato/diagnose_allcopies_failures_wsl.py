"""Get OpenStructure failure states for the small unscored all-copies tail."""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import score_static_docking_allcopies_wsl as scorer_module
from score_static_docking_wsl import initialize_worker


def read_csv(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def reason(scorer) -> str:
    value, _ = scorer_module._assigned_original(scorer)
    if value is not None:
        return "scored"
    try:
        issue, description = scorer.guess_model_ligand_unassigned_reason(0)
        return f"{issue}: {description}"
    except Exception as error:  # noqa: BLE001 - preserve OST diagnostic failures
        return f"unresolved:{type(error).__name__}:{error}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", type=Path,
                        default=Path("misato_output/misato_score_allcopies_v1.csv"))
    parser.add_argument("--out-csv", type=Path, required=True)
    parser.add_argument("--timeout-seconds", type=int, default=120)
    args = parser.parse_args()
    if args.out_csv.exists():
        parser.error("Output exists")
    root = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(root / "scoring" / "plb_bench"))
    initialize_worker()
    scorer_module._assigned_original = scorer_module._assigned
    results = []
    for row in read_csv(args.scores):
        if row["status"] != "partial_or_no_metrics":
            continue
        issues = []

        def capture(scorer, issues=issues):
            issues.append(reason(scorer))
            return scorer_module._assigned_original(scorer)

        scorer_module._assigned = capture
        rescored = scorer_module.score_target(row["target_id"], [row], str(root),
                                               args.timeout_seconds)[0]
        results.append({"target_id": row["target_id"], "method": row["method"],
                        "old_status": row["status"], "retry_status": rescored["status"],
                        "bisy_reason": issues[0] if issues else "not_attempted",
                        "lddt_reason": issues[1] if len(issues) > 1 else "not_attempted",
                        "retry_error": rescored["error"]})
    with args.out_csv.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["target_id", "method", "old_status",
                                                     "retry_status", "bisy_reason", "lddt_reason",
                                                     "retry_error"])
        writer.writeheader()
        writer.writerows(results)
    print(f"Diagnosed {len(results)} residual rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
