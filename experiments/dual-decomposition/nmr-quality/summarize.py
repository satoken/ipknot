#!/usr/bin/env python3
"""Write per-RNA comparisons and compact, machine-readable benchmark evidence."""
import argparse
import csv
import json
import math
from pathlib import Path
from statistics import mean, median


def read(path, mode="dd"):
    return {row["pdb"]: row for row in json.loads((path/"results.json").read_text())
            if row["mode"] == mode}


def f1(a, b):
    a, b = set(map(tuple, a)), set(map(tuple, b))
    return 2*len(a & b)/(len(a)+len(b)) if a or b else 1.


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--baseline", type=Path, required=True)
    p.add_argument("--candidate", type=Path, required=True)
    p.add_argument("--ilp", type=Path, required=True)
    p.add_argument("--long-ilp", type=Path)
    p.add_argument("--timing-baseline", type=Path, nargs="+", required=True)
    p.add_argument("--timing-candidate", type=Path, nargs="+", required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    baseline, candidate, ilp = read(args.baseline), read(args.candidate), read(args.ilp, "ilp")
    if args.long_ilp: ilp.update(read(args.long_ilp, "ilp"))
    timings = [[read(path) for path in paths] for paths in [args.timing_baseline, args.timing_candidate]]
    assert baseline.keys() == candidate.keys() == ilp.keys()
    rows = []
    for pdb, new in candidate.items():
        old, exact = baseline[pdb], ilp[pdb]
        assert old["input_sha256"] == new["input_sha256"] == exact["input_sha256"]
        assert new["objective"] >= old["objective"]-1e-8
        for family, expected in zip(timings, [old, new]):
            for repetition in family:
                assert repetition[pdb]["pairs"] == expected["pairs"]
                assert math.isclose(repetition[pdb]["objective"], expected["objective"], abs_tol=1e-8)
        compared = exact["status"] == "ok"
        rows.append(dict(pdb=pdb, baseline_objective=old["objective"], candidate_objective=new["objective"],
                         ilp_objective=exact.get("objective"), ilp_status=exact["status"],
                         baseline_ilp_f1=f1(old["pairs"],exact["pairs"]) if compared else None,
                         candidate_ilp_f1=f1(new["pairs"],exact["pairs"]) if compared else None,
                         baseline_reference_f1=old["accuracy"]["f1"], candidate_reference_f1=new["accuracy"]["f1"],
                         baseline_violations=old["violations"], candidate_violations=new["violations"],
                         ilp_violations=exact.get("violations"),
                         baseline_seconds=median(r[pdb]["seconds"] for r in timings[0]),
                         candidate_seconds=median(r[pdb]["seconds"] for r in timings[1]),
                         ilp_seconds=exact["seconds"],
                         polish_states=new["summaries"][-1].get("polish_states",0)))
    completed = [r for r in rows if r["ilp_status"] == "ok"]
    summary = dict(rnas=len(rows), ilp_completed=len(completed),
                   objective_improved=sum(r["candidate_objective"]>r["baseline_objective"]+1e-8 for r in rows),
                   baseline_ilp_matches=sum(r["baseline_ilp_f1"]==1 for r in completed),
                   candidate_ilp_matches=sum(r["candidate_ilp_f1"]==1 for r in completed),
                   candidate_ilp_objective_matches=sum(math.isclose(r["candidate_objective"],r["ilp_objective"],abs_tol=1e-7) for r in completed),
                   baseline_reference_f1=mean(r["baseline_reference_f1"] for r in rows),
                   candidate_reference_f1=mean(r["candidate_reference_f1"] for r in rows),
                   baseline_seconds=sum(r["baseline_seconds"] for r in rows),
                   candidate_seconds=sum(r["candidate_seconds"] for r in rows),
                   completed_ilp_seconds=sum(r["ilp_seconds"] for r in completed),
                   completed_candidate_seconds=sum(r["candidate_seconds"] for r in completed))
    args.output.mkdir(exist_ok=True, parents=True)
    with (args.output/"comparison.csv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader(); writer.writerows(rows)
    (args.output/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    evidence = dict(baseline=list(baseline.values()), candidate=list(candidate.values()), ilp=list(ilp.values()),
                    timing_baseline=[list(r.values()) for r in timings[0]],
                    timing_candidate=[list(r.values()) for r in timings[1]])
    (args.output/"results.json").write_text(json.dumps(evidence,indent=2)+"\n")
    print(json.dumps(summary,indent=2))


if __name__ == "__main__":
    main()
