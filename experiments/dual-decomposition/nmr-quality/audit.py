#!/usr/bin/env python3
"""Independently audit recorded DD incumbents, domains, rows and RNA planarity."""
import argparse
import json
import math
from pathlib import Path


def audit(path):
    problems, summaries, final = {}, {}, {}
    checked = 0
    for line in path.read_text().splitlines():
        event = json.loads(line)
        solve = event["solve_id"]
        if not event.get("constrained"):
            continue
        if event["event"] == "problem":
            assert event["noe_solver"] == "relaxed", "This auditor expects all rows in the relaxed trace"
            problems[solve] = event
            continue
        if event["event"] == "summary":
            summaries[solve] = event
            assert event["lower_bound"] <= event["upper_bound"] + 1e-6
            budget = problems[solve]["repair_state_budget"]
            assert not budget or event["repair_states"] <= budget
            assert event.get("polish_states", 0) <= problems[solve].get("nmr_polish_state_budget", 0)
        if event["event"] == "iteration" and event["lower_bound"] is not None:
            values, objective = event["recovered"], event["lower_bound"]
        elif event["event"] == "repair":
            values, objective = event["selected"], event["objective"]
            final[solve] = objective
        else:
            continue
        problem = problems[solve]
        assert "variable_domains" in problem, "Run with --trace-state"
        assert len(values) == len(problem["variable_domains"])
        score = 0
        for value, (coefficient, lower, upper, integer) in zip(values, problem["variable_domains"]):
            assert math.isfinite(value) and lower-1e-7 <= value <= upper+1e-7
            assert not integer or abs(value-round(value)) < 1e-7
            score += coefficient*value
        assert math.isclose(score, objective, abs_tol=1e-7), (score, objective)
        for bound, equality, terms in problem["rows"]:
            total = sum(values[col]*coefficient for col, coefficient in terms)
            assert abs(total-bound) < 1e-7 if equality else total <= bound+1e-7
        selected = [(i,j,level) for i,j,level,col in problem["pairs"] if values[col] > .5]
        for a, (i,j,level) in enumerate(selected):
            for k,l,other in selected[:a]:
                # Applications impose endpoint exclusion across all levels.
                assert not {i,j} & {k,l}
                if level == other:
                    assert not (i < k < j < l or k < i < l < j)
        checked += 1
    assert set(final) == set(problems), "Missing final state dump"
    for solve, objective in final.items():
        assert math.isclose(objective, summaries[solve]["lower_bound"], abs_tol=1e-7)
    return checked


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directories", nargs="+", type=Path)
    args = parser.parse_args()
    total = 0
    files = 0
    for directory in args.directories:
        for path in sorted(directory.glob("*-dd*/trace.jsonl")):
            total += audit(path)
            files += 1
    assert files, "No traces found"
    print(f"Audited {total} feasible incumbents in {files} traces.")
