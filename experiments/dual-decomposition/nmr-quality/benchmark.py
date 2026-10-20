#!/usr/bin/env python3
"""Compare NMR DD and ILP on identical saved posteriors and observations."""
import argparse
import concurrent.futures
import hashlib
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import threading
import time


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def metrics(reference, prediction):
    tp = len(reference & prediction)
    fp, fn = len(prediction - reference), len(reference - prediction)
    return dict(tp=tp, fp=fp, fn=fn,
                sensitivity=tp / len(reference) if reference else None,
                ppv=tp / len(prediction) if prediction else float(not reference),
                f1=2*tp / (2*tp+fp+fn) if 2*tp+fp+fn else 1.0)


def reference_pairs(structure):
    matching = dict(zip(")]}>", "([{<"))
    stacks = {c: [] for c in matching.values()}
    result = set()
    for i, c in enumerate(structure, 1):
        if c in stacks: stacks[c].append(i)
        elif c in matching: result.add((stacks[matching[c]].pop(), i))
        elif c != ".": raise ValueError(f"Unsupported bracket: {c}")
    assert not any(stacks.values())
    return result


def bpseq_pairs(path):
    result = set()
    for line in path.read_text().splitlines():
        if line and line[0].isdigit():
            left, base, right = line.split()
            if 0 < int(left) < int(right): result.add((int(left), int(right)))
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data", type=Path, help="Optional existing new_nmr directory")
    p.add_argument("--inputs", type=Path, default=Path(__file__).with_name("inputs.json"))
    p.add_argument("--posteriors", type=Path, default=Path(__file__).with_name("posteriors"))
    p.add_argument("--binary", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--modes", nargs="+", default=["dd", "ilp"])
    p.add_argument("--threshold", default=".25,.25")
    p.add_argument("--penalty", default=".1")
    p.add_argument("--counts", action="store_true", help="Also use HSQC counts (soft, exact, same penalty)")
    p.add_argument("--timeout", type=float, default=120)
    p.add_argument("--cpus", default="5,6,7,8")
    p.add_argument("--pdbs", nargs="*")
    p.add_argument("--trace-state", action="store_true")
    p.add_argument("--extra", nargs=argparse.REMAINDER, default=[])
    args = p.parse_args()
    if args.data:
        sys.path.insert(0, str(args.data.resolve()))
        import nmr_followup_experiments as data
        entries = data.read_inputs()
    else:
        entries = json.loads(args.inputs.read_text())
    if args.pdbs:
        entries = [e for e in entries if e["pdb"] in args.pdbs]
    args.output.mkdir(parents=True, exist_ok=True)
    binary = args.binary.resolve()
    manifest = dict(binary=str(binary), binary_sha256=digest(binary), inputs=entries,
                    arguments={k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
                    runner_sha256=digest(__file__))
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    cpus = list(map(int, args.cpus.split(",")))
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    def queue(cpu, assigned):
        results = []
        for entry in assigned:
            source = args.posteriors / f'L2_{entry["pdb"]}_linear/posterior.bpp'
            if not source.exists(): source = args.posteriors / (entry["pdb"] + ".bpp")
            posterior = args.output / (entry["pdb"] + ".bpp")
            posterior.write_text("".join(line + "\n" for line in source.read_text().splitlines()
                                         if line and not line.startswith("#")))
            for mode in args.modes:
                directory = args.output / (entry["pdb"] + "-" + mode)
                directory.mkdir(exist_ok=True)
                command = [str(binary), "--decoder", "ilp" if mode == "ilp" else "dd",
                           "-x", "-r", "0", "-n", "1", "-t", args.threshold, "--loglevel", "debug",
                           "--nmr-soft", "--nmr-stack-penalty", args.penalty,
                           "--nmr-bulge-mode", "all", "--coaxial-stacking",
                           "--nmr-coaxial-no-adjacent-bulge", "-B", str(directory / "prediction.bpseq")]
                for pattern in entry["patterns"]:
                    command += ["--stack-constraint", pattern]
                if args.counts:
                    command += ["--base-pairs", entry["counts"], "--nmr-count-penalty", args.penalty]
                if mode != "ilp":
                    command += ["--dd-trace", str(directory / "trace.jsonl")]
                    if args.trace_state: command += ["--dd-trace-state"]
                    if mode == "dd-more": command += ["--dd-constraint-states", "16384"]
                    if mode == "dd-periodic": command += ["--dd-constraint-recovery-every", "10"]
                    if mode == "dd-full": command += ["--dd-constraints", "full"]
                    command += args.extra
                command.append(str(posterior))
                target = directory / "result.json"
                signature = hashlib.sha256((json.dumps(command) + digest(binary) + digest(posterior)).encode()).hexdigest()
                if target.exists():
                    row = json.loads(target.read_text())
                    if row["signature"] != signature: raise ValueError("Incompatible saved run")
                    results.append(row)
                    continue
                (directory / "trace.jsonl").unlink(missing_ok=True)
                start = time.perf_counter()
                with (directory / "run.log").open("w") as log:
                    process = subprocess.Popen(["taskset", "-c", str(cpu), *command], env=env,
                                               stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
                    expired = threading.Event()
                    def expire():
                        if process.poll() is None:
                            expired.set()
                            try: os.killpg(process.pid, signal.SIGKILL)
                            except ProcessLookupError: pass
                    timer = threading.Timer(args.timeout, expire)
                    timer.start()
                    try: code = process.wait()
                    finally:
                        timer.cancel()
                        timer.join()
                    status = "timeout" if expired.is_set() else "ok" if code == 0 else "error"
                text = (directory / "run.log").read_text()
                row = dict(pdb=entry["pdb"], mode=mode, status=status, command=command,
                           seconds=time.perf_counter()-start, cpu=cpu, signature=signature,
                           input_sha256=digest(posterior), binary_sha256=digest(binary))
                if status == "ok":
                    pairs = bpseq_pairs(directory / "prediction.bpseq")
                    row.update(pairs=sorted(pairs),
                               accuracy=metrics(reference_pairs(entry["structure"]), pairs))
                    row["violations"] = len(re.findall(r"Soft NMR stack constraint \d+ violated", text))
                    objectives = re.findall(r"(?:DD|IP) objective: ([-+0-9.eE]+)", text)
                    if mode == "ilp":
                        row["objective"] = float(objectives[-1])
                    else:
                        events = [json.loads(line) for line in (directory / "trace.jsonl").read_text().splitlines()]
                        row["summaries"] = [e for e in events if e["event"] == "summary"]
                        row["objective"] = row["summaries"][-1]["lower_bound"]
                    # Automatic threshold selection compares models with different
                    # coefficients; its last subsolve is not the returned model.
                    if "auto" in args.threshold:
                        row["objective"] = None
                        row["violations"] = None
                else:
                    row["error"] = text[-3000:]
                target.write_text(json.dumps(row, indent=2) + "\n")
                results.append(row)
                print(json.dumps({k: row.get(k) for k in ("pdb", "mode", "status", "objective", "seconds", "violations")}), flush=True)
        return results
    with concurrent.futures.ThreadPoolExecutor(max_workers=len(cpus)) as pool:
        jobs = [pool.submit(queue, cpu, entries[i::len(cpus)]) for i, cpu in enumerate(cpus)]
        rows = [row for job in jobs for row in job.result()]
    assert digest(binary) == manifest["binary_sha256"]
    (args.output / "results.json").write_text(json.dumps(rows, indent=2) + "\n")


if __name__ == "__main__":
    main()
