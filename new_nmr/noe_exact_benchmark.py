#!/usr/bin/env python3
"""Compare ordered model fingerprints and archived exact predictions."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
import os
from pathlib import Path
import re
import random
import shlex
import signal
import subprocess
import sys
import time
from noe_build_helpers import ilp_decoder_args, legacy_ip_api, link_command

ROOT = Path(__file__).resolve().parent


def recorder(build, output, counts_only=False):
    obj = output.with_suffix(".o")
    flags = ["-DNMR_RECORD_COUNTS_ONLY"] if counts_only else []
    if legacy_ip_api(build):
        flags.append("-DNMR_RECORDER_LEGACY_IP_API")
    subprocess.run(["c++", "-std=c++17", "-O3", *flags, "-I", str(ROOT.parent / "src"),
                    "-c", str(ROOT.parent / "tests/nmr_model_recorder.cpp"), "-o", str(obj)], check=True)
    args = link_command(build)
    args = [arg for arg in args if not arg.startswith("-Wl,--dependency-file=") and
            arg != "CMakeFiles/ipknot.dir/src/ip.cpp.o"]
    args.insert(1, str(obj))
    args[args.index("-o") + 1] = str(output)
    subprocess.run(args, cwd=build, check=True)


def run(command, timeout, log):
    environment = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
                       PYTHONHOME="/tmp/ipknot-nmr-python/cpython-3.12.14-linux-aarch64-gnu",
                       PYTHONPATH="/tmp/ipknot-nmr-mxfold2-native-venv/lib/python3.12/site-packages")
    start = time.perf_counter()
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               text=True, env=environment, start_new_session=True)
    try:
        output, _ = process.communicate(timeout=timeout)
        status = "ok" if process.returncode == 0 else f"exit_{process.returncode}"
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGKILL)
        output, _ = process.communicate()
        status = "timeout"
    log.write_text(output)
    return dict(status=status, seconds=time.perf_counter() - start, output=output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--before-build", type=Path, default=Path("/tmp/ipknot-nmr-followup-final-build"))
    parser.add_argument("--after-build", type=Path, default=Path("/tmp/ipknot-noe-optimized-build"))
    parser.add_argument("--output", type=Path, default=ROOT / "noe_optimization_results/final")
    parser.add_argument("--mode", choices=("models", "predictions", "formulations", "fixtures", "random"), required=True)
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--case", help="Run only cases whose name matches this regular expression")
    parser.add_argument("--engine", choices=("lpc", "lpv"), help="Replace the archived probability engine")
    parser.add_argument("--repeats", type=int, default=1, help="Number of paired runs per case")
    args = parser.parse_args()
    if args.mode == "predictions" and args.engine:
        parser.error("--engine is not supported when comparing archived MXfold2 predictions")
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    source = ROOT / "nmr_followup_results"
    cases = []
    for level in (1, 2):
        for path in sorted((source / f"mxfold2_topology_level{level}").glob("*/*/result.json")):
            row = json.loads(path.read_text())
            if args.mode in ("models", "formulations") or row["status"] == "ok":
                cases.append((f"L{level}_{row['condition']}_{row['pdb']}", row, path))
    if args.mode == "fixtures":
        suite = json.loads(subprocess.check_output(["ctest", "--test-dir", str(args.after_build),
                                                   "--show-only=json-v1"], text=True))
        cases = [(test["name"], dict(command=test["command"], pdb=test["name"], condition="fixture"), None)
                 for test in suite["tests"] if Path(test["command"][0]).name == "ipknot"]
    if args.mode == "random":
        cases = []
        generator = random.Random(20261004)
        for seed in range(6):
            fasta = args.output / f"random_{seed}.fa"
            sequence = "".join(generator.choice("ACGUacgutN") for _ in range(12 + 8 * seed))
            fasta.write_text(f">random_{seed}\n{sequence}\n")
            for level in (1, 2):
                for mode in ("none", "fallback", "all"):
                    for shared in (False, True):
                        for neighbors in (False, True):
                            name = f"random_{seed}_L{level}_{mode}_shared{int(shared)}_neighbors{int(neighbors)}"
                            command = ["ipknot", "-r", "0", "-t", ",".join(["0.2"] * level),
                                       "--nmr-soft", "--nmr-bulge-mode", mode, "--coaxial-stacking",
                                       "--base-pairs", "AA=2,CC=1,UU=3,AG=1,AC=2,CU=1",
                                       "--nmr-count-mode", "lower-bound" if seed % 2 else "exact"]
                            for pattern in ("GC AU", "(GU)-G-U", "UU CC", "GC AU"):
                                command += ["--stack-constraint", pattern]
                            if shared: command += ["--nmr-allow-shared-stack-pairs"]
                            if not neighbors: command += ["--without-canonical-neighbor"]
                            command.append(str(fasta))
                            cases.append((name, dict(command=command, pdb=str(seed), condition=mode), None))
    # One additional completed, formerly timed-out bulge case.
    if args.mode == "predictions":
        path = source / "mxfold2_topology_bulge_extended_level2/stacks_s1_ball_c0/2N8V/result.json"
        cases.append(("L2_extended_bulge_2N8V", json.loads(path.read_text()), path))
    if args.case:
        cases = [case for case in cases if re.search(args.case, case[0])]
        if not cases:
            parser.error("--case did not match any cases")
    if args.mode != "predictions":
        binaries = [args.output / "before_recorder", args.output / "after_recorder"]
        recorder(args.before_build, binaries[0], args.mode == "formulations")
        recorder(args.after_build, binaries[1], args.mode == "formulations")
    else:
        binaries = [args.after_build / "ipknot"]
    manifest = dict(mode=args.mode, engine=args.engine, repeats=args.repeats,
                    command=shlex.join(sys.argv),
                    binaries={str(binary): hashlib.sha256(binary.read_bytes()).hexdigest() for binary in binaries},
                    objects={})
    for label, build in (("before", args.before_build), ("after", args.after_build)):
        for name in ("ip.cpp", "ipknot.cpp", "main.cpp", "fold.cpp", "linearpartition/LinearPartition.cpp"):
            obj = build / f"CMakeFiles/ipknot.dir/src/{name}.o"
            manifest["objects"][f"{label}_{name}"] = hashlib.sha256(obj.read_bytes()).hexdigest()
    (args.output / f"{args.mode}_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    results = []

    def compare(case):
        name, expected, path = case
        directory = args.output / args.mode / name
        directory.mkdir(parents=True, exist_ok=True)
        outputs = []
        for index, binary in enumerate(binaries):
            command = list(expected["command"])
            command[0] = str(binary)
            # DD is now the CLI default; these are ILP formulation comparisons.
            if "--decoder" in command:
                offset = command.index("--decoder")
                if ilp_decoder_args(binary):
                    command[offset + 1] = "ilp"
                else:
                    del command[offset:offset + 2]
            else:
                command[1:1] = ilp_decoder_args(binary)
            if "--bpseq-file" in command:
                command[command.index("--bpseq-file") + 1] = str(directory / f"{index}.bpseq")
            if args.engine:
                if "-e" in command:
                    command[command.index("-e") + 1] = args.engine
                else:
                    command[1:1] = ["-e", args.engine]
            outputs.append(run(command, 420, directory / f"{index}.log"))
        result = dict(case=name, pdb=expected["pdb"], condition=expected["condition"])
        if args.mode != "predictions":
            fingerprints = [re.findall(r"^MODEL .*$", output["output"], re.M) for output in outputs]
            result.update(equal=(all(output["status"] == "ok" for output in outputs) or
                                 args.mode == "fixtures" and outputs[0]["status"] == outputs[1]["status"]) and
                          fingerprints[0] == fingerprints[1], models=len(fingerprints[0]),
                          verification="shape_only" if args.mode == "formulations" else "ordered_model_hash",
                          fingerprints=fingerprints[0],
                          before_seconds=outputs[0]["seconds"], after_seconds=outputs[1]["seconds"],
                          before_build_seconds=sum(map(float, re.findall(r"built in ([0-9.]+)s", outputs[0]["output"]))),
                          after_build_seconds=sum(map(float, re.findall(r"built in ([0-9.]+)s", outputs[1]["output"]))))
        else:
            actual = outputs[0]
            pairs = []
            if actual["status"] == "ok":
                for line in (directory / "0.bpseq").read_text().splitlines():
                    fields = line.split()
                    if len(fields) == 3 and fields[0].isdigit() and fields[2].isdigit():
                        i, j = int(fields[0]), int(fields[2])
                        if i < j: pairs.append([i, j])
            structures = re.findall(r"^[.()\[\]{}<>]+$", actual["output"], re.M)
            reference_structures = re.findall(r"^[.()\[\]{}<>]+$", (path.parent / "run.log").read_text(), re.M)
            result.update(status=actual["status"], equal=actual["status"] == "ok" and
                          sorted(pairs) == expected["pairs"] and structures == reference_structures,
                          seconds=actual["seconds"], archived_seconds=expected["seconds"])
        return result

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        if args.repeats < 1:
            parser.error("--repeats must be positive")
        cases = [(f"{name}_r{repeat}" if args.repeats > 1 else name, row, path)
                 for name, row, path in cases for repeat in range(args.repeats)]
        futures = [pool.submit(compare, case) for case in cases]
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            print(f"{len(results)}/{len(cases)} {result['case']} equal={result['equal']}", flush=True)
            (args.output / f"{args.mode}.json").write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(dict(cases=len(results), equal=sum(r["equal"] for r in results),
                          models=sum(r.get("models", 0) for r in results)), indent=2))
    if any(not r["equal"] for r in results):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
