#!/usr/bin/env python3
"""Reproducible matrix-handoff and sparse-posterior time/RSS measurements."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shlex
import subprocess
from noe_build_helpers import legacy_ip_api, link_command

ROOT = Path(__file__).resolve().parent.parent


def compile_probe(source, output, build, object_name):
    obj = output.with_suffix(".o")
    # Matrix probes must use the ABI of the retained baseline's IP object.
    flags = []
    if object_name == "ip.cpp" and legacy_ip_api(build):
        flags = ["-include", str(ROOT / "tests/nmr_legacy_ip.h")]
    subprocess.run(["c++", "-std=c++17", "-O3", "-UNDEBUG", *flags, "-I", str(ROOT / "src"),
                    "-c", str(ROOT / "tests" / source), "-o", str(obj)], check=True)
    wanted = f"CMakeFiles/ipknot.dir/src/{object_name}.o"
    args = link_command(build)
    args = [arg for arg in args if not arg.startswith("-Wl,--dependency-file=") and
            (not arg.endswith(".o") or arg == wanted)]
    args.insert(1, str(obj))
    args[args.index("-o") + 1] = str(output)
    subprocess.run(args, cwd=build, check=True)


def measure(binary, arguments, log):
    result = subprocess.run(["/usr/bin/time", "-f", "RESOURCE %e %M", str(binary), *arguments],
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=60)
    log.write_text(result.stdout)
    if result.returncode:
        raise RuntimeError(f"Probe failed ({result.returncode}): {log}")
    resource = re.search(r"^RESOURCE ([0-9.]+) (\d+)$", result.stdout, re.M)
    return dict(seconds=float(resource[1]), max_rss_kib=int(resource[2]), output=result.stdout)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--before-build", type=Path, default=Path("/tmp/ipknot-nmr-followup-final-build"))
    parser.add_argument("--after-build", type=Path, default=Path("/tmp/ipknot-noe-optimized-build"))
    parser.add_argument("--output", type=Path, default=ROOT / "new_nmr/noe_optimization_results/sparse_performance")
    parser.add_argument("--mode", choices=("all", "matrix", "posterior"), default="all")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--lengths", default="2000,4000,8000")
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    results = []
    provenance = dict()
    for label, build in (("before", args.before_build), ("after", args.after_build)):
        for name in ("ip.cpp", "ipknot.cpp", "main.cpp", "fold.cpp", "linearpartition/LinearPartition.cpp"):
            path = build / f"CMakeFiles/ipknot.dir/src/{name}.o"
            provenance[f"{label}_{name}"] = hashlib.sha256(path.read_bytes()).hexdigest()
    (args.output / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")

    def save(record):
        results.append(record)
        print(json.dumps(record), flush=True)
        (args.output / "measurements.json").write_text(json.dumps(results, indent=2) + "\n")

    if args.mode in ("all", "matrix"):
        binaries = []
        for label, build in (("before", args.before_build), ("after", args.after_build)):
            binary = args.output / f"{label}_matrix"
            compile_probe("ip_matrix_transfer_test.cpp", binary, build, "ip.cpp")
            binaries.append((label, binary))
        for repeat in range(args.repeats):
            # Reverse paired order on alternate runs to reduce order bias.
            for label, binary in (binaries if repeat % 2 == 0 else reversed(binaries)):
                measured = measure(binary, ["10000", "20000", "1000"], args.output / f"{label}_matrix_{repeat}.log")
                save(dict(kind="matrix", version=label, repeat=repeat, columns=10000,
                          rows=20000, coefficients=10000000, **measured))
    if args.mode in ("all", "posterior"):
        binary = args.output / "sparse_posterior"
        compile_probe("linearpartition_sparse_benchmark.cpp", binary, args.after_build,
                      "linearpartition/LinearPartition.cpp")
        for engine in ("lpc", "lpv"):
            for length in map(int, args.lengths.split(",")):
                for repeat in range(args.repeats):
                    measured = measure(binary, [str(length), engine],
                                       args.output / f"{engine}_{length}_{repeat}.log")
                    fields = re.search(r"entries=(\d+) seconds=([0-9.]+)", measured["output"])
                    save(dict(kind="posterior", engine=engine, length=length, beam=100, repeat=repeat,
                              entries=int(fields[1]), posterior_seconds=float(fields[2]), **measured))


if __name__ == "__main__":
    main()
