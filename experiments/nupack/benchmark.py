"""Reproduce exact speed/identity and beam scaling; uses only Python's stdlib.

The reference receives the same correctness repairs (default parameter
alignment and single-pair pseudoknot marginals) before comparing performance.
"""
import argparse
import json
import os
from pathlib import Path
import random
import statistics
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--binary", type=Path, default=ROOT / "build/nupack_benchmark")
parser.add_argument("--reference", default="b0f965c")
parser.add_argument("--repeats", type=int, default=3)
parser.add_argument("--exact-lengths", type=int, nargs="+", default=[20, 40, 60])
parser.add_argument("--linear-lengths", type=int, nargs="+", default=[250, 500, 1000, 2000])
parser.add_argument("--beams", type=int, nargs="+", default=[25, 100])
parser.add_argument("--output", type=Path)
args = parser.parse_args()
args.binary = args.binary.resolve()


def run(binary, mode, seq, beam=100):
    results = []
    for _ in range(args.repeats):
        process = subprocess.run([str(binary), mode, str(beam), "dump"], input=seq + "\n",
                                 text=True, capture_output=True, check=True, timeout=300)
        results.append(json.loads(process.stdout))
    summary = results[0]
    summary["seconds"] = statistics.median(r["seconds"] for r in results)
    summary["rss_kib"] = max(r["rss_kib"] for r in results)
    assert all(r["log_z"] == summary["log_z"] and r["pairs"] == summary["pairs"] for r in results)
    return summary


def errors(reference, approximate):
    a = {(i, j): p for i, j, p in reference["pairs"]}
    b = {(i, j): p for i, j, p in approximate["pairs"]}
    differences = [abs(a.get(k, 0) - b.get(k, 0)) for k in a.keys() | b.keys()]
    n = reference["length"]
    return {"bpp_mae": sum(differences) / max(1, n * (n - 1) / 2),
            "bpp_max_error": max(differences, default=0),
            "log_z_error": approximate["log_z"] - reference["log_z"]}


report = {"reference": args.reference, "repeats": args.repeats, "seed": 20261010,
          "reference_correctness_repaired": True, "exact": [], "linear": []}
with tempfile.TemporaryDirectory(prefix="nupack-reference-") as directory:
    directory = Path(directory)
    for filename in ("nupack.cpp", "nupack.h", "dptable.h", "def_param.h"):
        source = subprocess.check_output(["git", "show", f"{args.reference}:src/nupack/{filename}"],
                                         cwd=ROOT, text=True)
        if filename == "nupack.cpp":
            old = "asymmetry_penalty[i] = *(v++)/100.0;"
            if old + "\n  max_asymmetry" not in source:
                source = source.replace(old, old + "\n  max_asymmetry = *(v++)/100.0;", 1)
            source = source.replace("Pg(i,a,d,e) += p;\n                Pg(b,c,f,j) += p;",
                                    "Pbg(i,e) += p;\n                Pbg(b,j) += p;")
            source = source.replace("Pg(i,i,e,f) += p;", "Pbg(i,f) += p;")
            source = source.replace("Pg(d,d,j,j) += p;", "Pbg(d,j) += p;")
        (directory / filename).write_text(source)
    reference = directory / "reference"
    subprocess.run([os.environ.get("CXX", "c++"), "-std=c++17", "-O3", "-DNDEBUG",
                    "-DNUPACK_REFERENCE", f"-I{directory}", str(ROOT / "experiments/nupack/driver.cpp"),
                    str(directory / "nupack.cpp"), "-o", str(reference)], check=True)
    for n in args.exact_lengths:
        rng = random.Random(20261010 + n)
        seq = "".join(rng.choices("ACGU", k=n))
        before, after = run(reference, "exact", seq), run(args.binary, "exact", seq)
        assert before["pairs"] == after["pairs"], "Exact float posteriors changed"
        assert before["log_z"] == after["log_z"], "Exact partition function changed"
        row = {"length": n, "before_seconds": before["seconds"], "after_seconds": after["seconds"],
               "speedup": before["seconds"] / after["seconds"], "bitwise_float_bpp_equal": True}
        row["approximation"] = []
        for beam in args.beams:
            approximate = run(args.binary, "linear", seq, beam)
            row["approximation"].append({"beam": beam, "seconds": approximate["seconds"],
                                          **errors(after, approximate)})
        report["exact"].append(row)
        print(json.dumps(row), flush=True)
    for n in args.linear_lengths:
        rng = random.Random(20261010 + n)
        seq = "".join(rng.choices("ACGU", k=n))
        for beam in args.beams:
            row = run(args.binary, "linear", seq, beam)
            row.pop("pairs")
            row["beam"] = beam
            report["linear"].append(row)
            print(json.dumps(row), flush=True)
if args.output:
    args.output.write_text(json.dumps(report, indent=2) + "\n")
