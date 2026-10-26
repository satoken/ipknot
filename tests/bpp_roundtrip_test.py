"""A saved initial posterior must preserve native decoder evidence exactly."""
from pathlib import Path
import json
import subprocess
import sys
import tempfile

binary = sys.argv[1]
fixture = Path(__file__).resolve().parents[1] / "examples/drz_Ppac_1_1.fa"


def pairs(path):
    result = set()
    for line in path.read_text().splitlines():
        f = line.split()
        if len(f) == 3 and f[0].isdigit() and f[2].isdigit():
            if 0 < int(f[0]) < int(f[2]):
                result.add((int(f[0]), int(f[2])))
    return result


with tempfile.TemporaryDirectory() as temporary:
    root = Path(temporary)
    raw, cache = root / "raw.bpp", root / "cache.bpp"
    direct, cached = root / "direct.bpseq", root / "cached.bpseq"
    common = ["-r", "0", "-t", "auto,auto"]
    first = subprocess.run([binary, "--decoder", "ilp", *common, "-B", str(direct), "--bpp", str(raw),
                            str(fixture)],
                           capture_output=True, text=True, timeout=20)
    assert first.returncode == 0, first.stdout + first.stderr
    cache.write_text("\n".join(line for line in raw.read_text().splitlines()
                              if not line.startswith("#")) + "\n")
    second = subprocess.run([binary, "--decoder", "ilp", *common, "-x", "-B", str(cached),
                             str(cache)],
                            capture_output=True, text=True, timeout=20)
    assert second.returncode == 0, second.stdout + second.stderr
    assert pairs(direct) == pairs(cached), "Cached BPP changed R0 prediction"
    # The cached-input path does not re-export BPP. Instead compare DD's
    # actual candidate coordinates and float-derived objective coefficients.
    graphs, predictions = [], []
    for name, source, input_flags in (("native", fixture, []), ("cached", cache, ["-x"])):
        trace, prediction = root / f"{name}.jsonl", root / f"{name}.bpseq"
        result = subprocess.run([binary, "--decoder", "dd", *common,
                                 *input_flags, "--dd-trace", str(trace),
                                 "--dd-trace-state", "-B", str(prediction), str(source)],
                                capture_output=True, text=True, timeout=20)
        assert result.returncode == 0, result.stdout + result.stderr
        graphs.append([event["pairs"] for event in map(json.loads, trace.read_text().splitlines())
                       if event["event"] == "problem"])
        predictions.append(pairs(prediction))
    assert graphs[0] and graphs[0] == graphs[1], "BPP export/import changed decoder evidence"
    assert predictions[0] == predictions[1], "Cached BPP changed DD R0 prediction"

print("Verified exact saved-posterior and prediction roundtrip without refinement")
