"""A saved initial posterior must preserve native decoder evidence exactly."""
from pathlib import Path
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
    direct_features, cached_features = root / "direct.tsv", root / "cached.tsv"
    common = ["-r", "0", "-t", "auto,auto", "--pk-h-max-stem", "200",
              "--pk-h-max-loop", "200"]
    first = subprocess.run([binary, "--decoder", "ilp", *common, "-B", str(direct), "--bpp", str(raw),
                            "--pk-feature-output", str(direct_features), str(fixture)],
                           capture_output=True, text=True, timeout=20)
    assert first.returncode == 0, first.stdout + first.stderr
    cache.write_text("\n".join(line for line in raw.read_text().splitlines()
                              if not line.startswith("#")) + "\n")
    second = subprocess.run([binary, "--decoder", "ilp", *common, "-x", "-B", str(cached),
                             "--pk-feature-output", str(cached_features), str(cache)],
                            capture_output=True, text=True, timeout=20)
    assert second.returncode == 0, second.stdout + second.stderr
    assert pairs(direct) == pairs(cached), "Cached BPP changed R0 prediction"
    assert any(line.startswith("B\t") for line in direct_features.read_text().splitlines()), \
        "Fixture does not exercise posterior-dependent crossing features"
    assert direct_features.read_bytes() == cached_features.read_bytes(), \
        "BPP export rounded native evidence and changed crossing features"

print("Verified exact posterior-feature and prediction roundtrip without refinement")
