"""Exercise learned CLI wiring, signed decoding, neutral export and validation."""
from pathlib import Path
import subprocess
import sys
import tempfile

binary = sys.argv[1]
names = ("bias anchor_support weak_support weak_min_support support_product "
         "support_gap weak_specificity anchor_specificity_product min_length "
         "weak_length outer_loop middle_loop").split()


def model(path, bias):
    path.write_text("IPKNOT_PK_LINEAR_V1\n" + "".join(
        f"{name} {bias if name == 'bias' else 0}\n" for name in names))


def pairs(path):
    return {(int(f[0]), int(f[2])) for line in path.read_text().splitlines()
            if len(f := line.split()) == 3 and f[0].isdigit()
            and f[2].isdigit() and 0 < int(f[0]) < int(f[2])}


with tempfile.TemporaryDirectory() as temporary:
    root = Path(temporary)
    bpp, output, weights = root / "pairs.bpp", root / "out.bpseq", root / "model.txt"
    expected = {(1, 10), (2, 9), (5, 14), (6, 13)}
    evidence = dict.fromkeys({(1, 10), (2, 9)}, .8)
    evidence.update(dict.fromkeys({(5, 14), (6, 13)}, .45))
    sequence = ["A"] * 14
    for i, j in expected:
        sequence[i - 1], sequence[j - 1] = "G", "C"
    bpp.write_text("\n".join(f"{i} {base}" + "".join(
        f" {j}:{p}" for (left, j), p in sorted(evidence.items()) if left == i)
        for i, base in enumerate(sequence, 1)) + "\n")

    def run(extra):
        return subprocess.run([binary, "--decoder", "ilp", "-x", "-t", "0.4,0.4", "-B", str(output),
                               "--loglevel", "info", *extra, str(bpp)],
                              capture_output=True, text=True, timeout=15)

    result = run([])
    assert result.returncode == 0 and pairs(output) == expected, result.stdout + result.stderr
    baseline = output.read_bytes()
    model(weights, -1)
    flags = ["--pk-learned-model", str(weights), "--pk-selection-weight", "0"]
    result = run(flags)
    assert result.returncode == 0 and pairs(output) == {(1, 10), (2, 9)}, result.stdout + result.stderr
    assert "0 additional integer variables" in result.stdout + result.stderr
    result = run([*flags, "--pk-crossing-simplify"])
    assert result.returncode == 0 and pairs(output) == {(1, 10), (2, 9)}, result.stdout + result.stderr
    assert "PK crossing simplification:" in result.stdout + result.stderr
    projected = [*flags, "--pk-h-formulation", "projected"]
    result = run(projected)
    assert result.returncode == 0 and pairs(output) == {(1, 10), (2, 9)}, result.stdout + result.stderr
    assert "PK projection:" in result.stdout + result.stderr
    assert "0 additional variables, 0 additional rows" in result.stdout + result.stderr
    result = run([*projected, "--pk-hybrid-shape", "--pk-h-intercept", "3"])
    assert result.returncode == 0 and pairs(output) == expected, result.stdout + result.stderr
    result = run([*projected, "--pk-learned-scale", "0"])
    assert result.returncode == 0 and output.read_bytes() == baseline, result.stdout + result.stderr
    for extra in [["--pk-crossing-hypograph"], ["--pk-crossing-normalize"], ["--pk-crossing-tight-bounds"],
                  ["--pk-crossing-normalize", "--pk-crossing-tight-bounds"]]:
        result = run([*flags, *extra])
        assert result.returncode == 0 and pairs(output) == {(1, 10), (2, 9)}, result.stdout + result.stderr
    # Both components contribute to one contact score, with independent scales.
    result = run([*flags, "--pk-hybrid-shape", "--pk-h-intercept", "3"])
    assert result.returncode == 0 and pairs(output) == expected, result.stdout + result.stderr
    result = run([*flags, "--pk-hybrid-shape", "--pk-h-intercept", "3",
                  "--pk-h-weight", "0"])
    assert result.returncode == 0 and pairs(output) == {(1, 10), (2, 9)}, result.stdout + result.stderr
    result = run([*flags, "--pk-hybrid-shape", "--pk-h-intercept", "-3",
                  "--pk-learned-scale", "0"])
    assert result.returncode == 0 and pairs(output) == {(1, 10), (2, 9)}, result.stdout + result.stderr
    result = run([*flags, "--pk-learned-scale", "0"])
    assert result.returncode == 0 and output.read_bytes() == baseline, result.stdout + result.stderr
    model(weights, 0)
    result = run(flags)
    assert result.returncode == 0 and output.read_bytes() == baseline, result.stdout + result.stderr
    features = root / "features.tsv"
    result = run(["--pk-feature-output", str(features)])
    assert result.returncode == 0 and output.read_bytes() == baseline, result.stdout + result.stderr
    rows = features.read_text().splitlines()
    assert any(line.startswith("G\t") for line in rows) and any(line.startswith("B\t") for line in rows)
    # Large corrections cannot introduce pairs discarded by the original cut.
    model(weights, 100)
    bpp.write_text(bpp.read_text().replace(":0.45", ":0.3"))
    result = run(flags)
    assert result.returncode == 0 and pairs(output) == {(1, 10), (2, 9)}, result.stdout + result.stderr
    result = run(projected)
    assert result.returncode == 0 and pairs(output) == {(1, 10), (2, 9)}, result.stdout + result.stderr
    for invalid in [[*flags, "--pk-learned-scale", "nan"],
                    [*flags, "--pk-learned-scale", "-1"],
                    [*flags, "--pk-energy-model", "dp"],
                    [*flags, "--pk-h-intercept", "1"],
                    ["--pk-feature-output", str(features), "--pk-h-formulation", "projected"],
                    [*flags, "--pk-h-allocation", "substems"],
                    ["--pk-hybrid-shape"],
                    ["--pk-crossing-normalize", "--pk-h-formulation", "projected"],
                    ["--pk-crossing-tight-bounds", "--pk-h-formulation", "projected"],
                    ["--pk-crossing-hypograph", "--pk-h-formulation", "projected"],
                    ["--pk-crossing-simplify", "--pk-h-formulation", "projected"]]:
        assert run(invalid).returncode != 0, invalid
    weights.write_text("IPKNOT_PK_LINEAR_V1\nbias 1\n")
    assert run(flags).returncode != 0, "incomplete model accepted"
    bpp.write_text(bpp.read_text().replace(":0.3", ":0.45"))
    for decoder in ('ilp','dd'):
        for formulation in ('crossing','projected'):
            for cost in ('log','cc'):
                weights.write_text('IPKNOT_PK_BOUNDED_V1\nbias 1\nsupport 0\ncompetition 0\nloop_cost .01\ncap 8\nloop_model '+cost+'\n')
                command=[binary,'--decoder',decoder,'-x','-r','0','-t','0.4,0.4','-B',str(output),
                         '--pk-learned-model',str(weights),'--pk-h-formulation',formulation,
                         '--pk-h-allocation','blocks',str(bpp)]
                result=subprocess.run(command,capture_output=True,text=True,timeout=15)
                assert result.returncode == 0 and pairs(output)=={(1,10),(2,9)}, result.stderr
                weights.write_text('IPKNOT_PK_BOUNDED_V1\nbias 0\nsupport 100\ncompetition 0\nloop_cost 0\ncap 8\nloop_model '+cost+'\n')
                result=subprocess.run(command,capture_output=True,text=True,timeout=15)
                assert result.returncode == 0 and pairs(output)==expected, result.stderr
                weights.write_text(weights.read_text().replace('cap 8','cap 0'))
                result=subprocess.run(command,capture_output=True,text=True,timeout=15)
                assert result.returncode == 0 and output.read_bytes()==baseline, result.stderr

print("Verified learned CLI, negative correction, zero compatibility, neutral feature export and candidate preservation")
