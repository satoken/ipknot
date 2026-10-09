"""End-to-end candidate retention, activation, table and CLI validation."""
from pathlib import Path
import subprocess
import sys
import tempfile

binary = sys.argv[1]


def run(extra, bpp, output):
    return subprocess.run([binary, "--decoder", "ilp", "-x", "-t", "0.4,0.4", "-B", str(output),
                           "--loglevel", "info", *extra, str(bpp)],
                          capture_output=True, text=True, timeout=15)


def pairs(path):
    result = set()
    for line in path.read_text().splitlines():
        fields = line.split()
        if len(fields) >= 3 and fields[0].isdigit() and fields[2].isdigit():
            i, j = int(fields[0]), int(fields[2])
            if i < j:
                result.add((i, j))
    return result


with tempfile.TemporaryDirectory() as temporary:
    root = Path(temporary)
    bpp = root / "fixture.bpp"
    output = root / "fixture.bpseq"
    expected = {(1, 10), (2, 9), (5, 14), (6, 13)}
    sequence = ["A"] * 14
    for i, j in expected:
        sequence[i - 1] = "G"
        sequence[j - 1] = "C"
    bpp.write_text("\n".join(f"{i} {base}" + "".join(f" {j}:0.3" for left, j in sorted(expected) if left == i)
                             for i, base in enumerate(sequence, 1)) + "\n")
    result = run([], bpp, output)
    assert result.returncode == 0 and pairs(output) == set(), result.stdout + result.stderr
    result = run(["--pk-h-intercept", "1"], bpp, output)
    assert result.returncode == 0 and pairs(output) == set(), "score resurrected a discarded candidate"
    retained = ["--pk-candidate-threshold", "0.1"]
    result = run([*retained, "--pk-h-intercept", "1"], bpp, output)
    assert result.returncode == 0 and pairs(output) == expected, result.stdout + result.stderr
    assert "PK score: 1" in result.stdout + result.stderr
    result = run([*retained, "--pk-h-intercept", "-1"], bpp, output)
    assert result.returncode == 0 and pairs(output) != expected, "negative cost was evaded"
    result = run([*retained, "--pk-h-intercept", "1", "--pk-h-formulation", "crossing"], bpp, output)
    assert result.returncode == 0 and pairs(output) == expected, result.stdout + result.stderr
    assert "PK score: 1" in result.stdout + result.stderr
    assert "0 additional integer variables" in result.stdout + result.stderr
    result = run([*retained, "--pk-h-intercept", "1", "--pk-h-formulation", "crossing", "--pk-h-weight", "0"], bpp, output)
    assert result.returncode == 0 and pairs(output) == set(), "zero crossing weight changed the original model"
    table = root / "scores.tsv"
    table.write_text("# stem1 stem2 loop1 loop2 loop3 score\n2 2 2 2 2 1\n")
    result = run([*retained, "--pk-h-table", str(table)], bpp, output)
    assert result.returncode == 0 and pairs(output) == expected, result.stdout + result.stderr
    # Automatic threshold selection must include the realized PK score.
    result = subprocess.run([binary, "--decoder", "ilp", "-x", "-t", "0.4_0.5,0.4_0.5", "-B", str(output),
                             "--loglevel", "info", *retained, "--pk-h-table", str(table), str(bpp)],
                            capture_output=True, text=True, timeout=15)
    assert result.returncode == 0 and pairs(output) == expected, result.stdout + result.stderr
    assert "Selected PK score=1" in result.stdout + result.stderr
    result = run([*retained, "--pk-h-table", str(table), "--pk-h-weight", "0"], bpp, output)
    assert result.returncode == 0 and pairs(output) == set(), "zero H weight did not disable scoring"
    for invalid in [["--pk-h-max-loop", "-1"], ["--pk-h-loop-mode", "bad"],
                    ["--pk-h-formulation", "bad"], ["--pk-candidate-threshold", "2"],
                    ["--pk-h-intercept", "nan"], ["--pk-h-weight", "-1"],
                    ["--pk-h-formulation", "projected", "--no-levelwise"],
                    ["--pk-h-formulation", "supported", "--no-levelwise"],
                    ["--pk-h-formulation", "crossing", "--no-levelwise"],
                    ["--pk-support-tolerance", "-0.1"]]:
        assert run(invalid, bpp, output).returncode != 0, invalid
    table.write_text("2 2 2 2 2 1\n2 2 2 2 2 -1\n")
    assert run(["--pk-h-table", str(table)], bpp, output).returncode != 0, "duplicate table row accepted"
    # Energy models default to the inexpensive fixed-block crossing decoder.
    energy_flags = [*retained, "--pk-energy-model", "dp", "--pk-energy-intercept", "0.5",
                    "--pk-energy-scale", "0.01"]
    result = run(energy_flags, bpp, output)
    assert result.returncode == 0 and pairs(output) == expected, result.stdout + result.stderr
    assert "PK blocks:" in result.stdout + result.stderr
    result = run([*retained, "--pk-energy-model", "dp", "--pk-energy-scale", "0.01"], bpp, output)
    assert result.returncode == 0 and pairs(output) == set(), "DP loop cost was not applied"
    result = run(["--pk-energy-model", "dp", "--pk-energy-scale", "0"], bpp, output)
    assert result.returncode == 0 and pairs(output) == set(), "zero energy weight changed baseline"
    for invalid in [["--pk-energy-model", "bad"], ["--pk-energy-scale", "0.01"],
                    ["--pk-energy-model", "dp", "--pk-energy-scale", "-0.1"],
                    ["--pk-energy-model", "cc", "--pk-energy-scale", "nan"],
                    ["--pk-energy-model", "dp", "--pk-energy-temperature", "-273.15"],
                    ["--pk-energy-model", "dp", "--pk-h-intercept", "1"],
                    ["--pk-h-allocation", "bad"], ["--pk-h-allocation", "blocks"],
                    ["--pk-energy-model", "dp", "--pk-energy-table", str(table)]]:
        assert run(invalid, bpp, output).returncode != 0, invalid
    # Projected coefficients leave model dimensions unchanged. They are a
    # look-ahead surrogate; the fixed underlying model still chooses a PK.
    result = run([*retained, "--pk-h-intercept", "1", "--pk-h-formulation", "projected"], bpp, output)
    assert result.returncode == 0 and pairs(output) == expected, result.stdout + result.stderr
    assert "0 additional variables, 0 additional rows" in result.stdout + result.stderr
    # Reranking does not change the fixed-threshold ILP optimum.
    result = run([*retained, "--pk-h-intercept", "1", "--pk-h-formulation", "rerank"], bpp, output)
    assert result.returncode == 0 and pairs(output) == set(), result.stdout + result.stderr
    # Its signed, realized score can change which automatic-threshold
    # candidate is chosen, without changing any candidate's ILP objective.
    result = subprocess.run([binary, "--decoder", "ilp", "-x", "-t", "0.1_0.4,0.1_0.4", "-B", str(output),
                             "--pk-h-intercept", "-1", "--pk-h-formulation", "rerank", str(bpp)],
                            capture_output=True, text=True, timeout=15)
    assert result.returncode == 0 and pairs(output) == set(), result.stdout + result.stderr
    # The motif budget is a visible error rather than silent score pruning.
    crowded = {(1, 11), (2, 10), (3, 9), (6, 16), (7, 15), (8, 14)}
    sequence = ["A"] * 16
    for i, j in crowded:
        sequence[i - 1] = "G"
        sequence[j - 1] = "C"
    bpp.write_text("\n".join(f"{i} {base}" + "".join(f" {j}:0.3" for left, j in sorted(crowded) if left == i)
                             for i, base in enumerate(sequence, 1)) + "\n")
    result = run([*retained, "--pk-h-intercept", "1", "--pk-h-max-motifs", "1"], bpp, output)
    assert result.returncode != 0 and "budget exceeded" in result.stdout + result.stderr, result.stdout + result.stderr
    # A-B offers the best projected shape, but BPP prefers A-C. Projection
    # can use the A-B bonus with the selected C partner; supported mode must
    # select a crossing partner compatible with that bonus.
    a = {(1, 10), (2, 9)}
    b = {(5, 14), (6, 13)}
    c = {(3, 12), (4, 11)}
    probability = {**dict.fromkeys(a, 0.55), **dict.fromkeys(b, 0.41), **dict.fromkeys(c, 0.75)}
    sequence = ["A"] * 14
    for i, j in a | b | c:
        sequence[i - 1] = "G"
        sequence[j - 1] = "C"
    bpp.write_text("\n".join(f"{i} {base}" + "".join(f" {j}:{p}" for (left, j), p in sorted(probability.items()) if left == i)
                             for i, base in enumerate(sequence, 1)) + "\n")
    table.write_text("2 2 2 2 2 1\n2 2 0 4 0 0.2\n")
    result = run(["--pk-h-table", str(table), "--pk-h-formulation", "projected"], bpp, output)
    assert result.returncode == 0 and pairs(output) == a | c, result.stdout + result.stderr
    result = run(["--pk-h-table", str(table), "--pk-h-formulation", "crossing"], bpp, output)
    assert result.returncode == 0 and pairs(output) == a | b, result.stdout + result.stderr
    result = run(["--pk-h-table", str(table), "--pk-h-formulation", "supported"], bpp, output)
    assert result.returncode == 0 and pairs(output) == a | b, result.stdout + result.stderr
    assert "tightened existing rows" in result.stdout + result.stderr
    assert "0 additional variables, 0 additional rows" in result.stdout + result.stderr
    result = run(["--pk-h-table", str(table), "--pk-h-formulation", "supported", "--pk-support-tolerance", "0.4"], bpp, output)
    assert result.returncode == 0 and pairs(output) == a | c, result.stdout + result.stderr
    result = run(["--pk-h-table", str(table), "--pk-h-formulation", "supported", "--pk-h-weight", "0"], bpp, output)
    assert result.returncode == 0 and pairs(output) == a | c, result.stdout + result.stderr

    # A supplied CC09 loop entropy changes the ILP while an unsupported
    # shape falls back to DP. The direct table has no CC06 assembly term.
    energy_pairs = {(1, 12), (2, 11), (3, 10), (5, 16), (6, 15), (7, 14)}
    sequence = ["A"] * 16
    for i, j in energy_pairs:
        sequence[i - 1], sequence[j - 1] = "G", "C"
    bpp.write_text("\n".join(f"{i} {base}" + "".join(f" {j}:0.3" for left, j in sorted(energy_pairs) if left == i)
                             for i, base in enumerate(sequence, 1)) + "\n")
    cc_table = root / "cc09.tsv"
    cc_table.write_text("Q 3 3 1 2 1 5\n")
    energy_flags = [*retained, "--pk-energy-model", "cc", "--pk-energy-intercept", "0.45",
                    "--pk-energy-scale", "0.02"]
    result = run(energy_flags, bpp, output)
    assert result.returncode == 0 and pairs(output) == set(), "missing CC parameter did not use DP"
    result = run([*energy_flags, "--pk-energy-table", str(cc_table)], bpp, output)
    assert result.returncode == 0 and pairs(output) == energy_pairs, result.stdout + result.stderr
    assert "PK score: 0.35" in result.stdout + result.stderr
    cc_table.write_text("Q 3 3 1 2 1 5\nQ 3 3 1 2 1 6\n")
    assert run([*energy_flags, "--pk-energy-table", str(cc_table)], bpp, output).returncode != 0

print("Verified candidate cut, signed scores, energy models, tables, automatic selection, projection, compatible witnesses and invalid options")
