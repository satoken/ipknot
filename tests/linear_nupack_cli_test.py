"""Exercise LinearNUPACK factory aliases, sparse export, refinement and alignment."""
from pathlib import Path
import math
import subprocess
import sys
import tempfile

binary = str(Path(sys.argv[1]).resolve())
fixture = Path(__file__).resolve().parents[1] / "examples/drz_Ppac_1_1.fa"


def run(*args):
    result = subprocess.run([binary, "--decoder", "dd", "-t", "0.2,0.1", *args],
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
    return result.stdout


with tempfile.TemporaryDirectory() as temporary:
    root = Path(temporary)
    outputs = []
    for model in ("LinearNUPACK", "lnupack", "linearnupack"):
        bpp = root / f"{model}.bpp"
        outputs.append(run("-e", model, "--beam-size", "25", "-r", "0", "--bpp", str(bpp), str(fixture)))
        probabilities = []
        for line in bpp.read_text().splitlines():
            fields = line.split()
            if len(fields) >= 3 and fields[0].isdigit():
                for field in fields[2:]:
                    probability = float(field.split(":")[-1])
                    assert math.isfinite(probability) and 0 <= probability <= 1
                    probabilities.append(probability)
        assert probabilities, bpp.read_text()
    assert outputs[0] == outputs[1] == outputs[2]
    ambiguous = root / "ambiguous.fa"
    ambiguous.write_text(">ambiguous\nGGNYMRSWPOKVXIDHCC\n")
    for model in ("NUPACK", "LinearNUPACK", "lpc", "lpv"):
        bpseq = root / f"{model}-ambiguous.bpseq"
        run("-e", model, "-r", "0", "-B", str(bpseq), str(ambiguous))
        rows = [line.split() for line in bpseq.read_text().splitlines()
                if line and not line.startswith("#")]
        assert len(rows) == 18
        assert "".join(row[1] for row in rows) == "GGNYMRSWPOKVXIDHCC"
        for row in rows:
            if row[1] not in "ACGU":
                assert row[2] == "0", row
    # Constrained posteriors are used in refinement, including sparse output.
    assert run("-e", "LinearNUPACK", "--beam-size", "25", "-r", "1", str(fixture))
    # Exercise the alignment factory and sparse AveragedModel path.
    alignment = root / "alignment.aln"
    alignment.write_text("CLUSTAL W\n\nseq1 GGGAAACCC\nseq2 GGGAAACCC\n\n")
    assert run("-e", "LinearNUPACK", "--beam-size", "0", "-r", "0", str(alignment))
    missing = subprocess.run([binary, "--decoder", "dd", "-e", "LinearNUPACK",
                              "-P", str(root / "missing.par"), str(fixture)],
                             capture_output=True, text=True, timeout=60)
    assert missing.returncode != 0
    assert "Cannot load LinearNUPACK parameters" in missing.stderr + missing.stdout
