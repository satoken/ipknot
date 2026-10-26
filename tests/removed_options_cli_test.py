"""Removed observation and experimental scoring options must fail explicitly."""
from pathlib import Path
import re
import subprocess
import sys

binary = sys.argv[1]
fixture = Path(__file__).resolve().parents[1] / 'examples/drz_Ppac_1_1.fa'
removed = [
    "dd-noe-solver", "dd-nmr-pair-beam", "dd-nmr-witnesses",
    "dd-nmr-pattern-beam", "dd-nmr-polish-states", "stack-constraint",
    "base-pairs", "without-canonical-neighbor", "coaxial-stacking",
    "nmr-coaxial-no-adjacent-bulge", "without-nmr-bulge", "nmr-bulge-mode",
    "nmr-soft", "nmr-count-penalty", "nmr-stack-penalty",
    "nmr-threshold-penalty-scale", "nmr-count-mode", "nmr-allow-shared-stack-pairs",
    "pk-level-penalty",
    "pk-h-intercept",
    "pk-h-weight",
    "pk-h-loop-penalty",
    "pk-h-stem-reward",
    "pk-h-coax-bonus",
    "pk-h-table",
    "pk-energy-model",
    "pk-energy-scale",
    "pk-energy-intercept",
    "pk-energy-temperature",
    "pk-energy-table",
    "pk-learned-model",
    "pk-learned-scale",
    "pk-core-width",
    "pk-best-partner",
    "pk-rank-model",
    "pk-rank-scale",
    "pk-rank-output",
    "pk-hybrid-shape",
    "pk-crossing-simplify",
    "pk-crossing-normalize",
    "pk-crossing-tight-bounds",
    "pk-crossing-hypograph",
    "pk-feature-output",
    "pk-h-max-stem",
    "pk-h-max-loop",
    "pk-h-max-motifs",
    "pk-h-loop-mode",
    "pk-h-formulation",
    "pk-h-allocation",
    "pk-support-tolerance",
    "pk-support-complete-stems",
    "pk-candidate-threshold",
    "pk-ensemble",
    "pk-ensemble-scale",
    "pk-ensemble-intercept",
    "pk-ensemble-temperature",
    "pk-ensemble-threshold",
    "pk-ensemble-output",
    "pk-selection-weight"
]

help_result = subprocess.run([binary, '--help'], capture_output=True, text=True, timeout=10)
assert help_result.returncode == 0, help_result.stdout + help_result.stderr
assert not re.search(r'--(?:pk-|nmr-|dd-nmr-|dd-noe-)', help_result.stdout + help_result.stderr)
for option in removed:
    result = subprocess.run([binary, f'--{option}=0', str(fixture)],
                            capture_output=True, text=True, timeout=10)
    text = result.stdout + result.stderr
    assert result.returncode != 0, (option, text)
    assert option in text and 'not exist' in text, (option, text)
print(f'{len(removed)} removed options explicitly rejected')
