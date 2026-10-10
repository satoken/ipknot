"""Independent metric and child-resource checks for the measurement harness."""
import gzip
import json
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest

from dataset import crossing_pairs, metrics
from report import failure


class EvaluationTests(unittest.TestCase):
    def test_crossings_against_brute_force(self):
        rng = random.Random(9)
        for n in range(2, 70):
            endpoints = list(range(1, n + 1))
            rng.shuffle(endpoints)
            pairs = [tuple(sorted(endpoints[i:i+2])) for i in range(0, n-1, 2)]
            expected = {p for p in pairs for q in pairs
                        if p[0] < q[0] < p[1] < q[1] or q[0] < p[0] < q[1] < p[1]}
            self.assertEqual(crossing_pairs(pairs), expected)

    def test_metrics(self):
        score = metrics([(1, 8), (2, 7), (3, 6)], [(1, 8), (2, 7), (4, 5)])
        self.assertEqual((score['tp'], score['fp'], score['fn']), (2, 1, 1))
        self.assertAlmostEqual(score['f1'], 2/3)
        self.assertEqual(metrics([], [])['f1'], 1)
        self.assertEqual(metrics([(1, 5)], [])['f1'], 0)

    def test_failure_classification(self):
        self.assertEqual(failure(dict(exit_code=137, timeout=True, stderr='')), 'timeout')
        self.assertEqual(failure(dict(exit_code=134, timeout=False, stderr='std::bad_alloc')), 'memory_limit')
        self.assertEqual(failure(dict(exit_code=1, timeout=False, stderr='unexpected')), 'other_error')

    def test_common_sets_and_audit(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            (root / 'final').mkdir()
            engines = ['nupack', 'lnupack', 'lpc', 'lpv']
            refs = [dict(id=i, member=i, sequence='GGAAACC', pairs=[[1,7],[2,6]],
                         dataset='Rfam14.5', length=7, size='short', pk=False,
                         input_sha256=i) for i in ['a', 'b', 'c', 'd']]
            with gzip.open(root / 'references.json.gz', 'wt') as f:
                json.dump(refs, f)
            (root / 'final-manifest.json').write_text(json.dumps(dict(sample_ids=['a','b','c','d'], engines=engines)))
            for ref in refs:
                for engine in engines:
                    row = dict(id=ref['id'], mode=engine, exit_code=0, timeout=False,
                               stderr='', input_sha256=ref['input_sha256'], wall_seconds=1,
                               user_seconds=.5, system_seconds=.1, rss_kib=100)
                    if (ref['id'], engine) == ('b', 'nupack'):
                        row.update(exit_code=137, timeout=True)
                    elif (ref['id'], engine) == ('c', 'lnupack'):
                        row.update(exit_code=134, stderr='std::bad_alloc')
                    else:
                        row.update(pairs=ref['pairs'], accuracy=metrics(ref['pairs'], ref['pairs']),
                                   crossing_accuracy=metrics([], []))
                    (root / 'final' / f"{ref['id']}-{engine}.json").write_text(json.dumps(row))
            command = [sys.executable, str(Path(__file__).with_name('report.py')),
                       '--root', str(root), '--output', str(root / 'out'), '--no-bootstrap']
            subprocess.run(command, check=True, capture_output=True, text=True)
            result = json.loads((root / 'out/summary.json').read_text())
            self.assertTrue(result['complete'])
            self.assertEqual(result['common_linear_n'], 3)
            self.assertEqual(result['common_all_four_n'], 2)
            self.assertEqual(result['groups']['all_attempts']['nupack']['macro_f1_failures_zero'], .75)
            # A wrong saved score must invalidate the report instead of being trusted.
            p = root / 'final/a-lpc.json'
            row = json.loads(p.read_text()); row['accuracy']['tp'] = 500
            p.write_text(json.dumps(row))
            self.assertNotEqual(subprocess.run(command, capture_output=True).returncode, 0)

    def test_launcher(self):
        launcher = os.environ.get('NUPACK_BENCHMARK_LAUNCHER')
        if not launcher:
            self.skipTest('set NUPACK_BENCHMARK_LAUNCHER to check wait4/time/memory limits')
        with tempfile.TemporaryDirectory() as d:
            result = Path(d) / 'resources'
            run = subprocess.run([launcher, str(result), '1', '1', '/bin/sleep', '5'])
            row = json.loads(result.read_text())
            self.assertEqual(run.returncode, 137)
            self.assertTrue(row['timeout'])
            self.assertTrue(.9 <= row['wall_seconds'] < 3)
            run = subprocess.run([launcher, str(result), '5', '1', sys.executable,
                                  '-c', 'bytearray(2 * 1024**3)'], capture_output=True)
            self.assertNotEqual(run.returncode, 0)
            self.assertIn(b'MemoryError', run.stderr)
            self.assertFalse(json.loads(result.read_text())['timeout'])


if __name__ == '__main__':
    unittest.main()
