"""Time existing BPP fixtures through the CLI, including process and file I/O."""
import argparse
import hashlib
import json
from pathlib import Path
import statistics
import subprocess
import tempfile
import time

parser = argparse.ArgumentParser()
parser.add_argument('--baseline', type=Path, required=True)
parser.add_argument('--candidate', type=Path, required=True)
parser.add_argument('--input-dir', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
parser.add_argument('--intra', action='store_true', help='Compare within-DP threads, including wavefront')
args = parser.parse_args()
files = sorted(args.input_dir.glob('*.bpp'))
if not files:
    parser.error('No .bpp files in input directory')
results = []
with tempfile.TemporaryDirectory() as temporary:
    trace = Path(temporary) / 'trace.jsonl'
    structure = Path(temporary) / 'output.bpseq'
    for dp in (('beam', 'improved-beam', 'nussinov', 'wavefront') if args.intra else
               ('beam', 'improved-beam', 'nussinov')):
        reference = None
        for revision, binary, threads in [('baseline', args.baseline, 1),
                                          ('candidate', args.candidate, 1),
                                          ('candidate', args.candidate, 4 if args.intra else 2)]:
            timings = []
            for repeat in range(6):
                start = time.perf_counter()
                summaries = []
                for file in files:
                    engine = 'nussinov' if revision == 'baseline' and dp == 'wavefront' else dp
                    command = [str(binary), '-x', '-r', '0', '-t', '.2,.1', '--dd-dp', engine,
                               '--dd-trace', str(trace), '--loglevel', 'error', '-B', str(structure)]
                    if revision == 'candidate':
                        command += ['--dd-dp-threads' if args.intra else '--dd-threads', str(threads)]
                    subprocess.run(command + [str(file)], capture_output=True, text=True,
                                   timeout=30, check=True)
                    events = [json.loads(line) for line in trace.read_text().splitlines()]
                    summary = next(event for event in reversed(events) if event['event'] == 'summary')
                    numeric_lines = '\n'.join(line for line in structure.read_text().splitlines()
                                             if line and line[0].isdigit())
                    summaries.append(dict(file=file.name, objective=summary['lower_bound'],
                                          upper_bound=summary['upper_bound'], iterations=summary['iterations'],
                                          structure_sha256=hashlib.sha256(numeric_lines.encode()).hexdigest()))
                if repeat:
                    timings.append(time.perf_counter() - start)
            if revision == 'candidate' and threads == 1:
                reference = summaries
            if threads > 1:
                assert reference == summaries, 'Parallel execution changed results'
            result = dict(dp=dp, revision=revision, threads=threads, files=len(files),
                          median_seconds=statistics.median(timings), summaries=summaries)
            if args.intra:
                result['thread_scope'] = 'dp'
            print({k: v for k, v in result.items() if k != 'summaries'}, flush=True)
            results.append(result)
args.output.write_text(json.dumps(results, indent=2) + '\n')
