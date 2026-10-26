"""Static accuracy/scaling figures. Requires matplotlib; uses all saved runs."""
import argparse
import collections
import csv
import math
from pathlib import Path
import random
import statistics

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

LABELS = {'nupack': 'NUPACK (r=0)', 'lnupack': 'LinearNUPACK (b=100,r=0)',
          'lpc': 'LPC (b=100,r=1)', 'lpv': 'LPV (b=100,r=1)'}
COLORS = {'nupack': '#8b5d9b', 'lnupack': '#d45c25', 'lpc': '#2477b5', 'lpv': '#32955c'}
BINS = [0, 50, 75, 100, 125, 150, 200, 300, 500, 1000, 2000, 4500]


def binned(rows, key, average=statistics.median):
    x, y = [], []
    for a, b in zip(BINS, BINS[1:]):
        group = [r for r in rows if a < r['length'] <= b]
        if group:
            x.append(statistics.median(r['length'] for r in group))
            y.append(average(r[key] for r in group))
    return x, y


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    args.output.mkdir(exist_ok=True, parents=True)
    with args.input.open() as f:
        rows = list(csv.DictReader(f))
    for r in rows:
        r['length'] = int(r['length'])
        for key in ['wall_seconds', 'rss_mib', 'f1']:
            r[key] = float(r[key]) if r[key] else None
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    for ax, common, title in zip(axes, ['common_linear', 'common_all_four'],
                                ['Three linear engines: common successful RNAs', 'Four engines: common successful RNAs']):
        engines = ['lnupack', 'lpc', 'lpv'] if common == 'common_linear' else list(LABELS)
        for engine in engines:
            selected = [r for r in rows if r['engine'] == engine and r[common] == 'True' and r['status'] == 'success']
            x, y = binned(selected, 'f1', statistics.mean)
            ax.plot(x, y, '.-', color=COLORS[engine], label=LABELS[engine])
        ax.set(xscale='log', ylim=(0, 1), xlabel='Sequence length (nt)', ylabel='Macro F1', title=title)
        ax.grid(alpha=.2)
        ax.legend(fontsize=8)
    fig.savefig(args.output / 'accuracy.png', dpi=180)
    fig.savefig(args.output / 'accuracy.svg')
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    for engine in LABELS:
        good = [r for r in rows if r['engine'] == engine and r['status'] == 'success']
        sample = random.Random(20261010).sample(good, min(1200, len(good)))
        for ax, key in zip(axes, ['wall_seconds', 'rss_mib']):
            ax.scatter([r['length'] for r in sample], [r[key] for r in sample],
                       s=3, alpha=.09, color=COLORS[engine], rasterized=True)
            x, y = binned(good, key)
            ax.plot(x, y, '.-', color=COLORS[engine], label=LABELS[engine])
    for ax, label in zip(axes, ['Child wall time (s)', 'Peak RSS (MiB)']):
        ax.set(xscale='log', yscale='log', xlabel='Sequence length (nt)', ylabel=label)
        ax.grid(alpha=.2)
        ax.legend(fontsize=8)
    fig.suptitle('Successful runs only; lines = medians within length bins\nExact NUPACK: 60 s / 8 GiB limit; linear engines: 1800 s / 8 GiB')
    fig.savefig(args.output / 'scaling.png', dpi=180)
    fig.savefig(args.output / 'scaling.svg')
    plt.close(fig)


if __name__ == '__main__':
    main()
