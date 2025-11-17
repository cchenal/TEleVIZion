#!/usr/bin/env python3
"""
Make a BED file of GC content windows colored by a Spectral colormap (hex colors).

Output columns (tab-separated): chr  start  end  name  itemRgb
- chr: sequence ID from FASTA
- start, end: 0-based half-open coordinates
- name: unique label (gc1, gc2, ...)
- itemRgb: hex color (e.g. #E41A1C) mapped to GC fraction

Usage:
  python gc_bed_hex.py input.fasta --window 100000 --out gc_windows.bed
"""

import argparse
import sys
import seaborn as sns
from matplotlib.colors import to_hex   # <-- use this instead of sns.utils.rgb2hex

def parse_args():
    p = argparse.ArgumentParser(description="Create GC-content colored BED (non-overlapping windows, hex colors).")
    p.add_argument("fasta", help="Input FASTA file (use '-' for stdin).")
    p.add_argument("--window", type=int, default=100000, help="Window size (bp). Default: 100000")
    p.add_argument("--out", default="-", help="Output BED path (default: stdout)")
    p.add_argument("--min-bases", type=int, default=1,
                   help="Min A/C/G/T bases required to compute GC%% (default: 1).")
    return p.parse_args()

def fasta_iter(path):
    """Simple FASTA reader yielding (seq_id, sequence_string)."""
    name = None
    chunks = []
    fh = sys.stdin if path == "-" else open(path, "rt")
    with fh:
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if name is not None:
            yield name, "".join(chunks)

def gc_fraction(seq, min_bases=1):
    """Return GC fraction in [0,1] ignoring non-ACGT; None if <min_bases informative bases."""
    gc = at = 0
    for b in seq:
        if b in "GgCc":
            gc += 1
        elif b in "AaTt":
            at += 1
    denom = gc + at
    if denom < min_bases:
        return None
    return gc / max(denom, 1)

def frac_to_hex(frac, cmap):
    """Map fraction in [0,1] to hex color. None → light gray."""
    if frac is None:
        return "#C8C8C8"  # light gray
    rgba = cmap(frac)
    return to_hex(rgba[:3])  # <-- fixed line

def main():
    args = parse_args()
    cmap = sns.color_palette("Spectral", as_cmap=True)  # or "Spectral_r" for reversed

    out_f = sys.stdout if args.out == "-" else open(args.out, "wt")
    with out_f as out:
        # Header line
        out.write("chr\tstart\tend\tname\titemRgb\tgc_content\n")

        uid = 0
        for chrom, seq in fasta_iter(args.fasta):
            L = len(seq)
            for start in range(0, L, args.window):  # NON-OVERLAPPING
                end = min(start + args.window, L)
                frac = gc_fraction(seq[start:end], min_bases=args.min_bases)
                color = frac_to_hex(frac, cmap)
                uid += 1
                out.write(f"{chrom}\t{start}\t{end}\tgc{uid}\t{color}\t{frac}\n")

if __name__ == "__main__":
    main()
