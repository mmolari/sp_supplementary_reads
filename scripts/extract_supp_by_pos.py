import numpy as np
import pandas as pd
import pysam
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Write csv with normalized number of supplementary reads.')
parser.add_argument('--sam', help='SAM input file name.')
parser.add_argument('--csv', help='CSV output file name.')
parser.add_argument('--header', help='CSV output file header.')
parser.add_argument('--step', help='Histogram step size (default = 10000)')
args = parser.parse_args()
step = 10000 if args.step is None else args.step
bins = 4643630 / step # use genome size instead

start_prim = 0
start_suppl = np.zeros(int(np.ceil(bins)))
with pysam.AlignmentFile(args.sam) as af:
    for read in af:
        if read.is_supplementary:
            start_suppl[read.reference_start // step] += 1
        elif not read.is_secondary:
            start_prim += 1

# Expected primaries per bin (faster with pandas but useful for different normalization)
exp_prim = start_prim / bins

df = pd.DataFrame({args.header : start_suppl/exp_prim})
df.to_csv(args.csv)
