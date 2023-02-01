import numpy as np
import pysam
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Plot histogram of read positions.')
parser.add_argument('--sam', help='SAM input file name.')
parser.add_argument('--png', help='PNG output file name.')
parser.add_argument('--step', help='Histogram step size (default = 10000)')
args = parser.parse_args()
step = 10000 if args.step is None else args.step
bins = np.arange(0, 4643630 + 2*step, step) # use genome size instead

start_prim = list()
start_supp = list()
with pysam.AlignmentFile(args.sam) as af:
    for read in af:
        if read.is_supplementary:
            start_supp.append(read.reference_start)
        elif not read.is_secondary:
            start_prim.append(read.reference_start)

plt.hist(start_prim,bins=bins,color='k',label='primary')
plt.hist(start_supp,bins=bins,color='r',label='supplementary')
plt.xlabel('reference position')
plt.ylabel('number of reads')
plt.title('Reference position of mapping start')
plt.legend()
plt.savefig(args.png)
