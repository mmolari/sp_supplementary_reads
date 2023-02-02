import pandas as pd
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Plot trajectories.')
parser.add_argument('--csv', help='CSV input file name')
parser.add_argument('--png', help='PNG output file name')
parser.add_argument('--n', help='Number of trajectories (default = 5)')
args = parser.parse_args()

n = 5 if args.n is None else int(args.n)
df = pd.read_csv(args.csv, index_col=0)

# Get times from csv
times = []
for col in df.columns:
    headers = col.rsplit('_',-1)
    for h in range(len(headers)):
        if headers[h] == 'time':
            times.append(headers[h-1])
for t in range(len(times)):
    times[t] = int(times[t])

# Create DataFrame with n positions of highest variability over time
diff = pd.DataFrame({'diff' : (df.max(axis=1)).sub(df.min(axis=1))})
df = pd.concat([df, diff], axis=1)
df = df[1:-1] # remove first and last position
df.to_csv("test.csv")
df = df.sort_values(by='diff', ascending=False)
traj = df.head(n)

# Plot trajectories
for row in traj.itertuples():
    y = []
    for t in range(len(times)):
        y.append(row[t+1])
    plt.plot(times,y,'-x',label='{}0 kbp'.format(row[0]))
plt.legend()
plt.xlabel('time')
plt.ylabel('supplementary mappings per expected primary')
plt.savefig(args.png)
