import numpy as np
import pysam
import pandas as pd
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Extract supplementary mappings.')
parser.add_argument('--sam', help='SAM input file name.')
parser.add_argument('--csv', help='CSV output file name.')
parser.add_argument('--png', help='PNG output file name for primary vs supplementary start positions.')
args = parser.parse_args()

def get_supp_query_list(sam_file):
    supplementary = []
    with pysam.AlignmentFile(args.sam) as af:
        for read in af:
            if read.is_supplementary:
                supplementary.append(read.query_name)
    return supplementary

def write_to_csv(sam_file, csv_file, query_list):
    df = []
    with pysam.AlignmentFile(sam_file) as af:
        for read in af:
            if read.query_name in query_list:
                mapping_type = 'primary'
                if read.is_secondary:
                    mapping_type = 'secondary'
                elif read.is_supplementary:
                    mapping_type = 'supplementary'
                df.append({
                    "query_name" : read.query_name,
                    "is_forward" : read.is_forward,
                    "type" : mapping_type,
                    "ref_start" : read.reference_start,
                    "ref_end" : read.reference_end,
                    "query_start" : read.query_alignment_start,
                    "query_end" : read.query_alignment_end,
                    "alignment_length" : read.query_alignment_length,
                    "query_length" : read.infer_read_length(),
                    "matches" : read.get_cigar_stats()[0][0],
                    "insertions" : read.get_cigar_stats()[0][1],
                    "deletions" : read.get_cigar_stats()[0][2],
                    "soft_clips" : read.get_cigar_stats()[0][4],
                    "hard_clips" : read.get_cigar_stats()[0][5],                
                })
    df = pd.DataFrame(df)
    df.to_csv(csv_file)
    return df

def get_query_info_from_name(query_list, df):
    query_ids = []
    for s in query_list:
        prim = 0
        sec = []
        supp = []
        for row in df.itertuples():
            if row.query_name == s:
                if row.type == "secondary":
                    sec.append(row.ref_start)
                elif row.type == "supplementary":
                    supp.append(row.ref_start)
                else:
                    prim = row.ref_start
        query_ids.append({
            "query_name" : s,
            "primary" : prim,
            "secondary" : sec,
            "supplementary" : supp
        })

def plot_prim_suppl_start(query_list):
    for q in query_list:
        for s in q["supplementary"]:
            p = q["primary"]
            color = "r" if (np.abs(s-p) > 1e5 and np.abs(s-p) < 4e6) else "k"
            plt.scatter(q["primary"],s,c=color,alpha=0.1)
    plt.xlabel("primary start")
    plt.ylabel("supplementary start")
    plt.savefig(args.png)

supplementary = get_supp_query_list(args.sam)
df = write_to_csv(args.sam, args.csv, supplementary)
plot_prim_suppl_start(get_query_info_from_name(supplementary, df))
