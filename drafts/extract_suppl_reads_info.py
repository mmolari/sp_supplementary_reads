import pysam
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Extract supplementary mappings.')
parser.add_argument('--sam', help='SAM file name.')
parser.add_argument('--csv', help='CSV file name.')
args = parser.parse_args()
csv_file = args.csv
sam_file = args.sam

supplementary = []
with pysam.AlignmentFile(sam_file) as af:
    for read in af:
        if read.is_supplementary:
            supplementary.append(read.query_name)

df = []
with pysam.AlignmentFile(sam_file) as af:
    for read in af:
        if read.query_name in supplementary:
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