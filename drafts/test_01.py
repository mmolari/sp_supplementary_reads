# %%

import pysam

# %%

sam_file = "../results/reads.sam"
# sam_file = "../results/reads.sorted.bam"

# open sam/bam file
with pysam.AlignmentFile(sam_file) as af:

    print("description = ", af.description)
    print("n. references = ", af.nreferences)
    print("references = ", af.references)
    print("lengths = ", af.lengths)

    for read in af:

        print("cigar string = ", read.cigarstring)
        # see https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.get_cigar_stats
        print("cigar stats = ", read.get_cigar_stats())

        print("in unmmapped = ", read.is_unmapped)
        print("is forward = ", read.is_forward)
        print("is secondary = ", read.is_secondary)
        print("is supplementary = ", read.is_supplementary)

        print("reference name = ", read.reference_name)
        print("reference start = ", read.reference_start)
        print("reference end = ", read.reference_end)
        print("reference length = ", read.reference_length)

        print("query name = ", read.query_name)
        print("query alignment start = ", read.query_alignment_start)
        print("query alignment end = ", read.query_alignment_end)
        print("query alignment length = ", read.query_alignment_length)
        print("query lenght = ", read.query_length)
        break
# %%

# 1) get a list of all query names that have a secondary mapping associated
# 2) get all mappings for all of these reads. For every read save in a dataframe:
#   - query id
#   - fwd/rev mapping
#   - primary/secondary/supplementary mapping
#   - query and reference start/end/length/alignment length
#   - total number of matches/insertions/deletions
# 3) save the results


# What we could do: look at interesting entries in this dataframe,
# check individual reads that have large secondary mappings,
# see how much sequence across the two mappings is aligned, and whether the region is disordered/duplicated
# (cross-check with IGV)

# Start by looking at histograms of density of secondary/primary mappings on the reference.
