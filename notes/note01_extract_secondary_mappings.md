# Extract supplementary mappings

After getting familiar with sam files, `pysam` and `pandas`, we can work on extracting information on supplementary mappings from the sam file.

The aim is to create a dataframe that contains information on primary and supplementary mappings of interesting reads, that we can later use to locate which regions of the genome are connected by these reads.

## parsing the sam file

One would need to parse the sam file twice. The first time one only considers supplementary mappings, and extract all of the query names of these mappings in a list.

Then one can parse this list again and consider all reads whose query name is in the previously-created list (i.e. all queries that have a supplementary mapping). For each of these mappings one would extract:
- the query name
- whether the read is primary, supplementary or secondary
- whether the mapping is forward or reverse
- the start/end of the alignment on the reference/query
- the alignment length (e.g. on the query) and the total query length (complete of unaligned/clipped regions, this can be done with the `infer_read_length` function)
- the total number of matches/insertions/deletions in the alignment

The results should be stored in a pandas dataframe with the above-specified columns, and then saved in `csv` format.

## creating a python script

We might need to perform this operation multiple times on different sam files. Ideally one would create a python script that can be used from the command line. The script should take as input a sam file and save the resulting dataframe in a file with the name specified by the user. The usage would be something like this:
```
python3 extract_suppl_reads_info.py
    --sam reads.sam
    --csv suppl_reads_info.csv
```
For convenient parsing of command line arguments one can use [argparse](https://docs.python.org/3/library/argparse.html).
