import pysam
import numpy as np
import pandas as pd
#from bin.helperClasses.geneIntervalTree.gene_interval_tree import GeneIntervalTree
from scripts.src.helperClasses.refFlat_repr import RefFlatParsed


N_CORES = 10
step = 500000

correct_genenames = {}
tested_genenames = {}


#bam file
infile_bam = pysam.AlignmentFile(snakemake.input["inbam"], "rb")

#parse refflat
refFlat = RefFlatParsed(snakemake.input["refflat"], infile_bam)

#outfile bam
outfile = pysam.AlignmentFile(snakemake.output["outbam"], "wb", template=infile_bam)

#for testing
correct_bam = pysam.AlignmentFile(snakemake.input["correctbam"], "rb")

for read in correct_bam:
    correct_genenames[read.query_name] = set() if not read.has_tag("gn") else set(read.get_tag("gn").split(","))



reads_list = [read for read in infile_bam]
blocks_list = [reads_list[i].get_blocks()[j] for i in range(len(reads_list)) for j in range(len(reads_list[i].get_blocks()))]
starts_list = [blocks_list[i][0] for i in range(len(blocks_list))]
ends_list = [blocks_list[i][1] for i in range(len(blocks_list))]

total = 0
for read in reads_list:
    total += len(read.get_blocks())

R = np.zeros(total, dtype='int64')
B = np.arange(total)
refs = []

ctr = 0
for i in range(len(reads_list)):
    for j in range(len(reads_list[i].get_blocks())):
        R[ctr] = i
        refs.append(infile_bam.getrname(reads_list[i].tid))
        ctr += 1

#group by ref for querying gene tree
refs = pd.DataFrame(data={"R": R, "B": B, "ref": refs, "start": starts_list, "end": ends_list}).groupby("ref")

for ref, group in refs:

    refFlat_intervals =  refFlat.as_intervals(ref)
    ncl = NCLS(refFlat_intervals.starts, refFlat_intervals.ends, refFlat_intervals.index.to_list())
    query_index, ncl_index = ncl.all_overlaps_both(group.start, group.end, group.index.to_list())

    overlaps = pd.DataFrame({"I1" : query_index, "I2" : ncl_index})
    merged = overlaps \
        .merge(group, left_on="I1",right_index=True) \
        .merge(refFlat_intervals, left_on="I2",right_index=True) \
        [["R","B","G"]]

    # how many distinct B's does an R have?
    merged["RB"] = merged[["R", "B"]].groupby("R").B.transform("nunique")
    # how many distinct B's does a G belong to in each R?
    merged["GB"] = merged.groupby(["R", "G"]).B.transform("nunique")

    res = merged[merged.RB==merged.GB]

    for read, grouped_by_read in res.groupby("R"):
        genes_for_read = grouped_by_read.G.to_list().sort()
        as_string = ','.join(genes_for_read)

        tested_genenames[reads_list[read].query_name] = as_string

    print("finished ref \n", ref)





exit()


