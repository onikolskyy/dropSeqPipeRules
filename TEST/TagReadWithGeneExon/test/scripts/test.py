import pysam
import numpy as np
import pandas as pd
from bin.helperClasses.geneIntervalTree.gene_interval_tree import GeneIntervalTree
from bin.funcs import *
import time


infile_bam = pysam.AlignmentFile(snakemake.input["inbam"], "rb")
gi_tree = GeneIntervalTree(snakemake.input["refflat"], infile_bam)
outfile = pysam.AlignmentFile(snakemake.output["outbam"], "wb", template=infile_bam)
correct_bam = pysam.AlignmentFile(snakemake.input["correctbam"], "rb")

correct_genenames = {}
tested_genenames = {}


# construct dataframe with reads and corresponding blocks
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

RB = pd.DataFrame(data={"R": R, "B": B, "ref": refs, "starts": starts_list, "ends": ends_list})

STEP_SIZE = 100000
i = 0
tot_reads = len(reads_list)
while i < tot_reads:
    lo = i
    hi = i+STEP_SIZE if i+STEP_SIZE < tot_reads else tot_reads
    chunk = RB[RB["R"].isin(range(tot_reads)[lo:hi])]
    res = gi_tree.get_overlaps(chunk)

    # how many distinct B's does an R have?
    res["RB"] = res[["R", "B"]].groupby("R").B.transform("nunique")
    # how many distinct B's does a G belong to in each R?
    res["GB"] = res.groupby(["R", "G"]).B.transform("nunique")

    tags = res\
        .merge(res, right_on=["R", "RB"], left_on=["R", "GB"], how="inner")[["R", "G_x","LF_x"]]\
        .rename(columns={"G_x": "G", "LF_x": "LF" })\
        .groupby(["R", "B"]).agg({"G": lambda x: set(x)})\
        .reset_index()\
        .groupby("R")\
        .agg({"G": lambda x: set(x), "LF": lambda x: [e for sublist in list(x) for e in sublist]})

    for index, row in tags.iterrows():
        read = reads_list[index]
        tested_genenames[read.query_name] = row["G"]

# test stats

tot_corr = 0
tot_wrong = 0

with open(snakemake.output[0], "w") as output:
    for query_name, genes in correct_genenames.items():
        if tested_genenames[query_name] == genes:
            tot_corr += 1
        else:
            output.write("%s\t%s" % ((",").join([g in genes]),(",").join([g in tested_genenames[query_name]])) )


exit()


