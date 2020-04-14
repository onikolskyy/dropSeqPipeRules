import pysam
import numpy as np
import pandas as pd
from bin.helperClasses.geneIntervalTree.gene_interval_tree import GeneIntervalTree
from bin.funcs import *
import time
import ray

N_CORES = 10
step = 500000
ray.init()

@ray.remote
class Worker:
    def __init__(self, t):
        self.t = t

    def work(self,query,ref):
        response = ray.get(self.t.get_overlaps.remote(query,ref))
        merged = query\
            .merge(response["overlaps"], on="B") \
            .merge(response["intervals"], right_index=True,left_on="index")

        # how many distinct B's does an R have?
        merged["RB"] = merged[["R", "B"]].groupby("R").B.transform("nunique")
        # how many distinct B's does a G belong to in each R?
        merged["GB"] = merged.groupby(["R", "G"]).B.transform("nunique")

        res  = merged[merged.RB==merged.GB]
        
        return res

infile_bam = pysam.AlignmentFile(snakemake.input["inbam"], "rb")
refFlat_repr = GeneIntervalTree(snakemake.input["refflat"], infile_bam)
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

refs = pd.DataFrame(data={"R": R, "B": B, "ref": refs, "starts": starts_list, "ends": ends_list}).groupby("ref")

workers = [Worker.remote(refFlat_repr) for _ in range(N_CORES)]
result_ids = []

for ref, reads_for_ref in refs:
    nrows = reads_for_ref.shape[0]
    while i < nrows:
        lo = i
        hi = i + step if i + step < nrows else nrows
        result_ids.append(workers[ctr % N_CORES].do_query.remote(reads_for_ref[lo:hi]))

results= ray.get(result_ids)

exit()


