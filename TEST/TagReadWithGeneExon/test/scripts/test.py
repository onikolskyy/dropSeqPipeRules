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

construction_total = 0
query_total = 0
genes_filtering_total = 0
genes_tagging_total = 0

for read in correct_bam:
    correct_genenames[read.query_name] = set() if not read.has_tag("gn") else set(read.get_tag("gn").split(","))

reads_list = [read for read in infile_bam]
print(len(reads_list))
blocks_list = [reads_list[i].get_blocks()[j] for i in range(len(reads_list)) for j in range(len(reads_list[i].get_blocks()))]
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

RB = pd.DataFrame(data={"R": R, "B": B, "ref": refs})

construction_end = time.time()
query_start = time.time()
query_starts = [blocks_list[i][0] for i in range(len(blocks_list))]
query_ends = [blocks_list[i][1] for i in range(len(blocks_list))]
query = pd.DataFrame({'starts': query_starts, 'ends': query_ends, 'ids': B})
res = gi_tree.get_overlaps(query)

merged = pd.merge(RB, res, on=["B","ref"])[["R","B","G"]]

query_end = time.time()
genes_filtering_start = time.time()

# how many distinct B's does an R have?
merged["RB"] = merged[["R", "B"]].groupby("R").B.transform("nunique")
# how many distinct B's does a G belong to in each R?
merged["GB"] = merged.groupby(["R", "G"]).B.transform("nunique")

tags = merged\
    .merge(merged, right_on=["R", "RB"], left_on=["R", "GB"], how="inner")[["R", "G_x"]]\
    .rename(columns={"G_x": "G"})\
    .groupby("R").agg({"G": lambda x: set(x)})

genes_filtering_end = time.time()
genes_tagging_start = time.time()

for index, row in tags.iterrows():
    read = reads_list[index]
    tested_genenames[read.query_name] = set([gi_tree.get_gene_by_index(ref, index).name for index in row["G"]])

genes_tagging_end = time.time()

print(" \"tagged\"", count_reads, "for", ref)
print("time elapsed: construction->", construction_end - construction_start, "; query->", query_end - query_start, "; gene filtering->", genes_filtering_end - genes_filtering_start, "tagging", genes_tagging_end - genes_filtering_start)

construction_total += construction_end - construction_start
query_total += query_end - query_start
genes_filtering_total += genes_filtering_end - genes_filtering_start
genes_tagging_total += genes_tagging_end - genes_filtering_start

ctr_wrong = 0
ctr_correct = 0
ctr = 0

for query_name, tagged_genes in correct_genenames.items():
    ctr+= 1
    set_to_test = set() if query_name not in tested_genenames else tested_genenames[query_name]
    if set_to_test == tagged_genes:
        ctr_correct+=1
    else:
        ctr_wrong+=1
        print("corrrect:", ctr_correct, "; wrong:", ctr_wrong)
    if (ctr % 100000 == 0): print("tested %i reads" % ctr)

print("corrrect:", ctr_correct, "; wrong:", ctr_wrong)
print("time elapsed TOTAL: construction->", construction_total, "; query->", query_total, "; gene filtering->", genes_filtering_total, "tagging", genes_tagging_total)

exit()


