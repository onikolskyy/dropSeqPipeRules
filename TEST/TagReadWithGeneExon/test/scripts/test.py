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

reads_dict = collections.defaultdict(lambda: {"reads_list": []})

correct_genenames = {}
tested_genenames = {}

for read in correct_bam:
    correct_genenames[read.tid] = set() if not read.has_tag("gn") else set(read.get_tag("gn").split(","))

for read in infile_bam:
    ref = infile_bam.getrname(read.tid)
    reads_dict[ref]["reads_list"].append(read)


print("start \"tagging\"")

for ref in reads_dict:
    if ref not in gi_tree.trees:
        print("scipping", ref)
        continue

    construction_start = time.time()

    count_reads = len(reads_dict[ref]["reads_list"])
    reads_list = reads_dict[ref]["reads_list"]
    blocks_list = [reads_list[i].get_blocks()[j] for i in range(len(reads_list)) for j in range(len(reads_list[i].get_blocks()))]
    total = 0
    for read in reads_list:
        total += len(read.get_blocks())

    R = np.zeros(total, dtype='int16')
    B = np.arange(total)

    ctr = 0
    for i in range(len(reads_list)):
        for j in range(len(reads_list[i].get_blocks())):
            R[ctr] = i
            ctr += 1

    RB = pd.DataFrame(data={"R": R, "B": B})

    construction_end = time.time()
    query_start = time.time()
    query_starts = [blocks_list[i][0] for i in range(len(blocks_list))]
    query_ends = [blocks_list[i][1] for i in range(len(blocks_list))]
    query = pd.DataFrame({'starts': query_starts, 'ends': query_ends, 'ids': B})
    BG = gi_tree.get_all_overlaps_by_ref(query, ref)

    RBG = pd.merge(RB, BG, on="B")

    query_end = time.time()
    genes_filtering_start = time.time()

    # how many distinct B's does an R have?
    RBG["RB"] = RBG[["R", "B"]].groupby("R").B.transform("nunique")
    # how many distinct B's does a G belong to in each R?
    RBG["GB"] = RBG.groupby(["R", "G"]).B.transform("nunique")

    print(RBG.sort(by=["R","B"]))

    print(RBG.merge(RBG, right_on=["R", "RB"], left_on=["R", "GB"], how="inner")[["R", "G_x"]])
    exit()

    tags = RBG\
        .merge(RBG, right_on=["R", "RB"], left_on=["R", "GB"], how="inner")[["R", "G_x"]]\
        .rename(columns={"G_x": "G"})\
        .groupby("R").agg({"G": lambda x: set(x)})

    print(tags)

    genes_filtering_end = time.time()
    genes_tagging_start = time.time()

    for index, row in tags.iterrows():
        read = reads_list[index]
        tested_genenames[read.tid] = set([gi_tree.get_gene_by_index(ref, index).name for index in row["G"]])

    genes_tagging_end = time.time()

    print(" \"tagged\"", count_reads, "for", ref)
    print("time elapsed: construction->", construction_end - construction_start, "; query->", query_end - query_start, "; gene filtering->", genes_filtering_end - genes_filtering_start, "tagging", genes_tagging_end - genes_filtering_start)

ctr_wrong = 0
ctr_correct = 0
ctr = 0

print(len(correct_genenames))

for pid, tagged_genes in correct_genenames.items():
    ctr+= 1
    set_to_test = set() if pid not in tested_genenames else tested_genenames[pid]
    if set_to_test == tagged_genes:
        ctr_correct+=1
    else:
        ctr_wrong+=1
        print("corrrect:", ctr_correct, "; wrong:", ctr_wrong)
    if (ctr % 100000 == 0): print("tested %i reads" % ctr)

print("corrrect:", ctr_correct, "; wrong:", ctr_wrong)
exit()


