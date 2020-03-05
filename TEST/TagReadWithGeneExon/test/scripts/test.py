import pysam
import numpy as np
import pandas as pd
from bin.helperClasses.geneIntervalTree.gene_interval_tree import GeneIntervalTree
from bin.funcs import *



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


for ref in reads_dict:
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

    query_starts = [blocks_list[i][0] for i in range(len(blocks_list))]
    query_ends = [blocks_list[i][1] for i in range(len(blocks_list))]
    query = pd.DataFrame({'starts': query_starts, 'ends': query_ends, 'ids': B})
    overlaps = gi_tree.get_all_overlaps_by_ref(query, ref)

    RBG = pd.merge(RB, overlaps, on="B")

    grouped_by_reads = RBG.groupby("R")
    for r, grouped_by_read in grouped_by_reads:
        gene_ids = set.intersection(*[set(grouped_by_block['G']) for b, grouped_by_block in grouped_by_read.groupby("B")])
        tested_genenames[reads_list[r].tid] = set([gi_tree.get_gene_by_id(ref, gene_id).name for gene_id in gene_ids])


ctr_wrong = 0
ctr_correct = 0

for pid, tagged_genes in correct_genenames:
    if len(tagged_genes) == 0:
        if pid not in tested_genenames or len(tested_genenames[pid]) == 0:
            ctr_correct += 1
            continue
        else:
            ctr_wrong += 1
    else:
        if tested_genenames[pid] == correct_genenames[pid]:
            ctr_correct += 1
        else:
            ctr_wrong += 1

print("corrrect:", ctr_correct, "; wrong:", ctr_wrong)
exit()


