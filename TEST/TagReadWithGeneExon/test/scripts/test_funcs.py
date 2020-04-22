import pysam
import numpy as np
import pandas as pd
from ncls import NCLS
from time import time
#from bin.helperClasses.geneIntervalTree.gene_interval_tree import GeneIntervalTree
from scripts.src.refFlat_repr import RefFlatParsed

def parse_correct(correct_bam):
    correct_genenames = {}
    correct_lfs = {}
    for read in correct_bam:
        correct_genenames[read.query_name] = read.get_tag("gn") if read.has_tag("gn") else ""
        correct_lfs[read.query_name] = read.get_tag("gf") if read.has_tag("gf") else ""

def make_bam_df(infile_bam):
    reads_list = [read for read in infile_bam]
    blocks_list = [reads_list[i].get_blocks()[j] for i in range(len(reads_list)) for j in
                   range(len(reads_list[i].get_blocks()))]

    # generate query data as list of one base long halfopen intervals (in 0-based coords)

    start = [i for b in range(len(blocks_list)) for i in range(blocks_list[b][0], blocks_list[b][1], 1)]
    end = [i + 1 for b in range(len(blocks_list)) for i in range(blocks_list[b][0], blocks_list[b][1], 1)]
    block = [b for b in range(len(blocks_list)) for i in range(blocks_list[b][0], blocks_list[b][1], 1)]
    read = [r for r in range(len(reads_list)) for b in reads_list[r].get_blocks() for i in range(b[0], b[1], 1)]
    ref = [infile_bam.getrname(reads_list[r].tid) \
           for r in range(len(reads_list)) for b in reads_list[r].get_blocks() for i in range(b[0], b[1], 1)]

    # group by ref for querying gene tree
    refs = pd.DataFrame(data={"read": read, "ref": ref, "start": start, "end": end, "block": block})

    return refs

def make_refflat(file_name,infile_bam):
    return RefFlatParsed(file_name, infile_bam)

def merge_overlaps_df(q,refFlat,ref):
    refFlat_intervals = refFlat.as_intervals(ref)
    if not isinstance(refFlat_intervals, pd.DataFrame):
        return

    q_filtered = q[q["ref"]==ref]

    ncl = NCLS(refFlat_intervals.start.to_numpy(), refFlat_intervals.end.to_numpy(), refFlat_intervals.index.to_numpy())
    query_index, ncl_index = ncl.all_overlaps_both(q_filtered.start.to_numpy(), q_filtered.end.to_numpy(), q_filtered.index.to_numpy())

    overlaps = pd.DataFrame({"I1": query_index, "I2": ncl_index})
    merged = (overlaps \
              .merge(q_filtered, left_on="I1", right_index=True) \
              .merge(refFlat_intervals[["gene", "LF"]], left_on="I2", right_index=True)) \
        [["read", "block", "start", "gene", "LF"]]

    # how many distinct B's does an R have?
    merged["RB"] = merged[["read", "block"]].groupby("read").block.transform("nunique")
    # how many distinct B's does a G belong to in each R?
    merged["GB"] = merged.groupby(["read", "gene"]).block.transform("nunique")

    return merged
