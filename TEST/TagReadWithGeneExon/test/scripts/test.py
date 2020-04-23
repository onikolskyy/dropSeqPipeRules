import pysam
import math
import numpy as np
import pandas as pd
from ncls import NCLS
from time import time
from multiprocessing import Pool, Array

#from bin.helperClasses.geneIntervalTree.gene_interval_tree import GeneIntervalTree
from src.refFlat_repr import RefFlatParsed


GLOB_REFFLAT = {}

def initializer(sharr, shape):
    GLOB_REFFLAT["sharr"] = sharr
    GLOB_REFFLAT["shape"] = shape

def worker(blocks_list):
    LFs = pd.DataFrame({"name": ["INTERGENIC", "INTRONIC", "UTR", "CODING"]})

    sharr_np = np.frombuffer(GLOB_REFFLAT["sharr"]).reshape(GLOB_REFFLAT["shape"])
    refFlat_intervals = pd.DataFrame(sharr_np,collumns=["gene","start","end","LF"])

    # init NCLS
    ncl = NCLS(refFlat_intervals.start.to_numpy(), refFlat_intervals.end.to_numpy(), refFlat_intervals.index.to_numpy())

    # unpack blocks and create pandas
    blocks_list_unpacked = [blocks_list[b][i] for b in range(len(blocks_list)) for i in range(len(blocks_list[b]))]
    ranges_start = [np.arange(blocks_list_unpacked[b][0], blocks_list_unpacked[b][1]) for b in
                    range(len(blocks_list_unpacked))]
    ranges_end = [np.arange(blocks_list_unpacked[b][0] + 1, blocks_list_unpacked[b][1] + 1) for b in
                  range(len(blocks_list_unpacked))]
    ranges_read = [
        np.full(blocks_list[r][b][1] - blocks_list[r][b][0], r)
        for r in range(len(blocks_list))
        for b in range(len(blocks_list[r]))
    ]

    ranges_block = [np.full(blocks_list_unpacked[b][1] - blocks_list_unpacked[b][0], b) for b in
                    range(len(blocks_list_unpacked))]

    start = np.concatenate(ranges_start)
    end = np.concatenate(ranges_end)
    block = np.concatenate(ranges_block)
    read = np.concatenate(ranges_read)

    reads = pd.DataFrame({"read": read, "start": start, "end": end, "block": block})

    # merge annotations and reads

    query_index, ncl_index = ncl.all_overlaps_both(reads.start.to_numpy(), reads.end.to_numpy(), reads.index.to_numpy())

    overlaps = pd.DataFrame({"I1": query_index, "I2": ncl_index})
    merged = (overlaps \
              .merge(reads, left_on="I1", right_index=True) \
              .merge(refFlat_intervals[["gene", "LF"]], left_on="I2", right_index=True)) \
        [["read", "block", "start", "gene", "LF"]]

    # how many distinct B's does an R have?
    merged["RB"] = merged[["read", "block"]].groupby("read").block.transform("nunique")
    # how many distinct B's does a G belong to in each R?
    merged["GB"] = merged.groupby(["read", "gene"]).block.transform("nunique")

    # split into reads with singl block and reads with multiple blocks, handle separately
    single_block = merged[merged.RB == 1]
    multi_block = merged[merged.GB != 1]

    # process single block
    single_block["maxLF"] = single_block[["read", "block", "start", "gene", "LF"]].groupby(
        ["read", "block", "start", "gene"]).transform(max)
    single_block = single_block[single_block["maxLF"] == single_block["LF"]][
        ["read", "gene", "LF"]].drop_duplicates().sort_values(["read", "gene"])
    single_block = single_block.merge(LFs, left_on="LF", right_index=True)[["read", "gene", "name"]]

    # process multi block
    # retain only genes which are overlapped by all blocks
    multi_block = multi_block[multi_block.RB == multi_block.GB]
    multi_block["maxLF"] = multi_block[["read", "block", "start", "gene", "LF"]].groupby(
        ["read", "block", "start", "gene"]).transform(max)
    multi_block = multi_block[multi_block["maxLF"] == multi_block["LF"]][["read", "gene", "LF"]].drop_duplicates()
    multi_block = multi_block.merge(LFs, left_on="LF", right_index=True)[["read", "gene", "name"]]

    # concat splitted data
    res = pd.concat([multi_block, single_block]).merge(LFs, right_index=True, left_on="LF")[["read", "gene"]]

    return res

####################################################################################



#bam file
infile_bam = pysam.AlignmentFile(snakemake.input["inbam"], "rb")
#parse refflat
refFlat = RefFlatParsed(snakemake.input["refflat"], infile_bam)
logfile = open(snakemake.output["out"],"a")

#################################################################################
# make reads data frame
#################################################################################
bam_header = infile_bam.header


# subsample for each ref
refs = [pair["SN"] for pair in bam_header["SQ"]]
reads_for_ref = {}
for ref in refs:
    reads_for_ref["ref"] = [read for read in infile_bam.fetch(ref)]

packed_blocks_for_ref = {}
for ref in refs:
    packed_blocks_for_ref[ref] = [reads_for_ref[ref][r].get_blocks() for r in range(len(reads_for_ref[ref]))]

for ref, packed_blocks in packed_blocks_for_ref.items():
    refFlat_intervals =  refFlat.as_intervals(ref).to_numpy()
    sharr = Array('c_int64', refFlat_intervals.shape[0]*refFlat_intervals.shape(1), False)
    sharr_np = np.frombuffer(sharr).reshape(refFlat_intervals.shape)
    np.copyto(sharr_np,refFlat_intervals)

    with Pool(processes=snakemake.threads, initializer=initializer, initargs=(sharr,refFlat_intervals.shape)) as pool:
        results = pool.map(worker, packed_blocks_for_ref[ref], chunksize=10**6)


