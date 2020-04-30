import itertools
import pysam
import math
import numpy as np
import pandas as pd
from ncls import NCLS
from time import time
from multiprocessing import Pool, RawArray

#from bin.helperClasses.geneIntervalTree.gene_interval_tree import GeneIntervalTree
from src.refFlat_repr import RefFlatParsed

GLOB_REFFLAT = {}


def initializer(glob_intervals, glob_intervals_shape):
    GLOB_REFFLAT["intervals"] = glob_intervals
    GLOB_REFFLAT["intervals_shape"] = glob_intervals_shape


def worker(blocks_list):
    LFs = pd.DataFrame({"name": ["INTERGENIC", "INTRONIC", "UTR", "CODING"]})

    # get intervals from shared memory
    sharr_np = np.frombuffer(GLOB_REFFLAT["intervals"],dtype=np.int64).reshape(GLOB_REFFLAT["intervals_shape"])
    reference = pd.DataFrame(sharr_np,columns=["gene","start","end","LF"])
    # init NCLS
    ncl = NCLS(
        np.ascontiguousarray(reference.start.to_numpy()),
        np.ascontiguousarray(reference.end.to_numpy()),
        np.ascontiguousarray(reference.index.to_numpy())
    )
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
              .merge(reference[["gene", "LF"]], left_on="I2", right_index=True)) \
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
    res = pd.concat([multi_block, single_block])

    return res

####################################################################################

#bam file
infile_bam = pysam.AlignmentFile(snakemake.input["inbam"], "rb")
#parse refflat
refFlat = RefFlatParsed(snakemake.input["refflat"], infile_bam)
outfile = pysam.AlignmentFile(snakemake.output[0], "wb", template=infile_bam)

#################################################################################
# make reads data frame
#################################################################################


def chunked_iterable(iterable, size):
    it = iter(iterable)
    while True:
        chunk = tuple(itertools.islice(it, size))
        if not chunk:
            break
        yield chunk


bam_header = infile_bam.header
reads_chunked = chunked_iterable(infile_bam,5*10**6)

chunk_ctr = 1
tot_reads = 0
for chunk in reads_chunked:
    reads_from_chunk = [r for r in chunk]
    print("chunk", chunk_ctr, len(reads_from_chunk))
    blocks_packed = [reads_from_chunk[i].get_blocks() for i in range(len(reads_from_chunk))]
    refs = pd.DataFrame({"ref": [infile_bam.getrname(reads_from_chunk[r].tid) for r in range(len(reads_from_chunk))]})
    for ref, group in refs.groupby("ref"):

        # get reference data
        reference_data = refFlat.as_intervals(ref)
        if reference_data is not None:
            gene, refFlat_intervals = reference_data[0], reference_data[1]
            # store packed blocks
            chunk_for_ref = [blocks_packed[i] for i in group.index.to_list()]
            iterator_over_chunk = chunked_iterable(chunk_for_ref, 500000)
            intervals = refFlat_intervals[["gene","start","end","LF"]].to_numpy()

            # store it in shared memory
            sharr = RawArray('l', intervals.shape[0]*intervals.shape[1])
            sharr_np = np.frombuffer(sharr, dtype=np.int64).reshape(intervals.shape)
            np.copyto(sharr_np, intervals)

            with Pool(processes=snakemake.threads, initializer=initializer,
                      initargs=(sharr, intervals.shape)) as pool:
                results = pool.map(worker, iterator_over_chunk)
    print("finished chunk",chunk_ctr)

    #TODO: write data to output bam!

    chunk_ctr+= 1
    tot_reads+= len(reads_from_chunk)
print("TOTAL:", len(reads_from_chunk))


