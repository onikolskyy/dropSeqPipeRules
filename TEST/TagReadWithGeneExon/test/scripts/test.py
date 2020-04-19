import pysam
import numpy as np
import pandas as pd
from ncls import NCLS
from time import time
#from bin.helperClasses.geneIntervalTree.gene_interval_tree import GeneIntervalTree
from src.refFlat_repr import RefFlatParsed

LFs = pd.DataFrame({"name" : ["INTERGENIC", "INTRONIC", "UTR", "CODING"]})

def find_LF_for_genes(df, is_multi_block=False):
    # what is maximum LF of a Base?
    df["maxLF"] = df[["read", "block", "start", "gene","LF"]].groupby(["read", "block", "start", "gene"]).transform(max)
    # reatain only max LF
    df = df[df.maxLF == df.LF]
    if not is_multi_block:
        # retain unique list of LF for each gene
        return df[["read","gene","LF"]].drop_duplicates()
    else:
        # filter out only those genes which are overlapped by all blocks of a read
        # retain a list of unique LFs for each gene
        return df[df.RB == df.GB][["read", "gene", "LF"]].drop_duplicates()



N_CORES = 10
step = 500000

correct_genenames = {}

correct_lfs = {}

tested_genenames = {}


#bam file
infile_bam = pysam.AlignmentFile(snakemake.input["inbam"], "rb")
print("start parsing refflat")
#parse refflat
refFlat = RefFlatParsed(snakemake.input["refflat"], infile_bam)
print("end parsing refflat")
#outfile bam
#outfile = pysam.AlignmentFile(snakemake.output["outbam"], "wb", template=infile_bam)
logfile = open(snakemake.output["out"],"a")
#for testing
correct_bam = pysam.AlignmentFile(snakemake.input["correctbam"], "rb")

for read in correct_bam:
    correct_genenames[read.query_name] =  read.get_tag("gn") if read.has_tag("gn") else ""
    correct_lfs[read.query_name] = read.get_tag("gf") if read.has_tag("gf") else ""

reads_list = [read for read in infile_bam]
blocks_list = [reads_list[i].get_blocks()[j] for i in range(len(reads_list)) for j in range(len(reads_list[i].get_blocks()))]

# generate query data as list of one base long halfopen intervals (in 0-based coords)

start = [i for b in range(len(blocks_list)) for i in range(blocks_list[b][0],blocks_list[b][1],1)]
end = [i+1 for b in range(len(blocks_list)) for i in range(blocks_list[b][0],blocks_list[b][1],1)]
block = [b for b in range(len(blocks_list)) for i in range(blocks_list[b][0],blocks_list[b][1],1)]
read = [r for r in range(len(reads_list)) for b in reads_list[r].get_blocks() for i in range(b[0],b[1],1)]
ref = [infile_bam.getrname(reads_list[r].tid)\
       for r in range(len(reads_list)) for b in reads_list[r].get_blocks() for i in range(b[0],b[1],1)]

# group by ref for querying gene tree
refs = pd.DataFrame(data={"read": read, "ref": ref, "start": start, "end": end, "block": block})
grouped = refs.groupby("ref")

wrg_single =0
wrg_multiple = 0

for ref, group in grouped:

    print("starting ref",ref)

    refFlat_intervals =  refFlat.as_intervals(ref)

    t_start = time()

    #if ref is not in refFlat, continue to next ref
    if not isinstance(refFlat_intervals,pd.DataFrame):
        continue

    ncl = NCLS(refFlat_intervals.start.to_numpy(), refFlat_intervals.end.to_numpy(), refFlat_intervals.index.to_numpy())
    query_index, ncl_index = ncl.all_overlaps_both(group.start.to_numpy(), group.end.to_numpy(), group.index.to_numpy())

    overlaps = pd.DataFrame({"I1" : query_index, "I2" : ncl_index})
    merged = ( overlaps \
        .merge(group, left_on="I1",right_index=True) \
        .merge(refFlat_intervals[["gene","LF"]], left_on="I2",right_index=True) )\
        [["read","block","start","gene","LF"]]

    print("merge complete \n")

    # how many distinct B's does an R have?
    merged["RB"] = merged[["read", "block"]].groupby("read").block.transform("nunique")
    # how many distinct B's does a G belong to in each R?
    merged["GB"] = merged.groupby(["read", "gene"]).block.transform("nunique")

    # split into reads with singl block and reads with multiple blocks, handle separately
    single_block = find_LF_for_genes( merged[merged.RB == 1] )
    multiple_blocks = find_LF_for_genes( merged[merged.GB != 1], True)

    single_block = single_block.merge(LFs, left_on="LF", right_index=True)
    multiple_blocks = multiple_blocks.merge(LFs, left_on="LF", right_index=True)

    single_block = single_block.sort_values(["read","gene"])
    multiple_blocks = multiple_blocks.sort_values(["read","gene"])

    print(single_block)
    print(multiple_blocks)

    for read, grouped_by_read in single_block.groupby("read"):
        genes_for_read = grouped_by_read.gene.to_list()
        genes_as_string = ','.join(genes_for_read)
        lf_as_string = ",".join(grouped_by_read.name.to_list())
        if not correct_genenames[reads_list[read].query_name]  == genes_as_string:
            print("SINGLE: correct:%s-->%s, wrong:%s-->%s \n"\
                  %(correct_genenames[reads_list[read].query_name],
                    correct_lfs[reads_list[read].query_name],
                    genes_as_string,
                    lf_as_string))
            wrg_single+=1

    for read, grouped_by_read in multiple_blocks.groupby("read"):
        genes_for_read = grouped_by_read.gene.to_list()
        genes_as_string = ','.join(genes_for_read)
        lf_as_string = ",".join(grouped_by_read.name.to_list())
        if not correct_genenames[reads_list[read].query_name] == genes_as_string:
            print("MULTIPLE: correct:%s-->%s, wrong:%s-->%s \n" \
                  % (correct_genenames[reads_list[read].query_name],
                     correct_lfs[reads_list[read].query_name],
                     genes_as_string,
                     lf_as_string))
            wrg_multiple += 1

    # concat splitted data
    res = pd.concat([multiple_blocks,single_block]).merge(LFs, right_index=True, left_on="LF")[["read","gene"]]
    t_end = time()

    logfile.write("ref %s took %s \n"%(ref,str(t_end-t_start)))
    print("ref %s took %s \n"%(ref,str(t_end-t_start)))

    # for read, grouped_by_read in res.groupby("read"):
    #     genes_for_read = grouped_by_read.gene.to_list()
    #     genes_for_read.sort()
    #     as_string = ','.join(genes_for_read)

    #    tested_genenames[reads_list[read].query_name] = as_string

    print("finished ref \n", ref)

ctr_correct = 0
ctr_wrong = 0
ctr_not_found = 0
# for qname, genenames in tested_genenames.items():
#     if qname not in correct_genenames:
#         ctr_not_found+=1
#     else:
#         if genenames == correct_genenames[qname]:
#             ctr_correct += 1
#         else:
#             ctr_wrong += 1
#             logfile.write("correct:%s; wrong: %s \n" %(correct_genenames[qname], genenames))
#             print("correct:", correct_genenames[qname], "; tested:", tested_genenames[qname])
print("------->single_wrg:%i, mult_wrg:%i"%(wrg_single,wrg_multiple))

#print("correct: %i; wrong: %i"%(ctr_correct,ctr_wrong))
#logfile.write("correct: %i; wrong: %i"%(ctr_correct,ctr_wrong))
logfile.close()
exit()


