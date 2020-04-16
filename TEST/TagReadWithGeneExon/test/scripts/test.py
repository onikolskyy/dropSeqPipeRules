import pysam
import numpy as np
import pandas as pd
from ncls import NCLS
from time import time
#from bin.helperClasses.geneIntervalTree.gene_interval_tree import GeneIntervalTree
from src.refFlat_repr import RefFlatParsed


N_CORES = 10
step = 500000

correct_genenames = {}
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


reads_list = [read for read in infile_bam]
blocks_list = [reads_list[i].get_blocks()[j] for i in range(len(reads_list)) for j in range(len(reads_list[i].get_blocks()))]
starts_list = [blocks_list[i][0] for i in range(len(blocks_list))]
ends_list = [blocks_list[i][1] for i in range(len(blocks_list))]

total = 0
for read in reads_list:
    total += len(read.get_blocks())

R = np.zeros(total, dtype='int64')
#B = np.arange(total)
refs = []

ctr = 0
for i in range(len(reads_list)):
    for j in range(len(reads_list[i].get_blocks())):
        R[ctr] = i
        refs.append(infile_bam.getrname(reads_list[i].tid))
        ctr += 1

#group by ref for querying gene tree
refs = pd.DataFrame(data={"R": R, "ref": refs, "start": starts_list, "end": ends_list}).groupby("ref")

for ref, group in refs:

    print("starting ref",ref)

    t_start = time()

    refFlat_intervals =  refFlat.as_intervals(ref)

    #if ref is not in refFlat, continue to next ref
    if not isinstance(refFlat_intervals,pd.DataFrame):
        continue

    ncl = NCLS(refFlat_intervals.start.to_numpy(), refFlat_intervals.end.to_numpy(), refFlat_intervals.index.to_numpy())
    query_index, ncl_index = ncl.all_overlaps_both(group.start.to_numpy(), group.end.to_numpy(), group.index.to_numpy())

    overlaps = pd.DataFrame({"I1" : query_index, "I2" : ncl_index})
    merged = overlaps \
        .merge(group, left_on="I1",right_index=True) \
        .merge(refFlat_intervals, left_on="I2",right_index=True) \
        [["R","B","G"]]

    # how many distinct B's does an R have?
    merged["RB"] = merged[["R", "B"]].groupby("R").B.transform("nunique")
    # how many distinct B's does a G belong to in each R?
    merged["GB"] = merged.groupby(["R", "G"]).B.transform("nunique")

    # split into reads with singl block and reads with multiple blocks
    single_block = merged[merged.RB == 1][["R","B","G"]]
    multiple_blocks = merged[merged.GB != 1]

    # filter out only those genes which are overlapped by all blocks of a read
    multiple_blocks_filterd = multiple_blocks[multiple_blocks.RB == multiple_blocks.GB][["R","B","G"]]
    # retain unique gene set for each read
    unique = multiple_blocks_filterd[["R","G"]].drop_duplicates()
    multiple_blocks_unique = multiple_blocks_filterd[multiple_blocks_filterd.index.isin(unique.index.to_list())]

    # concat splitted data
    res = pd.concat([multiple_blocks_unique,single_block])

    t_end = time()

    logfile.write("ref %s took %s"%(ref,str(t_end-t_start)))
    print("ref %s took %s"%(ref,str(t_end-t_start)))

    for read, grouped_by_read in res.groupby("R"):
        genes_for_read = grouped_by_read.G.to_list()
        genes_for_read.sort()
        as_string = ','.join(genes_for_read)

        tested_genenames[reads_list[read].query_name] = as_string

    print("finished ref \n", ref)

ctr_correct = 0
ctr_wrong = 0
ctr_not_found = 0
for qname, genenames in tested_genenames.items():
    if qname not in correct_genenames:
        ctr_not_found+=1
    else:
        if genenames == correct_genenames[qname]:
            ctr_correct += 1
        else:
            ctr_wrong += 1
            logfile.write("correct:" ,correct_genenames[qname],"wrong:", genenames, "\n" )
           # print("correct:", genenames, "; tested:", tested_genenames[qname])
print("correct: %i; wrong: %i"%(ctr_correct,ctr_wrong))
logfile.write("correct: %i; wrong: %i"%(ctr_correct,ctr_wrong))
logfile.close()
exit()


