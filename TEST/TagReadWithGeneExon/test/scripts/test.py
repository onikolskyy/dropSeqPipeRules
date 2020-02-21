import pysam
from .bin.helperClasses.geneIntervalTree.gene_interval_tree import GeneIntervalTree
from .bin.funcs import *

infile_bam = pysam.AlignmentFile(snakemake.input["inbam"], "rb")
gi_tree = GeneIntervalTree(snakemake.input["refflat"])
outfile = pysam.AlignmentFile(snakemake.output["outbam"], "wb", template=infile_bam)
correct_bam = pysam.AlignmentFile(snakemake.input["correctbam"], "rb")

reads_to_test = {}
correct_reads = {}

for read in correct_bam:
    correct_reads[read.query_name] = {}
    for tag in Tags:
        correct_reads[read.query_name][tag.name] = read.get_tag(tag.name)


#write to output file
for read in infile_bam:
    tag_read_with_functional_data(read, gi_tree)
    reads_to_test[read.query_name] = {}
    for tag in Tags:
        reads_to_test[read.query_name][tag.name] = read.get_tag(tag.name)
    outfile.write(read)

#test against correct bam file
ctr = 0

if not len(reads_to_test.keys()) == len(correct_reads.keys()):
    raise Exception("read ids not equal")

for read_id in reads_to_test.keys():
    if read_id in correct_reads.keys():
        correct_read = correct_reads[read_id]
        read_to_test = reads_to_test[read_id]
        for tag in Tags:
            if read_to_test[tag.name] == correct_read[tag.name]:
                pass
            else:
                raise Exception("Unequal reads: on tag %s - EXITING"%tag.name)
    else:
        raise Exception("id %s is present in test data set but is not present in correct data set"%read_id)
    ctr=ctr+1
    if ctr % 10000 == 0:
        print("tested %i reads"%ctr)



