import pysam
from bin.helperClasses.geneIntervalTree.gene_interval_tree import GeneIntervalTree
from bin.funcs import *

infile_bam = pysam.AlignmentFile(snakemake.input["inbam"], "rb")
gi_tree = GeneIntervalTree(snakemake.input["refflat"], infile_bam)
outfile = pysam.AlignmentFile(snakemake.output["outbam"], "wb", template=infile_bam)
correct_bam = pysam.AlignmentFile(snakemake.input["correctbam"], "rb")

reads_to_test = {}
correct_reads = {}

for read in correct_bam:
    correct_reads[read.query_name] = {}
    for tag_name, bam_tag in Tags.tags_dict.items():
        if read.has_tag(bam_tag):
            correct_reads[read.query_name][tag_name] = read.get_tag(bam_tag)
        else:
            correct_reads[read.query_name][tag_name] = ""
    correct_reads[read.query_name]["blocks"] = [b for b in read.get_blocks()]

ctr = 0
print("start tagging reads...")
#write to output file
for read in infile_bam:
    tag_read_with_functional_data(read, gi_tree)
    reads_to_test[read.query_name] = {}
    for tag_name, bam_tag in Tags.tags_dict.items():
        if read.has_tag(bam_tag):
            reads_to_test[read.query_name][tag_name] = read.get_tag(bam_tag)
        else:
            reads_to_test[read.query_name][tag_name] = ""
    # blocks (debugging)
    reads_to_test[read.query_name]["blocks"] = [b for b in read.get_blocks()]
    ctr+=1
    if ctr % 100000 == 0:
        print("tagged %ctr reads",ctr)

#test against correct bam file
ctr = 0

if not len(reads_to_test.keys()) == len(correct_reads.keys()):
    raise Exception("read ids not equal")

for read_id in reads_to_test.keys():
    ctr+=1
    if read_id in correct_reads.keys():
        correct_read = correct_reads[read_id]
        read_to_test = reads_to_test[read_id]

        if not correct_read["GENE_NAME_TAG"] == read_to_test["GENE_NAME_TAG"]:
            print("gene name mismatch for read_id", read_id)
            print("correct read blocks:", correct_reads[read_id]['blocks'] )
            print("tested read blocks:", reads_to_test[read_id]['blocks'] )
            print("correct read mapped to genes:", correct_read["GENE_NAME_TAG"])
            print("tested read mapped to genes:",  read_to_test["GENE_NAME_TAG"])
            exit()


