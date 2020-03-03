import pysam
from bin.helperClasses.geneIntervalTree.gene_interval_tree import GeneIntervalTree
from bin.funcs import *

infile_bam = pysam.AlignmentFile(snakemake.input["inbam"], "rb")
gi_tree = GeneIntervalTree(snakemake.input["refflat"], infile_bam)
outfile = pysam.AlignmentFile(snakemake.output["outbam"], "wb", template=infile_bam)
correct_bam = pysam.AlignmentFile(snakemake.input["correctbam"], "rb")

reads_to_test = {}
correct_reads = {}

ctr_correct = 0
ctr_wrong = 0
CTR_TEST = 0
for read in correct_bam:
    CTR_TEST+=1
    correct_reads[read.query_name] = {}
    for tag_name, bam_tag in Tags.tags_dict.items():
        if read.has_tag(bam_tag):
            correct_reads[read.query_name][tag_name] = read.get_tag(bam_tag)
        else:
            correct_reads[read.query_name][tag_name] = ""
    correct_reads[read.query_name]["blocks"] = [b for b in read.get_blocks()]
    # test if all genes a read is mapped to are only overlapped in coding section
    if not read.has_tag("gn"):
        ctr_correct+=1
        continue
    correct_genes = set(read.get_tag("gn").split(","))
    if "" in correct_genes:
        ctr_correct+=1
        continue
    else:
        ref = infile_bam.getrname(read.tid)
        filter_chrom = set.intersection(*[set(gi_tree.get_overlaps(block, ref)) for block in read.get_blocks()])
        if correct_genes == set([filtered.name for filtered in filter_chrom]):
            ctr_correct+=1
        else:
            print("correct genes",  correct_genes)
            print("filter_chrom",  filter_chrom)
print(ctr_correct)
print(ctr_wrong)
exit()

