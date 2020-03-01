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
    correct_genes = set([gi_tree.genes[g_id] for g_id in read.get_tag("gn").split(",")])
    if "" in correct_genes:
        ctr_correct+=1
        continue
    else:
        filter_transcripts  = defGenesOverlappedByFn(read.get_blocks(), gi_tree, "transcript")
        filter_exon  = defGenesOverlappedByFn(read.get_blocks(), gi_tree, "exon")
        filter_utr  = defGenesOverlappedByFn(read.get_blocks(), gi_tree, "utr")
        filter_coding = defGenesOverlappedByFn(read.get_blocks(), gi_tree, "coding")
        if correct_genes == set().union(*[filter_transcripts,filter_exons, filter_utr, filter_coding]):
            ctr_correct+=1
        else:
            print("correct genes", [gene.name for gene in correct_genes])
            print("filter_transcripts", [gene.name for gene in filter_transcripts])
            print("filter_exon", [gene.name for gene in filter_exons])
            print("filter_utr", [gene.name for gene in filter_utr])
            print("filter_coding", [gene.name for gene in filter_coding])
    if CTR_TEST == 10:
        break
print(ctr_correct)
print(ctr_wrong)
exit()

