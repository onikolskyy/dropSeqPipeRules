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
        correct_genes_with_coding_overlapped = getGenesWithOverlappedCoding(read.get_blocks(),
                                                                                    correct_genes)
        correct_genes_with_exon_overlapped = getGenesWithOverlappedExon(read.get_blocks(),
                                                                                    correct_genes)
        correct_genes_with_transcript_overlapped = getGenesWithOverlappedTranscript(read.get_blocks(), correct_genes)

        correct_genes_with_utr_overlapped = getGenesWitOverlappedUtr(read.get_blocks(), correct_genes)


        #if correct_genes == correct_genes_with_coding_overlapped.union(correct_genes_with_exon_overlapped):
        if correct_genes == set().union(*[correct_genes_with_coding_overlapped, correct_genes_with_exon_overlapped,
                                          correct_genes_with_transcript_overlapped, correct_genes_with_utr_overlapped]):
            ctr_correct+=1
        else:
            # print("_____________________________")
            # print("correct genes", [gene.name for gene in correct_genes])
            # print("correct_genes_with_transcript_overlapped", [gene.name for gene in correct_genes_with_transcript_overlapped])
            # print("correct_genes_with_coding_overlapped", [ gene.name for gene in correct_genes_with_coding_overlapped])
            # print("correct_genes_with_exon_overlapped", [ gene.name for gene in correct_genes_with_exon_overlapped])
            # dif = correct_genes - correct_genes_with_coding_overlapped.union(correct_genes_with_exon_overlapped)
            # for gene in dif:
            #     gene.verbose()
            #print("_____________________________")
            ctr_wrong += 1

    if CTR_TEST == 1000000:
        break()
print(ctr_correct)
print(ctr_wrong)
exit()

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
ctr_wrong = 1
for read_id in reads_to_test.keys():
    ctr+=1

    if (ctr_wrong >= 10):
        exit()
    if read_id in correct_reads.keys():
        correct_read = correct_reads[read_id]
        read_to_test = reads_to_test[read_id]

        if not correct_read["GENE_NAME_TAG"] == read_to_test["GENE_NAME_TAG"]:
            ctr_wrong+=1
            print("gene name mismatch for read_id", read_id)
            print("correct read blocks:", correct_reads[read_id]['blocks'] )
            print("tested read blocks:", reads_to_test[read_id]['blocks'] )
            print("correct read mapped to genes:", correct_read["GENE_NAME_TAG"])
            correct_genes = [gi_tree.genes[g_id] for g_id in correct_read["GENE_NAME_TAG"].split(",")]
            correct_reads_with_exon_overlapped = getGenesWithOverlappedExon(reads_to_test[read_id]['blocks'], correct_genes)
            correct_genes_with_coding_overlapped = getGenesWithOverlappedCoding(reads_to_test[read_id]['blocks'], correct_genes)
            print("correct_reads_with_exon_overlapped",[g.name for g  in correct_reads_with_exon_overlapped])
            print("correct_reads_with_coding_overlapped",[g.name for g  in correct_genes_with_coding_overlapped])
            print("_________________________")

            # for gene_id in correct_read["GENE_NAME_TAG"].split(","):
            #     gene = gi_tree.genes[gene_id]
            #     print(gene_id, "--> (", gene.start, ",", gene.end, ")", gene.chrom)
            #     for t_name, t in gene.transcripts.items():
            #         print("\t", t_name)
            #         print("\t ", t.transcription_start, t.transcription_end)
            #         print("\t ", t.coding_start, t.coding_end)
            #         print("\t ", t.exons)
            print('\n')
            print("tested read mapped to genes:",  read_to_test["GENE_NAME_TAG"])
            # for gene_id in read_to_test["GENE_NAME_TAG"].split(","):
            #     gene = gi_tree.genes[gene_id]
            #     print(gene_id, "--> (", gene.start, ",", gene.end, ")", gene.chrom)
            #     for t_name, t in gene.transcripts.items():
            #         print("\t", t_name)
            #         print("\t ", t.transcription_start, t.transcription_end)
            #         print("\t ", t.coding_start, t.coding_end)
            #         print("\t ", t.exons)
            print('\n')
            print("________________________________________________________-")



