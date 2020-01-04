from Bio import SeqIO
from collections import defaultdict
import pysam
import gzip
import sys

discard_secondary_alignements = snakemake.params['discard_secondary_alignements']
fastq_parser = SeqIO.parse(gzip.open(snakemake.input[0], "rt"), "fastq")
infile_bam = pysam.AlignmentFile(snakemake.input[1], "rb")
outfile = pysam.AlignmentFile(snakemake.output[0], "wb", template=infile_bam)

mapping=defaultdict(dict)
read_barcodes = defaultdict(lambda :{'XC':'','XM':''})

barcodes_struct = {
	'BC_start':snakemake.params['BC_start'],
	'BC_end':snakemake.params['BC_end'],
	'UMI_start':snakemake.params['UMI_start'],
	'UMI_end':snakemake.params['UMI_end']
	}

print("parsing witelist")

#construct bins
with open(snakemake.input[2],'r') as whitelist:
    for line in whitelist:
        if len(line.strip().split()) > 2:  # This means we didn't find any other linked barcode
            (reference, extended_ref, counts_ref, counts_ext) = line.strip().split()
            for barcode in extended_ref.split(','):
                mapping[barcode]["ref"]=reference
            mapping[reference]["ref"]=reference
        else:
            (reference, counts_ref) = line.strip().split()
            mapping[reference]["ref"] = reference

print("parsing whitelist done")
print("saving barcode and umi for reads")

ctr = 0

#write umi/bc to dict
for fastq_R1 in fastq_parser:
    bc = str(fastq_R1.seq)[barcodes_struct['BC_start']:barcodes_struct['BC_end']]      #barcode
    umi = str(fastq_R1.seq)[barcodes_struct['UMI_start']:barcodes_struct['UMI_end']]    #umi
    if(bc==''):
        sys.SystemExit('UMI empty for read {}.\n The barcode is: {}.\nWhole entry is:{}'.format(fastq_R1.id, fastq_R1.seq,fastq_R1))
    #repair bc
    map=mapping[bc]
    found="ref" in map
    true_bc=map["ref"] if found else bc
    #write
    read_barcodes[fastq_R1.id]['XC'] = true_bc
    read_barcodes[fastq_R1.id]['XM'] = umi
    ctr+=1

print("saved ", ctr, "reads")
print("start tagging bam")

ctr=0

#tag each read in bam with umi/bc from dict
for bam_read in infile_bam:
    if (discard_secondary_alignements & bam_read.is_secondary):
        continue
    if (bam_read.query_name) in read_barcodes:
        current_barcodes = read_barcodes.pop(bam_read.query_name)
        bam_read.set_tags([
            ('XC', current_barcodes['XC'], 'Z'),
            ('XM', current_barcodes['XM'], 'Z')])
        ctr+=1
    else:
        raise SystemExit('Read from mapped file is missing in reference fastq file!')
        os.remove(snakemake.output[0])

    outfile.write(bam_read)

print("tagged ", ctr, "reads")



