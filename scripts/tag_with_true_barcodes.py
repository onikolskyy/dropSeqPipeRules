from Bio import SeqIO
from sortedcontainers import SortedList
from collections import defaultdict
import pysam
import gzip
import sys


discard_secondary_alignements = snakemake.params['discard_secondary_alignements']
fastq_parser = SeqIO.parse(gzip.open(snakemake.input[0], "rt"), "fastq")
read_barcodes = defaultdict(lambda :{'XC':'','XM':''})
infile_bam = pysam.AlignmentFile(snakemake.input[1], "rb")
outfile = pysam.AlignmentFile(snakemake.output[0], "wb", template=infile_bam)

bins = []

def find_in_bins(bins, ref):
    for bin in bins:
        try:
            i = bin[1].index(ref)
            return i
        except:
            continue
    return -1

barcodes_struct = {
	'BC_start':snakemake.params['BC_start'],
	'BC_end':snakemake.params['BC_end'],
	'UMI_start':snakemake.params['UMI_start'],
	'UMI_end':snakemake.params['UMI_end']
	}

#construct bins
with open(snakemake.input[2],'r') as whitelist:
    for line in whitelist:
        if len(line.strip().split()) == 2:  # This means we didn't find any other linked barcode
            (reference, counts_ref) = line.strip().split()
            bins.append ([reference, SortedList()]) #Push an empty bin and tag it with reference barcode
        else:
            (reference, extended_ref, counts_ref, counts_ext) = line.strip().split()
            list = SortedList()
            for barcode in extended_ref.split(','):
                list.add(barcode)
            bins.append([reference, list])

#write umi/bc to dict
for fastq_R1 in fastq_parser:
    ref = str(fastq_R1.seq)[barcodes_struct['BC_start']:barcodes_struct['BC_end']]      #barcode
    umi = str(fastq_R1.seq)[barcodes_struct['UMI_start']:barcodes_struct['UMI_end']]    #umi
    if(ref==''):
        sys.SystemExit('UMI empty for read {}.\n The barcode is: {}.\nWhole entry is:{}'.format(fastq_R1.id, fastq_R1.seq,fastq_R1))

    #repair bc
    index = find_in_bins(bins, ref)
    true_ref = bins[index][0] if index >= 0 else ref
    #write
    read_barcodes[fastq_R1.id]['XC'] = ref
    read_barcodes[fastq_R1.id]['XM'] = umi

#tag each read in bam with umi/bc from dict
for bam_read in infile_bam:
    if (discard_secondary_alignements & bam_read.is_secondary):
        continue
    if (bam_read.query_name) in read_barcodes:
        current_barcodes = read_barcodes.pop(bam_read.query_name)
        bam_read.set_tags([
            ('XC', current_barcodes['XC'], 'Z'),
            ('XM', current_barcodes['XM'], 'Z')])
    else:
        raise SystemExit('Read from mapped file is missing in reference fastq file!')
        os.remove(snakemake.output[0])
    outfile.write(bam_read)



