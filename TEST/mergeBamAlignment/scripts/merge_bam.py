import pysam
import gzip
from collections import defaultdict
import os
import mmap


def extract_barcodes(seq, bc_len, umi_len):
    return seq[:bc_len], seq[bc_len:bc_len + umi_len]


discard_secondary_alignements = snakemake.params['discard_secondary_alignements']

infile_bam = pysam.AlignmentFile(snakemake.input["bam   "], "rb")
outfile = pysam.AlignmentFile(snakemake.output[0], "wb", template=infile_bam)

fastgz = os.open(snakemake.input["fastq"], os.O_RDONLY)
mm_fastqgz = mmap.mmap(fastgz, 0, prot=mmap.PROT_READ)
b_fastq = gzip.GzipFile(mode="r", fileobj=mm_fastqgz).read()
lines_fastq = b_fastq.decode().split("\n")

read_barcodes = defaultdict(lambda: {'XC': '', 'XM': ''})

line_ctr = 4

for line in lines_fastq:
    if line_ctr % 4 == 0:
        read_id = line.split(' ')[0][1:]
    elif line_ctr % 5 == 0:
        bc, umi = extract_barcodes(line, snakemake.params["cell_barcode_length"], snakemake.params["umi_barcode_length"])
        read_barcodes[read_id]["XC"] = bc
        read_barcodes[read_id]["XM"] = umi
    line_ctr = line_ctr + 1 if line_ctr < 7 else 4


for bam_read in infile_bam:
    if (discard_secondary_alignements & bam_read.is_secondary):
        continue
    if (bam_read.query_name) in read_barcodes:
        current_barcodes = read_barcodes.pop(bam_read.query_name)
        bam_read.set_tags([
            ('XC', current_barcodes['XC'], 'Z'),
            ('XM', current_barcodes['XM'], 'Z')])
    else:
        if (bam_read.query_name) not in read_barcodes:
            raise SystemExit('Read from mapped file is missing in reference fastq file!')
            os.remove(snakemake.output[0])
    outfile.write(bam_read)
