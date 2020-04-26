import mmap
import os
import sys
import pickle
import re
import gzip
import pysam
from collections import defaultdict, Counter

bases = ['T', 'G', 'A', 'C', 'N']

def save_obj(obj,name):
    with open(name, 'wb') as f:
        pickle.dump(obj,f,pickle.HIGHEST_PROTOCOL)

def extract_barcodes(seq, regex):
    pattern = re.compile(regex)
    match = pattern.match(seq)
    groupdict = match.groupdict()
    cell = ""
    umi = ""

    for k in sorted(list(groupdict)):
        span = match.span(k)
        if k.startswith("cell_"):
            cell += groupdict[k]
        if k.startswith("umi_"):
            umi += groupdict[k]
    return cell, umi


def subst_base(seq, base, i):
    b = seq[i]
    if i == 0:
        seq = base + seq[1:]
    elif i >= len(seq):
        seq = seq[:i] + base
    else:
        seq = seq[:i] + base + seq[i + 1:]
    return seq


def generate_possible(sequence):
    possible = set()
    for i in range(len(sequence)):
        for changed_base in filter(lambda b: b!=sequence[i], bases):
            seq = subst_base(sequence, changed_base, i)
            possible.add(seq)
    return possible


# mapping id -> (UMI,bc)
tags_for_id = defaultdict(lambda: {"UMI": "", "cellBC": ""})

# mapping barcode -> [id1,id2,..]
ids_for_barcode = defaultdict(lambda: set())

# mapping for csv
mapping = defaultdict(lambda: set())

# counter of reads for each bc
barcode_counts = Counter()

fastgz = os.open(snakemake.input["R1"], os.O_RDONLY)
mm_fastqgz = mmap.mmap(fastgz, 0, prot=mmap.PROT_READ)
b_fastq = gzip.GzipFile(mode="r", fileobj=mm_fastqgz).read()
lines_fastq = b_fastq.decode().split("\n")

line_ctr = 0
read_id = ""


regex = "(?P<cell_1>.{"+str(snakemake.params["cell_barcode_length"])+"})(?P<umi_1>.{"+str(snakemake.params["umi_barcode_length"])+"})"

for line in lines_fastq:
    if line_ctr % 4 == 0:
        read_id = line[1:]
    elif line_ctr % 2 == 0:
        bc, umi = extract_barcodes(line, regex)
        tags_for_id[read_id]["UMI"] = umi
        tags_for_id[read_id]["cellBC"] = bc
        ids_for_barcode[bc].add(read_id)
        barcode_counts[bc] += 1


sorted_counts = sorted(barcode_counts.values(), reverse=True)
threshold = sorted_counts[snakemake.params['num_cells']]
final_barcodes = set([
    x for x, y in barcode_counts.items() if y > threshold])


for cell, ids in ids_for_barcode.items():
    if cell in final_barcodes:
        pass
    else:
        possible = generate_possible(cell)
        match = possible.intersection(final_barcodes)
        if len(match) == 1:
            match = next(iter(match))
            for read_id in ids:
                tags_for_id[read_id]["BC"] = match
                mapping[match].add(cell)
        else:
            for read_id in ids:
                del tags_for_id[read_id]

infile_bam = pysam.AlignmentFile(snakemake.input["R2"], "rb")

out_bam = pysam.AlignmentFile(snakemake.output["repaired_bam"], "wb", template=infile_bam)

for read in infile_bam:
    if (snakemake.params["discard_secondary_alignements"] & read.is_secondary):
        continue
    if read.query_name in tags_for_id:
        read.set_tags([
            ('XC', tags_for_id["BC"], 'Z'),
            ('XM', tags_for_id["umi"], 'Z')])
    out_bam.write(read)

out_csv = open(snakemake.output["mapping"])

for final_barcode in final_barcodes:
        corrected_barcodes = ",".join(
            sorted(mapping[final_barcode]))
        corrected_barcode_counts = ",".join(
            map(str, [barcode_counts[x] for x
                      in sorted(mapping[final_barcode])]))
        out_csv.write("%s\t%s\t%s\t%s\n" % (
            final_barcode, corrected_barcodes, barcode_counts[final_barcode],
            corrected_barcode_counts))

out_csv.close()