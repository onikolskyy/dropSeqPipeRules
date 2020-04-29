import mmap
import os
import gzip
from multiprocessing import RawArray, Pool
from collections import defaultdict, Counter

bases = ['T', 'G', 'A', 'C', 'N']
BUFFER = {}

BARCODE_LENGTH = snakemake.params["cell_barcode_length"]
UMI_LENGTH = snakemake.params["umi_barcode_length"]


def extract_barcodes(seq, bc_len, umi_len):
    return seq[:bc_len], seq[bc_len:bc_len + umi_len]


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
        for changed_base in filter(lambda b: b != sequence[i], bases):
            seq = subst_base(sequence, changed_base, i)
            possible.add(seq)
    return possible


def init_reader(rawarr):
    print("reader init")
    BUFFER["buf"] = rawarr






# mapping id -> (UMI,bc)
tags_for_id = defaultdict(lambda: {"UMI": "", "cellBC": ""})

# mapping barcode -> [id1,id2,..]
ids_for_barcode = defaultdict(lambda: set())

# mapping for csv
mapping = defaultdict(lambda: set())

# counter of reads for each bc
barcode_counts = Counter()

# mmap and unzip fastq
fastgz = os.open(snakemake.input["fastq"], os.O_RDONLY)
mm_fastqgz = mmap.mmap(fastgz, 0, prot=mmap.PROT_READ)
b_fastq = gzip.GzipFile(mode="r", fileobj=mm_fastqgz).read()
print(type(b_fastq))
lines_fastq = len(b_fastq.decode().split("\n"))
num_reads = int(lines_fastq / 4)
print("ungzipped and mmaped")
# save to shared memory
raw = RawArray('b', len(b_fastq))
raw[:] = b_fastq
print("built rawarray")
print(raw)


def reader(startend):
    print("starting worker")
    bfastq = bytearray(raw)
    lines = bfastq.decode().split('\n')
    start = startend[0]
    end = startend[1]

    print("Worker work")

    tags_for_id = defaultdict(lambda: {"UMI": "", "cellBC": ""})

    line_ctr = 4

    for l in range(start, end):
        line = lines[l]
        if line_ctr % 4 == 0:
            read_id = line.split(' ')[0][1:]
        elif line_ctr % 5 == 0:
            bc, umi = extract_barcodes(line, BARCODE_LENGTH,
                                       UMI_LENGTH)
            tags_for_id[read_id]["UMI"] = umi
            tags_for_id[read_id]["cellBC"] = bc
            ids_for_barcode[bc].add(read_id)
            barcode_counts[bc] += 1
        line_ctr = line_ctr + 1 if line_ctr < 7 else 4

    return (tags_for_id, ids_for_barcode)


# make chunks
size_chunk = int(num_reads / snakemake.threads)
chunks = []
for i in range(snakemake.threads):
    if i == snakemake.threads - 1:
        chunks.append((4 * i * size_chunk, num_reads * 4))
    else:
        chunks.append((4 * (i) * size_chunk, 4 * (i + 1) * size_chunk))

print(chunks)

with Pool(processes=snakemake.threads, initializer=init_reader, initargs=raw) as pool:
    print("start pool")
    results = pool.map(reader, iterable=chunks)
    print(results)

exit()













sorted_counts = sorted(barcode_counts.values(), reverse=True)
threshold = sorted_counts[snakemake.params['num_cells']]
final_barcodes = set([
    x for x, y in barcode_counts.items() if y > threshold])

print("start mapping")
ctr = 0
for cell, ids in ids_for_barcode.items():
    if cell in final_barcodes:
        pass
    else:
        possible = generate_possible(cell)
        match = possible.intersection(final_barcodes)
        if len(match) == 1:
            match = next(iter(match))
            for read_id in ids:
                tags_for_id[read_id]["cellBC"] = match
                mapping[match].add(cell)
        else:
            for read_id in ids:
                del tags_for_id[read_id]
    ctr += 1
    if ctr % 10 ** 6 == 0:
        print("processed", ctr)

out_csv = open(snakemake.output["mapping"], "w")

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