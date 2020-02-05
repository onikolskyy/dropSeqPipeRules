import threading
import collections
import gzip
import time
import re

bases = ['T', 'G', 'A', 'C', 'N']


def worker(cells, index, mapping, whitelist):
    print("Worker to work on cells[%i:%i]" % (index[0], index[1]))
    for i in range(index[0], index[1]):
        find_match(cells[i], mapping, whitelist)


def find_match(cell, mapping, whitelist):
    if cell in whitelist:
            if cell not in mapping.keys():
                mapping[cell] = set()
    else:
        possible = generate_possible(cell)
        match = None

        for p in possible:
            if p in whitelist :
                if match is None:
                    mapping[p].add(cell)
                    match = (p,cell)
                else:
                    # don't allow mapping of bc to multiple true
                    mapping[match[0]].remove(cell)
                    break


def convert2string(b):
    if type(b) == str:
        return b
    else:
        return b.decode("utf-8")

def fastq_seq_iterate(infile):
    while 1:
        line1 = convert2string(infile.readline())
        line2 = convert2string(infile.readline())
        line3 = convert2string(infile.readline())
        line4 = convert2string(infile.readline())
        yield line2[:-1]


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
    possible=[]
    for i in range(len(sequence)):
        for changed_base in filter(lambda b: b!=sequence[i], bases):
            seq = subst_base(sequence, changed_base, i)
            possible.append(seq)
    return possible


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


cell_barcode_counts = collections.Counter()

ctr = 0
with gzip.open(snakemake.input[0], "rt") as infastq:
    while 1:
        line1 = convert2string(infastq.readline())
        if len(line1) == 0:
            break

        seq = convert2string(infastq.readline())
        line3 = convert2string(infastq.readline())
        line4 = convert2string(infastq.readline())

        cell, umi = extract_barcodes(seq, snakemake.params["regex"])
        cell_barcode_counts[cell] +=1
        ctr += 1
        if ctr % 500000 == 0:
            print ("parsed %i reads" % ctr)

counts = sorted(cell_barcode_counts.values(), reverse=True)

threshold = counts[snakemake.params['CELL_NUMBER']]
final_barcodes = set([
            x for x, y in cell_barcode_counts.items() if y > threshold])

#now we have true barcodes
#and can do work

mapping = collections.defaultdict(set)

cell_barcode_counts_list = list(cell_barcode_counts.keys())
n_cells_per_thread = round(len(cell_barcode_counts_list)/N_THREADS)

threads = []
sub_mappings  = []
ctr = 0

print("All cells : %i"%len(cell_barcode_counts_list))

for i in range(snakemake.params["N_THREADS"]):
    low = ctr
    high = ctr + n_cells_per_thread if i < N_THREADS-1 else len(cell_barcode_counts_list)
    sub_mapping = collections.defaultdict(set)
    sub_mappings.append(sub_mapping)
    threads.append(threading.Thread(args=(cell_barcode_counts_list, (low, high), sub_mapping, final_barcodes), target=worker))
    ctr = high

for thread in threads:
    thread.start()

for thread in threads:
    thread.join()

for sub_mapping in sub_mappings:
    for key in sub_mapping.keys():
        for cell in sub_mapping[key]:
            mapping[key].add(cell)

#write map:
out = opts["output"]
outfile = open(snakemake.output[0], "w")

for final_barcode in final_barcodes:
    corrected_barcodes = ",".join(
        sorted(mapping[final_barcode]))
    corrected_barcode_counts = ",".join(
        map(str, [cell_barcode_counts[x] for x
                  in sorted(mapping[final_barcode])]))
    outfile.write("%s\t%s\t%s\t%s\n" % (
        final_barcode, corrected_barcodes, cell_barcode_counts[final_barcode],
        corrected_barcode_counts))
outfile.close()







