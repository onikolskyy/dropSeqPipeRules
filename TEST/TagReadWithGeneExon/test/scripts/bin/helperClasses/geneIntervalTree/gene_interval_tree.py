import re
import collections
import numpy as np
#todo: correct importing
from ..transcript import Transcript
from ..gene import Gene
from ncls import NCLS
from .refflat_entries import RefflatEntries


#####################################################################################
# clone for picard's OverlapDetector in Python
####################################################################################


class GeneIntervalTree:
    def __init__(self, in_refflat, bam_file):
        self.genes = GeneIntervalTree.get_genes(in_refflat, bam_file)
        self.gene_ids = list(self.genes.keys())
        starts = np.array([self.genes[self.gene_ids[i]].start for i in range(len(self.gene_ids))])
        ends = np.array([self.genes[self.gene_ids[i]].end for i in range(len(self.gene_ids))])
        self.tree = NCLS(starts, ends, np.arange(0,len(self.gene_ids)))

    @staticmethod
    def get_genes(in_refflat, bam_file):
        bam_header = bam_file.header
        SQ = bam_header["SQ"]
        ctr = 0
        print("starting to parse refflat...")
        with open(in_refflat, "r") as refflat_file:
            refflat_line = refflat_file.readline()
            # mapping of parsed gene names to corresponding transcripts
            genes = collections.defaultdict(lambda : Gene())
            while refflat_line:
                parsed_entries = RefflatEntries.get_dict(re.split(r'\t+', refflat_line.rstrip('\t')))
                # find chromosome in header
                found = False
                for pair in bam_header["SQ"]:
                    if pair["SN"] == parsed_entries["chrom"]:
                        found = True
                        break

                if found:
                    genes[parsed_entries["gene_name"]].transcripts[parsed_entries["transcription_name"]] = Transcript(
                        parsed_entries["transcription_start"],
                        parsed_entries["transcription_end"],
                        parsed_entries["coding_start"],
                        parsed_entries["coding_end"],
                        [(parsed_entries["exon_starts"][i], parsed_entries["exon_ends"][i]) for i in range(len(parsed_entries["exon_starts"]))]
                    )
                    # todo : correctness (should I take transcriptions or coding?)
                    genes[parsed_entries["gene_name"]].start = min(genes[parsed_entries["gene_name"]].start, parsed_entries["transcription_start"]+1)
                    genes[parsed_entries["gene_name"]].end = max(genes[parsed_entries["gene_name"]].end, parsed_entries["transcription_end"])
                    genes[parsed_entries["gene_name"]].strand = parsed_entries["strand"]
                refflat_line = refflat_file.readline()
                ctr = ctr+1
                if ctr % 100000 == 0:
                    print("parsed %i lines"%ctr)
        return genes


    def get_overlaps(self, block):
        #overlap_tuples = filter(lambda o : block[0] <= o[0] and block[1] >= o[1], self.tree.find_overlap(block[0], block[1]))
        overlap_tuples = self.tree.find_overlap(block[0], block[1])
        overlaps = [overlap_tuple[2] for overlap_tuple in overlap_tuples]
        return [self.gene_ids[overlap] for overlap in overlaps]

