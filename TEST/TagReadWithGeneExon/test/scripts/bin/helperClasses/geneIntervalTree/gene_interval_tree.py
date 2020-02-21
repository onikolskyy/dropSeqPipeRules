import re
import collections
#todo: correct importing
from ..transcript import Transcript
from ..gene import Gene
from ncls import NCLS
from .refflat_entries import RefflatEntries


#####################################################################################
# clone for picard's OverlapDetector in Python
####################################################################################


class GeneIntervalTree:
    def __init__(self, in_refflat):
        self.genes = GeneIntervalTree.get_genes(in_refflat)
        genes_ids_list = [gene_name for gene_name in self.genes.keys()]
        starts = [self.genes[genes_ids_list[i]]["start"] for i in range(len(genes_ids_list))]
        ends = [self.genes[genes_ids_list[i]]["end"] for i in range(len(genes_ids_list))]
        self.tree = NCLS(starts, ends, genes_ids_list)

    @staticmethod
    def get_genes(in_refflat):
        with open(in_refflat, "r") as refflat_file:
            refflat_line = refflat_file.readline()
            # mapping of parsed gene names to corresponding transcripts
            genes = collections.defaultdict(Gene())
            while refflat_line:
                parsed_entries = RefflatEntries.get_dict(re.split(r'\t+', refflat_line.rstrip('\t')))
                genes[parsed_entries["gene_name"]].transcripts[parsed_entries["transcription_name"]] = Transcript(
                    parsed_entries["transcription_start"],
                    parsed_entries["transcription_end"],
                    parsed_entries["coding_start"],
                    parsed_entries["coding_end"],
                    [(parsed_entries["exon_starts"][i], parsed_entries["exon_ends"][i]) for i in range(len(parsed_entries["exon_starts"]))]
                )
                #todo : correctness (should I take transcriptions or coding?)
                genes[parsed_entries["gene_name"]].start = min(genes[parsed_entries["gene_name"]]['start'], parsed_entries["transcription_start"])
                genes[parsed_entries["gene_name"]].end = max(genes[parsed_entries["gene_name"]]['end'], parsed_entries["transcription_end"])
                genes[parsed_entries["gene_name"]].strand = parsed_entries["strand"]
                refflat_line = refflat_file.readline()
        return genes

    def get_overlaps(self, block):
        return self.tree.find_overlap(block[0], block[1])

