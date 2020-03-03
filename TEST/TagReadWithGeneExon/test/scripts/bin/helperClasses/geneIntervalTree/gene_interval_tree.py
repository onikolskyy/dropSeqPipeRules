import re
import collections
import numpy as np
#todo: correct importing
from ..transcript import Transcript
from ..gene import Gene
from ncls import NCLS
from .refflat_entries import RefflatEntries
import math


#####################################################################################
# clone for picard's OverlapDetector in Python
####################################################################################


class GeneIntervalTree:
    def __init__(self, in_refflat, bam_file):
        _genes = GeneIntervalTree.get_genes(in_refflat, bam_file)
        self.trees = {}
        for chrom, genes in _genes.items():
            self.trees[chrom] = {}
            ids = [gene_id for gene_id, gene in genes.items()]
            starts = np.array([genes[ids[i]].start-1 for i in range(len(genes))])
            ends = np.array([genes[ids[i]].end+1 for i in range(len(genes))])
            self.trees[chrom]["ncls"] = NCLS(starts, ends, np.arange(0, len(ids)))
            self.trees[chrom]["ids"] = ids
            self.trees[chrom]["genes"] = genes

    @staticmethod
    def get_genes(in_refflat, bam_file):
        bam_header = bam_file.header
        SQ = bam_header["SQ"]
        ctr = 0

        genes = {}
        parsed_mapping = {}

        with open(in_refflat, "r") as refflat_file:
            refflat_line = refflat_file.readline()
            # mapping of parsed gene names to corresponding transcripts

            while refflat_line:
                parsed_entries = RefflatEntries.get_dict(re.split(r'\t+', refflat_line.rstrip('\t')))
                # find chromosome in header
                found = False
                for pair in bam_header["SQ"]:
                    if pair["SN"] == parsed_entries["chrom"]:
                        found = True
                        break

                if found:
                    if parsed_entries["gene_name"] in parsed_mapping.keys():
                        parsed_gene = parsed_mapping[parsed_entries["gene_name"]]
                        if parsed_gene["chrom"] != parsed_entries["chrom"] or parsed_gene["strand"] != parsed_entries["strand"]:
                            parsed_gene["mismatch"] = True              # do not retain this gene
                        else:
                            parsed_gene["transcripts"][parsed_entries["transcription_name"]] = Transcript(
                                parsed_entries["transcription_start"],
                                parsed_entries["transcription_end"],
                                parsed_entries["coding_start"],
                                parsed_entries["coding_end"],
                                [(parsed_entries["exon_starts"][i], parsed_entries["exon_ends"][i]) for i in range(len(parsed_entries["exon_starts"]))]
                            )
                            parsed_gene["start"] = min(parsed_gene["start"], parsed_entries["transcription_start"])
                            parsed_gene["end"] = max(parsed_gene["end"], parsed_entries["transcription_end"])
                    else:
                        parsed_mapping[parsed_entries["gene_name"]] = {
                            "start" : parsed_entries["transcription_start"],
                            "end" : parsed_entries["transcription_end"],
                            "strand" : parsed_entries["strand"],
                            "chrom" : parsed_entries["chrom"],
                            "mismatch" : False,
                            "transcripts" : {}
                        }
                        parsed_mapping[parsed_entries["gene_name"]]["transcripts"][parsed_entries["transcription_name"]] = Transcript(
                                parsed_entries["transcription_start"],
                                parsed_entries["transcription_end"],
                                parsed_entries["coding_start"],
                                parsed_entries["coding_end"],
                                [(parsed_entries["exon_starts"][i], parsed_entries["exon_ends"][i]) for i in range(len(parsed_entries["exon_starts"]))]
                        )

                refflat_line = refflat_file.readline()
                ctr = ctr + 1
                if ctr % 100000 == 0:
                    print("parsed %i lines" % ctr)


            # save the parsed genes

            for gene_id, parsed_gene in parsed_mapping.items():
                if parsed_gene["mismatch"]:
                    continue
                gene = Gene()
                gene.name = gene_id
                gene.start = parsed_gene["start"]
                gene.end = parsed_gene["end"]
                gene.chrom = parsed_gene["chrom"]
                gene.strand = parsed_gene["strand"]
                for transcript_name, transcript in parsed_gene["transcripts"].items():
                    gene.transcripts[transcript_name] = transcript
                if gene.chrom in genes:
                    genes[gene.chrom][gene_id] = gene
                else:
                    genes[gene.chrom] = {}
                    genes[gene.chrom][gene_id] = gene
        return genes

    def get_overlaps_by_ref(self, block, ref):
        tree_obj = self.trees[ref]
        overlap_tuples = tree_obj["ncls"].find_overlap(block[0], block[1])
        overlaps = [overlap_tuple[2] for overlap_tuple in overlap_tuples]
        genes = tree_obj["genes"]
        gene_ids = tree_obj["ids"]
        result =  [genes[gene_ids[i]] for i in overlaps]
        return result

    def get_all_overlaps(self,block):
        result = {}
        for ref,tree in self.trees.items():
            overlapped = self.get_overlaps_by_ref(block,ref)
            result[ref] = [o.name for o in overlapped]

        return result