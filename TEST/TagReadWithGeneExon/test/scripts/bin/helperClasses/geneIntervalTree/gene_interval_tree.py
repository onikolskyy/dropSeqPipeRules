import re
import pandas as pd
import numpy as np
import ray
from ..transcript import Transcript
from ..locus_function import LocusFunctions
from ..gene import Gene
from ncls import NCLS
from .refflat_entries import RefflatEntries
from time import time
import math


#####################################################################################
# clone for picard's OverlapDetector in Python
####################################################################################
ray.init()

@ray.remote
class GeneIntervalTree:
    def __init__(self, in_refflat, bam_file):
        self.I_starts = []
        self.I_ends = []
        self.LF = []
        self.REFs = []
        self.genes = []

        self.data = {}
        _genes = GeneIntervalTree.get_genes(in_refflat, bam_file)
        self.gene_to_tree(_genes)

    @staticmethod
    def get_genes(in_refflat, bam_file):
        bam_header = bam_file.header
        SQ = bam_header["SQ"]
        ctr = 0

        genes = {}
        with open(in_refflat, "r") as refflat_file:
            # save the parsed genes

            parsed_mapping = {}

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
                        if parsed_gene["chrom"] != parsed_entries["chrom"] or parsed_gene["strand"] != parsed_entries[
                            "strand"]:
                            parsed_gene["mismatch"] = True  # do not retain this gene
                        else:
                            parsed_gene["transcripts"][parsed_entries["transcription_name"]] = Transcript(
                                parsed_entries["transcription_start"],
                                parsed_entries["transcription_end"],
                                parsed_entries["coding_start"],
                                parsed_entries["coding_end"],
                                [(parsed_entries["exon_starts"][i], parsed_entries["exon_ends"][i]) for i in
                                 range(len(parsed_entries["exon_starts"]))]
                            )
                            parsed_gene["start"] = min(parsed_gene["start"], parsed_entries["transcription_start"])
                            parsed_gene["end"] = max(parsed_gene["end"], parsed_entries["transcription_end"])
                    else:
                        parsed_mapping[parsed_entries["gene_name"]] = {
                            "start": parsed_entries["transcription_start"],
                            "end": parsed_entries["transcription_end"],
                            "strand": parsed_entries["strand"],
                            "chrom": parsed_entries["chrom"],
                            "mismatch": False,
                            "transcripts": {}
                        }
                        parsed_mapping[parsed_entries["gene_name"]]["transcripts"][
                            parsed_entries["transcription_name"]] = Transcript(
                            parsed_entries["transcription_start"],
                            parsed_entries["transcription_end"],
                            parsed_entries["coding_start"],
                            parsed_entries["coding_end"],
                            [(parsed_entries["exon_starts"][i], parsed_entries["exon_ends"][i]) for i in
                             range(len(parsed_entries["exon_starts"]))]
                        )

                refflat_line = refflat_file.readline()
        return parsed_mapping

    def add_Locus(self,start, end, lf, ref, gene):
        if start > end:
            return

        if start < 0 :
            start = 0

        if ref not in self.data:
            self.data[ref] = {
                "I_starts" : [],
                "I_ends" : [],
                "LF" : [],
                "genes" : [],
            }

        self.data[ref]["I_starts"].append(start)
        self.data[ref]["I_ends"].append(end)
        self.data[ref]["LF"].append(lf)
        self.data[ref]["REFs"].append(ref)
        self.data[ref]["genes"].append(gene)

    def gene_to_tree(self, genes_dict):
        for ref, genes in genes_dict.items():
            for gene_id, gene in genes.items():
                for transcript_name, transcript in gene.transcripts.items():
                    self.add_Locus(gene.start-1, transcript.transcription_start+1,LocusFunctions.INTERGENIC,ref, gene.name)
                    self.add_Locus(transcript.transcription_end-1,gene.end+1,LocusFunctions.INTERGENIC,ref,gene.name)  ##todo

                    for i in range(len(transcript.exons)):
                        exons = transcript.exons

                        if i < len(exons) - 1 and exons[i][0] > transcript.coding_start:
                            self.add_Locus(exons[i][1], exons[i + 1][0], LocusFunctions.INTRONIC, ref, gene.name)#9

                        if exons[i][0] < transcript.coding_start:
                            self.add_Locus(transcript.coding_start-1,exons[i][1]+1,LocusFunctions.CODING,ref,gene.name)
                            self.add_Locus(exons[i][0]-1,transcript.coding_start,LocusFunctions.CODING,ref,gene.name)
                            continue
                        if exons[i][1] > transcript.coding_end:
                            self.add_Locus(transcript.coding_end+1, exons[i][1]+1, LocusFunctions.UTR, ref, gene.name)
                            self.add_Locus(exons[i][0]-1, transcript.coding_end+1, LocusFunctions.CODING, ref, gene.name)
                            continue
                        self.add_Locus(exons[i][0]-1,exons[i][1]+1,LocusFunctions.CODING, ref, gene.name)

        self.data[ref]["ncl"] = NCLS(
            np.array(self.data[ref]["I_starts"]),
            np.array(self.data[ref]["I_ends"]),
            np.arange(len(self.data[ref]["I_starts"]))
        )

        self.data[ref]["intervals"] = pd.DataFrame({
         "start": self.data[ref]["I_starts"],
         "end": self.data[ref]["I_ends"],
         "LF": self.data[ref]["LF"],
         "G": self.data[ref]["genes"]
        # "index": np.arange(len(self.I_starts))
        })

    def get_overlaps(self, query, ref):
        t1 = time()
        overlaps = self.data[ref]["ncl"].all_overlaps_both(query["starts"].values,
                                               query["ends"].values,
                                               query["B"].values)
        t2= time()
        print("query took", t2-t1)

        return {
            "overlaps" : pd.DataFrame({"B":overlaps[0], "index":overlaps[1]}),
            "intervals" : self.data[ref]["intervals"]
        }





