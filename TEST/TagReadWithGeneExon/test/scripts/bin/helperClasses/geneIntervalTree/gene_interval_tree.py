import re
import pandas as pd
import numpy as np
#todo: correct importing
from ..transcript import Transcript
from ..locus_function import LocusFunctions
from ..gene import Gene
from ncls import NCLS
from .refflat_entries import RefflatEntries
import math


#####################################################################################
# clone for picard's OverlapDetector in Python
####################################################################################


class GeneIntervalTree:
    def __init__(self, in_refflat, bam_file):
        self.I_starts = []
        self.I_ends = []
        self.LF = []
        self.REFs = []
        self.genes = []

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

    def add_Locus(self,start, end, lf, ref, gene):
        if start > end:
            return

        if start < 0 :
            start = 0

        self.I_starts.append(start)
        self.I_ends.append(end)
        self.LF.append(lf)
        self.REFs.append(ref)
        self.genes.append(gene)

    def gene_to_tree(self, genes_dict):
        for chrom, genes in genes_dict.items():
            for gene_id, gene in genes:
                for transcript_name, transcript in gene.transcripts:
                    self.add_Locus(gene.start-1, transcript.transcription_start+1,LocusFunctions.INTERGENIC,chrom, gene.name)
                    self.add_Locus(transcript.transcription_end-1,gene.end+1,LocusFunctions.INTERGENIC,chrom,gene.name)  ##todo

                    for i in range(len(transcript.exons)):
                        exons = transcript.exons

                        if i < len(exons) - 1 and exons[i][0] > transcript.coding_start:
                            self.add_Locus(exons[i][1], exons[i + 1][0], LocusFunctions.INTRONIC, chrom, gene.name)#9

                        if exons[i][0] < transcript.coding_start:
                            self.add_Locus(transcript.coding_start-1,exons[i][1]+1,LocusFunctions.CODING,chrom,gene.name)
                            self.add_Locus(exons[i][0]-1,transcript.coding_start,LocusFunctions.CODING,chrom,gene.name)
                            continue
                        if exons[i][1] > transcript.coding_end:
                            self.add_Locus(transcript.coding_end+1, exons[1]+1, LocusFunctions.UTR, chrom, gene.name)
                            self.add_Locus(exons[0]-1, transcript.coding_end+1, LocusFunctions.CODING, chrom, gene.name)
                            continue
                        self.add_Locus(exons[i][0]-1,exons[i][1]+1,LocusFunctions.CODING, chrom, gene.name)


            self.intervals = pd.DataFrame({
             "start": self.I_starts,
             "end": self.I_ends,
             "ref": self.REFs,
             "LF": self.LF,
             "G": self.genes,
             "index": np.arange(len(self.I_starts))
            })

            self.tree = NCLS(self.I_starts,self.I_ends,np.arange(len(self.I_starts)))

    def get_overlaps(self, query):
        overlaps = self.tree.all_overlaps_both(query["starts"].values, query["ends"].values, query["ids"].values)
        return pd.DataFrame({"B":overlaps[0], "index":overlaps[1]}).merge(self.intervals, on="index")[["B","LF","G","ref"]]



