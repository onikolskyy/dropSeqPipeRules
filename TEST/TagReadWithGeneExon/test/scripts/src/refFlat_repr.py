import re
import pandas as pd

class RefflatEntries:
    entries = {
        "gene_name": {
            "index": 0,
            "convert_func": lambda gene_name: gene_name
        },
        "transcription_name": {
            "index": 1,
            "convert_func": lambda transcription_name: transcription_name
        },
        "ref": {
            "index": 2,
            "convert_func": lambda ref: ref
        },
        "strand": {
            "index": 3,
            "convert_func": lambda strand: -1 if strand == "-" else 1
        },
        "transcription_start": {
            "index": 4,
            "convert_func": lambda start: int(start)+1
        },
        "transcription_end": {
            "index": 5,
            "convert_func": lambda end: int(end)
        },
        "coding_start": {
            "index": 6,
            "convert_func": lambda start: int(start)
        },
        "coding_end": {
            "index": 7,
            "convert_func": lambda end: int(end)
        },
        "num_exons": {
            "index": 8,
            "convert_func": lambda num: int(num)
        },
        "exon_starts" : {
            "index": 9,
            "convert_func": lambda exon_starts: [int(exon_start) for exon_start in exon_starts.split(",")[:-1]] if (len(exon_starts) > 0 and not re.compile(" +\B").match(exon_starts)) else []
        },
        "exon_ends" : {
            "index": 10,
            "convert_func": lambda exon_ends: [int(exon_end) for exon_end in exon_ends.split(",")[:-1]] if (len(exon_ends) > 0 and not re.compile(" +\B").match(exon_ends)) else []
        },
    }

    @staticmethod
    def get_dict(entries_list):

        # handle error
        if len(entries_list) < 11:
            raise Exception("refflat file seems to be damaged: not enough columns on a row")

        res = {}

        for i in range(len(entries_list)):
            for entry_name in RefflatEntries.entries.keys():
                if RefflatEntries.entries[entry_name]["index"] == i:
                    convert_func = RefflatEntries.entries[entry_name]["convert_func"]
                    res[entry_name] = convert_func(entries_list[i])
                    break
        if not len(res["exon_starts"]) == len(res["exon_ends"]):
            raise Exception("Error while parsing refflat: length of exon starts does not match the length of exon ends...")


        return res


class RefFlatParsed:
    def __init__(self,in_refflat, bam_file):
        bam_header = bam_file.header
        self.parsed_mapping = {}

        with open(in_refflat, "r") as refflat_file:
            # save the parsed genes
            refflat_line = refflat_file.readline()
            # mapping of parsed gene names to corresponding transcripts

            while refflat_line:
                parsed_entries = RefflatEntries.get_dict(re.split(r'\t+', refflat_line.rstrip('\t')))
                # find ref in header
                found = False
                for pair in bam_header["SQ"]:
                    if pair["SN"] == parsed_entries["ref"]:
                        found = True
                        break

                if found:
                    if not parsed_entries["ref"] in self.parsed_mapping.keys():
                        self.parsed_mapping[parsed_entries["ref"]] = {}

                    parsed_mapping_for_ref = self.parsed_mapping[parsed_entries["ref"]]

                    if not parsed_entries["gene_name"] in parsed_mapping_for_ref:
                        parsed_mapping_for_ref[parsed_entries["gene_name"]] = {
                            "transcripts" : [],
                            "start" : float('inf'),
                            "end" : -float('inf')
                        }
                    parsed_mapping_for_gene = parsed_mapping_for_ref[parsed_entries["gene_name"]]

                    transcript = {
                        "coding_start": parsed_entries["coding_start"],
                        "coding_end" : parsed_entries["coding_end"],
                        "transcription_start" : parsed_entries["transcription_start"],
                        "transcription_end" : parsed_entries["transcription_end"],
                        "exons" : [(parsed_entries["exon_starts"][i], parsed_entries["exon_ends"][i])\
                                   for i in range(len(parsed_entries["exon_starts"]))]
                    }

                    parsed_mapping_for_gene["transcripts"].append(transcript)

                    parsed_mapping_for_gene["start"] = min(parsed_mapping_for_gene["start"],parsed_entries["transcription_start"] )
                    parsed_mapping_for_gene["end"] = max(parsed_mapping_for_gene["start"],parsed_entries["transcription_end"] )

                    # proceed to next line

                refflat_line = refflat_file.readline()
        print(self.parsed_mapping.keys())

    def as_intervals(self,ref):

        if ref not in self.parsed_mapping:
            return None

        parsed_mapping_for_ref = self.parsed_mapping[ref]

        start = []
        end = []
        gene = []
        LF = []

        ctr_coding = 0
        ctr_utr = 0

        for gene_name, parsed_gene in parsed_mapping_for_ref.items():
            for transcript in parsed_gene["transcripts"]:
                for i in range(len(transcript["exons"])):
                    exon = transcript["exons"][i]
                    #todo: can coding region be inside an exon?
                    if exon[0] < transcript["coding_start"] or exon[1] > transcript["coding_end"]:
                        # UTR
                        ctr_utr+=1
                        if exon[0] < transcript["coding_start"]:
                            # exon preceeds coding region
                            start.append(exon[0])
                            if exon[1] <= transcript["coding_start"]:
                                end.append(exon[1])
                            else:
                                end.append(transcript["coding_start"])
                            LF.append(2)
                            gene.append(gene_name)
                        else:
                            # exon exceeds coding region
                            end.append(exon[1])
                            if exon[0] >= transcript["coding_end"]:
                                start.append(exon[0])
                            else:
                                start.append(transcript["coding_end"])
                            LF.append(2)
                            gene.append(gene_name)
                    else:
                        ctr_coding +=1
                        #CODING
                        start.append(exon[0])
                        end.append(exon[1])
                        gene.append(gene_name)
                        LF.append(3)

                    #INTRONIC
                    if i < len(transcript["exons"])-1:
                        start.append(transcript["exons"][i][1])
                        end.append(transcript["exons"][i+1][1])
                        gene.append(gene_name)
                        LF.append(1)
            # Intergenic
            gene.append(gene_name)
            start.append(parsed_gene["start"])
            end.append(parsed_gene["end"])
            LF.append(0)

        print("--->coding: %i, UTR: %i"%(ctr_coding,ctr_utr))

        return pd.DataFrame({"gene": gene, "start": start, "end":end, "LF": LF})

