import re

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
        "chrom": {
            "index": 2,
            "convert_func": lambda chrom: chrom
        },
        "strand": {
            "index": 3,
            "convert_func": lambda strand: -1 if strand == "-" else 1
        },
        "transcription_start": {
            "index": 4,
            "convert_func": lambda start: int(start)
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
            "convert_func": lambda exon_starts: [int(exon_start)+1 for exon_start in exon_starts.split(",")[:-1]] if (len(exon_starts) > 0 and not re.compile(" +\B").match(exon_starts)) else []
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
