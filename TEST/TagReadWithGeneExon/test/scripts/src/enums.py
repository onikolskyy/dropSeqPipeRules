from enum import IntEnum

class Tags:
    tags_dict = {
        'GENE_FUNCTION_TAG' : "gf",
        'GENE_NAME_TAG' : "gn",
        'GENE_STRAND_TAG' : "gs",
        'READ_FUNCTION_TAG' : "XF"
    }

class LocusFunctions(IntEnum):
    CODING = 1
    INTERGENIC = 2
    INTRONIC = 3
    RIBOSOMAL = 4
    UTR = 5
