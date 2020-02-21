from .locus_function import LocusFunctions


class Transcript:
    def __init__(self, transcription_start, transcription_end, coding_start, coding_end, exons):
        self.transcription_start = transcription_start
        self.transcription_end = transcription_end
        self.coding_start = coding_start
        self.coding_end = coding_end
        self.exons = exons

    def assign_locus_function_for_range(self, start, locus_functions):
        low = max(start, self.transcription_start)
        high = min(self.transcription_end, start+len(locus_functions)-1)
        for i in range(low, high):
            if locus_functions[i - start] > LocusFunctions.CODING:
                if self.inExon(i):
                    if self.utr(i):
                        lf = LocusFunctions.UTR
                    else:
                        lf = LocusFunctions.CODING
                else:
                    lf = LocusFunctions.INTRONIC

                if lf > locus_functions[i - start]:
                    locus_functions[i - start] = lf

    def inExon(self, pos):
        for exon in self.exons:
            if exon[0] >= pos and exon[1] <= pos:
                return True
        return False

    def utr(self, pos):
        return pos < self.coding_start or pos > self.coding_end
