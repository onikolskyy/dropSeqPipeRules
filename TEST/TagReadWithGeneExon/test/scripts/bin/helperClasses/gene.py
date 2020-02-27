#this is a clone of Gene class from picard
import math
from .transcript import Transcript


class Gene:
    def __init__(self):
        self.transcripts = {}
        self.strand = -1
        self.start = math.inf
        self.end = -math.inf
        self.chrom = ""
        self.name = ""

    def add_transcript(self, name,
                       transcription_start, transcription_end, coding_start, coding_end, exons):
        self.transcripts[name] = Transcript(transcription_start, transcription_end, coding_start, coding_end, exons)

    def is_negative_strand(self):
        return self.strand < 0

    def verbose(self, very):
        print("strand:", self.strand, "start,end:" , self.start, self.end, "chrom", self.chrom)
        if not very:
            return
        else:
            for t_name, t in self.transcripts.items():
                print("TRANSCRIPT", t_name)
                t.verbose()