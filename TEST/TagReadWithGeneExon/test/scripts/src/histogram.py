infile_bam = pysam.AlignmentFile("/input/in.bam","rb")

blocks = [block[1]-block[0] for read in infile_bam for block in read.get_blocks()]
