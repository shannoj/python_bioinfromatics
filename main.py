from DNAtoolkit import *
import random

rndDNAstr = ''.join([random.choice(Nucleotides) for nuc in range(50)])

DNAseq = validateSeq(rndDNAstr)

print(countNucFrequency(rndDNAstr))