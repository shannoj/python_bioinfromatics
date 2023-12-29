from DNAtoolkit import *
import random



rndDNAstr = ''.join([random.choice(Nucleotides) for nuc in range(20)])

print(validateSeq(rndDNAstr))