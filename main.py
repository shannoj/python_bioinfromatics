from DNAtoolkit import *
import random

rndDNAstr = ''.join([random.choice(Nucleotides) for nuc in range(50)])

DNAseq = validateSeq(rndDNAstr)

RNAseq = transcription(DNAseq)

translation = translate_seq(DNAseq)

print(f'\nSequence: {DNAseq}\n')
print(f'\nRNA Sequence: {RNAseq}\n')

print(f'\nTranslation: {translation}\n')