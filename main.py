from DNAtoolkit import *
import random

rndDNAstr = ''.join([random.choice(Nucleotides) for nuc in range(50)])

DNAseq = validateSeq(rndDNAstr)

RNAseq = transcription(DNAseq)

translation = translate_seq(DNAseq)

frequency = codon_usage(DNAseq, "L")

reading_frames = gen_reading_frames(DNAseq)

print(f'\nSequence: {DNAseq}\n')
print(f'\nRNA Sequence: {RNAseq}\n')

print(f'\nTranslation: {translation}\n')

print(f'\nCodon Frequency Table: {frequency}\n')

for frame in reading_frames:
    print(f'\n{frame}\n')

test_rf = ['M', 'K', 'Q', 'T', 'S', 'P', 'Y', 'G', 'H', 'I', 'G', 'Y', 'E', 'W', 'R', '_','V', 'R']

print(f'\nProteins: {proteins_from_rf(test_rf)}\n')