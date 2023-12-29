import collections

Nucleotides = ["A", "T", "C", "G"]

def validateSeq(dna_seq):
    tmpseq = dna_seq.upper()
    for nuc in tmpseq:
        if nuc not in Nucleotides:
            return False
    return tmpseq

def countNucFrequency(seq):
    return dict(collections.Counter(seq))
