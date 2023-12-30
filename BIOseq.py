from BIOstruct import *
import collections
import random

class bio_seq:
    def __init__(self, seq="ATCG", seq_type="DNA", label="No Label"):
        self.seq = seq.upper()
        self.label = label
        self.seq_type = seq_type
        self.is_valid = self.__validate()
        assert self.is_valid, f"Provided Data is not a valid {self.seq_type} sequence"


    def __validate(self):
        """Check the sequence for valid alphabet characters"""
        return set(Nucleotides).issuperset(self.seq)
    
    def random_seq(self, length=10, seq_type="DNA"):
        """Generate a random DNA sequence"""
        seq = ''.join([random.choice(Nucleotides) for x in range(length)])
        self.__init__(seq, seq_type, "Randomly Generated Sequence")
    
    def get_seq_biotype(self):
        """Returns the biotype of the sequence"""
        return self.seq_type
    
    def countNucFrequency(self):
        return dict(collections.Counter(self.seq))
    
    def show_seq_info(self):
        return f"[Label]: {self.label}\n[Sequence]: {self.seq}\n[Type]: {self.seq_type}\n[Length]: {len(self.seq)}\n[Is Valid]: {self.is_valid}\n[Nucleotide Frequency]: {self.countNucFrequency()}"