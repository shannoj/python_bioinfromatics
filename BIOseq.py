from BIOstruct import *
import collections

class bio_seq:
    def __init__(self, seq="ATCG", seq_type="DNA", label="No Label"):
        self.seq = seq.upper()
        self.label = label
        self.seq_type = seq_type
        self.is_valid = self.validate()
        assert self.is_valid, f"Provided Data is not a valid {self.seq_type} sequence"


    def validate(self):
        """Check the sequence for valid alphabet characters"""
        return set(Nucleotides).issuperset(self.seq)
    
    def countNucFrequency(self):
        return dict(collections.Counter(self.seq))
    
    def show_seq_info(self):
        return f"[Label]: {self.label}\n[Sequence]: {self.seq}\n[Type]: {self.seq_type}\n[Length]: {len(self.seq)}\n[Is Valid]: {self.is_valid}\n[Nucleotide Frequency]: {self.countNucFrequency()}"