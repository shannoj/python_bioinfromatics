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
        return set(NUCLEOTIDE_BASE[self.seq_type]).issuperset(self.seq)
    
    def random_seq(self, length=10, seq_type="DNA"):
        """Generate a random DNA sequence"""
        seq = ''.join([random.choice(NUCLEOTIDE_BASE[self.seq_type]) for x in range(length)])
        self.__init__(seq, seq_type, "Randomly Generated Sequence")
    
    def get_seq_biotype(self):
        """Returns the biotype of the sequence"""
        return self.seq_type
    
    def countNucFrequency(self):
        return dict(collections.Counter(self.seq))
    
    def show_seq_info(self):
        if self.seq_type == "DNA":
            return f"[Label]: {self.label}\n[Sequence]: {self.seq}\n[Type]: {self.seq_type}\n[Length]: {len(self.seq)}\n[Is Valid]: {self.is_valid}\n[Nucleotide Frequency]: {self.countNucFrequency()}\n[GC Content]: {self.calculate_gc_content()}\n[Transcription]: {self.show_rna_transcription()}\n[Reverse Complement]: {self.generate_reverse_complement()}\n[Protein Translation]: {self.translate_seq()}\n[Reading Frames]: {self.gen_reading_frames()}\n[All Proteins]: {self.all_proteins_from_orfs()}"
        else:
            return f"[Label]: {self.label}\n[Sequence]: {self.seq}\n[Type]: {self.seq_type}\n[Length]: {len(self.seq)}\n[Is Valid]: {self.is_valid}\n[Nucleotide Frequency]: {self.countNucFrequency()}\n[GC Content]: {self.calculate_gc_content()}\n[Reverse Transcription]: {self.show_dna_transcription()}\n[Reverse Complement]: {self.generate_reverse_complement()}\n[Protein Translation]: {self.translate_seq()}\n[Reading Frames]: {self.gen_reading_frames()}\n[All Proteins]: {self.all_proteins_from_orfs()}"
    
    def generate_rna(self):
        """Convert DNA sequence into RNA sequence"""
        if self.seq_type == "DNA":
            self.seq_type = "RNA"
            return self.seq.replace("T", "U")
        if self.seq_type == "RNA":
            return "Already an RNA sequence"
        return "Not a DNA sequence"
    
    def show_rna_transcription(self):
        """Show RNA transcription"""
        if self.seq_type == "DNA":
            return self.seq.replace("T", "U")
        if self.seq_type == "RNA":
            return "Already an RNA sequence"
        return "Not a DNA sequence"
    
    def generate_dna(self):
        """Convert RNA sequence into DNA sequence"""
        if self.seq_type == "RNA":
            self.seq_type = "DNA"
            return self.seq.replace("U", "T")
        if self.seq_type == "DNA":
            return "Already a DNA sequence"
        return "Not an RNA sequence"
    
    def show_dna_transcription(self):
        """Show DNA transcription"""
        if self.seq_type == "RNA":
            return self.seq.replace("U", "T")
        if self.seq_type == "DNA":
            return "Already a DNA sequence"
        return "Not an RNA sequence"
    
    def generate_reverse_complement(self):
        """Generate reverse complement of the given DNA sequence"""
        if self.seq_type == "DNA":
            mapping = str.maketrans('ATCG', 'TAGC')
            return self.seq.translate(mapping)[::-1]
        else:
            mapping = str.maketrans('AUCG', 'UAGC')
            return self.seq.translate(mapping)[::-1]
    
    def calculate_gc_content(self):
        """Calculate the GC content of a DNA sequence"""
        return round((self.seq.count('C') + self.seq.count('G')) / len(self.seq) * 100)
    
    def calculate_gc_content_subsec(self, k=20):
        """Calculate the GC content of a DNA sequence subsection"""
        res = []
        for i in range(0, len(self.seq) - k + 1, k):
            subseq = self.seq[i:i + k]
            res.append(self.calculate_gc_content(subseq))
        return res
    
    def translate_seq(self, init_pos=0):
        """Translate the DNA sequence into an aminoacid sequence"""
        if self.seq_type == "DNA":
            return [DNA_Codons[self.seq[pos:pos + 3]] for pos in range(init_pos, len(self.seq) - 2, 3)]
        else:
            return [RNA_Codons[self.seq[pos:pos + 3]] for pos in range(init_pos, len(self.seq) - 2, 3)]
    
    def codon_usage(self, aminoacid):
        """Provide the frequency of each codon encoding a given aminoacid in a DNA sequence"""
        tmpList = []
        if self.seq_type == "DNA":
            for i in range(0, len(self.seq) - 2, 3):
                if DNA_Codons[self.seq[i:i + 3]] == aminoacid:
                    tmpList.append(self.seq[i:i + 3])
        else:
            for i in range(0, len(self.seq) - 2, 3):
                if RNA_Codons[self.seq[i:i + 3]] == aminoacid:
                    tmpList.append(self.seq[i:i + 3])         
        freqDict = dict(collections.Counter(tmpList))
        totalWight = sum(freqDict.values())
        for seq in freqDict:
            freqDict[seq] = round(freqDict[seq] / totalWight, 2)
        return freqDict
    
    def gen_reading_frames(self):
        """Generate the six reading frames of a DNA sequence, including the reverse complement"""
        if self.seq_type == "DNA":
            frames = []
            frames.append(self.translate_seq(0))
            frames.append(self.translate_seq(1))
            frames.append(self.translate_seq(2))
            reverse_complement = self.generate_reverse_complement()
            reverse_seq = bio_seq(reverse_complement, "DNA", "Reverse Complement")
            frames.append(reverse_seq.translate_seq(0))
            frames.append(reverse_seq.translate_seq(1))
            frames.append(reverse_seq.translate_seq(2))
            return frames
        return "Not a DNA sequence"
    
    def proteins_from_rf(self, aa_seq):
        """Compute all possible proteins in an aminoacid sequence and return a list of possible proteins"""
        if self.seq_type == "DNA":
            current_prot = []
            proteins = []
            for aa in aa_seq:
                if aa == "_":
                    if current_prot:
                        for p in current_prot:
                            proteins.append(p)
                        current_prot = []
                else:
                    if aa == "M":
                        current_prot.append("")
                    for i in range(len(current_prot)):
                        current_prot[i] += aa
            return proteins
        return "Not a DNA sequence"
    
    def all_proteins_from_orfs(self, startReadPos=0, endReadPos=0, ordered=False):
        """Compute all possible proteins for all open reading frames"""
        """Proteins must be in the same order of the sequence"""
        if self.seq_type == "DNA":
            if endReadPos > startReadPos:
                rfs = self.gen_reading_frames(self.seq[startReadPos: endReadPos])
            else:
                rfs = self.gen_reading_frames()
            res = []
            for rf in rfs:
                prots = self.proteins_from_rf(rf)
                for p in prots:
                    res.append(p)
            if ordered:
                return sorted(res, key=len, reverse=True)
            return res
        return "Not a DNA sequence"