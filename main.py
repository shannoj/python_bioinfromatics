from BIOseq import bio_seq
from BIOstruct import *

insulin_dna = bio_seq(insulin_dna_str, "DNA", "Insulin")

rna_insulin = insulin_dna.generate_rna()

rna_insulin_obj = bio_seq(rna_insulin, "RNA", "RNA Insulin")

print(rna_insulin_obj.show_seq_info())


