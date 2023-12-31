from BIOseq import bio_seq
from BIOstruct import *

insulin_dna = bio_seq(insulin_dna_str, "DNA", "Insulin")

print(insulin_dna.show_seq_info())
