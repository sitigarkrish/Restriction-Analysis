# Restriction Site Analysis of a DNA Sequence Using Biopython (Single Enzymes)

from Bio import SeqIO
from Bio.Restriction import EcoRI,RestrictionBatch,Analysis

record = SeqIO.read("BRCA1_cds.fasta", "fasta") # using BRCA1 Coding Sequence
sequence = record.seq

rb = RestrictionBatch([EcoRI])

analysis = Analysis(rb, sequence)

print("Coding Sequence of BRCA1 : ",sequence)
print("Length of the Coding Sequence : ",len(sequence))
print("Recognition Position for EcoRI : ",analysis.full())
print("Recognition Site : " ,EcoRI.site)
print("*"*30)

# Restriction Site Analysis of a DNA Sequence Using Biopython (Multiple Enzymes)

from Bio.Restriction import EcoRI,BamHI,HindIII,RestrictionBatch,Analysis
rb_multi = RestrictionBatch([EcoRI,BamHI,HindIII])
analysis_multi = Analysis(rb_multi,sequence)

result=analysis_multi.full()

for enzyme,positions in result.items():
    print("\nEnzyme : ",enzyme)
    print("Recognition Site : ",enzyme.site)
    if positions:
        print("Recognition Positions : ",positions)

    else:
        print("No Restriction Site Found")
