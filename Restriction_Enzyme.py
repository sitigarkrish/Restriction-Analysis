# Restriction site analysis of a Dna Sequence Using BioPython using single enzyme

from Bio import SeqIO
from Bio.Restriction import EcoRI,RestrictionBatch,Analysis

record = SeqIO.read("Brca_cds.fasta", "fasta") # using BRCA1 Coding Sequence
sequence = record.seq

rb = RestrictionBatch([EcoRI])

analysis = Analysis(rb, sequence)

print("Coding Sequence of BRCA1 : ",sequence)
print("Length Of Coding Sequence : ",len(sequence))
print("Rcognition position for EcoRI : ",analysis.full())
print("Recognition Site : " ,EcoRI.site)
print("*"*30)

# Restriction site analysis of a Dna Sequence Using BioPython using Multiple enzyme

from Bio.Restriction import EcoRI,BamHI,HindIII,RestrictionBatch,Analysis
rb_multi = RestrictionBatch([EcoRI,BamHI,HindIII])
analysis_multi = Analysis(rb_multi,sequence)

result=analysis_multi.full()

for enzyme,positions in result.items():
    print("\nEnzyme : ",enzyme)
    print("Recognition Site : ",enzyme.site)
    if positions:
        print("Recognition Postitions : ",positions)

    else:
        print("no restriction site found")
