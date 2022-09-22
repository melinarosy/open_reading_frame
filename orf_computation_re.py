from revc_computation import rev_complement
from rna_computation import dna2rna
from subs_computation import solve_subs
from prot_computation import rna2protein
import re

def dna2proteins_set(dna):
    orf_pattern = r'(AUG([ACGU]{3})*?)(UAA|UAG|UGA)'
    revc = rev_complement(dna)
    frames = [dna, revc] # double-stranded DNA. that is why we have both DNA and Reverse DNA
    proteins = set()
    for frame in frames:
        rna = dna2rna(frame)
        for pos in range(len(rna)):
            sub_rna = rna[pos:]
            match = re.match(orf_pattern, sub_rna)
            if match != None:
                prot = rna2protein(match.group(0))
                proteins.add(prot)
    return proteins
