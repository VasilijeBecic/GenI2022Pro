# -*- coding: utf-8 -*-
"""
Created on Tue May 24 20:41:40 2022

@author: vbeci
"""
from seedandextend import seed_and_extend
from seedandextend import reverse_complement
from globalalignment import Aligner
from fm import FMIndex
from Bio import SeqIO
from Bio.Seq import Seq

### Main program
# Main program
def main():
    # TESTS
    referenceGenome = 'AAGAAGTCAGGGAGCAAGCAGAGTCAGGGAGCAAGCCACCAC'
    read = 'AGTCAGGGAGCAAGC'
    reversedRead = reverse_complement(read)
    seek_length = 10
    margin = 2
    
    '''
    for match in range(0,3):
        for mismatch in range(-3,-1):
            for gap in range(-7, -4):
                if gap == -5:
                    continue
                aligner = Aligner(match, mismatch, gap)
    '''
    
    
    aligner = Aligner(1, -3, -7)
    fmIndex = FMIndex(referenceGenome)
    results = seed_and_extend(referenceGenome, read, seek_length, margin, aligner, fmIndex)
    print(results)

if __name__ == "__main__":
    main()

