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
import sys
import pandas as pd


def readFASTA(fasta_file):
  arr = []
  for seq_record in SeqIO.parse(fasta_file, "fasta"):
    arr.append(seq_record.seq)
  return arr[0]


def readFASTQ(fastq_file):
  arr = []
  for seq_record in SeqIO.parse(fastq_file, "fastq"):
    arr.append(seq_record.seq)
  return arr


# fasta = 'reference'
# fastq = 'reads'

# print(readFASTA(fasta + '.fasta'))
# print(readFASTQ(fastq + '.fastq')[0])


def store_to_csv(data, exportFileName):
    df = pd.DataFrame(data, columns =['Position', 'alignmentScore', 'transcript'])
    df.to_csv(exportFileName)


# Main program
def main():
    # To pass argv in Spyder -> Run>Configuration per file>Command line options. i.e: 2 3 a
    script, fasta, fastq, match, mismatch, gap, seedlength, margin = sys.argv
    
    # fasta += '.fasta'
    # fastq += '.fastq'
    print('Read files')
    referenceGenome = readFASTA(fasta)
    reads = readFASTQ(fastq)
    print('Files read')
    print()
    
    # for read in reads:
        #
     
    print('Begin fm index init.')
    fmIndex = FMIndex(referenceGenome)
    print('End fm index init.')
    print()
    
    print('Begin aligner init.')
    aligner = Aligner(match, mismatch, gap)
    print('End aligner init.')
    print()
    
    i = 1
    print('Begin for read ' + str(i))
    results = aligner.seed_and_extend(referenceGenome, reads[0].seq, seedlength, margin, aligner, fmIndex)
    print('Finish for read ' + str(i))
    
    print()
    store_to_csv(results, 'result.csv')
    
    '''
    # TESTS
    referenceGenome = 'AAGAAGTCAGGGAGCAAGCAGAGTCAGGGAGCAAGCCACCAC'
    read = 'AGTCAGGGAGCAAGC'
    reversedRead = reverse_complement(read)
    seek_length = 10
    margin = 2

    
    
    aligner = Aligner(1, -3, -7)
    fmIndex = FMIndex(referenceGenome)
    results = seed_and_extend(referenceGenome, read, seek_length, margin, aligner, fmIndex)
    print(results)
    '''

if __name__ == "__main__":
    main()


    '''
    for match in range(0,3):
        for mismatch in range(-3,-1):
            for gap in range(-7, -4):
                if gap == -5:
                    continue
                aligner = Aligner(match, mismatch, gap)
    '''

