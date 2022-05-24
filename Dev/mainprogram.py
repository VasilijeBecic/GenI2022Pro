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
from datetime import datetime


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

def print_curr_datetime(message):
    now = datetime.now()
    dtString = now.strftime("%d/%m/%Y %H:%M:%S")
    print(str(message) + " = ", dtString)


def store_to_csv(data, exportFileName):
    df = pd.DataFrame(data, columns =['Position', 'alignmentScore', 'transcript'])
    df.to_csv(exportFileName)


# Main program
def main():
    # To pass argv in Spyder -> Run>Configuration per file>Command line options. i.e: 2 3 a
    script, fasta, fastq, match, mismatch, gap, seedlength, margin = sys.argv
    
    match = int(match)
    mismatch= int(mismatch)
    gap = int(gap)
    seedlength = int(seedlength)
    margin = int(margin)
    
    print_curr_datetime("ANALYSIS START")
    
    
    # fasta += '.fasta'
    # fastq += '.fastq'
    print('Read files')
    referenceGenome = readFASTA(fasta)
    reads = readFASTQ(fastq)
    print('Files read')
    print()

        
     
    print('Begin fm index init.')
    fmIndex = FMIndex(referenceGenome)
    print('End fm index init.')
    print()
    
    print('Begin aligner init.')
    aligner = Aligner(match, mismatch, gap)
    print('End aligner init.')
    print()
    
    i = 1
    print("Start reads")
    allresults = []
    for read in reads:
        print_curr_datetime("Start datetime for read " + str(i))
        
        results = seed_and_extend(referenceGenome, read, seedlength, margin, aligner, fmIndex)
        results.append(i)
        allresults.append(results)
        
        print_curr_datetime("End datetime for read " + str(i))
        
        
        print_curr_datetime("Start datetime for reversed read " + str(i))
        
        read = reverse_complement(read)
        results = seed_and_extend(referenceGenome, read, seedlength, margin, aligner, fmIndex)
        results.append(i)
        allresults.append(results)
        
        print_curr_datetime("End datetime for reversed read " + str(i))
        
        i += 1
    
    
    print("Finish reads and all")
    store_to_csv(allresults, 'result.csv')
    
    print_curr_datetime("ANALYSIS FINISHED")
    
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

