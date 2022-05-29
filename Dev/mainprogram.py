# -*- coding: utf-8 -*-
"""
Created on Tue May 24 20:41:40 2022

@author: vbeci
"""
from seedandextend import seed_and_extend
from seedandextend import reverse_complement
from globalalignment import Aligner
from bwt import print_curr_datetime
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
    arr.append([seq_record.id, seq_record.seq])
  return arr


# fasta = 'reference'
# fastq = 'reads'

# print(readFASTA(fasta + '.fasta'))
# print(readFASTQ(fastq + '.fastq')[0])


def store_to_csv(data, exportFileName):
    df = pd.DataFrame(data, columns =['read_id', 'is_rev_comp', 'position','alignment_score','transcript'])
    df.to_csv(exportFileName)


# Main program
# def main():
def main(fasta, fastq, match, mismatch, gap, seedlength, margin):
    # To pass argv in Spyder -> Run>Configuration per file>Command line options. i.e: 2 3 a
    # script, fasta, fastq, match, mismatch, gap, seedlength, margin = sys.argv
    
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
    # fileName = 'results/result_read'

    step = 20
    bestResults = []
    for read in reads:
        if i % step == 0:
            print_curr_datetime("Start datetime for read " + str(i))
        
        results = seed_and_extend(referenceGenome, read[0], False, read[1], seedlength, margin, aligner, fmIndex)
        if (len(results) > 0):
            bestResults.append(results[0]) # For all
            
        # store_to_csv(results, fileName + str(i)+'.csv')
        
        if i % step == 0:
            print_curr_datetime("End datetime for read " + str(i))
        
        
        if i % step == 0:
            print_curr_datetime("Start datetime for reversed read " + str(i))
        
        reverseCompRead = reverse_complement(read[1])
        
        results = seed_and_extend(referenceGenome, read[0], True, reverseCompRead, seedlength, margin, aligner, fmIndex)
        if (len(results) > 0):
            bestResults.append(results[0]) # For all
            
        # store_to_csv(results, fileName + 'reversed' + str(i)+'.csv')
        
        if i % step == 0:
            print_curr_datetime("End datetime for reversed read " + str(i))
        
        i += 1
    
    
    store_to_csv(bestResults, 'results/results_ours/our_results' + str(match) + str(mismatch) + str(gap) + '.csv')
    
    print("Finish reads and all")
    
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
    fasta = 'example_human_reference.fasta'
    fastq = 'example_human_Illumina.pe_1.fastq'
    seed = 10
    margin = 2
    
    for match in range(1,3):
        for mismatch in range(-3,-1):
            for gap in range(-7, -4):
                if gap == -6:
                    continue
                main(fasta, fastq, match, mismatch, gap, seed, margin)
    # main()


    '''
    for match in range(0,3):
        for mismatch in range(-3,-1):
            for gap in range(-7, -4):
                if gap == -6:
                    continue
                aligner = Aligner(match, mismatch, gap)
    '''

