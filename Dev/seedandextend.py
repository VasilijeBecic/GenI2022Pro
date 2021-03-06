# -*- coding: utf-8 -*-
"""
"""
from fm import FMIndex
from globalalignment import Aligner


def reverse_complement(read):
    ''' read '''
    reversedRead = []
    i = len(read) - 1
    while i >= 0:
        if read[i] == 'A':
            reversedRead.append('T')
        elif read[i] == 'T':
            reversedRead.append('A')
        elif read[i] == 'C':
            reversedRead.append('G')
        elif read[i] == 'G':
            reversedRead.append('C')
        i -= 1
    return reversedRead


def seed_and_extend(referenceGenome, readId, isReverseComplement, read, seedLength, margin, aligner, fmIndex):
    ''' referenceGenome, readId, read, seedLength, margin, aligner, fmIndex'''
    results = []
    seed = read[0:seedLength]
    
    seedPositions = fmIndex.query(seed)
    
    i = 0
    for position in seedPositions:
        # print(position)
        start = position + seedLength
        end = start + len(read) - seedLength + margin # end is not counted in slicing
        if end >= len(referenceGenome):
            continue
            
        alignedRef = referenceGenome[start:end]
        alignedRead = read[seedLength:]
        
        D, alignmentScore = aligner.global_alignment(alignedRef, alignedRead)
        alignment, transcript = aligner.traceback(alignedRef, alignedRead, D)
        
        # print('alignment completed for ' + str(i))
        
        i += 1
        
        results.append((readId, isReverseComplement, start - seedLength, alignmentScore, transcript))
    
    # print('seed_and_extend completed')
    results.sort(key=lambda x: x[3], reverse=True)
    return results
    



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
print(reversedRead)

# Expected for referenceGenome = 'AAGAAGTCAGGGAGCAAGCAGAGTCAGGGAGCAAGCCACCAC'
# read = 'AGTCAGGGAGCAAGC'
# found seed positions: 4, 21
# [(4, -9, 'MMMMMDD'), (21, -9, 'MMMMDMD')]
# reversed read 
# ['G', 'C', 'T', 'T', 'G', 'C', 'T', 'C', 'C', 'C', 'T', 'G', 'A', 'C', 'T']
'''