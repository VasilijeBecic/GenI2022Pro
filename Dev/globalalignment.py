# -*- coding: utf-8 -*-
"""
Created on Sun May 22 20:34:09 2022

@author: vbeci
"""

import numpy


'''
    def scoring_matrix(self, a, b):
        if a == b: return 1
        if a == '_' or b == '_' : return -7
        maxb, minb = max(a, b), min(a, b)
        if minb == 'A' and maxb == 'G': return -1
        if minb == 'C' and maxb == 'T': return -1
        return -2
'''

class Aligner:
    match = 0
    mismatch = 0
    gap = 0
    
    def __init__(self, match, mismatch, gap):
        ''' Match, mismatch, gap '''
        self.match = match
        self.mismatch = mismatch
        self.gap = gap
        

    def scoring_matrix(self, a, b):
        if a == b: return self.match
        if a == '_' or b == '_' : return self.gap
        return self.mismatch


    def global_alignment(self, x, y):
        ''' Accepts x, y '''
        D = numpy.zeros((len(x) + 1, len(y) + 1), dtype=int)
        
        for i in range(1, len(x) + 1):
            D[i,0] = D[i-1,0] + self.scoring_matrix(x[i-1], '_')  
        for j in range(1, len(y)+1):
            D[0,j] = D[0,j-1] + self.scoring_matrix('_', y[j-1])
        
        for i in range(1, len(x) + 1):
            for j in range(1, len(y) + 1):
                D[i,j] = max(D[i-1,j]   + self.scoring_matrix(x[i-1], '_'),
                             D[i,j-1]   + self.scoring_matrix('_', y[j-1]), 
                             D[i-1,j-1] + self.scoring_matrix(x[i-1], y[j-1]))
                
        # function returns table and global alignment score
        #alignment score is in cell (n,m) of the matrix
        return D, D[len(x),len(y)] 


    def traceback(self, x, y, V):
        ''' Accepts x, y, V '''
        # initializing starting position cell(n,m)
        i=len(x)
        j=len(y)
        
        # initializing strings we use to represent alignments in x, y, edit transcript and global alignment
        ax, ay, am, tr = '', '', '', ''
        
        # exit condition is when we reach cell (0,0)
        while i > 0 or j > 0:
            
            # calculating diagonal, horizontal and vertical scores for current cell
            d, v, h = -100, -100, -100
            
            if i > 0 and j > 0:
                delta = 1 if x[i-1] == y[j-1] else 0
                d = V[i-1,j-1] + self.scoring_matrix(x[i-1], y[j-1])  # diagonal movement   
            if i > 0: v = V[i-1,j] + self.scoring_matrix(x[i-1], '_')  # vertical movement
            if j > 0: h = V[i,j-1] + self.scoring_matrix('_', y[j-1])  # horizontal movement
                
            # backtracing to next (previous) cell
            if d >= v and d >= h:
                ax += x[i-1]
                ay += y[j-1]
                if delta == 1:
                    tr += 'M'
                    am += '|'
                else:
                    tr += 'R'
                    am += ' '
                i -= 1
                j -= 1
            elif v >= h:
                ax += x[i-1]
                ay += '_'
                tr += 'D'
                am += ' '
                i -= 1
            else:
                ay += y[j-1]
                ax += '_'
                tr += 'I'
                am += ' '
                j -= 1
                
        alignment='\n'.join([ax[::-1], am[::-1], ay[::-1]])
        return alignment, tr[::-1]


'''
# TESTS

x = 'TACGTCAGC'
y = 'TATGTCATGC'


match = 2 # 0 1 2
mismatch = -3 # -3 -2
gap = -7 # -5, -7

aligner = Aligner(match, mismatch, gap)

D, alignmentScore = aligner.global_alignment(x, y)
alignment, transcript = aligner.traceback(x, y, D)

print(alignment)
print(transcript)
print(D)
print(alignmentScore)
print(transcript)

# Expected for match = 2, mismatch = -3, gap = -7: 

TACGTCA_GC
|| |||| ||
TATGTCATGC
MMRMMMMIMM
[[  0  -7 -14 -21 -28 -35 -42 -49 -56 -63 -70]
 [ -7   2  -5 -12 -19 -26 -33 -40 -47 -54 -61]
 [-14  -5   4  -3 -10 -17 -24 -31 -38 -45 -52]
 [-21 -12  -3   1  -6 -13 -15 -22 -29 -36 -43]
 [-28 -19 -10  -6   3  -4 -11 -18 -25 -27 -34]
 [-35 -26 -17  -8  -4   5  -2  -9 -16 -23 -30]
 [-42 -33 -24 -15 -11  -2   7   0  -7 -14 -21]
 [-49 -40 -31 -22 -18  -9   0   9   2  -5 -12]
 [-56 -47 -38 -29 -20 -16  -7   2   6   4  -3]
 [-63 -54 -45 -36 -27 -23 -14  -5  -1   3   6]]
6
MMRMMMMIMM
'''