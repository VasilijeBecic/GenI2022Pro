# -*- coding: utf-8 -*-
"""
Created on Sun May 22 20:34:09 2022

@author: vbeci
"""

import numpy


'''
    def scoringMatrix(self, a, b):
        if a == b: return self.match
        if a == '_' or b == '_' : return self.gap
        return self.mismatch
'''

class Aligner:
    match = 0
    mismatch = 0
    gap = 0
    
    def __init__(self, match, mismatch, gap):
        self.match = match
        self.mismatch = mismatch
        self.gap = gap
        

    def scoringMatrix(self, a, b):
        if a == b: return 1
        if a == '_' or b == '_' : return -7
        maxb, minb = max(a, b), min(a, b)
        if minb == 'A' and maxb == 'G': return -1
        if minb == 'C' and maxb == 'T': return -1
        return -2


    def globalAlignment(self, x, y, scoringMatrix):
        D = numpy.zeros((len(x) + 1, len(y) + 1), dtype=int)
        
        for i in range(1, len(x) + 1):
            D[i,0] = D[i-1,0] + scoringMatrix(x[i-1], '_')  
        for j in range(1, len(y)+1):
            D[0,j] = D[0,j-1] + scoringMatrix('_', y[j-1])
        
        for i in range(1, len(x) + 1):
            for j in range(1, len(y) + 1):
                D[i,j] = max(D[i-1,j]   + scoringMatrix(x[i-1], '_'),
                             D[i,j-1]   + scoringMatrix('_', y[j-1]), 
                             D[i-1,j-1] + scoringMatrix(x[i-1], y[j-1]))
                
        # function returns table and global alignment score
        #alignment score is in cell (n,m) of the matrix
        return D, D[len(x),len(y)] 


    def traceback(self, x, y, V, scoringMatrix):
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
                d = V[i-1,j-1] + scoringMatrix(x[i-1], y[j-1])  # diagonal movement   
            if i > 0: v = V[i-1,j] + scoringMatrix(x[i-1], '_')  # vertical movement
            if j > 0: h = V[i,j-1] + scoringMatrix('_', y[j-1])  # horizontal movement
                
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




# TESTS
x = 'TACGTCAGC'
y = 'TATGTCATGC'
match = 2 # 0 1 2
mismatch = -3 # -3 -2
gap = -7 # -5, -7

aligner = Aligner(match, mismatch, gap)

D, alignmentScore = aligner.globalAlignment(x, y, aligner.scoringMatrix)
alignment, transcript = aligner.traceback(x, y, D,  aligner.scoringMatrix)

print(alignment)
print(transcript)
print(D)
print(alignmentScore)
print(transcript)
# Expected: 
'''
TACGTCA_GC
|| |||| ||
TATGTCATGC
MMRMMMMIMM
[[  0  -7 -14 -21 -28 -35 -42 -49 -56 -63 -70]
 [ -7   1  -6 -13 -20 -27 -34 -41 -48 -55 -62]
 [-14  -6   2  -5 -12 -19 -26 -33 -40 -47 -54]
 [-21 -13  -5   1  -6 -13 -18 -25 -32 -39 -46]
 [-28 -20 -12  -6   2  -5 -12 -19 -26 -31 -38]
 [-35 -27 -19 -11  -5   3  -4 -11 -18 -25 -32]
 [-42 -34 -26 -18 -12  -4   4  -3 -10 -17 -24]
 [-49 -41 -33 -25 -19 -11  -3   5  -2  -9 -16]
 [-56 -48 -40 -32 -24 -18 -10  -2   3  -1  -8]
 [-63 -55 -47 -39 -31 -25 -17  -9  -3   1   0]]
0
MMRMMMMIMM
'''