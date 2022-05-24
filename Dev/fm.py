# -*- coding: utf-8 -*-
"""
Created on Sun May 22 17:39:34 2022

@author: vbeci
"""

from bwt import BWT

class FMIndex:
    
    bwt = {}
    tally = dict() # Key = character, Value = occurences of character (len of BWT(T))
    c = dict()  # Key = character, Value = (value for first row occurence)
    
    def __init__(self, text):
        text = text + '$'
        self.clear_parameters()
        self.bwt = BWT() # We only care about L and SA
        self.bwt.transform(text)
        self.init_tally() # Initialize tally
        self.init_c() # Initialize c
        
    
    def clear_parameters(self):
        self.bwt = {}
        self.tally = dict()
        self.c = dict()
        
    
    def init_tally(self):
        empty = []
        occurences = dict()
        for i in range(len(self.bwt.transformedText)):
            empty.append(0)
            occurences[self.bwt.transformedText[i]] = 0
            i += 1
            
        for char in self.bwt.firstCol.keys():
            self.tally[char] = empty.copy()
            
        for i in range(len(self.bwt.transformedText)):
            character = self.bwt.transformedText[i]
            occurences[character] += 1
            
            for char in self.tally.keys():
                self.tally[char][i] = occurences[char]
                

    def init_c(self):
        ''' Initializes c dictionary, where key = character 
        and value = tuple (index, value for first row occurence)'''
        i = 0
        for char in self.bwt.firstCol.keys():
            self.c[char] = (i, self.bwt.firstCol[char][0])
            i += 1
            
    
    def query(self, pattern):
        ''' Returns list of found positions for pattern in the text. 
        Rows are counted from 0 to n, not from 1.'''
        positions = []
        start = 1
        end = len(self.bwt.transformedText) - 1
        k = len(pattern) - 1
        while k > -1:
            ch = pattern[k]
            # start = C[ch], end = C[ch + 1] with break if ch + 1 > all keys
            if k == len(pattern) -1:
                start = self.c[ch][1]
                if self.c[ch][0] == len(self.c) -1:
                    break
                # Find C[ch + 1]
                index = self.c[ch][0] + 1
                newKey = ''
                for key in self.c:
                    if self.c[key][0] == index:
                        newKey = key
                        break
                end = self.c[newKey][1] - 1
            else:
                start = self.c[ch][1] + self.tally[ch][start - 1]
                end = self.c[ch][1] + self.tally[ch][end] - 1
                
            k -= 1
        for i in range(start, end + 1):
            positions.append(self.bwt.saIndexes[i])
        return sorted(positions)
        
    

# TESTS
# text = 'abaaba' # Example - Expected
# text = 'GATGCGAGAGATG'
# text = 'BANANA'
# fmIndex = FMIndex(text)

# print('Tally:')
# print(fmIndex.tally)
# print()
# Expected {'$': [0, 0, 0, 0, 1, 1, 1], 'a': [1, 1, 1, 2, 2, 3, 4], 'b': [0, 1, 2, 2, 2, 2, 2]}

# print('C')
# print(fmIndex.c)
# print()
# Expected {'$': (0, 0), 'a': (1, 1), 'b': (2, 5)}

# pattern = 'aba'
# pattern = 'GAGA'
# pattern = 'ANA'
# print('Query for pattern: ' + pattern)
# print(fmIndex.query(pattern))
# print()
# Expected
# For text='GAGA' -> [5,7]
# For text='ANA' -> [1,3]




