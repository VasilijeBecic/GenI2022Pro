# -*- coding: utf-8 -*-
"""
Created on Sat May 21 15:09:48 2022

@author: vbeci
"""

"""Traditional implementation directly from the classes"""
class BurrowsWheelerFMIndex:
    
    def rotations(self, t):
        """ Return list of rotations of input string t """
        tt = t * 2
        return [tt[i:i+len(t)] for i in range(0, len(t)) ]


    def bwm(self, t):
        """ Return lexicographically sorted list of t’s rotations """
        return sorted(self.rotations(t))


    def bwtViaBwm(self, t):
        """ Given T, returns BWT(T) by creating BWM """
        return ''.join(map(lambda x: x[-1], self.bwm(t)))
    
    
    def suffixArray(self, s):
        """ Given T return suffix array SA(T).  We use Python's sorted function here for simplicity, but we can do better."""
        satups = sorted([(s[i:], i) for i in range(len(s))])
        # Extract and return just the offsets
        # print(satups)
        return map(lambda x: x[1], satups)


    def bwtViaSa(self, t):
        """ Given T, returns BWT(T) by way of the suffix array. """
        bw = []
        sa = self.suffixArray(t)
        for si in sa:
            if si == 0: bw.append('$')
            else: bw.append(t[si-1])
        return ''.join(bw) # return string-ized version of list bw
    
    
bwFM = BurrowsWheelerFMIndex()

t = 'BANANA$'
print(bwFM.rotations(t))