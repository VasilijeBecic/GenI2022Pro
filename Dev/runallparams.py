# -*- coding: utf-8 -*-
"""
Created on Sun May 29 17:10:24 2022

@author: vbeci
"""
import sys


sys.runfile('mainprogram.py', args='example_human_reference.fasta example_human_Illumina.pe_1.fastq 1 -3 -7 10 2')


'''
for match in range(0,3):
        for mismatch in range(-3,-1):
            for gap in range(-7, -4):
                if gap == -6:
                    continue
                '''