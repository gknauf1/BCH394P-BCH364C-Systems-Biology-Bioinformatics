# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 11:09:54 2020

@author: Gregory Knauf
"""
#import os
#os.chdir('C:\\Users\\Gregory Knauf\\Documents\\Nb seq transfer')

import pandas as pd

#Combined CDRs

sequence = "TU_combined_CDRs.txt"
total_length = 0
aa = {}

file = open(sequence, "r")
for raw_line in file:
    line = raw_line.rstrip("\r\n")  #removes return keys from text lines
    line = line.replace('"', '')    #removes quote marks
    length = len(line)              #function to find length of string
    for amino in line:
        if amino in aa:             #new has_key function replacement
            aa[amino] += 1
        else:
            aa[amino] = 1
    total_length += length
file.close()

for i in aa.keys():
    fraction = 100.0 * aa[i] / total_length
    print("The nucleotide {0} occurs {1} times, or {2}%".format(i, aa[i], fraction))    

#make frequency dict

for k, v in aa.items():
    aa[k] = 100.0 * aa[k] / total_length

#create combined CDR freq datateble

aa_freq_table = pd.DataFrame.from_dict(aa, orient = 'index')
aa_freq_table = aa_freq_table.rename(columns = {0:'Combined'})

#CDR1

sequence = "TU_CDR1.txt"
total_length1 = 0
aa1 = {}

file = open(sequence, "r")
for raw_line in file:
    line = raw_line.rstrip("\r\n")  #removes return keys from text lines
    line = line.replace('"', '')    #removes quote marks
    length = len(line)              #function to find length of string
    for amino in line:
        if amino in aa1:             #new has_key function replacement
            aa1[amino] += 1
        else:
            aa1[amino] = 1
    total_length1 += length
file.close()

for i in aa1.keys():
    fraction = 100.0 * aa1[i] / total_length1
    print("The nucleotide {0} occurs {1} times, or {2}%".format(i, aa1[i], fraction))    

#make frequency dict

for k, v in aa1.items():
    aa1[k] = 100.0 * aa1[k] / total_length1

#create CDR1 freq datateble

aa1_freq_table = pd.DataFrame.from_dict(aa1, orient = 'index')
aa1_freq_table = aa1_freq_table.rename(columns = {0:'CDR1'})


#CDR2

sequence = "TU_CDR2.txt"
total_length2 = 0
aa2 = {}

file = open(sequence, "r")
for raw_line in file:
    line = raw_line.rstrip("\r\n")  #removes return keys from text lines
    line = line.replace('"', '')    #removes quote marks
    length = len(line)              #function to find length of string
    for amino in line:
        if amino in aa2:             #new has_key function replacement
            aa2[amino] += 1
        else:
            aa2[amino] = 1
    total_length2 += length
file.close()

for i in aa2.keys():
    fraction = 100.0 * aa2[i] / total_length2
    print("The nucleotide {0} occurs {1} times, or {2}%".format(i, aa2[i], fraction))    

#make frequency dict

for k, v in aa2.items():
    aa2[k] = 100.0 * aa2[k] / total_length2

#create combined CDR freq datateble

aa2_freq_table = pd.DataFrame.from_dict(aa2, orient = 'index')
aa2_freq_table = aa2_freq_table.rename(columns = {0:'CDR2'})



#CDR3

sequence = "TU_CDR3.txt"
total_length3 = 0
aa3 = {}

file = open(sequence, "r")
for raw_line in file:
    line = raw_line.rstrip("\r\n")  #removes return keys from text lines
    line = line.replace('"', '')    #removes quote marks
    length = len(line)              #function to find length of string
    for amino in line:
        if amino in aa3:             #new has_key function replacement
            aa3[amino] += 1
        else:
            aa3[amino] = 1
    total_length3 += length
file.close()

for i in aa3.keys():
    fraction = 100.0 * aa3[i] / total_length3
    print("The nucleotide {0} occurs {1} times, or {2}%".format(i, aa3[i], fraction))    

#make frequency dict

for k, v in aa3.items():
    aa3[k] = 100.0 * aa3[k] / total_length3

#create combined CDR freq datateble

aa3_freq_table = pd.DataFrame.from_dict(aa3, orient = 'index')
aa3_freq_table = aa3_freq_table.rename(columns = {0:'CDR3'})


#combine frequency tables

TU_Amino_Acid_Freq = pd.concat([aa_freq_table, aa1_freq_table, aa2_freq_table, aa3_freq_table], axis = 1)

TU_Amino_Acid_Freq.to_csv('TU_Amino_Acid_Freq.csv')
