#!/usr/bin/env python3

import re

def Levenshtein_dist(s, t):
    # Hjelmqvist algorithm from wikipedia,
    # translated to python
     
    # degenerate cases
    if s == t:
        return 0
    if len(s) == 0:
        return len(t)
    if len(t) == 0:
        return len(s)

    # create two work vectors of integer distances
    v0 = []
    v1 = []

    # initialize v0 (the previous row of distances)
    #   and v1 (the current row of distances)
    # the initial state of v0 is the edit distance for an empty s
    # the distance is just the number of characters to delete from t
    # we initialize v1 to 0.
    for i in range(len(t)+1):
        v0.append(i)
        v1.append(0)

    # now do the algorithm
    for i in range(len(s)):
        # calculate v1 (current row distances) from the previous row v0

        # first element of v1 is A[i+1][0]
        # edit distance is delete (i+1) chars from s to match empty t
        v1[0] = i + 1

        # now fill in the rest of the row
        for j in range(len(t)):
            if s[i] == t[j]:
                cost = 0
            else:
                cost = 1
            v1[j + 1] = min(v1[j] + 1,
                            v0[j + 1] + 1,
                            v0[j] + cost)

        # copy v1 (current row) to v0 (previous row) for next iteration
        v0 = v1.copy()

    return v1[len(t)]



# read in peptides
with open('Combined_CDR_AbTU.txt', 'r') as f:
    lines = f.readlines()


# make list of peptides (without newline)
peptides = []
imax = 10000 # maximum number of peptides to analyze
i = 0
for line in lines:
    p = line.strip() # remove newline
    p = line.replace('"', '')
    p = p.split('*')[0] # remove everything after first stop codon
    
    if len(p) > 22: # exclude peptides of length 3 or less
        peptides.append(p)
    i += 1
    if i >= imax:
        break

i = 0
j = 0 
   
z=open("Combined_CDR_AbTU_dis10_2.csv","w")    
for i in range(len(peptides)):
    for j in range(len(peptides)):
        p1 = peptides[i]
        p2 = peptides[j]
        p1 = p1.rstrip("\r\n")
        p2 = p2.rstrip("\r\n")
        p1 = p1.replace('"', '')
        p2 = p2.replace('"', '')
        d_full = Levenshtein_dist(p1, p2)
       
        
        pr = (str(p1), str(p2), str(d_full))
        z.write(','.join(pr))
        z.write("\n")
        
f.close()
        
z.close()        


        
#cat distances.txt | awk '{if ($8<.3){print $0}}'

        

'C:\\Users\\Gregory Knauf\\Documents\\Python Scripts\\Peptides\\Scripts\\full.csv','w'