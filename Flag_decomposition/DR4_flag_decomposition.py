#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 11:28:10 2023

@author: scottjoffre


Flag decomposition for 4FGL-DR4, DR3
"""

#GOAL: Want to easily input a file, and then output the types of flags for each of the sources in a new column (according to DR4)

#Notes
#1 & number -> if even gives 0, if odd gives 1

# x<<y  means x*2**y
#so i<<= 1 means i = i + 1*2**1
#i += 1 -> i = i + 1

#Imports
import pandas as pd
import numpy as np


#input .csv file of interest
df = pd.read_csv('4fgl-DR4.csv')

flags = df['Flags']
#Add 2 new columns 
df['flag_decomp_number'] = '' #give number of the flag for easy counting/comparision to Table4. 
df['Flag_Decomposition'] = '' #will literally write out flag for easy readability



#Thanks to nullptr for this function as provided here: https://stackoverflow.com/questions/30226094/how-do-i-decompose-a-number-into-powers-of-2
def decomp(x):
    powers = []
    i = 1
    while i <= x:
        if i & x:
            powers.append(i)
        i <<= 1 # really i + i*2**1 = i + i*2
    return powers



for j in range(len(flags)):
    mine = flags[j]
    
    #Decompose into the different powers of 2 which added together give the Flag number reported in the DR4
    mydecomp = decomp(mine)
    
    #iterate through each individual value and convert to the actual flag
    n_flag_list = []
    #get which power of 2 this value corresponds to
    power_val = np.log10(mydecomp)/np.log10(2)
    #Since it is 2^flag-1 then the flag = value + 1
    n_flag = power_val + 1 #this adds 1 to each value in the list
    n_flag_list.append(n_flag)


    #Now add the flags to the ['Flag_Decomposition] list 
    # As well as add the number of the flag itself as reported. See Table 4 from 4FGL-DR4(https://arxiv.org/pdf/2307.12546.pdf) or DR3 (https://arxiv.org/pdf/2201.11184.pdf)  
    #Additional notes can be found here: https://heasarc.gsfc.nasa.gov/w3browse/fermi/fermilpsc.html
    if 1 in n_flag_list[0]:
        df.loc[j,'Flag_Decomposition'] += 'TS<25 with other model,' # note: was at TS>35
        df.loc[j,'flag_decomp_number'] += '1,'

    if 2 in n_flag_list[0]:
        df.loc[j,'Flag_Decomposition'] += 'Moved beyond 95% error ellipse,'
        df.loc[j,'flag_decomp_number'] += '2,'
    
    if 3 in n_flag_list[0]:
        df.loc[j,'Flag_Decomposition'] += 'Flux changed with other model or analysis,' # flux or energy flux change by >3sigma
        df.loc[j,'flag_decomp_number'] += '3,'

    if 4 in n_flag_list[0]:
        df.loc[j,'Flag_Decomposition'] += 'Src/bkg < 10%,'
        df.loc[j,'flag_decomp_number'] += '4,'

        
        
    if 5 in n_flag_list[0]:
        df.loc[j,'Flag_Decomposition'] += 'Confused,'
        df.loc[j,'flag_decomp_number'] += '5,'

        
    if 6 in n_flag_list[0]:
        df.loc[j,'Flag_Decomposition'] += 'Interstellar gas clump,'
        df.loc[j,'flag_decomp_number'] += '6,'

        
    if 9 in n_flag_list[0]:
        df.loc[j,'Flag_Decomposition'] += 'Localization flag from pointlike'
        df.loc[j,'flag_decomp_number'] += '9,'

        
    if 10 in n_flag_list[0]:
        df.loc[j,'Flag_Decomposition'] += 'Bad spectral fit quality'
        df.loc[j,'flag_decomp_number'] += '10,'

        
    if 12 in n_flag_list[0]:
        df.loc[j,'Flag_Decomposition'] += 'Highly curved spectrum' # LogParabola beta fixed to 1 or PLEC_Index fixed to 0
        df.loc[j,'flag_decomp_number'] += '12,'

        
    if 13 in n_flag_list[0]:
        df.loc[j,'Flag_Decomposition'] += 'TS<25 at 12yr'
        df.loc[j,'flag_decomp_number'] += '13,'

        
    if 14 in n_flag_list[0]:
        df.loc[j,'Flag_Decomposition'] += 'Soft Galactic Unassociated (SGU)'
        df.loc[j,'flag_decomp_number'] += '14,'

#Save the newly updated .csv file
df.to_csv('4FGLdr4_flags_decomp.csv') #example with DR4 input


#Break up (and possible save as .csv files) for all sources containing one of the individual flags listed. 
df1 = df[df['flag_decomp_number'].str.contains('1,')]
print(len(df1))

df2 = df[df['flag_decomp_number'].str.contains('2,')]
print(len(df2))

df3 = df[df['flag_decomp_number'].str.contains('3,')]
print(len(df3))

df4 = df[df['flag_decomp_number'].str.contains('4,')]
print(len(df4))

df5 = df[df['flag_decomp_number'].str.contains('5,')]
print(len(df5))

df6 = df[df['flag_decomp_number'].str.contains('6,')]
print(len(df6))

df9 = df[df['flag_decomp_number'].str.contains('9,')]
print(len(df9))

df10 = df[df['flag_decomp_number'].str.contains('10,')]
print(len(df10))

df12 = df[df['flag_decomp_number'].str.contains('12,')]
print(len(df12))

df13 = df[df['flag_decomp_number'].str.contains('13,')]
print(len(df13))

df14 = df[df['flag_decomp_number'].str.contains('14,')]
print(len(df14))
