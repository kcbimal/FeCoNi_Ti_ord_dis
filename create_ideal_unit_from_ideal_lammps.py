# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 14:07:13 2023

@author: biknb
"""


import os, sys
import pandas as pd

sys = 'CoTi'
alat = [2.91, 2.92, 2.93, 2.94, 2.95, 2.96, 2.97, 2.98, 2.99, 3.00, 3.01, 3.02, 3.03, 3.04, 3.05, 3.06, 3.07, 3.08, 3.09, 3.10] #you can input here a list or single volume
temp = [10, 50, 90] #CoTi
# temp = [250, 290, 330]  #NiTi   

for x in alat:
    for y in temp:
        df = pd.read_csv('C:\\Users\\biknb\\Downloads\\bkc\\Phonons_0.5\\ord\\'+sys+'\\LAMMPS\\'+sys+'_'+str(x)+'_'+str(y)+'K' +'\\ideal_'+sys+'_B2_'+str(y)+'K' +'.txt', sep=' ', skiprows=9, skipinitialspace=True, header=None)
        
        df.to_csv('C:\\Users\\biknb\\Downloads\\bkc\\Phonons_0.5\\ord\\'+sys+'\\LAMMPS\\'+sys+'_'+str(x)+'_'+str(y)+'K' +'\\ideal_'+sys+'_B2_'+str(y)+'K' +'_unit.txt', header=None, index=None, sep=' ')#, float_format='%.4f')

