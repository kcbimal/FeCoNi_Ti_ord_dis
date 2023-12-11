# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 11:13:57 2022

@author: bkc/Munoz
"""
# run phonopy first to get pdos (using terminal and use cmd: phonopy -p -s mesh.conf)

import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
from scipy.integrate import quad
import scipy.integrate as integrate

# s = 's=0.1'        # p/dos_sigma value = 0.1, 0.9
syst = 'ord'        # ord, dis
sys = 'FeTi'
t = 700
kb = 8.617333262e-5 # eV/K

vols = [2.85, 2.86, 2.87, 2.88, 2.89, 2.90, 2.91, 2.92, 2.93, 2.94, 2.95, 2.96, 2.97, 2.98, 2.99, 3.00] #you can input here a list or s aingle volume
temp = [700]  

for x in vols:
    for y in temp:
                                                                                                                                                  
        #loading projected DOS
        # df = pd.read_csv("C:\\Users\\biknb\\Downloads\\bkc\\Phonons_0.5\\Phonopy\\"+root+"_"+str(x)+"-"+syst+"\\projected_dos_"+s+".dat", engine='python', sep="        ", skiprows=1, names = ["f", "fe", "v"])
        df = pd.read_csv('C:\\Users\\biknb\\Downloads\\bkc\\Phonons_0.5\\'+sys+'\\'+syst+'\\Phonopy\\'+str(y)+'\\'+sys+'_'+str(x)+'_'+str(y)+'K\\projected_dos.dat', engine='python', sep='        ', skiprows=1, names = ["f", "fe", "v"])
        #getting total DOS usng Fe + V, converting THz to eV for freq.
        dos = pd.Series(df.fe.add(df.v).values, index=df.f.mul(0.0041).values)
        # print(dos)        #(1/0.0245998) ord
        # dos.plot()        #(1/0.024599449497920003) dis
        w = dos.mul(0.05*0.0041)     #delta_f = 0.05   &    1THz = 0.0041ev
        # print("DOS_sum_"+str(y)+"_"+str(x)+"_"+syst+" = ", w.sum())
            
        #normalized DOS using 1/w.sum() as normalization factor
        dos_n = pd.Series(df.fe.add(df.v).mul(1/w.sum()).values, index=df.f.mul(0.0041).values)    
        
        # # Area under the curve
        # rise = dos.iloc[len(dos)-1] - dos.iloc[0]
        # run = dos.index[len(dos)-1] - dos.index[0]
        # integral1 = rise*run
        # print("Area under the curve = ", integral1)
        
    
        
        #vibrational entropy
        def BE(ei, T):              #Bose-Einstein distribution
            return (1/(math.exp((ei)/(kb*T)) - 1))
        
        n = pd.Series([BE(e, t)for e in dos.index[1:]], index=dos.index[1:])      #n(w); 
        # print(n) 
        # n.plot()                                                    #BE distribution
        # n.to_csv("C:\\Users\\biknb\\Downloads\\n.dat")
        
        N = n.add(1).mul(np.log(n.add(1))).subtract(n.mul(np.log(n)))
        # print(N)
        # N.plot()
        # N.to_csv("C:\\Users\\biknb\\Downloads\\N.dat")
        
        #integrand for vibrational entropy    g(w)*n(w)
        I = dos_n.mul(N).dropna()    #w = dos
        # print(I)
        I.plot()
        # I.to_csv("C:\\Users\\biknb\\Downloads\\I.dat")
        
        z = 3*I.mul(0.05*0.0041)
        z = z.sum()
        print("S_vib_"+str(y)+"K_"+syst+"      "+str(x)+"      ", z)