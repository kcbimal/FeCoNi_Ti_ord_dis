import numpy as np
from math import sqrt, isclose
import h5py
import os, sys
import pandas as pd
#make sure to check "POSCAR" and "band.conf" file before run
#           !!!!!!!!!!!!!!!!!!!!!! remove '' in line 12 !!!!!!!!!!!!!!!!!!


# Remember to modify POSCAR file to correct lattice parameter
syst = 'ord'
sys = 'FeTi'
root = 'latt_param_tmpK'
alat = latt_param * 10 # Lattice parameter times number of unit cells (so the size of the box)
ndisps = 1     #500    # Number of dispersions to generate (equal to number of yalm files)
n = 2000 # Number of atoms
init_no = 0    #1500 for 2000 steps to generate ndisps = 500


# in_path = 'C:\\Users\\biknb\\Downloads\\bkc\\Phonons_0.5\\'+syst+'\'+sys+'\\LAMMPS/'+sys+'_'+root+'/'
in_path = '/home/bkc/fall2023/'+syst+'/'+sys+'/LAMMPS/'+sys+'_'+root+'/'
in_filename = 'ideal_'+sys+'_B2_tmpK_unit.txt' # Modify pos file
ide_lat = np.genfromtxt(in_path+in_filename)
ide_lat = ide_lat[:n, 1:]

# in_filename = '432_ideal'
# ide_lat = np.genfromtxt(in_path+in_filename)
# ide_lat = ide_lat[:432, 1:]

# print(len(ide_lat))
    
a_val = alat / ((n / 2.)**(1./3.)) #Side length of each cube in bcc

ideal_distances = np.zeros((n, 3, n))
ideal_dist_sca = np.zeros((n, n))

for i in range(n):
    for j in range(n):
        if i != j:
            ideal_distances[i, :, j] = ide_lat[j, :] - ide_lat[i, :] #Ideal distance vector from i to j

            #Apply mic:
            for d in range(3):
                if ideal_distances[i, d, j] > (alat / 2):
                    ideal_distances[i, d, j] -= alat

                elif ideal_distances[i, d, j] <= (-alat / 2):
                    ideal_distances[i, d, j] += alat

            ideal_dist_sca[i, j] = np.linalg.norm(ideal_distances[i, :, j]) #Distance scalar

neighbors = [] #List of dictionaries holding a list of the IDs for all types of nn

#Distances to each type of neighbor:
first = (sqrt(3) * a_val) / 2
second = a_val
third = sqrt(2) * a_val
fourth = (sqrt(11) * a_val) / 2
fifth = sqrt(3) * a_val
nn_dist = [first, second, third, fourth, fifth]

for i in range(n):
    i_neigh = [[i]] 
    for prox in range(5):
        i_neigh.append([]) #Empty placeholder lists for neighbors of i
    
    #Add each list of nn to the placeholder dictionary, append the dictionary to the list:
    for j in range(n):
        if i != j:
            for prox in range(1,6):
                if isclose(ideal_dist_sca[i, j], nn_dist[prox-1], rel_tol=.001):
                    i_neigh[prox].append(j)

    neighbors.append(i_neigh)

fc_mat = np.zeros((n, n, 3, 3))
base_mat = np.zeros((6, 3, 3))
nn_matrix = .5 * a_val * np.array([[0.,0.,0.],
                                   [1.,1.,1.],
                                   [2.,0.,0.],
                                   [0.,2.,2.],
                                   [3.,1.,1.],
                                   [2.,2.,2.]])

fc_arr = np.array([[0, 0, 14, 14],
                   [1, 1, 2, 2],
                   [3, 4, 14, 14],
                   [5, 6, 14, 7],
                   [8, 9, 10, 11],
                   [12, 12, 13, 13]])

# For B2
# fc_in_path = 'C:\\Users\\biknb\\Downloads\\bkc\\Phonons_0.5\\'+syst+'\'+sys+'\\Force_constants/'+sys+'_'+root+'/'
fc_in_path = '/home/bkc/fall2023/'+syst+'/'+sys+'/Force_constants/'+sys+'_'+root+'/'
fc_filename = 'fc_'+sys+'_'+root+'.csv'

fc_per_ts_df = pd.read_csv(fc_in_path+fc_filename, skiprows=2, header=None)  #takes avg. from .csv file

for m in range(init_no, init_no+ndisps):
# for m in left:

    fc = [fc_per_ts_df.iloc[m, 0], fc_per_ts_df.iloc[m, 2], fc_per_ts_df.iloc[m, 3], fc_per_ts_df.iloc[m, 4], fc_per_ts_df.iloc[m, 5], fc_per_ts_df.iloc[m, 8], 
          fc_per_ts_df.iloc[m, 9], fc_per_ts_df.iloc[m, 10], fc_per_ts_df.iloc[m, 14], fc_per_ts_df.iloc[m, 15], fc_per_ts_df.iloc[m, 16], fc_per_ts_df.iloc[m, 17], 
          fc_per_ts_df.iloc[m, 18], fc_per_ts_df.iloc[m, 19], 0.0]
    
    fc_2 = [fc_per_ts_df.iloc[m, 1], fc_per_ts_df.iloc[m, 2], fc_per_ts_df.iloc[m, 3], fc_per_ts_df.iloc[m, 6], fc_per_ts_df.iloc[m, 7], fc_per_ts_df.iloc[m, 11], 
          fc_per_ts_df.iloc[m, 12], fc_per_ts_df.iloc[m, 13], fc_per_ts_df.iloc[m, 14], fc_per_ts_df.iloc[m, 15], fc_per_ts_df.iloc[m, 16], fc_per_ts_df.iloc[m, 17], 
          fc_per_ts_df.iloc[m, 20], fc_per_ts_df.iloc[m, 21], 0.0]
    
    #### 1ST SPECIES ####
    
    for prox in range(6):
        base_mat[prox, :, :] = np.array([[fc[fc_arr[prox, 0]], fc[fc_arr[prox, 2]], fc[fc_arr[prox, 2]]],
                                         [fc[fc_arr[prox, 2]], fc[fc_arr[prox, 1]], fc[fc_arr[prox, 3]]],
                                         [fc[fc_arr[prox, 2]], fc[fc_arr[prox, 3]], fc[fc_arr[prox, 1]]]])
    
    for i in range(1000):   #1000 for nats=2000
        for prox in range(6):
            for j in neighbors[i][prox]:
                fc_mat[i, j, :, :] = base_mat[prox, :, :]
                if prox in [2, 3, 4]:
                    for r in range(1,3):
                        if isclose(abs(ideal_distances[i, r, j]), nn_matrix[prox, 0], rel_tol=.001):
                            fc_mat[i, j, [r, 0], :] = fc_mat[i, j, [0, r], :]
                            fc_mat[i, j, :, [r, 0]] = fc_mat[i, j, :, [0, r]]
                            
                for r in range(3):
                    if ideal_distances[i, r, j] < 0:
                        fc_mat[i, j, r, :] = -fc_mat[i, j, r, :]
                        fc_mat[i, j, :, r] = -fc_mat[i, j, :, r]
                        
    #### 2ND SPECIES ####
                        
    for prox in range(6):
        base_mat[prox, :, :] = np.array([[fc_2[fc_arr[prox, 0]], fc_2[fc_arr[prox, 2]], fc_2[fc_arr[prox, 2]]],
                                         [fc_2[fc_arr[prox, 2]], fc_2[fc_arr[prox, 1]], fc_2[fc_arr[prox, 3]]],
                                         [fc_2[fc_arr[prox, 2]], fc_2[fc_arr[prox, 3]], fc_2[fc_arr[prox, 1]]]])
    
    for i in range(1000,2000):   #1000,2000 for nats=2000 for Fe,V
        for prox in range(6):
            for j in neighbors[i][prox]:
                fc_mat[i, j, :, :] = base_mat[prox, :, :]
                if prox in [2, 3, 4]:
                    for r in range(1,3):
                        if isclose(abs(ideal_distances[i, r, j]), nn_matrix[prox, 0], rel_tol=.001):
                            fc_mat[i, j, [r, 0], :] = fc_mat[i, j, [0, r], :]
                            fc_mat[i, j, :, [r, 0]] = fc_mat[i, j, :, [0, r]]
                            
                for r in range(3):
                    if ideal_distances[i, r, j] < 0:
                        fc_mat[i, j, r, :] = -fc_mat[i, j, r, :]
                        fc_mat[i, j, :, r] = -fc_mat[i, j, :, r]
                        
    
                        
    np.around(fc_mat, decimals=5)
    # out_path = 'C:\\Users\\biknb\\Downloads\\bkc\\Phonons_0.5\\'+syst+'\'+sys+'\\Phonopy/'+sys+'_'+root+'/'
    out_path = '/home/bkc/fall2023/'+syst+'/'+sys+'/Phonopy/'+sys+'_'+root+'/'
    out_filename = 'FORCE_CONSTANTS'
    with open(out_path+out_filename, 'w') as file:
        print(str(n) + ' ' + str(n), file=file)
        for i in range(n):
            for j in range(n):
                print(str(i+1) + ' ' + str(j+1), file=file)
                np.savetxt(file, fc_mat[i, j, :, :])
                
    
    os.chdir(out_path)
    command = 'phonopy -p -s band.conf'
    # print('here 1', m)
    os.system(command)
    # print('here 2', m)
    command = 'mv band.yaml band_' + str(m) + '.yaml'
    os.system(command)
    command = 'mv band.pdf band_' + str(m) + '.pdf'
    os.system(command)
    command = 'rm phonopy.yaml'
    os.system(command)
    print(root, m, )

sys.exit()



