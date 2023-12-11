# -*- coding: utf-8 -*-
"""
@author: cesar                        
"""             
                      
# Code specific for ordered composition 0.5 (B2)

import numpy as np
import os
import shutil
import pandas as pd
import fileinput

sys = 'FeTi' #'CoTi', 'NiTi'
tmp = '300'
pos_dir = "C:\\Users\\biknb\\Downloads\\bkc\\Phonons_0.5\\"
positions_file = os.mkdir(pos_dir + 'positions_file')

parent_dir = "C:\\Users\\biknb\\Downloads\\bkc\\Phonons_0.5\\positions_file\\"


lattice_parameters = []
vanadium_composition = 0.5
lattice_value_for = 2.8

vanadium_atoms_unit_cell = int(vanadium_composition * 16)

for i in range(39):
    lattice_parameters.append(lattice_value_for)
    lattice_value_for += 0.01

range_unit_cell = np.arange(1, 17)

# Unit Cell for 0.5 composition 
random_int_array = np.array([9, 10, 11, 12, 13, 14, 15, 16])

df_random_numbers = pd.DataFrame(random_int_array)
df_random_numbers.to_csv(parent_dir +'Random array for atomic configuration.csv')

for lattice_parameter in lattice_parameters:

    # Unit cell length for FeV
    unit_cell_length = 2*lattice_parameter
    
    # Defining unit cell
    basis = np.array([[1.0, 0.0, 0.0], 
                      [0.0, 1.0, 0.0],
                      [0.0, 0.0, 1.0]])*unit_cell_length
    
    base_atoms = np.array([[0.0, 0.0, 0.0], 
                           [0.5, 0.0, 0.0], 
                           [0.0, 0.5, 0.0],
                           [0.5, 0.5, 0.0],
                           [0.0, 0.0, 0.5],
                           [0.0, 0.5, 0.5],
                           [0.5, 0.0, 0.5],
                           [0.5, 0.5, 0.5],
                           [0.75, 0.25, 0.25],
                           [0.25, 0.75, 0.25],
                           [0.25, 0.25, 0.75],
                           [0.75, 0.75, 0.75],
                           [0.25, 0.25, 0.25],
                           [0.75, 0.25, 0.75],
                           [0.25, 0.75, 0.75],
                           [0.75, 0.75, 0.25]])*unit_cell_length
    
    # Size of the system cell in unit cells
    # assuming a cubic cell, starting at the origin
    
    system_size = 5
    
    # Generate atom positions
    
    positions = []
    for i in range(system_size):
        for j in range(system_size):
            for k in range(system_size):
                base_position = np.array([i, j, k])
                cart_position = np.inner(basis.T, base_position)
                for atom in base_atoms:
                    positions.append(cart_position + atom)
    
    sub_directory = f'{sys}_{lattice_parameter:.2f}'#'Simulation_Lattice_{}'.format(round(lattice_parameter, 2))
    # parent_dir = os.getcwd()
    
    new_directory = os.path.join(parent_dir, sub_directory)
    os.makedirs(new_directory)
    
    os.chdir(new_directory)
    
    # Copies potentials and input files to each of the directories
    
    potential_library_dest = os.path.join(new_directory, ""+sys+".library.meam")
    shutil.copy("C:\\Users\\biknb\\Downloads\\bkc\\Phonons_0.5\\"+sys+".library.meam", potential_library_dest)
    
    potential_FeV_dest = os.path.join(new_directory, ""+sys+".meam")
    shutil.copy("C:\\Users\\biknb\\Downloads\\bkc\\Phonons_0.5\\"+sys+".meam", potential_FeV_dest)
    
    input_file_dest = os.path.join(new_directory, ""+sys+".in")
    shutil.copy("C:\\Users\\biknb\\Downloads\\bkc\\Phonons_0.5\\"+sys+".in", 
                input_file_dest)

    
    
    with open('atoms_positions.data', 'w') as fdata:
            #First line included as comment
            fdata.write('Crystalline Al atoms - written by EnCodeVentor tutorial\n\n')
            
            #-----Header-----#
            # Specify number of atoms and atom types
            fdata.write('# lattice parameter {} \n \n'.format(round(lattice_parameter, 2)))
            fdata.write('{} atoms\n'.format(len(positions)))
            fdata.write('{} atom types\n'.format(2))
            
            # Specify box dimensions
            fdata.write('{} {} xlo xhi\n'.format(0.0, system_size*unit_cell_length))
            fdata.write('{} {} ylo yhi\n'.format(0.0, system_size*unit_cell_length))
            fdata.write('{} {} zlo zhi\n'.format(0.0, system_size*unit_cell_length))
            fdata.write('\n')
            
            # Atoms section
            fdata.write('Atoms\n\n')
            
            # Write each position
            
            unit_cell_index = 1
            
            for i,pos in enumerate(positions):
                if unit_cell_index in random_int_array:
                    fdata.write('{} 2 {} {} {}\n'.format(i+1, *pos))
                    unit_cell_index += 1
                else:
                    fdata.write('{} 1 {} {} {}\n'.format(i+1, *pos))
                    unit_cell_index += 1
                if unit_cell_index > 16:
                    unit_cell_index = 1  
    
    #os.popen("lammps-shell.exe -in FeV_external_positions_input_file.in")
    
    os.chdir(parent_dir)
 
    
#change 'tmp' in .in file 
root_directory = "C:\\Users\\biknb\\Downloads\\bkc\\Phonons_0.5\\positions_file\\"
search_string, replacement_string = 'tmp', tmp

for dirpath, dirnames, filenames in os.walk(root_directory):
    for filename in filenames:
        if filename.endswith('.in'):
            file_path = os.path.join(dirpath, filename)
            with fileinput.FileInput(file_path, inplace=True, backup='.bak') as file:
                for line in file:
                    print(line.replace(search_string, replacement_string), end='')
 
