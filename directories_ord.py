#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 18:23:13 2021

@author: biknb
"""

import os, shutil
sys = 'CoTi'         #'FeTi'    #'NiTi'
system = 'Co Ti'     #'Fe Ti'   #'Ni Ti'

#check line "36" and just change line "33" and don't change anything below that line
#check line 53, 81
#This would be the location of the files. potential1/2 are the fe.meam and meam.library files
syst = 'ord'
parent_dir = "C:\\Users\\biknb\\Downloads\\bkc\\Phonons_0.5\\"+syst+"\\"
par_dir = "C:\\Users\\biknb\\Downloads\\bkc\\Phonons_0.5\\"
l_dir = "C:\\Users\\biknb\\Downloads\\bkc\\Phonons_0.5\\"+syst+"\\"+sys+"\\"
eam = par_dir + sys+'.meam'
libmeam = par_dir+ sys + '.library.meam'
lammps_in_file = par_dir  + sys + '.in'
ljob = par_dir + 'lammps_multirun.sh'


#Create the LAMMPS directory with the V and T folders & input files
lammps_dir = os.mkdir(l_dir + 'LAMMPS')
lammps_path = "C:\\Users\\biknb\\Downloads\\bkc\\Phonons_0.5\\"+syst+"\\"+sys+"\\LAMMPS\\"
p_path = "C:\\Users\\biknb\\Downloads\\bkc\\Phonons_0.5\\positions_file\\"
shutil.copy(ljob, lammps_path )

vols = [2.91, 2.92, 2.93, 2.94, 2.95, 2.96, 2.97, 2.98, 2.99, 3.00, 3.01, 3.02, 3.03, 3.04, 3.05, 3.06, 3.07, 3.08, 3.09, 3.10] #you can input here a list or single volume
temp = [10, 50, 90] #CoTi
# temp = [250, 290, 330]  #NiTi

################ don't change anything below this line ####################### 

with open(lammps_in_file,'r') as msg:
            data = msg.read()
for x in vols:
    for y in temp:
        dir_name = lammps_path + sys + '_' + str(x) + '_' + str(y) + 'K'
        dirs = os.mkdir(dir_name)
        shutil.copy(eam, dir_name)
        shutil.copy(p_path + 'FeTi_' + str(x) + '\\atoms_positions.data', dir_name)
        shutil.copy(libmeam, dir_name)
        
        rename = data.replace("latt_par", str(x)).replace("tmp", str(y)).replace("sys", sys) 
        with open(parent_dir + sys + '_' + str(x) + '_' + str(y) +'K.in','w') as out_msg:
            out_msg.write(rename)
        shutil.move(parent_dir + sys + '_' + str(x) + '_' + str(y) +'K.in', dir_name)
        
        
# #create the force constants directory that uses the .csv files from previous calculation
force_constants_dir = os.mkdir(l_dir + 'Force_constants')
force_constants_path = "C:\\Users\\biknb\\Downloads\\bkc\\Phonons_0.5\\"+syst+"\\"+sys+"\\Force_constants\\"
bvk_py_file ="C:\\Users\\biknb\\Downloads\\bkc\\Phonons_0.5\\BvK Class BCC np ord.py"
multijob_file =  par_dir + 'hjob.sh'
hjob = par_dir + 'multirun_hjob.sh'
shutil.copy(hjob, force_constants_path)


with open(bvk_py_file,'r') as msg:
            data_f = msg.read()
                
for x in vols:
    for y in temp:
        dir_name_f = force_constants_path + sys + '_' + str(x) + '_' + str(y) + 'K'
        dirs_f = os.mkdir(dir_name_f)
        #shutil.copy(bvk_py_file, dir_name_f)
        shutil.copy(multijob_file, dir_name_f)
        
        rename_bvk = data_f.replace("latt_par", str(x)).replace("tmp", str(y))
        # Replace 'Fe Ti' with 'sys'
        rename_bvk = rename_bvk.replace('FeTi', sys)
        
        # Replace 'ord' with 'syst'
        rename_bvk = rename_bvk.replace('ord', syst)
        
        with open(force_constants_path +'BvK' + '_' + str(x) + '_' + str(y) + '.py','w') as out_msg:
              out_msg.write(rename_bvk)
        shutil.move(force_constants_path +'BvK' +'_' + str(x) + '_' + str(y) + '.py', dir_name_f)
 
        
#creates the phonopy directory (generates the FORCE_CONSTANTS file needed to compute phonon dispersions)
phonopy_dir = os.mkdir(l_dir + 'Phonopy')
phonopy_path = "C:\\Users\\biknb\\Downloads\\bkc\\Phonons_0.5\\"+syst+"\\"+sys+"\\Phonopy\\"
shutil.copy(hjob, phonopy_path)
disp_gen_file = par_dir + 'dispersion_generator_'+syst+'.py'
bandconf_file = par_dir + 'band.conf'
poscar_file   = par_dir + 'POSCAR'
multirun_file =  par_dir + 'hjob.sh'
multi_file =  par_dir + 'mjob.sh'
meshconf_file     = par_dir + 'mesh.conf'

# with open(disp_gen_file,'r') as msg:
#             data_f = msg.read()
# for x in vols:
#     for y in temp:
#         dir_name_disp = phonopy_path + sys + '_' + str(x) + '_' + str(y) + 'K'
#         dirs_dis = os.mkdir(dir_name_disp)
#         shutil.copy(bandconf_file, dir_name_disp)
#         shutil.copy(multirun_file, dir_name_disp)
#         rename_disp = data_f.replace("latt_param", str(x)).replace("tmp", str(y))
        
#         with open(phonopy_path +'disp_gen' + '_' + str(x) + '_' + str(y) + '.py','w') as out_msg:
#               out_msg.write(rename_disp)
#         shutil.copy(phonopy_path +'disp_gen' +'_' + str(x) + '_' + str(y) + '.py', dir_name_disp)
        #shutil.move(disp_gen_file, phonopy_path)

with open(disp_gen_file, 'r') as msg:
    data_f = msg.read()

for x in vols:
    for y in temp:
        dir_name_disp = phonopy_path + sys + '_' + str(x) + '_' + str(y) + 'K'
        dirs_dis = os.mkdir(dir_name_disp)
        shutil.copy(bandconf_file, dir_name_disp)
        shutil.copy(multirun_file, dir_name_disp)
        shutil.copy(multi_file, dir_name_disp)
        
        # Replace "latt_param" with the current value of 'x' and "tmp" with the current value of 'y'
        rename_disp = data_f.replace("latt_param", str(x)).replace("tmp", str(y))
        
        # Replace 'FeTi' with 'sys'
        rename_disp = rename_disp.replace('FeTi', sys)
        
        # Replace 'ord' with 'syst'
        rename_disp = rename_disp.replace('ord', syst)
        
        with open(phonopy_path + 'disp_gen' + '_' + str(x) + '_' + str(y) + '.py', 'w') as out_msg:
            out_msg.write(rename_disp)
        
        shutil.copy(phonopy_path + 'disp_gen' + '_' + str(x) + '_' + str(y) + '.py', dir_name_disp)

        
# #change 'Fe V' in sposcar file "system", if not needed, comment all below this line
# with open(sposcar_file,'r') as msg:
#             data_s = msg.read()
# for x in vols:
#     for y in temp:
#         dir_name_disp = phonopy_path + sys + '_' + str(x) + '_' + str(y) + 'K'
#         shutil.copy(sposcar_file, dir_name_disp)
#         rename_sys = data_s.replace('system', system)
#         with open(phonopy_path +'SPOSCAR','w') as out_msg:
#               out_msg.write(rename_sys)
#         shutil.copy(phonopy_path +'SPOSCAR', dir_name_disp)

#change 'FeV' in mesh file "sys", if not needed, comment all below this line
with open(meshconf_file,'r') as msg:
            data_m = msg.read()
for x in vols:
    for y in temp:
        dir_name_disp = phonopy_path + sys + '_' + str(x) + '_' + str(y) + 'K'
        shutil.copy(meshconf_file, dir_name_disp)
        rename_sys = data_m.replace('sys', sys)
        with open(phonopy_path +'mesh.conf','w') as out_msg:
              out_msg.write(rename_sys)
        shutil.copy(phonopy_path +'mesh.conf', dir_name_disp)

#change lattice parameter in poscar file "latt_param", if not needed, comment all below this line
# with open(poscar_file,'r') as msg:
#             data_p = msg.read()
# for x in vols:
#     for y in temp:
#         dir_name_disp = phonopy_path + sys + '_' + str(x) + '_' + str(y) + 'K'
#         shutil.copy(poscar_file, dir_name_disp)
#         rename_lat = data_p.replace('latt_param', str(x))
#         with open(phonopy_path +'POSCAR','w') as out_msg:
#               out_msg.write(rename_lat)
#         shutil.copy(phonopy_path +'POSCAR', dir_name_disp)
with open(poscar_file, 'r') as msg: 
    data_p = msg.read()

for x in vols:
    for y in temp:
        dir_name_disp = phonopy_path + sys + '_' + str(x) + '_' + str(y) + 'K'
        shutil.copy(poscar_file, dir_name_disp)
        
        # Replace 'latt_param' with the current value of 'x'
        rename_lat = data_p.replace('latt_param', str(x))
        
        # Replace 'Fe Ti' with 'system'
        rename_lat = rename_lat.replace('Fe Ti', system)
        
        with open(phonopy_path + 'POSCAR', 'w') as out_msg:
            out_msg.write(rename_lat)
        
        shutil.copy(phonopy_path + 'POSCAR', dir_name_disp)
     
        

        