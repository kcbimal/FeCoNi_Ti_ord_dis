# FeCoNi_Ti_ord_dis
################## Creating Directories #####################
0:Create dir: "Phonons_0.5" with all files in it

1: a) Create dir: "ord" and "dis" inside dir. "Phonons_0.5"  
                i) with sub-dir  "FeTi", "CoTi", "NiTi" (empty) for each
   b) Create dir: "positions_file" inside dir. "Phonons_0.5" 


################ Generating Atomic Positions ###################
2: a)Go to dir. "positions_file" and run 
     "cesar_atomic_position_generator_Fe=0.5.py" 

      Rename dir "**2.90" to "**2.9" manually
      Rename dir "**3.00" to "**3.0" manually
      Rename dir "**3.10" to "**3.1" manually
or,
  b) if there is already a "positions_file" directory exists and the 
     lattice param. is within our required range, no need to run step 2. a)

   NOTE: JUST CREATE POSITIONS: AS THEY ARE TEMP. INDEPENDENT FOR ANY GIVEN  
         LATTICE PARAM.
         THESE SAME POSITIONS CAN BE USED FOR ALL STRUCTURES 
         (ord/dis/300K/../1300K, etc) AND TEMPERATURES.
         It's BECAUSE "Step 4" MAKES THE POSITIONS CHANGES ACCORDINGLY.


################ Stay in dir. "Phonons_0.5" and; ############### 
3: Run 
	a) "directories_ord.py" for ordered 
            (directories_outdated.py is ignored)
        b) "directories_dis.py" for dis-ordered 
            (directories_outdated.py is ignored)

4: To generate atomic positions: 

         a) Run "get_disordered_str_from_SPOSCAR_V2_all.py" 
            for disordered structure
	 b) Run "get_ord_atom_position_from_SPOSCAR_all.py" 
            for ordered structure
  

################## Running Simulation in PC ####################
5: a) Run "LAMMPS" in PC using ./lammps_multirun.sh
   b) Run "create_ideal_unit_from_ideal_lammps.py" 
      to get the ideal unit file 
      (for  ord structure only) 
      (dis-ord don't need.It reads "dump_***.0" instead of ideal_unit***)


################## Running Simulation in Jakar ##################
6: Open Jakar and run using multirun_hjob.sh: 

	a) HELD FC's code "BvK" 
 
	b) Phonopy for "Phonopy_python" and "Phonopy_mesh" 
 

################## Running Simulation in PC ####################
7: Run "vib_entropy_from_pdos_bkc.py" for entropy in PC
       (for all lattice but single temp at a time)
