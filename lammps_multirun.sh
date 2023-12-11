for f in */ ; 
	do 
	
	cd ${f} ; 
	 mpiexec -np 2 lmp_mpi -in *.in -e screen
         # ~/lammps-stable_3Mar2020/src/lmp_mpi -in *.in
	
	cd ../ ; 
done
