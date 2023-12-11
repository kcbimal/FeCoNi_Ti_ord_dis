for f in */ ; 
	do 
	
	cd ${f} ; 
	sbatch hjob.sh
	
	cd ../ ; 
done