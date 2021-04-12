#!/bin/bash

#================Dream network parameters=================
network='dream4'
numofnode=(10 100)
totalnet=(5)
networknum=(1 2 3 4 5)
discretization=('efd' 'ewd' 'kmeans')

path="$PWD"

results='results'
source='sourcefiles'
inputfiles='input-files'
code='log_regression'
groundtruth='groundtruth'
scores='scores'

mkdir ./$results
cd ./$results

mkdir ./$scores

for nn in "${numofnode[@]}"
do 
	folder=$nn'-gene'
	
	mkdir ./$folder
	cd ./$folder
	echo $folder
	#==============================PEPN-GRN-v3 code=======================================

		for disc in "${discretization[@]}"
		do
		mkdir ./$disc
		cd ./$disc
		echo $disc

				cp $path/$source/*.m .
				cp $path/$folder/*.mat .				
				cp $path/$folder/$inputfiles/$disc/*.mat .
	
				matlab -nodisplay -nosplash -nodesktop -r "run $code($nn,$totalnet).m; quit;"				

				#===================== ROC and PR plots generated for each network =============
				cp $path/$folder/$groundtruth/*.mat .				
				for netnum in "${networknum[@]}"
				do
				matlab -nodisplay -nosplash -nodesktop -r "run rankwise_roc_pr_plot($nn,$netnum).m; quit;"

				cp 'all_edges_'$nn'_'$netnum'.txt' $path/$results/$scores/'all_edges_'$nn'_'$netnum'_'$disc'.txt'

				done
				#=========================================================
		cd ..		
		done
done
cd ..

cd ./$scores
cp $path/$source/dream4_scoring.py .
python3 dream4_scoring.py

cd ..
cd ..


