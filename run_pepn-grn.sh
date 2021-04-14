#!/bin/bash

#Set variant to 1 for PEPN-GRN_v1, 2 for PEPN-GRN_v2 and 3 for PEPN-GRN_v3 implementation

variant=1

#================Dream network parameters=================
network='dream4'
numofnode=(10 100)
totalnet=(5)
networknum=(1 2 3 4 5)
discretization=('efd' 'ewd' 'kmeans')
disclevel='3bin'

path="$PWD"

results='PEPN-GRN_v'$variant'-results'
source='source_v'$variant
code='code_v'$variant
log_regression='log_regression'


mkdir ./$results
cd ./$results



for nn in "${numofnode[@]}"
do 
	if  [ $variant = "3" ]
	then folder=$nn'-gene'
		mkdir ./$folder
	fi


	for netnum in "${networknum[@]}"
	do
	currentNet=$network'_'$nn'_'$netnum

	mkdir ./$currentNet
	cd ./$currentNet
	echo $currentNet
		#==============================PEPN-GRN code=======================================

		for disc in "${discretization[@]}"
		do
		mkdir ./$disc
		cd ./$disc
		echo $disc
				#gdtruth='gdTruth_'$currentNet
				cp $path/$source/*.m .
				cp $path/multi-bin-disc-dream4-data-repository/$currentNet/*.{mat,tsv} .
				cp $path/multi-bin-disc-dream4-data-repository/$currentNet/$disc/*.mat .
				cp $path/multi-bin-disc-dream4-data-repository/'genenames_'$nn'gene.mat' .
				
				if [ $variant = "1" ]
				then matlab -nodisplay -nosplash -nodesktop -r "run $code($nn,$netnum,'$disc','$disclevel',$variant).m; run rankwise_roc_pr_plot($nn,$netnum).m; quit;"
				elif [ $variant = "2" ]
				then matlab -nodisplay -nosplash -nodesktop -r "run $code($nn,$netnum,'$disc','$disclevel',$variant).m; run rankwise_roc_pr_plot($nn,$netnum).m; quit;"
				elif [ $variant = "3" ]
				then matlab -nodisplay -nosplash -nodesktop -r "run $code($nn,$netnum,'$disc','$disclevel',$variant).m; run edge_features($nn,$netnum).m; quit;"

					if [ ! -d ../../$folder/$disc ]; then
						mkdir ../../$folder/$disc
					fi

					cp 'pos_act_edges_'$nn'_'$netnum'.mat' ../../$folder/$disc/
					cp 'neg_act_edges_'$nn'_'$netnum'.mat' ../../$folder/$disc/
					cp 'pos_inh_edges_'$nn'_'$netnum'.mat' ../../$folder/$disc/
					cp 'neg_inh_edges_'$nn'_'$netnum'.mat' ../../$folder/$disc/	

				fi
							
		cd ..		
		done

	cd ..			
	done
		
	if [ $variant = "3" ]
	cd ./$folder
		then for disc in "${discretization[@]}"
			do
			cd ./$disc
			echo $disc
			cp $path/multi-bin-disc-dream4-data-repository/'genenames_'$nn'gene.mat' .
			cp $path/$source/$log_regression/*.m .

			matlab -nodisplay -nosplash -nodesktop -r "run $log_regression($nn,$totalnet).m; quit;"				

				for netnum in "${networknum[@]}"
				do
				currentNet=$network'_'$nn'_'$netnum
				cp $path/multi-bin-disc-dream4-data-repository/$currentNet/'all_pos_edge_'$nn'_'$netnum'.mat' .
				cp $path/multi-bin-disc-dream4-data-repository/$currentNet/'all_neg_edge_'$nn'_'$netnum'.mat' .
	
				matlab -nodisplay -nosplash -nodesktop -r "run rankwise_roc_pr_plot($nn,$netnum).m; quit;"
				done
			cd ..
			done

		cd ..
	fi			
done
	
cd ..


