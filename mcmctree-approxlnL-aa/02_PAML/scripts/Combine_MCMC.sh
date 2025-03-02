#!/bin/bash

curr_dir=$( pwd )
dat=$1
dirname=$2
seqchains=$3
clock=$4
traces=$6
trace_name=$7
# Num samples specified in MCMCtree (calculated with arg5)
mcmctree_numsamp=$(( $5 + 2 ))
printf "\nNumber of samples specified in the control file: "$5"\n"
printf "==========================================================\n"
#printf "mcmctree_numsapm="$mcmctree_numsamp"\n"

mkdir -p $dirname
if [[ $traces =~ [Yy] ]]
then
if [[ ! -d mcmcf4traces_$trace_name ]]
then
mkdir -p mcmcf4traces_$trace_name
fi
fi
count=0
for i in $seqchains
do
	count=$(( count + 1 ))	
	if [[ ! -f $dat/$i/mcmc.txt ]]
	then 
		printf "Sorry, no samples for run"$i" under "$clock" ...\n"
		printf $clock"\trun"$i"\n" >> "Not_collected_samples.tsv"
	else 
		printf "Parsing dat for run"$i" under "$clock" ... ... \n"
		end=$( wc -l $dat/$i/mcmc.txt | awk '{print $1}' )
		if [ $count -eq 1 ]
		then
			begin=1
			# if [[ $mcmctree_numsamp =~ ^[0-9]+$ ]]
			# then
			# 	echo mcmctree_numsamp is valid
			# fi
			# if [[ $end =~ ^[0-9]+$ ]]
			# then
			# 	echo end is valid
			# fi
			# echo $end is equal to $mcmctree_numsamp "?"
			if [ ! $end -eq $mcmctree_numsamp ]
			then
				count2=0
				while IFS= read -r line || [ -n "$line" ]; do
					count2=$(( count2 + 1 ))
				done < $dat/$i/mcmc.txt
				if [ $end -eq $count2 ]
				then
					printf "   [[ One additional line will be removed to match analyses in R ]]\n"
					end=$(( end - 1 ))
					sed -n ''${begin}','${end}'p' $dat/$i/mcmc.txt > $dirname/mcmc.txt
					sed -n ''${begin}','${end}'p' $dat/$i/mcmc.txt > mcmcf4traces_$trace_name/mcmc_$i.txt
				else
					sed -n ''${begin}','${end}'p' $dat/$i/mcmc.txt > $dirname/mcmc.txt
					sed -n ''${begin}','${end}'p' $dat/$i/mcmc.txt > mcmcf4traces_$trace_name/mcmc_$i.txt
				fi
			else
				printf "   [[ You collected all the samples specified in your control file for this chain! ]]\n"
				sed -n ''${begin}','${end}'p' $dat/$i/mcmc.txt > $dirname/mcmc.txt
				sed -n ''${begin}','${end}'p' $dat/$i/mcmc.txt > mcmcf4traces_$trace_name/mcmc_$i.txt
			fi
		else 
			begin=2 
			if [ ! $end -eq $mcmctree_numsamp ]
			then
				count2=0
				while IFS= read -r line || [ -n "$line" ]; do
					count2=$(( count2 + 1 ))
				done < $dat/$i/mcmc.txt
				if [ $end -eq $count2 ]
				then
					printf "   [[ One additional line will be removed to match analyses in R ]]\n"
					end=$(( end - 1 ))
					sed -n ''${begin}','${end}'p' $dat/$i/mcmc.txt >> $dirname/mcmc.txt
					sed -n '1,'${end}'p' $dat/$i/mcmc.txt > mcmcf4traces_$trace_name/mcmc_$i.txt
				else
					sed -n ''${begin}','${end}'p' $dat/$i/mcmc.txt >> $dirname/mcmc.txt
					sed -n '1,'${end}'p' $dat/$i/mcmc.txt > mcmcf4traces_$trace_name/mcmc_$i.txt
				fi
			else
				printf "   [[ You collected all the samples specified in your control file for this chain! ]]\n"
				sed -n ''${begin}','${end}'p' $dat/$i/mcmc.txt >> $dirname/mcmc.txt
				sed -n '1,'${end}'p' $dat/$i/mcmc.txt > mcmcf4traces_$trace_name/mcmc_$i.txt
			fi
		fi
	fi
done  

# NOTE:
# After this script, you need to copy dummy aln, ctl file, and tree file 
# to generate the FigTree file using the option -1 !