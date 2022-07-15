#!/bin/bash
#$ -cwd                     # Run the code from the current directory
#$ -V                       # Export environment to job
#$ -j y                     # Merge the standard output and standard error
#$ -l h_rt=240:00:00        # Limit each task to 10 days
#$ -l h_vmem=5G             # Request 1GB RAM
#$ -l node_type=nxv         # Request fastest nodes
#$ -t 1-3

##################################################################################################
## Creating architecture to run a test in MCMCtree to check if ST are OK		   				##
##																								##
## Contact Sandra Alvarez-Carretero for any doubts about this script: sandra.ac93@gmail.com		##
##################################################################################################
# Create a variable with the pathway to the home directory
home_dir=$( pwd )
if [[ $SGE_TASK_ID -eq 2 ]]
then 
$SGE_TASK_ID=$( echo 4 )
fi
# RAxML is exported in my PATH, so there is no need for a `bin` variable here

# Get the date in order to append it to the log file
#today=`date '+%d%m%y'`

# Print the time for every line in the log file
# At the same time, use "3>&1>" and "2>&1" in order to print the help manual if required
# 2>&1 --> stderr ("2") to the same log that stdout ("1")
# 3>&1 --> Creates a new file descriptor and redirects it to "1" (stdout)
# Without printing help manual: exec &> >(while read line; do echo "$(date +'[%d-%m-%Y | %H:%M:%S]') $line" >> log_$today.log; done;)
exec 3>&1> >(while read line; do echo "$line" >> log.raxml.data$SGE_TASK_ID.txt; done;) 2>&1

#############################################
## START SETTING FILES/FOLDER ARCHITECTURE ##
#############################################
echo The analyses will be carried out in the directory $home_dir/$SGE_TASK_ID/
printf "\n"
# Move to the home directory
cd $home_dir/$SGE_TASK_ID
# Get file names and set vars
seq=$( ls *aln )
name=$( echo $seq |sed 's/\.aln//' )

#################
## Run RAxML   ##
#################
printf "\nRunning RAxML...\n"
# NOTE: 
#   -f a  This runs a rapid bllootstrap analysis. The random seed is set by 
#         `-x 12345` and I set it to `-# 500` iterations 
raxmlHPC -f a -m GTRGAMMA -p 12345 -# 100 -x 12345 -# 500 -s $seq -o Duck,Chicken -n $name

printf "\n"
echo RAxML finished"!"
printf "\n"
