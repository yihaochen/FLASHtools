#!/bin/bash

#SBATCH -J 0801L430							# job name
#SBATCH -o Group_L430.%j.out				# output and error file name (%j expands to jobID)
#SBATCH -p skx-normal							# queue (partition) -- normal, development, etc.
#SBATCH -N 8									# total # of nodes
#SBATCH -n 384									# total # of mpi tasks
#SBATCH -t 48:00:00								# run time (hh:mm:ss) 
#SBATCH -A TG-AST140042							# Project number
#SBATCH --mail-user=ychen@astro.wisc.edu
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finishes

module restore intel18

python ${HOME}/set_filenumbers.py

# Make a copy of current flash.par file
cp flash.par flash.par.${SLURM_JOB_ID}.$(date "+%Y%m%d.%H%M%S")

if [ -f .dump_restart ]; then
    rm .dump_restart
fi

# Set up a timer to stop FLASH before timeout
# echo "touch .dump_restart" | at $(date -d "+2 days -2 minutes" "+%H:%M %b %d")
ibrun ./flash4             # run the MPI executable
