#!/bin/bash

#SBATCH -J 0801L430_dev						# job name
#SBATCH -d singleton
#SBATCH -o Group_L430.%j.out				# output and error file name (%j expands to jobID)
#SBATCH -p skx-dev								# queue (partition) -- normal, development, etc.
#SBATCH -N 4									# total # of nodes
#SBATCH -n 192									# total # of mpi tasks
#SBATCH -t 02:00:00								# run time (hh:mm:ss) 
#SBATCH -A TG-AST170026							# Project number
##SBATCH --mail-user=ychen@astro.wisc.edu
##SBATCH --mail-type=begin  # email me when the job starts
##SBATCH --mail-type=end    # email me when the job finishes

module restore

if [ -f .dump_restart ]; then
    rm .dump_restart
fi

python ${HOME}/set_filenumbers.py 7200

#env | /bin/grep SLURM

# Make a copy of current flash.par file
cp flash.par flash.par.${SLURM_JOB_ID}.$(date "+%Y%m%d.%H%M%S")

# Set up a timer to stop FLASH before timeout
#echo "touch .dump_restart" | at $(date -d "+2 hours -1 minutes" "+%H:%M %b %d")

ibrun ./flash4             # run the MPI executable
