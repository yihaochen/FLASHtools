#!/bin/sh
#SBATCH -J Jet_h0_10Myr
#SBATCH --partition=astro3				# default "univ" if not specified
#SBATCH --time=7-00:00:00				# run time in days-hh:mm:ss
#SBATCH --nodes=4						# require nodes
#SBATCH --ntasks-per-node=20			# default 16 if this line not specified
#SBATCH --mem-per-cpu=6000				# RAM per CPU core, in MB (default 4 GB/core)
#SBATCH --error=MHD_Jet_10Myr.%j.out
#SBATCH --output=MHD_Jet_10Myr.%j.out
#Make sure to change the above two lines to reflect your appropriate
# file locations for standard error and output                                                      
#SBATCH --mail-user=ychen@astro.wisc.edu
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finishes

#Now list your executable command (or a string of them).
if [ -f .dump_restart ]; then
	rm .dump_restart
fi
cp flash.par flash.par.$(date '+%Y%m%d.%H%M%S')

# This fixes the FLASH using only half of cores on astro3 partition
# when compiled with mvapich2
export MV2_CPU_BINDING_LEVEL="core"

mpirun ./flash4
