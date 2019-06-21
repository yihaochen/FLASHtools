#!/bin/bash

##############################################################################
# Modified from "An rclone backup script" by Jared Males (jaredmales@gmail.com)
# 
# Copyright (C) 2018 Jared Males <jaredmales@gmail.com>
#
# This script is licensed under the terms of the MIT license.
# https://opensource.org/licenses/MIT
##############################################################################

#### rclone sync options

# Directories to be managed
DIRS="
20180801_L430_rc10_beta07
20180802_L438_rc10_beta07
20180803_L446_rc10_beta07
20180906_L430_rc30_beta07
20180907_L446_rc30_beta07
20180824_L438_rc30_beta07
20181202_L438_rc100_beta07
20181203_L446_rc100_beta07
"

# Scratch variable (not set in cron job)
SCRATCH=/scratch/03551/tg828509

# Path to the rclone executable
RCLONE=${HOME}/.local/bin/rclone

#---- Edit this to the desired destination
DEST=gdrive:2018_production_runs

#---- This is the path to a file with a list of exclude rules
EXCLUDEFILE=${HOME}/.rclone_cron/excludes 

#---- Name of exclude file
# NOTE: you need "v1.39-036-g2030dc13β" or later for this to work.
#EXIFPRESENT=.rclone-ignore

#---- The bandwidth time table
BWLIMIT="8.9M"

#---- Minimum age for identifying checkpoint files (minutes)
MINAGE=10

#---- Number of transfers to do in parallel.
TRANSFERS=1

#---- Location of cron log
CRONLOG=${HOME}/.rclone_cron/rclone_cron.log

#---- Location of rclone log
LOGS='-v --stats-one-line --log-file='$CRONLOG


###################################################
## Locking Boilerplate from https://gist.github.com/przemoc/571091
## Included under MIT License:
###################################################

## Copyright (C) 2009 Przemyslaw Pawelczyk <przemoc@gmail.com>
##
## This script is licensed under the terms of the MIT license.
## https://opensource.org/licenses/MIT
#
# Lockable script boilerplate

### HEADER ###

LOCKFILE="${HOME}/.rclone_cron/`basename $0`.lock"
LOCKFD=99

# PRIVATE
_lock()             { flock -$1 $LOCKFD; }
_no_more_locking()  { _lock u; _lock xn && rm -f $LOCKFILE; }
_prepare_locking()  { eval "exec $LOCKFD>\"$LOCKFILE\""; trap _no_more_locking EXIT; }

# ON START
_prepare_locking

# PUBLIC
exlock_now()        { _lock xn; }  # obtain an exclusive lock immediately or fail
exlock()            { _lock x; }   # obtain an exclusive lock
shlock()            { _lock s; }   # obtain a shared lock
unlock()            { _lock u; }   # drop a lock

###################################################
# End of locking code from Pawelczyk
###################################################


#make a log entry if we exit because locked
exit_on_lock()      { echo $(date)' | rclone_cron.sh already running.' >> $CRONLOG; exit 1; }


#Now check for lock
exlock_now || exit_on_lock
#We now have the lock.

#Log startup
echo $(date)' | starting rclone_cron.sh......' >> $CRONLOG

#Loop through the directories!
for dir in $DIRS
do
	echo $(date)' | checking '${dir} >> $CRONLOG
	# Find the modification time of the last checkpoint file; only consider checkpoint files older than 30 minutes to avoid writing in progress
	LASTCHKTIME=$(find ${SCRATCH}/${dir}/Group_L4??_hdf5_chk_???? -maxdepth 0 -type f -mmin +$MINAGE -exec stat --format '%Y' "{}" \; | sort -nr | head -1)

	# Convert it into minutes from now
	MMIN=$((($(date +%s) - $LASTCHKTIME)/60))

	# Find all checkpoint files older than the last one (by at least 10 minutes)
	chkfiles=$(find ${SCRATCH}/${dir}/Group_L4??_hdf5_chk_???? -maxdepth 0 -type f -mmin +$(($MMIN+10)) 2>/dev/null)
	if [ ! -z "$chkfiles" ]
	then
		for f in $chkfiles
		do
			# Age of the checkpoint file (in hours)
			CHKAGE=$((($(date +%s) - $(stat --format '%Y' $f))/3600))
			HDLINK=$(stat --format '%h' $f)
			# Move the checkpoint file if older than 10 days and not hard linked
			if [[ $CHKAGE -gt 240 ]] && [[ $HDLINK -eq 1 ]]
			then
				CHKMOVED=true
				echo $(date)' | moving '$f' to data/' >> $CRONLOG
				mv -n $f -t ${SCRATCH}/${dir}/data/
			fi
		done

	fi

	files=$(find ${SCRATCH}/${dir}/Group_L4??_hdf5_{plt_cnt,part}_???? -maxdepth 0 -type f -mmin +$MMIN 2>/dev/null)
	if [[ ! -z "$files" ]] || [[ ! -z "$CHKMOVED" ]]
	then

		# Copy log files to log/
		logfiles=$(find ${SCRATCH}/${dir}/ -maxdepth 1 -name "Group_L4??.???" -o -name "nozzleVec.dat" -o -name "Group_L4??.*.out" -type f 2>/dev/null)
		cp -a ${logfiles} -t ${SCRATCH}/${dir}/log/

		# Move plt_cnt and part files to data/
		for f in $files
		do
			mv -n $f -t ${SCRATCH}/${dir}/data/
		done

		# Now do the rclone copy
		$RCLONE copy ${SCRATCH}/${dir}/log/ ${DEST}/${dir}/log/ ${LOGS}
		$RCLONE copy ${SCRATCH}/${dir}/data/ ${DEST}/${dir}/data/ --transfers=$TRANSFERS --bwlimit "${BWLIMIT}" --exclude-from ${EXCLUDEFILE} ${LOGS}
	fi
done


#log success
echo $(date)' | completed rclone_cron.sh.' >> $CRONLOG
echo '================================================================================' >> $CRONLOG
echo '' >> $CRONLOG

#release the lock
unlock

exit
