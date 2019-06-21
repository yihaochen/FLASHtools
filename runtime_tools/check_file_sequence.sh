#!/bin/bash

for FTYPE in chk plt_cnt part
do
    LASTFILE=$(ls ./data/*_hdf5_${FTYPE}_???? | sort | tail -1)
    FNAME=${LASTFILE:7:-4}
    NMAX=${LASTFILE: -4}
    
    for i in $(seq -w 0000 $NMAX)
    do
        if [ ! -f ./data/${FNAME}${i} ]
        then
            echo "${FNAME}${i}"
        fi
    done
done
