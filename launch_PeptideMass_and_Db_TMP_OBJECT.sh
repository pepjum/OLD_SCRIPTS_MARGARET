#!/bin/bash

while getopts "f:" OPTION
do
	case $OPTION in
	f) sdir=$OPTARG
    ;;
	?)
	usage
	exit
	esac
done

EXP=$sdir'files_to_TRAINING_OBJECT_PeptideMass.txt'

while IFS='' read -r line || [[ -n "$line" ]]; do        # go through the file to get every line

    printf "$line\n"
    Rscript /home/margaret/data/pepe/scripts/PeptideMASS_and_DB_TRAINING_OBJECT.R $line &

done < "${EXP}"
