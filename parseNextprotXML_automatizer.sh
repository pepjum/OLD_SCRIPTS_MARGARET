#!/bin/bash

day=`date +"%d/%m/%Y"`
hour=`date +"%H:%M"`
printf "Started at $day $hour \n"

sPath="`dirname \"$0\"`"
sPath="`( cd \"$sPath\" && pwd )`"

#printf $sPath\n

usage()
{
cat << EOF

OPTIONS:
   -f   Directory containing Sample file (required)
   -p   procesors
EOF
}

#Defaults --

while getopts "f:p:" OPTION
do
	case $OPTION in
	f) sdir=$OPTARG
	;;
    p) cores=$OPTARG
    ;;
	?)
	usage
	exit
	;;
	esac
done

if [[ -z $sdir ]]
then
     usage
     exit 1
fi

directory=$sdir"xml/"
cd $directory
#Files=$(find . -maxdepth 1 -iname "*.xml")
#printf '%s\n' "${Files[@]}" > $directory"file_xml_list.txt"

#
# #fi
export LANG=C
export LC_ALL=C

EXP=$directory"file_xml_list.txt"

while IFS='' read -r line || [[ -n "$line" ]]; do
              
              outnamefileR=$line
              printf "Running Rscript.....$sPath/parseNEXTPROT_XML_AGO20.R "$directory" "$outnamefileR"  \n" &
              Rscript $sPath/parseNEXTPROT_XML_AGO20.R "$directory" "$outnamefileR" &
              NPROC=$(($NPROC+1))
              if [ "$NPROC" -ge "$cores" ]; then
                wait
                NPROC=0
              fi
done < "${EXP}"


wait

#rm -rf $sdir/TMP

day=`date +"%d/%m/%Y"`
hour=`date +"%H:%M"`
printf "finished at $day $hour !\n"
