printf "\n\n===================================================================\n PepSPECTRAST. An R version of SpectraST tool\n===================================================================\n\n"


day=`date +"%d/%m/%Y"`
hour=`date +"%H:%M"`
printf "Started at $day $hour \n"

sPath="`dirname \"$0\"`"
sPath="`( cd \"$sPath\" && pwd )`"

usage()
{
cat << EOF

OPTIONS:
   -f   Project containing mgfs file (required)
   -d   Spectral database
   -r   Range for mz filtering  (default: 0.5)
   -s   Spread value ("Intensity value of each peak for transfering to their neigbours peaks") (default: 0)
   -b   Bining value for peaks (default: 1)
   -p   Number of processors used by R scripts (default: 1)
EOF
}

range="0.5"
spread="0"
bining="1"
proc="1"


while getopts "f:d:r:s:b:p:" OPTION
do
	case $OPTION in
	f) sdir=$OPTARG
	;;
    d) dbase=$OPTARG
    ;;
    r) range=$OPTARG
    ;;
    s) spread=$OPTARG
    ;;
    b) bining=$OPTARG
    ;;
	p) proc=$OPTARG
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

if [[ -z $dbase ]]
then
     usage
     exit 1
fi


printf "launching searches\n"

#fi
export LANG=C
export LC_ALL=C

EXP=$MGFFOLDER'mgfFolderList.txt'
find $MGFFOLDER -mindepth 1 -maxdepth 1 -type d -print > $EXP   #descomentar luego


if [[ ! -z $sdir ]]; then

  while IFS='' read -r line || [[ -n "$line" ]]; do        # go through the file to get every line

               outnamefileR=$line
               #printf "Running Rscript..... $sPath/PepSpectraSt_SEP18.R -f $sdir/TMP/chip/$line.chrom -c $sdir/TMP/input/$line.chrom -d "$odir" -n "$normalization" -w "$decimating" -l "$length" -s "$scales" -g "$levth" -i "$threshold" -a "$area" -k "$foldchange" -o "$outnamefileR" -v "$annotfile" \n" &
               Rscript $sPath/PepSpectraSt_SEP18.R -f "$line" -d "$dbase" -r "$range" -s "$spread" -b "$bining"  &
               NPROC=$(($NPROC+1))
               if [ "$NPROC" -ge "$proc" ]; then
                 wait
                 NPROC=0
               fi

  done < "${EXP}"
fi
