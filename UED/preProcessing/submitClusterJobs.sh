#!/bin/bash 

MAINDIR=/reg/neh/home/khegazy/analysis/CHD/UED
FILETORUN=preProcessing.exe

#if [ -z "$1" ]; then
#  echo "ERROR SUBMITTING JOBS!!!   Must give name of runList sub directory (i.e. UVpump)!"
#  exit
#fi

DATE="20161212"
if [ ! -z "$2" ]; then
  DATE=${2}
fi

#make clean; make
OUTPUTDIR=${MAINDIR}/preProcessing/output/
for file in runLists/runList_Date-${DATE}*
do

  OUTPUTFILENAME=${file:9:-4}
  echo ${OUTPUTFILENAME}

  bsub -q psanaq -o${OUTPUTDIR}${OUTPUTFILENAME}".log" ./${FILETORUN} ${file}

done


#cd ${OUTPUTDIR}/run
#cp /reg/neh/home/khegazy/analysis/UED/${MAINDIR}/${FILETORUN} .
#cp /reg/neh/home/khegazy/analysis/UED/${MAINDIR}/${EXECUTABLE} .
#sleep 5

#sed -i 's/CLUSTER=0/CLUSTER=1/g' ${FILETORUN}

#if [ ${1} -eq 1 ]; then
#  sed -i 's/DATAMC=0/DATAMC=1/g' ${FILETORUN}
#fi

#sleep 1
#stat ${FILETORUN}
#stat ${EXECUTABLEi}
#job=1
#until [ $job -gt $NJOBS ]; do
#  sleep 5
#  bsub -q psanaq -o "../logs/output"${job}".log" ./${FILETORUN} ${1} "job" ${job}
#  job=$((job + 1))
#done
