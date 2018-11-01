#!/bin/bash

resolution=$1
nsteps=$2
nruns=$3
qtype=$4

for i in `seq 1 $nruns`
do 
		currDir='run.'$i 
		rm -r $currDir
		mkdir $currDir
		cd $currDir

		sed -e s/RESOLUTION/$resolution/g -e s/RUNNUMBER/$i/g -e s/NSTEPS/$nsteps/g /netapp/sali/shruthi/optrep/code/scripts/tfiih/job_sample_$qtype'.sh' > job_sample.sh

		qsub job_sample.sh 

		cd ../

done


