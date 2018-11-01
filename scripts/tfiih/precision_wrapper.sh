#!/bin/bash

num_cores=$1

resolution=$2

numsteps=$3

sed -e s/NUMCORES/$num_cores/g -e s/RESOLUTION/$resolution/g -e s/NSTEPS/$numsteps/g ~/optrep/code/scripts/tfiih/job_precision.sh > job_precision.sh
qsub job_precision.sh 

