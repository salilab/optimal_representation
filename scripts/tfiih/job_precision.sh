#$ -S /bin/bash
#$ -N pTFIIH_RESOLUTION_NSTEPS
#$ -o ./
#$ -e ./
#$ -r n
#$ -j n
#$ -l netappsali=5G
#$ -l h_rt=20:00:00
#$ -l arch=linux-x64
#$ -R yes
#$ -cwd
#$ -q lab.q
#$ -l hostname="i*"
#$ -t 1-NUMCORES

hostname
# load MPI and Sali modules
module load openmpi-1.6-nodlopen
module load sali-libraries

j=$(( $SGE_TASK_ID - 1 ))

export PYTHONDONTWRITEBYTECODE=1

# IMP stuff
export IMP=/netapp/sali/shruthi/imp-clean/build/setup_environment.sh

date

$IMP python -B /netapp/sali/shruthi/imp-clean/imp/modules/optrep/pyext/src/estimate_sampling_precision.py -n NUMCORES -cn $j -pl ccl1 kin28 kin28 tfb3 rad3 rad3 ssl2 ssl2 tfb1 tfb2 tfb2 tfb4 tfb4 tfb5 ssl1 ssl1 -dl ccl1_2 kin28_1 kin28_3 tfb3_2 rad3_1 rad3_3 ssl2_1 ssl2_3 tfb1_2 tfb2_1 tfb2_3 tfb4_1 tfb4_3 tfb5_2 ssl1_1 ssl1_3 -rd ./ -tf  /netapp/sali/shruthi/optrep/input/tfiih/tfiih_topology.txt -gs 3.0 -xs 1.0 -lc 15.0 -o ../bead_precisions_sub 

date
