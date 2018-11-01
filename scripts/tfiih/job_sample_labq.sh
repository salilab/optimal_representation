#$ -S /bin/bash
#$ -N sTFIIH_RESOLUTION_NSTEPS_rnRUNNUMBER
#$ -o /scrapp/shruthi/tfiih/OUT_ERR
#$ -e /scrapp/shruthi/tfiih/OUT_ERR
#$ -r n
#$ -j n
#$ -l netappsali=5G
#$ -l h_rt=40:00:00
#$ -l arch=linux-x64
#$ -R yes
#$ -cwd
#$ -pe ompi 6
#$ -q lab.q
##$ -p 0
##$ -l hostname="id*|iq*|ih*|io*"
#$ -l hostname="i*"

hostname
# load MPI and Sali modules
module load openmpi-1.6-nodlopen
module load sali-libraries

# IMP stuff
export IMP=/netapp/sali/shruthi/imp-clean/build/setup_environment.sh

date

mpirun -np $NSLOTS $IMP python -B /netapp/sali/shruthi/optrep/code/scripts/tfiih/sample_multires.py ../bead_map_RESOLUTION.txt NSTEPS

date
