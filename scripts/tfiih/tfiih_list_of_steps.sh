## Iteration 1 (Coarse-grain regions without known structure to 10 residues per bead).

# 1. Create beadmap
~/imp-clean/build/setup_environment.sh python ~/imp-clean/imp/modules/optrep/pyext/src/set_next_beadmap.py -u create -r 10 -tf 
~/optrep/input/tfiih/tfiih_topology.txt

# 2. Sample models in 10 residues/bead representation 
sh ~/optrep/code/scripts/tfiih/job_sample_wrapper.sh 10 3500 24 labq

# 3. Get good-scoring models from sampling
# (use total score to get a few thousand models)
date ; ~/imp-clean/build/setup_environment.sh python ~/imp-clean/imp/modules/optrep/pyext/src/select_good_scoring_models.py -rd ./ 
-rp run. -cl Total_Score -kl Total_Score -agl -99999.9 -aul -110 -mlt 0.0 -mut 0.0 ; date
# (takes 20 mins)

# 4. Get bead-wise sampling precision
sh ~/optrep/code/scripts/tfiih/precision_wrapper.sh 50 10 3500

# 5. Concatenate the beadmaps from different cores 
for i in `seq 0 49`; do cat bead_precisions_sub.$i >> bead_precisions_10.txt; done 

## Iteration 2 (Coarse-grain imprecise beads to 30 residues per bead).
# 6. Get updated beadmap 
~/imp-clean/build/setup_environment.sh python ~/imp-clean/imp/modules/optrep/pyext/src/set_next_beadmap.py -u update -r 30 -bmf r10/bead_map_10.txt -pf r10/bead_precisions_10.txt

# 7. Sample models
sh ~/optrep/code/scripts/tfiih/job_sample_wrapper.sh 30 3500 24 labq 

# 8. Get good-scoring models from sampling
(use total score to get a few thousand models)
date ; ~/imp-clean/build/setup_environment.sh python ~/imp-clean/imp/modules/optrep/pyext/src/select_good_scoring_models.py -rd ./ -rp run. -cl Total_Score -kl Total_Score -agl -99999.9 -aul -144 -mlt 0.0 -mut 0.0 ; date

# 9. Get bead-wise sampling precision
sh ~/optrep/code/scripts/tfiih/precision_wrapper.sh 50 30 3500

# 10. Concatenate the beadmaps from different cores 
for i in `seq 0 49`; do cat bead_precisions_sub.$i >> bead_precisions_30.txt; done 

## Iteration 3 (Coarse-grain imprecise beads to 50 residues per bead).
 ~/imp-clean/build/setup_environment.sh python ~/imp-clean/imp/modules/optrep/pyext/src/set_next_beadmap.py -u update -r 50 -bmf r30/bead_map_30.txt -pf r30/bead_precisions_30.txt

# Perform complete sampling
sh ~/optrep/code/scripts/tfiih/job_sample_wrapper.sh 50 7000 50 labq 
