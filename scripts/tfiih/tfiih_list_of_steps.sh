
~/imp-clean/build/setup_environment.sh python ~/imp-clean/imp/modules/optrep/pyext/src/set_next_beadmap.py -u create -r 10 -tf ~/optrep/input/tfiih/tfiih_topology.txt

Sample initially

For full runs:
7000 steps at 50 runs

sh ~/optrep/code/scripts/tfiih/job_sample_wrapper.sh 30 7000 50 labq

For quick run (to get ropt):
3500 steps at 24 runs

sh ~/optrep/code/scripts/tfiih/job_sample_wrapper.sh 10 3500 24 labq

Trying lab.q: it takes a couple of hours! 

3. Get gsms:
(use total score to get a few thousand models)
date ; ~/imp-clean/build/setup_environment.sh python ~/imp-clean/imp/modules/optrep/pyext/src/select_good_scoring_models.py -rd ./ -rp run. -cl Total_Score -kl Total_Score -agl -99999.9 -aul -110 -mlt 0.0 -mut 0.0 ; date

(takes 20 mins)

4. Get sampling precision of bm domains 
Job_precision.sh 
Which contains only domains subject to change (beadmap)

Cat … 
5. Get next level beadmap (r=?)
~/imp-clean/build/setup_environment.sh python ~/imp-clean/imp/modules/optrep/pyext/src/set_next_beadmap.py -u update -r 30 -bmf r10/bead_map_10.txt -pf r10/bead_precisions_10.txt

6. Sample again
sh ~/optrep/code/scripts/tfiih/job_sample_wrapper.sh 30 3500 24 labq 

7. Get gsms:
(use total score to get a few thousand models)
date ; ~/imp-clean/build/setup_environment.sh python ~/imp-clean/imp/modules/optrep/pyext/src/select_good_scoring_models.py -rd ./ -rp run. -cl Total_Score -kl Total_Score -agl -99999.9 -aul -144 -mlt 0.0 -mut 0.0 ; date

8. Job_precision.sh

9. ~/imp-clean/build/setup_environment.sh python ~/imp-clean/imp/modules/optrep/pyext/src/set_next_beadmap.py -u update -r 50 -bmf r30/bead_map_30.txt -pf r30/bead_precisions_30.txt

10. Final sampling
sh ~/optrep/code/scripts/tfiih/job_sample_wrapper.sh 50 7000 50 labq 
