#!/usr/bin/env python

import os,sys
import subprocess

import string
import numpy

def get_models_from_file(fil):
    
    lst = []
    
    fl = open(fil,'r')
    for ln in fl.readlines():
        lst.append(int(ln.strip()))
    fl.close()
    return lst 

def get_model_identity_from_file(model_id_file):
    
    model_ids = {}
    # runid,replicaid and frameid tuple for each model
    
    fl = open(model_id_file,'r')
    for ln in fl.readlines():
        if ln.startswith('Model'):
            continue

        fields=ln.strip().split()
    
        model_ids[int(fields[0])] = (int(fields[1]),int(fields[2]),int(fields[3]))
   

    fl.close()
    
    return model_ids
    
    
class Violations(object):

    def __init__(self, threshold):

        self.violation_threshold  = threshold 
        self.violation_counts = {}  # dictionary with a key per restraint
    
    def get_number_violated_restraints(self, frame_out,rst_keyword):
        num_violated = 0

        stat_lines = frame_out.strip().split('\n')

        for ln in stat_lines:

            if not ln.startswith(rst_keyword):
                continue

            [rst,value] = ln.strip().split()
           
            if float(value) > self.violation_threshold: 
                num_violated += 1
                if rst not in self.violation_counts:
                    self.violation_counts[rst] = 1
                else:
                    self.violation_counts[rst] += 1
        return num_violated

cluster_models_file = sys.argv[1]
model_id_file = sys.argv[2]

model_indices = get_models_from_file(cluster_models_file)
model_ids = get_model_identity_from_file(model_id_file)

Analysis = Violations(35.0)
keyword = "SigmoidalCrossLinkMS_Distance_"

for mdl in model_indices:
    (run,rep,frame) = model_ids[mdl]
   
    print mdl

    stat_file_line_command = subprocess.Popen(["/home/shruthi/imp-clean/build/setup_environment.sh","python","/home/shruthi/imp-clean/imp/modules/pmi/pyext/src/process_output.py","-f","run."+str(run)+"/output/stat."+str(rep)
+".out","-n",str(frame+1)],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    frame_out,frame_err = stat_file_line_command.communicate()

    num_violated_in_model = Analysis.get_number_violated_restraints(frame_out,keyword)

num_violated_in_all_models = 0

for rst in Analysis.violation_counts:
    if Analysis.violation_counts[rst] == len(model_indices): # violated in all models
        num_violated_in_all_models+=1
        
print "Number of crosslinks violated in all models",num_violated_in_all_models

