#!/usr/bin/env python

import IMP
import IMP.core
import IMP.algebra
import IMP.atom
import IMP.container
import os,sys,string
import IMP.rmf
import RMF
import restraints
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.crosslinking
import IMP.pmi.representation as representation
import IMP.pmi.tools as tools
import IMP.pmi.samplers as samplers
import IMP.pmi.output as output
import IMP.optrep
import IMP.optrep.pmi_representation_builder
from numpy.random import rand as nrrand
from numpy import array

def parse_move_sizes(move_size_file):
   x=[]
   y=[]

   msf = open(move_size_file,'r')
   for ln in msf.readlines():
        fields = ln.strip().split()    
        x.append(float(fields[0]))
        y.append(float(fields[1]))

   msf.close()

   return (x,y)
   
def get_centroid_anchors(anchor_residues,anchors_file):
    
    with open(anchors_file) as data:
            D = data.readlines()
    dcoords={}
    for d in D:
            d=d.strip().split('|')
            if len(d)==6: dcoords[int(d[1])] = IMP.algebra.Vector3D(float(d[2]),float(d[3]),float(d[4]))
    
    centroid = IMP.algebra.Vector3D(0.,0.,0.)
    
    for res in anchor_residues:
        centroid += dcoords[res]    
    
    centroid = centroid/float(len(anchor_residues))
    
    #print "Centroid:",centroid

    data.close()
    
    return centroid
    
        
def shuffle_configuration_no_translation(dof, bounding_box_length=300.):
    "shuffle configuration, used to restart the optimization"
    "it only works if rigid bodies were initialized"
    if len(dof.get_rigid_bodies())==0:
        print("MultipleStates: rigid bodies were not intialized")
    hbbl=bounding_box_length/2
    ub = IMP.algebra.Vector3D(-hbbl,-hbbl,-hbbl)
    lb = IMP.algebra.Vector3D( hbbl, hbbl, hbbl)
    bb = IMP.algebra.BoundingBox3D(ub, lb)
    for rb in dof.get_rigid_bodies():
        translation = (rb.get_x(), rb.get_y(), rb.get_z())
        rotation = IMP.algebra.get_random_rotation_3d()
        transformation = IMP.algebra.Transformation3D(rotation, translation)
        rb.set_reference_frame(IMP.algebra.ReferenceFrame3D(transformation))
    for fb in dof.get_flexible_beads():
        translation = IMP.algebra.get_random_vector_in(bb)
        IMP.core.XYZ(fb).set_coordinates(translation)
            
def randomize_coords(coords):
  
  randomized = 1.*(nrrand(3)-0.5)+array(coords)

  return IMP.algebra.Vector3D(randomized[0],randomized[1],randomized[2])
    
bead_map_file = sys.argv[1]

nframes = int(sys.argv[2]) #7000

input_dir = os.path.expanduser('~')+'/optrep/input/tfiih'
topology_file = os.path.join(input_dir,'tfiih_topology.txt')
fasta_file =  os.path.join(input_dir,'tfiih_human.seqs')
move_sizes_file = os.path.join(input_dir,'move_sizes.txt')

xlinkAmplitude=17  

rbmaxtrans=1.5
rbmaxrot=0.04

resolutions_list,fbmaxtrans_list=parse_move_sizes(move_sizes_file) # 1.0

outputobjects=[]

# Setup System and add a State
m = IMP.Model()
s = IMP.pmi.topology.System(m)
st = s.create_state()

mols= IMP.optrep.pmi_representation_builder.add_representation(st,input_dir,fasta_file,topology_file,bead_map_file)

root_hier = s.build()
	      
### Uncomment this for verbose output of the representation
#IMP.atom.show_with_representations(root_hier)

### output to RMF

#fname = 'test_tfiih.rmf'
#rh = RMF.create_rmf_file(fname)
#IMP.rmf.add_hierarchy(rh, root_hier)
#IMP.rmf.save_frame(rh)

######################################################
##Set up Degrees of Freedom
######################################################
dof = IMP.pmi.dof.DegreesOfFreedom(m)

mols_dict={"kinase_dimer":[],"tfb2_tfb5_dimer":[]}

for mol in mols:
    if "ccl1" in mol.get_name() or "kin28" in mol.get_name():
        mols_dict["kinase_dimer"].append(mol)
    elif "tfb2" in mol.get_name() or "tfb5" in mol.get_name():    
        mols_dict["tfb2_tfb5_dimer"].append(mol)
        
    else:
        mols_dict[mol.get_name()]=mol
        
## --- set rigid bodies
## Sub component 1: Kinase

Kinase = [0,8,14,16,24,39,45,49,51,53,62,63,70,72,73,75,86,89,93,98,101,103,110,115,117,118,132,133,137,139,145,146,154,159,164,168,170]
Tfb3 = [3,9,20,21,36,59,67,91,94,97,121,126,156,167]
Rad3 = [4,7,18,23,26,30,31,33,37,41,43,47,52,57,60,65,84,85,99,102,108,112,119,120,134,147,148,149,152,155,157,161,162,163,165,166,169,171]
Tfiihcore = [1,2,5,6,10,11,12,13,15,17,19,22,25,26,28,29,31,32,34,35,36,37,38,40,42,44,46,48,50,52,54,55,56,59,60,61,66,67,68,69,71,74,76,77,
    78,79,80,81,82,83,88,90,91,95,96,100,104,105,106,107,109,111,112,113,114,116,122,123,124,125,127,128,129,130,131,135,136,138,140,141,143,
    144,149,150,151,153,157,158,160,162,165,166,171]

emxk,emyk,emzk = get_centroid_anchors(Kinase,os.path.join(input_dir,'tfiih_apo_multifit_moved.txt'))

kinase_unstructured = [mol.get_non_atomic_residues() for mol in mols_dict["kinase_dimer"]]
rbmovers, rb = dof.create_rigid_body(mols_dict["kinase_dimer"], nonrigid_parts=kinase_unstructured,max_trans=rbmaxtrans,max_rot=rbmaxrot,nonrigid_max_trans=-1,
fb_resolutions_list=resolutions_list,fb_move_sizes_list = fbmaxtrans_list)
rb.set_coordinates(randomize_coords((emxk,emyk,emzk)))

# Sub component 2: Tfb3
emxf,emyf,emzf = get_centroid_anchors(Tfb3,os.path.join(input_dir,'tfiih_apo_multifit_moved.txt'))

rbmovers, rb = dof.create_rigid_body(mols_dict["tfb3"], nonrigid_parts=mols_dict["tfb3"].get_non_atomic_residues(),max_trans=rbmaxtrans,max_rot=rbmaxrot,nonrigid_max_trans=-1,
fb_resolutions_list=resolutions_list,fb_move_sizes_list = fbmaxtrans_list)
rb.set_coordinates(randomize_coords((emxf,emyf,emzf)))

# Sub component 2: Rad3
emxr,emyr,emzr = get_centroid_anchors(Rad3,os.path.join(input_dir,'tfiih_apo_multifit_moved.txt'))
rbmovers, rb = dof.create_rigid_body(mols_dict["rad3"], nonrigid_parts=mols_dict["rad3"].get_non_atomic_residues(),max_trans=rbmaxtrans,max_rot=rbmaxrot,nonrigid_max_trans=-1,
fb_resolutions_list=resolutions_list,fb_move_sizes_list = fbmaxtrans_list)
rb.set_coordinates(randomize_coords((emxr,emyr,emzr)))

# Sub component 3: Ssl2 + Tfiihcore
emxt,emyt,emzt = get_centroid_anchors(Tfiihcore,os.path.join(input_dir,'tfiih_apo_multifit_moved.txt'))

tfb2_tfb5_unstructured= [mol.get_non_atomic_residues() for mol in mols_dict["tfb2_tfb5_dimer"]]
rbmovers, rb = dof.create_rigid_body(mols_dict["tfb2_tfb5_dimer"], nonrigid_parts=tfb2_tfb5_unstructured,max_trans=rbmaxtrans,max_rot=rbmaxrot,nonrigid_max_trans=-1,
fb_resolutions_list=resolutions_list,fb_move_sizes_list = fbmaxtrans_list)
rb.set_coordinates(randomize_coords((emxt,emyt,emzt)))

rbmovers, rb = dof.create_rigid_body(mols_dict["tfb1"], nonrigid_parts=mols_dict["tfb1"].get_non_atomic_residues(),max_trans=rbmaxtrans,max_rot=rbmaxrot,nonrigid_max_trans=-1,
fb_resolutions_list=resolutions_list,fb_move_sizes_list = fbmaxtrans_list)
rb.set_coordinates(randomize_coords((emxt,emyt,emzt)))

rbmovers, rb = dof.create_rigid_body(mols_dict["tfb4"], nonrigid_parts=mols_dict["tfb4"].get_non_atomic_residues(),max_trans=rbmaxtrans,max_rot=rbmaxrot,nonrigid_max_trans=-1,
fb_resolutions_list=resolutions_list,fb_move_sizes_list = fbmaxtrans_list)
rb.set_coordinates(randomize_coords((emxt,emyt,emzt)))

rbmovers, rb = dof.create_rigid_body(mols_dict["ssl1"], nonrigid_parts=mols_dict["ssl1"].get_non_atomic_residues(),max_trans=rbmaxtrans,max_rot=rbmaxrot,nonrigid_max_trans=-1,
fb_resolutions_list=resolutions_list,fb_move_sizes_list = fbmaxtrans_list)
rb.set_coordinates(randomize_coords((emxt,emyt,emzt)))

rbmovers, rb = dof.create_rigid_body(mols_dict["ssl2"], nonrigid_parts=mols_dict["ssl2"].get_non_atomic_residues(),max_trans=rbmaxtrans,max_rot=rbmaxrot,nonrigid_max_trans=-1,
fb_resolutions_list=resolutions_list,fb_move_sizes_list = fbmaxtrans_list)
rb.set_coordinates(randomize_coords((emxt,emyt,emzt)))

#--- re-orient initial orientation only
shuffle_configuration_no_translation(dof)

######################################################
##restraint setup
###################
crs = []
for mol in mols:
    cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(mol,resolution=30)
    cr.add_to_model()
    outputobjects.append(cr)
    crs.append(cr)
    
#EV restraint
ev=IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects=mols, resolution=30)
ev.add_to_model()
outputobjects.append(ev)

#cross-link restraint
xl=IMP.pmi.restraints.crosslinking.SigmoidalCrossLinkMS(root_hier,  
                     os.path.join(input_dir,'Ranish_Kornberg_thiih_xlinks_human.txt'),
                     inflection=25.0,slope=5.0,amplitude=xlinkAmplitude,
                     linear_slope=0.05,resolution=1)
xl.add_to_model()
xl.set_weight(0.1)
outputobjects.append(xl)

#GaussianEMRestraint

 
# Sub component 1: Kinase
em1=restraints.GaussianEMRestraint(root_hier,  
       os.path.join(input_dir,'tfiih_apo_multifit_moved.txt'),segment_anchors=Kinase,
	   segment_parts=['ccl1','kin28'],resolution=30)
em1.set_label('kinase_em')
em1.add_to_model()
outputobjects.append(em1)

# Sub component 2: Tfb3
em2=restraints.GaussianEMRestraint(root_hier,
        os.path.join(input_dir,'tfiih_apo_multifit_moved.txt'),segment_anchors=Tfb3,
        segment_parts=['tfb3'],resolution=30)
em2.set_label('tfb3_em')
em2.add_to_model()
outputobjects.append(em2)

#Sub component 2: Rad3
em3=restraints.GaussianEMRestraint(root_hier,
        os.path.join(input_dir,'tfiih_apo_multifit_moved.txt'),segment_anchors=Rad3,
        segment_parts=['rad3'],resolution=30)
em3.set_label('rad3_em')
em3.add_to_model()
outputobjects.append(em3)

## Sub component 3:Ssl2+Tfiihcore
em4=restraints.GaussianEMRestraint(root_hier,
        os.path.join(input_dir,'tfiih_apo_multifit_moved.txt'),segment_anchors=Tfiihcore,
        segment_parts=['ssl2','tfb1','tfb2','tfb4','tfb5','ssl1'],resolution=30)
em4.set_label('tfiicore_em')
em4.add_to_model()
outputobjects.append(em4)

#####################################################
## Monte Carlo
#####################################################
print "starting sampling"

# Run replica exchange Monte Carlo sampling
rex=IMP.pmi.macros.ReplicaExchange0(m,
                                    root_hier=root_hier,                          # pass the root hierarchy
                                    monte_carlo_temperature = 1.0,
                                    replica_exchange_minimum_temperature = 1.0,
                                    replica_exchange_maximum_temperature = 2.0,
                                    num_sample_rounds = 1,
                                    number_of_best_scoring_models = 0,
                                    monte_carlo_sample_objects=dof.get_movers(),  # pass MC movers
                                    global_output_directory='output/',
                                    output_objects=outputobjects,
                                    monte_carlo_steps=10,
                                    number_of_frames=nframes
                                    )

rex.execute_macro()
