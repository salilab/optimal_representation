import IMP
import RMF
import IMP.atom
import IMP.rmf
import IMP.pmi
import IMP.pmi.tools
import IMP.pmi.topology
import IMP.pmi.dof
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.basic
import IMP.pmi.io.crosslink
import IMP.pmi.restraints.crosslinking
import os,sys,math,numpy
import IMP.optrep.pmi_representation_builder

def recenter(mol):
    "recenter the system using the coordinate of a chain specified by its chain_id"
    m=mol.mdl 
    sel = IMP.atom.Selection(mol.get_hierarchy(),resolution=IMP.atom.ALL_RESOLUTIONS)
    ps = sel.get_selected_particles()

    rb=IMP.atom.create_rigid_body(ps)
    rbcoord=rb.get_coordinates()

    rot=IMP.algebra.get_identity_rotation_3d()
    tmptrans=IMP.algebra.Transformation3D(rot,rbcoord)
    trans=tmptrans.get_inverse()
    IMP.core.transform(rb,trans)
    IMP.core.RigidBody.teardown_particle(rb)
    m.remove_particle(rb.get_particle_index())
    # change this back!Now rigid member coords need to be optimized as no longer in rigid body
    for p in ps:
        IMP.core.XYZ(p).set_coordinates_are_optimized(True)

    return(mol)

def get_radius_center_residue(mol):

    sel = IMP.atom.Selection(mol.get_hierarchy(),resolution=IMP.atom.ALL_RESOLUTIONS)
    ps = sel.get_selected_particles()
    
    max_dist=0.0
    min_dist=100000.0

    for p in ps:
        crd= IMP.core.XYZ(p).get_coordinates()

        dist=IMP.algebra.get_distance(crd,IMP.algebra.Vector3D(0.,0.,0.))

        if dist<min_dist:
            min_dist=dist
            center_res=p

        if dist>max_dist:
            max_dist=dist

    return(max_dist,center_res)

def get_ligand_protein(topology_file):
    tf=open(topology_file)
    for ln in tf.readlines():
        fields=ln.strip().split()
        if "bm" in fields[7]:
            ligand_protein = fields[0] 
            return ligand_protein # there is only 1! 
    tf.close()
    
    return "" 

###################### SYSTEM SETUP #####################
# Parameters to tweak
tgt = sys.argv[1]
topology_file=sys.argv[2]
beadmap_file = sys.argv[3]

move_sizes_file=sys.argv[4]
xlink_file=sys.argv[5]
xlink_length = float(sys.argv[6])

steps_per_run = int(sys.argv[7])

ligand_max_translation=float(sys.argv[8]) #TODO needed?
 
ev_weight=float(sys.argv[9]) #TODO needed?

# Note that the receptor and ligand ordering are different from that of the previous case
#ligand_protein="A"
ligand_protein = get_ligand_protein(topology_file)

input_dir=os.path.expanduser('~')+'/optrep/input/'+tgt

fasta_file=input_dir+'/'+tgt+'.fasta'

# Move parameters
RB_MAX_TRANS = 2.0
RB_MAX_ROT = 0.1

# Setup System and add a State
mdl = IMP.Model()
s = IMP.pmi.topology.System(mdl)
st = s.create_state()

mols= IMP.optrep.pmi_representation_builder.add_representation(st,input_dir,fasta_file,topology_file,beadmap_file)
# since objects are passed by reference I dont need to return state.

###  Once you call build(), anything without representation is destroyed.
###  You can still use handles like molecule[a:b], molecule.get_atomic_residues() or molecule.get_non_atomic_residues()
###  However these functions will only return BUILT representations
root_hier = s.build()
	      
### Uncomment this for verbose output of the representation
IMP.atom.show_with_representations(root_hier)

for mol in mols:
    recenter(mol)
    
for mol in mols:
    if ligand_protein in mol.get_name():
        ligand_mol=mol
    else:
        receptor_mol=mol
        
(receptor_radius,receptor_center_residue)=get_radius_center_residue(receptor_mol)
(ligand_radius,ligand_center_residue)=get_radius_center_residue(ligand_mol)

print receptor_center_residue
print ligand_center_residue

print "radius",receptor_radius +ligand_radius+5.0

## Setup degrees of freedom
##  The DOF functions automatically select all resolutions
##  Objects passed to nonrigid_parts move with the frame but also have their own independent movers.
dof = IMP.pmi.dof.DegreesOfFreedom(mdl)
dof.create_rigid_body(ligand_mol,max_trans=RB_MAX_TRANS,max_rot=RB_MAX_ROT,name="ligand")

######################### RESTRAINTS #####################
output_objects = [] # keep a list of functions that need to be reported
display_restraints = [] # display as springs in RMF
 
# Connectivity keeps things connected along the backbone (ignores if inside same rigid body)
crs = []
for mol in mols:
    cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(mol)
    cr.add_to_model()
    output_objects.append(cr)
    crs.append(cr)

# Excluded volume - automatically more efficient due to rigid bodies
#evr = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects = mols)
evr = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects = ligand_mol,other_objects=receptor_mol)
evr.set_weight(ev_weight) #TODO
evr.add_to_model()
output_objects.append(evr)

# Crosslink restraint
# Not using the proxl database loader for now
kw = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
kw.set_protein1_key("prot1")
kw.set_protein2_key("prot2")
kw.set_residue1_key("res1")
kw.set_residue2_key("res2")
xldb = IMP.pmi.io.crosslink.CrossLinkDataBase(kw)
xldb.create_set_from_file(xlink_file)

xlr = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(root_hier=root_hier,CrossLinkDataBase=xldb,length=xlink_length,label="XLDSS",filelabel='dss',resolution=1,slope=0.05)
xlr.add_to_model()
output_objects.append(xlr)
display_restraints.append(xlr)
            
# add distance to point restraint from receptor to ligand, to keep ligand floating away #TODO do we need this?
#ligand_center_residue_index=IMP.atom.Residue(ligand_center_residue).get_index()
if IMP.atom.Residue.get_is_setup(ligand_center_residue):
    ligand_center_residue_index=IMP.atom.Residue(ligand_center_residue).get_index()
else:
    ligand_center_residue_index=IMP.atom.Fragment(ligand_center_residue).get_residue_indexes()[0]

receptor_center_residue_coords=IMP.core.XYZ(receptor_center_residue).get_coordinates()
dsr= IMP.pmi.restraints.basic.DistanceToPointRestraint(root_hier=root_hier,tuple_selection=(ligand_center_residue_index,ligand_center_residue_index,ligand_protein,0),anchor_point=receptor_center_residue_coords,radius=receptor_radius + ligand_radius + 3.0,weight=3.0)
dsr.add_to_model()
output_objects.append(dsr)

######################## SAMPLING #####################
## First shuffle the system 
IMP.pmi.tools.shuffle_configuration(ligand_mol,max_translation=ligand_max_translation) #TODO do we need this?

nframes = steps_per_run

# Run replica exchange Monte Carlo sampling
rex=IMP.pmi.macros.ReplicaExchange0(mdl,
                                    root_hier=root_hier,                          # pass the root hierarchy
                                    crosslink_restraints=display_restraints,                     # will display like XLs
                                    monte_carlo_temperature = 1.0,
                                    replica_exchange_minimum_temperature = 1.0,
                                    replica_exchange_maximum_temperature = 2.5,
                                    num_sample_rounds = 1,
                                    number_of_best_scoring_models = 0,
                                    monte_carlo_sample_objects=dof.get_movers(),  # pass MC movers
                                    global_output_directory='output/',
                                    output_objects=output_objects,
                                    monte_carlo_steps=10,
                                    number_of_frames=nframes
				                    )

rex.execute_macro()
