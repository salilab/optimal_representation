import IMP.isd
import IMP.atom
import IMP.pmi.tools as tools

class GaussianEMRestraint():
    """Old version of the GaussianEMRestraint (from salilab/pmi@c522e0abc),
       updated to work with current IMP."""

    def __init__(self,root_hier,map_anchors_fn,segment_anchors=None,segment_parts=None,rigid=True,resolution=None):
        #segment parts should be a list containing protein names or a tuple with protein names and residue ranges:
        # [("ABC",1,100),"CYT","CDR"]

        self.root_hier=root_hier
        self.m=self.root_hier.get_model()
        self.model_anchors=[]
        
        for seg in segment_parts:
            self.model_anchors.extend(IMP.atom.Selection(root_hier,state_index=0,molecule=seg,copy_index=0,resolution=resolution).get_selected_particles()) #TODO

        #print "segment parts",segment_parts
        #for p in self.model_anchors:
        #    print(p.get_name())

        with open(map_anchors_fn) as data:
            D = data.readlines()
        dcoords={}
        for d in D:
            d=d.strip().split('|')
            if len(d)==6: dcoords[int(d[1])] = IMP.algebra.Vector3D(float(d[2]),float(d[3]),float(d[4]))

        # parameters
        self.model_sigmas=[15.0]*len(self.model_anchors)

        self.model_sigmas=[]
        self.model_weights=[]
        for p in self.model_anchors: self.model_sigmas.append(IMP.core.XYZR(p).get_radius())
        for p in self.model_weights: self.model_weights.append(len(IMP.atom.Fragment(p).get_residue_indexes()))

        #self.model_sigmas=[float(anch.get_as_xyzr().get_radius()) for anch in self.segment_parts]
        self.model_weights=[1.0]*len(self.model_anchors)
        self.density_sigmas=[15.0]*len(segment_anchors)
        self.density_weights=[5.0]*len(segment_anchors)
        self.sigmamaxtrans=0.1
        self.sigmamin=1.
        self.sigmamax=100.0
        self.sigmainit=10.0
        self.cutoff_dist_for_container=10.0
        self.rigid=rigid
        self.segment_anchors=segment_anchors
        self.segment_parts=segment_parts
        self.tabexp=True


        self.density_anchors=[]
        for d in dcoords:
            if self.segment_anchors==[]:
                p=IMP.Particle(self.m)
                self.density_anchors.append(p)
                IMP.core.XYZR.setup_particle(p,\
                                         IMP.algebra.Sphere3D(d,\
                                         self.density_sigmas[nd]*1.5))
            else:
                if d in self.segment_anchors:
                    p=IMP.Particle(self.m)
                    self.density_anchors.append(p)
                    IMP.core.XYZR.setup_particle(p,\
                                         IMP.algebra.Sphere3D(dcoords[d],\
                                         self.density_sigmas[0]*1.5))

        for np,p in enumerate(self.model_anchors):
            self.model_sigmas[np]=IMP.core.XYZR(p).get_radius()/1.5

        self.sigmaglobal=tools.SetupNuisance(self.m,self.sigmainit,
                 self.sigmamin,self.sigmamax,True).get_particle()
        #print('setting up restraint')

        self.gaussianEM_restraint=IMP.isd.GaussianAnchorEMRestraint(
            self.model_anchors,self.model_sigmas,self.model_weights,
            self.density_anchors,self.density_sigmas,self.density_weights,
            self.sigmaglobal.get_particle(),self.cutoff_dist_for_container,
            self.rigid,self.tabexp)
        #print('done setup')
        self.rs = IMP.RestraintSet(self.m, 'GaussianEMRestraint')
        self.rs.add_restraint(self.gaussianEM_restraint)

    def set_label(self,label):
        self.label=label

    def add_to_model(self):
        IMP.pmi.tools.add_restraint_to_model(self.m, self.rs)

    def get_particles_to_sample(self):
        ps={}
        ps["Nuisances_GaussianEMRestraint_sigma_"+self.label]=([self.sigmaglobal],self.sigmamaxtrans)
        return ps

    def get_hierarchy(self):
        return self.prot

    def get_restraint_set(self):
        return self.rs

    def get_output(self):
        self.m.update()
        output={}
        score=self.rs.unprotected_evaluate(None)
        output["_TotalScore"]=str(score)
        output["GaussianEMRestraint_"+self.label]=str(self.rs.unprotected_evaluate(None))
        output["GaussianEMRestraint_sigma_"+self.label]=str(self.sigmaglobal.get_scale())
        return output
