
import IMP
import RMF
import IMP.atom
import IMP.rmf
import IMP.pmi
import IMP.pmi.topology
import IMP.pmi.dof
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.basic
import IMP.pmi.io.crosslink
import IMP.pmi.restraints.crosslinking
import os,sys,string

statFile = sys.argv[1]

mover_type = sys.argv[2] 
# rigid or floppy 

# read the stat file
sf =open(statFile,'r')

for i,ln in enumerate(sf.readlines()):
        lndict = eval(ln.strip())
	
	if i==0: # first line
	    # get the keys corresponding to these fields
	    keysAsked = []		

	    for k in lndict:
                if mover_type == "rigid" and lndict[k].startswith("MonteCarlo_Acceptance_P"):
                    keysAsked.append(k)
                
                if mover_type == "floppy" and lndict[k].startswith("MonteCarlo_Acceptance_BallMover"):
                    keysAsked.append(k)
            
	    print "Keys",keysAsked
 
        else:
            avgAcceptance =0.0
    
            for k in keysAsked:
                avgAcceptance = avgAcceptance + float(lndict[k])

            avgAcceptance = avgAcceptance/float(len(keysAsked))
    
            print avgAcceptance	

sf.close()


