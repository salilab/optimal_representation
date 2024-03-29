## binary_complexes
Each sample binary complex has a [folder](binary_complexes/) which contains:

- `config files`: which specify the complete set of parameters for running sampling and analysis iteratively. The parameters include the name
of the script to be used for sampling, the inputs directory where files are stored, a topology file (below), sampling parameters, scoring function 
parameters, criteria for selecting good-scoring models for analysis, the names of proteins/domains whose representation is to be optimized, 
the set of incremental bead sizes to perform coarse-graining, a file containing Monte Carlo maximum translations for each bead size, and finally 
parameters for the representation optimizer including a scale, a linear cutoff, and the grid size used to estimate bead-wise sampling 
precision.\
Note that the standard values of parameters related to the representation optimizer are discussed in the publication (Materials and Methods).\
There is one config file per experiment (one for finding the optimal representation `r*`, one for finding an approximately optimal representation `r'` for all complexes, one for 
finding an optimal representation for complex 7CEI when the component protein structures are known, and one for finding optimal representations given sparse and dense input crosslinks for complex 2IDO).

- `inputs for modeling`: crosslinks files, FASTA and PDB files. In the case of complex 1AVX, Modeller files are also present since missing residues were added to both proteins in the complex. 

- `topology file`: a file with lines in the following format for fixing the representation for each protein domain.\
`*Protein domain protein_chain fastakey start_residue end_residue pdb pdb_chain resolution color*`\
The resolution is fixed by giving a number or specified as "bm" (stands for beadmap) if it is to be optimized. Note that one of the proteins here has a fixed resolution while the representation of the other protein is optimized. 

The domain is used to separate parts with known structure (rigid bodies) from parts of the sequence with no known structure.
Missing residues in the domains with a PDB structure are not modeled. If there is a significant stretch of missing residues in such domains, it is better to represent it as a new domain whose resolution is to be optimized (resolution "bm"). 

## tfiih
The TFIIH [folder](tfiih/) contains all modeling and analysis inputs for modeling TFIIH, including component PDB files, crosslinks file, EM map, FASTA file, move sizes for different sized beads, file containing EM anchor points, and so on. See the [TFIIH repository](https://salilab.org/tfiih) and paper for details on these.\
Additionally the folder contains a topology file similar in format to the ones used for the sample binary complexes. 
