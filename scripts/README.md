## binary_complexes
Scripts for running sampling are provided. 
- `sample_rigid_to_snake.py` is the sampling script when one protein's structure is known and its representation is kept fixed ("rigid"), and the representation of the other protein with unknown structure ("snake") is optimized with non-uniform resolution, flexible beads.

- `sample_rigid_to_rigid.py` is a similar script to be used when both subunit structures are known, but one protein's representation is fixed while the 
other protein's representation is optimized. 

These scripts are called by the master scripts in the `utilities` folder of the `optrep` module (`incremental_coarse_grain_local.py` or
`incremental_coarse_grain.py`). To see how to use these scripts, refer to the master scripts which use the config files to pass command line 
parameters to these sampling scripts. 

There is also a script to test the Monte Carlo acceptance values for different runs (useful for generating approximate move sizes for different
sized beads). 

## tfiih

### Sampling


### Calculating precision


### Analysis
- `sampcon.sh` is used to get the final cluster models and densities at the optimal representation. 

- `get_xlink_violations.py` is used for 
- `chimera_densities.py` is used for visualizing density files in the results.
