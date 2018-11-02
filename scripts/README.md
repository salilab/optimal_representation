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
These scripts are for performing limited sampling at every coarse-graining iteration during the representation optimization, and finally complete sampling with the optimal representation.

- `job_sample_wrapper.sh` is the master script that takes as argument the resolution (i.e. coarse-grained bead size: 10,30, or 50 in this case), number of steps to run sampling for, number of independent sampling runs, and type of job script depending on the queue used. It calls the cluster job scripts `job_sample_all.sh` or `job_sample_labq.sh` depending on the queue.

The job scripts inturn call `sample_multires.py` (which uses `restraints.py`): these two are adapted versions of the TFIIH scripts to work with PMI2.

Note that one cannot set non-uniform resolutions in PMI1, hence the older scripts had to be adapted to PMI2.

### Calculating precision
These scripts are for calculating the bead-wise sampling precision at every coarse-graining iteration during the representation optimization.

- `precision_wrapper.sh` is the master script that takes as argument the number of cores for calculating precision (usually 4), resolution (i.e. coarse-grained bead size: 10,30, or 50 in this case), and number of steps that sampling was run for (7000 here). It calls the cluster job script `job_precision.sh` which in turn calls `estimate_sampling_precision.py` from the `optrep` module. 

### Visualization
- `chimera_densities.py` is used for visualizing density files in the results.
