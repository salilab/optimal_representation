This repository pertains to the benchmarks for testing the method to determine the optimal coarse-grained bead representation for a system given input data, a scoring function, and a sampling scheme.

It contains input data, scripts for modeling, and results including bead models and localization probability density maps. It uses [IMP](https://integrativemodeling.org) (*Integrative Modeling Platform*).

## Prerequisites 
Note that the results were obtained using code and scripts in the [optrep] (https://github.com/salilab/optrep) module of IMP. 

## Folder structure:
1) [inputs](inputs/) : contains all the input data used for modeling such as PDB files, crosslink files, configuration files, and so on.
2) [results](results/) : contains the optimal representations obtained, as well as final representative models. 
4) [scripts](scripts/) : scripts that are *specific to the benchmark* are included in this repository. For more general scripts and code, the user is referred to the [optrep] (https://github.com/salilab/optrep) module. 

## Information
_Author(s)_: Shruthi Viswanath 

_Maintainer_: shruthivis

_Date_: October 31st, 2018 

_License_: [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)
This work is licensed under the Creative Commons Attribution-ShareAlike 4.0
International License.

_Publications_:
- S. Viswanath and A. Sali, Optimizing model representation for integrative structure determination of macromolecular assemblies, submitted. 
