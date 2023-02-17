## Pan-Arctic Behavioural and Life-history Simulator, PASCAL (v3.10)
#### An improved behavioural and life-history simulation model for boreal Atlantic copepod, Calanus finmarchicus (3rd updgrade to v.3.10 | majour structural overhaul)
This repository contains the third version of PASCAL, initiated in 2021 & completed in 2023. PASCAL v.3.10 contains a species-specific, high-resolution behavioural and life-history simulation model for the boreal Atlantic copepod, _C. finmarchicus_. Here, the previous strategy-oriented modelling framework is changed to a much versatile Individual Based Modelling (IBM) framework. R codes are shifted to much efficient FORTRAN. Most biological functions are revised and updated. For the first time in a history of copepod IBM, a male-female (bi-sexual) breeding system was built and no artifical optimization is performed as a result. The model is self-adaptive and open-ended.  

The outputs from this model have already been analyzed and published once at:

_Bandara, K., Varpe, Ø., Maps, F., Ji, R., Eiane, K., & Tverberg, V. (2021). Timing of Calanus finmarchicus diapause in stochastic environments. Ecological Modelling, 460, 109739._

_Bandara, K., Varpe, Ø., Maps, F., Ji, R., Eiane, K., & Tverberg, V. (2023). Raw Data: Timing of Calanus finmarchicus diapause in stochastic environments (v.1.00.000) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.7646734_

##### Contents
This repository contains several files: (i) the PASCAL v3.00 hull files for _C. finmarchicus_ (ii) PASCAL v3.00 module files, which are collective subroutines called from the hull file and (iii) model launchcode (build & execution) in GNU Linux environment. The hull file is the main file from which all other subroutines are called. 

##### Warnings
The environment files (2D arrays) are needed to execute the model. These files are not included here. Check the corresponding Zenodo repository for example files.
The HPC framework has been redarcted due to an ongoing publication. We will add this as soon as we can in a separate repository (this README file will be updated when the HPC framework is made Open Access)

###### Runtime
The model is set-up to run for 100 generations (calendar years). With the HPC framework built-in, it takes about 11-13 hours to complete execution in a CORE i9 7920X 24-Thread, 2.9-3.5 GHz gaming rig and 6-8 hours in a RYZEN 9 5950X 16-core, 4.5 GHz gaming rig. For single threaded operations (no HPC), please dial down the no. of individuals in the model (no huge impact on runtime, nonehtless). 

##### Lisence
CC BY 4.0
