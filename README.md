## Pan-Arctic Behavioural and Life-history Simulator, PASCAL (v3.10)
#### An improved behavioural and life-history simulation model for boreal Atlantic copepod, Calanus finmarchicus (3rd updgrade to v.3.10 | majour structural overhaul)
This repository contains the third version of PASCAL, initiated in 2021 & completed in 2023. PASCAL v.3.10 contains a species-specific, high-resolution behavioural and life-history simulation model for the boreal Atlantic copepod, _C. finmarchicus_. Here, the previous strategy-oriented modelling framework is changed to a much versatile Individual Based Modelling (IBM) framework. R codes are shifted to much efficient FORTRAN. Most biological functions are revised and updated. For the first time in a history of copepod IBM, a male-female (bi-sexual) breeding system was built and no artifical optimization is performed as a result. The model is self-adaptive and open-ended.  

The outputs from this model have already been analyzed and published once at:

_Bandara, K., Varpe, Ø., Maps, F., Ji, R., Eiane, K., & Tverberg, V. (2021). Timing of Calanus finmarchicus diapause in stochastic environments. Ecological Modelling, 460, 109739._

_Bandara, K., Varpe, Ø., Maps, F., Ji, R., Eiane, K., & Tverberg, V. (2023). Raw Data: Timing of Calanus finmarchicus diapause in stochastic environments (v.1.00.000) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.7646734_

##### Contents
This repository contains several main files: (i) the PASCAL v2.00 hull files for the 3 species and (ii) PASCAL v2.00 simulation drive for all species, (iii) model environment generation files and (iv). two functions to generate body size trajectories for simulated evolution. The hull file is the main file from which all other functions (model environment, size trajectories, stage-specific behavioural & life-history simulations) are called. Similar to v.1.00, the additional functions file contains three (3) small functions to be loaded besides the two main files.

##### Warnings
The environment files (2D arrays) are needed to execute the model. These files are not included here. Check the corresponding Zenodo repository for example files.
The HPC framework has been redarcted due to an ongoing publication. We will add this as soon as we can in a separate repository (this README file will be updated when the HPC framework is made Open Access)

###### Runtime
The model is set-up to run for 400 generations (calendar years). With the HPC framework built-in, it takes about 42-48 hours to complete execution in a CORE i9 7920X 24-Thread, 2.9-3.5 GHz gaming rig. If you have the PushoverR app, the model will send an update to your phone once the simulation has been completed. For single threaded operations (no HPC), please dial down the no. of individuals in the model. However, the impact of the simulation population size on the emergent ecological properties has not been assessed. 

##### Lisence
CC BY 4.0
