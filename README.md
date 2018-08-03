# FAME-RSA
Home of the FAME RSA Analysis for the PSU CAN Lab

**rootdir** = `/gpfs/group/nad12/default/nad12/FAME8/RSA`

## preprocess

Subdirectory holding preprocessing scripts. In order to run our multivariate analyses, we need to grab the unsmoothed, non-normalized functional images.

`wildcard_preprocess.m` = batch preprocessing script from [KyleSPMToolbox](https://github.com/kkurkela/KyleSPMToolbox). The scripts reads in the raw data that @kkurkela transfered over to his scratch directory and runs a full preprocessing pipeline on it.

`wildcard_parameters12.m` = the parameters used for preprocessing. Realignment (SPM12 defaults) --> Slicetiming --> Coregistration (SPM12 defaults) --> Segmentation (SPM12 defaults) --> Normalization (SPM12 defaults) --> Smoothing (6mm FWHM).

## glm

Subdirectory holding general linear modeling (glm) scripts. In order to run our multivariate analyses, we need to run Least Square- All (LSA) single trial models in order to get an estimate of the BOLD activation pattern elicited for each trial at encoding and retreival. 
 
`SpecifyRetrievalModel.m` = a heavily modified version of the `SpecifyModel.m` script from [KyleSPMToolbox](https://github.com/kkurkela/KyleSPMToolbox), this scripts directly creates a LSA single trial model for _Retrieval_. Each trial is input as a seperate trial type, given a descriptive name detailing all of the information for that trial. For example:  
 
`imagename-backpack(1)_visualcategory-backpack_response-remember_trialtype-target_enctype-scrambled`  

Given these trial type names, we should be able to determine which trial type this trial would fall into (A RecHit backpack from the scrambled condition).  
 
`SpecifyEncodingModel.m` = a heavily modified version of the `SpecifyModel.m` script from [KyleSPMToolbox](https://github.com/kkurkela/KyleSPMToolbox), this scripts directly creates a LSA single trial model for _Encoding_. Each trial is input as a seperate trial type, given a descriptive name detailing all of the information for that trial. For example: 
 
`imagename-backpack(1)_visualcategory-backpack_response-pleasent_gistPosition-1_enctype-scrambled`
 
Given these trial type names, we should be able to determine which trial type this trial would fall into (e.g., a backpack from the scrambled condition that was presented first and was labeled as pleasent).    
 
`SpecifyGistEncModel2.m` = a heavily modified version of the `SpecifyModel.m` script from [KyleSPMToolbox](https://github.com/kkurkela/KyleSPMToolbox), this scripts created a special "gist" model at encoding. The gist model consists of a trial type for each visualcategory, in order to estimate the "gist" neural pattern from encoding for, for example, backpacks. We are left with 93 regressors; once for each visual category in the experiment.

`EstimateModel.m` = a slighly modified version of the `EstimateModel.m` script from [KyleSPMToolbox](https://github.com/kkurkela/KyleSPMToolbox), this script takes the multiple condition files from the previous three scripts and estimates them.
