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

## ers

Subdirectory holding the Encoding-Retrieval Similarity (ERS) scripts. These scripts run the actual multivariate analyses we are interested in.  
 
`run_ers_searchlight.m` = CoSMoMVPA based script for running an "Item-Level" ERS searchlight analysis on the FAME data. The script reads in the spm_T*.nii images from each subject's Encoding and Retrieval LSA glms and calculates the mean of the correlation values **on the diagnol** of the corrlation matrix created between each trial's encoding and retrieval pattern for a specified trial type. What results is a whole brain `*.nii` file for each of (RecHits, FamHits, Miss) for each subject. These `*.nii` files are then submitted to a second level, random effects analysis.  
 
`run_ers_gist_searchlight.m` = CoSMoMVPA based script for running a "Gist-Level" ERS searchlight analysis on the FAME data. The script reads in the spm_T*.nii images from each subject's Retrieval LSA glm and Encoding "gist" glm and calculates the mean correlation between Retrieval trials and their **corresponding encoding "gist" category**. What results is a whole brain `*.nii` file for each of (RecHits, FamHits, Miss, RecRelFA, FamRelFA, RelCR) for each subject. These `*.nii` files are then submitted to a second level, random effects analysis.  
 
`run_ers_global_searchlight.m` = CoSMoMVPA based script for running an "Global-Level" ERS searchlight analysis on the FAME data. The script reads in the spm_T*.nii images from each subject's Encoding and Retrieval LSA glms and calculates the mean of the correlation values of each Retrieval trial and **all Encoding trials**. What results is a whole brain `*.nii` file for _every possible trial type_ for each subject. These `*.nii` files are then submitted to a second level, random effects analysis.  

`run_ers_roi.m` = an older script, designed to run the "Item-level" ERS searchlight on a specified ROI (instead of a full-brain searchlight).  
 
`kyles_cosmo_ers_measure.m` and `correlation_summary_measure.m` = "measure" functions @kkurkela wrote, based heavily on the default measure functions in the CoSMoMVPA software package. These functions calculate a specified mean. See the "measure" function description on the CoSMoMVPA website.
