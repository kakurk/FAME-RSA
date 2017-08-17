# FAME-RSA-mean
Home of the FAME RSA-Mean Analysis for the PSU CAN Lab

**rootdir** = `/gpfs/group/nad12/default/nad12/FAME8/RSA`

This repository creates 3 models of the FAME data:

1. A **visual_category** model
2. A **trial_type** model
3. A **memory** model

## Step 1: Model Specification

### Scripts 

- `SpecifyModel_visualcategory.m`
- `SpecifyModel_trialtype.m`
- `SpecifyModel_memory.m`

### Directories  
- `rootdir/models/001_visualcategory/FAMEret8/subID/`
- `rootdir/models/002_trialtype/FAMEret8/subID/`
- `rootdir/models/003_memory/FAMEret8/subID/`

### Files  

Run multiple_condition files (see the [spm manual](http://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf) for more information)

- `rootdir/models/001_visualcategory/FAMEret8/subID/Run#_multiple_conditions.mat`
- `rootdir/models/002_trialtype/FAMEret8/subID/Run#_multiple_conditions.mat`
- `rootdir/models/003_memory/FAMEret8/subID/Run#_multiple_conditions.mat`

## Step 2: Model Estimation

### Scripts 

- `EstimateModel.m`

### Directories  
- `rootdir/models/001_visualcategory/FAMEret8/subID/`
- `rootdir/models/002_trialtype/FAMEret8/subID/`
- `rootdir/models/003_memory/FAMEret8/subID/`

### Files Created 

1. `SPM.mat`: the file that defines the current model. See information on this data structure [here](http://people.duke.edu/~njs28/spmdatastructure.htm) and [here](http://andysbrainblog.blogspot.com/2013/10/whats-in-spmmat-file.html)
2. `beta_*.nii`: 3-d brain image for each regressor in the model
3. `job.mat`: file that defines the SPM job used to estimate the model. You can look at this file in the SPM batch editor
4. `mask.nii`: 3-d brain image that defines where in 3-d space the model was esimated
5. `ResMS.nii`: Residual error left over from the model
6. `RPV.nii`: no idea

## Step 3: Single Trial Model Creation

### Scripts 

- `generate_single_trial.m`

### Directories  
- `rootdir/models/001_visualcategory/SingleTrialModel/subID/`
- `rootdir/models/002_trialtype/SingleTrialModel/subID/`
- `rootdir/models/003_memory/SingleTrialModel/subID/`

### Files Created 

Same as Step 2: Model Estimate, but instead there is a *single regressor for each trial*

## Step 4: Run RSA

### Scripts 

- `run_rsa.m`

### Directories  
- `rootdir/models/001_visualcategory/RSA_Results/subID/`
- `rootdir/models/002_trialtype/RSA_Results/subID/`
- `rootdir/models/003_memory/RSA_Results/subID/`

### Files Created 

1. `subID_roiID_rho_matrix.csv`: the MVPA correlation matrix of every combination of trials. Each column and row represents a trial. The r-value in each cell is the correlation of the voxel pattern at this ROI between those two trials. The higher the r value, the more similar the pattern is.
2. `subID_roiID_rho_matrix.fig`: same as above, just instead a MATLAB figure for visualization. Note: the figure may be hidden be default, so you may need to [make it visible](https://www.mathworks.com/help/matlab/ref/figure-properties.html) in order to see it.
3. `subID_roiID_z_matrix.csv: same as above, just z-transformed the R values
4. `subID_roiID_trialtypeRSAmatrix.csv`: an MVPA correlation matrix that represents **the mean** correlation between trials **within each trial type of interest**
5. `subID_roiID_trialtypeRSAmatrix.fig`: same as above, just instead a MATLAB figure for visualization. Note: the figure may be hidden be default, so you may need to [make it visible](https://www.mathworks.com/help/matlab/ref/figure-properties.html) in order to see it.

## Step 5: Compile RSA

### Scripts 

- `compile_rsa.m`

### Directories  
- `rootdir/models/001_visualcategory/RSA_Results/`
- `rootdir/models/002_trialtype/RSA_Results`
- `rootdir/models/003_memory/RSA_Results`

### Files Created 

1. `subID_roiID_rho_matrix.csv`: the MVPA correlation matrix of every combination of trials. Each column and row represents a trial. The r-value in each cell is the correlation of the voxel pattern at this ROI between those two trials. The higher the r value, the more similar the pattern is.
2. `subID_roiID_rho_matrix.csv`: the MVPA correlation matrix of every combination of trials. Each column and row represents a trial. The r-value in each cell is the correlation of the voxel pattern at this ROI between those two trials. The higher the r value, the more similar the pattern is.
