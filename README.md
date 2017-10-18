# FAME-RSA
Home of the FAME RSA Analysis for the PSU CAN Lab

**rootdir** = `/gpfs/group/nad12/default/nad12/FAME8/RSA`

## Step 1: Model Specification

### Scripts

- `SpecifyModel.m`

### Directories  

- `rootdir/model/`

### Files  

Run multiple_condition files (see the [spm manual](http://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf) for more information)

- `rootdir/model/subID/Run#_multiple_conditions.mat`

## Step 2: Model Estimation

### Scripts

- `EstimateModel.m`

### Directories  

- `rootdir/model/`

### Files Created

1. `SPM.mat`: the file that defines the current model. See information on this data structure [here](http://people.duke.edu/~njs28/spmdatastructure.htm) and [here](http://andysbrainblog.blogspot.com/2013/10/whats-in-spmmat-file.html)
2. `beta_*.nii`: 3-d brain image for each regressor in the model
3. `job.mat`: file that defines the SPM job used to estimate the model. You can look at this file in the SPM batch editor
4. `mask.nii`: 3-d brain image that defines where in 3-d space the model was esimated
5. `ResMS.nii`: Residual error left over from the model
6. `RPV.nii`:

## Step 3: Pattern Similarity Analysis (PSA)

### Scripts

- `run_rsa.m`

### Directories  

- `rootdir/model/subID`

### Files Created

1. `subID_roiID_rho_matrix.csv`: the MVPA correlation matrix of every combination of trials. Each column and row represents a trial. The r-value in each cell is the correlation of the voxel pattern at this ROI between those two trials. The higher the r value, the more similar the pattern is.  
2. `subID_roiID_rho_matrix.fig`: same as above, just instead a MATLAB figure for visualization.  

## Step 3: Group Level Statistics

### Scripts

- `group_statistics.m`

### Directories

- `rootdir/model/subID`

### Files Created

1.
