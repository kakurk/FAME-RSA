function run_rsa(subject, ROI, targetDSMID)
% ROI-based MVPA analysis for a single subject
%
% Load single-trial beta images from each subject, apply ROI mask, calculate 
% correlations between all trials, correlate neural dissimilarity matrix
% with hypothesized target dissimilarity matrices
%
% Written by Kyle Kurkela, kyleakurkela@gmail.com
% August, 2017

%% Pre analysis

% Add CoSMoMVPA to the MATLAB search path
addpath(genpath('/storage/home/kak53/CoSMoMVPA'))

% add the functions subfolder to the MATLAB search path
path = fileparts(mfilename('fullpath'));
addpath([path filesep 'functions'])

% turn cosmo warnings off
cosmo_warning('off')

%% Set analysis parameters

% Parameters:
%   rois                 = cell array of rois mask filenames. Assumes that
%                          this ROI is in the roi_path directory.
%   roi_path             = directory that holds the ROIs
%   studypath            = directory that holds the Single Trial SPM model.
roi_path             = '/gpfs/group/nad12/default/nad12/FAME8/RSA/ROIs';
study_path           = fullfile('/gpfs/group/nad12/default/nad12/FAME8/RSA/models', 'SingleTrialModel');

%% Directories

% path to save results into
parentDir   = fileparts(study_path);
out_path    = fullfile(parentDir, 'RSA_Results');

% This subject's:
%   data_path   = fullpath to this subject's Single Trial Model directory
%   spm_path    = fullpath to this subject's SPM.mat file. Note: the
%                 :beta appended to the end tells cosmo to pull the beta 
%                 information from the SPM.mat file.
%   output_path = fullpath to this subject's RSA output directory
data_path   = fullfile(study_path, subject);
spm_path    = fullfile(data_path, 'SPM.mat:beta');
output_path = fullfile(out_path, subject);

% create the output path if it doesn't already exist
if ~exist(output_path, 'dir')
    mkdir(output_path)
end

%% Pattern Similarity Matrices

% full path to ROI mask
mask_fn  = fullfile(roi_path, [ROI '.nii']);

% load beta images, utilizing cosmo_frmi_dataset's ability to extract
% infortmation from this subject's SPM.mat
ds  = cosmo_fmri_dataset(spm_path, 'mask', mask_fn);

% Create a targets field, which is required by
% cosmo_dissimilarity_matrix_measure. Each trial is a different
% target.
ds.sa.targets = (1:size(ds.samples, 1))';

% cosmo_remove_useless_data removes the NaNs from the data                                 
ds = cosmo_remove_useless_data(ds);

% cosmo check to make sure data in right format
cosmo_check_dataset(ds);

% compute correlation values between all trials, resulting
% in a nTrials x nTrials matrix, where each cell of the matrix represents
% the correlation between the voxel patterns for each pair of
% trials to do this, we are going to use cosmo's
% cosmo_dissimilarity_matrix_measure, which has some nice
% data organiziations features. NOTE: the output of this function
% is a **dissimilarity** matrix of 1-r, which has values between -2 and 2
%ds_dsm = cosmo_dissimilarity_matrix_measure(ds);

% There are two things to note about the outout of
% cosmo_dissimilarity_matrix measure:
%   1. It is a dissimilarity measure, 1-r
%   2. It is in vector form, an arbitrary data format designed to
%      save space and computation time
%
% We want to convert it back to its matrix form for nice, human
% readable visualization AND convert it back to similarity for 
% sanity's sake. The default dissimilarity measure in 
% cosmo_dissimilarity_matrix_measure is:
%   1 - r
% So, in order to get the similarity matrix we are looking for, we
% need to do (dsm - 1) * -1 AND convert it to matrix form
%rho = (cosmo_squareform(ds_dsm.samples) - 1) * -1;

%% Display Pattern Similarity
% display the resulting rho matrices

% % visualize the rho matrix using imagesc
% imagesc(rho);
% 
% % set axis labels
% %   set axis labels by figuring out the half way mark for each
% %   session
% labelPositions = [];
% for sess = 1:6
%     firstID = find(~cellfun(@isempty, regexp(ds.sa.labels, ['Sn\(' num2str(sess) '\).*'])), 1, 'first');
%     lastID  = find(~cellfun(@isempty, regexp(ds.sa.labels, ['Sn\(' num2str(sess) '\).*'])), 1, 'last');
%     curlabelPosition = firstID + ceil((lastID - firstID)/2);
%     labelPositions = horzcat(labelPositions, curlabelPosition);
% end
% set(gca, 'xtick', labelPositions, 'xticklabel', {'Sn(1)' 'Sn(2)' 'Sn(3)' 'Sn(4)' 'Sn(5)' 'Sn(6)'})
% set(gca, 'ytick', labelPositions, 'yticklabel', {'Sn(1)' 'Sn(2)' 'Sn(3)' 'Sn(4)' 'Sn(5)' 'Sn(6)'})
% 
% % title
% desc=sprintf(['Pattern Similarity among all trials for subject %s '...
%                 'in roi ''%s'''], subjectID, regexprep(roi_label, '_', ' '));
% title(desc)
% 
% % colorbar
% colorbar('EastOutside');
% 
% %% Write rho matrix to Excel
% filename = ['sub-' subjectID, '_roi-' roi_label '_psa-matix.xlsx'];
% xlswrite(fullfile(output_path, filename), rho)
% 
% %% Save the MATLAB figure
% filename = ['sub-' subjectID, '_roi-' roi_label '_psa-matrix.fig'];
% saveas(gcf, fullfile(output_path, filename))

%% Compare Neural Pattern Similarity to Hypothesized Target DSM using a searchlight

% memory trial types
tarFilt       = ~cellfun(@isempty, strfind(ds.sa.labels, 'trialtype-target'));
relLureFilt   = ~cellfun(@isempty, strfind(ds.sa.labels, 'trialtype-relatedLure'));
unrelLureFilt = ~cellfun(@isempty, strfind(ds.sa.labels, 'trialtype-unrealtedLure'));

% response types
rememberFilt  = ~cellfun(@isempty, strfind(ds.sa.labels, 'response-remember'));
familarFilt   = ~cellfun(@isempty, strfind(ds.sa.labels, 'response-familar'));
newFilt       = ~cellfun(@isempty, strfind(ds.sa.labels, 'response-new'));
nrFilt        = ~cellfun(@isempty, strfind(ds.sa.labels, 'response-nr'));

%%% Target DSMs

switch targetDSMID
    
    case 1
        % Hypothesis tested: Rec Hits and Fam Hits are represented
        % similarily AND are represented differently then all other trial
        % types
        RecHitsFilt  = rememberFilt & tarFilt;
        FamHitsFilt  = familarFilt  & tarFilt;
        combinedFilt = RecHitsFilt | FamHitsFilt;
        
    case 2
        % Hypothesis tested: Rec FAs and Fam FAs are represented
        % similarily AND are represented differently then all other trial
        % types
        RecFAsFilt   = rememberFilt & (relLureFilt | unrelLureFilt);
        FamFAsFilt   = familarFilt  & (relLureFilt | unrelLureFilt);
        combinedFilt = RecFAsFilt | FamFAsFilt;
        
    case 3
        % Hypothesis tested: Rec Hits and Rec FAs are represented
        % similarily AND are represented differently then all other trial
        % types
        RecHitsFilt  = rememberFilt & tarFilt;
        RecFAsFilt   = rememberFilt & (relLureFilt | unrelLureFilt);
        combinedFilt = RecHitsFilt | RecFAsFilt;
        
    case 4
        % Hypothesis tested: Rec Hits and Fam Hits are represented
        % similarily AND are represented differently then all other trial
        % types
        FamHitsFilt  = familarFilt  & tarFilt;
        FamFAsFilt   = familarFilt  & (relLureFilt | unrelLureFilt);
        combinedFilt = FamHitsFilt | FamFAsFilt;
        
end

% initalize target_dsm
target_dsm   = kron(combinedFilt, combinedFilt');

% exclude within run correlations
within_run_dsm = false(length(ds.sa.chunks));
for iChunk = unique(ds.sa.chunks)'
    boolean_vector = ds.sa.chunks == iChunk;
    within_run_dsm = within_run_dsm | kron(boolean_vector, boolean_vector');
end
target_dsm(within_run_dsm) = NaN;

% force diagnol to be zeros
target_dsm(logical(eye(size(target_dsm)))) = 0;

%%% Compare neural DSM to Target DSM
args.target_dsm = target_dsm;
args.type       = 'Spearman';
args.center     = 1;

ds_results = cosmo_target_dsm_corr_measure(ds, args);

%%% Save results to file
correlation = ds_results.samples;
subject     = {subject};
ROI         = {ROI};
results     = table(subject, ROI, targetDSMID, correlation);

writetable(results, fullfile(output_path, sprintf('sub-%s_roi-%s_targetdsm-%d.csv', subject{1}, ROI{1}, targetDSMID)))

end
