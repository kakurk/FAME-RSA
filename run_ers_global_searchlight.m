function run_ers_global_searchlight(iteration, varargin)
% ROI-based RSA analysis for a single subject and single ROI
%
% Load single-trial beta images from each subject, apply ROI mask, calculate 
% correlations between trial patterns, take the mean across trial types
%
% Written by Kyle Kurkela, kyleakurkela@gmail.com
% August, 2017

%% Pre analysis

% Add CoSMoMVPA to the MATLAB search path
addpath(genpath('/gpfs/group/nad12/default/nad12/CoSMoMVPA-master/CoSMoMVPA-master'))

% add the functions subfolder to the MATLAB search path
path = fileparts(mfilename('fullpath'));
addpath([path filesep 'functions'])

% turn cosmo warnings off
cosmo_warning('off')

% subject ids
% subjects = {'67o136','67o153','67o178','69o144','69o277','70o118','70o316','71o152','71o193','72o164','73o165','75o320','76o120','76o162','78o113','79o108','79o117','79o279','80o121','80o128','81o125','81o312','83o197'};
subjects     = {'18y404','18y566','20y297','20y396','20y415','20y439','20y441','20y444','20y455','21y299','21y437','21y521','21y534','22y422','23y452','23y546','25y543','67o136','67o153','67o178','69o144','69o277','70o118','70o316','71o152','71o193','72o164','73o165','75o320','76o120','76o162','78o113','79o108','79o117','79o279','80o121','80o128','81o125','81o312','83o197'};

% memory trial types
memtrialtypes = {'RecHits' 'FamHits' 'Misses' 'RelRecFAs' 'RelFamFAs' 'RelCRs' 'UnrelRecFAs' 'UnrelFamFAs' 'UnrelCRs'};

% all possible combinations of subjects and trial types
Combos(length(subjects) * length(memtrialtypes)).subject      = cell(1);
Combos(length(subjects) * length(memtrialtypes)).memtrialtype = cell(1);

sm = 0;
for s = 1:length(subjects)
    for m = 1:length(memtrialtypes)
        sm = sm + 1;
        Combos(sm).subject      = subjects{s};
        Combos(sm).memtrialtype = memtrialtypes{m};
    end
end

% the current subject to be run on this iteration
subject      = Combos(iteration).subject;
memtrialtype = Combos(iteration).memtrialtype;

%% Set analysis parameters

% GLM models path. Full path to the directory containing the FAME general
% linear models (GLMs).
glm_path = '/gpfs/group/nad12/default/nad12/FAME8/RSA/models';

% Encoding Single Trial Model (STM) path. Full path to this subject's 
% single trial encoding model.
encoding_STM_path  = fullfile(glm_path, 'SingleTrialEncodingModel', subject, 'SPM.mat:beta');

% Retrieval Single Trial Model (STM) path. Full path to this subject's
% single trial retrieval model.
retrieval_STM_path = fullfile(glm_path, 'SingleTrialModel', subject, 'SPM.mat:beta');

% Output path. Directory where we are going to save the results. For now,
% we will put it in `glm_path` in a subject subfolder
output_path    = fullfile(glm_path, 'ERS_results_full');

% create the output path if it doesn't already exist
if ~exist(output_path, 'dir')
    mkdir(output_path)
end

%% Pattern Similarity Matrices

% load beta images, utilizing cosmo_fmri_dataset's ability to extract
% information from the SPM.mat
if isempty(varargin)
    ds_enc  = cosmo_fmri_dataset(encoding_STM_path, 'Mask', '-all');
    ds_ret  = cosmo_fmri_dataset(retrieval_STM_path, 'Mask', '-all');
else
    ds_enc = varargin{1};
    ds_ret = varargin{2};
end

% logical filters for identifying trial types
target_Filt     = ~cellfun(@isempty, strfind(ds_ret.sa.labels, 'trialtype-target'));
rel_lure_Filt    = ~cellfun(@isempty, strfind(ds_ret.sa.labels, 'trialtype-relatedLure'));
unrel_lure_Filt = ~cellfun(@isempty, strfind(ds_ret.sa.labels, 'trialtype-unrelatedLure'));

remember_response_Filt = ~cellfun(@isempty, regexp(ds_ret.sa.labels, 'remember')); % only remember-targets
familiar_response_Filt = ~cellfun(@isempty, regexp(ds_ret.sa.labels, 'familiar')); % only familiar-targets
new_response_Filt      = ~cellfun(@isempty, regexp(ds_ret.sa.labels, 'new')); % only new-targets

% add a chunks field to both datasets
ds_enc.sa.chunks = ones(size(ds_enc.samples, 1), 1);
ds_ret.sa.chunks = repmat(2, size(ds_ret.samples, 1), 1);

% stack the datasets
ds_stacked = cosmo_stack({ds_enc ds_ret});

% cosmo_remove_useless_data removes the NaNs from the data                                 
ds_stacked = cosmo_remove_useless_data(ds_stacked);

% cosmo check to make sure data in right format
cosmo_check_dataset(ds_stacked);

% create a template matrix for @kyles_cosmo_ers_measure. We want the
% to take the mean of **entire columns** of the correlation matrices. The
% entire column of a correlation matrix corresponds to the "global" value.
template_matrix = zeros(length(ds_enc.sa.labels), length(ds_ret.sa.labels));

switch memtrialtype
    case 'RecHits'
        template_matrix(:, remember_response_Filt & target_Filt) = 1;
    case 'FamHits'
        template_matrix(:, familiar_response_Filt & target_Filt) = 1;
    case 'Misses'
        template_matrix(:, new_response_Filt & target_Filt) = 1;
    case 'RelRecFAs'
        template_matrix(:, remember_response_Filt & rel_lure_Filt) = 1;
    case 'RelFamFAs'
        template_matrix(:, familiar_response_Filt & rel_lure_Filt) = 1;
    case 'RelCRs'
        template_matrix(:, new_response_Filt & rel_lure_Filt) = 1;
    case 'UnrelRecFAs'
        template_matrix(:, remember_response_Filt & unrel_lure_Filt) = 1;
    case 'UnrelFamFAs'
        template_matrix(:, familiar_response_Filt & unrel_lure_Filt) = 1;
    case 'UnrelCRs'
        template_matrix(:, new_response_Filt & unrel_lure_Filt) = 1;
end

% weight the template matrix
template_matrix = template_matrix/length(find(template_matrix));

% Use correlation_summary_measure.
% This measure returns a single summary value of a correlation matrix,
% defined as the weighted sum of the template matrix.
measure=@correlation_summary_measure;
args.template       = template_matrix;
args.corr_type      = 'spearman';
args.post_corr_func = @atanh;

% define spherical neighborhood with radius of 3 voxels
radius=3; % voxels
nbrhood=cosmo_spherical_neighborhood(ds_stacked,'radius',radius);

% Run the searchlight with a 3 voxel radius
corr_results=cosmo_searchlight(ds_stacked, nbrhood, measure, args);

% Define output location
outfilename = sprintf('sub-%s_memtrialtype-%s_global-searchight.nii', subject, memtrialtype);
outputfile  = fullfile(output_path, outfilename);

% Store results to disc
cosmo_map2fmri(corr_results, outputfile);

end
