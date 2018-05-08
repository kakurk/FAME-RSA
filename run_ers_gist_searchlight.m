function run_ers_gist_searchlight(iteration, varargin)
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
memtrialtypes = {'RecHits' 'FamHits' 'Misses' 'RelRecFAs' 'RelFamFAs' 'RelCRs'};

% Combinations
Combos(length(subjects) * length(memtrialtypes)).subject = cell(1);
Combos(length(subjects) * length(memtrialtypes)).memtrialtype = cell(1);

sm = 0;
for s = 1:length(subjects)
    for m = 1:length(memtrialtypes)
        sm = sm + 1;
        Combos(sm).subject = subjects{s};
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
encoding_STM_path  = fullfile(glm_path, 'GistEncodingModel2', subject, 'SPM.mat:beta');

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

% slice the retrieval dataset. We want to drop unrelated lures, since they
% do not have a "gist" pattern at encoding, by definition
ds_ret = cosmo_slice(ds_ret, cellfun(@isempty, strfind(ds_ret.sa.labels, 'trialtype-unrelatedLure'))); % NOT unrelatedLure

% slice up the retrieval dataset so that it only contains the memory trial
% type of interest for this particular iteration.
switch memtrialtype
    case 'RecHits'
        filter = ~cellfun(@isempty, strfind(ds_ret.sa.labels, 'trialtype-target'));
        filter = filter & ~cellfun(@isempty, strfind(ds_ret.sa.labels, 'response-remember'));
    case 'FamHits'
        filter = ~cellfun(@isempty, strfind(ds_ret.sa.labels, 'trialtype-target'));
        filter = filter & ~cellfun(@isempty, strfind(ds_ret.sa.labels, 'response-familiar'));
    case 'Misses'
        filter = ~cellfun(@isempty, strfind(ds_ret.sa.labels, 'trialtype-target'));
        filter = filter & ~cellfun(@isempty, strfind(ds_ret.sa.labels, 'response-new'));
    case 'RelRecFAs'
        filter = ~cellfun(@isempty, strfind(ds_ret.sa.labels, 'trialtype-relatedLure'));
        filter = filter & ~cellfun(@isempty, strfind(ds_ret.sa.labels, 'response-remember'));
    case 'RelFamFAs'
        filter = ~cellfun(@isempty, strfind(ds_ret.sa.labels, 'trialtype-relatedLure'));
        filter = filter & ~cellfun(@isempty, strfind(ds_ret.sa.labels, 'response-familiar'));
    case 'RelCRs'
        filter = ~cellfun(@isempty, strfind(ds_ret.sa.labels, 'trialtype-relatedLure'));
        filter = filter & ~cellfun(@isempty, strfind(ds_ret.sa.labels, 'response-new'));
end
ds_ret = cosmo_slice(ds_ret, filter);

% add a chunks field to both datasets
ds_enc.sa.chunks = ones(size(ds_enc.samples, 1), 1);
ds_ret.sa.chunks = repmat(2, size(ds_ret.samples, 1), 1);

% match trials. Each sample is a trial pattern. We want to sort the trials
% so that they match up nicely.

    % encoding
    ds_enc.sa.visual_categories = cellfun(@(x) x{:}, regexp(ds_enc.sa.labels, '(?<=visual_category-).*(?=*)', 'match'), 'UniformOutput', false);
    [~, I]                      = sort(ds_enc.sa.visual_categories);
    ds_enc.sa.visual_categories = ds_enc.sa.visual_categories(I);
    ds_enc.sa.chunks            = ds_enc.sa.chunks(I);
    ds_enc.sa.fname             = ds_enc.sa.fname(I);
    ds_enc.sa.labels            = ds_enc.sa.labels(I);
    ds_enc.samples              = ds_enc.samples(I, :);

    % retrieval
    ds_ret.sa.visual_categories = cellfun(@(x) x{:}, regexp(ds_ret.sa.labels, '(?<=visualcategory-).*(?=_response)', 'match'), 'UniformOutput', false);
    [~, I]                      = sort(ds_ret.sa.visual_categories);
    ds_ret.sa.visual_categories = ds_ret.sa.visual_categories(I);
    ds_ret.sa.chunks            = ds_ret.sa.chunks(I);
    ds_ret.sa.fname             = ds_ret.sa.fname(I);
    ds_ret.sa.labels            = ds_ret.sa.labels(I);
    ds_ret.samples              = ds_ret.samples(I, :);

% stack the datasets
ds_stacked = cosmo_stack({ds_enc ds_ret});

% cosmo_remove_useless_data removes the NaNs from the data                                 
ds_stacked = cosmo_remove_useless_data(ds_stacked);

% cosmo check to make sure data in right format
cosmo_check_dataset(ds_stacked);

% create a template matrix for kyles measure function. See subfunction
% below.
template_matrix = create_custom_template_matrix(ds_enc, ds_ret);
template_matrix = template_matrix / length(find(template_matrix)); 

% Use kyle's correlation_summary_measure.
% This measure returns a summary value which is the sum of the weighted 
% correlation matrix. Also need to define a post correlation function
% and the type of correlation to use. See correlation_summary_measure.m
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
outfilename = sprintf('sub-%s_memtrialtype-%s_gist-searchight.nii', subject, memtrialtype);
outputfile  = fullfile(output_path, outfilename);

% Store results to disc
cosmo_map2fmri(corr_results, outputfile);

%% subfunctions

function contrast_matrix = create_custom_template_matrix(ds_enc, ds_ret)
% custom template matrix
    
    % extract the encoding labels
    enc_labels     = ds_enc.sa.labels;

    % extract the encoding categories
    enc_categories = cellfun(@(x) x{:}, regexp(enc_labels, '(?<=visual_category-).*(?=*)', 'match'), 'UniformOutput', false);

    % extract the retrieval labels
    ret_labels     = ds_ret.sa.labels;

    % extract the retrieval categories
    ret_categories = cellfun(@(x) x{:}, regexp(ret_labels, '(?<=visualcategory-).*(?=_response)', 'match'), 'UniformOutput', false);

    contrast_matrix = zeros(length(enc_labels), length(ret_labels));

    for e = 1:length(enc_categories)

        currentCategory = enc_categories{e};

        logicvector1    = strcmp(currentCategory, enc_categories);
        logicvector2    = strcmp(currentCategory, ret_categories);

        tmpmatrix       = kron(logicvector1, logicvector2');
        contrast_matrix = tmpmatrix + contrast_matrix;

    end

    figure; imagesc(contrast_matrix);

end

end