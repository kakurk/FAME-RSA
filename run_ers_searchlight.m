function run_ers(iteration, varargin)
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

% the current subject to be run on this iteration
subject     = subjects{iteration};

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
output_path    = fullfile(glm_path, 'ERS_results');

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

% slice the retrieval dataset. We want to start by **only looking at
% targets**
ds_ret = cosmo_slice(ds_ret, ~cellfun(@isempty, strfind(ds_ret.sa.labels, 'trialtype-target'))); % only targets

% we need to match encoding with retrieval
[ds_enc.sa.targets, ds_ret.sa.targets] = match_trials(ds_enc.sa.labels, ds_ret.sa.labels);

% slice the encoding dataset. Remove all trials that did NOT have a match
% at retrieval (i.e., were not brought to retrieval).
ds_enc = cosmo_slice(ds_enc, ~isnan(ds_enc.sa.targets));

% the number of features and number of samples must be the same across both
% datasets at this point
assert(all(size(ds_enc.samples) == size(ds_ret.samples)))

% add a chunks field to both datasets
ds_enc.sa.chunks = ones(size(ds_enc.samples, 1), 1);
ds_ret.sa.chunks = repmat(2, size(ds_ret.samples, 1), 1);

% stack the datasets
ds_stacked = cosmo_stack({ds_enc ds_ret});

% cosmo_remove_useless_data removes the NaNs from the data                                 
ds_stacked = cosmo_remove_useless_data(ds_stacked);

% the number of unque targets MUST be equal to the number of targets/2
assert(length(unique(ds_stacked.sa.targets)) == length(ds_stacked.sa.targets) / 2)

% cosmo check to make sure data in right format
cosmo_check_dataset(ds_stacked);

% create a template matrix for @kyles_cosmo_ers_measure
template_matrix = diag(ones(1, 270)) / 270;

% Use cosmo_correlation_measure.
% This measure returns, by default, a split-half correlation measure
% based on the difference of mean correlations for matching and
% non-matching conditions (a la Haxby 2001).
measure=@kyles_cosmo_ers_measure;
args.template = template_matrix;

% define spherical neighborhood with radius of 3 voxels
radius=3; % voxels
nbrhood=cosmo_spherical_neighborhood(ds_stacked,'radius',radius);

% Run the searchlight with a 3 voxel radius
corr_results=cosmo_searchlight(ds_stacked, nbrhood, measure, args);

% Define output location
outfilename = sprintf('sub-%s_ers_searchight.nii', subject);
outputfile  = fullfile(output_path, outfilename);

% Store results to disc
cosmo_map2fmri(corr_results, outputfile);

%% subfunctions

function [trialmatch1, trialmatch2] = match_trials(triallabels1, triallabels2)
   % match_trials  match trials in cellstr1 with those in cellstring 2
   % and create two output arrays that match the two. `cellstr1match`
   % and `cellstr2match` are the same length as cellstr1 and cellstr2.

   trialmatch1 = NaN(length(triallabels1), 1);
   trialmatch2 = NaN(length(triallabels2) ,1);

   c = 0;
   
   for t1 = 1:length(triallabels1)
       
       % filename to search for
       iTrial_imagefilename = cellfun(@(x) x, regexp(triallabels1{t1}, '(?<=imagename-).*(?=_visualcategory)', 'match'), 'UniformOutput', false);

       % where it is in the first dataset
       triallabels1Filt = ~cellfun(@isempty, strfind(triallabels1, strcat('-', iTrial_imagefilename, '_')));
       
       % where it is in the second dataset
       triallabels2Filt = ~cellfun(@isempty, strfind(triallabels2, strcat('-', iTrial_imagefilename, '_')));
       
       % label them both. If there isn't a match in the second cellstring,
       % continue the loop
       if ~isempty(find(triallabels2Filt, 1))
            c = c + 1;
            trialmatch1(triallabels1Filt) = c;
            trialmatch2(triallabels2Filt) = c;
            assert(trialmatch1(triallabels1Filt) == trialmatch2(triallabels2Filt))
       end
       
   end


end


end
