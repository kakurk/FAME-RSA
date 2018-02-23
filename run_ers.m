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
subjects     = {'18y404','18y566','20y297','20y396','20y415','20y439','20y441','20y444','20y455','21y299','21y437','21y521','21y534','22y422','23y452','23y546','25y543'};

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
output_path    = fullfile(glm_path, 'ERS_results', subject);

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
% trials. to do this, we are going to use cosmo's
% cosmo_dissimilarity_matrix_measure, which has some nice
% data organiziations features. NOTE: the output of this function
% is a **dissimilarity** matrix of 1-r, which has values between -2 and 2.
ds_dsm = cosmo_dissimilarity_matrix_measure(ds);

% There are two things to note about the output of
% cosmo_dissimilarity_matrix_measure:
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
% need to do (dsm - 1) * -1 AND convert it to matrix form:
rho = (cosmo_squareform(ds_dsm.samples) - 1) * -1;

%% Display Pattern Similarity
% display the resulting rho matrices

% visualize the rho matrix using imagesc. NOTE: the trials are in
% **chronological order**
imagesc(rho);

% set axis labels
%   set axis labels by figuring out the half way mark for each
%   session
labelPositions = zeros(1, 6);
for sess = 1:6
    firstID = find(~cellfun(@isempty, regexp(ds.sa.labels, ['Sn\(' num2str(sess) '\).*'])), 1, 'first');
    lastID  = find(~cellfun(@isempty, regexp(ds.sa.labels, ['Sn\(' num2str(sess) '\).*'])), 1, 'last');
    curlabelPosition = firstID + ceil((lastID - firstID)/2);
    labelPositions(sess) = curlabelPosition;
end
set(gca, 'xtick', labelPositions, 'xticklabel', {'Sn(1)' 'Sn(2)' 'Sn(3)' 'Sn(4)' 'Sn(5)' 'Sn(6)'})
set(gca, 'ytick', labelPositions, 'yticklabel', {'Sn(1)' 'Sn(2)' 'Sn(3)' 'Sn(4)' 'Sn(5)' 'Sn(6)'})

% title
desc=sprintf(['Pattern Similarity among all trials for subject %s '...
                'in roi ''%s'''], subject, regexprep(ROI, '_', ' '));
title(desc)

% colorbar
colorbar('EastOutside');

% Write rho matrix to Excel
filename = ['sub-' subject, '_roi-' ROI '_psa-matix.xlsx'];
xlswrite(fullfile(output_path, filename), rho)

% Save this MATLAB figure
filename = ['sub-' subject, '_roi-' ROI '_psa-matrix.fig'];
saveas(gcf, fullfile(output_path, filename))

%% Calculate averaged trial type matrices
% Take advantage of the nice trial labeling to parse trials into trial
% types.

% memory trial types as boolean vectors
tarFilt       = ~cellfun(@isempty, strfind(ds.sa.labels, 'trialtype-target'));
relLureFilt   = ~cellfun(@isempty, strfind(ds.sa.labels, 'trialtype-relatedLure'));
% unrelLureFilt = ~cellfun(@isempty, strfind(ds.sa.labels, 'trialtype-unrealtedLure'));

% response types as boolean vectors
rememberFilt  = ~cellfun(@isempty, strfind(ds.sa.labels, 'response-remember'));
familarFilt   = ~cellfun(@isempty, strfind(ds.sa.labels, 'response-familiar'));
%newFilt       = ~cellfun(@isempty, strfind(ds.sa.labels, 'response-new'));
%nrFilt        = ~cellfun(@isempty, strfind(ds.sa.labels, 'response-nr'));

% behavior trial types as boolean vectors
RecHitsFilt  = rememberFilt & tarFilt;
FamHitsFilt  = familarFilt  & tarFilt;
RecFAsFilt   = rememberFilt & relLureFilt;
FamFAsFilt   = familarFilt  & relLureFilt;

% cell array of the behavior trial types boolean vectors
AllTrialTypes = {RecHitsFilt, FamHitsFilt, RecFAsFilt, FamFAsFilt};

% creating the average pattern similarity square matrix, see function
% below
trial_type_matrix = create_average_pattern_similarity_square_matrix(rho, AllTrialTypes);

% visualize the trial_type_matrix matrix using imagesc
imagesc(trial_type_matrix);

% set axis labels
set(gca, 'xtick', 1:4, 'xticklabel', {'RecHits' 'FamHits' 'RecFAs' 'FamFAs'})
set(gca, 'ytick', 1:4, 'yticklabel', {'RecHits' 'FamHits' 'RecFAs' 'FamFAs'})

% title
desc=sprintf(['Average Pattern Similarity among select trials types for subject %s '...
                'in roi ''%s'''], subject, regexprep(ROI, '_', ' '));
title(desc)

% colorbar
colorbar('EastOutside');

% excel spreadsheet of the means correlation values
filename = ['sub-' subject, '_roi-' ROI '_trial-type-psa-matix.xlsx'];
xlswrite(fullfile(output_path, filename), trial_type_matrix)

% matlab .fig
filename = ['sub-' subject, '_roi-' ROI '_trial-type-psa-matrix.fig'];
saveas(gcf, fullfile(output_path, filename))

%% Create a tidyverse formatted table for final statistical analysis

% Create TrialTypeCombo, a column that defines the trial-type similarity
% comparison
TrialTypeCombo = nchoosek({'RecHits' 'FamHits' 'RecFAs' 'FamFAs'}, 2);
TrialTypeCombo = strcat(TrialTypeCombo(:, 1), {'-'}, TrialTypeCombo(:, 2));
withins = strcat({'RecHits' 'FamHits' 'RecFAs' 'FamFAs'}', {'-'}, {'RecHits' 'FamHits' 'RecFAs' 'FamFAs'}');
TrialTypeCombo = vertcat(TrialTypeCombo, withins);

% Create correlation, self explanatory
trial_type_matrix_without_the_daignol = trial_type_matrix;
trial_type_matrix_without_the_daignol(logical(eye(4))) = 0;
correlation = vertcat(squareform(trial_type_matrix_without_the_daignol)', diag(trial_type_matrix));

% create subjectid and roiid columns
subjectid   = repmat({subject}, length(TrialTypeCombo), 1);
roiid       = repmat({ROI}, length(TrialTypeCombo), 1);

% create the stats table
stats_table = table(subjectid, roiid, TrialTypeCombo, correlation);

% write the stats table
filename = sprintf('sub-%s_roi-%s_statistics-table.csv', subject, ROI);
writetable(stats_table, fullfile(output_path, filename))

%% subfunctions

function matrix = create_average_pattern_similarity_square_matrix(neural_pattern, trial_type_logical_vectors)


    %% On diagnol

    on_diagnol_averages = zeros(1, length(trial_type_logical_vectors));

    for on = 1:length(trial_type_logical_vectors)

        trial_type = trial_type_logical_vectors{on};

        on_diagnol_averages(on) = average_neural_pattern_similarity(neural_pattern, trial_type, trial_type);

    end

    %% Off diagnol

    off_diagnol_possibilites = nchoosek(1:length(trial_type_logical_vectors), 2);
    off_diagnol_averages     = zeros(1, length(off_diagnol_possibilites));

    for off = 1:length(off_diagnol_possibilites)

        trial_type1 = trial_type_logical_vectors{off_diagnol_possibilites(off,1)};
        trial_type2 = trial_type_logical_vectors{off_diagnol_possibilites(off,2)};

        off_diagnol_averages(off) = average_neural_pattern_similarity(neural_pattern, trial_type1, trial_type2);

    end

    %% create matrix

    diagnol    = diag(on_diagnol_averages);
    offdiagnol = squareform(off_diagnol_averages);
    matrix     = diagnol + offdiagnol;

    %% subfunction

    function average = average_neural_pattern_similarity(neural_pattern, logicalvector1, logicalvector2)
        % average neural pattern similarity for between given trial types, denoted
        % by logical vectors

        % create a contrast matrix that weights the lower triangle of the
        % neural pattern matrix by the number of unique trial pairs.
        lowerTriangle  = tril(kron(logicalvector1, logicalvector2'), -1);
        ContrastMatrix = lowerTriangle / length(find(lowerTriangle));

        % Perform element multiplication to get a weighted matrix
        WeightedMatrix = neural_pattern .* ContrastMatrix;

        % Sum the elements of the WeightedMatrix to get the trial type
        % average
        average        = sum(WeightedMatrix(:)); 

    end
end

end
