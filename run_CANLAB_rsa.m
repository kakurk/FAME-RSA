%% ROI-based MVPA analysis for the CAN Lab
%
% Load single-trial beta images from each subject, apply ROI mask, calculate 
% correlations between all trials, take the mean within and across trial
% types of interest, and then mean across subjects.
%
% Written by Kyle Kurkela, kyleakurkela@gmail.com
% May, 2017

%% Pre analysis 

% Add CoSMoMVPA to the MATLAB search path
if isunix % if we are on Hammer, a unix system
    addpath(genpath('/gpfs/group/n/nad12/RSA/Scripts/CoSMoMVPA-master'))
else % if not on unix, assume we are on Anvil
    addpath(genpath('S:\nad12\CoSMoMVPA-master'))
end

% add the functions subfolder to the MATLAB search path
path = fileparts(mfilename('fullpath'));
addpath([path filesep 'functions'])

% turn cosmo warnings off
cosmo_warning('off')

%% Set analysis parameters

% Parameters:
%   subjects             = cell array of subject IDs.
%   rois                 = cell array of rois mask filenames. Assumes that
%                          this ROI is in the studypath directory.
%   studypath            = directory that holds the Single Trial SPM model.
%   trialtypesOfInterest = cell array of trial types of interest for this
%                          analysis.
subjects             = { '18y404' '18y566'  '20y297' };
rois                 = { 'rLTG_left' };
study_path           = '/gpfs/group/n/nad12/RSA/Analysis_ret/SingleTrialModel';
trialtypesOfInterest = { 'RecHits' 'FamHits' 'RecFAs' 'FamFAs' };

%% Routine

% path to save results into
output_path = fullfile(study_path, 'RSA_Results');

% create the output path if it doesn't already exist
if ~exist(output_path, 'dir')
    mkdir(output_path)
end

% initalizing cell arrays for z_all, rho_all, and trial_labels
z_all        = cell(1,length(rois));
rho_all      = cell(1,length(rois));
trial_labels = cell(1, length(subjects));

for ss = 1:length(subjects)
 
    % Edit the SPM.mat file to use paths here on Hammer
    if isunix % only execute if we are on a Unix system (i.e., Hammer)
        spm_changepath(fullfile(study_path, subjects{ss}, 'SPM.mat'), 'S:\nad12\FAME8', '/gpfs/group/n/nad12/RSA')
        spm_changepath(fullfile(study_path, subjects{ss}, 'SPM.mat'), '\', '/')
    end
    
    % This subject's:
    %   data_path = fullpath to this subject's Single Trial Model directory
    %   spm_path  = fullpath to this subject's SPM.mat file. Note: the
    %   :beta appended to the end tells cosmo to pull the beta information
    %   from the SPM.mat file.
    data_path  = fullfile(study_path, subjects{ss});
    spm_path   = fullfile(data_path, 'SPM.mat:beta');

    for rr = 1:length(rois)

        % current ROI label
        roi_label = rois{rr};

        % full path to ROI mask
        mask_fn  = fullfile(study_path, [roi_label '.nii']);

        % load beta images, utilizing cosmo_frmi_dataset's ability to extract
        % infortmation from this subject's SPM.mat
        ds_all  = cosmo_fmri_dataset(spm_path, 'mask', mask_fn);

        % Record the trial labels for this subject
        if isempty(trial_labels{ss})

            % labels pulled from the SPM.mat by cosmo
            trial_labels{ss} = ds_all.sa.labels;

            % Extract just the trial type label from the trial_labels list
            trial_labels{ss} = regexp(trial_labels{ss}, ' \w*\_', 'match');
            trial_labels{ss} = unNest_cell_array(trial_labels{ss});

            % Remove the space and the underscore
            trial_labels{ss} = regexp(trial_labels{ss}, '[a-zA-Z]*', 'match');
            trial_labels{ss} = unNest_cell_array(trial_labels{ss});
        end

        % cosmo_remove_useless_data removes the NaNs from the data                                 
        ds_all = cosmo_remove_useless_data(ds_all);

        % cosmo check to make sure data in right format
        cosmo_check_dataset(ds_all);

        % get the samples. The samples are the beta voxels extracted from the
        % current ROI from ALL of the beta images, stored in a nBetaImages x
        % nVoxels matrix
        all_ds_samples = ds_all.samples;

        % compute correlation values between all trials, resulting
        % in a nTrials x nTrials matrix, where each cell of the matrix represents
        % the correlation between the voxel patterns for each pair of trials
        rho = cosmo_corr(all_ds_samples');

        % Correlations are limited between -1 and +1, thus they cannot be normally
        % distributed. To make these correlations more 'normal', apply a Fisher
        % transformation and store this in a variable 'z'
        z = atanh(rho);

        %% Display Results
        % display the resulting rho matrices

        % create a new figure
        figure

        % visualize the rho matrix using imagesc
        imagesc(rho);

        % set axis labels
        %   set axis labels by figuring out the half way mark for each
        %   session
        labelPositions = [];
        for sess = 1:6
            firstID = find(~cellfun(@isempty, regexp(ds_all.sa.labels, ['Sn\(' num2str(sess) '\).*'])), 1, 'first');
            lastID  = find(~cellfun(@isempty, regexp(ds_all.sa.labels, ['Sn\(' num2str(sess) '\).*'])), 1, 'last');
            curlabelPosition = firstID + ceil((lastID - firstID)/2);
            labelPositions = horzcat(labelPositions, curlabelPosition);
        end
        set(gca, 'xtick', labelPositions, 'xticklabel', {'Sn(1)' 'Sn(2)' 'Sn(3)' 'Sn(4)' 'Sn(5)' 'Sn(6)'})
        set(gca, 'ytick', labelPositions, 'yticklabel', {'Sn(1)' 'Sn(2)' 'Sn(3)' 'Sn(4)' 'Sn(5)' 'Sn(6)'})

        % title
        desc=sprintf(['Average correlations among trials for subject %s '...
                        'in mask ''%s'''], subjects{ss}, rois{rr});
        title(desc)

        % colorbar
        colorbar('EastOutside');

        %% Write rho matrix to Excel
        filename = [subjects{ss}, '_' roi_label '_rho_matix.xlsx'];
        xlswrite(fullfile(output_path, filename), rho)

        %% Write z matrix to Excel
        filename = [subjects{ss}, '_' roi_label '_z_matrix.xlsx'];
        xlswrite(fullfile(output_path, filename), z)

        %% Store result in a cell array for later calculations
        z_all{rr}   = cat(3, z_all{rr}, z); % a nTrials x nTrials x nSubjects 3-D Matrix
        rho_all{rr} = cat(3, rho_all{rr}, rho); % a nTrials x nTrials x nSubjects 3-D Matrix
        
        %% Save the MATLAB figure
        filename = [subjects{ss}, '_' roi_label '_rho_matrix.png'];
        saveas(gcf, fullfile(output_path, filename))
        
    end

end

%% Within and Between Trial Type Correlations
% Now that we calculated a nTrial x nTrials RSA matrix for each subject, we
% want to combine trials by our trial types of interest to create a
% nTrialTypesOfInterest x nTrialTypesOfInterest RSA matrix which averages
% the within and between trial types correlations.

% A cell array that will hold the AverageRSAmatrices
AllSubjectsAverageRSAmatrix = cell(1, length(rois));

for s = 1:length(subjects)
    
    for r = 1:length(rois)

        % The current subjects rho matrix
        subjectsRhoMatrix = rho_all{r}(:,:,s);

        % Grab just the lower off-diagonal of the current rho matrix. See
        % tril documentation
        lowerTriangle      = tril(subjectsRhoMatrix);

        % set the diagnol (i.e., the identity correlations) to 0
        lowerTriangle(lowerTriangle == 1) = 0;
        
        % initalize an empty nTrialTypesOfInterest x nTrialTypesOfInterest cell matrix
        % note, the matrix is slightly bigger to accomondate row and column labels
        trialtypeRSAmatrix   = cell(length(trialtypesOfInterest) + 1);
        
        % add column and row labels
        trialtypeRSAmatrix(1,2:end)   = trialtypesOfInterest;
        trialtypeRSAmatrix(2:end, 1)  = trialtypesOfInterest;
        
        % initalize an empty vector
        vector = zeros(1, length(trialtypesOfInterest));
        
        % calculate the mean 'identity' correlations (i.e., the
        % correlations within the same trial type)
        for k = 1:length(trialtypesOfInterest)
            curTT = trialtypesOfInterest{k};
            vector(k)  = mean(extractCorrelations(lowerTriangle, trial_labels{s}, curTT, curTT));
        end
        
        % Put the within trial type mean correlations on the diagnol of an
        % empty square matrix. See diag documentation.
        diagnol = diag(vector);
        
        % initalize an empty vector for the off-diagnol correlations.
        % nchoosek determines that number of off-diagnol correlations. See
        % nchoosek documentation.
        vector = zeros(1, nchoosek(length(trialtypesOfInterest), 2));
        
        % all of the possible off-diagnol combinations. See combnk
        % documentation.
        offdiagnoal_combinations = sortrows(combnk(1:length(trialtypesOfInterest),2));
        
        % calculate the mean off-diagnol correlations (i.e., the
        % correaltions between trial types)
        for cp = 1:length(offdiagnoal_combinations)
            firstinpair  = trialtypesOfInterest{offdiagnoal_combinations(cp, 1)};
            secondinpair = trialtypesOfInterest{offdiagnoal_combinations(cp, 2)};
            vector(cp)   = mean(extractCorrelations(lowerTriangle, trial_labels{s}, firstinpair, secondinpair));
        end
        
        % put the off-diagnol mean correlations into a square matrix. See
        % squareform documenation.
        offdiagnol = squareform(vector);
        
        % Combine the diagnol and and the off-diagnol matrices
        AverageRSAmatrix = diagnol + offdiagnol;
        
        % Save this AverageRSAmatrix in a 3-D matrix for the next section
        AllSubjectsAverageRSAmatrix{rr} = cat(3, AllSubjectsAverageRSAmatrix{rr}, AverageRSAmatrix);
        
        % Add this new square correlation matrix to the empty cell array we
        % initalized further up in this section
        trialtypeRSAmatrix(2:end, 2:end) = num2cell(AverageRSAmatrix);
        
        %% Display Results
        % display the Average RSA matrices

        % create a new figure
        figure

        % visualize the rho matrix using imagesc
        imagesc(AverageRSAmatrix);

        % set axis labels
        set(gca, 'xtick', 1:length(trialtypesOfInterest), 'xticklabel', trialtypesOfInterest)
        set(gca, 'ytick', 1:length(trialtypesOfInterest), 'yticklabel', trialtypesOfInterest)

        % title
        desc=sprintf(['Average correlations among trials types for subject %s '...
                        'in mask ''%s'''], subjects{s}, rois{r});
        title(desc)

        % colorbar
        colorbar('EastOutside');
        
        %% Save Trial Type Averaged rho matrix
        filename = [subjects{s}, '_' rois{r} '_trialtypeRSAmatrix.xlsx'];
        xlswrite(fullfile(output_path, filename), AverageRSAmatrix)
        
        %% Save the MATLAB figure
        filename = [subjects{s}, '_' rois{r} '_trialtypeRSAmatrix.png'];
        saveas(gcf, fullfile(output_path, filename))
        
    end
    
end

%% Across Subject's correlations
% Now that we have averaged within trial types, now we will average across
% subjects

for r = 1:length(rois)
    
    % Calculate mean across subjects
    AcrossSubjectsTrialTypeMatrix = mean(AllSubjectsAverageRSAmatrix{rr}, 3);
    
    %% Display Results
    
    % create a new figure
    figure
    
    % visualize the rho matrix using imagesc
    imagesc(AcrossSubjectsTrialTypeMatrix);
    
    % set axis labels
    set(gca, 'xtick', 1:length(trialtypesOfInterest), 'xticklabel', trialtypesOfInterest)
    set(gca, 'ytick', 1:length(trialtypesOfInterest), 'yticklabel', trialtypesOfInterest)
    
    % title
    desc=sprintf(['Average correlations among trials types across subjects'...
                    'in mask ''%s'''], rois{r});
    title(desc)
    
    % colorbar
    colorbar('EastOutside');
    
    %% Save Group Results
    filename = [rois{r} '_averagetrialtypeRSAmatrix.xlsx'];
    xlswrite(fullfile(output_path, filename), AverageRSAmatrix)
    
    %% Save the MATLAB figure
    filename = [rois{r} '_averagetrialtypeRSAmatrix.png'];
    saveas(gcf, fullfile(output_path, filename))
    
end