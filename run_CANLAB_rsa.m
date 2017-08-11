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
    addpath(genpath('/gpfs/group/nad12/default/nad12/CoSMoMVPA-master'))
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
%                          this ROI is in the roi_path directory.
%   roi_path             = directory that holds the ROIs
%   studypath            = directory that holds the Single Trial SPM model.
%   trialtypesOfInterest = cell array of trial types of interest for this
%                          analysis.
subjects             = { '18y404'  '20y297'  '20y415'  '20y441'  '20y455' ... 
                         '21y437'  '21y534'  '23y452'  '25y543'  '18y566' ... 
                         '20y396'  '20y439'  '20y444'  '21y299'  '21y521' ...
                         '22y422'  '23y546'};
%subjects             = { '18y404' }; 
rois                 = { 'rHC_bilat' 'rLTG_bilat' 'rPHG_bilat' 'roccip_bilat' 'rSMA_bilat'};
study_path           = fullfile(root, 'SingleTrialModel'); % do NOT put the "/" at the end of the path name
roi_path             = '/gpfs/group/nad12/default/nad12/FAME8/RSA/ROIs';
changepath_flag      = false;
if m == 3
    trialtypesOfInterest = { 'RecHits' 'FamHits' 'RecFAs' 'FamFAs' };
elseif m == 2
    trialtypesOfInterest = { 'Target' 'RelLure' 'UnrelLure' };
elseif m == 1
    trialtypesOfInterest = {'airplane','backpack','balloons','baseballglove','beachball',...
                            'bed','beetle','bell','bench','bicycle',...
                            'binoculars','birdhouse','blender','bonsai','bouquet',...
                            'bowtie','broccoli','butterfly','cake','camera',...
                            'candles','cat','chandelier','clock','coffeemaker',...
                            'computer','cooler','couch','cow','crib',...
                            'cupcake','dartboard','dog','earrings','fan',...
                            'fish','flashlight','flipflop','float','fork',...
                            'frame','gargoyle','gift','glass','gravyboat',...
                            'grill','guitar','gun','hamburger','hat',...
                            'helmet','horse','hotairballoon','icecreamcone','jackolantern',...
                            'key','kite','knife','ladder','lamp',...
                            'leaf','lightbulb','lipstick','luggage','mask',...
                            'menorah','mittens','mixer','mobile','motorcycle',...
                            'mower','muffin','mug','pants','patioset',...
                            'pen','picnicbasket','piggybank','pillow','pinecone',...
                            'pizza','playground','pooltable','pot','purse',...
                            'ring','rockingchair','rollerskate','scissors','scooter',...
                            'shades','shed','sheep','shell','shoe',...
                            'shovel','sled','socks','spoon','stapler',...
                            'stove','stroller','suitcase','sundae','swimsuit',...
                            'table','teapot','teddybear','television','tent',...
                            'toaster','toilet','toiletscrub','toothbrush','trafficlight',...
                            'trashcan','turtle','umbrella','urn','vacuum',...
                            'wagon','watercan','wheelbarrow','whistle'};
end

%% Routine

% path to save results into
parentDir   = fileparts(study_path);
out_path    = fullfile(parentDir, 'RSA_Results');

% initalizing cell arrays for z_all, rho_all, and trial_labels
rho_all      = cell(1,length(rois));
trial_labels = cell(1, length(subjects));

% create a new figure
figure('Visible', 'off');

for ss = 1:length(subjects)
 
    % Edit the SPM.mat file to use paths here on Hammer
    if isunix && changepath_flag % only execute if we are on a Unix system (i.e., Hammer) && flag is set to true
        spm_changepath(fullfile(study_path, subjects{ss}, 'SPM.mat'), 'S:\nad12\FAME8', '/gpfs/group/nad12/default/nad12/FAME8')
        spm_changepath(fullfile(study_path, subjects{ss}, 'SPM.mat'), '\', '/')
    end
    
    % This subject's:
    %   data_path   = fullpath to this subject's Single Trial Model directory
    %   spm_path    = fullpath to this subject's SPM.mat file. Note: the
    %                 :beta appended to the end tells cosmo to pull the beta 
    %                 information from the SPM.mat file.
    %   output_path = fullpath to this subject's RSA output directory
    data_path   = fullfile(study_path, subjects{ss});
    spm_path    = fullfile(data_path, 'SPM.mat:beta');
    output_path = fullfile(out_path, subjects{ss});
    
    % create the output path if it doesn't already exist
    if ~exist(output_path, 'dir')
        mkdir(output_path)
    end

    for r = 1:length(rois)

        % current ROI label
        roi_label = rois{r};

        % full path to ROI mask
        mask_fn  = fullfile(roi_path, [roi_label '.nii']);

        % load beta images, utilizing cosmo_frmi_dataset's ability to extract
        % infortmation from this subject's SPM.mat
        ds_all  = cosmo_fmri_dataset(spm_path, 'mask', mask_fn);

        % Create a targets field, which is required by
        % cosmo_dissimilarity_matrix_measure. Each trial is a different
        % target.
        ds_all.sa.targets = (1:size(ds_all.samples, 1))';
        
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

        % compute correlation values between all trials, resulting
        % in a nTrials x nTrials matrix, where each cell of the matrix represents
        % the correlation between the voxel patterns for each pair of
        % trials
        % to do this, we are going to use cosmo's
        % cosmo_dissimilarity_matrix_measure, which has some nice
        % data organiziations features. NOTE: the output of this function
        % is a dissimilarity matrix
        ds_dsm = cosmo_dissimilarity_matrix_measure(ds_all);

        
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
        rho = (cosmo_squareform(ds_dsm.samples) - 1) * -1;
        
        % Correlations are limited between -1 and +1, thus they cannot be normally
        % distributed. To make these correlations more 'normal', apply a Fisher
        % transformation and store this in a variable 'z'
        z   = atanh(ds_dsm.samples);

        %% Display Results
        % display the resulting rho matrices

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
                        'in mask ''%s'''], subjects{ss}, regexprep(rois{r}, '_', ' '));
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
        rho_all{r} = cat(3, rho_all{r}, rho); % a nTrials x nTrials x nSubjects 3-D Matrix
        
        %% Save the MATLAB figure
        filename = [subjects{ss}, '_' roi_label '_rho_matrix.fig'];
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
    
    % This subject's
    %   output_path = fullpath to this subject's RSA output directory
    output_path = fullfile(out_path, subjects{s});
    
    for r = 1:length(rois)

        % Step 1: Create a ds_dsm.sa.targets1 and ds_dsm.sa.targets2 fields
        % that turns the numeric target representations --> labels
        
        ds_dsm.sa.target1labels = cell(length(ds_all.sa.labels), 1);
        ds_dsm.sa.target2labels = cell(length(ds_all.sa.labels), 1);
        
        for i = unique(ds_dsm.sa.targets1)'
           filter = ds_dsm.sa.targets1 == i;
           ds_dsm.sa.target1labels(filter) = ds_all.sa.labels(i);
        end
        
        for i = unique(ds_dsm.sa.targets2)'
           filter = ds_dsm.sa.targets2 == i;
           ds_dsm.sa.target2labels(filter) = ds_all.sa.labels(i);
        end
        
        % Step 2: Calculate the Within-Trial Type Mean Correlations
        
        % initalize an empty vector
        vector = zeros(1, length(trialtypesOfInterest));
        
        % calculate the mean 'identity' correlations (i.e., the
        % correlations within the same trial type)
        for k = 1:length(trialtypesOfInterest)
            curTT = trialtypesOfInterest{k};
            vector(k)  = mean(extractCorrelations(ds_dsm, curTT, curTT));
        end
        
        % Put the within trial type mean correlations on the diagnol of an
        % empty square matrix. See diag documentation.
        diagnol = diag(vector);
        
        % Step 3: Calculate the Between-Trial Type Mean Correlations
        
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
            vector(cp)   = mean(extractCorrelations(ds_dsm, firstinpair, secondinpair));
        end
        
        % put the off-diagnol mean correlations into a square matrix. See
        % squareform documenation.
        offdiagnol = squareform(vector);
        
        % Combine the diagnol and and the off-diagnol matrices
        AverageRSAmatrix = diagnol + offdiagnol;
        
        % Save this AverageRSAmatrix in a 3-D matrix for the next section
        AllSubjectsAverageRSAmatrix{r} = cat(3, AllSubjectsAverageRSAmatrix{r}, AverageRSAmatrix);
        
        % Add this new square correlation matrix to the empty cell array we
        % initalized further up in this section
        trialtypeRSAmatrix(2:end, 2:end) = num2cell(AverageRSAmatrix);
        
        %% Display Results
        % display the Average RSA matrices

        % visualize the rho matrix using imagesc
        imagesc(AverageRSAmatrix);

        % set axis labels
        set(gca, 'xtick', 1:length(trialtypesOfInterest), 'xticklabel', trialtypesOfInterest, 'XTickLabelRotation', 90)
        set(gca, 'ytick', 1:length(trialtypesOfInterest), 'yticklabel', trialtypesOfInterest)

        % title
        desc=sprintf(['Average correlations among trials types for subject %s '...
                        'in mask ''%s'''], subjects{s}, regexprep(rois{r}, '_', ' '));
        title(desc)

        % colorbar
        colorbar('EastOutside');
        
        %% Save Trial Type Averaged rho matrix
        filename = [subjects{s}, '_' rois{r} '_trialtypeRSAmatrix.xlsx'];
        xlswrite(fullfile(output_path, filename), AverageRSAmatrix)
        
        %% Save the MATLAB figure
        filename = [subjects{s}, '_' rois{r} '_trialtypeRSAmatrix.fig'];
        saveas(gcf, fullfile(output_path, filename))
        
    end
    
end

%% Across Subject's correlations
% Now that we have averaged within trial types, now we will average across
% subjects

ax_handles = zeros(length(rois), 1);
fg_handles = zeros(length(rois), 1);
col_limits = zeros(length(rois), 2);

for r = 1:length(rois)
    
    % Calculate mean across subjects
    AcrossSubjectsTrialTypeMatrix = mean(AllSubjectsAverageRSAmatrix{r}, 3);
    
    %% Display Results
    
    % create a new figure
    figure('Visible', 'off')
    
    % current axis handle
    ax_handles(r) = gca;
    fg_handles(r) = gcf;
    
    % visualize the rho matrix using imagesc
    imagesc(AcrossSubjectsTrialTypeMatrix);
    
    % set axis labels
    set(gca, 'xtick', 1:length(trialtypesOfInterest), 'xticklabel', trialtypesOfInterest, 'XTickLabelRotation', 90)
    set(gca, 'ytick', 1:length(trialtypesOfInterest), 'yticklabel', trialtypesOfInterest)
    
    % title
    desc=sprintf(['Average correlations among trials types across subjects'...
                    'in mask ''%s'''], rois{r});
    title(desc)
    
    % colorbar
    colorbar('EastOutside');
    
    % colorbar limits
    col_limits(r, :) = get(gca, 'clim');
    
    %% Save Group Results
    filename = [rois{r} '_averagetrialtypeRSAmatrix.xlsx'];
    xlswrite(fullfile(out_path, filename), AverageRSAmatrix)
    
end

% Save the MATLAB figure for each ROI, after colorbar correction

for r = 1:length(rois)
    
    set(ax_handles(r), 'clim', [min(col_limits(:)), max(col_limits(:))])
    
    filename = [rois{r} '_averagetrialtypeRSAmatrix.fig'];
    saveas(fg_handles(r), fullfile(out_path, filename))
    
end