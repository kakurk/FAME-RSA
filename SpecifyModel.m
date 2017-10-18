function [] = SpecifyModel()
% SpecifyModel      function designed to build multiple conditions files
%                   for later use with SPM's matlabbatch system
%
%   This function takes no input. All relevenat variables are defined 
%   within the body of this function.
%
%   Assumes that the behavioral data are organized as follows:
%
%   /StudyDir/BehavDir/s001/ConcatenatedBehavioralData.csv
%   /StudyDir/BehavDir/s002/ConcatenatedBehavioralData.csv
%   /StudyDir/BehavDir/s003/ConcatenatedBehavioralData.csv
%
%   Where ConcatenatedBehavioralData.csv is a comma seperated text file
%   where all functional runs are concatenated, with a column that
%   indicates which rows (i.e., trials) belong to which functional run.
%
%   See also EstimateModel.m

%% 
%==========================================================================
%                           User Input
%==========================================================================

% User Input Step 1: Directories

% Please specify the name of the current analysis, the directory the
% current analysis is in, and the directoy which houses the behavioral
% data.

Analysis.name             = 'SingleTrialModel';
Analysis.directory        = fullfile('/gpfs/group/nad12/default/nad12/FAME8/RSA/models', Analysis.name);
Analysis.behav.directory  = '/gpfs/group/nad12/default/nad12/FAME8/Behav';


% User Input Step 2: Subjects

% Please list the subjects to model in a 1 x N cell array.

Subjects       = { '18y404'  '20y297'  '20y415'  '20y441'  '20y455' ... 
                   '21y437'  '21y534'  '23y452'  '25y543'  '18y566' ... 
                   '20y396'  '20y439'  '20y444'  '21y299'  '21y521' ...
                   '22y422'  '23y546' }';


% User Input Step 3: Model Specifics

% - Behavioral File Regular Expression

Analysis.behav.regexp         = '.*ret.xls$';

%% Routine

% Clean up and print update to the command window

clc
fprintf('Model: %s\n\n', Analysis.name)
fprintf('Model Directory: \n')
disp(Analysis.directory)
fprintf('\n')
fprintf('Behavioral Data Directory: \n')
disp(Analysis.behav.directory)
fprintf('\n')

% for each subject...

for indexS = 1:length(Subjects)
    
    %-- Build Path to this Subject's Behavioral Data
    % This section builds a path to the current subjects behavioral
    % file using spm_select's recursive function. See spm_select for more
    % details

    curSubj.name      = Subjects{indexS};
    curSubj.behavdir  = spm_select('FPListRec', Analysis.behav.directory, 'dir', ['.*' curSubj.name '.*']);
    curSubj.behavFile = spm_select('FPListRec', curSubj.behavdir, Analysis.behav.regexp);

    %-- Read in this Subject's Behavioral Data
    % This section of the code reads in the subjects behavioral data using
    % the readtable command. See readtable for more details

    fprintf('Reading in Subject %s ''s Behav Data ...\n\n\n\n', curSubj.name)
    BehavData     = readtable(curSubj.behavFile);
    
    % Clean up variable names
    BehavData.Properties.VariableNames = regexprep(regexprep(BehavData.Properties.VariableNames, '_', ''), '^x', '');
    
    Number.OfRows = height(BehavData);

    %-- Build path to this subjects analysis directory
    % This section builds a path to this subjects analysis directory,
    % and creates that directory if it does not already exist.

    curSubj.directory = fullfile(Analysis.directory, curSubj.name);
    if isdir(curSubj.directory)
    else
        mkdir(curSubj.directory)
    end

    %-- Initalize the counter cell array
    % The counter cell array will keep track of how many trials occur in
    % each trial type in each functional run

    Number.OfRuns = max(unique(BehavData.runID));

    %-- Build the multiple conditions *.mat file for each run

    fprintf('Sorting Behavioral Data...\n\n\n\n')

    for curRun = 1:Number.OfRuns
        
        %-- Initialize the names, onsets, durations, and pmods structure arrays
        % This section inilaized the names, onsets, and durations
        % structure arrays, which will be filled in with the
        % approrpiate information in a nested for loop.

        % the raw onset column for this run divided by 1000 to put it in
        % seconds
        onsets    = num2cell(BehavData.RAWONSET(BehavData.runID == curRun)/1000)';
        
        % zero for each trial, a stick function
        number_of_trials_in_this_run = length(find(BehavData.runID == curRun));
        durations                    = num2cell(zeros(1, number_of_trials_in_this_run));
       
        % loop over all the trials in this run
        this_run_idxs = find(BehavData.runID == curRun);
        names         = cell(1, length(this_run_idxs)); % initialize
        
        for i = 1:length(this_run_idxs)
            
            % the current trial index
            iIDX = this_run_idxs(i);
            
            %%% pull variables to add to the trial file name
            
            % visual category
            visual_category = regexp(BehavData.image(iIDX), '(?<=\\)[a-z]+', 'match');
            visual_category = strtrim(visual_category{1}{:});
            
            % response
            if BehavData.response(iIDX) == 0
                response        = 'nr';
            elseif BehavData.response(iIDX) == 28
                response        = 'remember';
            elseif BehavData.response(iIDX) == 29
                response        = 'familiar';
            elseif BehavData.response(iIDX) == 30
                response        = 'new';
            end
            
            % trialtype
            if BehavData.type(iIDX) == 1
                trialtype        = 'target';
            elseif BehavData.type(iIDX) == 3
                trialtype        = 'relatedLure';
            elseif BehavData.type(iIDX) == 4
                trialtype        = 'unrelatedLure';
            end
            
            % enctype
            if BehavData.encType(iIDX) == 1
                enctype        = 'blocked';
            elseif BehavData.encType(iIDX) == 2
                enctype        = 'scrambled';
            end
            
            % informative, unique, BIDS style trial name
            names{i} = sprintf('visualcategory-%s_response-%s_trialtype-%s_enctype-%s', visual_category, response, trialtype, enctype);
        
        end
        
        %-- Save the Multiple Conditions *.mat file
        % Save the names, onsets, durations, and pmod variables in a .mat
        % file to be uploaded back into MATLAB/SPM at a later date for
        % model estimation

        matfilename = fullfile(curSubj.directory, ['Run', num2str(curRun, '%03d'), '_multiple_conditions.mat']);
        fprintf('Saving Subject %s''s Run %d multiple conditions file...\n\n\n', curSubj.name, curRun)
        pause(3)
        save(matfilename, 'names', 'onsets', 'durations');

    end
end

disp('All Finished!!')

end