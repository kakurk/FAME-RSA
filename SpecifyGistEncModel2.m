function [] = SpecifyGistEncModel2(~)
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

Analysis.name             = 'GistEncodingModel2';
Analysis.directory        = fullfile('/gpfs/group/nad12/default/nad12/FAME8/RSA/models/unsmoothed/', Analysis.name);
Analysis.behav.directory  = '/gpfs/group/nad12/default/nad12/FAME8/Behav';
 

% User Input Step 2: Subjects

% Please list the subjects to model in a 1 x N cell array.

Subjects       = { '18y404' '20y297'  '20y415'  '20y441'  '20y455' ... 
                    '21y437'  '21y534'  '23y452'  '25y543'  '18y566' ... 
                    '20y396'  '20y439'  '20y444'  '21y299'  '21y521' ...
                    '22y422'  '23y546' '67o136' '67o153' '67o178' '69o144' '69o277' '70o118' '70o316' '71o152' ...
                   '71o193' '72o164' '73o165' '75o320' '76o120' '76o162' '78o113' '79o108' ...
                   '79o117' '79o279' '80o121' '80o128' '81o125' '81o312' '83o197'}';

% User Input Step 3: Model Specifics

% - Behavioral File Regular Expression

Analysis.behav.regexp         = '.*\wencDM.xls$';
Number.OfTrialTypes           = 93;
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
    Number.OfRows = height(BehavData);
    
    % Clean up variable names
    BehavData.Properties.VariableNames = regexprep(regexprep(BehavData.Properties.VariableNames, '_', ''), '^x', '');
    
  

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
    counter       = zeros(Number.OfRuns, Number.OfTrialTypes);
   
    %-- Build the multiple conditions *.mat file for each run

    fprintf('Sorting Behavioral Data...\n\n\n\n')

    for curRun = 1:Number.OfRuns
        
        fprintf('Sorting Run Data...\n\n\n\n')
        disp(curRun)
        
        %-- Initialize the names, onsets, durations, and pmods structure arrays
        % This section inilaized the names, onsets, and durations
        % structure arrays, which will be filled in with the
        % approrpiate information in a nested for loop.

        % the raw onset column for this run divided by 1000 to put it in
        % seconds
%         Trialonsets    = num2cell(BehavData.RAWONSET(BehavData.runID == curRun)/1000)'; %#ok<*NASGU>
        
        % zero for each trial, a stick function

       
        % loop over all the trials in this run
        names     = cell(1, Number.OfTrialTypes); % initalizing TT names
        onsets    = cell(1, Number.OfTrialTypes); % initalizing TT onset vector
        durations = cell(1, Number.OfTrialTypes); % intializing TT durations vector

        run_idxs = find(BehavData.runID == curRun);
        
        
        for curTrial=  run_idxs'
            
            % the current trial index
          
            rawonset    = BehavData.RAWONSET(curTrial); 
            curcat        = BehavData.catID(curTrial);        % type
            %%% pull variables to add to the trial file name
           
            indexTT=0;
            for Categories = 1:93
                indexTT=indexTT+1;
                if curcat == Categories
                    visual_category = regexp(BehavData.image(curTrial), '(?<=\\)[a-z]+', 'match');
                    visual_category = strtrim(visual_category{1}{:});
                    
                    counter(curRun,indexTT)                     = counter(curRun,indexTT)+1;
                    names{indexTT}                              = sprintf('visual_category-%s', visual_category);
                    onsets{indexTT}(counter(curRun,indexTT))    = rawonset/1000;
                    durations{indexTT}(counter(curRun,indexTT)) = 0;
                end
                
            end
            
        end
        
        emptyFilt = cellfun(@isempty, names);
        names(emptyFilt) = [];
        onsets(emptyFilt) = [];
        durations(emptyFilt) = [];
            
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

