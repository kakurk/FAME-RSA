function wildcard_preprocess(iteration)
% wildcard_preprocess   a function designed to preprocess a group of
%                       subjects using a custom designed preprocessing 
%                       routine in SPM12 or SPM8
%
% This script is designed to batch multiple subjects through a preprocessing
% pipeline designed to collect ALL AVAIABLE FUNCTIONAL RUNS USING WILDCARDS,
% and PREPROCESSING THEM ALLTOGETHER, realightning the images to the first
% image of the first run of the first collected functional run (caps added 
% for empahsis).
%
% Script assumes each subject has a high resolution anatomical in its 
% 'anat' folder. For example:
%
%   /StudyDir/AnatDir/s001/t1image.nii
%
% Also assumes that each functional run for each subject is in its own
% sub-directory. For example:
%
%   /StudyDir/FuncDir/s001/run1/fmridata.nii
%   /StudyDir/FuncDir/s001/run2/fmridata.nii
%   /StudyDir/FuncDir/s001/run3/fmridata.nii
%
% See also: wildcard_parameters8, wildcard_parameters12

%% User Input

addpath(genpath('/gpfs/group/nad12/default/nad12/spm12'))

% User Input Step 1: The subjects cell array
% List the subjects to preprocess in a cell array. Note, these IDs MUST
% match the name of the subject's subdirectory
subjects = {'18y404','18y566','20y297','20y396','20y415','20y439','20y441','20y444','20y455','20y461','21y299','21y300','21y437','21y521','21y534','22y422','23y452','23y546','25y543','67o136','67o153','67o178','69o144','69o277','70o118','71o152','71o193','72o164','73o165','76o103','76o120','76o162','78o113','79o108','79o117','79o279','80o121','80o128','81o125','83o197'}; % List of Subject IDs to batch through

% User Input Step 2: The Flag
% Set the flag to 1 to look at the parameters interactively in the GUI
% and 2 to actually run the parameters through SPM 12
flag     = 2;

% User Input 3: Wildcards
% Please specify a regular expression (google regular expressions) that
% will select the run directories, the raw functional series and the raw 
% anatomical image respectively.
regularexpr.runs = 'run.*';
regularexpr.func = '^run.*\.img'; % \.nii
regularexpr.anat = '^T1.*\.img';  % \.nii

% User Input 4: Directories
% Please specify the paths to the directories that hold the functional
% images, the anatomical images, and the desired location of the psfiles.
directories.func    = '/storage/home/kak53/scratch/Func_enc';
directories.anat    = '/storage/home/kak53/scratch/Anat_enc';
directories.psfiles = '/storage/home/kak53/scratch/psfiles';
    
%% Routine

spm('defaults', 'FMRI'); % load SPM default options
spm_jobman('initcfg')    % Configure the SPM job manger

for csub = subjects(iteration) % for each subject...
    
    % create the path to this subjects' functional folder
    subject_funcfolder = fullfile(directories.func, csub{:});

    % select all of the available run folders
    runs               = cellstr(spm_select('FPList', subject_funcfolder, 'dir', regularexpr.runs));
    
    % set batch parameters
    matlabbatch        = wildcard_parameters12(runs, csub{:}, regularexpr, directories);
    
    if flag == 1
        
        spm_jobman('interactive', matlabbatch)
        pause
        
    elseif flag == 2
        
        % configure spm graphics window. Ensures a .ps file is saved during preprocessing
        spm_figure('GetWin','Graphics');
        
        % make the psfilesdir the working directory. Ensures .ps file is saved in this directory
        cd(directories.psfiles)
        
        % run preprocessing
        spm_jobman('run', matlabbatch);
        
        % Rename the ps file from "spm_CurrentDate.ps" to "SubjectID.ps"
        temp = date;
        date_rearranged = [temp(end-3:end) temp(4:6) temp(1:2)];
        movefile(['spm_' date_rearranged '.ps'],sprintf('%s.ps',csub{:}))
        
    end

end

end