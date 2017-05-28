function [betaDescr] = gather_betas(SPMmat)

% Initalize SPM; load this subjects SPM.mat
SPM    = [];
load(SPMmat)

% Initalize beta logical vector
betalogicFilt = false(1, length(SPM.xX.name));

% Initalize a variable to hold the beta descriptions
betaDesc

%-- Idenfity betas

% for each session..
for s = 1:length(SPM.Sess)
    
    % This session's task regressors
    taskRegs = vertcat(SPM.Sess(s).U.name);
    
    for r = 1:length(taskRegs)
        
        % current regressor
        curRegressor = taskRegs{r};

        % find where the current regressor is located in the design
        logicFilt = regexp(SPM.xX.name, ['Sn\(' num2str(s) '\).' curRegressor '\*.*']);
        name      = regexp(SPM.xX.name, ['Sn\(' num2str(s) '\).' curRegressor '\*.*']);
        logicFilt = ~cellfun(@isempty, logicFilt);

        % merge this location with betalogicFilt
        betalogicFilt = betalogicFilt | logicFilt;

    end
    
end

%--Identified beta volume names

% grab vBetas
Vbetas = SPM.Vbeta(betalogicFilt);

% names of the identified beta images
filenames = vertcat(Vbetas.fname);

% append to create full path filenames
fp_files = cellstr(strcat(SPM.swd, filesep, filenames));

%-- Concatenate Betas

% set mattlabatch
matlabbatch = set_3D_to_4D(fp_files, 'concatenate_task_betas.nii');

% run spm_jobman
spm_jobman('run', matlabbatch)

end