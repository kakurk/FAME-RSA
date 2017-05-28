function contrasts(SPMmat)

% Initalize SPM; load this subjects SPM.mat
SPM    = [];
load(SPMmat)

% Initalize a matrix of zeros to store the contrast vectors
contrast_vectors = zeros(length(SPM.xX.name));

% For each regressor in this SPM model..
for r = 1:length(SPM.xX.name)
    
    % Current Regressor
    curRegressor = SPM.xX.name{r};
    
    %-- Create Contrast Vector
    
    % find where the current regressor is located
    logicFilt = strcmp(curRegressor. SPM.xX.name);
    
    % create the contrast vector for this regressor
    contrast_vectors(r, logicFilt) = 1;
    
end

% Set the contrast manager parameters
matlabbatch = set_conmanger(SPMmat, SPM.xX.name, cont_vector);

% Run through spm_jobman
spm_jobman('run', matlabbatch)

%% Subfunctions

function matlabbatch = set_conmanger(fullpath2SPM, cont_names, con_vectors)
% set_conmanager  set SPM's contrast manager parameters
%
%   matlabbatch = set_conmanger(fullpath2SPM, cont_names, con_vectors)
%
%       fullpath2SPM = string, full path to the SPM.mat file
%       cont_names   = cell of strings of the names of the contrasts
%       con_vectors  = a matrix of the contrast vectors, with each row
%                      corresponding to a contrast vector

    number_of_contrasts = size(con_vectors, 1); % rows = contrasts

    matlabbatch{1}.spm.stats.con.spmmat = {fullpath2SPM}; % path to SPM.mat    
    for c = 1:number_of_contrasts
        matlabbatch{1}.spm.stats.con.consess{c}.tcon.name    = cont_names{c};
        matlabbatch{1}.spm.stats.con.consess{c}.tcon.weights = con_vectors(c,:);
        matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none';
    end
    matlabbatch{1}.spm.stats.con.delete = 1;
    
end
       

end