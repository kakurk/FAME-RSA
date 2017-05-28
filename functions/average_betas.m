function average_betas(subject, studypath, cond, express)

% Initalize SPM; load this subjects SPM.mat
SPM    = [];
SPMmat = fullfile(studypath, subject, 'SPM.mat');
load(SPMmat)

% Find all of the betas for this condition using the regular expression
matches = regexp(SPM.xX.name, express);
matches = find(~cellfun('isempty', matches));

% Gather full paths to all of the matched betas into a cell array
betas_to_average = cell(1, length(matches));
count            = 0;
for k = matches
    count = count + 1;
    betas_to_average{count} = fullfile(SPM.swd, SPM.Vbeta(k).fname);
end

% Set the IMcalc parameters to average all matched betas
matlabbatch = set_imcalc(SPM, betas_to_average, cond);

% Run through spm_jobman
spm_jobman('initcfg')
spm_jobman('run', matlabbatch)


% set_imcalc subfunction

function matlabbatch = set_imcalc(SPM, betas_to_average, cond)
    matlabbatch{1}.spm.util.imcalc.input          = betas_to_average';
    matlabbatch{1}.spm.util.imcalc.output         = ['average_beta_' cond '.nii'];
    matlabbatch{1}.spm.util.imcalc.outdir         = {SPM.swd};

    for i = 1:length(betas_to_average)
        if i == 1
            expression = strcat('(i', num2str(i), '+ ');
        elseif i == length(betas_to_average)
            expression = strcat(expression, 'i', num2str(i), ')');
        else
            expression = strcat(expression, 'i', num2str(i), '+');
        end
    end

    expression = strcat(expression, '/', num2str(length(betas_to_average)));

    matlabbatch{1}.spm.util.imcalc.expression     = expression;
    matlabbatch{1}.spm.util.imcalc.var            = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx   = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask   = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype  = 4;
end
       

end