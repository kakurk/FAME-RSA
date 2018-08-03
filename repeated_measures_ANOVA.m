function repeated_measures_ANOVA(A)

rootdir = '/gpfs/group/nad12/default/nad12/FAME8/RSA/models/unsmoothed/ERS_results_full_unsmoothed/';

yas = {'18y404','18y566','20y297','20y396','20y415','20y439','20y441','20y444','20y455','21y299','21y437','21y521','21y534','22y422','23y546','25y543'}; % '23y452'
oas = {'67o136','67o153','67o178','71o152','71o193','72o164','73o165','76o120','76o162','78o113','79o108','79o117','79o279','80o121','80o128','83o197'};

% '69o144'  '69o277' '70o118' 
% '70o316' '75o320' '81o125','81o312'

all = horzcat(yas, oas);

switch A
	case 'ers'
		conds = [2 3 1];
	case 'gist'
		conds = [2 3 1 4 5 6];
	case 'global'
		conds = [2 3 1 X Y Z AA BB CC];
end


for a = 1:3

    if a == 1
        anovadir = fullfile(rootdir, 'ants-normal', 'ANOVAs', A, 'YAs');
        loop = yas;
    elseif a == 2
        anovadir = fullfile(rootdir, 'ants-normal', 'ANOVAs', A, 'OAs');
        loop = oas;
    else
        anovadir = fullfile(rootdir, 'ants-normal', 'ANOVAs', A, 'All');
        loop = all;
    end
    
    if ~exist(anovadir, 'dir')
        mkdir(anovadir)
    end
    
    matlabbatch{1}.spm.stats.factorial_design.dir = {anovadir};

    for s = 1:length(loop)
        filter = sprintf('^asub-%s.*%s', loop{s}, A);
        matlabbatch{1}.spm.stats.factorial_design.des.anovaw.fsubject(s).scans = cellstr(spm_select('FPList', rootdir, filter));
        matlabbatch{1}.spm.stats.factorial_design.des.anovaw.fsubject(s).conds = conds;
    end
    matlabbatch{1}.spm.stats.factorial_design.des.anovaw.dept = 1;
    matlabbatch{1}.spm.stats.factorial_design.des.anovaw.variance = 1;
    matlabbatch{1}.spm.stats.factorial_design.des.anovaw.gmsca = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.anovaw.ancova = 0;
    matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

    spm_jobman('run', matlabbatch)

end

end
