% Normalize the ERS results images

addpath(genpath('/gpfs/group/nad12/default/nad12/spm12'))

anatdir    = '/storage/home/kak53/scratch/Anat_enc';
rsadir     = '/gpfs/group/nad12/default/nad12/FAME8/RSA/models/unsmoothed/ERS_results_full_unsmoothed';

ya_subjects = {'18y404'  '20y297'  '20y415'  '20y441'  '20y455' ...
               '21y437'  '21y534'  '23y452'  '25y543'  '18y566' ...
               '20y396'  '20y439'  '20y444'  '21y299'  '21y521' ...
               '22y422'  '23y546'};

oa_subjects = {'67o136' '67o153' '67o178' '69o144' '69o277' '70o118' '70o316' '71o152' ...
               '71o193' '72o164' '73o165' '75o320' '76o120' '76o162' '78o113' '79o108' ...
               '79o117' '79o279' '80o121' '80o128' '81o125' '81o312' '83o197'};
           
subjects    = horzcat(ya_subjects, oa_subjects);

nrun = length(subjects); % enter the number of runs here
jobfile = {'/gpfs/group/nad12/default/nad12/FAME8/RSA/updated-scripts/normalize_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(2, nrun);
for crun = 1:nrun
    anat = spm_select('FPListRec', fullfile(anatdir, subjects{crun}), '^T1MprageAxial.img');
    if size(anat, 1) == 2
        anat = anat(1,:);
    end
    fprintf('%s\n\n', anat)
    inputs{1, crun} = cellstr(anat); % Segment: Volumes - cfg_files
    ers = spm_select('FPListRec', rsadir, sprintf('^sub-%s.*gist', subjects{crun}));
    inputs{2, crun} = cellstr(ers); % Normalise: Write: Images to Write - cfg_files
end

% remove empty

inputs = inputs(:, ~cellfun(@isempty, cellfun(@char, inputs(1,:), 'UniformOutput', false)));
jobs   = jobs(~cellfun(@isempty, cellfun(@char, inputs(1,:), 'UniformOutput', false)));

inputs = inputs(:, ~cellfun(@isempty, cellfun(@char, inputs(2,:), 'UniformOutput', false)));
jobs   = jobs(~cellfun(@isempty, cellfun(@char, inputs(2,:), 'UniformOutput', false)));

spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});
