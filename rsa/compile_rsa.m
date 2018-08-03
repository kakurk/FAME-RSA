function compile_rsa(modeldir, trialtypesOfInterest)
% ROI-based MVPA analysis for the CAN Lab
%
% Take the mean AverageRSAmatrix across subjects.
%
% Written by Kyle Kurkela, kyleakurkela@gmail.com
% August, 2017

%% Pre analysis

% Add CoSMoMVPA to the MATLAB search path
if isunix % if we are on Hammer, a unix system
    addpath(genpath('/gpfs/group/nad12/default/nad12/CoSMoMVPA-master'))
else % if not on unix, assume we are on Anvil
    addpath(genpath('S:\nad12\CoSMoMVPA-master'))
end

addpath(genpath('/gpfs/group/nad12/default/nad12/spm12'))

% add the functions subfolder to the MATLAB search path
path = fileparts(mfilename('fullpath'));
addpath([path filesep 'functions'])

% turn cosmo warnings off
cosmo_warning('off')

%% Set analysis parameters

rois                 = { 'rHC_bilat' 'rLTG_bilat' 'rPHG_bilat' 'roccip_bilat' 'rSMA_bilat'};

%% Across Subjects correlations
% Now that we have averaged within trial types, now we will average across
% subjects

ax_handles = zeros(length(rois), 1);
fg_handles = zeros(length(rois), 1);
col_limits = zeros(length(rois), 2);

for r = 1:length(rois)
    
    files = cellstr(spm_select('FPListRec', modeldir, [ '.*[0-9]{3}_'  rois{r} '.*trialtypeRSAmatrix.*\.csv']));
    
    AverageRSAmatrices = cellfun(@csvread, files, 'UniformOutput', false);
    
    AllSubjectsAverageRSAmatrix = cat(3, AverageRSAmatrices{:});
        
    % Calculate mean across subjects
    AcrossSubjectsTrialTypeMatrix = mean(AllSubjectsAverageRSAmatrix, 3);
    
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
                    'in mask ''%s'''], regexprep(rois{r}, '_', ' '));
    title(desc)
    
    % colorbar
    colorbar('EastOutside');
    
    % colorbar limits
    col_limits(r, :) = get(gca, 'clim');
    
    %% Save Group Results
    filename = [rois{r} '_averagetrialtypeRSAmatrix.xlsx'];
    xlswrite(fullfile(modeldir, filename), AcrossSubjectsTrialTypeMatrix)
    
end

% Save the MATLAB figure for each ROI, after colorbar correction

for r = 1:length(rois)
    
    set(ax_handles(r), 'clim', [min(col_limits(:)), max(col_limits(:))])
    
    filename = [rois{r} '_averagetrialtypeRSAmatrix.fig'];
    saveas(fg_handles(r), fullfile(modeldir, filename))
    
end

end
