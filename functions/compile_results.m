function compile_results()
% function for compiling the RSA results into a single excel sheet

% parameters
rois                 = { 'rHC_bilat' 'rLTG_bilat' 'rPHG_bilat' 'roccip_bilat' 'rSMA_bilat'};
rsa_results_dir      = '/gpfs/group/nad12/default/nad12/FAME8/RSA/models/003_memory/SingleTrialModel/RSA_Results';
trialtypesOfInterest = { 'RecHits' 'FamHits' 'RecFAs' 'FamFAs' };

% routine
for r = 1:length(rois)
    
    % current ROI name, as a string
    curROI = rois{r};
    
    % regular expression for selecting the averagetrialtypeRSAmatrices
    regExp = ['^[0-9]{2}.*' curROI '.*trialtypeRSAmatrix\.csv'];
    
    % selecting averagetrialtypeRSAmatrices
    RSA_matrix_filenames = cellstr(spm_select('FPList', rsa_results_dir, regExp));
    
    % reading in the averagetrialtypeRSAmatrices
    RSA_matrices         = cellfun(@dlmread, RSA_matrix_filenames, 'UniformOutput', false);
    
    % select only the lower triangular part of the matrix
    RSA_lwrTrian         = cellfun(@tril, RSA_matrices, 'UniformOutput', false);
    correlation          = cellfun(@extract_nonzero, RSA_lwrTrian, 'UniformOutput', false);  
    
    % Figure out all possible combinations of the trial types, select on
    % the unique combinations, join the combinations into a
    % TRIALTYPE-TRIALTYPE naming format, and repeat it to match the length
    % of correlation
    TrialCombination     = allcomb(trialtypesOfInterest, trialtypesOfInterest);
    TrialCombination     = TrialCombination([1:4, 6:8, 11:12, 16], :);
    TrialCombination     = strcat(TrialCombination(:,1), repmat({'-'}, size(TrialCombination,1), 1), TrialCombination(:,2));
    TrialCombination     = repmat(TrialCombination, size(correlation, 1), 1);
    
    % extract the subject IDs from the RSA_matrix_filenames
    subjectID            = regexp(RSA_matrix_filenames, '[0-9]{2}[oy][0-9]{3}', 'match');
    subjectID            = unNest_cell_array(subjectID);
    subjectID            = repelem(subjectID, 10);
    
    % Vertically concatenate all of the correlation values
    correlation          = vertcat(correlation{:});
    
    % create a table
    TABLE                = table(subjectID, TrialCombination, correlation);
    
    % save the table
    writetable(TABLE, fullfile(rsa_results_dir, ['compiled_' curROI '_trialtypeRSAmatrix.csv']))

    
end
    
%%% Subfunction

function out = extract_nonzero(in)

    out = in(in~=0);

end
  
% function out = grab_unique(in)
%     
%     % paramaters
%     totalNofCombos  = size(in, 1);
%     totalNofUniques = length(unique(in(:)));
%      
%     ind = P;
%     out = in(P);
%     
% end

end