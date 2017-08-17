function compile_results(rois, rsa_results_dir, trialtypesOfInterest)
% compile_results  function for compiling the RSA results into a single 
% excel sheet
% 
% Input Arguments:
%
%   rois                 = cell array of strings of the IDs of the ROIs to 
%                          be compiled
%   rsa_results_dir      = string of the path to the RSA_Results directory
%   trialtypesOfInterest = cell array of strings, indiciating the IDs of
%                          the trial types of interest that you want to 
%                          compile
%
% Output Arguments:
%   
%   None
%
% Detailed Description:
%
%   Combinds all of the information needed to run an ANOVA on the
%   correlation coefficients created during the MVPA analyses. Creates an
%   excel sheet called "compiled_roiID_trialtypeRSAmatrix.csv" in the
%   rsa_results_dir. The columns of this excel sheet are subjectID (self

%   explanatory), TrialCombination (description of the unique trial type 
%   combination), and correlation (the MEAN correlation coefficeint "r" of
%   the individual trial patterns)
%
% Example Input Parameters:
%   rois                 = { 'rHC_bilat' 'rLTG_bilat' 'rPHG_bilat' 'roccip_bilat' 'rSMA_bilat' };
%   rsa_results_dir      = '/gpfs/group/nad12/default/nad12/FAME8/RSA/models/003_memory/RSA_Results';
%   trialtypesOfInterest = { 'RecHits' 'FamHits' 'RecFAs' 'FamFAs' };
%
% See also: compile_rsa, run_rsa, spm_select
%% Routine

% for each ROI...
for r = 1:length(rois)
    
    % current ROI name, as a string
    curROI = rois{r};
    
    % regular expression for selecting the averagetrialtypeRSAmatrices. In
    % English: "all files that start with 2 numbers, have `curROI`
    % somewhere in the middle of the file name, and ends with 
    % `trialtypeRSAmatrix.csv`"
    regExp = ['^[0-9]{2}y.*' curROI '.*trialtypeRSAmatrix\.csv'];
    
    % select all of the averagetrialtypeRSAmatrices using spm_select. Put
    % the full paths into a cell string
    RSA_matrix_filenames = cellstr(spm_select('FPListRec', rsa_results_dir, regExp));
    
    % reading in the averagetrialtypeRSAmatrices using `cellfun`
    RSA_matrices         = cellfun(@dlmread, RSA_matrix_filenames, 'UniformOutput', false);
    
    % select only the lower triangular part of these matrices using `tril`
    % and `cellfun`
    RSA_lwrTrian         = cellfun(@tril, RSA_matrices, 'UniformOutput', false);
    
    % extract the non-zero correlations using the `extract_nonzero` function
    % defined at the bottom of this function
    correlation          = cellfun(@extract_nonzero, RSA_lwrTrian, 'UniformOutput', false);  
    
    % Figure out all possible combinations of the trial types, select only
    % the unique combinations, join the combinations into a
    % TrialType1-TrialType2 naming format, and repeat it to match the length
    % of `correlation`
    
    % Find all possible combinations of the trialtypesOfInterest using the
    % function `allcomb`
    TrialCombination     = allcomb(trialtypesOfInterest, trialtypesOfInterest);
    
    % Select ONLY the unique TrialCombinations
    TrialCombination     = grab_unique(RSA_lwrTrian{1}, TrialCombination);
    
    % Record the number of unique trial combinations
    nUniqueCombos        = length(TrialCombination);
    
    % create the TrialType1-TrialType2 labels
    TrialCombination     = strcat(TrialCombination(:, 1), repmat({'-'}, size(TrialCombination, 1), 1), TrialCombination(:, 2));
    
    % repeate the TrialType1-TrialType2 labels `size(correlation)` times
    TrialCombination     = repmat(TrialCombination, size(correlation, 1), 1);
    
    % extract the subject IDs from the RSA_matrix_filenames
    
    % use `regexp` to grab the strings the match the pattern "two numbers
    % followed by an "o" OR a "y" followed by 3 more numbers"
    subjectID            = regexp(RSA_matrix_filenames, '[0-9]{2}[oy][0-9]{3}', 'match');
    
    % un-nest the resulting cell array, using Kyle's `unNest_cell_array`
    % function
    subjectID            = unNest_cell_array(subjectID);
    
    % repeate each element in the subjectID vector `nUniqueCombos` times, 
    % one for each unique combination of the `trialtypesOfInterest`
    subjectID            = repelem(subjectID, nUniqueCombos);
    
    % Vertically concatenate all of the correlation values
    correlation          = vertcat(correlation{:});
    
    % create a table
    TABLE                = table(subjectID, TrialCombination, correlation);
    
    % write the table to a csv (comma seperated value) file, basically an
    % excel sheet
    writetable(TABLE, fullfile(rsa_results_dir, ['compiled_' curROI '_trialtypeRSAmatrix.csv']))

end
    
%%% Subfunction

function out = extract_nonzero(in)

    out = in(in~=0);

end
  
function out = grab_unique(lwrTri, combos)

    vec = find(lwrTri ~= 0);
    out = combos(vec, :); %#ok<FNDSB>
    
end

end