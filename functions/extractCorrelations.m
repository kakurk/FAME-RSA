function [vector] = extractCorrelations(ds_dsm, trialtype1, trialtype2)
% extract correlations from a RSA matrix matching these two trial types
% ds_dsm is a cosmo dataset output from the cosmo_dissimilarity_matrix_measure
% function, with the following fields:
%
%   ds_dsm
%       .samples = a squareform vector of dissimilarity values
%       .sa
%           .targets1 = a double vector of the target IDs for this
%                       correlation pair
%           .targets2 = a double vector of the target IDs for this
%                       correlation pair
%           .target1labels = cell array of labels for this correlation pair
%           .target2labels = cell array of labels for this correlation pair
%
%   trialtype1 = name of the first trial type of interest, e.g. "RecHit"
%   trialtype2 = name of the second trial type of interest, e.g., "FamHit"
%
%   vector = a vector of all of the trialtype1-trialtype2 correlations
%            found in the in input ds_dsm
%
% See also: squareform, cosmo_dissimilarity_matrix_measure

% create two filters, on for the first target label and one for the second
% target label and combined the two to find all trial pairs that are off
% the trialtype1-trialtype2 variety
filter1 = cellfun(@isempty, regexp(ds_dsm.sa.target1labels, trialtype1));
filter2 = cellfun(@isempty, regexp(ds_dsm.sa.target2labels, trialtype2));
vector  = ds_dsm.samples(filter1 & filter2);

% 1-r --> r
vector  = (vector - 1) * -1;

end