function  ds_sa = correlation_summary_measure(ds, varargin)
% Computes a ers correlation measure
%
% d=cosmo_correlation_measure(ds[, args])
%
% Inputs:
%  ds             dataset structure with fields .samples, .sa.targets and
%                 .sa.chunks
%  args           optional struct with the following optional fields:
%    .template    QxM matrix for Q classes in each chunk. This matrix
%                 weights the correlations across the two halves.
%                 If ds.sa.targets has only one unique value, it must be
%                 set to the scalar value 1; otherwise it should
%                 have a mean of zero. If omitted, it has positive values
%                 of (1/Q) on the diagonal and (-1/(Q*(Q-1)) off the
%                 diagonal.
%                 (Note: this can be used to test for representational
%                 similarity matching)
%    .corr_type   Type of correlation: 'Pearson','Spearman','Kendall'.
%                 The default is 'Pearson'.
%    .post_corr_func  Operation performed after correlation. (default:
%                     @atanh)
%
% Output:
%    ds_sa        Struct with fields:
%      .samples   Scalar indicating how well the template matrix
%                 correlates with the correlation matrix from the two
%                 halves (averaged over partitions). By default:
%                 - this value is based on Fisher-transformed correlation
%                   values, not raw correlation values
%                 - this is the average of the (Fisher-transformed)
%                   on-diagonal minus the average of the
%                   (Fisher-transformed) off-diagonal elements of the
%                   correlation matrix based on the two halves of the data.
%      .sa        Struct with field:
%        .labels  if output=='corr'
%
% Example:
%
% Notes:
%   - by default the post_corr_func is set to @atanh. This is equivalent to
%     a Fisher transformation from r (correlation) to z (standard z-score).
%     The underlying math is z=atanh(r)=.5*log((1+r)./log(1-r)).
%     The rationale is to make data more normally distributed under the
%     null hypothesis.
%     Fisher-transformed correlations can be transformed back to
%     their original correlation values using 'tanh', which is the inverse
%     of 'atanh'.
%   - To disable the (by default used) Fisher-transformation, set the
%     'post_corr_func' option to [].
%
% References
%   - Haxby, J. V. et al (2001). Distributed and overlapping
%     representations of faces and objects in ventral temporal cortex.
%     Science 293, 2425?2430

%% Parse Arguments

nchunk1samples = length(find(ds.sa.chunks == 1));
nchunk2samples = length(find(ds.sa.chunks == 2));

if isempty(varargin) 
    % defaults
    template       = ones(nchunk1samples, nchunk2samples) / (nchunk1samples * nchunk2samples);
    corr_type      = 'spearman';
    post_corr_func = @atanh;
else
    template       = varargin{1}.template;
    corr_type      = varargin{1}.corr_type;
    post_corr_func = varargin{1}.post_corr_func;
end

%% Check Assumptions

% assume that there are two and only two chunks
assert(all(ismember(unique(ds.sa.chunks), [1;2])), 'There must be two and only two chunks.')

% assume that the template matrix is the same size as expected size of the
% correlation matrix
assert(all(ismember(size(template), [nchunk1samples, nchunk2samples])), 'The size of the template martix must be the same size as the correlation matrix.')

%% Correlate and Summarize

% split data into "halves"
ds_half1 = cosmo_slice(ds, ds.sa.chunks == 1);
ds_half2 = cosmo_slice(ds, ds.sa.chunks == 2);

% correlate
craw = cosmo_corr(ds_half1.samples', ds_half2.samples', corr_type);

% post corr function
cpost = post_corr_func(craw);

% summarize
summary = cpost .* template;
summary = sum(summary(:));

%% Output

ds_sa.samples   = summary;
ds_sa.sa.labels = 'correlation_summary';

