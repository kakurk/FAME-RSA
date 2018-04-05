function statstable = extract_top_n(images, masks, n, summaryfunction)
% extract_top_n.  extract the top n voxels from a given region, summarize
% them, and return a tidyverse formatted statistics table for analysis with
% a different software package.
%
%   images = cell string of full paths to statistics maps that you want to
%            extract from. Example = 
%
%       images = {'/gpfs/group/nad12/default/nad12/FAME8/RSA/models/ERS_results_diagnol/sub-18y404_ers_searchight.nii'
%                 '/gpfs/group/nad12/default/nad12/FAME8/RSA/models/ERS_results_diagnol/sub-18y566_ers_searchight.nii'
%                 ...
%                 '/gpfs/group/nad12/default/nad12/FAME8/RSA/models/ERS_results_diagnol/sub-25y543_ers_searchight.nii'}
%
%   masks = cell string of full paths to the binary masks detailing the
%          search areas. Example = 
%
%       masks = {'/gpfs/group/nad12/default/nad12/FAME8/RSA/ROIs/HC_left.nii'
%                '/gpfs/group/nad12/default/nad12/FAME8/RSA/ROIs/HC_right.nii'}
%
%
%   n = a double set to the number of top voxels to extract. Example = 
%
%       n = 10; % top 10 voxels
%
%   summaryfunction = a function handle of a function to summarize the top
%                     n voxels. Example = 
%
%       summaryfunction = @mean
%
%
%   statstable = MATLAB table formatted in tidyverse format. Example = 
%
%
%                                                    images                                                                               masks                                n     summaryfunction    summaryvalues
%     ____________________________________________________________________________________________________    _____________________________________________________________    __    _______________    _____________
% 
%     '/gpfs/group/nad12/default/nad12/FAME8/RSA/models/ERS_results_diagnol/sub-18y404_ers_searchight.nii'    '/gpfs/group/nad12/default/nad12/FAME8/RSA/ROIs/HC_left.nii'     10    'mean'             0.044364     
%     '/gpfs/group/nad12/default/nad12/FAME8/RSA/models/ERS_results_diagnol/sub-18y566_ers_searchight.nii'    '/gpfs/group/nad12/default/nad12/FAME8/RSA/ROIs/HC_left.nii'     10    'mean'             0.054113     
%     '/gpfs/group/nad12/default/nad12/FAME8/RSA/models/ERS_results_diagnol/sub-20y297_ers_searchight.nii'    '/gpfs/group/nad12/default/nad12/FAME8/RSA/ROIs/HC_left.nii'     10    'mean'             0.051168     
%     '/gpfs/group/nad12/default/nad12/FAME8/RSA/models/ERS_results_diagnol/sub-20y396_ers_searchight.nii'    '/gpfs/group/nad12/default/nad12/FAME8/RSA/ROIs/HC_left.nii'     10    'mean'             0.025874 
%
%
% Written by Kyle Kurkela, kyleakurkela@gmail.com 4/4/2018

%% routine

%-- Sanity Checks

% spm must be on the searchpath
assert(~isempty(which('spm')), 'spm must be on the MATLAB seach path')

% 'vol' the ROI mask(s). see spm_vol.m
Vmasks = spm_vol(masks);

% check that the input masks are binary
for v = 1:length(Vmasks)
    
    % read in the current ROI
    R = spm_read_vols(Vmasks{v});
    
    % assert that this mask is binary
    isbinary = all(ismember(unique(R), [0;1]));
    assert(isbinary, 'One of the input masks is not binary. All input musts MUST be binary.')
    
end

% 'vol' the imput images. see spm_vol.m
Vimages = spm_vol(images);

% check that each of the input images are in the same space
mat = Vimages{1}.mat;

for v = 2:length(Vimages)
     areInSameSpace = all(all(Vimages{v}.mat == mat));
     mat = Vimages{v}.mat;
end

assert(areInSameSpace, 'All of the input images are NOT in the same space')

%-- Extraction!

% initalize
summaryvalues = [];

% for each roi mask
for v = 1:length(Vmasks)
    
    % read in current ROI
    [R, XYZmm] = spm_read_vols(Vmasks{v});
    
    % inclusive mask (i.e., only look at where the image is == 1). Mask
    % MUST BE BINARY!!
    XYZmm = XYZmm(:, R == 1);
    
    % convert mm coordinates --> vox coordinates using the images space
    XYZ = Vimages{1}.mat\[XYZmm;ones(1,size(XYZmm,2))];
    
    % extract all voxels in this mask for each image, creating an nimages
    % x nvoxels-in-this-mask matrix
    extracted_data = spm_get_data(images, XYZ);
    
    % find the top n voxels for each subject
    topnvoxelsvalues = maxn(extracted_data, n, 2);

    % summarize across the top n voxels using the summarizefunction
    summaryvalues  = vertcat(summaryvalues, summaryfunction(topnvoxelsvalues, 2, 'omitnan')); %#ok<AGROW>
    
end

%-- Build the Statistics Table

number_of_images = length(images);
number_of_masks  = length(masks);
combinations_of_masks_and_images = length(images) * length(masks);

% repeat the image names nmasks times.
images = repmat(images, number_of_masks, 1);

% repeat **each element** of masks nimages times.
masks  = repelem(masks, number_of_images);

% repeat n nimages x nmasks times
n      = repmat(n, combinations_of_masks_and_images, 1);

% repeat summaryfunction nimages x nmasks times. make a string
summaryfunction = repmat({func2str(summaryfunction)}, combinations_of_masks_and_images, 1);

% actually build the table
statstable = table(images, masks, n, summaryfunction, summaryvalues);

%% subfunction

function outmatrix = maxn(matrix, n, dim)
    
    outmatrix = [];
    
    for i = 1:n
        
        % find the max
        [Y, I] = max(matrix, [], dim);
        
        % build the output matrix
        outmatrix = cat(dim, outmatrix, Y);
        
        % remove the current max from the matrix
        matrix(:,I) = [];
        
    end
end

end