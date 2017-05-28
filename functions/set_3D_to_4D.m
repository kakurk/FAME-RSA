function matlabbatch = set_3D_to_4D(volumes, outputname)
    % set_3D_to_4D  set the matlabbatch structure for 3-D --> 4-D
    %               concatenation.

    matlabbatch{1}.spm.util.cat.vols  = volumes;
    matlabbatch{1}.spm.util.cat.name  = outputname;
    matlabbatch{1}.spm.util.cat.dtype = 0; % same as input
end