function B95_reslice_masks_to_functional(mean_img, mask1, mask2)

fileset = [mean_img; {[mask1 ',1']; [mask2 ',1']}];

matlabbatch{1}.spm.spatial.realign.write.data = fileset;
matlabbatch{1}.spm.spatial.realign.write.roptions.which = [1 0];
matlabbatch{1}.spm.spatial.realign.write.roptions.interp = 4;
matlabbatch{1}.spm.spatial.realign.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.realign.write.roptions.mask = 1;
matlabbatch{1}.spm.spatial.realign.write.roptions.prefix = 'r';

spm('defaults', 'FMRI');
spm_jobman('run', matlabbatch);

end