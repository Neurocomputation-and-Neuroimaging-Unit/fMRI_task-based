function B6_normalization_run(data_dir, struct_dir, sj, runs, vox_size, currPrefix)

spm('defaults','fmri');
spm_jobman('initcfg');

warning off

f1 = spm_select('List', struct_dir, '^y.*\.nii$');
numVols = size(f1,1);
parameter_file = cellstr([repmat([struct_dir filesep], numVols, 1) f1 ]);


% 
% V=spm_vol([repmat([data_dir filesep], numVols2, 1) f2 ',:']);
% for i=1:size(V,1)
%     functional_imgs{i} = [repmat([data_dir filesep], numVols2, 1) f2 repmat([',' int2str(i)] , numVols2, 1)];
% end

for r = 1:size(runs,2)

    f2 = spm_select('List', data_dir, ['^' currPrefix runs{sj,r}])
    numVols2 = size(f2,1);
%     data_dir = [sj_dir filesep sbj_level_folder filesep 'mean_bin_' num2str(b)];
%     f3 = spm_select('List', data_dir, '^res_accuracy_minus_chance.nii');
%     numVols = size(f3,1);
%     Images{b} = cellstr([repmat([data_dir filesep], numVols, 1) f3 repmat(',1', numVols, 1)]);
    V=spm_vol([repmat([data_dir filesep], numVols2, 1) f2 ',:']);
    for i=1:size(V,1)
        functional_imgs{i} = [repmat([data_dir filesep], numVols2, 1) f2 repmat([',' int2str(i)] , numVols2, 1)];
    end
% end
% functional_imgs'

%--------------------------------------------------------------------------
% --------------------------- Normalize  ----------------------------
%--------------------------------------------------------------------------
% for j = 1:size(functional_imgs,1)
    %%
    matlabbatch{1}.spm.spatial.normalise.write.subj.def = parameter_file;
    
    %%
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = functional_imgs';
    
    %%
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70;78 76 85];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = vox_size;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    
    
    %  save([tgt_dir, 'job_normalize.mat'], 'jobs');
    fprintf(['Normalizing ', num2str(sj), '\n']);
    spm_jobman('run', matlabbatch)
    clear jobs
end
