function B5_coregister_est(currPrefix, func_dir, struct_dir, sj, filter_struct, runs)

warning off

f1 = spm_select('List', struct_dir, filter_struct);
numVols = size(f1,1);
structural = cellstr([repmat([struct_dir filesep], numVols, 1) f1 repmat(',1', numVols, 1)]);

f2 = spm_select('List', func_dir, ['^mean.*\.nii']);
numVols = size(f2,1);
mean_img   = cellstr([repmat([func_dir filesep], numVols, 1) f2 repmat(',1', numVols, 1)]);

for r = 1:size(runs,2)
%     data_dir = [sj_dir filesep sbj_level_folder filesep 'mean_bin_' num2str(b)];
%     f3 = spm_select('List', data_dir, '^res_accuracy_minus_chance.nii');
%     numVols = size(f3,1);
%     Images{b} = cellstr([repmat([data_dir filesep], numVols, 1) f3 repmat(',1', numVols, 1)]);
    f3 = spm_select('ExtFPList', func_dir, ['^' currPrefix runs{sj,r}],Inf);
    numVols = size(f3,1);
%     Images{r,:}=cellstr([repmat([func_dir filesep], numVols, 1) f3 repmat(',1', numVols, 1)]);
%     Images{r,:}=cellstr([repmat([func_dir filesep], numVols, 1) f3]);
    Images{r,:}=cellstr(spm_select('ExtFPList', func_dir,['^' currPrefix runs{sj,r}],Inf));
end

%--------------------------------------------------------------------------
%---------------------------- Coregister (Estimate) -----------------------
%--------------------------------------------------------------------------
for i = 1:size(Images,1)
    matlabbatch{1}.spm.spatial.coreg.estimate.ref = structural;
    matlabbatch{1}.spm.spatial.coreg.estimate.source = mean_img;
    matlabbatch{1}.spm.spatial.coreg.estimate.other = Images{i,:};
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    
    % run job
    spm_jobman('run', matlabbatch)
    clear jobs
end
