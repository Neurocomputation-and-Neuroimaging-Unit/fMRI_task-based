function B9_smoothing_run(run_dir, SJ, filter,kernel_size)

warning off

fileset = {};
% select the files
f = spm_select('List', run_dir, filter);
% number of volumes
numVols = size(f,1);
% create SPM style file list for model specification
V=spm_vol([repmat([run_dir filesep], numVols, 1) f ',:']);
for i=1:size(V,1)
    fileset{i} = [repmat([run_dir filesep], numVols, 1) f repmat([',' int2str(i)] , numVols, 1)];
end
clear f;
% create prefix
aa=num2str(unique(kernel_size));
if length(aa)>1
    aa=num2str(kernel_size);
end

matlabbatch{1}.spm.spatial.smooth.data = fileset';
matlabbatch{1}.spm.spatial.smooth.fwhm = kernel_size;
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = ['s' aa(~isspace(aa))];

inputs = cell(0, 1);

% save the job variable to disc
fprintf(['Smoothing ', SJ, '\n'])
spm('defaults', 'FMRI');
spm_jobman('serial', matlabbatch, '', inputs{:});
%spm_jobman('serial', matlabbatch, inputs{:});
clear jobs

