function B3_slice_time_correction(SJ,run, run_dir, filter,n_slices,slice_order,refslice,TR)

f = spm_select('List', run_dir, filter);
numVols = size(f,1);
% create SPM style file list for model specification
%V=spm_vol([repmat([run_dir filesep], 1, 1) f(1,:) ',:']);
V=spm_vol([run_dir filesep f(1,:)]);
files={};
for i=1:(size(V,1))
    files{i} = [run_dir filesep strtrim(f(1,:)) ',' int2str(i)];
end

matlabbatch{1}.spm.temporal.st.scans = {files'};
matlabbatch{1}.spm.temporal.st.nslices = n_slices;
matlabbatch{1}.spm.temporal.st.tr = TR;
matlabbatch{1}.spm.temporal.st.ta = TR-(TR/n_slices);
matlabbatch{1}.spm.temporal.st.so = slice_order;
matlabbatch{1}.spm.temporal.st.refslice = refslice;
matlabbatch{1}.spm.temporal.st.prefix = 'a';

display(['Slice time correction: ' SJ ', ' run ])
fprintf('\n');
inputs = cell(0, 1);

% save the job variable to disc
spm('defaults', 'FMRI');
spm_jobman('Run', matlabbatch, '', inputs{:});
%spm_jobman('serial', matlabbatch, inputs{:});
clear jobs


