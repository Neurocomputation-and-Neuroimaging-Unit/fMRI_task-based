function B8_compcorr_run(data_dir, SJ, filter_data, numcomp, wm_mask_file, csf_mask_file)

spm('defaults','fmri');
spm_jobman('initcfg');

warning off

f = spm_select('List', data_dir, filter_data);
rest_file= [data_dir filesep f];

y_CompCor_PC(rest_file,{csf_mask_file, wm_mask_file},[rest_file(1:end-4) '_CompCorPCs.txt'], numcomp, 1, [], 2, 1);

display([SJ ' is done'])
