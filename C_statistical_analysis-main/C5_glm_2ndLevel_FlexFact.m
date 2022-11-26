function C5_glm_2ndLevel_FlexFact(src_dir, SJin, outputfolder, dir_1st, con_images, runs)

% This function performs a second level statistical analysis by performing
% a Flexible Factorial Design over the specified con_images

% target directory that will contain the created job.mat file
tgt_dir      = [src_dir filesep outputfolder];

% Create tgt_dir
if ~exist(tgt_dir, 'dir')
    mkdir(tgt_dir)
end

matlabbatch{1}.spm.stats.factorial_design.dir = {tgt_dir};
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).name = 'Subject';
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).variance = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).name = 'Type';
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).dept = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).variance = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).ancova = 0;
% matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).name = 'Frequency';
% matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).dept = 1;
% matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).variance = 0;
% matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).gmsca = 0;
% matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).ancova = 0;
for sj = 1:numel(SJin)
    % create SPM style file list for model specification.
    fNames = [];
    for con = con_images
%         for dir = 1:length(dir_1st)
            dir
            if con < 10
                filt            = [['con_000' num2str(con)] '*\.nii$'];
            else
                filt            = [['con_00'  num2str(con)] '*\.nii$'];
            end
            temp_dir      = fullfile(src_dir, SJin{sj}, dir_1st)
            cd(temp_dir)
            f             = spm_select('List', temp_dir, filt);
%             f             = spm_select('List', temp_dir, 's5wres_accuracy_minus_chance.nii');
            fs            = cellstr([repmat([temp_dir filesep], 1, 1) f repmat(',1', 1, 1)]);
            fNames          = [fNames; fs];
%         end
    end
        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(sj).scans = fNames;

        matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(sj).conds = [1:length(con_images)]'; %% AAACHTUNG
%     matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(sj).conds = [(1:length(con_images))'];
end
matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{1}.fmain.fnum = 1; % 1 for subject factor
matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{2}.fmain.fnum = 2;
% matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{3}.fmain.fnum = 3;
% matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{2}.inter.fnums = [2 3];
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

%%  Model Estimation
matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(tgt_dir, 'SPM.mat')};
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

spm('defaults', 'FMRI');
spm_jobman('run', matlabbatch);

% clear job variable
clear matlabbatch

