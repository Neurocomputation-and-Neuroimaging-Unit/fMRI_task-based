%%%% batch script for second level AFTER decoding of any type

clear all
close all

%SPM-path
SPM_path  = 'C:\Users\saraw\Desktop\BA\EXPRA2019_HIVR\Toolboxes\spm12';

%data source directory
src_dir      = 'C:\Users\saraw\Desktop\BIDS\fmri';
dec_folder = 'DEC';

addpath(genpath('C:\Users\saraw\Desktop\BA\EXPRA2019_HIVR\Toolboxes\hMRI-toolbox-0.4.0'));
addpath('C:\Users\saraw\Desktop\BA\EXPRA2019_HIVR\Toolboxes\spm12');

%subject identifiers if all subjects are to be included
%%%% to do: SJs to sub
cd(src_dir) 
pb=dir('sub*');
for i=1:length(pb)
    SJs(1,i)={pb(i).name};
end

display('Subjects found:')
SJs % analysis for these subjects 
display('Subjects to exclude:')
excludeSJ = [] % delete those sjs
zip_files = dir(fullfile(src_dir, '**', '*.gz'));

if ~isempty(zip_files) 
	for z = 1:size(zip_files, 1)
		gunzip([zip_files(z).folder filesep zip_files(z).name]);
        delete([zip_files(z).folder filesep zip_files(z).name]);
	end
end

start_prefix = 's5w';

for sb = 1:numel(SJs)
    if ismember(sb, excludeSJ)
        continue;
    else
        cd([src_dir filesep SJs{sb} filesep dec_folder]);
        dd = dir('mean*');
        for d = 1:size(dd,1)
            cd([src_dir filesep SJs{sb} filesep dec_folder filesep dd(d).name])
            rd = dir([start_prefix 'res*.nii']);
            folds(sb, d) = {rd.folder};
            decs(sb, d) = {rd.name};
        end
    end
end

%% INPUT

outputfolder = '2nd_level_test';
filt = start_prefix;
explicit_mask={'C:\Users\saraw\Desktop\BA\EXPRA2019_HIVR\brain_mask.nii'}; %%%%%%%%%% find

%%

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
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).name = 'Bin';
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).dept = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).variance = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).ancova = 0;

for sj = 1:numel(SJs)
    % create SPM style file list for model specification.
    fNames = [];
%     for con = con_images
        for ds = 1:size(decs, 2)
            cd(folds{sj,ds});
            f             = spm_select('List', folds{sj,ds}, decs{sj,ds});
            fs            = cellstr([repmat([folds{sj,ds} filesep], 1, 1) f repmat(',1', 1, 1)]);
            fNames          = [fNames; fs];
        end
%     end
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(sj).scans = cellstr(fNames);

    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(sj).conds = [1:size(decs,2)]'; %% AAACHTUNG
end
matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{1}.fmain.fnum = 1; % 1 for subject factor
matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{2}.fmain.fnum = 2;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;


matlabbatch{1}.spm.stats.factorial_design.masking.em = explicit_mask;


%%  Model Estimation
matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(tgt_dir, 'SPM.mat')};
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

spm('defaults', 'FMRI');
spm_jobman('run', matlabbatch);

% clear job variable
clear matlabbatch











