function C2_glm_1stLevel(SJs, sj, subj_dir, outputfolder_1st, prefix_func, tr, fmri_t, fmri_t0, hpf, runs, condnames, onsets, duration, hm, cc)
% similar to C1_glm_1stLevel

% set SPM defaults
spm('defaults','fmri')
global defaults;
global UFp;
spm_jobman('initcfg');

tgt_dir=fullfile(subj_dir, outputfolder_1st);
name_add='';
if hm
    name_add=[name_add 'hm'];
end
if cc
    name_add=[name_add 'cc'];
end
if ~isempty(name_add)
    tgt_dir=[tgt_dir '_' name_add];
end

if ~exist(tgt_dir, 'dir')
    mkdir(tgt_dir)
end

% output directory
matlabbatch{1}.spm.stats.fmri_spec.dir = cellstr(tgt_dir);
% timing parameters
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = tr;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = fmri_t;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = fmri_t0;

%% allocate the data per session (session is run in this case, confusing, i know)
for r = 1:size(runs,2)
    run_dir = fullfile(subj_dir, 'func');
    if exist(run_dir)==7
%         run_dir
        f = spm_select('List',run_dir, ['^' prefix_func runs{sj,r}]);
        V=spm_vol([run_dir filesep f(1,:)]);
        files={};
        for i=1:(size(V,1))
            files{i} = [run_dir filesep strtrim(f(1,:)) ',' int2str(i)];
        end
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).scans = files';
        
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi = {''}; % multiple conditions
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress = struct('name', {}, 'val', {}); %regressors
        
        %%% multiple regressors
        cov=1;
        multi_reg={''};
        %head motion
        if hm
            
            n=1;
            k = strfind(runs{sj, r}, 'd.nii')
            f1=spm_select('List', run_dir, ['^rp_' currPrefix(n:end) runs{sj, r}(1:k) '.txt']);
            while isempty(f1)
                n=n+1;
                f1=spm_select('List', run_dir, ['^rp_' currPrefix(n:end) runs{sj, r}(1:k) '.txt']);
            end
            
            if ~isempty(f1)
                multi_reg(cov,1)={[run_dir filesep f1]};
                cov=cov+1;
            else
                display('################################################################################')
                display(['############### ' SJs{sj} ', ' runs{r} ': Headmotion Parameters do not exsist #############'])
                display('################################################################################')
            end
        end
        
        %CompCorr
        if cc
            k = strfind(runs{sj, r}, 'd.nii')
            f2=spm_select('List',run_dir,['^' prefix_func runs{sj, r}(1:k) '_CompCorPCs.txt']);
            n=1;
            while isempty(f2)
                f2=spm_select('List',run_dir,['^' prefix_func(n:end) runs{sj, r}(1:k) '_CompCorPCs.txt']);
                n=n+1;
            end
            if ~isempty(f2)
                multi_reg(cov,1)={[run_dir filesep f2]};
                cov=cov+1;
            end
        end
        
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi_reg = multi_reg; % multiple regressors
        
        matlabbatch{1}.spm.stats.fmri_spec.sess(r).hpf = hpf; % high pass filter
        
        for c = 1:length(condnames)
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(c).name = condnames{c};
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(c).onset = onsets{sj,r,c};
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(c).duration = duration;            

             %single events für motion/antwort
%             if c == 10
%                 matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(c).duration = 0;
%             end
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(c).tmod = 0; % no time modulation
            
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(c).pmod = struct('name', {}, 'param', {}, 'poly', {}); % no parametric modulations
            
            matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(c).orth = 1; % orthogonalise modulations
        end
    end
    clear V f run_dir
end


%% additional model parameters
% no factorial design
matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
% model hrf and first temporal derivative
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
% model interactions
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
% global normalization
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
% masking
matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;   % masking threshold, proportion of globals
matlabbatch{1}.spm.stats.fmri_spec.mask = {''}; % explicit mask
% autocorrelation modelling (whitening filter)
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

%%  Model Estimation 
matlabbatch{2}.spm.stats.fmri_est.spmmat = {[tgt_dir filesep 'SPM.mat']};
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

%% create the model
fprintf(['Creating GLM\n'])
spm('defaults', 'FMRI');
spm_jobman('run', matlabbatch);

% clear job variable
clear matlabbatch
