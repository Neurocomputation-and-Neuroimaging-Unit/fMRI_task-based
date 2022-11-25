function D2_FIR(src_dir, onsets, SJs, TR, WMdelay, currPrefix, runs, excludeSJ, fir_out)
% This script implements a FIR model as commonly used in Time-resolved
% Decoding.
% Each time-bin of a WM-delay phase is modelled individually
% The WM-delay is of 12seconds duration, which results in 6 bins (TR=2s)
% The underlying directory-structure is explained in the Tutorial
%
% Copyright: Timo Torsten Schmidt, Freie Universität Berlin

% enter data directory here
data_dir  = src_dir;

%% Loop over Subjects
for sj = 1:numel(SJs)

    if ismember(sj, excludeSJ)
        continue;
    else
        % Composition of the Subject directory
        sj_dir = fullfile(data_dir, SJs{sj});
        
        % SPM defaults
        spm('defaults','fmri');
        global defaults;
        global UFp;
        spm_jobman('initcfg');
        % OUTPUT Directory (as subdirectory of the SJ directory)
        tgt_dir = fullfile(sj_dir, ['FIR_' fir_out]);
        if ~exist(tgt_dir, 'dir')
            mkdir(tgt_dir)
        end
    
        %********************************************************
        %               Setting FIR parameter
        % *******************************************************
        % Output directory
        jobs{1}.stats{1}.fmri_spec.dir = cellstr(tgt_dir);
        % timing parameters
        jobs{1}.stats{1}.fmri_spec.timing.units     = 'secs';
        jobs{1}.stats{1}.fmri_spec.timing.RT        = TR;
        jobs{1}.stats{1}.fmri_spec.timing.fmri_t    = 16;
        jobs{1}.stats{1}.fmri_spec.timing.fmri_t0   = 1;
        % FIR-SPECIFICATION
        jobs{1}.stats{1}.fmri_spec.bases.fir.length = WMdelay; % 9 seconds
        jobs{1}.stats{1}.fmri_spec.bases.fir.order  = WMdelay/TR;  % 6 time-bins
        % Other Specifications
        jobs{1}.stats{1}.fmri_spec.fact              = struct('name', {}, 'levels', {});
        jobs{1}.stats{1}.fmri_spec.volt              = 1;
        jobs{1}.stats{1}.fmri_spec.global            = 'None';
        jobs{1}.stats{1}.fmri_spec.mask              = {''};
        jobs{1}.stats{1}.fmri_spec.cvi               = 'None';
    
        % ********************************************************
        %            Specify the Desing/Conditions/Onsets
        % ********************************************************
        % CONDITON NAMES
        condnames = {'WM1', 'WM2', 'WM3', 'WM4'};
    
        % TODO
        % get all runIDs from file names in 'func' folder
        nifti_dir = fullfile(sj_dir, 'func');
        runIDs = 1:4;
    
        % Loop over Sessions, as the logfiles (with onset data) are saved
        % separately per session
        for s = 1:size(runs,2)
            s
            runID = runIDs(s);
            % high-pass cut-off (might also be reasonable at 300 for this design
            jobs{1}.stats{1}.fmri_spec.sess(s).hpf     = 128;
            % Allocation of Data (EPIs/Images) for the current Session
            % TODO get rid of hardcoded file names here
            %filename = fullfile(nifti_dir, sprintf('%s_task-pain_run-%i_bold.nii', SJs{sj}, runID));
            filename = fullfile(nifti_dir, [currPrefix runs{sj,runID}]);
            % use 'expand' to read 4d nifti file
            f = spm_select('expand', filename); %%%%%%%%%%%%%%%%%%%% for 4d ??????
            %jobs{1}.stats{1}.fmri_spec.sess(s).scans   = cellstr([f, repmat(',1', size(f,1),1)]);
            jobs{1}.stats{1}.fmri_spec.sess(s).scans   = cellstr(f);
    
            % The consecutive codes does not have to be understood neccessarily
            % It is highly specific to the way the onsets were coded in the
            % Logfiles
            % load logfiles - there is one per run
            % TODO actual file names are probably different, also probably tsv
            % files
            % WM1
            onsets2 = onsets{sj,s,1};
            jobs{1}.stats{1}.fmri_spec.sess(s).cond(1).name     = condnames{1};
            jobs{1}.stats{1}.fmri_spec.sess(s).cond(1).onset    = onsets2+3;
            jobs{1}.stats{1}.fmri_spec.sess(s).cond(1).duration = 0;
            jobs{1}.stats{1}.fmri_spec.sess(s).cond(1).tmod     = 0;
            jobs{1}.stats{1}.fmri_spec.sess(s).cond(1).pmod     = struct('name', {}, 'param', {}, 'poly', {});
            clear onsets2;
            % WM2
            onsets2 = onsets{sj,s,2};
            jobs{1}.stats{1}.fmri_spec.sess(s).cond(2).name     = condnames{2};
            jobs{1}.stats{1}.fmri_spec.sess(s).cond(2).onset    = onsets2+3;
            jobs{1}.stats{1}.fmri_spec.sess(s).cond(2).duration = 0;
            jobs{1}.stats{1}.fmri_spec.sess(s).cond(2).tmod     = 0;
            jobs{1}.stats{1}.fmri_spec.sess(s).cond(2).pmod     = struct('name', {}, 'param', {}, 'poly', {});
            clear onsets2;
            % WM3
            onsets2 = onsets{sj,s,3};
            jobs{1}.stats{1}.fmri_spec.sess(s).cond(3).name     = condnames{3};
            jobs{1}.stats{1}.fmri_spec.sess(s).cond(3).onset    = onsets2+3;
            jobs{1}.stats{1}.fmri_spec.sess(s).cond(3).duration = 0;
            jobs{1}.stats{1}.fmri_spec.sess(s).cond(3).tmod     = 0;
            jobs{1}.stats{1}.fmri_spec.sess(s).cond(3).pmod     = struct('name', {}, 'param', {}, 'poly', {});
            clear onsets2;
            % WM4
            onsets2 = onsets{sj,s,4};
            jobs{1}.stats{1}.fmri_spec.sess(s).cond(4).name     = condnames{4};
            jobs{1}.stats{1}.fmri_spec.sess(s).cond(4).onset    = onsets2+3;
            jobs{1}.stats{1}.fmri_spec.sess(s).cond(4).duration = 0;
            jobs{1}.stats{1}.fmri_spec.sess(s).cond(4).tmod     = 0;
            jobs{1}.stats{1}.fmri_spec.sess(s).cond(4).pmod     = struct('name', {}, 'param', {}, 'poly', {});
            clear onsets2;
        end
    
        % Create model
        fprintf(['Creating GLM\n'])
        spm_jobman('run', jobs);
        clear jobs
    
        %  Model Estimation
        load(fullfile(tgt_dir, 'SPM.mat'));
        fprintf(['Estimating GLM \n']);
        cd(tgt_dir);
        SPM = spm_spm(SPM);
        clear SPM;
    end
end
end