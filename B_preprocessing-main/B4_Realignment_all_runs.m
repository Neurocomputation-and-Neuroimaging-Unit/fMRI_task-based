function B4_Realignment_all_runs(runPath, run_files)

warning off
spm_figure('GetWin','Graphics');

fileset = {};

% LOOP over sessions to assemble SPM-style file-set
for r = 1:length(run_files)
        
%     if exist([sj_dir filesep 'func']) == 7
    
        funcPath = runPath;
        % select the files
        f = run_files{r};
        % number of volumes
        numVols = size(f,1);
        % create SPM style file list for model specification
        V=spm_vol([repmat([funcPath filesep], numVols, 1) f ',:']);
        for i=1:size(V,1)
            files{i,1} = [repmat([funcPath filesep], numVols, 1) f repmat([',' int2str(i)] , numVols, 1)];
        end
        fileset{r}=files;
        clear f V files;
        
%     end
end

if ~isempty(fileset)
    matlabbatch{1}.spatial{1}.realign{1}.estwrite.data             = fileset;
    matlabbatch{1}.spatial{1}.realign{1}.estwrite.eoptions.quality = 1; %%% for testing, set this lower, for actual preprocessing, leave it at 1!!
    matlabbatch{1}.spatial{1}.realign{1}.estwrite.eoptions.sep     = 4;
    matlabbatch{1}.spatial{1}.realign{1}.estwrite.eoptions.fwhm    = 5;
    matlabbatch{1}.spatial{1}.realign{1}.estwrite.eoptions.rtm     = 1;
    matlabbatch{1}.spatial{1}.realign{1}.estwrite.eoptions.interp  = 2;
    matlabbatch{1}.spatial{1}.realign{1}.estwrite.eoptions.wrap    = [0 0 0];
    matlabbatch{1}.spatial{1}.realign{1}.estwrite.eoptions.weight  = {''};
    matlabbatch{1}.spatial{1}.realign{1}.estwrite.roptions.which   = [2 1];
    matlabbatch{1}.spatial{1}.realign{1}.estwrite.roptions.interp  = 4;
    matlabbatch{1}.spatial{1}.realign{1}.estwrite.roptions.wrap    = [0 0 0];
    matlabbatch{1}.spatial{1}.realign{1}.estwrite.roptions.mask    = 1;
    matlabbatch{1}.spatial{1}.realign{1}.estwrite.roptions.prefix    = 'r';
    
    % save the job variable to disc
    spm('defaults', 'FMRI');
    spm_jobman('run', matlabbatch);
    clear matlabbatch
end

