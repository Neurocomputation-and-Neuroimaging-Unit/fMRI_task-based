function C3_contrast_1stLevel(subj_dir, analysisfolder, cnames, cvecs, del, n_hm, n_cc, runs)

%% Allocate SPM.mat file
matlabbatch{1}.spm.stats.con.spmmat = {fullfile(subj_dir, analysisfolder, 'SPM.mat')};
    
% Number of Runs
% rundir=dir([subj_dir filesep 'func']);
numRuns=size(runs,2);       
% Number of contrasts
numCons  = numel(cnames);

% Cycle over contrast specifications
for c = 1:numCons
    convec=cvecs{c};
    if n_hm>0
        convec=[convec zeros(1,n_hm)];
    end
    if n_cc>0
        convec=[convec zeros(1,n_cc)];
    end
    
    % Allocate t-contrast structure
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.name    = cnames{c};       
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.weights  = repmat(convec, 1, numRuns);         
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none';             
end 

% Delete existing contrasts (1=yes)
matlabbatch{1}.spm.stats.con.delete = del;
    
% Run the job
fprintf(['Computing 1st Level Contrasts\n'])
spm('defaults', 'FMRI');
spm_jobman('run', matlabbatch);

% clear job variable
clear matlabbatch