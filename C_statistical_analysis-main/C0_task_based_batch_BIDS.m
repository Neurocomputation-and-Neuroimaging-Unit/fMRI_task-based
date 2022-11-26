%% decoding BATCH %%

% required toolboxes:
% the decoding toolbox TDT
%  https://doi.org/10.3389/fninf.2014.00088

%#####################################################
%#################### INPUT ##########################
%#####################################################


%SPM-path
SPM_path  = 'C:\Users\saraw\Desktop\BA\EXPRA2019_HIVR\Toolboxes\spm12';

%data source directory
src_dir      = 'C:\Users\saraw\Desktop\BIDS\test3_pilotdataIMACU';

addpath(genpath('C:\Users\saraw\Desktop\BIDS')); %%%%%%% your current working directory %%%%%%%%%%%%%% eg: 'F:\GForce\Dataanalysis\D_decoding_BIDS-main'
addpath(genpath('C:\Users\saraw\Desktop\BA\EXPRA2019_HIVR\Toolboxes\hMRI-toolbox-0.4.0'));
% addpath('C:\Users\saraw\Desktop\BA\EXPRA2019_HIVR\Toolboxes\decoding_toolbox');
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
excludeSJ = [4] % delete those sjs

zip_files = dir(fullfile(src_dir, '**', ['sub-', '*.gz']));
if ~isempty(zip_files) 
	for z = 1:size(zip_files, 1)
		gunzip([zip_files(z).folder filesep zip_files(z).name]);
        delete([zip_files(z).folder filesep zip_files(z).name]);
	end
end

%session & run identifiers
sessNum = 0;
if exist([src_dir filesep SJs{1} filesep 'ses-1'])==7 %################################### skript can't (yet) process multiple sessions per subject #################################
    cd([src_dir filesep SJs{1}])
    sd = dir('ses*');
    sessNum = length(sd);
    for sess = sessNum
        sessions(1, sess) = {sd(sess).name};
    end
    for sb = 1:numel(SJs)
        cd([src_dir filesep SJs{sb} filesep sessions{1} filesep 'func']);
        rd = dir('sub*.nii');
        for r = 1:length(rd)
            runs(sb, r) = {rd(r).name};
        end
    end
else
    for sb = 1:numel(SJs)
        cd([src_dir filesep SJs{sb} filesep 'func']);
        rd = dir('sub*.nii');
        for r = 1:length(rd)
            runs(sb, r) = {rd(r).name};
        end
    end
end

%anatomy identifier
ana=['anat'];

nifti_files = dir(fullfile(src_dir, '**', ['sub-', '*bold.nii'])); %look for all functional nifti files
anat_files = dir(fullfile(src_dir, '**', ['sub-', '*T1w.nii'])); %look for all anat nifti files

%now we get the data from the json file 
json_files = (dir(fullfile(src_dir, '**', ['task', '*json']))); %extract all json files, althoguh they should have the same info 
%because for some datasets (flicker and Ganzfeld) we have json files for
%each nift file and these are named differently, we have to check if the
%first command returns an empty structure. If yes, it means the json files
%have a different naming, starting with subject 
if isequal(size(json_files), [0, 1]) 
    json_files = (dir(fullfile(src_dir, '**', ['sub-', '*bold.json'])));
end

json_file = [json_files(1).folder, filesep, json_files(1).name]; %we select the first json file to extract metadata from 
TR_json = get_metadata_val(json_file,'RepetitionTime') / 1000; % repetition time in sec
slice_timing = get_metadata_val(json_file,'SliceTiming'); %extract slice timing
n_slices_json = height(slice_timing); %compute number of slices from slice timing
[~,y]= sort(slice_timing); %compute slice order 
slice_order = y';

%now get the same info from nifti header  
nifti_file_metadata = [nifti_files(1).folder, filesep, nifti_files(1).name]; 
info = niftiinfo(nifti_file_metadata);
TR_nifti = info.PixelDimensions(4); 
n_slices_nifti = info.ImageSize(3);
vox_size=repmat(info.PixelDimensions(1),1,3);

%compare json and nifti header 
if round(TR_nifti, 4) ~= round(TR_json, 4)
    warning ("TR does not match between json file and nifti") 
end 
if n_slices_json ~= n_slices_nifti
    warning ("Number of slices does not match between json file and nifti") 
end 

%%% I suggest using the json values at least for TR since we know it is
%%% missing in the nifti headers for some datasets
n_slices = n_slices_json; % number of slices
TR=TR_json; % repetition time in sec.

%% selection of analysis steps (1-5) to be performed
analysis_switch = [4]; % 1 2 3 4 5
start_prefix='s8wrba'; %eg. s8wra
%'' when bids and unpreprocessed; 
%rba/arb if already slicetime-corrected and realigned
%'' if doing normilazation of accuracy maps
%w if aready normalized
% watch order of prefixes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# step 1: extract onsets reeeealy fast if you havn't done that yet
 % condition names
logDir=''; % in case you keep them somewhere differnt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# step 2:  1st level glm for NOT-normalized data
% must include variable 'onsets' (cell with runs x conditions) including onset-times in sec

beta_dir = '1st_level_something'; % folder that will contain the created job.mat file and SPM file
% condnames = {'WM1', 'WM2', 'WM3', 'WM4', 'hit', 'miss', 'catch1', 'catch2'}; 
condnames = {'WM1', 'WM2', 'WM3', 'WM4'};
duration = 3; % epoch duration; for single events set 0
tr = TR;
fmri_t = n_slices; % Microtime resolution; If you have performed slice-timing correction, change this parameter to match the number of slices specified there; otherwise, set default 16
fmri_t0	= slice_order(round(length(slice_order)/2));
hpf      = 128; % High-pass filter cut-off; default 128 sec
% include multiple regressors (1=yes)
% if 1 (yes), 'hm' and/or 'cc' will be appended to outputfolder_1st
hm=0;   % head motion parameters from realignment (step 4 in B0_preprocessing)
cc=0;   % CompCorr WM and CSF principal components (step 8 in B0_preprocessing)
% do you want to normalize the realigned data too, to use masks in MNI-space for ROI-based decoding?
% if yes (1) nomalization will be initiated before the construction of the glm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# step 3:  contrasts 
analysisfolder = beta_dir;
cnames = condnames;
%very basic, adjust yourself, maybe model response regressor?
% cvecs  = {[ 1  0  0  0  0  0  0  0 ], ...   % 1
%           [ 0  1  0  0  0  0  0  0 ], ...   % 2
%           [ 0  0  1  0  0  0  0  0 ], ...   % 3
%           [ 0  0  0  1  0  0  0  0 ], ...   % 4
%           [ 0  0  0  0  1  0  0  0 ], ...   % 5
%           [ 0  0  0  0  0  1  0  0 ], ...   % 6
%           [ 0  0  0  0  0  0  1  0 ], ...   % 7
%           [ 0  0  0  0  0  0  0  1 ]};      % 8
cvecs = {[ 1  0  0  0 ], ...
         [ 0  1  0  0 ], ...
         [ 0  0  1  0 ], ...
         [ 0  0  0  1 ]};
del=1; % Delete existing contrasts (1=yes)
% were multiple regressors included in 1st level (step 1)?
n_hm=0;   % number of head motion parameters from realignment (step 4 in B0_preprocessing)
n_cc=0;   % number of CompCorr WM and CSF principal components (step 8 in B0_preprocessing)
    % if >0, "zeros" will be appended in design matrix
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# step 4:  2nd level
outputfolder_2nd = '2nd_level_something2';
dir_1st   =  beta_dir;
cnames_2nd = cnames;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# step 5:  2nd level flexfact
con_images=1:4; % bei 2xN: 1:9,b 2xA: 1:10
outputfolder_2nd_b = '2nd_level_FlexFact_something'; % folder that will contain the created SPM file
dir_1st_b   = beta_dir; % Name of corresponding first Level Analysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

currPrefix = start_prefix;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = analysis_switch

    switch n

        %% extract onsets
		case 1

			for s = 1:numel(SJs)

                if ismember(s, excludeSJ)
                    continue;
                else
                    subj_dir = fullfile(src_dir, SJs{s});
    			    cd(subj_dir)

    			    oruns=dir([SJs{s} '*.tsv']); %do they have names??
                    %oruns=dir(['fmri_log_*.tsv']);

    			    if ~isempty(runs)
        			    for r=1:length(runs)
            			    this_log = tdfread([subj_dir filesep oruns(r).name]);
%                             this_log = readtable([src_dir filesep SJs{s} filesep 'func' filesep 'fmri_log_run' num2str(r) '.txt']);

            			    % cond 1: WM level 1
            			    onsets{s,r,1} = this_log.Trial_Onset(this_log.Condition_n==1 & this_log.WM_delay==9000)/1000; %%% change to struct with onsets.irgendwas{
            			    % cond 2: WM level 2
            			    onsets{s,r,2} = this_log.Trial_Onset(this_log.Condition_n==2 & this_log.WM_delay==9000)/1000;
            			    % cond 3: WM level 3
            			    onsets{s,r,3} = this_log.Trial_Onset(this_log.Condition_n==3 & this_log.WM_delay==9000)/1000;
            			    % cond 4: WM level 4
            			    onsets{s,r,4} = this_log.Trial_Onset(this_log.Condition_n==4 & this_log.WM_delay==9000)/1000;

            			    % cond 5: hits
            			    onsets{s,r,5} = this_log.Trial_Onset(this_log.Hit==1 & this_log.WM_delay==9000)/1000;
            			    % cond 6: misses
            			    onsets{s,r,6} = this_log.Trial_Onset(this_log.Hit==0 & this_log.WM_delay==9000)/1000;

            			    % cond 7: catch trials or whatever you want to read out
            			    rando = randperm(12);
            			    catchs = this_log.Trial_Onset(this_log.WM_delay~=9000)/1000;
            			    onsets{s,r,7} = sort(catchs(rando([1:6])));
            			    onsets{s,r,8} = sort(catchs(rando([7:12])));

        			    end 
    			    end %end if
                end
			end

		%% glm
        case 2

            for sj = 1:numel(SJs)
                if ismember(sj, excludeSJ)
                    continue;
                else                
                    display(['Step 3, 1st level glm: ' SJs{sj} ])
                    subj_dir = fullfile(src_dir, SJs{sj});

                    %try
                        C2_glm_1stLevel(SJs, sj, subj_dir, beta_dir, currPrefix, tr, fmri_t, fmri_t0, hpf, runs, condnames, onsets, duration, hm, cc);
                %catch
                    %display('###################################################################')
                    %display(['################## ' SJs{sj} ', ERROR GLM ###################'])
                    %display('###################################################################')
                end
            end


		%% contrasts
		case 3

			for sj = 1:numel(SJs)

                if ismember(sj, excludeSJ)
                    continue;
                else

				    SJ_dir = [src_dir filesep SJs{sj}];
                    C3_contrast_1stLevel(SJ_dir, analysisfolder, cnames, cvecs, del, n_hm, n_cc, runs);

                end

            end

	    %% 2nd level
		case 4 

            SJin=SJs;
            SJin(excludeSJ)=[];
            C4_glm_2ndLevel_OneSampleTTest(src_dir, SJin, outputfolder_2nd, dir_1st, cnames_2nd);
			
        %% 2nd level flex fact
        case 5  

            SJin=SJs;
            SJin(excludeSJ)=[];
            C5_glm_2ndLevel_FlexFact(src_dir, SJin, outputfolder_2nd_b, dir_1st_b, con_images, runs);
    end
end
            