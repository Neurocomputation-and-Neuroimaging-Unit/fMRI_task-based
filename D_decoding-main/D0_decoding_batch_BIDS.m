%% decoding BATCH %%

% required toolboxes:
% the decoding toolbox TDT
% https://doi.org/10.3389/fninf.2014.00088

%#####################################################
%#################### INPUT ##########################
%#####################################################


%SPM-path
SPM_path  = 'C:\Users\...\Toolboxes\spm12';

%data source directory
src_dir      = 'C:\Users\saraw\Desktop\BIDS\fmri';

addpath(genpath('C:\Users\saraw\Desktop\BIDS')); %%%%%%% your current working directory %%%%%%%%%%%%%% eg: 'F:\GForce\Dataanalysis\D_decoding_BIDS-main'
addpath(genpath('C:\Users\...\Toolboxes\hMRI-toolbox-0.4.0'));
addpath('C:\Users\...\Toolboxes\decoding_toolbox');
addpath('C:\Users\...\Toolboxes\spm12');

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
excludeSJ = [1] % delete those sjs

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
start_prefix='r'; 
%'' when bids and unpreprocessed; 
%ra/ar if already slicetime-corrected and realigned
%'' if doing normilazation of accuracy maps
%w if aready normalized
% watch order of prefixes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# step 1: extract onsets reeeealy fast if you havn't done that yet
 % condition names
logDir=''; % in case you keep them somewhere differnt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# step 2.a:  compute FIR if needed
WMdelay = 9; %in secs
condnamesFIR = {'WM1', 'WM2', 'WM3', 'WM4'};
refslice = slice_order(round(length(slice_order)/2));
fir_out = [num2str(WMdelay/TR) '_bins'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# step 3: 1st level glm for NOT-normalized data
% must include variable 'onsets' (cell with sj x runs x conditions) including onset-times in sec

condnames = {'WM1', 'WM2', 'WM3', 'WM4', 'hit', 'miss', 'catch1', 'catch2'}; 
duration = 3; % epoch duration; for single events set 0
tr = TR;
fmri_t = n_slices; % Microtime resolution; If you have performed slice-timing correction, change this parameter to match the number of slices specified there; otherwise, set default 16
fmri_t0	= refslice;
hpf      = 128; % High-pass filter cut-off; default 128 sec
% include multiple regressors (1=yes)
% if 1 (yes), 'hm' and/or 'cc' will be appended to outputfolder_1st
hm=0;   % head motion parameters from realignment (step 4 in B0_preprocessing)
cc=0;   % CompCorr WM and CSF principal components (step 8 in B0_preprocessing)
% do you want to normalize the realigned data too, to use masks in MNI-space for ROI-based decoding?
% if yes (1) nomalization will be initiated before the construction of the glm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# step 4:  decode 

betas = 'C:\Users\saraw\Desktop\BIDS\fmri\sub-001\TRSLD_FIR';

dec_type = "searchlight"; % searchlight or roi
dec_task = "regression"; % regression (SVR) and classification (SVM) are implemented so far

labelnames 	= 	{'WM1 bin 1', 'WM1 bin 2', 'WM1 bin 3', 'WM1 bin 4', 'WM1 bin 5', 'WM1 bin 6'; ...
                 'WM2 bin 1', 'WM2 bin 2', 'WM2 bin 3', 'WM2 bin 4', 'WM2 bin 5', 'WM2 bin 6'; ...
                 'WM3 bin 1', 'WM3 bin 2', 'WM3 bin 3', 'WM3 bin 4', 'WM3 bin 5', 'WM3 bin 6'; ...
                 'WM4 bin 1', 'WM4 bin 2', 'WM4 bin 3', 'WM4 bin 4', 'WM4 bin 5', 'WM4 bin 6'}; 

Nbins = WMdelay/TR; %default
con_pairs = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];

outputDir4 = 'D4';
sbj_level_folder = 'DEC';

rad = 3; % radius of searchlight

if dec_type == "roi" 
	cd(''); % dir that includes all roi-masks (binary!!)
	masks = dir(['*.nii']); %or whatever you use as an identifier
else
    cd(src_dir);
    masks = dir('brain_mask.nii');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# step 5:  average over conditions
outputDir5 = 'mean_';
% (9s at TR = 1,5 equals 6 Intervals; 12s at TR = 2 equals 6 Intervals)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# step 6:  normalize

% WARNING % 
% --- 
% T1 MUST HAVE BEEN SEGMENTED AND FUNCTIONAL DATA MUST HAVE
% BEEN REALIGNED (see preprocessing Batch B0 for these steps) 
% --- 
% WARNING %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%# step 7:  smooth
s_kernel = [5 5 5];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


currPrefix=start_prefix;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = analysis_switch

    switch n
	
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

		%% FIR
		case 2
		
            D2_FIR(src_dir, onsets, SJs, TR, WMdelay, currPrefix, runs, excludeSJ, fir_out)

		%% glm
        case 3
		
			beta_dir = ['1st_level_D0_' currPrefix]; % folder that will contain the created job.mat file and SPM file
            
            for sj = 1:numel(SJs)
                if ismember(sj, excludeSJ)
                    continue;
                else                
                    display(['Step 3, 1st level glm: ' SJs{sj} ])
                    subj_dir = fullfile(src_dir, SJs{sj});
				    
                    %try
                        D3_glm_1stLevel(SJs, sj, subj_dir, beta_dir, currPrefix, tr, fmri_t, fmri_t0, hpf, runs, condnames, onsets, duration, hm, cc);
                %catch
                    %display('###################################################################')
                    %display(['################## ' SJs{sj} ', ERROR GLM ###################'])
                    %display('###################################################################')
                end
            end

            
		%% decoding
		case 4
		
			for sj = 1:numel(SJs)

                if ismember(sj, excludeSJ)
                    continue;
                else
				
				    SJ_dir = [src_dir filesep SJs{sj}];
                    betas = [SJ_dir filesep ['FIR_' fir_out]];
                    dec_folder = [SJ_dir filesep sbj_level_folder];
                    if ~exist(dec_folder, 'dir')
                        mkdir(dec_folder);
                    end
                       
                    for b = 1:Nbins


                        if dec_task == "classification"

                            for cp = 1:size(con_pairs,1) % maybe adjust to allow restriction of used betas beforehand

                                label1 = labelnames(con_pairs(cp,1),b);
                                label2 = labelnames(con_pairs(cp,2),b);
                                these_labelnames = [label1, label2];

                                display(['Step 4, Decoding: ' SJs{sj} ', Conditions: ' cell2mat(label1) ', ' cell2mat(label2)])
                                beta_path    = fullfile(betas);
                                output_path  = fullfile(dec_folder, [outputDir4 '_' cell2mat(label1) '_' cell2mat(label2)]);
                                D4class_Decoding_batch(beta_path, output_path, these_labelnames, dec_type, rad, masks); %% updated D4, no further readjustments within

                                clear label1 label2 these_labelnames

                            end

                        elseif dec_task == "regression"

                            beta_path    = fullfile(betas);
                            output_path  = fullfile(dec_folder, [outputDir4 '_bin_' num2str(b)]);
                            these_labelnames = labelnames(:,b)';
                            D4regress_Decoding_batch(beta_path, output_path, these_labelnames, dec_type, rad, masks); %% updated D4, no further readjustments within
                            clear these_labelnames

                        end

                    end

                end
			
            end
            
	    %% average
		case 5 %% ADJUST 
		
			for sj = 1:numel(SJs)

                if ismember(sj, excludeSJ)
                    continue;
                else

                    sj_dir = [src_dir filesep SJs{sj}];
                    dec_folder = [SJ_dir filesep sbj_level_folder];
    
                    for b = 1:Nbins %bins
    
                        bin_folders = spm_select('List', dec_folder, 'dir', [' bin ' num2str(b)]);
                        
    
                        matlabbatch{1}.spm.util.imcalc.input = {strcat(dec_folder, filesep, bin_folders(1,:), filesep, 'res_accuracy_minus_chance.nii')
                                                                strcat(dec_folder, filesep, bin_folders(2,:), filesep, 'res_accuracy_minus_chance.nii')
                                                                strcat(dec_folder, filesep, bin_folders(3,:), filesep, 'res_accuracy_minus_chance.nii')
                                                                strcat(dec_folder, filesep, bin_folders(4,:), filesep, 'res_accuracy_minus_chance.nii')
                                                                strcat(dec_folder, filesep, bin_folders(5,:), filesep, 'res_accuracy_minus_chance.nii')
                                                                strcat(dec_folder, filesep, bin_folders(6,:), filesep, 'res_accuracy_minus_chance.nii')};
    
                        matlabbatch{1}.spm.util.imcalc.expression = '(i1 + i2 + i3 + i4 + i5 + i6) / 6';
    
        
                        matlabbatch{1}.spm.util.imcalc.output = 'res_accuracy_minus_chance.nii';
        
        %               matlabbatch{1}.spm.util.imcalc.output = 's5wres_accuracy_minus_chance.nii';
                                
               	        if ~exist(strcat(dec_folder, filesep, outputDir5, ['_bin_' num2str(b)]), 'dir')
                            mkdir(strcat(dec_folder, filesep, outputDir5, ['_bin_' num2str(b)]));
                        end
                                
                        matlabbatch{1}.spm.util.imcalc.outdir = {strcat(dec_folder, filesep, outputDir5, ['_bin_' num2str(b)])};
        
                        matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
                        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
                        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
                        matlabbatch{1}.spm.util.imcalc.options.interp = 1;
                        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
        
                        spm('defaults', 'FMRI');
                        spm_jobman('run', matlabbatch);
    
                    end
				    
                end
            end

        case 6 %% norm 

            for sj = 1:numel(SJs)

                if ismember(sj, excludeSJ)
                    continue;
                else

                    sj_dir = [src_dir filesep SJs{sj}];
                    func_dir = [sj_dir filesep 'func'];
                    struct_dir = [sj_dir filesep 'anat'];

                    for b = 1:Nbins

                        data_dir = [sj_dir filesep sbj_level_folder filesep 'mean_bin_' num2str(b)];
                        f3 = spm_select('List', data_dir, '^res_accuracy_minus_chance.nii');
                        if isempty(f3)
                            data_dir = [sj_dir filesep sbj_level_folder filesep [outputDir4 '_bin_' num2str(b)]];
                            f3 = spm_select('List', data_dir, '^res_zcorr.nii');
                        end
                        numVols = size(f3,1);
                        Images{b} = cellstr([repmat([data_dir filesep], numVols, 1) f3 repmat(',1', numVols, 1)]);

                    end
                    
                    D6a_coregister_est(func_dir, struct_dir, Images);
                    D6b_segmentation(struct_dir, SJs{sj}, SPM_path, '^s.*\.nii');
                    D6c_normalization(Images, struct_dir, vox_size);
                        
                end

            end
            currPrefix = 'w';

        case 7 %% smooth

            for sj = 1:numel(SJs)

                if ismember(sj, excludeSJ)
                    continue;
                else

                    sj_dir = [src_dir filesep SJs{sj}];
                    struct_dir = [sj_dir filesep 'anat'];

                    for b = 1:Nbins
                        
                        if dec_task == "classification"
                            data_dir = [sj_dir filesep sbj_level_folder filesep 'mean_bin_' num2str(b)];
                            D7_smoothing_run(data_dir, SJ{sj}, ['^' currPrefix 'res_accuracy_minus_chance.nii'], s_kernel);
                        elseif dec_task == "regression"
                            data_dir = [sj_dir filesep sbj_level_folder filesep [outputDir4 '_bin_' num2str(b)]];
                            D7_smoothing_run(data_dir, SJ{sj}, ['^' currPrefix 'res_zscore.nii'], s_kernel);
                        end
                        
                    end

                end

            end

    end
end
