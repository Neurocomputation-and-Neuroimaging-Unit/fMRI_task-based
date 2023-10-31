% function B0_preprocessing_batch

%### step A, structure data before running the preprocessing
%### --> convert DICOM images to 4D NIFTI image IN BIDS FORMAT (ESSENTIAL)
%### --> subject folders (e.g. 'sub-002') including functional (e.g. 'func') and anatomy (e.g. 'anat') folders
%### --> optional instance for sessions
%############################################################################################

% required toolboxes (only download and add what you need):
% SPM12 ( http://www.fil.ion.ucl.ac.uk/spm/ )
% Rest_plus (1.8) ( http://restfmri.net/forum/index.php?q=rest )
% hMRI (https://www.cbs.mpg.de/departments/neurophysics/software/hmri-toolbox)
% Step 7  Scrubbing: BRAMILA toolbox ( https://users.aalto.fi/~eglerean/bramila.html )
% Step 8  CompCorr: dPABI toolbox ( http://rfmri.org/dpabi )
%#####################################################
%#################### INPUT ##########################
%#####################################################

clc
clear all
close all

addpath('C:\Users\saraw\Desktop\BIDS\task-based_fMRI_preprocessing-main')
% addpath(genpath('C:\Users\...\Toolboxes\bramila-master'))
% addpath(genpath('C:\Users\...\Toolboxes\REST_V1.8'))
% addpath(genpath('C:\Users\...\Toolboxes\RESTplus_v1.24')) fucks up normalization somehow
addpath(genpath('C:\Users\...\Toolboxes\hMRI-toolbox-0.4.0'))
% addpath(genpath('C:\Users\saraw\...\Toolboxes\DPABI_V5.1_201201'))
addpath('C:\Users\...\Toolboxes\spm12')

%SPM-path
SPM_path  = 'C:\Users\...\Toolboxes\spm12';

%data source directory
src_dir      = 'C:\Users\saraw\Desktop\BIDS\fmri';

%subject identifiers
cd(src_dir) 
pb=dir('sub*');
for i=1:length(pb)
    SJs(1,i)={pb(i).name};
end

excludeSJ = []; % keep empty if preprocessing is planned for all subjects, other wise just put them in by number: [4 6 9]
excludeRun = []; % same thing

%unzip
zip_files = dir(fullfile(src_dir, '**', ['sub-', '*.gz']));
if ~isempty(zip_files) 
	for z = 1:size(zip_files, 1)
		gunzip([zip_files(z).folder filesep zip_files(z).name]);
        delete([zip_files(z).folder filesep zip_files(z).name]);
	end
end

%session & run identifiers
sessNum = 0;
if exist([src_dir filesep SJs{1} filesep 'ses-1'])==7
    cd([src_dir filesep SJs{1}])
    sd = dir('ses*')
    sessNum = length(sd);
    for sess = 1:sessNum
        sessions(1, sess) = {sd(sess).name};
    end
    for sb = 1:numel(SJs)
        cd([src_dir filesep SJs{sb} filesep sessions{1} filesep 'func']);
        rd = dir('sub*.nii')
        for r = 1:length(rd)
			if ismember(r, excludeRuns)
                continue;
            else
            	runs(sb, r) = {rd(r).name};
			end
        end
    end
else
    for sb = 1:numel(SJs)
        cd([src_dir filesep SJs{sb} filesep 'func']);
        rd = dir('sub*.nii')
        for r = 1:length(rd)
			if ismember(r, excludeRuns)
                continue;
            else
            	runs(sb, r) = {rd(r).name};
			end
        end
    end
end

%anatomy identifier
ana=['anat'];

nifti_files = dir(fullfile(src_dir, '**', ['sub-', '*bold*.nii'])); %look for all functional nifti files
anat_files = dir(fullfile(src_dir, '**', ['sub-', '*T1w.nii'])); %look for all anat nifti files

%anatomical masks (for comp corr)
wm_mask=['C:\Users\...\wm_mask_eroded.nii']; %white matter mask file
csf_mask=['C:\Users\...\csf_mask_eroded.nii']; %csf mask file
full_brain_mask=['C:\Users\...\full_brain_mask.nii']; %full brain mask file

% selection of analysis steps (1-12) to be performed
analysis_switch = [7 10]; %1 4 3 5 7 8 9
start_prefix=''; %'gf4d'; %if totally raw data, then keep empty, otherwise add prefix, e.g. 's8wra'

%now we get the data from the json file 
json_files = (dir(fullfile(src_dir, '**', ['task', '*json']))); %extract all json files, althoguh they should have the same info 
%to check if the first command returns an empty structure. If yes, it means the json files
%have a different naming, starting with subject 
if isequal(size(json_files), [0, 1]) 
    json_files = (dir(fullfile(src_dir, '**', ['sub-', '*bold*.json'])));
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

%compare json and nifti header 
if round(TR_nifti, 4) ~= round(TR_json, 4)
    warning ("TR does not match between json file and nifti") 
end 
if n_slices_json ~= n_slices_nifti
    warning ("Number of slices does not match between json file and nifti") 
end 

TR = TR_json;
n_slices = n_slices_json; 


%# ------ NOTE: spike removal (e.g. "artrepair") should be performed as first step

%# step 1  Segmentation
%# ------ Create nuisance masks on your own or take the provided ones
%# step 2 --> remove first x scans                       --> prefix: x(number of cut volumes)
x=0;
%# step 3 --> slice time correction                      --> prefix: a
%  for interleaved slice order: do slice time correction, then realignment
%  otherwise do first realignment, then slice time correction (in analysis_switch 4 before 3)
% n_slices = 37; % number of slices
% slice_order=[1:n_slices];
refslice=slice_order(round(length(slice_order)/2)); % reference slice
% TR=2; % repetition time in sec.
%# step 4  Realignment                                --> prefix: r
% realign over all runs
%# step 5  Coregister (estimate) mean-epi 2 anatomy (DEFAULT)
corrPrefix = ''; % so if you perform differnt kinds of preprocessing, there will be 
% multiple 'mean...nii' files -> check which one you want and if theres anything 
% in between the 'mean' and 'sub-00...' then put that there (probably nothing or an a?)
%# step 5b  Coregister (estimate & resclice) mean-epi 2 anatomy --> prefix c
Co_er = 0; %default: 0, if 1, then estimate & reslice
%# step 6  Normalization                              --> prefix: w
vox_size=[2 2 2]; % preferred voxel size after Normailzation (in mm)
%vox_size=repmat(info.PixelDimensions(1),1,3); % Voxel size from JSON-file (not changed)
%# step 7  Scrubbing: calculate, interpolate outliers --> prefix: m(scrub_thresh)
scrub_thresh=0.4; % threshhold FD for scrubbing
%# step 8 Calculate WM and CSF Nuisance Signal
numComp = 5; % number of principle components
%# step 9 Smoothing                                   --> prefix: s
kernel_size=[8 8 8]; %FWHM kernel size
%# step 10 Detrending                                 --> prefix: d

%#####################################################
%#################### INPUT end ######################
%#####################################################

%%
currPrefix=start_prefix;

for n=analysis_switch
    
    switch n
        
        case 1 %% Segmentation
            warning off
            for sj = 1:numel(SJs)
                if ismember(sj, excludeSJ)
                    continue;
                elseif exist([src_dir filesep SJs{sj} filesep ana])==7
                    display(['Step 1, segmentation: ' SJs{sj}])
                    struct_dir = fullfile(src_dir, SJs{sj}, ana);
                    B1_segmentation(struct_dir, SJs{sj}, SPM_path, '^s.*\.nii');
                elseif sessNum > 0 && exist([src_dir filesep SJs{sj} filesep 'ses-1' filesep ana])==7
                    for ses = 1:sessNum
                        sesPath = [src_dir filesep SJs{sj} filesep 'ses-' num2str(ses) filesep ana];
                        display(['Step 1, segmentation: ' SJs{sj}])
                        struct_dir = fullfile(sesPath);
                        B1_segmentation(struct_dir, SJs{sj}, SPM_path, '^s.*\.nii');
                    end
                else
                    display('###########################################################')
                    display(['############### ' SJs{sj} ', ' ana ' does not exsist ###########'])
                end
            end
            
        case 2 %% Delete first X scans
            if x>0
                for sj = 1:numel(SJs)
                    if ismember(sj, excludeSJ)
                        continue;
                    else
                        for r = 1:size(runs, 2)
                            if exist([src_dir filesep SJs{sj} filesep runs{sj, r}])
                                display(['Step 2, delete first ' num2str(x) ' volumes: ' SJs{sj} ', ' runs{sj, r}])
                                run_dir = fullfile(src_dir, SJs{sj}, 'func');
                                B2_delete_scans(run_dir, ['^' currPrefix runs{sj, r}],x);
                            elseif sessNum > 0 && exist([run_dir filesep SJs{sj} filesep 'ses-1' filesep runs{sj, r}])
                                for ses = 1:sessNum
                                    display(['Step 2, delete first ' num2str(x) ' volumes: ' SJs{sj} ', ' runs{sj, r}])
                                    sesPath = [src_dir filesep SJs{sj} filesep 'ses-' num2str(ses)];
                                    run_dir = fullfile(sesPath, 'func');
                                    B2_delete_scans(run_dir, ['^' currPrefix runs{sj, r}],x);
                                end
                            else
                                display('###########################################################')
                                display(['############### ' SJs{sj} ', ' runs{sj, r} ' does not exsist ###########'])
                            end
                        end
                    end
                end
                currPrefix=['x' num2str(x) currPrefix];
            end
            
        case 3 %% Slice time correction
            for sj = 1:numel(SJs)
                if ismember(sj, excludeSJ)
                    continue;
                else
                    for r = 1:size(runs, 2)
                        if exist([src_dir filesep SJs{sj} filesep 'func' filesep runs{sj, r}])
                            display(['Step 3, slice time correction: ' SJs{sj} ', ' runs{sj, r}])
                            funcPath = [src_dir filesep SJs{sj} filesep 'func'];
                            run_dir = fullfile(funcPath);
                            B3_slice_time_correction(SJs{sj},runs{sj, r}, run_dir, ['^' currPrefix runs{sj, r}],n_slices,slice_order,refslice,TR);
                        elseif sessNum > 0 && exist([src_dir filesep SJs{sj} filesep 'ses-1' filesep 'func' filesep runs{sj, r}])
                            for ses = 1:sessNum
                                display(['Step 3, slice time correction: ' SJs{sj} ', ' runs{sj, r}])
                                sesPath = [src_dir filesep SJs{sj} filesep 'ses-' num2str(ses) filesep 'func'];
                                run_dir = fullfile(sesPath);
                                B3_slice_time_correction(SJs{sj},runs{sj, r},run_dir, ['^' currPrefix runs{sj, r}],n_slices,slice_order,refslice,TR);
                            end
                        else
                            display('###########################################################')
                            display(['############### ' SJs{sj} ', ' runs{sj, r} ' does not exsist ###########'])
                        end
                    end
                end
            end
            currPrefix=['a' currPrefix];
            
        case 4 %% Realignment
            for sj = 1:numel(SJs)
                if ismember(sj, excludeSJ)
                    continue;
                else
                    sj_dir = [src_dir filesep SJs{sj}];
                        if exist([src_dir filesep SJs{sj} filesep 'func' filesep runs{sj, r}])
                            display(['Step 4, realignment: ' SJs{sj} ', ' runs{sj, r}])
                            funcPath = [src_dir filesep SJs{sj} filesep 'func'];
                            run_dir = fullfile(funcPath);
                            for r = 1:size(runs, 2)
                                run_files{r} = spm_select('List',run_dir,['^' currPrefix runs{sj, r}]);
                            end
                            B4_Realignment_all_runs(run_dir, run_files);
                            
                        elseif sessNum > 0 && exist([src_dir filesep SJs{sj} filesep 'ses-1' filesep 'func' filesep runs{sj, r}])
                                for ses = 1:sessNum
                                    display(['Step 4, realignment: ' SJs{sj} ', ' runs{sj, r}])
                                    sesPath = [src_dir filesep SJs{sj} filesep 'ses-' num2str(ses) filesep 'func'];
                                    run_dir = fullfile(sesPath);
                                    for r = 1:size(runs, 2)
                                        run_files{r} = spm_select('List',run_dir,['^' currPrefix runs{sj, r}]);
                                    end                                    
                                    B4_Realignment_all_runs(sesPath, run_files);
								end
                        end
                        else
                            display('###########################################################')
                            display(['############### ' SJs{sj} ', ' runs{sj, r} ' does not exsist ###########'])
                        end
                    %end
                end
            end
            currPrefix=['r' currPrefix]; %%%%%%%%%%%%%% fix
            
        case 5 %% Coregister (estimate) mean-epi 2 anatomy %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             warning off
            for sj = 1:numel(SJs)
                if ismember(sj, excludeSJ)
                    continue;
                else
%                     for r = 1:size(runs, 2)
                        if exist([src_dir filesep SJs{sj} filesep 'func'])
                            
                            funcPath = [src_dir filesep SJs{sj}];
                            func_dir        = fullfile(funcPath, 'func');
                            struct_dir      = fullfile(funcPath, ana);
							if Co_er ~= 1
								display(['Step 5, coregistration (estimate): ' SJs{sj}])
                            	B5_coregister_est(currPrefix, func_dir, struct_dir, sj, '^s.*\.nii', runs, corrPrefix);
							else
								display(['Step 5b, coregistration (estimate & reslice): ' SJs{sj}])
								B5b_coregister_est_re(currPrefix, func_dir, struct_dir, SJs{sj}, '^s.*\.nii', runs);
							end
                        elseif sessNum > 0 && exist([src_dir filesep SJs{sj} filesep 'ses-1' filesep 'func'])
                                for ses = 1:sessNum
                                    display(['Step 5, coregistration: ' SJs{sj}])
                                    sesPath = [src_dir filesep SJs{sj} filesep 'ses-' num2str(ses)];
                                    func_dir        = fullfile(sesPath, 'func');
                                    struct_dir      = fullfile(sesPath, ana);
									if Co_er ~= 1
										display(['Step 5, coregistration (estimate): ' SJs{sj}])
                                    	B5_coregister_est(currPrefix, func_dir, struct_dir, sj, '^s.*\.nii', runs, corrPrefix);
									else
										display(['Step 5b, coregistration (estimate & reslice): ' SJs{sj}])
										B5b_coregister_est_re(currPrefix, func_dir, struct_dir, SJs{sj}, '^s.*\.nii', runs);
									end
                                end
                        else
                            display('###########################################################')
                            display(['############### ' SJs{sj} 's functional data do not exsist ###########'])
                        end
%                     end
                end
            end
            
        case 6 %% Normalization
            for sj = 1:numel(SJs)
                if ismember(sj, excludeSJ)
                    continue;
                else
                    
%                     for r = 1:size(runs, 2)
                        if exist([src_dir filesep SJs{sj} filesep 'func'])
                            display(['Step 6, normalization: ' SJs{sj}])
                            funcPath = [src_dir filesep SJs{sj}];
                            struct_dir = fullfile(funcPath, ana);
                            data_dir = fullfile(funcPath, 'func');
                            B6_normalization_run(data_dir, struct_dir, sj, runs, vox_size, currPrefix);
                        elseif sessNum > 0 && exist([src_dir filesep SJs{sj} filesep 'ses-1' filesep 'func'])
                                for ses = 1:sessNum
                                    display(['Step 6, normalization: ' SJs{sj}])
                                    sesPath = [src_dir filesep SJs{sj} filesep 'ses-' num2str(ses)];
                                    func_dir        = fullfile(sesPath, 'func');
                                    struct_dir      = fullfile(sesPath, ana);
                                    B6_normalization_run(func_dir, struct_dir, sj, runs, vox_size, currPrefix);
                                end
                        else
                            display('###########################################################')
                            display(['############### ' SJs{sj} ', ' runs{sj, r} ' does not exsist ###########'])
                        end
%                     end
                end
            end
            currPrefix=['w' currPrefix];
            
        case 7 %% Scrubbing: calculate outliers
            scrub_prefix=['m' num2str(scrub_thresh)];
            for sj = 1:numel(SJs)
                if ismember(sj, excludeSJ)
                    continue;
                else
                    for r = 1:size(runs, 2)
                        if exist([src_dir filesep SJs{sj} filesep 'func' filesep runs{sj, r}])
                            display(['Step 7, scrubbing: ' SJs{sj} ', ' runs{sj, r}])
                            funcPath = [src_dir filesep SJs{sj}];
                            data_dir = fullfile(funcPath, 'func');
                            %estimate and save motion statistics
                            n=1;
                            k = strfind(runs{sj, r}, 'd.nii'); %% use prefix without file extension
                            f=spm_select('List', data_dir, ['^rp_' currPrefix(n:end) runs{sj, r}(1:k) '.txt']);
                            while isempty(f)
                                n=n+1;
                                f=spm_select('List', data_dir, ['^rp_' currPrefix(n:end) runs{sj, r}(1:k) '.txt']);
                            end
                            cfg.motionparam=[data_dir filesep f];
                            cfg.prepro_suite = 'spm';
                            [fwd,rms]=bramila_framewiseDisplacement(cfg);
                            outliers=fwd>scrub_thresh;
                            percent_out=(sum(outliers)/length(outliers))*100;
                            disp(['outliers for ' num2str(SJs{sj}) ', ' runs{r} ': ' num2str(percent_out) '%']);
                            save([data_dir filesep scrub_prefix currPrefix runs{sj, r}(1:k) '_FWDstat.mat'],'fwd','rms','outliers','percent_out','scrub_thresh','cfg')
                            %srub outliers by replacing them with average of nearest neighbors
                            B7_scrub_data(data_dir, ['^' currPrefix runs{sj, r}], outliers,  scrub_prefix);
                            all_percent_out(sj,r)=percent_out;
                            all_rp{sj,r}=load(cfg.motionparam);

                        elseif exist([src_dir filesep SJs{sj} filesep 'ses-1' filesep 'func' filesep runs{sj, r}])
                            for ses = 1:sessNum
                                display(['Step 7, scrubbing: ' SJs{sj} ', ' runs{sj, r}])
                                funcPath = [src_dir filesep SJs{sj} filesep 'ses-' num2str(ses)];
                                data_dir = fullfile(funcPath, 'func');
                                %estimate and save motion statistics
                                n=1;
								k = strfind(runs{sj, r}, 'd.nii'); %% use prefix without file extension
                                f=spm_select('List', data_dir, ['^rp_' currPrefix(n:end) runs{sj, r}(1:k) '.*\.txt']);
                                while isempty(f)
                                    n=n+1;
                                    f=spm_select('List', data_dir, ['^rp_' currPrefix(n:end) runs{sj, r}(1:k) '.*\.txt']);
                                end
                                cfg.motionparam=[data_dir filesep f];
                                cfg.prepro_suite = 'spm';
                                [fwd,rms]=bramila_framewiseDisplacement(cfg);
                                outliers=fwd>scrub_thresh;
                                percent_out=(sum(outliers)/length(outliers))*100;
                                disp(['outliers for ' num2str(SJs{sj}) ', ' runs{sj, r} ': ' num2str(percent_out) '%']);
                                save([data_dir filesep scrub_prefix currPrefix runs{sj, r}(1:k) '_FWDstat.mat'],'fwd','rms','outliers','percent_out','scrub_thresh','cfg')
                                %srub outliers by replacing them with average of nearest neighbors
                                B7_scrub_data(data_dir, ['^' currPrefix runs{sj, r}], outliers,  scrub_prefix);
                                all_percent_out(sj,r)=percent_out;
                                all_rp{sj,r}=load(cfg.motionparam);
                            end
                        else
                            display('###########################################################')
                            display(['############### ' SJs{sj} ', ' runs{sj, r} ' does not exsist ###########'])
                        end
                    end
                end
            end
            currPrefix=[scrub_prefix currPrefix];
            save([src_dir filesep 'all_MOTIONstat_' currPrefix '.mat'],'SJs','runs','scrub_thresh','all_percent_out','all_rp')
            
        case 8 %% CompCorr
            for sj = 1:numel(SJs)
                if ismember(sj, excludeSJ)
                    continue;
                else
                    for r = 1:size(runs, 2)
                        if exist([src_dir filesep SJs{sj} filesep 'func' filesep runs{sj, r}])
                            display(['Step 8, CompCorr: ' SJs{sj} ', ' runs{sj, r}])
                            funcPath = [src_dir filesep SJs{sj}];
                            data_dir = fullfile(funcPath, 'func');
                            B8_compcorr_run(data_dir, SJs{sj}, ['^' currPrefix runs{sj, r}], numComp, wm_mask, csf_mask);
                        elseif exist([src_dir filesep SJs{sj} filesep 'ses-1' filesep 'func' filesep runs{sj, r}])
                            for ses = 1:sessNum
                                    display(['Step 8, CompCorr: ' SJs{sj} ', ' runs{sj, r}])
                                    sesPath = [src_dir filesep SJs{sj} filesep 'ses-' num2str(ses) filesep 'func'];
                                    data_dir = fullfile(sesPath);
                                    B8_compcorr_run(data_dir, SJs{sj}, ['^' currPrefix runs{sj, r}], numComp, wm_mask, csf_mask);
                            end
                        else
                            display('###########################################################')
                            display(['############### ' SJs{sj} ', ' runs{sj, r} ' does not exsist ###########'])
                        end
                    end
                end
            end
            
        case 9 %% Smoothing
            for sj = 1:numel(SJs)
                if ismember(sj, excludeSJ)
                    continue;
                else
                    for r = 1:size(runs, 2)
                        if exist([src_dir filesep SJs{sj} filesep 'func' filesep runs{sj, r}])
                            display(['Step 9, smoothing: ' SJs{sj} ', ' runs{sj, r}])
                            funcPath = [src_dir filesep SJs{sj}];
                            run_dir = fullfile(funcPath, 'func');
                            B9_smoothing_run(run_dir, SJs{sj}, ['^' currPrefix runs{sj, r}],kernel_size);
                            display([SJs{sj} ', ' runs{r} ' is done'])
                        elseif exist([src_dir filesep SJs{sj} filesep 'ses-1' filesep 'func' filesep runs{sj, r}])
                            for ses = 1:sessNum
                                    display(['Step 9, smoothing: ' SJs{sj} ', ' runs{sj, r}])
                                    sesPath = [src_dir filesep SJs{sj} filesep 'ses-' num2str(ses) filesep 'func'];
                                    run_dir = fullfile(sesPath);
                                    B9_smoothing_run(run_dir, SJs{sj}, ['^' currPrefix runs{sj, r}],kernel_size);
                            end
                        else
                            display('###########################################################')
                            display(['############### ' SJs{sj} ', ' runs{sj, r} ' does not exsist ###########'])
                        end
                    end
                end
            end
            currPrefix=['s' num2str(unique(kernel_size)) currPrefix];
            
        case 10 %% Detrending
            for sj = 1:numel(SJs)
                if ismember(sj, excludeSJ)
                    continue;
                else
                    display(['Step 10, Detrending: ' SJs{sj}])
                    for r = 1:size(runs, 2)
                        if exist([src_dir filesep SJs{sj} filesep 'func' filesep runs{sj, r}])
                            funcPath = [src_dir filesep SJs{sj}];
                            run_dir = fullfile(funcPath, 'func');
                            f = spm_select('List',run_dir, ['^' currPrefix runs{sj, r}]);
                            V=spm_vol([run_dir filesep f(1,:)]);
                            files={};
                            for i=1:(size(V,1))
                                files{i} = [run_dir filesep strtrim(f(1,:)) ',' int2str(i)];
                            end
                            fileset{r}=char(files);

                        elseif exist([src_dir filesep SJs{sj} filesep 'ses-1' filesep 'func' filesep runs{sj, r}])
                            for ses = 1:sessNum
                                sessPath = [src_dir filesep SJs{sj} filesep 'ses-' num2str(ses)];
                                run_dir = fullfile(sessPath, 'func');
                                f = spm_select('List',run_dir, ['^' currPrefix runs{sj, r}]);
                                V=spm_vol([run_dir filesep f(1,:)]);
                                files={};
                                for i=1:(size(V,1))
                                    files{i} = [run_dir filesep strtrim(f(1,:)) ',' int2str(i)];
                                end
                                fileset{r}=char(files);
                            end
                        else
                            display('###########################################################')
                            display(['############### ' SJs{sj} ', ' runs{sj, r} ' does not exsist ###########'])
                        end
                    end
                    if exist('fileset')
                        B10_detrending_lmgs(fileset);
                        clear fileset files
                    end
                    display([SJs{sj} ', is done'])
                end
            end
        
        otherwise
            display('######################################################################################')
            display(['############################## Case ' num2str(n) ' does not exsist ##############################'])
            display('######################################################################################')
    end
end
%%
% adjust the following to delete files with certain prefixes from the subject
% folders
% if ismember(99,analysis_switch)
%     for sj = sbjs 1:numel(SJs)
%         for runnr = 1:numel(runs)
%             data_dir = fullfile(src_dir, SJs{sj}, ['run0', num2str(runnr)]);
%             cd(data_dir)
%             delete('r*.mat')
%             delete('f*.mat')
%             delete('wraf4d*.mat')
%             delete('ds8wraf4d*.mat')
%             delete('ds8m0.4wraf4d*.mat')
%             delete('m*.mat')
%             delete('a*.mat')
% %             delete('o*.mat')
% %             delete('d*.mat')
% %             delete('s*.mat')
% %             delete('w*.mat')
%         end
%     end
% end


