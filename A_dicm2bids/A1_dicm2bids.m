function A1_dicm2bids(dcmDir, niiDir, varargin)
% This function calls dicm2nii.m to convert dicom into nii under niiDir. 
% Then moves nii files into BIDS style folders. The unknown series will stay in
% subject folder.
%
% To teach the converter for those unknown series, do the following:
% 
%   tbl = table('Size', [0 3], ...
%     'VariableTypes', {'cellstr' 'categorical' 'categorical'}, ...
%     'VariableNames', {'Name' 'Type' 'Modality'});
%   tbl = getpref('dicm2nii_gui_para', 'ModalityTable', tbl);
% 
% This will load the existing modality table or initilaize an empty one. Then
% you can append a row like:
%   tbl(end+1,:) = {'EPI_correction_PA' 'fmap' 'fieldmap'};
% 
% The three components are {ProtocalName subFolder modality} for the series.
% You can also edit the table by
%  open tbl
% After done, save the table by 
%  setpref('dicm2nii_gui_para', 'ModalityTable', tbl);

% 210530 XiangruiLi at gmail.com wrote it for CCBBI dataset

% types = {'anat' 'dwi' 'fmap' 'func' 'perf'};
% modality = {'FLAIR' 'PD' 'PDmap' 'T1map' 'T1rho' 'T1w' 'T2map' ...
%     'T2star' 'T2w' 'dwi' 'fieldmap' 'magnitude1' 'magnitude2' ...
%     'phase1' 'phase2' 'phasediff' 'task-motor' 'task-rest'};
% suffix = {'_sbref' '_bold'};

% perf/sub-_acq-_rec-_dir-_run-_asl/m0scan/aslcontext/asllabeling.jpg
% fmap/sub-_ses-_acq-_ce-_dir-_run-_epi
% func/sub-_ses-_task-_acq-_ce-_rec-_dir-_run-_echo-_part-_bold.nii[.gz]
% dwi/sub-_ses-_acq-_dir-_run-_part-_dwi

if nargin<2 || isempty(niiDir), niiDir = './'; end
if ~endsWith(niiDir, {'\' '/'}), niiDir = [niiDir filesep]; end
save_json = getpref('dicm2nii_gui_para', 'save_json', false);
setpref('dicm2nii_gui_para', 'save_json', true);
dicm2nii(dcmDir, niiDir); 
setpref('dicm2nii_gui_para', 'save_json', save_json); % restore it

load([niiDir 'dcmHeaders.mat'], 'h');
NIIs = fieldnames(h);
s = h.(NIIs{1});

if nargin > 2 %% allow batch-skript with number of subject in case of many subjects

    subj = varargin(1);
    subj = subj{1,1};
    
    try
        if subj < 10
            subj = ['00' num2str(subj)];
        elseif subj < 100
            subj = ['0' num2str(subj)];
        else
            subj = num2str(subj);
        end        
    end

    fSubj = [niiDir 'sub-' subj filesep];

    if nargin > 3 %% allow manual naming of tasks in batch skript

        task = varargin(2);
        task = task{1,1};

        if nargin > 4
            
            sess = varargin(3);
            sess = sess{1,1};

            if ischar(sess)
                fSubj = [niiDir 'sub-' subj filesep 'sub-' subj '_ses-' sess];
            end

            if nargin > 5
                
                ses = varargin(4);
                ses = ses{1,1};

                if nargin > 6

                    skip_sbref = varargin(5);
                    skip_sbref = skip_sbref{1,1};

                end

            end

        end

    end

else

    subj = regexprep(s.PatientName, '[^0-9a-zA-Z]', ''); 
    fSubj = [niiDir 'sub-' subj filesep];

end

if ~exist(fSubj, 'dir'), mkdir(fSubj); end
% csvFile = [subj '_' s.AcquisitionDate(3:8) '_seriesInfo.csv'];
% copyfile(['/mnt/fmri/downloads/' PI_info(s) '/' csvFile], [fSubj '.']);

tbl = table('Size', [0 3], ...
    'VariableTypes', {'cellstr' 'categorical' 'categorical'}, ...
    'VariableNames', {'Name' 'Type' 'Modality'});
tbl = getpref('dicm2nii_gui_para', 'ModalityTable', tbl);

N = numel(NIIs);
types = repmat({''}, N, 1); % empty mean to skip
names = repmat({''}, N, 1);
try
    if sess(1) == 1 || sess(2) > 1
        sessions = repmat({num2str(ses)}, N, 1);
    end
catch
    warning('no number of sessions assigned! check for possible errors.')
end
PEdir = repmat({''}, N, 1);
fNum = 0;
for i = 1:N
    s = h.(NIIs{i});
    if ~isfield(s, 'UnwarpDirection') % suppose axial slices
    elseif s.UnwarpDirection == "-y", PEdir{i} = 'dir-AP';
    elseif s.UnwarpDirection ==  "y", PEdir{i} = 'dir-PA';
    elseif s.UnwarpDirection ==  "x", PEdir{i} = 'dir-LR';
    elseif s.UnwarpDirection == "-x", PEdir{i} = 'dir-RL';
    end
    seqContains = @(p)any(cellfun(@(c)contains(s.SequenceName,c), p));
    j = find(strcmp(s.ProtocolName, tbl.Name), 1, 'last');
    if ~isempty(j) % user-provided tbl takes precedence
        names{i} = char(tbl.Modality(j));
        types{i} = char(tbl.Type(j));
        if types{i} == "skip", types{i}  = ''; end
    elseif s.isDTI
        types{i} = 'dwi';
        if endsWith(s.SeriesDescription, '_SBRef')
            names{i} = [names{i} '_sbref'];
        else
            names{i} = [names{i} '_dwi'];
        end
    elseif seqContains({'*fl3d1_ns' '*fl2d1' 'ABCD3d1_32ns'}) || ...
            contains(s.ImageType, 'DERIVED\')
        continue;
    elseif contains(s.ImageType, '\ASL\')
        types{i} = 'perf'; names{i} = 'asl';
    elseif seqContains({'epfid2d'})
        prot = s.SeriesDescription;
        types{i} = 'func';
        fNum = fNum+1;
        if contains(prot, 'rest', 'IgnoreCase', true)
            names{i} = 'task-rest';
        else
            c = regexpi(prot, '(.*)?run[-_]*(\d*)', 'tokens', 'once'); %%% verbessern, sodass task und run funktionieren und nur dingens, hier, wie heißts, _ greifen
            if isempty(c)
                if nargin > 3 %% as task might be reqired but not found in protocol
                    c = [task '_run-' num2str(fNum)];
                    warning('Used user-input for task-value and number of file for run-value as protocol was not readable! Please check input and output data for possible mistakes.');
                else
                    c = ['_run-' num2str(fNum)];
                    warning('Used number of file for run-value as protocol was not readable! Please check input and output data for possible mistakes.');
                end
%                 c = '';
%                 c = regexprep(prot, '[^0-9a-zA-Z]', '');
            else
                c = regexprep(c, '[^0-9a-zA-Z]', '');
                c = sprintf('%s_run-%s', c{:});
            end
            names{i} = ['task-' c]; 
        end
        if sess(2) > 1

        end
        ec = regexp(NIIs{i}, '_e(\d+)', 'tokens', 'once');
        if ~isempty(ec)
            names{i} = [names{i} '_echo-' ec{1}];
        end
        if endsWith(s.SeriesDescription, '_SBRef')
            if skip_sbref
                types{i} = 'other';
                continue;
            else
                names{i} = [names{i} '_sbref'];
            end
        else
            names{i} = [names{i} '_bold'];                               
        end
    elseif seqContains({'epse'}) % need to verify
        types{i} = 'fmap'; names{i} = 'fieldmap';
    elseif seqContains({'fm2d'})
        types{i} = 'fmap';
        if endsWith(NIIs{i}, '_e1')
            names{i} = 'magnitude1';
        elseif endsWith(NIIs{i}, '_e2')
            names{i} = 'magnitude2';
        elseif contains(s.ImageType, '\P\')
            names{i} = 'phasediff';
        else
            names{i} = 'fieldmap';
        end
    elseif seqContains({'tir2d'})
        types{i} = 'anat'; names{i} = 'FLAIR';
    elseif seqContains({'spc' 'tse2d'})
        types{i} = 'anat'; names{i} = 'T2w';
    elseif seqContains({'tfl3d'})
        types{i} = 'anat'; names{i} = 'T1w';
    else
        warning('Unknown SequenceName %s for %s', s.SequenceName, NIIs{i});
    end
end

for i = 1:N*(numel(unique(PEdir(types == "func"))) > 1)
    if types{i} ~= "func", continue; end
    if ~isempty(PEdir{i})
        names{i} = regexprep(names{i}, 'task-\w*?(?=_)', ['$0_' PEdir{i}]);
    else
        names{i} = regexprep(names{i}, 'task-', task);
    end
    nam = dir([niiDir NIIs{i} '.nii*']);
    hdr = nii_tool('hdr', [niiDir nam(1).name]);
    if hdr.dim(5)<=20 && ~endsWith(names{i}, '_sbref') % arbituray thre for task
        names{i} = regexprep(names{i}, '_bold$', '_epi');
        types{i} = 'fmap';
    end
    if i>1 && endsWith(names{i-1}, '_sbref') && ...
         strcmp(regexprep(names{i}, '_epi$', ''), regexprep(names{i-1}, '_sbref$', ''))
        types{i-1} = 'fmap';
    end
end

for i = 1:N*(numel(unique(PEdir(types == "dwi"))) > 1)
    if types{i} ~= "dwi", continue; end
    if ~isempty(PEdir{i})
        names{i} = [PEdir{i} names{i}];
    end
    nam = dir([niiDir NIIs{i} '.nii*']);
    hdr = nii_tool('hdr', [niiDir nam(1).name]);
    if hdr.dim(5)<6 && ~endsWith(names{i}, '_sbref') % arbituray thre for dwi
        types{i} = 'fmap';
    end
    if i>1 && endsWith(names{i-1}, '_sbref') && ...
         strcmp(regexprep(names{i}, '_dwi$', ''), regexprep(names{i-1}, '_sbref$', ''))
        types{i-1} = 'fmap';
    end
end

for i = 1:N*(numel(unique(PEdir(types == "fmap"))) > 1)
    if types{i} == "fmap" && ~isempty(PEdir{i})
        names{i} = [PEdir{i} '_' names{i}];
    end
end

for i = 1:N-1 % deal with name conflict
    if isempty(types{i}), continue; end
    ind = find(strcmp(strcat(types,names), [types{i} names{i}]));
    ind = ind(ind>i);
    if isempty(ind)
        continue; 
    end
%     if types{i} == "func"       %%%%%%%%%%% verbessern! bis dahin lassen...
%         curr = dir([niiDir NIIs{i} '.nii*']);
%         last = dir([niiDir NIIs{ind(end)} '.nii*']);
%         if curr.bytes < last.bytes-(last.bytes/8) % arbitrary               %%%%%%%%%%% changed
%             types{i} = ''; continue; % stopped/incomplete run, skip it
%         end
%     end
    names{ind} = [names{ind} '_' num2str(ind)]; % append nii name to avoid overwrite (PROBLEM: nii-names might be indentical, too!) OR append i to avoid overwrite
    warning('Identical names detected! Please check input and output data for possible mistakes.');
end

for i = 1:N
    if isempty(types{i}), continue; end
    if sess(1) == 1
        f = [fSubj 'ses-' num2str(ses) filesep types{i} filesep];
        if ~exist(f, 'dir'), mkdir(f); end
        dst = [f 'sub-' subj '_' names{i}];
    elseif sess(2) > 1
        f = [fSubj types{i} '/'];
        if ~exist(f, 'dir'), mkdir(f); end
        dst = [f 'sub-' subj '_ses-' num2str(ses) '_' names{i}];
    else
        f = [fSubj types{i} '/'];
        if ~exist(f, 'dir'), mkdir(f); end
        dst = [f 'sub-' subj '_' names{i}];
    end
    nam = dir([niiDir NIIs{i} '.*']);
    for j = 1:numel(nam)
        ext = strrep(nam(j).name, NIIs{i}, '');
        movefile([niiDir nam(j).name], [dst ext]);
    end
end

nams = dir(niiDir);
nams([nams.isdir]) = [];
for i = 1:numel(nams)
    movefile([niiDir nams(i).name], [fSubj '.']);
end
