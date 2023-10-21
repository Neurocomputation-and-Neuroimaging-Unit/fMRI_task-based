% Proviedes input for dicm2bids

%% input
%% you'll be needing the dicom2nifti toolbox ( https://de.mathworks.com/matlabcentral/fileexchange/42997-xiangruili-dicm2nii )

cd('C:\Users\saraw\Desktop\BIDS')

sourceF = 'C:\Users\saraw\Desktop\BIDS\raw'; %this is the folder that contains your raw data (folders eg. ccnb_12345)
targetF = 'C:\Users\saraw\Desktop\BIDS\test';

addpath('C:\...\Toolboxes\dicm2nii-master')

sess = [0 1]; % do wou want an extra instance for session folders? (yes = 1, no = 0) 
    % How many sessions are there per subject? (eg. 1 or 4) NOT RUNS!
    % so if i had 4 sessions per subject and want an extra folder per session: [1 4]
    % if i had a single session and dont want an extra folder: [0 1]
    % if i had 370ish sessions and wanted no extra folder per session: [0 370] (oficially a bad decision though)
    % this script can NOT (yet?) handle different amounts of sessions for differnt subjects! 

if ~exist(targetF, 'dir')
    mkdir(targetF);
    
end

cd(sourceF);
% pb = dir('*loe_*');
% pb = [pb dir('*pha_*')];
pb = [pb dir('*ccnb_*')]; % this one to find our standard-style exported data
% pb = dir('**/*ccnb_*'); % this one if you'd also like to search in sub-folders

task = 'task'; % the only thing you get to able to personalize about the file-names

%% convert (no need to change anything here)

if sess(1) == 1 && sess(2) == length(pb) %% 1 sj and a lot of seesions to store in diff folders
    fSubj = [targetF '/sub-001/'];
    if ~exist(fSubj, 'dir'), mkdir(fSubj); end
    for ses = 1:sess(2)
        fSess = [fSubj 'ses-' num2str(ses) '/'];
        if ~exist(fSess, 'dir'), mkdir(fSess); end
        dcmDir = [pb(ses).folder filesep pb(ses).name];
        A1_dicm2bids(dcmDir, targetF, 1, task, sess, ses);    %%%% übergib session: hier halten wir uns an session-ordner & benennung! irgendwie...
    end

elseif sess(1) == 1 % more sessions and subjects + extra instance for session folder
    if mod(length(pb), sess(2)) ~= 0
        error('Unexpected number of input or sessions!')
    end
    numSJ = length(pb)/sess(2);
    for subj = 1:numSJ
        if subj < 10
            sub = ['00' num2str(subj)];
        elseif subj < 100
            sub = ['0' num2str(subj)];
        else
            sub = num2str(subj);
        end
        fSubj = [targetF '/sub-' sub '/'];
        if ~exist(fSubj, 'dir'), mkdir(fSubj); end
        for ses = 1:sess(2)
            fSess = [fSubj 'ses-' num2str(ses) '/'];
            if ~exist(fSess, 'dir'), mkdir(fSess); end
            dcmDir = [pb(ses+(subj-1)*(sess(2))).folder filesep pb(ses+(subj-1)*(sess(2))).name];
            A1_dicm2bids(dcmDir, targetF, 1, task, sess, ses);    %%%% übergib session: hier halten wir uns an session-ordner & benennung! irgendwie...
        end
    end

elseif sess(2) > 1 % no session-folder-instance
    if sess(2) == length(pb) % 1 subject
        fSubj = [targetF '/sub-001/'];
        if ~exist(fSubj, 'dir'), mkdir(fSubj); end
        for ses = 1:sess(2)
            dcmDir = [pb(ses).folder filesep pb(ses).name];
            A1_dicm2bids(dcmDir, targetF, 1, task, sess, ses);    %%%% übergib session: hier halten wir uns an session-ordner & benennung! irgendwie...
        end
    else % more that 1 session and subject
        if mod(length(pb), sess(2)) ~= 0
        error('Unexpected number of input or sessions!')
        end
        numSJ = length(pb)/sess(2);
        for subj = 1:numSJ
            if subj < 10
                sub = ['00' num2str(subj)];
            elseif subj < 100
                sub = ['0' num2str(subj)];
            else
                sub = num2str(subj);
            end
            fSubj = [targetF '/sub-' sub '/'];
            if ~exist(fSubj, 'dir'), mkdir(fSubj); end
            for ses = 1:sess(2)
                dcmDir = [pb(ses+(subj-1)*(sess(2))).folder filesep pb(ses+(subj-1)*(sess(2))).name];
                A1_dicm2bids(dcmDir, targetF, subj, task, sess, ses);    %%%% übergib session: hier halten wir uns an session-ordner & benennung! irgendwie...
            end
        end
    end
else % no sessions at all, only subjects!
    for subj = 1:length(pb) %% subj must be numeric
        if subj < 10
                sub = ['00' num2str(subj)];
            elseif subj < 100
                sub = ['0' num2str(subj)];
            else
                sub = num2str(subj);
        end
        fSubj = [targetF '/sub-' sub '/'];
        if ~exist(fSubj, 'dir'), mkdir(fSubj); end
        dcmDir = [pb(subj).folder filesep pb(subj).name];
        A1_dicm2bids(dcmDir, targetF, subj, task, sess);
    end
end




