% Proviedes input for dicm2bids

%% input
%% you'll be needing the dicom2nifti toolbox ( https://de.mathworks.com/matlabcentral/fileexchange/42997-xiangruili-dicm2nii )

cd('C:\Users\saraw\Desktop\BIDS')

sourceF = 'C:\Users\saraw\Desktop\SEBs\biStim\raw\'; %this is the folder that contains your raw data (folders eg. ccnb_12345)
targetF = 'C:\Users\saraw\Desktop\BIDS\test\';  % C:\Users\saraw\Desktop\SEBs\multi_stim

addpath('C:\Users\saraw\Desktop\BIDS\convert\fMRI-BIDS-main\dicm2nii-master')

sess = [0 1]; % do wou want an extra instance for session folders? (yes = 1, no = 0) 
    % How many sessions are there per subject? (eg. 1 or 4) NOT RUNS!

individ_ses_length = 0; % set to 0 if every subject has the same amount of sessions,
                        % if you have differnt amount of sessions, the put
                        % a 1, ignore the second part of the sess variable
sub_offset = 0;% 1999; %if you don't start with sub-001

skip_sbref = 1;

if ~exist(targetF, 'dir')
    mkdir(targetF);
end

cd(sourceF);
% pb = dir('*loe_*');
% pb = [pb dir('*pha_*')];
% pb = dir('*raw_*');
pb = dir('*ccnb*'); % 

task = 'test';

%% convert

if individ_ses_length == 0

    if sess(1) == 1 && sess(2) == length(pb) %% 1 sj and a lot of seesions to store in diff folders
        fSubj = [targetF '/sub-001/'];
        if ~exist(fSubj, 'dir'), mkdir(fSubj); end
        for ses = 1:sess(2)
            fSess = [fSubj 'ses-' num2str(ses) '/'];
            if ~exist(fSess, 'dir'), mkdir(fSess); end
            dcmDir = [sourceF filesep pb(ses).name];
            A1_dicm2bids(dcmDir, targetF, 1, task, sess, ses);    %%%% übergib session: hier halten wir uns an session-ordner & benennung! irgendwie...
        end
    
    elseif sess(1) == 1
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
                dcmDir = [sourceF filesep pb(ses+(subj-1)*(sess(2))).name];
                A1_dicm2bids(dcmDir, targetF, 1, task, sess, ses);    %%%% übergib session: hier halten wir uns an session-ordner & benennung! irgendwie...
            end
        end
    
    elseif sess(2) > 1
        if sess(2) == length(pb)
            fSubj = [targetF '/sub-001/'];
            if ~exist(fSubj, 'dir'), mkdir(fSubj); end
            for ses = 1:sess(2)
                dcmDir = [sourceF filesep pb(ses).name];
                A1_dicm2bids(dcmDir, targetF, 1, task, sess, ses);    %%%% übergib session: hier halten wir uns an session-ordner & benennung! irgendwie...
            end
        else
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
                    dcmDir = [sourceF filesep pb(ses+(subj-1)*(sess(2))).name];
                    A1_dicm2bids(dcmDir, targetF, subj, task, sess, ses);    %%%% übergib session: hier halten wir uns an session-ordner & benennung! irgendwie...
                end
            end
        end
    else
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
            dcmDir = [sourceF filesep pb(subj).name];
            A1_dicm2bids(dcmDir, targetF, subj, task, sess, 1, skip_sbref);
        end
    end

else

    for subj = [1:length(pb)] %% subj must be numeric

        sj = subj+sub_offset;

        cd(pb(subj));
        fd = dir('*ses-*'); % assumed that b vor r vor t

        if subj < 10
            sub = ['00' num2str(sj)];
        elseif subj < 100
            sub = ['0' num2str(sj)];
        else
            sub = num2str(sj);
        end

        fSubj = [targetF '/sub-' sub '/'];
        if ~exist(fSubj, 'dir'), mkdir(fSubj); end

        for s = 1:length(fd)

            this_name = fd(s).name
            ses_name = this_name(strfind(this_name, 'ses-')+4:end);

            fSes = [targetF filesep 'sub-' sub filesep 'sub-' sub '_ses-' ses_name filesep];
            if ~exist(fSes, 'dir'), mkdir(fSes); end

            dcmDir = [sourceF filesep pb(subj).name filesep fd(s).name];
            A1_dicm2bids(dcmDir, targetF, sub, task, ses_name, 0, 1);
        end
    end


end


