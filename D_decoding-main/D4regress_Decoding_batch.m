function D4regress_Decoding_batch(beta_dir, output_dir, labelnames, dec_type, rad, masks)
% This Batch Script first specifies what features should be decoded and
% executes the Decoding analysis.
% The resulting accuracy maps should later be normalized and smoothed

% CREATE THE OUTPUT DIRECTORY IF IT DOES NOT EXIST YET
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% CONFIGURE THE DECODING TOOLBOX

clear cfg
cfg = decoding_defaults;
cfg.software = 'SPM12';

% Specify where the results should be saved
cfg.results.overwrite     = 1;
cfg.results.dir           = output_dir;

% DATA SCALING
cfg.scale.method          = 'z'; %z-score for all voxels over all samples (2 conditions x 4 runs =8)
cfg.scale.estimation      = 'all'; % only training, only test data, or all

% SEARCHLIGHT SPECIFICATIONS
cfg.analysis              = dec_type;   % or 'roi'
cfg.searchlight.unit      = 'voxels'; % or mm 
cfg.searchlight.radius    = rad; % 4 voxel or mm
cfg.searchlight.spherical = 0;  % only useful for mm
% The amount of information you want to have printed on the screen
% 0: no, 1: normal, output, 2: all)
cfg.verbose               = 0;  

% Method and model parameters 
%cfg.decoding.method = 'classification';
%cfg.decoding.train.classification.model_parameters = '-s 0 -t 0 -c 1 -b 0 -q'; 
cfg.decoding.method = 'regression';

% OUTPUTS SPECIFICATION
cfg.results.output = {'corr', 'zcorr'};

% DISPLAY:
cfg.plot_selected_voxels  = 0; % 0: no plotting, 1: every step, 2: every second step, 100: every hundredth step...

% This is by default set to 1, but if you repeat the same design again and again, it can get annoying...
cfg.plot_design           = 0;

% Set the filename of your brain mask (or your ROI masks as cell matrix) 
% for searchlight or wholebrain e.g. 'c:\exp\glm\model_button\mask.img' OR 
% for ROI e.g. {'c:\exp\roi\roimaskleft.img', 'c:\exp\roi\roimaskright.img'}
% You can also use a mask file with multiple masks inside that are
% separated by different integer values (a "multi-mask")

% cfg.files.mask = fullfile(beta_dir, 'mask.nii');
if dec_type == 'roi'
    cfg.files.mask = masks;
else
    cfg.files.mask = [beta_dir filesep 'mask.nii'];
end
% cfg.files.mask = 'C:\Users\saraw\Desktop\Masks\MASKS2\thatmask\rCONJ_att_CUT_allOfSomato.nii';

% cfg.basic_checks.DoubleFilenameEntriesOk = 1; %%%% ?

% Decoding DESIGN
% labels = [1 2 3 4];
% labels = [-1.5*ones(4,1); 
%     -0.5*ones(4,1);
%     0.5*ones(4,1);
%     1.5*ones(4,1)];
labels = [-1.5; 
    -0.5;
    0.5;
    1.5];

% cfg.files.chunk = [1:4, 1:4, 1:4, 1:4]';
% cfg.files.chunk = [1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4]';

% The following function extracts all beta names and corresponding run
% numbers from the SPM.mat
regressor_names = design_from_spm(beta_dir);

% Extract all information for the cfg.files structure
cfg = decoding_describe_data(cfg,labelnames,labels,regressor_names,beta_dir); %cfg = decoding_describe_data(cfg,{labelname1,labelname2,labelname3,labelname4},[1 2 3 4],regressor_names,beta_dir); 


% This creates the leave-one-run-out cross validation design:
cfg.design = make_design_cv(cfg);
display_design(cfg);

% Run decoding
results = decoding(cfg);
