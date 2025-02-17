
clear batch;
batch.filename='/data/conn_batchexample_simple/conn_test.mat';          % new project
if conn_existfile(batch.filename)&&~isequal(questdlg('Project already exists, Overwrite?','warning','Yes','No','No'),'Yes'), return; end

%% Define CONN Setup info
batch.parallel.N = 2;
batch.Setup.isnew = 1;                                                    % note: this will overwrite/reset any prior info on this project (WARNING!!! remove this line if you are editing this script to continue processing/analyzing an existing CONN project)
batch.Setup.nsubjects=2;
batch.Setup.RT=2;                                                       % TR in seconds
filepath = '/data/conn_batchexample_simple/data/'; 
batch.Setup.functionals{1}{1} = [filepath 'sub-0001/func/sub-0001_run-01_bold.nii'];          % note: list functionals{subject}{session} (in this example 2 subjects, each with 3 functional runs)
batch.Setup.functionals{1}{2} = [filepath 'sub-0001/func/sub-0001_run-02_bold.nii'];
batch.Setup.functionals{1}{3} = [filepath 'sub-0001/func/sub-0001_run-03_bold.nii'];
batch.Setup.structurals{1}    = [filepath 'sub-0001/anat/sub-0001_run-01_T1w.nii'];
batch.Setup.functionals{2}{1} = [filepath 'sub-0002/func/sub-0002_run-01_bold.nii'];
batch.Setup.functionals{2}{2} = [filepath 'sub-0002/func/sub-0002_run-02_bold.nii'];
batch.Setup.functionals{2}{3} = [filepath 'sub-0002/func/sub-0002_run-03_bold.nii'];
batch.Setup.structurals{2}    = [filepath 'sub-0002/anat/sub-0002_run-01_T1w.nii'];

batch.Setup.conditions.names = {'rest'};
batch.Setup.conditions.onsets{1}{1} = {0, 0, 0};                          % note: onsets{condition}{subject}{session} (use brackets for multiple values, e.g. onsets={ [0,30,60], [15,45], [0,30,60] };
batch.Setup.conditions.onsets{1}{2} = {0, 0, 0}; 
batch.Setup.conditions.durations{1}{1} = {inf, inf, inf};
batch.Setup.conditions.durations{1}{2} = {inf, inf, inf};

batch.Setup.preprocessing.steps = {'default_mni'};                        % name of preprocessing pipeline (or list of individual preprocessing steps)
batch.Setup.preprocessing.sliceorder = 'interleaved (Siemens)';           % slice acquisition order
batch.Setup.preprocessing.art_thresholds = [3, 0.5];                      % user "conservative" outlier identification thresholds
batch.Setup.done = 1;                                                     % note: this will run the 'Setup' step

%% Define CONN Denoising info
batch.Denoising.filter = [0.01, 0.1];                                     % frequency filter (band-pass values, in Hz)
batch.Denoising.done = 1;                                                 % note: this will run the 'Denoising' step

%% Define CONN Analysis info
batch.Analysis.analysis_name = 'SBC_01';                                  % 1st-level SBC/RRC analysis name
%batch.Analysis.sources = {};                                              % (defaults to all ROIs)
batch.Analysis.sources = {'networks'};                                   % select all ROIs from the "networks" ROI
%batch.Analysis.sources = {'networks.DefaultMode.MPFC','networks.DefaultMode.PCC'}; % select explicit list of ROIs
batch.Analysis.done = 1;                                                  % note: this will run the '1st-level Analysis' step for the analysis defined above

%% Run batch commands
conn_batch(batch);
