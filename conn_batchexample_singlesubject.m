function conn_batchexample_singlesubject
% This example prompts the user to select one subject functional series
% (one resting state session) and one anatomical volume, and performs all 
% standard SPM preprocessing steps (realignment,coregistration,
% normalization,smoothing, etc.) and connectivity analyses using all default atlas ROIs
%

%% batch preprocessing for single-subject single-session data 

% Selects functional / anatomical volumes
FUNCTIONAL_FILE=cellstr(spm_select([1,Inf],'\.img$|\.nii$','Select functional volumes'));
if isempty(FUNCTIONAL_FILE),return;end
STRUCTURAL_FILE=cellstr(spm_select(1,'\.img$|\.nii$','Select structural volume'));
if isempty(STRUCTURAL_FILE),return;end
% Selects TR (seconds)
TR=conn_menu_inputdlg({'Enter Repetition-Time (in seconds)'},'TR',1,{'2'});
TR=str2num(TR{1});

%% CONN New experiment
batch.filename=fullfile(pwd,'conn_singlesubject01.mat');
if ~isempty(dir(batch.filename)), 
    Ransw=questdlg('conn_singlesubject01 project already exists, Overwrite?','warning','Yes','No','No');
    if ~strcmp(Ransw,'Yes'), return; end; 
end
%% CONN Setup
batch.Setup.nsubjects=1;
batch.Setup.functionals{1}{1}=FUNCTIONAL_FILE;
batch.Setup.structurals{1}=STRUCTURAL_FILE;
batch.Setup.RT=TR;
batch.Setup.rois.names={'atlas'};
batch.Setup.rois.files{1}=fullfile(fileparts(which('conn')),'rois','atlas.nii');
batch.Setup.conditions.names={'rest'};       
batch.Setup.conditions.onsets{1}{1}{1}=0;
batch.Setup.conditions.durations{1}{1}{1}=inf; 
batch.Setup.preprocessing.steps={}; % user-ask 
batch.Setup.isnew=1;
batch.Setup.done=1;
%% CONN Denoising
batch.Denoising.filter=[0.01, 0.1];          % frequency filter (band-pass values, in Hz)
batch.Denoising.done=1;
%% CONN Analysis
batch.Analysis.analysis_number=1;       % Sequential number identifying each set of independent first-level analyses
batch.Analysis.measure=1;               % connectivity measure used {1 = 'correlation (bivariate)', 2 = 'correlation (semipartial)', 3 = 'regression (bivariate)', 4 = 'regression (multivariate)';
batch.Analysis.weight=2;                % within-condition weight used {1 = 'none', 2 = 'hrf', 3 = 'hanning';
batch.Analysis.sources={};              % (defaults to all ROIs)
batch.Analysis.done=1;

conn_batch(batch);

%% CONN Display
% launches conn gui to explore results
conn
conn('load',fullfile(pwd,'conn_singlesubject01.mat'));
conn gui_results


end