% CONN_BATCH_HUMANCONNECTOMEPROJECT batch processing script for the Human Connectome Project resting-state dataset (HCP; http://www.humanconnectome.org/; tested on Q6 497-subjects release)
% 
% This script assumes that the dataset is already available in your system (set the CONNECTOMEpath variable in this script to the appropriate folder)
%
% The script will create a new conn_HCP project and:
%   Load the preprocessed "clean" ICA+FIX functional series (four sessions)
%   Apply structural segmentation, functional ART scrubbing, and functional smoothing
%   Apply default aCompCor denoising and band-pass filtering
%   Compute seed-to-voxel and ROI-to-ROI bivariate correlation connectivity measures for all CONN default ROIs (atlas + dmn)
%
% Default settings (see code for details):
%    Run each individual subject in a separate parallel stream (edit RUNPARALLEL and NJOBS variables in this script to modify these settings)
%    Assumes user has write permission into connectome data folders (edit COPYFILES variable in this script when users have only read-permissions)
%    Run all subjects in dataset (edit NSUBJECTS in this script to process instead only a subset of subjects)
%


% note: before running this script it is recommended to test it first with just a few subjects to make sure everything is working
% as expected. To do so, set NSUBJECTS=4; in the DEFAULT SETTINGS section of this script
%

%% DEFAULT SETTINGS: EDIT THE LINES BELOW (minimally set CONNECTOMEpath to the actual location of your connectome data)
TARGETpath=pwd;                                                             % target folder for conn project (default current folder)
CONNECTOMEpath='/projectnb/connectomedb/data/%s/MNINonLinear';              % source of connectome dataset (%s stands for numeric subject-specific folders)
RUNPARALLEL=true;                                                           % run in parallel using computer cluster
NSUBJECTS=[];                                                               % number of subjects to include in your project (leave empty for all subjects)
NJOBS=[];                                                                   % number of parallel jobs to submit (leave empty for one job per subject)
COPYFILES=false;                                                            % true/false: set to true if you do not have write-permissions into connectome data folders
                                                                            %   This will create a local copy (in same folder as your conn project) of the structural/functional data
                                                                            %   where any post-processed files will also be stored.
OVERWRITE=false;                                                            % overwrites files if they exist in target folder (unzipped files and/or files in local-copy folder)
                                                                            %   Set to false if you have already unzipped / copied-to-local-folder your data and would like to skip this step


%% FINDS STRUCTURAL/FUNCTIONAL/REALIGNMENT FILES
clear FUNCTIONAL_FILE* REALIGNMENT_FILE* STRUCTURAL_FILE;
subs=dir(regexprep(CONNECTOMEpath,'%s.*$','*')); 
subs=subs([subs.isdir]>0);
subs={subs.name};
subs=subs(cellfun(@(s)all(s>='0'&s<='9'),subs));
if isempty(NSUBJECTS), NSUBJECTS=numel(subs); 
else subs=subs(1:NSUBJECTS);
end
if isempty(NJOBS), NJOBS=NSUBJECTS; end
NJOBS=min(NSUBJECTS,NJOBS);

for n=1:numel(subs)
    fprintf('Locating subject %s files\n',subs{n});
    
    t1=fullfile(sprintf(CONNECTOMEpath,subs{n}),'T1w_restore_brain.nii.gz');                                        % STRUCTURAL VOLUME
    f1=fullfile(sprintf(CONNECTOMEpath,subs{n}),'Results','rfMRI_REST1_LR','rfMRI_REST1_LR_hp2000_clean.nii.gz');   % FUNCTIONAL VOLUME (1/4)
    f2=fullfile(sprintf(CONNECTOMEpath,subs{n}),'Results','rfMRI_REST1_RL','rfMRI_REST1_RL_hp2000_clean.nii.gz');   % FUNCTIONAL VOLUME (2/4)
    f3=fullfile(sprintf(CONNECTOMEpath,subs{n}),'Results','rfMRI_REST2_LR','rfMRI_REST2_LR_hp2000_clean.nii.gz');   % FUNCTIONAL VOLUME (3/4)
    f4=fullfile(sprintf(CONNECTOMEpath,subs{n}),'Results','rfMRI_REST2_RL','rfMRI_REST2_RL_hp2000_clean.nii.gz');   % FUNCTIONAL VOLUME (4/4)
    r1=fullfile(sprintf(CONNECTOMEpath,subs{n}),'Results','rfMRI_REST1_LR','Movement_Regressors.txt');              % REALIGNMENT FILE (1/4)
    r2=fullfile(sprintf(CONNECTOMEpath,subs{n}),'Results','rfMRI_REST1_RL','Movement_Regressors.txt');              % REALIGNMENT FILE (2/4)
    r3=fullfile(sprintf(CONNECTOMEpath,subs{n}),'Results','rfMRI_REST2_LR','Movement_Regressors.txt');              % REALIGNMENT FILE (3/4)
    r4=fullfile(sprintf(CONNECTOMEpath,subs{n}),'Results','rfMRI_REST2_RL','Movement_Regressors.txt');              % REALIGNMENT FILE (4/4)
    if isempty(dir(t1)), error('file %s not found',t1); end
    if isempty(dir(f1)), error('file %s not found',f1); end
    if isempty(dir(f2)), error('file %s not found',f2); end
    if isempty(dir(r1)), warning('file %s not found',r1); end
    if isempty(dir(r2)), warning('file %s not found',r2); end
    
    if COPYFILES
        fprintf('Copying files to local folder\n');
        [ok,nill]=mkdir(TARGETpath,'LocalCopyDataFiles');
        [ok,nill]=mkdir(fullfile(TARGETpath,'LocalCopyDataFiles'),subs{n});
        t1b=fullfile(TARGETpath,'LocalCopyDataFiles',subs{n},'structural.nii.gz');  if OVERWRITE||isempty(dir(t1b)), [ok,nill]=system(sprintf('cp ''%s'' ''%s''',t1,t1b)); end; t1=t1b;
        f1b=fullfile(TARGETpath,'LocalCopyDataFiles',subs{n},'functional1.nii.gz'); if OVERWRITE||isempty(dir(f1b)), [ok,nill]=system(sprintf('cp ''%s'' ''%s''',f1,f1b)); end; f1=f1b;
        f2b=fullfile(TARGETpath,'LocalCopyDataFiles',subs{n},'functional2.nii.gz'); if OVERWRITE||isempty(dir(f2b)), [ok,nill]=system(sprintf('cp ''%s'' ''%s''',f2,f2b)); end; f2=f2b;
        f3b=fullfile(TARGETpath,'LocalCopyDataFiles',subs{n},'functional3.nii.gz'); if OVERWRITE||isempty(dir(f3b)), [ok,nill]=system(sprintf('cp ''%s'' ''%s''',f3,f3b)); end; f3=f3b;
        f4b=fullfile(TARGETpath,'LocalCopyDataFiles',subs{n},'functional4.nii.gz'); if OVERWRITE||isempty(dir(f4b)), [ok,nill]=system(sprintf('cp ''%s'' ''%s''',f4,f4b)); end; f4=f4b;
        r1b=fullfile(TARGETpath,'LocalCopyDataFiles',subs{n},'Movement1.txt'); if OVERWRITE||isempty(dir(r1b)), [ok,nill]=system(sprintf('cp ''%s'' ''%s''',r1,r1b)); end; r1=r1b;
        r2b=fullfile(TARGETpath,'LocalCopyDataFiles',subs{n},'Movement2.txt'); if OVERWRITE||isempty(dir(r2b)), [ok,nill]=system(sprintf('cp ''%s'' ''%s''',r2,r2b)); end; r2=r2b;
        r3b=fullfile(TARGETpath,'LocalCopyDataFiles',subs{n},'Movement3.txt'); if OVERWRITE||isempty(dir(r3b)), [ok,nill]=system(sprintf('cp ''%s'' ''%s''',r3,r3b)); end; r3=r3b;
        r4b=fullfile(TARGETpath,'LocalCopyDataFiles',subs{n},'Movement4.txt'); if OVERWRITE||isempty(dir(r4b)), [ok,nill]=system(sprintf('cp ''%s'' ''%s''',r4,r4b)); end; r4=r4b;
    end
    fprintf('Unzipping files\n');
    if OVERWRITE||isempty(dir(regexprep(t1,'\.gz$',''))), gunzip(t1); end; t1=regexprep(t1,'\.gz$','');
    if OVERWRITE||isempty(dir(regexprep(f1,'\.gz$',''))), gunzip(f1); end; f1=regexprep(f1,'\.gz$','');
    if OVERWRITE||isempty(dir(regexprep(f2,'\.gz$',''))), gunzip(f2); end; f2=regexprep(f2,'\.gz$','');
    if OVERWRITE||isempty(dir(regexprep(f3,'\.gz$',''))), gunzip(f3); end; f3=regexprep(f3,'\.gz$','');
    if OVERWRITE||isempty(dir(regexprep(f4,'\.gz$',''))), gunzip(f4); end; f4=regexprep(f4,'\.gz$','');
    r1b=regexprep(r1,'\.txt$','\.deg.txt'); if OVERWRITE||isempty(dir(r1b)), [ok,nill]=system(sprintf('cp ''%s'' ''%s''',r1,r1b)); end; r1=r1b; % note: angles in degrees
    r2b=regexprep(r2,'\.txt$','\.deg.txt'); if OVERWRITE||isempty(dir(r2b)), [ok,nill]=system(sprintf('cp ''%s'' ''%s''',r2,r2b)); end; r2=r2b; 
    r3b=regexprep(r3,'\.txt$','\.deg.txt'); if OVERWRITE||isempty(dir(r3b)), [ok,nill]=system(sprintf('cp ''%s'' ''%s''',r3,r3b)); end; r3=r3b; 
    r4b=regexprep(r4,'\.txt$','\.deg.txt'); if OVERWRITE||isempty(dir(r4b)), [ok,nill]=system(sprintf('cp ''%s'' ''%s''',r4,r4b)); end; r4=r4b; 
    
    STRUCTURAL_FILE{n,1}=t1;
    FUNCTIONAL_FILE(n,1:4)={f1,f2,f3,f4};
    REALIGNMENT_FILE(n,1:4)={r1,r2,r3,r4};
end
nsessions=4;
fprintf('%d subjects, %d sessions\n',NSUBJECTS,nsessions);



%% CREATES CONN BATCH STRUCTURE
clear batch;
batch.filename=fullfile(TARGETpath,'conn_HCP.mat');
if RUNPARALLEL
    batch.parallel.N=NJOBS;                             % number of parallel processing batch jobs
end                                                     % note: use default parallel profile (defined in GUI Tools.GridSettings)

% CONN Setup                                           
batch.Setup.isnew=1;
batch.Setup.nsubjects=NSUBJECTS;
batch.Setup.RT=0.72;                                    % TR (seconds)

batch.Setup.conditions.names={'rest'};                  % single condition (aggregate across all sessions)
for ncond=1,for nsub=1:NSUBJECTS,for nses=1:nsessions,      batch.Setup.conditions.onsets{ncond}{nsub}{nses}=0; batch.Setup.conditions.durations{ncond}{nsub}{nses}=inf;end;end;end     % rest condition (all sessions)

batch.Setup.functionals=repmat({{}},[NSUBJECTS,1]);     % Point to functional volumes for each subject/session
for nsub=1:NSUBJECTS,for nses=1:nsessions,                  batch.Setup.functionals{nsub}{nses}=FUNCTIONAL_FILE(nsub,nses);end; end 

batch.Setup.structurals=STRUCTURAL_FILE;                % Point to anatomical volumes for each subject

batch.Setup.voxelresolution=1;                          % default 2mm isotropic voxels analysis space

batch.Setup.covariates.names={'realignment'};
batch.Setup.covariates.files{1}=repmat({{}},[NSUBJECTS,1]);      
for nsub=1:NSUBJECTS,for nses=1:nsessions,                  batch.Setup.covariates.files{1}{nsub}{nses}=REALIGNMENT_FILE(nsub,nses);end; end 

batch.Setup.analyses=[1,2];                             % seed-to-voxel and ROI-to-ROI pipelines
batch.Setup.overwrite='Yes';                            
batch.Setup.done=1;

batch.Setup.preprocessing.steps={'structural_segment','functional_art','functional_smooth'};  % Run additional preprocessing steps: segmentation,ART,smoothing
batch.Setup.preprocessing.fwhm=6;                       % smoothing fwhm (mm)
conn_batch(batch);

clear batch;
% CONN Denoising                                    
batch.Denoising.filter=[0.01, 0.10];                    % frequency filter (band-pass values, in Hz)
batch.Denoising.detrending=2;
batch.Denoising.done=1;                                 % use default denoising step (CompCor, motion regression, scrubbing, detrending)
batch.Denoising.overwrite='Yes';
conn_batch(batch);

clear batch;
% CONN Analysis                                         % Default options (uses all ROIs in conn/rois/ as connectivity sources); see conn_batch for additional options 
batch.Analysis.done=1;
batch.Analysis.overwrite='Yes';


%% RUNS CONN BATCH STRUCTURE
conn_batch(batch);

