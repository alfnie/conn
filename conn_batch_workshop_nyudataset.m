function conn_batch_workshop_nyudataset(varargin)
%
% batch processing script for the NYU_CSC_TestRetest dataset (published in Shehzad et al., 2009, The Resting Brain: Unconstrained yet Reliable. Cerebral Cortex. doi:10.1093/cercor/bhn256)
% 
% Steps:
% 1. Run conn_batch_workshop_nyudataset. The script will:
%       a) Download the entire dataset from: https://www.nitrc.org/projects/nyu_trt/ (6 NYU_TRT_session*.tar.gz files)
%       b) Decompress the NYU_TRT_session*.tar.gz files into NYU_TRT_session* folders
%          Decompress the *.nii.gz files
%       c) Preprocessing of the anatomical and functional volumes
%         (normalization & segmentation of anatomical volumes; realignment,
%         coregistration, normalization, outlier detection, and smooting of the 
%         functional volumes)
%       d) Estimate first-level seed-to-voxel connectivity maps for each of 
%         the default seeds (located in the conn/rois folder), separately 
%         for each subject and for each of the three test-retest sessions.
%
% Optionally, to download only a portion of the data, specify the desired data subsets (1a,1b,2a,2b,3a,3b) as additional arguments. For example:
%   conn_batch_workshop_nyudataset_parallel 1a;           % downloads and processes only 13-subjects one-session data 
%   conn_batch_workshop_nyudataset_parallel 1a 1b;        % downloads and processes only 25-subjects one-session data 
%   conn_batch_workshop_nyudataset_parallel 1a 2a 3a;     % downloads and processes only 13-subjects three-session data 
%       
%


OVERWRITE=false;    % skip downloading/uncompressing already-existing files/folders
DOWNLOADFILES=true; % set to false if already manually-downloaded the dataset (or a portion of it)
UNZIPFILES=true;    % set to false if already manually-unzipped the dataset
DOPARALLEL=false;   % set to true/false to run in parallel or locally
PARALLELPROFILE=[];
    
iopt=strcmp(varargin,'-overwrite'); if any(iopt), varargin=varargin(~iopt); OVERWRITE=true; end
iopt=strcmp(varargin,'-download'); if any(iopt), varargin=varargin(~iopt); DOWNLOADFILES=true; end
iopt=strcmp(varargin,'-unzip'); if any(iopt), varargin=varargin(~iopt); UNZIPFILES=true; end
iopt=strcmp(varargin,'-donotoverwrite'); if any(iopt), varargin=varargin(~iopt); OVERWRITE=false; end
iopt=strcmp(varargin,'-donotdownload'); if any(iopt), varargin=varargin(~iopt); DOWNLOADFILES=false; end
iopt=strcmp(varargin,'-donotunzip'); if any(iopt), varargin=varargin(~iopt); UNZIPFILES=false; end
iopt=strcmp(varargin,'-parallel'); if any(iopt), varargin=varargin(~iopt); DOPARALLEL=true; end

if ~nargin||isempty(varargin), 
    data={'1a','1b','2a','2b','3a','3b'}; % define the subsets to be downloaded (all data by default) 
    try
        answ=conn_questdlg({'Preparing to download and analyze the publicly available NYU test-retest dataset (https://www.nitrc.org/projects/nyu_trt/)','This procedure will create a new CONN project with the fully analyzed NYU test-retest dataset','The full dataset requires approximately 200Gb of hard-drive space and up to 30 hours to finish. You may download/process the full dataset or only a portion of it',' ','Please select dataset to download/use:'},'Sample data download','Full dataset','Partial dataset','Already downloaded dataset','Cancel','Full dataset');
        if isempty(answ)||isequal(answ,'Cancel'), return; end
        if isequal(answ,'Partial dataset')
            datastr={'Subset A (13 subjects), session 1','Subset B (12 subjects), session 1','Subset A (13 subjects), session 2','Subset B (12 subjects), session 2','Subset A (13 subjects), session 3','Subset B (12 subjects), session 3'};
            s = listdlg('PromptString','Select data subset(s) to download','ListSize',[300 100],...
                      'SelectionMode','multiple',...
                      'ListString',datastr,...
                      'InitialValue',1:6);
            if isempty(s), return; end
            data=data(s);
        end
        if isequal(answ,'Already downloaded dataset')
            str='Select folder containing NYU_TRT_session* files/folders';
            disp(str);
            cwd=uigetdir(pwd,str);
            if isequal(cwd,0), return; end
            folders=dir(fullfile(cwd,'NYU_TRT_session*'));
            if isempty(folders), return; end
            folders={folders.name};
            folders=regexp(folders,'^NYU_TRT_session([123][ab]).*','tokens','once');
            data=unique([folders{:}]);
            DOWNLOADFILES=false;
        else
            cwd=pwd;
            answ=conn_questdlg(sprintf('Files will be downloaded to %s',cwd),'','Continue','Modify','Cancel','Continue');
            if isempty(answ)||strcmp(answ,'Cancel'), return; end
            if strcmp(answ,'Modify'), cwd=uigetdir(pwd,'Select target folder to store the new CONN project and data'); end
        end
        if isequal(cwd,0), return; end
        cd(cwd);
        pr0='Run locally on this computer';
        pr1=conn_jobmanager('getprofile');
        pr={pr0, pr1};
        %pr=conn_jobmanager('profiles');
        %pr=[{pr0} {pr1} pr(~ismember(pr,pr1))];
        answ=conn_questdlg('','Parallelization',pr{:},pr0);
        if isempty(answ), return; end
        if isequal(answ,pr0), DOPARALLEL=false; 
        else DOPARALLEL=true; PARALLELPROFILE=answ;
        end
    end
else data=varargin;
end




%% DOWNLOADS/LOCATES FUNCTIONAL&ANATOMICAL DATA
%% DOWNLOAD *.tar.gz files
if DOWNLOADFILES
    for n=1:numel(data), 
        filename=['NYU_TRT_session',data{n},'.tar.gz'];
        if OVERWRITE||~exist(filename,'file')
            fprintf('Downloading %s (file %d/%d). This process may take several minutes. Please wait...\n',filename,n,numel(data));
            [fname,ok]=urlwrite(sprintf('https://www.nitrc.org/frs/download.php/%d/%s',1070+find(ismember({'1a','1b','2a','2b','3a','3b'},data{n})),filename),filename);
        end
    end
end
NSUBJECTS=0;
if any(ismember(data,{'1a','2a','3a'})), NSUBJECTS=NSUBJECTS+13; end
if any(ismember(data,{'1b','2b','3b'})), NSUBJECTS=NSUBJECTS+12; end

%% UNTAR .tar.gz files
if UNZIPFILES
    for n=1:numel(data), 
        filename=['NYU_TRT_session',data{n},'.tar.gz'];
        a=dir(filename);
        for n1=1:length(a),
            [a_path,a_name,a_ext]=fileparts(a(n1).name);[nill,a_name2,a_ext2]=fileparts(a_name);
            dirname=fullfile(a_path,a_name2);
            if ~isdir(dirname),
                disp(['extracting contents from file ',a(n1).name]);
                untar(a(n1).name,dirname);
            end
        end
        %% UNZIP .nii.gz files
        filename=['NYU_TRT_session',data{n}];
        a=strvcat(conn_dir(fullfile(pwd,filename,'lfo.nii.gz')),conn_dir(fullfile(pwd,filename,'mprage_anonymized.nii.gz')));
        for n1=1:size(a,1),
            [a_path,a_name,a_ext]=fileparts(a(n1,:));
            if isempty(dir(fullfile(a_path,a_name))),
                disp(['unzipping file ',a(n1,:)]);
                gunzip(deblank(a(n1,:)));
            end
        end
    end
end

%% FIND functional/structural files
% note: this will look for all data in these folders, irrespective of the specific download subsets entered as command-line arguments
cwd=pwd;
FUNCTIONAL_FILE={};
STRUCTURAL_FILE={};
for n=1:numel(data),
    filename=['NYU_TRT_session',data{n}];
    tFUNCTIONAL_FILE=cellstr(conn_dir(fullfile(pwd,filename,'lfo.nii')));
    tSTRUCTURAL_FILE=cellstr(conn_dir(fullfile(pwd,filename,'mprage_anonymized.nii')));
    FUNCTIONAL_FILE=[FUNCTIONAL_FILE;tFUNCTIONAL_FILE(:)];
    STRUCTURAL_FILE=[STRUCTURAL_FILE;tSTRUCTURAL_FILE(:)];
end
if ~NSUBJECTS, NSUBJECTS=length(STRUCTURAL_FILE); end
if rem(length(FUNCTIONAL_FILE),NSUBJECTS),error('mismatch number of functional files %d', length(FUNCTIONAL_FILE));end
if rem(length(STRUCTURAL_FILE),NSUBJECTS),error('mismatch number of anatomical files %d', length(FUNCTIONAL_FILE));end
nsessions=length(FUNCTIONAL_FILE)/NSUBJECTS;
FUNCTIONAL_FILE=reshape(FUNCTIONAL_FILE,[NSUBJECTS,nsessions]);
STRUCTURAL_FILE={STRUCTURAL_FILE{1:NSUBJECTS}};
disp([num2str(size(FUNCTIONAL_FILE,1)),' subjects']);
disp([num2str(size(FUNCTIONAL_FILE,2)),' sessions']);
TR=2; % Repetition time = 2 seconds


%% CONN-SPECIFIC SECTION: RUNS PREPROCESSING/SETUP/DENOISING/ANALYSIS STEPS
%% Prepares batch structure
clear batch;
batch.filename=fullfile(cwd,'conn_NYU.mat');            % New conn_*.mat experiment name
if DOPARALLEL, 
    batch.parallel.N=NSUBJECTS; 
    if ~isempty(PARALLELPROFILE), batch.parallel.profile=PARALLELPROFILE; end
end

%% SETUP & PREPROCESSING step (using default values for most parameters, see help conn_batch to define non-default values)
% CONN Setup                                            % Default options (uses all ROIs in conn/rois/ directory); see conn_batch for additional options 
% CONN Setup.preprocessing                               (realignment/coregistration/segmentation/normalization/smoothing)
batch.Setup.isnew=1;
batch.Setup.nsubjects=NSUBJECTS;
batch.Setup.RT=TR;                                        % TR (seconds)
batch.Setup.functionals=repmat({{}},[NSUBJECTS,1]);       % Point to functional volumes for each subject/session
for nsub=1:NSUBJECTS,
    for nses=1:nsessions,
        batch.Setup.functionals{nsub}{nses}=FUNCTIONAL_FILE{nsub,nses}; 
    end
end %note: each subject's data is defined by three sessions and one single (4d) file per session
batch.Setup.structurals=STRUCTURAL_FILE;                  % Point to anatomical volumes for each subject
nconditions=nsessions;                                  % treats each session as a different condition (comment the following three lines and lines 84-86 below if you do not wish to analyze between-session differences)
if nconditions==1
    batch.Setup.conditions.names={'rest'};
    for ncond=1,for nsub=1:NSUBJECTS,for nses=1:nsessions,              batch.Setup.conditions.onsets{ncond}{nsub}{nses}=0; batch.Setup.conditions.durations{ncond}{nsub}{nses}=inf;end;end;end     % rest condition (all sessions)
else
    batch.Setup.conditions.names=[{'rest'}, arrayfun(@(n)sprintf('Session%d',n),1:nconditions,'uni',0)];
    for ncond=1,for nsub=1:NSUBJECTS,for nses=1:nsessions,              batch.Setup.conditions.onsets{ncond}{nsub}{nses}=0; batch.Setup.conditions.durations{ncond}{nsub}{nses}=inf;end;end;end     % rest condition (all sessions)
    for ncond=1:nconditions,for nsub=1:NSUBJECTS,for nses=1:nsessions,  batch.Setup.conditions.onsets{1+ncond}{nsub}{nses}=[];batch.Setup.conditions.durations{1+ncond}{nsub}{nses}=[]; end;end;end
    for ncond=1:nconditions,for nsub=1:NSUBJECTS,for nses=ncond,        batch.Setup.conditions.onsets{1+ncond}{nsub}{nses}=0; batch.Setup.conditions.durations{1+ncond}{nsub}{nses}=inf;end;end;end % session-specific conditions
end
batch.Setup.preprocessing.steps='default_mni';
batch.Setup.preprocessing.sliceorder='interleaved (Siemens)';
batch.Setup.done=1;
batch.Setup.overwrite='Yes';                            

% uncomment the following 3 lines if you prefer to run one step at a time:
% conn_batch(batch); % runs Preprocessing and Setup steps only
% clear batch;
% batch.filename=fullfile(cwd,'conn_NYU.mat');            % Existing conn_*.mat experiment name

%% DENOISING step
% CONN Denoising                                    % Default options (uses White Matter+CSF+realignment+scrubbing+conditions as confound regressors); see conn_batch for additional options 
batch.Denoising.filter=[0.01, 0.1];                 % frequency filter (band-pass values, in Hz)
batch.Denoising.done=1;
batch.Denoising.overwrite='Yes';

% uncomment the following 3 lines if you prefer to run one step at a time:
% conn_batch(batch); % runs Denoising step only
% clear batch;
% batch.filename=fullfile(cwd,'conn_NYU.mat');            % Existing conn_*.mat experiment name

%% FIRST-LEVEL ANALYSIS step
% CONN Analysis                                     % Default options (uses all ROIs in conn/rois/ as connectivity sources); see conn_batch for additional options 
batch.Analysis.done=1;
batch.Analysis.overwrite='Yes';

%% Run all analyses
conn_batch(batch);

%% CONN Display
% launches conn gui to explore results
conn
conn('load',fullfile(cwd,'conn_NYU.mat'));
conn gui_results

