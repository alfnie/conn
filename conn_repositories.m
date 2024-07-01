function varargout=conn_repositories(varargin)
% CONN_REPOSITORIES downloads and process data from public repositories
% 
% conn_repositories
%   launches GUI to download available datasets
% conn_repositories <projectname> <dataset1> <dataset2> ...
%   creates project <projectname> including the selected datasets
% conn_repositories ... [-overwrite -donotoverwrite] 
%   overwrites dataset in target directory if it already exists (default -donotoverwrite)
% conn_repositories ... [-create -donotcreate] 
%   creates a new CONN project with the downloaded data (default -create)
% conn_repositories ... [-process -donotprocess] 
%   runs preprocessing/denoising/first-level analysis steps on downloaded data (default -process)
% 
%
varargout={};

OVERWRITE=false;    % skip downloading/uncompressing already-existing files/folders
DOWNLOADFILES=true; % set to false if already manually-downloaded the dataset (or a portion of it)
UNZIPFILES=true;    % set to false if already manually-unzipped the dataset
DOCREATE=true;      % set to false to skip creating a CONN project (download only)
DOPROCESS=true;     % set to false to skip preprocessing/denoising/first-level steps
DOPARALLEL=false;   % set to true/false to run in parallel or locally
PARALLELPROFILE=[];
DODISPLAY=true;     % set to false to skip loading/displaying conn project after finishing

iopt=strcmp(varargin,'-overwrite'); if any(iopt), varargin=varargin(~iopt); OVERWRITE=true; end
iopt=strcmp(varargin,'-download'); if any(iopt), varargin=varargin(~iopt); DOWNLOADFILES=true; end
iopt=strcmp(varargin,'-unzip'); if any(iopt), varargin=varargin(~iopt); UNZIPFILES=true; end
iopt=strcmp(varargin,'-display'); if any(iopt), varargin=varargin(~iopt); DODISPLAY=true; end
iopt=strcmp(varargin,'-create'); if any(iopt), varargin=varargin(~iopt); DOCREATE=true; end
iopt=strcmp(varargin,'-process'); if any(iopt), varargin=varargin(~iopt); DOPROCESS=true; end
iopt=strcmp(varargin,'-parallel'); if any(iopt), varargin=varargin(~iopt); DOPARALLEL=true; end
iopt=strcmp(varargin,'-donotoverwrite'); if any(iopt), varargin=varargin(~iopt); OVERWRITE=false; end
iopt=strcmp(varargin,'-donotdownload'); if any(iopt), varargin=varargin(~iopt); DOWNLOADFILES=false; end
iopt=strcmp(varargin,'-donotunzip'); if any(iopt), varargin=varargin(~iopt); UNZIPFILES=false; end
iopt=strcmp(varargin,'-donotdisplay'); if any(iopt), varargin=varargin(~iopt); DODISPLAY=false; end
iopt=strcmp(varargin,'-donotcreate'); if any(iopt), varargin=varargin(~iopt); DOCREATE=false; end
iopt=strcmp(varargin,'-donotprocess'); if any(iopt), varargin=varargin(~iopt); DOPROCESS=false; end
iopt=strcmp(varargin,'-donotparallel'); if any(iopt), varargin=varargin(~iopt); DOPARALLEL=false; end

if isempty(varargin), projectname='';  % fist parameter: project name
else projectname=varargin{1}; varargin=varargin(2:end); 
end
if isempty(varargin), % other parameters: names of selected datasets
    repofiles=spm_jsonread(fullfile(fileparts(which(mfilename)),'conn_repositories.json'));

    data={repofiles.name}; % define the subsets to be downloaded (all data by default)
    datastr=cellfun(@(a,b)sprintf('<HTML>%s (%s) <small>%s</small></HTML>',a,regexp(b,'n\s*=\s*\d+','match','once'),b),{repofiles.name},{repofiles.description},'uni',0);
    s = listdlg('PromptString','Select dataset(s) from the 1000 Functional Connectomes Project to download','ListSize',[800 400],...
        'SelectionMode','multiple',...
        'ListString',datastr,...
        'InitialValue',1);
    if isempty(s), return; end
    data=data(s);
    repofiles=repofiles(s);

    DATA=0; NSUBJECTS=0; for n=1:numel(data), NSUBJECTS=NSUBJECTS+repofiles(n).subjects; DATA=DATA+ceil(repofiles(n).duration/repofiles(n).tr)*repofiles(n).subjects; end
    if ~DOCREATE
        cwd=conn_projectmanager('pwd');
        ok=false; while ~ok
            answ=conn_questdlg({sprintf('Preparing to download %d datasets (%d participants) from the 1000 Functional Connectomes Project',numel(data),NSUBJECTS),'(Mennes, M., Biswal, B. B., Castellanos, F. X., & Milham, M. P. (2013). Making data sharing work: the FCP/INDI experience. Neuroimage, 82, 683-691)',' ',sprintf('Processing the selected datasets may require approximately %d Gbs of hard-drive space as well as %d hours to finish.',ceil(3*DATA/180),ceil(DATA/180/10)),' ','Data will be downloaded to the following folder:',cwd},'','Continue (download only)','Modify target folder','Cancel','Continue (download only)');
            if isempty(answ)||strcmp(answ,'Cancel'), return; end
            if strcmp(answ,'Modify target folder'),
                if conn_projectmanager('inserver'),
                    cwd=conn_menu_inputdlg('Select target folder to store the downloaded data','data folder',1,{cwd});
                    if isempty(cwd), return; end
                    cwd=conn_server('util_remotefile',char(cwd));
                else
                    cwd=uigetdir(pwd,'Select target folder to store the downloaded data');
                    if isequal(cwd,0), return; end
                end
            else
                ok=true;
                DOPROCESS=strcmp(answ,'Continue (download and process)');
            end
        end
        PROJECTNAME='conn_FCP.mat';
    elseif ~isempty(projectname)
        [cwd,tfilename,tfileext]=fileparts(projectname);
        PROJECTNAME=[tfilename,tfileext];
        answ=conn_questdlg({sprintf('Preparing to download and process %d datasets (%d participants) from the 1000 Functional Connectomes Project',numel(data),NSUBJECTS),'(Mennes, M., Biswal, B. B., Castellanos, F. X., & Milham, M. P. (2013). Making data sharing work: the FCP/INDI experience. Neuroimage, 82, 683-691)',' ',sprintf('Processing the selected datasets may require approximately %d Gbs of hard-drive space as well as %d hours to finish.',ceil(3*DATA/180),ceil(DATA/180/10)),' ','A new ',projectname,' project will be created'},'','Continue (download and process)','Continue (download only)','Cancel','Continue (download and process)');
        if isempty(answ)||strcmp(answ,'Cancel'), return; end
        DOPROCESS=strcmp(answ,'Continue (download and process)');
    else
        cwd=conn_projectmanager('pwd');
        ok=false; while ~ok
            answ=conn_questdlg({sprintf('Preparing to download and process %d datasets (%d participants) from the 1000 Functional Connectomes Project',numel(data),NSUBJECTS),'(Mennes, M., Biswal, B. B., Castellanos, F. X., & Milham, M. P. (2013). Making data sharing work: the FCP/INDI experience. Neuroimage, 82, 683-691)',' ',sprintf('Processing the selected datasets may require approximately %d Gbs of hard-drive space as well as %d hours to finish.',ceil(3*DATA/180),ceil(DATA/180/10)),' ','A new conn_FCP.mat project will be created in the following folder:',cwd},'','Continue (download and process)','Continue (download only)','Modify target folder','Cancel','Continue (download and process)');
            if isempty(answ)||strcmp(answ,'Cancel'), return; end
            if strcmp(answ,'Modify target folder'),
                if conn_projectmanager('inserver'),
                    cwd=conn_menu_inputdlg('Select target folder to store the new conn_FCP.mat CONN project and data','data folder',1,{cwd});
                    if isempty(cwd), return; end
                    cwd=conn_server('util_remotefile',char(cwd));
                else
                    cwd=uigetdir(pwd,'Select target folder to store the new conn_FCP.mat CONN project and data');
                    if isequal(cwd,0), return; end
                end
            else
                ok=true;
                DOPROCESS=strcmp(answ,'Continue (download and process)');
            end
        end
        PROJECTNAME='conn_FCP.mat';
    end

    conn_fileutils('cd',char(cwd));
    if DOPROCESS
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
    if conn_projectmanager('inserver'), % if connected to remote server
        opts={};
        if OVERWRITE, opts{end+1}='-overwrite'; else opts{end+1}='-donotoverwrite'; end
        if DOWNLOADFILES, opts{end+1}='-download'; else opts{end+1}='-donotdownload'; end
        if UNZIPFILES, opts{end+1}='-unzip'; else opts{end+1}='-donotunzip'; end
        if 0&&DODISPLAY, opts{end+1}='-display'; else opts{end+1}='-donotdisplay'; end
        if DOCREATE, opts{end+1}='-create'; else opts{end+1}='-donotcreate'; end
        if DOPROCESS, opts{end+1}='-process'; else opts{end+1}='-donotprocess'; end
        if DOPARALLEL, opts{end+1}='-parallel'; else opts{end+1}='-donotparallel'; end
        hmsg=conn_msgbox({'Process running remotely',' ','CONN will resume automatically when this process finishes','Please wait...'},'');
        conn_server('run','conn_process','conn','repositories',fullfile(cwd,PROJECTNAME),data{:},opts{:});
        if ~isempty(hmsg)&&ishandle(hmsg), delete(hmsg); end
        if DODISPLAY&&DOCREATE
            conn('load',fullfile(cwd,PROJECTNAME));
            conn gui_setup
        end
        return
    end
else
    data=varargin
    repofiles=spm_jsonread(fullfile(fileparts(which(mfilename)),'conn_repositories.json'));
    alldata={repofiles.name}; % define the subsets to be downloaded (all data by default)
    [ok,s]=ismember(data,alldata);
    assert(all(ok), 'unrecognized datasets %s', sprintf('%s ',data{~ok}));
    repofiles=repofiles(s);
    DODISPLAY=false;
end

%% DOWNLOADS/LOCATES FUNCTIONAL&ANATOMICAL DATA
%% DOWNLOAD *.tar files
if DOWNLOADFILES
    if DODISPLAY, hmsg=conn_msgbox('Downloading dataset(s). Please wait...',''); end
    for n=1:numel(data),
        %[ok,nill]=mkdir(repofiles(n).name);
        repofiles_link=repofiles(n).link;
        repofiles_file=repofiles(n).file;
        if ~iscell(repofiles_link), repofiles_link={repofiles_link}; end
        if ~iscell(repofiles_file), repofiles_file={repofiles_file}; end
        assert(numel(repofiles_link)==numel(repofiles_file));
        for nfile=1:numel(repofiles_link)
            filename=fullfile(cwd,repofiles_file{nfile});
            if OVERWRITE||~exist(filename,'file')
                fprintf('Downloading %s (repo %d/%d file %d/%d). This process may take several minutes. Please wait...\n',filename,n,numel(data),nfile,numel(repofiles_link));
                [fname,ok]=conn_loadgdfile(repofiles_link{nfile},filename,true);
            end
        end
    end
    if DODISPLAY, if ishandle(hmsg), delete(hmsg); end; end
end
NSUBJECTS=0;
for n=1:numel(data), NSUBJECTS=NSUBJECTS+repofiles(n).subjects; end

%% UNTAR .tar.gz files
if UNZIPFILES
    if DODISPLAY, hmsg=conn_msgbox('Uncompressing data files. Please wait...',''); end
    for n=1:numel(data),
        repofiles_link=repofiles(n).link;
        repofiles_file=repofiles(n).file;
        if ~iscell(repofiles_link), repofiles_link={repofiles_link}; end
        if ~iscell(repofiles_file), repofiles_file={repofiles_file}; end
        assert(numel(repofiles_link)==numel(repofiles_file));
        for nfile=1:numel(repofiles_link)
            filename=fullfile(cwd,repofiles_file{nfile});
            a=dir(filename);
            if ~isempty(a),
                [a_path,a_name,a_ext]=fileparts(fullfile(cwd,a.name));[nill,a_name2,a_ext2]=fileparts(a_name); % note: .tar or .tar.gz files
                dirname=fullfile(a_path,a_name2);
                if ~isdir(dirname),
                    disp(['extracting contents from file ',a.name]);
                    untar(a.name,dirname);
                end
                %% UNZIP .nii.gz files
                a=strvcat(conn_dir(fullfile(dirname,repofiles(n).func)),conn_dir(fullfile(dirname,repofiles(n).anat)));
                for n1=1:size(a,1),
                    [a_path,a_name,a_ext]=fileparts(a(n1,:));
                    if isempty(dir(fullfile(a_path,a_name))),
                        disp(['unzipping file ',a(n1,:)]);
                        gunzip(deblank(a(n1,:)));
                    end
                end
            end
        end
    end
    if DODISPLAY, if ishandle(hmsg), delete(hmsg); end; end
end

%% FIND functional/structural files
cwd=pwd;
sFUNCTIONAL_FILE={};
sSTRUCTURAL_FILE={};
SUB_UIDS=[];
DIR_NAMES={};
SUB_DEMODATASET={};
SUB_DEMOGRAPHICS={};
SUB_DATASET=[];
for n=1:numel(data),
    repofiles_link=repofiles(n).link;
    repofiles_file=repofiles(n).file;
    if ~iscell(repofiles_link), repofiles_link={repofiles_link}; end
    if ~iscell(repofiles_file), repofiles_file={repofiles_file}; end
    assert(numel(repofiles_link)==numel(repofiles_file));
    for nfile=1:numel(repofiles_link)
        filename=fullfile(cwd,repofiles_file{nfile});
        a=dir(filename);
        if ~isempty(a),
            [a_path,a_name,a_ext]=fileparts(fullfile(cwd,a.name));[nill,a_name2,a_ext2]=fileparts(a_name); % note: .tar or .tar.gz files
            dirname=fullfile(a_path,a_name2);
            dirsubs=dir(fullfile(dirname,'sub*'));
            dirsubs=dirsubs([dirsubs.isdir]>0);
            for ndirsubs=1:numel(dirsubs)
                ename=[data{n},'.',dirsubs(ndirsubs).name];
                [ok,NSUB]=ismember(ename,DIR_NAMES);
                if ~ok,
                    DIR_NAMES{end+1}=ename;
                    NSUB=numel(DIR_NAMES);
                end
                tFUNCTIONAL_FILE=conn_dir(fullfile(dirname,dirsubs(ndirsubs).name,regexprep(repofiles(n).func,'\.gz$','')),'-cell');
                tSTRUCTURAL_FILE=conn_dir(fullfile(dirname,dirsubs(ndirsubs).name,regexprep(repofiles(n).anat,'\.gz$','')),'-cell');
                for n1=1:numel(tFUNCTIONAL_FILE)
                    jsonfile=conn_prepend('',tFUNCTIONAL_FILE{n1},'.json');
                    if ~conn_existfile(jsonfile), conn_fileutils('spm_jsonwrite',jsonfile,struct('RepetitionTime',repofiles(n).tr,'SliceOrder',repofiles(n).sliceorder),struct('indent',' ')); end
                end
                if ~isempty(tFUNCTIONAL_FILE)&&~isempty(tSTRUCTURAL_FILE) % note: skip subjects without functional or anatomical data
                    sFUNCTIONAL_FILE=[sFUNCTIONAL_FILE;{reshape(tFUNCTIONAL_FILE,1,[])}];
                    sSTRUCTURAL_FILE=[sSTRUCTURAL_FILE;{reshape(tSTRUCTURAL_FILE,1,[])}];
                    SUB_UIDS=[SUB_UIDS;NSUB];
                    SUB_DATASET=[SUB_DATASET;n];
                end
            end
            dirdemographics=dir(fullfile(dirname,'*demographics.txt'));
            if ~isempty(dirdemographics), 
                SUB_DEMODATASET{n}=data{n}; 
                SUB_DEMOGRAPHICS{n}=conn_loadcsvfile(fullfile(dirname,dirdemographics(1).name),false); 
            end
        end
    end
end
[uSUB_UIDS,uSUB_idx]=unique(SUB_UIDS);
NSUBJECTS=numel(uSUB_UIDS);
FUNCTIONAL_FILE={};
STRUCTURAL_FILE={};
tIDs={};
tDATASET=[];
for nsub=1:NSUBJECTS,
    idx=find(SUB_UIDS==uSUB_UIDS(nsub));
    FUNCTIONAL_FILE{nsub}=cat(2,sFUNCTIONAL_FILE{idx});
    STRUCTURAL_FILE{nsub}=cat(2,sSTRUCTURAL_FILE{idx});
    tIDs{nsub}=DIR_NAMES{uSUB_UIDS(nsub)}; % note: extended IDs (e.g. 'NewYork_a.sub02503') is used to sort subjects
    tDATASET(nsub)=SUB_DATASET(uSUB_idx(nsub));
end
COVS=struct;
for n=1:numel(SUB_DEMOGRAPHICS), % matches covariates field sub##### id's to directory names
    if ~isempty(SUB_DEMOGRAPHICS{n})&&isstruct(SUB_DEMOGRAPHICS{n})
        fnames=fieldnames(SUB_DEMOGRAPHICS{n});
        NF=1:numel(fnames);
        COVS_ID=[];
        for nf=NF(NF>0),
            try, if iscell(SUB_DEMOGRAPHICS{n}.(fnames{nf})) && all(cellfun(@(x)~isempty(regexp(x,'^sub\d+$')),SUB_DEMOGRAPHICS{n}.(fnames{nf}))), COVS_ID=cellfun(@(x)[SUB_DEMODATASET{n},'.',x],SUB_DEMOGRAPHICS{n}.(fnames{nf}),'uni',0); NF(nf)=0; end; end
        end
        if ~isempty(COVS_ID)
            [ok,idx]=ismember(tIDs,COVS_ID); % match COVS_ID from file tu tIDs from directory names
            if any(ok),
                if ~isfield(COVS,'ID'), COVS.ID=nan(NSUBJECTS,1); end
                COVS.ID(ok)=str2double(regexprep(COVS_ID(idx(ok)),'^.*sub',''));
                for nf=NF(NF>0),
                    try, if iscell(SUB_DEMOGRAPHICS{n}.(fnames{nf})) && all(ismember(lower(SUB_DEMOGRAPHICS{n}.(fnames{nf})),{'f','m'})),
                            if ~isfield(COVS,'MALE'), COVS.MALE=nan(NSUBJECTS,1); end
                            if ~isfield(COVS,'FEMALE'), COVS.FEMALE=nan(NSUBJECTS,1); end
                            COVS.MALE(ok)=cellfun(@(x)strcmpi(x,'m'),SUB_DEMOGRAPHICS{n}.(fnames{nf})(idx(ok))); COVS.FEMALE(ok)=cellfun(@(x)strcmpi(x,'f'),SUB_DEMOGRAPHICS{n}.(fnames{nf})(idx(ok))); NF(nf)=0;
                    end; end
                try, if iscell(SUB_DEMOGRAPHICS{n}.(fnames{nf})) && all(cellfun(@(x)~isempty(regexp(x,'^r|^l')),lower(SUB_DEMOGRAPHICS{n}.(fnames{nf})))),
                        if ~isfield(COVS,'HANDRIGHT'), COVS.HANDRIGHT=nan(NSUBJECTS,1); end
                        if ~isfield(COVS,'HANDLEFT'), COVS.HANDLEFT=nan(NSUBJECTS,1); end
                        COVS.HANDRIGHT(ok)=cellfun(@(x)~strcmpi(x,'l'),SUB_DEMOGRAPHICS{n}.(fnames{nf})(idx(ok)));
                        COVS.HANDLEFT(ok)=cellfun(@(x)~strcmpi(x,'r'),SUB_DEMOGRAPHICS{n}.(fnames{nf})(idx(ok))); NF(nf)=0;
                end; end
            try, if ~iscell(SUB_DEMOGRAPHICS{n}.(fnames{nf})),
                    if ~isfield(COVS,'AGE'), COVS.AGE=nan(NSUBJECTS,1); end
                    COVS.AGE(ok)=SUB_DEMOGRAPHICS{n}.(fnames{nf})(idx(ok)); COVS.AGE(COVS.AGE==9999)=nan; NF(nf)=0;
            end; end
                end
            end
        end
    end
end

nsessions=cellfun('length',FUNCTIONAL_FILE);
disp([num2str(numel(FUNCTIONAL_FILE)),' subjects']);
if all(nsessions==mean(nsessions)), disp([num2str(mean(nsessions)),' runs']);
else disp(['between ',num2str(min(nsessions)),' and ',num2str(max(nsessions)),' runs']);
end

if ~DOCREATE,
    fprintf('Data downloaded to %s\n',cwd);
    return;
end

if DODISPLAY, 
    if DOPROCESS, hmsg=conn_msgbox('Importing and processing data. Please wait...',''); 
    else hmsg=conn_msgbox('Importing data. Please wait...',''); 
    end
end

%% CONN-SPECIFIC SECTION: RUNS PREPROCESSING/SETUP/DENOISING/ANALYSIS STEPS
%% Prepares batch structure
clear batch;
batch.filename=fullfile(cwd,PROJECTNAME);            % New conn_*.mat experiment name
if DOPARALLEL,
    batch.parallel.N=1;
    batch.parallel.immediatereturn=true;
    if ~isempty(PARALLELPROFILE), batch.parallel.profile=PARALLELPROFILE; end
end

%% SETUP step (using default values for most parameters, see help conn_batch to define non-default values)
% CONN Setup                                            % Default options (uses all ROIs in conn/rois/ directory); see conn_batch for additional options
% CONN Setup.preprocessing                               (realignment/coregistration/segmentation/normalization/smoothing)
batch.Setup.isnew=1;
batch.Setup.nsubjects=NSUBJECTS;
batch.Setup.RT=NaN;                                       % TR (seconds)
batch.Setup.functionals=FUNCTIONAL_FILE;                  % Point to functional volumes for each subject/session
batch.Setup.structurals=STRUCTURAL_FILE;                  % Point to anatomical volumes for each subject
nconditions=max(nsessions);                               % treats each session as a different condition (comment the following three lines and lines 84-86 below if you do not wish to analyze between-session differences)
if nconditions==1
    batch.Setup.conditions.names={'rest'};
    for ncond=1,for nsub=1:NSUBJECTS,for nses=1:nsessions(nsub),              batch.Setup.conditions.onsets{ncond}{nsub}{nses}=0; batch.Setup.conditions.durations{ncond}{nsub}{nses}=inf;end;end;end     % rest condition (all sessions)
else
    batch.Setup.conditions.names=[{'rest'}, arrayfun(@(n)sprintf('Session%d',n),1:nconditions,'uni',0)];
    for ncond=1,for nsub=1:NSUBJECTS,for nses=1:nsessions(nsub),              batch.Setup.conditions.onsets{ncond}{nsub}{nses}=0; batch.Setup.conditions.durations{ncond}{nsub}{nses}=inf;end;end;end     % rest condition (all sessions)
    for ncond=1:nconditions,for nsub=1:NSUBJECTS,for nses=1:nsessions(nsub),  batch.Setup.conditions.onsets{1+ncond}{nsub}{nses}=[];batch.Setup.conditions.durations{1+ncond}{nsub}{nses}=[]; end;end;end
    for ncond=1:nconditions,for nsub=1:NSUBJECTS,if nsessions(nsub)>=ncond,   batch.Setup.conditions.onsets{1+ncond}{nsub}{ncond}=0; batch.Setup.conditions.durations{1+ncond}{nsub}{ncond}=inf;end;end;end % session-specific conditions
end
batch.Setup.subjects.effects={};
batch.Setup.subjects.effect_names={};
batch.Setup.subjects.add=true;
uDATASET=unique(tDATASET);
for nd=1:numel(uDATASET),
    batch.Setup.subjects.effects{end+1}=tDATASET(:)==uDATASET(nd);
    batch.Setup.subjects.effect_names{end+1}=sprintf('SITE_%s',data{uDATASET(nd)});
end
fCOVS=fieldnames(COVS);
if ~isempty(fCOVS)
    for nfCOVS=1:numel(fCOVS), 
        batch.Setup.subjects.effects{end+1}=COVS.(fCOVS{nfCOVS}); 
        batch.Setup.subjects.effect_names{end+1}=fCOVS{nfCOVS};
    end
end

if DOPROCESS
    %% PREPROCESSING step (using default values for most parameters, see help conn_batch to define non-default values)
    batch.Setup.preprocessing.steps='default_mni';    % default preprocessing pipeline
    batch.Setup.preprocessing.sliceorder='BIDS';      % info in .json files
    batch.Setup.preprocessing.art_thresholds=[3 0.5]; % conservative settings
    %batch.Setup.preprocessing.rmask=0;                % skip masking during realignment
    batch.Setup.done=1;
    batch.Setup.overwrite='Yes';

    %% DENOISING step
    % CONN Denoising                                    % Default options (uses White Matter+CSF+realignment+scrubbing+conditions as confound regressors); see conn_batch for additional options
    batch.Denoising.done=1;
    batch.Denoising.overwrite='Yes';

    %% QC plots
    batch.QA.plots=[1,2,4,5,7,9,11,12,13,31];           % Default Quality Control plots

    %% FIRST-LEVEL ANALYSIS step
    % CONN Analysis                                     % Default options (uses all network ROIs in conn/rois/networks.nii as connectivity sources); see conn_batch for additional options
    batch.Analysis.sources={'networks'};
    batch.Analysis.done=1;
    batch.Analysis.overwrite='Yes';
end

%% Run all analyses
conn_batch(batch);
if DODISPLAY, if ishandle(hmsg), delete(hmsg); end; end

%% CONN Display
% launches conn gui to explore results
if DOPARALLEL
    fprintf('Process running in parallel. Check its status on Tools -> HPC options -> Active/pending jobs"\n')
    conn
    conn gui_setup
elseif DODISPLAY
    conn
    conn('load',fullfile(cwd,PROJECTNAME));
    conn gui_setup
end
end

