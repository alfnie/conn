function varargout=conn_module(option,varargin)
% CONN_MODULE provides access to independent CONN modules 
%
% conn_module(module_name, ...) runs individual CONN's module "module_name" on user-defined data
%
% Current module names: PREP, GLM, EL, FL
%
%    PREP : runs CONN preprocessing pipeline on user-defined data (see www.conn-toolbox.org/resources/prep for details)
%
%       basic syntax: conn_module preprocessing
%       advanced syntax: conn_module('PREP', fieldname1, fieldvalue1, fieldname2, fieldvalue2, ...)
%
%                 Input data is specified with field name/value pairs as defined in batch.Setup documentation
%                 Preprocessing options are specified with field name/value pairs as defined in batch.Setup.preprocessing documentation
%                    functionals       : list of functional data files { { Sub1Ses1, Sub1Ses2, ...}, {Sub2Ses1, Sub2Ses2, ...}, ...}
%                    structurals       : list of structural data files { Sub1, Sub2, ...}
%                    steps             : list of preprocessing steps (tpye "conn_module PREP steps" for a list of valid preprocessing step names)
%
%                 See "doc conn_batch" for a complete list of these options
%                 See Nieto-Castanon, 2020 for details about these preprocessing steps and pipelines (www.conn-toolbox.org/fmri-methods)
%
%       e.g. conn_module('PREP',...
%             'steps','default_mni',...
%             'functionals',{'./func.nii'},...
%             'structurals',{'./anat.nii'},...
%             'RT',2,...
%             'sliceorder','interleaved (Siemens)');
%            runs default MNI-space preprocessing pipeline on the specified functional/structural data
%
%       alternative syntax: conn_module('PREP',optionsfile) 
%          input data and preprocessing options defined in .cfg (see conn_loadcfgfile/conn_savecfgfile) or .json (see spm_jsonread/spm_jsonwrite) structure text file 
%
%       alternative syntax: conn_module preprocessing steps
%          returns the list of valid preprocessing steps
%
%       See https://web.conn-toolbox.org/resources/conn-extensions/prep for additional options
%
%    GLM : runs CONN second-level analyses on user-defined data (see www.conn-toolbox.org/resources/glm for details)
%
%       basic syntax: conn_module GLM
%       advanced syntax: conn_module('GLM', fieldname1, fieldvalue1, fieldname2, fieldvalue2, ...)
%               with the following field name/value pairs
%                  data            : list of nifti files entered into second-level analysis (Nsubjects x Nmeasures) defining one or multiple outcome / dependent measures 
%                                       note: when entering multiple files per subject (e.g. repeated measures) enter first all files (one per subject) for measure#1, followed by all files for measure#2, etc.
%                                       note: nifti files may contani 3d volume-level data, fsaverage surface-level data, or ROI-to-ROI data (see conn_surf_write and conn_mtx_write to create surface/matrix nifti files)
%                                       note: enter a single 4d file or a filename with wildcards (e.g. vol_*.nii) to simplify the specification of multiple files 
%                                    alternatively list of SPM.mat files containing first-level analyses (Nsubjects x 1, or Nsubjects x Nmeasures)
%                                    alternatively list of folder names containing SPM.mat first-level analyses (Nsubjects x 1, or Nsubjects x Nmeasures)
%                  design_matrix   : design matrix (Nsubjects x Neffects) defining different explanatory / independent measures or subject-effects 
%                                       enter one row for each subject
%                                       each row should contain one value/number per modeled effect/covariate
%                  contrast_between: between-subjects contrast vector/matrix (Nc1 x Neffects) 
%                  contrast_within : within-subjects contrast vector/matrix (Nc2 x Nmeasures)
%                  contrast_names  : (optional, only when entering SPM.mat files in #data field) list of contrast names to select from first-level analysis files (Nmeasures x 1)
%                  data_labels     : (optional) labels of columns of data matrix
%                  design_labels   : (optional) labels of columns of design matrix
%                  mask            : (optional) analysis-mask file
%                  analysistype    : (optional) analysis type 1: include both parametric and non-parametric stats; 2: include only parametric stats (Random Field Theory assumptions); 3: include only non-parametric stats (permutation/randomization analyses)
%                  folder          : (optional) folder where analysis are stored; default current folder
%                  design_file     : (optional) (alternative to design_matrix field) file containing design_matrix data (.csv/.tsv/.txt/.json/.mat with Nsubjects x Neffects matrix of explanatory / independent measures)
%                  design          : (optional) (alternative to design_matrix field) transpose of design_matrix (Neffects x Nsubjects); enter one row for each modeled effect (across subjects); each row should contain one value/number per subject
%
%         eg: conn_module('GLM', ...
%            'design_matrix',[1; 1; 1; 1] ,...
%            'data',{'subject1.img'; 'subject2.img'; 'subject3.img'; 'subject4.img'} );
%             performs a one-sample t-test and stores the analysis results in the current folder
%
%         eg: conn_module('GLM', ...
%            'design_matrix',[1 0; 1 0; 0 1; 0 1; 0 1],...
%            'data', {'subject1_group1.img'; 'subject2_group1.img'; 'subject1_group2.img'; 'subject2_group2.img'; 'subject3_group2.img'},...
%            'contrast_between',[1 -1]);
%             performs a two-sample t-test and stores the analysis results in the current folder
%
%         eg: conn_module('GLM', ...
%            'design_matrix', [1; 1; 1; 1],...
%            'data', {'subject1_time1.img', subject1_time2.img'; 'subject2_time1.img', subject2_time2.img'; 'subject3_time1.img', subject3_time2.img'; 'subject4_time1.img', subject4_time2.img'},...
%            'contrast_beetween',1,...
%            'contrast_within',[1 -1]);
%             performs a paired t-test and stores the analysis results in the current folder
%
%       alternative syntax: conn_module('GLM',optionsfile) 
%          input data and GLM options defined in .cfg (see conn_loadcfgfile/conn_savecfgfile) or .json (see spm_jsonread/spm_jsonwrite) structure text file 
%
%       alternative syntax: spmfolder=conn_module('GLM',...) 
%          skips results display step (only computes second-level analysis, and returns folder where results are stored)
%          use conn_display(spmfolder) syntax to then launch the results explorer window on previously computed analyses
%
%       See also CONN_DISPLAY for displaying GLM results
%       See Nieto-Castanon, 2020 for details about General Linear Model analyses (www.conn-toolbox.org/fmri-methods)
%       See https://web.conn-toolbox.org/resources/conn-extensions/glm for additional options
%
%    EL : runs EvLab (evlab.mit.edu) fMRI pipeline for subject-centric task-activation analyses
%       See https://web.conn-toolbox.org/resources/conn-extensions/el for descriptiona and options
% 
%    FL : runs FrankLab (sites.bu.edu/guentherlab) fMRI pipeline for group-centric task-activation analyses (see www.conn-toolbox.org/resources/franklab for details)
%       See https://web.conn-toolbox.org/resources/conn-extensions/fl for description and options
% 
%    Additional functionality: conn_module('get',...)
%          conn_module('get','structurals');             outputs current structural files (e.g. output of structural preprocessing steps)
%          conn_module('get','functionals' [,setlabel]); outputs current functional files (e.g. output of functional preprocessing steps)
%          conn_module('get','l1covariates' [,covname]); outputs first-level covariate files (e.g. other potential outputs of functional preprocessing)
%          conn_module('get','l2covariates' [,covname]); outputs second-level covariate values (e.g. other potential outputs of functional preprocessing)
%          conn_module('get','masks' [,roiname]);        outputs Grey Matter/White Matter/CSF files (e.g. other potential outputs of functional preprocessing)
%          conn_module('get','rois' [,roiname]);         outputs ROI files
%    Additional functionality: conn_module('set',...)
%          conn_module('set','l1covariates',files,covname [,add]);
%          conn_module('set','l2covariates',values,covname [,covdescrip ,add]);
%          conn_module('set','masks',files [,roiname]);
%          conn_module('set','rois',files ,roiname);
%
%    Note: before using conn_module functionality with externally defined data it is recommended to close CONN's gui in order to avoid potentially loosing any unsaved changes
%    

persistent defaults modules modulespath cpm_lastjob;

if isempty(defaults), defaults=struct('mat_format','-v7.3'); end
if isempty(modules), 
    LETOVERLOAD=false; % true: allows module functions to be overloaded by same-name functions within folders higher in Matlab path
                       % (set to 'true' in multiuser environments where you would like a common CONN version but possibly user-specific versions of module functions)
    modulespath=conn_dir(fullfile(fileparts(which(mfilename)),'modules','*'),'^[^\.].*','-dir','-cell','-R'); 
    modules=regexprep(modulespath,'^.*[\\\/]',''); 
    for n=1:numel(modules), 
        wmodules=fileparts(which(modules{n}));
        if LETOVERLOAD&&isempty(wmodules), addpath(modulespath{n});
        elseif ~LETOVERLOAD&&~isequal(wmodules,modulespath{n}), 
            if ~isempty(wmodules), fprintf('warning: there is an existing version of module %s in your path (%s), overloading with version in %s. Please remove original version from your Matlab path to stop seeing this message in the future\n',modules{n},wmodules,modulespath{n}); end
            addpath(modulespath{n});
        end
    end
end
if ~nargin, help(mfilename); return; end

varargout={[]};
switch(lower(option))
    
    case 'get'
        names={};files={};other={};
        switch(lower(varargin{1}))
            case 'structurals',
                files=conn('get','Setup.structural');
                for nsub=1:numel(files),
                    if conn('get','Setup.structural_sessionspecific'),
                        for nses=1:numel(files{nsub}),
                            files{nsub}{nses}=files{nsub}{nses}{1};
                        end
                    else
                        files{nsub}=files{nsub}{1}{1};
                    end
                end
            case 'functionals',
                if numel(varargin)>1, nset=varargin{2}; 
                else nset=0;
                end
                files=conn('get','Setup.functional');
                for nsub=1:numel(files),
                    for nses=1:numel(files{nsub}),
                        files{nsub}{nses}=conn_get_functional(nsub,nses,nset); 
                        %if ~nset, files{nsub}{nses}=conn_get_functional(nsub,nses,nset); 
                        %else files{nsub}{nses}=files{nsub}{nses}{1};
                        %end
                    end
                end
            case {'spm','dicom','bids'}
                files=conn('get',['Setup.',(varargin{1})]);
                for nsub=1:numel(files),
                    files{nsub}=files{nsub}{1};
                end
            case 'l1covariates',
                data=conn('get','Setup.l1covariates');
                files={};
                for nsub=1:numel(data.files),
                    for ncov=1:numel(data.names)-1,
                        for nses=1:numel(data.files{nsub}{ncov}),
                            files{ncov}{nsub}{nses}=char(data.files{nsub}{ncov}{nses}{1});
                            if isequal(files{ncov}{nsub}{nses},'[raw values]'), files{ncov}{nsub}{nses}=data.files{nsub}{ncov}{nses}{3}; end
                        end
                    end
                end
                for ncov=1:numel(data.names)-1,
                    names{ncov}=data.names{ncov};
                end
                if numel(varargin)>1,
                    icov=find(ismember(names,cellstr(varargin{2})));
                    if numel(icov)==1,
                        names=names{icov};
                        files=files{icov};
                    else
                        names=names(icov);
                        files=files(icov);
                    end
                end
            case 'l2covariates',
                data=conn('get','Setup.l2covariates');
                files=[];
                for nsub=1:numel(data.values),
                    for ncov=1:numel(data.values{nsub}),
                        files(nsub,ncov)=data.values{nsub}{ncov};
                    end
                end
                for ncov=1:numel(data.names)-1,
                    names{ncov}=data.names{ncov};
                    other{ncov}=data.descrip{ncov};
                end
                if numel(varargin)>1,
                    [ok,idx]=ismember(names,cellstr(varargin{2}));
                    icov=find(ok);
                    nmatch=numel(icov);
                    if ~nmatch, icov=find(cellfun('length',regexp(names,varargin{2}))); end
                    if nmatch==1,
                        names=names{icov};
                        other=other{icov};
                        files=files(:,icov);
                    else
                        [ok,idx]=ismember(cellstr(varargin{2}),names(icov));
                        if all(ok), icov=icov(idx); end
                        names=names(icov);
                        other=other(ncov);
                        files=files(:,icov);
                    end
                end
            case {'masks','rois'}
                data=conn('get','Setup.rois');
                files={};
                if strcmpi(varargin{1},'rois'), maxrois=inf; else maxrois=3; end
                for nsub=1:numel(data.files),
                    for nroi=1:min(maxrois, numel(data.files{nsub})),
                        for nses=1:numel(data.files{nsub}{nroi}),
                            files{nroi}{nsub}{nses}=data.files{nsub}{nroi}{nses}{1};
                        end
                    end
                end
                for nroi=1:min(maxrois, numel(data.names)-1),
                    names{nroi}=data.names{nroi};
                end
                if numel(varargin)>1,
                    iroi=find(ismember(names,cellstr(varargin{2})));
                    if numel(iroi)==1,
                        names=names{iroi};
                        files=files{iroi};
                    else
                        names=names(iroi);
                        files=files(iroi);
                    end
                end
            otherwise
                files=conn('get',varargin{:});
        end
        varargout={files,names,other};

    case 'set'
        switch(lower(varargin{1}))
            case 'structurals',
                conn_batch('Setup.structurals',varargin{2});
            case 'functionals',
                conn_batch('Setup.functionals',varargin{2});
            case {'spm','dicom','bids'}
                files=varargin{2};
                data.files={};
                for nsub=1:numel(files),
                    data.files{nsub}=conn_file(files{nsub});
                end
                conn('set',['Setup.',(varargin{1})],data.files);
            case 'l1covariates',
                files=varargin{2};
                names=varargin{3};
                if numel(varargin)>3, add=varargin{4};
                else add=false;
                end
                if ~iscell(names),
                    names={names};
                    files={files};
                end
                if add
                    [files0,names0]=conn_module('get','l1covariates');
                    if ischar(names0), names0={names0}; files0={files0}; end
                    [ok,idx]=ismember(names,names0);
                    files0(idx(ok))=files(ok);
                    files=[files0(:);reshape(files(~ok),[],1)];names=[names0(:);reshape(names(~ok),[],1)];
                end
                data.names={};
                data.files={};
                for ncov=1:numel(files),
                    for nsub=1:numel(files{ncov}),
                        for nses=1:numel(files{ncov}{nsub}),
                            if isempty(files{ncov}{nsub}{nses})&&isnumeric(files{ncov}{nsub}{nses}), data.files{nsub}{ncov}{nses}={'[raw values]',[],files{ncov}{nsub}{nses}};
                            else data.files{nsub}{ncov}{nses}=conn_file(files{ncov}{nsub}{nses});
                            end
                        end
                    end
                end
                for ncov=1:numel(names),
                    data.names{ncov}=names{ncov};
                end
                data.names{numel(names)+1}=' ';
                conn('set','Setup.l1covariates',data);
            case 'l2covariates',
                files=varargin{2};
                names=varargin{3};
                if numel(varargin)>3, descrip=varargin{4}; 
                else descrip={}; 
                end
                if numel(varargin)>4, add=varargin{5};
                else add=false;
                end
                if ~iscell(names),
                    names={names};
                end
                if add
                    [files0,names0]=conn_module('get','l2covariates');
                    if ischar(names0), names0={names0}; end
                    [ok,idx]=ismember(names,names0);
                    files0(:,idx(ok))=files(:,ok);
                    files=[files0,files(:,~ok)];names=[names0(:);reshape(names(~ok),[],1)];
                end
                data_prev=conn('get','Setup.l2covariates');
                data.names={};
                data.values={};
                data.descrip={};
                for ncov=1:size(files,2),
                    for nsub=1:size(files,1),
                        data.values{nsub}{ncov}=files(nsub,ncov);
                    end
                end
                for ncov=1:numel(names),
                    data.names{ncov}=names{ncov};
                    if numel(descrip)<ncov, data.descrip{ncov}='';
                    else data.descrip{ncov}=descrip{ncov};
                    end
                end
                try
                    [ok,idx]=ismember(data.names,data_prev.names);
                    ok(cellfun('length',descrip)>0)=false;
                    data.descrip(ok)=data_prev.descrip(idx(ok));
                end
                data.names{numel(names)+1}=' ';
                conn('set','Setup.l2covariates',data);
            case {'masks','rois'}
                files=varargin{2};
                names=varargin{3};
                if ~iscell(names),
                    names={names};
                    files={files};
                end
                data.names={};
                data.files={};
                for nroi=1:min(3, numel(files)),
                    for nsub=1:numel(files{nroi}),
                        for nses=1:numel(files{nroi}{nsub}),
                            data.files{nsub}{nroi}{nses}=conn_file(files{nroi}{nsub}{nses});
                        end
                    end
                end
                for nroi=1:min(3, numel(names)),
                    data.names{nroi}=names{nroi};
                end
                data.names{numel(names)+1}=' ';
                conn('set','Setup.rois.files',data.files);
                conn('set','Setup.rois.names',data.names);
            otherwise
                if nargout>0, [varargout{1:nargout}]=conn('set',varargin{:});
                else conn('set',varargin{:});
                end
        end
        
    case {'prep','preprocessing'}
        options=struct;
        if isempty(varargin)
            answ=conn_menu_inputdlg('Enter number of subjects','conn_module PREP',1,{'1'});
            if isempty(answ), return; end
            options.nsubjects=str2num(answ{1});
            options.functionals=cell(options.nsubjects,1);
            nsessions=ones(1,options.nsubjects);
            for nsubject=1:options.nsubjects,
                temp=conn_menu_inputdlg(['Subject ',num2str(nsubject),': Enter number of runs/sessions'],'conn_module PREP',1,{num2str(nsessions(min(length(nsessions),nsubject)))});
                if isempty(temp), return; end
                temp=str2num(temp{1});
                nsessions(nsubject)=temp;
                options.functionals{nsubject}=cell(nsessions(min(length(nsessions),nsubject)),1);
                for nsession=1:nsessions(min(length(nsessions),nsubject)),
                    options.functionals{nsubject}{nsession}=cellstr(spm_select(Inf,'\.img$|\.nii$',['SUBJECT ',num2str(nsubject),' SESSION ',num2str(nsession),' functional volumes'],options.functionals{nsubject}{nsession}));
                    if isempty(options.functionals{nsubject}{nsession}{1}),return;end
                end
            end
            options.structurals=cell(options.nsubjects,1);
            for nsubject=1:options.nsubjects,
                options.structurals{nsubject}=cellstr(spm_select(1,'\.img$|\.nii$',['SUBJECT ',num2str(nsubject),' structural volume'],options.structurals{nsubject}));
                if isempty(options.structurals{nsubject}{1}),return;end
            end
            options.steps='';
            options.multiplesteps=1;
        elseif isstruct(varargin{1})
            % syntax: conn_module('PREP',struct('data',...),...)
            n=1;
            options=varargin{n};
        elseif nargin==2&&isequal(varargin{1},'steps')
            if nargout>0, varargout{1}=conn_setup_preproc('steps');
            else
                disp('List of valid preprocessing steps:');
                disp(char(conn_setup_preproc('steps')));
            end
            return
        elseif nargin==2
            % syntax: conn_module PREP optionsfile.cfg
            n=1;
            if ~isempty(varargin{n}), 
                filename=varargin{n};
                if ischar(filename)
                    if ~conn_existfile(filename),
                        if conn_existfile(fullfile(fileparts(which(mfilename)),filename)), filename=fullfile(fileparts(which(mfilename)),filename);
                        else error('file %s not found',filename);
                        end
                    end
                    fprintf('loading file %s\n',filename);
                    [nill,nill,fext]=fileparts(filename);
                    if strcmp(fext,'.json')
                        options=conn_jsonread(filename);
                    else
                        options=conn_loadcfgfile(filename,options);
                        if isfield(options,'functionals')
                            if ischar(options.functionals), options.functionals=cellstr(options.functionals); end
                            if numel(options.functionals)>1, options.functionals={options.functionals}; end; % note: single-subject data in each .cfg file
                        end
                        if isfield(options,'structurals')
                            if ischar(options.structurals), options.structurals=cellstr(options.structurals); end
                            if numel(options.structurals)>1, options.structurals={options.structurals}; end; % note: single-subject data in each .cfg file
                        end
                    end
                else
                    options=conn_loadcfgfile(filename,options);
                end
            end
        else
            % syntax: conn_module PREP fieldname1 fieldvalue1 ...
            n=0;
        end
        for n=n+1:2:numel(varargin)-1
            fieldname=regexp(varargin{n},'\.','split');
            fieldvalue=varargin{n+1};
            options=setfield(options,fieldname{:},fieldvalue);
        end
        Batch=struct;
        if isfield(options,'functionals')&&~isfield(options,'nsubjects'), 
            if iscell(options.functionals), options.nsubjects=numel(options.functionals); 
            else options.nsubjects=1;
            end
        end
        if isfield(options,'structurals')&&~isfield(options,'nsubjects'), 
            if iscell(options.structurals), options.nsubjects=numel(options.structurals); 
            else options.nsubjects=1;
            end
        end
        if isfield(options,'parallel')&&~isfield(options,'filename'),
            str=fullfile(pwd,sprintf('CONN_module_%s%s.mat',datestr(now,'dd-mmm-yyyy-HHMMSSFFF'),char(floor('0'+10*rand(1,8)))));
            Batch.filename=str;
            Batch.Setup.isnew=true;
        elseif isfield(options,'nsubjects')&&~isfield(options,'filename')
            conn('init');
            conn('set','filename','');
        end
        if isfield(options,'functionals')
            if ~iscell(options.functionals), options.functionals={options.functionals}; end
            for nsub=1:numel(options.functionals),
                if ~iscell(options.functionals{nsub}), options.functionals{nsub}={options.functionals{nsub}}; end
            end
        end
        if isfield(options,'structurals')
            if ~iscell(options.structurals), options.structurals={options.structurals}; end
            for nsub=1:numel(options.structurals),
                if ~iscell(options.structurals{nsub}), options.structurals{nsub}={options.structurals{nsub}}; end
            end
        end
        Fields={'functionals','structurals','secondarydatasets','nsubjects','RT','covariates','unwarp_functionals','vdm_functionals','fmap_functionals','coregsource_functionals','masks','rois','localcopy','localcopy_reduce','isnew'}; % send these to Setup
        for n=1:numel(Fields)
            if isfield(options,Fields{n}),
                Batch.Setup.(Fields{n})=options.(Fields{n});
                options=rmfield(options,Fields{n});
            elseif isfield(options,lower(Fields{n})),
                Batch.Setup.(Fields{n})=options.(lower(Fields{n}));
                options=rmfield(options,lower(Fields{n}));
            end
        end
        Fields=fieldnames(options);
        for n=1:numel(Fields)
            if any(strcmp(Fields{n},{'filename','parallel','subjects'})), % leave these in root
                Batch.(Fields{n})=options.(Fields{n});
                options=rmfield(options,Fields{n});
            else,                                                         % send others to Setup.preprocessing
                Batch.Setup.preprocessing.(Fields{n})=options.(Fields{n});
                options=rmfield(options,Fields{n});
            end
        end
        conn_batch(Batch);
        
    case {'cpm','pm'} % PREDICTIVE MODELS
        switch(lower(varargin{1}))
            case {'create','predict'}
                filename_project=varargin{2}; 
                options=struct('predictor',[],'outcome',[],'covariate',[],'fit',false,'null',false,'parallel',false,'nnull',100,'folder','');
                validfields=fieldnames(options);
                for n=3:2:numel(varargin)-1
                    fieldname=lower(varargin{n});
                    fieldvalue=varargin{n+1};
                    assert(isfield(options,fieldname),'incorrect field %s (valid field names: %s)',fieldname,sprintf('%s ',validfields{:}));
                    options=setfield(options,fieldname,fieldvalue);
                end
                if ischar(options.fit), options.fit=str2num(options.fit); end
                if ischar(options.null), options.null=str2num(options.null); end
                if ischar(options.parallel), options.parallel=str2num(options.parallel); end
                if ischar(options.nnull), options.nnull=str2num(options.nnull); end
                if isempty(options.folder), options.folder=pwd; end
                data_type='unknown';
                do_gui=(isempty(options.predictor)|isempty(options.outcome))&~options.parallel;
                if isempty(options.predictor)||ischar(options.predictor)
                    if isempty(options.predictor), 
                        disp('Select file containing PREDICTOR MEASURES'); 
                        [tfilename,tpathname]=uigetfile({ '*',  'All Files (*)'; '*.mtx.nii','matrix NIFTI files (*.mtx.nii)'; '*.surf.nii','surface NIFTI files (*.surf.nii)'; '*.nii','volume NIFTI files (*.nii)'},'Select file with PREDICTOR MEASURES');
                        if ~ischar(tfilename)||isempty(tfilename), return; end
                        options.predictor=fullfile(tpathname,tfilename);
                    end
                    if ~isempty(regexp(options.predictor,'\.mtx\.nii$'))
                        [data_predictor,info{1:3}]=conn_mtx_read(options.predictor); % info: names,coords,samples
                        data_predictor=reshape(data_predictor,[],size(data_predictor,3))'; % samples x connections
                        data_type='mtx';
                    elseif ~isempty(regexp(options.predictor,'\.surf\.nii$'))
                        data_predictor=conn_surf_read(options.predictor);
                        data_predictor=data_predictor'; % samples x vertices
                        data_type='surf';
                    else
                        [data_predictor,info]=conn_vol_read(options.predictor);
                        data_type='vol';
                        data_predictor=reshape(data_predictor,[],size(data_predictor,4))'; % samples x voxels
                    end
                else 
                    data_predictor=options.predictor;
                end
                if strcmp(lower(varargin{1}),'predict')
                    load(conn_prepend('',filename_project,'.mat'),'model');
                    if isempty(options.covariate), Yfit = model.predict(data_predictor);
                    else Yfit = model.predict(data_predictor,options.covariate);
                    end
                    if nargout>0, varargout={Yfit}
                    else disp(Yfit);
                    end
                    return
                end
                if isempty(options.outcome)||ischar(options.outcome)
                    if isempty(options.outcome)
                        disp('Select file containing OUTCOME MEASURE'); 
                        [tfilename,tpathname]=uigetfile({'*',  'All Files (*)';'*.mat','MATLAB file (*.mat)'; '*.csv','CSV file (*.csv)'},'Select file with OUTCOME MEASURE');
                        if ~ischar(tfilename)||isempty(tfilename), return; end
                        options.outcome=fullfile(tpathname,tfilename);
                    end
                    if ~isempty(regexp(options.outcome,'\.mat$'))
                        data=load(options.outcome);
                    else
                        data=conn_loadtextfile(options.outcome)
                    end
                    if isfield(data,'B')&&isfield(data,'Bnames') % for .mat files that contain a "B" field ([samples x variables] values) and a "Bnames" field (variable names)
                        fnames=data.Bnames;
                        idx=1:numel(fnames);
                        if numel(idx)>1
                            disp('Select outcome measure');
                            nselect=listdlg('liststring',fnames(idx),'selectionmode','single','promptstring',{'Select Outcome Measure'},'ListSize',[500 200]);
                            if isempty(nselect), return; end
                            idx=idx(nselect);
                        end
                        options.outcome=data.B(:,idx);
                    else % for .csv and for all other .mat files
                        fnames=fieldnames(data);
                        idx=find(cellfun(@(v)isnumeric(data.(v)),fnames));
                        if numel(idx)>1
                            disp('Select outcome measure');
                            nselect=listdlg('liststring',fnames(idx),'selectionmode','single','promptstring',{'Select Outcome Measure'},'ListSize',[500 200]);
                            if isempty(nselect), return; end
                            idx=idx(nselect);
                        end
                        options.outcome=data.(fnames{idx});
                    end
                end
                if size(options.outcome,1)~=size(data_predictor,1)&&size(options.outcome,2)==size(data_predictor,1), options.outcome=options.outcome'; end
                if ~isempty(options.covariate)&&ischar(options.covariate)
                    if isempty(options.covariate)
                        disp('Select file containing COVARIATE MEASURE'); 
                        [tfilename,tpathname]=uigetfile({ '*',  'All Files (*)'; '*.mat','MATLAB file (*.mat)'; '*.csv','CSV file (*.csv)'},'Select file with COVARIATE MEASURE');
                        if ~ischar(tfilename)||isempty(tfilename), return; end
                        options.covariate=fullfile(tpathname,tfilename);
                    end
                    if ~isempty(regexp(options.covariate,'\.mat$'))
                        data=load(options.covariate);
                    else
                        data=conn_loadtextfile(options.covariate)
                    end
                    if isfield(data,'B')&&isfield(data,'Bnames') % for .mat files that contain a "B" field ([samples x variables] values) and a "Bnames" field (variable names)
                        fnames=data.Bnames;
                        idx=1:numel(fnames);
                        if numel(idx)>1
                            disp('Select covariate measure');
                            nselect=listdlg('liststring',fnames(idx),'selectionmode','multiple','promptstring',{'Select Covariate Measure'},'ListSize',[500 200]);
                            if isempty(nselect), return; end
                            idx=idx(nselect);
                        end
                        options.covariate=data.B(:,idx);
                    else
                        fnames=fieldnames(data);
                        idx=find(cellfun(@(v)isnumeric(data.(v)),fnames));
                        if numel(idx)>1
                            disp('Select covariate measure');
                            nselect=listdlg('liststring',fnames(idx),'selectionmode','multiple','promptstring',{'Select Covariate Measure(s)'},'ListSize',[500 200]);
                            if isempty(nselect), return; end
                            idx=idx(nselect);
                        end
                        options.covariate=[]; for n=1:numel(idx), options.covariate=[options.covariate, data.(fnames{idx(n)})]; end
                    end
                end
                if do_gui
                    thfig=conn_figure('units','norm','position',[.4,.4,.2,.3],'color',1*[1 1 1],'name','CPM options','numbertitle','off','menubar','none');
                    uicontrol('style','text','units','norm','position',[.05,.7,.9,.10],'string','Other options (select all that apply)','backgroundcolor',1*[1 1 1]);
                    ht1=uicontrol('style','checkbox','units','norm','position',[.2,.55,.6,.10],'string','Compute fit values','value',options.fit,'backgroundcolor',1*[1 1 1],'tooltipstring','Computes fit values in training dataset using nested crossvalidation');
                    ht2=uicontrol('style','checkbox','units','norm','position',[.2,.45,.6,.10],'string','Compute p values','value',options.null,'backgroundcolor',1*[1 1 1],'tooltipstring','Computes null-hypothesis distribution of MSE scores using permutation analyses');
                    ht3=uicontrol('style','checkbox','units','norm','position',[.2,.35,.6,.10],'string','Run in background','value',options.parallel,'backgroundcolor',1*[1 1 1],'tooltipstring','Runs remotely or in background using default HPC settings in CONN');
                    uicontrol('style','pushbutton','string','OK','units','norm','position',[.1,.01,.38,.10],'callback','uiresume');
                    uicontrol('style','pushbutton','string','Cancel','units','norm','position',[.51,.01,.38,.10],'callback','delete(gcbf)');
                    uiwait(thfig);
                    if ~ishandle(thfig), return; end
                    options.fit=get(ht1,'value');
                    options.null=get(ht2,'value');
                    options.parallel=get(ht3,'value');
                    delete(thfig);
                    drawnow;
                end
                if options.parallel
                    cpm_lastjob = conn('submit',mfilename,option,varargin{:},'predictor',options.predictor,'outcome',options.outcome,'covariate',options.covariate,'fit',options.fit,'null',options.null,'parallel',false,'folder',options.folder);
                    varargout={cpm_lastjob};
                    return
                end                
                if options.fit
                    [model,Yfit]=conn_clusterregress(data_predictor,options.outcome,'covariates',options.covariate);
                    fprintf('Model %s: MSE = %f\n',filename_project,model.fit.MSE_std);
                else
                    [model]=conn_clusterregress(data_predictor,options.outcome,'covariates',options.covariate);
                    fprintf('Model %s: MSE = %f\n',filename_project,model.error.MSE_std);
                    Yfit = [];
                end
                save(fullfile(options.folder,conn_prepend('',filename_project,'.mat')), 'model', 'Yfit', 'options');
                switch(data_type)
                    case 'mtx', conn_mtx_write(fullfile(options.folder,conn_prepend('',filename_project,'.mtx.nii')),reshape(model.parameters.B_std,sqrt(numel(model.parameters.B_std))*[1 1]),info{1:2});
                    case 'surf',conn_surf_write(fullfile(options.folder,conn_prepend('',filename_project,'.surf.nii')),model.parameters.B_std);
                    case 'vol', conn_vol_write(fullfile(options.folder,conn_prepend('',filename_project,'.nii')),reshape(model.parameters.B_std,info.dim(1:3)),info);
                end
                if options.null
                    rand('state',0);
                    clear ModelNull;
                    mse=[];
                    for n=1:options.nnull,
                        randidx=randperm(size(data_predictor,1));
                        ModelNull(n)=conn_clusterregress(data_predictor,options.outcome(randidx,:),'covariates',options.covariate);
                        mse(n)=ModelNull(n).error.MSE_std;
                        fprintf('Model NULL %d: MSE = %f\n',n,mse(n));
                    end
                    p = mean(mse<=model.error.MSE_std); 
                    fprintf('Model %s: MSE = %f ; p = %f\n',filename_project,model.error.MSE_std,p);
                    save(fullfile(options.folder,conn_prepend('',filename_project,'.null.mat')),'ModelNull', 'p');
                end

            case {'parallel.status','parallel_status'}
                if ~isempty(cpm_lastjob)
                    conn_jobmanager('statusjob', cpm_lastjob, [], true, true);
                end

            case {'display.scatter','display_scatter'}
                filename_project=varargin{2}; 
                load(conn_prepend('',filename_project,'.mat'), 'model','options');
                assert(isfield(model,'fit'),'Fit values not computed yet. Please run CPM using the syntax conn_module(''CPM'',...,''fit'', true) to compute fit values for the training set');
                figure;
                conn_menu_plotscatter( options.outcome, model.fit.Yfit, varargin{3:end});
                xlabel 'Outcome measure';
                ylabel 'Predicted values';

            case {'display.forward','display_forward'}
                filename_project=varargin{2}; 
                load(conn_prepend('',filename_project,'.mat'), 'options');
                if (iscell(options.predictor)||ischar(options.predictor))&&~all(conn_existfile(options.predictor)), 
                    fprintf('warning: unable to find data file %s\n',options.predictor);
                    options.predictor=regexprep(options.predictor,'.*[\/\\]','');
                    assert(all(conn_existfile(options.predictor)),'unable to find data file %s',options.predictor);
                end                
                conn_module('glm','data',options.predictor,'design_matrix',[ones(size(options.outcome,1),1),options.outcome],'contrast_between',[zeros(size(options.outcome,2),1), eye(size(options.outcome,2))],'folder',filename_project);

            case 'display'
                filename_project=varargin{2}; 
                dispoptions=varargin(3:end);
                if conn_existfile(conn_prepend('',filename_project,'.mtx.nii'))
                    if isempty(dispoptions), dispoptions={'top100'};end
                    conn_mtx_braindisplay(conn_prepend('',filename_project,'.mtx.nii'), dispoptions{:});
                elseif conn_existfile(conn_prepend('',filename_project,'.surf.nii'))
                    conn_mesh_display(conn_prepend('',filename_project,'.surf.nii'),dispoptions{:});
                elseif conn_existfile(conn_prepend('',filename_project,'.vol.nii'))
                    conn_mesh_display(conn_prepend('',filename_project,'.nii'),dispoptions{:});
                end
        end

    case 'glm'        
        % loads .cfg files
        options=struct;
        if nargin==1 
            % syntax: conn_module GLM
            conn_module_glminternal;
            return            
        elseif nargin==2||(nargin>2&&rem(numel(varargin),2)>0&&(ischar(varargin{1})||isstruct(varargin{1})||isempty(varargin{1})))
            % syntax: conn_module GLM optionsfile.cfg
            % syntax: conn_module GLM optionsfile.cfg fieldname1 fieldvalue1 ...
            n=1;
            if ~isempty(varargin{n}), 
                filename=varargin{n};
                if ischar(filename)
                    if ~conn_existfile(filename),
                        if conn_existfile(fullfile(fileparts(which(mfilename)),filename)), filename=fullfile(fileparts(which(mfilename)),filename);
                        else error('file %s not found',filename);
                        end
                    end
                    fprintf('loading file %s\n',filename);
                    [nill,nill,fext]=fileparts(filename);
                    if strcmp(fext,'.json')
                        options=conn_jsonread(filename);
                    else
                        options=conn_loadcfgfile(filename,options);
                    end
                else
                    options=conn_loadcfgfile(filename,options);
                end
            end
        elseif nargin>1&&~isempty(varargin{1})&&isnumeric(varargin{1}) 
            % syntax: conn_module GLM X,Y,... (for back-compatibility)
            if nargout>0, conn_module_glminternal(varargin{:});
            else [varargout{1:nargout}]=conn_module_glminternal(varargin{:});
            end
            return            
        else 
            % syntax: conn_module GLM fieldname1 fieldvalue1 ...
            n=0;
        end
        for n=n+1:2:numel(varargin)-1
            fieldname=regexp(varargin{n},'\.','split');
            fieldvalue=varargin{n+1};
            options=setfield(options,fieldname{:},fieldvalue);
        end
        fprintf('%s options:\n',mfilename);
        disp(options);
        options0=options;
        options0.arguments=varargin;
        
        % interprets info
        
        if isfield(options,'design')&&~isfield(options,'design_matrix'), 
            options.design_matrix=options.design'; 
            options=rmfield(options,'design');
        end
        if isfield(options,'design_file')&&~isfield(options,'design_matrix'), 
            fprintf('loading design information from %s\n',options.design_file);
            x=conn_loadtextfile(options.design_file);
            if isstruct(x)
                enames=fieldnames(x);
                X=[]; teffectnames={};
                for n=1:numel(enames),
                    t=x.(enames{n});
                    X=[X t];
                    if size(t,2)>1, teffectnames=[teffectnames arrayfun(@(m)[enames{n},num2str(m)],1:size(t,2),'uni',0)];
                    else teffectnames=[teffectnames enames(n)];
                    end
                end
                if ~isfield(options,'design_labels')||isempty(options.design_labels),design_labels=teffectnames; end
            else X=x;
            end
            options=rmfield(options,'design_file');
        end
        assert(isfield(options,'design_matrix'),'missing #design or #design_matrix information');
        X = options.design_matrix;
        options=rmfield(options,'design_matrix');
        assert(isfield(options,'data'),'missing #data information');
        Y = options.data;
        options=rmfield(options,'data');
        if ischar(Y), Y=cellstr(Y); end
        if numel(Y)==1&&any(Y{1}=='*'), Y=conn_dir(Y{1},'-ls'); end
        if numel(Y)==1, try, Y=conn_expandframe(Y); end; end
        assert(~rem(numel(Y),size(X,1)),'mismatch number of subjects in #design (%d subjects) and total number of entries in #data (%d)',size(X,1),numel(Y));
        Y = reshape(Y, size(X,1),[]);
        Yisdir=conn_existfile(Y,2);
        if any(Yisdir), Y(Yisdir)=cellfun(@(x)fullfile(x,'SPM.mat'),Y(Yisdir),'uni',0); end
        if isfield(options,'contrast_between'),
            contrast_between=options.contrast_between;
            options=rmfield(options,'contrast_between');
        else contrast_between=[];
        end
        if isfield(options,'contrast_within'),
            contrast_within=options.contrast_within;
            options=rmfield(options,'contrast_within');
        else contrast_within=[];
        end
        if isfield(options,'folder'),
            spmfolder=char(options.folder);
            options=rmfield(options,'folder');
        else spmfolder=pwd;
        end
        if isfield(options,'explicitmask')&&~isfield(options,'mask')
            options.mask=options.explicitmask;
            options=rmfield(options,'explicitmask');
        end
        if isfield(options,'mask'),
            maskfile=char(options.mask);
            options=rmfield(options,'mask');
        else maskfile='';
        end
        if isfield(options,'analysistype'),
            analysistype=char(options.analysistype);
            options=rmfield(options,'analysistype');
        else analysistype=[]; % 1:all; 2:param only; 3:nonparam only
        end
        if isfield(options,'design_labels')&&~isempty(options.design_labels),
            effectnames=options.design_labels;
            options=rmfield(options,'design_labels');
        else effectnames={};
        end
        if isfield(options,'data_labels')&&~isempty(options.data_labels),
            contrastnames=options.data_labels;
            options=rmfield(options,'data_labels');
        elseif isfield(options,'contrast_names')&&~isempty(options.contrast_names)
            contrastnames=options.contrast_names;
        else contrastnames={};
        end
        if isfield(options,'subjectids')&&~isempty(options.subjectids),
            subjectids=options.subjectids;
            options=rmfield(options,'subjectids');
        else subjectids={};
        end
        [nill,nill,ext]=cellfun(@spm_fileparts,Y,'uni',0);
        ext=unique(ext);
        if isfield(options,'conditions')&&~isfield(options,'contrast_names'), options.contrast_names=options.conditions; end % back-compatibility
        
        if any(strcmp(ext,'.mat'))
            assert(numel(ext)==1,'mixed format files not supported (%s)',sprintf('%s ',ext{:}));
            assert(isfield(options,'contrast_names'),'missing #contrast_names information');
            cnames=options.contrast_names;
            newY=cell(size(Y,1),numel(cnames));
            for n1=1:size(Y,1),
                missing=1:numel(cnames);
                for n2=1:size(Y,2)
                    data=conn_loadmatfile(Y{n1,n2},'SPM');
                    if isfield(data.SPM.xCon,'name')
                        names={data.SPM.xCon.name};
                        [ok,idx]=ismember(cnames(missing),names); % find in list of contrast names
                        if numel(cnames)==size(Y,2), ok(missing~=n2)=false; end % if number of SPM.mat files matches number of contrast names, assume the two are paired
                        if any(ok),
                            for n3=reshape(find(ok),1,[])
                                filename=data.SPM.xCon(idx(n3)).Vcon.fname;
                                if isempty(fileparts(filename)), filename=fullfile(data.SPM.swd,filename); end
                                newY{n1,missing(n3)}=filename;
                            end
                            missing(ok)=[];
                        end
                        if isempty(missing), break; end
                    end
                end
                assert(isempty(missing)|any(isnan(X(n1,:)))|all(X(n1,:)==0),'unable to find contrast_names %s in %s',sprintf('%s ',cnames{missing}),sprintf('%s ',Y{n1,:}));
            end
            Y=newY;
            options=rmfield(options,'contrast_names');
            assert(isempty(contrast_within)||size(Y,2)==size(contrast_within,2),'mismatch number of columns in #contrast_within (%d) and number of measures in #contrast_names (%d)',size(contrast_within,2),size(Y,2));
        else
            assert(isempty(contrast_within)||size(Y,2)==size(contrast_within,2),'mismatch number of columns in #contrast_within (%d) and number of measures in #data (%d)',size(contrast_within,2),size(Y,2));
        end
        if isempty(analysistype)
            if any(cellfun('length',regexp(Y,'\.mtx\.nii(,\d+)?$|\.surf\.nii(,\d+)?$'))), analysistype=3;
            else analysistype=1;
            end
        end
        tnames=fieldnames(options);
        %asssert(isempty(tnames), fprintf('unable to interpret the following information fields: %s\n',sprintf('%s ',tnames{:})));
        assert(isempty(contrast_between)||size(X,2)==size(contrast_between,2),'mismatch number of columns in #contrast_between (%d) and number of effects in #design (%d)',size(contrast_between,2),size(X,2));
        
        % runs analysis
        if ~nargout, 
            conn_module_glminternal(X,Y,contrast_between,contrast_within,spmfolder,effectnames,contrastnames,analysistype,maskfile,subjectids);
        else 
            SPM=conn_module_glminternal(X,Y,contrast_between,contrast_within,spmfolder,effectnames,contrastnames,analysistype,maskfile,subjectids);
            varargout={spmfolder,SPM};
        end

    case 'default'
        if isempty(varargin), varargout={defaults}; 
        elseif isfield(defaults,varargin{1})
            if numel(varargin)>1, defaults.(varargin{1})=varargin{2};
            else varargout={defaults.(varargin{1})};
            end
        else error('unrecognized default %s',varargin{1});
        end
        
    otherwise
        if ~isempty(which(sprintf('conn_module_%s',option))), % conn_module_* functions
            fh=eval(sprintf('@conn_module_%s',option));
            if nargout, [varargout{1:nargout}]=feval(fh,varargin{:});
            else feval(fh,varargin{:});
            end
        elseif ismember(regexprep(lower(option),'_.*$',''),modules) % conn/modules/[root]/[root_*] functions
            if ismember(lower(option),modules), fh=eval(sprintf('@%s',lower(option)));
            else fh=eval(sprintf('@%s',option));
            end
            if nargout, [varargout{1:nargout}]=feval(fh,varargin{:});
            else feval(fh,varargin{:});
            end
        else
            disp(sprintf('unrecognized conn_module option %s or function conn_module_%s or module %s',option,option,option));
        end
end
end
