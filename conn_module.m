function varargout=conn_module(option,varargin)
% CONN_MODULE provides access to independent CONN modules 
%
% conn_module(module_name, ...) runs individual CONN's module "module_name" on user-defined data
% Current module names:
%
%    PREPROCESSING : runs CONN preprocessing pipeline on user-defined data
%
%       basic syntax: conn_module preprocessing
%       advanced syntax: conn_module('preprocessing', fieldname1, fieldvalue1, fieldname2, fieldvalue2, ...)
%
%                 Input data is specified with field name/value pairs as defined in batch.Setup documentation
%                 Preprocessing options are specified with field name/value pairs as defined in batch.Setup.preprocessing documentation
%                    functionals       : list of functional data files { { Sub1Ses1, Sub1Ses2, ...}, {Sub2Ses1, Sub2Ses2, ...}, ...}
%                    structurals       : list of structural data files { Sub1, Sub2, ...}
%                    steps             : list of preprocessing steps (tpye "conn_module preprocessing steps" for a list of valid preprocessing step names)
%
%                 See "doc conn_batch" for a complete list of these options
%                 See Nieto-Castanon, 2020 for details about these preprocessing steps and pipelines (www.conn-toolbox.org/fmri-methods)
%
%       e.g. conn_module('preprocessing',...
%             'steps','default_mni',...
%             'functionals',{'./func.nii'},...
%             'structurals',{'./anat.nii'},...
%             'RT',2,...
%             'sliceorder','interleaved (Siemens)');
%            runs default MNI-space preprocessing pipeline on the specified functional/structural data
%
%       alternative syntax: conn_module('preprocessing',optionsfile) 
%          input data and preprocessing options defined in .cfg (see conn_loadcfgfile/conn_savecfgfile) or .json (see spm_jsonread/spm_jsonwrite) structure text file 
%
%       alternative syntax: conn_module preprocessing steps
%          returns the list of valid preprocessing steps
%
%    GLM : runs CONN second-level analyses on user-defined data
%
%       basic syntax: conn_module glm
%       advanced syntax: conn_module('glm', fieldname1, fieldvalue1, fieldname2, fieldvalue2, ...)
%               with the following field name/value pairs
%                  data            : list of nifti files entered into second-level analysis (Nsubjects x Nmeasures) defining one or multiple outcome / dependent measures 
%                                       note: when entering multiple files per subject (e.g. repeated measures) enter first all files (one per subject) for measure#1, followed by all files for measure#2, etc.
%                                       note: nifti files may contani 3d volume-level data, fsaverage surface-level data, or ROI-to-ROI data (see conn_surf_write and conn_mtx_write to create surface/matrix nifti files)
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
%         eg: conn_module('glm', ...
%            'design_matrix',[1; 1; 1; 1] ,...
%            'data',{'subject1.img'; 'subject2.img'; 'subject3.img'; 'subject4.img'} );
%             performs a one-sample t-test and stores the analysis results in the current folder
%
%         eg: conn_module('glm', ...
%            'design_matrix',[1 0; 1 0; 0 1; 0 1; 0 1],...
%            'data', {'subject1_group1.img'; 'subject2_group1.img'; 'subject1_group2.img'; 'subject2_group2.img'; 'subject3_group2.img'},...
%            'contrast_between',[1 -1]);
%             performs a two-sample t-test and stores the analysis results in the current folder
%
%         eg: conn_module('glm', ...
%            'design_matrix', [1; 1; 1; 1],...
%            'data', {'subject1_time1.img', subject1_time2.img'; 'subject2_time1.img', subject2_time2.img'; 'subject3_time1.img', subject3_time2.img'; 'subject4_time1.img', subject4_time2.img'},...
%            'contrast_beetween',1,...
%            'contrast_within',[1 -1]);
%             performs a paired t-test and stores the analysis results in the current folder
%
%       alternative syntax: conn_module('glm',optionsfile) 
%          input data and GLM options defined in .cfg (see conn_loadcfgfile/conn_savecfgfile) or .json (see spm_jsonread/spm_jsonwrite) structure text file 
%
%       alternative syntax: spmfolder=conn_module('glm',...) 
%          skips results display step (only computes second-level analysis, and returns folder where results are stored)
%          use conn_display(spmfolder) syntax to then launch the results explorer window on previously computed analyses
%
%       See also CONN_DISPLAY for displaying GLM results
%       See Nieto-Castanon, 2020 for details about General Linear Model analyses (www.conn-toolbox.org/fmri-methods)
%
%    Additional functionality: conn_module('get',...)
%          conn_module('get','structurals');             outputs current structural files (e.g. output of structural preprocessing steps)
%          conn_module('get','functionals' [,setlabel]); outputs current functional files (e.g. output of functional preprocessing steps)
%          conn_module('get','l1covariates' [,covname]); outputs first-level covariate files (e.g. other potential outputs of functional preprocessing)
%          conn_module('get','l2covariates' [,covname]); outputs second-level covariate values (e.g. other potential outputs of functional preprocessing)
%          conn_module('get','masks');                   outputs Grey Matter/White Matter/CSF files (e.g. other potential outputs of functional preprocessing)
%          conn_module('get','masks',roiname);           outputs roiname files (e.g. other potential outputs of functional preprocessing)
%    Additional functionality: conn_module('set',...)
%          conn_module('set','l1covariates',files,covname [,add]);
%          conn_module('set','l2covariates',values,covname [,covdescrip ,add]);
%          conn_module('set','masks',files [,roiname]);
%
%    Note: before using conn_module functionality with externally defined data it is recommended to close CONN's gui in order to avoid potentially loosing any unsaved changes
%    

persistent defaults;

if isempty(defaults), defaults=struct('mat_format','-v7.3'); end
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
                        for nses=1:nsessall,
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
        
    case 'preprocessing'
        options=struct;
        if isempty(varargin)
            answ=inputdlg('Enter number of subjects','conn_module preprocessing',1,{'1'});
            if isempty(answ), return; end
            options.nsubjects=str2num(answ{1});
            options.functionals=cell(options.nsubjects,1);
            nsessions=ones(1,options.nsubjects);
            for nsubject=1:options.nsubjects,
                temp=inputdlg(['Subject ',num2str(nsubject),': Enter number of runs/sessions'],'conn_module preprocessing',1,{num2str(nsessions(min(length(nsessions),nsubject)))});
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
            % syntax: conn_module('preprocessing',struct('data',...),...)
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
            % syntax: conn_module preprocessing optionsfile.cfg
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
                        options=spm_jsonread(filename);
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
            % syntax: conn_module preprocessing fieldname1 fieldvalue1 ...
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
        
    case 'glm'
        % loads .cfg files
        options=struct;
        if nargin==1 
            % syntax: conn_module glm
            conn_module_glminternal;
            return            
        elseif nargin==2||(nargin>2&&rem(numel(varargin),2)>0&&(ischar(varargin{1})||isstruct(varargin{1})||isempty(varargin{1})))
            % syntax: conn_module glm optionsfile.cfg
            % syntax: conn_module glm optionsfile.cfg fieldname1 fieldvalue1 ...
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
                        options=spm_jsonread(filename);
                    else
                        options=conn_loadcfgfile(filename,options);
                    end
                else
                    options=conn_loadcfgfile(filename,options);
                end
            end
        elseif nargin>1&&~isempty(varargin{1})&&isnumeric(varargin{1}) 
            % syntax: conn_module glm X,Y,... (for back-compatibility)
            if nargout>0, conn_module_glminternal(varargin{:});
            else [varargout{1:nargout}]=conn_module_glminternal(varargin{:});
            end
            return            
        else %if nargin>2&&ischar(varargin{1})&&(ismember(varargin{1},{'data','design','design_matrix','contrast_between','contrast_within','contrast_names','mask','data_labels','design_labels','analysistype','folder'})||isempty(dir(varargin{1}))) 
            % syntax: conn_module glm fieldname1 fieldvalue1 ...
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
                    data=load(Y{n1,n2},'SPM');
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
        if ~isempty(which(sprintf('conn_module_%s',option))),
            fh=eval(sprintf('@conn_module_%s',option));
            if nargout, [varargout{1:nargout}]=feval(fh,varargin{:});
            else feval(fh,varargin{:});
            end
        else
            disp(sprintf('unrecognized option %s or conn_module_%s function',option,option));
        end
end
end
