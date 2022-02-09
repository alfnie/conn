function varargout=el(option,varargin)
% EVLAB tools
%
% INITIALIZATION SYNTAX:
%
%   conn module EL init                  : initializes EVLAB tools
%   el root.subjects <subjects_folder>   : defines root directory where subject folders are stored
%                                    Altnernatively, this is defined by the symbolic link .../conn/modules/el/root.subjects if it exists
%                                    (syntax "el root.subjects <subjects_folder> all" may be used to create this symbolic link)
%   el root.tasks <tasks_folder>         : defines root directory where experimental design files are stored
%                                    Altnernatively, this is defined by the symbolic link .../conn/modules/el/root.tasks if it exists
%                                    (syntax "el root.tasks <tasks_folder> all" may be used to create this symbolic link)
%   el root.pipelines <pipelines_folder> : defines root directory where preprocessing and model pipeline files are stored
%                                    By default this points to .../conn/modules/el, where several default pipelines already exist
%   el remote on                         : work remotely (default 'off') (0/'off': when working with data stored locally on your computer; 1/'on': when working with data stored in a remote server -to enable this functionality on a remote server run on the remote server the command "conn remotely setup")
%
% PREPROCESSING SYNTAX:
%
%   el('preprocessing',subjectID [, pipelineID]) imports (e.g. from DICOM) and preprocesses the functional&structural data for subject 'subjectID'
%      subjectID                : subject folder name (e.g. '408_FED_20160617a_3T2')
%      pipelineID               : (optional) 
%                                   .cfg file defining preprocessing steps (e.g. .../modules/el/pipeline_preproc_DefaultMNI.cfg)
%                                   Alternatively, a pipelineID is just a shortcut to a .cfg file located in conn/modules/el/ and named "pipeline_preproc_<pipelineID>.cfg"; e.g. 'DefaultMNI' is a shortcut to the .cfg file named ..../modules/el/pipeline_preproc_DefaultMNI.cfg)
%                                   Alternatively, if unspecified or left empty:
%                                       'DefaultMNI_PlusStructural' : if the data contains structural/anatomical image(s)
%                                       'DefaultMNI'                : if the data does not contain any structural/anatomical image
%      data.cfg                 : file located in the subject folder, and defining the location and source of functional and anatomical data (e.g. example_datafiles_DicomFunctionals.cfg)
%                                   note: the original DICOM or NIFTI data may or may not be located within the subject folder; when preprocessing EVLAB will create a new folder within the subject 
%                                   folder and with name equal to the pipelineID (or name of preprocessing_pipeline .cfg file), which will contain the imported NIFTI functional/anatomical data
%                                   as well as all the associated files resulting from running all of the preprocessing steps
%
%      e.g.
%           >> el preprocessing 408_FED_20160617a_3T2
%
%
%   el('preprocessing.append',subjectID, pipelineID, pipelineID_append) runs additional preprocessing steps to the (already previously preprocessed) functional&anatomical data
%
%      e.g.
%           >> el preprocessing.append 408_FED_20160617a_3T2 DefaultMNI_PlusStructural OnlyDenoise
%
%
%   submitID = el('submit','preprocessing',subjectID [, preprocessing_pipeline]) runs import&preprocessing steps on remote node
%   submitID = el('submit','preprocessing.append',subjectID, preprocessing_pipeline, preprocessing_pipeline_append) runs additional preprocessing steps on remote node
%
%      e.g.  
%           >> el submit preprocessing 408_FED_20160617a_3T2
%      e.g.  
%           >> el submit preprocessing.append 408_FED_20160617a_3T2 DefaultMNI_PlusStructural OnlyDenoise
%
%   el('preprocessing.qa',subjectID, pipelineID) creates QA plots on already-preprocessed dataset
%   el('preprocessing.qa.plot',subjectID, pipelineID) displays already-created QA plots
%
%
% MODEL SYNTAX:
%
%   el('model',subjectID, pipelineID, experimentID [, modelOPTIONS]) runs first-level GLM model estimation step
%      subjectID                : subject folder name (e.g. '408_FED_20160617a_3T2')
%      pipelineID               : preprocessing pipeline that will be used as source of data in this analysis
%      experimentID             : .cat file defining the experimental design associated with each functional run
%                                   Alternatively, a experimentID is just a shortcut to a .cat file located in the subject directory and named "*_<experimentID>.cat"
%      modelOPTIONS             : (optional) 
%                                   .cfg file defining additional first-level model estimation options
%                                   Alternatively, a modelOPTIONS is just a shortcut to a .cfg file located in conn/modules/el/ and named "pipeline_model_<modelOPTIONS>.cfg"; e.g. 'Default' is a shortcut to the .cfg file named ..../modules/el/pipeline_model_Default.cfg)
%                                   Alternatively, if unspecified or left empty, it takes the value 'Default' (conn/modules/el/pipeline_model_Default.cfg)
%      e.g.
%           >> el model 408_FED_20160617a_3T2 DefaultMNI_PlusStructural langlocSN
%
%   submitID = el('submit','model',...) runs first-level GLM model estimation step on remote node
%
%      e.g.  
%           >> el submit preprocessing 408_FED_20160617a_3T2
%
%   el('model.plot',subjectID, pipelineID, experimentID) displays first-level effect-size estimates
%   el('model.stats',subjectID, pipelineID, experimentID) displays first-level statistics
%   el('model.qa',subjectID, pipelineID, experimentID) creates QA plots on first-level GLM analyses
%   el('model.qa.plot',subjectID, pipelineID) displays already-created QA plots
%
%
% CONFIGURATION OPTIONS:
%
%   foldername = el('default','folder_dicoms' [,foldername])
%      foldername               : subdirectory name within folder_subjects where DICOM files may be found. Default '*dicoms'
%   list = el('default','dicom_disregard_functional' [,list])
%      list                     : list of keywords defining SeriesDescription values of DICOM sessions to be excluded from the list of valid functional runs
%   list = el('default','dicom_isstructural'[,list])
%      list                     : list of keywords defining SeriesDescription values of DICOM sessions to be interpreted as anatomical/structural scans
% 
%

persistent defaults;
if isempty(defaults), 
    defaults=struct(...
        'folder_subjects',fullfile(fileparts(which(mfilename)),'root.subjects') ,...
        'folder_tasks',fullfile(fileparts(which(mfilename)),'root.tasks') ,...
        'folder_pipelines',fileparts(which(mfilename)) ,...
        'folder_dicoms','*dicoms' ,...
        'dicom_isstructural',{{'^T1_MPRAGE_1iso'}} ,...
        'dicom_disregard_functional',{{'^localizer','^AAScout','^AAHScout','^MoCoSeries','^T1_MPRAGE_1iso','^DIFFUSION_HighRes'}},...
        'create_model_cfg_files',false,... % 1/0 specifying whether "el model *" commands create an associated .cfg file (for backwards compatibility) [false]
        'isremote', 0); 
end

conn_module('evlab17','init','silent');
fileout=[];
varargout=cell(1,nargout);

if defaults.isremote&&~(~isempty(regexp(lower(char(option)),'plots?$'))||ismember(lower(char(option)),{'root.subjects','root.tasks','root.pipelines','remote','init','initforce','default','model.stats'})); % run these locally
    [hmsg,hstat]=conn_msgbox({'Process running remotely','Please wait...',' ',' '},[],[],true);
    if ~isempty(hmsg), [varargout{1:nargout}]=conn_server('run_withwaitbar',hstat,mfilename,option,varargin{:}); 
    else [varargout{1:nargout}]=conn_server('run',mfilename,option,varargin{:}); 
    end
    if ~isempty(hmsg)&&ishandle(hmsg), delete(hmsg); end
    return
end

switch(lower(option))
    case 'root.subjects'
        if numel(varargin)>1&&isequal(varargin{2},'all')&&~defaults.isremote, conn_fileutils('linkdir',varargin{1},fullfile(fileparts(which(mfilename)),'root.subjects')); end
        if numel(varargin)>0, 
            defaults.folder_subjects=varargin{1};
            if defaults.isremote, conn_server('run',mfilename,'root.subjects',defaults.folder_subjects); end
        else varargout={defaults.folder_subjects};
        end
    case 'root.tasks'
        if numel(varargin)>1&&isequal(varargin{2},'all')&&~defaults.isremote, conn_fileutils('linkdir',varargin{1},fullfile(fileparts(which(mfilename)),'root.tasks')); end
        if numel(varargin)>0, 
            defaults.folder_tasks=varargin{1};
            if defaults.isremote, conn_server('run',mfilename,'root.tasks',defaults.folder_tasks); end
        else varargout={defaults.folder_tasks};
        end
        
    case 'root.pipelines'
        if numel(varargin)>0, 
            defaults.folder_pipelines=varargin{1};
            if defaults.isremote, conn_server('run',mfilename,'root.pipelines',defaults.folder_pipelines); end
        else varargout={defaults.folder_pipelines};
        end
        
    case 'remote'
        if nargin>1&&~isempty(varargin{1}), 
            defaults.isremote=varargin{1}; 
            if ischar(defaults.isremote)
                switch(lower(defaults.isremote))
                    case 'on',  defaults.isremote=1;
                    case 'off', defaults.isremote=0;
                    otherwise, defaults.isremote=str2num(defaults.isremote); 
                end
            end
            if defaults.isremote&&~conn_server('isconnected')
                fprintf('Starting new remote connection to server\n');
                conn remotely on;
            end
            if defaults.isremote, fprintf('working with remote projects now\n');
            else fprintf('working with local projects now\n');
            end
            if defaults.isremote&&conn_server('isconnected')
                conn_server('run','conn_module','el','init');
                defaults=conn_server('run',mfilename,'default');
                defaults.isremote=true;
                fprintf('ROOT.SUBJECTS folder set to %s\n',defaults.folder_subjects);
                fprintf('ROOT.TASKS folder set to %s\n',defaults.folder_tasks);
                fprintf('ROOT.PIPELINES folder set to %s\n',defaults.folder_pipelines);
            end
            if nargout>0, varargout={defaults.isremote}; end
        else
            if ~nargout,
                if defaults.isremote, fprintf('working with remote projects (remote=1)\n');
                else fprintf('working with local projects (remote=0)\n');
                end
            else varargout={defaults.isremote};
            end
        end
        
    case 'default'
        if isempty(varargin), varargout={defaults}; 
        elseif isfield(defaults,varargin{1})
            if numel(varargin)>1, 
                if defaults.isremote, conn_server('run',mfilename,option,varargin{:}); end
                defaults.(varargin{1})=varargin{2};
            else varargout={defaults.(varargin{1})};
            end
        else
            if numel(varargin)>1, conn_module('evlab17','default',varargin{:});
            else varargout={conn_module('evlab17','default',varargin{1})};
            end
        end
        
    case 'submit'
        if ~nargout, conn('submit',mfilename,varargin{:}); % e.g. el submit preprocessing 408_FED_20160617a_3T2
        else [varargout{1:nargout}]=conn('submit',mfilename,varargin{:});
        end
        
    case 'preprocessing'
        assert(numel(varargin)>=1,'incorrect usage >> el preprocessing subject_id [, pipeline_id]');
        % adapted from msieg preprocess_PL2017
        % el('preprocessing',subject_id [, preprocessing_pipeline])
        % e.g. el preprocessing 408_FED_20160617a_3T2
        
        pwd0=pwd;
        [subject,data_config_file,subject_path]=el_readsubject(varargin{1},defaults);
        assert(conn_existfile(data_config_file),'unable to find data configuration file %s',data_config_file);

        info=conn_loadcfgfile(data_config_file);
        if isfield(info,'structurals'), struct_run=info.structurals;
        else struct_run=[];
        end
        if numel(varargin)<2||isempty(varargin{2})
            if isempty(struct_run), preproc_config_file = fullfile(defaults.folder_pipelines,'pipeline_preproc_DefaultMNI.cfg');
            else preproc_config_file = fullfile(defaults.folder_pipelines,'pipeline_preproc_DefaultMNI_PlusStructural.cfg');
            end
            assert(conn_existfile(preproc_config_file),'unable to find preprocessing pipeline %s',preproc_config_file);
        else
            preproc_config_file = el_readpreprocconfig(varargin{2},defaults);
            assert(conn_existfile(preproc_config_file),'unable to find preprocessing pipeline %s',varargin{2});
        end
        
        %run preproc
        if strcmpi(option,'preprocessing.main')
            [ok,msg]=mkdir(fullfile(subject_path,'nii'));
            cd(fullfile(subject_path,'nii'));
            fileout=conn_module('evlab17','run_preproc',data_config_file,preproc_config_file,[],...
                'alt_functionals_path', fullfile(subject_path,'func'), ...
                'alt_structurals_path', fullfile(subject_path,'anat'), ...
                'alt_rois_path', fullfile(subject_path,'rois'),...
                varargin(3:end));
        else % copies data to pipeline-specific folder before preprocessing
            [nill,tname]=fileparts(preproc_config_file);
            tname=regexprep(tname,'^pipeline_preproc_','');
            dataset=fullfile(subject_path,[tname,'.mat']);
            [ok,msg]=mkdir(fileparts(dataset));
            cd(fileparts(dataset));
            fileout=conn_module('evlab17','run_preproc',data_config_file,[],...
                'dataset',dataset, ...
                'alt_functionals_path', fullfile(subject_path,'func'), ...
                'alt_structurals_path', fullfile(subject_path,'anat'), ...
                'alt_rois_path', fullfile(subject_path,'rois'));
            conn_module('evlab17','save',dataset);
            conn_module('evlab17','update'); % import
            conn_importvol2bids(2);
            conn save;
            fileout=conn_module('evlab17','run_preproc',preproc_config_file,[],...
                'dataset',dataset,...
                'alt_functionals_path', fullfile(subject_path,'func'), ...
                'alt_structurals_path', fullfile(subject_path,'anat'), ...
                'alt_rois_path', fullfile(subject_path,'rois'),...
                varargin(3:end));
        end
        cd(pwd0);
        varargout{1}=fileout;
    
    case 'preprocessing.append'
        assert(numel(varargin)>=3,'incorrect usage >> el preprocessing.append subject_id, original_pipeline_id, added_pipeline_id');
        % el('preprocessing.append',subject_id, preprocessing_pipeline_original, preprocessing_pipeline_append)
        
        pwd0=pwd;
        [subject,data_config_file,subject_path]=el_readsubject(varargin{1},defaults);
        dataset=el_readpipeline(varargin{2},subject,subject_path,defaults);
        assert(conn_existfile(dataset),'unable to find dataset %s',dataset);
        conn_module('evlab17','load',dataset);
        preproc_config_file = el_readpreprocconfig(varargin{3},defaults);
        assert(conn_existfile(preproc_config_file),'unable to find preprocessing pipeline %s',varargin{3});
        if isempty(preproc_config_file), opts={[]};
        else opts={preproc_config_file, []};
        end
        %run preproc
        cd(fileparts(dataset));
        fileout=conn_module('evlab17','run_preproc',opts{:},...
            'dataset',dataset, ...
            'alt_functionals_path', fullfile(subject_path,'func'), ...
            'alt_structurals_path', fullfile(subject_path,'anat'), ...
            'alt_rois_path', fullfile(subject_path,'rois'), ...
            varargin(4:end));
        cd(pwd0);
        varargout{1}=fileout; 
    
    case {'data.cfg.dicom','createdatacfg'} % subject id
        assert(numel(varargin)>=1,'incorrect usage >> el createdatacfg subject_id');
        pwd0=pwd;
        [subject,data_config_file,subject_path]=el_readsubject(varargin{1},defaults);
        assert(~conn_existfile(data_config_file),'data configuration file %s already exist. Please delete this file before proceeding',data_config_file);
        
        subject_path_dicoms = fullfile(subject_path,defaults.folder_dicoms);
        if any(subject_path_dicoms=='*'), subject_path_dicoms = conn_dir(subject_path_dicoms ,'-dir','-R','-cell','-sort'); end
        assert(~isempty(subject_path_dicoms),'unable to find dicom folder %s\n',fullfile(subject_path,defaults.folder_dicoms));
        if iscell(subject_path_dicoms), subject_path_dicoms=subject_path_dicoms{1}; end
        data_config_file='';
        func_runs=[];
        struct_run=[];
        
        Series=conn_dcmdir(fullfile(subject_path_dicoms,'*-1.dcm'),false);
        idx=find(cellfun('length',{Series.SeriesDescription})>0&cellfun('length',{Series.SeriesNumber})>0);
        SeriesDescription={Series(idx).SeriesDescription};
        SeriesNumber=[Series(idx).SeriesNumber];
        struct_run = SeriesNumber(find(cellfun(@(x)~isempty(regexp(char(x),strjoin(defaults.dicom_isstructural,'|'))),SeriesDescription)>0,1,'first')); % keep first structural
        func_runs = setdiff(SeriesNumber(find(cellfun(@(x)isempty(regexp(char(x),strjoin(defaults.dicom_disregard_functional,'|'))),SeriesDescription)>0)),struct_run); % keep all functionals
        assert(~isempty(func_runs),'unable to find any functional runs in %s\n',subject_path_dicoms);
        try,
            fid=1; %fopen(fullfile(subject_path_dicoms,'runs.csv'),'wt');
            fprintf(fid,'SeriesNumber,SeriesDescription\n');
            for n=1:numel(SeriesNumber), fprintf(fid,'%s,%s\n',num2str(SeriesNumber(n)),SeriesDescription{n}); end
            fclose(fid);
            %fprintf('DICOM series information stored in %s\n',fullfile(subject_path_dicoms,'runs.csv'))
        end
        
        all_runs=[struct_run func_runs];
        fid=fopen(fullfile(subject_path,'data.cfg'),'wt');
        fprintf(fid,'\n#dicoms\n');
        for n=1:numel(all_runs), fprintf(fid,'%s\n',fullfile(subject_path_dicoms,['*-' num2str(all_runs(n)) '-1.dcm'])); end
        fprintf(fid,'\n#functionals\n');
        fprintf(fid,'%s\n',num2str(func_runs));
        if ~isempty(struct_run)
            fprintf(fid,'\n#structurals\n');
            fprintf(fid,'%s\n',num2str(struct_run));
        end
        fprintf(fid,'\n#RT nan\n'); % note: will read RT from .json files
        fclose(fid);
        data_config_file=fullfile(subject_path,'data.cfg');
        fprintf('Created data configuration file %s\n',data_config_file);
        %             else % run numbers
        %                 if ischar(varargin{2})&&conn_existfile(varargin{2}), data_config_file=varargin{2};
        %                 elseif isstruct(varargin{2}), data_config_file=varargin{2};
        %                 elseif ischar(varargin{2})&&~isempty(str2num(varargin{2})), func_runs=str2num(varargin{2});
        %                 else error('unrecognized data descriptor argument');
        %                 end
        %             end
        varargout{1}=data_config_file;
            
    case {'data.cfg.nifti','data.cfg'} % subject id
        assert(numel(varargin)>=1,'incorrect usage >> el data.cfg subject_id');
        pwd0=pwd;
        [subject,data_config_file,subject_path]=el_readsubject(varargin{1},defaults);
        assert(~conn_existfile(data_config_file),'data configuration file %s already exist. Please delete this file before proceeding',data_config_file);
        
        subject_path_func = fullfile(subject_path,'func');
        assert(conn_existfile(subject_path_func,2),'unable to find func folder %s\n',subject_path_func);
        subject_path_anat = fullfile(subject_path,'anat');
        assert(conn_existfile(subject_path_anat,2),'unable to find anat folder %s\n',subject_path_anat);
        subject_path_fmap = fullfile(subject_path,'fmap');
        
        func=conn_dir(fullfile(subject_path_func,'*.nii'),'-ls');
        anat=conn_dir(fullfile(subject_path_anat,'*.nii'),'-ls');
        if conn_existfile(subject_path_fmap,2),fmap=conn_dir(fullfile(subject_path_fmap,'*.nii'),'-ls'); else fmap={}; end
        assert(~isempty(func),'unable to find any functional files in %s\n',subject_path_func);

        fid=fopen(fullfile(subject_path,'data.cfg'),'wt');
        fprintf(fid,'\n#functionals\n');
        for n=1:numel(func), fprintf(fid,'%s\n',func{n});end
        if ~isempty(anat)
            fprintf(fid,'\n#structurals\n');
            for n=1:numel(anat), fprintf(fid,'%s\n',anat{n});end
        end
        if ~isempty(fmap)
            fprintf(fid,'\n#fmap_functionals\n');
            for n=1:numel(fmap), fprintf(fid,'%s\n',fmap{n});end
        end
        fprintf(fid,'\n#RT nan\n'); % note: will read RT from .json files
        fclose(fid);
        data_config_file=fullfile(subject_path,'data.cfg');
        fprintf('Created data configuration file %s\n',data_config_file);
        varargout{1}=data_config_file;
            
    case {'preprocessing.qa','model.qa','qa.preprocessing','qa.model'}
        assert(numel(varargin)>=1,'incorrect usage >> el preprocessing.qa subject_id [,pipeline_id ,model_id]');
        [subject,data_config_file,subject_path]=el_readsubject(varargin{1},defaults);
        if numel(varargin)<2, pipeline=''; else pipeline=varargin{2}; end
        dataset=el_readpipeline(pipeline,subject,subject_path,defaults);
        assert(conn_existfile(dataset),'unable to find dataset %s',dataset);
        if ~isempty(regexp(lower(option),'model')), 
            if numel(varargin)<3, expt=''; else expt=varargin{3}; end
            expt=el_readexpt(expt,dataset);
            opts={'qa_plist',['firstlevel_',expt]};
        else opts={'qa_plist','preprocessing'};
        end
        conn_module('evlab17','run_qa','dataset',dataset,opts{:},varargin{4:end});
        
    case {'qa.plot','preprocessing.qa.plot','model.qa.plot'}
        assert(numel(varargin)>=1,'incorrect usage >> el qa.plot subject_id [,pipeline_id]');
        [subject,data_config_file,subject_path]=el_readsubject(varargin{1},defaults);
        if numel(varargin)<2, pipeline=''; else pipeline=varargin{2}; end
        dataset=el_readpipeline(pipeline,subject,subject_path,defaults);
        assert(conn_existfile(dataset),'unable to find dataset %s',dataset);
        conn_module('evlab17','qaplots',dataset,varargin{4:end});
        
    case {'model','firstlevel'}
        assert(numel(varargin)>=2,'incorrect usage >> el model subject_id, pipeline_id, model_id [, modeloptions]');
        pwd0=pwd;
        [subject,data_config_file,subject_path]=el_readsubject(varargin{1},defaults);
        if numel(varargin)<2, pipeline=''; else pipeline=varargin{2}; end
        dataset=el_readpipeline(pipeline,subject,subject_path,defaults);
        if numel(varargin)<3, expt=''; else expt=varargin{3}; end
        expt=el_readexpt(expt,dataset);
        if numel(varargin)<4, model_config_file=''; else model_config_file=varargin{4}; end
        try, if isempty(model_config_file)&&conn_existfile(fullfile(subject_path,[expt,'.cfg'])), model_config_file=fullfile(subject_path,[expt,'.cfg']); end; end
        model_config_file=el_readmodelconfig(model_config_file,defaults);
        assert(conn_existfile(model_config_file),'unable to find model estimation options %s',model_config_file);
        
        % adapted from msieg firstlevel_PL2017
        contrasts_file=fullfile(el_readtaskfolder(defaults),[expt,'.con']);
        if ~conn_existfile(contrasts_file), contrasts_file=fullfile(el_readtaskfolder(defaults),[expt,'.txt']); end
        if ~conn_existfile(contrasts_file), contrasts_file=fullfile(el_readtaskfolder(defaults),[expt,'.csv']); end
        if ~conn_existfile(contrasts_file), contrasts_file=fullfile(el_readtaskfolder(defaults),[expt,'.tsv']); end
        if conn_existfile(contrasts_file), % contrast definitions
            fprintf('contrast file %s found\n',contrasts_file); 
            str=conn_fileutils('fileread',contrasts_file);
            str=reshape(regexp(str,'\n','split'),1,[]);
            cons=reshape(str(cellfun('length',str)>0),1,[]);
        else % back-compatibility
            all_contrasts_files=fullfile(el_readtaskfolder(defaults),'contrasts_by_expt.txt'); % single-file, all expt contrasts
            assert(conn_existfile(all_contrasts_files), 'no contrast-file found\n');
            fprintf('contrast file %s found\n',all_contrasts_files); 
            %if ~conn_existfile(all_contrasts_files), all_contrasts_files=fullfile(defaults.folder_subjects,'..','ANALYSIS','contrasts_by_expt.txt'); end
            str=conn_fileutils('fileread',all_contrasts_files);
            str=reshape(regexp(str,'\n','split'),1,[]);
            emptyspaces=cellfun('length',str)==0;
            idx=find([true emptyspaces(1:end-1) true]&[~emptyspaces true]);
            con=strmatch(expt,str(idx(1:end-1)),'exact');
            assert(numel(con)==1,'found %d matches to %s in %s',numel(con),expt,all_contrasts_files);
            cons=str(idx(con)+1:idx(con+1)-3);
        end
        cat_file=conn_dir(fullfile(subject_path,[expt,'.cat']),'-ls');
        if isempty(cat_file), cat_file=conn_dir(fullfile(subject_path,['*_',expt,'.cat']),'-ls'); end
        assert(numel(cat_file)==1,'%d %s files found',numel(cat_file),fullfile(subject_path,['*_',expt,'.cat']));
        cat_info=conn_loadcfgfile(char(cat_file),struct('path',el_readtaskfolder(defaults)));
        DOSAVECFG=defaults.create_model_cfg_files;
        opts=struct('dataset',char(dataset),'design',char(cat_file),'model_folder','root','model_name',expt,'contrasts',{cons});
        if DOSAVECFG
            model_definition_file=conn_prepend('',char(cat_file),'.cfg');
            conn_savecfgfile(model_definition_file,opts);            
            conn_module('evlab17','run_model',model_definition_file,model_config_file,varargin{5:end});
        else
            conn_module('evlab17','run_model',opts,model_config_file,varargin{5:end});
        end
        cd(pwd0);

    case {'model.plot', 'model.stats'}
        assert(numel(varargin)>=3,'incorrect usage >> el model.plot subject_id, pipeline_id, and model_id')
        [subject,data_config_file,subject_path]=el_readsubject(varargin{1},defaults);
        if numel(varargin)<2, pipeline=''; else pipeline=varargin{2}; end
        dataset=el_readpipeline(pipeline,subject,subject_path,defaults);
        if numel(varargin)<3, expt=''; else expt=varargin{3}; end
        expt=el_readexpt(expt,dataset);
        opts={dataset,expt,...
            varargin{4:end}};
        assert(conn_existfile(dataset),'file %s not found',dataset);
        if strcmpi(option,'model.stats'), conn_module('evlab17','modelplots',opts{:},'stats');
        else conn_module('evlab17','modelplots',opts{:});
        end
        
    otherwise
        [varargout{1:nargout}]=conn_module('evlab17',option,varargin{:});
end

end

function [subject,data_config_file,subject_path]=el_readsubject(subject,defaults)
subject=char(subject); % subject id
if ~isempty(regexp(subject,'\.cfg')) % old format
    data_config_file=subject;
    subject_path=filepath(subject);
else
    subject_path=fullfile(defaults.folder_subjects,subject);
    if defaults.isremote, subject_path=fullfile('/CONNSERVER',subject_path); end
    assert(conn_existfile(subject_path,2),'unable to find directory %s',subject_path);
    data_config_file=fullfile(subject_path,'data.cfg');
end
end

function preproc_config_file = el_readpreprocconfig(preproc_config_file,defaults);
preproc_config_file=conn_prepend('',char(preproc_config_file),'.cfg');
if isempty(fileparts(preproc_config_file)),
    preproc_config_file=fullfile(defaults.folder_pipelines,preproc_config_file);
    if ~conn_existfile(preproc_config_file)&&conn_existfile(conn_prepend('pipeline_preproc_',preproc_config_file,'.cfg')), preproc_config_file=conn_prepend('pipeline_preproc_',preproc_config_file,'.cfg'); end
end
end

function model_config_file = el_readmodelconfig(model_config_file,defaults);
model_config_file=char(model_config_file);
if isempty(model_config_file)
    model_config_file = fullfile(defaults.folder_pipelines,'pipeline_model_Default.cfg');
else
    model_config_file=conn_prepend('',model_config_file,'.cfg');
    if isempty(fileparts(model_config_file)),
        model_config_file=fullfile(defaults.folder_pipelines,model_config_file);
        if ~conn_existfile(model_config_file)&&conn_existfile(conn_prepend('pipeline_model_',model_config_file,'.cfg')), model_config_file=conn_prepend('pipeline_model_',model_config_file,'.cfg'); end
    end
end
end

function dataset=el_readpipeline(dataset,subject,subject_path,defaults);
dataset=char(dataset);
if isempty(dataset)
    files=conn_dir(fullfile(subject_path,'*.mat'),'-ls');
    assert(~isempty(files), 'unable to find any *.mat files in %s\n',subject_path);
    if numel(files)>1, [nill,tnames]=cellfun(@fileparts,files,'uni',0); [nill,idx]=sort(tnames); files=files(idx); end
    dataset=files{end};
elseif ~isempty(regexp(dataset,'\.mat$')) % .mat file
    dataset=conn_fullfile(dataset);
else % pipeline id
    dataset=fullfile(defaults.folder_subjects,subject,[dataset,'.mat']);
    if defaults.isremote, dataset=fullfile('/CONNSERVER',dataset); end
end
end

function expt=el_readexpt(expt,dataset);
expt=char(expt);
if isempty(expt),
    expt=conn_dir(fullfile(conn_prepend('',dataset,''),'results','firstlevel','*'),'^[^\.]*','-dir','-cell','-R');
    assert(~isempty(expt),'unable to find first-level analysis results in %s',fullfile(conn_prepend('',dataset,''),'results','firstlevel'));
    expt=expt{end};
end
end

function fpath = el_readtaskfolder(defaults)
fpath = defaults.folder_tasks;
if defaults.isremote, fpath=fullfile('/CONNSERVER',fpath); end
end




