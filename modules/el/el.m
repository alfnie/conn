function varargout=el(option,varargin)
% EVLAB tools
%
% INITIALIZATION SYNTAX:
%
%   el root.subjects <subjects_folder>   : defines root directory where subject folders are stored
%                                    Altnernatively, this is defined by the symbolic link conn/modules/el/root.subjects if it exists
%                                    (syntax "el root.subjects <subjects_folder> all" may be used to create this symbolic link)
%   el root.tasks <tasks_folder>         : defines root directory where experimental design files are stored
%                                    Altnernatively, this is defined by the symbolic link conn/modules/el/root.tasks if it exists
%                                    (syntax "el root.tasks <tasks_folder> all" may be used to create this symbolic link)
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
%   el('model',subjectID, pipelineID, designID [, modelOPTIONS]) runs first-level GLM model estimation step
%      subjectID                : subject folder name (e.g. '408_FED_20160617a_3T2')
%      pipelineID               : preprocessing pipeline that will be used as source of data in this analysis
%      designID                  : .cfg file defining the experimental design associated with each functional run
%                                   Alternatively, a designID is just a shortcut to a .cfg file located in the subject directory and named "*_<designID>.cfg"
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
%   el('model.plot',subjectID, pipelineID, designID) displays first-level effect-size estimates
%   el('model.qa',subjectID, pipelineID, designID) creates QA plots on first-level GLM analyses
%   el('qa.plot',subjectID, pipelineID) displays already-created QA plots
%
%
% CONFIGURATION OPTIONS:
%
%   foldername = el('default','folder_subjects' [,foldername])
%      foldername               : default root.subjects directory where subject folders may be found. By default this is defined as conn/modules/el/root.subjects
%   foldername = el('default','folder_tasks' [,foldername])
%      foldername               : default root.tasks directory where task files may be found. By default this is defined as conn/modules/el/root.tasks
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
        'folder_dicoms','*dicoms' ,...
        'dicom_isstructural',{{'^T1_MPRAGE_1iso'}} ,...
        'dicom_disregard_functional',{{'^localizer','^AAScout','^AAHScout','^MoCoSeries','^T1_MPRAGE_1iso','^DIFFUSION_HighRes'}} );
end

conn_module('evlab17','init','silent');
fileout=[];
varargout=cell(1,nargout);

switch(lower(option))
    case 'root.subjects'
        if numel(varargin)>1&&isequal(varargin{2},'all'), conn_fileutils('linkdir',varargin{1},fullfile(fileparts(which(mfilename)),'root.subjects')); end
        if numel(varargin)>0, defaults.folder_subjects=varargin{1};
        else varargout={defaults.folder_subjects};
        end
    case 'root.tasks'
        if numel(varargin)>1&&isequal(varargin{2},'all'), conn_fileutils('linkdir',varargin{1},fullfile(fileparts(which(mfilename)),'root.tasks')); end
        if numel(varargin)>0, defaults.folder_tasks=varargin{1};
        else varargout={defaults.folder_tasks};
        end
        
    case 'default'
        if isempty(varargin), varargout={defaults}; 
        elseif isfield(defaults,varargin{1})
            if numel(varargin)>1, defaults.(varargin{1})=varargin{2};
            else varargout={defaults.(varargin{1})};
            end
        else
            if numel(varargin)>1, conn_module('evlab17','default',varargin{:});
            else varargout={conn_module('evlab17','default',varargin{1})};
            end
        end
    case 'submit'
        if ~nargout, conn('submit',@el,varargin{:}); % e.g. el submit preprocessing 408_FED_20160617a_3T2
        else [varargout{1:nargout}]=conn('submit',@el,varargin{:});
        end
            
    case {'preprocessing','preprocessing.test'}
        % adapted from msieg preprocess_PL2017
        % el('preprocessing',subject_id [, preprocessing_pipeline_file])
        % e.g. el preprocessing 408_FED_20160617a_3T2
        
        pwd0=pwd;
        subject=char(varargin{1}); % subject id
        subject_path=fullfile(defaults.folder_subjects,subject);
        data_config_file=fullfile(subject_path,'data.cfg');
        %subject_keys=regexp(subject,'_','split');
        %subject_path_dicoms = fullfile(subject_path,strjoin([subject_keys(2:4),{'dicoms'}],'_')); % 268_FED_20170929a_3T2_PL2017_unsmoothed/FED_20170929a_3T2_dicoms, 268_KAN_EvDB_20150317a_PL2017/KAN_EvDB_20150317a_dicoms, 230_KAN_prodsemphon_12_PL2017/KAN_prodsemphon_12_dicoms, 183_POLY01_20160420_3T1/POLY01_20160420_3T1_dicoms

        if strcmpi(option,'preprocessing.test')||~conn_existfile(data_config_file)
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
            fprintf('Ready to run preprocessing\n');
            varargout{1}=data_config_file;
        else
            info=conn_loadcfgfile(data_config_file);
            if isfield(info,'structurals'), struct_run=info.structurals;
            else struct_run=[];
            end
        end
        if strcmpi(option,'preprocessing')
            if numel(varargin)<2||isempty(varargin{2})
                if isempty(struct_run), preproc_config_file = fullfile(fileparts(which(mfilename)),'pipeline_preproc_DefaultMNI.cfg');
                else preproc_config_file = fullfile(fileparts(which(mfilename)),'pipeline_preproc_DefaultMNI_PlusStructural.cfg');
                end
                assert(conn_existfile(preproc_config_file),'unable to find preprocessing pipeline %s',preproc_config_file);
            else
                preproc_config_file = varargin{2};
                preproc_config_file=conn_prepend('',preproc_config_file,'.cfg');
                if isempty(fileparts(preproc_config_file)), 
                    preproc_config_file=fullfile(fileparts(which(mfilename)),preproc_config_file); 
                    if ~conn_existfile(preproc_config_file), preproc_config_file=conn_prepend('pipeline_preproc_',preproc_config_file,'.cfg'); end
                end
                assert(conn_existfile(preproc_config_file),'unable to find preprocessing pipeline %s',varargin{2});
            end
            
            %run preproc
            if strcmpi(option,'preprocessing.main') 
                [ok,msg]=mkdir(fullfile(subject_path,'nii'));
                cd(fullfile(subject_path,'nii'));
                fileout=conn_module('evlab17','run_preproc',data_config_file,preproc_config_file,[],varargin(3:end));
            else % copies data to pipeline-specific folder before preprocessing
                [nill,tname]=fileparts(preproc_config_file);
                tname=regexprep(tname,'^pipeline_preproc_','');
                dataset=fullfile(subject_path,[tname,'.mat']);
                [ok,msg]=mkdir(fileparts(dataset));
                cd(fileparts(dataset));
                fileout=conn_module('evlab17','run_preproc',data_config_file,[],'dataset',dataset);
                conn_module('evlab17','save',dataset);
                conn_module('evlab17','update'); % import
                conn_importvol2bids(2);
                conn save;
                fileout=conn_module('evlab17','run_preproc',preproc_config_file,[],'dataset',dataset,varargin(3:end));
            end
            cd(pwd0);
            varargout{1}=fileout;
        end
    
    case 'preprocessing.append'
        % el('preprocessing.append',subject_id, preprocessing_pipeline_original, preprocessing_pipeline_append)
        
        pwd0=pwd;
        subject=char(varargin{1}); % subject id
        if isempty(regexp(subject,'\.mat$'))
            if numel(varargin)>=2&&~isempty(varargin{2}),
                pipeline_id=char(varargin{2}); % pipeline id
                dataset=fullfile(defaults.folder_subjects,subject,[pipeline_id,'.mat']);
            else
                subject_path=fullfile(defaults.folder_subjects,subject);
                files=conn_dir(fullfile(subject_path,'*.mat'),'-ls');
                assert(~isempty(files), 'unable to find any *.mat files in %s\n',subject_path);
                if numel(files)>1, [nill,tnames]=cellfun(@fileparts,files,'uni',0); [nill,idx]=sort(tnames); files=files(idx); end
                dataset=files{end};
            end
        else dataset=conn_fullfile(subject); 
        end
        conn_module('evlab17','load',dataset);
        preproc_config_file = char(varargin{3}); % preprocessing pipeline
        preproc_config_file=conn_prepend('',preproc_config_file,'.cfg');
        if isempty(fileparts(preproc_config_file)),
            preproc_config_file=fullfile(fileparts(which(mfilename)),preproc_config_file);
            if ~conn_existfile(preproc_config_file), preproc_config_file=conn_prepend('pipeline_preproc_',preproc_config_file,'.cfg'); end
        end
        assert(conn_existfile(preproc_config_file),'unable to find preprocessing pipeline %s',varargin{2});
        if isempty(preproc_config_file), opts={[]};
        else opts={preproc_config_file, []};
        end
        %run preproc
        cd(fileparts(dataset));
        fileout=conn_module('evlab17','run_preproc',opts{:},'dataset',dataset, varargin(4:end));
        cd(pwd0);
        varargout{1}=fileout; 
    
    case {'preprocessing.qa','model.qa','qa.preprocessing','qa.model'}
        subject=char(varargin{1}); % subject id
        if isempty(regexp(subject,'\.mat$'))
            if numel(varargin)>=2&&~isempty(varargin{2}),
                pipeline_id=char(varargin{2}); % pipeline id
                dataset=fullfile(defaults.folder_subjects,subject,[pipeline_id,'.mat']);
            else
                subject_path=fullfile(defaults.folder_subjects,subject);
                files=conn_dir(fullfile(subject_path,'*.mat'),'-ls');
                assert(~isempty(files), 'unable to find any *.mat files in %s\n',subject_path);
                if numel(files)>1, [nill,tnames]=cellfun(@fileparts,files,'uni',0); [nill,idx]=sort(tnames); files=files(idx); end
                dataset=files{end};
            end
        else dataset=conn_fullfile(subject); 
        end
        if ~isempty(regexp(lower(option),'model')), 
            if numel(varargin)<3||isempty(varargin{3}), 
                expt=conn_dir(fullfile(fileparts(dataset),'results','firstlevel','*'),'^[^\.]*','-dir','-cell','-R');
                assert(~isempty(expt),'unable to find first-level analysis results in %s',fullfile(fileparts(dataset),'results','firstlevel'));
                expt=expt{end};
            else expt=char(varargin{3}); % design
            end
            opts={'qa_plist',['firstlevel_',expt]};
        else opts={'qa_plist','preprocessing'};
        end
        conn_module('evlab17','run_qa','dataset',dataset,opts{:},varargin(3:end));
        
    case {'qa.plot','preprocessing.qa.plot','model.qa.plot'}
        subject=char(varargin{1}); % subject id
        if isempty(regexp(subject,'\.mat$'))
            if numel(varargin)>=2&&~isempty(varargin{2}),
                pipeline_id=char(varargin{2}); % pipeline id
                dataset=fullfile(defaults.folder_subjects,subject,[pipeline_id,'.mat']);
            else
                subject_path=fullfile(defaults.folder_subjects,subject);
                files=conn_dir(fullfile(subject_path,'*.mat'),'-ls');
                assert(~isempty(files), 'unable to find any *.mat files in %s\n',subject_path);
                if numel(files)>1, [nill,tnames]=cellfun(@fileparts,files,'uni',0); [nill,idx]=sort(tnames); files=files(idx); end
                dataset=files{end};
            end
        else dataset=conn_fullfile(subject); 
        end
        conn_module('evlab17','qaplots',dataset,varargin(3:end));
        
    case {'model','firstlevel'}
        pwd0=pwd;
        subject=char(varargin{1}); % subject id
        if isempty(regexp(subject,'\.mat$'))
            if numel(varargin)>=2&&~isempty(varargin{2}),
                pipeline_id=char(varargin{2}); % pipeline id
                dataset=fullfile(defaults.folder_subjects,subject,[pipeline_id,'.mat']);
            else
                subject_path=fullfile(defaults.folder_subjects,subject);
                files=conn_dir(fullfile(subject_path,'*.mat'),'-ls');
                assert(~isempty(files), 'unable to find any *.mat files in %s\n',subject_path);
                if numel(files)>1, [nill,tnames]=cellfun(@fileparts,files,'uni',0); [nill,idx]=sort(tnames); files=files(idx); end
                dataset=files{end};
            end
        else dataset=conn_fullfile(subject); 
        end
        subject_path=fileparts(dataset);
        expt=char(varargin{3}); % design
        if numel(varargin)<4||isempty(varargin{4})
            model_config_file = fullfile(fileparts(which(mfilename)),'pipeline_model_Default.cfg');
            assert(conn_existfile(model_config_file),'unable to find model estimation options %s',model_config_file);
        else
            model_config_file = varargin{4};
            model_config_file=conn_prepend('',model_config_file,'.cfg');
            if isempty(fileparts(model_config_file)),
                model_config_file=fullfile(fileparts(which(mfilename)),model_config_file);
                if ~conn_existfile(model_config_file), model_config_file=conn_prepend('pipeline_model_',model_config_file,'.cfg'); end
            end
            assert(conn_existfile(model_config_file),'unable to find preprocessing pipeline %s',varargin{4});
        end
        
        % adapted from msieg firstlevel_PL2017
        contrasts_file=fullfile(subject_path,['contrasts_',expt,'.txt']);
        if ~conn_existfile(contrasts_file), contrasts_file=fullfile(subject_path,['contrast_',expt,'.txt']); end
        if conn_existfile(contrasts_file), % contrast definitions
            str=conn_fileutils('fileread',all_contrasts_files);
            str=reshape(regexp(str,'\n','split'),1,[]);
            cons=reshape(str(cellfun('length',str)>0),1,[]);
        else % back-compatibility
            all_contrasts_files=fullfile(defaults.folder_tasks,'contrasts_by_expt.txt'); % single-file, all expt contrasts
            if ~conn_existfile(all_contrasts_files), all_contrasts_files=fullfile(fileparts(defaults.folder_subjects),'ANALYSIS','contrasts_by_expt.txt'); end
            str=conn_fileutils('fileread',all_contrasts_files);
            str=reshape(regexp(str,'\n','split'),1,[]);
            emptyspaces=cellfun('length',str)==0;
            idx=find([true emptyspaces(1:end-1) true]&[~emptyspaces true]);
            con=strmatch(expt,str(idx(1:end-1)),'exact');
            assert(numel(con)==1,'found %d matches to %s in %s',numel(con),expt,all_contrasts_files);
            cons=str(idx(con)+1:idx(con+1)-3);
        end
        cat_file=conn_dir(fullfile(subject_path,['*_',expt,'.cat']),'-ls');
        assert(numel(cat_file)==1,'%d %s files found',numel(cat_file),fullfile(subject_path,['*_',expt,'.cat']));
        cat_info=conn_loadcfgfile(char(cat_file),struct('path',defaults.folder_tasks));
        DOSAVECFG=false;
        opts=struct('dataset',char(dataset),'design',char(cat_file),'model_name',expt,'contrasts',{cons});
        if DOSAVECFG
            model_definition_file=conn_prepend('',char(cat_file),'.cfg');
            conn_savecfgfile(model_definition_file,opts);            
            conn_module('evlab17','run_model',model_definition_file,model_config_file,varargin{5:end});
        else
            conn_module('evlab17','run_model',opts,model_config_file,varargin{5:end});
        end
        cd(pwd0);

    case {'model.plot'}
        assert(numel(varargin)>=3,'incorrect usage: please specify subject_id, pipeline_id, and firstlevel_id')
        subject=char(varargin{1}); % subject id
        if isempty(regexp(subject,'\.mat$'))
            if numel(varargin)>=2&&~isempty(varargin{2}),
                pipeline_id=char(varargin{2}); % pipeline id
                dataset=fullfile(defaults.folder_subjects,subject,[pipeline_id,'.mat']);
            else
                subject_path=fullfile(defaults.folder_subjects,subject);
                files=conn_dir(fullfile(subject_path,'*.mat'),'-ls');
                assert(~isempty(files), 'unable to find any *.mat files in %s\n',subject_path);
                if numel(files)>1, [nill,tnames]=cellfun(@fileparts,files,'uni',0); [nill,idx]=sort(tnames); files=files(idx); end
                dataset=files{end};
            end
        else dataset=conn_fullfile(subject); 
        end
        expt=char(varargin{3}); % design
        opts={dataset,expt,...
            varargin{4:end}};
        assert(conn_existfile(dataset),'file %s not found',dataset);
        conn_module('evlab17','modelplots',opts{:});
        
    case 'model.qa'
        if ~isempty(design_id), opts={'qa_plist',['firstlevel_',design_id],opts{:}};
        else opts={'qa_plist','preprocessing',opts{:}};
        end
        
    otherwise
        [varargout{1:nargout}]=conn_module('evlab17',option,varargin{:});
end

end
