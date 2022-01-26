
function varargout=fl(STEPS,varargin)
% FL Frank Lab fMRI preprocessing/QA/analysis procedures
%
% BASIC COMMANDS:
%
% ************************************************************************************
%
% fl('INIT')
% initializes all FL toolboxes and performs full check of software-versions
% (this is done automatically everytime that other fl commands are run, but users may use this syntax in their own scripts
% to make sure that their functions use the same software-versions as the FL functions)
% 
% ************************************************************************************
%
% fl('IMPORT',subject_id)
%   imports subject data from dicoms/niftis
%   Inputs:
%       subject_id : string of the form [Experiment_code][Subject_number][OptionalSubject_code] identifying individual subject/session
%   Input config files:
%       [subject_id].cfg file in $FLDATA/[Experiment_code]/config/FL/IMPORT detailing original functional/structural files and first-level design information (see subject_info in HELP FL)
%   Output in $FLDATA/[Experiment_code]/sub-[subject_id]
%
%   fl('IMPORT.DICOM2NII',subject_id) 
%      performs only the dicom-to-nifti conversion of subject /scanner/dicom folder
%
%   EXAMPLE: fl IMPORT XMP01
%
% ************************************************************************************
%
% fl('PREPROCESSING',subject_id,pipeline_id)
%   preprocessing of subject data
%   Inputs:
%       subject_id : string of the form [Experiment_code][Subject_number][OptionalSubject_code] identifying individual subject/session
%       pipeline_id : string identifying an individual preprocessing pipeline
%   Input config files:
%       preprocessing_[pipeline_id].cfg file in $FLDATA/[Experiment_code]/config/FL/PREPROCESSING detailing data preprocessing steps and options (see pipeline_info in HELP FL_INTERNAL)
%       (optional) preprocessing_[subject_id]_[pipeline_id].cfg file in $FLDATA/[Experiment_code]/config/FL/PREPROCESSING containing any subject-specific preprocessing steps and options
%   Output in $FLDATA/[Experiment_code]/derivatives/FL/[pipeline_id]/sub-[subject_id]
%
%   fl('PREPROCESSING',subject_id,pipeline_id,'parallel',1)
%      submits preprocessing job to cluster (instead of running it locally)
%   fl('PREPROCESSING.REPORT',subject_id,pipeline_id) 
%      checks the status of submitted job
%   fl('PREPROCESSING.REPORT.GUI',subject_id,pipeline_id) 
%      checks the status of submitted job
%   fl('PREPROCESSING.DELETE',subject_id,pipeline_id) 
%      cancels&deletes submitted job
%   fl('PREPROCESSING.APPEND',subject_id,pipeline_id,extra_pipeline_id) 
%      runs additional preprocessing steps on an existing dataset (results stored in same preprocessing pipeline)
%   fl('PREPROCESSING.BRANCH',subject_id,pipeline_id,extra_pipeline_id) 
%      same as preprocessing.append but creating a new preprocessing pipeline "extra_pipeline_id" with the results
%
%   EXAMPLE: fl PREPROCESSING XMP01 mnispace
%   EXAMPLE: fl PREPROCESSING XMP* mnispace
%
% ************************************************************************************
%
% fl('FIRSTLEVEL',subject_id,pipeline_id,model_id)
%   first-level analysis of subject data
%   Inputs:
%       subject_id : string of the form [Experiment_code][Subject_number][OptionalSubject_code] identifying individual subject/session
%       pipeline_id : string identifying an individual preprocessing pipeline
%       model_id : string identifying an individual first-level analysis
%   Input config files:
%       firstlevel_[model_id].cfg file in $FLDATA/[Experiment_code]/config/FL/FIRSTLEVEL detailing first level estimation options (see design_info in HELP FL)
%       (optional) firstlevel_[subject_id]_[model_id].cfg file in $FLDATA/[Experiment_code]/config/FL/FIRSTLEVEL detailing design information for each subject (if exists, this information will take precedence over information in the original [subject_id].cfg file used during the IMPORT step; see subject_info #design field in HELP FL)
%   Output in $FLDATA/[Experiment_code]/derivatives/FL/[pipeline_id]/sub-[subject_id]/results/firstlevel/[model_id]
%
% fl('FIRSTLEVEL.PLOT',subject_id,pipeline_id,model_id [,contrast_name])
%   displays first-level contrast estimation of subject data
% fl('FIRSTLEVEL.CONTRAST',subject_id,pipeline_id,model_id [,contrast_name])
%   runs only contrast-estimation step (skips first-level analysis estimation step)
%
%   EXAMPLE: fl FIRSTLEVEL XMP01 mnispace speechrate
%   EXAMPLE: fl FIRSTLEVEL XMP* mnispace speechrate
%
% ************************************************************************************
%
% fl('SECONDLEVEL',experiment_id,pipeline_id,model_id,results_id)
%   second-level analysis of subject data
%   Inputs:
%       experiment_id : string identifying an individual experiment ([Experiment_code] portion of subject_id)
%       pipeline_id : string identifying an individual preprocessing pipeline
%       model_id : string identifying an individual first-level analysis
%   Input config files:
%       secondlevel_[results_id].cfg file in $FLDATA/[Experiment_code]/config/FL/SECONDLEVEL detailing second level estimation options
%   Output in $FLDATA/[Experiment_code]/derivatives/FL/[pipeline_id]/results/secondlevel/[results_id]
%
% fl('SECONDLEVEL.PLOT',experiment_id,pipeline_id,model_id,results_id)
%   displays second-level results (previously computed) 
%
% fl('SECONDLEVEL.ROI',experiment_id,pipeline_id,model_id,results_id,roi_id)
%   displays second-level results (previously computed) aggregated within user-defined ROIs
%   Inputs:
%       roi_id : string identifying an individual ROI or atlas file 
%   Input config files:
%       roi_[roi_id].cfg file in $FLDATA/[Experiment_code]/config/FL/IMPORT detailing roi/atlas information
%
%   EXAMPLE: fl SECONDLEVEL XMP mnispace speechrate groupcomparison
%
% ************************************************************************************
%
% fl('QA.CREATE',subject_id,pipeline_id [,model_id] [,arg_name,arg_value,...])
%   creates QA plots
%   Inputs:
%       subject_id : string of the form [Experiment_code][Subject_number][OptionalSubject_code] identifying individual subject/session
%       pipeline_id : string identifying an individual preprocessing pipeline
%       (optional) model_id : string identifying an individual first-level analysis
%       (optional) 'qa_plist',labels : list of QA plots to create (see "help conn_qaplots" for additional details)
%       (optional) 'qa_set',label : functional dataset (label or index) for preprocessing qa plots
%
%   EXAMPLE: fl QA.CREATE XMP01 mnispace
%
% ************************************************************************************
%
% fl('QA.PLOT',subject_id,pipeline_id)
% fl('QA.PLOTS',Experiment_code,pipeline_id)
%   displays QA plots
%   Inputs:
%       subject_id : string of the form [Experiment_code][Subject_number][OptionalSubject_code] identifying individual subject/session
%       pipeline_id : string identifying an individual preprocessing pipeline
%
%   EXAMPLE: fl QA.PLOTS XMP mnispace
%
% ************************************************************************************
%
% fl('CONVERT.DICOM',subject_id)
%   converts dicom files to nifti format
%   Inputs:
%       dicom files in $FLDATA/[Experiment_code]/sub-[subject_id]/scanner folder
%   Output nifti files run-#.nii in $FLDATA/[Experiment_code]/sub-[subject_id]/scanner folder
%   See conn_dcmconvert_*.log file in same output folder for details about each DICOM session
%
%   EXAMPLE: fl CONVERT.DICOM XMP01
%
% ************************************************************************************
%
% fl('OPEN',subject_id)
% fl('OPEN',subject_id,pipeline_id)
%   opens CONN GUI with subject data
%   Inputs:
%       subject_id : string of the form [Experiment_code][Subject_number][OptionalSubject_code] identifying individual subject/session
%       pipeline_id : string identifying an individual preprocessing pipeline
%
%   EXAMPLE: fl OPEN XMP mnispace
%
% ************************************************************************************
%
% fl('LIST',Experiment_code)
%   lists all subjects within this experiment
% fl('LIST',Experiment_code,pipeline_id)
%   lists all subjects within this preprocessing pipeline
% fl('LIST',Experiment_code,pipeline_id,model_id)
%   lists all subjects within this preprocessing pipeline
% fl('LIST',subject_id)
%   lists all pipelines derived from this subject data
% fl('LIST',subject_id,pipeline_id)
%   lists all first-level analyses derived from this preprocessed dataset
% fl('LIST',subject_id,pipeline_id,model_id)
%   lists all effect-names and contrast-names modeled in this first-level analysis
%   Inputs:
%       subject_id : string of the form [Experiment_code][Subject_number][OptionalSubject_code] identifying individual subject/session
%       pipeline_id : string identifying an individual preprocessing pipeline
%       model_id : string identifying an individual first-level analysis
%
%   EXAMPLE: fl LIST XMP
%
% ************************************************************************************
%
% fl('FIXPERMISSIONS',subject_id)
% fl('FIXPERMISSIONS',subject_id,pipeline_id)
%   fix possible permission issues (grants read&write permissions to users within your same group) in dataset files/folders 
%   Inputs:
%       subject_id : string of the form [Experiment_code][Subject_number][OptionalSubject_code] identifying individual subject/session
%       pipeline_id : string identifying an individual preprocessing pipeline
%
%   EXAMPLE: fl FIXPERMISSIONS XMP*
%
% ************************************************************************************
%
% fl('ROOT',rootfolder) 
% defines the location of $FLDATA (root folder for all experiments) as rootfolder
%
% ************************************************************************************
%
% see also: FL_INTERNAL

persistent rootfolder isremote;

if isempty(rootfolder), rootfolder=fullfile(fileparts(which(mfilename)),'REPOSITORY'); end
if isempty(isremote), isremote=false; end
fileout=[];
varargout=cell(1,nargout);
if isremote, OUTPUT_FOLDER=fullfile('/CONNSERVER',rootfolder); 
else OUTPUT_FOLDER=rootfolder; 
end

if isremote&&~isempty(varargin)&&isempty(regexp(lower(char(STEPS)),'^secondlevel'))&&~(~isempty(regexp(lower(char(STEPS)),'plots?$'))||ismember(lower(char(STEPS)),{'root','remote','list','init','initforce','open','preprocessing.report','preprocessing.report.gui','preprocessing.delete','parallel.report','parallel.report.gui','parallel.delete','report','report.gui','delete'})); % run these locally
    [hmsg,hstat]=conn_msgbox({'Process running remotely','Please wait...',' ',' '},[],[],true);
    if ~isempty(hmsg), [varargout{1:nargout}]=conn_server('run_withwaitbar',hstat,mfilename,STEPS,varargin{:}); 
    else [varargout{1:nargout}]=conn_server('run',mfilename,STEPS,varargin{:}); 
    end
    if ~isempty(hmsg)&&ishandle(hmsg), delete(hmsg); end
    return
end

if ~isempty(varargin)&&isempty(regexp(lower(char(STEPS)),'^secondlevel'))&&~ismember(lower(char(STEPS)),{'init','initforce','root','remote','qa.plots','qaplots'}); %,'list'}) % exceptions to subject-expansion
    subject_id=varargin{1};
    if ischar(subject_id)&&any(ismember(subject_id,'*?'))
        project_id=regexprep(subject_id,'[\d\*]+.*$','');
        dataset=fullfile(OUTPUT_FOLDER,project_id,sprintf('sub-%s',subject_id));
        subject_info=fullfile(OUTPUT_FOLDER,project_id,'config','FL','IMPORT',conn_prepend('',subject_id,'.cfg'));
        f1=conn_dir(dataset,'-cell','-R','-dir'); f1=f1(cellfun('length',regexp(f1,'\.qlog$'))==0);
        f2=conn_dir(subject_info,'-cell','-R'); 
        if 1, 
            [nill,subject_id2]=cellfun(@fileparts,f2,'uni',0);
            subject_id2=conn_prepend('',subject_id2,'');
            if numel(f1)~=numel(f2), 
                fprintf('Warning: mismatched number of subject folders in %s vs. number of subject config files in %s. Subject-expansion will derive list of subjects from list of existing config files\n',dataset,subject_info); 
                [nill,subject_id1]=cellfun(@fileparts,f1,'uni',0);
                subject_id1=conn_prepend('',regexprep(subject_id1,'^sub-',''),'');
                subject_id=intersect(subject_id1,subject_id2);
            else
                subject_id=subject_id2;
            end
        else
            if numel(f1)~=numel(f2), fprintf('Warning: mismatched number of subject folders in %s vs. number of subject config files in %s. Subject-expansion will derive list of subjects from list of existing subject folders\n',dataset,subject_info); end
            subject_id=regexprep(f1,['^.*sub-(',project_id,'\d+.*?)$'],'$1');
        end
    end
    if iscell(subject_id)
        for n=1:numel(subject_id)
            fprintf('PROCESSING SUBJECT %s\n',subject_id{n});
            fl(STEPS,subject_id{n},varargin{2:end});
        end
        return
    end
end

STEPS=char(STEPS);
switch(lower(STEPS))
    case {'init','initforce'}
        conn_module('evlab17',STEPS,varargin{:});
        
    case {'fixpermissions','setpermissions'}
        assert(numel(varargin)>=1,'incorrect usage: please specify subject_id')
        subject_id=varargin{1};
        project_id=regexprep(subject_id,'[\d\*]+.*$','');
        if numel(varargin)>=2, 
            pipeline_id=varargin{2}; 
            dataset=fullfile(OUTPUT_FOLDER,project_id,'derivatives','FL',pipeline_id,sprintf('sub-%s',subject_id));
            conn_setpermissions;
            conn_fixpermissions(dataset,[],true);
            [ok,msg]=system(sprintf('chmod %s ''%s''','g+rw',fullfile(OUTPUT_FOLDER,project_id,'derivatives','FL',pipeline_id)));
            [ok,msg]=system(sprintf('chmod %s ''%s''','g+rw',fullfile(OUTPUT_FOLDER,project_id,'derivatives','FL')));
            [ok,msg]=system(sprintf('chmod %s ''%s''','g+rw',fullfile(OUTPUT_FOLDER,project_id,'derivatives')));
            [ok,msg]=system(sprintf('chmod %s ''%s''','g+rw',fullfile(OUTPUT_FOLDER,project_id)));
        else
            dataset=fullfile(OUTPUT_FOLDER,project_id,sprintf('sub-%s',subject_id));
            conn_setpermissions;
            conn_fixpermissions(dataset,[],true);
            [ok,msg]=system(sprintf('chmod %s ''%s''','g+rw',fullfile(OUTPUT_FOLDER,project_id)));
        end
        
    case 'root'
        if nargin>1&&~isempty(varargin{1}), 
            rootfolder=varargin{1}; 
            if isremote, conn_server('run',mfilename,'root',rootfolder); end
            if nargout>0, varargout={rootfolder}; end
            fprintf('ROOT folder set to %s\n',rootfolder);
        else
            if ~nargout, disp(rootfolder);
            else varargout={rootfolder};
            end
        end
        
    case 'remote'
        if nargin>1&&~isempty(varargin{1}), 
            isremote=varargin{1}; 
            if ischar(isremote), isremote=str2num(isremote); end
            if isremote&&~conn_server('isconnected')
                fprintf('Starting new remote connection to server\n');
                conn remotely on;
            end
            if isremote, fprintf('working with remote projects now\n');
            else fprintf('working with local projects now\n');
            end
            if isremote&&conn_server('isconnected')
                conn_server('run','conn_module','fl','init');
                rootfolder=conn_server('run',mfilename,'root');
                fprintf('ROOT folder set to %s\n',rootfolder);
            end
            if nargout>0, varargout={isremote}; end
        else
            if ~nargout,
                if isremote, fprintf('working with remote projects (remote=1)\n');
                else fprintf('working with local projects (remote=0)\n');
                end
            else varargout={isremote};
            end
        end
        
    case 'open'
        assert(numel(varargin)>=1,'incorrect usage: please specify subject_id')
        subject_id=varargin{1};
        
        project_id=regexprep(subject_id,'[\d\*]+.*$','');
        if numel(varargin)>1, 
            pipeline_id=varargin{2}; 
            dataset=fullfile(OUTPUT_FOLDER,project_id,'derivatives','FL',pipeline_id,sprintf('sub-%s',subject_id));
        else
            dataset=fullfile(OUTPUT_FOLDER,project_id,sprintf('sub-%s',subject_id));
        end
        conn_module evlab17 init silent;
        conn_module('evlab17','load',dataset);
        conn;
        conn('load',conn_prepend('',dataset,'.mat'));
        conn gui_setup;
        
    case {'preprocessing.report','preprocessing.report.gui','preprocessing.delete','parallel.report','parallel.report.gui','parallel.delete','report','report.gui','delete'}
        assert(numel(varargin)>=2,'incorrect usage: please specify subject_id and pipeline_id')
        subject_id=varargin{1};
        pipeline_id=varargin{2};
        project_id=regexprep(subject_id,'\d+.*$','');
        dataset=fullfile(OUTPUT_FOLDER,project_id,'derivatives','FL',pipeline_id,sprintf('sub-%s',subject_id));
        if ~isremote, [ok,msg]=system('qstat -u $USER');if ~ok, disp(msg); end; end
        conn_module evlab17 init silent;
        conn_module('evlab17','load',dataset);
        if ~isempty(regexp(lower(STEPS),'report.gui')), conn_jobmanager;
        else conn_jobmanager(regexprep(lower(STEPS),'^preprocessing\.|^parallel\.',''));
        end
        
    case {'convert.dicom'}  % convert DICOM to NIFTI
        assert(numel(varargin)>=1,'incorrect usage: please specify subject_id')
        subject_id=varargin{1};
        project_id=regexprep(subject_id,'[\d\*]+.*$','');
        dataset=fullfile(OUTPUT_FOLDER,project_id,sprintf('sub-%s',subject_id));
        conn_dcmconvert(fullfile(dataset,'scanner','*'),'folderout',fullfile(dataset,'scanner'));
        %if strcmpi(STEPS,'import.dicom'), fl('import.func',varargin{:}); end
        
    case {'import','import.func','import.anat','import.dicoms','import.dicom'}  % import ACE01
        assert(numel(varargin)>=1,'incorrect usage: please specify subject_id')
        subject_id=varargin{1};
        project_id=regexprep(subject_id,'[\d\*]+.*$','');
        dataset=fullfile(OUTPUT_FOLDER,project_id,sprintf('sub-%s',subject_id));
        subject_info=fullfile(OUTPUT_FOLDER,project_id,'config','FL','IMPORT',conn_prepend('',subject_id,'.cfg'));
        issubject_info=conn_existfile(subject_info);
        opts={};
        switch(lower(STEPS))
            case {'import.dicoms','import.dicom'}
                if conn_existfile(fullfile(dataset,'scanner','dicom'),2)
                    opts={'PREPROC.IMPORT',dataset,...
                        'dicoms',fullfile('dicom','*'),...
                        'dicoms_path','scanner',...
                        'path',fullfile(OUTPUT_FOLDER,project_id,'config','FL','IMPORT'),...
                        varargin{2:end}};
                    if issubject_info, opts=[opts,{'subject_info',subject_info}]; end
                end
            case {'import.func','import.anat'}
                if conn_existfile(fullfile(dataset,'scanner','anat'),2)&&conn_existfile(fullfile(dataset,'scanner','func'),2)
                    opts={'PREPROC.IMPORT',dataset,...
                        'structurals',fullfile('anat','*.nii'),...
                        'functionals',fullfile('func','*.nii'),...
                        'structurals_path','scanner',...
                        'functionals_path','scanner',...
                        'rois_path','scanner',...
                        'dicoms_path','scanner',...
                        'path',fullfile(OUTPUT_FOLDER,project_id,'config','FL','IMPORT'),...
                        varargin{2:end}};
                    if issubject_info, opts=[opts,{'subject_info',subject_info}]; end
                end
            otherwise
                if issubject_info,
                    opts={'PREPROC.IMPORT',dataset,...
                        'subject_info',subject_info,...
                        'structurals_path','scanner',...
                        'functionals_path','scanner',...
                        'rois_path','scanner',...
                        'dicoms_path','scanner',...
                        'path',fullfile(OUTPUT_FOLDER,project_id,'config','FL','IMPORT'),...
                        'isnew',1,...
                        varargin{2:end}};
                end
        end
        assert(~isempty(opts), 'Could not find neither a subject info file in %s or anat/func/dicom folders in %s',subject_info,dataset);
        fprintf('Project id %s\n',project_id);
        fprintf('Subject id %s\n',subject_id);
        fprintf('Output subject folder %s\n',dataset);
        if issubject_info, fprintf('Subject information read from %s\n',subject_info);
        else fprintf('Subject information implicitly defined (file %s not found)\n',subject_info);
        end
        fprintf('FL_INTERNAL command syntax:\n'); for n=1:numel(opts), fprintf('   '); disp(opts{n}); end
        fl_internal(opts{:});
        fprintf('Done\n');
        
    case {'preproc','preprocessing'}  % preproc ACE01 try01
        assert(numel(varargin)>=2,'incorrect usage: please specify subject_id and pipeline_id')
        subject_id=varargin{1};
        pipeline_id=varargin{2};
        [subject_idpath,subject_idname,subject_idext]=fileparts(subject_id);
        subject_id=[subject_idname,subject_idext];
        project_id=regexprep(subject_id,'[\d\*]+.*$','');
        dataset=fullfile(OUTPUT_FOLDER,project_id,'derivatives','FL',pipeline_id,sprintf('sub-%s',subject_id));
        datainput=fullfile(OUTPUT_FOLDER,project_id,subject_idpath,conn_prepend('',sprintf('sub-%s',subject_id),'.mat'));
        pipeline_info=fullfile(OUTPUT_FOLDER,project_id,'config','FL','PREPROCESSING',conn_prepend('preprocessing_',pipeline_id,'.cfg'));
        pipeline_info_opt=fullfile(OUTPUT_FOLDER,project_id,'config','FL','PREPROCESSING',conn_prepend('preprocessing_',[subject_id,'_',pipeline_id],'.cfg'));
        if conn_existfile(pipeline_info_opt), 
            if ~conn_existfile(pipeline_info), pipeline_info=struct; 
            else pipeline_info=conn_loadcfgfile(pipeline_info); 
            end
            pipeline_info=conn_loadcfgfile(pipeline_info_opt,pipeline_info);
        end
        opts=varargin(3:end);
        if rem(numel(opts),2), opts=[opts,{[]}]; end
        opts={'PREPROC',dataset,....
            'subject_info',datainput,...
            'pipeline_info',pipeline_info,...
            opts{:}};
        fprintf('Project id %s\n',project_id);
        fprintf('Subject id %s\n',subject_id);
        fprintf('Pipeline id %s\n',pipeline_id);
        fprintf('Output subject folder %s\n',dataset);
        fprintf('Input subject folder %s\n',datainput);
        if ~isstruct(pipeline_info)&&~conn_existfile(pipeline_info), 
            tpipeline_info=fullfile(fileparts(which(mfilename)),conn_prepend('preprocessing_',pipeline_id,'.cfg'));
            if ~isempty(tpipeline_info), 
                fprintf('Warning: no %s or %s pipeline information files found. Copying pipeline definition from %s to %s\n',pipeline_info,pipeline_info_opt,tpipeline_info,pipeline_info); 
                if ispc, [ok,msg]=system(sprintf('copy "%s" "%s"',tpipeline_info,pipeline_info));
                else [ok,msg]=system(sprintf('cp ''%s'' ''%s''',tpipeline_info,pipeline_info));
                end
            else error('No %s or %s pipeline information files found',pipeline_info,pipeline_info_opt); 
            end
            %fprintf('No pipeline info found in %s. Preprocessing pipeline %s will be only initialized (no steps run)\n',pipeline_info);
            %pipeline_info=struct;
        elseif isstruct(pipeline_info), fprintf('Preprocessing pipeline info read from %s\n',pipeline_info_opt);
        else fprintf('Preprocessing pipeline info read from %s\n',pipeline_info);
        end
        fprintf('FL_INTERNAL command syntax:\n'); for n=1:numel(opts), fprintf('   '); disp(opts{n}); end
        assert(conn_existfile(datainput,1),'file %s not found',datainput);
        fl_internal(opts{:});
        fprintf('Done\n');
                
    case {'preproc.branch','preprocessing.branch'}  % preproc.branch ACE01 try01 try02
        assert(numel(varargin)>=3,'incorrect usage: please specify subject_id, pipeline_id, and AdditionalStepsPipeline_id')
        subject_id=varargin{1};
        pipeline_id=varargin{2};
        pipeline_idadd=varargin{3};
        opts=varargin(4:end);
        if rem(numel(opts),2), opts=[opts,{[]}]; end
        fl('preprocessing',fullfile('derivatives','FL',pipeline_id,subject_id),pipeline_idadd,opts{:});
        
    case {'preproc.append','preprocessing.append'}  % preproc.append ACE01 try01 try02
        assert(numel(varargin)>=3,'incorrect usage: please specify subject_id, pipeline_id, and AdditionalStepsPipeline_id')
        subject_id=varargin{1};
        pipeline_id=varargin{2};
        pipeline_idadd=varargin{3};
        project_id=regexprep(subject_id,'[\d\*]+.*$','');
        dataset=fullfile(OUTPUT_FOLDER,project_id,'derivatives','FL',pipeline_id,sprintf('sub-%s',subject_id));
        pipeline_info=fullfile(OUTPUT_FOLDER,project_id,'config','FL','PREPROCESSING',conn_prepend('preprocessing_',pipeline_idadd,'.cfg'));
        pipeline_info_opt=fullfile(OUTPUT_FOLDER,project_id,'config','FL','PREPROCESSING',conn_prepend('preprocessing_',[subject_id,'_',pipeline_idadd],'.cfg'));
        if conn_existfile(pipeline_info_opt), 
            if ~conn_existfile(pipeline_info), pipeline_info=struct; 
            else pipeline_info=conn_loadcfgfile(pipeline_info); 
            end
            pipeline_info=conn_loadcfgfile(pipeline_info_opt,pipeline_info);
        end
        opts=varargin(4:end);
        if rem(numel(opts),2), opts=[opts,{[]}]; end
        opts={'PREPROC.APPEND',dataset,....
            'pipeline_info',pipeline_info,...
            opts{:}};
        fprintf('Project id %s\n',project_id);
        fprintf('Subject id %s\n',subject_id);
        fprintf('Pipeline id %s\n',pipeline_id);
        fprintf('Output subject folder %s\n',dataset);
        if ~isstruct(pipeline_info)&&~conn_existfile(pipeline_info), error('No %s or %s additional pipeline information files found',pipeline_info,pipeline_info_opt); 
            %fprintf('No pipeline info found in %s. Preprocessing pipeline %s will be only initialized (no steps run)\n',pipeline_info);
            %pipeline_info=struct;
        elseif isstruct(pipeline_info), fprintf('Additional preprocessing pipeline info read from %s\n',pipeline_info_opt);
        else fprintf('Additional preprocessing pipeline info read from %s\n',pipeline_info);
        end
        fprintf('FL_INTERNAL command syntax:\n'); for n=1:numel(opts), fprintf('   '); disp(opts{n}); end
        assert(conn_existfile(dataset,1),'file %s not found',dataset);
        fl_internal(opts{:});
        fprintf('Done\n');
        
    case {'model','firstlevel','firstlevel.contrast'}  % firstlevel ACE01 try01 Speech_Raw
        assert(numel(varargin)>=3,'incorrect usage: please specify subject_id, pipeline_id, and firstlevel_id')
        subject_id=varargin{1};
        pipeline_id=varargin{2};
        design_id=varargin{3};
        project_id=regexprep(subject_id,'[\d\*]+.*$','');
        dataset=fullfile(OUTPUT_FOLDER,project_id,'derivatives','FL',pipeline_id,sprintf('sub-%s',subject_id));
        analysis_info=fullfile(OUTPUT_FOLDER,project_id,'config','FL','FIRSTLEVEL',conn_prepend('firstlevel_',design_id,'.cfg'));
        analysis_info_opt=fullfile(OUTPUT_FOLDER,project_id,'config','FL','FIRSTLEVEL',conn_prepend('firstlevel_',[subject_id,'_',design_id],'.cfg'));
        subject_info=fullfile(OUTPUT_FOLDER,project_id,'config','FL','IMPORT',conn_prepend('',subject_id,'.cfg'));
        isanalysis_info_opt=conn_existfile(analysis_info_opt);
        design_info=struct;
        design_info.design=struct('path',fullfile(OUTPUT_FOLDER,project_id,'config','FL','DESIGN'));
        if conn_existfile(analysis_info), design_info=conn_loadcfgfile(analysis_info,design_info); end
        if isanalysis_info_opt, 
            design_info=conn_loadcfgfile(analysis_info_opt,design_info);
        elseif conn_existfile(subject_info), % note: if no subject-specific (e.g. config/EXAMPLE01_speech.cfg) file is found, it will re-read the info in config/EXAMPLE01.cfg looking for design info there
            design_info_temp=conn_loadcfgfile(subject_info,design_info);
            if isfield(design_info_temp,'design'), design_info.design=design_info_temp.design; end
            if isfield(design_info_temp,'files'), design_info.files=design_info_temp.files; end
            if isfield(design_info_temp,'runs'), design_info.runs=design_info_temp.runs; end
            if isfield(design_info_temp,'path'), design_info.path=design_info_temp.path; end
        end
        if strcmpi(STEPS,'firstlevel.contrast'), 
            if isfield(design_info,'design')&&isfield(design_info.design,'files'), design_info.design=rmfield(design_info.design,'files');  end
            if isfield(design_info,'files'), design_info=rmfield(design_info,'files');  end
        end
        opts={'MODEL',dataset,...
            'design_info',design_info,...
            'model_name',design_id,...
            varargin{4:end}};
        fprintf('Project id %s\n',project_id);
        fprintf('Subject id %s\n',subject_id);
        fprintf('First-level analysis id %s\n',design_id);
        fprintf('Subject folder %s\n',dataset);
        fprintf('First-level model estimation info read from %s\n',analysis_info);
        if isanalysis_info_opt, fprintf('First-level design info read from %s\n',design_info_opt);
        else fprintf('First-level design info read from %s (file %s not found)\n',subject_info);
        end
        fprintf('FL_INTERNAL command syntax:\n'); for n=1:numel(opts), fprintf('   '); disp(opts{n}); end
        assert(conn_existfile(dataset,1),'file %s not found',dataset);
        assert(conn_existfile(analysis_info),'file %s not found',analysis_info);
        assert(isanalysis_info_opt|conn_existfile(subject_info),'files %s or %s not found',analysis_info_opt,subject_info);
        fl_internal(opts{:});
        fprintf('Done\n');
        
    case {'firstlevel.plot','firstlevel.plots'}
        assert(numel(varargin)>=3,'incorrect usage: please specify subject_id, pipeline_id, and firstlevel_id')
        subject_id=varargin{1};
        pipeline_id=varargin{2};
        design_id=varargin{3};
        project_id=regexprep(subject_id,'[\d\*]+.*$','');
        dataset=fullfile(OUTPUT_FOLDER,project_id,'derivatives','FL',pipeline_id,sprintf('sub-%s',subject_id));
        opts={dataset,design_id,...
            varargin{4:end}};
        fprintf('Project id %s\n',project_id);
        fprintf('Subject id %s\n',subject_id);
        fprintf('First-level analysis id %s\n',design_id);
        fprintf('Experiment folder %s\n',dataset);
        assert(conn_existfile(dataset,1),'file %s not found',dataset);
        conn_module evlab17 init silent;
        conn_module('evlab17','modelplots',opts{:});
       
    case {'secondlevel'}  % secondlevel ACE try01 Speech_Raw Group_Comparison
        assert(numel(varargin)>=4,'incorrect usage: please specify experiment_id, pipeline_id, firstlevel_id, and results_id')
        project_id=varargin{1};
        pipeline_id=varargin{2};
        design_id=varargin{3};
        results_id=varargin{4};
        results_path=fullfile(OUTPUT_FOLDER,project_id,'derivatives','FL',pipeline_id,'results','secondlevel',design_id,results_id);
        results_info=fullfile(OUTPUT_FOLDER,project_id,'config','FL','SECONDLEVEL',conn_prepend('secondlevel_',results_id,'.cfg'));
        options=varargin(5:end);
        
        fprintf('Project id %s\n',project_id);
        fprintf('First-level analysis id %s\n',design_id);
        fprintf('Second-level results id %s\n',results_id);
        fprintf('Second-level results folder %s\n',results_path);
        assert(conn_existfile(results_info),'file %s not found',results_info);
        assert(conn_existfile(fullfile(OUTPUT_FOLDER,project_id,'derivatives','FL',pipeline_id),1),'preprocessing path %s not found',fullfile(OUTPUT_FOLDER,project_id,'derivatives','FL',pipeline_id));
        
        results_cfg=conn_loadcfgfile(results_info);
        if ~isfield(results_cfg,'subjects')&&isfield(results_cfg,'data'), 
            fprintf('warning: missing #subjects field in %s. Using #data field instead\n',results_info);
            dataset=results_cfg.data;
            subjects={};
        elseif ~isfield(results_cfg,'subjects'), 
            fprintf('#subjects field not found in %s. Assuming all subjects\n',results_info)
            f1=conn_dir(fullfile(OUTPUT_FOLDER,project_id,'derivatives','FL',pipeline_id,'sub-*'),'-cell','-R','-dir');f1=f1(cellfun('length',regexp(f1,'\.qlog$'))==0);
            assert(~isempty(f1),'no first-level results found in %s. Check syntax or re-run first-level analyses',fullfile(OUTPUT_FOLDER,project_id,'derivatives','FL',pipeline_id,'sub-*')); 
            dataset=cellfun(@(x)fullfile(x,'results','firstlevel',design_id,'SPM.mat'),f1,'uni',0);
            subjects=regexprep(f1,'^.*[\\\/]','');
        else
            if ischar(results_cfg.subjects), results_cfg.subjects=cellstr(results_cfg.subjects); end
            if numel(results_cfg.subjects)==1, f1=conn_dir(fullfile(OUTPUT_FOLDER,project_id,'derivatives','FL',pipeline_id,['sub-',regexprep(results_cfg.subjects{1},'^sub-','')]),'-cell','-R','-dir');f1=f1(cellfun('length',regexp(f1,'\.qlog$'))==0);
            else f1=cellfun(@(x)fullfile(OUTPUT_FOLDER,project_id,'derivatives','FL',pipeline_id,['sub-',regexprep(x,'^sub-','')]),results_cfg.subjects,'uni',0);
            end
            assert(~isempty(f1),'no first-level results found in %s. Check syntax or re-run first-level analyses',fullfile(OUTPUT_FOLDER,project_id,'derivatives','FL',pipeline_id,['sub-',regexprep(results_cfg.subjects{1},'^sub-','')])); 
            dataset=cellfun(@(x)fullfile(x,'results','firstlevel',design_id,'SPM.mat'),f1,'uni',0);
            subjects=regexprep(f1,'^.*[\\\/]','');
        end
        okdataset=conn_existfile(dataset,1);
        assert(~isempty(dataset), 'unable to find any first-level files');
        if any(~okdataset), fprintf('Warning: unable to find first-level files in:\n'); disp(char(dataset(~okdataset))); end
        dataset=dataset(okdataset>0);
        subjects=subjects(okdataset>0);
        fprintf('%d subjects: %s\n',numel(dataset),sprintf('%s ',subjects{:}));
        assert(numel(dataset)>1,'insufficient subjects. Unable to perform second-level analyses');
        
        assert(isfield(results_cfg,'design'),'missing #design field in %s',results_info);
        X=results_cfg.design;
        X_labels={};
        if iscell(X)
            newX=[];
            pdata_name=[];
            pdata_idx=[];
            X_labels=X;
            if conn_existfile(fullfile(OUTPUT_FOLDER,project_id,'participants.tsv')),pdata_name=fullfile(OUTPUT_FOLDER,project_id,'participants.tsv');
            elseif conn_existfile(fullfile(OUTPUT_FOLDER,project_id,'config','FL','participants.tsv')),pdata_name=fullfile(OUTPUT_FOLDER,project_id,'config','FL','participants.tsv');
            elseif conn_existfile(fullfile(OUTPUT_FOLDER,project_id,'config','FL','SECONDLEVEL','participants.tsv')),pdata_name=fullfile(OUTPUT_FOLDER,project_id,'config','FL','SECONDLEVEL','participants.tsv');
            end
            if ~isempty(pdata_name), pdata=conn_loadtextfile(pdata_name); end
            for n=1:numel(X)
                if strcmpi(X{n},'AllSubjects'), x=ones(numel(dataset),1); 
                elseif isempty(regexprep(X{n},'[\d\s\t\.-+\/]','')), x=str2num(X{n}); X_labels{n}=sprintf('unknown%d',n);
                else
                    if isempty(pdata_idx)
                        assert(~isempty(pdata_name),'unable to find behavioral data in %s',fullfile(OUTPUT_FOLDER,project_id,'participants.tsv'));
                        assert(isfield(pdata,'participant_id'),'unable to find field "%s" in %s','participant_id',pdata_name);
                        [ok,idx]=ismember(subjects,pdata.participant_id);
                        assert(all(ok),'unable to find subject id "%s" in %s',sprintf('%s ',subjects{~ok}),pdata_name);
                        pdata_idx=idx;
                    end
                    if ~isempty(regexp(X{n},'\*')), xvars=regexp(X{n},'\s*\*\s*','split'); % interactions (A*B*...)
                    else xvars={X{n}};
                    end
                    x=1;
                    for nxvars=1:numel(xvars)
                        assert(isfield(pdata,xvars{nxvars}),'unable to find field "%s" in %s',xvars{nxvars},pdata_name);
                        tx=pdata.(xvars{nxvars});
                        tx=tx(pdata_idx);
                        if iscell(tx)
                            [ux,nill,ix]=unique(regexprep(tx,'^\s+|\s+$',''));
                            fprintf('Warning: converting field %s levels (%s) to numeric values (%s)\n',xvars{nxvars}, sprintf('%s ',ux{:}),mat2str(1:numel(ux)));
                            tx=ix;
                        end
                        x=x.*tx;
                    end
                end
                assert(numel(x)==numel(dataset),'unexpected number of elements in %s (found %d, expected %d)',X{n},numel(x),numel(dataset));
                newX=cat(1,newX,x(:)');
            end
            X=newX;
        end
        fprintf('Design matrix: %s\n',mat2str(X));
        
        conn_module evlab17 init silent;
        conn_module('evlab17','run_results', ...
            results_info,...
            'data',dataset,...
            'design',X,...
            'design_labels',X_labels,...
            'subjectids',subjects,...
            'folder',results_path,...
            options{:});        
        fprintf('Done\n');
        
    case {'secondlevel.plot'}  % secondlevel ACE try01 Speech_Raw Group_Comparison
        assert(numel(varargin)>=4,'incorrect usage: please specify experiment_id, pipeline_id, firstlevel_id, and results_id')
        project_id=varargin{1};
        pipeline_id=varargin{2};
        design_id=varargin{3};
        results_id=varargin{4};
        options=varargin(5:end);
        results_path=fullfile(OUTPUT_FOLDER,project_id,'derivatives','FL',pipeline_id,'results','secondlevel',design_id,results_id);
        conn_module evlab17 init silent;
        fh=conn_module('evlab17','resultsplots',results_path,options{:});
        varargout={fh};
        
    case {'secondlevel.roi','secondlevel.rois','secondlevel.wroi','secondlevel.wrois'}  % secondlevel.roi ACE try01 Speech_Raw Group_Comparison SpeechROIs
        assert(numel(varargin)>=5,'incorrect usage: please specify experiment_id, pipeline_id, firstlevel_id, secondlevel_id , rois_id')
        project_id=varargin{1};
        pipeline_id=varargin{2};
        design_id=varargin{3};
        results_id=varargin{4};
        roi_id=varargin{5};
        results_path=fullfile(OUTPUT_FOLDER,project_id,'derivatives','FL',pipeline_id,'results','secondlevel',design_id,results_id);
        options=varargin(6:end);
        if ~conn_existfile(fullfile(results_path,'SPM.mat'))
            fprintf('Second-level results file %s does not exist\n',fullfile(results_path,'SPM.mat'));
            fl('secondlevel',project_id,pipeline_id,design_id,results_id,'analysistype',3); % note: non-param only
        end
        roi_keys={@(roiid,subid)fullfile(OUTPUT_FOLDER,project_id,'derivatives','FL',pipeline_id,sprintf('sub-%s',subid),'roi',sprintf('sub-%s_roi-%s.nii',subid,roiid)),...
                  @(roiid,subid)fullfile(OUTPUT_FOLDER,project_id,sprintf('sub-%s',subid),'roi',sprintf('sub-%s_roi-%s.nii',subid,roiid))}; % note: change to line below if ROIs can be modified by preprocessing pipeline; % note: not supported yet session-specific rois
        troi_id=[];
        roi_file={};
        if ~iscell(roi_id), roi_id={roi_id}; end
        for n=1:numel(roi_id)
            roi_file{n}=fullfile(OUTPUT_FOLDER,project_id,'config','FL','IMPORT',conn_prepend('roi_',roi_id{n},'.cfg')); % cfg file
            if conn_existfile(roi_file{n}), troi_id=roi_id{n};
            elseif conn_existfile(roi_id{n}), roi_file{n}=roi_id{n}; % roi file
            else error('unable to find configuration file %s',roi_file{n}); 
            end
        end
        if ~isempty(troi_id), options=[options,{'roi_id',troi_id}]; end % last cfg file
        if ismember(lower(char(STEPS)),{'secondlevel.wroi','secondlevel.wrois'}), options=[options,{'summary_measure','weighted mean'}]; end
        conn_module evlab17 init silent;
        conn_module('evlab17','run_roiresults', ...
            roi_file{:},...
            [],...
            'folder',results_path,...
            'roi_keys',roi_keys,...
            options{:});
        
    case {'secondlevel.roi.plot','secondlevel.rois.plot','secondlevel.wroi.plot','secondlevel.wrois.plot'}  % secondlevel.roi.plot ACE try01 Speech_Raw Group_Comparison SpeechROIs
        assert(numel(varargin)>=4,'incorrect usage: please specify experiment_id, pipeline_id, firstlevel_id, secondlevel_id , rois_id')
        project_id=varargin{1};
        pipeline_id=varargin{2};
        design_id=varargin{3};
        results_id=varargin{4};
        results_path=fullfile(OUTPUT_FOLDER,project_id,'derivatives','FL',pipeline_id,'results','secondlevel',design_id,results_id);
        if numel(varargin)>=5, 
            roi_id=varargin{5};
            rex_file=fullfile(results_path,sprintf('REX_%s.mat',roi_id));
        else
            [tfilename,tpathname]=conn_fileutils('uigetfile',{'REX*.mat','REX file (REX*.mat)'; '*',  'All Files (*)'},'Select ROI results file');
            if ~ischar(tfilename)||isempty(tfilename), return; end
            rex_file=fullfile(tpathname,tfilename);
        end
        conn_module evlab17 init silent;
        conn_module('evlab17','roiresultsplots',rex_file);
        
    case {'qa.create','qacreate'}  % qacreate ACE01 try01    OR   qacreate ACE01 try01 Speech_Raw
        assert(numel(varargin)>=2,'incorrect usage: please specify subject_id and pipeline_id [, firstlevel_id]')
        subject_id=varargin{1};
        pipeline_id=varargin{2};
        opts=varargin(3:end);
        if rem(numel(opts),2), design_id=opts{1}; opts=opts(2:end);
        else design_id='';
        end
        project_id=regexprep(subject_id,'[\d\*]+.*$','');
        if ~isempty(design_id), opts={'qa_plist',conn_prepend('firstlevel_',design_id),opts{:}};
        else opts={'qa_plist','preprocessing',opts{:}};
        end
        dataset=fullfile(OUTPUT_FOLDER,project_id,'derivatives','FL',pipeline_id,sprintf('sub-%s',subject_id));
        opts={'QA.CREATE',dataset,opts{:}};
        fprintf('Project id %s\n',project_id);
        fprintf('Subject id %s\n',subject_id);
        if ~isempty(design_id), fprintf('First-level analysis id %s\n',design_id); end
        fprintf('Subject folder %s\n',dataset);
        fprintf('FL_INTERNAL command syntax:\n'); for n=1:numel(opts), fprintf('   '); disp(opts{n}); end
        assert(conn_existfile(dataset,1),'file %s not found',dataset);
        fl_internal(opts{:});
        fprintf('Done\n');
        
    case {'qa.plots','qaplots'}  % qaplots ACE try01
        assert(numel(varargin)>=2,'incorrect usage: please specify experiment_id and pipeline_id')
        experiment_id=varargin{1};
        pipeline_id=varargin{2};
        opts=varargin(3:end);
        if rem(numel(opts),2), design_id=opts{1}; opts=opts(2:end);
        else design_id='';
        end
        project_id=regexprep(experiment_id,'\d+.*$','');
        dataset=fullfile(OUTPUT_FOLDER,project_id,'derivatives','FL',pipeline_id);
        opts={'QA.PLOTS',dataset,...
            opts{:}};
        fprintf('Project id %s\n',project_id);
        fprintf('Experiment id %s\n',experiment_id);
        fprintf('Experiment folder %s\n',dataset);
        assert(conn_existfile(dataset,1),'file %s not found',dataset);        
        conn_qaplotsexplore(dataset,opts{:},'flfolders');
                
    case {'qa.plot','qaplot'}  % qaplot ACE01 try01
        assert(numel(varargin)>=2,'incorrect usage: please specify subject_id and pipeline_id')
        subject_id=varargin{1};
        pipeline_id=varargin{2};
        opts=varargin(3:end);
        if rem(numel(opts),2), design_id=opts{1}; opts=opts(2:end);
        else design_id='';
        end
        project_id=regexprep(subject_id,'\d+.*$','');
        dataset=fullfile(OUTPUT_FOLDER,project_id,'derivatives','FL',pipeline_id,sprintf('sub-%s',subject_id));
        opts={'QA.PLOT',dataset,...
            opts{:}};
        fprintf('Project id %s\n',project_id);
        fprintf('Subject id %s\n',subject_id);
        fprintf('Experiment folder %s\n',dataset);
        fprintf('FL_INTERNAL command syntax:\n'); for n=1:numel(opts), fprintf('   '); disp(opts{n}); end
        assert(conn_existfile(dataset,1),'file %s not found',dataset);
        fl_internal(opts{:});
                
    case 'list'
        if isempty(varargin)
            project_id=conn_dir(fullfile(OUTPUT_FOLDER,'*'),'-cell','-sort','-R','-dir');
            [nill,project_id]=cellfun(@fileparts,project_id,'uni',0);
            project_id=project_id(cellfun('length',regexprep(project_id,'^\.+',''))>0);
            if isempty(project_id), fprintf('no experiments created\n')
            else
                if ~nargout, disp(char(project_id));
                else varargout={project_id};
                end
            end
        else 
            subject_id=varargin{1};
            project_id=regexprep(subject_id,'[\d\*]+.*$','');
            if isequal(project_id,subject_id) % exp
                if numel(varargin)>1 % exp pipeline_id
                    pipeline_id=varargin{2};
                else % exp
                    subject_info=fullfile(OUTPUT_FOLDER,project_id,'config','FL','IMPORT',conn_prepend('',project_id,'*.cfg'));
                    f2=conn_dir(subject_info,'-cell');
                    subject_id=regexprep(f2,'^.*[\\\/]config[\\\/]FL[\\\/]IMPORT[\\\/]',''); %[nill,subject_id]=cellfun(@fileparts,f2,'uni',0);
                    subject_id=conn_prepend('',subject_id,'');
                    if isempty(subject_id), fprintf('no subjects defined (no %s files found)\n',subject_info)
                    else
                        if ~nargout, disp(char(subject_id));
                        else varargout={subject_id};
                        end
                    end
                end
            else % subj
                if numel(varargin)>2, % subj_id pipeline_id model_id
                    pipeline_id=varargin{2};
                    model_id=varargin{3};
                    files=fullfile(OUTPUT_FOLDER,project_id,'derivatives','FL',pipeline_id,sprintf('sub-%s',subject_id),'results','firstlevel',model_id,'SPM.mat');
                    if ~conn_existfile(files)
                        fprintf('This first-level analysis has not been run yet');
                    else
                        conn_loadmatfile(files,'SPM');
                        SPMcolnames=SPM.xX.name;
                        S=regexp(SPMcolnames,'^Sn\((\d+)\)\s+(.*?)(\*bf\(1\))?$','tokens','once');
                        SPMidxvalid=find(cellfun(@(s)isequal(size(s),[1,3]),S));
                        S=cat(1,S{SPMidxvalid});
                        SPMconditions=unique(regexprep(S(:,2),'\^1$',''));
                        disp(char(SPMconditions)); % effect-names
                        disp(' ');
                        if ~nargout, disp(char({SPM.xCon.name})); % contrast-names
                        else varargout={{SPM.xCon.name}};
                        end
                    end
                elseif numel(varargin)>1, % subj_id pipeline_id
                    pipeline_id=varargin{2};
                    files=conn_dir(fullfile(OUTPUT_FOLDER,project_id,'derivatives','FL',pipeline_id,'SPM.mat'),'-cell','-sort');
                    [files,nill]=cellfun(@fileparts,files,'uni',0);
                    firstlevel_id=regexprep(files,'^.*[\\\/]derivatives[\\\/]FL[\\\/]',''); %[nill,firstlevel_id]=cellfun(@fileparts,files,'uni',0);
                    firstlevel_id=unique(firstlevel_id);
                    if isempty(firstlevel_id), fprintf('no first-level analyses run\n')
                    else
                        if ~nargout, disp(char(firstlevel_id));
                        else varargout={firstlevel_id};
                        end
                    end
                else % subj_id
                    pipeline_id=conn_dir(fullfile(OUTPUT_FOLDER,project_id,'derivatives','FL','*'),'-cell','-sort','-R','-dir');
                    [nill,pipeline_id]=cellfun(@fileparts,pipeline_id,'uni',0);
                    pipeline_id=pipeline_id(cellfun('length',regexprep(pipeline_id,'^\.+',''))>0);
                    if isempty(pipeline_id), fprintf('no preprocessing pipelines run\n')
                    else
                        if ~nargout, disp(char(pipeline_id));
                        else varargout={pipeline_id};
                        end
                    end
                end
            end
        end
        
    otherwise
        disp(sprintf('unrecognized option %s',STEPS));
end
end



