function conn_x=conn_updatefolders(conn_x,doextended)

% backward compatibility check of CONN_x structure
if nargin<2||isempty(doextended), doextended=true; end
if nargin<1||isempty(conn_x),
    global CONN_x;
    if ~isfield(CONN_x,'opt'), CONN_x.opt=[]; end
    if ~isfield(CONN_x.opt,'fmt1'), CONN_x.opt.fmt1='%d'; end
    if ~isfield(CONN_x,'pobj'), CONN_x.pobj=conn_projectmanager('null'); end
    if ~isfield(CONN_x,'ispending'), CONN_x.ispending=0; end
    [path,name,ext]=fileparts(CONN_x.filename);if isempty(path),path=pwd;end
    if ~isfield(CONN_x,'folders'),CONN_x.folders=struct('rois',[],'data',[],'preprocessing',[],'bids',[],'qa',[],'bookmarks',[],'firstlevel',[],'firstlevel_vv','firstlevel_dyn','secondlevel',[]);end
    if 1
        CONN_x.folders.rois=fullfile(fileparts(which(mfilename)),'rois');
        CONN_x.folders.data=fullfile(path,name,'data'); 
        CONN_x.folders.bids=fullfile(path,name,'data','BIDS'); 
        CONN_x.folders.qa=fullfile(path,name,'results','qa'); 
        CONN_x.folders.methods=fullfile(path,name,'results','methods'); 
        CONN_x.folders.bookmarks=fullfile(path,name,'results','bookmarks'); 
        CONN_x.folders.preprocessing=fullfile(path,name,'results','preprocessing'); 
        CONN_x.folders.firstlevel=fullfile(path,name,'results','firstlevel'); 
        CONN_x.folders.secondlevel=fullfile(path,name,'results','secondlevel');
        if CONN_x.pobj.holdsdata&&(~isfield(CONN_x,'isready')||numel(CONN_x.isready)<2||CONN_x.isready(2))
            if ~conn_existfile(CONN_x.folders.data,true), conn_fileutils('mkdir',path,name);conn_fileutils('mkdir',fullfile(path,name),'data'); end
            if ~conn_existfile(CONN_x.folders.bids,true), conn_fileutils('mkdir',path,name);conn_fileutils('mkdir',fullfile(path,name),'data'); conn_fileutils('mkdir',fullfile(path,name,'data'),'BIDS'); end
            if ~conn_existfile(CONN_x.folders.qa,true), conn_fileutils('mkdir',path,name);conn_fileutils('mkdir',fullfile(path,name),'results'); conn_fileutils('mkdir',fullfile(path,name,'results'),'qa'); end
            if ~conn_existfile(CONN_x.folders.methods,true), conn_fileutils('mkdir',path,name);conn_fileutils('mkdir',fullfile(path,name),'results'); conn_fileutils('mkdir',fullfile(path,name,'results'),'methods'); end
            if ~conn_existfile(CONN_x.folders.bookmarks,true), conn_fileutils('mkdir',path,name);conn_fileutils('mkdir',fullfile(path,name),'results'); conn_fileutils('mkdir',fullfile(path,name,'results'),'bookmarks'); end
            if ~conn_existfile(CONN_x.folders.preprocessing,true), conn_fileutils('mkdir',path,name);conn_fileutils('mkdir',fullfile(path,name),'results'); conn_fileutils('mkdir',fullfile(path,name,'results'),'preprocessing'); end
            if ~conn_existfile(CONN_x.folders.firstlevel,true), conn_fileutils('mkdir',path,name);conn_fileutils('mkdir',fullfile(path,name),'results'); conn_fileutils('mkdir',fullfile(path,name,'results'),'firstlevel'); end
            if ~conn_existfile(CONN_x.folders.secondlevel,true), conn_fileutils('mkdir',path,name);conn_fileutils('mkdir',fullfile(path,name),'results'); conn_fileutils('mkdir',fullfile(path,name,'results'),'secondlevel'); end
        end
        if conn('checkver','15.h',CONN_x.ver)
            if isfield(CONN_x,'vvAnalyses')&&~isfield(CONN_x.vvAnalyses,'name'), 
                [nill,CONN_x.vvAnalyses.name]=fileparts(CONN_x.folders.firstlevel_vv); 
            end
            if isfield(CONN_x,'dynAnalyses')&&~isfield(CONN_x.dynAnalyses,'name'), 
                [nill,CONN_x.dynAnalyses.name]=fileparts(CONN_x.folders.firstlevel_dyn); 
            end
            CONN_x.folders.firstlevel_vv=fullfile(path,name,'results','firstlevel');
            CONN_x.folders.firstlevel_dyn=fullfile(path,name,'results','firstlevel');
        else
            if isfield(CONN_x,'vvAnalyses')&&~isfield(CONN_x.vvAnalyses,'name'), CONN_x.vvAnalyses.name=''; end
            if isfield(CONN_x,'dynAnalyses')&&~isfield(CONN_x.dynAnalyses,'name'), CONN_x.dynAnalyses.name=''; end
            CONN_x.folders.firstlevel_vv=CONN_x.folders.firstlevel;
            CONN_x.folders.firstlevel_dyn=CONN_x.folders.preprocessing;
        end
    end
    if ~doextended, return; end
    if ~isfield(CONN_x,'isready'), CONN_x.isready=[1 1 1 1]; end
    if ~isfield(CONN_x.Setup,'steps'), CONN_x.Setup.steps=[1,1,0,0]; end 
    if numel(CONN_x.Setup.steps)~=4, CONN_x.Setup.steps=[CONN_x.Setup.steps(1:min(numel(CONN_x.Setup.steps),4)) zeros(1,max(0,4-numel(CONN_x.Setup.steps)))]; end
    if ~isfield(CONN_x.Setup,'spatialresolution'), CONN_x.Setup.spatialresolution=2; end    
    if ~isfield(CONN_x.Setup,'outputfiles'), CONN_x.Setup.outputfiles=[0,0,0,0,0,0]; end
    if numel(CONN_x.Setup.outputfiles)<6, CONN_x.Setup.outputfiles=[CONN_x.Setup.outputfiles,zeros(1,6-numel(CONN_x.Setup.outputfiles))]; end
    if ~isfield(CONN_x.Setup,'analysismask'), CONN_x.Setup.analysismask=1; end    
    if ~isfield(CONN_x.Setup,'analysisunits'), CONN_x.Setup.analysisunits=1; end    
    if ~isfield(CONN_x.Setup,'secondlevelanalyses'), CONN_x.Setup.secondlevelanalyses=1; end    
    if ~isfield(CONN_x.Setup,'explicitmask'), CONN_x.Setup.explicitmask=fullfile(fileparts(which(mfilename)),'utils','surf','mask.volume.brainmask.nii'); end
    if ~isempty(CONN_x.Setup.explicitmask)&&ischar(CONN_x.Setup.explicitmask)
        try, CONN_x.Setup.explicitmask=conn_file(CONN_x.Setup.explicitmask);
        catch, CONN_x.Setup.explicitmask={CONN_x.Setup.explicitmask, [], []};
        end
    end
    if ~isfield(CONN_x.Setup,'roifunctional')&&~isfield(CONN_x.Setup,'secondarydataset'), 
        if ~isfield(CONN_x.Setup,'roiextract'), CONN_x.Setup.roiextract=2; end
        if ~isfield(CONN_x.Setup,'roiextract_rule'), CONN_x.Setup.roiextract_rule={}; end
        if ~isfield(CONN_x.Setup,'roiextract_functional'), CONN_x.Setup.roiextract_functional={}; end
        CONN_x.Setup.roifunctional.roiextract=CONN_x.Setup.roiextract; 
        CONN_x.Setup.roifunctional.roiextract_rule=CONN_x.Setup.roiextract_rule;
        CONN_x.Setup.roifunctional.roiextract_functional=CONN_x.Setup.roiextract_functional;
        CONN_x.Setup=rmfield(CONN_x.Setup,{'roiextract','roiextract_rule','roiextract_functional'});
    end
    if ~isfield(CONN_x.Setup,'secondarydataset')&&isfield(CONN_x.Setup,'roifunctional'), 
        [CONN_x.Setup.roifunctional.functionals_type]=deal(CONN_x.Setup.roifunctional.roiextract);
        [CONN_x.Setup.roifunctional.functionals_explicit]=deal(CONN_x.Setup.roifunctional.roiextract_functional);
        [CONN_x.Setup.roifunctional.functionals_rule]=deal(CONN_x.Setup.roifunctional.roiextract_rule);
        CONN_x.Setup.secondarydataset=rmfield(CONN_x.Setup.roifunctional,{'roiextract','roiextract_functional','roiextract_rule'});
        CONN_x.Setup=rmfield(CONN_x.Setup,'roifunctional');
    end
    if ~isfield(CONN_x.Setup.secondarydataset,'label'), 
        [CONN_x.Setup.secondarydataset.label]=deal(''); 
        [CONN_x.Setup.secondarydataset([CONN_x.Setup.secondarydataset.functionals_type]==2).label]=deal('non-smoothed data'); 
    end
    if ~isfield(CONN_x.Setup,'dicom'), CONN_x.Setup.dicom=repmat({{[],[],[]}},1,CONN_x.Setup.nsubjects); end
    if ~isfield(CONN_x.Setup,'bids'), CONN_x.Setup.bids={[],[],[]}; end
    if ~isfield(CONN_x.Setup,'unwarp_functional'), CONN_x.Setup.unwarp_functional={}; end
    if ~isfield(CONN_x.Setup,'coregsource_functional'), CONN_x.Setup.coregsource_functional={}; end
    if ~isfield(CONN_x.Setup,'erosion')
        if ~isfield(CONN_x.Setup,'cwthreshold'), cwthreshold=[.5 1]; 
        else cwthreshold=CONN_x.Setup.cwthreshold; 
        end
        ithr=[3,1,2];
        if isfield(CONN_x.Setup,'cwthreshold')&&numel(cwthreshold)>0, THRall=cwthreshold(1:2:end);
        else THRall=.5;  % default threshold for white/CSF probability maps before erosion
        end
        if isfield(CONN_x.Setup,'cwthreshold')&&numel(cwthreshold)>1, ERODEall=cwthreshold(2:2:end);
        else ERODEall=1; % default erosing level for white/csf masks (voxels) (set to 0 for no erosion)
        end
        for nroi=1:3 % fills-out defaults in a back-compatible manner if these values are not defined
            if nroi==1,
                if numel(ERODEall)>=3, ERODE=ERODEall(3); else ERODE=0; end
                if numel(THRall)>=3, THR=THRall(3); else THR=.5; end
            else
                ERODE=ERODEall(max(1,min(numel(ERODEall),nroi-1)));
                THR=THRall(max(1,min(numel(THRall),nroi-1)));
            end
            cwthreshold(2*ithr(nroi)-1)=THR;
            cwthreshold(2*ithr(nroi))=ERODE;
        end
        cwthreshold=reshape(cwthreshold,2,3);
        CONN_x.Setup.erosion=struct('binary_threshold',cwthreshold(1,ithr),'erosion_steps',cwthreshold(2,ithr),'erosion_neighb',[1 1 1]);
    end
    if ~isfield(CONN_x.Setup.erosion,'binary_threshold_type'), CONN_x.Setup.erosion.binary_threshold_type=[1 1 1]; end
    if ~isfield(CONN_x.Setup.erosion,'exclude_grey_matter'), CONN_x.Setup.erosion.exclude_grey_matter=[nan nan nan]; end
    if ~isfield(CONN_x.Setup,'acquisitiontype'), 
        if isfield(CONN_x.Setup,'sparseacquisition'), CONN_x.Setup.acquisitiontype=1+(CONN_x.Setup.sparseacquisition>0);
        else CONN_x.Setup.acquisitiontype=1; end
    end    
    if ~isfield(CONN_x.Setup.conditions,'model'), CONN_x.Setup.conditions.model=cell(1,numel(CONN_x.Setup.conditions.names)-1); end
    if ~isfield(CONN_x.Setup.conditions,'param'), CONN_x.Setup.conditions.param=zeros(1,numel(CONN_x.Setup.conditions.names)-1); end
    if ~isfield(CONN_x.Setup.conditions,'filter'), CONN_x.Setup.conditions.filter=cell(1,numel(CONN_x.Setup.conditions.names)-1); end
    if ~isfield(CONN_x.Setup.conditions,'allnames'),CONN_x.Setup.conditions.allnames=CONN_x.Setup.conditions.names(1:end-1); end
    if ~isfield(CONN_x.Setup.conditions,'missingdata'),CONN_x.Setup.conditions.missingdata=0; end
    if ~isfield(CONN_x.Setup.rois,'unsmoothedvolumes'), CONN_x.Setup.rois.unsmoothedvolumes=ones(1,numel(CONN_x.Setup.rois.dimensions)); end
    if ~isfield(CONN_x.Setup.rois,'sessionspecific'),
        CONN_x.Setup.rois.sessionspecific=zeros(1,numel(CONN_x.Setup.rois.names)-1);
        for nsub=1:CONN_x.Setup.nsubjects
            for nroi=1:numel(CONN_x.Setup.rois.names)-1
                CONN_x.Setup.rois.files{nsub}{nroi}=repmat({CONN_x.Setup.rois.files{nsub}{nroi}},[1,CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsub))]);
            end
        end
    end
    if ~isfield(CONN_x.Setup.rois,'subjectspecific'),
        CONN_x.Setup.rois.subjectspecific=zeros(1,numel(CONN_x.Setup.rois.names)-1);
        for nroi=1:numel(CONN_x.Setup.rois.names)-1
            temp=CONN_x.Setup.rois.files{1}{nroi}{1};
            for nsub=1:CONN_x.Setup.nsubjects
                for nses=1:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsub))
                    if ~isequal(CONN_x.Setup.rois.files{nsub}{nroi}{nses},temp), CONN_x.Setup.rois.subjectspecific(nroi)=1; break; end
                end
            end
        end
    end
    if ~isfield(CONN_x.Setup.rois,'weighted'), CONN_x.Setup.rois.weighted=cellfun(@(x)isequal(x,0),CONN_x.Setup.rois.dimensions); end
    if ~isfield(CONN_x.Setup,'structural_sessionspecific'),
        CONN_x.Setup.structural_sessionspecific=0;
        for nsub=1:CONN_x.Setup.nsubjects
            CONN_x.Setup.structural{nsub}=repmat({CONN_x.Setup.structural{nsub}},[1,CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsub))]);
        end
    end
    if ~isfield(CONN_x.Setup.l2covariates,'descrip'), CONN_x.Setup.l2covariates.descrip=cell(1,numel(CONN_x.Setup.l2covariates.names)-1); end
    if ~conn('checkver','18.a',CONN_x.ver), CONN_x.Setup.l1covariates.names=regexprep(CONN_x.Setup.l1covariates.names,'^QA_','QC_'); CONN_x.Setup.l2covariates.names=regexprep(CONN_x.Setup.l2covariates.names,'^QA_','QC_'); end
    if ~isfield(CONN_x.Preproc,'regbp'), CONN_x.Preproc.regbp=1; end
    if ~isfield(CONN_x.Preproc,'despiking'), CONN_x.Preproc.despiking=0; end
    if ~isfield(CONN_x.Preproc,'detrending'), CONN_x.Preproc.detrending=0; end
    if isfield(CONN_x.Preproc,'variables')&&~isfield(CONN_x.Preproc.variables,'power'), CONN_x.Preproc.variables.power=repmat({1},1,length(CONN_x.Preproc.variables.names)); end
    if isfield(CONN_x.Preproc,'confounds')&&~isfield(CONN_x.Preproc.confounds,'power'), CONN_x.Preproc.confounds.power=repmat({1},1,length(CONN_x.Preproc.confounds.names)); end
    if isfield(CONN_x.Preproc,'variables')&&~isfield(CONN_x.Preproc.variables,'filter'), CONN_x.Preproc.variables.filter=repmat({0},1,length(CONN_x.Preproc.variables.names)); end
    if isfield(CONN_x.Preproc,'confounds')&&~isfield(CONN_x.Preproc.confounds,'filter'), CONN_x.Preproc.confounds.filter=repmat({0},1,length(CONN_x.Preproc.confounds.names)); end
    if ~isfield(CONN_x,'Analysis'), CONN_x.Analysis=1; end
    if ~isempty(CONN_x.Analyses)&&~isfield(CONN_x,'Analysis_variables')
        CONN_x.Analysis_variables=CONN_x.Analyses(1).variables;
        CONN_x.Analyses=rmfield(CONN_x.Analyses,'variables');
        if ~isempty(CONN_x.Analysis_variables)&&isstruct(CONN_x.Analysis_variables)&&~isfield(CONN_x.Analysis_variables,'fbands'), CONN_x.Analysis_variables.fbands=repmat({1},size(CONN_x.Analysis_variables.names)); end
    end
    for ianalysis=1:numel(CONN_x.Analyses)
        if ~isfield(CONN_x.Analyses(ianalysis),'name'),CONN_x.Analyses(ianalysis).name=['SBC_',num2str(ianalysis,'%02d')]; end
        if ~isfield(CONN_x.Analyses(ianalysis),'sourcenames'),CONN_x.Analyses(ianalysis).sourcenames={}; end
        if ~isfield(CONN_x.Analyses(ianalysis),'modulation')||isempty(CONN_x.Analyses(ianalysis).modulation),CONN_x.Analyses(ianalysis).modulation=0; end
        if ~isfield(CONN_x.Analyses(ianalysis),'measure')||isempty(CONN_x.Analyses(ianalysis).measure),CONN_x.Analyses(ianalysis).measure=1; end
        if ~isfield(CONN_x.Analyses(ianalysis),'weight')||isempty(CONN_x.Analyses(ianalysis).weight),CONN_x.Analyses(ianalysis).weight=2; end
        if ~isfield(CONN_x.Analyses(ianalysis),'type')||isempty(CONN_x.Analyses(ianalysis).type),CONN_x.Analyses(ianalysis).type=3; end
        if ~isfield(CONN_x.Analyses(ianalysis),'conditions'),CONN_x.Analyses(ianalysis).conditions={''}; end
        %if ~isempty(CONN_x.Analyses(ianalysis).variables)&&isstruct(CONN_x.Analyses(ianalysis).variables)&&~isfield(CONN_x.Analyses(ianalysis).variables,'fbands'), CONN_x.Analyses(ianalysis).variables.fbands=repmat({1},size(CONN_x.Analyses(ianalysis).variables.names)); end
        if ~isempty(CONN_x.Analyses(ianalysis).regressors)&&isstruct(CONN_x.Analyses(ianalysis).regressors)&&~isfield(CONN_x.Analyses(ianalysis).regressors,'fbands'), CONN_x.Analyses(ianalysis).regressors.fbands=repmat({1},size(CONN_x.Analyses(ianalysis).regressors.names)); end
    end    
    if ~isfield(CONN_x,'vvAnalyses'), CONN_x.vvAnalyses=struct('name','','measurenames',{{}},'variables', conn_v2v('measures'),'regressors', conn_v2v('measures'),'measures',{{}},'mask',[],'options',''); end
    if ~isfield(CONN_x,'vvAnalysis'), CONN_x.vvAnalysis=1; end
    for ianalysis=1:numel(CONN_x.vvAnalyses)
        if ~isfield(CONN_x.vvAnalyses(ianalysis),'name')||isempty(CONN_x.vvAnalyses(ianalysis).name), CONN_x.vvAnalyses(ianalysis).name=''; end
        if ~isfield(CONN_x.vvAnalyses(ianalysis),'mask'), CONN_x.vvAnalyses(ianalysis).mask=[]; end
        if ~isfield(CONN_x.vvAnalyses(ianalysis),'options')||isempty(CONN_x.vvAnalyses(ianalysis).options), CONN_x.vvAnalyses(ianalysis).options=''; end
        if ~isfield(CONN_x.vvAnalyses(ianalysis),'measures'), CONN_x.vvAnalyses(ianalysis).measures={}; end
        if isfield(CONN_x.vvAnalyses(ianalysis),'variables')&&~isfield(CONN_x.vvAnalyses(ianalysis).variables,'norm'), CONN_x.vvAnalyses(ianalysis).variables.norm=repmat({1},1,numel(CONN_x.vvAnalyses(ianalysis).variables.names)); end
        if isfield(CONN_x.vvAnalyses(ianalysis),'regressors')&&~isfield(CONN_x.vvAnalyses(ianalysis).regressors,'norm'), CONN_x.vvAnalyses(ianalysis).regressors.norm=repmat({1},1,numel(CONN_x.vvAnalyses(ianalysis).regressors.names)); end
        if isfield(CONN_x.vvAnalyses(ianalysis),'variables')&&~isfield(CONN_x.vvAnalyses(ianalysis).variables,'alt_names'), CONN_x.vvAnalyses(ianalysis).variables.alt_names=repmat({{}},1,numel(CONN_x.vvAnalyses(ianalysis).variables.names)); end
        if isfield(CONN_x.vvAnalyses(ianalysis),'regressors')&&~isfield(CONN_x.vvAnalyses(ianalysis).regressors,'alt_names'), CONN_x.vvAnalyses(ianalysis).regressors.alt_names=repmat({{}},1,numel(CONN_x.vvAnalyses(ianalysis).regressors.names)); end
        %if isfield(CONN_x.vvAnalyses(ianalysis),'variables')&&~isfield(CONN_x.vvAnalyses(ianalysis).variables,'mask'), CONN_x.vvAnalyses(ianalysis).variables.mask=repmat({{}},1,numel(CONN_x.vvAnalyses(ianalysis).variables.names)); end
        %if isfield(CONN_x.vvAnalyses(ianalysis),'regressors')&&~isfield(CONN_x.vvAnalyses(ianalysis).regressors,'mask'), CONN_x.vvAnalyses(ianalysis).regressors.mask=repmat({{}},1,numel(CONN_x.vvAnalyses(ianalysis).regressors.names)); end
    end
    if ~isfield(CONN_x,'dynAnalyses'), CONN_x.dynAnalyses=struct('name','','regressors', struct('names',{{}}),'variables', struct('names',{{}}),'Ncomponents',[20],'condition',[1],'window',10,'output',[1 1 0]); end
    if ~isfield(CONN_x,'dynAnalysis'), CONN_x.dynAnalysis=1; end
    if isfield(CONN_x.dynAnalyses,'analyses'), CONN_x.dynAnalyses=rmfield(CONN_x.dynAnalyses,'analyses'); end
    for ianalysis=1:numel(CONN_x.dynAnalyses)
        if ~isfield(CONN_x.dynAnalyses(ianalysis),'name'), CONN_x.dynAnalyses(ianalysis).name=''; end
        if ~isfield(CONN_x.dynAnalyses(ianalysis),'condition')||isempty(CONN_x.dynAnalyses(ianalysis).condition), CONN_x.dynAnalyses(ianalysis).condition=[]; end
        if ~isfield(CONN_x.dynAnalyses(ianalysis),'window')||isempty(CONN_x.dynAnalyses(ianalysis).window),
            if isfield(CONN_x.dynAnalyses(ianalysis),'filter')&&~isempty(CONN_x.dynAnalyses(ianalysis).filter), CONN_x.dynAnalyses(ianalysis).window=round(1/CONN_x.dynAnalyses(ianalysis).filter);
            else CONN_x.dynAnalyses(ianalysis).window=30;
            end
        end
    end
    if ~isfield(CONN_x,'Results')||~isfield(CONN_x.Results,'saved')||isempty(CONN_x.Results.saved), CONN_x.Results.saved=struct('names',{{}},'labels',{{}},'descrip',{{}},'nsubjecteffects',{{}},'csubjecteffects',{{}},'nconditions',{{}},'cconditions',{{}}); end
    if ~isfield(CONN_x.Results.saved,'names'), CONN_x.Results.saved.names={}; end
    if ~isfield(CONN_x,'isready')||numel(CONN_x.isready)<2||CONN_x.isready(2)
        for ianalysis=1:length(CONN_x.Analyses),
            if ~conn_existfile(fullfile(CONN_x.folders.firstlevel,CONN_x.Analyses(ianalysis).name),true), conn_fileutils('mkdir',CONN_x.folders.firstlevel,CONN_x.Analyses(ianalysis).name); end;
            filesourcenames=fullfile(CONN_x.folders.firstlevel,CONN_x.Analyses(ianalysis).name,'_list_sources.mat');
            filesourcenames=conn_projectmanager('projectfile',filesourcenames,CONN_x.pobj,'.mat');
            if conn_existfile(filesourcenames),
                sourcenames={}; conn_loadmatfile(filesourcenames,'sourcenames');
                if numel(sourcenames)>=numel(CONN_x.Analyses(ianalysis).sourcenames), %~isempty(sourcenames)
                    if ~isequal(CONN_x.Analyses(ianalysis).sourcenames,sourcenames), fprintf('warning: mismatch list_sources info. Using automatic recovery from %s\n',filesourcenames); end
                    CONN_x.Analyses(ianalysis).sourcenames=sourcenames;
                end
            end
        end
        for ianalysis=1:length(CONN_x.vvAnalyses),
            if ~conn_existfile(fullfile(CONN_x.folders.firstlevel_vv,CONN_x.vvAnalyses(ianalysis).name),true), conn_fileutils('mkdir',CONN_x.folders.firstlevel_vv,CONN_x.vvAnalyses(ianalysis).name); end;
            filemeasurenames=fullfile(CONN_x.folders.firstlevel_vv,CONN_x.vvAnalyses(ianalysis).name,'_list_measures.mat');
            filemeasurenames=conn_projectmanager('projectfile',filemeasurenames,CONN_x.pobj,'.mat');
            if conn_existfile(filemeasurenames),
                measurenames={}; conn_loadmatfile(filemeasurenames,'measurenames');
                if numel(measurenames)>=numel(CONN_x.vvAnalyses(ianalysis).measurenames), %~isempty(measurenames)
                    if ~isequal(CONN_x.vvAnalyses(ianalysis).measurenames,measurenames), fprintf('warning: mismatch list_measures info. Using automatic recovery from %s\n',filemeasurenames); end
                    CONN_x.vvAnalyses(ianalysis).measurenames=measurenames;
                end
            end
        end
        for ianalysis=1:length(CONN_x.dynAnalyses),
            if ~conn_existfile(fullfile(CONN_x.folders.firstlevel_dyn,CONN_x.dynAnalyses(ianalysis).name),true), conn_fileutils('mkdir',CONN_x.folders.firstlevel_dyn,CONN_x.dynAnalyses(ianalysis).name); end;
        end
        fileconditionnames=fullfile(CONN_x.folders.preprocessing,'_list_conditions.mat');
        fileconditionnames=conn_projectmanager('projectfile',fileconditionnames,CONN_x.pobj,'.mat');
        if conn_existfile(fileconditionnames),
            allnames={}; conn_loadmatfile(fileconditionnames,'allnames');
            if numel(allnames)>=numel(CONN_x.Setup.conditions.allnames), %~isempty(allnames)
                if ~isequal(CONN_x.Setup.conditions.allnames,allnames), fprintf('warning: mismatch list_conditions info. Using automatic recovery from %s\n',fileconditionnames); end
                CONN_x.Setup.conditions.allnames=allnames;
            end
        end
        fileresultsnames=fullfile(CONN_x.folders.secondlevel,'_list_results.mat');
        if conn_existfile(fileresultsnames),
            results=struct; conn_loadmatfile(fileresultsnames,'results');
            if ~isequal(CONN_x.Results.saved,results), fprintf('warning: mismatch list_results info. Using automatic recovery from %s\n',fileresultsnames); end
            if ~isempty(results)
                CONN_x.Results.saved=results;
            end
        end
        if ~isfield(CONN_x.Results.saved,'labels'), CONN_x.Results.saved.labels=repmat({''},size(CONN_x.Results.saved.names)); end
        if ~isfield(CONN_x.Results.saved,'nsubjecteffects'), CONN_x.Results.saved.nsubjecteffects=cell(size(CONN_x.Results.saved.names)); end
        if ~isfield(CONN_x.Results.saved,'csubjecteffects'), CONN_x.Results.saved.csubjecteffects=cell(size(CONN_x.Results.saved.names)); end
        if ~isfield(CONN_x.Results.saved,'nconditions'), CONN_x.Results.saved.nconditions=cell(size(CONN_x.Results.saved.names)); end
        if ~isfield(CONN_x.Results.saved,'cconditions'), CONN_x.Results.saved.cconditions=cell(size(CONN_x.Results.saved.names)); end
        if ~isfield(CONN_x.Results.saved,'descrip'), CONN_x.Results.saved.descrip=repmat({''},size(CONN_x.Results.saved.names)); end
    end

else
    if ~isfield(conn_x,'opt'), conn_x.opt=[]; end
    if ~isfield(conn_x.opt,'fmt1'), conn_x.opt.fmt1='%d'; end
    if ~isfield(conn_x,'pobj'), conn_x.pobj=conn_projectmanager('null'); end
    if ~isfield(conn_x,'ispending'), conn_x.ispending=0; end
    [path,name,ext]=fileparts(conn_x.filename);if isempty(path),path=pwd;end
    if ~isfield(conn_x,'folders'),conn_x.folders=struct('rois',[],'data',[],'preprocessing',[],'bids',[],'qa',[],'bookmarks',[],'firstlevel',[],'firstlevel_vv','firstlevel_dyn','secondlevel',[]);end
    if 1
        conn_x.folders.rois=fullfile(fileparts(which(mfilename)),'rois');
        conn_x.folders.data=fullfile(path,name,'data'); 
        conn_x.folders.bids=fullfile(path,name,'data','BIDS'); 
        conn_x.folders.qa=fullfile(path,name,'results','qa'); 
        conn_x.folders.methods=fullfile(path,name,'results','methods'); 
        conn_x.folders.bookmarks=fullfile(path,name,'results','bookmarks'); 
        conn_x.folders.preprocessing=fullfile(path,name,'results','preprocessing'); 
        conn_x.folders.firstlevel=fullfile(path,name,'results','firstlevel'); 
        conn_x.folders.secondlevel=fullfile(path,name,'results','secondlevel');
        if conn_x.pobj.holdsdata&&(~isfield(conn_x,'isready')||numel(conn_x.isready)<2||conn_x.isready(2))
            if ~conn_existfile(conn_x.folders.data,true), conn_fileutils('mkdir',path,name);conn_fileutils('mkdir',fullfile(path,name),'data'); end
            if ~conn_existfile(conn_x.folders.bids,true), conn_fileutils('mkdir',path,name);conn_fileutils('mkdir',fullfile(path,name),'data'); conn_fileutils('mkdir',fullfile(path,name,'data'),'BIDS'); end
            if ~conn_existfile(conn_x.folders.qa,true), conn_fileutils('mkdir',path,name);conn_fileutils('mkdir',fullfile(path,name),'results'); conn_fileutils('mkdir',fullfile(path,name,'results'),'qa'); end
            if ~conn_existfile(conn_x.folders.methods,true), conn_fileutils('mkdir',path,name);conn_fileutils('mkdir',fullfile(path,name),'results'); conn_fileutils('mkdir',fullfile(path,name,'results'),'methods'); end
            if ~conn_existfile(conn_x.folders.bookmarks,true), conn_fileutils('mkdir',path,name);conn_fileutils('mkdir',fullfile(path,name),'results'); conn_fileutils('mkdir',fullfile(path,name,'results'),'bookmarks'); end
            if ~conn_existfile(conn_x.folders.preprocessing,true), conn_fileutils('mkdir',path,name);conn_fileutils('mkdir',fullfile(path,name),'results'); conn_fileutils('mkdir',fullfile(path,name,'results'),'preprocessing'); end
            if ~conn_existfile(conn_x.folders.firstlevel,true), conn_fileutils('mkdir',path,name);conn_fileutils('mkdir',fullfile(path,name),'results'); conn_fileutils('mkdir',fullfile(path,name,'results'),'firstlevel'); end
            if ~conn_existfile(conn_x.folders.secondlevel,true), conn_fileutils('mkdir',path,name);conn_fileutils('mkdir',fullfile(path,name),'results'); conn_fileutils('mkdir',fullfile(path,name,'results'),'secondlevel'); end
        end
        if conn('checkver','15.h',conn_x.ver)
            if isfield(conn_x,'vvAnalyses')&&~isfield(conn_x.vvAnalyses,'name'), [nill,conn_x.vvAnalyses.name]=fileparts(conn_x.folders.firstlevel_vv); end
            if isfield(conn_x,'dynAnalyses')&&~isfield(conn_x.dynAnalyses,'name'), [nill,conn_x.dynAnalyses.name]=fileparts(conn_x.folders.firstlevel_dyn); end
            conn_x.folders.firstlevel_vv=fullfile(path,name,'results','firstlevel'); 
            conn_x.folders.firstlevel_dyn=fullfile(path,name,'results','firstlevel');
        else
            conn_x.folders.firstlevel_vv=conn_x.folders.firstlevel;
            conn_x.folders.firstlevel_dyn=conn_x.folders.preprocessing;
        end
    end
    if ~doextended, return; end
    if ~isfield(conn_x,'isready'), conn_x.isready=[1 1 1 1]; end
    if ~isfield(conn_x.Setup,'steps'), conn_x.Setup.steps=[1,1,0,0]; end
    if numel(conn_x.Setup.steps)~=4, conn_x.Setup.steps=[conn_x.Setup.steps(1:min(numel(conn_x.Setup.steps),4)) zeros(1,max(0,4-numel(conn_x.Setup.steps)))]; end
    if ~isfield(conn_x.Setup,'spatialresolution'), conn_x.Setup.spatialresolution=2; end
    if ~isfield(conn_x.Setup,'outputfiles'), conn_x.Setup.outputfiles=[0,0,0,0,0,0]; end
    if numel(conn_x.Setup.outputfiles)<6, conn_x.Setup.outputfiles=[conn_x.Setup.outputfiles,zeros(1,6-numel(conn_x.Setup.outputfiles))]; end
    if ~isfield(conn_x.Setup,'analysismask'), conn_x.Setup.analysismask=1; end    
    if ~isfield(conn_x.Setup,'analysisunits'), conn_x.Setup.analysisunits=1; end    
    if ~isfield(conn_x.Setup,'secondlevelanalyses'), conn_x.Setup.secondlevelanalyses=1; end    
    if ~isfield(conn_x.Setup,'explicitmask'), conn_x.Setup.explicitmask=fullfile(fileparts(which(mfilename)),'utils','surf','mask.volume.brainmask.nii'); end
    if ~isempty(conn_x.Setup.explicitmask)&&ischar(conn_x.Setup.explicitmask)
        try, conn_x.Setup.explicitmask=conn_file(conn_x.Setup.explicitmask);
        catch, conn_x.Setup.explicitmask={conn_x.Setup.explicitmask, [], []};
        end
    end
    if ~isfield(conn_x.Setup,'roifunctional')&&~isfield(conn_x.Setup,'secondarydataset'), 
        if ~isfield(conn_x.Setup,'roiextract'), conn_x.Setup.roiextract=2; end
        if ~isfield(conn_x.Setup,'roiextract_rule'), conn_x.Setup.roiextract_rule={}; end
        if ~isfield(conn_x.Setup,'roiextract_functional'), conn_x.Setup.roiextract_functional={}; end
        conn_x.Setup.roifunctional.roiextract=conn_x.Setup.roiextract; 
        conn_x.Setup.roifunctional.roiextract_rule=conn_x.Setup.roiextract_rule;
        conn_x.Setup.roifunctional.roiextract_functional=conn_x.Setup.roiextract_functional;
        conn_x.Setup=rmfield(conn_x.Setup,{'roiextract','roiextract_rule','roiextract_functional'});
    end
    if ~isfield(conn_x.Setup,'secondarydataset')&&isfield(conn_x.Setup,'roifunctional'), 
        [conn_x.Setup.roifunctional.functionals_type]=deal(conn_x.Setup.roifunctional.roiextract);
        [conn_x.Setup.roifunctional.functionals_explicit]=deal(conn_x.Setup.roifunctional.roiextract_functional);
        [conn_x.Setup.roifunctional.functionals_rule]=deal(conn_x.Setup.roifunctional.roiextract_rule);
        conn_x.Setup.secondarydataset=rmfield(conn_x.Setup.roifunctional,{'roiextract','roiextract_functional','roiextract_rule'});
        conn_x.Setup=rmfield(conn_x.Setup,'roifunctional');
    end
    if ~isfield(conn_x.Setup.secondarydataset,'label'), 
        [conn_x.Setup.secondarydataset.label]=deal(''); 
        [conn_x.Setup.secondarydataset([conn_x.Setup.secondarydataset.functionals_type]==2).label]=deal('non-smoothed data'); 
    end
    if ~isfield(conn_x.Setup,'dicom'), conn_x.Setup.dicom=repmat({{[],[],[]}},1,conn_x.Setup.nsubjects); end
    if ~isfield(conn_x.Setup,'bids'), conn_x.Setup.bids={[],[],[]}; end
    if ~isfield(conn_x.Setup,'unwarp_functional'), conn_x.Setup.unwarp_functional={}; end
    if ~isfield(conn_x.Setup,'coregsource_functional'), conn_x.Setup.coregsource_functional={}; end
    if ~isfield(conn_x.Setup,'erosion')
        if ~isfield(conn_x.Setup,'cwthreshold'), cwthreshold=[.5 1]; 
        else cwthreshold=conn_x.Setup.cwthreshold; 
        end
        ithr=[3,1,2];
        if isfield(conn_x.Setup,'cwthreshold')&&numel(cwthreshold)>0, THRall=cwthreshold(1:2:end);
        else THRall=.5;  % default threshold for white/CSF probability maps before erosion
        end
        if isfield(conn_x.Setup,'cwthreshold')&&numel(cwthreshold)>1, ERODEall=cwthreshold(2:2:end);
        else ERODEall=1; % default erosing level for white/csf masks (voxels) (set to 0 for no erosion)
        end
        for nroi=1:3 % fills-out defaults in a back-compatible manner if these values are not defined
            if nroi==1,
                if numel(ERODEall)>=3, ERODE=ERODEall(3); else ERODE=0; end
                if numel(THRall)>=3, THR=THRall(3); else THR=.5; end
            else
                ERODE=ERODEall(max(1,min(numel(ERODEall),nroi-1)));
                THR=THRall(max(1,min(numel(THRall),nroi-1)));
            end
            cwthreshold(2*ithr(nroi)-1)=THR;
            cwthreshold(2*ithr(nroi))=ERODE;
        end
        cwthreshold=reshape(cwthreshold,2,3);
        conn_x.Setup.erosion=struct('binary_threshold',cwthreshold(1,ithr),'erosion_steps',cwthreshold(2,ithr),'erosion_neighb',[1 1 1]);
    end
    if ~isfield(conn_x.Setup.erosion,'binary_threshold_type'), conn_x.Setup.erosion.binary_threshold_type=[1 1 1]; end
    if ~isfield(conn_x.Setup.erosion,'exclude_grey_matter'), conn_x.Setup.erosion.exclude_grey_matter=[nan nan nan]; end
    if ~isfield(conn_x.Setup,'acquisitiontype'), 
        if isfield(conn_x.Setup,'sparseacquisition'), conn_x.Setup.acquisitiontype=1+(conn_x.Setup.sparseacquisition>0);
        else conn_x.Setup.acquisitiontype=1; end
    end    
    if ~isfield(conn_x.Setup.conditions,'model'), conn_x.Setup.conditions.model=cell(1,numel(conn_x.Setup.conditions.names)-1); end
    if ~isfield(conn_x.Setup.conditions,'param'), conn_x.Setup.conditions.param=zeros(1,numel(conn_x.Setup.conditions.names)-1); end
    if ~isfield(conn_x.Setup.conditions,'filter'), conn_x.Setup.conditions.filter=cell(1,numel(conn_x.Setup.conditions.names)-1); end
    if ~isfield(conn_x.Setup.conditions,'allnames'),conn_x.Setup.conditions.allnames=conn_x.Setup.conditions.names(1:end-1); end
    if ~isfield(conn_x.Setup.conditions,'missingdata'),conn_x.Setup.conditions.missingdata=0; end
    if ~isfield(conn_x.Setup.rois,'unsmoothedvolumes'), conn_x.Setup.rois.unsmoothedvolumes=ones(1,numel(conn_x.Setup.rois.dimensions)); end
    if ~isfield(conn_x.Setup.rois,'sessionspecific'),
        conn_x.Setup.rois.sessionspecific=zeros(1,numel(conn_x.Setup.rois.names)-1);
        for nsub=1:conn_x.Setup.nsubjects
            for nroi=1:numel(conn_x.Setup.rois.names)-1
                conn_x.Setup.rois.files{nsub}{nroi}=repmat({conn_x.Setup.rois.files{nsub}{nroi}},[1,conn_x.Setup.nsessions(min(numel(conn_x.Setup.nsessions),nsub))]);
            end
        end
    end
    if ~isfield(conn_x.Setup.rois,'subjectspecific'),
        conn_x.Setup.rois.subjectspecific=zeros(1,numel(conn_x.Setup.rois.names)-1);
        for nroi=1:numel(conn_x.Setup.rois.names)-1
            temp=conn_x.Setup.rois.files{1}{nroi}{1};
            for nsub=1:conn_x.Setup.nsubjects
                for nses=1:conn_x.Setup.nsessions(min(numel(conn_x.Setup.nsessions),nsub))
                    if ~isequal(conn_x.Setup.rois.files{nsub}{nroi}{nses},temp), conn_x.Setup.rois.subjectspecific(nroi)=1; break; end
                end
            end
        end
    end
    if ~isfield(conn_x.Setup.rois,'weighted'), conn_x.Setup.rois.weighted=cellfun(@(x)isequal(x,0),conn_x.Setup.rois.dimensions); end
    if ~isfield(conn_x.Setup,'structural_sessionspecific'),
        conn_x.Setup.structural_sessionspecific=0;
        for nsub=1:conn_x.Setup.nsubjects
            conn_x.Setup.structural{nsub}=repmat({conn_x.Setup.structural{nsub}},[1,conn_x.Setup.nsessions(min(numel(conn_x.Setup.nsessions),nsub))]);
        end
    end
    if ~isfield(conn_x.Setup.l2covariates,'descrip'), conn_x.Setup.l2covariates.descrip=cell(1,numel(conn_x.Setup.l2covariates.names)-1); end
    if ~conn('checkver','18.a',conn_x.ver), conn_x.Setup.l1covariates.names=regexprep(conn_x.Setup.l1covariates.names,'^QA_','QC_'); conn_x.Setup.l2covariates.names=regexprep(conn_x.Setup.l2covariates.names,'^QA_','QC_'); end
    if ~isfield(conn_x.Preproc,'regbp'), conn_x.Preproc.regbp=1; end
    if ~isfield(conn_x.Preproc,'despiking'), conn_x.Preproc.despiking=0; end
    if ~isfield(conn_x.Preproc,'detrending'), conn_x.Preproc.detrending=0; end
    if isfield(conn_x.Preproc,'variables')&&~isfield(conn_x.Preproc.variables,'power'), conn_x.Preproc.variables.power=repmat({1},1,length(conn_x.Preproc.variables.names)); end
    if isfield(conn_x.Preproc,'confounds')&&~isfield(conn_x.Preproc.confounds,'power'), conn_x.Preproc.confounds.power=repmat({1},1,length(conn_x.Preproc.confounds.names)); end
    if isfield(conn_x.Preproc,'variables')&&~isfield(conn_x.Preproc.variables,'filter'), conn_x.Preproc.variables.filter=repmat({0},1,length(conn_x.Preproc.variables.names)); end
    if isfield(conn_x.Preproc,'confounds')&&~isfield(conn_x.Preproc.confounds,'filter'), conn_x.Preproc.confounds.filter=repmat({0},1,length(conn_x.Preproc.confounds.names)); end
    if ~isfield(conn_x,'Analysis'), conn_x.Analysis=1; end
    if ~isempty(conn_x.Analyses)&&~isfield(conn_x,'Analysis_variables')
        conn_x.Analysis_variables=conn_x.Analyses(1).variables;
        conn_x.Analyses=rmfield(conn_x.Analyses,'variables');
        if ~isempty(conn_x.Analysis_variables)&&isstruct(conn_x.Analysis_variables)&&~isfield(conn_x.Analysis_variables,'fbands'), conn_x.Analysis_variables.fbands=repmat({1},size(conn_x.Analysis_variables.names)); end
    end
    for ianalysis=1:numel(conn_x.Analyses)
        if ~isfield(conn_x.Analyses(ianalysis),'name'),conn_x.Analyses(ianalysis).name=['SBC_',num2str(ianalysis,'%02d')]; end
        if ~isfield(conn_x.Analyses(ianalysis),'sourcenames'),conn_x.Analyses(ianalysis).sourcenames={}; end
        if ~isfield(conn_x.Analyses(ianalysis),'modulation')||isempty(conn_x.Analyses(ianalysis).modulation),conn_x.Analyses(ianalysis).modulation=0; end
        if ~isfield(conn_x.Analyses(ianalysis),'measure'),conn_x.Analyses(ianalysis).measure=1; end
        if ~isfield(conn_x.Analyses(ianalysis),'weight'),conn_x.Analyses(ianalysis).weight=2; end
        if ~isfield(conn_x.Analyses(ianalysis),'type'),conn_x.Analyses(ianalysis).type=3; end
        if ~isfield(conn_x.Analyses(ianalysis),'conditions'),conn_x.Analyses(ianalysis).conditions={''}; end
        %if ~isempty(conn_x.Analyses(ianalysis).variables)&&isstruct(conn_x.Analyses(ianalysis).variables)&&~isfield(conn_x.Analyses(ianalysis).variables,'fbands'), conn_x.Analyses(ianalysis).variables.fbands=repmat({1},size(conn_x.Analyses(ianalysis).variables.names)); end
        if ~isempty(conn_x.Analyses(ianalysis).regressors)&&isstruct(conn_x.Analyses(ianalysis).regressors)&&~isfield(conn_x.Analyses(ianalysis).regressors,'fbands'), conn_x.Analyses(ianalysis).regressors.fbands=repmat({1},size(conn_x.Analyses(ianalysis).regressors.names)); end
    end    
    if ~isfield(conn_x,'vvAnalyses'), conn_x.vvAnalyses=struct('name','','measurenames',{{}},'variables', conn_v2v('measures'),'regressors', conn_v2v('measures'),'measures',{{}},'mask',[],'options',''); end
    if ~isfield(conn_x,'vvAnalysis'), conn_x.vvAnalysis=1; end
    for ianalysis=1:numel(conn_x.vvAnalyses)
        if ~isfield(conn_x.vvAnalyses(ianalysis),'name')||isempty(conn_x.vvAnalyses(ianalysis).name), conn_x.vvAnalyses(ianalysis).name=''; end
        if ~isfield(conn_x.vvAnalyses(ianalysis),'mask'), conn_x.vvAnalyses(ianalysis).mask=[]; end
        if ~isfield(conn_x.vvAnalyses(ianalysis),'options')||isempty(conn_x.vvAnalyses(ianalysis).options), conn_x.vvAnalyses(ianalysis).options=''; end
        if ~isfield(conn_x.vvAnalyses(ianalysis),'measures'), conn_x.vvAnalyses(ianalysis).measures={}; end
        if isfield(conn_x.vvAnalyses(ianalysis),'variables')&&~isfield(conn_x.vvAnalyses(ianalysis).variables,'norm'), conn_x.vvAnalyses(ianalysis).variables.norm=repmat({1},1,numel(conn_x.vvAnalyses(ianalysis).variables.names)); end
        if isfield(conn_x.vvAnalyses(ianalysis),'regressors')&&~isfield(conn_x.vvAnalyses(ianalysis).regressors,'norm'), conn_x.vvAnalyses(ianalysis).regressors.norm=repmat({1},1,numel(conn_x.vvAnalyses(ianalysis).regressors.names)); end
        if isfield(conn_x.vvAnalyses(ianalysis),'variables')&&~isfield(conn_x.vvAnalyses(ianalysis).variables,'alt_names'), conn_x.vvAnalyses(ianalysis).variables.alt_names=repmat({{}},1,numel(conn_x.vvAnalyses(ianalysis).variables.names)); end
        if isfield(conn_x.vvAnalyses(ianalysis),'regressors')&&~isfield(conn_x.vvAnalyses(ianalysis).regressors,'alt_names'), conn_x.vvAnalyses(ianalysis).regressors.alt_names=repmat({{}},1,numel(conn_x.vvAnalyses(ianalysis).regressors.names)); end
        %if isfield(conn_x.vvAnalyses(ianalysis),'variables')&&~isfield(conn_x.vvAnalyses(ianalysis).variables,'mask'), conn_x.vvAnalyses(ianalysis).variables.mask=repmat({{}},1,numel(conn_x.vvAnalyses(ianalysis).variables.names)); end
        %if isfield(conn_x.vvAnalyses(ianalysis),'regressors')&&~isfield(conn_x.vvAnalyses(ianalysis).regressors,'mask'), conn_x.vvAnalyses(ianalysis).regressors.mask=repmat({{}},1,numel(conn_x.vvAnalyses(ianalysis).regressors.names)); end
    end
    if ~isfield(conn_x,'dynAnalyses'), conn_x.dynAnalyses=struct('name','','regressors', struct('names',{{}}),'variables', struct('names',{{}}),'Ncomponents',[20],'condition',[1],'window',10,'output',[1 1 0]); end
    if ~isfield(conn_x,'dynAnalysis'), conn_x.dynAnalysis=1; end
    if isfield(conn_x.dynAnalyses,'analyses'), conn_x.dynAnalyses=rmfield(conn_x.dynAnalyses,'analyses'); end
    for ianalysis=1:numel(conn_x.dynAnalyses)
        if ~isfield(conn_x.dynAnalyses(ianalysis),'name'), conn_x.dynAnalyses(ianalysis).name=''; end
        if ~isfield(conn_x.dynAnalyses(ianalysis),'condition')||isempty(conn_x.dynAnalyses(ianalysis).condition), conn_x.dynAnalyses(ianalysis).condition=[]; end
        if ~isfield(conn_x.dynAnalyses(ianalysis),'window')||isempty(conn_x.dynAnalyses(ianalysis).window),
            if isfield(conn_x.dynAnalyses(ianalysis),'filter')&&~isempty(conn_x.dynAnalyses(ianalysis).filter), conn_x.dynAnalyses(ianalysis).window=round(1/conn_x.dynAnalyses(ianalysis).filter);
            else conn_x.dynAnalyses(ianalysis).window=30;
            end
        end
    end
    if ~isfield(conn_x,'Results')||~isfield(conn_x.Results,'saved')||isempty(conn_x.Results.saved), conn_x.Results.saved=struct('names',{{}},'labels',{{}},'descrip',{{}},'nsubjecteffects',{{}},'csubjecteffects',{{}},'nconditions',{{}},'cconditions',{{}}); end
    if ~isfield(conn_x.Results.saved,'names'), conn_x.Results.saved.names={}; end
    if ~isfield(conn_x,'isready')||numel(conn_x.isready)<2||conn_x.isready(2)
        for ianalysis=1:length(conn_x.Analyses),
            if ~conn_existfile(fullfile(conn_x.folders.firstlevel,conn_x.Analyses(ianalysis).name),true), conn_fileutils('mkdir',conn_x.folders.firstlevel,conn_x.Analyses(ianalysis).name); end;
            filesourcenames=fullfile(conn_x.folders.firstlevel,conn_x.Analyses(ianalysis).name,'_list_sources.mat');
            filesourcenames=conn_projectmanager('projectfile',filesourcenames,conn_x.pobj,'.mat');
            if conn_existfile(filesourcenames),
                sourcenames={}; conn_loadmatfile(filesourcenames,'sourcenames');
                if numel(sourcenames)>=numel(conn_x.Analyses(ianalysis).sourcenames), %~isempty(sourcenames)
                    if ~isequal(conn_x.Analyses(ianalysis).sourcenames,sourcenames), fprintf('warning: mismatch list_sources info. Using automatic recovery from %s\n',filesourcenames); end
                    conn_x.Analyses(ianalysis).sourcenames=sourcenames;
                end
            end
        end
        for ianalysis=1:length(conn_x.vvAnalyses),
            if ~conn_existfile(fullfile(conn_x.folders.firstlevel_vv,conn_x.vvAnalyses(ianalysis).name),true), conn_fileutils('mkdir',conn_x.folders.firstlevel_vv,conn_x.vvAnalyses(ianalysis).name); end;
            filemeasurenames=fullfile(conn_x.folders.firstlevel_vv,conn_x.vvAnalyses(ianalysis).name,'_list_measures.mat');
            filemeasurenames=conn_projectmanager('projectfile',filemeasurenames,conn_x.pobj,'.mat');
            if conn_existfile(filemeasurenames),
                measurenames={}; conn_loadmatfile(filemeasurenames,'measurenames');
                if numel(measurenames)>=numel(conn_x.vvAnalyses(ianalysis).measurenames), %~isempty(measurenames)
                    if ~isequal(conn_x.vvAnalyses(ianalysis).measurenames,measurenames), fprintf('warning: mismatch list_measures info. Using automatic recovery from %s\n',filemeasurenames); end
                    conn_x.vvAnalyses(ianalysis).measurenames=measurenames;
                end
            end
        end
        for ianalysis=1:length(conn_x.dynAnalyses),
            if ~conn_existfile(fullfile(conn_x.folders.firstlevel_dyn,conn_x.dynAnalyses(ianalysis).name),true), conn_fileutils('mkdir',conn_x.folders.firstlevel_dyn,conn_x.dynAnalyses(ianalysis).name); end;
        end
        fileconditionnames=fullfile(conn_x.folders.preprocessing,'_list_conditions.mat');
        fileconditionnames=conn_projectmanager('projectfile',fileconditionnames,conn_x.pobj,'.mat');
        if conn_existfile(fileconditionnames),
            allnames={}; conn_loadmatfile(fileconditionnames,'allnames');
            if numel(allnames)>=numel(conn_x.Setup.conditions.allnames), %~isempty(allnames)
                if ~isequal(conn_x.Setup.conditions.allnames,allnames), fprintf('warning: mismatch list_conditions info. Using automatic recovery from %s\n',fileconditionnames); end
                conn_x.Setup.conditions.allnames=allnames;
            end
        end
        fileresultsnames=fullfile(conn_x.folders.secondlevel,'_list_results.mat');
        if conn_existfile(fileresultsnames),
            results={}; conn_loadmatfile(fileresultsnames,'results');
            if ~isequal(conn_x.Results.saved,results), fprintf('warning: mismatch list_results info. Using automatic recovery from %s\n',fileresultsnames); end
            if ~isempty(results)
                conn_x.Results.saved=results;
            end
        end
        if ~isfield(conn_x.Results.saved,'labels'), conn_x.Results.saved.labels=repmat({''},size(conn_x.Results.saved.names)); end
        if ~isfield(conn_x.Results.saved,'nsubjecteffects'), conn_x.Results.saved.nsubjecteffects=cell(size(conn_x.Results.saved.names)); end
        if ~isfield(conn_x.Results.saved,'csubjecteffects'), conn_x.Results.saved.csubjecteffects=cell(size(conn_x.Results.saved.names)); end
        if ~isfield(conn_x.Results.saved,'nconditions'), conn_x.Results.saved.nconditions=cell(size(conn_x.Results.saved.names)); end
        if ~isfield(conn_x.Results.saved,'cconditions'), conn_x.Results.saved.cconditions=cell(size(conn_x.Results.saved.names)); end
        if ~isfield(conn_x.Results.saved,'descrip'), conn_x.Results.saved.descrip=repmat({''},size(conn_x.Results.saved.names)); end
    end
end
