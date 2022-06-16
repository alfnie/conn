function varargout=conn_process(options,varargin)

global CONN_x
if nargin<1, options=[]; end
if isequal(options,'aminserver') % if running from server
    if ~isfield(CONN_x,'gui'), conn_x_gui=0; else conn_x_gui=CONN_x.gui; end
    CONN_x.gui=varargin{1};
    skiploadsave=varargin{2};
    try
        if isnumeric(CONN_x.gui), CONN_x.gui=0;
        elseif isstruct(CONN_x.gui)
            if isfield(CONN_x.gui,'display'), CONN_x.gui.display=0; end
            if ~isfield(CONN_x.gui,'overwrite'), CONN_x.gui.overwrite='no'; end
        end
    end
    [varargout{1:nargout}]=conn_process(varargin{3:end});
    CONN_x.gui=conn_x_gui;
    if ~skiploadsave, 
        CONN_x.gui=1;
        conn save; 
    end
    return
elseif conn_projectmanager('inserver')&&isnumeric(options)&&nnz(~ismember(options,[1.5 5 9 9.1 9.2 9.3 17 17.5 19 34])) % if running from client connected to server; note: list of processes which may be run from client
    hmsg=[];
    if isfield(CONN_x,'gui')&&(isnumeric(CONN_x.gui)&&CONN_x.gui || isfield(CONN_x.gui,'display')&&CONN_x.gui.display),
        info=conn_remotely('info');
        if isfield(info,'host')&&~isempty(info.host), tnameserver=info.host;
        elseif isfield(info,'remote_ip')&&~isempty(info.remote_ip), tnameserver=info.remote_ip;
        else tnameserver='CONN server';
        end
        [hmsg,hstat]=conn_msgbox({sprintf('Process running remotely (%s)',tnameserver),'Please wait...',' ',' '},[],[],true);
    end
    if ~isfield(CONN_x,'filename')||isempty(CONN_x.filename), skiploadsave=true; else skiploadsave=false; end
    if ~isfield(CONN_x,'gui'), conn_x_gui=0; else conn_x_gui=CONN_x.gui; end
    if ~skiploadsave, conn save; end % note: save+push+rload
    if ~isempty(hmsg), [varargout{1:nargout}]=conn_server('run_withwaitbar',hstat,'conn_process','aminserver',conn_x_gui,skiploadsave,options,varargin{:});
    else [varargout{1:nargout}]=conn_server('run','conn_process','aminserver',conn_x_gui,skiploadsave,options,varargin{:});
    end
    if ~skiploadsave, conn load; end % note: (rload+rsave)+pull+load
    if ~isempty(hmsg)&&ishandle(hmsg), delete(hmsg); end
    return
end

if ischar(options),
    [optionsnow,options]=strtok(options,';');
    while ~isempty(optionsnow),
        %conn_disp('fprintf','conn_process %s started at %s\n',optionsnow,datestr(now));
        switch(lower(optionsnow)),
            case 'checkerrors',     conn_disp(['CONN: CHECKING PROJECT DEFINITIONS']); conn_process(0);
            case 'all',             conn_disp(['CONN: RUNNING ALL ANALYSIS STEPS']); conn_process([0:4,4.5,5,1.5,2,[6:11 13 14 15]]);
            case 'all_voxel',       conn_disp(['CONN: RUNNING ALL ANALYSIS STEPS (voxel-based analyses only)']); conn_process([0:3,5,1.5,2,6,9:10]);
            case 'all_roi',         conn_disp(['CONN: RUNNING ALL ANALYSIS STEPS (roi-based analyses only)']); conn_process([0:2,4,4.5,5,1.5,2,7,9,11,15]);
            case 'all_vv',          conn_disp(['CONN: RUNNING ALL ANALYSIS STEPS (voxel-to-voxel analyses only)']); conn_process([0:3,5,1.5,2,6,8:9,13]);
            case 'setup',           conn_disp(['CONN: RUNNING SETUP STEP']); conn_process([0:4,4.5,5]);
            case 'setup_roi',       conn_disp(['CONN: RUNNING SETUP STEP (roi-based analyses only)']); conn_process([0:2,4,4.5,5]);
            case 'setup_voxel',     conn_disp(['CONN: RUNNING SETUP STEP (roi-based analyses only)']); conn_process([0:3,5]);
            case 'setup_rois',      conn_disp(['CONN: RUNNING SETUP STEP (roi-based analyses only)']); conn_process([4,4.5,5]); % propagate rois
            case 'setup_covariates',conn_disp(['CONN: RUNNING SETUP STEP (covariate-setup only)']); conn_process([2 5]);        % propagate covariates
            case 'setup_conditions',conn_disp(['CONN: RUNNING SETUP STEP (condition-setup only)']); conn_process([1.5 2 5]);    % propagate conditions
            case 'setup_readyaggregate',conn_disp(['CONN: RUNNING SETUP STEP (condition-setup only)']); conn_process([0 2]);
            case 'setup_conditionsdecomposition',conn_disp(['CONN: RUNNING SETUP STEP (condition decomposition -setup only)']); conn_process([1.5]);
            case 'setup_skipchecks', conn_disp(['CONN: RUNNING SETUP STEP (skipping project integrity checks)']); conn_process([1:4,4.5,5]);
            case 'setup_skipmasks', conn_disp(['CONN: RUNNING SETUP STEP (skipping Grey/White/Mask processing)']); conn_process([2:4,4.5,5]);
            case 'setup_update',    conn_disp(['CONN: RUNNING SETUP STEP (update)']); conn_process([4.5 5]);
            case 'setup_preprocessing', conn_disp(['CONN: RUNNING SETUP.PREPROCESSING STEP']); conn_process(35,varargin{:});
            case 'setup_updatedenoising',    conn_disp(['CONN: RUNNING SETUP STEP (update)']); conn_process([5]);
            case {'preprocessing','denoising'},   conn_disp(['CONN: RUNNING DENOISING STEP']); conn_process([1.5,2,5:9],varargin{:});
            case {'preprocessing_gui','denoising_gui'}, conn_disp(['CONN: RUNNING DENOISING STEP']); conn_process([1.5,2,6:9],varargin{:});
            case {'preprocessing_roi','denoising_roi'}, conn_disp(['CONN: RUNNING DENOISING STEP (roi-based analyses only)']); conn_process([1.5,2,5,7,9],varargin{:});
            case 'denoising_finish', conn_process(9,varargin{:});
            case 'analyses',        conn_disp(['CONN: RUNNING ANALYSIS STEP']); conn_process([9,10,11,13,15],varargin{:});
            case 'analyses_nocombine', conn_disp(['CONN: RUNNING ANALYSIS STEP']); conn_process([9,10,11,13],varargin{:});
            case 'analyses_seed',   conn_disp(['CONN: RUNNING ANALYSIS STEP (ROI-to-ROI or seed-to-voxel analyses)']); conn_process([9.1,10],varargin{:});
            case 'analyses_seedsetup',   conn_disp(['CONN: RUNNING ANALYSIS STEP (ROI-to-ROI or seed-to-voxel analyses setup)']); conn_process([9.1],varargin{:});
            case 'analyses_vv',     conn_disp(['CONN: RUNNING ANALYSIS STEP (voxel-to-voxel analyses)']); conn_process([9.2,13],varargin{:});
            case 'analyses_vvmvpa',     conn_disp(['CONN: RUNNING ANALYSIS STEP (voxel-to-voxel MVPA analyses)']); conn_process([9.2,13.1],varargin{:});
            case 'analyses_vvsetup',     conn_disp(['CONN: RUNNING ANALYSIS STEP (voxel-to-voxel MVPA analyses setup)']); conn_process([9.2],varargin{:});
            case 'analyses_roi',    conn_disp(['CONN: RUNNING ANALYSIS STEP (roi-based analyses only)']); conn_process([9.1,11,15],varargin{:});
            case 'analyses_seedandroi', conn_disp(['CONN: RUNNING ANALYSIS STEP (seed-to-voxel and roi-to-roi analyses)']); conn_process([9.1,10,11,15],varargin{:});
            case 'analyses_dyn',    conn_disp(['CONN: RUNNING DYNAMIC CONNECTIVITY STEP']); conn_process([9.3,14],varargin{:});
            case 'analyses_dyn_step1', conn_disp(['CONN: RUNNING DYNAMIC CONNECTIVITY STEP']); conn_process([9.3,14.1],varargin{:});
            case 'analyses_dyn_step2', conn_disp(['CONN: RUNNING DYNAMIC CONNECTIVITY STEP']); conn_process([9.3,14.2],varargin{:});
            case 'analyses_dynsetup',    conn_disp(['CONN: RUNNING DYNAMIC CONNECTIVITY STEP (setup)']); conn_process([9.3],varargin{:});
            case 'analyses_gui',    conn_disp(['CONN: RUNNING ANALYSIS STEP']); conn_process([10,11,13,15],varargin{:});
            case 'analyses_gui_seedandroi',conn_disp(['CONN: RUNNING ANALYSIS STEP (ROI-to-ROI or seed-to-voxel analyses)']); conn_process([10,11,15],varargin{:});
            case 'analyses_gui_vv', conn_disp(['CONN: RUNNING ANALYSIS STEP (voxel-to-voxel analyses)']); conn_process([13],varargin{:});
            case 'analyses_gui_vvmvpa', conn_disp(['CONN: RUNNING ANALYSIS STEP (voxel-to-voxel MVPA analyses)']); conn_process([13.1],varargin{:});
            case 'analyses_gui_dyn',conn_disp(['CONN: RUNNING DYNAMIC CONNECTIVITY STEP']); conn_process(14,varargin{:});
            case 'analyses_gui_dyn_step1',conn_disp(['CONN: RUNNING DYNAMIC CONNECTIVITY STEP']); conn_process(14.1,varargin{:});
            case 'analyses_gui_dyn_step2',conn_disp(['CONN: RUNNING DYNAMIC CONNECTIVITY STEP']); conn_process(14.2,varargin{:});
            case 'results_nonparametric',conn_process(20,varargin{:});
            case 'results_nonparametricroi',conn_process(21,varargin{:});
            case 'results_nonparametric_collapse',conn_process(22,varargin{:});
            case 'extract_connectome',conn_process(31,varargin{:});
            case 'prepare_results', conn_process(15,varargin{:});
            case 'prepare_results_roi',conn_process(15,varargin{:});
            case 'results',         conn_process(16:18,varargin{:});
            case 'results_voxel',   [varargout{1:nargout}]=conn_process(16,varargin{:});
            case 'results_roi',     [varargout{1:nargout}]=conn_process(17,varargin{:});
            case 'results_roi_seed',[varargout{1:nargout}]=conn_process(17.5,varargin{:});
            case 'results_connectome',[varargout{1:nargout}]=conn_process(18,varargin{:});
            case 'postmerge',       conn_process([4.5,5,9,15],varargin{:});
            case 'maskserode',      conn_process(33,varargin{:});
            case 'qaplots',         conn_process(32,varargin{:});
            case 'nyudataset',      conn_process(34,varargin{:});
            case 'update',          conn gui_setup_saveas; conn_process all; conn save;
            case 'conn',            conn(varargin{:});
            case 'batch',           conn_batch(varargin{:});
            case 'fcn',             if ischar(varargin{1}), fh=eval(sprintf('@%s',varargin{1})); else fh=varargin{1}; end
                                    feval(fh,varargin{2:end});
            case 'spmbatch',
                spm_jobman('initcfg');
                try, spm_get_defaults('mat.format','-v7.3'); end
                warning('off','MATLAB:RandStream:ActivatingLegacyGenerators');
                job_id=spm_jobman('run',varargin{:});
                warning('on','MATLAB:RandStream:ActivatingLegacyGenerators');
            case 'multiplesteps'
                arglist=varargin{1};
                for n=1:numel(arglist),
                    conn_process(arglist{n}{:});
                end
            otherwise,
                if all(ismember(options,'0123456789.,()[]+- ')), conn_process(str2num(options));
                else disp(sprintf('conn_process: unrecognized option %s',options));
                end
        end
        %conn_disp('fprintf','conn_process %s ended at %s\n',optionsnow,datestr(now));
        [optionsnow,options]=strtok(options,';');
    end        
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checks valid/complete setup fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(options==0),
    if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'subjects'), validsubjects=CONN_x.gui.subjects; else validsubjects=1:CONN_x.Setup.nsubjects; end
    if isfield(CONN_x,'pobj')&&isstruct(CONN_x.pobj)&&isfield(CONN_x.pobj,'subjects'), validsubjects=CONN_x.pobj.subjects; end
    h=conn_waitbar(0,['Step ',num2str(sum(options<=0)),'/',num2str(length(options)),': Checking data completeness']);
    CHECKEMPTYROIS=false;
    CHECKFILEPATHS=isequal(validsubjects,1:CONN_x.Setup.nsubjects);
    ERR={};
	if isempty(CONN_x.filename)||isempty(dir(CONN_x.filename)),ERR{end+1}=['WARNING: Project not saved or empty project information']; end
	nl1covariates=length(CONN_x.Setup.l1covariates.names)-1;
	nconditions=length(CONN_x.Setup.conditions.names)-1;
    if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'conditions')&&~isempty(CONN_x.gui.conditions), validconditions=CONN_x.gui.conditions; else validconditions=1:length(CONN_x.Setup.conditions.names)-1; end
    nrois=length(CONN_x.Setup.rois.names)-1;
    if nconditions<1, ERR{end+1}=['ERROR: No conditions defined. You must define at least one condition (use onset=0; duration=inf; for encompassing the entire scanning session)']; end
	for nsub=validsubjects,
		nsess=CONN_x.Setup.nsessions(min(length(CONN_x.Setup.nsessions),nsub));
        if nsess<1,ERR{end+1}=['ERROR: Subject ',num2str(nsub),' has no functional data assigned (number of sessions = 0); increase the number of sessions in Setup->Basic or remove this subject']; end
        if numel(CONN_x.Setup.conditions.values)<nsub, ERR{end+1}=['ERROR: Subject ',num2str(nsub),' has no experimental conditions defined']; end
        if numel(CONN_x.Setup.conditions.values)>=nsub&&numel(CONN_x.Setup.conditions.values{nsub})<nconditions, ERR{end+1}=['ERROR: Missing experimental conditions on Subject ',num2str(nsub)]; end
        for nroi=4:nrois,%1:nrois,
            if ~CONN_x.Setup.rois.sessionspecific(nroi), nsesstemp=1; else nsesstemp=nsess; end
            if ~CONN_x.Setup.rois.subjectspecific(nroi), nsubstemp=1; else nsubstemp=nsub; end
            for nses=1:nsesstemp,
                if numel(CONN_x.Setup.rois.files)<nsubstemp||numel(CONN_x.Setup.rois.files{nsubstemp})<nroi||numel(CONN_x.Setup.rois.files{nsubstemp}{nroi})<nses||numel(CONN_x.Setup.rois.files{nsubstemp}{nroi}{nses})<3||isempty(CONN_x.Setup.rois.files{nsubstemp}{nroi}{nses}{1}),
                    ERR{end+1}=['ERROR: Subject ',num2str(nsubstemp),' Session ',num2str(nses),' ROI ',num2str(nroi), ' file has not been defined'];
                elseif CHECKEMPTYROIS&&isstruct(CONN_x.Setup.rois.files{nsubstemp}{nroi}{nses}{3})
                    try
                        a=spm_vol(deblank(CONN_x.Setup.rois.files{nsubstemp}{nroi}{nses}{1}));
                        if ~nnz(spm_read_vols(a)>0)
                            ERR{end+1}=['ERROR: Subject ',num2str(nsubstemp),' Session ',num2str(nses),' ROI ',num2str(nroi), ' no voxels with >0 values in file ',deblank(CONN_x.Setup.rois.files{nsubstemp}{nroi}{nses}{1})];
                        end
                    catch
                        conn_disp(['warning: Subject ',num2str(nsubstemp),' Session ',num2str(nses),' ROI ',num2str(nroi), ' reading file ',deblank(CONN_x.Setup.rois.files{nsubstemp}{nroi}{nses}{1})]);
                    end
                end
            end
        end
        okconditions=false(1,nconditions);
        oksessions=false(1,nsess);
		for nses=1:nsess,
            if length(CONN_x.Setup.nscans)<nsub, ERR{end+1}=['ERROR: Subject ',num2str(nsub),' has no functional data defined']; 
            elseif length(CONN_x.Setup.nscans{nsub})<nses, ERR{end+1}=['ERROR: Subject ',num2str(nsub),' has no functional data defined for session ',num2str(nses)]; 
            else,
                nscans=CONN_x.Setup.nscans{nsub}{nses};
                if isempty(nscans)||nscans<1,ERR{end+1}=['ERROR: Subject ',num2str(nsub),' Session ',num2str(nses), ' functional data has not been defined']; 
                elseif nscans<2, ERR{end+1}=['ERROR: Subject ',num2str(nsub),' Session ',num2str(nses), ' functional data has only ',num2str(nscans),' time-points (scans) defined']; 
                elseif nscans<16, conn_disp(['warning: Subject ',num2str(nsub),' Session ',num2str(nses), ' functional data has only ',num2str(nscans),' time-points (scans) defined']); 
                end
                for nl1covariate=1:nl1covariates,
                    if numel(CONN_x.Setup.l1covariates.files)<nsub||numel(CONN_x.Setup.l1covariates.files{nsub})<nl1covariate||numel(CONN_x.Setup.l1covariates.files{nsub}{nl1covariate})<nses||numel(CONN_x.Setup.l1covariates.files{nsub}{nl1covariate}{nses})<1
                        ERR{end+1}=['ERROR: Missing first-level covariate ',num2str(nl1covariate),' for Subject ',num2str(nsub),' Session ',num2str(nses)];
                    else
                        filename=CONN_x.Setup.l1covariates.files{nsub}{nl1covariate}{nses}{1};
                        names=CONN_x.Setup.l1covariates.names{nl1covariate};
                        if isempty(filename)||~ischar(filename), ERR{end+1}=['ERROR: Subject ',num2str(nsub),' Session ',num2str(nses), ' first-level covariate ',names,' not defined'];
                        else
                            switch(filename),
                                case '[raw values]',
                                    data=CONN_x.Setup.l1covariates.files{nsub}{nl1covariate}{nses}{3};
                                otherwise,
                                    data=conn_loadtextfile(filename,false);
                                    %if isstruct(data), tempnames=fieldnames(data); data=data.(tempnames{1}); end
                            end
                            if size(data,1)~=nscans, 
                                if isempty(data)&&isequal(filename,'[raw values]'),
                                else ERR{end+1}=['ERROR: Subject ',num2str(nsub),' Session ',num2str(nses), ' first-level covariate ',names,' mismatched dimensions (',num2str(size(data,1)),' rows, while functional data has ',num2str(nscans),' scans; the number of rows of a first-level covariate should equal the number of scans for this subject/session)']; 
                                end
                            end
                        end
                    end
                end
                for ncondition=validconditions,
                    if numel(CONN_x.Setup.conditions.values)<nsub||numel(CONN_x.Setup.conditions.values{nsub})<ncondition||numel(CONN_x.Setup.conditions.values{nsub}{ncondition})<nses||numel(CONN_x.Setup.conditions.values{nsub}{ncondition}{nses})<2, 
                        if strcmp(CONN_x.Setup.conditions.names{ncondition},'rest')
                            conn_disp(['note: Subject ',num2str(nsub),' condition ',CONN_x.Setup.conditions.names{ncondition},' has not been defined for session ',num2str(nses),'. Assuming this condition is present in this session']);
                            CONN_x.Setup.conditions.values{nsub}{ncondition}{nses}={0,inf};
                        else
                            conn_disp(['note: Subject ',num2str(nsub),' condition ',CONN_x.Setup.conditions.names{ncondition},' has not been defined for session ',num2str(nses),'. Assuming this condition is not present in this session']);
                            CONN_x.Setup.conditions.values{nsub}{ncondition}{nses}={[],[]};
                        end
                    elseif ~isempty(CONN_x.Setup.conditions.model{ncondition})
                        okconditions(ncondition)=true;
                    else
                        oktemp=~isempty(CONN_x.Setup.conditions.values{nsub}{ncondition}{nses}{1})&&~isempty(CONN_x.Setup.conditions.values{nsub}{ncondition}{nses}{2});
                        okconditions(ncondition)=okconditions(ncondition)|oktemp;
                        oksessions(nses)=oksessions(nses)|oktemp;
                    end
                    if any(CONN_x.Setup.conditions.values{nsub}{ncondition}{nses}{1}>conn_get_rt(nsub,nses)*CONN_x.Setup.nscans{nsub}{nses})
                        ERR{end+1}=['ERROR: Subject ',num2str(nsub),' Session ',num2str(nses), ' condition ',CONN_x.Setup.conditions.names{ncondition},' contains onset times beyond the end of the scanning session']; 
                    end
                end
            end
        end
        if any(~okconditions(validconditions))
            for ncondition=validconditions(~okconditions(validconditions));
                if CONN_x.Setup.conditions.missingdata, conn_disp(['note: Subject ',num2str(nsub),' does not have any scan associated with condition ',CONN_x.Setup.conditions.names{ncondition}]); 
                else ERR{end+1}=['ERROR: Subject ',num2str(nsub),' does not have any scan associated with condition ',CONN_x.Setup.conditions.names{ncondition},'. If this is expected/correct please select ''Allow missing data'' in Setup.Conditions to avoid this error message'];
                end
            end
        end
        if any(~oksessions)&&isequal(validconditions(:)',1:length(CONN_x.Setup.conditions.names)-1)
            for nses=reshape(find(~oksessions),1,[]);
                if CONN_x.Setup.conditions.missingdata, conn_disp(['note: Subject ',num2str(nsub),' does not have any condition associated with data from session ',num2str(nses)]); 
                else ERR{end+1}=['ERROR: Subject ',num2str(nsub),' does not have any condition associated with data from session ',num2str(nses),'. If this is expected/correct please select ''Allow missing data'' in Setup.Conditions to avoid this error message'];
                end
            end
        end
        CONN_x.Setup.functional{nsub}=CONN_x.Setup.functional{nsub}(1:min(numel(CONN_x.Setup.functional{nsub}),nsess)); 
        CONN_x.Setup.structural{nsub}=CONN_x.Setup.structural{nsub}(1:min(numel(CONN_x.Setup.structural{nsub}),nsess)); 
        try, for nl1covariate=1:nl1covariates,CONN_x.Setup.l1covariates.files{nsub}{nl1covariate}=CONN_x.Setup.l1covariates.files{nsub}{nl1covariate}(1:min(numel(CONN_x.Setup.l1covariates.files{nsub}{nl1covariate}),nsess));end; end
        try, for nroi=1:nrois,CONN_x.Setup.rois.files{nsub}{nroi}=CONN_x.Setup.rois.files{nsub}{nroi}(1:min(numel(CONN_x.Setup.rois.files{nsub}{nroi}),nsess));end; end
    end
    if CONN_x.Setup.spatialresolution==4
        for nsub=validsubjects,
            nsess=CONN_x.Setup.nsessions(min(length(CONN_x.Setup.nsessions),nsub));
            if ~CONN_x.Setup.structural_sessionspecific, nsesstemp=1; else nsesstemp=nsess; end
            for nses=1:nsesstemp,
                if ~conn_checkFSfiles(CONN_x.Setup.structural{nsub}{nses}{3})
                    ERR{end+1}=['ERROR: Subject ',num2str(nsub),' Session ',num2str(nses),'. No Freesurfer-generated surface files found. See Setup->Structural for more details'];
                end
            end
        end
    end
    for ncondition=validconditions,
        if numel(CONN_x.Setup.conditions.model)<ncondition, CONN_x.Setup.conditions.model{ncondition}=[]; end
        if numel(CONN_x.Setup.conditions.param)<ncondition, CONN_x.Setup.conditions.param(ncondition)=0; end
        if numel(CONN_x.Setup.conditions.filter)<ncondition, CONN_x.Setup.conditions.filter{ncondition}=[]; end
        if CONN_x.Setup.conditions.param(ncondition)>numel(CONN_x.Setup.l1covariates.names)-1
            ERR{end+1}=['ERROR: Condition ',CONN_x.Setup.conditions.names{ncondition},' temporal-modulation links to non-existing first-level covariate ',num2str(CONN_x.Setup.conditions.param(ncondition))];
        end
    end
    % other warnings or common issues
    conn_checkdistributionfiles;
    if isnan(spm_type('float32'))
        ERR{end+1}=['ERROR: spm_type overloaded by SPM2 version in folder ',fileparts(which('spm_type')),'. Update this file to a more recent version or deprecate this folder in the Matlab path to somewhere below SPM folder'];
    end
    if ~isempty(ERR),
        conn_disp('SOME ERRORS FOUND!: Please revise the Setup information');
        if isfield(CONN_x,'gui')&&(isnumeric(CONN_x.gui)&&CONN_x.gui || isfield(CONN_x.gui,'display')&&CONN_x.gui.display),
            error(['<nodetailsflag>',sprintf('%s\n',ERR{:})]);
        else
            conn_disp(strvcat(ERR));
            error(sprintf('%s\n',ERR{:}));
        end
        return;
    end
    if CHECKFILEPATHS,
        try, conn_updatefilepaths; catch, ERR{end+1}=['ERROR: Files not found. Check that all your functional/structural/ROIs/first-level covariate files point to existing files.']; end
    end
    % check coregistration of reference files
    CONN_x.Setup.normalized=1;
    if ~ismember(CONN_x.Setup.spatialresolution,[1 4])
        ok=1;
        for nsub=validsubjects,
            switch(CONN_x.Setup.spatialresolution)
                case 2, v=CONN_x.Setup.structural{nsub}{1}{3}.dim;
                case 3, v=CONN_x.Setup.functional{nsub}{1}{3}.dim;
                otherwise, error('Invalid value in batch.Setup.voxelresolution');
            end
            if nsub==validsubjects(1), vref=v;
            elseif any(v~=vref), ok=0; 
            end
        end
        if ~ok,
            switch(CONN_x.Setup.spatialresolution)
                case 2, vname='structural';
                case 3, vname='functional';
            end
            txt={sprintf('Warning: %s volumes do not have the same dimensions',vname),'Multi-subject voxel-level analyses will NOT be performed.','Choose ''same as mask'' option in ''analysis space'' field in Setup.Options if you wish','to resample the data to a common voxel-resolution instead'};
            if isfield(CONN_x,'gui')&&(isnumeric(CONN_x.gui)&&CONN_x.gui || isfield(CONN_x.gui,'display')&&CONN_x.gui.display), 
                  answ= conn_questdlg(txt,'','Continue', 'Stop','Stop');
                  if strcmp(answ,'Stop'), return; end
                  CONN_x.Setup.normalized=0;
            else
                conn_disp(char(txt));
                CONN_x.Setup.normalized=0;
            end
        end
    end
	conn_waitbar('close',h);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Segments structural data if appropriate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(options==1),
    if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'subjects'), validsubjects=CONN_x.gui.subjects; else validsubjects=1:CONN_x.Setup.nsubjects; end
    if isfield(CONN_x,'pobj')&&isstruct(CONN_x.pobj)&&isfield(CONN_x.pobj,'subjects'), validsubjects=CONN_x.pobj.subjects; end
    h=conn_waitbar(0,['Step ',num2str(sum(options<=1)),'/',num2str(length(options)),': Segmentation']);
	[path,name,ext]=fileparts(CONN_x.filename);
	for nsub=validsubjects,
		nsess=CONN_x.Setup.nsessions(min(length(CONN_x.Setup.nsessions),nsub));
        if ~CONN_x.Setup.structural_sessionspecific, nsesstemp=1; else nsesstemp=nsess; end
        dosegm=0;for n1=1:3, for nses=1:nsesstemp, if length(CONN_x.Setup.rois.files{nsub})<n1||numel(CONN_x.Setup.rois.files{nsub}{n1})<nses||isempty(CONN_x.Setup.rois.files{nsub}{n1}{nses})||~iscell(CONN_x.Setup.rois.files{nsub}{n1}{nses})||isempty(CONN_x.Setup.rois.files{nsub}{n1}{nses}{1}),dosegm=1;end;end;end
        if dosegm,
            dosegm=0;
            file=cell(1,nsesstemp);
            for nses=1:nsesstemp
                file{nses}=deblank(CONN_x.Setup.structural{nsub}{nses}{1});
                if ~conn_existfile(conn_prepend('c1',file{nses}))||~conn_existfile(conn_prepend('c2',file{nses}))||~conn_existfile(conn_prepend('c3',file{nses})), dosegm=1; end
            end
            if dosegm
                conn_setup_preproc('run_structural_segment','subjects',nsub);
            else
                for n1=1:3,
                    for nses=1:nsess
                        ffile=conn_prepend(['c',num2str(n1)],file{min(nsesstemp,nses)});
                        CONN_x.Setup.rois.files{nsub}{n1}{nses}=conn_file(ffile);
                        %[V,str,icon]=conn_getinfo(ffile);
                        %CONN_x.Setup.rois.files{nsub}{n1}={ffile,str,icon};
                    end
                end
            end
		end
		conn_waitbar(sum(validsubjects<=nsub)/numel(validsubjects),h,sprintf('Subject %d',nsub));
	end
	conn_waitbar('close',h);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates non-parametric temporal/frequency decomposition conditions if necessary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(options==1.5),
	nconditions=length(CONN_x.Setup.conditions.names)-1;
    if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'conditions')&&~isempty(CONN_x.gui.conditions), validconditions=CONN_x.gui.conditions; else validconditions=1:length(CONN_x.Setup.conditions.names)-1; end
    if length(CONN_x.Setup.conditions.model)<nconditions, CONN_x.Setup.conditions.model=[CONN_x.Setup.conditions.model, cell(1,nconditions-length(CONN_x.Setup.conditions.model))]; end
    if length(CONN_x.Setup.conditions.param)<nconditions, CONN_x.Setup.conditions.param=[CONN_x.Setup.conditions.param, zeros(1,nconditions-length(CONN_x.Setup.conditions.param))]; end
    if length(CONN_x.Setup.conditions.filter)<nconditions, CONN_x.Setup.conditions.filter=[CONN_x.Setup.conditions.filter, cell(1,nconditions-length(CONN_x.Setup.conditions.filter))]; end
    maxrt=nan;
	h=conn_waitbar(0,['Step ',num2str(sum(options<=1.5)),'/',num2str(length(options)),': Expanding conditions']);
    for ncondition=validconditions,
        % Frequency-band non-parametric modulation
        if ncondition>=length(CONN_x.Setup.conditions.names)
        elseif numel(CONN_x.Setup.conditions.filter{ncondition})==1, 
            nbands=CONN_x.Setup.conditions.filter{ncondition};
            % adds new conditions
            newcond=[];
            for nparam=1:nbands+2,
                fband=[];
                model=[];
                if nparam==nbands+1
                    condname=[CONN_x.Setup.conditions.names{ncondition},' x Frequency Average'];
                    model=[{'avg'},CONN_x.Setup.conditions.names(newcond(1:nbands))];
                elseif nparam==nbands+2
                    condname=[CONN_x.Setup.conditions.names{ncondition},' x Frequency Variability'];
                    model=[{'std'},CONN_x.Setup.conditions.names(newcond(1:nbands))];
                else
                    ffilter=CONN_x.Preproc.filter;
                    if any(isinf(ffilter))&&isnan(maxrt), maxrt=max(conn_get_rt); end
                    ffilter(isinf(ffilter))=1/maxrt/2;
                    fband=ffilter(1)+(ffilter(2)-ffilter(1))/nbands*[nparam-1,nparam];
                    condname=[CONN_x.Setup.conditions.names{ncondition},' x FrequencyBand',num2str(nparam)];
                end
                idx=strmatch(condname,CONN_x.Setup.conditions.names(1:nconditions),'exact');
                if ~isempty(idx),
                    icond=idx(1);
                else,
                    icond=numel(CONN_x.Setup.conditions.names);
                end
                if icond==numel(CONN_x.Setup.conditions.names), CONN_x.Setup.conditions.names{icond+1}=' '; end
                CONN_x.Setup.conditions.names{icond}=condname;
                for n1=1:length(CONN_x.Setup.conditions.values), CONN_x.Setup.conditions.values{n1}{icond}=CONN_x.Setup.conditions.values{n1}{ncondition}; end
                CONN_x.Setup.conditions.model{icond}=model;
                CONN_x.Setup.conditions.param(icond)=CONN_x.Setup.conditions.param(ncondition);
                CONN_x.Setup.conditions.filter{icond}=fband;
                newcond=[newcond icond];
            end
            if ~isempty(newcond)&&isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'conditions'),
                CONN_x.gui.conditions=union(CONN_x.gui.conditions,newcond);
            end
            % removes other matching conditions
            condname=[CONN_x.Setup.conditions.names{ncondition},' x Frequency'];
            nconditions0=numel(CONN_x.Setup.conditions.names);
            idx=setdiff(strmatch(condname,CONN_x.Setup.conditions.names(1:nconditions0-1)),newcond); % other cond
            if ~isempty(idx)
                idx=setdiff(1:nconditions0,idx); % keep these
                CONN_x.Setup.conditions.names=CONN_x.Setup.conditions.names(idx);
                idx=setdiff(idx,nconditions0);
                for n1=1:length(CONN_x.Setup.conditions.values), CONN_x.Setup.conditions.values{n1}=CONN_x.Setup.conditions.values{n1}(idx); end
                CONN_x.Setup.conditions.model=CONN_x.Setup.conditions.model(idx);
                CONN_x.Setup.conditions.param=CONN_x.Setup.conditions.param(idx);
                CONN_x.Setup.conditions.filter=CONN_x.Setup.conditions.filter(idx);
                if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'conditions'),
                    [ok,idx2]=ismember(CONN_x.gui.conditions,idx);
                    CONN_x.gui.conditions=idx2(ok);
                end
            end
        % Temporal sliding-window non-parametric modulation
        elseif numel(CONN_x.Setup.conditions.filter{ncondition})>2, 
            duration=CONN_x.Setup.conditions.filter{ncondition}(1);
            ronsets=CONN_x.Setup.conditions.filter{ncondition}(2:end);
            nsteps=numel(ronsets);
            % adds new conditions
            newcond=[];
            for nparam=1:nsteps+2,
                model=[];
                if nparam==nsteps+1
                    condname=[CONN_x.Setup.conditions.names{ncondition},' x Temporal Average'];
                    model=[{'avg'},CONN_x.Setup.conditions.names(newcond(1:nsteps))];
                elseif nparam==nsteps+2
                    condname=[CONN_x.Setup.conditions.names{ncondition},' x Temporal Variability'];
                    model=[{'std'},CONN_x.Setup.conditions.names(newcond(1:nsteps))];
                else
                    condname=[CONN_x.Setup.conditions.names{ncondition},' x Time',num2str(nparam)];
                end
                idx=strmatch(condname,CONN_x.Setup.conditions.names(1:nconditions),'exact');
                if ~isempty(idx),
                    icond=idx(1);
                else,
                    icond=numel(CONN_x.Setup.conditions.names);
                end
                if icond==numel(CONN_x.Setup.conditions.names), CONN_x.Setup.conditions.names{icond+1}=' '; end
                CONN_x.Setup.conditions.names{icond}=condname;
                if nparam>nsteps, 
                    for n1=1:length(CONN_x.Setup.conditions.values), CONN_x.Setup.conditions.values{n1}{icond}=CONN_x.Setup.conditions.values{n1}{ncondition}; end
                else
                    for n1=1:length(CONN_x.Setup.conditions.values),
                        for n2=1:length(CONN_x.Setup.conditions.values{n1}{ncondition}),
                            CONN_x.Setup.conditions.values{n1}{icond}{n2}{1}=CONN_x.Setup.conditions.values{n1}{ncondition}{n2}{1}+ronsets(nparam);
                            CONN_x.Setup.conditions.values{n1}{icond}{n2}{2}=duration;
                        end
                    end
                end
                CONN_x.Setup.conditions.model{icond}=model;
                CONN_x.Setup.conditions.param(icond)=CONN_x.Setup.conditions.param(ncondition);
                CONN_x.Setup.conditions.filter{icond}=[];
                newcond=[newcond icond];
            end
            if ~isempty(newcond)&&isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'conditions'),
                CONN_x.gui.conditions=union(CONN_x.gui.conditions,newcond);
            end
            % removes other matching conditions
            condname1=[CONN_x.Setup.conditions.names{ncondition},' x Time'];
            condname2=[CONN_x.Setup.conditions.names{ncondition},' x Temporal'];
            nconditions0=numel(CONN_x.Setup.conditions.names);
            idx=setdiff(union(strmatch(condname1,CONN_x.Setup.conditions.names(1:nconditions0-1)),strmatch(condname2,CONN_x.Setup.conditions.names(1:nconditions0-1))),newcond); % other cond
            if ~isempty(idx)
                idx=setdiff(1:nconditions0,idx); % keep these
                CONN_x.Setup.conditions.names=CONN_x.Setup.conditions.names(idx);
                idx=setdiff(idx,nconditions0);
                for n1=1:length(CONN_x.Setup.conditions.values), CONN_x.Setup.conditions.values{n1}=CONN_x.Setup.conditions.values{n1}(idx); end
                CONN_x.Setup.conditions.model=CONN_x.Setup.conditions.model(idx);
                CONN_x.Setup.conditions.param=CONN_x.Setup.conditions.param(idx);
                CONN_x.Setup.conditions.filter=CONN_x.Setup.conditions.filter(idx);
                if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'conditions'),
                    [ok,idx2]=ismember(CONN_x.gui.conditions,idx);
                    CONN_x.gui.conditions=idx2(ok);
                end
            end
        end
		conn_waitbar(mean(validconditions<=ncondition),h,sprintf('Condition %s',CONN_x.Setup.conditions.names{ncondition}));
	end
	conn_waitbar('close',h);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates COV_Subject###_Session###.mat files (first-level temporal covariate values)
% Creates COND_Subject###_Session###.mat files (temporal samples included in each condition)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(options==2),
	[path,name,ext]=fileparts(CONN_x.filename);
	%[ok,nill]=mkdir(path,name);
	%[ok,nill]=mkdir(fullfile(path,name),'data');
	%filepath=fullfile(path,name,'data');
    filepath=CONN_x.folders.data;
    try, conn_fileutils('mkdir',filepath); end
    if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'subjects'), validsubjects=CONN_x.gui.subjects; else validsubjects=1:CONN_x.Setup.nsubjects; end
    if isfield(CONN_x,'pobj')&&isstruct(CONN_x.pobj)&&isfield(CONN_x.pobj,'subjects'), validsubjects=CONN_x.pobj.subjects; end
	h=conn_waitbar(0,['Step ',num2str(sum(options<=2)),'/',num2str(length(options)),': Importing conditions/covariates']);
    REDO='Yes';%[];
	nl1covariates=length(CONN_x.Setup.l1covariates.names)-1;
	nconditions=length(CONN_x.Setup.conditions.names)-1;
	if length(CONN_x.Setup.nsessions)==1, N=CONN_x.Setup.nsessions*numel(validsubjects); else N=sum(CONN_x.Setup.nsessions(validsubjects)); end
	N=N*(nl1covariates+nconditions);
    maxdims=[];n=0;
    CROPCONDITIONSAMPLES=false; % new behavior to allow temporal modulation analyses (Preprocessing condition-specific files include zero-value weights)
%     uniqueconditions=[];
%     for n1=1:nconditions
%         ok=true;
%         for n2=1:n1-1
%             ok=ok&~all(arrayfun(@(n)isequal(CONN_x.Setup.conditions.values{n}{n1},CONN_x.Setup.conditions.values{n}{n2}),1:numel(CONN_x.Setup.conditions.values)));
%         end
%         if ok, uniqueconditions=[uniqueconditions, n1]; end
%     end
	for nsub=validsubjects,
		nsess=CONN_x.Setup.nsessions(min(length(CONN_x.Setup.nsessions),nsub));
		for nses=1:nsess,
            filename=fullfile(filepath,['COND_Subject',num2str(nsub,'%03d'),'_Session',num2str(nses,'%03d'),'.mat']);
            if isempty(REDO)&&conn_existfile(filename),
                if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'overwrite'), REDO=CONN_x.gui.overwrite; else, REDO=conn_questdlg('Overwrite existing subject results?','','Yes', 'No', 'Yes');end;
            end
            if strcmp(lower(REDO),'yes')||~conn_existfile(filename),
                clear data names;
                clear samples weights;
                nscans=CONN_x.Setup.nscans;
                for nl1covariate=1:nl1covariates,
                    filename=CONN_x.Setup.l1covariates.files{nsub}{nl1covariate}{nses}{1};
                    switch(filename),
                        case '[raw values]',
                            data{nl1covariate}=CONN_x.Setup.l1covariates.files{nsub}{nl1covariate}{nses}{3};
                        otherwise,
                            data{nl1covariate}=conn_loadtextfile(filename,false);
                            %if isstruct(data{nl1covariate}), tempnames=fieldnames(data{nl1covariate}); data{nl1covariate}=data{nl1covariate}.(tempnames{1}); end
                    end
                    names{nl1covariate}=CONN_x.Setup.l1covariates.names{nl1covariate};
                    n=n+1;
                    conn_waitbar(n/N,h,sprintf('Subject %d Session %d',nsub,nses));
                end
                crop=CROPCONDITIONSAMPLES;
                RT=[];
                for ncondition=1:nconditions,
                    if ~isempty(CONN_x.Setup.conditions.model{ncondition})
                        data{nl1covariates+ncondition}=[];
                        names{nl1covariates+ncondition}='';
                        samples{ncondition}=[];
                        weights{ncondition}=cell(1,5);
                    else
                        onset=CONN_x.Setup.conditions.values{nsub}{ncondition}{nses}{1};
                        durat=CONN_x.Setup.conditions.values{nsub}{ncondition}{nses}{2};
                        if isempty(RT), RT=conn_get_rt(nsub,nses); end
                        rt=RT/10;
                        offs=ceil(100/rt);
                        hrf=spm_hrf(rt);
                        x=zeros(offs+ceil(CONN_x.Setup.nscans{nsub}{nses}*RT/rt),1);
                        for n1=1:length(onset),
                            tdurat=max(rt,min(offs*rt+RT*CONN_x.Setup.nscans{nsub}{nses}-onset(n1),durat(min(length(durat),n1))));
                            in=offs+round(1+onset(n1)/rt+(0:tdurat/rt-1));
                            x(in(in>0))=1;
                        end
                        if CONN_x.Setup.acquisitiontype==1,
                            x=convn(x,hrf);
                        end
                        x=mean(reshape(x(offs+(1:10*CONN_x.Setup.nscans{nsub}{nses})),[10,CONN_x.Setup.nscans{nsub}{nses}]),1)';%x=x(1+10*(0:CONN_x.Setup.nscans{nsub}{nses}-1));
                        condname=['Effect of ',CONN_x.Setup.conditions.names{ncondition}];
                        data{nl1covariates+ncondition}=x;
                        names{nl1covariates+ncondition}=condname;
                        
                        idx1=find(x>0); idx2=[0;find(diff(idx1)>1);length(idx1)]; % crop to within-condition samples
                        if CONN_x.Setup.conditions.param(ncondition)>0
                            xw=sum(data{CONN_x.Setup.conditions.param(ncondition)},2); %note: assuming covariates already in BOLD-signal time-frame (e.g. realignment params, dynamic temporal states, etc.)
                            if 0,%CONN_x.Setup.acquisitiontype==1,
                                hrf=spm_hrf(conn_get_rt(nsub,nses));
                                xw=convn(xw,hrf,'same');
                            end
                            xw=xw.*x;
                        else
                            xw=x;%ones(size(x,1),1);
                        end
                        tsamples=[]; tweights=cell(1,5);
                        for n1=1:length(idx2)-1,
                            idx=idx1(idx2(n1)+1:idx2(n1+1));
                            tsamples=cat(1,tsamples,idx(:));
                            tweights{1}=cat(1,tweights{1},x(idx(:)));                       % HRF  %.*xw(idx(:)))
                            tweights{2}=cat(1,tweights{2},conn_hanning(length(idx)));       % HANNING %.*xw(idx(:)));
                            tweights{3}=cat(1,tweights{3},xw(idx(:)));                      % HRF.*MODUL
                            tweights{4}=cat(1,tweights{4},idx(:));                          % IDX
                            tweights{5}=cat(1,tweights{5},(1:numel(idx))');                 % JDX
                        end
                        if ~crop
                            for n1=1:5, tweights0=zeros(size(x,1),1); tweights0(tsamples)=tweights{n1}; tweights{n1}=tweights0; end
                            tsamples=(1:size(x,1))';
                        end
                        samples{ncondition}=tsamples;
                        weights{ncondition}=tweights;
                    end
                    n=n+1;
                    conn_waitbar(n/N,h,sprintf('Subject %d Session %d',nsub,nses));
                end
                filename=fullfile(filepath,['COV_Subject',num2str(nsub,'%03d'),'_Session',num2str(nses,'%03d'),'.mat']);
                %names=names([1:nl1covariates,nl1covariates+uniqueconditions]);
                %data=data([1:nl1covariates,nl1covariates+uniqueconditions]);
                [nill,idx]=unique(regexprep(names,' x FrequencyBand\d+$| x Time\d+$',''),'first');
                idx=sort(idx);
                idx=idx(cellfun('length',names(idx))>0);
                names=names(idx);
                data=data(idx);
                save(filename,'data','names');
                % 			if str2num(version('-release'))>=14, save(filename,'-V6','data','names');
                % 			else, save(filename,'data','names'); end
                for n1=1:length(data), if length(maxdims)<n1, maxdims(n1)=size(data{n1},2); else, maxdims(n1)=max(maxdims(n1),size(data{n1},2)); end; end
                names={CONN_x.Setup.conditions.names{1:nconditions}};
                filename=fullfile(filepath,['COND_Subject',num2str(nsub,'%03d'),'_Session',num2str(nses,'%03d'),'.mat']);
                save(filename,'samples','weights','names','crop');
                % 			if str2num(version('-release'))>=14, save(filename,'-V6','samples','weights','names');
                % 			else, save(filename,'samples','weights','names'); end
            end
		end
	end
	conn_waitbar('close',h);
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates DATA_Subject###_Session###.mat files (percentage signal change data at analysis mask)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(options==3) && any(CONN_x.Setup.steps([2,3])) && ~(isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'steps')&&~any(CONN_x.gui.steps([2,3])))
	[path,name,ext]=fileparts(CONN_x.filename);
    filepath=CONN_x.folders.data;
    try, conn_fileutils('mkdir',filepath); end
    if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'subjects'), validsubjects=CONN_x.gui.subjects; else validsubjects=1:CONN_x.Setup.nsubjects; end
    if isfield(CONN_x,'pobj')&&isstruct(CONN_x.pobj)&&isfield(CONN_x.pobj,'subjects'), validsubjects=CONN_x.pobj.subjects; end
	h=conn_waitbar(0,['Step ',num2str(sum(options<=3)),'/',num2str(length(options)),': Importing functional data']);
    REDO='Yes';filename=fullfile(filepath,['DATA_Subject',num2str(validsubjects(1),'%03d'),'_Session',num2str(1,'%03d'),'.mat']);
    if conn_existfile(filename),if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'overwrite'), REDO=CONN_x.gui.overwrite; else, REDO=conn_questdlg('Overwrite existing subject results?','','Yes', 'No', 'Yes');end; end
	if length(CONN_x.Setup.nsessions)==1, N=CONN_x.Setup.nsessions*numel(validsubjects); else N=sum(CONN_x.Setup.nsessions(validsubjects)); end
	n=0;
    if CONN_x.Setup.analysismask==1,%CONN_x.Setup.normalized,
        if isfield(CONN_x.Setup,'explicitmask'), filename=CONN_x.Setup.explicitmask{1};
        elseif CONN_x.Setup.spatialresolution~=4, filename=fullfile(fileparts(which(mfilename)),'utils','surf','mask.volume.brainmask.nii');
        else filename=fullfile(fileparts(which(mfilename)),'utils','surf','mask.surface.brainmask.nii');
        end
        Vmask=spm_vol(filename);
    else, Vmask=[]; end
	for nsub=validsubjects,
        switch CONN_x.Setup.spatialresolution
            case 1, 
                if isfield(CONN_x.Setup,'explicitmask'), filename=CONN_x.Setup.explicitmask{1};
                else filename=fullfile(fileparts(which(mfilename)),'utils','surf','mask.volume.brainmask.nii');
                end
            case 2, filename=deblank(CONN_x.Setup.structural{nsub}{1}{1});
            case 3, filename=deblank(CONN_x.Setup.functional{nsub}{1}{1}(1,:));
            case 4, filename=deblank(CONN_x.Setup.structural{nsub}{1}{1});
        end
        Vref=spm_vol(filename);
        Vref=Vref(1);
        CONN_x.Setup.spatialresolutionvolume{nsub}=Vref;

%         if ~CONN_x.Setup.normalized,
%             filename=deblank(CONN_x.Setup.structural{nsub}{1});
%             Vref=spm_vol(filename);
%         end
		nsess=CONN_x.Setup.nsessions(min(length(CONN_x.Setup.nsessions),nsub));
		sfile=[];%conn_prepend('',CONN_x.Setup.structural{nsub}{1},'_seg_inv_sn.mat');
		for nses=1:nsess,
			filename=fullfile(filepath,['DATA_Subject',num2str(nsub,'%03d'),'_Session',num2str(nses,'%03d'),'.mat']);
            if strcmp(lower(REDO),'yes')||~conn_existfile(filename),
                warning off;Vsource=spm_vol(CONN_x.Setup.functional{nsub}{nses}{1});warning on;
                CONN_x.Setup.nscans{nsub}{nses}=prod(size(Vsource));
                if CONN_x.Setup.analysismask==2&&nses==1,%~CONN_x.Setup.normalized&&(isempty(Vmask)), % computes analysis mask
                    Vmask=Vsource(1);
                    Vmask.pinfo=[1;0;0];
                    Vmask.fname=conn_prepend('mask_',CONN_x.Setup.functional{nsub}{1}{1}(1,:),'.nii');
                    Vmask.dt=[spm_type('uint8') spm_platform('bigend')];
                    a=ones(Vmask.dim(1:3));
                    [gridx,gridy,gridz]=ndgrid(1:Vmask.dim(1),1:Vmask.dim(2),1:Vmask.dim(3));xyz=Vmask.mat*[gridx(:),gridy(:),gridz(:),ones(numel(gridx),1)]';
                    for nsest=1:nsess,
                        if nsest==nses,Vsourcet=Vsource;else,warning off;Vsourcet=spm_vol(CONN_x.Setup.functional{nsub}{nsest}{1});warning on;end
                        for nt=1:numel(Vsourcet),
                            b=reshape(spm_get_data(Vsourcet(nt),pinv(Vsourcet(nt).mat)*xyz),Vmask.dim(1:3)); %b=spm_read_vols(Vsourcet(nt));
                            mb=mean(b(b>mean(b(~isnan(b)))/8));
                            a=a&(b>0.80*mb);
                        end
                    end
                    try,
                        Vmask=spm_write_vol(Vmask,a);
                    catch,
                        [temp_path,temp_filename,temp_ext]=fileparts(Vmask.fname);
                        Vmask.fname=fullfile(pwd,[temp_filename,temp_ext]);
                        Vmask=spm_write_vol(Vmask,a);
                    end
                end
                [filename,cache]=conn_tempcache(filename);
                V=conn_create_vol(filename,Vsource,[],Vref,sfile,Vmask,CONN_x.Setup.analysisunits==1,CONN_x.Setup.spatialresolution==4,0); %CONN_x.Setup.surfacesmoothing);
                conn_tempcache(cache,'matc');
            end
			n=n+1;
			conn_waitbar(n/N,h,sprintf('Subject %d Session %d',nsub,nses));
		end
	end
	conn_waitbar('close',h);
	clear data;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates ROI_Subject###_Session###.mat files (activation timecourses for each roi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(options==4) && any(CONN_x.Setup.steps([1,2,3,4])) && ~(isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'steps')&&~any(CONN_x.gui.steps([1,2,3]))),
	[path,name,ext]=fileparts(CONN_x.filename);
    filepath=CONN_x.folders.data;
    try, conn_fileutils('mkdir',filepath); end
    if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'subjects'), validsubjects=CONN_x.gui.subjects; else validsubjects=1:CONN_x.Setup.nsubjects; end
    if isfield(CONN_x,'pobj')&&isstruct(CONN_x.pobj)&&isfield(CONN_x.pobj,'subjects'), validsubjects=CONN_x.pobj.subjects; end
	h=conn_waitbar(0,['Step ',num2str(sum(options<=4)),'/',num2str(length(options)),': Importing ROI data']);
    REDO='Yes';filename=fullfile(filepath,['ROI_Subject',num2str(validsubjects(1),'%03d'),'_Session',num2str(1,'%03d'),'.mat']);
    if conn_existfile(filename),if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'overwrite'), REDO=CONN_x.gui.overwrite; else, REDO=conn_questdlg('Overwrite existing subject results?','','Yes', 'No', 'Yes');end; end
    USEEXPLICITMASK=true;
	nrois=length(CONN_x.Setup.rois.names)-1;
    allrois=1:nrois;
    if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'rois'), validrois=CONN_x.gui.rois; else validrois=allrois; end %find(arrayfun(@(n)~isempty(CONN_x.Setup.rois.files{1}{n}{1}{1}),1:nrois));
    refpial=[];
    %Vmask_ref_vol=fullfile(fileparts(which(mfilename)),'utils','surf','referenceGM.nii');
    if ~USEEXPLICITMASK, Vmask_ref_vol=fullfile(fileparts(which(mfilename)),'utils','surf','referenceGM.nii');
    elseif CONN_x.Setup.analysismask==1&&~isempty(regexp(CONN_x.Setup.explicitmask{1},'utils[\\\/]surf[\\\/]mask.volume.brainmask.nii$')), Vmask_ref_vol=fullfile(fileparts(which(mfilename)),'utils','surf','referenceGM.nii');
    else Vmask_ref_vol='';
    end
    Vmask_ref_surf=fullfile(fileparts(which(mfilename)),'utils','surf','mask.surface.brainmask.nii');
	if length(CONN_x.Setup.nsessions)==1, N=CONN_x.Setup.nsessions*numel(validsubjects); else N=sum(CONN_x.Setup.nsessions(validsubjects)); end
	N=N*nrois;
	n=0;
    all_tnames={};
	for nsub=validsubjects,
		nsess=CONN_x.Setup.nsessions(min(length(CONN_x.Setup.nsessions),nsub));
        clear Vmask; 
        for nroi=allrois, 
            if (nroi>3&&~CONN_x.Setup.rois.sessionspecific(nroi))||(nroi<=3&&~CONN_x.Setup.structural_sessionspecific), nsesstemp=1; else nsesstemp=nsess; end
            if nroi>3&&~CONN_x.Setup.rois.subjectspecific(nroi), nsubstemp=1; else nsubstemp=nsub; end
            for nses=1:nsesstemp,
                Vmask{nroi}{nses}=CONN_x.Setup.rois.files{nsubstemp}{nroi}{nses}{1};
            end
        end
        FORCERECOMPUTE=false;   % force recomputing eroded file
        SKIPRECOMPUTE=false;    % skip recomputing eroded file if it already exists
        QASTEP2=[0,1];          % QA volume variables: 0 before erosion; 1 after erosion
        DETRENDBEFOREPCA=true;
        if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'overwrite')&&strcmp(CONN_x.gui.overwrite,'No'), SKIPRECOMPUTE=true; end
        for nroi=intersect(validrois,1:3),
            THR=CONN_x.Setup.erosion.binary_threshold(nroi);
            THRTYPE=CONN_x.Setup.erosion.binary_threshold_type(nroi);
            THRGM=CONN_x.Setup.erosion.exclude_grey_matter(nroi);
            ERODE=CONN_x.Setup.erosion.erosion_steps(nroi);
            NEIGHB=CONN_x.Setup.erosion.erosion_neighb(nroi);
            if ~CONN_x.Setup.structural_sessionspecific, nsesstemp=1; else nsesstemp=nsess; end
            for nses=1:nsesstemp
                Vmask{nroi}{nses}=conn_prepend('e',CONN_x.Setup.rois.files{nsub}{nroi}{nses}{1});
                if FORCERECOMPUTE||(~SKIPRECOMPUTE&&strcmp(lower(REDO),'yes'))||~conn_existfile(Vmask{nroi}{nses}),
                    [nill,nill,ext,fnum]=spm_fileparts(Vmask{nroi}{nses});
                    switch(ext),
                        case {'.img','.nii','.hdr'},
                            V0=spm_vol(CONN_x.Setup.rois.files{nsub}{nroi}{nses}{1}); % mask
                            assert(numel(V0)==1,'%s ROI should contain a single 3d image (found %d)',CONN_x.Setup.rois.names{nroi},numel(V0));
                            X0=spm_read_vols(V0);
                            if nroi>1&&~isempty(THRGM)&&~isnan(THRGM)
                                V2=spm_vol(CONN_x.Setup.rois.files{nsub}{1}{nses}{1}); % grey matter exclusion mask
                                [tx,ty,tz]=ndgrid(1:V0(1).dim(1),1:V0(1).dim(2),1:V0(1).dim(3));
                                txyz=V0(1).mat*[tx(:) ty(:) tz(:) ones(numel(tx),1)]';
                                X2=spm_get_data(V2,pinv(V2(1).mat)*txyz);
                                X0(X2>THRGM)=nan;
                            end
                            if THRTYPE==1, tTHR=THR;
                                idx1=find(X0(:)>tTHR);
                            else
                                sX0=sort(X0(~isnan(X0)&X0~=0));
                                tTHR=sX0(max(1,min(numel(sX0),round(THR*numel(sX0)))));
                                idx1=find(X0(:)>=tTHR);
                            end
                            Nstep1=numel(idx1);
                            if isempty(idx1), error(sprintf('No suprathreshold voxels in ROI file %s (this typically indicates a problem during normalization/segmentation for this subject; please try re-running normalization/segmentation step after manually re-aligning the structural volumes to better match the default template orientation)',CONN_x.Setup.rois.files{nsub}{nroi}{nses}{1}));
                            end
                            if ERODE>0
                                finished=false;
                                if rem(ERODE,1), tERODE=1;
                                else             tERODE=ERODE;
                                end
                                while ~finished
                                    [idxx,idxy,idxz]=ind2sub(size(X0),idx1);
                                    idxt=find(idxx>tERODE&idxx<size(X0,1)+1-tERODE&idxy>tERODE&idxy<size(X0,2)+1-tERODE&idxz>tERODE&idxz<size(X0,3)+1-tERODE);
                                    [idxdx,idxdy,idxdz]=ndgrid(-tERODE:tERODE);
                                    idxd=idxdx+size(X0,1)*idxdy+size(X0,1)*size(X0,2)*idxdz;
                                    tidx1=idx1(idxt);
                                    valt=zeros(size(tidx1));
                                    for n1=1:numel(idxd), valt=valt+(X0(tidx1+idxd(n1))<tTHR); end %for n1=1:length(idxt), valt(n1)=sum(sum(sum(X0(idxx(idxt(n1))+(-tERODE:tERODE),idxy(idxt(n1))+(-tERODE:tERODE),idxz(idxt(n1))+(-tERODE:tERODE))<tTHR,3),2),1); end
                                    if rem(ERODE,1),
                                        svalt=sort(valt);
                                        svalt=svalt(min(numel(svalt),round(ERODE*Nstep1)));
                                        finished=svalt>0;
                                    else svalt=NEIGHB;
                                        finished=true;
                                    end
                                    idxt=idxt(valt<=svalt);
                                    idx1=idx1(idxt);
                                    tERODE=tERODE+1;
                                end
                                Nstep2=numel(idx1);
%                                 [idxx,idxy,idxz]=ind2sub(size(X0),idx1);
%                                 idxt=find(idxx>ERODE&idxx<size(X0,1)+1-ERODE&idxy>ERODE&idxy<size(X0,2)+1-ERODE&idxz>ERODE&idxz<size(X0,3)+1-ERODE);
%                                 for n1=1:length(idxt), if (sum(sum(sum(X0(idxx(idxt(n1))+(-ERODE:ERODE),idxy(idxt(n1))+(-ERODE:ERODE),idxz(idxt(n1))+(-ERODE:ERODE))<tTHR,3),2),1))>NEIGHB, idxt(n1)=0; end; end
%                                 idxt=idxt(idxt>0);
%                                 idx1=idx1(idxt);
%                                 Nstep2=numel(idx1);
                            else Nstep2=Nstep1;
                            end
                            X1=zeros(size(X0));X1(idx1)=1;
                            V0.fname=regexprep(Vmask{nroi}{nses},',\d+$','');
                            spm_write_vol(V0,X1);
                            if ~Nstep2, error(sprintf('No suprathreshold voxels in ROI file %s after erosion (original file contained %d suprathreshold voxels)',CONN_x.Setup.rois.files{nsub}{nroi}{nses}{1},Nstep1));
                            end
                            for qastep2=QASTEP2,
                                roiname=regexprep(CONN_x.Setup.rois.names{nroi},'\s','');
                                if qastep2, roiname=[roiname,'_eroded']; end
                                if ~CONN_x.Setup.structural_sessionspecific, qa_name=['QC_',roiname,'_vol'];
                                else qa_name=['QC_',roiname,'_vol_session',num2str(nses)];
                                end
                                qa_icov=find(strcmp(qa_name,CONN_x.Setup.l2covariates.names(1:end-1)),1);
                                if isempty(qa_icov),
                                    qa_icov=numel(CONN_x.Setup.l2covariates.names);
                                    CONN_x.Setup.l2covariates.names{qa_icov}=qa_name;
                                    CONN_x.Setup.l2covariates.descrip{qa_icov}=['CONN Quality Assurance: # of voxels in ',roiname];
                                    CONN_x.Setup.l2covariates.names{qa_icov+1}=' ';
                                    for tnsub=1:CONN_x.Setup.nsubjects, CONN_x.Setup.l2covariates.values{tnsub}{qa_icov}=nan; end
                                end
                                if qastep2, CONN_x.Setup.l2covariates.values{nsub}{qa_icov}=Nstep2;
                                else CONN_x.Setup.l2covariates.values{nsub}{qa_icov}=Nstep1;
                                end
                            end
                        otherwise,
                            Vmask{nroi}{nses}=CONN_x.Setup.rois.files{nsub}{nroi}{nses}{1};
                    end
                end
            end
        end
        
		for nses=1:nsess,
            sampledatachange=false;
            anychange=false;
			filename=fullfile(filepath,['ROI_Subject',num2str(nsub,'%03d'),'_Session',num2str(nses,'%03d'),'.mat']);
            if conn_existfile(filename),
                old=load(filename);%,'data','names','source','xyz','sampledata','samplexyz');
            else, old=[]; end
            
            filename=fullfile(filepath,['COV_Subject',num2str(nsub,'%03d'),'_Session',num2str(nses,'%03d'),'.mat']);
            covariates=load(filename);
            defaultcov=~cellfun('length',regexp(covariates.names,'^QA_|^QC_')); % default subset of covariates (includes everything but the ones listed here) : note:default copied from step=5 code
            covariates=cat(2,covariates.data{defaultcov});
            if DETRENDBEFOREPCA, covariates=cat(2,linspace(-1,1,size(covariates,1))',covariates); end
            source=fullfile(filepath,['DATA_Subject',num2str(nsub,'%03d'),'_Session',num2str(nses,'%03d'),'.mat']);
            Vsource=CONN_x.Setup.functional{nsub}{nses}{1};
                clear VsourceUnsmoothed;
                for nalt=1:numel(CONN_x.Setup.secondarydataset)
                    VsourceUnsmoothed{nalt}=conn_get_functional(nsub,nses,nalt,true);
                end
            issurface=conn_surf_dimscheck(CONN_x.Setup.functional{nsub}{nses}{3}); % note: when sampling gray-matter data for Denoising plots use Primary Dataset if surface-level
            if issurface, Vmask_ref=Vmask_ref_surf;
            else Vmask_ref=Vmask_ref_vol;
            end
            if isempty(Vmask_ref)
                if ~CONN_x.Setup.structural_sessionspecific, nsesstemp=1; else nsesstemp=nses; end
                Vmask_ref=Vmask{1}{min(nses,nsesstemp)};
                if nsess==1, conn_disp('warning: 1000-node network sample for histogram displays created using subject-specific grey-matter mask voxels. Quality Control QC-FC correlation plots may be inaccurate (if grey-matter masks are different across subjects). Enter a subject-independent Grey Matter ROI in order to avoid this warning'); end
                if ~issurface&&CONN_x.Setup.rois.unsmoothedvolumes(1), Vsourcethis=VsourceUnsmoothed{CONN_x.Setup.rois.unsmoothedvolumes(1)}; else, Vsourcethis=Vsource; end
            else
                Vsourcethis=Vsource; %if ~issurface&&CONN_x.Setup.rois.unsmoothedvolumes(1), Vsourcethis=VsourceUnsmoothed{CONN_x.Setup.rois.unsmoothedvolumes(1)}; else, Vsourcethis=Vsource; end
            end
            if CONN_x.Setup.analysisunits==1, scalinglevel='roi'; else scalinglevel='none'; end
            if any(ismember(validrois,[0,1]))
                if strcmp(lower(REDO),'yes')||isempty(old)||~isfield(old,'sampledata')||size(old.sampledata,2)<1000,%||~isfield(old,'samplexyz'),
                    try,
                        [sampledata,samplexyz]=conn_rex(Vsourcethis,Vmask_ref,'dims',1000,'roi_threshold',.25,'level','subsetvoxels','scaling',scalinglevel,'select_clusters',0,'output_type','none');
                        if issurface, samplexyz=num2cell(conn_surf_coords(cell2mat(samplexyz)),1); end
                        %[sampledata,samplexyz]=conn_rex(Vsource,Vmask{1}{min(nses,numel(Vmask{1}))},'dims',512,'level','subsetvoxels','scaling',scalinglevel,'select_clusters',0,'output_type','none');
                    catch
                        error('Error extracting Grey Matter BOLD signal for subject %d session %d',nsub,nses);
                    end
                    sampledatachange=true;
                else sampledata=old.sampledata; if isfield(old,'samplexyz'), samplexyz=old.samplexyz; else samplexyz=[]; end; 
                end
                if isempty(sampledata), error('Empty grey-matter mask in subject %d',nsub); end
            elseif isfield(old,'sampledata'), sampledata=old.sampledata; if isfield(old,'samplexyz'), samplexyz=old.samplexyz; else samplexyz=[]; end;
            else sampledata=zeros(CONN_x.Setup.nscans{nsub}{nses},0); samplexyz=[];
            end
            clear data names xyz;
            nroi1=1;
            for nroi=allrois,
                if ~isfield(CONN_x.Setup.rois,'mask') || length(CONN_x.Setup.rois.mask)<nroi, CONN_x.Setup.rois.mask(nroi)=0; end
                if ~isfield(CONN_x.Setup.rois,'subjectspecific') || length(CONN_x.Setup.rois.subjectspecific)<nroi, CONN_x.Setup.rois.subjectspecific(nroi)=1; end
                if ~isfield(CONN_x.Setup.rois,'sessionspecific') || length(CONN_x.Setup.rois.sessionspecific)<nroi, CONN_x.Setup.rois.sessionspecific(nroi)=0; end
                if ~isfield(CONN_x.Setup.rois,'multiplelabels') || length(CONN_x.Setup.rois.multiplelabels)<nroi, CONN_x.Setup.rois.multiplelabels(nroi)=0; end
                if ~isfield(CONN_x.Setup.rois,'regresscovariates') || length(CONN_x.Setup.rois.regresscovariates)<nroi, CONN_x.Setup.rois.regresscovariates(nroi)=double(CONN_x.Setup.rois.dimensions{nroi}>1); end
                if ~isfield(CONN_x.Setup.rois,'unsmoothedvolumes') || length(CONN_x.Setup.rois.unsmoothedvolumes)<nroi, CONN_x.Setup.rois.unsmoothedvolumes(nroi)=1; end
                if ~isfield(CONN_x.Setup.rois,'weighted') || length(CONN_x.Setup.rois.weighted)<nroi, CONN_x.Setup.rois.weighted(nroi)=double(CONN_x.Setup.rois.dimensions{nroi}==0); end
                if (nroi>3&&~CONN_x.Setup.rois.sessionspecific(nroi))||(nroi<=3&&~CONN_x.Setup.structural_sessionspecific), nsesstemp=1; else nsesstemp=nsess; end
                if nroi>3&&~CONN_x.Setup.rois.subjectspecific(nroi), nsubstemp=1; else nsubstemp=nsub; end
                
                failed=0;
                isredo=strcmp(lower(REDO),'yes')&&ismember(nroi,validrois);
                if ~isredo&&~isempty(old),
                    if CONN_x.Setup.rois.multiplelabels(nroi), 
                        [roi_path_dir,roi_path_name]=fileparts(Vmask{nroi}{min(nses,nsesstemp)});
                        idx=strmatch([CONN_x.Setup.rois.names{nroi},'.'],old.names);
                        %idx=strmatch([roi_path_name,'.'],old.names);
                        if isempty(idx), failed=1;
                        else
                            try, tnames=all_tnames{nsubstemp,min(nses,nsesstemp),nroi};
                            catch, tnames={};
                            end
                            if isempty(tnames)
                                if 1,%isempty(dir(fullfile(roi_path_dir,[roi_path_name,'.txt']))),
                                    if CONN_x.Setup.rois.mask(nroi), mask=CONN_x.Setup.rois.files{nsubstemp}{1}{min(nses,nsesstemp)}{1}; else, mask=''; end
                                    if CONN_x.Setup.rois.multiplelabels(nroi)&&numel(CONN_x.Setup.rois.files{nsubstemp}{nroi}{min(nses,nsesstemp)}{3})<=1, level='clusters'; else, level='rois'; end
                                    if CONN_x.Setup.rois.regresscovariates(nroi), entercovariates=covariates; else, entercovariates=[]; end
                                    if 0,%CONN_x.Setup.rois.dimensions{nroi}>1, [nill,tnames]=conn_rex(Vmask{nroi}{min(nses,nsesstemp)},Vmask{nroi}{min(nses,nsesstemp)},'summary_measure','eigenvariate','disregard_zeros',0,'conjunction_mask',mask,'level',level,'scaling',scalinglevel,'select_clusters',0,'covariates',entercovariates,'output_type','none');
                                    else [nill,tnames]=conn_rex(Vmask{nroi}{min(nses,nsesstemp)},Vmask{nroi}{min(nses,nsesstemp)},'summary_measure','mean','disregard_zeros',0,'conjunction_mask',mask,'level',level,'scaling',scalinglevel,'select_clusters',0,'covariates',entercovariates,'output_type','none');
                                    end
                                    for n1=1:numel(tnames)
                                        if ~isempty(strmatch(roi_path_name,tnames{n1})), tnames{n1}=[tnames{n1}(numel(roi_path_name)+2:end)]; end
                                    end
                                else
                                    tnames=textread(fullfile(roi_path_dir,[roi_path_name,'.txt']),'%s','delimiter','\n');
                                end
                            end
                            all_tnames{nsubstemp,min(nses,nsesstemp),nroi}=tnames;
                            clear tidx; for n1=1:length(tnames),
                                tidx{n1}=find(strcmp([CONN_x.Setup.rois.names{nroi},'.',tnames{n1}],old.names));
                                if isempty(tidx{n1}), failed=1; else, tidx{n1}=tidx{n1}(1); end
                            end
                            if ~failed,
                                for n1=1:length(tnames),
                                    data{nroi1}=old.data{tidx{n1}};
                                    names{nroi1}=old.names{tidx{n1}};
                                    xyz{nroi1}=old.xyz{tidx{n1}};
                                    nroi1=nroi1+1;
                                end
                            end
                        end
                    else
                        idx=find(strcmp(CONN_x.Setup.rois.names{nroi},old.names));
                        if isempty(idx), failed=1;
                        else, 
                            idx=idx(1);
                            data{nroi1}=old.data{idx};
                            names{nroi1}=old.names{idx};
                            xyz{nroi1}=old.xyz{idx};
                            nroi1=nroi1+1;
                        end
                    end
                end
                if isredo||isempty(old)||failed
                    anychange=true;
                    %disp(nroi);
                    if CONN_x.Setup.rois.mask(nroi), 
                        if ~CONN_x.Setup.structural_sessionspecific, nsesstemp1=1; else nsesstemp1=nsess; end
                        mask=conn_prepend('e',CONN_x.Setup.rois.files{nsub}{1}{min(nses,nsesstemp1)}{1});
                        %mask=CONN_x.Setup.rois.files{nsub}{1}{min(nses,nsesstemp1)}{1};
                    else mask='';
                    end
                    if CONN_x.Setup.rois.multiplelabels(nroi)&&numel(CONN_x.Setup.rois.files{nsubstemp}{nroi}{min(nses,nsesstemp)}{3})<=1, level='clusters'; else, level='rois'; end
                    if CONN_x.Setup.rois.regresscovariates(nroi), entercovariates=covariates; else, entercovariates=[]; end
                    if CONN_x.Setup.rois.unsmoothedvolumes(nroi), Vsourcethis=VsourceUnsmoothed{CONN_x.Setup.rois.unsmoothedvolumes(nroi)}; else, Vsourcethis=Vsource; end
                    if conn_surf_dimscheck(CONN_x.Setup.rois.files{nsubstemp}{nroi}{min(nses,nsesstemp)}{3})&&~conn_surf_dimscheck(Vsourcethis), fsanatomical=CONN_x.Setup.structural{nsub}{min(nses,nsesstemp)}{1}; else fsanatomical=''; end
                    if isfield(CONN_x.Setup,'outputfiles')&&numel(CONN_x.Setup.outputfiles)>=6&&CONN_x.Setup.outputfiles(6), outputtype='saverex'; filenamerex=['REX_Subject',num2str(nsub,'%03d'),'_Session',num2str(nses,'%03d'),'_ROI',num2str(nroi),'.mat'];
                    else outputtype='none'; filenamerex='';
                    end
                    if CONN_x.Setup.rois.dimensions{nroi}>1&&isfield(CONN_x.Setup,'extractSVD')&&CONN_x.Setup.extractSVD % SVD instead of PCA
                        if CONN_x.Setup.rois.weighted(nroi), [data{nroi1},namesroi{nroi},params]=conn_rex(Vsourcethis,Vmask{nroi}{min(nses,nsesstemp)},'summary_measure','weighted eigenvariate','PCA',0,'dims',CONN_x.Setup.rois.dimensions{nroi},'conjunction_mask',mask,'level',level,'scaling','none','select_clusters',0,'covariates',entercovariates,'fsanatomical',fsanatomical,'output_type',outputtype,'output_rex',filenamerex,'output_folder',filepath);
                        else [data{nroi1},namesroi{nroi},params]=conn_rex(Vsourcethis,Vmask{nroi}{min(nses,nsesstemp)},'summary_measure','eigenvariate','PCA',0,'dims',CONN_x.Setup.rois.dimensions{nroi},'conjunction_mask',mask,'level',level,'scaling',scalinglevel,'select_clusters',0,'covariates',entercovariates,'fsanatomical',fsanatomical,'output_type',outputtype,'output_rex',filenamerex,'output_folder',filepath);
                        end
                        if size(data{nroi1},2)<CONN_x.Setup.rois.dimensions{nroi}, 
                            conn_disp('fprintf','WARNING: not enough voxels or scans to extract %d dimensions from %s @ %s\n',CONN_x.Setup.rois.dimensions{nroi},Vsourcethis,Vmask{nroi}{min(nses,nsesstemp)});
                            data{nroi1}=[data{nroi1} zeros(size(data{nroi1},1),CONN_x.Setup.rois.dimensions{nroi}-size(data{nroi1},2))]; 
                        end
                    else
                        if CONN_x.Setup.rois.dimensions{nroi}>1,        % average&pca
                            if CONN_x.Setup.rois.weighted(nroi), [data{nroi1},namesroi{nroi},params]=conn_rex(Vsourcethis,Vmask{nroi}{min(nses,nsesstemp)},'summary_measure','weighted eigenvariate','dims',CONN_x.Setup.rois.dimensions{nroi},'conjunction_mask',mask,'level',level,'scaling','none','select_clusters',0,'covariates',entercovariates,'fsanatomical',fsanatomical,'output_type',outputtype,'output_rex',filenamerex,'output_folder',filepath);
                            else [data{nroi1},namesroi{nroi},params]=conn_rex(Vsourcethis,Vmask{nroi}{min(nses,nsesstemp)},'summary_measure','eigenvariate','dims',CONN_x.Setup.rois.dimensions{nroi},'conjunction_mask',mask,'level',level,'scaling',scalinglevel,'select_clusters',0,'covariates',entercovariates,'fsanatomical',fsanatomical,'output_type',outputtype,'output_rex',filenamerex,'output_folder',filepath);
                            end
                            if size(data{nroi1},2)<CONN_x.Setup.rois.dimensions{nroi},
                                conn_disp('fprintf','WARNING: not enough voxels or scans to extract %d dimensions from %s @ %s\n',CONN_x.Setup.rois.dimensions{nroi},Vsourcethis,Vmask{nroi}{min(nses,nsesstemp)});
                                data{nroi1}=[data{nroi1} zeros(size(data{nroi1},1),CONN_x.Setup.rois.dimensions{nroi}-size(data{nroi1},2))];
                            end
                        elseif CONN_x.Setup.rois.dimensions{nroi}==0||CONN_x.Setup.rois.weighted(nroi),   % weighted sum
                            [data{nroi1},namesroi{nroi},params]=conn_rex(Vsourcethis,Vmask{nroi}{min(nses,nsesstemp)},'summary_measure','weighted sum','conjunction_mask',mask,'level',level,'scaling','none','select_clusters',0,'covariates',entercovariates,'fsanatomical',fsanatomical,'output_type',outputtype,'output_rex',filenamerex,'output_folder',filepath); 
                        else                                            % average
                            [data{nroi1},namesroi{nroi},params]=conn_rex(Vsourcethis,Vmask{nroi}{min(nses,nsesstemp)},'summary_measure','mean','conjunction_mask',mask,'level',level,'scaling',scalinglevel,'select_clusters',0,'covariates',entercovariates,'fsanatomical',fsanatomical,'output_type',outputtype,'output_rex',filenamerex,'output_folder',filepath); 
                        end
                    end
                    [data{nroi1},ok]=conn_nan(data{nroi1});
                    data{nroi1}=detrend(data{nroi1},'constant');
                    if ~isempty(fsanatomical)&&isempty(refpial), refpial=conn_surf_readsurf(fullfile(fileparts(which(mfilename)),'utils','surf','pial.surf')); end
                    if CONN_x.Setup.rois.multiplelabels(nroi),
                        datatemp=data{nroi1};
                        [Vmaskfilepath,Vmaskfilename,Vmaskfileext]=fileparts(Vmask{nroi}{min(nses,nsesstemp)});
                        for n1=1:max(1,CONN_x.Setup.rois.dimensions{nroi}):size(datatemp,2),
                            data{nroi1}=datatemp(:,n1+(0:max(0,CONN_x.Setup.rois.dimensions{nroi}-1))); 
                            names{nroi1}=namesroi{nroi}{n1};
                            if ~isempty(strmatch(Vmaskfilename,names{nroi1})), names{nroi1}=[CONN_x.Setup.rois.names{nroi},names{nroi1}(numel(Vmaskfilename)+1:end)]; end
                            if any(strcmp(names{nroi1},names(1:nroi1-1))), error('duplicated ROI name %s',names{nroi1}); end
                            Z=params.ROIinfo.basis{params.ROIinfo.trans{n1}{1}}{params.ROIinfo.trans{n1}{2}}(params.ROIinfo.trans{n1}{4},params.ROIinfo.trans{n1}{3});
                            XYZ=params.ROIinfo.voxels{params.ROIinfo.trans{n1}{1}}{params.ROIinfo.trans{n1}{2}}(params.ROIinfo.trans{n1}{4},:);
                            if ~isempty(fsanatomical)
                                sroi=CONN_x.Setup.rois.files{nsubstemp}{nroi}{min(nses,nsesstemp)}{3}(1).dim;
                                psroi=prod(sroi);
                                tXYZ=1+(XYZ-1)*[1;sroi(1);sroi(1)*sroi(2)];
                                them=tXYZ<=psroi/2;
                                XYZ(them,:)=refpial(1).vertices(tXYZ(them),:);
                                them=tXYZ>psroi/2;
                                XYZ(them,:)=refpial(2).vertices(tXYZ(them)-psroi/2,:);
                            end
                            xyz{nroi1}=mean(XYZ(Z==max(Z),:),1);
                            nroi1=nroi1+1;
                        end
                    else,
                        names{nroi1}=CONN_x.Setup.rois.names{nroi};
                        if isempty(params.ROIinfo.trans)
                            error('Empty ROI mask in Subject %d Session %d ROI %d = %s',nsub,nses,nroi,names{nroi1});
                        end
                        n1=1;
                        Z=params.ROIinfo.basis{params.ROIinfo.trans{n1}{1}}{params.ROIinfo.trans{n1}{2}}(params.ROIinfo.trans{n1}{4},params.ROIinfo.trans{n1}{3});
                        XYZ=params.ROIinfo.voxels{params.ROIinfo.trans{n1}{1}}{params.ROIinfo.trans{n1}{2}}(params.ROIinfo.trans{n1}{4},:);
                        if ~isempty(fsanatomical)
                            sroi=CONN_x.Setup.rois.files{nsubstemp}{nroi}{min(nses,nsesstemp)}{3}(1).dim;
                            psroi=prod(sroi);
                            tXYZ=1+(XYZ-1)*[1;sroi(1);sroi(1)*sroi(2)];
                            them=tXYZ<=psroi/2;
                            XYZ(them,:)=refpial(1).vertices(tXYZ(them),:);
                            them=tXYZ>psroi/2;
                            XYZ(them,:)=refpial(2).vertices(tXYZ(them)-psroi/2,:);
                        end
                        xyz{nroi1}=mean(XYZ(Z==max(Z),:),1);
                        nroi1=nroi1+1;
                    end
%                     if isfield(CONN_x.Setup,'outputfiles')&&numel(CONN_x.Setup.outputfiles)>=6&&CONN_x.Setup.outputfiles(6),
%                         filename=fullfile(filepath,['REX_Subject',num2str(nsub,'%03d'),'_Session',num2str(nses,'%03d'),'_ROI',num2str(nroi),'.mat']);
%                         load(fullfile(filepath,'REX.mat'),'params');save(filename,'params');
%                     end
                end
                
                n=n+1;
                conn_waitbar(n/N,h,sprintf('Subject %d Session %d',nsub,nses));
            end
            if anychange||isempty(old)||~isequal(old.names,names)
                filename=fullfile(filepath,['ROI_Subject',num2str(nsub,'%03d'),'_Session',num2str(nses,'%03d'),'.mat']);
                save(filename,'data','names','source','xyz','sampledata','samplexyz');
            elseif sampledatachange,
                filename=fullfile(filepath,['ROI_Subject',num2str(nsub,'%03d'),'_Session',num2str(nses,'%03d'),'.mat']);
                save(filename,'sampledata','samplexyz','-append');
            end
%             if str2num(version('-release'))>=14, save(filename,'-V6','data','names','source','xyz','sampledata');
%             else, save(filename,'data','names','source','xyz','sampledata');end
		end
    end
	conn_waitbar('close',h);
end    
    
% double-checks same ROIs across all subjects
if any(options==4.5) && any(CONN_x.Setup.steps([1,2,3,4])) && ~(isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'steps')&&~any(CONN_x.gui.steps([1,2,3]))),
    % eliminates ROIs not present in all subjects/sessions
    filepath=CONN_x.folders.data;
    validsubjects=1:CONN_x.Setup.nsubjects; %if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'subjects'), validsubjects=CONN_x.gui.subjects; else validsubjects=1:CONN_x.Setup.nsubjects; end
    if isfield(CONN_x,'pobj')&&isstruct(CONN_x.pobj)&&isfield(CONN_x.pobj,'subjects'), validsubjects=CONN_x.pobj.subjects; conn_projectmanager('addstep',4.5); end
	h=conn_waitbar(0,['Step ',num2str(sum(options<=4.5)),'/',num2str(length(options)),': Checking ROI data consistency across subjects']);
    n=0;N=numel(validsubjects);
    anyissue=false;
    ROInamesall={};
	for nsub=validsubjects,
		nsess=CONN_x.Setup.nsessions(min(length(CONN_x.Setup.nsessions),nsub));
		for nses=1:nsess,
			filename=fullfile(filepath,['ROI_Subject',num2str(nsub,'%03d'),'_Session',num2str(nses,'%03d'),'.mat']);
            if conn_existfile(filename)
                load(filename,'names');
                if nsub==validsubjects(1)&&nses==1,ROInamesall=names;
                else
                    if ~anyissue, anyissue=~isequal(ROInamesall,names); end
                    temp=zeros(1,length(ROInamesall));
                    for n1=1:length(ROInamesall),
                        idx=find(strcmp(ROInamesall{n1},names));
                        temp(n1)=~isempty(idx);
                    end
                    if ~all(temp),
                        idx=find(~temp);
                        for n1=1:length(idx), conn_disp(['Warning: ROI ',ROInamesall{idx(n1)},' missing from Subject ',num2str(nsub),' Session ',num2str(nses)]); end
                        anyissue=true;
                    end
                    ROInamesall={ROInamesall{find(temp)}};
                end
            end
        end
        n=n+.5;
        conn_waitbar(n/N,h,sprintf('Subject %d',nsub));
    end
    if anyissue
        if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'overwrite'), tmp='Yes'; else tmp=[]; end
        for nsub=validsubjects,
            nsess=CONN_x.Setup.nsessions(min(length(CONN_x.Setup.nsessions),nsub));
            for nses=1:nsess,
                filename=fullfile(filepath,['ROI_Subject',num2str(nsub,'%03d'),'_Session',num2str(nses,'%03d'),'.mat']);
                if conn_existfile(filename)
                    load(filename,'data','names','source','xyz','sampledata','samplexyz');
                    temp=zeros(1,length(names));
                    for n1=1:length(names),
                        idx=find(strcmp(names{n1},ROInamesall));
                        temp(n1)=~isempty(idx);
                    end
                    if ~all(temp),
                        if isempty(tmp), tmp=conn_questdlg('Remove regions with incomplete data?','','Yes', 'No', 'Yes'); end
                        if strcmp(tmp,'Yes'),
                            data={data{find(temp)}};
                            names={names{find(temp)}};
                            xyz={xyz{find(temp)}};
                            save(filename,'data','names','source','xyz','sampledata','samplexyz');
                        end
                    end
                end
            end
            n=n+.5;
            conn_waitbar(n/N,h,sprintf('Subject %d',nsub));
        end
    else
        n=n+.5*numel(validsubjects);
        conn_waitbar(n/N,h,sprintf('Subject %d',nsub));
    end
	conn_waitbar('close',h);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% updates CONN_x.Preproc structure with a subset of ROIs plus all covariates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(options==5),
	[path,name,ext]=fileparts(CONN_x.filename);
    filepath=CONN_x.folders.data;
    validsubjects=1:CONN_x.Setup.nsubjects; %if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'subjects'), validsubjects=CONN_x.gui.subjects; else validsubjects=1:CONN_x.Setup.nsubjects; end
    if isfield(CONN_x,'pobj')&&isstruct(CONN_x.pobj)&&isfield(CONN_x.pobj,'subjects'), validsubjects=CONN_x.pobj.subjects; conn_projectmanager('addstep',5); end
	h=conn_waitbar(0,['Step ',num2str(sum(options<=5)),'/',num2str(length(options)),': Updating Denoising variables']);
	filename1=fullfile(filepath,['ROI_Subject',num2str(validsubjects(1),'%03d'),'_Session',num2str(1,'%03d'),'.mat']);
	filename2=fullfile(filepath,['COV_Subject',num2str(validsubjects(1),'%03d'),'_Session',num2str(1,'%03d'),'.mat']);
    if ~conn_existfile(filename1), conn_disp(['Not ready to process step conn_process_5']); return; end
	x1=conn_loadmatfile(filename1);
	x2=conn_loadmatfile(filename2);
	CONN_x.Preproc.variables.names=cat(2,x1.names,x2.names);
	CONN_x.Preproc.variables.types=cat(2,repmat({'roi'},[1,length(x1.names)]),repmat({'cov'},[1,length(x2.names)]));
	%CONN_x.Preproc.variables.deriv=cat(2,repmat({0},[1,length(x1.names)]),repmat({1},[1,length(x2.names)]));
    %[CONN_x.Preproc.variables.deriv{strcmp(CONN_x.Preproc.variables.names,'scrubbing')>0}]=deal(0);
	CONN_x.Preproc.variables.deriv=cat(2,repmat({0},[1,length(x1.names)]),cellfun(@(x)~isempty(x),regexpi(x2.names,'^effect of|realign|movement|motion'),'uni',0));
	CONN_x.Preproc.variables.power=cat(2,repmat({1},[1,length(x1.names)]),repmat({1},[1,length(x2.names)]));
	CONN_x.Preproc.variables.filter=cat(2,repmat({0},[1,length(x1.names)]),repmat({0},[1,length(x2.names)]));
	CONN_x.Preproc.variables.dimensions={};
	for n1=1:length(x1.names), CONN_x.Preproc.variables.dimensions{end+1}=size(x1.data{n1},2)*ones(1,2); end
	for n1=1:length(x2.names), CONN_x.Preproc.variables.dimensions{end+1}=size(x2.data{n1},2)*ones(1,2); end
    N1=length(x1.names);N2=length(x2.names);
	for nsub=validsubjects,
		nsess=CONN_x.Setup.nsessions(min(length(CONN_x.Setup.nsessions),nsub));
		for nses=1:nsess,
            filename2=fullfile(filepath,['COV_Subject',num2str(nsub,'%03d'),'_Session',num2str(nses,'%03d'),'.mat']);
            if conn_existfile(filename2)
                x2=conn_loadmatfile(filename2,'data','names');
                for n1=1:N2, CONN_x.Preproc.variables.dimensions{N1+n1}=max(CONN_x.Preproc.variables.dimensions{N1+n1},size(x2.data{n1},2)*ones(1,2)); end
            end
        end
    end
    conn_waitbar(1/2,h);

    defaultcov=~cellfun('length',regexp(x2.names,'^QA_|^QC_')); % default subset of covariates (includes everything but the ones listed here)
    defaultroi={'White Matter','CSF','MotionMask'}; % default subset of rois (includes only those present among the ones listed here)
    if isfield(CONN_x.Preproc.confounds,'names') && ~isempty(CONN_x.Preproc.confounds.names), initial=CONN_x.Preproc.confounds.names; dims=CONN_x.Preproc.confounds.dimensions; ders=CONN_x.Preproc.confounds.deriv; pows=CONN_x.Preproc.confounds.power; filts=CONN_x.Preproc.confounds.filter; 
    else, initial={defaultroi{:},x2.names{defaultcov}};dims={5,5}; ders={}; pows={}; filts={}; end
	CONN_x.Preproc.confounds.names={};
	CONN_x.Preproc.confounds.types={};
	CONN_x.Preproc.confounds.power={};
	CONN_x.Preproc.confounds.deriv={};
	CONN_x.Preproc.confounds.filter={};
	CONN_x.Preproc.confounds.dimensions={};
	for n1=1:length(initial), 
        idx=strmatch(initial{n1},CONN_x.Preproc.variables.names,'exact'); 
        if isempty(idx),
            idx=strmatch(initial{n1},CONN_x.Preproc.variables.names); % allows partial-name matches
            if numel(idx)~=1, idx=[]; end
        end
        if ~isempty(idx),
			CONN_x.Preproc.confounds.names{end+1}=CONN_x.Preproc.variables.names{idx};
			CONN_x.Preproc.confounds.types{end+1}=CONN_x.Preproc.variables.types{idx};
			if length(ders)>=n1&&~isempty(ders{n1}), CONN_x.Preproc.confounds.deriv{end+1}=ders{n1}; 
            else CONN_x.Preproc.confounds.deriv{end+1}=CONN_x.Preproc.variables.deriv{idx};
            end
			if length(pows)>=n1&&~isempty(pows{n1}), CONN_x.Preproc.confounds.power{end+1}=pows{n1}; 
            else CONN_x.Preproc.confounds.power{end+1}=CONN_x.Preproc.variables.power{idx};
            end
			if length(filts)>=n1&&~isempty(filts{n1}), CONN_x.Preproc.confounds.filter{end+1}=filts{n1}; 
            else CONN_x.Preproc.confounds.filter{end+1}=CONN_x.Preproc.variables.filter{idx};
            end
% 			if length(dims)>=n1&&~isempty(dims{n1}), CONN_x.Preproc.confounds.dimensions{end+1}=[min(dims{n1}(1),CONN_x.Preproc.variables.dimensions{idx}(1)),CONN_x.Preproc.variables.dimensions{idx}(1)]; 
%             else CONN_x.Preproc.confounds.dimensions{end+1}=CONN_x.Preproc.variables.dimensions{idx}; 
%             end
			if length(dims)>=n1&&~isempty(dims{n1}), CONN_x.Preproc.confounds.dimensions{end+1}=[dims{n1}(1),CONN_x.Preproc.variables.dimensions{idx}(1)]; 
            else CONN_x.Preproc.confounds.dimensions{end+1}=[inf,CONN_x.Preproc.variables.dimensions{idx}(1)]; 
            end
		end
	end
    conn_waitbar(1,h);
	if ~isfield(CONN_x.Preproc,'filter')||isempty(CONN_x.Preproc.filter), 
        maxrt=max(conn_get_rt);
        CONN_x.Preproc.filter=[0,1/(2*maxrt)]; 
    end
	conn_waitbar('close',h);
    CONN_x.isready(2)=1;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates DATA_Subject###_Condition###.mat files (whole-brain data after removal of confounding effects & filtering)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(options==6) && any(CONN_x.Setup.steps([2,3])) && ~(isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'steps')&&~any(CONN_x.gui.steps([2,3]))),
    warning('off','MATLAB:DELETE:FileNotFound');
	[path,name,ext]=fileparts(CONN_x.filename);
    filepath=CONN_x.folders.data;
	filepathresults=CONN_x.folders.preprocessing;
    try, conn_fileutils('mkdir',filepathresults); end
    if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'subjects'), validsubjects=CONN_x.gui.subjects; else validsubjects=1:CONN_x.Setup.nsubjects; end
    if isfield(CONN_x,'pobj')&&isstruct(CONN_x.pobj)&&isfield(CONN_x.pobj,'subjects'), validsubjects=CONN_x.pobj.subjects; end
	nconditions=length(CONN_x.Setup.conditions.names)-1;
    if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'conditions')&&~isempty(CONN_x.gui.conditions), validconditions=CONN_x.gui.conditions; else validconditions=1:length(CONN_x.Setup.conditions.names)-1; end
    icondition=[];isnewcondition=[];for ncondition=validconditions,[icondition(ncondition),isnewcondition(ncondition)]=conn_conditionnames(CONN_x.Setup.conditions.names{ncondition},'+'); end
    validconditions=validconditions(cellfun('length',CONN_x.Setup.conditions.model(validconditions))==0); 
	h=conn_waitbar(0,['Step ',num2str(sum(options<=6)),'/',num2str(length(options)),': Denoising functional data']);
    REDO=[];%'Yes';filename=fullfile(filepathresults,['DATA_Subject',num2str(1,'%03d'),'_Condition',num2str(1,'%03d'),'.mat']);
    %if ~isempty(dir(filename)),if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'overwrite'), REDO=CONN_x.gui.overwrite; else, REDO=conn_questdlg('Overwrite existing subject results?','','Yes', 'No', 'Yes');end; end
    DONEWOUTPUTCONFCORR=1; % create confound-corrected timeseries option (0: old format -one file per condition-; 1: new format -one file per session-; 2: both)
    redone_files=0;
    reportedsettings=false;
	N=numel(validsubjects); %sum(CONN_x.Setup.nsessions); if length(CONN_x.Setup.nsessions)==1, N=N*CONN_x.Setup.nsubjects; end
	n=0;
	for nsub=validsubjects,
        missingdata=arrayfun(@(n)isempty(dir(fullfile(filepathresults,['DATA_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(n),'%03d'),'.mat']))),validconditions);
        if isempty(REDO)&&~all(missingdata),
            if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'overwrite'), REDO=CONN_x.gui.overwrite; else, REDO=conn_questdlg('Overwrite existing subject results?','','Yes', 'No', 'Yes');end;
        end
        if strcmp(lower(REDO),'yes')||any(missingdata),
            if ~reportedsettings
                reportedsettings=true;
                try,
                    t12={'non-filtered','band-pass filtered'};
                    tnames=CONN_x.Preproc.confounds.names;
                    tnames=cellfun(@(name,dim,deriv,power)sprintf('%s (%dP)',name,min(dim)*(power*(1+deriv))),CONN_x.Preproc.confounds.names,CONN_x.Preproc.confounds.dimensions,CONN_x.Preproc.confounds.deriv,CONN_x.Preproc.confounds.power,'uni',0);
                    for nreg=1:numel(tnames),
                        switch(CONN_x.Preproc.confounds.power{nreg}), case 1, t13='no polynomial expansion'; case 2, t13='add quadratic effects'; case 3, t13='add cubic effects'; otherwise, t13=''; end
                        switch(CONN_x.Preproc.confounds.deriv{nreg}), case 0, t14='no temporal expansion'; case 1, t14='added first-order derivative'; case 2, t14='added second-order derivative'; otherwise, t14=''; end
                        switch(CONN_x.Preproc.confounds.filter{nreg}), case 0, t15='non-filtered'; otherwise, t15='band-pass filtered'; end
                        conn_disp('fprintf','  Regression : %s (%d dimensions; %s; %s; %s)\n',tnames{nreg},min(CONN_x.Preproc.confounds.dimensions{nreg}),t13,t14,t15);
                    end
                    switch(CONN_x.Preproc.regbp), case 2, t16='Simultaneous regression and filtering (Simult)'; otherwise, t16='Regression followed by filtering (RegBP)'; end
                    conn_disp('fprintf','  Order : %s\n',t16);
                    conn_disp('fprintf','  Detrending : %s\n',mat2str(CONN_x.Preproc.detrending));
                    conn_disp('fprintf','  Despiking : %s\n',mat2str(CONN_x.Preproc.despiking));
                    conn_disp('fprintf','  Bandpass : %s\n',mat2str(CONN_x.Preproc.filter));
                end
            end
            if strcmp(lower(REDO),'yes'), missingdata(:)=true; end
            nsess=CONN_x.Setup.nsessions(min(length(CONN_x.Setup.nsessions),nsub));
            clear Y X iX X1 X2 C Xnames;
            clear Youtnorm0 cachenorm0 Voutputfiles;
            for nses=1:nsess, % loads all ROI COV COND data for this subject 
                filename=fullfile(filepath,['DATA_Subject',num2str(nsub,'%03d'),'_Session',num2str(nses,'%03d'),'.mat']);
                if ~conn_existfile(filename), conn_disp(['Not ready to process step conn_process_6']); conn_waitbar('close',h); return; end
                Y{nses}=conn_vol(filename);
                if Y{nses}.size.Nt==0||sum(Y{nses}.size.Nv)==0, error(sprintf('Empty data (samples=%d,voxels=%d) in file %s',Y{nses}.size.Nt,sum(Y{nses}.size.Nv),filename)); end
                filename=fullfile(filepath,['ROI_Subject',num2str(nsub,'%03d'),'_Session',num2str(nses,'%03d'),'.mat']);
                X1{nses}=load(filename);
                filename=fullfile(filepath,['COV_Subject',num2str(nsub,'%03d'),'_Session',num2str(nses,'%03d'),'.mat']);
                X2{nses}=load(filename);
                filename=fullfile(filepath,['COND_Subject',num2str(nsub,'%03d'),'_Session',num2str(nses,'%03d'),'.mat']);
                C{nses}=load(filename);
                if ~isequal(CONN_x.Setup.conditions.names(1:end-1),C{nses}.names), error(['Incorrect conditions in file ',filename,'. Re-run previous step']); end
                confounds=CONN_x.Preproc.confounds;
                nfilter=find(cellfun(@(x)max(x),CONN_x.Preproc.confounds.filter));
                if isfield(CONN_x.Preproc,'detrending')&&CONN_x.Preproc.detrending, 
                    confounds.types{end+1}='detrend'; 
                    if CONN_x.Preproc.detrending>=2, confounds.types{end+1}='detrend2'; end
                    if CONN_x.Preproc.detrending>=3, confounds.types{end+1}='detrend3'; end
                end
                [X{nses},ifilter,nill,nill,Xnames{nses}]=conn_designmatrix(confounds,X1{nses},X2{nses},{nfilter});
                Xnames{nses}=regexprep(Xnames{nses},'^.*$',sprintf('Session %d: $0',nses));
                Xconstant{nses}=cellfun('length',regexp(Xnames{nses},'constant term$'))>0;
                RT=conn_get_rt(nsub,nses);
                if isfield(CONN_x.Preproc,'regbp')&&CONN_x.Preproc.regbp==2,
                    X{nses}(:,~Xconstant{nses})=conn_filter(RT,CONN_x.Preproc.filter,X{nses}(:,~Xconstant{nses}));
                elseif nnz(ifilter{1})
                    X{nses}(:,find(ifilter{1}))=conn_filter(max(RT),CONN_x.Preproc.filter,X{nses}(:,find(ifilter{1})));
                end
                if size(X{nses},1)~=CONN_x.Setup.nscans{nsub}{nses}, error('Wrong dimensions'); end
                try, iX{nses}=pinv(X{nses});
                catch, iX{nses}=pinv(X{nses}'*X{nses})*X{nses}';
                end
                filename=fullfile(filepathresults,['NORMS_Subject',num2str(nsub,'%03d'),'_Session',num2str(nses,'%03d'),'.mat']);
                [filename,cachenorm0(nses)]=conn_tempcache(filename); 
                Youtnorm0{nses}=Y{1};
                Youtnorm0{nses}.size.Nt=1;
                Youtnorm0{nses}.fname=filename;
                Youtnorm0{nses}=conn_init_vol(Youtnorm0{nses});
                if isfield(CONN_x.Setup,'outputfiles')&&numel(CONN_x.Setup.outputfiles)>=2&&CONN_x.Setup.outputfiles(2)&&DONEWOUTPUTCONFCORR>=1,
                    try,
                        fileout=conn_prepend('d',cellstr(conn_get_functional(nsub,nses)));
                        if numel(fileout)>1, fileout=fileout(1); end
                        Voutputfiles{nses}=struct('fname',regexprep(char(fileout),',\d+$',''),...
                            'mat',Y{nses}.matdim.mat,...
                            'dim',Y{nses}.matdim.dim(1:3),...
                            'n',[1,1],...
                            'pinfo',[1;0;0],...
                            'dt',[spm_type('float32'),spm_platform('bigend')],...
                            'descrip','denoised');
                        Voutputfiles{nses}=repmat(Voutputfiles{nses},[Y{nses}.size.Nt,1]);for nt=1:Y{nses}.size.Nt,Voutputfiles{nses}(nt).n=[nt,1];end
                        spm_unlink(fileout{:});
                        Voutputfiles{nses}=spm_create_vol(Voutputfiles{nses}); 
                    catch
                        Voutputfiles{nses}=[];
                        try, str=lasterror; conn_disp('Warning: unable to create denoised files; error message:'); conn_disp(str.message); end
                    end
                end
                if isfield(CONN_x.Preproc,'regbp')&&CONN_x.Preproc.regbp==2, dof2=max(0,CONN_x.Setup.nscans{nsub}{nses}*(min(1/(2*max(RT)),CONN_x.Preproc.filter(2))-max(0,CONN_x.Preproc.filter(1)))/(1/(2*max(RT)))+0-size(X{nses},2));
                elseif nnz(ifilter{1}), dof2=max(0,(CONN_x.Setup.nscans{nsub}{nses}-size(X{nses},2)+nnz(ifilter{1}))*(min(1/(2*max(RT)),CONN_x.Preproc.filter(2))-max(0,CONN_x.Preproc.filter(1)))/(1/(2*max(RT)))+0-nnz(ifilter{1}));
                else dof2=max(0,(CONN_x.Setup.nscans{nsub}{nses}-size(X{nses},2))*(min(1/(2*max(RT)),CONN_x.Preproc.filter(2))-max(0,CONN_x.Preproc.filter(1)))/(1/(2*max(RT)))+0);
                end
                edof_name=sprintf('QC_DOF_session%d',nses);
                edof_icov=find(strcmp(edof_name,CONN_x.Setup.l2covariates.names(1:end-1)),1);
                if isempty(edof_icov),
                    edof_icov=numel(CONN_x.Setup.l2covariates.names);
                    CONN_x.Setup.l2covariates.names{edof_icov}=edof_name;
                    CONN_x.Setup.l2covariates.descrip{edof_icov}='CONN Quality Assurance: Effective degrees of freedom (after denoising)';
                    CONN_x.Setup.l2covariates.names{edof_icov+1}=' ';
                    for tnsub=1:CONN_x.Setup.nsubjects, CONN_x.Setup.l2covariates.values{tnsub}{edof_icov}=nan; end
                end
                CONN_x.Setup.l2covariates.values{nsub}{edof_icov}=dof2;
            end
            clear nsamples time0 dataroi conditionsweights;
            crop=0;
            for ncondition=validconditions(missingdata), % computes number of samples per condition
                nsamples{ncondition}=0; for nses=1:nsess, nsamples{ncondition}=nsamples{ncondition}+length(C{nses}.samples{ncondition}); end
                % 			dataroi{ncondition}=cell(1,length(X1{1}.data));
                conditionsweights{ncondition}=cell(1,length(C{1}.weights{ncondition}));
                for nses=1:nsess,
                    if ~isfield(C{nses},'crop')||C{nses}.crop>0, crop=1; end
                    for nweight=1:length(C{nses}.weights{ncondition}),
                        conditionsweights{ncondition}{nweight}=cat(1,conditionsweights{ncondition}{nweight},C{nses}.weights{ncondition}{nweight});
                    end
                end
            end
            clear Yout cache Youtnorm cachenorm;
            softdone=false;
            softlinkoverwrite=[];
            softlinkcache=[];
            for ncondition=validconditions(missingdata),
                filename=fullfile(filepathresults,['DATA_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'.mat']);
                [filename,cache(ncondition)]=conn_tempcache(filename); 
                Yout{ncondition}=Y{1}; Yout{ncondition}.fname=filename;
                Yout{ncondition}.size.Nt=nsamples{ncondition};
                if ~Yout{ncondition}.size.Nt, error(['No samples in file ',filename]); end
                Yout{ncondition}.conditionsweights=conditionsweights{ncondition};
                Yout{ncondition}.crop=crop;
                wx=conditionsweights{ncondition}{1};
                emptycondition=~nnz(~isnan(wx)&wx~=0);
                filename=fullfile(filepathresults,['NORMS_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'.mat']);
                [filename,cachenorm(ncondition)]=conn_tempcache(filename); 
                Youtnorm{ncondition}=Yout{ncondition};
                Youtnorm{ncondition}.size.Nt=2;
                Youtnorm{ncondition}.fname=filename;
                Youtnorm{ncondition}.EmptyData=emptycondition;
                if crop||numel(CONN_x.Setup.conditions.filter{ncondition})==2, softlink=[]; softlinkoverwrite(ncondition)=true;
                else
                    softlink=fullfile(filepathresults,['DATA_Subject',num2str(nsub,'%03d'),'_Condition',num2str(0,'%03d'),'.mat']);
                    softlinkoverwrite(ncondition)=strcmp(lower(REDO),'yes')||~conn_existfile([softlink 'c']); % avoids data duplication
                    if softlinkoverwrite(ncondition), [softlink,softlinkcache]=conn_tempcache(softlink); end
                    softdone=true;
                end
                Yout{ncondition}=conn_init_vol(Yout{ncondition},[],[],softlink,softlinkoverwrite(ncondition));
                Youtnorm{ncondition}=conn_init_vol(Youtnorm{ncondition});
                redone_files=redone_files+1;
            end
            if any(softlinkoverwrite)
                B=[]; RT=[]; 
                for slice=1:Y{1}.matdim.dim(3),
                    Bb=[];
                    if 1, % analyses per slice (all sessions together, faster but requires more memory)
                        clear y ypre;
                        for nses=1:nsess,
                            [y{nses},idx]=conn_get_slice(Y{nses},slice);
                            if size(y{nses},1)~=CONN_x.Setup.nscans{nsub}{nses}, error('Wrong dimensions'); end
                            if slice==1&&rank(X{nses})>=size(y{nses},1), conn_disp(['Warning: Over-determined model (no degrees of freedom for subject ',num2str(nsub),' session ',num2str(nses),'). Please consider reducing the number, dimensions, or covariates order of the confounds or disregarding this subject/session']); end
                            if isfield(CONN_x.Preproc,'despiking')&&CONN_x.Preproc.despiking==1,
                                my=repmat(median(y{nses},1),[size(y{nses},1),1]);
                                sy=repmat(4*median(abs(y{nses}-my)),[size(y{nses},1),1]);
                                y{nses}=my+sy.*tanh((y{nses}-my)./max(eps,sy));
                            end
                            b=iX{nses}*y{nses};
                            y{nses}=y{nses}-X{nses}*b;
                            if isfield(CONN_x.Preproc,'despiking')&&CONN_x.Preproc.despiking==2,
                                my=repmat(median(y{nses},1),[size(y{nses},1),1]);
                                sy=repmat(4*median(abs(y{nses}-my)),[size(y{nses},1),1]);
                                y{nses}=my+sy.*tanh((y{nses}-my)./max(eps,sy));
                            end
                            ypre{nses}=y{nses};
                            if numel(RT)<nses, RT(nses)=conn_get_rt(nsub,nses); end
                            y{nses}=conn_filter(RT(nses),CONN_x.Preproc.filter,y{nses});
                            if isfield(CONN_x.Setup,'outputfiles')&&numel(CONN_x.Setup.outputfiles)>=1&&CONN_x.Setup.outputfiles(1),
                                Bb=cat(1,Bb,b);
                            end
                            norm_pre=sqrt(max(0,mean(ypre{nses}.^2,1)));
                            conn_write_slice(Youtnorm0{nses},norm_pre,slice);
                            if DONEWOUTPUTCONFCORR>=1&&isfield(CONN_x.Setup,'outputfiles')&&numel(CONN_x.Setup.outputfiles)>=2&&CONN_x.Setup.outputfiles(2)&&~isempty(Voutputfiles{nses}),
                                t=zeros(Y{1}.matdim.dim(1:2));
                                for nt=1:size(y{nses},1),
                                    t(idx)=y{nses}(nt,:) + X{nses}(nt,Xconstant{nses})*b(Xconstant{nses},:);
                                    Voutputfiles{nses}(nt)=spm_write_plane(Voutputfiles{nses}(nt),t,slice);
                                end
                            end
                        end
                        for ncondition=validconditions(missingdata),
                            ytemp=[];
                            for nses=1:nsess,
                                ytemp0=ypre{nses};
                                ytemp=cat(1,ytemp,ytemp0(C{nses}.samples{ncondition},:));
                            end
                            norm_pre=sqrt(max(0,sum(repmat(Yout{ncondition}.conditionsweights{1},1,size(ytemp,2)).*ytemp.^2,1)/max(eps,sum(Yout{ncondition}.conditionsweights{1}))));
                            ytemp=[];
                            for nses=1:nsess,
                                if numel(CONN_x.Setup.conditions.filter{ncondition})==2, ytemp0=conn_filter(RT(nses),CONN_x.Setup.conditions.filter{ncondition},y{nses});
                                else ytemp0=y{nses};
                                end
                                ytemp=cat(1,ytemp,ytemp0(C{nses}.samples{ncondition},:));
                            end
                            norm_post=sqrt(max(0,sum(repmat(Yout{ncondition}.conditionsweights{1},1,size(ytemp,2)).*ytemp.^2,1)/max(eps,sum(Yout{ncondition}.conditionsweights{1}))));
                            conn_write_slice(Yout{ncondition},ytemp,slice);
                            if ~Youtnorm{ncondition}.EmptyData, conn_write_slice(Youtnorm{ncondition},cat(1,norm_pre,norm_post),slice); end
                        end
                    else, % analyses per slice/session (slower but requires less memory; obsolete now)
                        for nses=1:nsess,
                            [y,idx]=conn_get_slice(Y{nses},slice);
                            if size(y,1)~=CONN_x.Setup.nscans{nsub}{nses}, error('Wrong dimensions'); end
                            if slice==1&&rank(X{nses})>=size(y,1), conn_disp(['Warning: Over-determined model (no degrees of freedom for subject ',num2str(nsub),' session ',num2str(nses),'). Please consider reducing the number, dimensions, or covariates order of the confounds or disregarding this subject/session']); end
                            if isfield(CONN_x.Preproc,'despiking')&&CONN_x.Preproc.despiking==1,
                                my=repmat(median(y,1),[size(y,1),1]);
                                sy=repmat(4*median(abs(y-my)),[size(y,1),1]);
                                y=my+sy.*tanh((y-my)./max(eps,sy));
                            end
                            b=iX{nses}*y;
                            y=y-X{nses}*b;
                            if isfield(CONN_x.Preproc,'despiking')&&CONN_x.Preproc.despiking==2,
                                my=repmat(median(y,1),[size(y,1),1]);
                                sy=repmat(4*median(abs(y-my)),[size(y,1),1]);
                                y=my+sy.*tanh((y-my)./max(eps,sy));
                            end
                            y=conn_filter(conn_get_rt(nsub,nses),CONN_x.Preproc.filter,y);
                            if isfield(CONN_x.Setup,'outputfiles')&&numel(CONN_x.Setup.outputfiles)>=1&&CONN_x.Setup.outputfiles(1),
                                Bb=cat(1,Bb,b);
                            end
                            for ncondition=validconditions(missingdata),
                                if nses==1, time0{ncondition}=1; end
                                if numel(CONN_x.Setup.conditions.filter{ncondition})==2, ytemp0=conn_filter(conn_get_rt(nsub,nses),CONN_x.Setup.conditions.filter{ncondition},y);
                                else ytemp0=y;
                                end
                                conn_write_slice(Yout{ncondition},ytemp0(C{nses}.samples{ncondition},:),slice,time0{ncondition});
                                time0{ncondition}=time0{ncondition}+length(C{nses}.samples{ncondition});
                            end
                        end
                    end
                    if isfield(CONN_x.Setup,'outputfiles')&&numel(CONN_x.Setup.outputfiles)>=1&&CONN_x.Setup.outputfiles(1),
                        if 0, %%% 0 to avoid memory errors
                            B=cat(2,B,Bb);
                        else
                            if slice==1
                                V=struct('mat',Y{1}.matdim.mat,'dim',Y{1}.matdim.dim,'pinfo',[1;0;0],'fname',fullfile(filepathresults,['BETA_denoising_Subject',num2str(nsub,'%03d'),'.nii']),...
                                    'dt',[spm_type('float32') spm_platform('bigend')]);
                                try, delete(V.fname); end
                                %V=CONN_x.Setup.structural{nsub}{3}; V=V(1);
                                %if isfield(V,'dt'), V.dt=[spm_type('float32') spm_platform('bigend')];
                                %elseif length(V.dim)>3, V.dim(4)=spm_type('float32'); end
                                %V.fname=fullfile(filepathresults,['BETA_Subject',num2str(nsub,'%03d'),'.nii']);
                                V=repmat(V,[size(Bb,1),1]);for nh=1:numel(V),V(nh).n=[nh,1];end
                                V=spm_create_vol(V);
                                try
                                    BETAnames=cat(2,Xnames{:});
                                    save(fullfile(filepathresults,'_list_BETA_denoising.mat'),'BETAnames','X','Xnames');
                                    fileout=fullfile(filepathresults,'_list_BETA_denoising.txt');
                                    fh=fopen(fileout,'wt');
                                    for n1=1:length(BETAnames),fprintf(fh,'%s\n',BETAnames{n1});end
                                    fclose(fh);
                                end
                                %                         for n1=1:size(Bb,1)
                                %                             V(n1)=V(1);
                                %                             V(n1).fname=fullfile(filepathresults,['BETA_Subject',num2str(nsub),'_Regressor',num2str(n1,'%04d'),'.nii']);
                                %                             V(n1)=spm_create_vol(V(n1));
                                %                         end
                            end
                            t=nan+zeros(Y{1}.matdim.dim(1:2));
                            for n1=1:size(Bb,1),
                                t(idx)=Bb(n1,:);
                                V(n1)=spm_write_plane(V(n1),t,slice);
                            end
                        end
                    end
                    conn_waitbar((n+slice/Y{1}.matdim.dim(3))/N,h,sprintf('Subject %d',nsub));
                end
                
                if 0,%%% 0 to avoid memory errors
                    if isfield(CONN_x.Setup,'outputfiles')&&numel(CONN_x.Setup.outputfiles)>=1&&CONN_x.Setup.outputfiles(1),
                        t=nan+zeros(Y{1}.matdim.dim);
                        V=struct('mat',Y{1}.matdim.mat,'dim',Y{1}.matdim.dim,'pinfo',[1;0;0],'fname','',...
                            'dt',[spm_type('float32') spm_platform('bigend')]);
                        %V=CONN_x.Setup.structural{nsub}{3}; V=V(1);
                        %if isfield(V,'dt'), V.dt=[spm_type('float32') spm_platform('bigend')];
                        %elseif length(V.dim)>3, V.dim(4)=spm_type('float32'); end
                        for n1=1:size(B,1),
                            t(Y{1}.voxels)=B(n1,:);
                            V.fname=fullfile(filepathresults,['BETA_Subject',num2str(nsub,'%03d'),'_Regressor',num2str(n1,'%04d'),'.nii']);
                            spm_write_vol(V,t);
                        end
                    end
                end
            end
            if ~isempty(softlinkcache)
                conn_tempcache(softlinkcache,'matc');
            end
            if any(softlinkoverwrite)
                bstd=0;
                for nses=1:nsess,
                    conn_tempcache(cachenorm0(nses),'matc');
                    filename=fullfile(filepathresults,['NORMS_Subject',num2str(nsub,'%03d'),'_Session',num2str(nses,'%03d'),'.mat']);
                    [nill,outvals]=conn_matc2nii(filename,0);
                    bstd=bstd+outvals{1};
                end
                bstd=bstd/nsess;
                bstd_name='QC_BOLDstd';
                bstd_icov=find(strcmp(bstd_name,CONN_x.Setup.l2covariates.names(1:end-1)),1);
                if isempty(bstd_icov),
                    bstd_icov=numel(CONN_x.Setup.l2covariates.names);
                    CONN_x.Setup.l2covariates.names{bstd_icov}=bstd_name;
                    CONN_x.Setup.l2covariates.descrip{bstd_icov}='CONN Quality Assurance: BOLD signal standard deviation (after denoising)';
                    CONN_x.Setup.l2covariates.names{bstd_icov+1}=' ';
                    for tnsub=1:CONN_x.Setup.nsubjects, CONN_x.Setup.l2covariates.values{tnsub}{bstd_icov}=nan; end
                end
                CONN_x.Setup.l2covariates.values{nsub}{bstd_icov}=bstd;
            end
            for ncondition=validconditions(missingdata),
                conn_tempcache(cachenorm(ncondition),'matc');
                conn_tempcache(cache(ncondition),'matc');
            end
            if DONEWOUTPUTCONFCORR~=1&&isfield(CONN_x.Setup,'outputfiles')&&numel(CONN_x.Setup.outputfiles)>=2&&CONN_x.Setup.outputfiles(2),
                filename={}; for ncondition=validconditions(missingdata),filename{ncondition}=fullfile(filepathresults,['DATA_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'.mat']); end
                conn_matc2nii(filename(validconditions(missingdata)),0);
            end
        end
        
        n=n+1;
        conn_waitbar(n/N,h,sprintf('Subject %d',nsub));
    end
    conn_disp('fprintf','      processed %d files\n',redone_files);
    conn_waitbar('close',h);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates ROI_Subject###_Condition###.mat files (roi data after removal of confounding effects & filtering)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(options==7) && any(CONN_x.Setup.steps([1,2,4])) && ~(isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'steps')&&~any(CONN_x.gui.steps([1,2]))),
	[path,name,ext]=fileparts(CONN_x.filename);
    filepath=CONN_x.folders.data;
	filepathresults=CONN_x.folders.preprocessing;
    try, conn_fileutils('mkdir',filepathresults); end
    if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'subjects'), validsubjects=CONN_x.gui.subjects; else validsubjects=1:CONN_x.Setup.nsubjects; end
    if isfield(CONN_x,'pobj')&&isstruct(CONN_x.pobj)&&isfield(CONN_x.pobj,'subjects'), validsubjects=CONN_x.pobj.subjects; end
	nconditions=length(CONN_x.Setup.conditions.names)-1;
    runallconditions=true;
    if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'conditions')&&~isempty(CONN_x.gui.conditions), validconditions=CONN_x.gui.conditions; runallconditions=isempty(setdiff(1:length(CONN_x.Setup.conditions.names)-1,validconditions)); else validconditions=1:length(CONN_x.Setup.conditions.names)-1; end
    icondition=[];isnewcondition=[];for ncondition=validconditions,[icondition(ncondition),isnewcondition(ncondition)]=conn_conditionnames(CONN_x.Setup.conditions.names{ncondition},'+'); end
    validconditions=validconditions(cellfun('length',CONN_x.Setup.conditions.model(validconditions))==0);     
	h=conn_waitbar(0,['Step ',num2str(sum(options<=7)),'/',num2str(length(options)),': Denoising ROI data']);
    REDO=[];%filename=fullfile(filepathresults,['ROI_Subject',num2str(1,'%03d'),'_Condition',num2str(icondition(validconditions(1)),'%03d'),'.mat']);
    %if ~isempty(dir(filename)),if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'overwrite'), REDO=CONN_x.gui.overwrite; else, REDO=conn_questdlg('Overwrite existing subject results?','','Yes', 'No', 'Yes');end; end
    PREPROCESSCOVARIATES=false; % set to true if 1st-level covariate should be denoised before moving to first-level analyses
    NUMBEROFFREQBANDS=8;
    reportedsettings=false;
    redone_files=0;
    maxrt=nan;
	N=numel(validsubjects); %sum(CONN_x.Setup.nsessions); if length(CONN_x.Setup.nsessions)==1, N=N*CONN_x.Setup.nsubjects; end
	n=0;
	for nsub=validsubjects,
        missingdata=arrayfun(@(n)isempty(dir(fullfile(filepathresults,['ROI_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(n),'%03d'),'.mat']))),validconditions);
        if isempty(REDO)&&~any(missingdata),
            if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'overwrite'), REDO=CONN_x.gui.overwrite; else, REDO=conn_questdlg('Overwrite existing subject results?','','Yes', 'No', 'Yes');end;
        end
        
        nsess=CONN_x.Setup.nsessions(min(length(CONN_x.Setup.nsessions),nsub));
        clear Y X iX X1 X2 C;
        for nses=1:nsess, % loads all ROI COV COND data for this subject
            filename=fullfile(filepath,['ROI_Subject',num2str(nsub,'%03d'),'_Session',num2str(nses,'%03d'),'.mat']);
            if ~conn_existfile(filename), conn_disp(['Not ready to process step conn_process_7']); conn_waitbar('close',h); return; end
            X1{nses}=load(filename);
            filename=fullfile(filepath,['COV_Subject',num2str(nsub,'%03d'),'_Session',num2str(nses,'%03d'),'.mat']);
            X2{nses}=load(filename);
            filename=fullfile(filepath,['COND_Subject',num2str(nsub,'%03d'),'_Session',num2str(nses,'%03d'),'.mat']);
            C{nses}=load(filename);
            if ~isequal(CONN_x.Setup.conditions.names(1:end-1),C{nses}.names), error(['Incorrect conditions in file ',filename,'. Re-run previous step']); end
            confounds=CONN_x.Preproc.confounds;
            nfilter=find(cellfun(@(x)max(x),CONN_x.Preproc.confounds.filter));
            if isfield(CONN_x.Preproc,'detrending')&&CONN_x.Preproc.detrending,
                confounds.types{end+1}='detrend';
                if CONN_x.Preproc.detrending>=2, confounds.types{end+1}='detrend2'; end
                if CONN_x.Preproc.detrending>=3, confounds.types{end+1}='detrend3'; end
            end
            [X{nses},ifilter]=conn_designmatrix(confounds,X1{nses},X2{nses},{nfilter});
            if isfield(CONN_x.Preproc,'regbp')&&CONN_x.Preproc.regbp==2,
                X{nses}=conn_filter(conn_get_rt(nsub,nses),CONN_x.Preproc.filter,X{nses});
            elseif nnz(ifilter{1})
                X{nses}(:,find(ifilter{1}))=conn_filter(max(conn_get_rt(nsub,nses)),CONN_x.Preproc.filter,X{nses}(:,find(ifilter{1})));
            end
            if size(X{nses},1)~=CONN_x.Setup.nscans{nsub}{nses}, error('Wrong dimensions'); end
            try, iX{nses}=pinv(X{nses});
            catch, iX{nses}=pinv(X{nses}'*X{nses})*X{nses}';
            end
        end
        
        if strcmp(lower(REDO),'yes')||any(missingdata), redo=true;
        else
            names=cat(2,X1{1}.names,X2{1}.names);
            filename=fullfile(filepathresults,['ROI_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(validconditions(1)),'%03d'),'.mat']);
            namestomatch=load(filename,'names');
            redo=~isequal(namestomatch.names,names);
        end
        if redo
            if ~reportedsettings
                reportedsettings=true;
                try,
                    t12={'non-filtered','band-pass filtered'};
                    tnames=CONN_x.Preproc.confounds.names;
                    tnames=cellfun(@(name,dim,deriv,power)sprintf('%s (%dP)',name,min(dim)*(power*(1+deriv))),CONN_x.Preproc.confounds.names,CONN_x.Preproc.confounds.dimensions,CONN_x.Preproc.confounds.deriv,CONN_x.Preproc.confounds.power,'uni',0);
                    for nreg=1:numel(tnames),
                        switch(CONN_x.Preproc.confounds.power{nreg}), case 1, t13='no polynomial expansion'; case 2, t13='add quadratic effects'; case 3, t13='add cubic effects'; otherwise, t13=''; end
                        switch(CONN_x.Preproc.confounds.deriv{nreg}), case 0, t14='no temporal expansion'; case 1, t14='added first-order derivative'; case 2, t14='added second-order derivative'; otherwise, t14=''; end
                        switch(CONN_x.Preproc.confounds.filter{nreg}), case 0, t15='non-filtered'; otherwise, t15='band-pass filtered'; end
                        conn_disp('fprintf','  Regression : %s (%d dimensions; %s; %s; %s)\n',tnames{nreg},min(CONN_x.Preproc.confounds.dimensions{nreg}),t13,t14,t15);
                    end
                    switch(CONN_x.Preproc.regbp), case 2, t16='Simultaneous regression and filtering (Simult)'; otherwise, t16='Regression followed by filtering (RegBP)'; end
                    conn_disp('fprintf','  Order : %s\n',t16);
                    conn_disp('fprintf','  Detrending : %s\n',mat2str(CONN_x.Preproc.detrending));
                    conn_disp('fprintf','  Despiking : %s\n',mat2str(CONN_x.Preproc.despiking));
                    conn_disp('fprintf','  Bandpass : %s\n',mat2str(CONN_x.Preproc.filter));
                end
            end
            clear nsamples time0 dataroi d1dataroi d2dataroi dbdataroi conditionsweights;
            crop=0;
            dataroiall=cell(1,length(X1{1}.data)); % no covariates there +length(X2{1}.data));
            dataroiall_sessions=[];
            for ncondition=validconditions, % computes number of samples per condition
                nsamples{ncondition}=0; for nses=1:nsess, nsamples{ncondition}=nsamples{ncondition}+length(C{nses}.samples{ncondition}); end
                dataroi{ncondition}=cell(1,length(X1{1}.data)+length(X2{1}.data));
                d1dataroi{ncondition}=cell(1,length(X1{1}.data)+length(X2{1}.data));
                d2dataroi{ncondition}=cell(1,length(X1{1}.data)+length(X2{1}.data));
                fbdataroi{ncondition}=repmat({cell(1,NUMBEROFFREQBANDS)},1,length(X1{1}.data)+length(X2{1}.data));
                conditionsweights{ncondition}=cell(1,length(C{1}.weights{ncondition}));
                for nses=1:nsess,
                    if ~isfield(C{nses},'crop')||C{nses}.crop>0, crop=1; end
                    for nweight=1:length(C{nses}.weights{ncondition}),
                        conditionsweights{ncondition}{nweight}=cat(1,conditionsweights{ncondition}{nweight},C{nses}.weights{ncondition}{nweight});
                    end
                end
            end
            
            RT=[];
            for nroi=1:length(X1{1}.data),
                for nses=1:nsess,
                    y=X1{nses}.data{nroi};
                    if size(y,1)~=CONN_x.Setup.nscans{nsub}{nses}, error('Wrong dimensions'); end
                    if isfield(CONN_x.Preproc,'despiking')&&CONN_x.Preproc.despiking==1,
                        my=repmat(median(y,1),[size(y,1),1]);
                        sy=repmat(4*median(abs(y-my)),[size(y,1),1]);
                        y=my+sy.*tanh((y-my)./max(eps,sy));
                    end
                    b=iX{nses}*y;
                    y=y-X{nses}*b;
                    if isfield(CONN_x.Preproc,'despiking')&&CONN_x.Preproc.despiking==2,
                        my=repmat(median(y,1),[size(y,1),1]);
                        sy=repmat(4*median(abs(y-my)),[size(y,1),1]);
                        y=my+sy.*tanh((y-my)./max(eps,sy));
                    end
                    if numel(RT)<nses, RT(nses)=conn_get_rt(nsub,nses); end
                    y=conn_filter(RT(nses),CONN_x.Preproc.filter,y);
                    ffilter=CONN_x.Preproc.filter;
                    if any(isinf(ffilter))&&isnan(maxrt), maxrt=max(conn_get_rt); end
                    ffilter(isinf(ffilter))=1/maxrt/2;
                    fby={};
                    for nfb=1:NUMBEROFFREQBANDS
                        for nfb2=1:nfb
                            fby{nfb}(:,:,nfb2)=conn_filter(RT(nses),ffilter(1)+(ffilter(2)-ffilter(1))/nfb*[nfb2-1,nfb2],y);
                        end
                    end
                    y0=y;
                    if ~isempty(dataroiall{nroi})&&~isempty(y)&&size(dataroiall{nroi},2)~=size(y,2), error('WARNING: incompatible dimensions in ROI %s subject %d session %d (observed %d expected %d)',X1{1}.names{nroi}, nsub, nses, size(y,2), size(dataroiall{nroi},2)); end
                    dataroiall{nroi}=cat(1,dataroiall{nroi},y);
                    if nroi==1, dataroiall_sessions=cat(1,dataroiall_sessions,nses+zeros(size(y,1),1)); end
                    for ncondition=validconditions,
                        y=y0;
                        if ~isempty(y)&&numel(CONN_x.Setup.conditions.filter{ncondition})==2, y=conn_filter(RT(nses),CONN_x.Setup.conditions.filter{ncondition},y); end
                        if ~isempty(y), d1y=convn(cat(1,y(1,:),y,y(end,:)),[1;0;-1]/2,'valid');
                        else d1y=y; end
                        if ~isempty(d1y), d2y=convn(cat(1,d1y(1,:),d1y,d1y(end,:)),[1;0;-1]/2,'valid');
                        else d2y=d1y; end
                        dataroi{ncondition}{nroi}=cat(1,dataroi{ncondition}{nroi},y(C{nses}.samples{ncondition},:));
                        d1dataroi{ncondition}{nroi}=cat(1,d1dataroi{ncondition}{nroi},d1y(C{nses}.samples{ncondition},:));
                        d2dataroi{ncondition}{nroi}=cat(1,d2dataroi{ncondition}{nroi},d2y(C{nses}.samples{ncondition},:));
                        for nfb=1:NUMBEROFFREQBANDS
                            fbdataroi{ncondition}{nroi}{nfb}=cat(1,fbdataroi{ncondition}{nroi}{nfb},fby{nfb}(C{nses}.samples{ncondition},:,:));
                        end
                    end
                end
            end
            
            for ncov=1:length(X2{1}.data),
                for nses=1:nsess,
                    y=X2{nses}.data{ncov};
                    if size(y,1)~=CONN_x.Setup.nscans{nsub}{nses}, error('Wrong dimensions'); end
                    if numel(RT)<nses, RT(nses)=conn_get_rt(nsub,nses); end
                    if PREPROCESSCOVARIATES
                        b=iX{nses}*y;
                        y=y-X{nses}*b;
                        y=conn_filter(RT(nses),CONN_x.Preproc.filter,y);
                    end
                    ffilter=CONN_x.Preproc.filter;
                    if any(isinf(ffilter))&&isnan(maxrt), maxrt=max(conn_get_rt); end
                    ffilter(isinf(ffilter))=1/maxrt/2;
                    fby={};
                    for nfb=1:NUMBEROFFREQBANDS
                        for nfb2=1:nfb
                            fby{nfb}(:,:,nfb2)=conn_filter(RT(nses),ffilter(1)+(ffilter(2)-ffilter(1))/nfb*[nfb2-1,nfb2],y);
                        end
                    end
                    y0=y;
                    fby0=fby;
                    for ncondition=validconditions,
                        y=y0;
                        fby=fby0;
                        if ~isempty(y)&&numel(CONN_x.Setup.conditions.filter{ncondition})==2, y=conn_filter(RT(nses),CONN_x.Setup.conditions.filter{ncondition},y); end
                        if ~isempty(y), d1y=convn(cat(1,y(1,:),y,y(end,:)),[1;0;-1]/2,'valid');
                        else d1y=y; end
                        if ~isempty(d1y), d2y=convn(cat(1,d1y(1,:),d1y,d1y(end,:)),[1;0;-1]/2,'valid');
                        else d2y=d1y; end
                        if size(y,2)>size(dataroi{ncondition}{length(X1{1}.data)+ncov},2),%&&size(dataroi{ncondition}{length(X1{1}.data)+ncov},2)>0, 
                            dataroi{ncondition}{length(X1{1}.data)+ncov}=cat(2,dataroi{ncondition}{length(X1{1}.data)+ncov}, zeros(size(dataroi{ncondition}{length(X1{1}.data)+ncov},1),size(y,2)-size(dataroi{ncondition}{length(X1{1}.data)+ncov},2))); 
                            d1dataroi{ncondition}{length(X1{1}.data)+ncov}=cat(2,d1dataroi{ncondition}{length(X1{1}.data)+ncov}, zeros(size(d1dataroi{ncondition}{length(X1{1}.data)+ncov},1),size(y,2)-size(d1dataroi{ncondition}{length(X1{1}.data)+ncov},2))); 
                            d2dataroi{ncondition}{length(X1{1}.data)+ncov}=cat(2,d2dataroi{ncondition}{length(X1{1}.data)+ncov}, zeros(size(d2dataroi{ncondition}{length(X1{1}.data)+ncov},1),size(y,2)-size(d2dataroi{ncondition}{length(X1{1}.data)+ncov},2))); 
                            for nfb=1:NUMBEROFFREQBANDS
                                fbdataroi{ncondition}{length(X1{1}.data)+ncov}{nfb}=cat(2,fbdataroi{ncondition}{length(X1{1}.data)+ncov}{nfb}, zeros([size(fbdataroi{ncondition}{length(X1{1}.data)+ncov}{nfb},1),size(y,2)-size(fbdataroi{ncondition}{length(X1{1}.data)+ncov}{nfb},2),nfb])); 
                            end
                        end
                        if size(y,2)<size(dataroi{ncondition}{length(X1{1}.data)+ncov},2), 
                            y=cat(2,y,zeros(size(y,1),size(dataroi{ncondition}{length(X1{1}.data)+ncov},2)-size(y,2))); 
                            d1y=cat(2,d1y,zeros(size(d1y,1),size(d1dataroi{ncondition}{length(X1{1}.data)+ncov},2)-size(d1y,2))); 
                            d2y=cat(2,d2y,zeros(size(d2y,1),size(d2dataroi{ncondition}{length(X1{1}.data)+ncov},2)-size(d2y,2))); 
                            for nfb=1:NUMBEROFFREQBANDS
                                fby{nfb}=cat(2,fby{nfb},zeros([size(fby{nfb},1),size(fbdataroi{ncondition}{length(X1{1}.data)+ncov}{nfb},2)-size(fby{nfb},2),nfb]));
                            end
                        end
                        dataroi{ncondition}{length(X1{1}.data)+ncov}=cat(1,dataroi{ncondition}{length(X1{1}.data)+ncov},y(C{nses}.samples{ncondition},:));
                        d1dataroi{ncondition}{length(X1{1}.data)+ncov}=cat(1,d1dataroi{ncondition}{length(X1{1}.data)+ncov},d1y(C{nses}.samples{ncondition},:));
                        d2dataroi{ncondition}{length(X1{1}.data)+ncov}=cat(1,d2dataroi{ncondition}{length(X1{1}.data)+ncov},d2y(C{nses}.samples{ncondition},:));
                        for nfb=1:NUMBEROFFREQBANDS
                            fbdataroi{ncondition}{length(X1{1}.data)+ncov}{nfb}=cat(1,fbdataroi{ncondition}{length(X1{1}.data)+ncov}{nfb},fby{nfb}(C{nses}.samples{ncondition},:,:));
                        end
                    end
                end
            end  
            conn_waitbar((n+1/2)/N,h,sprintf('Subject %d',nsub));
            
            for ncondition=validconditions,
                data=dataroi{ncondition};
                d1data=d1dataroi{ncondition};
                d2data=d2dataroi{ncondition};
                fbdata=fbdataroi{ncondition};
                names=cat(2,X1{1}.names,X2{1}.names);
                xyz=cat(2,X1{1}.xyz,repmat({[nan,nan,nan]},[1,length(X2{1}.data)]));
                %voxels=X1{1}.voxels;
                %weights=X1{1}.weights;
                source=X1{1}.source;
                conditionname=C{1}.names{ncondition};
                conditionweights=conditionsweights{ncondition};
                filename=fullfile(filepathresults,['ROI_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'.mat']);
                save(filename,'data','d1data','d2data','fbdata','names','xyz','source','conditionweights','conditionname','crop');
                redone_files=redone_files+1;
%                 if str2num(version('-release'))>=14, save(filename,'-V6','data','names','xyz','source','conditionweights','conditionname');
%                 else, save(filename,'data','names','xyz','source','conditionweights','conditionname');end
            end
            filename=fullfile(filepathresults,['ROI_Subject',num2str(nsub,'%03d'),'_Condition',num2str(0,'%03d'),'.mat']);
            names=X1{1}.names;
            xyz=X1{1}.xyz;
            source=X1{1}.source;
            data=dataroiall;
            data_sessions=dataroiall_sessions;
            conditionsnames=C{1}.names;
            conditionsnamesvalid=C{1}.names(validconditions);
            if ~runallconditions&&conn_existfile(filename)
                try % fill with previous info
                    temp=load(filename,'conditionsweights','conditionsnames');
                    for ncondition=1:numel(temp.conditionsweights)
                        if ~ismember(ncondition,validconditions), conditionsweights{ncondition}=temp.conditionsweights{ncondition}; end % assert(isequal(conditionsnames{ncondition}),temp.conditionsnames{ncondition}); end
                    end
                end
            end
            save(filename,'data','names','xyz','source','crop','data_sessions',   'conditionsweights','conditionsnames','conditionsnamesvalid');
        end
        
        n=n+1;
        conn_waitbar(n/N,h,sprintf('Subject %d',nsub));
    end
    conn_disp('fprintf','      processed %d files\n',redone_files);
    conn_waitbar('close',h);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates vvPC_Subject###_Condition###.mat files (voxel-to-voxel first-level SVD decomposition)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(options==8) && any(CONN_x.Setup.steps([3])) && ~(isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'steps')&&~any(CONN_x.gui.steps([3]))),
    [path,name,ext]=fileparts(CONN_x.filename);
    filepath=CONN_x.folders.preprocessing;
    try, conn_fileutils('mkdir',filepath); end
    if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'subjects'), validsubjects=CONN_x.gui.subjects; else validsubjects=1:CONN_x.Setup.nsubjects; end
    if isfield(CONN_x,'pobj')&&isstruct(CONN_x.pobj)&&isfield(CONN_x.pobj,'subjects'), validsubjects=CONN_x.pobj.subjects; end
    nconditions=length(CONN_x.Setup.conditions.names)-1;
    if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'conditions')&&~isempty(CONN_x.gui.conditions), validconditions=CONN_x.gui.conditions; else validconditions=1:length(CONN_x.Setup.conditions.names)-1; end
    SKIPTIMEFREQ=true; % set to false if you wish to allow V2V analyses of sliding-window and frequency-bank conditions
    MAXDIMS=256; % maximum number of dimensions to keep in first-level dimensionality reduction
    if SKIPTIMEFREQ&&isequal(validconditions,1:length(CONN_x.Setup.conditions.names)-1), validconditions=validconditions(~cellfun('length',regexp(CONN_x.Setup.conditions.names(validconditions),' x FrequencyBand\d+$| x Time\d+$'))); end
    icondition=[];isnewcondition=[];for ncondition=validconditions,[icondition(ncondition),isnewcondition(ncondition)]=conn_conditionnames(CONN_x.Setup.conditions.names{ncondition},'+'); end
    validconditions=validconditions(cellfun('length',CONN_x.Setup.conditions.model(validconditions))==0); 
    h=conn_waitbar(0,['Step ',num2str(sum(options<=8)),'/',num2str(length(options)),': preprocessing voxel-to-voxel covariance']);
    REDO=[];
% 	filename=fullfile(fileparts(which(mfilename)),'utils','surf','mask.volume.brainmask.nii');
% 	Vmask=spm_vol(filename);
%     THR_MASK=.25;
    missingdata=arrayfun(@(n)isempty(dir(fullfile(filepath,['DATA_Subject',num2str(validsubjects(1),'%03d'),'_Condition',num2str(icondition(n),'%03d'),'.mat']))),validconditions);
    if any(missingdata), conn_disp(['Not ready to process step conn_process_8']); return; end
    filename=fullfile(filepath,['DATA_Subject',num2str(validsubjects(1),'%03d'),'_Condition',num2str(icondition(validconditions(1)),'%03d'),'.mat']);
    Y=conn_vol(filename);
    N=1.1*numel(validsubjects)*numel(validconditions)*Y.matdim.dim(3);n=0; 
        
    for nsub=validsubjects,
        maxrt=[];
        for ncondition=validconditions,
            %filename=fullfile(filepath,['vvPCcov_SubjectA',num2str(nsub,'%03d'),'_SubjectB',num2str(nsub,'%03d'),'_ConditionA',num2str(ncondition,'%03d'),'_ConditionB',num2str(ncondition,'%03d'),'.mat']); 
            filename=fullfile(filepath,['vvPC_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'.mat']);
            if isempty(REDO)&&conn_existfile(filename),
                if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'overwrite'), REDO=CONN_x.gui.overwrite; else, REDO=conn_questdlg('Overwrite existing subject results?','','Yes', 'No', 'Yes');end;
            end
            if strcmp(lower(REDO),'yes')||(~conn_existfile(filename)),
                filename=fullfile(filepath,['DATA_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'.mat']);
                if ~conn_existfile(filename), conn_disp(['Not ready to process step conn_process_8']); return; end
                Y=conn_vol(filename);
                if isfield(Y,'issurface')&&Y.issurface, issurface=true; else issurface=false; end
                if isfield(Y,'conditionsweights')
                    clear X1;
                    X1.conditionweights=Y.conditionsweights;
                else
                    filename=fullfile(filepath,['ROI_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'.mat']);
                    X1=load(filename,'conditionweights');
                    Y.conditionsweights=X1.conditionweights;
                end
                wx=X1.conditionweights{1};
                emptycondition=~nnz(~isnan(wx)&wx~=0);

                if isempty(maxrt), maxrt=max(conn_get_rt(nsub)); end
                DOF=max(0,Y.size.Nt*(min(1/(2*maxrt),CONN_x.Preproc.filter(2))-max(0,CONN_x.Preproc.filter(1)))/(1/(2*maxrt))+1);
                Cy=0;
                Cidx=Y.voxels; %Cidx=[];
                for slice=1:Y.matdim.dim(3),
                    [y,idx]=conn_get_slice(Y,slice);
                    %[xyzx,xyzy]=ind2sub(Y.matdim.dim(1:2),idx);xyz=Y.matdim.mat*[xyzx,xyzy,zeros(size(xyzx))+slice,ones(size(xyzx))]';
                    %z=spm_get_data(Vmask,pinv(Vmask(1).mat)*xyz);
                    %idxt=find(z>THR_MASK&~any(isnan(y),1));
                    %idxt=find(~any(isnan(y),1));
                    idxt=1:size(y,2); y(:,any(isnan(y),1))=0;
                    if ~isempty(idxt),
                        %Cidx=cat(1,Cidx,idx(idxt)+(slice-1)*Y.matdim.dim(1)*Y.matdim.dim(2));
                        y=conn_wdemean(y(:,idxt),wx);
                        y=y.*repmat(wx,[1,size(y,2)]);
                        y=y.*repmat(sqrt(1./max(eps,sum(abs(y).^2,1))),[size(y,1),1]);
                        Cy=Cy+y*y';
                    end
                    n=n+.1;
                    conn_waitbar(n/N,h,sprintf('Subject %d Condition %d',nsub,ncondition));
                end
                [Q1,D]=svd(Cy); % Q1*D*Q1 = Ctt (time-by-time covariance matrix; note:unit variance timeseries)
                D=diag(D);
                DIMS=min([MAXDIMS,ceil(DOF),sum(D>1e-8)]);
                Q1orig=Q1; Dorig=D;
                Q1=Q1(:,1:DIMS);
                D=D(1:DIMS);
                filename_D=fullfile(filepath,['vvPCeig_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'.mat']);
                save(filename_D,'D','Q1','Dorig');
                try
                    temp1=[.7:.1:.9];
                    temp2=[0;cumsum(Dorig)/sum(Dorig)];
                    temp2(temp2==temp2(end)&(1:numel(temp2)<numel(temp2)))=[];
                    temp3=ceil(interp1(temp2,0:numel(temp2)-1,temp1));
                    fileout=fullfile(filepath,['vvPCeig_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'.log']);
                    fh=fopen(fileout,'wt');
                    fprintf(fh,'%d components computed\n',numel(D));
                    for n1=1:length(temp1),fprintf(fh,'retain %d components for a %d%% MSE-fit to the Voxel-to-Voxel correlation matrix\n',temp3(n1),round(temp1(n1)*100));end
                    fclose(fh);
                end
%                 filename_A=fullfile(filepath,['vvPC_Subject',num2str(nsub,'%03d'),'_Condition',num2str(ncondition,'%03d'),'.nii']); % spatial base
%                 Yout_A=struct('fname',filename_A,'mat',Y.matdim.mat,'dim',Y.matdim.dim,'pinfo',[1;0;0],'n',[1,1],'dt',[spm_type('float32'),spm_platform('bigend')],'descrip','conn SVD (spatial base)');
%                 Yout_A=repmat(Yout_A,[DIMS,1]);for nh=1:DIMS,Yout_A(nh).n=[nh,1];end
%                 Yout_A=spm_create_vol(Yout_A);
                filename_A=fullfile(filepath,['vvPC_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'.mat']); % spatial base
                [filename_A,cache]=conn_tempcache(filename_A);
                Yout_A=Y; Yout_A.fname=filename_A;Yout_A.size.Nt=DIMS;Yout_A.DOF=DOF; Yout_A.EmptyData=emptycondition;  %Yout_A.BASE.Q1=Q1;Yout_A.BASE.D=D;
                Yout_A=conn_init_vol(Yout_A,Cidx);
                if ~emptycondition
                    e1=0;e2=0;dX_prev=[];
                    L=zeros([0,3,3]);
                    for slice=1:Y.matdim.dim(3),
                        [y,idx]=conn_get_slice(Y,slice);
                        %[xyzx,xyzy]=ind2sub(Y.matdim.dim(1:2),idx);xyz=Y.matdim.mat*[xyzx,xyzy,zeros(size(xyzx))+slice,ones(size(xyzx))]';
                        %z=spm_get_data(Vmask,pinv(Vmask(1).mat)*xyz);
                        %idxt=find(z>THR_MASK&~any(isnan(y),1));
                        %idxt=find(~any(isnan(y),1));
                        idxt=1:size(y,2); y(:,any(isnan(y),1))=0;
                        if ~isempty(idxt),
                            y=conn_wdemean(y(:,idxt),wx);
                            y=y.*repmat(wx,[1,size(y,2)]);
                            y=y.*repmat(sqrt(1./max(eps,sum(abs(y).^2,1))),[size(y,1),1]);
                            Q2=Q1'*y;
                            e1=e1+Q1orig'*sum(y,2); % for GCOR computation
                            e2=e2+size(y,2);
                            %e=(sum(sum(Q2,2).^2,1)-size(Q2,2))/(size(Q2,2)^2-size(Q2,2)); %gcor
                            %e=sum(abs(y).^2,1);
                            conn_write_slice(Yout_A,Q2,slice); % Y = Q1*sqrtm(D)*R1'; (SVD decomposition)
                                                               % Q2 = Q1'*Y = sqrtm(D)*R1'(components saved in vvPC_Subject*.matc)
                                                               % Cxx = Y'*Y = R1*D*R1' = Q2'*Q2 (voxel-to-voxel correlation matrix)
                                                               % Ctt = Y*Y' = Q1*D*Q1 (time-by-time correlation matrix)
                                                               % Y = Q1*Q2 (time-by-voxel timeseries; unit norm)
                                                               % D = Q2*Q2' (norms saved in vvPCeig)
                            if ~issurface % computes L matrix for smoothness estimation
                                x0=zeros([DIMS,Y.matdim.dim(1:2)]);
                                x0(:,idx)=Q2;
                                dx=(x0(:,[2:end end],:)~=0).*(x0~=0).*(x0(:,[2:end end],:)-x0);
                                dy=(x0(:,:,[2:end end])~=0).*(x0~=0).*(x0(:,:,[2:end end])-x0);
                                if isempty(dX_prev), dz=zeros([size(x0,1),size(x0,2),size(x0,3)]); 
                                else dz=(dX_prev~=0).*(x0~=0).*(dX_prev-x0); 
                                end
                                %dx=(x0(:,2:end,1:end-1)~=0).*(x0(:,1:end-1,1:end-1)~=0).*(x0(:,2:end,1:end-1)-x0(:,1:end-1,1:end-1));
                                %dy=(x0(:,1:end-1,2:end)~=0).*(x0(:,1:end-1,1:end-1)~=0).*(x0(:,1:end-1,2:end)-x0(:,1:end-1,1:end-1));
                                %if isempty(dX_prev), dz=zeros([size(x0,1),size(x0,2)-1,size(x0,3)-1]); 
                                %else dz=(dX_prev(:,1:end-1,1:end-1)~=0).*(x0(:,1:end-1,1:end-1)~=0).*(dX_prev(:,1:end-1,1:end-1)-x0(:,1:end-1,1:end-1)); 
                                %end
                                dx=dx(:,idx)';
                                dy=dy(:,idx)';
                                dz=dz(:,idx)';
                                tL=zeros([numel(idx),3,3]);
                                tL(:,1,1) = sum(dx.*dx,2); % sum over components
                                tL(:,1,2) = sum(dx.*dy,2);
                                tL(:,1,3) = sum(dx.*dz,2);
                                tL(:,2,2) = sum(dy.*dy,2);
                                tL(:,2,3) = sum(dy.*dz,2);
                                tL(:,3,3) = sum(dz.*dz,2);
                                L=cat(1,L,tL);
                                dX_prev=x0;
                            end
                        end
                        n=n+1;
                        conn_waitbar(n/N,h,sprintf('Subject %d Condition %d',nsub,ncondition));
                    end
                    gcor=(sum(e1.^2,1)-e2)/(e2^2-e2); % note: removing diagonal 1's to avoid positive-only constrain (gcor_full=sum(e1.^2,1)/(e2^2);)
                    gcor_name=['QC_GCOR_',CONN_x.Setup.conditions.names{ncondition}];
                    gcor_icov=find(strcmp(gcor_name,CONN_x.Setup.l2covariates.names(1:end-1)),1);
                    if isempty(gcor_icov),
                        gcor_icov=numel(CONN_x.Setup.l2covariates.names);
                        CONN_x.Setup.l2covariates.names{gcor_icov}=gcor_name;
                        CONN_x.Setup.l2covariates.descrip{gcor_icov}=['CONN Quality Assurance: Global Correlation @ ',CONN_x.Setup.conditions.names{ncondition}];
                        CONN_x.Setup.l2covariates.names{gcor_icov+1}=' ';
                        for tnsub=1:CONN_x.Setup.nsubjects, CONN_x.Setup.l2covariates.values{tnsub}{gcor_icov}=nan; end
                    end
                    CONN_x.Setup.l2covariates.values{nsub}{gcor_icov}=gcor; 
                    if ~issurface % note: only used in volume-level analyses (disregard for surfaces)
                        resel_xyz = [L(:,1,1) L(:,2,2)  L(:,3,3)]; % note: sum(Dorig) = size(L,1)
                        resel_img = L(:,1,1).*L(:,2,2).*L(:,3,3) + ...
                            L(:,1,2).*L(:,2,3).*L(:,1,3)*2 - ...
                            L(:,1,1).*L(:,2,3).*L(:,2,3) - ...
                            L(:,1,2).*L(:,1,2).*L(:,3,3) - ...
                            L(:,1,3).*L(:,2,2).*L(:,1,3);
                        resel_img(resel_img<0) = 0;
                        % Convert det(Lambda) and diag(Lambda) to units of resels
                        resel_img = sqrt(resel_img/(4*log(2))^3);
                        resel_xyz = sqrt(resel_xyz/(4*log(2)));
                        i = any(isnan(resel_xyz),2)|any(resel_xyz==0,2);
                        resel_img = mean(resel_img(~i,:),1);
                        resel_xyz = mean(resel_xyz(~i,:),1);
                        L = mean(L(~i,:,:),1); % note: change this and keep original L for future voxel-specific RPV computations
                        RESEL = resel_img^(1/3)*(resel_xyz/prod(resel_xyz)^(1/3));
                        FWHM  = full(sparse(1,1:3,1./RESEL,1,3));
                        FWHM(isnan(FWHM)) = 0;
                        FWHM(~FWHM) = 1;
                        filename_D=fullfile(filepath,['vvPCeig_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'.mat']);
                        save(filename_D,'L','resel_img','resel_xyz','FWHM','-append');
                    end
                    
                    if 0 % QC_LCOR computation (this implementation takes too long, preferable slice-serial computations)
                        if issurface, localsupport=25; % kernel size for LCOR 
                        else localsupport=25./sqrt(sum(abs(Y.matdim.mat(1:3,1:3)).^2,1));
                        end
                        e3=0;
                        for ndim=1:Yout_A.size.Nt, %DIMS
                            [q2,idxq2]=conn_get_time(Yout_A,ndim);
                            q3=conn_conv(q2,localsupport,[],issurface); %%%
                            e3=e3+q2(:).*q3(:); % for LCOR computation
                        end
                        lcor=sum(e3)/e2;
                        lcor_name=['QC_LCOR_',CONN_x.Setup.conditions.names{ncondition}];
                        lcor_icov=find(strcmp(lcor_name,CONN_x.Setup.l2covariates.names(1:end-1)),1);
                        if isempty(lcor_icov),
                            lcor_icov=numel(CONN_x.Setup.l2covariates.names);
                            CONN_x.Setup.l2covariates.names{lcor_icov}=lcor_name;
                            CONN_x.Setup.l2covariates.descrip{lcor_icov}=['CONN Quality Assurance: Local Correlation @ ',CONN_x.Setup.conditions.names{ncondition}];
                            CONN_x.Setup.l2covariates.names{lcor_icov+1}=' ';
                            for tnsub=1:CONN_x.Setup.nsubjects, CONN_x.Setup.l2covariates.values{tnsub}{lcor_icov}=nan; end
                        end
                        CONN_x.Setup.l2covariates.values{nsub}{lcor_icov}=lcor;
                    end
                end
                conn_tempcache(cache,'matc');
            else, 
                n=n+1.1*Y.matdim.dim(3); 
                conn_waitbar(n/N,h,sprintf('Subject %d Condition %d',nsub,ncondition));
            end
        end
    end
    conn_waitbar('close',h);
    
    FWHM=0; % FWHM (mm) for spheres in post-hoc connectome display (from voxel-to-voxel analyses) 
             % (set to 0 for using individual voxels instead of spheres)
    REDO=[];
    filepath=CONN_x.folders.preprocessing;
    nconditions=length(CONN_x.Setup.conditions.names)-1;
    if FWHM>0
        h=conn_waitbar(0,['Step ',num2str(sum(options<=8)),'/',num2str(length(options)),': smoothing voxel-to-voxel covariance']);
        for ncondition=validconditions,
            for nsub=validsubjects
                filename=fullfile(filepath,['svvPC_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'.mat']);
                if isempty(REDO)&&conn_existfile(filename),
                    if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'overwrite'), REDO=CONN_x.gui.overwrite; else, REDO=conn_questdlg('Overwrite existing subject results?','','Yes', 'No', 'Yes');end;
                end
                if strcmp(lower(REDO),'yes')||(~conn_existfile(filename)),
                    filename_B1=fullfile(filepath,['vvPC_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'.mat']);
                    Y1=conn_vol(filename_B1);
                    [X,idx]=conn_get_volume(Y1);
                    fwhm=FWHM./sqrt(sum(abs(Y1.matdim.mat(1:3,1:3)).^2,1));
                    X2=X;
                    x0=zeros(Y1.matdim.dim);
                    x0(idx)=1;
                    x0=conn_conv(x0,fwhm);
                    for nt=1:Y1.size.Nt,
                        x=zeros(Y1.matdim.dim);
                        x(idx)=X(nt,:);
                        x=conn_conv(x,fwhm)./max(eps,x0);
                        X2(nt,:)=x(idx);
                    end
                    Y2=Y1;
                    Y2.fname=fullfile(filepath,['svvPC_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'.mat']);
                    Y2=conn_init_vol(Y2,[],X2);
                    %conn_waitbar(((ncondition-1)*CONN_x.Setup.nsubjects+nsub)/nconditions/CONN_x.Setup.nsubjects,h);
                end
            end
        end
        conn_waitbar('close',h);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% updates CONN_x.Analyses structure with all ROIs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(floor(options)==9),
    if any(options==9.1), h=conn_waitbar(0,['Step ',num2str(sum(floor(options)<=9)),'/',num2str(length(options)),': Checking/Updating S2V/R2R Analysis variables']);
    elseif any(options==9.2), h=conn_waitbar(0,['Step ',num2str(sum(floor(options)<=9)),'/',num2str(length(options)),': Checking/Updating V2V Analysis variables']);
    elseif any(options==9.3), h=conn_waitbar(0,['Step ',num2str(sum(floor(options)<=9)),'/',num2str(length(options)),': Checking/Updating dyn-ICA Analysis variables']);
    else h=conn_waitbar(0,['Step ',num2str(sum(floor(options)<=9)),'/',num2str(length(options)),': Updating Analysis variables']);
    end
    nconditions=length(CONN_x.Setup.conditions.names)-1;
    if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'conditions')&&~isempty(CONN_x.gui.conditions), validconditions=CONN_x.gui.conditions; else validconditions=1:length(CONN_x.Setup.conditions.names)-1; end
    validconditions=validconditions(cellfun('length',CONN_x.Setup.conditions.model(validconditions))==0); 
    icondition=[];isnewcondition=[];for ncondition=validconditions,[icondition(ncondition),isnewcondition(ncondition)]=conn_conditionnames(CONN_x.Setup.conditions.names{ncondition}); end
    if any(isnewcondition(validconditions)), error(['Some conditions have not been processed yet. Re-run previous step']); end
    [path,name,ext]=fileparts(CONN_x.filename);
    filepath=CONN_x.folders.preprocessing;
    if any(CONN_x.Setup.steps([1,2,4])) && any(options==9.1|options==9) && ~(isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'steps')&&~any(CONN_x.gui.steps([1,2,4])))
        if nargin>1&&~isempty(varargin{1}),analyses=varargin{1}; % selected analysis only
        else analyses=1:length(CONN_x.Analyses);
        end
        if ischar(analyses)||iscell(analyses), analyses=find(ismember({CONN_x.Analyses.name},analyses)); end
        analyses(analyses<1|analyses>numel(CONN_x.Analyses))=[];
        validsubjects=1:CONN_x.Setup.nsubjects; %if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'subjects'), validsubjects=CONN_x.gui.subjects; else validsubjects=1:CONN_x.Setup.nsubjects; end
        if isfield(CONN_x,'pobj')&&isstruct(CONN_x.pobj)&&isfield(CONN_x.pobj,'subjects'), validsubjects=CONN_x.pobj.subjects; conn_projectmanager('addstep',9.1,analyses); end
        missingdata=~conn_existfile(conn_fullfile(filepath,arrayfun(@(n)['ROI_Subject',num2str(validsubjects(1),'%03d'),'_Condition',num2str(icondition(n),'%03d'),'.mat'],validconditions,'uni',0)));
        if any(missingdata), conn_disp(['Not ready to process step conn_process_9']); return; end
        filename1=fullfile(filepath,['ROI_Subject',num2str(validsubjects(1),'%03d'),'_Condition',num2str(icondition(validconditions(1)),'%03d'),'.mat']);
        x1=conn_fcnutils('loadnamesandsize',filename1); % x1=conn_loadmatfile(filename1);        
        CONN_x.Analysis_variables.names={};
        CONN_x.Analysis_variables.types={};
        CONN_x.Analysis_variables.deriv={};
        CONN_x.Analysis_variables.fbands={};
        CONN_x.Analysis_variables.dimensions={};
        for n1=1:length(x1.names),
            idx=strmatch(x1.names{n1},CONN_x.Preproc.confounds.names,'exact');
            %if isempty(idx),
            %    idx=strmatch(x1.names{n1},CONN_x.Preproc.confounds.names);  % allows partial-name matches
            %    if numel(idx)~=1, idx=[]; end
            %end
            if isempty(idx),
                CONN_x.Analysis_variables.names{end+1}=x1.names{n1};
                CONN_x.Analysis_variables.types{end+1}='roi';
                CONN_x.Analysis_variables.deriv{end+1}=0;
                CONN_x.Analysis_variables.fbands{end+1}=1;
                CONN_x.Analysis_variables.dimensions{end+1}=[x1.dimensions{n1},x1.dimensions{n1}];
            end
        end
        analysisbak=CONN_x.Analysis;
        for ianalysis=analyses,
            CONN_x.Analysis=ianalysis;
            if isempty(CONN_x.Analyses(ianalysis).name),CONN_x.Analyses(ianalysis).name=['SBC_',num2str(ianalysis,'%02d')]; end;
            if ~conn_existfile(fullfile(CONN_x.folders.firstlevel,CONN_x.Analyses(ianalysis).name),2), try, conn_fileutils('mkdir',CONN_x.folders.firstlevel,CONN_x.Analyses(ianalysis).name); end; end;
            if isfield(CONN_x.Analyses(ianalysis).regressors,'names') && ~isempty(CONN_x.Analyses(ianalysis).regressors.names), initial=CONN_x.Analyses(ianalysis).regressors.names; dims=CONN_x.Analyses(ianalysis).regressors.dimensions; ders=CONN_x.Analyses(ianalysis).regressors.deriv; fbands=CONN_x.Analyses(ianalysis).regressors.fbands;
            else, initial=CONN_x.Analysis_variables.names; %(strcmp(CONN_x.Analysis_variables.types,'roi')); 
                initial=initial(~cellfun('length',regexp(initial,'^Grey Matter$|^White Matter$|^CSF Matter$|^QA_|^QC_|^Effect of'))); dims={}; ders={}; fbands={}; end
            CONN_x.Analyses(ianalysis).regressors.names={};
            CONN_x.Analyses(ianalysis).regressors.types={};
            CONN_x.Analyses(ianalysis).regressors.deriv={};
            CONN_x.Analyses(ianalysis).regressors.fbands={};
            CONN_x.Analyses(ianalysis).regressors.dimensions={};
            for n1=1:length(initial),
                idx=strmatch(initial{n1},CONN_x.Analysis_variables.names,'exact');
                if isempty(idx),
                    idx=strmatch(initial{n1},CONN_x.Analysis_variables.names); % allows partial-name matches
                    %if numel(idx)~=1, idx=[]; end % allows multiple partial-name matches
                end
                for idx=idx(:)' %if ~isempty(idx),%&&~strcmp(initial{n1},'Grey Matter')&&~strcmp(initial{n1},'White Matter')&&~strcmp(initial{n1},'CSF'),
                    CONN_x.Analyses(ianalysis).regressors.names{end+1}=CONN_x.Analysis_variables.names{idx};
                    CONN_x.Analyses(ianalysis).regressors.types{end+1}=CONN_x.Analysis_variables.types{idx};
                    if length(ders)>=n1&&~isempty(ders{n1}), CONN_x.Analyses(ianalysis).regressors.deriv{end+1}=ders{n1}; else, CONN_x.Analyses(ianalysis).regressors.deriv{end+1}=CONN_x.Analysis_variables.deriv{idx};end
                    if length(fbands)>=n1&&~isempty(fbands{n1}), CONN_x.Analyses(ianalysis).regressors.fbands{end+1}=fbands{n1}; else, CONN_x.Analyses(ianalysis).regressors.fbands{end+1}=CONN_x.Analysis_variables.fbands{idx};end
                    if length(dims)>=n1&&~isempty(dims{n1}), CONN_x.Analyses(ianalysis).regressors.dimensions{end+1}=[min(dims{n1}(1),CONN_x.Analysis_variables.dimensions{idx}(1)),CONN_x.Analysis_variables.dimensions{idx}(1)];
                    else, CONN_x.Analyses(ianalysis).regressors.dimensions{end+1}=CONN_x.Analysis_variables.dimensions{idx}; end
                end
            end
            if ~isfield(CONN_x.Analyses(ianalysis),'modulation') || isempty(CONN_x.Analyses(ianalysis).modulation), CONN_x.Analyses(ianalysis).modulation=0; end
            if ~isfield(CONN_x.Analyses(ianalysis),'measure') || isempty(CONN_x.Analyses(ianalysis).measure), CONN_x.Analyses(ianalysis).measure=1; end
            if ~isfield(CONN_x.Analyses(ianalysis),'weight') || isempty(CONN_x.Analyses(ianalysis).weight), CONN_x.Analyses(ianalysis).weight=2; end
            if ~isfield(CONN_x.Analyses(ianalysis),'type') || isempty(CONN_x.Analyses(ianalysis).type), CONN_x.Analyses(ianalysis).type=3; end
        end
        CONN_x.Analysis=analysisbak;
    end
    conn_waitbar(1/3,h);
    if any(CONN_x.Setup.steps([3])) && any(options==9.2|options==9) && ~(isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'steps')&&~any(CONN_x.gui.steps([3])))
        if nargin>1&&~isempty(varargin{1}),analyses=varargin{1}; % selected analysis only
        else analyses=1:length(CONN_x.vvAnalyses);
        end
        if ischar(analyses)||iscell(analyses), analyses=find(ismember({CONN_x.vvAnalyses.name},analyses)); end
        analyses(analyses<1|analyses>numel(CONN_x.vvAnalyses))=[];
        validsubjects=1:CONN_x.Setup.nsubjects; %if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'subjects'), validsubjects=CONN_x.gui.subjects; else validsubjects=1:CONN_x.Setup.nsubjects; end
        if isfield(CONN_x,'pobj')&&isstruct(CONN_x.pobj)&&isfield(CONN_x.pobj,'subjects'), validsubjects=CONN_x.pobj.subjects; conn_projectmanager('addstep',9.2,analyses); end
        analysisbak=CONN_x.vvAnalysis;
        for ianalysis=analyses,
            CONN_x.vvAnalysis=ianalysis;
            if ~conn_existfile(fullfile(CONN_x.folders.firstlevel_vv,CONN_x.vvAnalyses(ianalysis).name),2), try, conn_fileutils('mkdir',CONN_x.folders.firstlevel_vv,CONN_x.vvAnalyses(ianalysis).name); end; end;
            
            CONN_x.vvAnalyses(ianalysis).variables=conn_v2v('measures');
            if isfield(CONN_x.vvAnalyses(ianalysis).regressors,'names') && ~isempty(CONN_x.vvAnalyses(ianalysis).regressors.names),
                initial=CONN_x.vvAnalyses(ianalysis).regressors;
            else % default analyses
                initial=CONN_x.vvAnalyses(ianalysis).variables;
                idx=find(strcmp(initial.names,'group-MVPA'),1);
                if ~isempty(idx)
                    optionsnames=fieldnames(initial);
                    for n2=1:numel(optionsnames), initial.(optionsnames{n2})=initial.(optionsnames{n2})(idx); end
                    initial.dimensions_out{1}=max(1,min(20,round(CONN_x.Setup.nsubjects/5)));
                end
            end
            CONN_x.vvAnalyses(ianalysis).regressors=conn_v2v('empty');
            for n1=1:numel(initial.names)
                idx=strmatch(initial.names{n1},CONN_x.vvAnalyses(ianalysis).variables.names,'exact');
                if isempty(idx)&&isfield(CONN_x.vvAnalyses(ianalysis).variables,'alt_names'),
                    for n=1:numel(CONN_x.vvAnalyses(ianalysis).variables.alt_names), % allows alternative-name matches
                        if ~isempty(CONN_x.vvAnalyses(ianalysis).variables.alt_names{n})
                            idx=strmatch(initial.names{n1},CONN_x.vvAnalyses(ianalysis).variables.alt_names{n},'exact');
                            if numel(idx)==1, idx=n; break;
                            else idx=[];
                            end
                        end
                    end
                end
                if isempty(idx),
                    idx=strmatch(initial.names{n1},CONN_x.vvAnalyses(ianalysis).variables.names); % allows partial-name matches
                    if numel(idx)~=1, idx=[]; conn_disp(['warning: no match for voxel-to-voxel analysis named ',initial.names{n1},'. Skipping analysis...']); end
                end
                optionsnames=fieldnames(CONN_x.vvAnalyses(ianalysis).variables);
                for n2=1:numel(optionsnames),
                    if ~isempty(idx)&&(~isfield(initial,optionsnames{n2})||numel(initial.(optionsnames{n2}))<n1||isempty(initial.(optionsnames{n2}){n1})), CONN_x.vvAnalyses(ianalysis).regressors.(optionsnames{n2}){n1}=CONN_x.vvAnalyses(ianalysis).variables.(optionsnames{n2}){idx(1)};
                    else CONN_x.vvAnalyses(ianalysis).regressors.(optionsnames{n2}){n1}=initial.(optionsnames{n2}){n1};
                    end
                end
                %if isempty(idx), CONN_x.vvAnalyses(ianalysis).regressors.measuretype{n1}=1; end
            end
        end
        CONN_x.vvAnalysis=analysisbak;
    end
    conn_waitbar(2/3,h);
    if any(CONN_x.Setup.steps([4])) && any(options==9.3|options==9) && ~(isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'steps')&&~any(CONN_x.gui.steps([4])))
        if nargin>1&&~isempty(varargin{1}),analyses=varargin{1}; % selected analysis only
        else analyses=1:length(CONN_x.dynAnalyses);
        end
        if ischar(analyses)||iscell(analyses), analyses=find(ismember({CONN_x.dynAnalyses.name},analyses)); end
        analyses(analyses<1|analyses>numel(CONN_x.dynAnalyses))=[];
        validsubjects=1:CONN_x.Setup.nsubjects; %if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'subjects'), validsubjects=CONN_x.gui.subjects; else validsubjects=1:CONN_x.Setup.nsubjects; end
        if isfield(CONN_x,'pobj')&&isstruct(CONN_x.pobj)&&isfield(CONN_x.pobj,'subjects'), validsubjects=CONN_x.pobj.subjects; conn_projectmanager('addstep',9.3,analyses); end
        %CONN_x.dynAnalyses.variables.names=CONN_x.Analyses(analyses(1)).variables.names;
        filename1=fullfile(filepath,['ROI_Subject',num2str(validsubjects(1),'%03d'),'_Condition',num2str(0,'%03d'),'.mat']);
        x1=conn_loadmatfile(filename1,'names'); %,'conditionsnames');
        analysisbak=CONN_x.dynAnalysis;
        for ianalysis=analyses,
            CONN_x.dynAnalysis=ianalysis;
            if ~conn_existfile(fullfile(CONN_x.folders.firstlevel_dyn,CONN_x.dynAnalyses(ianalysis).name),2), try, conn_fileutils('mkdir',CONN_x.folders.firstlevel_dyn,CONN_x.dynAnalyses(ianalysis).name); end; end;
            CONN_x.dynAnalyses(CONN_x.dynAnalysis).variables.names={};
            for n1=1:length(x1.names),
                idx=strmatch(x1.names{n1},CONN_x.Preproc.confounds.names,'exact');
                %if isempty(idx),
                %    idx=strmatch(x1.names{n1},CONN_x.Preproc.confounds.names);  % allows partial-name matches
                %    if numel(idx)~=1, idx=[]; end
                %end
                if isempty(idx),
                    CONN_x.dynAnalyses(CONN_x.dynAnalysis).variables.names{end+1}=x1.names{n1};
                end
            end
            if isfield(CONN_x.dynAnalyses(CONN_x.dynAnalysis).regressors,'names') && ~isempty(CONN_x.dynAnalyses(CONN_x.dynAnalysis).regressors.names), initial=CONN_x.dynAnalyses(CONN_x.dynAnalysis).regressors.names;
            else, initial=CONN_x.dynAnalyses(CONN_x.dynAnalysis).variables.names; initial=initial(~cellfun('length',regexp(initial,'^Grey Matter$|^White Matter$|^CSF Matter$|^QA_|^QC_'))); end
            CONN_x.dynAnalyses(CONN_x.dynAnalysis).regressors.names={};
            for n1=1:length(initial),
                idx=strmatch(initial{n1},CONN_x.dynAnalyses(CONN_x.dynAnalysis).variables.names,'exact');
                if isempty(idx),
                    idx=strmatch(initial{n1},CONN_x.dynAnalyses(CONN_x.dynAnalysis).variables.names); % allows partial-name matches
                    %if numel(idx)~=1, idx=[]; end
                end
                for idx=idx(:)' %if ~isempty(idx),%&&~strcmp(initial{n1},'Grey Matter')&&~strcmp(initial{n1},'White Matter')&&~strcmp(initial{n1},'CSF'),
                    CONN_x.dynAnalyses(CONN_x.dynAnalysis).regressors.names{end+1}=CONN_x.dynAnalyses(CONN_x.dynAnalysis).variables.names{idx};
                end
            end
        end
        CONN_x.dynAnalysis=analysisbak;
    end
    conn_waitbar(3/3,h);
    conn_waitbar('close',h);
    CONN_x.isready(3)=1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates resultsDATA_Subject###_Condition###_Source###.mat files (first-level ROI-to-voxel analysis)
% Creates BETA_Subject###_Condition###_Source###.mat files (first-level seed-to-voxel analyses)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(options==10) && any(CONN_x.Setup.steps([2])) && ~(isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'steps')&&~any(CONN_x.gui.steps([2]))),
    warning('off','MATLAB:DELETE:FileNotFound');
    [path,name,ext]=fileparts(CONN_x.filename);
    filepath=CONN_x.folders.preprocessing;
    if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'subjects'), validsubjects=CONN_x.gui.subjects; else validsubjects=1:CONN_x.Setup.nsubjects; end
    if isfield(CONN_x,'pobj')&&isstruct(CONN_x.pobj)&&isfield(CONN_x.pobj,'subjects'), validsubjects=CONN_x.pobj.subjects; end
    nconditions=length(CONN_x.Setup.conditions.names)-1;
    if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'conditions')&&~isempty(CONN_x.gui.conditions), validconditions=CONN_x.gui.conditions; else validconditions=1:length(CONN_x.Setup.conditions.names)-1; end
    icondition=[];isnewcondition=[];for ncondition=1:nconditions,[icondition(ncondition),isnewcondition(ncondition)]=conn_conditionnames(CONN_x.Setup.conditions.names{ncondition}); end
    secondaryconditions=validconditions(cellfun('length',CONN_x.Setup.conditions.model(validconditions))>0);
    validconditions=validconditions(cellfun('length',CONN_x.Setup.conditions.model(validconditions))==0); 
    if any(isnewcondition(validconditions)), error(['Some conditions have not been processed yet. Re-run previous step']); end
    missingdata=arrayfun(@(n)isempty(dir(fullfile(filepath,['ROI_Subject',num2str(validsubjects(1),'%03d'),'_Condition',num2str(icondition(n),'%03d'),'.mat']))),validconditions);
    if any(missingdata), conn_disp(['Not ready to process step conn_process_10']); return; end
    if ~isempty(validconditions), referenceconditions=validconditions;
    else
        evaluatefunction=CONN_x.Setup.conditions.model{secondaryconditions(1)}{1};
        if isequal(evaluatefunction,'lin'), primaryconditionsnames=CONN_x.Setup.conditions.model{secondaryconditions(1)}(3:end);
        else primaryconditionsnames=CONN_x.Setup.conditions.model{secondaryconditions(1)}(2:end);
        end
        referenceconditions=find(ismember(CONN_x.Setup.conditions.names(1:end-1),primaryconditionsnames));
    end
    filename=fullfile(filepath,['ROI_Subject',num2str(validsubjects(1),'%03d'),'_Condition',num2str(icondition(referenceconditions(1)),'%03d'),'.mat']);
    X1=load(filename);
    if nargin>1&&~isempty(varargin{1}),analyses=varargin{1}; % selected analysis only
    else, analyses=1:length(CONN_x.Analyses); 
    end
    if ischar(analyses)||iscell(analyses), analyses=find(ismember({CONN_x.Analyses.name},analyses)); end
    analyses(analyses<=0|analyses>numel(CONN_x.Analyses))=[];
    doanalyses=false(1,numel(analyses));
    for nanalyses=1:numel(analyses),if analyses(nanalyses)>0&&any(CONN_x.Analyses(analyses(nanalyses)).type==[2,3]), doanalyses(nanalyses)=true; end; end
    analyses=analyses(doanalyses);
    h=conn_waitbar(0,['Step ',num2str(sum(options<=10)),'/',num2str(length(options)),': ROI-to-voxel first-level analyses']);
    REDO=[];
    analysisbak=CONN_x.Analysis;
    for nanalyses=1:length(analyses),
        redone_files=0;
        ianalysis=analyses(nanalyses);
        CONN_x.Analysis=ianalysis;
        names={};
%         ianalysis=CONN_x.Analysis;
        filepathresults=fullfile(CONN_x.folders.firstlevel,CONN_x.Analyses(ianalysis).name);
        conn_disp('fprintf','      first-level data in %s\n',filepathresults);
        [X,nill,names]=conn_designmatrix(CONN_x.Analyses(ianalysis).regressors,X1,[]);
        nrois=size(X,2)-1;
        iroi=[];isnew=[];for nroi=1:nrois,[iroi(nroi),isnew(nroi)]=conn_sourcenames(names{nroi},'+');end
        isnew0=isnew;
        
        filename=fullfile(filepath,['DATA_Subject',num2str(validsubjects(1),'%03d'),'_Condition',num2str(icondition(referenceconditions(1)),'%03d'),'.mat']);
        Y=conn_vol(filename);
        if nanalyses==1, N=numel(validsubjects)*(numel(validconditions)+numel(secondaryconditions))*Y.matdim.dim(3)*length(analyses);n=0;end

        ConditionWeights={};
        if ~ischar(CONN_x.Analyses(ianalysis).modulation)&&CONN_x.Analyses(ianalysis).modulation>0
            if isempty(CONN_x.Analyses(ianalysis).conditions), validconditions2=validconditions;
            else validconditions2=find(ismember(CONN_x.Setup.conditions.names(1:end-1),CONN_x.Analyses(ianalysis).conditions));
            end
            for ncondition=union(validconditions,validconditions2),
                for nsub=validsubjects,
                    filename=fullfile(filepath,['ROI_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'.mat']);
                    X1=load(filename,'conditionweights');
                    for n1=1:numel(X1.conditionweights)
                        ConditionWeights{nsub,n1}(:,ncondition)=X1.conditionweights{n1};
                    end
                end
            end
        end
        for nsub=validsubjects,
            touched=false(length(CONN_x.Setup.conditions.names)-1,nrois);
            maxrt=[];
            for ncondition=validconditions,
                filename=fullfile(filepath,['DATA_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'.mat']);
                Y=conn_vol(filename);
                filename=fullfile(filepath,['ROI_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'.mat']);
                X1=load(filename);
                [X,nill,names]=conn_designmatrix(CONN_x.Analyses(ianalysis).regressors,X1,[]);
                if nrois~=size(X,2)-1,
                    error('Incorrect number of ROI components for subject %d condition %d. Please re-run Setup&Denoising step (all subjects/seeds)',nsub,ncondition);
                end
                clear Yout cache;
                isnew=isnew0;
                for nroi=1:nrois,
                    filename=fullfile(filepathresults,['BETA_Subject',num2str(nsub,CONN_x.opt.fmt1),'_Condition',num2str(icondition(ncondition),'%03d'),'_Source',num2str(iroi(nroi),'%03d'),'.nii']);
                    %filename=fullfile(filepathresults,['resultsDATA_Subject',num2str(nsub,'%03d'),'_Condition',num2str(ncondition,'%03d'),'_Source',num2str(iroi(nroi),'%03d'),'.mat']);
                    isnew(nroi)=isnew(nroi)|~conn_existfile(filename);
                end
                if any(isnew),
                    switch(CONN_x.Analyses(ianalysis).measure),
                        case {2,4}, %partial
                            isnew=ones(size(isnew));
                    end
                end
                if isempty(REDO)&&~all(isnew),
                    if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'overwrite'), REDO=CONN_x.gui.overwrite; else, REDO=conn_questdlg('Overwrite existing subject results?','','Yes', 'No', 'Yes');end;
                end
                if strcmp(lower(REDO),'no'),
                    idxredo=find(isnew);
                    X=cat(2,X(:,1),X(:,1+idxredo));
                else,
                    idxredo=1:nrois;
                end
                if length(idxredo)>0,
                    nX=size(X,2);
                    wx=ones(size(X,1),1);
                    switch(CONN_x.Analyses(ianalysis).weight),
                        case 1, wx=double(X1.conditionweights{1}>0);                            % none
                        case 2, wx=X1.conditionweights{1};                                      % hrf
                        case 3, wx=X1.conditionweights{2};                                      % hanning
                        case 4, wx=X1.conditionweights{3};                                      % condition-specific temporal weights
                    end
                    if ~(ischar(CONN_x.Analyses(ianalysis).modulation)||CONN_x.Analyses(ianalysis).modulation>0) % Weighted GLM
                        if any(wx<0), conn_disp('Warning: temporal weights have negative values; cropping at 0'); end
                        wx=max(0,wx);
                        X=cat(2,X(:,1),conn_wdemean(X(:,2:end),wx));
                        X=X.*repmat(wx,[1,size(X,2)]);
%                         if CONN_x.Analyses(ianalysis).measure<3, 
%                             conn_disp('Warning: correlation measure not recommended for PPI analyses'); 
%                         end 
                    else % parametric modulation
                        if ~ischar(CONN_x.Analyses(ianalysis).modulation), % PPI
                            %wx=X1.conditionweights{3};   %PPI
                            wx=ConditionWeights{nsub,3}(:,[setdiff(validconditions2,ncondition) ncondition]); %gPPI
                        elseif ~isempty(regexp(CONN_x.Analyses(ianalysis).modulation,'^(.*\/|.*\\)?Dynamic factor \d+$')), % Dyn FC
                            CONN_x.Analyses(ianalysis).modulation=regexprep(CONN_x.Analyses(ianalysis).modulation,'\\|\/',['\',filesep]);
                            [name1,name2]=fileparts(CONN_x.Analyses(ianalysis).modulation);
                            [ok,value1]=ismember(name1,{CONN_x.dynAnalyses.name}); if ~ok, if numel(CONN_x.dynAnalyses)==1, value1=1; else error('Analysis name %s not found',name1); end; end
                            filename=fullfile(CONN_x.folders.firstlevel_dyn,CONN_x.dynAnalyses(value1).name,['dyn_Subject',num2str(nsub,'%03d'),'.mat']);
                            xmod=load(filename);
                            [ok,idx]=ismember(name2,xmod.names);
                            if ok, wx=xmod.data(:,[setdiff(1:size(xmod.data,2),idx),idx]);
                            else error('Temporal factor not found');
                            end
                            wx=conn_bsxfun(@times,wx,X1.conditionweights{1}>0);
                        else % Physio-Physiological interaction
                            idx=find(strcmp(CONN_x.Analyses(ianalysis).modulation,X1.names));
                            if numel(idx)==1, wx=X1.data{idx};
                            elseif isempty(idx), error('Covariate not found. Please re-run ''Denoising'' step');
                            else,
                                idx=find(cellfun(@(x)all(isnan(x)),X1.xyz));
                                idx=idx(strcmp(CONN_x.Analyses(ianalysis).modulation,X1.names(idx)));
                                if numel(idx)==1, wx=X1.data{idx};
                                else error('Covariate not found');
                                end
                            end
                            wx=conn_bsxfun(@times,wx,X1.conditionweights{1}>0);
                        end
                        %if size(wx,2)>1, conn_disp('Warning: multivariate interaction term not supported. Summing interaction term across multiple components'); end
                        inter=wx;
                        X=[X(:,1) detrend([X(:,2:end) reshape(repmat(permute(inter,[1 3 2]),[1,size(X,2),1]),size(X,1),[]) reshape(conn_bsxfun(@times,X,permute(inter,[1 3 2])),size(X,1),[])],'constant')];
                    end
                    nVars=size(X,2)/nX;
                    if isempty(maxrt), maxrt=max(conn_get_rt(nsub)); end
                    if ischar(CONN_x.Analyses(ianalysis).modulation)||CONN_x.Analyses(ianalysis).modulation>0 % parametric modulation
                        switch(CONN_x.Analyses(ianalysis).measure),
                            case {1,3}, %bivariate
                                Xtemp=permute(reshape(X,[size(X,1),nX,nVars]),[1,3,2]);
                                iX=sparse(nVars*nX,nVars*nX);
                                for nXtemp=1:size(Xtemp,3), iXtemp=pinv(Xtemp(:,:,nXtemp)'*Xtemp(:,:,nXtemp)); iX(nXtemp:nX:end,nXtemp:nX:end)=iXtemp; end
                                %iX=pinv((X'*X).*kron(ones(nVars),eye(numel(idxredo)+1)));
                                DOF=max(0,Y.size.Nt*(min(1/(2*maxrt),CONN_x.Preproc.filter(2))-max(0,CONN_x.Preproc.filter(1)))/(1/(2*maxrt))-nVars);
                            case {2,4}, %partial
                                iX=pinv(X'*X);
                                DOF=max(0,Y.size.Nt*(min(1/(2*maxrt),CONN_x.Preproc.filter(2))-max(0,CONN_x.Preproc.filter(1)))/(1/(2*maxrt))-rank(X)+1);
                        end
                        r=sqrt(diag(iX));
                        if 0&&~ischar(CONN_x.Analyses(ianalysis).modulation) % gPPI absolute values (physiological+PPI)
                            iX=iX(1:nX,:)+iX((nVars-1)*nX+1:end,:);
                        else
                            iX=iX((nVars-1)*nX+1:end,:);
                        end
                        r=r((nVars-1)*nX+1:end);
                    else % standard functional connectivity
                        switch(CONN_x.Analyses(ianalysis).measure),
                            case {1,3}, %bivariate
                                iX=diag(1./max(eps,sum(X.^2,1)));
                                %iX=pinv(diag(diag(X'*X)));
                                DOF=max(0,Y.size.Nt*(min(1/(2*maxrt),CONN_x.Preproc.filter(2))-max(0,CONN_x.Preproc.filter(1)))/(1/(2*maxrt))-1);
                            case {2,4}, %partial
                                iX=pinv(X'*X);
                                DOF=max(0,Y.size.Nt*(min(1/(2*maxrt),CONN_x.Preproc.filter(2))-max(0,CONN_x.Preproc.filter(1)))/(1/(2*maxrt))-rank(X)+1);
                        end
                        r=sqrt(diag(iX));
                    end
                    filename=fullfile(filepathresults,['se_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'.nii']);
                    SEout=struct('mat',Y.matdim.mat,'dim',Y.matdim.dim,'fname',filename,'dt',[spm_type('float32') spm_platform('bigend')]);
                    conn_fileutils('deletefile',SEout.fname);
                    SEout=spm_create_vol(SEout);
%                     filename=fullfile(filepathresults,['seDATA_Subject',num2str(nsub,'%03d'),'_Condition',num2str(ncondition,'%03d'),'.mat']);
%                     SEout=Y; SEout.fname=filename;
%                     SEout.size.Nt=1;%CONN_x.Setup.nsubjects;
%                     SEout.DOF=DOF;
%                     SEout=conn_init_vol(SEout);
                    emptycondition=~nnz(~isnan(wx)&wx~=0);
                    emptycondition_roi=emptycondition | all(iX*X'==0,2);
                    for nroi=1:length(idxredo),
                        touched(ncondition,idxredo(nroi))=true;
                        filename=fullfile(filepathresults,['BETA_Subject',num2str(nsub,CONN_x.opt.fmt1),'_Condition',num2str(icondition(ncondition),'%03d'),'_Source',num2str(iroi(idxredo(nroi)),'%03d'),'.nii']);
                        [filename,cache(nroi)]=conn_tempcache(filename);
                        Yout{nroi}=struct('mat',Y.matdim.mat,'dim',Y.matdim.dim,'fname',filename,'pinfo',[1;0;0],'n',[1,1],'dt',[spm_type('float32') spm_platform('bigend')],'descrip','');
                        if emptycondition_roi(1+nroi), Yout{nroi}.descrip='CONNlabel:MissingData'; end
                        conn_fileutils('deletefile',Yout{nroi}.fname);
                        Yout{nroi}=spm_create_vol(Yout{nroi});
                        redone_files=redone_files+1;
%                         filename=fullfile(filepathresults,['resultsDATA_Subject',num2str(nsub,'%03d'),'_Condition',num2str(ncondition,'%03d'),'_Source',num2str(iroi(idxredo(nroi)),'%03d'),'.mat']);
%                         Yout{nroi}=Y; Yout{nroi}.fname=filename;
%                         Yout{nroi}.size.Nt=1;%CONN_x.Setup.nsubjects;
%                         Yout{nroi}=conn_init_vol(Yout{nroi});
                    end
                    if ~emptycondition
                        for slice=1:Y.matdim.dim(3),
                            [y,idx]=conn_get_slice(Y,slice);
                            wx=ones(size(y,1),1);
                            switch(CONN_x.Analyses(ianalysis).weight),
                                case 1, wx=double(X1.conditionweights{1}>0);
                                case 2, wx=X1.conditionweights{1};
                                case 3, wx=X1.conditionweights{2};
                                case 4, wx=X1.conditionweights{3};
                            end
                            if ~(ischar(CONN_x.Analyses(ianalysis).modulation)||CONN_x.Analyses(ianalysis).modulation>0),
                                %if any(wx<0), conn_disp('Warning: temporal weights have negative values; cropping at 0'); end
                                wx=max(0,wx);
                                y=conn_wdemean(y,wx);
                                y=y.*repmat(wx,[1,size(y,2)]);
                            else
                                y=detrend(y,'constant');
                            end
                            B=full(real(iX*(X'*y)));
                            e=sqrt(sum(abs(y).^2,1));
                            switch(CONN_x.Analyses(ianalysis).measure),
                                case {1,2}, %correlation
                                    B=B./max(eps,r*e);
                                    B=atanh(max(eps-1,min(1-eps,B)));
                                    e(:)=1./sqrt(max(eps,DOF-3));
                                case {3,4,5}, %regression
                                    e=e/max(eps,DOF);
                            end
                            %conn_disp([slice,nsub])
                            %                         conn_write_slice(SEout,e,slice);%nsub);
                            %                         for nroi=1:length(idxredo);
                            %                             conn_write_slice(Yout{nroi},B(1+nroi,:),slice);%nsub);
                            %                         end
                            t=zeros(Y.matdim.dim(1:2));
                            t(idx)=e;
                            SEout=spm_write_plane(SEout,t,slice);
                            for nroi=1:length(idxredo);
                                t(idx)=B(1+nroi,:);
                                Yout{nroi}=spm_write_plane(Yout{nroi},t,slice);
                            end
                            n=n+1;
                            conn_waitbar(n/N,h,sprintf('Subject %d Condition %d',nsub,ncondition));
                        end
                    else
                        for slice=1:Y.matdim.dim(3),
                            t=zeros(Y.matdim.dim(1:2));
                            for nroi=1:length(idxredo)
                                Yout{nroi}=spm_write_plane(Yout{nroi},t,slice);
                            end
                            n=n+1;
                            conn_waitbar(n/N,h,sprintf('Subject %d Condition %d',nsub,ncondition));
                        end
                    end
                    for nroi=1:length(idxredo)
                        conn_tempcache(cache(nroi),'nii');
                    end
                    if CONN_x.Analyses(ianalysis).measure<3&&isfield(CONN_x.Setup,'outputfiles')&&numel(CONN_x.Setup.outputfiles)>=3&&any(CONN_x.Setup.outputfiles(3:min(5,numel(CONN_x.Setup.outputfiles)))),
                        for nroi=1:length(idxredo),
                            filename=fullfile(filepathresults,['BETA_Subject',num2str(nsub,CONN_x.opt.fmt1),'_Condition',num2str(icondition(ncondition),'%03d'),'_Source',num2str(iroi(idxredo(nroi)),'%03d'),'.nii']);
                            Va=spm_vol(filename);
                            t=spm_read_vols(Va);
                            t(~t)=nan;
                            if numel(CONN_x.Setup.outputfiles)>=3&&CONN_x.Setup.outputfiles(3),
                                filename=fullfile(filepathresults,['corr_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'_Source',num2str(iroi(idxredo(nroi)),'%03d'),'.nii']);
                                V=struct('mat',Y.matdim.mat,'dim',Y.matdim.dim,'fname',filename,'pinfo',[1;0;0],'n',[1,1],'dt',[spm_type('float32') spm_platform('bigend')],'descrip','');
                                if emptycondition_roi(1+nroi), V.descrip='CONNlabel:MissingData'; end
                                spm_write_vol(V,tanh(t));
                            end
                            if numel(CONN_x.Setup.outputfiles)>=4&&any(CONN_x.Setup.outputfiles(4:min(5,numel(CONN_x.Setup.outputfiles))))
                                p=spm_Ncdf(t*sqrt(max(0,DOF-3)));
                                p=2*min(p,1-p);
                            end
                            if numel(CONN_x.Setup.outputfiles)>=4&&CONN_x.Setup.outputfiles(4),
                                filename=fullfile(filepathresults,['p_corr_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'_Source',num2str(iroi(idxredo(nroi)),'%03d'),'.nii']);
                                V=struct('mat',Y.matdim.mat,'dim',Y.matdim.dim,'fname',filename,'pinfo',[1;0;0],'n',[1,1],'dt',[spm_type('float32') spm_platform('bigend')]);
                                spm_write_vol(V,p);
                            end
                            if numel(CONN_x.Setup.outputfiles)>=5&&CONN_x.Setup.outputfiles(5),
                                filename=fullfile(filepathresults,['pFDR_corr_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'_Source',num2str(iroi(idxredo(nroi)),'%03d'),'.nii']);
                                V=struct('mat',Y.matdim.mat,'dim',Y.matdim.dim,'fname',filename,'pinfo',[1;0;0],'n',[1,1],'dt',[spm_type('float32') spm_platform('bigend')]);
                                p(:)=conn_fdr(p(:));
                                spm_write_vol(V,p);
                            end
                        end
                    end
                else
                    n=n+Y.matdim.dim(3);
                    conn_waitbar(n/N,h,sprintf('Subject %d Condition %d',nsub,ncondition));
                end
            end
            for ncondition=secondaryconditions
                evaluatefunction=CONN_x.Setup.conditions.model{ncondition}{1};
                if isequal(evaluatefunction,'lin'), primaryconditionsnames=CONN_x.Setup.conditions.model{ncondition}(3:end);
                else primaryconditionsnames=CONN_x.Setup.conditions.model{ncondition}(2:end);
                end
                if isempty(primaryconditionsnames), error('Unable to compute secondary condition %s. Undefined input conditions'); end
                primaryconditions=find(ismember(CONN_x.Setup.conditions.names(1:end-1),primaryconditionsnames));
                if numel(primaryconditions)~=numel(primaryconditionsnames), error('Unable to compute secondary condition %s. Could not find a match for primary condition in %s',CONN_x.Setup.conditions.names{ncondition},sprintf('%s ',primaryconditionsnames{:})); end
                if isempty(REDO)&&~all(any(touched(primaryconditions,:),1)),
                    if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'overwrite'), REDO=CONN_x.gui.overwrite; else, REDO=conn_questdlg('Overwrite existing subject results?','','Yes', 'No', 'Yes');end;
                end
                if strcmp(lower(REDO),'no'),
                    idxredo=find(any(touched(primaryconditions,:),1)|arrayfun(@(n)~conn_existfile(fullfile(filepathresults,['BETA_Subject',num2str(nsub,CONN_x.opt.fmt1),'_Condition',num2str(icondition(ncondition),'%03d'),'_Source',num2str(iroi(n),'%03d'),'.nii'])),1:nrois));
                else,
                    idxredo=1:nrois;
                end
                for nroi=1:length(idxredo),
                    fileout=fullfile(filepathresults,['BETA_Subject',num2str(nsub,CONN_x.opt.fmt1),'_Condition',num2str(icondition(ncondition),'%03d'),'_Source',num2str(iroi(idxredo(nroi)),'%03d'),'.nii']);
                    filein={};
                    for ncondition2=primaryconditions
                        filein{end+1}=fullfile(filepathresults,['BETA_Subject',num2str(nsub,CONN_x.opt.fmt1),'_Condition',num2str(icondition(ncondition2),'%03d'),'_Source',num2str(iroi(idxredo(nroi)),'%03d'),'.nii']);
                    end
                    [ok,evaluatefunction] = conn_imcalc(filein,fileout,evaluatefunction,CONN_x.Setup.conditions.model{ncondition}{2},cat(2,CONN_x.Setup.l2covariates.values{nsub}{:}),CONN_x.Setup.l2covariates.names(1:end-1));
                    redone_files=redone_files+ok;
                end
                n=n+Y.matdim.dim(3);
                conn_waitbar(n/N,h,sprintf('Subject %d Condition %d',nsub,ncondition));
            end
        end
        if ~isempty(names), CONN_x.Analyses(ianalysis).sources=names; end
        try
            fileout=fullfile(filepathresults,'_list_sources.txt');
            fh=fopen(fileout,'wt');
            for n1=1:length(CONN_x.Analyses(ianalysis).sourcenames),fprintf(fh,'Source%03d = %s\n',n1,CONN_x.Analyses(ianalysis).sourcenames{n1});end
            fclose(fh);
        end
        try
            fileout=fullfile(filepathresults,'_list_conditions.txt');
            fh=fopen(fileout,'wt');
            for ncondition=1:nconditions,
                fprintf(fh,'Condition%03d = %s\n',icondition(ncondition),CONN_x.Setup.conditions.names{ncondition});
            end
            fclose(fh);
        end
        conn_disp('fprintf','      processed %d files\n',redone_files);
    end
    CONN_x.Analysis=analysisbak;
    conn_waitbar('close',h);
    CONN_x.isready(4)=1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates resultsROI_Subject###_Condition###.mat files (first-level ROI-to-ROI analyses)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(options==11) && any(CONN_x.Setup.steps([1])) && ~(isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'steps')&&~any(CONN_x.gui.steps([1]))),
    [path,name,ext]=fileparts(CONN_x.filename);
    filepath=CONN_x.folders.preprocessing;
    if nargin>1&&~isempty(varargin{1}),analyses=varargin{1}; % selected analysis only
    else, analyses=1:length(CONN_x.Analyses); 
    end
    if ischar(analyses)||iscell(analyses), analyses=find(ismember({CONN_x.Analyses.name},analyses)); end
    analyses(analyses<=0|analyses>numel(CONN_x.Analyses))=[];
    doanalyses=false(1,numel(analyses));
    for nanalyses=1:numel(analyses),if analyses(nanalyses)>0&&any(CONN_x.Analyses(analyses(nanalyses)).type==[1,3]), doanalyses(nanalyses)=true; end; end
    analyses=analyses(doanalyses);
    h=conn_waitbar(0,['Step ',num2str(sum(options<=11)),'/',num2str(length(options)),': ROI-to-ROI first-level analyses']);
    analysisbak=CONN_x.Analysis;
    if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'subjects'), validsubjects=CONN_x.gui.subjects; else validsubjects=1:CONN_x.Setup.nsubjects; end
    if isfield(CONN_x,'pobj')&&isstruct(CONN_x.pobj)&&isfield(CONN_x.pobj,'subjects'), validsubjects=CONN_x.pobj.subjects; end
    nconditions=length(CONN_x.Setup.conditions.names)-1;
    if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'conditions')&&~isempty(CONN_x.gui.conditions), validconditions=CONN_x.gui.conditions; else validconditions=1:length(CONN_x.Setup.conditions.names)-1; end
    icondition=[];isnewcondition=[];for ncondition=1:nconditions,[icondition(ncondition),isnewcondition(ncondition)]=conn_conditionnames(CONN_x.Setup.conditions.names{ncondition}); end
    secondaryconditions=validconditions(cellfun('length',CONN_x.Setup.conditions.model(validconditions))>0);
    validconditions=validconditions(cellfun('length',CONN_x.Setup.conditions.model(validconditions))==0); 
    if any(isnewcondition(validconditions)), error(['Some conditions have not been processed yet. Re-run previous step']); end
    N=numel(validsubjects)*(numel(validconditions)+numel(secondaryconditions))*100*length(analyses); nrois2bak=100;
    DOREDUCED=true; % (compute RRC square matrices only) set to false for back-compatibility
    for nanalyses=1:length(analyses),
        ianalysis=analyses(nanalyses);
        CONN_x.Analysis=ianalysis;
        names={};
%         ianalysis=CONN_x.Analysis;
        filepathresults=fullfile(CONN_x.folders.firstlevel,CONN_x.Analyses(ianalysis).name);
        conn_disp('fprintf','      first-level data in %s\n',filepathresults);
        REDO=[];%filename=fullfile(filepathresults,['resultsROI_Subject',num2str(validsubjects(1),'%03d'),'_Condition',num2str(icondition(validconditions(1)),'%03d'),'.mat']);
        %if ~isempty(dir(filename)),if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'overwrite'), REDO=CONN_x.gui.overwrite; else, REDO=conn_questdlg('Overwrite existing subject results?','','Yes', 'No', 'Yes');end; end
        if nanalyses==1, n=0; end
        Ni=[];
        redone_files=0;
        
        ConditionWeights={};
        if ~ischar(CONN_x.Analyses(ianalysis).modulation)&&CONN_x.Analyses(ianalysis).modulation>0
            if isempty(CONN_x.Analyses(ianalysis).conditions), validconditions2=validconditions;
            else validconditions2=find(ismember(CONN_x.Setup.conditions.names(1:end-1),CONN_x.Analyses(ianalysis).conditions));
            end
            for ncondition=union(validconditions,validconditions2),
                for nsub=validsubjects,
                    filename=fullfile(filepath,['ROI_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'.mat']);
                    if ~conn_existfile(filename), conn_disp(['Not ready to process step conn_process_10']); conn_waitbar('close',h);return; end
                    X1=load(filename,'conditionweights');
                    for n1=1:numel(X1.conditionweights)
                        ConditionWeights{nsub,n1}(:,ncondition)=X1.conditionweights{n1};
                    end
                end
            end
        end
        for nsub=validsubjects,
            touched=false(length(CONN_x.Setup.conditions.names)-1,1);
            maxrt=[];
            for ncondition=validconditions,
                filename=fullfile(filepathresults,['resultsROI_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'.mat']);
                if isempty(REDO)&&conn_existfile(filename),
                    if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'overwrite'), REDO=CONN_x.gui.overwrite; else, REDO=conn_questdlg('Overwrite existing subject results?','','Yes', 'No', 'Yes');end;
                end
                tfilename=fullfile(filepath,['ROI_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'.mat']);
                if isempty(dir(tfilename)), conn_disp(['Not ready to process step conn_process_10']); conn_waitbar('close',h);return; end
                X1=load(tfilename);
                [X,nill,names,xyz]=conn_designmatrix(CONN_x.Analyses(ianalysis).regressors,X1,[]);
                nrois=size(X,2)-1;
                if strcmp(lower(REDO),'yes')||~conn_existfile(filename), redo=true;
                else namestomatch=load(filename,'names'); redo=~isequal(namestomatch.names,names); 
                end
                if redo,
                    touched(ncondition,1)=true;
                    if DOREDUCED&CONN_x.Analyses(ianalysis).type==1
                        X2=X; names2=names; xyz2=xyz;
                    else
                        [X2,nill,names2,xyz2]=conn_designmatrix({CONN_x.Analysis_variables,CONN_x.Analyses(ianalysis).regressors},X1,[]);
                    end
                    nrois2=size(X2,2)-1;
                    [nill,idxroi1roi2]=ismember(names,names2);
                    %idxroi1roi2=zeros(1,nrois);
                    %for n1=1:nrois,temp=strmatch(names{n1},names2,'exact'); idxroi1roi2(n1)=temp(1);end
                    X2=cat(2,X2(:,1),X2(:,1+idxroi1roi2),X2(:,1+setdiff(1:nrois2,idxroi1roi2)));
                    names2=cat(2,{names2{idxroi1roi2}},{names2{setdiff(1:nrois2,idxroi1roi2)}});
                    xyz=cat(2,{xyz2{idxroi1roi2}},{xyz2{setdiff(1:nrois2,idxroi1roi2)}});
%                     for n1=1:nrois,temp=strmatch(names{n1},names2,'exact'); if isempty(temp), idxroi1roi2(n1)=nrois2+n1; else, idxroi1roi2(n1)=temp(1);end; end
%                     X2t=[X2(:,2:end) X(:,2:end)]; names2t=[names2 names]; xyz2t=[xyz2 xyz];
%                     X2=cat(2,X2(:,1),X2t(:,idxroi1roi2),X2(:,1+setdiff(1:nrois2,idxroi1roi2)));
%                     names2=cat(2,{names2t{idxroi1roi2}},{names2{setdiff(1:nrois2,idxroi1roi2)}});
%                     xyz=cat(2,{xyz2t{idxroi1roi2}},{xyz2{setdiff(1:nrois2,idxroi1roi2)}});
                    
                    if isempty(Ni), N=N+(sum(validsubjects>=nsub)+numel(validsubjects)*((numel(validconditions)+numel(secondaryconditions))-1))*(nrois2-nrois2bak); Ni=nrois2; end
                    Z=zeros(nrois,nrois2);%+diag(nan+zeros(nrois,1));
                    SE=zeros(1,nrois2);
                    nX=size(X,2);
                    wx=ones(size(X,1),1);
                    switch(CONN_x.Analyses(ianalysis).weight),
                        case 1, wx=double(X1.conditionweights{1}>0);
                        case 2, wx=X1.conditionweights{1};
                        case 3, wx=X1.conditionweights{2};
                        case 4, wx=X1.conditionweights{3};
                    end
                    if ~(ischar(CONN_x.Analyses(ianalysis).modulation)||CONN_x.Analyses(ianalysis).modulation>0)
                        if any(wx<0), conn_disp('Warning: temporal weights have negative values; cropping at 0'); end
                        wx=max(0,wx);
                        X=cat(2,X(:,1),conn_wdemean(X(:,2:end),wx));
                        X=X.*repmat(wx,[1,size(X,2)]);
                        X2=cat(2,X2(:,1),conn_wdemean(X2(:,2:end),wx));
                        X2=X2.*repmat(wx,[1,size(X2,2)]);
                    else
                        if ~ischar(CONN_x.Analyses(ianalysis).modulation),
                            %wx=X1.conditionweights{3};   % PPI
                            wx=ConditionWeights{nsub,3}(:,[setdiff(validconditions2,ncondition) ncondition]); %gPPI
                        elseif ~isempty(regexp(CONN_x.Analyses(ianalysis).modulation,'^(.*\/|.*\\)?Dynamic factor \d+$')),
                            CONN_x.Analyses(ianalysis).modulation=regexprep(CONN_x.Analyses(ianalysis).modulation,'\\|\/',['\',filesep]);
                            [name1,name2]=fileparts(CONN_x.Analyses(ianalysis).modulation);
                            [ok,value1]=ismember(name1,{CONN_x.dynAnalyses.name}); if ~ok, if numel(CONN_x.dynAnalyses)==1, value1=1; else error('Analysis name %s not found',name1); end; end
                            filename=fullfile(CONN_x.folders.firstlevel_dyn,CONN_x.dynAnalyses(value1).name,['dyn_Subject',num2str(nsub,'%03d'),'.mat']);
                            xmod=load(filename);
                            [ok,idx]=ismember(name2,xmod.names);
                            if ok, wx=xmod.data(:,[setdiff(1:size(xmod.data,2),idx),idx]);
                            else error('Temporal factor not found');
                            end
                            if ~isempty(wx), wx=conn_bsxfun(@times,wx,X1.conditionweights{1}>0); end
                        else
                            idx=find(strcmp(CONN_x.Analyses(ianalysis).modulation,X1.names));
                            if numel(idx)==1, wx=X1.data{idx};
                            elseif isempty(idx), error('Covariate not found. Please re-run ''Denoising'' step');
                            else,
                                idx=find(cellfun(@(x)all(isnan(x)),X1.xyz));
                                idx=idx(strcmp(CONN_x.Analyses(ianalysis).modulation,X1.names(idx)));
                                if numel(idx)==1, wx=X1.data{idx};
                                else error('Covariate not found');
                                end
                            end
                            if ~isempty(wx), wx=conn_bsxfun(@times,wx,X1.conditionweights{1}>0); end
                        end
                        %if size(wx,2)>1, conn_disp('Warning: multivariate interaction term not supported. Summing interaction term across multiple components'); end
                        inter=wx;
                        X2=detrend(X2,'constant');
                        %X=[X(:,1) detrend([X(:,2:end) reshape(repmat(permute(inter,[1 3 2]),[1,size(X,2),1]),size(X,1),[]) reshape(conn_bsxfun(@times,X,permute(inter,[1 3 2])),size(X,1),[])],'constant')];
                        %X=[X(:,1) detrend([X(:,2:end) repmat(inter,[1,size(X,2)]) conn_bsxfun(@times,X,inter)],'constant')];
                    end
%                     X=cat(2,X(:,1),conn_wdemean(X(:,2:end),wx));
%                     X=X.*repmat(wx,[1,size(X,2)]);
%                     X2=cat(2,X2(:,1),conn_wdemean(X2(:,2:end),wx));
%                     X2=X2.*repmat(wx,[1,size(X2,2)]);                   

                    if isempty(maxrt), maxrt=max(conn_get_rt(nsub)); end
                    switch(CONN_x.Analyses(ianalysis).measure),
                        case {1,3}, %bivariate
                            if ischar(CONN_x.Analyses(ianalysis).modulation)||CONN_x.Analyses(ianalysis).modulation>0, DOF=max(0,size(X,1)*(min(1/(2*maxrt),CONN_x.Preproc.filter(2))-max(0,CONN_x.Preproc.filter(1)))/(1/(2*maxrt))-2-size(inter,2));
                            else DOF=max(0,size(X,1)*(min(1/(2*maxrt),CONN_x.Preproc.filter(2))-max(0,CONN_x.Preproc.filter(1)))/(1/(2*maxrt))-1);
                            end
                        case {2,4}, %partial
                            DOF=max(0,size(X,1)*(min(1/(2*maxrt),CONN_x.Preproc.filter(2))-max(0,CONN_x.Preproc.filter(1)))/(1/(2*maxrt))-rank(X)+1);
                    end
                    emptycondition=~nnz(~isnan(wx)&wx~=0);
                    if emptycondition, 
                        Z(:)=nan; 
                        n=n+nrois2;
                        conn_waitbar(n/N,h,sprintf('Subject %d Condition %d',nsub,ncondition));
                    else
                        if ismember(CONN_x.Analyses(ianalysis).measure,[1 3])&&(ischar(CONN_x.Analyses(ianalysis).modulation)||CONN_x.Analyses(ianalysis).modulation>0), % speed-up for bivariate & parametric modulation
                            clear iXpers Xpers;
                            idxrois=[1 1+(1:nrois)];
                            x=[X(:,1) detrend([X(:,idxrois(2:end)) reshape(repmat(permute(inter,[1 3 2]),[1,numel(idxrois),1]),size(X,1),[]) reshape(conn_bsxfun(@times,X(:,idxrois),permute(inter,[1 3 2])),size(X,1),[])],'constant')];
                            nX=numel(idxrois);
                            nVars=size(x,2)/nX;
                            Xtemp=permute(reshape(x,[size(x,1),nX,nVars]),[1,3,2]);
                            iXpers=sparse(nVars*nX,nVars*nX);
                            for nXtemp=1:nX, iXpers(nXtemp:nX:end,nXtemp:nX:end)=pinv(Xtemp(:,:,nXtemp)'*Xtemp(:,:,nXtemp)); end
                            Xpers=x;
                        end
                        for nroi=1:nrois2,
                            y=X2(:,1+nroi);
                            if ischar(CONN_x.Analyses(ianalysis).modulation)||CONN_x.Analyses(ianalysis).modulation>0 % parametric modulation
                                idxrois=[1 1+setdiff(1:nrois,nroi)];
                                nX=numel(idxrois);
                                switch(CONN_x.Analyses(ianalysis).measure),
                                    case {1,3}, %bivariate
                                        idxkeep=conn_bsxfun(@plus,idxrois(:),(1+nrois)*(0:size(Xpers,2)/(1+nrois)-1));
                                        x=Xpers(:,idxkeep);
                                        iX=iXpers(idxkeep,idxkeep);
                                        %iX=sparse(nVars*nX,nVars*nX);
                                        %for nXtemp=1:size(Xtemp,3), iXtemp=pinv(Xtemp(:,:,nXtemp)'*Xtemp(:,:,nXtemp)); iX(nXtemp:nX:end,nXtemp:nX:end)=iXtemp; end
                                        %iX=pinv((x'*x).*kron(ones(3),eye(numel(idxrois))));
                                    case {2,4}, %partial
                                        x=[X(:,1) detrend([X(:,idxrois(2:end)) reshape(repmat(permute(inter,[1 3 2]),[1,numel(idxrois),1]),size(X,1),[]) reshape(conn_bsxfun(@times,X(:,idxrois),permute(inter,[1 3 2])),size(X,1),[])],'constant')];
                                        iX=pinv(x'*x);
                                end
                                nVars=size(x,2)/nX;
                                r=sqrt(diag(iX));
                                iX=iX((nVars-1)*nX+1:end,:);
                                r=r((nVars-1)*nX+1:end);
                            else % standard functional connectivity
                                x=cat(2,X(:,1),X(:,1+setdiff(1:nrois,nroi)));
                                switch(CONN_x.Analyses(ianalysis).measure),
                                    case {1,3}, %bivariate
                                        iX=diag(1./max(eps,sum(x.^2,1)));
                                        %iX=pinv(diag(diag(x'*x)));
                                    case {2,4}, %partial
                                        iX=pinv(x'*x);
                                end
                                r=sqrt(diag(iX));
                            end
                            B=full(real(iX*(x'*y)));
                            e=sqrt(sum(abs(y).^2,1));
                            switch(CONN_x.Analyses(ianalysis).measure),
                                case {1,2}, %correlation
                                    B=B./max(eps,r*e);
                                    B(~isnan(B))=atanh(max(eps-1,min(1-eps,B(~isnan(B)))));
                                    SE(nroi)=1./sqrt(max(eps,DOF-3));
                                case {3,4,5},
                                    SE(nroi)=e/max(eps,DOF);
                            end
                            Z(setdiff(1:nrois,nroi),nroi)=B(2:end);
                            if nroi<=nrois,Z(nroi,nroi)=nan;end
                            n=n+1;
                        end
                        conn_waitbar(n/N,h,sprintf('Subject %d Condition %d',nsub,ncondition));
                    end

%                     xyz={}; %note: this assumes constant number of dimensions per subject for analysis regressors
%                     for n1=1:length(CONN_x.Analysis_variables.names),
%                         for n2=1:CONN_x.Analysis_variables.deriv{n1}+1,
%                             for n3=1:CONN_x.Analysis_variables.dimensions{n1}(1),
%                                 idx=strmatch(CONN_x.Analysis_variables.names{n1},X1.names,'exact');
%                                 if isempty(idx), xyz{end+1}=''; else, xyz{end+1}=X1.xyz{idx}; end
%                             end
%                         end;
%                     end
%                     xyz=cat(2,{xyz{idxroi1roi2}},{xyz{setdiff(1:nrois2,idxroi1roi2)}});
                    regressors=CONN_x.Analyses(ianalysis).regressors;
                    filename=fullfile(filepathresults,['resultsROI_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'.mat']);
                    save(filename,'Z','regressors','names','names2','xyz','SE','DOF');
                    redone_files=redone_files+1;
                else
                    if isempty(Ni), n=n+nrois2bak; 
                    else n=n+Ni;
                    end
                    conn_waitbar(n/N,h,sprintf('Subject %d Condition %d',nsub,ncondition));
                end
            end
            for ncondition=secondaryconditions
                evaluatefunction=CONN_x.Setup.conditions.model{ncondition}{1};
                if isequal(evaluatefunction,'lin'), primaryconditionsnames=CONN_x.Setup.conditions.model{ncondition}(3:end);
                else primaryconditionsnames=CONN_x.Setup.conditions.model{ncondition}(2:end);
                end
                if isempty(primaryconditionsnames), error('Unable to compute secondary condition %s. Undefined input conditions'); end
                primaryconditions=find(ismember(CONN_x.Setup.conditions.names(1:end-1),primaryconditionsnames));
                if numel(primaryconditions)~=numel(primaryconditionsnames), error('Unable to compute secondary condition %s. Could not find a match for primary condition in %s',CONN_x.Setup.conditions.names{ncondition},sprintf('%s ',primaryconditionsnames{:})); end
                fileout=fullfile(filepathresults,['resultsROI_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'.mat']);
                if isempty(REDO)&&~any(touched(primaryconditions,:),1)&&conn_existfile(fileout),
                    if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'overwrite'), REDO=CONN_x.gui.overwrite; else, REDO=conn_questdlg('Overwrite existing subject results?','','Yes', 'No', 'Yes');end;
                end
                if strcmp(lower(REDO),'yes')||any(touched(primaryconditions,:),1)||~conn_existfile(fileout), redo=true;
                else redo=false;
                end
                if redo
                    filein={};
                    for ncondition2=primaryconditions
                        filein{end+1}=fullfile(filepathresults,['resultsROI_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition2),'%03d'),'.mat']);
                    end
                    [ok,evaluatefunction] = conn_imcalc(filein,fileout,evaluatefunction,CONN_x.Setup.conditions.model{ncondition}{2},cat(2,CONN_x.Setup.l2covariates.values{nsub}{:}),CONN_x.Setup.l2covariates.names(1:end-1));
                    redone_files=redone_files+ok;
                end
                if isempty(Ni), n=n+nrois2bak;
                else n=n+Ni;
                end
                conn_waitbar(n/N,h,sprintf('Subject %d Condition %d',nsub,ncondition));
            end            
        end
        if ~isempty(names), CONN_x.Analyses(ianalysis).sources=names; end
        try
            fileout=fullfile(filepathresults,'_list_conditions.txt');
            fh=fopen(fileout,'wt');
            for ncondition=1:nconditions,
                fprintf(fh,'Condition%03d = %s\n',icondition(ncondition),CONN_x.Setup.conditions.names{ncondition});
            end
            fclose(fh);
        end
        conn_disp('fprintf','      processed %d files\n',redone_files);
    end
    CONN_x.Analysis=analysisbak;
    conn_waitbar('close',h);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates BETA_Subject###_Condition###_Measure###.nii files (first-level voxel-to-voxel analyses)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(options==13|options==13.1) && any(CONN_x.Setup.steps([3])) && ~(isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'steps')&&~any(CONN_x.gui.steps([3]))),
    warning('off','MATLAB:DELETE:FileNotFound');
    [path,name,ext]=fileparts(CONN_x.filename);
    filepath=CONN_x.folders.preprocessing;
    if nargin>1&&~isempty(varargin{1}),analyses=varargin{1}; % selected analysis only
    else, analyses=1:length(CONN_x.vvAnalyses); 
    end
    if ischar(analyses)||iscell(analyses), analyses=find(ismember({CONN_x.vvAnalyses.name},analyses)); end
    analyses(analyses<=0|analyses>numel(CONN_x.vvAnalyses))=[];
    analysisbak=CONN_x.vvAnalysis;
    if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'subjects'), validsubjects=CONN_x.gui.subjects; else validsubjects=1:CONN_x.Setup.nsubjects; end
    if isfield(CONN_x,'pobj')&&isstruct(CONN_x.pobj)&&isfield(CONN_x.pobj,'subjects'), validsubjects=CONN_x.pobj.subjects; end
    nconditions=length(CONN_x.Setup.conditions.names)-1;
    if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'conditions')&&~isempty(CONN_x.gui.conditions), validconditions=CONN_x.gui.conditions; else validconditions=1:length(CONN_x.Setup.conditions.names)-1; end
    icondition=[];isnewcondition=[];for ncondition=1:nconditions,[icondition(ncondition),isnewcondition(ncondition)]=conn_conditionnames(CONN_x.Setup.conditions.names{ncondition}); end
    secondaryconditions=validconditions(cellfun('length',CONN_x.Setup.conditions.model(validconditions))>0);
    validconditions=validconditions(cellfun('length',CONN_x.Setup.conditions.model(validconditions))==0); 
    if ~isempty(validconditions), referenceconditions=validconditions;
    else
        evaluatefunction=CONN_x.Setup.conditions.model{secondaryconditions(1)}{1};
        if isequal(evaluatefunction,'lin'), primaryconditionsnames=CONN_x.Setup.conditions.model{secondaryconditions(1)}(3:end);
        else primaryconditionsnames=CONN_x.Setup.conditions.model{secondaryconditions(1)}(2:end);
        end
        referenceconditions=find(ismember(CONN_x.Setup.conditions.names(1:end-1),primaryconditionsnames));
    end
    if any(isnewcondition(referenceconditions)), error(['Some conditions have not been processed yet. Please re-run previous step']); end
    missingdata=arrayfun(@(n)isempty(dir(fullfile(filepath,['vvPC_Subject',num2str(validsubjects(1),'%03d'),'_Condition',num2str(icondition(n),'%03d'),'.mat']))),validconditions);
    if any(missingdata), conn_disp('fprintf','Warning: some conditions have not been processed yet (%s), skipping these conditions. Please re-run previous step to avoid this warning\n',sprintf('%s ',CONN_x.Setup.conditions.names{validconditions(missingdata)})); validconditions=validconditions(~missingdata); end
    missingdata=arrayfun(@(n)isempty(dir(fullfile(filepath,['vvPC_Subject',num2str(n,'%03d'),'_Condition',num2str(icondition(referenceconditions(1)),'%03d'),'.mat']))),validsubjects);
    if any(missingdata), conn_disp(['Some subjects have not been processed yet. Please re-run previous step']); return; end
    filename_B1=fullfile(filepath,['vvPC_Subject',num2str(validsubjects(1),'%03d'),'_Condition',num2str(icondition(referenceconditions(1)),'%03d'),'.mat']);
    Y1=conn_vol(filename_B1);
    if isfield(Y1,'issurface')&&Y1.issurface, issurface=true; else issurface=false; end
    LEFT2RIGHT=[];RIGHT2LEFT=[];
    REPLACERESULTS=true; 
    dogrouplevel=any(options==13|options==13.1); if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'grouplevel')&&~ismember(CONN_x.gui.grouplevel,[0 1]),dogrouplevel=false; end 
    dosubjectlevel=any(options==13|options==13.2); if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'grouplevel')&&~ismember(CONN_x.gui.grouplevel,[0 2]),dosubjectlevel=false; end 
    
    REDO=[]; 
    for nanalyses=1:length(analyses),
        redone_files=0;
        ianalysis=analyses(nanalyses);
        CONN_x.vvAnalysis=ianalysis;
        filepathresults=fullfile(CONN_x.folders.firstlevel_vv,CONN_x.vvAnalyses(ianalysis).name);
        
        idx=find(cellfun(@(a,b)(isequal(a,1)|isequal(a,2))&isequal(b,3),CONN_x.vvAnalyses(ianalysis).regressors.deriv,CONN_x.vvAnalyses(ianalysis).regressors.dimensions_out));
        if ~isempty(idx), CONN_x.vvAnalyses(ianalysis).regressors.dimensions_out(idx)=repmat({2+~issurface},1,numel(idx)); end
        measures=CONN_x.vvAnalyses(ianalysis).regressors;
        [nill,idx]=sort(cat(1,measures.localsupport{:}));mfieldnames=fieldnames(measures);for nfieldnames=1:numel(mfieldnames), measures.(mfieldnames{nfieldnames})=measures.(mfieldnames{nfieldnames})(idx); end % make analyses with same kernel size consecutive for speed improvemnet (using conn_v2v local storage)
        measuretypes=cat(1,measures.measuretype{:});
        measuretypes_1=ismember(measuretypes,[1,5]);     %v2v
        measuretypes_2=ismember(measuretypes,[2,3,4]);   %groupanalyses
        measuretypes_3=ismember(measuretypes,[6,7,8]);   %other
        measuretypes=1*measuretypes_1+2*measuretypes_2+3*measuretypes_3;
        [measuretypes,idx]=sort(measuretypes);mfieldnames=fieldnames(measures);for nfieldnames=1:numel(mfieldnames), measures.(mfieldnames{nfieldnames})=measures.(mfieldnames{nfieldnames})(idx); end
        nmeasures=numel(measures.names);
        nmeasures1=nnz(measuretypes_1);
        nmeasures2=nnz(measuretypes_2);
        nmeasures3=nnz(measuretypes_3);
        imeasure=[];isnew=[];for nmeasure=1:nmeasures,[imeasure(nmeasure),isnew(nmeasure)]=conn_v2v('match_measures',measures,nmeasure,'+');end
        isnew0=isnew;
        
        if nmeasures1>0&&any(options==13) % voxel-to-voxel measures
            N=0;n=0;
            h=conn_waitbar(0,['Step ',num2str(sum(options<=13)),'/',num2str(length(options)),': Voxel-to-voxel first-level analysis ',num2str(nanalyses),'/',num2str(length(analyses))]);
            conn_disp('fprintf','      first-level data in %s\n',filepathresults);
            for nsub=validsubjects,
                touched=false(length(CONN_x.Setup.conditions.names)-1,1);
                for ncondition=validconditions,
                    isnew=isnew0;
                    for nmeasure=1:nmeasures1,
                        filename=fullfile(filepathresults,['BETA_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'_Measure',num2str(imeasure(nmeasure),'%03d'),'_Component',num2str(1,'%03d'),'.nii']);
                        isnew(nmeasure)=isnew(nmeasure)|~conn_existfile(filename);
                    end
                    if isempty(REDO)&&~all(isnew(1:nmeasures1)),
                        if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'overwrite'), REDO=CONN_x.gui.overwrite; else, REDO=conn_questdlg('Overwrite existing subject results?','','Yes', 'No', 'Yes');end;
                    end
                    if strcmp(lower(REDO),'yes')||any(isnew(1:nmeasures1)),
                        touched(ncondition)=true;
                        %filename_B1=fullfile(filepath,['vvPC_Subject',num2str(nsub,'%03d'),'_Condition',num2str(ncondition,'%03d'),'.nii']);
                        %Y1=spm_vol(filename_B1);
                        filename_B1=fullfile(filepath,['vvPC_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'.mat']);
                        Y1=conn_vol(filename_B1);
                        filename_D1=fullfile(filepath,['vvPCeig_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'.mat']);
                        D1=load(filename_D1,'D');
                        %if ~N, N=(CONN_x.Setup.nsubjects-nsub+1)*nconditions*numel(Y1); numelY1bak=numel(Y1);
                        %else, N=N-numelY1bak+numel(Y1);end
                        if ~N, N=(numel(validsubjects)-sum(validsubjects<=nsub)+1)*(numel(validconditions)+numel(secondaryconditions))*Y1.size.Nt; numelY1bak=Y1.size.Nt;
                        else, N=N-numelY1bak+Y1.size.Nt;end
                        if isfield(Y1,'EmptyData')&&Y1.EmptyData, emptycondition=1; else emptycondition=0; end
                        
                        params=cell(1,nmeasures1);
                        for nmeasure=1:nmeasures1,
                            if strcmp(lower(REDO),'yes')||isnew(nmeasure)
                                if isfield(Y1,'issurface')&&Y1.issurface, issurface=true; else issurface=false; end
                                params{nmeasure}=conn_v2v('compute_start',measures,nmeasure,Y1.matdim.mat,issurface);
                                redone_files=redone_files+1;
                            end
                        end
                        %[gridx,gridy,gridz]=ndgrid(1:Y1(1).dim(1),1:Y1(1).dim(2),1:Y1(1).dim(3));xyz=[gridx(:),gridy(:),gridz(:),ones(numel(gridx),1)]';
                        for ndim=1:Y1.size.Nt,%numel(Y1),
                            doneread=0;
                            for nmeasure=1:nmeasures1,
                                if (strcmp(lower(REDO),'yes')||isnew(nmeasure))&&ndim<=params{nmeasure}.dimensions_in,
                                    if ~doneread
                                        [y1,idxy1]=conn_get_time(Y1,ndim);
                                        %y1=reshape(spm_get_data(Y1(ndim),xyz),Y1(1).dim); %y1=spm_read_vols(Y1(ndim));
                                        doneread=1;
                                    end
                                    params{nmeasure}=conn_v2v('compute_step',params{nmeasure},y1,D1.D(ndim),D1.D,numel(idxy1));
                                end
                            end
                            n=n+1;
                            conn_waitbar(n/N,h,sprintf('Subject %d Condition %d',nsub,ncondition));
                        end
                        for nmeasure=1:nmeasures1,
                            if strcmp(lower(REDO),'yes')||isnew(nmeasure)
                                d=conn_v2v('compute_end',params{nmeasure});
                                if iscell(d)
                                    dsum=0;
                                    for nout=numel(d):-1:1
                                        filename=fullfile(filepathresults,['BETA_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'_Measure',num2str(imeasure(nmeasure),'%03d'),'_Component',num2str(nout,'%03d'),'.nii']);
                                        Yout=struct('fname',filename,'mat',Y1.matdim.mat,'dim',Y1.matdim.dim,'n',[1,1],'pinfo',[1;0;0],'dt',[spm_type('float32'),spm_platform('bigend')],'descrip','');
                                        if emptycondition, Yout.descrip='CONNlabel:MissingData'; end
                                        %Yout=Y1(1);Yout.fname=filename;
                                        if isequal(d{nout},0), d{nout}=zeros(Yout.dim); end
                                        spm_write_vol(Yout,d{nout});
                                        dsum=dsum+abs(d{nout}).^2;
                                    end
                                    if numel(d)>1
                                        dsum=sqrt(dsum);
                                        filename=fullfile(filepathresults,['BETA_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'_Measure',num2str(imeasure(nmeasure),'%03d'),'_Component',num2str(0,'%03d'),'.nii']);
                                        Yout=struct('fname',filename,'mat',Y1.matdim.mat,'dim',Y1.matdim.dim,'n',[1,1],'pinfo',[1;0;0],'dt',[spm_type('float32'),spm_platform('bigend')],'descrip','');
                                        if emptycondition, Yout.descrip='CONNlabel:MissingData'; end
                                        %Yout=Y1(1);Yout.fname=filename;
                                        if isequal(dsum,0), dsum=zeros(Yout.dim); end
                                        spm_write_vol(Yout,dsum);
                                    end
                                else
                                    filename=fullfile(filepathresults,['BETA_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'_Measure',num2str(imeasure(nmeasure),'%03d'),'_Component',num2str(1,'%03d'),'.nii']);
                                    Yout=struct('fname',filename,'mat',Y1.matdim.mat,'dim',Y1.matdim.dim,'n',[1,1],'pinfo',[1;0;0],'dt',[spm_type('float32'),spm_platform('bigend')],'descrip','');
                                    if emptycondition, Yout.descrip='CONNlabel:MissingData'; end
                                    %Yout=Y1(1);Yout.fname=filename;
                                    if isequal(d,0), d=zeros(Yout.dim); end
                                    spm_write_vol(Yout,d);
                                end
                            end
                        end
                    end
                end

                for ncondition=secondaryconditions
                    evaluatefunction=CONN_x.Setup.conditions.model{ncondition}{1};
                    if isequal(evaluatefunction,'lin'), primaryconditionsnames=CONN_x.Setup.conditions.model{ncondition}(3:end);
                    else primaryconditionsnames=CONN_x.Setup.conditions.model{ncondition}(2:end);
                    end
                    if isempty(primaryconditionsnames), error('Unable to compute secondary condition %s. Undefined input conditions'); end
                    primaryconditions=find(ismember(CONN_x.Setup.conditions.names(1:end-1),primaryconditionsnames));
                    if numel(primaryconditions)~=numel(primaryconditionsnames), error('Unable to compute secondary condition %s. Could not find a match for primary condition in %s',CONN_x.Setup.conditions.names{ncondition},sprintf('%s ',primaryconditionsnames{:})); end
                    
                    fileout=fullfile(filepathresults,['BETA_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'_Measure',num2str(imeasure(nmeasure),'%03d'),'_Component',num2str(1,'%03d'),'.nii']);
                    if isempty(REDO)&&~any(touched(primaryconditions,:),1)&&conn_existfile(fileout),
                        if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'overwrite'), REDO=CONN_x.gui.overwrite; else, REDO=conn_questdlg('Overwrite existing subject results?','','Yes', 'No', 'Yes');end;
                    end
                    if strcmp(lower(REDO),'yes')||any(touched(primaryconditions,:),1)||~conn_existfile(fileout), redo=true;
                    else redo=false;
                    end
                    if redo
                        filein={};
                        for ncondition2=primaryconditions
                            filein{end+1}=fullfile(filepathresults,['BETA_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition2),'%03d'),'_Measure',num2str(imeasure(nmeasure),'%03d'),'_Component',num2str(1,'%03d'),'.nii']);
                        end
                        [ok,evaluatefunction] = conn_imcalc(filein,fileout,evaluatefunction,CONN_x.Setup.conditions.model{ncondition}{2},cat(2,CONN_x.Setup.l2covariates.values{nsub}{:}),CONN_x.Setup.l2covariates.names(1:end-1));
                        redone_files=redone_files+ok;
                    end
                    if ~N, N=(numel(validsubjects)-sum(validsubjects<=nsub)+1)*(numel(validconditions)+numel(secondaryconditions))*1; numelY1bak=1; end
                    n=n+numelY1bak;
                    conn_waitbar(n/N,h,sprintf('Subject %d Condition %d',nsub,ncondition));
                end
            end
            conn_waitbar('close',h);
        end
        
        if nmeasures3>0&&any(options==13) % other measures
            N=0;n=0;
            h=conn_waitbar(0,['Step ',num2str(sum(options<=13)),'/',num2str(length(options)),': Voxel-to-voxel first-level analysis ',num2str(nanalyses),'/',num2str(length(analyses))]);
            conn_disp('fprintf','      first-level data in %s\n',filepathresults);
            for nsub=validsubjects,
                touched=false(length(CONN_x.Setup.conditions.names)-1,1);
                for ncondition=validconditions,
                    isnew=isnew0;
                    for nmeasure=nmeasures1+nmeasures2+(1:nmeasures3),
                        filename=fullfile(filepathresults,['BETA_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'_Measure',num2str(imeasure(nmeasure),'%03d'),'_Component',num2str(1,'%03d'),'.nii']);
                        isnew(nmeasure)=isnew(nmeasure)|~conn_existfile(filename);
                    end
                    if isempty(REDO)&&~all(isnew(nmeasures1+nmeasures2+(1:nmeasures3))),
                        if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'overwrite'), REDO=CONN_x.gui.overwrite; else, REDO=conn_questdlg('Overwrite existing subject results?','','Yes', 'No', 'Yes');end;
                    end
                    if strcmp(lower(REDO),'yes')||any(isnew(nmeasures1+nmeasures2+(1:nmeasures3))),
                        touched(ncondition)=true;
                        filename_B1=fullfile(filepath,['NORMS_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'.mat']);
                        if ~conn_existfile(filename_B1), conn_msgbox('Norm information not available for all subjects. Please re-run Denoising step (s2v pipeline)','',2); return; end
                        Y1=conn_vol(filename_B1);
                        if ~N, N=(numel(validsubjects)-sum(validsubjects<=nsub)+1)*(numel(validconditions)+numel(secondaryconditions)); end
                        if isfield(Y1,'EmptyData')&&Y1.EmptyData, emptycondition=1; else emptycondition=0; end
                        if ~emptycondition, 
                            [y1,idxy1]=conn_get_time(Y1,1);
                            [y2,idxy1]=conn_get_time(Y1,2);
                        end
                        params=cell(1,nmeasures1);
                        for nmeasure=nmeasures1+nmeasures2+(1:nmeasures3),
                            if strcmp(lower(REDO),'yes')||isnew(nmeasure)
                                if isfield(Y1,'issurface')&&Y1.issurface, issurface=true; else issurface=false; end
                                params{nmeasure}=conn_v2v('compute_start',measures,nmeasure,Y1.matdim.mat,issurface);
                            end
                        end
                        
                        for nmeasure=nmeasures1+nmeasures2+(1:nmeasures3),
                            if strcmp(lower(REDO),'yes')||isnew(nmeasure)
                                filename=fullfile(filepathresults,['BETA_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'_Measure',num2str(imeasure(nmeasure),'%03d'),'_Component',num2str(1,'%03d'),'.nii']);
                                Yout=struct('fname',filename,'mat',Y1.matdim.mat,'dim',Y1.matdim.dim,'n',[1,1],'pinfo',[1;0;0],'dt',[spm_type('float32'),spm_platform('bigend')],'descrip','');
                                if emptycondition, 
                                    Yout.descrip='CONNlabel:MissingData'; 
                                    d=zeros(Y1.matdim.dim);
                                    params{nmeasure}.m=d;
                                else
                                    if measures.measuretype{nmeasure}==6 % ALFF
                                        params{nmeasure}.m=y2;
                                    elseif measures.measuretype{nmeasure}==7 % fALFF
                                        params{nmeasure}.m=y2./max(eps,y1);
                                    elseif measures.measuretype{nmeasure}==8 % IHC
                                        filename=fullfile(filepath,['DATA_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'.mat']);
                                        Y=conn_vol(filename);
                                        if isfield(Y,'issurface')&&Y.issurface,
                                            y3=zeros([numel(y2)/2,2]);
                                            if isempty(LEFT2RIGHT),load(fullfile(fileparts(which(mfilename)),'utils','surf','lhrh.mat'),'LEFT2RIGHT','RIGHT2LEFT'); end
                                            for nt=1:Y.size.Nt,
                                                if Y.conditionsweights{1}(nt)~=0
                                                    [ytemp,idx]=conn_get_time(Y,nt);
                                                    ytemp=reshape(ytemp,[],2);
                                                    y3=y3+Y.conditionsweights{1}(nt)*[ytemp(:,1).*ytemp(RIGHT2LEFT,2), ytemp(:,2).*ytemp(LEFT2RIGHT,1)];
                                                end
                                            end
                                            y3=y3/max(eps,sum(Y.conditionsweights{1}));
                                            ty2=reshape(y2,[],2);
                                            y3=reshape(y3./max(eps,ty2.*[ty2(RIGHT2LEFT,2), ty2(LEFT2RIGHT,1)]),size(y2));
                                        else
                                            y3=zeros(size(y1));
                                            for slice=1:Y.matdim.dim(3),
                                                [ytemp,idx]=conn_get_slice(Y,slice);
                                                idx_involume=idx+Y.matdim.dim(1)*Y.matdim.dim(2)*(slice-1);
                                                t=zeros(Y.matdim.dim(1:2));
                                                t(idx)=1:numel(idx);
                                                t=flipud(t);
                                                flipped=t(idx);
                                                flippedok=flipped>0 & flipped~=reshape(1:numel(idx),size(flipped));
                                                if isempty(LEFT2RIGHT)&&(nnz(Y.matdim.mat(1:3,1:3)-diag(diag(Y.matdim.mat(1:3,1:3))))>0||Y.matdim.mat(1,:)*[Y.matdim.dim(1)+1;0;0;2]~=0), LEFT2RIGHT=false; conn_disp('Warning: data first-dimension does not represent x spatial-coordinates; homotopic interhemispheric correlations may be inaccurate'); end
                                                rnorm_post=zeros(size(flipped));
                                                rnorm_post(flippedok)=(sum(repmat(Y.conditionsweights{1},1,nnz(flippedok)).*ytemp(:,flippedok).*ytemp(:,flipped(flippedok)),1)/max(eps,sum(Y.conditionsweights{1})))./reshape(max(eps,y2(idx_involume(flippedok)).*y2(idx_involume(flipped(flippedok)))),1,[]);
                                                y3(idx_involume)=rnorm_post;
                                            end
                                        end
                                        params{nmeasure}.m=atanh(max(-1,min(1,y3)));
                                    end
                                    d=conn_v2v('compute_end',params{nmeasure});
                                end
                                spm_write_vol(Yout,d);
                                redone_files=redone_files+1;
                            end
                        end                        
                        n=n+1;
                        conn_waitbar(n/N,h,sprintf('Subject %d Condition %d',nsub,ncondition));
                    end
                end
                for ncondition=secondaryconditions
                    evaluatefunction=CONN_x.Setup.conditions.model{ncondition}{1};
                    if isequal(evaluatefunction,'lin'), primaryconditionsnames=CONN_x.Setup.conditions.model{ncondition}(3:end);
                    else primaryconditionsnames=CONN_x.Setup.conditions.model{ncondition}(2:end);
                    end
                    if isempty(primaryconditionsnames), error('Unable to compute secondary condition %s. Undefined input conditions'); end
                    primaryconditions=find(ismember(CONN_x.Setup.conditions.names(1:end-1),primaryconditionsnames));
                    if numel(primaryconditions)~=numel(primaryconditionsnames), error('Unable to compute secondary condition %s. Could not find a match for primary condition in %s',CONN_x.Setup.conditions.names{ncondition},sprintf('%s ',primaryconditionsnames{:})); end
                    
                    fileout=fullfile(filepathresults,['BETA_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'_Measure',num2str(imeasure(nmeasure),'%03d'),'_Component',num2str(1,'%03d'),'.nii']);
                    if isempty(REDO)&&~any(touched(primaryconditions,:),1)&&conn_existfile(fileout),
                        if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'overwrite'), REDO=CONN_x.gui.overwrite; else, REDO=conn_questdlg('Overwrite existing subject results?','','Yes', 'No', 'Yes');end;
                    end
                    if strcmp(lower(REDO),'yes')||any(touched(primaryconditions,:),1)||~conn_existfile(fileout), redo=true;
                    else redo=false;
                    end
                    if redo
                        filein={};
                        for ncondition2=primaryconditions
                            filein{end+1}=fullfile(filepathresults,['BETA_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition2),'%03d'),'_Measure',num2str(imeasure(nmeasure),'%03d'),'_Component',num2str(1,'%03d'),'.nii']);
                        end
                        [ok,evaluatefunction] = conn_imcalc(filein,fileout,evaluatefunction,CONN_x.Setup.conditions.model{ncondition}{2},cat(2,CONN_x.Setup.l2covariates.values{nsub}{:}),CONN_x.Setup.l2covariates.names(1:end-1));
                        redone_files=redone_files+ok;
                    end
                    if ~N, N=(numel(validsubjects)-sum(validsubjects<=nsub)+1)*(numel(validconditions)+numel(secondaryconditions))*1; end
                    n=n+1;
                    conn_waitbar(n/N,h,sprintf('Subject %d Condition %d',nsub,ncondition));
                end
            end
            conn_waitbar('close',h);
        end
        
        if nmeasures2>0, %&& isequal(validsubjects,1:CONN_x.Setup.nsubjects),% MVPA,ICA,PCA
            conn_disp(['Step ',num2str(sum(options<=13)),'/',num2str(length(options)),': Voxel-to-voxel group-level analysis ',num2str(nanalyses),'/',num2str(length(analyses))]);
            conn_disp('fprintf','      first-level data in %s\n',filepathresults);
        end
        if nmeasures2>0 && ~isequal(validsubjects,1:CONN_x.Setup.nsubjects),
            if isfield(CONN_x,'pobj')&&isstruct(CONN_x.pobj)&&isfield(CONN_x.pobj,'subjects'),
                if isfield(CONN_x.pobj,'partition')&&isequal(CONN_x.pobj.partition,[1 1])
                    conn_disp('fprintf','NOTE: single-job, subset of %d subjects only\n',numel(validsubjects));
                elseif 0,%~dogrouplevel
                    conn_disp('NOTE: performing subject-level analyses only');
                else
                    conn_disp('WARNING: group-level factorization analyses parallelization not yet available (run locally or as a single job instead). Skipping these analyses');
                    nmeasures2=0;
                end
            end
        end
        if nmeasures2>0, %&& isequal(validsubjects,1:CONN_x.Setup.nsubjects),% MVPA,ICA,PCA
            isnew=isnew0;
            %         MAXDIMS=256;
            %         h=conn_waitbar(0,['Step ',num2str(sum(options<=13)),'/',num2str(length(options)),': connectome-level MVPA analyses']);
            N=(numel(validconditions)+numel(secondaryconditions))*nmeasures2;
            n=0;
            for ncondition=validconditions,
                for nmeasure=nmeasures1+(1:nmeasures2),
                    filename=fullfile(filepathresults,['BETA_Subject',num2str(CONN_x.Setup.nsubjects,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'_Measure',num2str(imeasure(nmeasure),'%03d'),'_Component',num2str(1,'%03d'),'.nii']);
                    isnew(nmeasure)=isnew(nmeasure)|~conn_existfile(filename);
                end
            end
            if isempty(REDO)&&~all(isnew(nmeasures1+(1:nmeasures2))),
                if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'overwrite'), REDO=CONN_x.gui.overwrite; else, REDO=conn_questdlg('Overwrite existing subject results?','','Yes', 'No', 'Yes');end;
            end
            if strcmp(lower(REDO),'yes')||any(isnew(nmeasures1+(1:nmeasures2))),
                %             filename_D1=fullfile(filepath,['vvPCcov_SubjectA',num2str(1,'%03d'),'_SubjectB',num2str(1,'%03d'),'_ConditionA',num2str(ncondition,'%03d'),'_ConditionB',num2str(ncondition,'%03d'),'.mat']);
                %             D1=load(filename_D1,'D');
                %                 filename_B1=fullfile(filepath,['vvPC_Subject',num2str(1,'%03d'),'_Condition',num2str(ncondition,'%03d'),'.nii']);
                %                 Y1=spm_vol(filename_B1);
                %                 sA=prod(Y1(1).dim);
                %                 [gridx,gridy,gridz]=ndgrid(1:Y1(1).dim(1),1:Y1(1).dim(2),1:Y1(1).dim(3));xyz=[gridx(:),gridy(:),gridz(:),ones(numel(gridx),1)]';
                filename_B1=fullfile(filepath,['vvPC_Subject',num2str(validsubjects(1),'%03d'),'_Condition',num2str(icondition(validconditions(1)),'%03d'),'.mat']);
                Y1=conn_vol(filename_B1);
                sA=prod(Y1.matdim.dim);
                [gridx,gridy,gridz]=ind2sub(Y1(1).matdim.dim,Y1(1).voxels);xyz=Y1(1).matdim.mat*[gridx(:),gridy(:),gridz(:),ones(numel(gridx),1)]'; clear gridx gridy gridz;
                
                Y1Nt=zeros(CONN_x.Setup.nsubjects,nconditions);
                Y1MD=ones(CONN_x.Setup.nsubjects,nconditions);
                Y1weights=cell(CONN_x.Setup.nsubjects,nconditions);
                for nsub=validsubjects
                    for ncondition=validconditions,
                        filename=fullfile(filepath,['vvPC_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'.mat']);
                        Yout=conn_vol(filename);
                        if ~isequal(Y1.matdim.dim,Yout.matdim.dim) || ~isequal(Y1.voxels,Yout.voxels), error('Unequal analysis voxels across subjects. Modify ''spatial resolution'' and/or ''analysis mask'' fields in Setup->Options to match analysis voxels across subjects'); end
                        Y1Nt(nsub,ncondition)=Yout.size.Nt;
                        Y1weights{nsub,ncondition}=Yout.conditionsweights;
                        if ~isfield(Yout,'EmptyData')||Yout.EmptyData==0, Y1MD(nsub,ncondition)=0; end
                    end
                end
                if isfield(CONN_x.vvAnalyses(ianalysis),'mask')&&iscell(CONN_x.vvAnalyses(ianalysis).mask)&&~isempty(CONN_x.vvAnalyses(ianalysis).mask{1})
                    VmaskV=spm_vol(CONN_x.vvAnalyses(ianalysis).mask{1});
                    Vmask=spm_get_data(VmaskV,pinv(VmaskV.mat)*xyz)>0;
                else Vmask=[];
                end
                DOGICA3=true;
                ICAMETHOD='';
                if isfield(CONN_x.vvAnalyses(ianalysis),'options')&&ischar(CONN_x.vvAnalyses(ianalysis).options)&&~isempty(CONN_x.vvAnalyses(ianalysis).options)
                    if ~isempty(regexpi(CONN_x.vvAnalyses(ianalysis).options,'gica1')), DOGICA3=false; end
                    if ~isempty(regexpi(CONN_x.vvAnalyses(ianalysis).options,'tanh|gauss|pow3')), ICAMETHOD=char(regexp(lower(CONN_x.vvAnalyses(ianalysis).options),'tanh|gauss|pow3','match','once')); end
                end
                
                NdimsIn=0; for nmeasure=nmeasures1+(1:nmeasures2), NdimsIn=max(NdimsIn,min(measures.dimensions_in{nmeasure},max(Y1Nt(:)))); end
                if 1,%dogrouplevel
                    h2=conn_waitbar(0,['computing subject x subject covariance (voxel-to-voxel analysis ',num2str(nanalyses),'/',num2str(length(analyses)),')']);
                    MAXMEM=2e9; % max allowed C matrix memory usage for direct storage (default 2Gb)
                    ismtxC=8*(numel(validsubjects)*numel(validconditions)*NdimsIn)^2>MAXMEM;
                    if ismtxC
                        filename=fullfile(filepathresults,['TEMPORAL1_Measure',num2str(imeasure(nmeasures1+1),'%03d'),'.mtx']);
                        C=conn_mtx('init',[NdimsIn*numel(validsubjects)*numel(validconditions) NdimsIn numel(validsubjects) numel(validconditions)],filename);
                        ismtxC=true;
                    else C=0;
                    end
                    for slice=1:Y1.size.Ns,
                        for isub=1:numel(validsubjects)
                            nsub=validsubjects(isub);
                            for ivalidcondition=1:numel(validconditions),
                                ncondition=validconditions(ivalidcondition);
                                Y1.fname=fullfile(filepath,['vvPC_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'.mat']);
                                Y1.size.Nt=Y1Nt(nsub,ncondition);
                                temp=conn_get_slice(Y1,slice);
                                if ~isempty(Vmask), temp=temp(:,Vmask(sum(Y1.size.Nv(1:slice-1))+(1:Y1.size.Nv(slice)))); end
                                if isub==1&&ivalidcondition==1, x=zeros([size(temp,2),NdimsIn,numel(validsubjects),numel(validconditions)]); end
                                if ~isempty(temp)
                                    x(:,1:min(size(temp,1),NdimsIn),isub,ivalidcondition)=temp(1:min(size(temp,1),NdimsIn),:)';
                                end
                            end
                        end
                        conn_waitbar(slice/Y1.size.Ns,h2,sprintf('Slice %d',slice));
                        x=x(:,:);
                        if ismtxC
                            blockcolumns=conn_mtx('getblockcolumns',C);
                            for nt=1:numel(blockcolumns), conn_mtx('addtoblock',C,nt,x'*x(:,blockcolumns{nt})); end
                        else C=C+x'*x;
                        end
                    end
                    conn_waitbar('close',h2);%close(h2);
                    if ~ismtxC, C=reshape(C,[NdimsIn,numel(validsubjects),numel(validconditions),NdimsIn,numel(validsubjects),numel(validconditions)]); end
                    filename=fullfile(filepathresults,['TEMPORAL1_Measure',num2str(imeasure(nmeasures1+1),'%03d'),'.mat']);
                    save(filename,'C');
                end
                
                for nmeasure=nmeasures1+(1:nmeasures2),
                    if strcmp(lower(REDO),'yes')||isnew(nmeasure)

                        thisNdimsIn=min(measures.dimensions_in{nmeasure},max(Y1Nt(:)));
                        NdimsOut=min(thisNdimsIn*numel(validsubjects)*numel(validconditions),measures.dimensions_out{nmeasure});
                        if 1,%dogrouplevel
                            filename=fullfile(filepathresults,['TEMPORAL1_Measure',num2str(imeasure(nmeasures1+1),'%03d'),'.mat']);
                            load(filename,'C');
                            if NdimsIn>thisNdimsIn, % placeholder (to do: break down MVPA into separable group- and subject- level processes)
                                if ismtxC,
                                    C=conn_mtx('zerocolumns',C,thisNdimsIn+1:NdimsIn);
                                    C=conn_mtx('zerorows',C,conn_bsxfun(@plus,NdimsIn*(0:numel(validsubjects)*numel(validconditions)-1)',thisNdimsIn+1:NdimsIn));
                                else
                                    C(thisNdimsIn+1:NdimsIn,:,:,:,:,:)=0;
                                    C(:,:,:,thisNdimsIn+1:NdimsIn,:,:)=0;
                                end
                            end
                        end
                        
                        %                     filename=fullfile(filepathresults,['TEMPORAL2_Measure',num2str(imeasure(nmeasure),'%03d'),'.mat']);
                        %                     Yout=Y1; Yout.fname=filename;Yout.size.Nt=CONN_x.Setup.nsubjects*nconditions*NdimsOut;
                        %                     Yout=conn_init_vol(Yout);
                        %                     filename=fullfile(filepathresults,['TEMPORAL3_Measure',num2str(imeasure(nmeasure),'%03d'),'.mat']);
                        %                     Dout=Y1; Dout.fname=filename;Dout.size.Nt=NdimsOut;
                        %                     Dout=conn_init_vol(Dout);
                        
                        if measures.measuretype{nmeasure}==2 % group-MVPA
                            if 1,%dogrouplevel,
                                h2=conn_waitbar(0,['computing MVPA components']);
                                if ~ismtxC, C=permute(C,[1,4,2,3,5,6]); end
                                clear filesout filesoutCov cache;
                                for nsub=1:CONN_x.Setup.nsubjects
                                    for ncondition=validconditions,
                                        for ndim=1:NdimsOut
                                            filename=fullfile(filepathresults,['BETA_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'_Measure',num2str(imeasure(nmeasure),'%03d'),'_Component',num2str(ndim,'%03d'),'.nii']);
                                            conn_fileutils('deletefile',filename); % delete all conditions to avoid mix-up of different models
                                        end
                                    end
                                    for ncondition=validconditions,
                                        for ndim=1:NdimsOut
                                            filename=fullfile(filepathresults,['BETA_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'_Measure',num2str(imeasure(nmeasure),'%03d'),'_Component',num2str(ndim,'%03d'),'.nii']);
                                            [filename, cache(nsub,ncondition,ndim)]=conn_tempcache(filename);
                                            temp=struct('fname',filename,'mat',Y1.matdim.mat,'dim',Y1.matdim.dim,'n',[1,1],'pinfo',[1;0;0],'dt',[spm_type('float32'),spm_platform('bigend')],'descrip','');
                                            if Y1MD(nsub,ncondition)||~ismember(nsub,validsubjects), temp.descrip='CONNlabel:MissingData'; end
                                            filesout(nsub,ncondition,ndim)=spm_create_vol(temp);
                                            redone_files=redone_files+1;
                                        end
                                    end
                                end
                                for ndim=1:NdimsOut
                                    filename=fullfile(filepathresults,['PCAcov_Measure',num2str(imeasure(nmeasure),'%03d'),'_Component',num2str(ndim,'%03d'),'.nii']);
                                    conn_fileutils('deletefile',filename);
                                    filesoutCov(ndim)=spm_create_vol(struct('fname',filename,'mat',Y1.matdim.mat,'dim',Y1.matdim.dim,'n',[1,1],'pinfo',[1;0;0],'dt',[spm_type('float32'),spm_platform('bigend')],'descrip',mfilename));
                                end
                                maxvoxels=max(1,floor(MAXMEM/(8*(numel(validsubjects)*numel(validconditions))^2)));
                                if maxvoxels>=NdimsOut, Qglobal=zeros([numel(validsubjects)*numel(validconditions)*[1 1], NdimsOut]); 
                                else
                                    filenameglobal=fullfile(filepathresults,'TEMPORAL2.mtx');
                                    Qglobal=conn_mtx('init',[numel(validsubjects)*numel(validconditions)*[1 1], NdimsOut],filenameglobal);
                                end
                                for slice=1:Y1.size.Ns,
                                    for isub=1:numel(validsubjects)
                                        nsub=validsubjects(isub);
                                        for ivalidcondition=1:numel(validconditions)
                                            ncondition=validconditions(ivalidcondition);
                                            Y1.fname=fullfile(filepath,['vvPC_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'.mat']);
                                            Y1.size.Nt=Y1Nt(nsub,ncondition);
                                            [temp,idx]=conn_get_slice(Y1,slice);
                                            if isub==1&&ivalidcondition==1, x=zeros([size(temp,2),NdimsIn,numel(validsubjects),numel(validconditions)]); end
                                            if ~isempty(temp)
                                                x(:,1:min(size(temp,1),NdimsIn),isub,ivalidcondition)=temp(1:min(size(temp,1),NdimsIn),:)';
                                            end
                                        end
                                    end
                                    xEig=zeros([numel(validsubjects)*numel(validconditions),NdimsOut,Y1.size.Nv(slice)]);
                                    dCov=zeros([NdimsOut,Y1.size.Nv(slice)]);
                                    for nvoxelbase=1:maxvoxels:Y1.size.Nv(slice), % blocks of voxels (when single-slice c data does not fit in memory)
                                        ivox=nvoxelbase:min(Y1.size.Nv(slice),nvoxelbase+maxvoxels-1);
                                        c=zeros([numel(ivox),repmat([numel(validsubjects),numel(validconditions)],[1,2])]);
                                        for nsub2=1:numel(validsubjects)
                                            for ivalidcondition2=1:numel(validconditions)
                                                ncondition2=validconditions(ivalidcondition2);
                                                if ismtxC,
                                                    Cblock=permute(reshape(conn_mtx('getblock',C,[nsub2,ivalidcondition2]),[NdimsIn numel(validsubjects) numel(validconditions) NdimsIn]),[1,4,2,3]);
                                                end
                                                for nsub=1:numel(validsubjects)
                                                    for ivalidcondition=1:numel(validconditions)
                                                        ncondition=validconditions(ivalidcondition);
                                                        if ismtxC, c(:,nsub,ivalidcondition,nsub2,ivalidcondition2)=sum((x(ivox,:,nsub,ivalidcondition)*Cblock(:,:,nsub,ivalidcondition)).*x(ivox,:,nsub2,ivalidcondition2),2);
                                                        else c(:,nsub,ivalidcondition,nsub2,ivalidcondition2)=sum((x(ivox,:,nsub,ivalidcondition)*C(:,:,nsub,ivalidcondition,nsub2,ivalidcondition2)).*x(ivox,:,nsub2,ivalidcondition2),2);
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                        %nvox0=sum(Y1.size.Nv(1:slice-1));
                                        for nvox=1:numel(ivox), %Y1.size.Nv(slice),
                                            c1=reshape(c(nvox,:),numel(validsubjects)*numel(validconditions)*[1,1]);
                                            if measures.norm{nmeasure}
                                                c1=conn_bsxfun(@minus,c1,mean(c1,1));
                                                c1=conn_bsxfun(@minus,c1,mean(c1,2));
                                            end
                                            try, [Q,D]=svd(c1);
                                            catch, [Q,D]=svds(c1,NdimsOut);
                                            end
                                            d=diag(D);
                                            Q=Q(:,1:NdimsOut);
                                            %dr=sqrt(abs(d(1:NdimsOut)));
                                            %dr=1./sqrt(mean(abs(Q-repmat(mean(Q,1),size(Q,1),1)).^2,1))';
                                            dr=sqrt(size(Q,1))+zeros(NdimsOut,1); % note: scaled to mean(q^2)=1
                                            doflip=mean(sign(Q(1:end-1+rem(size(Q,1),2),:)),1)<0;
                                            dr(doflip)=-dr(doflip);
                                            Q=conn_bsxfun(@times,Q,dr');
                                            xEig(:,:,ivox(nvox))=Q;
                                            dCov(:,ivox(nvox))=cumsum(d(1:NdimsOut))./max(eps,sum(d)); % note: this d corresponds to D^2 in manuscript (as c1 above is R*R')
                                            if maxvoxels>=NdimsOut, for ndim=1:NdimsOut, Qglobal(:,:,ndim)=Qglobal(:,:,ndim)+Q(:,ndim)*Q(:,ndim)'; end; 
                                            else for ndim=1:NdimsOut, conn_mtx('addtoblock',Qglobal,ndim,Q(:,ndim)*Q(:,ndim)'); end; 
                                            end
                                            %conn_write_voxel(Yout,Q,nvox0+nvox);
                                        end
                                    end
                                    for nsub=1:CONN_x.Setup.nsubjects
                                        isub=find(validsubjects==nsub,1);
                                        for ivalidcondition=1:numel(validconditions)
                                            ncondition=validconditions(ivalidcondition);
                                            for ndim=1:NdimsOut
                                                t=zeros(Y1.matdim.dim(1:2));
                                                if ~isempty(isub)
                                                    t(idx)=squeeze(xEig(sub2ind([numel(validsubjects),numel(validconditions)],isub,ivalidcondition),ndim,:));
                                                end
                                                filesout(nsub,ncondition,ndim)=spm_write_plane(filesout(nsub,ncondition,ndim),t,slice);
                                            end
                                        end
                                        conn_waitbar(.75*(slice+nsub/CONN_x.Setup.nsubjects-1)/Y1.size.Ns,h2,sprintf('Slice %d Subject %d',slice,nsub));
                                    end
                                    for ndim=1:NdimsOut
                                        t=zeros(Y1.matdim.dim(1:2));
                                        t(idx)=squeeze(dCov(ndim,:));
                                        filesoutCov(ndim)=spm_write_plane(filesoutCov(ndim),t,slice);
                                    end
                                    %                         conn_write_slice(Yout,reshape(xEig,[CONN_x.Setup.nsubjects*nconditions*NdimsOut,Y1.size.Nv(slice)]),slice);
                                    %                         conn_write_slice(Dout,reshape(dCov,[NdimsOut,Y1.size.Nv(slice)]),slice);
                                end
                                for ndim=1:NdimsOut, % flip orthogonal basis signs for global consistency
                                    if maxvoxels>=NdimsOut, c1=Qglobal(:,:,ndim);
                                    else c1=conn_mtx('getblock',Qglobal,ndim);
                                    end
                                    try, [Q,D]=svd(c1);
                                    catch, [Q,D]=svds(c1,1);
                                    end
                                    Q=Q(:,1);
                                    if mean(Q)<0, Q=-Q; end
                                    conn_signflip(filesout(validsubjects,validconditions,ndim),Q);
                                    conn_waitbar(.75+.25*ndim/NdimsOut,h2,sprintf('Component %d',ndim));
                                end
                                for nsub=1:CONN_x.Setup.nsubjects
                                    for ncondition=validconditions,
                                        for ndim=1:NdimsOut
                                            conn_tempcache(cache(nsub,ncondition,ndim),'nii');
                                        end
                                    end
                                end
                                conn_waitbar('close',h2);%close(h2);
                            end
                            
                        elseif measures.measuretype{nmeasure}==3||measures.measuretype{nmeasure}==4 % group-ICA, group-PCA
                            clear filesout cache;
                            for nsub=1:CONN_x.Setup.nsubjects % note: delete/reset these files even if only computing group-level analyses
                                for ncondition=validconditions,
                                    for ndim=1:NdimsOut
                                        filename=fullfile(filepathresults,['BETA_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'_Measure',num2str(imeasure(nmeasure),'%03d'),'_Component',num2str(ndim,'%03d'),'.nii']);
                                        conn_fileutils('deletefile',filename); % delete all conditions to avoid mix-up of different models
                                    end
                                end
                                for ncondition=validconditions,
                                    for ndim=1:NdimsOut
                                        filename=fullfile(filepathresults,['BETA_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'_Measure',num2str(imeasure(nmeasure),'%03d'),'_Component',num2str(ndim,'%03d'),'.nii']);
                                        [filename, cache(nsub,ncondition,ndim)]=conn_tempcache(filename);
                                        temp=struct('fname',filename,'mat',Y1.matdim.mat,'dim',Y1.matdim.dim,'n',[1,1],'pinfo',[1;0;0],'dt',[spm_type('float32'),spm_platform('bigend')],'descrip','');
                                        if Y1MD(nsub,ncondition)||~ismember(nsub,validsubjects), temp.descrip='CONNlabel:MissingData'; end
                                        filesout(nsub,ncondition,ndim)=temp;
                                    end
                                end
                            end
                            if 1,%dogrouplevel,
                                h2=conn_waitbar(0,['computing group-level dimensionality reduction']);
                                if ismtxC,
                                    randstate=randn('state');
                                    randn('seed',0);
                                    v0=randn(NdimsIn*numel(validsubjects)*numel(validconditions),1);
                                    randn('state',randstate);
                                    fC=conn_mtx('multiplicationhandle',C);
                                    [Q0,D]=eigs(fC,NdimsIn*numel(validsubjects)*numel(validconditions),NdimsOut,'LM',struct('issym',true,'isreal',true,'v0',v0));
                                    D=diag(D); % note: no info on remaining eigenvalues (fix sum(eig(C))=trace(C))
                                else
                                    C=reshape(C,[NdimsIn*numel(validsubjects)*numel(validconditions),NdimsIn*numel(validsubjects)*numel(validconditions)]);
                                    [Q0,D]=svd(C);
                                    Q0=Q0(:,1:NdimsOut);
                                    D=diag(D);
                                end
                                Q=permute(reshape(Q0,[NdimsIn,numel(validsubjects),numel(validconditions),NdimsOut]),[4,1,2,3]);
                                conn_waitbar(1,h2);
                                conn_waitbar('close',h2);%close(h2);
                                h2=conn_waitbar(0,['computing group-PCA components']);
                                Y=0;
                                for isub=1:numel(validsubjects)
                                    nsub=validsubjects(isub);
                                    for ivalidcondition=1:numel(validconditions)
                                        ncondition=validconditions(ivalidcondition);
                                        Y1.fname=fullfile(filepath,['vvPC_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'.mat']);
                                        Y1.size.Nt=Y1Nt(nsub,ncondition);
                                        y1=conn_get_volume(Y1);
                                        if ~isempty(Vmask), y1=y1(:,Vmask); end
                                        y1=y1(1:min(NdimsIn,size(y1,1)),:);
                                        Y=Y+Q(:,1:size(y1,1),isub,ivalidcondition)*y1;
                                        conn_waitbar((isub+ivalidcondition/numel(validconditions)-1)/numel(validsubjects),h2,sprintf('Condition %d Subject %d',ivalidcondition,nsub));
                                    end
                                end %note: sum(Y.^2,2)==d(1:NdimsOut)
                                idx=find(sum(Y.^3,2)<0);
                                Y(idx,:)=-Y(idx,:);
                                Q0(:,idx)=-Q0(:,idx);
                                %Q(idx,:)=-Q(idx,:); note: not used
                                clear temp_out;
                                for ndim=1:NdimsOut
                                    filename=fullfile(filepathresults,['ICAPCAcov_Measure',num2str(imeasure(nmeasure),'%03d'),'_Component',num2str(ndim,'%03d'),'.nii']);
                                    conn_fileutils('deletefile',filename);
                                    filesoutCov=struct('fname',filename,'mat',Y1.matdim.mat,'dim',Y1.matdim.dim,'n',[1,1],'pinfo',[1;0;0],'dt',[spm_type('float32'),spm_platform('bigend')],'descrip',mfilename);
                                    t=nan(Y1.matdim.dim);
                                    if ~isempty(Vmask), t(Y1.voxels(Vmask))=Y(ndim,:);
                                    else t(Y1.voxels)=Y(ndim,:);
                                    end
                                    temp_out(ndim)=spm_write_vol(filesoutCov,t);
                                end
                                filename=fullfile(filepathresults,['PCA.ROIs.nii']);
                                spm_file_merge(temp_out,filename);
                                fh=fopen(conn_prepend('',filename,'.txt'),'wt');
                                for ndim=1:NdimsOut, fprintf(fh,'%d\n',ndim);end
                                fclose(fh);
                                filename=fullfile(filepathresults,['ICAPCAcov_Measure',num2str(imeasure(nmeasure),'%03d'),'.mat']);
                                d=D.^2; d=cumsum(d)/sum(d); d=d(1:NdimsOut);
                                save(filename,'d');
                                conn_waitbar('close',h2);%close(h2);
                                if measures.measuretype{nmeasure}==3 % group-ICA
                                    [Y,W]=conn_ica(Y,[],'method',ICAMETHOD);
                                    wW=diag(1./sqrt(max(eps,mean(abs(Y).^2,2))));
                                    W=wW*W;
                                    Y=wW*Y;
                                    idx=find(sum(Y.^3,2)<0);
                                    Y(idx,:)=-Y(idx,:);
                                    W(idx,:)=-W(idx,:);
                                    idx=[];
                                    if issurface % sorting for surfaces not implemented, reverts to default (not sorting)
                                    else
                                        mask=fullfile(fileparts(which('conn')),'utils','surf','sortICA.nii'); % v1:gray/(other=white+csf+bone+soft+air)
                                        if conn_existfile(mask)
                                            maskV=spm_vol(mask);
                                            if ~isempty(Vmask), v=spm_get_data(maskV,pinv(maskV(1).mat)*xyz(:,Vmask));
                                            else v=spm_get_data(maskV,pinv(maskV(1).mat)*xyz);
                                            end
                                            vY=Y.^2*v'./(sum(Y.^2,2)*sum(v,2)');
                                            [maxvY,idxmaxvY]=max(abs(vY),[],2);
                                            [nill,idx]=sort(maxvY,'descend');
                                            filename=fullfile(filepathresults,['ICA_Measure',num2str(imeasure(nmeasure),'%03d'),'_ComponentsSort.log']);
                                            fh=fopen(filename,'wt');
                                            fprintf(fh,'ICA components by sortICA volumes match:\n');
                                            fprintf(fh,'%s\n',mat2str(vY));
                                            for nt=1:NdimsOut, fprintf(fh,'ICA %d: best match sortICA %d (%f)\n',nt,idxmaxvY(idx(nt)),maxvY(idx(nt))); end
                                            fclose(fh);
                                        end
                                    end
                                    if isempty(idx)
                                        idx=1:NdimsOut;
                                        %sY=D;sY=sY(1:NdimsOut);
                                        %[nill,idx]=sort(W.^2*sY,'descend');
                                    end
                                    Y=Y(idx,:); W=W(idx,:);
                                    if DOGICA3, WDW=W*diag(1./D(1:NdimsOut))*pinv(W); 
                                    else WDW=[];
                                    end
                                    ICAPCA='ICA';
                                else % group-PCA
                                    W=diag(1./sqrt(max(eps,mean(abs(Y).^2,2))));
                                    Y=W*Y;
                                    %W=eye(NdimsOut);
                                    WDW=[];
                                    ICAPCA='PCA';
                                end
                                clear temp_out;
                                for ndim=1:NdimsOut
                                    filename=fullfile(filepathresults,[ICAPCA,'_Measure',num2str(imeasure(nmeasure),'%03d'),'_Component',num2str(ndim,'%03d'),'.nii']);
                                    conn_fileutils('deletefile',filename);
                                    temp=struct('fname',filename,'mat',Y1.matdim.mat,'dim',Y1.matdim.dim,'n',[1,1],'pinfo',[1;0;0],'dt',[spm_type('float32'),spm_platform('bigend')],'descrip',mfilename);
                                    t=nan(Y1.matdim.dim);
                                    tY=Y(ndim,:);
                                    if ~isempty(Vmask), t(Y1.voxels(Vmask))=tY;
                                    else t(Y1.voxels)=tY;
                                    end
                                    temp_out(ndim)=spm_write_vol(temp,t);
                                end
                                filename=fullfile(filepathresults,[ICAPCA,'.Maps.nii']);
                                spm_file_merge(temp_out,filename);
                                clear temp_out;
                                for ndim=1:NdimsOut
                                    filename=fullfile(filepathresults,[ICAPCA,'_Measure',num2str(imeasure(nmeasure),'%03d'),'_Component',num2str(ndim,'%03d'),'.nii']);
                                    conn_fileutils('deletefile',filename);
                                    temp=struct('fname',filename,'mat',Y1.matdim.mat,'dim',Y1.matdim.dim,'n',[1,1],'pinfo',[1;0;0],'dt',[spm_type('float32'),spm_platform('bigend')],'descrip',mfilename);
                                    t=nan(Y1.matdim.dim);
                                    if ~isempty(WDW), tY=WDW(ndim,:)*Y;
                                    else tY=Y(ndim,:);
                                    end
                                    if ~isempty(Vmask), t(Y1.voxels(Vmask))=tY;
                                    else t(Y1.voxels)=tY;
                                    end
                                    temp_out(ndim)=spm_write_vol(temp,t);
                                end
                                filename=fullfile(filepathresults,[ICAPCA,'.ROIs.nii']);
                                spm_file_merge(temp_out,filename);
                                fh=fopen(conn_prepend('',filename,'.txt'),'wt');
                                for ndim=1:NdimsOut, fprintf(fh,'%d\n',ndim);end
                                fclose(fh);
                                if ~isempty(Vmask), Ybak=Y; Y=zeros(size(Y,1),numel(Y1.voxels)); Y(:,Vmask)=Ybak; end
                            else %%% note: obsolete section below, to be removed in future release
                                if measures.measuretype{nmeasure}==3, ICAPCA='ICA';
                                else ICAPCA='PCA';
                                end
                                WDW=[]; 
                                for ndim=1:NdimsOut
                                    filename=fullfile(filepathresults,[ICAPCA,'_Measure',num2str(imeasure(nmeasure),'%03d'),'_Component',num2str(ndim,'%03d'),'.nii']);
                                    temp=spm_vol(filename);
                                    t=spm_read_vols(temp);
                                    if ndim==1, Y=zeros(NdimsOut,numel(Y1.voxels)); end
                                    Y(ndim,:)=t(Y1.voxels);
                                end
                            end
                            
                            if 1,% dosubjectlevel,
                                h2=conn_waitbar(0,['computing subject-level components']);
                                %WD=W*diag(D(1:NdimsOut));
                                %WQ=W*Q0';
                                data={};
                                for nsub=1:CONN_x.Setup.nsubjects
                                    isub=find(validsubjects==nsub,1);
                                    for ivalidcondition=1:numel(validconditions)
                                        ncondition=validconditions(ivalidcondition);
                                        if ~isempty(isub)
                                            Y1.fname=fullfile(filepath,['vvPC_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'.mat']);
                                            Y1.size.Nt=Y1Nt(nsub,ncondition);
                                            y1=conn_get_volume(Y1);
                                            y1=y1(1:min(NdimsIn,size(y1,1)),:);
                                            % placeholder: equivalent formulation1 (faster)
                                            %wqc=WD*Q0(sub2ind([NdimsIn,numel(validsubjects),numel(validconditions)],1:NdimsIn,isub+zeros(1,NdimsIn),ivalidcondition+zeros(1,NdimsIn)),:)';
                                            % placeholder: equivalent formulation2
                                            %if ismtxC, wqc=WQ*conn_mtx('getblock',C,[nsub,ivalidcondition]);
                                            %else wqc=WQ*C(:,sub2ind([NdimsIn,numel(validsubjects),numel(validconditions)],1:NdimsIn,isub+zeros(1,NdimsIn),ivalidcondition+zeros(1,NdimsIn)));
                                            %end
                                            if isempty(WDW), % (GICA1) wqc'*y2=y1
                                                wqc=Y*y1'/size(y1,2);
                                                y2=pinv(wqc*wqc')*wqc(:,1:size(y1,1))*y1;
                                            else % (GICA3)
                                                wqc=WDW*Y*y1'; %(W*diag(1./D(1:NdimsOut))*pinv(W))*Y*y1';
                                                y2=wqc(:,1:size(y1,1))*y1 *(numel(validconditions)*CONN_x.Setup.nsubjects);
                                            end
                                            filename_D1=fullfile(filepath,['vvPCeig_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'.mat']);
                                            D1=load(filename_D1,'Q1');
                                            data{nsub}{ivalidcondition}=D1.Q1(:,1:size(wqc,2))*wqc';
                                        else data{nsub}{ivalidcondition}=[];
                                        end
                                        for ndim=1:NdimsOut
                                            filename=fullfile(filepathresults,['BETA_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'_Measure',num2str(imeasure(nmeasure),'%03d'),'_Component',num2str(ndim,'%03d'),'.nii']);
                                            t=zeros(Y1.matdim.dim);
                                            if ~isempty(isub), t(Y1.voxels)=y2(ndim,:); end
                                            spm_write_vol(filesout(nsub,ncondition,ndim),t);
                                            redone_files=redone_files+1;
                                        end
                                        conn_waitbar((nsub+ivalidcondition/numel(validconditions)-1)/CONN_x.Setup.nsubjects,h2,sprintf('Condition %d Subject %d',ivalidcondition,nsub));
                                    end
                                end
                                if measures.measuretype{nmeasure}==3, ICAPCA='ICA';
                                else ICAPCA='PCA';
                                end
                                filename=fullfile(filepathresults,[ICAPCA,'.Timeseries.mat']);
                                conditions=CONN_x.Setup.conditions.names(validconditions);
                                weights=Y1weights(:,validconditions);
                                save(filename,'data','conditions','weights');
                                for nsub=1:CONN_x.Setup.nsubjects
                                    for ncondition=validconditions,
                                        for ndim=1:NdimsOut
                                            conn_tempcache(cache(nsub,ncondition,ndim),'nii');
                                        end
                                    end
                                end
                                conn_waitbar('close',h2);%close(h2);
                                
                                if 1, % saves new second-level covariates
                                    
                                    iremove=reshape(find(cellfun('length',regexp(CONN_x.Setup.l2covariates.names(1:end-1),['^_(Variability|Average|Frequency) ',ICAPCA,'\d+ ',CONN_x.vvAnalyses(ianalysis).name,' @ .*?$']))),1,[]);
                                    %iremove=[];
                                    %for ivalidcondition=1:numel(validconditions),
                                    %    iremove=[iremove,reshape(find(cellfun('length',regexp(CONN_x.Setup.l2covariates.names(1:end-1),['^_(Variability|Average|Frequency) ',ICAPCA,'\d+ ',CONN_x.vvAnalyses(ianalysis).name,' @ ',conditions{ivalidcondition},'$']))),1,[])];
                                    %end
                                    if ~isempty(iremove)
                                        ikeep=setdiff(1:numel(CONN_x.Setup.l2covariates.names),iremove);
                                        ikeep0=setdiff(ikeep,numel(CONN_x.Setup.l2covariates.names));
                                        CONN_x.Setup.l2covariates.names=CONN_x.Setup.l2covariates.names(ikeep);
                                        CONN_x.Setup.l2covariates.descrip=CONN_x.Setup.l2covariates.descrip(ikeep0);
                                        for nsub=1:CONN_x.Setup.nsubjects,
                                            CONN_x.Setup.l2covariates.values{nsub}=CONN_x.Setup.l2covariates.values{nsub}(ikeep0);
                                        end
                                    end
                                    new_names={};
                                    new_values={};
                                    maxrt=[]; for nsub=validsubjects, maxrt(nsub)=max(conn_get_rt(nsub)); end
                                    for ncomp=1:NdimsOut
                                        for ivalidcondition=1:numel(validconditions),
                                            VARt=nan(CONN_x.Setup.nsubjects,1);
                                            AVGt=nan(CONN_x.Setup.nsubjects,1);
                                            ZCt=nan(CONN_x.Setup.nsubjects,1);
                                            for nsub=validsubjects
                                                mask=weights{nsub,ivalidcondition}{1}>0;
                                                b=data{nsub}{ivalidcondition}(mask,ncomp)*sqrt(nnz(mask));
                                                remove=weights{nsub,ivalidcondition}{5}(mask)==1;
                                                VARt(nsub)=std(b,0,1);
                                                AVGt(nsub)=mean(b,1);
                                                b=b-mean(b);
                                                b(remove)=0;
                                                ZCt(nsub)=mean((b(1:end-1,:)<=0&b(2:end,:)>0)|(b(1:end-1,:)>=0&b(2:end,:)<0))/maxrt(nsub)/2;
                                            end
                                            if isempty(CONN_x.vvAnalyses(ianalysis).name)
                                                new_names{end+1}=sprintf('_Variability %s%02d @ %s',ICAPCA,ncomp,conditions{ivalidcondition});
                                                new_names{end+1}=sprintf('_Average %s%02d @ %s',ICAPCA,ncomp,conditions{ivalidcondition});
                                                new_names{end+1}=sprintf('_Frequency %s%02d @ %s (%s)',ICAPCA,ncomp,conditions{ivalidcondition});
                                            else
                                                new_names{end+1}=sprintf('_Variability %s%02d %s @ %s',ICAPCA,ncomp,CONN_x.vvAnalyses(ianalysis).name,conditions{ivalidcondition});
                                                new_names{end+1}=sprintf('_Average %s%02d %s @ %s',ICAPCA,ncomp,CONN_x.vvAnalyses(ianalysis).name,conditions{ivalidcondition});
                                                new_names{end+1}=sprintf('_Frequency %s%02d %s @ %s',ICAPCA,ncomp,CONN_x.vvAnalyses(ianalysis).name,conditions{ivalidcondition});
                                            end
                                            new_values{end+1}=VARt;
                                            new_values{end+1}=AVGt;
                                            new_values{end+1}=ZCt;
                                        end
                                    end
                                    for nvar=1:numel(new_names)
                                        icov=strmatch(new_names{nvar},CONN_x.Setup.l2covariates.names,'exact');
                                        if isempty(icov), icov=numel(CONN_x.Setup.l2covariates.names); CONN_x.Setup.l2covariates.names{end+1}=' '; end
                                        CONN_x.Setup.l2covariates.names{icov}=new_names{nvar};
                                        CONN_x.Setup.l2covariates.descrip{icov}='';
                                        for nsub=1:CONN_x.Setup.nsubjects,
                                            CONN_x.Setup.l2covariates.values{nsub}{icov}=new_values{nvar}(nsub);
                                        end
                                    end
                                end
                            end
                        end
                    end
                    n=n+1;
                    %                 conn_waitbar(n/N,h);
                end
                try
                    filename=fullfile(filepathresults,['TEMPORAL1_Measure',num2str(imeasure(nmeasures1+1),'%03d'),'.mat']);
                    spm_unlink(filename);
                end
            else
                n=n+nmeasures2;
                %             conn_waitbar(n/N,h);
            end
            %         conn_waitbar('close',h);
        end
        
        if REPLACERESULTS,%~isfield(CONN_x.vvAnalyses(ianalysis),'measures')||isempty(CONN_x.vvAnalyses(ianalysis).measures), 
            CONN_x.vvAnalyses(ianalysis).measures=conn_v2v('list_extended',CONN_x.vvAnalyses(ianalysis).regressors);
        else % remove from list of displayed results other analyses by the same name but otherwise keep existing results
            newmeasures=conn_v2v('list_extended',CONN_x.vvAnalyses(ianalysis).regressors);
            ko=[]; for n1=1:numel(newmeasures), ko=[ko conn_v2v('match_existing',newmeasures{n1})]; end
            if any(ko), CONN_x.vvAnalyses(ianalysis).measures=CONN_x.vvAnalyses(ianalysis).measures(setdiff(1:numel(CONN_x.vvAnalyses(ianalysis).measures),ko)); end
            CONN_x.vvAnalyses(ianalysis).measures=[CONN_x.vvAnalyses(ianalysis).measures newmeasures];
        end
        try
            fileout=fullfile(filepathresults,'_list_measures.txt');
            fh=fopen(fileout,'wt');
            for n1=1:length(CONN_x.vvAnalyses(ianalysis).measurenames),fprintf(fh,'Measure%03d = %s\n',n1,CONN_x.vvAnalyses(ianalysis).measurenames{n1});end
            fclose(fh);
        end
        try
            fileout=fullfile(filepathresults,'_list_conditions.txt');
            fh=fopen(fileout,'wt');
            for ncondition=1:nconditions,
                fprintf(fh,'Condition%03d = %s\n',icondition(ncondition),CONN_x.Setup.conditions.names{ncondition});
            end
            fclose(fh);
        end
        conn_disp('fprintf','      processed %d files\n',redone_files);
    end
    CONN_x.vvAnalysis=analysisbak;
    CONN_x.isready(4)=1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates dyn_Subject###.mat files (temporal components for dynamic connectivity analyses)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(floor(options)==14) && any(CONN_x.Setup.steps([4])) && ~(isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'steps')&&numel(CONN_x.gui.steps)>3&&~any(CONN_x.gui.steps([4]))),
    [path,name,ext]=fileparts(CONN_x.filename);
    if nargin>1&&~isempty(varargin{1}),analyses=varargin{1}; % selected analysis only
    else, analyses=1:length(CONN_x.dynAnalyses); 
    end
    if ischar(analyses)||iscell(analyses), analyses=find(ismember({CONN_x.dynAnalyses.name},analyses)); end
    analyses(analyses<=0|analyses>numel(CONN_x.dynAnalyses))=[];
    analysisbak=CONN_x.dynAnalysis;
    filepath=CONN_x.folders.preprocessing;
	nconditions=length(CONN_x.Setup.conditions.names)-1;
    if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'subjects'), validsubjects=CONN_x.gui.subjects; else validsubjects=1:CONN_x.Setup.nsubjects; end
    if isfield(CONN_x,'pobj')&&isstruct(CONN_x.pobj)&&isfield(CONN_x.pobj,'subjects'), validsubjects=CONN_x.pobj.subjects; end
    if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'conditions')&&~isempty(CONN_x.gui.conditions), validconditions=CONN_x.gui.conditions; else validconditions=1:length(CONN_x.Setup.conditions.names)-1; end
    icondition=[];isnewcondition=[];for ncondition=validconditions,[icondition(ncondition),isnewcondition(ncondition)]=conn_conditionnames(CONN_x.Setup.conditions.names{ncondition},'+'); end
    %h=conn_waitbar(0,['Step ',num2str(sum(options<=14)),'/',num2str(length(options)),': dynamic temporal-components estimation']);
    if any(arrayfun(@(n)isempty(dir(fullfile(filepath,['ROI_Subject',num2str(n,'%03d'),'_Condition',num2str(0,'%03d'),'.mat']))),validsubjects)), error(['Data not ready yet. Re-run Denoising step']); end
    dogrouplevel=any(options==14|options==14.1); if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'grouplevel')&&~ismember(CONN_x.gui.grouplevel,[0 1]),dogrouplevel=false; end 
    dosubjectlevel=any(options==14|options==14.2); if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'grouplevel')&&~ismember(CONN_x.gui.grouplevel,[0 2]),dosubjectlevel=false; end 
    for nanalyses=1:length(analyses),
        ianalysis=analyses(nanalyses);
        CONN_x.dynAnalysis=ianalysis;
        filepathresults=fullfile(CONN_x.folders.firstlevel_dyn,CONN_x.dynAnalyses(ianalysis).name);
        conn_disp('fprintf','      first-level data in %s\n',filepathresults);
    
        Ncomponents=CONN_x.dynAnalyses(CONN_x.dynAnalysis).Ncomponents; % number of components to estimate
        selectedcondition=CONN_x.dynAnalyses(CONN_x.dynAnalysis).condition;    % selected condition (determines span of BOLD timeseries to use in these analyses)
        roinames=CONN_x.dynAnalyses(CONN_x.dynAnalysis).regressors.names;  % cell-array of ROI names (empty to use all ROIs)
        window=CONN_x.dynAnalyses(CONN_x.dynAnalysis).window;
        roixyzs={};
        
%         do=true;
%         if ~isequal(validsubjects,1:CONN_x.Setup.nsubjects),
%             if isfield(CONN_x,'pobj')&&isstruct(CONN_x.pobj)&&isfield(CONN_x.pobj,'subjects'),
%                 if isfield(CONN_x.pobj,'partition')&&isequal(CONN_x.pobj.partition,[1 1])
%                     conn_disp('NOTE: single-job, subset of subjects only');
%                 elseif subjectlevelonly
%                     conn_disp('NOTE: performing subject-level analyses only');
%                 else
%                     subjectlevelonly=true;
%                     %                     else
%                     %                         conn_disp('WARNING: Dynamic factor analysis parallelization not yet available (run locally or as a single job instead). Skipping these analyses');
%                     %                         do=false;
%                 end
%             end
%         end
        REDO='Yes'; filename=fullfile(filepathresults,['dyn_Subject',num2str(validsubjects(1),'%03d'),'.mat']);
        if conn_existfile(filename),if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'overwrite'), REDO=CONN_x.gui.overwrite; else, REDO=conn_questdlg('Overwrite existing subject results?','','Yes', 'No', 'Yes');end; end
        if dogrouplevel && (~conn_existfile(filename)||strcmp(lower(REDO),'yes'))
            do=true;
            if ~isequal(validsubjects,1:CONN_x.Setup.nsubjects),
                if isfield(CONN_x,'pobj')&&isstruct(CONN_x.pobj)&&isfield(CONN_x.pobj,'subjects'),
                    if isfield(CONN_x.pobj,'partition')&&isequal(CONN_x.pobj.partition,[1 1])
                        conn_disp('NOTE: single-job, subset of subjects only');
                    else
                        conn_disp('WARNING: Dynamic factor analysis parallelization not yet available (run locally or as a single job instead). Skipping these analyses');
                        do=false;
                        dosubjectlevel=false;
                        newanalyses=[];
                    end
                end
            end
            if do
                REDO='Yes';
                DEMEAN=false; % true=disregard between-subject variability in average connectivity; set to true for back-compatibility
                X=[];
                Y=[];
                Xfilter=[];
                IDX_subject=[];
                IDX_session=[];
                COND_names=CONN_x.Setup.conditions.names(validconditions);
                COND_weights=cell(size(COND_names));
                if 1 % extract temporal measures only across selected conditions
                    COND_names_all=COND_names;
                    COND_weights_all=COND_weights;
                else % extract temporal measures across all conditions
                    COND_names_all=CONN_x.Setup.conditions.names(1:nconditions);
                    COND_weights_all=cell(size(COND_names_all));
                end
                for nsub=validsubjects,
                    filename=fullfile(filepath,['ROI_Subject',num2str(nsub,'%03d'),'_Condition',num2str(0,'%03d'),'.mat']);
                    if ~conn_existfile(filename), conn_disp(['Not ready to process step conn_process_135']); return; end %conn_waitbar('close',h);return; end
                    X1=load(filename); 
                    if isempty(roinames), roinames=X1.names(~ismember(X1.names,{'White Matter','CSF','MotionMask'})); end
                    [ok,idx]=ismember(roinames,X1.names);
                    if ~all(ok), error('Missing ROI data in subject %d',nsub); end
                    y=cat(2,X1.data{idx});
                    if isempty(roixyzs), roixyzs=X1.xyz(idx); end
                    if nsub~=validsubjects(1)&&size(y,2)~=size(Y,2), error('Incorrect number of ROI timeseries in subject %d',nsub); end
                    if ~isempty(selectedcondition)
                        matchedcondition=find(strcmp(X1.conditionsnames,CONN_x.Setup.conditions.names{selectedcondition}),1);
                        if isempty(matchedcondition), error('Condition name %s not found in %s. Please re-run Denoising step',CONN_x.Setup.conditions.names{selectedcondition},filename); end
                        if numel(X1.conditionsweights)>=matchedcondition&&~isempty(X1.conditionsweights{matchedcondition})
                            wx=X1.conditionsweights{matchedcondition}{1}; % hrf
                            wx=max(0,wx);
                            y=conn_wdemean(y,wx);
                            y=y.*repmat(wx,[1,size(y,2)]);
                        end
                    end
                    y=conn_bsxfun(@rdivide,y,max(eps,sqrt(mean(abs(y).^2,1))));
                    if any(all(y==0,1)), conn_disp(sprintf('Warning: %d null ROI timeseries in subject %d',sum(all(y==0,1)),nsub)); end
                    Y=cat(1,Y,y);
                    X=[X zeros(size(X,1),1); zeros(size(y,1),size(X,2)) ones(size(y,1),1)]; % subject-effects
                    xfilter=[];
                    tr=inf;for nses=1:max(X1.data_sessions), tr=min(tr,conn_get_rt(nsub,nses)); end
                    if ~isempty(window)&&window>tr,
                        nsess=max(X1.data_sessions);
                        xfilter=[];
                        for nses=1:nsess
                            tr=conn_get_rt(nsub,nses);
                            nt=CONN_x.Setup.nscans{nsub}{nses};
                            nh=min(nt,1+2*floor(window/tr));
                            temp=spm_convmtx(conn_hanning(nh),nt);
                            temp=temp(floor(nh/2)+(1:nt),:);
                            [i,j,v]=find(temp);
                            xf=sparse(i,j,v,nt,nt);
                            xfilter=[xfilter sparse(size(xfilter,1),size(xf,2)); sparse(size(xf,1),size(xfilter,2)) xf]; % within-session temporal smoothing kernel
                        end
                    end
                    Xfilter=[Xfilter sparse(size(Xfilter,1),size(xfilter,2)); sparse(size(xfilter,1),size(Xfilter,2)) xfilter];
                    IDX_subject=cat(1,IDX_subject,nsub+zeros(size(y,1),1));
                    IDX_session=cat(1,IDX_session,X1.data_sessions);
                    ok=cellfun(@(x)numel(x{1})==numel(X1.data_sessions),X1.conditionsweights,'ErrorHandler',@(varargin)false);
                    %if nsub==validsubjects(1), COND_names=X1.conditionsnames(ok); COND_weights=cell(1,nnz(ok)); end
                    ok1=find(ok&ismember(X1.conditionsnames(1:numel(ok)),COND_names));
                    [ok2,iok2]=ismember(COND_names,X1.conditionsnames(ok1));
                    for n=find(ok2), COND_weights{n}=cat(1,COND_weights{n},X1.conditionsweights{ok1(iok2(n))}{1}); end
                    COND_names=COND_names(ok2);
                    COND_weights=COND_weights(ok2);
                    ok1=find(ok&ismember(X1.conditionsnames(1:numel(ok)),COND_names_all));
                    [ok2,iok2]=ismember(COND_names_all,X1.conditionsnames(ok1));
                    for n=find(ok2), COND_weights_all{n}=cat(1,COND_weights_all{n},X1.conditionsweights{ok1(iok2(n))}{1}); end
                    COND_names_all=COND_names_all(ok2);
                    COND_weights_all=COND_weights_all(ok2);
                end
                drawnow;
                if DEMEAN, tX=X; else tX=[]; end
                [H,B,H0,B0]=conn_invPPI(Y,Ncomponents,tX,Xfilter,1);
                Ncomponents=size(H,2);
                CONN_x.dynAnalyses(CONN_x.dynAnalysis).sources=roinames;
                filename=fullfile(filepathresults,'dyn_Base.mat');
                names=arrayfun(@(n)sprintf('Dynamic factor %02d',n),1:Ncomponents,'uni',0);
                ROInames=roinames;
                ROIxyzs=roixyzs;
                if ~isempty(selectedcondition), selectedconditionname=CONN_x.Setup.conditions.names{selectedcondition};
                else selectedconditionname='';
                end
                save(filename,'X','Y','Xfilter','B','H','B0','H0','names','ROInames','ROIxyzs','IDX_subject','IDX_session','COND_names','COND_weights','selectedconditionname');
                %load(filename,'H','IDX_subject','IDX_session','COND_names','COND_weights','selectedconditionname')
                for nsub=1:CONN_x.Setup.nsubjects
                    filename=fullfile(filepathresults,['dyn_Subject',num2str(nsub,'%03d'),'.mat']);
                    thissub=IDX_subject==nsub;
                    nthissub=nnz(thissub);
                    data=H(thissub,:);
                    names=arrayfun(@(n)sprintf('Dynamic factor %02d',n),1:Ncomponents,'uni',0);
                    data_sessions=IDX_session(IDX_subject==nsub);
                    save(filename,'data','names','data_sessions');
                    filename=fullfile(filepathresults,['dyn_Subject',num2str(nsub,'%03d'),'.mtx.nii']);
                    C=0; 
                    for ncomp=1:size(H0,2), C=C+repmat(B0(:,:,ncomp),[1,1,nthissub]).*repmat(shiftdim(H0(thissub,ncomp),-2),[size(B0,1),size(B0,2),1]); end
                    for ncomp=1:size(H ,2), C=C+repmat( B(:,:,ncomp),[1,1,nthissub]).*repmat(shiftdim( H(thissub,ncomp),-2),[size(B ,1),size(B ,2),1]); end
                    %C=sum(permute(B0,[1,2,4,3]).*permute(H0(thissub,:),[3,4,1,2]),4)+sum(permute(B,[1,2,4,3]).*permute(H(thissub,:),[3,4,1,2]),4);
                    conn_mtx_write(filename,C,ROInames,ROIxyzs);
                end

                iremove=reshape(find(cellfun('length',regexp(CONN_x.Setup.l2covariates.names(1:end-1),['^_(Variability|Average|Frequency) Dynamic factor (.*_)?\d+ ',CONN_x.dynAnalyses(CONN_x.dynAnalysis).name,' @ .*$']))),1,[]);
                if ~isempty(iremove)
                    ikeep=setdiff(1:numel(CONN_x.Setup.l2covariates.names),iremove);
                    ikeep0=setdiff(ikeep,numel(CONN_x.Setup.l2covariates.names));
                    CONN_x.Setup.l2covariates.names=CONN_x.Setup.l2covariates.names(ikeep);
                    CONN_x.Setup.l2covariates.descrip=CONN_x.Setup.l2covariates.descrip(ikeep0);
                    for nsub=1:CONN_x.Setup.nsubjects,
                        CONN_x.Setup.l2covariates.values{nsub}=CONN_x.Setup.l2covariates.values{nsub}(ikeep0);
                    end
                end
                new_names={};
                new_values={};
                %                 for ncomp=1:Ncomponents
                %                     for n=1:numel(COND_names)
                %                         w=max(eps,accumarray(IDX_subject,max(0,COND_weights{n}),[CONN_x.Setup.nsubjects,1]));
                %                         H_avg_weighted=accumarray(IDX_subject,H(:,ncomp).*max(0,COND_weights{n}),[CONN_x.Setup.nsubjects,1],@sum,nan)./w;
                %                         if isempty(selectedcondition)||strcmp(CONN_x.Setup.conditions.names{selectedcondition},'rest'), new_names{end+1}=sprintf('Dynamic baseline factor %02d @ %s',ncomp,COND_names{n});
                %                         else new_names{end+1}=sprintf('Dynamic baseline factor %s_%02d @ %s',CONN_x.Setup.conditions.names{selectedcondition},ncomp,COND_names{n});
                %                         end
                %                         new_values{end+1}=H_avg_weighted;
                %                     end
                %                 end
                for ncomp=1:Ncomponents
                    %             H_std=accumarray(IDX_subject,H(:,ncomp),[CONN_x.Setup.nsubjects,1],@std);
                    %             new_names{end+1}=sprintf('Dynamic factor %02d',ncomp);
                    %             new_values{end+1}=H_std;
                    for n=1:numel(COND_names_all)
                        w=max(eps,accumarray(IDX_subject,max(0,COND_weights_all{n}),[CONN_x.Setup.nsubjects,1]));
                        H_std_weighted=sqrt(max(0,accumarray(IDX_subject,H(:,ncomp).^2.*max(0,COND_weights_all{n}),[CONN_x.Setup.nsubjects,1],@sum,nan)./w - (accumarray(IDX_subject,H(:,ncomp).*max(0,COND_weights_all{n}),[CONN_x.Setup.nsubjects,1],@sum,nan)./w).^2));
                        if isempty(CONN_x.dynAnalyses(CONN_x.dynAnalysis).name)
                            if isempty(selectedcondition)||strcmp(CONN_x.Setup.conditions.names{selectedcondition},'rest'), new_names{end+1}=sprintf('_Variability Dynamic factor %02d @ %s',ncomp,COND_names_all{n});
                            else new_names{end+1}=sprintf('_Variability Dynamic factor %s_%02d @ %s',CONN_x.Setup.conditions.names{selectedcondition},ncomp,COND_names_all{n});
                            end
                        else
                            if isempty(selectedcondition)||strcmp(CONN_x.Setup.conditions.names{selectedcondition},'rest'), new_names{end+1}=sprintf('_Variability Dynamic factor %02d %s @ %s',ncomp,CONN_x.dynAnalyses(CONN_x.dynAnalysis).name,COND_names_all{n});
                            else new_names{end+1}=sprintf('_Variability Dynamic factor %s_%02d %s @ %s',CONN_x.Setup.conditions.names{selectedcondition},ncomp,CONN_x.dynAnalyses(CONN_x.dynAnalysis).name,COND_names_all{n});
                            end
                        end
                        H_std_weighted(setdiff(1:CONN_x.Setup.nsubjects,validsubjects))=nan;
                        new_values{end+1}=H_std_weighted;
                    end
                end
                for ncomp=1:Ncomponents
                    for n=1:numel(COND_names_all)
                        w=max(eps,accumarray(IDX_subject,max(0,COND_weights_all{n}),[CONN_x.Setup.nsubjects,1]));
                        H_avg_weighted=accumarray(IDX_subject,H(:,ncomp).*max(0,COND_weights_all{n}),[CONN_x.Setup.nsubjects,1],@sum,nan)./w;
                        if isempty(CONN_x.dynAnalyses(CONN_x.dynAnalysis).name)
                            if isempty(selectedcondition)||strcmp(CONN_x.Setup.conditions.names{selectedcondition},'rest'), new_names{end+1}=sprintf('_Average Dynamic factor %02d @ %s',ncomp,COND_names_all{n});
                            else new_names{end+1}=sprintf('_Average Dynamic factor %s_%02d @ %s',CONN_x.Setup.conditions.names{selectedcondition},ncomp,COND_names_all{n});
                            end
                        else
                            if isempty(selectedcondition)||strcmp(CONN_x.Setup.conditions.names{selectedcondition},'rest'), new_names{end+1}=sprintf('_Average Dynamic factor %02d %s @ %s',ncomp,CONN_x.dynAnalyses(CONN_x.dynAnalysis).name,COND_names_all{n});
                            else new_names{end+1}=sprintf('_Average Dynamic factor %s_%02d %s @ %s',CONN_x.Setup.conditions.names{selectedcondition},ncomp,CONN_x.dynAnalyses(CONN_x.dynAnalysis).name,COND_names_all{n});
                            end
                        end
                        H_avg_weighted(setdiff(1:CONN_x.Setup.nsubjects,validsubjects))=nan;
                        new_values{end+1}=H_avg_weighted;
                    end
                end
                for n=1:numel(COND_names_all)
                    w=max(eps,accumarray(IDX_subject,max(0,COND_weights_all{n}),[CONN_x.Setup.nsubjects,1]));
                    H_std_weighted=0;
                    for ncomp=1:Ncomponents
                        H_std_weighted=H_std_weighted+(max(0,accumarray(IDX_subject,H(:,ncomp).^2.*max(0,COND_weights_all{n}),[CONN_x.Setup.nsubjects,1],@sum,nan)./w - (accumarray(IDX_subject,H(:,ncomp).*max(0,COND_weights_all{n}),[CONN_x.Setup.nsubjects,1],@sum,nan)./w).^2));
                    end
                    if isempty(CONN_x.dynAnalyses(CONN_x.dynAnalysis).name)
                        if isempty(selectedcondition)||strcmp(CONN_x.Setup.conditions.names{selectedcondition},'rest'), new_names{end+1}=sprintf('_Variability Dynamic Total @ %s',COND_names_all{n});
                        else new_names{end+1}=sprintf('_Variability Dynamic Total_%s @ %s',CONN_x.Setup.conditions.names{selectedcondition},COND_names_all{n});
                        end
                    else
                        if isempty(selectedcondition)||strcmp(CONN_x.Setup.conditions.names{selectedcondition},'rest'), new_names{end+1}=sprintf('_Variability Dynamic Total %s @ %s',CONN_x.dynAnalyses(CONN_x.dynAnalysis).name,COND_names_all{n});
                        else new_names{end+1}=sprintf('_Variability Dynamic Total_%s %s @ %s',CONN_x.Setup.conditions.names{selectedcondition},CONN_x.dynAnalyses(CONN_x.dynAnalysis).name,COND_names_all{n});
                        end
                    end
                    H_std_weighted(setdiff(1:CONN_x.Setup.nsubjects,validsubjects))=nan;
                    new_values{end+1}=sqrt(H_std_weighted);
                end
                allrt=[];
                for ncomp=1:Ncomponents
                    for n=1:numel(COND_names_all)
                        tempW=max(0,COND_weights_all{n});
                        w=max(eps,accumarray(IDX_subject,tempW,[CONN_x.Setup.nsubjects,1]));
                        H_avg_weighted=accumarray(IDX_subject,H(:,ncomp).*tempW,[CONN_x.Setup.nsubjects,1],@sum,nan)./w;
                        tempH=(H(:,ncomp)<=H_avg_weighted(IDX_subject,:)&H([2:end end],ncomp)>H_avg_weighted(IDX_subject,:))|(H(:,ncomp)>=H_avg_weighted(IDX_subject,:)&H([2:end end],ncomp)<H_avg_weighted(IDX_subject,:));
                        maskout=IDX_subject~=[IDX_subject(2:end);nan]|IDX_session~=[IDX_session(2:end);nan];
                        tempW(maskout)=0;
                        tempH(maskout)=0;
                        w=max(eps,accumarray(IDX_subject,tempW,[CONN_x.Setup.nsubjects,1]));
                        H_freq_weighted=max(0,accumarray(IDX_subject,tempH.*tempW,[CONN_x.Setup.nsubjects,1],@sum,nan))./w;
                        if isempty(allrt), allrt=conn_get_rt; end
                        H_freq_weighted=H_freq_weighted./reshape(allrt,CONN_x.Setup.nsubjects,1)/2;
                        if isempty(CONN_x.dynAnalyses(CONN_x.dynAnalysis).name)
                            if isempty(selectedcondition)||strcmp(CONN_x.Setup.conditions.names{selectedcondition},'rest'), new_names{end+1}=sprintf('_Frequency Dynamic factor %02d %s @ %s',ncomp,COND_names_all{n});
                            else new_names{end+1}=sprintf('_Frequency Dynamic factor %s_%02d @ %s',CONN_x.Setup.conditions.names{selectedcondition},ncomp,COND_names_all{n});
                            end
                        else
                            if isempty(selectedcondition)||strcmp(CONN_x.Setup.conditions.names{selectedcondition},'rest'), new_names{end+1}=sprintf('_Frequency Dynamic factor %02d %s @ %s',ncomp,CONN_x.dynAnalyses(CONN_x.dynAnalysis).name,COND_names_all{n});
                            else new_names{end+1}=sprintf('_Frequency Dynamic factor %s_%02d @ %s',CONN_x.Setup.conditions.names{selectedcondition},ncomp,CONN_x.dynAnalyses(CONN_x.dynAnalysis).name,COND_names_all{n});
                            end
                        end
                        H_freq_weighted(setdiff(1:CONN_x.Setup.nsubjects,validsubjects))=nan;
                        new_values{end+1}=H_freq_weighted;
                    end
                end
                for nvar=1:numel(new_names)
                    icov=strmatch(new_names{nvar},CONN_x.Setup.l2covariates.names,'exact');
                    if isempty(icov), icov=numel(CONN_x.Setup.l2covariates.names); CONN_x.Setup.l2covariates.names{end+1}=' '; end
                    CONN_x.Setup.l2covariates.names{icov}=new_names{nvar};
                    CONN_x.Setup.l2covariates.descrip{icov}='';
                    for n1=1:CONN_x.Setup.nsubjects, CONN_x.Setup.l2covariates.values{n1}{icov}=new_values{nvar}(n1); end
                end
                
                if 1,
                    if isempty(CONN_x.dynAnalyses(CONN_x.dynAnalysis).name)
                        if isempty(selectedcondition)||strcmp(CONN_x.Setup.conditions.names{selectedcondition},'rest'), name='QC_dynamicfactors';
                        else name=sprintf('QC_dynamicfactors_%s',CONN_x.Setup.conditions.names{selectedcondition});
                        end
                    else
                        if isempty(selectedcondition)||strcmp(CONN_x.Setup.conditions.names{selectedcondition},'rest'), name=sprintf('QC_dynamicfactors_%s',CONN_x.dynAnalyses(CONN_x.dynAnalysis).name);
                        else name=sprintf('QC_dynamicfactors_%s_%s',CONN_x.dynAnalyses(CONN_x.dynAnalysis).name,CONN_x.Setup.conditions.names{selectedcondition});
                        end
                    end
                    idx=strmatch(name,CONN_x.Setup.l1covariates.names,'exact');
                    if isempty(idx), idx=length(CONN_x.Setup.l1covariates.names); CONN_x.Setup.l1covariates.names{end+1}=' '; end
                    CONN_x.Setup.l1covariates.names{idx}=name;
                    for nsub=1:CONN_x.Setup.nsubjects
                        nsess=max(IDX_session(IDX_subject==nsub));
                        for nses=1:nsess
                            samples=IDX_subject==nsub&IDX_session==nses;
                            filename=fullfile(filepathresults,['dyn_Subject',num2str(nsub,CONN_x.opt.fmt1),'_Session',num2str(nses,'%03d'),'.txt']);
                            conn_savetextfile(filename,H(samples,:));
                            CONN_x.Setup.l1covariates.files{nsub}{idx}{nses}=conn_file(filename);
                        end
                    end
                end
                if CONN_x.dynAnalyses(CONN_x.dynAnalysis).output(3)
                    for ncomp=1:Ncomponents
                        if isempty(CONN_x.dynAnalyses(CONN_x.dynAnalysis).name)
                            if isempty(selectedcondition)||strcmp(CONN_x.Setup.conditions.names{selectedcondition},'rest'), name=sprintf('Dynamic factor %02d',ncomp);
                            else name=sprintf('Dynamic factor %s_%02d',CONN_x.Setup.conditions.names{selectedcondition},ncomp);
                            end
                        else
                            if isempty(selectedcondition)||strcmp(CONN_x.Setup.conditions.names{selectedcondition},'rest'), name=sprintf('Dynamic factor %02d %s',ncomp,CONN_x.dynAnalyses(CONN_x.dynAnalysis).name);
                            else name=sprintf('Dynamic factor %s_%02d %s',CONN_x.Setup.conditions.names{selectedcondition},ncomp,CONN_x.dynAnalyses(CONN_x.dynAnalysis).name);
                            end
                        end
                        idx=strmatch(name,CONN_x.Setup.l1covariates.names,'exact');
                        if isempty(idx), idx=length(CONN_x.Setup.l1covariates.names); CONN_x.Setup.l1covariates.names{end+1}=' '; end
                        CONN_x.Setup.l1covariates.names{idx}=name;
                        for nsub=1:CONN_x.Setup.nsubjects
                            nsess=max(IDX_session(IDX_subject==nsub));
                            for nses=1:nsess
                                samples=IDX_subject==nsub&IDX_session==nses;
                                CONN_x.Setup.l1covariates.files{nsub}{idx}{nses}={'[raw values]',[],H(samples,ncomp)};
                            end
                        end
                    end
                end
                
                %                 end
                %                 if 1, %CONN_x.dynAnalyses(CONN_x.dynAnalysis).output(1)
                newanalyses=[];
                for ncomp=1:Ncomponents
                    if isempty(selectedcondition)||strcmp(CONN_x.Setup.conditions.names{selectedcondition},'rest'), name=sprintf('Dynamic factor %02d',ncomp);
                    else name=sprintf('Dynamic factor %s_%02d',CONN_x.Setup.conditions.names{selectedcondition},ncomp);
                    end
                    name=fullfile(CONN_x.dynAnalyses(CONN_x.dynAnalysis).name,name);
                    names={CONN_x.Analyses.name};
                    tianalysis=find(strcmp(names,name));
                    if isempty(tianalysis), tianalysis=numel(CONN_x.Analyses)+1; end
                    newanalyses=[newanalyses tianalysis];
                    CONN_x.Analyses(tianalysis).name=name;
                    CONN_x.Analyses(tianalysis).modulation=fullfile(CONN_x.dynAnalyses(CONN_x.dynAnalysis).name,sprintf('Dynamic factor %02d',ncomp));
                    CONN_x.Analyses(tianalysis).measure=1;
                    CONN_x.Analyses(tianalysis).weight=2;
                    CONN_x.Analyses(tianalysis).type=1;
                    CONN_x.Analyses(tianalysis).regressors.names=CONN_x.dynAnalyses(CONN_x.dynAnalysis).regressors.names;
                    CONN_x.Analyses(tianalysis).regressors.dimensions=repmat({1},size(CONN_x.dynAnalyses(CONN_x.dynAnalysis).regressors.names));
                    CONN_x.Analyses(tianalysis).regressors.deriv=repmat({0},size(CONN_x.dynAnalyses(CONN_x.dynAnalysis).regressors.names));
                    CONN_x.Analyses(tianalysis).regressors.types=repmat({'roi'},size(CONN_x.dynAnalyses(CONN_x.dynAnalysis).regressors.names));
                    CONN_x.Analyses(tianalysis).regressors.fbands=repmat({1},size(CONN_x.dynAnalyses(CONN_x.dynAnalysis).regressors.names));
                    %CONN_x.Analysis_variables=CONN_x.Analyses(tianalysis).regressors;
                    CONN_x.Analyses(tianalysis).sourcenames={};
                    CONN_x.Analysis=tianalysis;
                    [ok,nill]=mkdir(CONN_x.folders.firstlevel,name);
                    if ispc, [ok,nill]=system(sprintf('del "%s"',fullfile(CONN_x.folders.firstlevel,name,'resultsROI*.*')));
                    else [ok,nill]=system(sprintf('rm ''%s''',fullfile(CONN_x.folders.firstlevel,name,'resultsROI*')));
                    end
                    if ispc, [ok,nill]=system(sprintf('del "%s"',fullfile(CONN_x.folders.firstlevel,name,'_list*.*')));
                    else [ok,nill]=system(sprintf('rm ''%s''',fullfile(CONN_x.folders.firstlevel,name,'_list*')));
                    end
                end
            end
        else
            filename=fullfile(filepathresults,'dyn_Base.mat');
            load(filename,'selectedconditionname','H'); %'X','Y','Xfilter','B','H','B0','H0','names','ROInames','IDX_subject','IDX_session','COND_names','COND_weights','selectedconditionname');
            Ncomponents=min(Ncomponents, size(H,2));
            if ~isempty(selectedconditionname),
                selectedcondition=find(strcmp(CONN_x.Setup.conditions.names(1:numel(CONN_x.Setup.conditions.names)-1),selectedconditionname));
                if numel(selectedcondition)~=1, error('Prior analyses have been run using a non-existing or deleted condition'); end
            else selectedcondition=[];
            end
            
            newanalyses=[];
            for ncomp=1:Ncomponents
                if isempty(selectedcondition)||strcmp(CONN_x.Setup.conditions.names{selectedcondition},'rest'), name=sprintf('Dynamic factor %02d',ncomp);
                else name=sprintf('Dynamic factor %s_%02d',CONN_x.Setup.conditions.names{selectedcondition},ncomp);
                end
                name=fullfile(CONN_x.dynAnalyses(CONN_x.dynAnalysis).name,name);
                names={CONN_x.Analyses.name};
                tianalysis=find(strcmp(names,name));
                if isempty(tianalysis), conn_disp('WARNING: Analysis %s not found. Skipping this analysis',name);
                else newanalyses=[newanalyses tianalysis];
                end
            end
        end
        
        if dosubjectlevel % subject-level backprojection
            if ~isequal(validsubjects,1:CONN_x.Setup.nsubjects)
                if ~(isfield(CONN_x,'pobj')&&isstruct(CONN_x.pobj)&&isfield(CONN_x.pobj,'subjects')), % not running in parallel: switch to all subjects
                    tempvalidsubjects=CONN_x.gui.subjects;
                    CONN_x.gui.subjects=1:CONN_x.Setup.nsubjects;
                    for n1=1:numel(newanalyses), conn_process('analyses_roi',newanalyses(n1)); end
                    CONN_x.gui.subjects=tempvalidsubjects;
                elseif isfield(CONN_x,'pobj')&&isstruct(CONN_x.pobj)&&isfield(CONN_x.pobj,'partition')&&isequal(CONN_x.pobj.partition,[1 1]) % running in parallel, single-job, not all subjects
                    tempvalidsubjects=CONN_x.pobj.subjects;
                    CONN_x.pobj.subjects=1:CONN_x.Setup.nsubjects;
                    for n1=1:numel(newanalyses), conn_process('analyses_roi',newanalyses(n1)); end
                    CONN_x.pobj.subjects=tempvalidsubjects;
                else % running in parallel, multiple-jobs (note: results may be incomplete if not all subjects were processed across the multiple jobs)
                    for n1=1:numel(newanalyses), conn_process('analyses_roi',newanalyses(n1)); end
                end
            else
                for n1=1:numel(newanalyses), conn_process('analyses_roi',newanalyses(n1)); end
            end
        end
        if ~isempty(newanalyses), CONN_x.Analysis=newanalyses(1); end
    end
    CONN_x.dynAnalysis=analysisbak;
    CONN_x.isready(4)=1;
    %conn_waitbar('close',h);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates swresultsROI_Subject###_Condition###.mat files (first-level sliding-window ROI-to-ROI analyses)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if any(options==13.5) && any(CONN_x.Setup.steps([1])) && ~(isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'steps')&&~any(CONN_x.gui.steps([1]))),
%     [path,name,ext]=fileparts(CONN_x.filename);
%     filepath=CONN_x.folders.preprocessing;
%     if nargin>1,analyses=varargin{1}; % selected analysis only
%     else, analyses=1:length(CONN_x.Analyses); end; 
%     doanalyses=false(1,numel(analyses));
%     for nanalyses=1:numel(analyses),if analyses(nanalyses)>0&&any(CONN_x.Analyses(analyses(nanalyses)).type==[1,3]), doanalyses(nanalyses)=true; end; end
%     analyses=analyses(doanalyses);
%     h=conn_waitbar(0,['Step ',num2str(sum(options<=13.5)),'/',num2str(length(options)),': ROI-to-ROI first-level analyses']);
%     analysisbak=CONN_x.Analysis;
%     nconditions=length(CONN_x.Setup.conditions.names)-1;
%     icondition=[];isnewcondition=[];for ncondition=1:nconditions,[icondition(ncondition),isnewcondition(ncondition)]=conn_conditionnames(CONN_x.Setup.conditions.names{ncondition}); end
%     if any(isnewcondition), error(['Some conditions have not been processed yet. Re-run previous step']); end
%     for nanalyses=1:length(analyses),
%         ianalysis=analyses(nanalyses);
%         CONN_x.Analysis=ianalysis;
% %         ianalysis=CONN_x.Analysis;
%         filepathresults=fullfile(CONN_x.folders.firstlevel,CONN_x.Analyses(ianalysis).name);
%         REDO='Yes';%filename=fullfile(filepathresults,['swresultsROI_Subject',num2str(1,'%03d'),'_Condition',num2str(icondition(validconditions(1)),'%03d'),'.mat']);
%         %if ~isempty(dir(filename)),if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'overwrite'), REDO=CONN_x.gui.overwrite; else, REDO=conn_questdlg('Overwrite existing subject results?','','Yes', 'No', 'Yes');end; end
%         if nanalyses==1, n=0; end
% 
%         CC=[];
%         XX=[];
%         for ncondition=1:nconditions,
%             for nsub=1:CONN_x.Setup.nsubjects,
%                 filename=fullfile(filepathresults,['swresultsROI_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'.mat']);
%                 if strcmp(lower(REDO),'yes')||isempty(dir(filename)),
%                     filename=fullfile(filepath,['ROI_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'.mat']);
%                     if isempty(dir(filename)), conn_disp(['Not ready to process step conn_process_10']); conn_waitbar('close',h);return; end
%                     X1=load(filename);
%                     if numel(X1.conditionweights)<5, 
%                         answ=conn_questdlg('Data has been processed using an old CONN version. Need to re-run Denoising.ROI step to be able to continue. This step should take only a few minutes. Do this now?','','Yes','No','Yes'); 
%                         if isequal(answ,'Yes')
%                             conn preprocessing_roi;
%                             X1=load(filename);
%                         else, conn_waitbar('close',h);return; 
%                         end
%                     end
%                     [X,nill,names,xyz]=conn_designmatrix(CONN_x.Analyses(ianalysis).regressors,X1,[]);
%                     nrois=size(X,2)-1;
%                     if nanalyses==1&&nsub==1 && ncondition==1,N=CONN_x.Setup.nsubjects*nconditions*nrois*length(analyses); nroisbak=nrois;
%                     elseif nrois~=nroisbak, error('Incorrect number of ROIs in subject %d condition %d (expected %d; found %d)',nsub,ncondition,nroisbak,nrois); 
%                     end
%                     wx=max(0,X1.conditionweights{1});
%                     
%                     x=X(:,2:end).*repmat(wx,1,size(x,2));
%                     x=x-(wx/sum(wx))*sum(x,1);
%                     x=bsxfun(@rdivide,x,sqrt(sum(x.^2,1)));
%                     
%                     %SL sliding window FWHM in samples
%                     SL=60/CONN_x.Setup.RT(min(numel(CONN_x.Setup.RT),nsub));
%                     wi=X1.conditionweights{4};
%                     ws=cumsum([1;diff(X1.conditionweights{5})~=1]);
%                     C=zeros(nrois*(nrois-1)/2,numel(wi));
%                     M=triu(true(nrois),1);
%                     for ni=1:numel(wi)
%                         wd=abs(wi-wi(ni));
%                         wl=(ws==ws(ni) & wd<=SL);
%                         wk=.5+cos(pi/SL*wd(wl))/2;
%                          t=x(wl,:).*repmat(wk,1,size(x,2));
%                          t=t-(wk/sum(wk))*sum(t,1);
%                         %t=x(wl,:);%.*repmat(wk,1,size(x,2));
%                         %t=t-repmat(mean(t,1),size(t,1),1);
%                         %%t=t./repmat(sqrt(sum(t.^2,1)),size(t,1),1);
%                         %%c=t'*t;%/numel(wk)/mean(wk.^2);
%                         c=t'*t/(1+sum(wk.^2));%/numel(wk)/mean(wk.^2);
%                         C(:,ni)=c(M);
%                         %imagesc(c);title(ni);colorbar; drawnow;
%                     end
%                     CC=[CC C];
%                     XX=[XX zeros(size(XX,1),size(C,2)); zeros(1,size(XX,2)) ones(1,size(C,2))];
%                     
%                     
% %                     [X2,nill,names2,xyz2]=conn_designmatrix({CONN_x.Analysis_variables,CONN_x.Analyses(ianalysis).regressors},X1,[]);
% %                     nrois2=size(X2,2)-1;
% %                     idxroi1roi2=zeros(1,nrois);
% %                     for n1=1:nrois,temp=strmatch(names{n1},names2,'exact'); idxroi1roi2(n1)=temp(1);end
% %                     X2=cat(2,X2(:,1),X2(:,1+idxroi1roi2),X2(:,1+setdiff(1:nrois2,idxroi1roi2)));
% %                     names2=cat(2,{names2{idxroi1roi2}},{names2{setdiff(1:nrois2,idxroi1roi2)}});
% %                     xyz=cat(2,{xyz2{idxroi1roi2}},{xyz2{setdiff(1:nrois2,idxroi1roi2)}});
% % %                     for n1=1:nrois,temp=strmatch(names{n1},names2,'exact'); if isempty(temp), idxroi1roi2(n1)=nrois2+n1; else, idxroi1roi2(n1)=temp(1);end; end
% % %                     X2t=[X2(:,2:end) X(:,2:end)]; names2t=[names2 names]; xyz2t=[xyz2 xyz];
% % %                     X2=cat(2,X2(:,1),X2t(:,idxroi1roi2),X2(:,1+setdiff(1:nrois2,idxroi1roi2)));
% % %                     names2=cat(2,{names2t{idxroi1roi2}},{names2{setdiff(1:nrois2,idxroi1roi2)}});
% % %                     xyz=cat(2,{xyz2t{idxroi1roi2}},{xyz2{setdiff(1:nrois2,idxroi1roi2)}});
% %                     
% %                     if nanalyses==1&&nsub==1 && ncondition==1,N=CONN_x.Setup.nsubjects*nconditions*nrois2*length(analyses); nrois2bak=nrois2;
% %                     elseif nsub==1 && ncondition==1, N=N-CONN_x.Setup.nsubjects*nconditions*nrois2bak+CONN_x.Setup.nsubjects*nconditions*nrois2; end
% %                     Z=zeros(nrois,nrois2);%+diag(nan+zeros(nrois,1));
% %                     SE=zeros(1,nrois2);
% %                     wx=ones(size(X,1),1);
% %                     switch(CONN_x.Analyses(ianalysis).weight),
% %                         case 1, if numel(X1.conditionweights)>2, wx=X1.conditionweights{3}; end
% %                         case 2, wx=X1.conditionweights{1};
% %                         case 3, wx=X1.conditionweights{2};
% %                     end
% %                     X=cat(2,X(:,1),conn_wdemean(X(:,2:end),wx));
% %                     X=X.*repmat(wx,[1,size(X,2)]);
% %                     X2=cat(2,X2(:,1),conn_wdemean(X2(:,2:end),wx));
% %                     X2=X2.*repmat(wx,[1,size(X2,2)]);                   
% %                     switch(CONN_x.Analyses(ianalysis).measure),
% %                         case {1,3}, %bivariate
% %                             DOF=max(0,size(X,1)*(min(1/(2*CONN_x.Setup.RT),CONN_x.Preproc.filter(2))-max(0,CONN_x.Preproc.filter(1)))/(1/(2*CONN_x.Setup.RT))-1);
% %                         case {2,4}, %partial
% %                             DOF=max(0,size(X,1)*(min(1/(2*CONN_x.Setup.RT),CONN_x.Preproc.filter(2))-max(0,CONN_x.Preproc.filter(1)))/(1/(2*CONN_x.Setup.RT))-rank(X)+1);
% %                     end
% %                     for nroi=1:nrois2,
% %                         y=X2(:,1+nroi);
% %                         switch(CONN_x.Analyses(ianalysis).measure),
% %                             case {1,3}, %bivariate
% %                                 x=cat(2,X(:,1),X(:,1+setdiff(1:nrois,nroi)));
% %                                 iX=diag(1./max(eps,sum(x.^2,1)));
% %                                 %iX=pinv(diag(diag(x'*x)));
% %                             case {2,4}, %partial
% %                                 x=cat(2,X(:,1),X(:,1+setdiff(1:nrois,nroi)));
% %                                 iX=pinv(x'*x);
% %                         end
% %                         B=iX*(x'*y);
% %                         e=sqrt(sum(abs(y).^2,1));
% %                         switch(CONN_x.Analyses(ianalysis).measure),
% %                             case {1,2}, %correlation
% %                                 r=sqrt(diag(iX));
% %                                 B=B./max(eps,r*e);
% %                                 B(~isnan(B))=atanh(max(eps-1,min(1-eps,B(~isnan(B)))));
% %                                 SE(nroi)=1./max(eps,sqrt(DOF-3));
% %                             case {3,4},
% %                                 SE(nroi)=e/max(eps,DOF);
% %                         end
% %                         Z(setdiff(1:nrois,nroi),nroi)=B(2:end);
% %                         if nroi<=nrois,Z(nroi,nroi)=nan;end
% %                         n=n+1;
% %                         conn_waitbar(n/N,h,sprintf('Subject %d Condition %d',nsub,ncondition));
% %                     end
% % %                     xyz={}; %note: this assumes constant number of dimensions per subject for analysis regressors
% % %                     for n1=1:length(CONN_x.Analysis_variables.names),
% % %                         for n2=1:CONN_x.Analysis_variables.deriv{n1}+1,
% % %                             for n3=1:CONN_x.Analysis_variables.dimensions{n1}(1),
% % %                                 idx=strmatch(CONN_x.Analysis_variables.names{n1},X1.names,'exact');
% % %                                 if isempty(idx), xyz{end+1}=''; else, xyz{end+1}=X1.xyz{idx}; end
% % %                             end
% % %                         end;
% % %                     end
% % %                     xyz=cat(2,{xyz{idxroi1roi2}},{xyz{setdiff(1:nrois2,idxroi1roi2)}});
% %                     regressors=CONN_x.Analyses(ianalysis).regressors;
% %                     filename=fullfile(filepathresults,['swresultsROI_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'.mat']);
% %                     save(filename,'Z','regressors','names','names2','xyz','SE','DOF');
%                 end
%             end
%         end
%     end
%     CONN_x.Analysis=analysisbak;
%     conn_waitbar('close',h);
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates resultsDATA_Condition###_Source###.mat files (combined first-level analyses)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if 0,%any(options==13) && any(CONN_x.Setup.steps([2])),
%     if ~isfield(CONN_x.Setup,'normalized'), CONN_x.Setup.normalized=1; end
%     if ~CONN_x.Setup.normalized,
%         conn_disp(['Not spatially-normalized data. Skipping ROI-to-voxel second-level analyses']);
%     else,
%         [path,name,ext]=fileparts(CONN_x.filename);
%         filepath=CONN_x.folders.preprocessing;
%         if nargin>1,analyses=varargin{1}; % selected analysis only
%         else, analyses=1:length(CONN_x.Analyses); end;
%         h=conn_waitbar(0,['Step ',num2str(sum(options<=13)),'/',num2str(length(options)),': Preparing second-level analyses']);
%         analysisbak=CONN_x.Analysis;
%         for nanalyses=1:length(analyses),
%             ianalysis=analyses(nanalyses);
%             CONN_x.Analysis=ianalysis;
%             %         ianalysis=CONN_x.Analysis;
%             filepathresults=fullfile(CONN_x.folders.firstlevel,CONN_x.Analyses(ianalysis).name);
%             nconditions=length(CONN_x.Setup.conditions.names)-1;
%             
%             filename=fullfile(filepath,['ROI_Subject',num2str(1,'%03d'),'_Condition',num2str(1,'%03d'),'.mat']);
%             if isempty(dir(filename)), conn_disp(['Not ready to process step conn_process_12']); conn_waitbar('close',h);return; end
%             X1=load(filename);
%             [X,nill,names]=conn_designmatrix(CONN_x.Analyses(ianalysis).regressors,X1,[]);
%             nrois=size(X,2)-1;
%             iroi=[];isnew=[];for nroi=1:nrois,[iroi(nroi),isnew(nroi)]=conn_sourcenames(names{nroi},'-');end
%             if any(isnew), error(['Non-existing ROI first-level data for ',names{find(isnew)},'. Please repeat first-level analyses']); return; end
%             if nanalyses==1, N=nconditions*nrois*length(analyses);n=0; nroisbak=nrois; 
%             else N=N-nconditions*nroisbak+nconditions*nrois; end
%             
%             for ncondition=1:nconditions,
%                 for nroi=1:nrois,
%                     clear Yin Yout;
%                     for nsub=1:CONN_x.Setup.nsubjects,
%                         filename=fullfile(filepathresults,['resultsDATA_Subject',num2str(nsub,'%03d'),'_Condition',num2str(ncondition,'%03d'),'_Source',num2str(iroi(nroi),'%03d'),'.mat']);
%                         Yin(nsub)=conn_vol(filename);
%                     end
%                     filename=fullfile(filepathresults,['resultsDATA_Condition',num2str(ncondition,'%03d'),'_Source',num2str(iroi(nroi),'%03d'),'.mat']);
%                     Yout=Yin(1); Yout.fname=filename;
% %                     Yout.size.Nt=CONN_x.Setup.nsubjects;
% %                     Yout=conn_init_vol(Yout);
%                     conn_write_combine(Yin,Yout);
%                     n=n+1;
%                     conn_waitbar(n/N,h);
%                 end
%                 clear Yin Yout;
%                 for nsub=1:CONN_x.Setup.nsubjects,
%                     filename=fullfile(filepathresults,['seDATA_Subject',num2str(nsub,'%03d'),'_Condition',num2str(ncondition,'%03d'),'.mat']);
%                     Yin(nsub)=conn_vol(filename);
%                 end
%                 filename=fullfile(filepathresults,['seDATA_Condition',num2str(ncondition,'%03d'),'.mat']);
%                 Yout=Yin(1); Yout.fname=filename;
%                 Yout.DOF=cat(2,Yin(:).DOF);
% %                 Yout.size.Nt=CONN_x.Setup.nsubjects;
% %                 Yout=conn_init_vol(Yout);
%                 conn_write_combine(Yin,Yout);
%             end
%             CONN_x.Analyses(ianalysis).sources=names;
%         end
%         CONN_x.Analysis=analysisbak;
%         conn_waitbar('close',h);
%         %CONN_x.Results.measure=CONN_x.Analyses(ianalysis).measure;
%     end
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates resultsROI_Condition###.mat files (combined first-level analyses)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(options==15) && any(CONN_x.Setup.steps([1])) && ~(isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'steps')&&~any(CONN_x.gui.steps([1]))),
    [path,name,ext]=fileparts(CONN_x.filename);
    filepath=CONN_x.folders.preprocessing;
    if nargin>1&&~isempty(varargin{1}),analyses=varargin{1}; % selected analysis only
    else analyses=1:length(CONN_x.Analyses); 
    end
    if ischar(analyses)||iscell(analyses), analyses=find(ismember({CONN_x.Analyses.name},analyses)); end
    doanalyses=false(1,numel(analyses));
    for nanalyses=1:numel(analyses),if analyses(nanalyses)>0&&any(CONN_x.Analyses(analyses(nanalyses)).type==[1,3]), doanalyses(nanalyses)=true; end; end
    analyses=analyses(doanalyses);
    validsubjects=1:CONN_x.Setup.nsubjects; %if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'subjects'), validsubjects=CONN_x.gui.subjects; else validsubjects=1:CONN_x.Setup.nsubjects; end
    if isfield(CONN_x,'pobj')&&isstruct(CONN_x.pobj)&&isfield(CONN_x.pobj,'subjects'), validsubjects=CONN_x.pobj.subjects; if ~isempty(analyses), conn_projectmanager('addstep',15,analyses); end; end
    DOREDUCED=true; % (compute RRC square matrices only) set to false for back-compatibility
    if isequal(validsubjects,1:CONN_x.Setup.nsubjects), 
        h=conn_waitbar(0,['Step ',num2str(sum(options<=15)),'/',num2str(length(options)),': Preparing second-level ROI analyses']);
        analysisbak=CONN_x.Analysis;
        nconditions=length(CONN_x.Setup.conditions.names)-1;
        if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'conditions')&&~isempty(CONN_x.gui.conditions), validconditions=CONN_x.gui.conditions; else validconditions=1:length(CONN_x.Setup.conditions.names)-1; end
        icondition=[];isnewcondition=[];for ncondition=1:nconditions,[icondition(ncondition),isnewcondition(ncondition)]=conn_conditionnames(CONN_x.Setup.conditions.names{ncondition}); end
        if any(isnewcondition(validconditions)), error(['Some conditions have not been processed yet. Re-run previous step']); end        
        secondaryconditions=validconditions(cellfun('length',CONN_x.Setup.conditions.model(validconditions))>0);
        validconditions=validconditions(cellfun('length',CONN_x.Setup.conditions.model(validconditions))==0);
        if ~isempty(validconditions), referenceconditions=validconditions;
        else
            evaluatefunction=CONN_x.Setup.conditions.model{secondaryconditions(1)}{1};
            if isequal(evaluatefunction,'lin'), primaryconditionsnames=CONN_x.Setup.conditions.model{secondaryconditions(1)}(3:end);
            else primaryconditionsnames=CONN_x.Setup.conditions.model{secondaryconditions(1)}(2:end);
            end
            referenceconditions=find(ismember(CONN_x.Setup.conditions.names(1:end-1),primaryconditionsnames));
        end        
        if any(isnewcondition(referenceconditions)), error(['Some conditions have not been processed yet. Re-run previous step']); end        
        for nanalyses=1:length(analyses),
            ianalysis=analyses(nanalyses);
            CONN_x.Analysis=ianalysis;
            names={};
            %         ianalysis=CONN_x.Analysis;
            filepathresults=fullfile(CONN_x.folders.firstlevel,CONN_x.Analyses(ianalysis).name);
            missingdata=arrayfun(@(n)isempty(dir(fullfile(filepath,['ROI_Subject',num2str(1,'%03d'),'_Condition',num2str(icondition(n),'%03d'),'.mat']))),validconditions);
            if any(missingdata), conn_disp(['Not ready to process step conn_process_15']); return; end
            filename=fullfile(filepath,['ROI_Subject',num2str(1,'%03d'),'_Condition',num2str(icondition(referenceconditions(1)),'%03d'),'.mat']);
            X1=load(filename);
            [X,nill,names,xyz]=conn_designmatrix(CONN_x.Analyses(ianalysis).regressors,X1,[]);
            nrois=size(X,2)-1;
            if DOREDUCED&CONN_x.Analyses(ianalysis).type==1
                X2=X; names2=names; xyz2=xyz;
            else
                [X2,nill,names2,xyz2]=conn_designmatrix({CONN_x.Analysis_variables,CONN_x.Analyses(ianalysis).regressors},X1,[]);
            end
            nrois2=size(X2,2)-1;
            [nill,idxroi1roi2]=ismember(names,names2);
            %idxroi1roi2=zeros(1,nrois);
            %for n1=1:nrois,temp=strmatch(names{n1},names2,'exact'); idxroi1roi2(n1)=temp(1);end
            X2=cat(2,X2(:,1),X2(:,1+idxroi1roi2),X2(:,1+setdiff(1:nrois2,idxroi1roi2)));
            names2=cat(2,{names2{idxroi1roi2}},{names2{setdiff(1:nrois2,idxroi1roi2)}});
            xyz=cat(2,{xyz2{idxroi1roi2}},{xyz2{setdiff(1:nrois2,idxroi1roi2)}});
            if nanalyses==1, n=0;N=CONN_x.Setup.nsubjects*(numel(validconditions)+numel(secondaryconditions))*nrois*length(analyses); nroisbak=nrois;
            else N=N-CONN_x.Setup.nsubjects*(numel(validconditions)+numel(secondaryconditions))*nroisbak+CONN_x.Setup.nsubjects*(numel(validconditions)+numel(secondaryconditions))*nrois; end
            
            names2checked=false;
            for ncondition=[validconditions,secondaryconditions],
                Z=zeros([nrois,nrois2,CONN_x.Setup.nsubjects]);
                xyz={};
                SE=zeros([CONN_x.Setup.nsubjects,nrois2]);
                DOF=zeros([1,CONN_x.Setup.nsubjects]);
                for nsub=1:CONN_x.Setup.nsubjects,
                    filename=fullfile(filepathresults,['resultsROI_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'.mat']);
                    t=load(filename,'Z','regressors','names','names2','xyz','SE','DOF');
                    IDX=[];IDX2=[];
                    for nroi=1:nrois,
                        idx=find(strcmp(names{nroi},t.names));
                        %idx=strmatch(names{nroi},t.names,'exact');
                        if isempty(idx), error(['Non-existing source ROI first-level data for ',names{nroi},' subject ',num2str(nsub),'. Please repeat first-level analyses']); return; end
                        IDX(nroi)=idx(1);
                        n=n+1;
                    end
                    for nroi=1:nrois2,
                        idx=find(strcmp(names2{nroi},t.names2));
                        %idx=strmatch(names2{nroi},t.names2,'exact');
                        %if isempty(idx)&&names2checked, error(['Non-existing target ROI first-level data for ',names2{nroi},' subject ',num2str(nsub),'. Please repeat first-level analyses']); return; 
                        if isempty(idx)&&~isempty(regexp(names2{nroi},'^QC_')), % skip warning
                        elseif isempty(idx)&&names2checked, %conn_disp(['Warning: Non-existing target ROI first-level data for ',names2{nroi},' subject ',num2str(nsub),'. Skipping this target ROI']); 
                        elseif ~isempty(idx), IDX2(nroi)=idx(1);
                        end
                    end
                    if nnz(IDX2)~=nrois2
                        if names2checked
                            Z=Z(:,IDX2>0,:);
                            SE=SE(:,IDX2>0);
                            names2=names2(IDX2>0);
                            nrois2=nnz(IDX2);
                            IDX2=IDX2(IDX2>0);
                        else
                            names2=names2(IDX2>0);
                            nrois2=nnz(IDX2);
                            IDX2=IDX2(IDX2>0);
                            Z=zeros([nrois,nrois2,CONN_x.Setup.nsubjects]);
                            SE=zeros([CONN_x.Setup.nsubjects,nrois2]);
                        end
                    end
                    names2checked=true;
                    Z(:,:,nsub)=t.Z(IDX,IDX2);
                    xyz={t.xyz{IDX2}};
                    SE(nsub,:)=t.SE(IDX2);
                    DOF(nsub)=t.DOF;
                    conn_waitbar(n/N,h,sprintf('Subject %d Condition %d',nsub,ncondition));
                end
                regressors=CONN_x.Analyses(ianalysis).regressors;
                filename=fullfile(filepathresults,['resultsROI_Condition',num2str(icondition(ncondition),'%03d'),'.mat']);
                save(filename,'Z','regressors','names','names2','xyz','SE','DOF');
            end
            if ~isempty(names), CONN_x.Analyses(ianalysis).sources=names; end
        end
        CONN_x.Analysis=analysisbak;
        conn_waitbar('close',h);
    end
    CONN_x.isready(4)=1;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates second-level SPM.mat results (second-level analysis results for SPM)
% CON_Subject### (source files for second-leve analyses)
% CON_Subject###_Source###.nii (source files for source-specific second-level analyses)
% between-subject contrasts of interest:
%    T-contrast: "connectivity results": as specified in "between subject contrast"
%    T-contrast: one per effect specified in "subject effects" (named as the corresponding second-level covariates)
%    F-contrast: "effects of interest": F- test for all effects in second-levle model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(options==16) && any(CONN_x.Setup.steps([2,3])) && ~(isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'steps')&&~any(CONN_x.gui.steps([2,3]))),
    FORCEORTH=true; % enforces orthogonal rows in within-subject contrast
    if nargin>1&&~isempty(varargin{1}), 
        if ischar(varargin{1})
            switch(lower(varargin{1})),
                case 'dosingle',dosinglecontrast=1;    % performs one specific analysis
                case 'doall',dosinglecontrast=0;       % one analysis per source
                case 'readsingle', dosinglecontrast=2; % performs one speficic analysis (if it does not exist already)
                otherwise, dosinglecontrast=0;
            end
            allsources=[];
        else
            dosinglecontrast=0;
            allsources=varargin{1};
        end
    else, dosinglecontrast=0; allsources=[]; end
    state=2;
    if nargin>2&&~isempty(varargin{2}), 
        switch(lower(varargin{2})),
            case 'seed-to-voxel', state=2;
            case 'voxel-to-voxel', state=3;
        end
    end
    if state==2,
        if nargin>3&&~isempty(varargin{3}),CONN_x.Analysis=varargin{3}; end % selected analysis only
        if ischar(CONN_x.Analysis), CONN_x.Analysis=find(strcmp({CONN_x.Analyses.name},CONN_x.Analysis)); end
    else
        if nargin>3&&~isempty(varargin{3}),CONN_x.vvAnalysis=varargin{3}; end % selected analysis only
        if ischar(CONN_x.vvAnalysis), CONN_x.vvAnalysis=find(strcmp({CONN_x.vvAnalyses.name},CONN_x.vvAnalysis)); end
    end
    if nargin>4&&~isempty(varargin{4}), % selected contrasts only
        ncontrast=varargin{4}; 
        if ischar(ncontrast), ncontrast=find(strcmp(CONN_x.Results.saved.names,ncontrast)); end
        assert(numel(ncontrast)==1,'incorrect number of contrasts');
        newneffects=CONN_x.Results.saved.nsubjecteffects{ncontrast};
        newceffects=CONN_x.Results.saved.csubjecteffects{ncontrast};
        newnconditions=CONN_x.Results.saved.nconditions{ncontrast};
        newcconditions=CONN_x.Results.saved.cconditions{ncontrast};
        [ok1,i1]=ismember(newneffects,CONN_x.Setup.l2covariates.names(1:end-1));
        [i1,idx1]=sort(i1);
        CONN_x.Results.xX.nsubjecteffects=i1;
        CONN_x.Results.xX.csubjecteffects=newceffects(:,idx1);
        CONN_x.Results.xX.nsubjecteffectsbyname=CONN_x.Setup.l2covariates.names(i1);
        [ok2,i2]=ismember(newnconditions,CONN_x.Setup.conditions.names(1:end-1));
        [i2,idx2]=sort(i2);
        CONN_x.Results.xX.nconditions=i2;
        CONN_x.Results.xX.cconditions=newcconditions(:,idx2);
        CONN_x.Results.xX.nconditionsbyname=CONN_x.Setup.conditions.names(i2);
    end 
    if ~isfield(CONN_x,'Results'), conn_disp(['Not ready to process step conn_process_16']); return; end
    if state==2
        filepathresults1=fullfile(CONN_x.folders.firstlevel,CONN_x.Analyses(CONN_x.Analysis).name);
    else
        filepathresults1=fullfile(CONN_x.folders.firstlevel_vv,CONN_x.vvAnalyses(CONN_x.vvAnalysis).name);
    end
    filepathresults2=CONN_x.folders.secondlevel;
    if state==2
        sources=CONN_x.Analyses(CONN_x.Analysis).sources; %CONN_x.Results.xX.sources;
        if dosinglecontrast==0
            if isempty(allsources), nsources=1:numel(sources);
            else nsources=allsources;
            end
            csources=ones(1,numel(nsources));
            CONN_x.Results.xX.nsources=nsources;
            CONN_x.Results.xX.csources=csources;
            CONN_x.Results.xX.nsourcesbyname=sources(nsources);
        else
            nsources=CONN_x.Results.xX.nsources;
            csources=CONN_x.Results.xX.csources;
        end
    else
        sources=CONN_x.vvAnalyses(CONN_x.vvAnalysis).measures; %CONN_x.Results.xX.measures;
        if dosinglecontrast==0
            if isempty(allsources), nsources=1:numel(sources);
            else nsources=allsources;
            end
            csources=ones(1,numel(nsources));
            CONN_x.Results.xX.nmeasures=nsources;
            CONN_x.Results.xX.cmeasures=csources;
            CONN_x.Results.xX.nmeasuresbyname=conn_v2v('cleartext',sources(nsources));
        else
            nsources=CONN_x.Results.xX.nmeasures;
            csources=CONN_x.Results.xX.cmeasures;
        end
    end
    nconditions=CONN_x.Results.xX.nconditions;
    icondition=[];isnewcondition=[];for ncondition=nconditions(:)',[icondition(ncondition),isnewcondition(ncondition)]=conn_conditionnames(CONN_x.Setup.conditions.names{ncondition}); end
    if any(isnewcondition), error(['Some conditions have not been processed yet. Re-run previous step']); end
    cconditions=CONN_x.Results.xX.cconditions;
    nsubjecteffects=CONN_x.Results.xX.nsubjecteffects;
    csubjecteffects=CONN_x.Results.xX.csubjecteffects;
    if ~isfield(CONN_x.Results.xX,'modeltype')||isempty(CONN_x.Results.xX.modeltype), CONN_x.Results.xX.modeltype=1; end; %note: obsolete 'modeltype=2' fixed-effect analyses option; remove in next release
    modeltype=CONN_x.Results.xX.modeltype;
%     if nargin>1&&varargin{1}>1, nsources=varargin{2}; csources=ones(1,length(nsources)); sources={sources{nsources}}; nsources=1:length(nsources); end
%     if nargin>1, dosinglecontrast=(varargin{1}==1); else, dosinglecontrast=length(nsources)~=1; end; %~(length(nsources)>1 & all(csources==1/length(nsources))); end
    if isfield(CONN_x.Results,'foldername')&&~isempty(CONN_x.Results.foldername),
        [ok,nill]=mkdir(filepathresults2,CONN_x.Results.foldername);
        filepathresults2=fullfile(filepathresults2,CONN_x.Results.foldername);
        CONN_x.Results.foldername=[];
    else,
        [foldername,foldername_back]=conn_resultsfolder('subjectsconditions',state,nsubjecteffects,csubjecteffects,nconditions,cconditions);
        for nfolderbak=1:numel(foldername_back), 
            if isdir(fullfile(filepathresults2,foldername_back{nfolderbak})), foldername=foldername_back{nfolderbak}; break; end % backwards-compatibility with existing results
        end
        [ok,nill]=mkdir(filepathresults2,foldername);
        if ok,filepathresults2=fullfile(filepathresults2,foldername);
        else,filepathresults2=uigetdir(filepathresults2,'Select a directory to write the results');end
    end
    if ~ischar(filepathresults2), return; end;
%     conn_disp(['Second-level results stored in ',filepathresults2]);

	MDok = conn_checkmissingdata(state);
    X=zeros(CONN_x.Setup.nsubjects,length(CONN_x.Setup.l2covariates.names)-1);
    for nsub=1:CONN_x.Setup.nsubjects,
        for ncovariate=1:length(CONN_x.Setup.l2covariates.names)-1;
            X(nsub,ncovariate)=CONN_x.Setup.l2covariates.values{nsub}{ncovariate};
        end
    end
    nsubjects=find(any(X(:,nsubjecteffects)~=0,2)&~any(isnan(X(:,nsubjecteffects)),2)&MDok);
    if isempty(nsubjects), error('Null design matrix (all subjects have missing data or all selected covariates have zero or NaN values). Please check your second-level covariates'); end
    
    clear SPMall;
    cwd=pwd;
    cd(filepathresults2);
%     REDUCENAMES=(state==2);
    if dosinglecontrast>0,
        [foldername,foldername_back]=conn_resultsfolder('sources',state,sources,nsources,csources);
        for nfolderbak=1:numel(foldername_back), 
            if isdir(fullfile(filepathresults2,foldername_back{nfolderbak})), foldername=foldername_back{nfolderbak}; break; end % backwards-compatibility with existing results
        end
        [ok,nill]=mkdir(filepathresults2,foldername);
        if ok,filepathresults3{1}=fullfile(filepathresults2,foldername);
        else,filepathresults3{1}=uigetdir(filepathresults2,'Select a directory to write the results');dosinglecontrast=1;end
        if dosinglecontrast==1&&isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'overwrite')&&strcmpi(CONN_x.gui.overwrite,'no')&&~isempty(dir(fullfile(filepathresults3{1},'SPM.mat'))), dosinglecontrast=2;
        elseif dosinglecontrast==2,
            if isempty(dir(fullfile(filepathresults3{1},'SPM.mat'))), dosinglecontrast=1;
            elseif ~(isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'overwrite')), REDO=conn_questdlg('','results explorer','Load existing analysis results', 'Recompute/overwrite results', 'Load existing analysis results'); if strcmp(lower(REDO),'recompute/overwrite results'),dosinglecontrast=1; end; end
            %elseif ~(isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'overwrite')), REDO=conn_questdlg('Re-estimate/Overwrite existing second-level results?','','Yes', 'No', 'No'); if strcmp(lower(REDO),'yes'),dosinglecontrast=1; end; end
        end
        if FORCEORTH&&numel(nconditions)>1&&~isequal(cconditions,eye(numel(nconditions))),cconditions=spm_orth(cconditions','norm')';end
        if FORCEORTH&&numel(nsources)>1&&~isequal(csources,eye(numel(nsources))),csources=spm_orth(csources','norm')';end
        if dosinglecontrast==1,
            h=conn_waitbar(0,['Step ',num2str(sum(options<=16)),'/',num2str(length(options)),': Functional data second-level analyses']);
            clear SE Y;
            altestsmooth=false;
            Zfiles={};
            Znames={};
            for nsub=1:length(nsubjects),%CONN_x.Setup.nsubjects,
                filename={};
                contrastname={};
                for n0=1:numel(nconditions)
                    for n1=1:numel(nsources)
                        nroi=nsources(n1);
                        roiname=sources{nroi};
                        ncondition=nconditions(n0);
                        % note: check here missing data
                        if state==2
                            [iroi,isnew]=conn_sourcenames(roiname,'-');
                            if isnew, error(['Non-existing seed-to-voxel first-level data for ',roiname,' subject ',num2str(nsubjects(nsub)),'. Please repeat first-level analyses']); return; end
                            tfilename=fullfile(filepathresults1,['BETA_Subject',num2str(nsubjects(nsub),CONN_x.opt.fmt1),'_Condition',num2str(icondition(ncondition),'%03d'),'_Source',num2str(iroi,'%03d'),'.nii']);
                            troiname=regexprep(roiname,{'^BA\.(\d+) \(([LR])\)\. .*','^\((-?\d+),(-?\d+),(-?\d+)\)$','^SLrois\.|^aal\.|^atlas\.|^networks\.','\s\(([LR])\)','([^\(\)]+)\(.+\)\s*$'},{'$1$2','($1 $2 $3)','',' ${lower($1)}','$1'});
                            if numel(nconditions)==1&&numel(nsources)>1, tcontrastname=troiname;
                            elseif numel(nconditions)>1&&numel(nsources)==1, tcontrastname=CONN_x.Setup.conditions.names{ncondition};
                            else tcontrastname=[CONN_x.Setup.conditions.names{ncondition},'_',troiname];
                            end
                        else
                            altestsmooth=altestsmooth|ismember(conn_v2v('fieldtext',roiname,1),{'2'});
                            %altestsmooth=altestsmooth|~isempty(strmatch('connectome',roiname)); 
                            [iroi,isnew,ncomp]=conn_v2v('match_extended',roiname);
                            %[iroi,isnew]=conn_v2v('match_measures',CONN_x.vvAnalyses.regressors,nroi,'-');
                            if isnew, error(['Non-existing voxel-to-voxel first-level data for ',roiname,' subject ',num2str(nsubjects(nsub)),'. Please repeat first-level analyses']); return; end
                            tfilename=fullfile(filepathresults1,['BETA_Subject',num2str(nsubjects(nsub),'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'_Measure',num2str(iroi,'%03d'),'_Component',num2str(ncomp,'%03d'),'.nii']);
                            if numel(nconditions)==1&&numel(nsources)>1, tcontrastname=conn_v2v('cleartext',roiname);
                            elseif numel(nconditions)>1&&numel(nsources)==1, tcontrastname=CONN_x.Setup.conditions.names{ncondition};
                            else tcontrastname=[CONN_x.Setup.conditions.names{ncondition},'_',conn_v2v('cleartext',roiname)];
                            end
                            
                        end
                        filename{n1,n0}=tfilename;
                        contrastname{n1,n0}=tcontrastname;
                    end
                end
                Zfiles(nsub,:)=filename(:)';
                Znames=contrastname(:)';
                Zcontr=kron(cconditions,csources);
                if numel(nconditions)>1&&~isequal(cconditions,eye(numel(nconditions))),%&&size(cconditions,1)==1
                    for n2=1:size(cconditions,1)
                        for n1=1:numel(nsources)
                            tcontrastname='';
                            b=0;
                            for n0=1:numel(nconditions)
                                a=spm_vol(filename{n1,n0});
                                b=b+spm_read_vols(a)*cconditions(n2,n0);
                                if cconditions(n2,n0)==1, tcontrastname=[tcontrastname,'+',contrastname{n1,n0}];
                                elseif cconditions(n2,n0)==-1, tcontrastname=[tcontrastname,'-',contrastname{n1,n0}];
                                elseif cconditions(n2,n0)~=0, tcontrastname=[tcontrastname,num2str(cconditions(n2,n0),'%+g'),'*',contrastname{n1,n0}];
                                end
                                %                             tcontrastname=[tcontrastname,char('+'+2*(sign(cconditions(n0))<0)),num2str(cconditions(n0)),contrastname{n1,n0}];
                            end
                            tfilename=fullfile(filepathresults3{1},['CON_Subject',num2str(nsubjects(nsub),'%03d'),'_compA',num2str(n2,'%03d'),'_',num2str(n1,'%03d'),'.nii']);
                            a.fname=tfilename;
                            spm_write_vol(a,b);
                            filename{n1,numel(nconditions)+n2}=tfilename;
                            contrastname{n1,numel(nconditions)+n2}=tcontrastname;
                        end
                    end
                    filename=filename(:,end-size(cconditions,1)+(1:size(cconditions,1)));
                    contrastname=contrastname(:,end-size(cconditions,1)+(1:size(cconditions,1)));
                    [nill,new_cconditions]=conn_mtxbase(cconditions);
                else
                    new_cconditions=cconditions;
                end
                if numel(nsources)>1&&~isequal(csources,eye(numel(nsources))),%&&size(csources,1)==1
                    for n2=1:size(csources,1)
                        for n0=1:size(filename,2)
                            tcontrastname='';
                            b=0;
                            for n1=1:numel(nsources)
                                a=spm_vol(filename{n1,n0});
                                b=b+spm_read_vols(a)*csources(n2,n1);
                                if csources(n2,n1)==1, tcontrastname=[tcontrastname,'+(',contrastname{n1,n0},')'];
                                elseif csources(n2,n1)==-1, tcontrastname=[tcontrastname,'-(',contrastname{n1,n0},')'];
                                elseif csources(n2,n1)~=0, tcontrastname=[tcontrastname,num2str(csources(n2,n1)),'*(',contrastname{n1,n0},')'];
                                end
                                %                             tcontrastname=[tcontrastname,char('+'+2*(sign(csources(n1))<0)),num2str(csources(n1)),contrastname{n1,n0}];
                            end
                            tfilename=fullfile(filepathresults3{1},['CON_Subject',num2str(nsubjects(nsub),'%03d'),'_compB',num2str(n2,'%03d'),'_',num2str(n0,'%03d'),'.nii']);
                            a.fname=tfilename;
                            spm_write_vol(a,b);
                            filename{numel(nsources)+n2,n0}=tfilename;
                            contrastname{numel(nsources)+n2,n0}=tcontrastname;
                        end
                    end
                    filename=filename(end-size(csources,1)+(1:size(csources,1)),:);
                    contrastname=contrastname(end-size(csources,1)+(1:size(csources,1)),:);
                    [nill,new_csources]=conn_mtxbase(csources);
                else
                    new_csources=csources;
                end
                filename=filename(:);
                contrastname=contrastname(:);
                contrast=kron(new_cconditions,new_csources);
                
                SPMall(1).xY.VY(nsub,:)=spm_vol(char(filename))';
                if nsub==1
                    SPMall(1).altestsmooth=altestsmooth;
                    SPMall(1).xX.name={};
                    for n01=1:numel(contrastname),for n00=1:numel(nsubjecteffects),SPMall(1).xX.name{n00,n01}=[CONN_x.Setup.l2covariates.names{nsubjecteffects(n00)},'_',contrastname{n01}]; end; end
                    SPMall(1).xX.name=SPMall(1).xX.name(:)';
                end
                for n1=1:numel(filename),SPMall(1).xY.VY(nsub,n1).fname=filename{n1}; end
                if modeltype==2&&nsub==1, SPMall(1).connvols.Y=Y;SPMall(1).connvols.SE=SE;SPMall(1).connvols.cconditions=cconditions;SPMall(1).connvols.csources=csources;end
                conn_waitbar(2/4*(nsub/length(nsubjects)),h,sprintf('Subject %d',nsub));
            end
            nrepeated=size(SPMall(1).xY.VY,2);
            SPMall(1).xX_multivariate.X=X(nsubjects,nsubjecteffects); 
            SPMall(1).xX_multivariate.C=csubjecteffects;
            SPMall(1).xX_multivariate.M=contrast;
            SPMall(1).xX_multivariate.Xnames=CONN_x.Setup.l2covariates.names(nsubjecteffects);
            SPMall(1).xX_multivariate.Ynames=contrastname;
            SPMall(1).xX_multivariate.Znames=Znames;
            SPMall(1).xX_multivariate.Zcontr=Zcontr;
            SPMall(1).xX_multivariate.Zfiles=Zfiles;
            SPMall(1).xX.SelectedSubjects=logical(full(sparse(nsubjects,1,1,CONN_x.Setup.nsubjects,1)));
            [SPMall(1).xX.isSurface,SPMall(1).xX.isMtx]=conn_surf_dimscheck(SPMall(1).xY.VY(1).dim); %,isequal(SPMall(1).xY.VY(1).dim,conn_surf_dims(8).*[1 1 2]);
            SPMall(1).xX.X=kron(eye(nrepeated),X(nsubjects,nsubjecteffects));%CONN_x.Results.xX.X;
            %SPMall(1).xX.name=repmat({CONN_x.Setup.l2covariates.names{nsubjecteffects}},[1,nrepeated]);%CONN_x.Results.xX.name;
            SPMall(1).xX.iH     = [];
            SPMall(1).xX.iC     = 1:size(SPMall(1).xX.X,2);
            SPMall(1).xX.iB     = [];
            SPMall(1).xX.iG     = [];
            SPMall(1).xGX       = [];
            
            if nrepeated>1
                xVi=struct('I',[repmat((1:size(SPMall(1).xY.VY,1))',[nrepeated,1]),reshape(repmat(1:nrepeated,[size(SPMall(1).xY.VY,1),1]),[],1)],'var',[0,1],'dep',[1,0]);
                SPMall(1).xVi=spm_non_sphericity(xVi);
            end
        end
    else
        h=conn_waitbar(0,['Step ',num2str(sum(options<=16)),'/',num2str(length(options)),': Functional data second-level analyses']);
        doskip=false(size(nsources));
        for n1=1:length(nsources),
            [foldername,foldername_back]=conn_resultsfolder('sources',state,sources(nsources(n1)),1,1);
            for nfolderbak=1:numel(foldername_back),
                if isdir(fullfile(filepathresults2,foldername_back{nfolderbak})), foldername=foldername_back{nfolderbak}; break; end % backwards-compatibility with existing results
            end
%             txttmp=sources{n1};
%             if state==3,txttmp=conn_v2v('pcleartext',txttmp); end
%             if REDUCENAMES&&length(txttmp)>2&&strcmp(txttmp(1:3),'BA.'),
%                 idxtxttmp=find(txttmp=='.');idxtxttmp2=find(txttmp=='_');
%                 if length(idxtxttmp)>1&&length(idxtxttmp2)>1&&idxtxttmp(end)-1<idxtxttmp2(end-1), 
%                     txttmp=txttmp([1:idxtxttmp(end)-1,idxtxttmp2(end-1):end]); 
%                 end;
%             end            
%             foldername=[txttmp,'.'];
%             foldername(foldername==' '|foldername==filesep|foldername=='*')='_';
%             if numel(foldername)>100, foldername=foldername(1:100); end
            [ok,nill]=mkdir(filepathresults2,foldername);
            if ok,filepathresults3{n1}=fullfile(filepathresults2,foldername);
            else,filepathresults3{n1}=uigetdir(filepathresults2,'Select a directory to write the results');end
            if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'overwrite')&&strcmpi(CONN_x.gui.overwrite,'no')&&conn_existfile(fullfile(filepathresults3{n1},'SPM.mat')), doskip(n1)=true; end
        end
        if any(doskip)
            nsources=nsources(~doskip);
            csources=csources(:,~doskip);
            filepathresults3=filepathresults3(~doskip);
            if isempty(nsources), SPMall=[]; end
        end
        clear SE Y;
        if FORCEORTH&&numel(nconditions)>1&&~isequal(cconditions,eye(numel(nconditions))),cconditions=spm_orth(cconditions','norm')';end
        Zfiles={};
        Znames={};
        Zcontr=cconditions;
        if ~isempty(nsources)
            for nsub=1:length(nsubjects),%CONN_x.Setup.nsubjects,
                for n1=1:length(nsources),
                    altestsmooth=false;
                    nroi=nsources(n1);
                    roiname=sources{nsources(n1)};
                    if state==2
                        [iroi,isnew]=conn_sourcenames(roiname,'-');
                        if isnew, error(['Non-existing ROI first-level data for ',roiname,' subject ',num2str(nsubjects(nsub)),'. Please repeat first-level analyses']); return; end
                        troiname=regexprep(roiname,{'^BA\.(\d+) \(([LR])\)\. .*','^\((-?\d+),(-?\d+),(-?\d+)\)$','^SLrois\.|^aal\.|^atlas\.|^networks\.','\s\(([LR])\)','([^\(\)]+)\(.+\)\s*$'},{'$1$2','($1 $2 $3)','',' ${lower($1)}','$1'});
                    else
                        altestsmooth=altestsmooth|ismember(conn_v2v('fieldtext',roiname,1),{'2'});
                        %altestsmooth=altestsmooth|~isempty(strmatch('connectome',roiname));
                        [iroi,isnew,ncomp]=conn_v2v('match_extended',roiname);
                        %                     [iroi,isnew]=conn_v2v('match_measures',CONN_x.vvAnalyses.regressors,nroi,'-');
                        if isnew, error(['Non-existing voxel-to-voxel first-level data for ',roiname,' subject ',num2str(nsubjects(nsub)),'. Please repeat first-level analyses']); return; end
                        troiname=conn_v2v('cleartext',roiname);
                    end
                    if modeltype==2,%||(length(nconditions)>1&&size(cconditions,1)==1),%||any(cconditions~=1),
                        b=0;
                        tcontrastname='';
                        for n0=1:length(nconditions),
                            ncondition=nconditions(n0);
                            if state==2
                                filename=fullfile(filepathresults1,['BETA_Subject',num2str(nsubjects(nsub),CONN_x.opt.fmt1),'_Condition',num2str(icondition(ncondition),'%03d'),'_Source',num2str(iroi,'%03d'),'.nii']);
                            else
                                filename=fullfile(filepathresults1,['BETA_Subject',num2str(nsubjects(nsub),'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'_Measure',num2str(iroi,'%03d'),'_Component',num2str(ncomp,'%03d'),'.nii']);
                            end
                            Zfiles{nsub,n0,n1}=filename;
                            Znames{n0}=CONN_x.Setup.conditions.names{ncondition};
                            a=spm_vol(filename);
                            b=b+cconditions(n0)*spm_read_vols(a);
                            if cconditions(n0)==1, tcontrastname=[tcontrastname,'+',CONN_x.Setup.conditions.names{ncondition},'_',troiname];
                            elseif cconditions(n0)==-1, tcontrastname=[tcontrastname,'-',CONN_x.Setup.conditions.names{ncondition},'_',troiname];
                            else tcontrastname=[tcontrastname,num2str(cconditions(n0)),'*',CONN_x.Setup.conditions.names{ncondition},'_',troiname];
                            end
                            %                             tcontrastname=[tcontrastname,char('+'+2*(sign(cconditions(n0))<0)),num2str(cconditions(n0)),CONN_x.Setup.conditions.names{ncondition},'_',roiname];
                            %                         if n0==1&&n1==1&&nsub==1
                            %                             [gridx,gridy,gridz]=ndgrid(1:a.dim(1),1:a.dim(2),1:a.dim(3));xyz=a.mat*[gridx(:),gridy(:),gridz(:),ones(numel(gridx),1)]'; adim=a.dim(1:3);
                            %                         end
                            %                         b=b+cconditions(n0)*reshape(spm_get_data(a,pinv(a.mat)*xyz),adim);%spm_read_vols(a);
                            %                         if state==2&&modeltype==2&&nsub==1,
                            %                             filename=fullfile(filepathresults1,['resultsDATA_Condition',num2str(ncondition,'%03d'),'_Source',num2str(iroi,'%03d'),'.mat']);
                            %                             Y(n1,n0)=conn_vol(filename);
                            %                             if n1==1,
                            %                                 filename=fullfile(filepathresults1,['seDATA_Condition',num2str(ncondition,'%03d'),'.mat']);
                            %                                 SE(n0)=conn_vol(filename);
                            %                             end
                            %                         end
                        end
                        filename=fullfile(filepathresults3{n1},['CON_Subject',num2str(nsubjects(nsub),'%03d'),'.nii']);
                        %filename=fullfile(filepathresults2,['CON_Subject',num2str(nsubjects(nsub)),'_Source',num2str(iroi,'%03d'),'.nii']);
                        a.fname=filename;
                        spm_write_vol(a,b);
                        filename={filename};
                        contrastname={tcontrastname};
                        contrast=1;
                    else
                        filename={};
                        contrastname={};
                        for n0=1:numel(nconditions)
                            ncondition=nconditions(n0);
                            if state==2
                                filename{end+1}=fullfile(filepathresults1,['BETA_Subject',num2str(nsubjects(nsub),CONN_x.opt.fmt1),'_Condition',num2str(icondition(ncondition),'%03d'),'_Source',num2str(iroi,'%03d'),'.nii']);
                                contrastname{end+1}=[CONN_x.Setup.conditions.names{ncondition},'_',troiname];
                            else
                                filename{end+1}=fullfile(filepathresults1,['BETA_Subject',num2str(nsubjects(nsub),'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'_Measure',num2str(iroi,'%03d'),'_Component',num2str(ncomp,'%03d'),'.nii']);
                                contrastname{end+1}=[CONN_x.Setup.conditions.names{ncondition},'_',troiname];
                            end
                            Zfiles{nsub,n0,n1}=filename{end};
                            Znames{n0}=CONN_x.Setup.conditions.names{ncondition};
                        end
                        if numel(nconditions)>1&&~isequal(cconditions,eye(numel(nconditions))),%&&size(cconditions,1)==1
                            for n2=1:size(cconditions,1)
                                for n1alt=1
                                    tcontrastname='';
                                    b=0;
                                    for n0=1:numel(nconditions)
                                        a=spm_vol(filename{n1alt,n0});
                                        b=b+spm_read_vols(a)*cconditions(n2,n0);
                                        if cconditions(n2,n0)==1, tcontrastname=[tcontrastname,'+',contrastname{n1alt,n0}];
                                        elseif cconditions(n2,n0)==-1, tcontrastname=[tcontrastname,'-',contrastname{n1alt,n0}];
                                        elseif cconditions(n2,n0)~=0, tcontrastname=[tcontrastname,num2str(cconditions(n2,n0),'%+g'),'*',contrastname{n1alt,n0}];
                                        end
                                        %                             tcontrastname=[tcontrastname,char('+'+2*(sign(cconditions(n0))<0)),num2str(cconditions(n0)),contrastname{n1alt,n0}];
                                    end
                                    tfilename=fullfile(filepathresults3{1},['CON_Subject',num2str(nsubjects(nsub),'%03d'),'_compA',num2str(n2,'%03d'),'_',num2str(n1alt,'%03d'),'.nii']);
                                    a.fname=tfilename;
                                    spm_write_vol(a,b);
                                    filename{n1alt,numel(nconditions)+n2}=tfilename;
                                    contrastname{n1alt,numel(nconditions)+n2}=tcontrastname;
                                end
                            end
                            filename=filename(:,end-size(cconditions,1)+(1:size(cconditions,1)));
                            contrastname=contrastname(:,end-size(cconditions,1)+(1:size(cconditions,1)));
                            [nill,contrast]=conn_mtxbase(cconditions);
                        else
                            contrast=cconditions;
                        end
                    end
                    SPMall(n1).xY.VY(nsub,:)=spm_vol(char(filename))';
                    SPMall(n1).altestsmooth=altestsmooth;
                    SPMall(n1).xX.name={};
                    for n01=1:numel(contrastname),for n00=1:numel(nsubjecteffects),SPMall(n1).xX.name{n00,n01}=[CONN_x.Setup.l2covariates.names{nsubjecteffects(n00)},'_',contrastname{n01}]; end; end
                    SPMall(n1).xX.name=SPMall(n1).xX.name(:)';
                    if state==2&&modeltype==2&&nsub==1, SPMall(n1).connvols.Y=Y(n1,:);SPMall(n1).connvols.SE=SE;SPMall(1).connvols.cconditions=cconditions;SPMall(1).connvols.csources=1;end
                end
                conn_waitbar(2/4*(nsub/length(nsubjects)),h,sprintf('Subject %d',nsub));
            end
        end
        for n1=1:length(nsources),
            nrepeated=size(SPMall(n1).xY.VY,2);
            SPMall(n1).xX_multivariate.X=X(nsubjects,nsubjecteffects);
            SPMall(n1).xX_multivariate.C=csubjecteffects;
            SPMall(n1).xX_multivariate.M=contrast;
            SPMall(n1).xX_multivariate.Xnames=CONN_x.Setup.l2covariates.names(nsubjecteffects);
            SPMall(n1).xX_multivariate.Ynames=contrastname;
            SPMall(n1).xX_multivariate.Znames=Znames;
            SPMall(n1).xX_multivariate.Zcontr=Zcontr;
            SPMall(n1).xX_multivariate.Zfiles=Zfiles(:,:,n1);
            SPMall(n1).xX.SelectedSubjects=logical(full(sparse(nsubjects,1,1,CONN_x.Setup.nsubjects,1)));
            [SPMall(n1).xX.isSurface,SPMall(n1).xX.isMtx]=conn_surf_dimscheck(SPMall(n1).xY.VY(1).dim); %isequal(SPMall(n1).xY.VY(1).dim,conn_surf_dims(8).*[1 1 2]);
            SPMall(n1).xX.X=kron(eye(nrepeated),X(nsubjects,nsubjecteffects));%CONN_x.Results.xX.X;
            %SPMall(n1).xX.name=repmat({CONN_x.Setup.l2covariates.names{nsubjecteffects}},[1,nrepeated]);%CONN_x.Results.xX.name;
            SPMall(n1).xX.iH     = [];
            SPMall(n1).xX.iC     = 1:size(SPMall(n1).xX.X,2);
            SPMall(n1).xX.iB     = [];
            SPMall(n1).xX.iG     = [];
            SPMall(n1).xGX       = [];
            if nrepeated>1
                xVi=struct('I',[repmat((1:size(SPMall(n1).xY.VY,1))',[nrepeated,1]),reshape(repmat(1:nrepeated,[size(SPMall(n1).xY.VY,1),1]),[],1)],'var',[0,1],'dep',[1,0]);
                SPMall(n1).xVi=spm_non_sphericity(xVi);
            end
        end
    end
    if dosinglecontrast==2,
        cd(filepathresults3{1});
        if isfield(CONN_x,'gui')&&(isfield(CONN_x.gui,'display_results')&&CONN_x.gui.display_results || isfield(CONN_x.gui,'display_contrast')&&~isempty(CONN_x.gui.display_contrast) || isfield(CONN_x.gui,'display_style')&&~isempty(CONN_x.gui.display_style) || isfield(CONN_x.gui,'display_options')&&~isempty(CONN_x.gui.display_options)),
            if isfield(CONN_x.gui,'display_contrast')&&~isempty(CONN_x.gui.display_contrast), ncon=CONN_x.gui.display_contrast; else ncon=1; end
            if isfield(CONN_x.gui,'display_style')&&~isempty(CONN_x.gui.display_style), style=CONN_x.gui.display_style; else style=[]; end
            if isfield(CONN_x.gui,'display_options')&&~isempty(CONN_x.gui.display_options), display_options=CONN_x.gui.display_options; else display_options={}; end
            fh=conn_display('SPM.mat',ncon,style);
            if ~isempty(display_options)
                if ~iscell(display_options),
                    conn_display(fh,display_options);
                elseif ~iscell(display_options{1})
                    conn_display(fh,display_options{:});
                else
                    for nfields=1:numel(display_options)
                        conn_display(fh,display_options{nfields}{:});
                    end
                end
            end
            varargout{1}=fh; 
        end
        cd(cwd);
    else,
        switch(modeltype),
            case 1, % random-effects
                for n1=1:length(SPMall),
                    cd(filepathresults3{n1});
                    SPM=SPMall(n1);
                    issurface=isfield(SPM.xX,'isSurface')&&SPM.xX.isSurface;
                    save('SPM.mat','SPM','-v7.3');
                    spm_unlink('mask.img','mask.hdr','mask.nii');
                    files=cat(1,dir('spmT_*.cluster.mat'),dir('nonparametric_p*.mat'));
                    if ~isempty(files)
                        files={files.name};
                        spm_unlink(files{:});
                    end
                    doboth=~issurface&&CONN_x.Setup.secondlevelanalyses==1;
                    if issurface||ismember(CONN_x.Setup.secondlevelanalyses,[1 3]) % nonparametric stats
                        if 0,% faster but requires more memory
                            Y=spm_read_vols(SPM.xY.VY);
                            mask=~any(isnan(Y),4)&any(diff(Y,1,4)~=0,4);
                            Y=permute(reshape(Y,[SPM.xY.VY(1).dim(1:3) size(SPM.xY.VY)]),[4,5,1,2,3]);
                            SPM.xX_multivariate.mask=mask;
                            [results_h,results_F,nill,SPM.xX_multivariate.dof,SPM.xX_multivariate.statsname]=conn_glm(SPM.xX_multivariate.X,Y(:,:,mask),SPM.xX_multivariate.C,SPM.xX_multivariate.M,[],false);
                            SPM.xX_multivariate.h=zeros([size(results_h,1),size(results_h,2),SPM.xY.VY(1).dim(1:3)]); SPM.xX_multivariate.h(:,:,mask)=results_h;
                            SPM.xX_multivariate.F=zeros([size(results_F,1),size(results_F,2),SPM.xY.VY(1).dim(1:3)]); SPM.xX_multivariate.F(:,:,mask)=results_F;
                        else
                            mask=ones(SPM.xY.VY(1).dim(1:3));
                            %[gridx,gridy]=ndgrid(1:SPM.xY.VY(1).dim(2),1:SPM.xY.VY(1).dim(3));
                            [gridx,gridy]=ndgrid(1:SPM.xY.VY(1).dim(1),1:SPM.xY.VY(1).dim(2));
                            xyz0=[gridx(:),gridy(:)]';
                            DOINCLUDEDOF=false; % note: placeholder for future voxel-specific dof values implementation
                            donefirst=false;
                            %for n2=1:SPM.xY.VY(1).dim(1)
                            for n2=1:SPM.xY.VY(1).dim(3)
                                %xyz=[n2+zeros(1,size(xyz0,2)); xyz0; ones(1,size(xyz0,2))];
                                xyz=[xyz0; n2+zeros(1,size(xyz0,2)); ones(1,size(xyz0,2))];
                                y=spm_get_data(SPM.xY.VY(:)',xyz);
                                maskthis=~any(isnan(y),1)&any(diff(y,1,1)~=0,1);
                                %mask(n2,:,:)=reshape(maskthis,[1 SPM.xY.VY(1).dim(2:3)]);
                                mask(:,:,n2)=reshape(maskthis,[SPM.xY.VY(1).dim(1:2)]);
                                if any(maskthis)
                                    fmaskthis=find(maskthis);
                                    %y=reshape(y,size(SPM.xY.VY,1),size(SPM.xY.VY,2),SPM.xY.VY(1).dim(2),SPM.xY.VY(1).dim(3));
                                    y=reshape(y,size(SPM.xY.VY,1),size(SPM.xY.VY,2),SPM.xY.VY(1).dim(1),SPM.xY.VY(1).dim(2));
                                    results_p=[];
                                    if ~donefirst
                                        donefirst=true;
                                        [results_h,results_F,results_p,results_dof,SPM.xX_multivariate.statsname]=conn_glm(SPM.xX_multivariate.X,y(:,:,maskthis),SPM.xX_multivariate.C,SPM.xX_multivariate.M);
                                        SPM.xX_multivariate.h=zeros([size(results_h,1),size(results_h,2),SPM.xY.VY(1).dim(1:3)]);
                                        SPM.xX_multivariate.F=zeros([size(results_F,1),size(results_F,2),SPM.xY.VY(1).dim(1:3)]);
                                        if DOINCLUDEDOF, SPM.xX_multivariate.dof=zeros([size(results_dof,1),size(results_dof,2),SPM.xY.VY(1).dim(1:3)]);
                                        else SPM.xX_multivariate.dof=results_dof;
                                        end
                                    else
                                        if DOINCLUDEDOF, [results_h,results_F,results_p,results_dof]=conn_glm(SPM.xX_multivariate.X,y(:,:,maskthis),SPM.xX_multivariate.C,SPM.xX_multivariate.M,[],false);
                                        else [results_h,results_F]=conn_glm(SPM.xX_multivariate.X,y(:,:,maskthis),SPM.xX_multivariate.C,SPM.xX_multivariate.M); 
                                        end
                                    end
                                    %SPM.xX_multivariate.h(:,:,n2,maskthis)=results_h;
                                    %SPM.xX_multivariate.F(:,:,n2,maskthis)=results_F;
                                    SPM.xX_multivariate.h(:,:,(n2-1)*prod(SPM.xY.VY(1).dim(1:2))+fmaskthis)=results_h;
                                    SPM.xX_multivariate.F(:,:,(n2-1)*prod(SPM.xY.VY(1).dim(1:2))+fmaskthis)=results_F;                                    
                                    if DOINCLUDEDOF, 
                                        if size(results_dof,3)==1&&nnz(maskthis)>1, results_dof=repmat(results_dof,[1,1,nnz(maskthis)]); end
                                        %SPM.xX_multivariate.dof(:,:,n2,maskthis)=results_dof; 
                                        SPM.xX_multivariate.dof(:,:,(n2-1)*prod(SPM.xY.VY(1).dim(1:2))+fmaskthis)=results_dof;
                                    end
                                end
                                %conn_waitbar(1/2+1/2*((n1-1+(.5+.5*~doboth)*(n2/SPM.xY.VY(1).dim(1)))/(length(SPMall))),h,sprintf('Analysis %d',n1));
                                conn_waitbar(1/2+1/2*((n1-1+(.5+.5*~doboth)*(n2/SPM.xY.VY(1).dim(3)))/(length(SPMall))),h,sprintf('Analysis %d',n1));
                            end
                        end
                        if size(SPM.xX_multivariate.F,1)==1&&size(SPM.xX_multivariate.F,2)==1
                            V=struct('mat',SPM.xY.VY(1).mat,'dim',SPM.xY.VY(1).dim,'fname','spmF_mv.nii','pinfo',[1;0;0],'n',[1,1],'dt',[spm_type('float32') spm_platform('bigend')]);
                            V=spm_write_vol(V,shiftdim(SPM.xX_multivariate.F,2));
                            try, spm_jsonwrite('spmF_mv.json',struct('dof',SPM.xX_multivariate.dof(:)','statsname',SPM.xX_multivariate.statsname)); end
                        end
                    end
                    if issurface
                        V=struct('mat',SPM.xY.VY(1).mat,'dim',SPM.xY.VY(1).dim,'fname','mask.img','pinfo',[1;0;0],'n',[1,1],'dt',[spm_type('uint8') spm_platform('bigend')]);
                        spm_write_vol(V,double(mask));
                        save('SPM.mat','SPM','-v7.3');
                        conn_disp('fprintf','\nSecond-level results saved in folder %s\n',pwd);
                        if isfield(CONN_x,'gui')&&(isfield(CONN_x.gui,'display_results')&&CONN_x.gui.display_results || isfield(CONN_x.gui,'display_contrast')&&~isempty(CONN_x.gui.display_contrast) || isfield(CONN_x.gui,'display_style')&&~isempty(CONN_x.gui.display_style) || isfield(CONN_x.gui,'display_options')&&~isempty(CONN_x.gui.display_options)),
                            if isfield(CONN_x.gui,'display_contrast')&&~isempty(CONN_x.gui.display_contrast), ncon=CONN_x.gui.display_contrast; else ncon=1; end
                            if isfield(CONN_x.gui,'display_style')&&~isempty(CONN_x.gui.display_style), style=CONN_x.gui.display_style; else style=[]; end
                            if isfield(CONN_x.gui,'display_options')&&~isempty(CONN_x.gui.display_options), display_options=CONN_x.gui.display_options; else display_options={}; end
                            fh=conn_display('SPM.mat',ncon,style);
                            if ~isempty(display_options)
                                if ~iscell(display_options), 
                                    conn_display(fh,display_options);
                                elseif ~iscell(display_options{1})
                                    conn_display(fh,display_options{:});
                                else
                                    for nfields=1:numel(display_options)
                                        conn_display(fh,display_options{nfields}{:});
                                    end
                                end
                            end
                        end
                    elseif ismember(CONN_x.Setup.secondlevelanalyses,[1 2]) % parametric stats
                        save('SPM.mat','SPM','-v7.3');
                        spm('Defaults','fmri');
                        try, spm_select('init'); end
                        try, spm_get_defaults('mat.format','-v7.3'); end
                        if isfield(CONN_x.Setup,'stats_ufp')&&~isempty(CONN_x.Setup.stats_ufp), spm_get_defaults('stats.fmri.ufp',CONN_x.Setup.stats_ufp); end
                        spm_unlink('mask.img','mask.hdr','mask.nii');
                        SPM=spm_spm(SPM);
                        c=kron(contrast,csubjecteffects); %CONN_x.Results.xX.C;
                        cname='connectivity result';
                        if size(c,1)==1, Statname='T'; else Statname='F'; end
                        if ~isfield(SPM.xX,'xKXs'), error('SPM analyses did not finish correctly'); end
                        SPM.xCon = spm_FcUtil('Set',cname,Statname,'c',c',SPM.xX.xKXs);
                        if 0,%adds other contrasts (one T-test per regressor, and one "effects of interest" F-test)
                            C=eye(length(nsubjecteffects));
                            if size(C,1)>1,
                                for n2=1:size(C,1),
                                    c=C(n2,:);
                                    cname=CONN_x.Setup.l2covariates.names{nsubjecteffects(n2)}; %CONN_x.Results.xX.name{n2};
                                    SPM.xCon(end+1) = spm_FcUtil('Set',cname,'T','c',c',SPM.xX.xKXs);
                                end
                                cname='effects of interest';
                                SPM.xCon(end+1) = spm_FcUtil('Set',cname,'F','c',C,SPM.xX.xKXs);
                            end
                        end
                        if isfield(SPM,'altestsmooth')&&SPM.altestsmooth, % modified smoothness estimation
                            SPM=conn_est_smoothness(SPM);
                            save('SPM.mat','SPM','-v7.3');
                        end
                        SPM=spm_contrasts(SPM,1:length(SPM.xCon));
                        SPM.xY.VY=SPM.xY.VY(:);
                        SPM.xsDes='';
                        save('SPM.mat','SPM','-v7.3');
                        %conn_cumdisp; 
                        conn_disp('fprintf','Second-level results saved in folder %s\n',pwd);
                        if state==2&&(~isempty(dir('con_0001.img'))||~isempty(dir('con_0001.nii'))),
                            switch(CONN_x.Analyses(CONN_x.Analysis).measure),
                                case {1,2}, %correlation
                                    if ~isempty(dir('con_0001.img')), V=spm_vol('con_0001.img');
                                    else V=spm_vol('con_0001.nii');
                                    end
                                    t=spm_read_vols(V);
                                    V.fname='corr_0001.img';
                                    spm_write_vol(V,tanh(t));
                            end
                        end
                        if isfield(CONN_x,'gui')&&(isfield(CONN_x.gui,'display_results')&&CONN_x.gui.display_results || isfield(CONN_x.gui,'display_contrast')&&~isempty(CONN_x.gui.display_contrast) || isfield(CONN_x.gui,'display_style')&&~isempty(CONN_x.gui.display_style) || isfield(CONN_x.gui,'display_options')&&~isempty(CONN_x.gui.display_options)),
                            if isfield(CONN_x.gui,'display_contrast')&&~isempty(CONN_x.gui.display_contrast), ncon=CONN_x.gui.display_contrast; else ncon=1; end
                            if isfield(CONN_x.gui,'display_style')&&~isempty(CONN_x.gui.display_style), style=CONN_x.gui.display_style; else style=[]; end
                            if isfield(CONN_x.gui,'display_options')&&~isempty(CONN_x.gui.display_options), display_options=CONN_x.gui.display_options; else display_options={}; end
                            fh=conn_display('SPM.mat',ncon,style);
                            if ~isempty(display_options)
                                if ~iscell(display_options), 
                                    conn_display(fh,display_options);
                                elseif ~iscell(display_options{1})
                                    conn_display(fh,display_options{:});
                                else
                                    for nfields=1:numel(display_options)
                                        conn_display(fh,display_options{nfields}{:});
                                    end
                                end
                            end
                        end
                    elseif ismember(CONN_x.Setup.secondlevelanalyses,[1 3]) % nonparametric stats
                        V=struct('mat',SPM.xY.VY(1).mat,'dim',SPM.xY.VY(1).dim,'fname','mask.img','pinfo',[1;0;0],'n',[1,1],'dt',[spm_type('uint8') spm_platform('bigend')]);
                        spm_write_vol(V,double(mask));
                        save('SPM.mat','SPM','-v7.3');
                        if state==2&&(~isempty(dir('con_0001.img'))||~isempty(dir('con_0001.nii'))),
                            switch(CONN_x.Analyses(CONN_x.Analysis).measure),
                                case {1,2}, %correlation
                                    if ~isempty(dir('con_0001.img')), V=spm_vol('con_0001.img');
                                    else V=spm_vol('con_0001.nii');
                                    end
                                    t=spm_read_vols(V);
                                    V.fname='corr_0001.img';
                                    spm_write_vol(V,tanh(t));
                            end
                        end
                        if isfield(CONN_x,'gui')&&(isfield(CONN_x.gui,'display_results')&&CONN_x.gui.display_results || isfield(CONN_x.gui,'display_contrast')&&~isempty(CONN_x.gui.display_contrast) || isfield(CONN_x.gui,'display_style')&&~isempty(CONN_x.gui.display_style) || isfield(CONN_x.gui,'display_options')&&~isempty(CONN_x.gui.display_options)),
                            if isfield(CONN_x.gui,'display_contrast')&&~isempty(CONN_x.gui.display_contrast), ncon=CONN_x.gui.display_contrast; else ncon=1; end
                            if isfield(CONN_x.gui,'display_style')&&~isempty(CONN_x.gui.display_style), style=CONN_x.gui.display_style; else style=[]; end
                            if isfield(CONN_x.gui,'display_options')&&~isempty(CONN_x.gui.display_options), display_options=CONN_x.gui.display_options; else display_options={}; end
                            fh=conn_display('SPM.mat',ncon,style);
                            if ~isempty(display_options)
                                if ~iscell(display_options), 
                                    conn_display(fh,display_options);
                                elseif ~iscell(display_options{1})
                                    conn_display(fh,display_options{:});
                                else
                                    for nfields=1:numel(display_options)
                                        conn_display(fh,display_options{nfields}{:});
                                    end
                                end
                            end
                        end
                    end
                    conn_waitbar(1/2+1/2*(n1/(length(SPMall))),h,sprintf('Analysis %d',n1));
                end
                conn_waitbar('close',h);
                cd(cwd);
            case 2, % fixed-effects
                for n1=1:length(SPMall),
                    cd(filepathresults3{n1});
                    csources=SPMall(n1).connvols.csources;
                    cconditions=SPMall(n1).connvols.cconditions;
                    Y=SPMall(n1).connvols.Y;
                    SE=SPMall(n1).connvols.SE;
                    xf=SPMall(n1).xX.X;
                    filename='CON_spmT_0001.mat';
                    Yout=Y(1); Yout.fname=filename;
                    Yout.size.Nt=1;
                    Yout.DOF=sum([SE(:).DOF]);
                    Yout=conn_init_vol(Yout);
                    for slice=1:Y(1).matdim.dim(3),
                        yf=0;
                        se.data=0;
                        se.dof=0;
                        for ncondition=1:size(Y,2),
                            for nsource=1:size(Y,1),
                                [temp,y.idx]=conn_get_slice(Y(nsource,ncondition),slice);
                                yf=yf+temp*csources(nsource)*cconditions(ncondition);
                            end;
                            [temp,nill]=conn_get_slice(SE(ncondition),slice);
                            se.data=se.data+sum(csources.^2)*(cconditions(ncondition)*temp).^2;
                            se.dof=se.dof+SE(ncondition).DOF;
                        end
                        se.data=sqrt(se.data);
                        [B,opt]=conn_glmunivariate('estimatefixed',xf,yf,se);
                        [H,F,p,dof,R]=conn_glmunivariate('evaluatefixed',opt,[],csubjecteffects); 
                        conn_write_slice(Yout,F,slice);
                        conn_waitbar(1/2+1/2*((n1-1+slice/Y(1).matdim.dim(3))/(length(SPMall))),h);
                    end
                    t=nan+zeros(Yout.matdim.dim);
                    y=conn_get_time(Yout,1,[],false);
                    %V=spm_vol(deblank(CONN_x.Setup.structural{nsub}{1})); V=V(1);
                    V=CONN_x.Setup.structural{1}{1}{3}; V=V(1);
                    t(Yout.voxels)=y;
                    V.fname='spmT_0001.img';
                    if isfield(V,'dt'), V.dt=[spm_type('float32') spm_platform('bigend')];
                    elseif length(V.dim)>3, V.dim(4)=spm_type('float32'); end
                    spm_write_vol(V,t);
                    clear SPM;
                    SPM.xCon(1).Vspm=V;
                    SPM.xCon(1).eidf=inf;
                    SPM.xX.erdf=Yout.DOF;
                    SPM.xVol.R=[];
                    SPM.xVol.S=[];
                    save('SPM.mat','SPM','-v7.3');
                    
                    if isfield(CONN_x,'gui')&&(isfield(CONN_x.gui,'display_results')&&CONN_x.gui.display_results),conn_display('SPM.mat',1); end
                    conn_waitbar(1/2+1/2*(n1/(length(SPMall))),h);
                end
                conn_waitbar('close',h);
                cd(cwd);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates second-level ROI.mat results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (any(floor(options)==17) && any(CONN_x.Setup.steps([1])) && ~(isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'steps')&&~any(CONN_x.gui.steps([1])))) || (any(options==18) && any(CONN_x.Setup.steps([3])) && ~(isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'steps')&&~any(CONN_x.gui.steps([3])))),
    dovv=any(options==18);
    filepathresults2=CONN_x.folders.secondlevel;
    if any(options==17)
        if nargin>1&&~isempty(varargin{1}),CONN_x.Analysis=varargin{1}; end % selected analysis only
        if ischar(CONN_x.Analysis), CONN_x.Analysis=find(strcmp({CONN_x.Analyses.name},CONN_x.Analysis)); end
        if nargin>2&&~isempty(varargin{2}), % selected contrasts only
            ncontrast=varargin{2};
            if ischar(ncontrast), ncontrast=find(strcmp(CONN_x.Results.saved.names,ncontrast)); end
            assert(numel(ncontrast)==1,'incorrect number of contrasts');
            newneffects=CONN_x.Results.saved.nsubjecteffects{ncontrast};
            newceffects=CONN_x.Results.saved.csubjecteffects{ncontrast};
            newnconditions=CONN_x.Results.saved.nconditions{ncontrast};
            newcconditions=CONN_x.Results.saved.cconditions{ncontrast};
            [ok1,i1]=ismember(newneffects,CONN_x.Setup.l2covariates.names(1:end-1));
            [i1,idx1]=sort(i1);
            CONN_x.Results.xX.nsubjecteffects=i1;
            CONN_x.Results.xX.csubjecteffects=newceffects(:,idx1);
            CONN_x.Results.xX.nsubjecteffectsbyname=CONN_x.Setup.l2covariates.names(i1);
            [ok2,i2]=ismember(newnconditions,CONN_x.Setup.conditions.names(1:end-1));
            [i2,idx2]=sort(i2);
            CONN_x.Results.xX.nconditions=i2;
            CONN_x.Results.xX.cconditions=newcconditions(:,idx2);
            CONN_x.Results.xX.nconditionsbyname=CONN_x.Setup.conditions.names(i2);
        end
        filepathresults1=fullfile(CONN_x.folders.firstlevel,CONN_x.Analyses(CONN_x.Analysis).name);
        sources=CONN_x.Analyses(CONN_x.Analysis).sources;
        dosinglecontrast=0;
        nsources=1:length(sources);
        csources=ones(1,length(nsources));
        hcw=0;
        domvpa=[];
        %if nargin<2, domvpa=-1;
        %else, domvpa=varargin{1}; 
        %end
    elseif any(options==17.5), % single between-sources contrast
        filepathresults1=fullfile(CONN_x.folders.firstlevel,CONN_x.Analyses(CONN_x.Analysis).name);
        sources=CONN_x.Analyses(CONN_x.Analysis).sources;
        dosinglecontrast=1;
        nsources=varargin{1};
        csources=varargin{2};
        domvpa=[];
        hcw=[];
    else % vv post hoc
        filepathresults1=varargin{1};
        filename=fullfile(filepathresults1,['resultsROI_Condition',num2str(icondition(1),'%03d'),'.mat']);
        names={}; conn_loadmatfile(filename,'names');
        sources=names;
        dosinglecontrast=0;
        nsources=1:length(sources);
        csources=ones(1,length(nsources));
        hcw=0;
        if nargin<3, domvpa=[];
        else, domvpa=varargin{2}; 
        end
    end
    nconditions=CONN_x.Results.xX.nconditions;
    cconditions=CONN_x.Results.xX.cconditions;
    nsubjecteffects=CONN_x.Results.xX.nsubjecteffects;
    csubjecteffects=CONN_x.Results.xX.csubjecteffects;
    icondition=[];isnewcondition=[];for ncondition=nconditions(:)',[icondition(ncondition),isnewcondition(ncondition)]=conn_conditionnames(CONN_x.Setup.conditions.names{ncondition}); end
    if any(isnewcondition), error(['Some conditions have not been processed yet. Re-run previous step']); end
    
    if isfield(CONN_x.Results,'foldername')&&~isempty(CONN_x.Results.foldername),
        try, conn_fileutils('mkdir',filepathresults2,CONN_x.Results.foldername); end
        %[ok,nill]=mkdir(filepathresults2,CONN_x.Results.foldername);
        filepathresults2=fullfile(filepathresults2,CONN_x.Results.foldername);
        CONN_x.Results.foldername=[];
    else,
        [foldername,foldername_back]=conn_resultsfolder('subjectsconditions',1,nsubjecteffects,csubjecteffects,nconditions,cconditions);
        for nfolderbak=1:numel(foldername_back), 
            if conn_fileutils('isdir',fullfile(filepathresults2,foldername_back{nfolderbak})), foldername=foldername_back{nfolderbak}; break; end % backwards-compatibility with existing results
        end
        if ~nargout||nargout>1,
            try, 
                conn_fileutils('mkdir',filepathresults2,foldername);
                filepathresults2=fullfile(filepathresults2,foldername);
            catch
                filepathresults2=uigetdir(filepathresults2,'Select a directory to write the results');
            end
        %             if isfield(CONN_x.Results,'foldername')&&~isempty(CONN_x.Results.foldername),
        %                 [ok,nill]=mkdir(filepathresults2,CONN_x.Results.foldername);
        %                 filepathresults2=fullfile(filepathresults2,CONN_x.Results.foldername);
        %             else,
        %                 filepathresults2=CONN_x.folders.secondlevel;
        %                 filepathresults2=uigetdir(filepathresults2,'Select a directory to write the results');
        %             end
            if ~ischar(filepathresults2), return; end;
        end
    end
    recompute=1;
    if conn_existfile(fullfile(filepathresults2,'ROI.mat'))
       recompute=0;
       if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'overwrite'), recompute=strcmpi(CONN_x.gui.overwrite,'yes'); 
       else REDO=conn_questdlg('','results explorer','Load existing analysis results', 'Recompute/overwrite results', 'Load existing analysis results'); if strcmp(lower(REDO),'recompute/overwrite results'),recompute=1; end; 
       end
       if recompute
           f=conn_dir(fullfile(filepathresults2,'nonparametricroi_*.mat'),'-R','-cell');
           if ~isempty(f), conn_fileutils('spm_unlink',f{:}); end
       end
    end
    if ~recompute
        ROI={}; conn_loadmatfile(fullfile(filepathresults2,'ROI.mat'),'ROI');
        varargout{1}=ROI;
        if nargout>1, varargout{2}=filepathresults2; end
    else
        if ~isempty(hcw), hcw=conn_waitbar(0,'Computing ROI-level results. Please wait...'); end
        X=zeros(CONN_x.Setup.nsubjects,length(CONN_x.Setup.l2covariates.names)-1);
        for nsub=1:CONN_x.Setup.nsubjects,
            for ncovariate=1:length(CONN_x.Setup.l2covariates.names)-1;
                X(nsub,ncovariate)=CONN_x.Setup.l2covariates.values{nsub}{ncovariate};
            end
        end
        se.data=0;
        se.dof=0;
        y=0;
        Y=repmat({0},[1,length(nsources)]);
        orig_conditions={}; orig_sources={};
        for n0=1:length(nconditions),
            ncondition=nconditions(n0);
            filename=fullfile(filepathresults1,['resultsROI_Condition',num2str(icondition(ncondition),'%03d'),'.mat']);
            if ~conn_existfile(filename),
                Ransw=conn_questdlg({'First-level ROI analyses have not completed.','Perform now?'},'warning','Yes','No','Yes');
                if strcmp(Ransw,'Yes'), conn_process('analyses_roi'); end
            end
            Z=[]; names={}; names2={}; xyz={}; conn_loadmatfile(filename,'Z','names','names2','xyz','-cache');%,'SE','DOF');
            orig_conditions{n0}=CONN_x.Setup.conditions.names{ncondition};
            iroi=[];clear ixyz; for nroi=1:length(sources),
                roiname=sources{nroi};
                tiroi=strmatch(roiname,names,'exact');
                if isempty(tiroi), error(['Non-existing ROI first-level data for ',roiname,'. Please repeat first-level analyses']); return; end
                iroi(nroi)=tiroi(1);
                ixyz{nroi}=xyz{tiroi(1)};
            end
            for n1=1:length(nsources),
                nroi=nsources(n1);
                orig_sources{n1}=sources{nroi};
                yt=permute(Z(iroi(nroi),:,:),[3,2,1]);
                if dosinglecontrast,
                    if n0==1&&n1==1, y=zeros([size(yt,1),size(yt,2),length(nsources),length(nconditions)]); end
                    y(:,:,n1,n0)=yt;% subjects x rois x nsources x nconditions
                else
                    if n0==1, Y{n1}=zeros([size(yt,1),size(yt,2),length(nconditions)]); end
                    Y{n1}(:,:,n0)=yt;% subjects x rois x nconditions
                end
                %y=y+bsxfun(@times,yt,shiftdim(csources(:,n1)*cconditions(:,n0)',-2)); % subjects x rois x sourcecontrasts x conditioncontrasts
                %Y{n1}=Y{n1}+bsxfun(@times,yt,shiftdim(cconditions(:,n0),-2)); % subjects x rois x conditioncontrasts
            end
            %         se.data=se.data+(cconditions(n0)*SE).^2;
            %         se.dof=se.dof+DOF;
        end
        if isequal(cconditions,'var'),
            y=std(y,1,4);
            for n1=1:length(nsources), Y{n1}=std(Y{n1},1,3); end
            cconditions=1;
        end
        %     se.data=sqrt(se.data);
        %if nargin>1, dosinglecontrast=(varargin{1}==1); else, dosinglecontrast=~(length(nsources)>1 & all(csources==1/length(nsources))); end
        clear SPMall filepathresults3;
        if dosinglecontrast,
            nsubjects=find(any(X(:,nsubjecteffects)~=0,2)&~any(isnan(X(:,nsubjecteffects)),2)&~any(any(all(isnan(y),2),3),4));
            ROIall.xX.X=X(nsubjects,nsubjecteffects);%CONN_x.Results.xX.X;
            ROIall.xX.name={CONN_x.Setup.l2covariates.names{nsubjecteffects}};%CONN_x.Results.xX.name;
            ROIall.xX.SelectedSubjects=logical(full(sparse(nsubjects,1,1,CONN_x.Setup.nsubjects,1)));
            ROIall.orig_sources=orig_sources;
            ROIall.orig_conditions=orig_conditions;
            %         if ~nargout, filepathresults3{n1}=filepathresults2; end
            ROIall.y=y(nsubjects,:,:,:);
            ROIall.c2=kron(cconditions,csources);
            %         ROIall.se=se;
            %         ROIall.se.data=sqrt(sum(csources.^2))*ROIall.se.data;
            %         ROIall.se.data=ROIall.se.data(nsubjects,:);
            %         ROIall.se.dof=ROIall.se.dof(nsubjects);
        else
            for n1=1:length(nsources),
                nsubjects=find(any(X(:,nsubjecteffects)~=0,2)&~any(isnan(X(:,nsubjecteffects)),2)&~any(any(all(isnan(Y{n1}),2),3),4));
                ROIall(n1).xX.X=X(nsubjects,nsubjecteffects);%CONN_x.Results.xX.X;
                ROIall(n1).xX.name={CONN_x.Setup.l2covariates.names{nsubjecteffects}};%CONN_x.Results.xX.name;
                ROIall(n1).xX.SelectedSubjects=logical(full(sparse(nsubjects,1,1,CONN_x.Setup.nsubjects,1)));
                nroi=nsources(n1);
                roiname=sources{nroi};
                %             if ~nargout,[ok,nill]=mkdir(filepathresults2,roiname); filepathresults3{n1}=fullfile(filepathresults2,roiname); end
                ROIall(n1).y=Y{n1}(nsubjects,:,:);
                ROIall(n1).c2=cconditions;
                %             ROIall(n1).se=se;
                %             ROIall(n1).se.data=ROIall(n1).se.data(nsubjects,:);
                %             ROIall(n1).se.dof=ROIall(n1).se.dof(nsubjects);
            end
        end
        clear ROIout;
        for n1=1:length(ROIall),
            ROI=ROIall(n1);
            ROI.names=sources;
            ROI.xyz=ixyz;
            ROI.names2=names2;
            ROI.xyz2=xyz;
            ROI.cname='connectivity result';
            ROI.c=csubjecteffects;%CONN_x.Results.xX.C; 
            ROI.ynames=orig_conditions;
            
            if ~isfield(CONN_x.Results.xX,'modeltype')||isempty(CONN_x.Results.xX.modeltype), CONN_x.Results.xX.modeltype=1; end
            if CONN_x.Results.xX.modeltype==1,
                [ROI.h,ROI.F,ROI.p,ROI.dof,ROI.statsname,ROI.B]=conn_glm(ROI.xX.X,permute(ROI.y(:,:,:),[1,3,2]),ROI.c,ROI.c2);
                ROI.B=reshape(permute(ROI.B,[2,1,3]),[size(ROI.y,3),size(ROI.y,4),size(ROI.B,1),size(ROI.B,3)]); % sources x conditions x subject-effects x targets
                ROI.h=reshape(permute(ROI.h,[2,1,3]),[size(ROI.h,1)*size(ROI.h,2),size(ROI.h,3)]);
                if size(ROI.h,1)>1, ROI.h=sqrt(sum(abs(ROI.h).^2,1)); end
                ROI.F=reshape(ROI.F,[size(ROI.F,1)*size(ROI.F,2),size(ROI.F,3)]);
                ROI.p=reshape(ROI.p,[size(ROI.p,1)*size(ROI.p,2),size(ROI.p,3)]);
                %[b,opt]=conn_glmunivariate('estimate',ROI.xX.X,ROI.y);
                %[ROI.h,ROI.F,ROI.p,ROI.dof,nill,ROI.statsname]=conn_glmunivariate('evaluate',opt,[],ROI.c);
            elseif CONN_x.Results.xX.modeltype==2,
                [b,opt]=conn_glmunivariate('estimatefixed',ROI.xX.X,ROI.y,ROI.se);
                [ROI.h,ROI.F,ROI.p,ROI.dof,nill,ROI.statsname]=conn_glmunivariate('evaluatefixed',opt,[],ROI.c);
            end
            if 0,%~isempty(domvpa),
                if any(domvpa<0), domvpa=1:size(ROI.y,2); end
                ndims=ceil(sqrt(size(ROI.y,1))/2);
                ndims=max(1,min(min(numel(domvpa),size(ROI.y,2)), ndims ));
                if ndims<numel(domvpa)
                    y=ROI.y(:,domvpa,:);
                    y(:,any(any(isnan(y),1),3),:)=[];
                    sy=[size(y),1,1];
                    y=reshape(permute(y,[1,3,2]),sy(1)*sy(3),sy(2));
                    [Q,D,R]=svd(y,0);
                    y=y*R(:,1:ndims);
                    ROI.MVPAy=permute(reshape(y,[sy(1),sy(3),ndims]),[1,3,2]);
                    d=D(1:size(D,1)+1:size(D,1)*min(size(D))).^2;
                    ROI.MVPApcacov=d(1:ndims)/sum(d);
                else
                    y=ROI.y(:,domvpa,:);
                    y=y(:,~any(any(isnan(y),1),3),:);
                    ROI.MVPAy=y;
                    ROI.MVPApcacov=[];
                end
                [ROI.MVPAh,ROI.MVPAF,ROI.MVPAp,ROI.MVPAdof,ROI.MVPAstatsname]=conn_glm(ROI.xX.X,ROI.MVPAy(:,:),ROI.c,kron(ROI.c2,eye(size(ROI.MVPAy,2))));
                if isequal(ROI.MVPAstatsname,'T'),
                    %ROI.MVPAstatsname='F'; ROI.MVPAdof=[1,ROI.MVPAdof]; ROI.MVPAF=ROI.MVPAF.^2;
                    ROI.MVPAp=2*min(ROI.MVPAp,1-ROI.MVPAp);
                end
            end
            ROIout(n1)=ROI;
            %         if ~(nargin>1), save(fullfile(filepathresults3{n1},'ROI.mat'),'ROI');
            %         elseif dosinglecontrast, varargout{1}=ROI;
            %         elseif n1<=length(nsources), varargout{1}(n1)=ROI; end
            if ~isempty(hcw), conn_waitbar(n1/(length(nsources)+1),hcw,sprintf('Source %d',n1)); end
        end
        varargout{1}=ROIout;
        if nargout>1, varargout{2}=filepathresults2; end
        if ~nargout||nargout>1,
            try
                ROI=ROIout;
                tinfo=whos('ROI');
                if tinfo.bytes>2e9, conn_savematfile(fullfile(filepathresults2,'ROI.mat'),'ROI','-v7.3');
                else conn_savematfile(fullfile(filepathresults2,'ROI.mat'),'ROI','-v7');
                end
                conn_disp(['ROI results saved in ',fullfile(filepathresults2,'ROI.mat')]);
            end
            try % adds summary info
                F=cat(1,ROI.F);
                p=cat(1,ROI.p);
                dof=cat(1,ROI.dof);
                statsname=ROI(1).statsname;
                if isequal(statsname,'T'), p=2*min(p,1-p); end
                names=ROI(1).names;
                F=F(:,1:size(F,1));
                p=p(:,1:size(F,1));
                P=nan(size(p));P(~isnan(p))=conn_fdr(p(~isnan(p)));
                ROI=ROI(1);
                
                summary.rois.names=ROI.names;
                summary.rois.xyz=ROI.xyz;
                summary.results.RRC_F=F;
                summary.results.RRC_p=p;
                summary.results.RRC_P=P;
                summary.results.RRC_names=names;               
                summary.design.contrast_within=ROI.c2;
                summary.design.contrast_between=ROI.c;
                summary.design.designmultivariateonly=1;
                summary.design.designmatrix=ROI.xX.X;
                summary.design.designmatrix_name=ROI.xX.name;
                try, summary.design.conditions=ROI.ynames;
                catch, summary.design.conditions={};
                end
                summary.design.subjects=find(ROI.xX.SelectedSubjects);
                summary.design.pwd=filepathresults2;
                if isfield(ROI,'ynames'),
                    summary.design.data=arrayfun(@(a,b)sprintf('subject%03d: %s',a,ROI.ynames{b}),repmat(reshape(find(ROI.xX.SelectedSubjects),[],1),[1,numel(ROI.ynames)]),repmat(1:numel(ROI.ynames),[nnz(ROI.xX.SelectedSubjects),1]),'uni',0);
                    summary.design.dataTitle=ROI.ynames;
                else
                    summary.design.data=arrayfun(@(a,b)sprintf('subject%03d: measure #%d',a,b),repmat(reshape(find(ROI.xX.SelectedSubjects),[],1),[1,size(ROI.c2,2)]),repmat(1:size(ROI.c2,2),[nnz(ROI.xX.SelectedSubjects),1]),'uni',0);
                    summary.design.dataTitle=arrayfun(@(a)sprintf('measure #%d',a),1:size(ROI.c2,2),'uni',0);
                end
                conn_savematfile(fullfile(filepathresults2,'ROI.mat'),'summary','-append');
            end
            try
                %if isfield(CONN_x,'gui')&&(isnumeric(CONN_x.gui)&&CONN_x.gui || isfield(CONN_x.gui,'display')&&CONN_x.gui.display || isfield(CONN_x.gui,'display_contrast')&&~isempty(CONN_x.gui.display_contrast) || isfield(CONN_x.gui,'display_style')&&~isempty(CONN_x.gui.display_style) || isfield(CONN_x.gui,'display_options')&&~isempty(CONN_x.gui.display_options)),
                if isfield(CONN_x,'gui')&&(isfield(CONN_x.gui,'display_results')&&CONN_x.gui.display_results || isfield(CONN_x.gui,'display_contrast')&&~isempty(CONN_x.gui.display_contrast) || isfield(CONN_x.gui,'display_style')&&~isempty(CONN_x.gui.display_style) || isfield(CONN_x.gui,'display_options')&&~isempty(CONN_x.gui.display_options)),
                    cwd=pwd;
                    try, cd(filepathresults2); end
                    if isfield(CONN_x.gui,'display_contrast')&&~isempty(CONN_x.gui.display_contrast), ncon=CONN_x.gui.display_contrast; else ncon=1; end
                    if isfield(CONN_x.gui,'display_style')&&~isempty(CONN_x.gui.display_style), style=CONN_x.gui.display_style; else style=[]; end
                    if isfield(CONN_x.gui,'display_options')&&~isempty(CONN_x.gui.display_options), display_options=CONN_x.gui.display_options; else display_options={}; end
                    fh=conn_display(fullfile(filepathresults2,'ROI.mat'),ncon,style);
                    if ~isempty(display_options)
                        if ~iscell(display_options),
                            conn_display(fh,display_options);
                        elseif ~iscell(display_options{1})
                            conn_display(fh,display_options{:});
                        else
                            for nfields=1:numel(display_options)
                                conn_display(fh,display_options{nfields}{:});
                            end
                        end
                    end
                    cd(cwd);
                end
            end
        end
        if ~isempty(hcw), conn_waitbar('close',hcw); end
    end
    %cd(cwd);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Displays second-level SPM.mat results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(options==19) && any(CONN_x.Setup.steps([2])) && ~(isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'steps')&&~any(CONN_x.gui.steps([2]))),
    conn_display;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% permutation/randomization anaylses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(options==20) % randomise voxel-level
    % 'X','Y','c','m','THR','THR_TYPE','SIDE','niters','simfilename','mask','groupingsamples'
    filename=varargin{1}; 
    if numel(varargin)>1, BlockDistributed=varargin{2}; 
    elseif isfield(CONN_x.pobj,'partition')&&numel(CONN_x.pobj.partition)>=1, BlockDistributed=CONN_x.pobj.partition(1); 
    else BlockDistributed=[];
    end
    load(filename);
    a=spm_vol(Y);
    y=[];
    if isempty(mask)
        y=spm_read_vols(a);
        mask=~any(isnan(y),4)&any(diff(y,1,4)~=0,4);
        y=reshape(y,size(y,1)*size(y,2)*size(y,3),size(y,4));
        y=y(mask,:);
    end
    fmask=find(mask);
    [i,j,k]=ind2sub(size(mask),fmask);
    xyz=[i(:) j(:) k(:)]';
    if isempty(y), y=spm_get_data(a,xyz)'; end
    y=permute(reshape(y,size(y,1),size(X,1),[]),[2,3,1]);
    if conn_surf_dimscheck(a), 
        surfparams=load(fullfile(fileparts(which(mfilename)),'utils','surf','surf_top.mat'),'A');
        adj=kron(speye(2),surfparams.A);
        adj=adj(fmask,fmask);
    else adj=xyz;
    end
    conn_randomise(X,y,c,m,THR,THR_TYPE,SIDE,[niters BlockDistributed],simfilename,[],adj,[],groupingsamples,false);
end
if any(options==21) % randomise ROI-to-ROI
    % 'X','Y','c','m','THR','THR_TYPE','SIDE','niters','simfilename'
    filename=varargin{1}; 
    if numel(varargin)>1, BlockDistributed=varargin{2}; 
    elseif isfield(CONN_x.pobj,'partition')&&numel(CONN_x.pobj.partition)>=1, BlockDistributed=CONN_x.pobj.partition(1); 
    else BlockDistributed=[];
    end
    load(filename);
    conn_randomise(X,Y,c,m,THR,THR_TYPE,SIDE,[niters BlockDistributed],simfilename,[],'matrix');
end
if any(options==22), % merges permutation/randomization results
    filename=varargin{1};
    newfilenames=conn_dir(conn_prepend('parallel_*_',filename),'-R');
    if ~isempty(newfilenames)
        newfilenames=cellstr(newfilenames);
        for n=1:numel(newfilenames)
            newfilename=newfilenames{n};
            if ~conn_existfile(filename)
                conn_disp('fprintf','merging permutation results file %s\n',newfilename);
                if ispc, [ok,msg]=system(sprintf('move "%s" "%s"',newfilename,filename));
                else [ok,msg]=system(sprintf('mv ''%s'' ''%s''',newfilename,filename));
                end
            elseif conn_existfile(newfilename)
                conn_disp('fprintf','merging permutation results file %s with %s\n',newfilename,filename);
                load(filename);
                new=load(newfilename);
                Hist_Voxel_stat=[Hist_Voxel_stat,new.Hist_Voxel_stat];
                Dist_Voxel_statmax=[Dist_Voxel_statmax,new.Dist_Voxel_statmax];
                Hist_Cluster_size=[Hist_Cluster_size,new.Hist_Cluster_size];
                Hist_Cluster_mass=[Hist_Cluster_mass,new.Hist_Cluster_mass];
                Hist_Cluster_score=[Hist_Cluster_score,new.Hist_Cluster_score];
                Dist_Cluster_sizemax=[Dist_Cluster_sizemax,new.Dist_Cluster_sizemax];
                Dist_Cluster_massmax=[Dist_Cluster_massmax,new.Dist_Cluster_massmax];
                Dist_Cluster_scoremax=[Dist_Cluster_scoremax,new.Dist_Cluster_scoremax];
                Hist_Seed_size=[Hist_Seed_size,new.Hist_Seed_size];
                Hist_Seed_mass=[Hist_Seed_mass,new.Hist_Seed_mass];
                Hist_Seed_score=[Hist_Seed_score,new.Hist_Seed_score];
                Dist_Seed_sizemax=[Dist_Seed_sizemax,new.Dist_Seed_sizemax];
                Dist_Seed_massmax=[Dist_Seed_massmax,new.Dist_Seed_massmax];
                Dist_Seed_scoremax=[Dist_Seed_scoremax,new.Dist_Seed_scoremax];
                Hist_Network_size=[Hist_Network_size,new.Hist_Network_size];
                Hist_Network_mass=[Hist_Network_mass,new.Hist_Network_mass];
                Hist_Network_score=[Hist_Network_score,new.Hist_Network_score];
                Dist_Network_sizemax=[Dist_Network_sizemax,new.Dist_Network_sizemax];
                Dist_Network_massmax=[Dist_Network_massmax,new.Dist_Network_massmax];
                Dist_Network_scoremax=[Dist_Network_scoremax,new.Dist_Network_scoremax];
                Pthr=[Pthr,new.Pthr];
                Pthr_type=[Pthr_type,new.Pthr_type];
                Pthr_side=[Pthr_side,new.Pthr_side];
                save(filename,'Pthr','Pthr_type','Pthr_side','Hist_Voxel_stat','Dist_Voxel_statmax','Hist_Cluster_size','Hist_Cluster_mass','Hist_Cluster_score','Dist_Cluster_sizemax','Dist_Cluster_massmax','Dist_Cluster_scoremax','Hist_Seed_size','Hist_Seed_mass','Hist_Seed_score','Dist_Seed_sizemax','Dist_Seed_massmax','Dist_Seed_scoremax','Hist_Network_size','Hist_Network_mass','Hist_Network_score','Dist_Network_sizemax','Dist_Network_massmax','Dist_Network_scoremax','-append');
                if ispc, [ok,msg]=system(sprintf('del "%s"',newfilename));
                else [ok,msg]=system(sprintf('rm ''%s''',newfilename));
                end
            else
                conn_disp('fprintf','WARNING: permutation results file %s not found. Statistics may be computed from fewer permutation/randomization samples than expected\n',newfilename);
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% post hoc ROI-to-ROI analyses from voxel-to-voxel analyses
% Creates resultsDATA_Condition###_Source###.mat files (combined first-level analyses)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(options==31) && any(CONN_x.Setup.steps([3])) && ~(isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'steps')&&~any(CONN_x.gui.steps([3]))),
    warning('off','MATLAB:DELETE:FileNotFound');
    [path,name,ext]=fileparts(CONN_x.filename);
    filepath=CONN_x.folders.preprocessing;
    %filepathresults=CONN_x.folders.firstlevel;
    nconditions=length(CONN_x.Setup.conditions.names)-1;
    if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'conditions')&&~isempty(CONN_x.gui.conditions), validconditions=CONN_x.gui.conditions; else validconditions=1:length(CONN_x.Setup.conditions.names)-1; end
    icondition=[];isnewcondition=[];for ncondition=1:nconditions,[icondition(ncondition),isnewcondition(ncondition)]=conn_conditionnames(CONN_x.Setup.conditions.names{ncondition}); end
    if any(isnewcondition(validconditions)), error(['Some conditions have not been processed yet. Re-run previous step']); end
    peaks=varargin{1};
    peaks2=varargin{2};
    folderout=varargin{3};
    optionPeaksSpheres=2;
    optionSourceROIs=1;
    optionTargetROIs=2;
    if isempty(folderout)||~ischar(folderout) % add new analysis to the current analysis list
        if ~isfield(CONN_x,'Analyses')||isempty(CONN_x.Analyses), n=1;tianalysis=1;
        else
            n=max(str2double(regexprep({CONN_x.Analyses.name},'PEAKS_(\d*)','$1')))+1;
            if isnan(n), n=1; end
            tianalysis=numel(CONN_x.Analyses)+1;
        end
        if ~ischar(folderout)
            hd=dialog('name','Post hoc ROI-to-ROI analysis options:','windowstyle','normal','resize','on','units','norm','position',[.3,.5,.3,.3],'color','w');
            set(hd,'userdata',struct(...
                'handles',[...
                uicontrol(hd,'style','text','units','norm','position',[.1,.8,.4,.1],'string','Analysis name','backgroundcolor','w'),...
                uicontrol(hd,'style','edit','units','norm','position',[.5,.8,.4,.1],'string',['PEAKS_',num2str(n,'%02d')]),...
                uicontrol(hd,'style','popupmenu','units','norm','position',[.1,.5,.8,.1],'string',{'Source-ROIs: suprathreshold peak voxels','Sources-ROIs: all peak voxels'},'value',optionSourceROIs,'callback','a=get(gcbf,''userdata'');set(a.handles(4),''value'',max(get(a.handles(3),''value''),get(a.handles(4),''value'')));'),...
                uicontrol(hd,'style','popupmenu','units','norm','position',[.1,.4,.8,.1],'string',{'Target-ROIs: suprathreshold peak voxels (ROI-to-ROI analyses)','Target-ROIs: all peak voxels (ROI-to-ROI analyses)','Target-ROIs: all voxels (seed-to-voxel analyses)'},'value',optionTargetROIs,'callback','a=get(gcbf,''userdata'');set(a.handles(4),''value'',max(get(a.handles(3),''value''),get(a.handles(4),''value'')));'),...
                uicontrol(hd,'style','popupmenu','units','norm','position',[.1,.6,.8,.1],'string',{'ROIs: individual peak voxels','ROIs: 10mm spheres around peak voxels'},'value',optionPeaksSpheres),...
                uicontrol(hd,'style','pushbutton','units','norm','position',[.33,.1,.33,.1],'string','OK','callback','uiresume(gcbf)'),...
                uicontrol(hd,'style','pushbutton','units','norm','position',[.66,.1,.33,.1],'string','Cancel','callback','set(gcbf,''userdata'',[]); uiresume(gcbf)')]));
            uiwait(hd);
            if ~ishandle(hd)||isempty(get(hd,'userdata')), delete(hd); return; end
            hdu=get(hd,'userdata');
            foldername=get(hdu.handles(2),'string');
            optionPeaksSpheres=get(hdu.handles(5),'value');
            optionSourceROIs=get(hdu.handles(3),'value');
            optionTargetROIs=get(hdu.handles(4),'value');
            delete(hd);
%             ok=0;
%             while ~ok,
%                 txt=inputdlg({'New analysis name:'},'conn',1,{['PEAKS_',num2str(n,'%02d')]});
%                 if isempty(txt), break; end
%                 txt{1}(txt{1}==' ')='_';
%                 [ok,nill]=mkdir(CONN_x.folders.firstlevel,txt{1});
%                 if ~ok, errordlg('Analysis name must be a valid name for a folder','conn'); end
%             end
%             foldername=txt{1};
            tianalysis=strmatch(foldername,{CONN_x.Analyses.name},'exact');
            if isempty(tianalysis), tianalysis=numel(CONN_x.Analyses)+1; else tianalysis=tianalysis(1); end
        else
            foldername=['PEAKS_',num2str(n,'%02d')];
        end
        [ok,nill]=mkdir(CONN_x.folders.firstlevel,foldername);
        folderout=fullfile(CONN_x.folders.firstlevel,foldername);
        addnew=true;
    else addnew=false;
    end
    if optionSourceROIs==2, peaks=peaks2;
    elseif optionTargetROIs==1, peaks2=peaks;
    end
    REDO=[];
    N=CONN_x.Setup.nsubjects*numel(validconditions);
    n=0;
    
    missingdata=arrayfun(@(n)isempty(dir(fullfile(folderout,['resultsROI_Condition',num2str(icondition(n),'%03d'),'.mat']))),validconditions);
    if isempty(REDO)&&~all(missingdata),
        if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'overwrite'), REDO=CONN_x.gui.overwrite; else, REDO=conn_questdlg('Overwrite existing results?','','Yes', 'No', 'Yes');end;
    else REDO='Yes'; 
    end
    if strcmp(lower(REDO),'yes'),
        if ~isequal(peaks2(1:size(peaks,1),:),peaks),error('peaks not sorted'); end
        xyz=num2cell(peaks,2)';
        names=cellfun(@(x)sprintf('(%d,%d,%d)',x(:)),xyz,'uni',0);
        xyz=num2cell(peaks2,2)';
        names2=cellfun(@(x)sprintf('(%d,%d,%d)',x(:)),xyz,'uni',0);
        if optionTargetROIs==3
            analysisbak=CONN_x.Analysis;
            CONN_x.Analysis=tianalysis;
            iroi=[];isnew=[];for nroi=1:numel(names),[iroi(nroi),isnew(nroi)]=conn_sourcenames(names{nroi},'+');end
        end
        h=conn_waitbar(0,'Extracting correlation matrix, please wait...');
        for ivalidcondition=1:numel(validconditions),
            ncondition=validconditions(ivalidcondition);
            Z=nan(size(peaks,1),size(peaks2,1),CONN_x.Setup.nsubjects);
            for nsub=1:CONN_x.Setup.nsubjects
                filename_B1=fullfile(filepath,['svvPC_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'.mat']); 
                filename_B10=fullfile(filepath,['vvPC_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'.mat']);
                if optionPeaksSpheres==1||isempty(dir(filename_B1)), 
                    if optionPeaksSpheres>1, conn_disp(['warning: non-existing file ',filename_B1,'. Using ROIs = individual peak voxels instead']); end
                    filename_B1=filename_B10;
                end
                Y1=conn_vol(filename_B1);
                if optionTargetROIs<3, % roi-to-roi
                    [x,idx]=conn_get_data(Y1,peaks2);
                    idxv=find(idx>0);
                    temp=zeros(size(peaks2,1));
                    temp(idxv,idxv)=x(:,idx(idxv))'*x(:,idx(idxv));
                    Z(:,:,nsub)=temp(1:size(peaks,1),:);
                else % seed-to-voxel
                    [x,idx]=conn_get_data(Y1,peaks);
                    idxv=find(idx>0);
                    Y10=conn_vol(filename_B10);
                    clear Yout;
                    for nroi=1:length(idxv);
                        filename=fullfile(folderout,['BETA_Subject',num2str(nsub,CONN_x.opt.fmt1),'_Condition',num2str(icondition(ncondition),'%03d'),'_Source',num2str(iroi(idxv(nroi)),'%03d'),'.nii']);
                        Yout{nroi}=struct('mat',Y10.matdim.mat,'dim',Y10.matdim.dim,'fname',filename,'pinfo',[1;0;0],'n',[1,1],'dt',[spm_type('float32') spm_platform('bigend')]);
                        conn_fileutils('deletefile',Yout{nroi}.fname);
                        Yout{nroi}=spm_create_vol(Yout{nroi});
                    end
                    for slice=1:Y10.matdim.dim(3),
                        [z,idxz]=conn_get_slice(Y10,slice);
                        temp=x(:,idx(idxv))'*z;
                        t=zeros(Y10.matdim.dim(1:2));
                        for nroi=1:length(idxv);
                            t(idxz)=temp(nroi,:);
                            Yout{nroi}=spm_write_plane(Yout{nroi},t,slice);
                        end
                    end
                end
                conn_waitbar(((ivalidcondition-1)*CONN_x.Setup.nsubjects+nsub)/numel(validconditions)/CONN_x.Setup.nsubjects,h,sprintf('Subject %d Condition %d',nsub,ncondition));
            end
            if optionTargetROIs<3
                filename=fullfile(folderout,['resultsROI_Condition',num2str(icondition(ncondition),'%03d'),'.mat']);
                save(filename,'Z','names','names2','xyz');
            end
        end
        conn_waitbar('close',h);
        if optionTargetROIs==3
            CONN_x.Analysis=analysisbak;
        end
        if addnew
            CONN_x.Analyses(tianalysis).name=foldername;
            CONN_x.Analysis=tianalysis;
            if optionTargetROIs<3
                CONN_x.Analyses(tianalysis).sourcenames={};
                CONN_x.Analyses(tianalysis).sources=names;
                CONN_x.Analyses(tianalysis).type=1;
                conn save;
                conn_msgbox({'ROI-to-ROI analyses created',['Go to ''second-level Results''-''>ROI-to-ROI'' and select analysis ',foldername]},'post-hoc analyses');
            else
                CONN_x.Analyses(tianalysis).sourcenames={};
                CONN_x.Analyses(tianalysis).sources=names;
                CONN_x.Analyses(tianalysis).type=2;
                conn save;
                conn_msgbox({'seed-to-voxel analyses created',['Go to ''second-level Results''-''>seed-to-voxel'' and select analysis ',foldername]},'post-hoc analyses');
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% creates QA plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(options==32)
    opts=varargin; % qafolder,procedures,SUBJECTS,validrois,validsets,nl2covariates,nl1contrasts
    if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'subjects'), validsubjects=CONN_x.gui.subjects; opts{3}=validsubjects; end
    if isfield(CONN_x,'pobj')&&isstruct(CONN_x.pobj)&&isfield(CONN_x.pobj,'subjects'), 
        validsubjects=CONN_x.pobj.subjects; opts{3}=validsubjects;
        if numel(opts)>=2&&any(ismember(opts{2},[13])) % plots that require all subjects in single procedure
            if ~(isfield(CONN_x.pobj,'partition')&&isequal(CONN_x.pobj.partition,[1 1]))
                conn_disp('WARNING: QA_DENOISE FC-QC plot parallelization not yet available (please run locally or as a single job instead). Skipping this plot');
                opts{2}(ismember(opts{2},[13]))=[];
            end
        end
    end    
    [varargout{1:nargout}]=conn_qaplots(opts{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% erodes masks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(options==33)
    opts=varargin; % validsubjects,validrois,REDO
    if isfield(CONN_x,'gui')&&isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'subjects'), validsubjects=CONN_x.gui.subjects; opts{1}=validsubjects; end
    if isfield(CONN_x,'pobj')&&isstruct(CONN_x.pobj)&&isfield(CONN_x.pobj,'subjects'),  validsubjects=CONN_x.pobj.subjects; opts{1}=validsubjects; end    
    varargout{1}=conn_maskserode(opts{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WORKSHOP DATASET
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(options==34)
    conn_batch_workshop_nyudataset(varargin{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONN_SETUP_PREPROC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(options==35)
    conn_setup_preproc(varargin{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% additional version-update step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if 0,%any(options==18),
% 	[path,name,ext]=fileparts(CONN_x.filename);
% % 	filepath=fullfile(path,name,'data');
%     filepath=CONN_x.folders.preprocessing;
%     filename=fullfile(filepath,['ROI_Subject',num2str(1,'%03d'),'_Condition',num2str(1,'%03d'),'.mat']);
%     X1=load(filename);
%     [X,nill,names]=conn_designmatrix(CONN_x.Analyses.regressors,X1,[]);
%     CONN_x.Analyses(CONN_x.Analysis).sources=names;
% %     CONN_x.Results.measure=CONN_x.Analyses.measure;
% end

return;
end

function [filename,cache]=conn_tempcache(in,ftype)
if isstruct(in)
    cache=in;
    if conn_existfile(cache.filename_cached), conn_process_movefile(cache.filename_cached,cache.filename_original); end
    switch(lower(ftype))
        case 'matc', 
            if conn_existfile(conn_prepend('',cache.filename_cached,'.matc')), 
                conn_process_movefile(conn_prepend('',cache.filename_cached,'.matc'),conn_prepend('',cache.filename_original,'.matc')); 
            end
    end
    filename=cache.filename_original;
else
    if 1 % same folder as target
        cache=struct('filename_original',in,'filename_cached',conn_prepend('cachetmp_',in));
    else % in conn_cache folder (e.g. scratch folder) %cache=struct('filename_original',in,'filename_cached',conn_cache('pull',in));
        [fpath,fname,fext]=fileparts(in);
        cache=struct('filename_original',in,'filename_cached',fullfile(conn_cache('getlocal'),['cachetmp_',fname,fext]));
    end
    filename=cache.filename_cached;
end
end

function conn_process_movefile(a,b,varargin)
if ispc, [ok,nill]=system(['move "',a,'" "',b,'"']);
else, [ok,nill]=system(['mv -f ''',a,''' ''',b,'''']);
end
if ~isequal(ok,0), error('Error moving file %s to %s, check target permissions',a,b); end
end


% from SPM5
function savefields(fnam,p)
if length(p)>1, error('Can''t save fields.'); end;
fn = fieldnames(p);
if numel(fn)==0, return; end;
for i=1:length(fn),
    eval([fn{i} '= p.' fn{i} ';']);
end;
if spm_matlab_version_chk('7') >= 0
    save(fnam,'-V6',fn{:});
else
    save(fnam,fn{:});
end
end
				
% function h=conn_waitbar(varargin)
% global CONN_x;
% h=[];
% if isfield(CONN_x,'gui')&&(isnumeric(CONN_x.gui)&&CONN_x.gui || isfield(CONN_x.gui,'display')&&CONN_x.gui.display),
%     if ischar(varargin{1}), h=varargin{2}; close(h);
%     else h=conn_timedwaitbar(varargin{:}); end
% else
%     if ischar(varargin{1}), conn_cumdisp;
%     elseif ~varargin{1}, 
%         conn_cumdisp; conn_disp(varargin{2}); 
%     else
%         str=[num2str(100*varargin{1},'%3.1f'),'% '];
%         if nargin>2,str=[str, ' (',varargin{3},')']; end
%         conn_cumdisp(str); 
%     end
% end
% end
% 
% function conn_cumdisp(txt)
% % CUMDISP persistent disp
% % cumdisp; initializes persistent display
% % cumdisp(text); displays persistent text
% %
% persistent oldtxt;
% 
% if nargin<1,
%     oldtxt=''; 
%     fprintf(1,'\n'); 
% else
%     fprintf(1,[repmat('\b',[1,length(oldtxt)]),'%s'],txt);
%     oldtxt=sprintf('%s',txt);
% end
% end


			
