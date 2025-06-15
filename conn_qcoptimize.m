
function varargout=conn_qcoptimize(varargin)

global CONN_x CONN_gui;

options=struct(...
    'removingoutliersubjects',true,...                                                  % 1/0 considers removing subjects
    'changingwmcsfcomponents',true,...                                                  % 1/0 considers changing the number of aCompCor components during denoising
    'usingglobalsignal',false,...                                                       % 1/0 considers using global signal regression (GSR)
    'exhaustivesearch',0,...                                                            % +1/-1/0 exhaustive search outlier subjects (+1 to sequentially find the optimal subject to remove; -1 to sequentially find the optijmal subject to keep)
    'minoutliers',[],...                                                                % consider removing at least this total number of subjects (leave empty to always remove subjects with "extreme outlier" scores, i.e. subjects with QC_OutlierScore>3 scores)
    'maxoutliers',[],...                                                                % consider removing up to this total number of subjects (leave empty to consider removing subjects up to QC_OutlierScore>1 scores)
    'maxwmcsfcomponents',10,...                                                         % considers using up to this number of aCompCor components (default 10)
    'maxwmcsfcomponentsdiff',10,...                                                      % considers only "similar" number of WM & CSF components (when difference abs(WM-CSF) is below maxwmcsfcomponentsdiff)
    'initwmcsfcomponents',[0 0],...                                                     % starting position for number of WM&CSF components
    'initoutliers',0,...                                                                % starting position for number of outlier subjects
    'initgsr',0,...                                                                     % starting position for global signal regression
    'initnpass',0,...                                                                   % starting pass index (0%3 for aCompCor, 1%3 for GSR, 2%3 for outliers)
    'foldername', 'qcoptimize',...                                                      % output QA-report name where results will be saved
    'scorevariable','QC_OutlierScore',...                                               % 2nd-level covariate name containing outlier score for each subject (larger numbers = worse)
    'outliersubjects',[],...                                                            % 2nd-level covariate name containing invalid subjects (1:subjects to exclude; 0:subjects to keep) (when outliersubjects is used, scorevariable is disregarded)
    'qcvariables',{{'QC_InvalidScans','QC_MeanMotion','QC_ProportionValidScans'}},...   % name of QC variables to use for QC-FC plots (only needed if DOSITES=true)
    'controlvariables',{{}},...                                                         % name of control variables to used in QC-FC correlations 
    'projectname',[],...                                                                % CONN project name (by default it will work on the currently-open CONN project)
    'numberofpasses',2);                                                                % number of line-search optimization passes (when optimizing both the number of subjects and the number of WMCSF components)
for n=1:2:nargin-1,
    assert(isfield(options,lower(varargin{n})),'unrecognize option %s',varargin{n});
    options.(lower(varargin{n}))=varargin{n+1};
end

if ~isempty(options.projectname)&&any(conn_server('util_isremotefile',options.projectname)), [varargout{1:nargout}]=conn_server('run',mfilename,varargin{:},'projectname',conn_server('util_localfile',options.projectname)); return;
elseif conn_projectmanager('inserver'), 
    [hmsg,hstat]=conn_msgbox({'Process running remotely','Please wait...',' ',' '},[],[],true);
    [varargout{1:nargout}]=conn_server('run_withwaitbar',hstat,mfilename,varargin{:}); 
    if ~isempty(hmsg)&&ishandle(hmsg), delete(hmsg); end
    return; 
end

if ~isempty(options.projectname)
    conn('load',options.projectname);
end

[filepath,filename,fileext] = fileparts(conn_projectmanager('projectfile')); % project folder
filepath = conn_fullfile(filepath,filename);
if options.removingoutliersubjects
    try
        [covs,covs_names]=conn_module('get','l2covariates');
        if ~isempty(options.outliersubjects), 
            idxscore=find(ismember(covs_names,options.outliersubjects)); 
            if isempty(options.minoutliers), options.minoutliers=sum(isnan(covs(:,idxscore))|covs(:,idxscore)>0); end % remove (at least) any subject marked as invalid
            if isempty(options.maxoutliers), options.maxoutliers=sum(isnan(covs(:,idxscore))|covs(:,idxscore)>0); end % remove (at most) any subject marked as invalid
            options.initoutliers=options.minoutliers; 
        else 
            idxscore=find(ismember(covs_names,options.scorevariable)); % list of subjects to keep is derived from ranks of options.scorevariable (for non-exhaustive search only)
            if isempty(options.minoutliers), options.minoutliers=sum(isnan(covs(:,idxscore))|covs(:,idxscore)>3); end % remove (at least) any subject with scores above 3 (extreme outliers)
            if isempty(options.maxoutliers), options.maxoutliers=sum(isnan(covs(:,idxscore))|covs(:,idxscore)>1); end % remove (at most) any subject with scores above 1
        end
        [nill,idxsubjects]=sort(covs(:,idxscore)); % note: nan's appear last
    end
end
if isempty(options.initoutliers), options.initoutliers=0; end
if options.changingwmcsfcomponents
    iWM=find(strcmp(CONN_x.Preproc.confounds.names,'White Matter'),1);
    iCSF=find(strcmp(CONN_x.Preproc.confounds.names,'CSF'),1);
    assert(~isempty(iWM),'unable to find White Matter component in list of confounds. Please add ''White Matter'' to the list of confounds in the Denoising tab and try again');
    assert(~isempty(iCSF),'unable to find CSF component in list of confounds. Please add ''CSF'' to the list of confounds in the Denoising tab and try again');
end
if options.usingglobalsignal
    iGM=find(strcmp(CONN_x.Preproc.confounds.names,'Grey Matter'),1);
    assert(~isempty(iGM),'unable to find Grey Matter component in list of confounds. Please add ''Grey Matter'' to the list of confounds in the Denoising tab and try again');
end


PAR=[]; % information for each trial (#WMcomponents #CSFcomponents #outliersubjects)
POV=[]; % minimum percent overlap measure among tested QC measures for each trial
if nnz([options.removingoutliersubjects,options.changingwmcsfcomponents,options.usingglobalsignal])>1, NPASS=options.initnpass+(0:3*options.numberofpasses-1);
else NPASS=options.initnpass+(0:2);
end
OPT_NDIM={options.initwmcsfcomponents(1) options.initwmcsfcomponents(2)}; % # of WM&CSF dimensions
OPT_NREMOVE=options.initoutliers; % # of outlier subjects to remove
OPT_NGLOBAL=options.initgsr; % # of GSR components
LIST_IKEEP=cell(1,CONN_x.Setup.nsubjects); % list of subjects to keep (for exhaustivesearch>0) or list of subjects to remove (for exhaustivesearch<0) (note:disregarded for exhaustivesearch=false)
for npass=NPASS    
    switch(rem(npass,3))
        case 0,
            if options.changingwmcsfcomponents, % 0%3 npass used to optimize aCompCor
                NREMOVE=OPT_NREMOVE;
                NDIM={0:options.maxwmcsfcomponents, 0:options.maxwmcsfcomponents};
                NGLOBAL=OPT_NGLOBAL;
            else continue; % on even npass if not optimizing aCompCor
            end
        case 1,  % 1%3 npass used to optimize GSR
            if options.usingglobalsignal, 
                NREMOVE=OPT_NREMOVE;
                NDIM=OPT_NDIM;
                NGLOBAL=0:1;
            else continue;
            end
        case 2,
            if options.removingoutliersubjects, % 2%3 used to optimize outliersubjects
                if options.exhaustivesearch<0, NREMOVE=options.maxoutliers:-1:options.minoutliers;
                else NREMOVE=options.minoutliers:options.maxoutliers;
                end
                NDIM=OPT_NDIM;
                NGLOBAL=OPT_NGLOBAL;
                %LIST_IKEEP=cell(1,CONN_x.Setup.nsubjects); % note: uncomment this line to reset index of outlier subjects in each pass
            else continue; % on odd npass if not optimizing outliersubjects
            end
    end % note: consider enabling several parameters simultaneoulsy for optimization over entire parameter space beyond linesearch (takes longer but produces always equal-or-better results)
    for nglobal=NGLOBAL,
        for ndim1=NDIM{1},
            for ndim2=NDIM{2},
                if options.changingwmcsfcomponents&&~isempty(options.maxwmcsfcomponentsdiff)&&abs(ndim1-ndim2)>options.maxwmcsfcomponentsdiff, continue; end
                if options.changingwmcsfcomponents,
                    CONN_x.Preproc.confounds.dimensions{iWM}(1)=ndim1;
                    CONN_x.Preproc.confounds.dimensions{iCSF}(1)=ndim2;
                end
                if options.usingglobalsignal,
                    CONN_x.Preproc.confounds.dimensions{iGM}(1)=nglobal;
                end
                %if (options.changingwmcsfcomponents||options.usingglobalsignal)&&conn_projectmanager('inserver'), conn save; end; % note: when working on a remote project
                addoptions={};
                if ~isempty(options.controlvariables), addoptions=[addoptions, {'QA.controlcovariates',options.controlvariables}]; end
                for nremove=NREMOVE,
                    tname=sprintf('QA_%s_%03d_%03d_%03d',options.foldername,ndim1,ndim2,nremove);
                    if options.removingoutliersubjects
                        if nremove>0&&options.exhaustivesearch&&isempty(LIST_IKEEP{nremove}) % updates list of subjects to keep/remove
                            if options.exhaustivesearch>0 % try removing one subject from the optimal result when removing nremove-1 subjects
                                if nremove==1, lastkeep=1:CONN_x.Setup.nsubjects;
                                else lastkeep=LIST_IKEEP{nremove-1}; 
                                end
                                if isempty(lastkeep) % initialize to OulierScores rank if starting exhaustive search from a number of outliers above 0
                                    lastkeep=sort(idxsubjects(1:end-nremove+1));
                                end
                            else % try adding one subject from the optimal result when removing nremove+1 subjects
                                lastkeep=LIST_IKEEP{nremove+1}; 
                                if isempty(lastkeep) % initialize to OulierScores rank if starting exhaustive search from a number of outliers above 0
                                    lastkeep=sort(idxsubjects(end-nremove:end));
                                end
                            end
                            tPOV=[];
                            for n1=1:numel(lastkeep)
                                keep=lastkeep; keep(n1)=[];
                                if options.exhaustivesearch<0, keep=setdiff(1:CONN_x.Setup.nsubjects,keep); end 
                                conn_batch(...
                                    'QA.plots',         {'QA_DENOISE FC-QC'},...
                                    'QA.l2covariates',  options.qcvariables,...
                                    'QA.foldername',    tname,...
                                    addoptions{:},...
                                    'subjects',keep);
                                tPOV(n1)=min(cellfun(@(name)conn_jsonread(fullfile(filepath,'results','qa',tname,sprintf('QA_DENOISE_QC-FC.measure%s.json',name)),'AfterDenoising_PercentOverlap'),options.qcvariables));
                                conn_disp('fprintf','exhaustive search: test %d_%d optimal overlap %.1f\n',nremove,n1,100*max(tPOV));
                            end
                            [tPOVmax,idx]=max(tPOV);
                            keep=lastkeep; keep(idx)=[];
                            LIST_IKEEP{nremove}=keep;
                            conn_disp('fprintf','CONSIDERING REMOVING THE FOLLOWING SUBJECTS (%d): %s\n',nremove,mat2str(setdiff(1:CONN_x.Setup.nsubjects,keep)));
                        end
                        if nremove>0
                            if options.exhaustivesearch>0, keep=LIST_IKEEP{nremove};
                            elseif options.exhaustivesearch<0, keep=setdiff(1:CONN_x.Setup.nsubjects,LIST_IKEEP{nremove});
                            else keep=sort(idxsubjects(1:end-nremove));
                            end
                            addoptions=[addoptions, {'subjects',keep}];
                        end
                    end
                    conn_batch(...
                        'QA.plots',         {'QA_DENOISE FC-QC'},...
                        'QA.l2covariates',  options.qcvariables,...
                        'QA.foldername',    tname,...
                        addoptions{:});
                    PAR=[PAR; nremove ndim1 ndim2 nglobal];
                    POV=[POV; min(cellfun(@(name)conn_jsonread(fullfile(filepath,'results','qa',tname,sprintf('QA_DENOISE_QC-FC.measure%s.json',name)),'AfterDenoising_PercentOverlap'),options.qcvariables))];

                    str=''; if options.removingoutliersubjects, str=[str sprintf(' #remove = %d ',PAR(end,1))]; end; if options.changingwmcsfcomponents, str=[str sprintf(' #WM = %d  #CSF = %d ',PAR(end,2),PAR(end,3))]; end; if options.usingglobalsignal, str=[str sprintf(' #GSR = %d ',PAR(end,4))]; end;
                    conn_disp('fprintf','%s  ====> Percent Match = %.1f%%\n',str,100*POV(end));
                    [nill,idx]=max(POV);
                    str=''; if options.removingoutliersubjects, str=[str sprintf(' #remove = %d ',PAR(idx,1))]; end; if options.changingwmcsfcomponents, str=[str sprintf(' #WM = %d  #CSF = %d ',PAR(idx,2),PAR(idx,3))]; end; if options.usingglobalsignal, str=[str sprintf(' #GSR = %d ',PAR(idx,4))]; end;
                    conn_disp('fprintf','OPTIMAL COMBINATION SO FAR: %s  ====> Percent Match = %.1f%%\n',str,100*POV(idx));
                end
            end
        end
    end
    [nill,idx]=max(POV);
    OPT_NREMOVE=PAR(idx,1);
    OPT_NDIM={PAR(idx,2) PAR(idx,3)}; 
    OPT_NGLOBAL=PAR(idx,4);
end
varargout={PAR,POV,LIST_IKEEP};
