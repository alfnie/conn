function fileout=evlab17_run_roiresults(varargin)
% EVLAB17_RUN_ROIRESULTS runs second-level ROI-based analyses
%
%  OPTIONS (entered in .cfg file with fieldnames preceeded by #, or entered as argument pairs to evlab17_run_results)
%
%      folder          : folder where original second-level analysis are stored (see evlab17_run_results)
%      rois            : list of ROI files (note: if roi_ss=1 enter first all files -one per subject- for ROI#1, followed by all files for ROI#2, etc.)
%      roi_measure     : (optional) measure computed within each ROI {'mean'},'eigenvariate','weighted mean','weighted eigenvariate','median','sum','weighted sum','count','max',min' 
%      roi_threshold   : (optional) absolute-threshold defining voxels within each ROI [0]
%      roi_labels      : (optional) list of ROIs to be included in analyses (by default when entering atlas files defining multiple ROIs all ROIs in the atlas are included in the analyses) ({})
%      roi_ss          : (optional) 1/0 indicating whether ROIs are subject-specific [0]
%      roi_id          : (optional) ROI name (note: output filename will be [folder]/REX_[roi_id].mat file)
%      roi_groups      : (optional) list of ROIs to group together as a single ROI in analyses (by default each ROI is analyzed separately) For multiple groups separate the list of ROI indexes/labels by the label of each new group-ROI (e.g. group1 roi1 roi2 group2 roi3 roi4 roi5) ({})
%
% note: if multiple .cfg files are entered, second-level analyses will include all of the ROIs defined in each .cfg file. 
%       and the output filename will use the last 'roi_id' field, if any is present within all .cfg files. 
%

evlab17_module init silent;
fileout=[];

% loads .cfg files
options={};ids={};
if nargin>1&&ischar(varargin{1})&&isempty(dir(varargin{1}))
    n=0;
else
    for n=1:nargin
        if isempty(varargin{n}), break; end
        filenames=varargin{n};
        if ~iscell(filenames), filenames={filenames}; end
        for n1=1:numel(filenames)
            filename=filenames{n1};
            if ischar(filename)
                if isempty(regexpi(filename,'<subjectids?>'))&&isempty(dir(filename)),
                    if ~isempty(dir(fullfile(fileparts(which(mfilename)),filename))), filename=fullfile(fileparts(which(mfilename)),filename);
                    else
                        fprintf('warning: file %s not found\n',filename);
                        filename=which(filename);
                    end
                end
                [nill,nill,fext]=fileparts(filename);
                if ~strcmp(fext,'.cfg'), filename=struct('rois',{cellstr(filename)}); % if not .cfg file assumes user is entering individual roi filename directly
                else fprintf('loading file %s\n',filename);
                end
            end
            options{end+1}=conn_loadcfgfile(filename);
        end
    end
end
if isempty(options), options={struct}; end
for n=n+1:2:nargin-1
    fieldname=regexp(varargin{n},'\.','split');
    fieldvalue=varargin{n+1};
    for n1=1:numel(options), options{n1}=setfield(options{n1},fieldname{:},fieldvalue); end
end
for n=1:numel(options),
    fprintf('%s options set #%d:\n',mfilename,n);
    disp(options{n}); 
end
options0=options;
options0_arguments=varargin;

% interprets info
results_name=sprintf('REX_%s.mat',datestr(now,'yyyy_mm_dd_HHMMSSFFF'));
roidata={}; roinames={}; SPMsources=[];
for n=1:numel(options)
    assert(isfield(options{n},'rois'),'missing #rois information in set %d',n);
    assert(isfield(options{n},'folder'),'missing #folder information in set %d',n);
    results_path=char(options{n}.folder);
    cROIfile=cellstr(options{n}.rois);

    fields={'select_clusters',0};
    names=fieldnames(options{n});
    listallowed={'rois','folder','roi_id','roi_keys','roi_groups','data','design','contast_between','contrast_within','contrast_names','mask'};
    for n1=1:numel(names), % all other options are passed directly to rex
        if ~ismember(names{n1},listallowed), 
            if strcmp(names{n1},'roi_measure'), tnames='summary_measure'; 
            elseif strcmp(names{n1},'roi_labels'), tnames='selected_clusters'; 
            else tnames=names{n1};
            end
            fields{end+1}=tnames; fields{end+1}=options{n}.(names{n1}); 
        end
    end
    if (isfield(options{n},'roi_ss')&&options{n}.roi_ss)||any(cellfun('length',regexpi(cROIfile,'<subjectids?>')))||any(~cellfun('length',regexp(cROIfile,'[\\\/\:]'))),
        if isempty(SPMsources)
            assert(conn_existfile(fullfile(results_path,'SPM.mat')),'unable to find file %s. Please re-run original second-level analyses',fullfile(results_path,'SPM.mat'));
            load(fullfile(results_path,'SPM.mat'),'SPM');
            if isfield(SPM,'xX_multivariate')&&isfield(SPM.xX_multivariate,'Zfiles')
                SPMsources=reshape(SPM.xX_multivariate.Zfiles,size(SPM.xX_multivariate.X,1),[]);
            else
                SPMsources=reshape({SPM.xY.VY(:).fname},size(SPM.xX.X,1),[]);
            end
        end
        Nscans=size(SPMsources,1);
        if any(cellfun('length',regexpi(cROIfile,'<subjectids?>')))||any(~cellfun('length',regexp(cROIfile,'[\\\/\:]'))),
            assert(isfield(SPM.xX,'SubjectIDs'),'unable to find SubjectID information in %s. Please re-run original second-level analyses',fullfile(results_path,'SPM.mat'));
            assert(numel(cellstr(SPM.xX.SubjectIDs))==Nscans,'unable to find SubjectID information in %s. Please re-run original second-level analyses',fullfile(results_path,'SPM.mat'));
            subjectids=regexprep(cellstr(SPM.xX.SubjectIDs),'^sub-','');
            cROIfile=repmat(cROIfile(:)',[Nscans,1]);
            for n1=1:size(cROIfile,2), 
                if ~isempty(regexpi(cROIfile{1,n1},'<subjectids?>')),
                    cROIfile(:,n1)=reshape(cellfun(@(x)regexprep(cROIfile{1,n1},'<subjectids?>',x,'ignorecase'),subjectids,'uni',0),[],1);
                elseif isfield(options{n},'roi_keys')
                    if iscell(options{n}.roi_keys)
                        for nr=1:numel(options{n}.roi_keys)
                            cROIfile(:,n1)=reshape(cellfun(@(x)options{n}.roi_keys{nr}(cROIfile{1,n1},x),subjectids,'uni',0),[],1);
                            if all(conn_existfile(cROIfile(:,n1))), break; end
                        end
                    else
                        cROIfile(:,n1)=reshape(cellfun(@(x)options{n}.roi_keys(cROIfile{1,n1},x),subjectids,'uni',0),[],1);
                    end
                else error('unable to match ROI id %s. Please try entering full paths to ROI files instead',cROIfile{1,n1});
                end
            end
        else subjectids={};
        end
        for n1=1:numel(cROIfile), 
            assert(conn_existfile(cROIfile{n1}),'unable to find ROI file %s',cROIfile{n1});
            [nill,nill,cROIfile{n1}]=conn_file(cROIfile{n1}); 
        end
        assert(~rem(numel(cROIfile),Nscans),'unexpected number of subject/scan-specific ROI files (%d). Expected a multiple of the number of subjects/scans in this analysis (%d)',numel(cROIfile),Nscans);
        cROIfile=reshape(cROIfile,Nscans,[]);
        tdata=[]; % conditions x ROIs x subjects
        for n1=1:Nscans
            fprintf('.');
            [ttdata,tnames,params]=conn_rex(char(SPMsources(n1,:)),char(cROIfile(n1,:)),'output_type','none','level','clusters_nospatial','gui',0,fields{:});
            ttdata(isnan(ttdata))=0; % missing-data
            if isfield(options{n},'roi_id'), vnames=1:numel(tnames); for n2=1:size(cROIfile,2), if isempty(vnames), break; end; [nill,tstr,nill]=fileparts(cROIfile{n1,n2}); match=strncmp(tstr,tnames(vnames),numel(tstr)); if any(match), tnames(vnames(match))=cellfun(@(x)[options{n}.roi_id,x(numel(tstr)+1:end)],tnames(vnames(match)),'uni',0); vnames(match)=[]; end; end
            elseif numel(subjectids)>=n1, tnames=regexprep(tnames,['_?',subjectids{n1}],''); 
            end
            if n1==1, tnames0=tnames;
            elseif ~isequal(tnames,tnames0), 
                tmissing=~ismember(tnames0,tnames);
                tnew=~ismember(tnames,tnames0);
                if any(tmissing), ttdata=cat(2,ttdata,zeros(size(ttdata,1),nnz(tmissing))); tnames=[tnames,tnames0(tmissing)]; fprintf('warning: missing ROIs %s in subject/scan #%d\n',sprintf('%s ',tnames0{tmissing}),n1); end
                if any(tnew), tdata=cat(2,tdata,zeros(size(tdata,1),nnz(tnew),size(tdata,3))); tnames0=[tnames0,tnames(tnew)]; fprintf('warning: new ROIs %s in subject/scan #%d\n',sprintf('%s ',tnames{tnew}),n1); end
                [ok,idx]=ismember(tnames0,tnames);
                ttdata=ttdata(:,idx);
                tnames=tnames(idx);
            end
            if n1==1, tdata=ttdata;
            elseif ~isequal(size(tdata(:,:,1)),size(ttdata)), error('dimensions mismatch: %d samples and %d regions for subject/scan #%d, %d samples and %d regions for subject/scan #1',size(ttdata,1),size(ttdata,2),n1,size(tdata,1),size(tdata,2)); 
            else tdata(:,:,n1)=ttdata;
            end
        end
        fprintf('\n');
        roidata{n}=reshape(permute(tdata,[3,1,2]),[],size(tdata,2)); % (subjects x conditions) x ROIs
        roinames{n}=tnames;
        SPMinfo=struct('SPM',SPM);
    else
        for n1=1:numel(cROIfile), [nill,nill,cROIfile{n1}]=conn_file(cROIfile{n1}); end
        [roidata{n},roinames{n},params]=conn_rex(fullfile(results_path,'SPM.mat'),char(cROIfile),'output_type','none','level','clusters_nospatial','gui',0,fields{:});
        roidata{n}(isnan(roidata{n}))=0; % missing-data
        if isfield(options{n},'roi_id'), vnames=1:numel(roinames); for n1=1:numel(cROIfile), if isempty(vnames), break; end; [nill,tstr,nill]=fileparts(cROIfile{n1}); match=strncmp(tstr,roinames(vnames),numel(tstr)); if any(match), roinames(vnames(match))=cellfun(@(x)[options{n}.roi_id,x(numel(tstr)+1:end)],roinames(vnames(match)),'uni',0); vnames(match)=[]; end; end; end
        SPMinfo=params.SPM;
    end
    if isfield(options{n},'roi_id'), results_name=sprintf('REX_%s.mat',options{n}.roi_id); end
    if isfield(options{n},'roi_groups'), 
        groups=cellstr(options{n}.roi_groups);
        groupname={'grouproi'};
        grouprois={[]};
        roinames=[cat(2,roinames{:})];
        roidata=cat(2,roidata{:});
        fprintf('available ROI labels: %s\n',sprintf('%s ',roinames{:}));
        for n1=1:numel(groups)
            idxc=find(strcmp(groups{n1},roinames));
            if numel(idxc)~=1, idxc=find(strncmp(groups{n1},roinames,numel(groups{n1}))); end
            if numel(idxc)~=1, idxc=find(strncmp(groups{n1},regexprep(roinames,'^.*\.',''),numel(groups{n1}))); end
            if numel(idxc)==1, grouprois{end}=[grouprois{end} idxc]; 
            elseif isempty(idxc), 
                fprintf('ROI-label %s not found. Assuming label of new group\n',groups{n1}); 
                groupname{end+1}=groups{n1}; 
                grouprois{end+1}=[];
            else fprintf('warning: multiple (%d) possible ROI-label matches found for %s in %s. Selecting first\n',numel(idxc),groups{n1},sprintf('%s ',roinames{idxc})); grouprois{end}=[grouprois{end} idxc(1)]; 
            end
        end
        newroidata=roidata(:,setdiff(1:numel(roinames),[grouprois{:}]));
        newroinames=roinames(setdiff(1:numel(roinames),[grouprois{:}]));
        for n1=1:numel(groupname)
            if ~isempty(grouprois{n1})
                temp=cat(2,roidata(:,grouprois{n1}));
                newroidata=[newroidata mean(temp,2)];
                newroinames=[newroinames groupname(n1)];
            end
        end
        roidata={newroidata};
        roinames={newroinames};
    end
end
params=struct(...
    'output_type','saverex',...
    'output_rex',results_name,...
    'output_folder',results_path,...
    's',[],...
    'spm_file',fullfile(results_path,'SPM.mat'),...
    'SPM',SPMinfo,...
    'ROIdata',cat(2,roidata{:}),...
    'ROInames',{cat(2,roinames{:})});
save(fullfile(results_path,results_name),'params');
fprintf('ROI-level statistics saved in %s\n',fullfile(results_path,results_name));
% save & display results
conn_rex('results',fullfile(results_path,results_name),'output_type','saverex'); 
end




