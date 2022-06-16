function varargout = conn_mvpacreateh(option,varargin)

% internal function
%
% creates effect-size H maps
% eig | scores | map
% 

global CONN_x;
if nargin<2||isempty(option), option=''; end
options=struct(...
    'vvAnalysis',[],... % analysis number
    'RESULTSfolder',[],... % folder containing second-level results
    'ncondition',[],...  % condition number
    'valideig',[],...    % number of eigenpatterns
    'ROIfile',[],...
    'dogui',true);
for n=1:2:numel(varargin)-1, assert(isfield(options,varargin{n}),'unrecognized option %s',varargin{n}); options.(varargin{n})=varargin{n+1}; end

if isempty(options.vvAnalysis), options.vvAnalysis=CONN_x.vvAnalysis; end
if ischar(options.vvAnalysis), options.vvAnalysis=strmatch(options.vvAnalysis,{CONN_x.vvAnalyses.name},'exact'); assert(~isempty(options.vvAnalysis),'unable to find analysis name matching vvAnalysis field'); end
if ischar(options.ncondition), options.ncondition=strmatch(options.ncondition,CONN_x.Setup.conditions.names(1:end-1),'exact'); assert(~isempty(options.ncondition),'unable to find condition name matching ncondition field'); end

filepathresults=fullfile(CONN_x.folders.firstlevel_vv,CONN_x.vvAnalyses(options.vvAnalysis).name);
nconditions=length(CONN_x.Setup.conditions.names)-1; %%%
icondition=[];isnewcondition=[];for ncondition=1:nconditions,[icondition(ncondition),isnewcondition(ncondition)]=conn_conditionnames(CONN_x.Setup.conditions.names{ncondition}); end
if nnz(~isnewcondition)==1&&isempty(options.ncondition), options.ncondition=find(~isnewcondition); end
if isempty(options.ncondition),
    assert(options.dogui,'please enter condition of interest in ''ncondition'' field');
    options.ncondition=find(~isnewcondition);
    answ=listdlg('name','','PromptString','Select condition','ListString',CONN_x.Setup.conditions.names(options.ncondition),'SelectionMode','single','ListSize',[300 200]);
    if isempty(answ), return; end
    options.ncondition=options.ncondition(answ);
end
outcomenames=CONN_x.vvAnalyses(options.vvAnalysis).measures;
shownsources=find(ismember(conn_v2v('fieldtext',outcomenames,1),{'2'}));
outcomeisource=[];for n1=1:length(outcomenames),
    [outcomeisource(n1),isnew,outcomencompsource(n1)]=conn_v2v('match_extended',outcomenames{n1});
    if isnew, error('Measure %s not found in global measures list. Please re-run first-level analyses',outcomenames{n1}); end
end

varargout=cell(1,nargout);

switch(lower(option))
    case {'eig'}
    case {'scores'} % conn_mvpacreateh('scores',...)
        assert(~isempty(options.RESULTSfolder),'please enter group-level analysis folder in RESULTSfolder field');
        load(fullfile(options.RESULTSfolder,'SPM.mat'),'SPM');
        assert(isequal(SPM.xX_multivariate.Znames,arrayfun(@(n)sprintf('group-MVPA_%d',n),1:numel(SPM.xX_multivariate.Znames),'uni',0)),'sorry, effect-size extraction procedure only implemented for individual conditions');
        assert(size(SPM.xX_multivariate.C,1)==1,'sorry, effect-size extraction procedure only implemented for univariate tests');
        if isempty(options.valideig), options.valideig=numel(SPM.xX_multivariate.Znames); end
        
        % creates H_scores files
        if options.dogui, htfig=conn_msgbox('Loading connectivity values. Please wait','',-1); end
        for nsub=1:CONN_x.Setup.nsubjects
            if SPM.xX.SelectedSubjects(nsub)
                z=0;
                for ncomp=1:options.valideig
                    filename=fullfile(filepathresults,['BETA_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(options.ncondition),'%03d'),'_Measure',num2str(outcomeisource(ncomp),'%03d'),'_Component',num2str(outcomencompsource(ncomp),'%03d'),'.nii']);
                    [data,v]=conn_vol_read(filename);
                    z=z+data.*shiftdim(SPM.xX_multivariate.h(1,ncomp,:,:,:),2);
                end
                filename=fullfile(options.RESULTSfolder,['ALPHA_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(options.ncondition),'%03d'),'.nii']);
                conn_vol_write(filename,z,v);
                fprintf('saved file %s\n',filename);
            end
        end
    case {'map','mapprojected'}
        assert(~isempty(options.RESULTSfolder),'please enter group-level analysis folder in RESULTSfolder field');
        assert(~isempty(options.ROIfile),'please enter ROI file ROIfile field');
        load(fullfile(options.RESULTSfolder,'SPM.mat'),'SPM');
        if strcmpi(option,'mapprojected')
            filename={};
            for nsub=1:CONN_x.Setup.nsubjects
                if SPM.xX.SelectedSubjects(nsub)
                    filename{end+1}=fullfile(options.RESULTSfolder,['ALPHA_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(options.ncondition),'%03d'),'.nii']);
                end
            end
            [z,zmat]=conn_vv2rr(options.ROIfile,'style','vv2rv','filepath',CONN_x.folders.preprocessing,'validconditions',options.ncondition,'validsubjects',find(SPM.xX.SelectedSubjects),'contrastsubjects',filename);
        else
            X=SPM.xX_multivariate.X;
            C=SPM.xX_multivariate.C;
            alpha=X*pinv(X'*X)*C'; % alpha for N=inf
            [z,zmat]=conn_vv2rr(options.ROIfile,'style','vv2rv','filepath',CONN_x.folders.preprocessing,'validconditions',options.ncondition,'validsubjects',find(SPM.xX.SelectedSubjects),'contrastsubjects',alpha);
        end            
        for n=1:size(z,1)
            filename=conn_prepend('',options.ROIfile,sprintf('_effectsize_condition%d_roi%d.nii',options.ncondition,n));
            conn_vol_write(filename,shiftdim(z(n,:,:,:),1),zmat);
            fprintf('saved file %s\n',filename);
        end
    otherwise 
        error('unrecognized option %s (valid options are ''scores'' ''eig'' or ''map'')',option);
end

% 
% 
% [y,name,info]=conn_rex(char(filename),ROIname,'level','clusters');
% if options.dogui&&ishandle(htfig), delete(htfig); end
% name=regexprep(name,'^results\.ROIs\.?','cluster '); % {'label',...}
% xyz=info.ROIinfo.voxels; % {[nx3],...}
% xyz=cat(2,xyz{:});
% y=reshape(y,[size(filename), size(y,2)]); % subjects x components x ROIs  (eigenpattern scores)
% 
% if isempty(options.contrasteig)&&~isempty(options.RESULTSfolder)
%     CiX=SPM.xX_multivariate.C*pinv(SPM.xX_multivariate.X);
%     for n=1:size(y,3), options.contrasteig(:,n)=CiX*y(:,:,n); end
% end
% 
% for nroi=1:size(y,3) % ROIs
%     Z=0;
%     for ncomp=1:size(y,2) % components
%         [z,zmat]=conn_vv2rr({xyz{nroi}'},'style','vv2rv','validconditions',options.ncondition,'contrastsubjects',y(:,ncomp,nroi)'); 
%         z=shiftdim(z,1);
%         filename=conn_prepend('',ROIname,sprintf('_cluster%d_eig%d.nii',nroi,ncomp));
%         conn_vol_write(filename,z,zmat);
%         fprintf('created file %s\n',filename);
%         if ~isempty(options.contrasteig), Z=Z+options.contrasteig(ncomp,min(size(options.contrasteig,2),nroi))*z; end
%     end
%     if ~isempty(options.contrasteig), 
%         filename=conn_prepend('',ROIname,sprintf('_cluster%d_effect.nii',nroi));
%         conn_vol_write(filename,Z,zmat);
%         fprintf('created file %s\n',filename);
%     end
% end 

%         fh=conn_mesh_display(filename);
%         fh('background',[1 1 1]);
%         fh('brain',4);
%         fh('colormap',parula);
%         fh('material',[]);
%         fh('view',[-1,0,0],[],-1); fh('print',1,conn_prepend('',filename,'_left.jpg'),'-nogui');
%         fh('view',[1,0,0],[1,0,.5],-1); fh('print',1,conn_prepend('',filename,'_leftmedial.jpg'),'-nogui');
%         fh('close');        


