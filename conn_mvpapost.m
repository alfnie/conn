function varargout = conn_mvpapost(ROIname,RESULTSfolder,varargin)

% internal function
%
% computes MVPA post-hoc effects
% 

global CONN_x;

options=struct(...
    'ncondition',[],...
    'valideig',[],...
    'contrasteig',[],... % neig x nrois
    'dogui',true);

if nargin<2, RESULTSfolder=[]; end    
for n=1:2:numel(varargin)-1, assert(isfield(options,varargin{n}),'unrecognized option %s',varargin{n}); options.(varargin{n})=varargin{n+1}; end

filepathresults=fullfile(CONN_x.folders.firstlevel_vv,CONN_x.vvAnalyses(CONN_x.vvAnalysis).name);
nconditions=length(CONN_x.Setup.conditions.names)-1; %%%
icondition=[];isnewcondition=[];for ncondition=1:nconditions,[icondition(ncondition),isnewcondition(ncondition)]=conn_conditionnames(CONN_x.Setup.conditions.names{ncondition}); end
if nnz(~isnewcondition)==1&&isempty(options.ncondition), options.ncondition=find(~isnewcondition); end
if isempty(options.ncondition),
    options.ncondition=find(~isnewcondition);
    answ=listdlg('name','','PromptString','Select condition','ListString',CONN_x.Setup.conditions.names(options.ncondition),'SelectionMode','single','ListSize',[300 200]);
    if isempty(answ), return; end
    options.ncondition=options.ncondition(answ);
end


outcomenames=CONN_x.vvAnalyses(CONN_x.vvAnalysis).measures;
shownsources=find(ismember(conn_v2v('fieldtext',outcomenames,1),{'2'}));
outcomeisource=[];for n1=1:length(outcomenames),
    [outcomeisource(n1),isnew,outcomencompsource(n1)]=conn_v2v('match_extended',outcomenames{n1});
    if isnew, error('Measure %s not found in global measures list. Please re-run first-level analyses',outcomenames{n1}); end
end

if isempty(options.contrasteig)&&~isempty(RESULTSfolder)
    load(fullfile(RESULTSfolder,'SPM.mat'),'SPM');
    %assert(isequal(SPM.xX_multivariate.Znames,arrayfun(@(n)sprintf('group-MVPA_%d',n),1:numel(SPM.xX_multivariate.Znames),'uni',0)),'sorry, effect-size extraction procedure only implemented for individual conditions');
    assert(size(SPM.xX_multivariate.C,1)==1,'sorry, effect-size extraction procedure only implemented for univariate tests');
    if isempty(options.valideig), options.valideig=1:numel(SPM.xX_multivariate.Znames); end
elseif isempty(options.valideig), options.valideig=1:numel(outcomeisource); 
end
assert(isempty(options.contrasteig)|numel(options.contrasteig)==numel(options.valideig),'mismatch between contrasteig and valideig fields');
varargout=cell(1,nargout);

% extract eigenpattern scores from each ROI
filename=cell(CONN_x.Setup.nsubjects,length(options.valideig));
for nsub=1:CONN_x.Setup.nsubjects
    for ncomp=1:length(options.valideig)
        filename{nsub,ncomp}=fullfile(filepathresults,['BETA_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(options.ncondition),'%03d'),'_Measure',num2str(outcomeisource(options.valideig(ncomp)),'%03d'),'_Component',num2str(outcomencompsource(options.valideig(ncomp)),'%03d'),'.nii']);
        %filename{nsub,ncomp}=fullfile(filepathresults,sprintf('BETA_Subject%03d_Condition%03d_Measure%03d_Component%03d.nii',nsub,icondition(options.ncondition),outcomeisource(options.valideig(ncomp)),outcomencompsource(options.valideig(ncomp))));
    end
end
if options.dogui, htfig=conn_msgbox('Loading connectivity values. Please wait','',-1); end
[y,name,info]=conn_rex(char(filename),ROIname,'level','clusters');
if options.dogui&&ishandle(htfig), delete(htfig); end
name=regexprep(name,'^results\.ROIs\.?','cluster '); % {'label',...}
xyz=info.ROIinfo.voxels; % {[nx3],...}
xyz=cat(2,xyz{:});
y=reshape(y,[size(filename), size(y,2)]); % subjects x components x ROIs  (eigenpattern scores)

if isempty(options.contrasteig)&&~isempty(RESULTSfolder)
    CiX=SPM.xX_multivariate.C*pinv(SPM.xX_multivariate.X);
    for n=1:size(y,3), options.contrasteig(:,n)=CiX*y(:,:,n); end
end

for nroi=1:size(y,3) % ROIs
    Z=0;
    for ncomp=1:size(y,2) % components
        [z,zmat]=conn_vv2rr({xyz{nroi}'},'style','vv2rv','validconditions',options.ncondition,'contrastsubjects',y(:,ncomp,nroi)'); 
        z=shiftdim(z,1);
        filename=conn_prepend('',ROIname,sprintf('_cluster%d_eig%d.nii',nroi,ncomp));
        conn_vol_write(filename,z,zmat);
        fprintf('created file %s\n',filename);
        if ~isempty(options.contrasteig), Z=Z+options.contrasteig(ncomp,min(size(options.contrasteig,2),nroi))*z; end
    end
    if ~isempty(options.contrasteig), 
        filename=conn_prepend('',ROIname,sprintf('_cluster%d_effect.nii',nroi));
        conn_vol_write(filename,Z,zmat);
        fprintf('created file %s\n',filename);
    end
end 

%         fh=conn_mesh_display(filename);
%         fh('background',[1 1 1]);
%         fh('brain',4);
%         fh('colormap',parula);
%         fh('material',[]);
%         fh('view',[-1,0,0],[],-1); fh('print',1,conn_prepend('',filename,'_left.jpg'),'-nogui');
%         fh('view',[1,0,0],[1,0,.5],-1); fh('print',1,conn_prepend('',filename,'_leftmedial.jpg'),'-nogui');
%         fh('close');        


