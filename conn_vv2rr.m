function [Z,zmat,xyz]=conn_vv2rr(ROI,varargin)
% internal function
%
% computes ROI-to-ROI matrix by averaging Voxel-to-Voxel correlations
% 
% conn_vv2rr(file, 'style','vv2rr','contrastsubjects',w1,'contrastconditions',w2)    : sum_(subjects&conditions) { w1(subject)*w2(condition)*  Q' * R[subject,condition] * Q }
% conn_vv2rr(file, 'style','vv2rv','contrastsubjects',w1,'contrastconditions',w2)    : sum_(subjects&conditions) { w1(subject)*w2(condition)*  Q' * R[subject,condition] }
% conn_vv2rr(file, 'style','vv2rv','contrastsubjects',W1,'contrastconditions',w2)    : sum_(subjects&conditions) { diag(W1(:,subject))*w2(condition)*  Q' * R[subject,condition] }
%
% Q is voxel weights vector/matrix in 'file' (e.g. input ROI file)
% R[subject,condition] is voxel-to-voxel correlation computed during denoising
% w(subject)
%

% note: call from within conn_process or conn_batch

global CONN_x;

options=struct(...
    'style','vv2rr',...             % vv2rr -> ROI-to-ROI matrix; vv2rv -> ROI-to-voxel matrix
    'isroi',true,...                % 1: treats input as ROI file (interprets 1:N 3d-file, and scales sum of values across voxels within each ROI to 1)
    'filepath',[],...               % [location of voxel-to-voxel data]
    'folderout',[],...              % [placeholder]
    'validsubjects',[],...          % selected subjects (default 1:N)
    'contrastsubjects',[],...       % weights for each subject (default ones(1,N)), alternatively one file for each subject (for voxel-specific weights)
    'validconditions',[],...        % selected conditions (default 1:N)
    'contrastconditions',[]);       % weights for each condition (default ones(1,N))

for n=1:2:numel(varargin)-1, assert(isfield(options,varargin{n}),'unrecognized option %s',varargin{n}); options.(varargin{n})=varargin{n+1}; end

nconditions=length(CONN_x.Setup.conditions.names)-1;
if isempty(options.style), style='vv2rr'; end
if isempty(options.filepath), options.filepath=CONN_x.folders.preprocessing; end
if isempty(options.validsubjects), options.validsubjects=1:CONN_x.Setup.nsubjects; end
if isempty(options.contrastsubjects), options.contrastsubjects=ones(1,numel(options.validsubjects))/numel(options.validsubjects); end
if isempty(options.validconditions), options.validconditions=1:nconditions; end
if isempty(options.contrastconditions), options.contrastconditions=ones(1,numel(options.validconditions))/numel(options.validconditions); end

icondition=[];isnewcondition=[];for ncondition=1:nconditions,[icondition(ncondition),isnewcondition(ncondition)]=conn_conditionnames(CONN_x.Setup.conditions.names{ncondition}); end
if any(isnewcondition(options.validconditions)), error(['Some conditions have not been processed yet. Re-run previous step']); end
% addnew=false;
if ischar(ROI), ROI=conn_fileutils('spm_vol',ROI); end
if ischar(options.contrastsubjects), options.contrastsubjects=cellstr(options.contrastsubjects); end
assert(isempty(options.contrastsubjects)||numel(options.contrastsubjects)==numel(options.validsubjects),'mismatched dimensions between options.contrastsubjects and validsubjects');

%h=conn_waitbar(0,'Extracting correlation matrix, please wait...');
lastmatdim=[];
Z=[];
xyz=[];
zmat=[];
for ivalidcondition=1:numel(options.validconditions),
    ncondition=options.validconditions(ivalidcondition);
    for ivalidsubject=1:numel(options.validsubjects)
        nsub=options.validsubjects(ivalidsubject);
        filename_B1=fullfile(options.filepath,['vvPC_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'.mat']);
        Y1=conn_vol(filename_B1);
        if ~isempty(zmat)&&~isequal(zmat,Y1.matdim.mat)&&numel(options.validsubjects)>1&&ivalidcondition==1, fprintf('warning, unequal orientation in vvPC files for subject %d\n',nsub); end
        zmat=Y1.matdim.mat;
        if isequal(lastmatdim, Y1.matdim)
        else
            if isstruct(ROI)||iscell(ROI), % 4d-file or 1:N 3d-file only; or {[3xN1], [3xN2], ...} (xyz coordinates of multiple ROIs)
                xyz=conn_convertcoordinates('idx2tal',1:prod(Y1.matdim.dim),Y1.matdim.mat,Y1.matdim.dim)';
                if isstruct(ROI)
                    W=conn_fileutils('spm_get_data',ROI,pinv(ROI(1).mat)*xyz);
                else
                    XYZ=cat(2,ROI{:});
                    if size(XYZ,1)==3, XYZ=[XYZ;ones(1,size(XYZ,2))]; end
                    i1=cell2mat(arrayfun(@(n)n+zeros(1,size(ROI{n},2)),1:numel(ROI),'uni',0));
                    XYZ=pinv(Y1.matdim.mat)*XYZ;
                    w=1; idx=1; for n=1:3, idx=idx+w*round(max(1,min(Y1.matdim.dim(n),XYZ(n,:)))-1); w=w*Y1.matdim.dim(n); end
                    W=zeros(numel(ROI),prod(Y1.matdim.dim));
                    W(i1+numel(ROI)*(idx-1))=1;
                end
                maxW=max(W(:));
                if options.isroi && maxW>1 && size(W,1)==1 && all(ismember(unique(W),0:maxW)), W=double(repmat(W,[maxW,1])==repmat((1:maxW)',[1,size(W,2)])); end
                w=W(:,Y1.voxels);
                if options.isroi, 
                    temp=w*xyz(1:3,Y1.voxels)';
                    sw=sum(w,2);
                    xyz=temp./repmat(max(eps,sw),1,3); 
                end
                lastmatdim=Y1.matdim;
            else W=ROI; % ROIs-by-voxels matrix (note: already in same space as vvPC_ files)
            end
            assert(isempty(W)|prod(Y1.matdim.dim)==size(W,2),'unequal number of voxels in ROI (%d) and VV (%d) files',size(W,2), prod(Y1.matdim.dim));
        end
        x=conn_get_volume(Y1);
        if isempty(W), w=[];
        else w=W(:,Y1.voxels);
        end
        if options.isroi&&~isempty(W), sw=sum(w,2); end
        if iscell(options.contrastsubjects)
            vol=spm_vol(options.contrastsubjects{ivalidsubject});
            assert(isequal(vol.dim,Y1.matdim.dim),'mismatched dimensions between contrastsubjects and vvPC files');
            assert(isequal(vol.mat,Y1.matdim.mat),'mismatched orientation between contrastsubjects and vvPC files');
            C=spm_read_vols(vol);
            C=C(Y1.voxels);
            if ~isempty(w), w=w.*repmat(C(:)',size(w,1),1); end
        elseif ~isempty(options.contrastsubjects), 
            C=options.contrastsubjects(ivalidsubject); 
            if ~isempty(w), w=w*C; end
        else C=1;
        end
        switch(options.style)
            case 'vv2rr'
                if isempty(Z), Z=zeros(size(W,1),size(W,1)); end
                temp=x*w';
                z=(temp'*temp);
                if options.isroi, Z=Z+options.contrastconditions(ivalidcondition)*z./max(eps,sw*sw'); % w*x'*x*w' / w*1*1'*w' (x: [Ncomp x Nvox], x'*x = R)
                else Z=Z+options.contrastconditions(ivalidcondition)*z; % w*x'*x*w' 
                end
            case 'vv2rv'
                if isempty(Z), Z=zeros([size(W,1),Y1.matdim.dim]); end
                temp=x*w';
                z=temp'*x;
                if options.isroi, Z(:,Y1.voxels)=Z(:,Y1.voxels)+diag(1./sw)*options.contrastconditions(ivalidcondition)*z;  % w*x'*x / w*1*1'
                else Z(:,Y1.voxels)=Z(:,Y1.voxels)+options.contrastconditions(ivalidcondition)*z;  % w*x'*x
                end
            case 'vv2vnorm'
                if isempty(Z), Z=zeros(Y1.matdim.dim); end
                nx=sum(x.^2,2);
                z=nx'*(x.^2);
                Z(Y1.voxels)=Z(Y1.voxels)+options.contrastsubjects(ivalidsubject)*options.contrastconditions(ivalidcondition)*z'; % sum-square of rows of x'*x
            %case 'vv2vnormr'
            %    if isempty(Z), Z=zeros([size(W,1),Y1.matdim.dim]); end
            %    nw=trace(w*w');
            %    temp=x*w';
            %    nx=temp*temp';
            %    z=sum(x.*(nx*x),1)/nw;
            %    Z(:,Y1.voxels)=Z(:,Y1.voxels)+options.contrastsubjects(ivalidsubject)*options.contrastconditions(ivalidcondition)*z; % sum-square of rows of x'*x
            otherwise
                error('unrecognized style %s (vv2rr or vv2rv options only)',options.style);
        end
        %conn_waitbar(((ivalidcondition-1)*CONN_x.Setup.nsubjects+nsub)/numel(options.validconditions)/CONN_x.Setup.nsubjects,h,sprintf('Subject %d Condition %d',nsub,ncondition));
        %n=n+1;
    end
%     if ~isempty(folderout)
%         filename=fullfile(folderout,['resultsROI_Condition',num2str(icondition(ncondition),'%03d'),'.mat']);
%         save(filename,'Z','names','names2','xyz'); %%%
%     end
end
%conn_waitbar('close',h);
% if addnew %%%
%     if ~isfield(CONN_x,'Analyses')||isempty(CONN_x.Analyses), tianalysis=1;
%     else tianalysis=numel(CONN_x.Analyses)+1;
%     end
%     CONN_x.Analyses(tianalysis).name=foldername;
%     CONN_x.Analysis=tianalysis;
%     if 1
%         CONN_x.Analyses(tianalysis).sourcenames={};
%         CONN_x.Analyses(tianalysis).sources=names;
%         CONN_x.Analyses(tianalysis).type=1;
%         conn save;
%         conn_msgbox({'ROI-to-ROI analyses created',['Go to ''second-level Results''-''>ROI-to-ROI'' and select analysis ',foldername]},'post-hoc analyses');
%     end
% end


end