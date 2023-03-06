function [Z,zmat,xyz]=conn_vv2rr(ROI,varargin)
% CONN_VV2RR computes ROI-to-ROI matrix by averaging Voxel-to-Voxel correlations across voxel-pairs 
%
% data = conn_vv2rr(ROIfile)
%   data : [Nrois x Nrois x Nsub x Ncond] average correlation matrix for each subject and condition 
% 
% data = conn_vv2rr(ROIfile, option1, value1, option2, value2, ...) uses alternative analysis options:
%   'style'               : 'vv2rr' outputs an ROI-to-ROI matrix (default 'vv2rr')  
%                           'vv2rv' outputs a seed-to-voxel image per ROI (i.e. aggregates only across one dimension)
%   'saveas'              : filename of output NIFTI file (default '', output not saved to file)
%   'validsubjects'       : list of subject indexes to include (default all subjects in current project) 
%   'contrastsubjects'    : contrast across subjects, output will be evaluated separately for each row in this contrast matrix (default eye(numel(validsubjects)), one contrast per subject )
%   'validconditions'     : list of condition indexes to include (default all conditions in current project) 
%   'contrastconditions'  : contrast across conditions, output will be evaluated separately for each row in this contrast matrix (default eye(numel(validconditions)), one contrast per condition )
%
% for option 'vv2rr'
%    output data size  : [Nrois x Nrois x Ncontrastsubjects x Ncontrastconditions] 
%    output datafile   : a single 3D .mtx.nii NIFTI file for option 'vv2rr' (see help conn_mtx_read) 
% for option 'vv2rv'
%    output data       : [Ndim1 x Ndim2 x Ndim3 x Nrois x Ncontrastsubjects x Ncontrastconditions] 
%    output datafile   : a single 4D .vol.nii NIFTI file named [filename].nii when specifying a single subject (or subject-contrast) and condition (or condition-contrast) (see help conn_vol_read) 
%    output datafile   : a separate 4D .vol.nii NIFTI file named [filename]_subjectcontrast###_conditioncontrast###.nii for each subject (or subject-contrast) and condition (or condition-contrast) (see help conn_vol_read) 
%

% notes:
% data = conn_vv2rr(ROIfile, 'style','vv2rr','contrastsubjects',w1,'contrastconditions',w2)           : data(:,:,i,j) = sum_(subjects&conditions) { w1(i,subject)*w2(j,condition)*  Q' * R[subject,condition] * Q }
% [data, vol] = conn_vv2rr(ROIfile, 'style','vv2rv','contrastsubjects',w1,'contrastconditions',w2)    : data(:,:,:,nrois,i,j) = sum_(subjects&conditions) { w1(i,subject)*w2(j,condition)*  Q' * R[subject,condition] }
%    Q is voxel weights vector/matrix in 'file' (e.g. input ROI file)
%    R[subject,condition] is voxel-to-voxel correlation matrix computed during denoising
%
%   'validslices'         : list of valid slices (in vvPC* files) [default all] 

% note: call from within conn_process or conn_batch

global CONN_x;

options=struct(...
    'style','vv2rr',...             % vv2rr -> ROI-to-ROI matrix; vv2rv -> ROI-to-voxel matrix
    'isroi',true,...                % 1: treats input as ROI file (interprets 1:N 3d-file as 4d file, 4d file ROI weights=abs(values), and scales sum of weights across voxels within each ROI to 1)
    'filepath',[],...               % [location of voxel-to-voxel data]
    'saveas',[],...                 % saves output as NIFTI file
    'hmap',[],...                   % [only for use in context of MVPA analyses, source folder of 2nd-level analysis]
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
if isempty(options.contrastsubjects), options.contrastsubjects=eye(numel(options.validsubjects)); end
if isempty(options.validconditions), options.validconditions=1:nconditions; end
if isempty(options.contrastconditions), options.contrastconditions=eye(numel(options.validconditions)); end
% if ~isempty(options.hmap)
%     load(fullfile(options.RESULTSfolder,'SPM.mat'),'SPM');
%     options.validsubjects=find(SPM.xX.SelectedSubjects~=0);
%     [nill,tnames]=fileparts(SPM.xX_multivariate.Zfiles);
%     regexprep(tnames,'.*BETA_Subject(\d+)_Condition(\d+)_Measure(\d+)_Component(\d+).*','$1,$2,$3,$4');
%     options.validconditions=[]; %BETA_Subject001_Condition001_Measure001_Component008.nii'
% end

icondition=[];isnewcondition=[];for ncondition=1:nconditions,[icondition(ncondition),isnewcondition(ncondition)]=conn_conditionnames(CONN_x.Setup.conditions.names{ncondition}); end
if any(isnewcondition(options.validconditions)), error(['Some conditions have not been processed yet. Re-run previous step']); end
% addnew=false;
ROInames={}; ROIxyz={}; ROIthrtype='';
if iscell(ROI)||ischar(ROI), 
    if iscell(ROI), [ROI,ROIthrtype,ROIthr]=deal(ROI{:}); end
    if ~isempty(options.saveas), try, [ROIxyz,ROInames]=conn_roicenters(ROI); end; end
    ROI=conn_fileutils('spm_vol',ROI); 
end
assert(size(options.contrastsubjects,2)==numel(options.validsubjects),'mismatched dimensions between options.contrastsubjects and options.validsubjects');
assert(size(options.contrastconditions,2)==numel(options.validconditions),'mismatched dimensions between options.contrastsubjects and options.validconditions');

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
                switch(lower(ROIthrtype))
                    case '>', W=double(W>ROIthr);
                    case '>=', W=double(W>=ROIthr);
                    case {'=','=='}, W=double(W==ROIthr);
                    case '<=', W=double(W<=ROIthr);
                    case '<', W=double(W<ROIthr);
                end
                maxW=max(W(:));
                if options.isroi && maxW>1 && size(W,1)==1 && all(ismember(unique(W),0:maxW)), W=double(repmat(W,[maxW,1])==repmat((1:maxW)',[1,size(W,2)])); end
                if options.isroi, W=abs(W); end
                w=W(:,Y1.voxels);
                if options.isroi, 
                    temp=w*xyz(1:3,Y1.voxels)';
                    sw=sum(abs(w),2);
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
        if options.isroi&&~isempty(W), sw=sum(abs(w),2); end
        if ~isempty(options.contrastsubjects), C=options.contrastsubjects(:,ivalidsubject); 
        else C=1;
        end
        if ~isempty(options.contrastconditions), D=options.contrastconditions(:,ivalidcondition); 
        else D=1;
        end
        % w  : [Nrois x Nvox] ROI definition weights 
        % sw : [Nrois x 1] ROI sizes
        % x  : [Ncomp x Nvox] first-level dimensionality reduction components (x'*x = R) 
        % Z  : (vv2rr) [Nrois x Nrois x Ncontrastsubjects x Ncontrastconditions]
        %      (vv2rv) [Nrois x Nvox_dim1 x Nvox_dim2 x Nvox_dim3 x Ncontrastsubjects x Ncontrastconditions]
        %      (vv2vnrom) [Nvox_dim1 x Nvox_dim2 x Nvox_dim3 x Ncontrastsubjects x Ncontrastconditions] 
        switch(options.style)
            case 'vv2rr'
                if isempty(Z), Z=zeros(size(W,1),size(W,1),size(C,1),size(D,1)); end
                temp=x*w';
                z=(temp'*temp);
                if options.isroi, Zt=z./max(eps,sw*sw'); % w*x'*x*w' / w*1*1'*w' (x: [Ncomp x Nvox], x'*x = R)
                else Zt=z; % w*x'*x*w' 
                end
                for n1=1:size(C,1), for n2=1:size(D,1), Z(:,:,n1,n2)=Z(:,:,n1,n2)+Zt*C(n1)*D(n2); end; end
            case 'vv2rv'
                if isempty(Z), Z=zeros([size(W,1),Y1.matdim.dim,size(C,1),size(D,1)]); end
                temp=x*w';
                z=temp'*x;
                Zt=zeros([size(W,1),Y1.matdim.dim]); 
                if options.isroi, Zt(:,Y1.voxels)=diag(1./sw)*z;  % w*x'*x / w*1*1'
                else Zt(:,Y1.voxels)=z;  % w*x'*x
                end
                for n1=1:size(C,1), for n2=1:size(D,1), Z(:,:,:,:,n1,n2)=Z(:,:,:,:,n1,n2)+Zt*C(n1)*D(n2); end; end
            case 'vv2vnorm'
                if isempty(Z), Z=zeros([Y1.matdim.dim, size(C,1), size(D,1)]); end
                nx=sum(x.^2,2);
                z=nx'*(x.^2);
                Zt=zeros(Y1.matdim.dim); 
                Zt(Y1.voxels)=z'; % sum-square of rows of x'*x
                for n1=1:size(C,1), for n2=1:size(D,1), Z(:,:,:,n1,n2)=Z(:,:,:,n1,n2)+Zt*C(n1)*D(n2); end; end
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
end
if strcmp(options.style,'vv2rv'), Z = permute(Z, [2,3,4,1,5,6]); end  % dim1 dim2 dim3 rois contrastsubjects contrastconditions 
if ~isempty(options.saveas)
    switch(options.style)
        case 'vv2rr'
            labels={}; for n1=1:size(C,1), for n2=1:size(D,1), labels{n1,n2}=sprintf('subject-contrast #%d condition-contrast #%d',n1,n2); end; end
            conn_mtx_write(options.saveas, Z(:,:,:), ROInames, ROIxyz, labels(:));
            if ~nargout, fprintf('output saved as %s\n',options.saveas); end
        case 'vv2rv'
            if size(C,1)==1&&size(D,1)==1
                conn_vol_write(options.saveas, Z, zmat);
                if ~nargout, fprintf('output saved as %s\n',options.saveas); end
            else
                for n1=1:size(C,1),
                    for n2=1:size(D,1),
                        filename=conn_prepend('',options.saveas, sprintf('_subjectcontrast%03d_conditioncontrast%03d.nii',n1,n2));
                        conn_vol_write(filename, Z(:,:,:,:,n1,n2), zmat);
                    end
                end
                if ~nargout, fprintf('output saved as %s\n',conn_prepend('',options.saveas,'_subjectcontrast*_conditioncontrast*.nii')); end
            end
        case 'vv2vnorm'
            if size(C,1)==1&&size(D,1)==1
                conn_vol_write(options.saveas, Z, zmat);
                if ~nargout, fprintf('output saved as %s\n',options.saveas); end
            else
                for n1=1:size(C,1),
                    for n2=1:size(D,1),
                        filename=conn_prepend('',options.saveas, sprintf('_subjectcontrast%03d_conditioncontrast%03d.nii',n1,n2));
                        conn_vol_write(filename, Z(:,:,:,n1,n2), zmat);
                    end
                end
                if ~nargout, fprintf('output saved as %s\n',conn_prepend('',options.saveas,'_subjectcontrast*_conditioncontrast*.nii')); end
            end
    end
end
    
%     if ~isempty(folderout)
%         filename=fullfile(folderout,['resultsROI_Condition',num2str(icondition(ncondition),'%03d'),'.mat']);
%         save(filename,'Z','names','names2','xyz'); %%%
%     end
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