function Z=conn_vv2rr(ROI,filepath,folderout,validconditions)
% computes ROI-to-ROI matrix from Voxel-to-Voxel SVD representation for each subject
global CONN_x;

nconditions=length(CONN_x.Setup.conditions.names)-1;
if nargin<2||isempty(filepath), filepath=CONN_x.folders.preprocessing; end
if nargin<3, folderout=[]; end
if nargin<4||isempty(validconditions), validconditions=1:length(CONN_x.Setup.conditions.names)-1; end
icondition=[];isnewcondition=[];for ncondition=1:nconditions,[icondition(ncondition),isnewcondition(ncondition)]=conn_conditionnames(CONN_x.Setup.conditions.names{ncondition}); end
if any(isnewcondition(validconditions)), error(['Some conditions have not been processed yet. Re-run previous step']); end
% addnew=false;
if iscell(ROI), ROI=char(ROI); end
if ischar(ROI), ROI=spm_vol(ROI); end

%h=conn_waitbar(0,'Extracting correlation matrix, please wait...');
for ivalidcondition=1:numel(validconditions),
    ncondition=validconditions(ivalidcondition);
    Z=zeros(size(ROI,1),size(ROI,1),CONN_x.Setup.nsubjects);
    for nsub=1:CONN_x.Setup.nsubjects
        filename_B1=fullfile(filepath,['vvPC_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'.mat']);
        Y1=conn_vol(filename_B1);
        if isstruct(ROI), 
            xyz=conn_convertcoordinates('idx2tal',1:prod(Y1.matdim.dim),Y1.matdim.mat,Y1.matdim.dim)';
            W=spm_get_data(ROI,pinv(ROI.mat)*xyz);
            if size(W,1)==1&&isequal(reshape(unique(W),[],1),(0:max(ROI(:)))'), W=double(repmat(W,[max(W(:)),1])==repmat((1:max(W(:)))',[1,size(W,2)])); end
        else W=ROI;
        end
        assert(prod(Y1.matdim.dim)==size(W,2),'unequal number of voxels in ROI (%d) and VV (%d) files',size(W,2), prod(Y1.matdim.dim));
        [x,idx]=conn_get_volume(Y1);
        w=W(:,Y1.voxels);
        Z(:,:,nsub)=((w*x')*(x*w'))./max(eps,w*w');
        %temp=temp./repmat(max(eps,sum(abs(W),2)'),size(x,1),1);
        %Z(:,:,nsub)=temp'*temp;
        %conn_waitbar(((ivalidcondition-1)*CONN_x.Setup.nsubjects+nsub)/numel(validconditions)/CONN_x.Setup.nsubjects,h,sprintf('Subject %d Condition %d',nsub,ncondition));
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