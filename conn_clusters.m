function [nlabels,labels]=conn_clusters(mask,A,aexp)
% labels connected components (bwlabel for volume/surface/matrix/arbitrary connectivity)
%

persistent surfparams;
if nargin<3, aexp=[]; end
if nargin<2, datatype='matrix';
elseif isempty(A), datatype='get_matrix';
elseif ischar(A), datatype=A; A=[];
else datatype='explicit';
end
   
if isempty(aexp), aexp=str2double(regexp(datatype,'\d+$','match','once')); end
datatype=regexprep(datatype,'\d+$','');

switch(datatype)
    case 'network' % uses NBS node-connectivity for ROI-to-ROI analyses
        [Nr1,Nr2]=size(mask);
        [i,j]=find(mask);
        M=sparse(i(i~=j),j(i~=j),1,max(Nr1,Nr2),max(Nr1,Nr2)); % handles non-square matrices
        M=M|M';
        valid=find(any(M,1));
        M=M(valid,valid)|speye(numel(valid));
        [p,q,r,s]=dmperm(M);
        n=diff(r);
        i=find(n>1);
        Nlabels=zeros(size(M,1),1); % node labels
        %Nnlabels=zeros(numel(i),1); % node count
        for i1=1:numel(i)
            Nlabels(p(r(i(i1)):r(i(i1)+1)-1))=i1;
            %Nnlabels(i1)=n(i(i1));
        end
        labels=zeros(Nr1,Nr2);
        rmask=mask(valid(valid<=Nr1),valid);
        rmask(1:size(rmask,1)+1:size(rmask,1)*size(rmask,1))=0;
        labels(valid(valid<=Nr1),valid)=conn_bsxfun(@times,Nlabels(valid<=Nr1),rmask); % edge labels
        nlabels=accumarray(Nlabels(valid<=Nr1),sum(rmask,2)); % edge count
    
    case 'matrix' % uses sorted ROIs one-dimensional lattice connectivity for ROI-to-ROI analyses
        [Nr1,Nr2]=size(mask);
        [i,j]=find(mask);
        M=sparse(i(i~=j),j(i~=j),1,max(Nr1,Nr2),max(Nr1,Nr2)); % handles non-square matrices
        M=M|M'; % forces symmetric neighborhoods
        if ~isempty(aexp)&&aexp>1, 
            [i1,i2]=ndgrid(-round(aexp)+1:round(aexp)-1,-round(aexp)+1:round(aexp)-1);
            M2=convn(full(double(M)),double(abs(i1)+abs(i2)<round(aexp)),'same')>0; 
        else M2=M;
        end
        M2(1:size(M,1)+1:size(M,1)^2)=0;
        [L,num]=spm_bwlabel(triu(double(full(M2))));
        labels=(L+L').*full(M);
        labels=labels(1:Nr1,1:Nr2);
        nlabels=accumarray(labels(labels>0),1);
    
    case 'get_matrix',  % returns default node-connectivity structure for ROI-to-ROI analyses
        [Nr1,Nr2]=size(mask);
        A=cat(1,repmat(sparse(1:Nr1*Nr1,repmat(1:Nr1,Nr1,1),1),1,Nr2),sparse(Nr1*(Nr2-Nr1),Nr1*Nr2));
        A=A|A'|sparse(1:size(A,1),1:size(A,1),1);
        nlabels=double(A);

    case 'none'
        idx=find(mask);
        labels=zeros(size(mask));
        labels(idx)=1:numel(idx);
        nlabels=ones(numel(idx),1);
        
    otherwise, % for surface, volume, or other arbitrary node-connectivity definitions
        if isempty(A) % if no info assume volume or surface
            if strcmp(datatype,'volume')
                [tx,ty,tz]=ndgrid(1:size(mask,1),1:size(mask,2),1:size(mask,3));
                %txyz=[tx(:) ty(:) tz(:)];
                %A=struct('xyz',txyz');
                A=struct('xyz',[tx(:) ty(:) tz(:)]');
            else
                if isempty(surfparams)
                    surfparams=load(fullfile(fileparts(which(mfilename)),'utils','surf','surf_top.mat'),'A');
                end
                A=surfparams.A;
            end
        end
        if (isstruct(A)&&isfield(A,'xyz'))||(~issparse(A)&&size(A,1)==3), % volume-case (entering xyz voxel coordinates)
            labels=zeros(size(mask));
            if isstruct(A), l=spm_clusters(A.xyz(:,mask));
            else l=spm_clusters(A(:,mask));
            end
            nlabels=accumarray(l(:),1);
            labels(mask)=l;
            return;
        end
        if (isstruct(A)&&isfield(A,'faces'))||(~issparse(A)&&size(A,2)==3), A=spm_mesh_adjacency(A); end % surface-case (entering triangular faces info)
        if ~isempty(aexp)&&aexp>1, A=A^aexp; end
        sizemask=size(mask);
        if size(mask,1)~=size(A,2), mask=reshape(mask,size(A,2),[]); end
        [n,n2]=size(mask);
        mask=mask>0;
        labels=zeros(size(mask));
        nlabels=zeros(n*n2,1);
        for nlabel=1:n*n2
            if ~nnz(mask), break; end
            [i,j]=find(mask,1);
            
            d=sparse(i,j,1,n,n2);
            e=sparse(i,j,1,n,n2);
            c=mask;c(i,j)=false;
            while 1
                d=c&(A*d);
                if ~nnz(d), break; end
                c(d)=false;
                e=e|d;
            end
            j=find(e);
            labels(j)=nlabel;
            mask(j)=false;
            nlabels(nlabel)=numel(j);
        end
        nlabels=nlabels(1:nlabel-1);
        if ~isequal(size(labels),sizemask), labels=reshape(labels,sizemask); end
end
end

