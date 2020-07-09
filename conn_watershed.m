function [I,P]=conn_watershed(h,varargin)
% CONN_WATERSHED watershed segmentation
%
% S = conn_watershed(h)

% alfnie@gmail.com 2019

params=struct(...
    'minh',-inf,...
    'adjacency',[],...            % connectivity criterion: a) adjacency=[] uses default connectivity criterion; b) adjacency='full' for 3^D-1 connectivity criterion in all cases; c) adjacency=sparse(...,N,N) uses explicit adjacency matrix; or d) adjacency=struct('size',...,'index',...) uses default connectivity criterion with voxels selected from n-dimensional matrix with size "size" and single-index positions "index"
    'removeborders',false);       % true/false: P output (identifying local peaks) does not include border voxels (voxels neighb to x|isnan(h(x)))
for n1=1:2:numel(varargin)-1, if ~isfield(params,lower(varargin{n1})), error('unknown option %s',lower(varargin{n1})); else params.(lower(varargin{n1}))=varargin{n1+1}; end; end

% resize
adj=params.adjacency;
if isstruct(adj)
    h0=h;
    h=zeros(adj.size(:)');
    h(adj.index)=h0;
    adj=[];
end
sh=size(h);
nh=numel(h);
% zero-padd
if isempty(adj)||isequal(adj,'full')
    for n=1:ndims(h)
        sh=size(h);
        h=cat(1,-inf([1,sh(2:end)]),h,-inf([1,sh(2:end)]));
        h=shiftdim(h,1);
    end
    sh=size(h);
    nh=numel(h);
    [dn{1:numel(sh)}]=ndgrid(1:3);
    neighb=reshape(sub2ind(sh,dn{:})-sum(cumprod(sh(1:end-1)))-2,[],1);
    if isempty(adj)&&numel(sh)>=3, neighb(all(abs(cat(4,dn{:})-2)==1,4))=[]; end % e.g. 8-connectivity in 2d, 18-connectivity in 3d
    neighb(~neighb)=[];
    adj=[];
end
mask=reshape(find(~isnan(h)&h>params.minh),[],1);
[nill,idx]=sort(h(mask),'descend');
idx=mask(idx);

I=zeros(sh);               % index to all clusters
P=false(sh);               % location of peaks
lastI=0;

%info=zeros(sh);
for j=idx(:)'
    if isempty(adj), t=I(j+neighb);
    else t=I(adj(:,j)>0);
    end
    i=t(t>0);
%info(j)=numel(i);        
    if isempty(i)                       % this voxel is a local peak => start a new cluster
        lastI=lastI+1;
        I(j)=lastI;
        P(j)=true; 
    elseif all(i==i(1))                 % this voxel is adjacent to an existing cluster => grow this cluster
        i=i(1);
        I(j)=i;
    else                                % this voxel is adjacent to several existing clusters => disregard
    end
end

% remove border peaks
if params.removeborders&&nargout>1
    idx=find(P);
    keep=true(numel(idx),1);
    if isempty(adj)
        for n=1:numel(neighb)
            keep=keep&h(idx+neighb(n))>params.minh;
        end
    else
        for n=1:numel(idx)
            keep(n)=all(h(adj(:,idx(n))>0)>params.minh);
        end
    end
    P(idx(~keep))=false;
end
% remove zero-padd
if isempty(adj)
    idx=cell(1,ndims(P)); for n=1:ndims(P),idx{n}=2:size(P,n)-1; end; 
    P=P(idx{:});
    I=I(idx{:});
end
% resize
if isstruct(params.adjacency)
    P=reshape(P(params.adjacency.index),size(h0));
    I=reshape(I(params.adjacency.index),size(h0));
end

end
