function [S,P,Q,I]=conn_tfce(h,varargin)
% CONN_TFCE threshold free cluster enhancement transformation (Smith & Nichols 2009)
%
% S = conn_tfce(h)
% where h is a n-d matrix characterizing an intensity function h(x) 
% computes S(x) = Integral between Hmin and h(x) of { extent(x,h)^E * h^H dh }
% where extent(x,h) is the size of a contiguous cluster of pixels/voxels including point x and with intensities above h
% Model default values: E=0.5, H=2, Hmin=0
%
% Alternative model parameters:
% S = conn_tfce(h [,param_name1, param_value1, param_name1, param_value2, ...]) uses alternative model parameters
%    'E'      : extent function exponent (default 0.5)
%    'H'      : intensity function exponent (default 2)
%    'Hmin'   : integration lower limit (default 0)
%
% Additional information:
% - integral estimated exactly, no dh discretization/approximation needed (equivalent to lim dh->0 results in other implementations)
% - connectivity criterion = 2 for vectors, 8 for images, 18 for 3-dimensional h matrices, and 3^D-2^D-1 for higher dimensions (note: for compatibility with spm_clusters)
% - properties:
%     conn_tfce(h,'E',0,'H',0) = h
%     conn_tfce(h,'E',E,'H',0) = conn_tfce(k*h,'E',E,'H',0) / k
%     conn_tfce(h,'E',0,'H',H) = h^(H+1)/(H+1)
%     conn_tfce(h,'E',E,'H',H) = conn_tfce(h^(H+1)/(H+1),'E',E,'H',0)
% - [S,P,Q] = conn_tfce(...) also returns 
%           P: local-peak voxels (local maxima/minima in output S(x) -and input h(x)-)
%           Q: maximum S(x) within same support region as x (dx | extent(x+dx,h) is a continuous function of h)
%

% alfnie@gmail.com 2019

params=struct(...
    'e',.5,...                    % E constant
    'h',2,...                     % H constant
    'hmin',0,...                  % Hmin costant
    'adjacency',[],...            % connectivity criterion: a) adjacency=[] or adjacency='volume' uses default connectivity criterion; b) adjacency='full' for 3^D-1 connectivity criterion in all cases; c) adjacency=sparse(...,N,N) uses explicit adjacency matrix; or d) adjacency=struct('size',...,'index',...) uses default connectivity criterion with voxels selected from n-dimensional matrix with size "size" and single-index positions "index"
    'removeborders',false);       % true/false: P output (identifying local peaks) does not include border voxels (voxels neighb to x|isnan(h(y)))
for n1=1:2:numel(varargin)-1, if ~isfield(params,lower(varargin{n1})), error('unknown option %s',lower(varargin{n1})); else params.(lower(varargin{n1}))=varargin{n1+1}; end; end

% resize
adj=params.adjacency;
if ischar(adj), adj=lower(adj); assert(ismember(adj,{'full','volume'}),'unrecognized adjacency value %s',adj); end
if isequal(adj,'volume'), adj=[]; end
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
        h=cat(1,nan([1,sh(2:end)]),h,nan([1,sh(2:end)]));
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
mask=reshape(find(~isnan(h)&h>params.hmin),[],1);
[nill,idx]=sort(h(mask),'descend');
idx=mask(idx);

%nh=2000;
%peaks=max(h(conn_bsxfun(@plus,mask,neighb')),[],2)<=h(mask);
Q=zeros(sh);               % local support
I=zeros(sh);               % index to all clusters
S=zeros(sh);               % TFCE integral between h(x) and h_max within each voxel
P=false(sh);               % location of peaks
C_open=zeros(nh,1);        % index to open clusters
C_extent=zeros(nh,1);      % number of voxels within each cluster at the threshold level of the last voxel sampled
C_inth=zeros(nh,1);        % indefinite integral of function h^H within each cluster at the last voxel sampled
C_score=zeros(nh,1);       % TFCE integral between h_min and h_max within each closed cluster
C_bak=zeros(nh,1);         % list of clusters merged (for baktrack)
lastI=0;

%info=zeros(sh);
for j=idx(:)'
    if isempty(adj), t=I(j+neighb);
    else t=I(adj(:,j)>0);
    end
    i=C_open(t(t>0));
%info(j)=numel(i);        
    if isempty(i)                       % this voxel is a local peak => start a new cluster
        lastI=lastI+1;
        S(j)=0;
        I(j)=lastI;
        P(j)=true; % local maxima
        C_score(lastI)=0;
        C_extent(lastI)=1;
        C_inth(lastI)=h(j).^(params.h+1)/(params.h+1);
        C_open(lastI)=lastI;
    elseif all(i==i(1))                 % this voxel is adjacent to an existing cluster => grow this cluster
        i=i(1);
        inth=h(j).^(params.h+1)/(params.h+1);
        ds=(C_extent(i).^params.e).*(C_inth(i)-inth);
        C_score(i)=C_score(i)+ds;
        S(j)=C_score(i);
        I(j)=i;
        C_extent(i)=C_extent(i)+1;
        C_inth(i)=inth;
    else                                % this voxel is adjacent to several existing clusters => merge these clusters
        i=find(sparse(i,1,1)); % i=unique(i)
        inth=h(j).^(params.h+1)/(params.h+1);
        ds=(C_extent(i).^params.e).*(C_inth(i)-inth);
        C_score(i)=C_score(i)+ds;
        lastI=lastI+1;
        S(j)=0;
        I(j)=lastI;
        %P(j)=true; % local minima
        C_score(lastI)=0;
        C_extent(lastI)=sum(C_extent(i))+1;
        C_inth(lastI)=inth;
        C_open(lastI)=lastI;        
        C_open(ismember(C_open(1:lastI),i))=lastI;
        C_bak(i)=lastI;
    end
end


% backpropagate
for i=lastI:-1:1
    if i<=nh&&C_bak(i), C_score(i)=C_score(i)+C_score(C_bak(i)); 
    else C_score(i)=C_score(i)+(C_extent(i).^params.e).*(C_inth(i)-params.hmin.^(params.h+1)/(params.h+1));
    end
end
Q(mask) = C_score(I(mask));
S(mask) = Q(mask)-S(mask);

% remove border peaks
if params.removeborders&&nargout>=2
    idx=find(P);
    keep=true(numel(idx),1);
    if isempty(adj)
        for n=1:numel(neighb)
            %keep=keep&h(idx+neighb(n))>params.hmin;
            keep=keep&~isnan(h(idx+neighb(n)));
        end
    else
        for n=1:numel(idx)
            %keep(n)=all(h(adj(:,idx(n))>0)>params.hmin);
            keep(n)=all(~isnan(h(adj(:,idx(n))>0)));
        end
    end
    P(idx(~keep))=false;
end
% remove non-peak clusters
%[nill,t]=ismember(1:lastI,I(P));
%I(I>0)=t(I(I>0));
% remove zero-padd
if isempty(adj)
    idx=cell(1,ndims(S)); for n=1:ndims(S),idx{n}=2:size(S,n)-1; end; 
    S=S(idx{:});
    if nargout>=2, P=P(idx{:}); end
    if nargout>=3, Q=Q(idx{:}); end
    if nargout>=4, I=I(idx{:}); end
end
% resize
if isstruct(params.adjacency)
    S=reshape(S(params.adjacency.index),size(h0));
    if nargout>=2, P=reshape(P(params.adjacency.index),size(h0)); end
    if nargout>=3, Q=reshape(Q(params.adjacency.index),size(h0)); end
    if nargout>=4, I=reshape(I(params.adjacency.index),size(h0)); end
end

end
