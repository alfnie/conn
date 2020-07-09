function x=conn_surf_morph(x,faces,op,N)
% conn_surf_morph morphological operations on the surface
%
% y = conn_surf_morph(x,faces,'dilate' [,N]) performs N steps of binary mask dilation on the surface
% y = conn_surf_morph(x,faces,'erode' [,N]) performs N steps of binary mask erosion on the surface
%     x      : [Mx1] input surface masks
%     faces  : [Mx3] triangular faces defining surface tessellation (if empty it assumes fsaverage conn_surf_sphere(8) tessellation)
%     y      : [Mx1] output surface mask
%
%

if isempty(faces), faces=conn_surf_sphere(8); end
if isstruct(faces), ref=faces;
else ref=struct('faces',faces);
end
if nargin<3||isempty(op), op='dilate'; end
if nargin<4||isempty(N), N=1; end

L=max(ref.faces(:));
A=sparse([ref.faces(:,1);ref.faces(:,2);ref.faces(:,3)],[ref.faces(:,2);ref.faces(:,3);ref.faces(:,1)],1,L,L);
A=A|A'|speye(L);

sX=size(x);
assert(~rem(numel(x),size(A,1)),'incorrect dimensions of input x (%d elements, expected a multiple of %d elements)',numel(x),size(A,1));
x=reshape(x,size(A,1),[]);

switch(lower(op))
    case {'dilate','dilation'}
        for n=1:N
            x=(A*double(x>0))>0;
        end
    case {'erode','erosion'}
        for n=1:N
            x=~(A*double(~x));
        end
    case 'dilaterois'
        nROIs=max(x(:));
        nvertices=accumarray(x(x>0),1);
        [nill,idx]=sort(nvertices,'descend'); %idx=1:nROIs;
        while N>0
            for n=1:nROIs
                mask1=x==idx(n);
                mask2=(A*double(mask1))>0;
                x(~x&mask2)=idx(n);
            end
            if ~nnz(~x), return; end
            N=N-1;
            fprintf('.');
        end
    case 'eroderois'
        nROIs=max(x(:));
        for n=1:nROIs
            mask1=x==n;
            mask2=mask1;
            for m=1:N
                mask2=~(A*double(~mask2));
            end
            x(mask1&~mask2)=0;
        end
    otherwise, error('unknown operation %s',op);
end
x=reshape(x,sX);