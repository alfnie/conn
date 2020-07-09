function Mref2xyz = conn_surf_interpmtx(ref,xyz,ref2xyz)
% conn_surf_interpmtx trilinear interpolation matrix
% 
% Mref2xyz = conn_surf_interpmtx(ref,xyz [,ref2xyz])
%   ref : reference surface: structure with fields "vertices" [Mx3] and "faces" [Px3] (e.g. output of conn_surf_sphere or conn_surf_readsurf)
%   xyz : [Nx3] matrix of coordinates of points within reference surface
%   ref2xyz : [Nx1] closest vertex in ref for each point in xyz (e.g. output of conn_surf_sphere)
%   Mref2xyz : [NxM] Trilinear interpolation sparse matrix such that Mref2xyz * val(ref) ~= val(xyz)
%

Norder=1; % distance of neighbor set
if nargin<3||isempty(ref2xyz),
    ref2xyz=zeros(size(xyz,1),1);
    %nxyz=sum(xyz.^2,2);
    nref=sum(ref.vertices.^2,2);
    for n1=1:size(xyz,1)
        match=ref.vertices*xyz(n1,:)';
        match=2*match-nref;%-nxyz(n1);
        [maxref2xyz,ref2xyz(n1)]=max(match); % for each point in xyz, index to closest point in ref
    end
end
A=sparse([ref.faces(:,1);ref.faces(:,2);ref.faces(:,3)],[ref.faces(:,2);ref.faces(:,3);ref.faces(:,1)],1,size(ref.vertices,1),size(ref.vertices,1));
A=A|A'|speye(size(ref.vertices,1));
rA=sparse(ref2xyz,1:numel(ref2xyz),1,size(ref.vertices,1),numel(ref2xyz));
for n=1:Norder, rA=A*rA; end
%rA=A(:,ref2xyz);
nrA=sum(rA,1);
unrA=accumarray(full(nrA'),1);
[i,j]=find(rA); % for each point in xyz, find closest point + its neighbors (5-6)
d=sum(abs(ref.vertices(i,:)-xyz(j,:)).^2,2); % distances
dj=[0 find(diff(j))' numel(j)];
ti=zeros(3*size(xyz,1),1);tj=ti;tw=ti;
for n=1:size(xyz,1)
    idx=dj(n)+1:dj(n+1);
    [sd,k]=sort(d(idx));
    ok=false; minttw_best=-inf; tti_best=ones(3,1); ttw_best=nan(3,1);
    for l=nchoosek(1:numel(idx),3)' % subsets of three ref points (to find face containing xyz; if it exists)
        tti=i(idx(k(l')));
        ttw=([xyz(n,:),1]/[ref.vertices(tti,:),ones(3,1)])'; % trilinear interp
        minttw=min(ttw);
        if minttw>=0; minttw_best=minttw; tti_best=tti; ttw_best=ttw; break;
        elseif minttw>minttw_best, minttw_best=minttw; tti_best=tti; ttw_best=ttw;
        end
    end
    ti(3*n-2:3*n)=tti_best;
    tj(3*n-2:3*n)=n+zeros(3,1);
    tw(3*n-2:3*n)=ttw_best;
end
Mref2xyz=sparse(tj,ti,tw,size(xyz,1),size(ref.vertices,1));
end