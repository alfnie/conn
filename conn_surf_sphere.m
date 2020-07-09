function [ref,sphere2xyz,xyz2sphere,Msphere2xyz,Mxyz2sphere]=conn_surf_sphere(resolution,xyz,issphere)
% ref=SURF_SPHERE(res) generates patch data (vertices/faces)
% for a unit sphere using a recursively-subdivided icosahedral grid.
%
% [ref,ref2xyz,xyz2ref]=SURF_SPHERE(res,xyz) returns additional
% coregistration indexes for secondary coordinates xyz 
% (Nx3 matrix of points in a unit sphere) such that:
%   xyz(xyz2ref,:) ~= ref.vertices
%   ref.vertices(ref2xyz,:) ~= xyz
%
% res vertices faces
% 1      12      20 (icosahedron)
% 2      42      80
% 3     162     320
% 4     642    1280
% 5    2562    5120
% 6   10242   20480
% 7   40962   81920
% 8  163842  327680
% n 2+5/2*4^n  5*4^n

if nargin<1||isempty(resolution), resolution=8; end
if nargin<2, xyz=[]; end
if nargin<3||isempty(issphere), issphere=true; end

if isstruct(xyz), xyz_faces=xyz.faces; xyz=xyz.vertices; 
else xyz_faces=[];
end

errorcheckxyz=false;
ang1=(0:2*pi/5:2*pi-pi/5)';
ang2=atan(1/2);
ref.vertices=[0,0,1; [cos(ang2)*cos(ang1), cos(ang2)*sin(ang1), sin(ang2)*ones(size(ang1))]; [cos(-ang2)*cos(ang1+pi/5), cos(-ang2)*sin(ang1+pi/5), sin(-ang2)*ones(size(ang1))]; 0,0,-1];
ref.faces=[1 2 3; 1 3 4; 1 4 5; 1 5 6; 1 6 2; 2 7 3; 3 8 4; 4 9 5; 5 10 6; 6 11 2; 3 7 8; 4 8 9; 5 9 10; 6 10 11; 2 11 7; 7 12 8; 8 12 9; 9 12 10; 10 12 11; 11 12 7];
ref.cdata=(1:20)';
if ~isempty(xyz)
    xyz=conn_bsxfun(@rdivide,xyz,sqrt(sum(abs(xyz).^2,2)));
    match=xyz*ref.vertices';
    [maxsphere2xyz,sphere2xyz]=max(match,[],2); % for each point in xyz, index to closest point in ref
    [maxxyz2sphere,xyz2sphere]=max(match',[],2);% for each point in ref, index to closest point in xyz
end

% resolution>1: sequential partition triangular faces (1->4)
for res=2:resolution
    n=size(ref.vertices,1);
    m=size(ref.faces,1);
    oldfaces=ref.faces;
    % for each edge, a new vertex
    edges=sort([ref.faces(:,1:2);ref.faces(:,[2,3]);ref.faces(:,[1,3])],2);
    [uedges,sedges,idx]=unique(edges,'rows');
    iedges=reshape(idx,m,3);
    newvertices=(ref.vertices(uedges(:,1),:)+ref.vertices(uedges(:,2),:))/2;
    if issphere, newvertices=newvertices./repmat(sqrt(sum(newvertices.^2,2)),1,3); end
    % divide each face in 4 faces
    newfaces=[ref.faces(:,1),n+iedges(:,1),n+iedges(:,3); ref.faces(:,2),n+iedges(:,2),n+iedges(:,1); ref.faces(:,3),n+iedges(:,3),n+iedges(:,2); n+iedges];
    
    ref.vertices=cat(1,ref.vertices,newvertices);
    ref.faces=newfaces;
    ref.cdata=repmat(ref.cdata,[4,1]);
    
    if ~errorcheckxyz&&~isempty(xyz), % registers secondary surface
        
        [m,i]=max([sum(newvertices.*xyz(xyz2sphere(uedges(:,1)),:),2),sum(newvertices.*xyz(xyz2sphere(uedges(:,2)),:),2)],[],2);
        xyz2sphere=[xyz2sphere;xyz2sphere(uedges((1:size(uedges,1))'+size(uedges,1)*(i-1)))];
        maxxyz2sphere=[maxxyz2sphere;m];

        % for each new-vertex, consider the four old-vertices surrounding it
        mask=sparse(oldfaces(:,1),oldfaces(:,2),1,n,n)|sparse(oldfaces(:,1),oldfaces(:,3),1,n,n)|sparse(oldfaces(:,2),oldfaces(:,3),1,n,n);
        mask=mask|mask';
        [i,j]=find(mask(:,uedges(:,1))&mask(:,uedges(:,2))); 
        assert(isequal(j,floor(1:.5:size(uedges,1)+.5)'), 'unexpected condition'); % [i,j]=find(...) finds always two non-zero elements per column
        toedges=reshape(i,2,[])';
        xuedges=[uedges,toedges];
        xuedges=uedges;
        
        % for each xyz subset closest to one given old vertex, 
        % reconsider the closest match with four of the new vertices.
        
        [oldVert_sorted,oldVert_idx]=sort(sphere2xyz);           % indices to old vertices (size xyz)
        oldVert_segments=find([true;diff(oldVert_sorted);true]); % indices to xyz (size segments1)
        [newVert_sorted,newVert_idx]=sort(xuedges(:));            % indices to old vertices (size 4*newvertices)
        newVert_idx=1+mod(newVert_idx-1,size(uedges,1));         % indices to newvertices (size 4*newvertices)
        newVert_segments=find([true;diff(newVert_sorted);true]); % indices to newvertices (size segments2)
        nsegments2=1;
        
        for nsegments1=1:numel(oldVert_segments)-1,
            %oldVert_set=find(sphere2xyz==oldvert);
            %newVert_set=find(any(uedges==oldvert,2));
            oldvert=oldVert_sorted(oldVert_segments(nsegments1));
            oldVert_set=oldVert_idx(oldVert_segments(nsegments1):oldVert_segments(nsegments1+1)-1);
            while newVert_sorted(newVert_segments(nsegments2))<oldvert, nsegments2=nsegments2+1; end
            if newVert_sorted(newVert_segments(nsegments2))~=oldvert, error('unknown condition'); end
            newVert_set=newVert_idx(newVert_segments(nsegments2):newVert_segments(nsegments2+1)-1);
            match=xyz(oldVert_set,:)*newvertices(newVert_set,:)';
            [mmatch,imatch]=max(match,[],2);
            newmax=mmatch>maxsphere2xyz(oldVert_set);
            sphere2xyz(oldVert_set(newmax))=n+newVert_set(imatch(newmax));
            maxsphere2xyz(oldVert_set(newmax))=mmatch(newmax);
            [mmatch,imatch]=max(match',[],2);
            newmax=mmatch>maxxyz2sphere(n+newVert_set);
            xyz2sphere(n+newVert_set(newmax))=oldVert_set(imatch(newmax));
            maxxyz2sphere(n+newVert_set(newmax))=mmatch(newmax);
        end
    end
end

if errorcheckxyz&&~isempty(xyz)
    sphere2xyz=zeros(size(xyz,1),1);
    for n1=1:size(xyz,1)
        match=ref.vertices*xyz(n1,:)';
        [maxsphere2xyz,sphere2xyz(n1)]=max(match); % for each point in xyz, index to closest point in ref
    end
    xyz2sphere=zeros(size(ref.vertices,1),1);
    for n1=1:size(ref.vertices,1)
        match=xyz*ref.vertices(n1,:)';
        [maxxyz2sphere,xyz2sphere(n1)]=max(match); % for each point in ref, index to closest point in xyz
    end
end

if resolution>8, conn_disp('warning: Non-existing FreeSurfer ico reference. Files will not be compatible with freesurfer ico space');
else % makes order of vertices compatible with freesurfer ico-reference files
    load conn_surf_sphere.mat
    ind2freesurfer=IND2FREESURFER{resolution};
    freesurfer2ind=FREESURFER2IND{resolution};
    ref.vertices=ref.vertices(ind2freesurfer,:);
    ref.faces=fliplr(freesurfer2ind(ref.faces));
    if ~isempty(xyz)
        sphere2xyz=freesurfer2ind(sphere2xyz);
        xyz2sphere=xyz2sphere(ind2freesurfer);
    end
end

if nargout>3&&~isempty(xyz)
    Msphere2xyz = conn_surf_interpmtx(ref,xyz,sphere2xyz);
end
if nargout>4&&~isempty(xyz)
    if ~isempty(xyz_faces), Mxyz2sphere = conn_surf_interpmtx(struct('vertices',xyz,'faces',xyz_faces),ref.vertices,xyz2sphere);
    else conn_disp('warning: this functionality needs additional xyz surface information; please enter a surface structure in xyz argument'); Mxyz2sphere=[]; 
    end
end

if ~nargout
    h=patch(ref);set(h,'edgecolor','k','facecolor',.75*[1 1 1]);axis equal; grid on
end
end