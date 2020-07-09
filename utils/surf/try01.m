ref=conn_surf_sphere(8);
N=size(ref.vertices,1);
Nvox=2;
A=double(spm_mesh_adjacency(ref.faces)|speye(N));
coords=ref.vertices';
dX=sparse(N,N);
dY=sparse(N,N);
I=[];
J=[];
V=[];
dX=sparse(N,N);
dY=sparse(N,N);
for i=1:N,
    if ~rem(i,1e3), disp(i); end
    j=sparse(i,1,1,N,1);
    for n=2:Nvox, j=A*j; end
    ij=find(j); %B*j);
    
    % base tangential to spherical surface
    xyz=coords(:,i);
    [nill,pivot]=max(abs(xyz));
    other=1:3;
    other(pivot)=[];
    base=zeros(3,2);
    base(pivot,:)=-xyz(other);
    base(other(1),1)=xyz(pivot);
    base(other(2),2)=xyz(pivot);
    base(:,1)=base(:,1)/sqrt(sum(abs(base(:,1)).^2));
    base(:,2)=base(:,2)-base(:,1)*(base(:,2)'*base(:,1));
    base(:,2)=base(:,2)/sqrt(sum(abs(base(:,2)).^2));
    xyz2=base'*coords(:,ij);
    temp=xyz2(1,:)>0; xyz2(1,temp)=xyz2(1,temp)/max(eps,sum(abs(xyz2(1,temp))));
    temp=xyz2(1,:)<0; xyz2(1,temp)=xyz2(1,temp)/max(eps,sum(abs(xyz2(1,temp))));
    temp=xyz2(2,:)>0; xyz2(2,temp)=xyz2(2,temp)/max(eps,sum(abs(xyz2(2,temp))));
    temp=xyz2(2,:)<0; xyz2(2,temp)=xyz2(2,temp)/max(eps,sum(abs(xyz2(2,temp))));
    I=[I,i+zeros(1,numel(ij))];
    J=[J,ij'];
    V=[V,xyz2];
    if ~rem(i,1e3)||i==N
        dX=dX+sparse(I,J,V(1,:),N,N);
        dY=dY+sparse(I,J,V(2,:),N,N);
        I=[];
        J=[];
        V=[];
    end
    %dX(i,ij)=xyz2(1,:);
    %dY(i,ij)=xyz2(2,:);
end
%save('surf_top.mat','A','dX','dY');
