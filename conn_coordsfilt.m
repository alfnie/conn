function [xyzout,idx]=conn_coordsfilt(xyz,DisMin)
% CONN_COORDSFILT filters sorted list of coordinates to remove points closer than DisMin from previous points
%
% xyz_output=conn_coordsfilt(xyz,DisMin)
%    xyz     : [3xN] x/y/z coordinates of input points
%    DisMin  : minimum distance between output points
%    xyz_out : [3xM] x/y/z coordinates of output points; (M<=N) so dist(xyz_out(:,i),xyz_out(:,j))>DisMin for any i<j
%

xyzout=[];
L=size(xyz,1);
N=size(xyz,2);
sources=zeros(N,L);
nsources=0;
idx=zeros(1,N);
for n1=1:N
    if DisMin>0,
        d=sqrt((sources(1:nsources,1)-xyz(1,n1)).^2+(sources(1:nsources,2)-xyz(2,n1)).^2+(sources(1:nsources,3)-xyz(3,n1)).^2); 
    end
    if DisMin==0||all(d>DisMin)
        nsources=nsources+1;
        sources(nsources,:)=xyz(:,n1)';
        idx(nsources)=n1;
    end
end
xyzout=sources(1:nsources,:)';
idx=idx(1:nsources);