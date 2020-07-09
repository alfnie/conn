function [c,ac,rac]=conn_network_clustering(C,M)
% CONN_NETWORK_CLUSTERING local clustering coefficient 
% c=conn_network_clustering(C);
%  C: NxN thresholded connectivity matrix (for N nodes)
%  c: Nx1 clustering coefficient of each node
%
% [c,avc,navc]=conn_network_clustering(C); also returns the average
% clustering coefficient network measure (avc) and the normalized average
% clustering coefficient (navc, divided by the equivalent measure on random graph with similar degree)
%

N=size(C,1);
C=C>0;
C(1:N+1:end)=false;
c=nan(N,1);
for n=1:N,
    d=C(n,:);
    nd=sum(d);
    if nd>1, c(n)=sum(sum(C(d,d)))/(nd*(nd-1)); end
end
if nargout>1,
    vc=~isnan(c);
    ac=sum(c(vc))/sum(vc);
end
if nargout>2,
    if nargin<2,M=20;end
    K=sum(C(:))/2;
    idx=find(triu(ones(N,N)-eye(N)));
    rac=zeros(1,M);
    for n=1:M,
        idx0=randperm(numel(idx));
        Cc=zeros(N,N);Cc(idx(idx0(1:K)))=1;Cc=Cc+Cc';
        [nill,rac(n)]=conn_network_clustering(Cc);
    end
    rac=ac/mean(rac);
end

end

function y=randperm(N)
[nill,y]=sort(rand(1,N));
end



