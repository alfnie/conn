function [D,ad,rad]=conn_network_mindist(C,M)
% CONN_NETWORK_MINDIST minimum-path distance 
% D=conn_network_mindist(C);
%  C: NxN thresholded connectivity matrix (for N nodes)
%  D: NxN minimum-path distance matrix (between each pair of nodes)
%
% [D,avd,navd]=conn_network_mindist(C); also returns the average path
% distance network measure (avd) and the normalized average path distance (navd, divided by the equivalent measure on random graph with similar degree)
%

N=size(C,1);
C=double(C>0);
C(1:N+1:end)=0;
X=sparse(1:N,1:N,true(1,N),N,N);
D=inf(N,N);
D(X)=0;

for n=1:N, 
    X=(C*X)>0;
    X=X&(D>n);
    if ~any(X(:)),break;end
    D(X)=n;
end

if nargout>1,
    vd=D>0&~isinf(D);
    ad=mean(D(vd));             % mean
% 	vd=D(:)>0;
% 	ad=1./mean(1./D(vd));      % harmonic mean
end
if nargout>2,
    if nargin<2,M=20;end
    K=sum(C(:))/2;
    idx=find(triu(ones(N,N)-eye(N)));
    rad=zeros(1,M);
    for n=1:M,
        idx0=randperm(numel(idx));
        Cc=zeros(N,N);Cc(idx(idx0(1:K)))=1;Cc=Cc+Cc';
        [nill,rad(n)]=conn_network_mindist(Cc);
    end
    rad=ad/mean(rad);
%     K=sum(C(:));
%     rad(2)=(log(N)-0.5772156649)/log(K/N)+1/2; % Fronczak et al. 2004
end

end

function y=randperm(N)
[nill,y]=sort(rand(1,N));
end
