function [E,D]=conn_network_centrality(C)
% CONN_NETWORK_CENTRALITY betweenness centrality 
% D=conn_network_centrality(C);
%  C: NxN thresholded connectivity matrix (for N nodes)
%  D: Nx1 betweenness centrality of each node
%

N=size(C,1);
C=double(C>0 | C'>0);
C(1:N+1:end)=0;
X=sparse(1:N,1:N,true(1,N),N,N);
D=inf(N,N);
D(X)=0;
E=zeros(N,N,N);
EN=full(C);
for n=1:N, 
    Y=(C*X)>0;
    Y=Y&(D>n);
    if ~any(Y(:)),break;end
    if n>1, 
        for m=1:N, 
            E(:,Y(:,m),m)=E(:,X(:,m),m)*C(X(:,m),Y(:,m)); 
            t=diag(EN(m,X(:,m)))*C(X(:,m),Y(:,m));
            E(X(:,m),Y(:,m),m)=t; 
            EN(m,Y(:,m))=sum(t,1);
        end; 
    end
    D(Y)=n;
    X=Y;
end
E=E./repmat(max(eps,shiftdim(EN,-1)),[N,1,1]);
E=sum(sum(E,2),3)/max(eps,(N-1)*(N-2));
