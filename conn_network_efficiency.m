function [E_global,E_local,K,e_global,e_local,k]=conn_network_efficiency(C,M)
% CONN_NETWORK_efficiency Global&local efficiency measures 
% [E_global,E_local,Cost,e_global,e_local,cost]=conn_network_efficiency(C);
%  C: NxN thresholded connectivity matrix (for N nodes)
%  E_global: Nx1 global efficiency of each node
%  E_local: Nx1 local efficiency of each node
%  Cost: Nx1 cost of each node
%  e_global: global efficiency of the network
%  e_local: local efficiency of the network
%  cost: cost of the network
%

if nargin>1,
    N=C(1);N0=round(C(2)*N*(N-1)/2);
%     N=size(C,1);
%     C=C>0;
%     C(1:N+1:end)=false;
%     N0=ceil(sum(C(:))/2);
%     t=ceil(sqrt(N));t2=reshape((1:t^2),[t,t]);Cl=zeros(t^2);for n1=1:t,for n2=1:t,Cl(t2(n1,n2),t2(n1,:))=1;Cl(t2(n1,n2),t2(:,n2))=1;end;end;Cl=(Cl+Cl')>0;Cl(1:t^2+1:end)=0;
    if M>0,
        % random network
        idx=find(triu(ones(N,N)-eye(N)));
        e_global=0;e_local=0;k=0;E_global=0;E_local=0;K=0;
        for n=1:M,
            N0=round(C(1+ceil((numel(C)-1)*rand))*N*(N-1)/2);
            idx0=randperm(numel(idx));
            Cc=zeros(N,N);Cc(idx(idx0(1:N0)))=1;Cc=Cc+Cc';
            [nill,nill,nill,t1,t2,t3]=conn_network_efficiency(Cc);
            e_global=e_global+t1;e_local=e_local+t2;k=k+t3;
        end
    else
        M=abs(M);
        % lattice
        idx=find(triu(ones(N,N)-eye(N)));
        [idxi,idxj]=ind2sub([N,N],idx);
        t=abs(idxi-idxj);
        t=min(t,N-t);
        e_global=0;e_local=0;k=0;E_global=0;E_local=0;K=0;
        for n=1:M,
            N0=round(C(1+ceil((numel(C)-1)*rand))*N*(N-1)/2);
            [nill,idx0]=sort(t+rand(size(t))); %randperm(numel(idx));
            Cc=zeros(N,N);Cc(idx(idx0(1:N0)))=1;Cc=Cc+Cc';
            [nill,nill,nill,t1,t2,t3]=conn_network_efficiency(Cc);
            e_global=e_global+t1;e_local=e_local+t2;k=k+t3;
        end
    end
    e_global=e_global/M;e_local=e_local/M;k=k/M;
    return;
end

N=size(C,1);
C=C>0;
C(1:N+1:end)=false;
D=conn_network_mindist(C);
iD=zeros(size(D));iD(D>0)=1./D(D>0);
E_global=zeros(N,1);E_local=E_global;K=zeros(N,1);
for n=1:N,
    d=find(C(n,:));
    nd=numel(d);
    if nd>1, 
%         E_local(n)=sum(sum(iD(d,d)))/(nd*(nd-1)); 
        d=conn_network_mindist(C(d,d));
        id=zeros(size(d));id(d>0)=1./d(d>0);
        E_local(n)=sum(sum(id))/(nd*(nd-1)); 
    end
    E_global(n)=sum(iD(n,:))/(N-1);
    K(n)=nd/(N-1);
end
e_global=mean(E_global);
e_local=mean(E_local);
k=mean(K);
%K=conn_network_centrality(C);
end

function y=randperm(N)
[nill,y]=sort(rand(1,N));
end
