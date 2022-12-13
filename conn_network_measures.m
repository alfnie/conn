function [measure_names,measure_values]=conn_network_measures(C)
% CONN_NETWORK_MEASURES Computes Graph-theory measures 
% [measure_names,measure_values]=conn_network_measures(C);
%  C: NxN thresholded connectivity matrix (for N nodes)
%

measure_names={'GlobalEfficiency','LocalEfficiency','BetweennessCentrality','ClosenessCentrality','EigenvectorCentrality','Eccentricity','Cost','AveragePathLength','ClusteringCoefficient','Degree'};
if ~nargin
    measure_values={}; 
    return
else
    N=size(C,1);
    C=C>0;
    C(1:N+1:end)=false;
    D=conn_network_mindist(C);
    iD=zeros(size(D));iD(D>0)=1./D(D>0);
    E_global=nan(N,1);E_local=E_global;PathLength=E_global;ClCoef=E_global;K=nan(N,1);Degree=K;Ecc=K;
    [EigV,EigD]=eig(double(C|C'));
    [nill,EigDidx]=max(real(diag(EigD))); EigV=EigV(:,EigDidx); EigV=max(0,EigV/mean(EigV));
    for n=1:N,
        d=find(C(n,:));
        nd=numel(d);
        Degree(n)=nd;
        connected=~isinf(D(n,:))&D(n,:)>0;
        if nnz(connected)>0, 
            PathLength(n)=mean(D(n,connected)); 
            Ecc(n)=max(D(n,connected)); 
        end
        if nd>1,
            ClCoef(n)=sum(sum(C(d,d)))/(nd*(nd-1)); 
            md=conn_network_mindist(C(d,d));
            id=zeros(size(md));id(md>0)=1./md(md>0);
            E_local(n)=sum(sum(id))/(nd*(nd-1));
        end
        E_global(n)=sum(iD(n,:))/(N-1);
        K(n)=nd/(N-1);
    end
    W=conn_network_centrality(C);
    PathLength(isinf(PathLength))=nan;
    measure_values={E_global,E_local,W,1./PathLength,EigV,Ecc,K,PathLength,ClCoef,Degree};
end
