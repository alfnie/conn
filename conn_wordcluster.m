function [T,idx,TT] = conn_wordcluster(names,nclusters)
% internal function
% groups names by Levenshtein distance
%
% [T,idx,TT] = conn_wordcluster(names,nclusters)
% 

if nargin<2, nclusters=[]; end

N=numel(names);
Y=conn_wordld(names,names);
I=tril(ones(N),-1);
Y(eye(size(Y))>0)=0;
Y=(Y+Y')/2;
Yb=Y(I>0)';
Z = conn_statslinkage(Yb, 'co');
idx=conn_statsoptimalleaforder(Z,Y); %reference: Bar-Joseph, Z., Gifford, D.K., and Jaakkola, T.S. (2001). Fast optimal leaf ordering for hierarchical clustering. Bioinformatics 17, Suppl 1:S22?9. PMID: 11472989
%[H,t,idx]=conn_statsdendrogram(Z,0,'labels',data.names2(i2),'orientation','left');
if isempty(nclusters) % automatic cutoff number of clusters
    score=(1:N-1)'/(N-1)-Z(:,3)/max(Z(:,3));
    [nill,idxmax]=max(score);
    nclusters=N-idxmax;
else nclusters=NCLUSTERS;
end
T = conn_statscluster(Z, 'maxclust', nclusters);
if nargout>2, 
    TT=[]; 
    for n1=1:min(N,nclusters), % variable number of clusters
        TT(:,n1) = conn_statscluster(Z, 'maxclust', n1); 
    end
end
end
