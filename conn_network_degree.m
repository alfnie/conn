function [k,avk]=conn_network_degree(C)
% CONN_NETWORK_DEGREE degree of each node 
% k=conn_network_degree(C);
%  C: NxN thresholded connectivity matrix (for N nodes)
%  k: Nx1 degree of each node
%

C=C>0;
k=sum(C,2)-diag(C);
if nargout>1,
    avk=mean(k);
end
