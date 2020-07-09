function q=conn_fdr(p,dim,thr)
% CONN_FDR False Discovery Rate
% Q=CONN_FDR(P); returns vector Q of estimated false discovery rates (set-level q-values) from 
% a vector P of multiple-test false positive levels (uncorrected p-values)
% Q=CONN_FDR(P,dim); where P is a matrix computes the FDR along the dimension dim of P
%

if nargin<2||isempty(dim), 
    if sum(size(p)>1)==1,dim=find(size(p)>1);
    else, dim=1; end
end
nd=length(size(p)); 
if dim~=1, p=permute(p,[dim,1:dim-1,dim+1:nd]); end

sp=size(p);
q=nan(sp);
N0=sp(1);
N2=prod(sp(2:end));
if nargin>2
    [sp,idx]=sort(p,1);
    spvalid=~isnan(sp);
    N1=sum(spvalid,1);
    i=sp<=conn_bsxfun(@rdivide,thr*(1:N0)',N1);
    sp(~i)=-inf;
    q=conn_bsxfun(@le,p,max(sp,[],1));
else
    [sp,idx]=sort(p,1);
    spvalid=~isnan(sp);
    N1=sum(spvalid,1);
    for n2=find(N1(:)'>0),
        n1=N1(n2);
        qt=min(1,n1*sp(spvalid(:,n2),n2)./(1:n1)');
        min1=nan;
        for n=n1:-1:1,
            min1=min(min1,qt(n));
            q(idx(n,n2),n2)=min1;
        end
    end
end
if dim~=1, q=ipermute(q,[dim,1:dim-1,dim+1:nd]); end

%for n=1:N,q(idx(n))=min(qt(n:end));end


