function cdf2 = conn_randomise_hist2cdf(pdf,xthr)

m=numel(xthr);
if ~m, cdf2=[]; return; end 
idx=find(pdf);
n=numel(idx);
cdf=full(cumsum(pdf(idx)));
cdf2=zeros(size(xthr));
mask=xthr>1;
[nill,idx3]=sort([idx(:)' xthr(:)'-1]);
mask=idx3>n;
mask=mask|[diff(mask)>0 false];
last=0;
for n1=idx3(mask), 
    if n1<=n, last=cdf(n1);
    else cdf2(n1-n)=last; 
    end
end
cdf2=cdf(end)-cdf2;
end
% resamples histogram to compute cdf
% cdf2(i)=sum(pdf(xthr(i):end));
% but faster for large pdf,xthr arrays, sparse pdf
