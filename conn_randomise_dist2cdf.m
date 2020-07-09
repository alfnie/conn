function cdf2 = conn_randomise_dist2cdf(x,xthr)

n=numel(x);
m=numel(xthr);
[nill,idx]=sort([xthr(:)' x(:)']);
idx(idx)=1:numel(idx);
cdf2=idx(1:m);
[nill,idx]=sort(cdf2);
cdf2(idx)=cdf2(idx)-(1:numel(idx));
cdf2=(n-cdf2)/n;
end
% computes cdf from samples
% cdf2(i) = mean(x>=xthr(i));
% but faster for large x,xthr arrays
