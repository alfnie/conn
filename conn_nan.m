function [x,ok]=conn_nan(x)
idx=find(any(isnan(x),1));
ok=isempty(idx);
if ~ok, % turn nan values into mean values for each column separately
    mx=conn_nanmean(x,1);
    idx2=find(isnan(mx));
    if ~isempty(idx2), % turn columns with all nan's to mean value across all columns
        mx(idx2)=conn_nanmean(mx,2);
        if any(isnan(mx)),mx(isnan(mx))=0;end; % turn matrix with all nan's to 0
    end
    for n1=1:length(idx),x(isnan(x(:,idx(n1))),idx(n1))=mx(idx(n1));end
end

function m = conn_nanmean(x,dim)
nans = isnan(x);
x(nans) = 0;
n = sum(~nans,dim);
n(n==0) = NaN;
m = sum(x,dim) ./ n;

