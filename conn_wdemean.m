function y=conn_wdemean(x,w)
w=max(0,w);
if numel(w)==1, w=repmat(w,size(x,1),1); end
y=conn_bsxfun(@minus,x,sum(x.*repmat(w,[1,size(x,2)]),1)./max(eps,sum(w)));