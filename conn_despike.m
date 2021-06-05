function yf = conn_despike(yf)
% internal function
if conn_server('util_isremotevar',yf), yf=conn_server('run_keep',mfilename,yf); return; end
my=repmat(median(yf,1),[size(yf,1),1]);
sy=repmat(4*median(abs(yf-my)),[size(yf,1),1]);
yf=my+sy.*tanh((yf-my)./max(eps,sy));

