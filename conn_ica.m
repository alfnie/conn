function [S,W]=conn_ica(X,N,varargin)
% fast-ICA algorithm (Aapo Hyvarinen "Fast and Robust Fixed-Point Algorithms for Independent Component Analysis")
%
% [S,W]=CONN_ICA(X,k)
%    X: [NxM] data matrix (each column is an N-dimensional observation vector)
%    k: number of independent components
%    S: independent components
%    W: mixing matrix
%        X = inv(W) * S
%        S = W * X
%
% [S,W]=CONN_ICA(X,k, parameter_name,parameter_value, ...)
%    'method' : choice of contrast function:
%               'tanh' (G1; default)
%               'gauss' (G2)
%               'pow3' (G3)
%    'whiten' : (1/0) apply pre-whitening to columns of X (default: 1)
%    'mu'     : (0-1) stabilization factor step size (default: 1)   
%    'maxiter': maximum number of itereations (default: 1e3)
%    'dodisp' : waitbar display (default: 1)
%

% alfnie@gmail.com

options=struct('method','tanh',...
               'docenter',true,...
               'dowhiten',true,...
               'mu',1,...
               'maxiter',1e3,...
               'miniter',1e1,...
               'maxchange',1e-8,...
               'decreasemuat',.75,...
               'rndseed',false,...
               'dodisp',1);
for n=1:2:numel(varargin), if isfield(options,lower(varargin{n})), options.(lower(varargin{n}))=varargin{n+1}; else error('unknown parameter %s',varargin{n}); end; end
if isempty(options.method), options.method='tanh'; end
if ~ismember(options.method,{'tanh','gauss','pow3'}), error('unknown method %s',options.method); end
[Nd,Ns]=size(X); % columns are samples
if nargin<2||isempty(N), N=Nd; end
if N>Nd, conn_disp('fprintf','conn_ica warning: number of components greater than dimensions of data; estimating %d components instead\n',Nd); N=Nd; end

if options.dodisp==1, hdl=conn_waitbar(0,'computing ICA decomposition. Please wait...'); end
if options.docenter
    mX=mean(X,2);
    X=conn_bsxfun(@minus,X,mX);
else mX=zeros(Nd,1);
end
if options.dowhiten
    [Q,D]=eig(X*X'/(Ns-1));
    iD=diag(1./sqrt(diag(D)));
    X=iD*Q'*X;
end
W=zeros(N,Nd);
P=0;
if ~options.rndseed, randn('seed',0); end
for n=1:N
    wbak=nan;
    w=randn(Nd,1);
    if n>1, w = w-P'*(P*w); end
    w = w/max(eps,sqrt(sum(w.^2)));
    if n<Nd, 
        for niter=1:options.maxiter
            wx=w'*X;
            switch(options.method)
                case 'tanh'
                    gwx=tanh(wx);
                    dgwx=Ns-sum(gwx.^2);
                case 'pow3'
                    gwx=(wx).^3;
                    dgwx=3*Ns;
                case 'gauss'
                    wx2=wx.^2;
                    expwx2=exp(-wx2/2);
                    gwx=wx.*expwx2;
                    dgwx=sum((1-wx2).*expwx2);
            end
            if options.mu==1,
                w = (X*gwx' - dgwx*w)/Ns;
            else
                beta = wx*gwx';
                w = (options.mu*(X*gwx') + ((1-options.mu)*beta - dgwx)*w)/Ns;
            end
            if n>1, w = w-P'*(P*w); end
            w = w/max(eps,sqrt(sum(w.^2)));
            if niter>options.miniter&&min(max(abs(wbak-w)),max(abs(wbak+w)))<options.maxchange, break; end
            wbak=w;
            if niter==ceil(options.maxiter*options.decreasemuat), options.mu=max(.1,options.mu/2); end
        end
        if niter==options.maxiter, conn_disp('conn_ica warning: maximum iteration reached');
            %else conn_disp('fprintf','component %d converged at %d iterations\n',n,niter);
        end
    end
    W(n,:)=w';
    P=W(1:n,:);
    
    if options.dodisp==1, conn_waitbar(n/N,hdl); end
end
S=W*X;
if options.dodisp==1, close(hdl(ishandle(hdl))); end

if options.dowhiten
    W=W*iD*Q';
end
if options.docenter
    S=conn_bsxfun(@plus,S,W*mX);
end
end



