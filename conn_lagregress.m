function [Y_fit,b,lag, fy, fymask]=conn_lagregress(X,Y,varargin)
% LAGREGRESS fits lagged regression model
%
% [Y_fit,B,lag]=conn_lagregress(X,Y)
%
% given Y[n,m] and X[n,k] 
% find parameters lag[m] and B[k,m] that best approximate
%       Y[n,m] = sum_k  X[n-lag[m],k]*B[k,m]
%
% [Y_fit,B,lag]=conn_lagregress(X,Y,'select',idx) limits lag terms to only selected columns of X
%       Y[n,m] = sum_k_in_idx  X[n-lag[m],k]*B[k,m] + sum_k_notin_idx X[n,k]*B[k,m]
% [Y_fit,B,lag]=conn_lagregress(..., 'maxdn',N)  defines maximum lag N in samples (default size(Y,1))
% [Y_fit,B,lag]=conn_lagregress(..., 'resdn',K)  defines sub-sample resolution (default 10)
%


params=struct(...
    'select',1:size(X,2),... % select regressors that have lag terms (default all)
    'maxdn',max(ceil(size(X,1)/2),size(Y,1)),...  % maximum lag in samples (lag = -maxdn:maxdn) (default Nt/2)
    'resdn',[],...           % sub-sample resolution (number of subsamples per sample) (default 10)
    'omit',[],...            % omit these regressors in final Y_fit estimation (default [])
    'interp','sym' ...  % 'sym', 'wrap', or 'linear' sub-sample interpolation and external extrapolation modes
    );
if nargin>2
    if isnumeric(varargin{1}), varargin=[{'maxdn'} varargin]; end
    for n1=1:2:numel(varargin)-1, if ~isfield(params,lower(varargin{n1})), error('unknown option %s',lower(varargin{n1})); else params.(lower(varargin{n1}))=varargin{n1+1}; end; end
end
if isempty(params.resdn), params.resdn=max(10,ceil(50/params.maxdn)); end

[Nnx,Nx]=size(X);
[Nn,Ny]=size(Y);
assert(Nnx>=Nn,'mismatched number of rows/timepoints in X and Y');
nX=sqrt(mean(abs(X.^2),1));
X=X./repmat(max(eps,nX),size(X,1),1); 
iswrap=strcmpi(params.interp,'wrap');
issym=strcmpi(params.interp,'sym');
if iswrap
    X=fft(X);
elseif issym
    X=[flipud(X(1:ceil(size(X,1)/2),:));X;flipud(X(ceil(size(X,1)/2)+1:end,:))];
    %X=fft([flipud(X); X; flipud(X)]); 
    %w=[0:ceil(Nn/2)-1 ceil(Nn/2)+zeros(1,3*Nn-2*ceil(Nn/2)) ceil(Nn/2)-1:-1:0]'/ceil(Nn/2);
    %X=fft( (1-w)*mean(X,1) + repmat(w,1,Nx).*interp1((0:Nn-1)',X,max(-Nn/2,min(2*Nn-1-Nn/2,-Nn:2*Nn-1))','linear','extrap') );
end
bbest=nan(Nx,Ny); Y_fit=nan(size(Y)); ilag=nan(1,Ny); mine=inf(1,Ny);
if nargout>3, fy=zeros(size(X,1),Ny); end
lag=0:1/params.resdn:abs(params.maxdn);
if numel(lag)>1, lag=[lag, -lag(2:end)]; end
f=mod(floor((size(X,1)-1)/2)+(0:size(X,1)-1),size(X,1))'-floor((size(X,1)-1)/2);
t=(0:size(X,1)-1)';
maskx=(1:Nn)+ceil((size(X,1)-Nn)/2);
for nlag=1:numel(lag)
    fx=X;
    if iswrap||issym
        fx(:,params.select)=fx(:,params.select).*repmat(exp(-1i*2*pi*f*lag(nlag)/size(fx,1)),1,numel(params.select));
        fx=real(ifft(fx));
    else
        fx(:,params.select)=interp1(t,fx(:,params.select),t-lag(nlag),'linear','extrap');
    end
    fxr=fx(maskx,:);
    fx2=real(fxr'*fxr);
    b=pinv(fx2)*real(fxr'*Y); 
    e=-sum(b.*(fx2*b),1);       % note:  e = sum(abs(Y-fxr*b).^2,1) = sum(abs(Y).^2,1)-sum(conj(b).*(fx2*b),1);
    ok=e<mine;
    if any(ok), 
        mine(ok)=e(ok); 
        bbest(:,ok)=b(:,ok); ilag(ok)=nlag; 
        if any(params.omit), b(params.omit,ok)=0; end
        Y_fit(:,ok)=fxr*b(:,ok); 
        if nargout>3, fy(:,ok)=fx*bbest(:,ok); end
    end
end

lag=lag(ilag);
if nargout>3 % clean Y timeseries: removed not-lag regressors (those not in 'select' or 'omit'), and de-lagged (Y[n+lag[m],m])
    %fy=zeros(size(X,1),Ny);
    fy(maskx,:)=Y;
    fixed=setdiff(1:Nx,[params.select,params.omit]);
    if ~isempty(fixed), fy=fy-X(:,fixed)*bbest(fixed,:); end
    if iswrap||issym
        fy=fft(fy);
        fy=fy.*exp(+1i*2*pi*f*lag/size(fy,1));
        fy=real(ifft(fy));
        fymask=1+mod(bsxfun(@plus,t,lag),size(X,1)); fymask=max(0,min(1,min(fymask-maskx(1)+1,maskx(end)+1-fymask)));
    else
        fy=interp1(t,fy,t+lag(nlag),'linear','extrap');
        fymask=1+mod(bsxfun(@plus,t,lag),size(X,1)); fymask=max(0,min(1,min(fymask-maskx(1)+1,maskx(end)+1-fymask)));
    end
    %fy=fy(maskx,:);
end
b=bbest./repmat(max(eps,nX'),1,Ny);

end


