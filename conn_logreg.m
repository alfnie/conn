function [B,Yfit,L]=conn_logreg(X,Y,varargin)
% LOGREG logistic regression
% [B,Yfit]=conn_logreg(X,Y,varargin)
%    X : predictor matrix (#samples x #predictors) (note: include constant term explicitly if needed)
%    Y : outcome matrix (0/1 values #samples x #outcomes)
% fits a model of the form:
%    Yfit = 1./(1+exp(-X*B))
% minimizing the negative log-likelihood function - (Y.*log(Yfit) + (1-Y).*log(1-Yfit))
%
% additional parameters: conn_logreg(...,'fieldname1',fieldvalue1,'fieldname2',fieldvalue2,...)
%   scaling             : scales predictor matrix to constant-norm columns (mean(X.^2,1)=1) (default true)
%   regularization      : regularization factor (ridge regression; cost function = negative log-likelihood + regularization/2*mean(B.^2)) (default 0.001)
%   stepsize            : initial gradient step size (default 1)
%   adaptive            : adaptive gradient step size factor (step-sizes change by multiplication factors of the form (1+adaptive)^k) (default 0.25)
%   display             : display cost function and step-size during optimization (default false)
%   [optimization stop criteria parameters] Run between "miniter" and "maxiter" iterations. Stop iterations if cost function changes less than "tolmin" for at least "toliter" iterations in a row 
%   miniter             : (default 100)
%   maxiter             : (default 5000)
%   tolmin              : (default 1e-7)
%   toliter             : (default 10)
%    
%

% alfnie@bu.edu 

% note: test against brute force search
% Lfun=@(B)mean(mean( -Y.*(X*B)+log(1+exp(X*B)) )); 
% B=fminsearch(Lfun,zeros(size(X,2),size(Y,2))); 
%

options=struct(...
    'stepsize', 1,... % step size
    'regularization', 0.001,... % regularization factor (ridge regression)
    'scaling',true,... % scales X
    'adaptive',0.25,...  % exponential increment/decrement of step size
    'fastadapt',10,... % number of iterations when step-size is allowed to increase fast
    'miniter', 100,... % minimum number of iterations
    'maxiter', 5000,... % maximum number of iterations
    'toliter',10,...
    'tolmin', 1e-7,...     % (between miniter and maxiter: stops iteration if loglikelihood change below this threshold for at least toliter iterations in a row)
    'display',false);
for n=1:2:numel(varargin)
    if isfield(options,lower(varargin{n}))
        options.(lower(varargin{n}))=varargin{n+1};
    else
        fnames=fieldnames(options);
        error('unable to match property %s (valid properties: %s)',varargin{n},sprintf('%s ',fnames{:}));
    end
end
[N1,Nx]=size(X);
[N,Ny]=size(Y);
assert(N==N1,'unequal number of samples (rows) in X and Y data');
assert(all(all(Y==0|Y==1)),'Y values should be 0/1 only');
assert(~nnz(isnan(X))&~nnz(isnan(Y)),'NaN values in X or Y');

B=zeros(Nx,Ny);
stepsize=options.stepsize; % learning rate
slowimprovement=0;
if options.scaling
    sX=sqrt(mean(X.^2,1)); % unit-norm scaling
    X=X./repmat(max(eps,sX),N,1);
end

for iter=1:options.maxiter
    [L,Yfit,dB]=conn_logreg_gradient(X,Y,B,options.regularization);
    if options.display, fprintf('Iteration %d : loss = %f   step-size = %f\n', iter, L, stepsize); end
    if options.adaptive<=0, maxfactor=0; minfactor=0; 
    elseif iter<=options.fastadapt, maxfactor=10; minfactor=-20; % fast stepsize change: from fast increase ((1+adaptive)^10*stepsize), to fast decrease (up to (1+adaptive)^-20*stepsize)
    else maxfactor=1; minfactor=-20; end  % default stepsize change: from slow increase ((1+adaptive)*stepsize), to fast decrease (up to (1+adaptive)^-20*stepsize)
    for factor=maxfactor:-1:minfactor 
        stepsize2=stepsize*((1+max(0,options.adaptive))^factor);
        B2=B+stepsize2*dB; % updates optimal B
        L2=conn_logreg_gradient(X,Y,B2,options.regularization);
        if L2<=L, break; end
    end
    if L2>L, break; end % early exit (no further L decrease found)
    if iter>options.miniter
        if L-L2<options.tolmin, % reached plateau
            slowimprovement=slowimprovement+1;
            if slowimprovement>=options.toliter, break; end
        else slowimprovement=0;
        end
    end
    B=B2;
    stepsize=stepsize2;
end
if iter==options.maxiter, fprintf('warning: exceeded maxiter iterations\n'); end
[L,Yfit]=conn_logreg_gradient(X,Y,B,options.regularization);
if options.scaling, B=B./repmat(max(eps,sX)',1,Ny); end % back to original units

end

function [L,Yfit,dB]=conn_logreg_gradient(X,Y,B,regularization)
z=X*B;
ez=exp(-abs(z));
L=mean(mean(max(z,0)+log1p(ez)-Y.*z)) + regularization/2*mean(B(:).^2); %  L=-mean(mean(Y.*log(max(1e-16,Yfit)) + (1-Y).*log(max(1e-16,1-Yfit)),1),2) + regularization/2*mean(B(:).^2);
if nargout>=2 %  Yfit=1./(1+exp(-z));
    Yfit=1./(1+ez);
    negz=z<0;
    Yfit(negz)=Yfit(negz).*ez(negz);
end    
if nargout>=3, dB=(X'*(Y-Yfit))/size(X,1)-regularization*B; end
end


