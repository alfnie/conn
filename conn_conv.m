function [Y,K]=conn_conv(X,s,compute,issurface)
% Y=conn_conv(X,s)
% X: n-dimensional matrix
% s: n-length vector: kernel width (FWHM) in pixels for each dimension of X
% Y: n-dimensional matrix same size as X: convolution of X with non-isotropic gaussian kernel
%

persistent surfparams;
if nargin<3||isempty(compute), compute=2; end
if nargin<4||isempty(issurface), issurface=false; end

sX=size(X);
nX=min(numel(s),numel(sX));

if issurface
    s=ceil((mean(s)/1.2922)^2); % FWHM(mm) to #iterations
    if isempty(surfparams)
        surfparams=load(fullfile(fileparts(which(mfilename)),'utils','surf','surf_top.mat'),'A');
        surfparams.A=surfparams.A*sparse(1:size(surfparams.A,1),1:size(surfparams.A,1),1./max(eps,sum(surfparams.A,1)));
    end
    Y=reshape(X,size(surfparams.A,1),[]);
    for n=1:mean(s)
        Y=surfparams.A*Y;
    end
    Y=reshape(X,sX);
else
    s=max(eps,s/sqrt(8*log(2))); % FWHM(vox) to sigma(vox)
    % if compute>1,
    %     idxnull=X==0;
    %     %X(idxnull)=nan;
    % end
    Y=X;
    if nargout>1, K=cell(1,nX); end
    for ndim=1:nX,
        n=min(round(3*s(ndim)),sX(ndim));
        k=exp(-.5*((-n:n)'/s(ndim)).^2);
        k=k/sum(k);
        if compute,
            if numel(k)>1
                Y(:,:)=convn( cat(1,flipud(Y(1:n,:)),Y(:,:),flipud(Y(end-n+1:end,:))), k, 'valid');
            end
            Y=permute(Y,[2:numel(sX),1]);
            %Y=shiftdim(Y,1);
        end
        if nargout>1, K{ndim}=k; end
    end
    if compute
        Y=permute(Y,[1+numel(sX)-nX:numel(sX),1:numel(sX)-nX]);
        %Y=shiftdim(Y,numel(sX)-nX);
    end
    % if compute>1
    %     Y(idxnull)=0;
    %     %Y(isnan(Y))=0;
    % end
end
