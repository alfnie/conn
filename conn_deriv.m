function varargout=conn_deriv(X,deriv,issurface)

persistent surfparams;
if nargin<3||isempty(issurface), issurface=false; end

sX=size(X);
if issurface
    if isempty(surfparams)
        surfparams=load(fullfile(fileparts(which(mfilename)),'utils','surf','surf_top.mat'),'dX','dY','A');
        surfparams.A(1:size(surfparams.A,1)+1:end)=0;
        surfparams.A=sparse(1:size(surfparams.A,1),1:size(surfparams.A,1),1./max(eps,sum(surfparams.A,2)))*surfparams.A;
    end
    X=reshape(X,size(surfparams.dX,1),[]);
    if deriv==1
        dx=surfparams.dX*X;
        dy=surfparams.dY*X;
        dx=reshape(dx,sX);
        dy=reshape(dy,sX);
        varargout={dx,dy};
    elseif deriv==2
        dx=X-surfparams.A*X;
        dx=reshape(dx,sX);
        varargout={dx};
    end
else
    if deriv==1
        dx=X([2:end,end],:,:)-X([1,1:end-1],:,:);
        dy=X(:,[2:end,end],:)-X(:,[1,1:end-1],:);
        dz=X(:,:,[2:end,end])-X(:,:,[1,1:end-1]);
        varargout={dx,dy,dz};
    elseif deriv==2
        dx=X([2:end,end],:,:)+X([1,1:end-1],:,:);
        dx=dx+X(:,[2:end,end],:)+X(:,[1,1:end-1],:);
        dx=dx+X(:,:,[2:end,end])+X(:,:,[1,1:end-1]);
        varargout={X-dx/6};
    end
end