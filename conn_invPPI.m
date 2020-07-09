function [H1,B1,H0,B0,dE,R]=conn_invPPI(X,Nk,X0,X1,DOPLOT,DOICA,FLIPSIGN,RNDSEED)
% EM iterative method for solving the main&interaction factors h in PPI equation:
% X_nj = X_ni*a_ij + sum_k (X_ni.*h_nk)*b_ikj + noise
%   where
%     i and j are indexes over ROIs
%     n is an index over scans/time
%     X_ni is a [Nscans x Nrois] matrix representing the ROI-level BOLD timeseries
%     k is an index over the unknown interaction terms
%     a_ij is a matrix characterizing the main connectivity effects
%     h_nj and b_ikj are matrices characterizing the interaction effects
%       (h: temporal modulation; b: ROI-to-ROI modulation)
%
% [H,B]=conn_invPPI(X,Nk)
%   X:  [Nscans x Nrois] matrix
%   Nk: desired number of components
% returns:
%   H:  [Nscans x Nk] interaction temporal terms
%   B:  [Nrois x Nrois x Nk] interaction ROI-to-ROI connectivity terms
%
% [H,B]=conn_invPPI(...,X0,X1)
%   X0: control for subject/session effects (by default X0 is fixed to ones(Nt,1); when entering multiple subjects/session use kron(eye(Nsubjects),ones(Nscans,1)) to remove main between-subject effects )
%   X1: temporal smoothing matrix
%

lambda=0;               % regularization term (optional)
stopcriterion=1e-4;     % converge criterion (minimum h param change between two consecutive iterations)
miniter=8;              % minimum number of iterations
maxiter=1e2;            % maximum number of iterations
[Nt,Nr]=size(X);
if nargin<2||isempty(Nk), Nk=4; end
if nargin<3||isempty(X0), X0=ones(Nt,1); end
if nargin<4||isempty(X1), X1=speye(Nt); end
if nargin<5||isempty(DOPLOT), DOPLOT=0; end
if nargin<6||isempty(DOICA), DOICA=1; end
if nargin<7||isempty(FLIPSIGN), FLIPSIGN=true; end
if nargin<7||isempty(RNDSEED), RNDSEED=false; end
REMOVEX0VAR=true;
Nk0=size(X0,2);
Nk=min(Nk,Nr*(Nr-1)/2);
Nk1=Nk;
Nk=Nk+Nk0;

mX1=mean(X1,1);
cYY=repmat(mean(mX1*X.^2,2),[2,Nr]);
hist_err=nan(Nk,maxiter*Nk);
hist_err(1)=cYY(1);
nhist_err=1;
B=zeros(Nr,Nk,Nr);
dE=nan(1,Nk);
if DOPLOT==1, hdl=conn_waitbar(0,'Computing dynamic decomposition. Please wait...');
elseif DOPLOT==2, hfig=figure('units','norm','position',[.3 .3 .4 .4],'color','w','name','estimation of temporal modulation factors','numbertitle','off','menubar','none');
end

err_old=nan;
idxnroi=1:Nr;
%if nk>Nk0&&alpha<1, idxnroi=randperm(Nr); idxnroi=idxnroi(1:ceil(numel(idxnroi)/2));end
idxnt=reshape(find(any(X~=0,2)),1,[]); 
X1X=X1*(X.^2);
X0b=X0~=0;
if ~RNDSEED, randn('seed',0); end
H=zeros(Nt,Nk);
H(idxnt,:)=randn(numel(idxnt),Nk);
if Nk0>0, H(:,1:Nk0)=X0; end

H_old=H;
for n=1:maxiter
    for nk=Nk0+1:Nk
        h=H(:,nk);
        if nk>1, h=h-H(:,1:nk-1)*(H(:,1:nk-1)\h); end
        h=h/max(eps,sqrt(mean(h.^2)));
        H(:,nk)=h;
    end
    % estimate B
    db1=zeros(Nk,Nr,Nr);
    db2=zeros(Nk,Nk,Nr);
    for nroi=1:Nr % (loop over sources)
        Y=X;            % target ROIs
        x=X(:,nroi);    % source ROI
        if n==1, cYY(1,nroi)=mean(mX1*Y.^2,2); end
        xy=sparse(1:Nt,1:Nt,x)*Y;
        cXY=H'*X1*xy;
        cXX=H'*sparse(1:Nt,1:Nt,X1*xy(:,nroi))*H;
        db1(:,:,nroi)=cXY;
        db2(:,:,nroi)=cXX;
    end
    for nroi1=1:Nr,
        for nroi2=1:Nr
            cXY=db1(:,nroi1,nroi2)+db1(:,nroi2,nroi1);
            cXX=db2(:,:,nroi1)+db2(:,:,nroi2);
            dcXX=cXX(1:Nk+1:end); maxdcXX=max(dcXX); dcXX(dcXX<1e-5*maxdcXX)=1e-5*maxdcXX;cXX(1:Nk+1:end)=dcXX; 
            B(nroi1,:,nroi2)=inv(cXX)*cXY; % ROI-to-ROI modulation
        end
    end
    % estimate H
    dh1=0; % [Nt Nk]
    err=0;
    for nroi=idxnroi, %1:Nr
        Y=X;            % target ROIs
        x=X(:,nroi);    % source ROI
        if REMOVEX0VAR
            b=B(:,:,nroi);
        else
            Y=Y-H(:,1:Nk0)*B(:,1:Nk0,nroi)';
            b=B(:,Nk0+1:end,nroi);
        end
        xy=sparse(1:Nt,1:Nt,x)*Y;
        dh1=dh1+xy*b;                               % cumulative term for computing temporal modulation
        if DOPLOT>1
            fitmse=cYY(1,nroi)-(X1*xy(:,nroi))'*sum((H*B(:,:,nroi)').^2,2)/Nr/Nt;           % current mse for each source ROI
            cYY(2,nroi)=fitmse;
            err=err+fitmse/Nr;
        end
    end
    dh1=X1*dh1;
    
    if REMOVEX0VAR, Nkt=Nk; else Nkt=Nk1; end
    bb=zeros([Nkt Nkt Nr]); % [Nk Nk Nr]
    for nroi=idxnroi, %1:Nr
        if REMOVEX0VAR, b=B(:,:,nroi);
        else b=B(:,Nk0+1:end,nroi);
        end
        bb(:,:,nroi)=b'*b;
    end
    bb=reshape(bb,[Nkt^2,Nr]);
    idxterms=true(1,Nk);
    for nt=idxnt,
        dh2=reshape(bb*X1X(nt,:)',[Nkt Nkt]);
        if REMOVEX0VAR
            if Nk0>0, idxterms(1:Nk0)=X0b(nt,:); end
            if any(dh1(nt,idxterms))
                ht=inv(Nr^2*lambda+dh2(idxterms,idxterms))*dh1(nt,idxterms)';
                H(nt,idxterms)=ht;               % temporal modulation
            else H(nt,idxterms)=0;
            end
        else
            if any(dh1(nt,:))
                ht=inv(Nr^2*lambda+dh2)*dh1(nt,:)';
                H(nt,Nk0+1:end)=ht;
            else H(nt,Nk0+1:end)=0;
            end
        end
    end
    if Nk0>0, H(:,1:Nk0)=X0; end
    
    nhist_err=nhist_err+1;
    if DOPLOT>1, hist_err(1,nhist_err)=err; end %err_last-err;
    stopnow=n>miniter & max(max(abs(H-H_old)))<stopcriterion;
    if DOPLOT>1&&(~rem(n,1)||stopnow) % plots
        subplot(211); plot(fliplr(hist_err(:,1:nhist_err)'),'.-'); set(gca,'xcolor',.75*[1 1 1],'ycolor',.75*[1 1 1]); xlabel('iterations'); ylabel('mse');
        subplot(223); plot(H); set(gca,'xcolor',.75*[1 1 1],'ycolor',.75*[1 1 1]); title('temporal modulation');
        subplot(224); imagesc(reshape(permute(B(:,:,:),[1,3,2]),Nr,[])); set(gca,'xcolor',.75*[1 1 1],'ycolor',.75*[1 1 1]); axis equal tight; title('ROI-to-ROI modulation'); colorbar;
        %disp(max(max(abs(H-H_old))));
        drawnow;
    end
    if DOPLOT==1, conn_waitbar(n/maxiter,hdl); end
    if stopnow, break; end
    H_old=H;
end

if DOPLOT==1, close(hdl(ishandle(hdl)));
elseif DOPLOT==2, close(hfig(ishandle(hfig)));
end
%rand('state',state);

% removes covariates-of-no-interest from returned matrices
H0=H(:,1:Nk0);
B0=permute(B(:,1:Nk0,:),[1 3 2]);
idx1=Nk0+1:Nk;
H1=H(:,idx1);
B1=permute(B(:,idx1,:),[1 3 2]);

if DOICA&&Nk1>1
    if DOICA==1, % spatial-ICA
        [S,W]=conn_ica(reshape(B1,[Nr*Nr Nk1])',[],'dodisp',DOPLOT,'rndseed',RNDSEED);
        B1=reshape(S',[Nr Nr Nk1]);
        H1=H1*pinv(W);
    else % temporal-ICA
        [S,W]=conn_ica(H1',[],'dodisp',DOPLOT,'rndseed',RNDSEED);
        H1=S';
        B1=reshape( reshape(B1,[Nr*Nr Nk1])*pinv(W), [Nr Nr Nk1]);
    end
end

% nk=sqrt(mean(H1.^2,1));
% H1=conn_bsxfun(@rdivide,H1,nk);
% B1=conn_bsxfun(@times,B1,shiftdim(nk,-1));

if FLIPSIGN
    signB=sum(sum(B1.^3,1),2)<0;
    B1(:,:,signB)=-B1(:,:,signB);
    H1(:,signB)=-H1(:,signB);
end

if nargout>5
    R=zeros(size(B1,1),size(B1,2),size(H1,1));
    for n=1:size(H1,1)
        r=0;
        for nk=1:size(H0,2)
            r=r+H0(n,nk)*B0(:,:,nk);
        end
        for nk=1:size(H1,2)
            r=r+H1(n,nk)*B1(:,:,nk);
        end
        R(:,:,n)=r;
    end
end






