function varargout=conn_clusterregress(X,Y,varargin)
% Ensemble Clustered Regression (ECR)
% Y(X) = A + X*B + noise
% for robust regression when size(X,2)>>size(X,1)
% with clustering over columns of B and split-half 
% resampling for variable selection and fitting
%
% [Model,Yfit]=conn_clusterregress(X,Y);
% Y                       : outcome variables [# samples x # outcomes]
% X                       : predictor variables [# samples x # of predictor variables]
% Model.parameters.A      : model bias term [1 x # outcomes] 
% Model.parameters.B      : model regressor coefficients [# of predictor variables x # outcomes] 
% Model.parameters.B_std  : model regressor standardized coefficients (scaled to unit standard deviation predictor and outcome measures) [# of predictor variables x # outcomes] 
% Model.fit.MSE_std       : model out-of-sample standardized Mean Square Error (scaled to unit standard deviation outcome measures)
% Model.predict           : predictor function for new samples | predict(X) = A+B*X with model parameters A&B approximates expected outcome data Y for a new sample X
% Yfit                    : outcome variable model fits for training data [# samples x # outcomes]
%
% [Model]=conn_clusterregress(X,Y);
% Use this option to skip the outer LOO Cross-validation loop needed to estimate Yfit (model fits for the training dataset)
% This will also output a faster but slightly less accurate estimate of out-of-sample MSE using partial estimates computed during training before bagging
%
% conn_clusterregress(X,Y,'paramName',paramValue,...);
% Use parameter name/value pairs to specify additional input arguments. Valid parameter names are:
%       pvariables          : proportion of predictor variables to consider for clustering (default 0.5)
%       ngroups             : number of clusters of predictor variables (default: min( size(X,1)/10, size(X,2)*pvariables )
%       ptrain              : proportion of samples used for variable selection (default 0.5)
%       nbagging            : number of base learners (default size(X,1); note: nbagging is always rounded to an integer multiple of size(X,1))
%       nouterloop          : outer loop number of iterations (adapt cluster means) (default 10)
%       ninnerloop          : inner loop number of iterations (adapt cluster sizes) (default 2)
%       covariates          : covariates matrix [# samples x # covariate variables] Z in a new model of the form: Y(X) = A + X*B + Z*C + noise (note: covariates are always included in all models; they are not part of feature selection)
%       display             : display info while fitting (default false)
%       warning             : display warnings while fitting (default true)
%
% example:
%  N=100;  % samples
%  R=1000; % predictors
%  P=0.1;  % percentage true effects
%  y=randn(2*N,1); % training+testing outcome
%  x=[.25*y+randn(2*N,ceil(R*P)),randn(2*N,R-ceil(R*P))]; % training+testing predictors
%  [A,B]=conn_clusterregress(x(1:N,:),y(1:N,:)); % training set
%  yfit=A+x*B;
%  subplot(121); bar(B); xlabel 'Predictors'; ylabel 'Regressor Coefficients'
%  subplot(222); plot(y(1:N,:),yfit(1:N,:),'ko','markerfacecolor','k'); grid on; xlabel target; ylabel prediction; title 'training set'
%  subplot(224); plot(y(N+1:end,:),yfit(N+1:end),'ko','markerfacecolor','k'); grid on; xlabel target; ylabel prediction; title 'testing set'
%  mse=mean(abs(y(N+1:end)-yfit(N+1:end)).^2); % testing set
%  disp(mse);

options=struct(...
    'pvariables',[],...       % proportion of predictor variables to keep (default 50%)
    'ngroups',[], ...         % maximum number of clusters among kept predictor variables (default: min( size(X,1)/10, size(X,2)*pvariables )
    'display',false,...       % display info while fitting
    'warning',false,...        % display warnings while fitting
    'ptrain',[],...           % proportion of samples used for variable selection (default 50%)
    'nbagging',[],...         % number of base learners (default size(X,1); note: nbagging is always rounded to an integer multiple of size(X,1))
    'skipcentering',false,... % skips centering step
    'covariates',[],...       % covariate matrix [# samples x # covariate variables]
    'nouterloop',10,...       % outer loop number of iterations (adapt cluster means)
    'ninnerloop',2 ...        % inner loop number of iterations (adapt cluster sizes)
    );
for n=1:2:numel(varargin)
    if isfield(options,lower(varargin{n}))
        options.(lower(varargin{n}))=varargin{n+1};
    else
        fnames=fieldnames(options);
        error('unable to match property %s (valid properties: %s)',varargin{n},sprintf('%s ',fnames{:}));
    end
end
Z=options.covariates;
if isempty(Z), Z=zeros(size(Y,1),0); end
if isempty(X), X=zeros(size(Y,1),0); end

% disregards predictors with constant values for all samples
% disregards samples with NaN in outcome variables or with all NaN predictor values
validS=~all(isnan(X),2) & all(~isnan(Y),2); % valid samples
if any(~validS), X=X(validS,:); Y=Y(validS,:); Z=Z(validS,:); if options.warning, fprintf('note: disregarded %d samples with NaN/missing values\n',nnz(~validS)); end; end
validX=std(X,1,1)>0; % valid predictors
if any(~validX), X=X(:,validX); if options.warning, fprintf('note: disregarded %d predictors with NaN/missing or constant values\n',nnz(~validX)); end; end

% checks input arguments
[N,Nx]=size(X);
[N2,Ny]=size(Y);
[N3,Nz]=size(Z);
assert(N2==N,'unequal number of samples in predictor and ouctome matrices');
assert(isempty(Z)|N3==N,'unequal number of samples in covariate matrix');
%obsolete (uncomment if using selectVars_v1) if Ny>1, disp('warning: expected vector output variable; variable selection will focus on first outcome variable only'); end
if isempty(options.pvariables), options.pvariables=min(1,max(1,N/Nx)); end % 50% (up to 100% when number of predictors is low / standard regression)
if isempty(options.ptrain), options.ptrain=min(.5,Nx/N); end % 50% (down to 0% when number of predictors is low / standard regression)
if isempty(options.nbagging), options.nbagging=N; end
options.nbagging=N*max(1,round(options.nbagging/N));
if isempty(options.ngroups), options.ngroups=ceil(N/10); end % minimum 10:1 samples:predictor ratio
options.ngroups=min(options.ngroups, ceil(Nx*options.pvariables));
Yfit=[];

if nargout>1 % use CV to estimate Yfit on training data
    Yfit=zeros([N,Ny]);
    Nfit=zeros(N,1);
    nsamples=min(N-1,ceil(N*.20)); % uses 80% of the samples for training; 20% of the samples for testing
    state=rand('seed'); rand('seed',0);
    for n1=1:4 % we want each sample included 4 times in a testing dataset
        samples=randperm(N); 
        for n2=0:floor(N/nsamples)
            test=samples(nsamples*n2+1:min(N,nsamples*(n2+1)));
            train=setdiff(1:N,test);
            if isempty(Z), 
                Model=conn_clusterregress(X(train,:),Y(train,:),varargin{:});
                Yfit(test,:)=Yfit(test,:) + Model.parameters.A + X(test,:)*Model.parameters.B;
            else 
                Model=conn_clusterregress(X(train,:),Y(train,:),varargin{:},'covariates',Z(train,:));
                Yfit(test,:)=Yfit(test,:) + Model.parameters.A + X(test,:)*Model.parameters.B + Z(test,:)*Model.parameters.C;
            end
            Nfit(test)=Nfit(test)+1;
        end
    end
    Yfit=Yfit./repmat(Nfit,1,Ny);
    clear Model;
    Model.fit.MSE=mean(mean(abs(Yfit-Y).^2,1),2);
    Model.fit.MSE_std=mean(mean(abs(Yfit-Y).^2,1)./var(Y,1,1),2);
    if any(~validS), temp=zeros(numel(validS),Ny); temp(validS,:)=Yfit; Yfit=temp; end
    Model.fit.Yfit=Yfit;
end
%     Model=conn_clusterregress(X,Y,varargin{:});
%     if any(~validX),
%         temp=zeros(numel(validX),Ny); temp(validX,:)=Model.parameters.B_std; Model.parameters.B_std=temp;
%         temp=zeros(numel(validX),Ny); temp(validX,:)=Model.parameters.B; Model.parameters.B=temp;
%     end
%     if isempty(Z), Model.predict=@(x,varargin)Model.parameters.A+x*Model.parameters.B;
%     else Model.predict=@(x,z,varargin)Model.parameters.A+x*Model.parameters.B+z*Model.parameters.C;
%     end
%     varargout={Model,Yfit};
%     return
% end

if options.skipcentering
    mX=zeros(1,Nx);
    mY=zeros(1,Ny);
    mZ=zeros(1,Nc);
    sX=ones(1,Nx);
    sY=ones(1,Ny);
    sZ=ones(1,Nz);
else % removes mean & scales to unit std
    mX=mean(X,1);
    mY=mean(Y,1);
    mZ=mean(Z,1);
    sX=std(X,0,1);
    sY=std(Y,0,1);
    sZ=std(Z,0,1);
    X=(X-repmat(mX,N,1))./repmat(max(eps,sX),N,1);
    Y=(Y-repmat(mY,N,1))./repmat(max(eps,sY),N,1);
    if ~isempty(Z), Z=(Z-repmat(mZ,N,1))./repmat(max(eps,sZ),N,1); end
end
if ~isempty(Z) % removes variability associated with covariates
    C=Z\Y;
    Y=Y-Z*C; 
end

B=zeros([Nx,Ny,options.ngroups]);
err=zeros(1,options.ngroups);
errnull=0;
samples=1:N; % placeholder for nested CV
ns=max(1,round(options.ptrain*N));
state=rand('seed'); rand('seed',0);
for nrepeat=1:options.nbagging
    if options.display, fprintf('.'); end
    idx=randperm(N);
    TRAIN1=samples(idx(1:ns));      % use these samples for variable selection
    TRAIN2=samples(idx(ns+1:N));    % use these samples for parameter estimates
    G=selectVars(X(TRAIN1,:),Y(TRAIN1,:),ceil(Nx*options.pvariables),options.ngroups,options.nouterloop,options.ninnerloop); % [Nx x ngroups] 0/1 group matrix (note: sum-normed columns)
    xg=X(TRAIN2,:)*G;
    syy=sum(sum(Y(TRAIN2,:).^2,1),2);
    for ngroups=1:options.ngroups % test different number of subgroups among kept predictor variables
        %b=xg(:,1:ngroups)\Y(TRAIN2,:);
        if ngroups<=options.ngroups,  % first n groups separtely
            tG=G(:,1:ngroups);
            txg=xg(:,1:ngroups);
        else % all groups as a single group
            tG=double(any(G,2)); tG=tG/nnz(tG);
            txg=X(TRAIN2,:)*tG; 
        end
        txgy=txg'*Y(TRAIN2,:);
        b=pinv(txg'*txg)*txgy;
        terr=(syy-sum(sum(b.*txgy,1),2))/(N-ns); % sum of MSE over outcome variables: sum(mean(abs(Y(TRAIN2,:)-X(TRAIN2,:)*Gb).^2,1),2)
        Gb=tG*b;
        B(:,:,ngroups)=B(:,:,ngroups)+Gb;
        err(ngroups)=err(ngroups)+terr;
    end
    errnull=errnull+syy/(N-ns);
end
rand('seed',state);
err=err/options.nbagging/Ny;
errnull=errnull/options.nbagging/Ny;
if options.display, fprintf('mse = %s : # groups %s \n',mat2str(err),mat2str(1:options.ngroups)); end
aic=log(err)*ns+2*[1:options.ngroups];
[nill,idx]=min(aic); % minimum AIC to determine ngroups
Model.parameters.B_std=B(:,:,idx)/options.nbagging;
Model.parameters.B=repmat(1./max(eps,sX'),1,Ny).*Model.parameters.B_std.*repmat(sY,Nx,1);
if isempty(Z), 
    Model.parameters.C=[];
    Model.parameters.A=mY-mX*Model.parameters.B;
else 
    Model.parameters.C=repmat(1./max(eps,sZ'),1,Ny).*C.*repmat(sY,Nz,1);; 
    Model.parameters.A=mY-mX*Model.parameters.B-mZ*Model.parameters.C;
end
Model.error.MSE_std=err(idx);
Model.error.NULL_std=errnull;
if Ny==1, Model.error.MSE=Model.error.MSE_std*sY.^2; end
if any(~validX), 
    temp=zeros(numel(validX),Ny); temp(validX,:)=Model.parameters.B_std; Model.parameters.B_std=temp; 
    temp=zeros(numel(validX),Ny); temp(validX,:)=Model.parameters.B; Model.parameters.B=temp; 
end
if isempty(Z), Model.predict=@(x,varargin)Model.parameters.A+x*Model.parameters.B;
else Model.predict=@(x,z,varargin)Model.parameters.A+x*Model.parameters.B+z*Model.parameters.C;
end
varargout={Model,Yfit};

end

function G=selectVars(x,y,Nvars,Ngroups,nouterloop,ninnerloop)
% keep "best" Nvars predictor variables, divides them into Ngroups groups, and sort them by largest-absolute-effect
N=size(x,2);
Ny=size(y,2);
b=x'*y;
[sb,idx]=sort(sum(abs(b).^2,2),'descend');
KEEP=idx(1:Nvars);
b=b(KEEP,:);
[sb,idx]=sort(b(:,1)); idx=idx(ceil(Nvars*(2*(1:Ngroups)-1)/(2*Ngroups))); 
%idx=ceil(Nvars*rand(Ngroups,1));
mb=b(idx,:)+eps*randn(Ngroups,1); % means in each cluster
ab=zeros(Ngroups,1);              % size-bias in each cluster
keps=NaN;

for nrepeat1=1:nouterloop
    %d=zeros(Nvars,Ngroups); for n=1:Ngroups, d(:,n)=sum(abs(b-mb(n+zeros(Nvars,1),:)).^2,2); end % note: use this syntax in older Matlab versions
    d=zeros(Nvars,Ngroups); for n=1:Ngroups, d(:,n)=sum(abs(b-mb(n,:)).^2,2); end
    if ninnerloop<1, [mind,T]=min(d,[],2);
    else
        for nrepeat2=1:ninnerloop
            [mind,T]=min(d+repmat(ab',Nvars,1),[],2);
            nb=accumarray(T,1,[Ngroups,1]); %disp(nb')
            if isnan(keps), keps=mean(mind)*Ngroups/Nvars; end
            ab=ab+keps*(nb-Nvars/Ngroups);
        end
    end
    G=sparse(T,1:Nvars,1,Ngroups,Nvars); mb=(G*b)./max(eps,sum(G,2));
end

[nill,idx]=sort(sum(mb.^2,2),'descend');
G=full(sparse(KEEP,T,1,N,Ngroups));
G=G./repmat(max(eps,sum(G,1)),N,1);
G=G(:,idx);
end

function G=selectVars_v3(x,y,Nvars,Ngroups,varargin)
% keep "best" Nvars predictor variables, divides them into Ngroups groups, and sort them by largest-absolute-effect
Nvars=min(1e6,Nvars); % memory-cutoff
M=size(x,1);
N=size(x,2);
Ny=size(y,2);
b0=x'*y;
b=repmat(b0,1,M).*repelem(x',1,Ny);
[sb,idx]=sort(sum(abs(b0).^2,2),'descend');
KEEP=idx(1:Nvars);
b=b(KEEP,:);
mb=b(ceil(Nvars*rand(Ngroups,1)),:)+eps*randn(Ngroups,1);

for nrepeat=1:10
    %d=zeros(Nvars,Ngroups); for n=1:Ngroups, d(:,n)=sum(abs(b-mb(n+zeros(Nvars,1),:)).^2,2); end % note: use this syntax in older Matlab versions
    %mb=zeros(Ngroups,Ny); for n=1:Ngroups, mb(n,:)=mean(b(T==n,:),1); end
    d=zeros(Nvars,Ngroups); for n=1:Ngroups, d(:,n)=sum(abs(b-mb(n,:)).^2,2); end
    T=minbalanced(d);
    %[nill,T]=min(d,[],2);
    G=sparse(T,1:Nvars,1,Ngroups,Nvars); mb=(G*b)./max(eps,sum(G,2));
end

[nill,idx]=sort(sum(mb.^2,2),'descend');
G=full(sparse(KEEP,T,1,N,Ngroups));
G=G./repmat(max(eps,sum(G,1)),N,1);
G=G(:,idx);
end

function G=selectVars_v2(x,y,Nvars,Ngroups,varargin)
% keep "best" Nvars predictor variables, divides them into Ngroups groups, and sort them by largest-absolute-effect
Nvars=min(1e4,Nvars); % memory-cutoff
N=size(x,2);
b=x'*y;
[sb,idx]=sort(sum(abs(b).^2,2),'descend');
KEEP=idx(1:Nvars);

I=tril(ones(Nvars),-1);
Y0=permute(sqrt(sum(abs(conn_bsxfun(@minus,b(KEEP,:),permute(b(KEEP,:),[3,2,1]))).^2,2)),[1,3,2]);
Y0(eye(size(Y0))>0)=0;
Y0=(Y0+Y0')/2;
Y0b=Y0(I>0)';
Z = conn_statslinkage(Y0b, 'co');
T = conn_statscluster(Z, 'maxclust', Ngroups);
mb=accumarray(T,b(KEEP,:),[Ngroups,1],@mean).^2;
%nb=accumarray(T,1,[Ngroups,1]);
[nill,idx]=sort(mb,'descend');
G=full(sparse(KEEP,T,1,N,Ngroups));
G=G./repmat(max(eps,sum(G,1)),N,1);
G=G(:,idx);
end

function G=selectVars_v1(x,y,Nvars,Ngroups,varargin) 
% keep "best" Nvars predictor variables, divides them into Ngroups groups, and sort them by largest-absolute-effect
N=size(x,2);
b=x'*y;
[sb,idx]=sort(sum(abs(b).^2,2),'descend');
KEEP=idx(1:Nvars);
[sb,idx]=sort(b(KEEP,1));
G=full(sparse(KEEP(idx),ceil(Ngroups*(1:Nvars)/Nvars),1,N,Ngroups));
G=G./repmat(max(eps,sum(G,1)),N,1);
mb=sum(abs(G'*b).^2,2);
[nill,idx]=sort(mb,'descend');
G=G(:,idx);
end

function T=minbalanced(d)
Nvars=size(d,1);
Ngroups=size(d,2);
[mind,idx1]=sort(d,2);
[nill,idx2]=sort(mind);
idx1=idx1(idx2,:);
T=zeros(Nvars,1);
SizeGroups=zeros(1,Ngroups);
maxSizeGroups=ceil(Nvars/Ngroups);
for n=1:Nvars,
    i=idx1(n,find(~isnan(idx1(n,:)),1));
    T(idx2(n))=i;
    SizeGroups(i)=SizeGroups(i)+1;
    if SizeGroups(i)==maxSizeGroups
        idx1(idx1==i)=nan;
    end
end
end
