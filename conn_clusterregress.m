function varargout=conn_clusterregress(X,Y,varargin)
% Ensemble Clustered Regression (ECR)
% Y(X) = A + X*B + noise
% for robust regression when size(X,2)>>size(X,1)
% with clustering over columns of B, split-half 
% resampling for variable selection and fitting,
% and bagging across multiple resampling splits.
%
% [Model,Yfit]=conn_clusterregress(X,Y);
% Y                       : outcome variables [# samples x # outcomes]
% X                       : predictor variables [# samples x # of predictor variables]
% Model.parameters.A      : model bias term [1 x # outcomes] 
% Model.parameters.B      : model regressor coefficients [# of predictor variables x # outcomes] 
% Model.parameters.B_std  : model regressor standardized coefficients (scaled to unit standard deviation predictor and outcome measures) [# of predictor variables x # outcomes] 
% Model.fit.MSE           : model out-of-sample Mean Square Error
% Model.fit.MSE_std       : model out-of-sample standardized Mean Square Error (scaled to unit standard deviation outcome measures)
% Model.predict           : predictor function for new samples | predict(X) = A+B*X for linear regression model and 1./(1+exp(-A-B*X)) for a logistic regression model, with model parameters A&B, which approximates expected outcome data Y for a new sample X
% Yfit                    : outcome variable model fits for training data [# samples x # outcomes]
%
% conn_clusterregress(X,Y,'paramName',paramValue,...);
% Use parameter name/value pairs to specify additional input arguments. Valid parameter names are:
%       covariates          : covariates matrix [# samples x # covariate variables] Z in a new model of the form: Y(X) = A + X*B + Z*C + noise (note: covariates are always included in all models; they are not part of feature selection)
%       display             : display info while fitting (default false)
%       warning             : display warnings while fitting (default true)
%       computeyfit         : computes Yfit values on training data using nested-crossvalidation (default true)
%                               set to false to skip the outer cross-validation loop needed to estimate Yfit (e.g. if model fits to the training dataset are unnecessary because we are embedding conn_clusterregress on our own nested crossvalidation loop)
%                               This will also output a faster but less accurate estimate of out-of-sample MSE using partial estimates computed during training before bagging (in Model.error.MSE)
%       probconvert         : converts outcome values to probabilities when appropriate (default true); when Y data is binary (0/1 values only) conn_clusterregress converts the clustered regression model outputs to values between 0 and 1 by performing a post-hoc logistic regression model between the fitted continuous values Yfit and the training binary values Y (default true) note: logreg model fits Y(X) = 1./(1+exp(-a-b*Yfit(X))) model estimating a&b values that minimize a negative log-likelihood cost function 
%       pvariables          : proportion of predictor variables to consider for clustering (default 1)
%       ngroups             : number of clusters of predictor variables (default: min( size(X,1)/10, size(X,2)*pvariables )
%       nbagging            : number of base learners (default size(X,1); note: nbagging is always rounded to an integer multiple of size(X,1))
%       nouterloop          : clustering algorithm outer loop number of iterations (adapting cluster means) (default 10)
%       ninnerloop          : clustering algorithm inner loop number of iterations (adapting cluster sizes) (default 2)
%
% example:
%  N=100;  % samples
%  R=1000; % predictors
%  P=0.1;  % percentage true effects
%  y=randn(10*N,1); % training+testing outcome
%  x=[.25*y+randn(10*N,ceil(R*P)),randn(10*N,R-ceil(R*P))]; % training+testing predictors
%  [model,yfit_training]=conn_clusterregress(x(1:N,:),y(1:N,:)); % training set
%  yfit_testing=model.parameters.A+x(N+1:end,:)*model.parameters.B;
%  subplot(121); plot(model.parameters.B,'.'); xlabel 'Predictors'; ylabel 'Regressor Coefficients'; title(sprintf('example with %d samples and %d predictors',N,R)); 
%  subplot(222); plot(y(1:N),yfit_training,'ko','markerfacecolor','k'); grid on; xlabel target; ylabel prediction; title 'training set'
%  subplot(224); plot(y(N+1:end),yfit_testing,'ko','markerfacecolor','k'); grid on; xlabel target; ylabel prediction; title 'testing set'
%  mse=mean(abs(y(N+1:end)-yfit_testing).^2); % testing set
%  disp(mse);

options=struct(...
    'pvariables',[],...       % proportion of predictor variables to consider for clustering (default 100%)
    'ngroups',[], ...         % maximum number of clusters among kept predictor variables (default: min( size(X,1)/10, size(X,2)*pvariables )
    'display',false,...       % display info while fitting
    'warning',false,...       % display warnings while fitting
    'computeyfit',true,...    % compute fit values for training data using nested cross-validation
    'probconvert',true,...     % fits binary output data using logreg post-hoc fit
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
if options.probconvert&&~all(Y(:)==0|Y(:)==1), options.probconvert=false; end
if isempty(options.pvariables), options.pvariables=1; end 
if isempty(options.nbagging), options.nbagging=N; end
options.nbagging=N*max(1,round(options.nbagging/N));
if isempty(options.ngroups), options.ngroups=ceil(N/10); end % minimum 10:1 samples:predictor ratio
options.ngroups=min(options.ngroups, ceil(Nx*options.pvariables));
Yfit=[];

if options.computeyfit||options.probconvert % use CV to estimate Yfit on training data
    K=5; % k-fold cv
    nsamples=ceil(N/K); % uses 80% of the samples for training; 20% of the samples for testing (5-fold crossvalidation)
    state=rng; rng(0,'twister');
    Yfit=zeros([N,Ny]);
    Nfit=zeros(N,1);
    for n1=1:4 % total 4*K folds (balanced sampling, each sample included 4 times in a testing dataset)
        samples=randperm(N); 
        for n2=0:K-1
            test=samples(nsamples*n2+1:min(N,nsamples*(n2+1)));
            train=setdiff(1:N,test);
            if isempty(Z), 
                Model=conn_clusterregress(X(train,:),Y(train,:),varargin{:},'computeyfit',false,'probconvert',false);
                Yfit(test,:)=Yfit(test,:) + Model.parameters.A + X(test,:)*Model.parameters.B;
            else 
                Model=conn_clusterregress(X(train,:),Y(train,:),varargin{:},'covariates',Z(train,:),'computeyfit',false,'probconvert',false);
                Yfit(test,:)=Yfit(test,:) + Model.parameters.A + X(test,:)*Model.parameters.B + Z(test,:)*Model.parameters.C;
            end
            Nfit(test)=Nfit(test)+1;
        end
    end
    Yfit=Yfit./repmat(Nfit,1,Ny);
    clear Model;
    Model.fit.MSE=mean(mean(abs(Yfit-Y).^2,1),2);
    Model.fit.MSE_std=mean(mean(abs(Yfit-Y).^2,1)./var(Y,1,1),2);
    if options.probconvert
        if any(~validS), Model.fit.YfitOriginal=nan(numel(validS),Ny); Model.fit.YfitOriginal(validS,:)=Yfit; 
        else Model.fit.YfitOriginal=Yfit; % keeps original fit values here (continuous linear regression model)
        end
        LOGREG_SCALE1=zeros(1,Ny);
        LOGREG_SCALE2=zeros(1,Ny);
        NLL=0; 
        for n1=1:Ny
            [t,yfit,nll]=conn_logreg([ones(N,1) Yfit(:,n1)],Y(:,n1),'regularization',1/N); % post-hoc logreg (note: post-hoc analysis performed outside of CV)
            Yfit(:,n1)=yfit;
            LOGREG_SCALE1(n1)=t(1);
            LOGREG_SCALE2(n1)=t(2);
            NLL=NLL+nll;
        end
        Model.fit.NegLogLikelihood=NLL/Ny; % Negative log-likelihood averaged across outcome measures
    end
    if any(~validS), temp=nan(numel(validS),Ny); temp(validS,:)=Yfit; Yfit=temp; end
    Model.fit.Yfit=Yfit;
    rng(state);
end

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

% divides (Nx*options.pvariables) predictor variables in ngroups groups (ngroups from 1 to options.ngroups); keep first n of these as predictors (n from 1 to ngroups)
%samples=1:N; % placeholder for nested CV
Ngroups=options.ngroups; while Ngroups(end)>1, Ngroups(end+1)=round(Ngroups(end)/2); end
ntests=sum(Ngroups);
B=zeros([Nx,Ny,ntests]);
err=zeros(1,ntests);
nparameters=zeros(1,ntests);
errnull=0;
[nill,samples]=sort(sum(Y.^2,2),'descend'); % sorts samples by normY to balance TRAIN1/TRAIN2 sets
state=rng; rng(0,'twister');
for nrepeat=1:options.nbagging
    if options.display, fprintf('.'); end
    idx=(1:2:N)+(rand(1,ceil(N/2))>.5);
    idx(idx>N)=[];
    ns=numel(idx);
    TRAIN1=samples(idx);
    TRAIN2=samples; TRAIN2(idx)=[];
    %idx=randperm(N);
    %TRAIN1=samples(idx(1:ns));      % use these samples for variable selection
    %TRAIN2=samples(idx(ns+1:N));    % use these samples for parameter estimates
    syy=sum(sum(Y(TRAIN2,:).^2,1),2);
    ntest=0;
    for ngroups=Ngroups % divide predictors in n groups
        if options.pvariables==1&&ngroups==1, G=ones(Nx,1)/Nx; 
        else G=selectVars(X(TRAIN1,:),Y(TRAIN1,:),ceil(Nx*options.pvariables),ngroups,options.nouterloop,options.ninnerloop); % [Nx x ngroups] 0/1 group matrix (note: sum-normed columns)
        end
        xg=X(TRAIN2,:)*G;
        for nkeep=1:ngroups % keep first m groups and use those as separate predictors
            tG=G(:,1:nkeep);
            txg=xg(:,1:nkeep);
            txgy=txg'*Y(TRAIN2,:);
            b=pinv(txg'*txg)*txgy;
            terr=(syy-sum(sum(b.*txgy,1),2))/(N-ns); % sum of MSE over outcome variables: sum(mean(abs(Y(TRAIN2,:)-X(TRAIN2,:)*Gb).^2,1),2)
            Gb=tG*b;
            ntest=ntest+1;
            B(:,:,ntest)=B(:,:,ntest)+Gb;
            err(ntest)=err(ntest)+terr;
            nparameters(ntest)=nkeep;
        end
    end
    errnull=errnull+syy/(N-ns);
end
rng(state);
err=err/options.nbagging/Ny;
errnull=errnull/options.nbagging/Ny;
if options.display, fprintf('mse = %s : # parameters %s \n',mat2str(err),mat2str(nparameters)); end
aic=log(err)*N/2+2*nparameters;
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
if options.probconvert
    Model.parameters.A=Model.parameters.A.*LOGREG_SCALE2+LOGREG_SCALE1;
    Model.parameters.B=Model.parameters.B.*LOGREG_SCALE2;
    if isempty(Z), 
        Model.predict=@(x,varargin)1./(1+exp(-Model.parameters.A-x*Model.parameters.B));
    else 
        Model.parameters.C=Model.parameters.C.*LOGREG_SCALE2;
        Model.predict=@(x,z,varargin)1./(1+exp(-Model.parameters.A-x*Model.parameters.B-z*Model.parameters.C));
    end
else
    if isempty(Z), Model.predict=@(x,varargin)Model.parameters.A+x*Model.parameters.B;
    else Model.predict=@(x,z,varargin)Model.parameters.A+x*Model.parameters.B+z*Model.parameters.C;
    end
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
ab=zeros(Ngroups,1);              % size-bias in each cluster (keeps clusters balanced)
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

