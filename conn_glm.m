function [h,F,p,dof,statsname,B,BB,EE,opt]=conn_glm(X,Y,C,M,D,opt,dop)
% CONN_GLM General Linear Model estimation and hypothesis testing.
%
%   [h,F,p,dof]=CONN_GLM(X,Y,C,M,D) estimates a linear model of the form Y = X*B + E
%   where Y is an observed matrix of response or output variables (rows are observations, columns are output variables)
%         X is an observed design or regressor matrix (rows are observations, columns are predictor variables)
%         B is a matrix of unknown regression parameters (rows are predictor variables, columns are output variables)
%         E is an matrix of unobserved multivariate normally distributed disturbances with zero mean and unknown covariance.
%   and tests a general null hypothesis of the form C*B*M' = D
%   where C is matrix or vector of "predictor" contrasts (rows are contrasts, columns are predictor variables, defaults to C=eye(size(X,2)) )
% 		  M is matrix or vector of "outcome" contrasts (rows are contrasts, columns are output variables, defaults to M=eye(size(Y,2)) )
% 		  D is matrix or vector of "baseline" values (rows are predictor contrasts, columns are outcome contrasts, defaults to D=0 )
%
%   CONN_GLM returns the following information:
%		  h:   a matrix of estimated contrast effect sizes (C*B*M'-D)
%		  F:   the test statistic(s) (T,F,or Chi2 value, depending on whether h is a scalar, a vector, or a matrix. See below)
%		  p:   p-value of the test(s)
%		  dof: degrees of freedom
%
%   Additional information:
%   By default CONN_GLM will use a T, F, or a Chi2 statistic for hypothesis testing depending on the size of h=C*B*M'. The default options are:
%                  when size(h)=[1,1]      -> T statistic (note: one-sided t-test)
%	                      			          Examples of use: one-sided two-sample t-test, linear regression
%                  when size(h)=[1,Ns]     -> F statistic (note: equivalent to two-sided t-test when Ns=1)
%   							  	 		  Examples of use: Hotelling's two sample t-square test, repeated measures ANOVA, multivariate regression
%                  when size(h)=[Nc,1]     -> F statistic (note: equivalent to two-sided t-test when Nc=1)
%					  			     		  Examples of use: ANOVA, ANCOVA, linear regression omnibus test
%                  when size(h)=[Nc,Ns]    -> Wilks' Lambda statistic
%								     		  Examples of use: MANOVA, MANCOVA, multivariate regression omnibus test, likelihood ratio test
%   The default option can be changed using the syntax CONN_GLM(X,Y,C,M,opt) where opt is one of the following character strings:
%               CONN_GLM(X,Y,C,M,'collapse_none') will perform a separate univariate test (T-statistic) on each of the elements of the matrix h=C*B*M'
%  			    CONN_GLM(X,Y,C,M,'collapse_outcomes') will perform a separate multivariate test (F-statistic) on each of the rows of the matrix h (collapsing across multiple outcome variables or outcome contrasts)
% 			    CONN_GLM(X,Y,C,M,'collapse_predictors') will perform a separate multivariate test (F-statistic) on each of the columns of the matrix h (collapsing across multiple predictor variables or predictor contrasts)
% 			    CONN_GLM(X,Y,C,M,'collapse_all_rao') will perform a single omnibus multivariate test (Wilks Lambda statistic, Rao's F approximation; no assumptions on form of M*E'*E*M covariance) on the matrix h
% 			    CONN_GLM(X,Y,C,M,'collapse_all_bartlett') will perform a single omnibus multivariate test (Wilks Lambda statistic, Bartlett's Chi2 approximation; no assumptions on form of M*E'*E*M covariance) on the matrix h
%               CONN_GLM(X,Y,C,M,'collapse_all_satterthwaite') will perform a single omnibus univariate test (F-statistic with conservative Satterthwaite dof correction; no assumptions on form of M*E'*E*M covariance) on the matrix h
%               CONN_GLM(X,Y,C,M,'collapse_all_sphericity') will perform a single omnibus univariate test (F-statistic assumming sphericity M*E'*E*M' = sigma*I) on the matrix h
% 			    CONN_GLM(X,Y,C,M,'collapse_all') same as 'collapse_all_rao'
%
% Example of use:
%   % MANOVA (three groups, two outcome variables)
%   % Data preparation
%    N1=10;N2=20;N3=30;
%    Y1=randn(N1,2)+repmat([0,0],[N1,1]); % data for group 1 (N1 samples, population mean = [0,0])
%    Y2=randn(N2,2)+repmat([0,1],[N2,1]); % data for group 2 (N2 samples, population mean = [0,1])
%    Y3=randn(N3,2)+repmat([1,0],[N3,1]); % data for group 2 (N3 samples, population mean = [1,0])
%    Y=cat(1,Y1,Y2,Y3);
%    X=[ones(N1,1),zeros(N1,2); zeros(N2,1),ones(N2,1),zeros(N2,1); zeros(N3,2),ones(N3,1)];
%   % Sample data analyses
%    [h,F,p,dof]=conn_glm(X,Y,[1,-1,0;0,1,-1]); disp(['Multivariate omnibus test of non-equality of means across the three groups:']);disp([' F(',num2str(dof(1)),',',num2str(dof(2)),') = ',num2str(F),'   p = ',num2str(p)]);
%    [h,F,p,dof]=conn_glm(X,Y,[1,-1,0]); disp(['Multivariate test of non-equality of means between groups 1 and 2:']);disp([' F(',num2str(dof(1)),',',num2str(dof(2)),') = ',num2str(F),'   p = ',num2str(p)]);
%    [h,F,p,dof]=conn_glm(X,Y,[-1,1,0],eye(2),0,'collapse_none'); disp(['Univariate one-sided test of non-equality of means between groups 1 and 2 on each outcome variable:']);disp([' T(',num2str(dof),') = ',num2str(F(:)'),'   p = ',num2str(p(:)')]);
%

% alfnie@gmail.com
% 04/03

[N1,Nx]=size(X);
[N2,Ns,Na]=size(Y);
if N1~=N2, error('wrong dimensions'); end
if nargin<3 || isempty(C), C=eye(Nx); end
if nargin<4 || isempty(M), M=speye(Ns,Ns); else, Ns=rank(M); end
if nargin>=5&&nargin<=6&&(isempty(D)||ischar(D))       % old syntax (X,Y,C,M, opt,dop)
    if nargin<5, D=[]; end
    if nargin<6||isempty(opt), opt=true; end
    [D,opt,dop]=deal(0,D,opt);
else                                        % new syntax (X,Y,C,M, D,opt,dop)
    if nargin<5, D=0; end
    if nargin<6, opt=[]; end
    if nargin<7||isempty(dop), dop=true; end
end
if ~isempty(opt),switch(lower(opt)),case {'collapse_none','aa','t'},opt='AA';case {'collapse_outcome','collapse_outcomes','ab','f','frow'},opt='AB';case {'collapse_predictor','collapse_predictors','ba','fcol'},opt='BA';case {'collapse_all','bb'},opt='BB';case {'collapse_all_rao','bb_rao'},opt='BB_RAO';case {'chi2','collapse_all_bartlett','bb_bartlett'},opt='BB_BARTLETT';case {'collapse_all_satterthwaite','bb_satterthwaite'},opt='BB_SATTERTHWAITE';case {'collapse_all_sphericity','bb_sphericity'},opt='BB_SPHERICITY';otherwise,error(['Unknown option ',opt]); end; end
BBDEFAULT='BB_RAO';         % default collapse_all procedure is Rao (BB_RAO, BB_BARTLETT, BB_SATTERTHWAITE, BB_SPHERICITY)
FIXDFSATTERTHWAITE=true;    % Satterthwaite procedure returns modified F statistic & original dof values (set to false to return original F statistic & modified dof values)

Nx=rank(X);
dofe=N1-Nx;
Nc0=rank(X*C');

iX=pinv(X'*X);
r=C*iX*C';
ir=pinv(r);
if isempty(opt), opt=[char(65+(size(C,1)>1)),char(65+(size(M,1)>1))]; 
else, opt=upper(opt); 
end
if Na>1, Na_h=[]; nas=find(~any(any(isnan(Y)|isinf(Y),1),2)); else nas=1; end
if isempty(nas), nas=(1:Na)'; end
univariate=(strcmp(opt,'AA')||strcmp(opt,'BA'))&&isequal(M,speye(Ns,Ns));
if strcmp(opt,'BB'),opt=BBDEFAULT;end
iXX=iX*X';

if Ns==1&&isequal(M,1) %(Ns=1,M=1)
    B=iXX*Y(:,:);
    E=Y(:,:)-X*B;
    EE=sum(abs(E).^2,1);
    h=roundeps(C*B-D);
    switch(opt),
        case 'AA',                          % h: [1,1]  T-test
            k=                sqrt(diag(r)*EE);
            F=			      real(h./max(eps,k))*sqrt(dofe);
        case 'AB',                          % h: [1,1] F-test
            k=                diag(r)*EE;
            F=			      real((abs(h).^2)./max(eps,k))*(dofe-Ns+1)/Ns;
        otherwise,                          % h: [Nc,1] F-test
            BB=sum(h.*(ir*h),1);
            F=                real(BB./max(eps,EE))*dofe/Nc0;
    end
    h=permute(h,[1 3 2]);
    F=permute(F,[1 3 2]);
    B=permute(B,[1 3 2]);
    dofs=Ns;
    doft=1;
else
    if univariate
        B=iXX*Y(:,:);
        E=Y(:,:)-X*B;
        EE=sum(abs(E).^2,1);
        h=roundeps(C*B-D);
        switch(opt),
            case 'AA',                          % h: [1,1]  T-test
                k=                sqrt(diag(r)*EE);
                F=			      real(h./max(eps,k))*sqrt(dofe);
            case 'BA',                          % h: [Nc,1] F-test
                dBB=              sum(conj(h).*(ir*h),1);     % between matrix
                F=                real(dBB./max(eps,EE))*dofe/Nc0;
        end
        F=reshape(F, size(F,1),Ns,[]);
        h=reshape(h, size(h,1),Ns,[]);
        B=reshape(B, size(B,1),size(Y,2),[]);
    else
        for na=nas(:)'
            B=iXX*Y(:,:,na);
            E=Y(:,:,na)-X*B;
            if univariate,EE=sparse(1:Ns,1:Ns,sum(abs(E).^2,1)); % univariate case within matrix
            else, if size(E,2)<size(E,1), EE=M*(E'*E)*M'; else, EE=E*M'; EE=EE'*EE; end; end          	% within matrix
            %EE=full(EE);
            h=roundeps(full(C*B*M'-D));
            dofs=Ns;
            
            switch(opt),
                case 'AA',                          % h: [1,1]  T-test
                    k=                sqrt(diag(r)*full(diag(EE)).');
                    F=			      real(h./max(eps,k))*sqrt(dofe);
                case 'AB',                          % h: [1,Ns] F-test
                    if dofe<Ns
                        F=            nan(size(h,1),1); % note: could use BB_SATTERTHWAITE here
                    else
                        F=			  real(sum((h*pinv(EE)).*conj(h),2)./max(eps,diag(r)))*(dofe-Ns+1)/Ns;
                    end
                case 'BA',                          % h: [Nc,1] F-test
                    dBB=sum(conj(h).*(ir*h),1).';                    % between matrix
                    F=                real(dBB./max(eps,full(diag(EE))))*dofe/Nc0;
                case 'BB_BARTLETT',                          % h: [Nc,Ns] Wilks Lambda-test (Ns,dof,Nc0)
                    if Ns==1 % rank-defficient EE case
                        dBB=sum(conj(h).*(ir*h),1).';                    % between matrix
                        F=            real(dBB./max(eps,full(diag(EE))))*dofe/Nc0; F=F(1);
                    elseif dofe<Ns
                        F=            nan;              % note: could use BB_SATTERTHWAITE here
                    else
                        BB=h'*ir*h;                    % between matrix
                        F=            -(dofe-1/2*(Ns-Nc0+1))*real(log(real(roundeps(det(EE/dofe))./roundeps(det((EE+BB)/dofe)))));
                        %F=    (dofe-1/2*(Ns-Nc0+1))*real(log(real(roundeps(det(eye(size(EE,2))+pinv(EE)*BB)))));  
                        %F=    (dofe-1/2*(Ns-Nc0+1))*real(log(real(roundeps(det(eye(size(EE,1))+EE\BB)))));
                        %F=    -(dofe-1/2*(Ns-Nc0+1))*real(sum(log(real(eig(full(EE)))))-sum(log(real(eig(full(EE+BB))))));
                    end
                case {'BB','BB_RAO'},          % h: [Nc,Ns] Wilks Lambda-test (Ns,dof,Nc0)
                    BB=h'*ir*h;                    % between matrix
                    l1=real(det((EE+BB)/dofe));
                    if l1>0, l=     max(0,real((det(EE/dofe))./l1));
                    else     l=     1;
                    end
                    if Ns^2+Nc0^2-5<=0, doft=1;
                    else doft=sqrt((Ns^2*Nc0^2-4)/(Ns^2+Nc0^2-5));
                    end
                    if dofe<Ns
                        F=            nan(size(l));     % note: could use BB_SATTERTHWAITE here
                    else
                        F=            ((dofe-1/2*(Ns-Nc0+1))*doft-(Ns*Nc0-2)/2)/(Ns*Nc0)*(1-l.^(1/doft))./max(eps,l.^(1/doft));
                    end
                case 'BB_SATTERTHWAITE',       % h: [Nc,Ns] F-test (Satterthwaite correction)
                    dofs=trace(EE)^2/sum(sum(EE.^2,1),2);  % tr(EE)^2/tr(EE*EE)
                    dBB=sum(conj(h).*(ir*h),1).';
                    F=                real(sum(dBB)./max(eps,full(sum(diag(EE)))))*dofe/Nc0;
                    if FIXDFSATTERTHWAITE&~(isnan(F)|F==0), F=spm_invFcdf(spm_Fcdf(F,dofs*Nc0,dofs*dofe),Nc0,dofe); end
                case 'BB_SPHERICITY',          % h: [Nc,Ns] F-test (sphericity assumption)
                    dBB=sum(conj(h).*(ir*h),1).';
                    F=                real(sum(dBB)./max(eps,full(sum(diag(EE)))))*dofe/Nc0;
            end
            if Na>1
                if isempty(Na_h), Na_h=nan(size(h,1),size(h,2),Na);Na_F=nan(size(F,1),size(F,2),Na); Na_B=nan(size(B,1),size(B,2),Na);Na_dofs=nan(size(dofs,1),size(dofs,2),Na);end
                Na_h(:,:,na)=h;
                Na_F(:,:,na)=F;
                Na_B(:,:,na)=B;
                Na_dofs(:,:,na)=dofs;
            end
        end
        if Na>1
            h=Na_h;
            F=Na_F;
            B=Na_B;
            dofs=Na_dofs;
        end
    end
end
if nargout>2
    p=nan(size(F));
    idxvalid=find(~(isnan(F)|F==0));
    switch(opt),
        case 'AA',                          % h: [1,1]  T-test
            dof=              dofe;
            statsname=        'T';
        case 'AB',                          % h: [1,Ns] F-test
            dof=              [Ns,dofe-Ns+1];
            statsname=        'F';
        case 'BA',                          % h: [Nc,1] F-test
            dof=              [Nc0,dofe];
            statsname=        'F';
        case 'BB_BARTLETT',                          % h: [Nc,Ns] Wilks Lambda-test (Ns,dof,Nc0)
            if Ns==1 % rank-defficient EE case
                dof=          [Nc0,dofe];
                statsname=    'F';
            else
                dof=          [Ns*Nc0];
                statsname=    'X';
            end
        case {'BB','BB_RAO'}
            dof=              [Ns*Nc0 (dofe-1/2*(Ns-Nc0+1))*doft-(Ns*Nc0-2)/2];
            statsname=        'F';
        case 'BB_SATTERTHWAITE',       % h: [Nc,Ns] F-test (Satterthwaite correction)
            if FIXDFSATTERTHWAITE, dof=[Nc0,dofe];
            else dof=         cat(2,dofs*Nc0,dofs*dofe);
            end
            statsname=        'F';
        case 'BB_SPHERICITY',          % h: [Nc,Ns] F-test (sphericity assumption)
            dof=              [Ns*Nc0,Ns*dofe];
            statsname=        'F';
    end
    if ~isempty(idxvalid)&&dop,
        switch(opt),
            case 'AA',                          % h: [1,1]  T-test
                p(idxvalid)=spm_Tcdf(-F(idxvalid),dofe);
            case 'AB',                          % h: [1,Ns] F-test
                p(idxvalid)=1-spm_Fcdf(F(idxvalid),Ns,dofe-Ns+1);
            case 'BA',                          % h: [Nc,1] F-test
                p(idxvalid)=1-spm_Fcdf(F(idxvalid),Nc0,dofe);
            case 'BB_BARTLETT',                 % h: [Nc,Ns] Wilks Lambda-test (Ns,dof,Nc0)
                if Ns==1 % rank-defficient EE case
                    p(idxvalid)=1-spm_Fcdf(F(idxvalid),Nc0,dofe);
                else
                    p(idxvalid)=1-spm_Xcdf(F(idxvalid),Ns*Nc0);
                end
            case {'BB','BB_RAO'}
                p(idxvalid)=1-spm_Fcdf(F(idxvalid),Ns*Nc0,(dofe-1/2*(Ns-Nc0+1))*doft-(Ns*Nc0-2)/2);
            case {'BB_SATTERTHWAITE','BB_SPHERICITY'}, % h: [Nc,Ns] F-test (Satterthwaite or sphericity assumption)
                if FIXDFSATTERTHWAITE, 
                    p(idxvalid)=1-spm_Fcdf(F(idxvalid),Nc0,dofe);
                else
                    p(idxvalid)=1-spm_Fcdf(F(idxvalid),dofs(idxvalid)*Nc0,dofs(idxvalid)*dofe);
                end
        end
    end
end
end



function y=roundeps(y)
y(abs(y)<eps)=0;
end
