function varargout=conn_glmunivariate(option,X,Y,C)

switch(lower(option)),
	case {'estimate','estimatefixed'}
        if any(conn_server('util_isremotevar',{X,Y})), [varargout{1:nargout}]=conn_server('run_keep',mfilename,option,X,Y); return; end
		[N1,Nx]=size(X);
		[N2,Ns]=size(Y);
		if N1~=N2, error('wrong dimensions'); end
		
		Nx=rank(X);
		dof=N1-Nx;

        XX=X'*X;
		iX=pinv(XX);
		XY=X'*Y;
		B=iX*XY; 
		EE0=sum(abs(Y).^2,1);
		EE=EE0-sum(XY.*B,1);

		varargout{1}=B;
        if strcmp(lower(option),'estimate'),
            varargout{2}=struct('B',B,'X',X,'Y',Y,'EE',EE,'EE0',EE0,'iX',iX,'dof',dof,'XX',XX);
        else,
            varargout{2}=struct('B',B,'X',X,'Y',Y,'EE',EE,'EE0',EE0,'iX',iX,'dof',C.dof,'XX',XX,'SE',C.data);
        end
	case {'evaluate','evaluatefixed','evaluater'},
        if any(conn_server('util_isremotevar',{X,C})), [varargout{1:nargout}]=conn_server('run',mfilename,option,X,Y,C); return; end % note: expands output
        %if conn_server('util_isremotevar',X), X=conn_server('run',X); end
        %if conn_server('util_isremotevar',C), C=conn_server('run',C); end
		[N1,Nx]=size(X.X);
		Nc0=rank(X.X*C');
		h=C*X.B;
		r=C*X.iX*C';
		BB=sum(conj(h).*(pinv(r)*h),1); 	% diag(between matrix)
		
        if strcmpi(option,'evaluate')||strcmpi(option,'evaluater'),
            if 0,%size(h,1)>1,
                F=                real(BB./max(eps,X.EE))*X.dof/max(eps,Nc0);
                idxvalid=find(~isnan(F)&~isinf(F));if isempty(idxvalid),p=F;else,p=nan+zeros(size(F));p(idxvalid)=1-spm_Fcdf(F(idxvalid),max(eps,Nc0),max(eps,X.dof));end
                statsname='F';
            else,
                k=                sqrt(diag(r)*X.EE);
                F=			      real(h./max(eps,k))*sqrt(X.dof);
                p=nan+zeros(size(F));idxvalid=find(~isnan(F)&~isinf(F));if ~isempty(idxvalid), p(idxvalid)=1-spm_Tcdf(F(idxvalid),X.dof); end
                statsname='T';
            end
        else,
            c=(C*X.iX*X.X');
            X.dof=sum(X.dof);
            F=(h./sqrt((abs(c).^2)*(abs(X.SE).^2)));
            idxvalid=find(~isnan(F)&~isinf(F));if isempty(idxvalid),p=F;else,p=nan+zeros(size(F));p(idxvalid)=1-spm_Ncdf(F(idxvalid));end
            statsname='Z';
        end
        
		dof=              [Nc0,X.dof];
        if strcmpi(option,'evaluater')
            varargout={(real(BB./max(eps,X.EE0)))};
        else
            varargout{1}=real(h);
            varargout{2}=F;
            varargout{3}=p;
            varargout{4}=dof;
            varargout{5}=(real(BB./max(eps,X.EE0)));
            varargout{6}=statsname;
            if size(h,1)==1, varargout{5}=sign(real(h)).*varargout{5}; end
            %varargout{4}=(X.EE0-X.EE)./max(eps,X.EE0);
        end
end
	