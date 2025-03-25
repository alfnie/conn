
function varargout=conn_designmatrix(confounds,X1,X2,nconfounds,selectonly)
if nargin<5||isempty(selectonly), selectonly=false; end
if nargin<4||isempty(nconfounds), nconfounds={}; end
if any(conn_server('util_isremotevar',{X1,X2})), [varargout{1:nargout}]=conn_server('run',mfilename,confounds,conn_server('util_cleanremotevar',X1),conn_server('util_cleanremotevar',X2),nconfounds,selectonly); return; end % note: returns expanded variables

select=[];X=[];names={};xyz={}; Xnames={};isfixed=[];
if iscell(confounds),
    temp=confounds{1};tempfields=fieldnames(temp);for n1=2:numel(confounds),for n2=1:numel(tempfields), temp.(tempfields{n2})=cat(2,temp.(tempfields{n2}),confounds{n1}.(tempfields{n2})); end; end; confounds=temp; 
    valid=ones(numel(confounds.names),1);
    for n1=1:numel(confounds.names),
        if valid(n1)
            %idx=intersect(strmatch(confounds.names{n1},confounds.names,'exact'),strmatch(confounds.types{n1},confounds.types,'exact'));
            idx=intersect(find(strcmp(confounds.names{n1},confounds.names)),find(strcmp(confounds.types{n1},confounds.types)));
            confounds.dimensions{n1}=max(cat(1,confounds.dimensions{idx}),[],1);
            confounds.deriv{n1}=max(cat(1,confounds.deriv{idx}),[],1);
            if isfield(confounds,'fbands')
                confounds.fbands{n1}=max(cat(1,confounds.fbands{idx}),[],1);
            end
            if isfield(confounds,'power')
                confounds.power{n1}=max(cat(1,confounds.power{idx}),[],1);
            end
            valid(idx(2:end))=0;
        end
    end
    if any(~valid)
        idx=find(valid);
        for n2=1:numel(tempfields), confounds.(tempfields{n2})=confounds.(tempfields{n2})(idx); end
    end
end
if nargin>3&&~isempty(nconfounds), select=cell(size(nconfounds)); end


for n1=1:length(confounds.names), 
    x=[];
	switch(confounds.types{n1}),
		case 'roi',
			%idx=strmatch(confounds.names{n1},X1.names,'exact');
			idx=find(strcmp(confounds.names{n1},X1.names));
			if isempty(idx)||length(idx)>1, error('Mismatch info: Please re-run Setup step'); end
            if isfield(X1,'fbdata')&&isfield(confounds,'fbands')&&numel(confounds.fbands)>=n1&&confounds.fbands{n1}>1
                x=X1.fbdata{idx}{confounds.fbands{n1}}(:,1:min(confounds.dimensions{n1}(1),size(X1.fbdata{idx}{confounds.fbands{n1}},2)),:);
                x=cat(3,X1.data{idx}(:,1:min(confounds.dimensions{n1}(1),size(X1.data{idx},2))),x);
                sx=[size(x,2),size(x,3)];
                x=x(:,:);
            else
                x=X1.data{idx}(:,1:min(confounds.dimensions{n1}(1),size(X1.data{idx},2)));
                sx=[size(x,2),size(x,3)];
            end
            if isfield(X1,'d1data'), d1x=X1.d1data{idx}(:,1:min(confounds.dimensions{n1}(1),size(X1.d1data{idx},2)));
            elseif ~isempty(x), d1x=convn(cat(1,x(1,:),x,x(end,:)),[1;0;-1],'valid'); 
            else d1x=x; end
            if isfield(X1,'d2data'), d2x=X1.d2data{idx}(:,1:min(confounds.dimensions{n1}(1),size(X1.d2data{idx},2)));
            elseif ~isempty(d1x), d2x=convn(cat(1,d1x(1,:),d1x,d1x(end,:)),[1;0;-1],'valid'); 
            else d2x=d1x; end
            dx={d1x,d2x};
		case 'cov',
			%idx=strmatch(confounds.names{n1},X2.names,'exact');
			idx=find(strcmp(confounds.names{n1},X2.names));
			if isempty(idx)||length(idx)>1, error('Mismatch info: Please re-run Setup step'); end
            if isfield(X2,'fbdata')&&isfield(confounds,'fbands')&&numel(confounds.fbands)>=n1&&confounds.fbands{n1}>1
                x=X2.fbdata{idx}{confounds.fbands{n1}}(:,1:min(confounds.dimensions{n1}(1),size(X2.fbdata{idx}{confounds.fbands{n1}},2)),:);
                sx=[size(x,2),size(x,3)];
                x=x(:,:);
            else
                x=X2.data{idx}(:,1:min(confounds.dimensions{n1}(1),size(X2.data{idx},2)));
                sx=[size(x,2),size(x,3)];
            end
            if isfield(X2,'d1data'), d1x=X2.d1data{idx}(:,1:min(confounds.dimensions{n1}(1),size(X2.d1data{idx},2)));
            elseif ~isempty(x), d1x=convn(cat(1,x(1,:),x,x(end,:)),[1;0;-1],'valid'); 
            else d1x=x; end
            if isfield(X2,'d2data'), d2x=X2.d2data{idx}(:,1:min(confounds.dimensions{n1}(1),size(X2.d2data{idx},2)));
            elseif ~isempty(d1x), d2x=convn(cat(1,d1x(1,:),d1x,d1x(end,:)),[1;0;-1],'valid'); 
            else d2x=d1x; end
            dx={d1x,d2x};
    end
    if ~isfield(confounds,'fixed')||numel(confounds.fixed)<n1, confounds.fixed{n1}=0; end
    if ~isempty(x),
        X=cat(2,X,x);
        if all(sx==1),          names{end+1}=[confounds.names{n1}];
        elseif sx(1)==1,        names{end+1}=[confounds.names{n1}]; for n0=1:size(x,2)-1,names{end+1}=[confounds.names{n1},'_Freq',num2str(n0)];end
        elseif sx(2)==1,        for n0=1:size(x,2),names{end+1}=[confounds.names{n1},'_Dim',num2str(n0)];end
        else                    for n0=1:sx(1),names{end+1}=[confounds.names{n1},'_Dim',num2str(n0)];end; for n02=1:sx(2)-1,for n01=1:sx(1), names{end+1}=[confounds.names{n1},'_Dim',num2str(n01),'_Freq',num2str(n02)];end; end
        end
        %for n0=1:size(x,2),names{end+1}=[confounds.names{n1},'_',num2str(1),'_',num2str(n0)];end
        if strcmp(confounds.types{n1},'roi')&&isfield(X1,'xyz'), for n0=1:size(x,2),xyz{end+1}=X1.xyz{idx}+(n0-1)*(confounds.deriv{n1}+1);end; 
        else for n0=1:size(x,2),xyz{end+1}=nan(1,3);end; 
        end
        if nargin>3&&~isempty(nconfounds), for n0=1:length(nconfounds),
                if any(n1==nconfounds{n0}), select{n0}=cat(2,select{n0},ones(1,size(x,2)));
                else select{n0}=cat(2,select{n0},zeros(1,size(x,2))); end
            end; end
        isfixed=[isfixed, repmat(confounds.fixed{n1},1,size(x,2))];
        if isfield(confounds,'power')
            for n2=2:confounds.power{n1},
                tx=x(:,:);
                tx=conn_bsxfun(@minus,tx,mean(tx,1));
                tx=conn_bsxfun(@rdivide,tx,max(eps,max(abs(tx),[],1)));
                tx=tx.^n2;
                %if ~isempty(X), tx=tx-X*(X\tx); end
                X=cat(2,X,tx);
                if all(sx==1),          names{end+1}=[confounds.names{n1},'_Pow',num2str(n2)];
                elseif sx(1)==1,        names{end+1}=[confounds.names{n1},'_Pow',num2str(n2)]; for n0=1:size(x,2)-1,names{end+1}=[confounds.names{n1},'_Pow',num2str(n2),'_Freq',num2str(n0)];end
                elseif sx(2)==1,        for n0=1:size(x,2),names{end+1}=[confounds.names{n1},'_Pow',num2str(n2),'_Dim',num2str(n0)];end
                else                    for n0=1:sx(1),names{end+1}=[confounds.names{n1},'_Pow',num2str(n2),'_Dim',num2str(n0)];end; for n02=1:sx(2)-1,for n01=1:sx(1), names{end+1}=[confounds.names{n1},'_Pow',num2str(n2),'_Dim',num2str(n01),'_Freq',num2str(n02)];end; end
                end
                %for n0=1:size(x,2),names{end+1}=[confounds.names{n1},'_',num2str(n2+1),'_',num2str(n0)];end
                if strcmp(confounds.types{n1},'roi')&&isfield(X1,'xyz'),for n0=1:size(x,2),xyz{end+1}=X1.xyz{idx}+(n0-1)*(confounds.power{n1})+n2;end;
                else for n0=1:size(x,2),xyz{end+1}=nan(1,3);end;
                end
                if nargin>3&&~isempty(nconfounds), for n0=1:length(nconfounds),
                        if any(n1==nconfounds{n0}), select{n0}=cat(2,select{n0},ones(1,size(x,2)));
                        else select{n0}=cat(2,select{n0},zeros(1,size(x,2))); end
                    end; end
                isfixed=[isfixed, repmat(confounds.fixed{n1},1,size(tx,2))];
            end
        end        
        for n2=1:min(confounds.deriv{n1},numel(dx)),
            %x=convn(cat(1,x(1,:),x,x(end,:)),[1;0;-1],'valid');
            x=dx{n2};
            sx=[size(x,2),size(x,3)];
            X=cat(2,X,x(:,:));
            if all(sx==1),          names{end+1}=[confounds.names{n1},'_Der',num2str(n2)];
            elseif sx(1)==1,        names{end+1}=[confounds.names{n1},'_Der',num2str(n2)]; for n0=1:size(x,2)-1,names{end+1}=[confounds.names{n1},'_Der',num2str(n2),'_Freq',num2str(n0)];end
            elseif sx(2)==1,        for n0=1:size(x,2),names{end+1}=[confounds.names{n1},'_Der',num2str(n2),'_Dim',num2str(n0)];end
            else                    for n0=1:sx(1),names{end+1}=[confounds.names{n1},'_Der',num2str(n2),'_Dim',num2str(n0)];end; for n02=1:sx(2)-1,for n01=1:sx(1), names{end+1}=[confounds.names{n1},'_Der',num2str(n2),'_Dim',num2str(n01),'_Freq',num2str(n02)];end; end
            end
            %for n0=1:size(x,2),names{end+1}=[confounds.names{n1},'_',num2str(n2+1),'_',num2str(n0)];end
            if strcmp(confounds.types{n1},'roi')&&isfield(X1,'xyz'),for n0=1:size(x,2),xyz{end+1}=X1.xyz{idx}+(n0-1)*(confounds.deriv{n1}+1)+n2;end;
            else for n0=1:size(x,2),xyz{end+1}=nan(1,3);end; 
            end
            if nargin>3&&~isempty(nconfounds), for n0=1:length(nconfounds),
                    if any(n1==nconfounds{n0}), select{n0}=cat(2,select{n0},ones(1,size(x,2)));
                    else select{n0}=cat(2,select{n0},zeros(1,size(x,2))); end
                end; end
            isfixed=[isfixed, repmat(confounds.fixed{n1},1,size(x(:,:),2))];
        end
        if isfield(confounds,'power')&&confounds.deriv{n1}>0
            for n2a=1:min(confounds.deriv{n1},numel(dx)),
                x=dx{n2a};
                sx=[size(x,2),size(x,3)];
                for n2=2:confounds.power{n1},
                    tx=x(:,:);
                    tx=conn_bsxfun(@minus,tx,mean(tx,1));
                    tx=conn_bsxfun(@rdivide,tx,max(eps,max(abs(tx),[],1)));
                    tx=tx.^n2;
                    %if ~isempty(X), tx=tx-X*(X\tx); end
                    X=cat(2,X,tx);
                    if all(sx==1),          names{end+1}=[confounds.names{n1},'_Der',num2str(n2a),'_Pow',num2str(n2)];
                    elseif sx(1)==1,        names{end+1}=[confounds.names{n1},'_Der',num2str(n2a),'_Pow',num2str(n2)]; for n0=1:size(x,2)-1,names{end+1}=[confounds.names{n1},'_Der',num2str(n2a),'_Pow',num2str(n2),'_Freq',num2str(n0)];end
                    elseif sx(2)==1,        for n0=1:size(x,2),names{end+1}=[confounds.names{n1},'_Der',num2str(n2a),'_Pow',num2str(n2),'_Dim',num2str(n0)];end
                    else                    for n0=1:sx(1),names{end+1}=[confounds.names{n1},'_Der',num2str(n2a),'_Pow',num2str(n2),'_Dim',num2str(n0)];end; for n02=1:sx(2)-1,for n01=1:sx(1), names{end+1}=[confounds.names{n1},'_Der',num2str(n2a),'_Pow',num2str(n2),'_Dim',num2str(n01),'_Freq',num2str(n02)];end; end
                    end
                    %for n0=1:size(x,2),names{end+1}=[confounds.names{n1},'_',num2str(n2+1),'_',num2str(n0)];end
                    if strcmp(confounds.types{n1},'roi')&&isfield(X1,'xyz'),for n0=1:size(x,2),xyz{end+1}=X1.xyz{idx}+(n0-1)*(confounds.power{n1})+n2;end;
                    else for n0=1:size(x,2),xyz{end+1}=nan(1,3);end;
                    end
                    if nargin>3&&~isempty(nconfounds), for n0=1:length(nconfounds),
                            if any(n1==nconfounds{n0}), select{n0}=cat(2,select{n0},ones(1,size(x,2)));
                            else select{n0}=cat(2,select{n0},zeros(1,size(x,2))); end
                        end; end
                    isfixed=[isfixed, repmat(confounds.fixed{n1},1,size(tx,2))];
                end
            end
        end
    end
end
Xnames=names;
N=size(X,1); if ~N, N=size(X1.data{1},1); end
if any(strcmp(confounds.types,'detrend'))
    X=[X,linspace(-1,1,N)'];
    Xnames{end+1}='trend linear term';
    if nargin>3&&~isempty(nconfounds), for n0=1:length(nconfounds),select{n0}=[select{n0},0];end; end
    isfixed(end+1)=0;
end
if any(strcmp(confounds.types,'detrend2'))
    X=[X,linspace(-1,1,N)'.^2];
    Xnames{end+1}='trend quadratic term';
    if nargin>3&&~isempty(nconfounds), for n0=1:length(nconfounds),select{n0}=[select{n0},0];end; end
    isfixed(end+1)=0;
end
if any(strcmp(confounds.types,'detrend3'))
    X=[X,linspace(-1,1,N)'.^3];
    Xnames{end+1}='trend cubic term';
    if nargin>3&&~isempty(nconfounds), for n0=1:length(nconfounds),select{n0}=[select{n0},0];end; end
    isfixed(end+1)=0;
end

%Xnames(isfixed>0)=regexprep(Xnames(isfixed>0),'.*','$0_fixed');
%Xnames(isfixed==0)=regexprep(Xnames(isfixed==0),'(_fixed)+$','');

X=[ones(N,1),X];
Xnames=[{'constant term'}, Xnames];
isfixed=[false,isfixed>0];
if nargin>3&&~isempty(nconfounds), for n0=1:length(nconfounds),select{n0}=[0,select{n0}];end; end

%X=[X(:,~isfixed), X(:,isfixed)];
%Xnames=[Xnames(~isfixed), Xnames(isfixed)];
%if nargin>3&&~isempty(nconfounds), for n0=1:length(nconfounds),select{n0}=[select{n0}(~isfixed) select{n0}(isfixed)];end; end

if selectonly, varargout={select};
else varargout={X,select,names,xyz,Xnames,isfixed};
end
end

