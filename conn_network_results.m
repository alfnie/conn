function ss=conn_network_results(filename,ask)
%

% define design matrix
ss=[];
posstr=1;

if nargin<1||isempty(filename),
    global CONN_x
    if ~isempty(CONN_x)&&isfield(CONN_x,'folders')&&isfield(CONN_x.folders,'firstlevel'),pathname=fullfile(CONN_x.folders.firstlevel,CONN_x.Analyses(1).name);
    else pathname=pwd; end
    [filename, pathname] = uigetfile('.mat', 'Select a network analysis file',fullfile(pathname,'*.network'));
    filename=fullfile(pathname,filename);
end
if iscell(filename), results=filename{1}; ss=filename{2}; filename=[];
else load(filename,'-mat'); end
%if ~isempty(ss),ss.askn=2;ss.ask='all';end
if nargin>1,ss.ask=ask;end
if isfield(results(1),'measure_names'),ss.n=numel(results(1).measures.(results(1).measure_names{1}));
elseif isfield(results(1),'pathlength'),ss.n=numel(results(1).pathlength);
else ss.n=numel(results(1).GlobalEfficiency); 
end
%results=struct('Z',ZZ,'Z_thr',ZZ_thr,'thr',thr,'rois',{names},'pathlength',net_rd,'clustering',net_rc,'norm_pathlength',net_d,'norm_clustering',net_c,'degree',net_k,'pathlength_roi',NET_d,'clustering_roi',NET_c,'degree_roi',NET_k);

if ~isfield(ss,'ask')||isempty(ss.ask), 
    ss.ask='missing'; ss.askn=1;
else
    types={'none','missing','all'};typesn=[0,1,2];
    if isnumeric(ss.ask), sstype=ss.ask;
    else sstype=strmatch(lower(ss.ask),lower(types),'exact'); end
    ss.ask=types{sstype};
    ss.askn=typesn(sstype);
end

if ss.askn>1||((~isfield(ss,'model')||isempty(ss.model))&&(~isfield(ss,'X')||isempty(ss.X))), 
    if ~isfield(ss,'model')||isempty(ss.model), ss.model=1; end
    if ss.askn,
        if isnumeric(posstr)&&isempty(findobj(0,'tag','Interactive')), spm('CreateIntWin'); end;
        ss.model=spm_input('Select model',posstr,'m','one-sample t-test|two-sample t-test|multiple regression',[],ss.model);posstr='+1';
    end
end

if ss.askn>1||~isfield(ss,'X')||isempty(ss.X),
    switch(ss.model),
        case 1,
            ss.X=ones(ss.n,1);
            default_xname={'Group'};
            default_ccontrast={1};default_cname={'one-sided t-test'};
        case 2,
           if isnumeric(posstr)&&isempty(findobj(0,'tag','Interactive')), spm('CreateIntWin'); end;
           if ~isfield(ss,'X')||isempty(ss.X), ss.X=[]; end
            done=0;
            if isempty(ss.X), idx1=[]; else idx1=find(ss.X(:,1)); end
            while ~done,
                idx1=spm_input('Subjects in group 1?',posstr,'n',mat2str(idx1(:)'),[1,inf],ss.n);posstr='+1';
                idx2=spm_input('Subjects in group 2?',posstr,'n',mat2str(setdiff(1:ss.n,idx1(:)')),[1,inf],ss.n);posstr='+1';
                if numel(idx1)+numel(idx2)==ss.n&&isempty(intersect(idx1,idx2)),done=1;end
            end
            ss.X=zeros(ss.n,2);ss.X(idx1,1)=1;ss.X(idx2,2)=1;
            default_xname={'Group 1','Group 2'};
            default_ccontrast={};
        case 3,
            if isnumeric(posstr)&&isempty(findobj(0,'tag','Interactive')), spm('CreateIntWin'); end;
            if ~isfield(ss,'X')||isempty(ss.X), ss.X=[]; end
            ss.X=spm_input('Regressor matrix?',posstr,'r',mat2str(ss.X),[ss.n,inf]);posstr='+1';
            default_xname=cellstr([repmat('Regressor #',[size(ss.X,2),1]),num2str((1:size(ss.X,2))')]);
            default_ccontrast={};
        otherwise,
            default_xname={};
            default_ccontrast={};
    end
else
    default_xname=repmat({'regressor-name'},[size(ss.X,2),1]);
    default_ccontrast={};
end

if ss.askn>1||~isfield(ss,'Xname')||numel(ss.Xname)~=size(ss.X,2), 
    if ~isfield(ss,'Xname')||numel(ss.Xname)~=size(ss.X,2), ss.Xname=default_xname; end
    for nc=1:size(ss.X,2),
        if numel(ss.Xname)<nc,ss.Xname{nc}=['Regressor #',num2str(nc)];end
        if ss.askn,
            if isnumeric(posstr)&&isempty(findobj(0,'tag','Interactive')), spm('CreateIntWin'); end;
            ss.Xname{nc}=spm_input(['Regressor #',num2str(nc),'name?'],posstr,'s',ss.Xname{nc});posstr='+1';
        end
    end
end

if ~isfield(ss,'C2'), ss.C2=eye(numel(results)); end 

if ~isfield(ss,'C'), ss.C={}; end 
if ~isfield(ss,'Cname'), ss.Cname={}; end 
if ~isempty(default_ccontrast),
    for nc=1:numel(default_ccontrast),
        ss.C{nc}=default_ccontrast{nc};
        ss.Cname{nc}=default_cname{nc};
    end
end

% define contrast
if ss.askn>1||isempty(ss.C)
    Ic=1;
    if numel(ss.C)<Ic,ss.C{Ic}=[];end
    ss.C{Ic}=spm_input('Contrast vector/matrix?',posstr,'r',mat2str(ss.C{Ic}),[inf,size(ss.X,2)]);posstr='+1';
    if numel(ss.Cname)<Ic,ss.Cname{Ic}='contrast name';end
    ss.Cname{Ic}=spm_input('Contrast name?',posstr,'s',ss.Cname{Ic});posstr='+1';
else Ic=1; end

ss.askn=1;ss.ask='missing';
ss.estimate={}; ss.evaluate={}; 

% estimate model
if isfield(results(1),'type')&&results(1).type==1,
    if isfield(results(1),'measure_names')
        vars=fieldnames(results(1).measures);
    else
        vars={'GlobalEfficiency','LocalEfficiency','Cost','GlobalEfficiency_roi','LocalEfficiency_roi','Cost_roi'};
    end
else
    vars={'pathlength','clustering','norm_pathlength','norm_clustering','degree','pathlength_roi','clustering_roi','degree_roi','norm_pathlength_roi','norm_clustering_roi'};
end

X=ss.X;
C=ss.C{Ic};
Nh=size(C,1);

for nvar=1:numel(vars),
    Y=[];
    for n1=1:numel(results)
        if isfield(results(n1),'measure_names')
            Y=cat(3,Y,results(n1).measures.(vars{nvar}));
        else
            Y=cat(3,Y,results(n1).(vars{nvar}));
        end
    end
%     [h,F,p,dof,statsname]=conn_glm(X,permute(Y,[1,3,2]),C,ss.C2);
%     glmresults.h=permute(h,[2,3,1]);
%     glmresults.F=permute(F,[2,3,1]);
%     glmresults.p=permute(p,[2,3,1]);
%     glmresults.dof=dof;
%     glmresults.statsname=statsname;
    
%     Y=0; for n1=1:numel(results), 
%         if isfield(results(n1),'measure_names')
%             Y=Y+results(n1).measures.(vars{nvar})*ss.C2(:,n1)';
%         else
%             Y=Y+results(n1).(vars{nvar})*ss.C2(:,n1)';
%         end
%     end
    glmresults=struct('h',nan(1+0*Nh,size(Y,2)),'F',nan(1,size(Y,2)),'p',nan(1,size(Y,2)),'dof',nan(1,size(Y,2)));
    for nout=1:size(Y,2),
        idxsubjects=find(~any(isnan(X),2)&~any(isnan(Y(:,nout,:)),3));
        x=X(idxsubjects,:);
        y=Y(idxsubjects,nout,:);
        if ~any(isnan(y)),
            [h,F,p,dof,statsname]=conn_glm(x,permute(y,[1,3,2]),C,ss.C2);
            if size(h,2)>1, h=sqrt(sum(abs(h).^2,2)); end
            if size(h,1)>1, h=sqrt(sum(abs(h).^2,1)); end
            glmresults.h(:,nout)=h;
            glmresults.F(nout)=F;
            glmresults.p(nout)=p;
            if numel(dof)>size(glmresults.dof,1), glmresults.dof=cat(1,glmresults.dof,nan(numel(dof)-size(glmresults.dof,1),size(glmresults.dof,2))); end
            glmresults.dof(:,nout)=dof';
            glmresults.statsname=statsname;
        end
%         N=numel(idxsubjects);
%         ix=pinv(x'*x);
%         dof=N-rank(x);
%         r=C*ix*C';
%         Nc0=rank(x*C');
%         if ~any(isnan(y)),
%             b=ix*(x'*y);
%             e=y-x*b;
%             ee=sum(abs(e).^2,1);
%             
%             h=C*b;
%             if Nh==1,% T-stats
%                 k=sqrt(r*ee);
%                 F=real(h./max(eps,k))*sqrt(dof);
%                 if isnan(F)||dof<=0,F=nan;p=nan;else p=1-spm_Tcdf(F,dof);end
%                 edof=dof;
%             else % F-stats
%                 bb=h'*pinv(r)*h;
%                 F=real(bb./ee)*dof/Nc0;
%                 if isnan(F)||dof<=0||Nc0<=0,F=nan;p=nan;else p=1-spm_Fcdf(F,Nc0,dof); end
%                 edof=[Nc0,dof];
%             end
%             glmresults.h(:,nout)=h;
%             glmresults.F(nout)=F;
%             glmresults.p(nout)=p;
%             glmresults.dof(nout)=dof;
%         end
    end
    ss.evaluate{Ic}.(vars{nvar})=glmresults;
end
if ~isempty(filename), save(filename,'ss','-append'); end

% plot results
%vars={'pathlength','clustering','norm_pathlength','norm_clustering','degree','pathlength_roi','clustering_roi','degree_roi','norm_pathlength_roi','norm_clustering_roi'};
if ~isempty(filename),
    figure;
    for nvar=1:numel(vars),
        %figure;
        subplot(floor(sqrt(numel(vars))),ceil(numel(vars)/floor(sqrt(numel(vars)))),nvar);
        var=ss.evaluate{1}.(vars{nvar});
        F=var.F;
        %p=2*min(var.p,1-var.p);
        P=conn_fdr(p);
        alpha=.05;
        
        [nill,idx]=find(p<alpha); [nill,idxsort]=sort(p(idx));idx=idx(idxsort);
        if isempty(idx), idx=1:numel(F); end
        %[nill,idx]=sort(p);
        hb=bar(F(idx));set(hb,'facecolor',.75*[1,1,1],'edgecolor','none');hold on;
        hb=bar(F(idx).*(p(idx)<alpha)); set(hb,'facecolor',1*[1,1,0],'edgecolor','none');
        hb=bar(F(idx).*(P(idx)<alpha)); set(hb,'facecolor',1*[1,0,0],'edgecolor','none');
        hold off;
        set(gca,'xlim',[.5,numel(idx)+.5]);
        set(gca,'fontsize',12,'xtick',[]);
        %set(gca,'units','norm','position',[.1,.5,.8,.4]);
        ht=title(vars{nvar});set(ht,'interpreter','none','fontsize',16);
        ht=ylabel('T');set(ht,'interpreter','none','fontsize',14);
        miny=0;%get(gca,'ylim')*[1;0];
        str={'right','left'};for n1=1:numel(idx),if p(idx(n1))<alpha, ht=text(n1,miny,[' ',results(1).rois{idx(n1)}(1:min(10,numel(results(1).rois{idx(n1)}))),' ']);set(ht,'rotation',90,'horizontalalignment',str{round((3+sign(F(idx(n1))))/2)},'fontsize',9,'fontweight','bold');end; end
        set(gcf,'color','w');
    end
end
    
