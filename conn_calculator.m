function ok=conn_calculator(icovariates2)
% CONN_CALCULATOR
%
% Calculator allows you to explore and analyze individual and simple combinations 
% of subject-level measures/variables
% 
% Select the measure that you want to explore in the 'Measures' list
% 
% Select AllSubjects in the 'Subject effects' list to simply display 
% the values of your measure, or select any combination of subject-effects 
% there to perform a second-level general linear analysis (GLM) of your measure 
% of interest
% 
% You may add new measures to both of these lists by:
%
%  1) Manually defining them in Setup->Covariates->Second-level (e.g. behavioral or 
%     demographic measures)
% 
%  2) Aggregating first-level covariates (e.g. movement parameters) in the 
%     Setup->Covariates->First-level tab, by clicking on the 'compute summary 
%     measures' button
% 
%  3) Importing subject-level ROI-to-ROI connectivity values in the Second-level 
%     Results tab, by clicking on the 'Import values' button (e.g. connectivity 
%     between MPFC and PCC for each subject)
% 
%  4) Importing subject-level seed-to-voxel connectivity values in the seed-to-voxel 
%     results explorer window, by clicking on the 'Import values' button (e.g. average 
%     connectivity for each subject between MPFC and a cluster of interest)
% 
%  5) Importing subject-level voxel-to-voxel measures in the voxel-to-voxel results
%     explorer window, by clicking on the 'Import values' button (e.g. average intrinsic 
%     connectivity for each subject at a cluster of interest)''}
%

%$

global CONN_x CONN_gui CONN_h;

if ~isfield(CONN_gui,'font_offset'), font_offset=0; else font_offset=CONN_gui.font_offset; end
if ~isfield(CONN_gui,'backgroundcolor'), backgroundcolor=[0 0 0]; else backgroundcolor=CONN_gui.backgroundcolor; end
if CONN_x.Setup.nsubjects==1, conn_msgbox({'Single-subject second-level analyses not supported (only population-level inferences via random-effect analyses available)','Please add more subjects before proceeding'},'',2); end
X=zeros(CONN_x.Setup.nsubjects,length(CONN_x.Setup.l2covariates.names)-1);
for nsub=1:CONN_x.Setup.nsubjects,
    for ncovariate=1:length(CONN_x.Setup.l2covariates.names)-1;
        X(nsub,ncovariate)=CONN_x.Setup.l2covariates.values{nsub}{ncovariate};
    end
end
if ~nargin||isempty(icovariates2), 
    icovariates1all=find(cellfun(@(x)isempty(regexp(x,'^Dynamic |^_')),CONN_x.Setup.l2covariates.names(1:end-1)));
    %icovariates1all=1:numel(CONN_x.Setup.l2covariates.names)-1;
%     try
%         [nill,i]=sort(CONN_x.Setup.l2covariates.names(icovariates1all));
%         icovariates1all=icovariates1all(i);
%     end
    %icovariates2=1:numel(CONN_x.Setup.l2covariates.names)-1; 
    icovariates2=icovariates1all(cellfun(@(x)isempty(regexp(x,'^Dynamic |^_')),CONN_x.Setup.l2covariates.names(icovariates1all)));
    icovariates1=icovariates1all(cellfun(@(x)isempty(regexp(x,'^Dynamic |^_|^QA_|^QC_')),CONN_x.Setup.l2covariates.names(icovariates1all)));
    fulltype=1; 
    showall=false;
else
    icovariates1all=setdiff(1:numel(CONN_x.Setup.l2covariates.names)-1,icovariates2); 
%     try
%         [nill,i]=sort(CONN_x.Setup.l2covariates.names(icovariates1all));
%         icovariates1all=icovariates1all(i);
%     end
    icovariates1=icovariates1all(cellfun(@(x)isempty(regexp(x,'^Dynamic |^_|^QA_|^QC_')),CONN_x.Setup.l2covariates.names(icovariates1all)));
    cnames=regexp(CONN_x.Setup.l2covariates.names(icovariates2), '@.*$','match','once');
    cnames=regexprep(cnames,'^@\s*|^_(\S+) Dynamic Total.*','');
    cnames2=regexp(CONN_x.Setup.l2covariates.names(icovariates2),'^_(\S+) (ICA|PCA|Dynamic factor\s*)(\d+) ','match','once');
    cnames2=regexprep(cnames2,{'^_(\S+) (ICA|PCA|Dynamic factor\s*)(\d+) ','^ x (Variability|Frequency|Average)'},{' x $1',' x Temporal $1'});
    cnames=cellfun(@(a,b)[a b],cnames,cnames2,'uni',0);
    mnames=regexprep(CONN_x.Setup.l2covariates.names(icovariates2), '@.*$','');
    mnames=regexprep(mnames,{'^_(\S+) (ICA|PCA|Dynamic factor\s*)0*(\d+) .*','^Dynamic\s*|^_|^QA_|^QC_','Dynamic factor'},{'_$2_$3 ','','Circuit'});
    [ucnames,nill,icnames]=uniquestable(cnames);
    [umnames,nill,imnames]=unique(mnames);
    [ucnames,nill,idx]=conn_sortfilenames(ucnames); icnames=idx(icnames);
    [umnames,nill,idx]=conn_sortfilenames(umnames); imnames=idx(imnames);
    matchnames=full(sparse(imnames,icnames,1:numel(cnames)));
    cvalid=find(all(matchnames>0,1)); mvalid=find(all(matchnames>0,2));
    matchnames=matchnames(mvalid,cvalid);
    ucnames=ucnames(cvalid); umnames=umnames(mvalid);
    fulltype=0;
    showall=false;
end
X1=X(:,icovariates1);
X2=X(:,icovariates2);

if fulltype, 
    fs0=[.05 .25];
    fs=[.47 .55];
    fs0b=[fs0(1)+fs(1)+.05 fs0(2)];
    fsb=[1-fs(1)-fs0(1)-.05 fs(2)];
    %thfig=figure('units','norm','position',[.15,.3,.8,.5],'name','conn second-level calculator','numbertitle','off','menubar','none','color',[1 1 1]);
    conn_menu('frame',[fs0 fs],'Group analyses (2nd-level)');
    conn_menu('pushbutton',[fs0 0 0]+[fs fs].*[.90,1.05,.10,.05],'','Help','','conn(''gui_help'',''help'',''conn_calculator.m'');'); 
    ht1=conn_menu('listbox',[fs0 0 0]+[fs fs].*[.02,.22,.30,.63],'Subject effects',conn_strexpand(CONN_x.Setup.l2covariates.names(icovariates1),CONN_x.Setup.l2covariates.descrip(icovariates1)),'Select subject effect(s) characterizing second-level analysis model / independent variables',@(varargin)conn_calculator_update);
    ht3=conn_menu('edit',[fs0 0 0]+[fs fs].*[.02,.10,.30,.05],'Between-subjects contrast','1','<HTML>Define desired contrast across selected subject-effects<br/> - enter contrast vector/matrix with as many elements/columns as subject-effects selected <br/> - use the list below to see a list of standard contrasts for the selected subject-effects <br/> - enter multiple rows separated by <b>;</b> (semicolon) for OR conjunction (multivariate test) of several contrasts</HTML>',@(varargin)conn_calculator_update);
    ht2=conn_menu('listbox',[fs0 0 0]+[fs fs].*[.33,.22,.65,.63],'Dependent measure(s)',conn_strexpand(CONN_x.Setup.l2covariates.names(icovariates2),CONN_x.Setup.l2covariates.descrip(icovariates2)),'Select outcome measure(s) / dependent variables',@(varargin)conn_calculator_update);
    ht4=conn_menu('edit',[fs0 0 0]+[fs fs].*[.33,.10,.65,.05],'Between-measures contrast','1','<HTML>Define desired contrast across selected outcome measures <br/> - enter contrast vector/matrix with as many elements/columns as outcome measures selected <br/> - use the list below to see a list of standard contrasts for the selected outcome measures<br/> - enter multiple rows separated by <b>;</b> (semicolon) for OR conjunction (multivariate test) of several contrasts</HTML>',@(varargin)conn_calculator_update);
    hc1b=uicontextmenu;
    if showall, uimenu(hc1b,'Label','Hide secondary variables','callback',@(varargin)conn_calculator_update_hideshow);
    else uimenu(hc1b,'Label','Show secondary variables','callback',@(varargin)conn_calculator_update_hideshow);
    end
    set(ht1,'uicontextmenu',hc1b);
    set(ht1,'max',2,'value',1);
    [nill,i]=max(icovariates2);
    set(ht2,'max',2,'value',i);
    ht5=conn_menu('textedit2',[fs0 0 0]-[0 .15 0 0]+[fs fs].*[.02,0,.96,.10],'Statistics','','Effect size and statistics for selected second-level analysis');
else
    fs0=[.01 .335];
    fs=[.605 .505];
    fs0b=[fs0(1)+fs(1)+.05 fs0(2)];
    fsb=[1-fs(1)-fs0(1)-2*.05 fs(2)];
    %thfig=figure('units','norm','position',[.15,.3,.8,.5],'name','conn second-level calculator','numbertitle','off','menubar','none','color',[1 1 1]);
    conn_menu('frame',[fs0 fs+0*[0 .06]],'Group analyses (2nd-level)');%Second-level design');
    ht1=conn_menu('listbox',[fs0 0 0]+[fs fs].*[.02,.30,.30,.45],'Subject effects',conn_strexpand(CONN_x.Setup.l2covariates.names(icovariates1),CONN_x.Setup.l2covariates.descrip(icovariates1)),'Select subject effect(s) characterizing second-level analysis model / independent variables',@(varargin)conn_calculator_update);
    ht3=conn_menu('edit',[fs0 0 0]+[fs fs].*[.02,.15,.30,.05],'Between-subjects contrast','1','<HTML>Define desired contrast across selected subject-effects<br/> - enter contrast vector/matrix with as many elements/columns as subject-effects selected <br/> - use the list below to see a list of standard contrasts for the selected subject-effects <br/> - enter multiple rows separated by <b>;</b> (semicolon) for OR conjunction (multivariate test) of several contrasts</HTML>',@(varargin)conn_calculator_update);
    ht2c=conn_menu('listbox',[fs0 0 0]+[fs fs].*[.34,.30,.30,.45],'Conditions',ucnames,'Select condition(s)',@(varargin)conn_calculator_update);
    ht4c=conn_menu('edit',[fs0 0 0]+[fs fs].*[.34,.15,.30,.05],'Between-conditions contrast','1','<HTML>Define desired contrast across selected conditions <br/> - enter contrast vector/matrix with as many elements/columns as selected conditions <br/> - use the list below to see a list of standard contrasts for the selected condition<br/> - enter multiple rows separated by <b>;</b> (semicolon) for OR conjunction (multivariate test) of several contrasts</HTML>',@(varargin)conn_calculator_update);
    ht2=conn_menu('listbox',[fs0 0 0]+[fs fs].*[.66,.30,.30,.45],'Measures',umnames,'<HTML>Select outcome measure(s) <br/> - measures represent properties of timecourses for each subject/condition</HTML>',@(varargin)conn_calculator_update);
    ht4=conn_menu('edit',[fs0 0 0]+[fs fs].*[.66,.15,.30,.05],'Between-measures contrast','1','<HTML>Define desired contrast across selected measures <br/> - enter contrast vector/matrix with as many elements/columns as measures selected <br/> - use the list below to see a list of standard contrasts for the selected measures<br/> - enter multiple rows separated by <b>;</b> (semicolon) for OR conjunction (multivariate test) of several contrasts</HTML>',@(varargin)conn_calculator_update);
    hc1b=uicontextmenu;
    if showall, uimenu(hc1b,'Label','Hide secondary variables','callback',@(varargin)conn_calculator_update_hideshow('hide'));
    else uimenu(hc1b,'Label','Show secondary variables','callback',@(varargin)conn_calculator_update_hideshow('show'));
    end
    set(ht1,'uicontextmenu',hc1b);
    set(ht1,'max',2,'value',1);
    set(ht2c,'max',2,'value',1);
    set(ht2,'max',2,'value',1);
    ht5=conn_menu('textedit2',[fs0 0 0]-[0 .15 0 0]+[fs fs].*[.22,.0,.76,.10],'Analysis results','','Effect size and statistics for selected second-level analysis');
end
handleplot=[];
% try
% if domacbugfix,
%     set(findobj(thfig,'style','popupmenu','-or','style','edit','-or','style','pushbutton','-or','style','togglebutton'),'foregroundcolor',get(0,'defaultuicontrolforegroundcolor'),'backgroundcolor',get(0,'defaultuicontrolbackgroundcolor'));
% end
% end
defaultvars={1,2,0,0;2,1,0,0;2,2,0,0};
[plot_sortplot, plot_spacing, plot_labels, plot_shownull]=deal(defaultvars{1,:});
changedvars=false;
conn_calculator_update;


    function conn_calculator_update(newsortplot,newspacing,newlabels,newshownull,newchangedvars)
        if nargin>=1&&~isempty(newsortplot), plot_sortplot=newsortplot; changedvars=true; end
        if nargin>=2&&~isempty(newspacing), plot_spacing=newspacing; changedvars=true; end
        if nargin>=3&&~isempty(newlabels), plot_labels=newlabels; changedvars=true; end
        if nargin>=4&&~isempty(newshownull), plot_shownull=newshownull; changedvars=true; end
        if nargin>=5&&~isempty(newchangedvars), changedvars=newchangedvars; end
        ncov1=get(ht1,'value');
        str=get(ht3,'string');
        ccov1=conn_contraststr2num(str);
        if size(ccov1,2)~=numel(ncov1), ccov1=eye(numel(ncov1)); end
        set(ht3,'string',mat2str(ccov1));
        if fulltype
            ncov2=get(ht2,'value');
            str=get(ht4,'string');
            ccov2=conn_contraststr2num(str);
            if size(ccov2,2)~=numel(ncov2), ccov2=eye(numel(ncov2)); end
            set(ht4,'string',mat2str(ccov2));
        else
            ncov2m=get(ht2,'value');
            str=get(ht4,'string');
            ccov2m=conn_contraststr2num(str);
            if size(ccov2m,2)~=numel(ncov2m), ccov2m=eye(numel(ncov2m)); end
            set(ht4,'string',mat2str(ccov2m));
            ncov2c=get(ht2c,'value');
            str=get(ht4c,'string');
            ccov2c=conn_contraststr2num(str);
            if size(ccov2c,2)~=numel(ncov2c), ccov2c=eye(numel(ncov2c)); end
            set(ht4c,'string',mat2str(ccov2c));
            ncov2=matchnames(ncov2m,ncov2c);
            ccov2=kron(ccov2c,ccov2m);
        end
        nsubjects=find(any(X1(:,ncov1)~=0,2)&~any(isnan(X1(:,ncov1)),2)&~any(isnan(X2(:,ncov2)),2));
        [h,F,p,dof,statsname]=conn_glm(X1(nsubjects,ncov1),X2(nsubjects,ncov2),ccov1,ccov2);
        if F>1e10, F=inf; p=nan; end
        if isempty(h), str='';
        else
            str=sprintf('beta = %s   %s(%s) = %.2f   p = %.6f',mat2str(h,max(0,max(ceil(log10(1e-10+abs(h(:))))))+2),statsname,deblank(sprintf('%d ',dof)),F,p);
            if isequal(statsname,'T'), str=[str sprintf('  (two-sided p = %.6f)',2*min(p,1-p))]; end
        end
        set(ht5,'string',str);
        if fulltype, tstr='predictor variables'; else tstr='subject effects'; end
        conn_contrasthelp(ht3, tstr, get(ht1,'string'),get(ht1,'value'),all(ismember(X1(:,get(ht1,'value')),[0 1]),1)+2*all(ismember(X1(:,get(ht1,'value')),[-1 1]),1));
        if fulltype
            conn_contrasthelp(ht4,'outcome variables',get(ht2,'string'),get(ht2,'value'),all(ismember(X2(:,get(ht2,'value')),[0 1]),1)+2*all(ismember(X2(:,get(ht2,'value')),[-1 1]),1));
        else
            conn_contrasthelp(ht4c,'conditions',ucnames,ncov2c,[]);
            conn_contrasthelp(ht4,'measures',umnames,ncov2m,[]);
        end
        x=X1(nsubjects,ncov1)*ccov1';
        y=X2(nsubjects,ncov2)*ccov2';
        x0=X1(nsubjects,ncov1)-X1(nsubjects,ncov1)*ccov1'*pinv(ccov1*ccov1')*ccov1;
        y_fit=X1(nsubjects,ncov1)*(pinv(X1(nsubjects,ncov1))*X2(nsubjects,ncov2))*ccov2';
        y_fit0=x0*(pinv(x0)*X2(nsubjects,ncov2))*ccov2';
        hasconst=numel(unique(all(X1(nsubjects,ncov1)==1,1)))==2;
        hasconst0=numel(unique(all(x0==1,1)))==2;
        if ~changedvars
            nux=numel(unique(x));
            if nux<=1, [plot_sortplot, plot_spacing, plot_labels, plot_shownull]=deal(defaultvars{1,:});
            elseif nux<=2, [plot_sortplot, plot_spacing, plot_labels, plot_shownull]=deal(defaultvars{2,:});
            else [plot_sortplot, plot_spacing, plot_labels, plot_shownull]=deal(defaultvars{3,:});
            end
        end
        switch(plot_sortplot)
            case 1, idx=1:numel(nsubjects);  cr=idx(:); xlabelstr='Subjects';
            case 2, cr1=x; [cr,idx]=sortrows([cr1, sum(X1(nsubjects,ncov1),2)]); cr=cr1(idx,:); xlabelstr='Predictor / Subject effects';
            case 3, cr1=y_fit; [cr,idx]=sortrows([cr1, sum(X1(nsubjects,ncov1),2)]); cr=cr1(idx,:); xlabelstr='Outcome (fitted)';
            case 4, cr1=y; [cr,idx]=sortrows([cr1 sum(X1(nsubjects,ncov1),2)]); cr=cr1(idx,:); xlabelstr='Outcome (observed)';
        end
        switch(plot_spacing)
            case 1, xplot=(1:numel(idx))'; if plot_sortplot~=1, xlabelstr=sprintf('Subjects (sorted by %s)',xlabelstr); end
            case 2, xplot=sum(cr,2); %sum(cr(:,any(cr~=1,1)),2);
        end
        switch(plot_labels)
            case 0, xlabels=[];
            case 1, xlabels=nsubjects;
            case 2, xlabels=sum(x,2); xlabelstr='Predictor';
            case 3, xlabels=sum(y_fit,2); xlabelstr='Outcome (fitted)';
            case 4, xlabels=sum(y,2); xlabelstr='Outcome (observed)';
        end
        if any(ishandle(handleplot)),delete(handleplot(ishandle(handleplot))); end
        str2=CONN_x.Setup.l2covariates.names(icovariates2);
        %disp(char(str2(ncov2)));
        %disp(ccov2);
        %str2=get(ht2,'string');
        handleplot=[];
        %figure;
        for n=1:size(ccov2,1)
            handle=axes('tag','conn_calculator_axes');
            cla;
            hold on;
            plot(repmat(xplot(:)',2,1),[y(idx,n),y_fit(idx,n)]','k:','color',.5*[1 1 1]);
            plot(xplot,y(idx,n),'ro','markerfacecolor','r','markeredgecolor',[1 .25 .25],'markersize',max(1,2+font_offset));
            if plot_shownull, 
                plot(xplot,y_fit0(idx,n),'ko','color',.5*[1 1 1],'markerfacecolor',.5*[1 1 1],'markeredgecolor',.5*[1 1 1],'markersize',max(1,(1+font_offset)/2)); 
                %plot(repmat(xplot(:)',2,1),[y(idx,n),y_fit0(idx,n)]','k:','color',.5*[1 1 1]);
            end
            plot(xplot,y_fit(idx,n),'bo','markerfacecolor','b','markeredgecolor',[.25 .25 1],'markersize',max(1,(1+font_offset)/2));
            hold off;
            axis tight;
            handle=gca;
            set(handle,'units','norm','position',[fs0b 0 0]+[fsb fsb].*[.1 .2+.75-.75*n/size(ccov2,1) .8 .6/size(ccov2,1)],'xcolor',.7*[1 1 1],'ycolor',.7*[1 1 1],...
                'color',CONN_gui.backgroundcolor,'box','off',...
                'xlim',get(gca,'xlim')+max(1e-4,abs(diff(get(gca,'xlim'))))*[-.05 .05],'ylim',get(gca,'ylim')+max(1e-4,abs(diff(get(gca,'ylim'))))*[-.05 .05]);
            [uxplot,ixplot]=unique(xplot);
            if n<size(ccov2,1), set(handle,'xtick',[]);
            elseif ~isempty(xlabels)
                if numel(uxplot)<=30, set(handle,'xtick',uxplot,'xticklabel',arrayfun(@num2str,xlabels(idx(ixplot)),'uni',0));
                else set(handle,'xtick',uxplot(round(linspace(1,numel(uxplot),30))),'xticklabel',arrayfun(@num2str,xlabels(idx(ixplot(round(linspace(1,numel(uxplot),30))))),'uni',0));
                end
            end
            handleplot=[handleplot,handle];
            s1=get(handle,'xtick'); s2=get(handle,'xticklabel');
            if ~isempty(xlabels)&&~isempty(s1)
                hl1=text(s1,get(handle,'ylim')*[1.10;-.10]+zeros(size(s1)),s2,'color',.7*[1 1 1],'rotation',-90);
                set(handle,'xticklabel',[]);
                handleplot=[handleplot,hl1(:)'];
            end
            hc1=uicontextmenu;
            hc2=[uimenu(hc1,'Label','x-axis: Subjects','callback',@(varargin)conn_calculator_update(1)),...
                 uimenu(hc1,'Label','x-axis: Predictor values','callback',@(varargin)conn_calculator_update(2)),...
                 uimenu(hc1,'Label','x-axis: Outcome (fitted) values','callback',@(varargin)conn_calculator_update(3)),...
                 uimenu(hc1,'Label','x-axis: Outcome (observed) values','callback',@(varargin)conn_calculator_update(4))];
            set(hc2,'checked','off'); if plot_sortplot>=1&&plot_sortplot<=numel(hc2), set(hc2(plot_sortplot),'checked','on'); end
            hc2=[uimenu(hc1,'Label','x-axis spacing: Rank','callback',@(varargin)conn_calculator_update([],1),'separator','on'),...
                 uimenu(hc1,'Label','x-axis spacing: Proportional','callback',@(varargin)conn_calculator_update([],2))];
            set(hc2,'checked','off'); if plot_spacing>=1&&plot_spacing<=numel(hc2), set(hc2(plot_spacing),'checked','on'); end
            hc2=[uimenu(hc1,'Label','x-axis labels: Default','callback',@(varargin)conn_calculator_update([],[],0),'separator','on'),...
                 uimenu(hc1,'Label','x-axis labels: Subject IDs','callback',@(varargin)conn_calculator_update([],[],1)),...
                 uimenu(hc1,'Label','x-axis labels: Predictor values','callback',@(varargin)conn_calculator_update([],[],2)),...
                 uimenu(hc1,'Label','x-axis labels: Outcome (fitted) values','callback',@(varargin)conn_calculator_update([],[],3)),...
                 uimenu(hc1,'Label','x-axis labels: Outcome (observed) values','callback',@(varargin)conn_calculator_update([],[],4))];
            set(hc2,'checked','off'); if plot_labels+1>=1&&plot_labels+1<=numel(hc2), set(hc2(plot_labels+1),'checked','on'); end
%             uimenu(hc1,'Label','x-axis labels: predictor','callback',@(varargin)conn_calculator_update([],[],2));
%             uimenu(hc1,'Label','x-axis labels: outcome (fitted)','callback',@(varargin)conn_calculator_update([],[],3));
%             uimenu(hc1,'Label','x-axis labels: outcome (observed)','callback',@(varargin)conn_calculator_update([],[],4));
            hc2=[uimenu(hc1,'Label','Hide null-hypothesis','callback',@(varargin)conn_calculator_update([],[],[],0),'separator','on'),...
                 uimenu(hc1,'Label','Show null-hypothesis','callback',@(varargin)conn_calculator_update([],[],[],1))];
            set(hc2,'checked','off'); if plot_shownull+1>=1&&plot_shownull+1<=numel(hc2), set(hc2(plot_shownull+1),'checked','on'); end
            uimenu(hc1,'Label','Restore default display options','callback',@(varargin)conn_calculator_update([],[],[],[],0),'separator','on');
            uimenu(hc1,'Label','Print (high-res)','callback',@(varargin)conn_print);
            set(handle,'uicontextmenu',hc1);
            set(handle,'buttondownfcn',@(varargin)set(hc1,'position',get(0,'pointerlocation'),'visible','on'));
            handleplot=[handleplot,hc1];
            if all(ismember(ccov2(n,:),[0 1]))&&sum(ccov2(n,:)==1)==1, handle=title(str2{ncov2(find(ccov2(n,:)==1))}); 
            else  handle=title(['contrast ',mat2str(ccov2(n,:))]);
            end
            set(handle,'interpreter','none','color',.7*[1 1 1],'fontsize',8+font_offset);
            handleplot=[handleplot,handle];
            if n==size(ccov2,1)
                handle=xlabel(xlabelstr);
                handleplot=[handleplot,handle];
            end
            tstr={};
            if hasconst, r2=max(0,1-sum(abs(y(idx,n)-y_fit(idx,n)).^2)/max(eps,sum(abs(detrend(y(idx,n),'constant')).^2)));   tstr{end+1}=sprintf('%sR^2 = %0.2f','\color[rgb]{.25 .25 1}',r2); end
            if hasconst0&&plot_shownull, r2=max(0,1-sum(abs(y(idx,n)-y_fit0(idx,n)).^2)/max(eps,sum(abs(detrend(y(idx,n),'constant')).^2))); tstr{end+1}=sprintf('%sR^2 = %0.2f','\color[rgb]{.5 .5 .5}',r2); end
            if hasconst
                handle=text(get(gca,'xlim')*[.95;.05],get(gca,'ylim')*[0;1],tstr,'fontsize',8+font_offset,'horizontalalignment','left','verticalalignment','top');
                handleplot=[handleplot,handle];
            end
            set(gca,'fontsize',7+font_offset);
        end
        handle=axes('units','norm','position',[fs0b 0 0]+[fsb fsb].*[.1 0 .8 .10],'visible','off');
        handleplot=[handleplot,handle];
        if plot_shownull, handle=text(0,0,{'red dots: observed values','blue dots: fitted values','gray dots: fitted values under null hypothesis (between-subjects contrast = 0)'},'color',.7*[1 1 1],'fontsize',8+font_offset,'horizontalalignment','left');
        else handle=text(0,0,{'red dots: observed values','blue dots: fitted values'},'color',.7*[1 1 1],'fontsize',8+font_offset,'horizontalalignment','left');
        end
        set(gca,'xlim',[0 1]);
        handleplot=[handleplot,handle];
        try, if isfield(CONN_h,'menus')&&isfield(CONN_h.menus,'waiticonObj'), CONN_h.menus.waiticonObj.stop; end; end
    end

    function conn_calculator_update_hideshow(varargin)
        showall=~showall;
        if showall, set(get(hc1b,'children'),'label','Hide secondary variables');
            icovariates1=icovariates1all;
        else
            set(get(hc1b,'children'),'Label','Show secondary variables');
            icovariates1=icovariates1all(cellfun(@(x)isempty(regexp(x,'^Dynamic |^_|^QA_|^QC_')),CONN_x.Setup.l2covariates.names(icovariates1all)));
        end
        set(ht1,'string',conn_strexpand(CONN_x.Setup.l2covariates.names(icovariates1),CONN_x.Setup.l2covariates.descrip(icovariates1)),'value',unique(min(numel(icovariates1), get(ht1,'value'))));
        X1=X(:,icovariates1);
        conn_calculator_update;
    end
end

function [a,b,c]=uniquestable(x)
[a,b,c]=unique(x);
[nill,idx]=sort(accumarray(c(:),(1:numel(c))',[],@min));
iidx=idx; iidx(idx)=1:numel(idx);
a=a(idx);
b=b(idx);
c=iidx(c);
end

function str=conn_strexpand(varargin)
global CONN_gui;
if nargin<1, str={}; return; end
str=varargin{1};
changed=false(size(str));
if isfield(CONN_gui,'isjava')&&CONN_gui.isjava,
    for n1=2:nargin
        for n2=1:min(numel(str),numel(varargin{n1}))
            if ~isempty(varargin{n1}{n2}), changed(n2)=true; str{n2}=[str{n2} ' <i>(' varargin{n1}{n2} ')</i>']; end
        end
    end
    for n2=find(changed(:))'
        str{n2}=['<HTML>' str{n2} '</HTML>'];
    end
end
end



