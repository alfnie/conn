function conn_convertl12l2covariate(varargin)
% CONN_CONVERTL12L2COVARIATE
%  computes summary measures (creates second-level covariates summarizing properties of first-level covariates)
%  

global CONN_x CONN_gui;
if ~isfield(CONN_gui,'font_offset'), conn_font_init; end

filepath=CONN_x.folders.data;
filename=fullfile(filepath,['COND_Subject',num2str(1,'%03d'),'_Session',num2str(1,'%03d'),'.mat']);
if ~conn_existfile(filename), uiwait(warndlg(['Not ready to compute summary measures yet. Please run Setup step first (minimally: conn_process setup_conditions)'],'')); return; end

if nargin<1||isempty(varargin{1}), 
    ncovariates=listdlg('liststring',CONN_x.Setup.l1covariates.names(1:end-1),'selectionmode','multiple','promptstring',{'Select first-level covariate(s)'},'ListSize',[200 100]);
    if isempty(ncovariates), return; end
else
    ncovariates=varargin;
    if iscell(ncovariates)
        if ischar(ncovariates{1})
            [ok,ncovariates]=ismember(ncovariates,CONN_x.Setup.l1covariates.names(1:end-1));
        else
            ncovariates=cell2mat(ncovariates);
        end
    end
end
ncovariates0=length(CONN_x.Setup.l1covariates.names);
ncovariates=setdiff(ncovariates,[0 ncovariates0]);
nconditions=numel(CONN_x.Setup.conditions.names)-1;
X=repmat({},CONN_x.Setup.nsubjects,numel(ncovariates));
ndims=nan(CONN_x.Setup.nsubjects,numel(ncovariates));
equaldims=true(CONN_x.Setup.nsubjects,numel(ncovariates));

transforms={'raw values',         @(x)x,         @(w)w,     'x_k(n)';
           'absolute values',     @(x)abs(x),    @(w)w,     '|x_k(n)|';
           'squared values',      @(x)abs(x).^2, @(w)w,     'x_k(n)^2';
           'other function of raw values',1,@(w)w,'f(x_k(n))';
           'scan-to-scan differences',          @(x)diff(x,1,1),                @(w)~convn(w<=0,[1;1],'valid'),     'x_k(n)-x_k(n-1)';
           'absolute scan-to-scan differences', @(x)abs(diff(x,1,1)),           @(w)~convn(w<=0,[1;1],'valid'),     '|x_k(n)-x_k(n-1)|';
           'squared scan-to-scan differences',  @(x)abs(diff(x,1,1)).^2,        @(w)~convn(w<=0,[1;1],'valid'),     '|x_k(n)-x_k(n-1)|^2';
           'other function of scan-to-scan differences',2,@(w)~convn(w<=0,[1;1],'valid'),'f(x_k(n)-x_k(n-1))'};
       
measures={ 'average',           @(x,w,dim)sum(x.*w,dim)./max(eps,sum(w,dim)),   @(x,dim,k1,k2,k3)['\frac{1}{',k3,'}\sum_{',dim,'=',k1,'}^{',k2,'} ',x];
           'sum',               @(x,w,dim)sum(x.*(w>0),dim),                    @(x,dim,k1,k2,k3)['\sum_{',dim,'=',k1,'}^{',k2,'} ',x,' '];
           'standard deviation',@(x,w,dim)sqrt(max(0,sum(x.^2.*w,dim)./max(eps,sum(w,dim))-(sum(x.*w,dim)./max(eps,sum(w,dim))).^2)), @(x,dim,k1,k2,k3)['\frac{1}{',k3,'}\sum_{',dim,'=',k1,'}^{',k2,'} (',x,')^2 - (\frac{1}{',k3,'}\sum_{',dim,'=',k1,'}^{',k2,'} ',x,')^2'];
           'minimum',           @(x,w,dim)rem(min(x+(w<=0).*1e20,[],dim),1e20), @(x,dim,k1,k2,k3)['\min_',dim,'\{ ',x,' \}'];
           'maximum',           @(x,w,dim)rem(max(x-(w<=0).*1e20,[],dim),1e20), @(x,dim,k1,k2,k3)['\max_',dim,'\{ ',x,' \}'];
           'any/union',         @(x,w,dim)any(x&(w>0),dim),                     @(x,dim,k1,k2,k3)['\bigcup_',dim,'\{ ',x,' \}'];
           'all/intersection',  @(x,w,dim)all(x&(w>0),dim),                     @(x,dim,k1,k2,k3)['\bigcap_',dim,'\{ ',x,' \}'];
           'user-defined',      @(x,w,dim)userdefined(x,w,dim),                 @(x,dim,k1,k2,k3)['g'+(k2=='N'),'(\{ ',x,' | ',k1,'\le ',dim,'\le ',k2,'\} )'];
           'weighted sum',      @(x,w,dim)sum(x.*w,dim),                        @(x,dim,k1,k2,k3)['\sum_{',dim,'=',k1,'}^{',k2,'} w_',dim,'\cdot ',x,' '];
           'do not aggregate',  @(x,w,dim)x,                                    @(x,dim,k1,k2,k3)x };
       
measures_step1=measures(1:end-2,:);
if numel(ncovariates)>1, covname=sprintf('%s ',CONN_x.Setup.l1covariates.names{ncovariates}); covnames1='Consider covariates:'; covnames2='first-level covariate'; 
else covname=CONN_x.Setup.l1covariates.names{ncovariates}; covnames1=sprintf('Consider covariate %s:',covname); covnames2=sprintf('%s covariate',regexprep(covname,'[^a-zA-Z0-9]',' ')); 
end
thfig=conn_dialog('units','norm','position',[.3,.3,.3,.6],'windowstyle','normal','name',['Summarize 1st-level covariate ',covname],'color','w','resize','on');
bg=.9*[1 1 1];
uicontrol(thfig,'style','frame','units','norm','position',[0,.4,1,.6],'backgroundcolor',bg,'foregroundcolor',bg);
uicontrol(thfig,'style','text','units','norm','position',[.1,.90,.8,.05],'string',covnames1,'horizontalalignment','left','backgroundcolor',bg,'fontsize',9+CONN_gui.font_offset,'fontweight','bold');
ht0a=uicontrol(thfig,'style','popupmenu','units','norm','position',[.1,.85,.8,.05],'string',transforms(:,1),'fontsize',8+CONN_gui.font_offset,'horizontalalignment','left','backgroundcolor',bg,'tooltipstring',...
    '<HTML>Define the original measure that you would like summarized across timepoints and/or dimensions<br/> - e.g. <i>raw values</i> uses the original first-level covariate raw values : x<br/> - e.g. <i>absolute values</i> uses the absolute values of the original first-level covariate : |x|<br/> - e.g. <i>scan-to-scan differences</i> uses the framewise differences between the original first-level covariate values : x(n) - x(n-1)</HTML>','callback',@conn_convertl12l2covariate_update);
hs1a=uicontrol(thfig,'style','text','units','norm','position',[.1,.78,.8,.05],'string','Summarize across timepoints:','horizontalalignment','left','backgroundcolor',bg,'fontsize',9+CONN_gui.font_offset,'fontweight','bold');
ht1a=uicontrol(thfig,'style','popupmenu','units','norm','position',[.1,.73,.8,.05],'string',measures_step1(:,1),'fontsize',8+CONN_gui.font_offset,'horizontalalignment','left','backgroundcolor',bg,'tooltipstring',...
    '<HTML>Define how to summarize the measure of interest above across multiple timepoints/scans and multiple sessions<br/> - e.g. <i>average</i> will compute the average of this covariate across all timepoints/scans</HTML>','callback',@conn_convertl12l2covariate_update);
hs1b=uicontrol(thfig,'style','text','units','norm','position',[.1,.67,.8,.05],'string','Summarize across dimensions:','horizontalalignment','left','backgroundcolor',bg,'fontsize',9+CONN_gui.font_offset,'fontweight','bold');
ht1b=uicontrol(thfig,'style','popupmenu','units','norm','position',[.1,.62,.8,.05],'string',measures(:,1),'fontsize',8+CONN_gui.font_offset,'horizontalalignment','left','backgroundcolor',bg,'tooltipstring',...
    '<HTML>Define how to summarize the measure of interest above across multiple dimensions/components (for covariates that contain more than one timeseries)<br/> - e.g. <i>average</i> will compute the average of this covariate dimensions/components<br/> - e.g. <i>weighted sum</i> will compute the weighted average of this covariate dimensions/components (using a user-defined weight for each component)<br/> - e.g. <i>do not aggregate</i> will create a separate measure / second-level covariate for each dimension/component of this first-level covariate</HTML>','callback',@conn_convertl12l2covariate_update);
ht1c=uicontrol(thfig,'style','checkbox','units','norm','position',[.1,.56,.8,.05],'string','summarize across dimensions first','fontsize',8+CONN_gui.font_offset,'backgroundcolor',bg,'horizontalalignment','left','tooltipstring',...
    '<HTML>Performs aggregation across covariate dimensions first, followed by aggregation across timepoints<br/> - by default (when this options is unchecked) first-level covariates are first aggregated across multiple timepoints, followed by aggregation across multiple dimensions</HTML>','value',0,'callback',@conn_convertl12l2covariate_update);
ht2=uicontrol(thfig,'style','listbox','units','norm','position',[.1,.41,.8,.10],'max',2,'string',CONN_x.Setup.conditions.names(1:end-1),'value',1:numel(CONN_x.Setup.conditions.names)-1,'fontsize',8+CONN_gui.font_offset,'backgroundcolor',bg,'horizontalalignment','left','tooltipstring','Compute summary measures for each of the selected conditions','callback',@conn_convertl12l2covariate_update);
ht2a=uicontrol(thfig,'style','checkbox','units','norm','position',[.1,.51,.8,.05],'string','condition-specific measures','fontsize',8+CONN_gui.font_offset,'backgroundcolor',bg,'horizontalalignment','left','tooltipstring','Computes summary measures separately for each condition','value',0,'callback',@conn_convertl12l2covariate_update);
uicontrol(thfig,'style','text','units','norm','position',[.1,.325,.8,.05],'string','Summary measure:','horizontalalignment','left','backgroundcolor','w','fontsize',9+CONN_gui.font_offset,'fontweight','bold');
ha=axes('units','norm','position',[.1 .225 .8 .10],'visible','off','parent',thfig);
ht3=text(0,0,' ','horizontalalignment','center','backgroundcolor','w','fontsize',9+CONN_gui.font_offset,'interpreter','latex','parent',ha);
set(ha,'xlim',[-1 1],'ylim',[-1 1],'visible','off');
ha=axes('units','norm','position',[.1 .125 .8 .10],'visible','off','parent',thfig);
text(0,0,{['$$x_k(n)$$ : ',covnames2],'$$n$$ : timepoints $$1\le n\le N$$','$$k$$ : dimensions $$1 \le k \le K$$'},'horizontalalignment','center','backgroundcolor','w','fontsize',9+CONN_gui.font_offset,'interpreter','latex','parent',ha);
set(ha,'xlim',[-1 1],'ylim',[-1 1],'visible','off');

uicontrol(thfig,'style','pushbutton','string','Ok','units','norm','position',[.1,.025,.26,.07],'callback',@(varargin)conn_convertl12l2covariate_ok('ok'),'fontsize',8+CONN_gui.font_offset,'tooltipstring','Compute summary measures and import them as new 2nd-level covariates');
uicontrol(thfig,'style','pushbutton','string','Apply','units','norm','position',[.37,.025,.26,.07],'callback',@(varargin)conn_convertl12l2covariate_ok('apply'),'fontsize',8+CONN_gui.font_offset,'fontweight','bold','tooltipstring','Compute summary measures and import them as new 2nd-level covariates');
uicontrol(thfig,'style','pushbutton','string','Cancel','units','norm','position',[.64,.025,.26,.07],'callback','delete(gcbf)','fontsize',8+CONN_gui.font_offset);
conn_convertl12l2covariate_update;

    function conn_convertl12l2covariate_ok(opt,varargin)
        if nargin<1||isempty(opt), opt='ok'; end
        ok=ishandle(thfig);
        if ok,
            newtransform=get(ht0a,'value');
            newmeasure1=get(ht1a,'value');
            newmeasure2=get(ht1b,'value');
            f2first=get(ht1c,'value');
            nconditions=get(ht2,'value');
            if ~get(ht2a,'value'), nconditions=0; end
            %delete(thfig);
            f0x=transforms{newtransform,2};
            f0w=transforms{newtransform,3};
            f1=measures{newmeasure1,2};
            f2=measures{newmeasure2,2};
            measurename0=transforms{newtransform,1};
            measurename1=measures{newmeasure1,1};
            measurename2=measures{newmeasure2,1};
            donotaggregate=newmeasure2==size(measures,1);
            if donotaggregate||isequal(measurename1,measurename2), measurename=measurename1;
            else measurename=[measurename1,'/',measurename2];
            end
            if strcmp(measurename2,'weighted sum'), WEIGHTS=[];
            else WEIGHTS=1;
            end
            userdefined('clear');
            
            Z=cell(numel(nconditions),numel(ncovariates));
            [samples,weights,names]=deal({});
            ht=conn_waitbar(0,'Loading covariate/condition info. Please wait...',false);
            for nsub=1:CONN_x.Setup.nsubjects
                nsess=CONN_x.Setup.nsessions(min(length(CONN_x.Setup.nsessions),nsub));
                Y=cell(numel(nconditions),numel(ncovariates));
                W=cell(numel(nconditions),numel(ncovariates));
                for nses=1:nsess
                    filename=fullfile(filepath,['COND_Subject',num2str(nsub,'%03d'),'_Session',num2str(nses,'%03d'),'.mat']);
                    if ~conn_existfile(filename), conn_disp(['Not ready. Please run Setup step first (minimally: conn_process setup_conditions)']); return; end
                    conn_loadmatfile(filename,'samples','weights','names');
                    if ~isequal(CONN_x.Setup.conditions.names(1:end-1),names), conn_waitbar('close',ht); conn_msgbox(['Incorrect conditions in file ',filename,'. Please re-run Setup step first (minimally: conn_process setup_conditions)'],'error summarizing first-level covariate',2); return; end
                    for icondition=1:numel(nconditions)
                        for icovariate=1:numel(ncovariates)
                            nl1covariate=ncovariates(icovariate);
                            if 1,%~isempty(CONN_x.Setup.l1covariates.files{nsub}{nl1covariate}{nses}{3})
                                x=conn_get_l1covariate(nsub,nl1covariate,nses); %CONN_x.Setup.l1covariates.files{nsub}{nl1covariate}{nses}{3};
                                if nconditions(icondition),
                                    y=x(samples{nconditions(icondition)},:);
                                    w=weights{nconditions(icondition)}{3};
                                else
                                    y=x;
                                    w=ones(size(x,1),1);
                                end
                                if isnumeric(f0x)
                                    switch(f0x)
                                        case 1,
                                            answ=conn_menu_inputdlg('function f(x) of raw values','',1,{'sin(x)'});
                                            if isempty(answ), conn_waitbar('close',ht); return; end
                                            g0x=str2num(['@(x)(',answ{1},')']);
                                            f0x=@(x)g0x(x);
                                        case 2,
                                            answ=conn_menu_inputdlg('function f(x) of scan-to-scan differences','',1,{'sin(x)'});
                                            if isempty(answ), conn_waitbar('close',ht); return; end
                                            g0x=str2num(['@(x)(',answ{1},')']);
                                            f0x=@(x)g0x(diff(x,1,1));
                                    end
                                    
                                end
                                y=f0x(y);
                                w=f0w(w);
                                if size(y,2)~=1,
                                    if f2first
                                        while isempty(WEIGHTS)
                                            wtemp=ones(1,size(y,2));
                                            answ=conn_menu_inputdlg(sprintf('Weights (%d values)',size(y,2)),'',1,{mat2str(wtemp)});
                                            if isempty(answ), conn_waitbar('close',ht); return; end
                                            WEIGHTS=str2num(answ{1});
                                            if numel(WEIGHTS)~=size(y,2), WEIGHTS=[]; end
                                        end
                                        Wtemp=WEIGHTS(repmat(1+mod(0:size(y,2)-1,numel(WEIGHTS)),size(y,1),1));
                                        y = f2(y,Wtemp,2); % aggregate across dimensions
                                        w = repmat(w,1,size(y,2));
                                    elseif 1,%donotaggregate
                                        w=repmat(w,1,size(y,2));
                                    else
                                        w=repmat(w,size(y,2),1);
                                        y=y(:);
                                    end
                                end
                                if ~isempty(y)
                                    if ~isempty(Y{icondition,icovariate}) && size(Y{icondition,icovariate},2)~=size(y,2), conn_waitbar('close',ht); conn_msgbox({'Covariate does not contain the same number of dimensions across subjects/sessions','Summarize across dimensions first to avoid this issue'},'error summarizing first-level covariate',2); return; end
                                    Y{icondition,icovariate}=cat(1,Y{icondition,icovariate},y);
                                    W{icondition,icovariate}=cat(1,W{icondition,icovariate},w);
                                end
                            end
                        end
                    end
                end
                for icondition=1:numel(nconditions)
                    for icovariate=1:numel(ncovariates)
                        temp=f1(Y{icondition,icovariate},W{icondition,icovariate},1); % aggregate across timepoints
                        if ~f2first
                            while isempty(WEIGHTS)
                                wtemp=ones(1,size(temp,2));
                                answ=conn_menu_inputdlg(sprintf('Weights (%d values)',size(temp,2)),'',1,{mat2str(wtemp)});
                                if isempty(answ), conn_waitbar('close',ht); return; end
                                WEIGHTS=str2num(answ{1});
                                if numel(WEIGHTS)~=size(temp,2), WEIGHTS=[]; end
                            end
                            Wtemp=WEIGHTS(repmat(1+mod(0:size(temp,2)-1,numel(WEIGHTS)),size(temp,1),1));
                            temp=f2(temp,Wtemp,2); % aggregate across dimensions
                        end
                        if ~isempty(temp), Z{icondition,icovariate}(nsub,:)=temp; else  Z{icondition,icovariate}(nsub,:)=0; end
                    end
                end
                conn_waitbar(nsub/CONN_x.Setup.nsubjects,ht);
            end
            conn_waitbar('close',ht);
            names=cell(numel(nconditions),numel(ncovariates));
            for icondition=1:numel(nconditions)
                for icovariate=1:numel(ncovariates)
                    if nconditions(icondition), names{icondition,icovariate}=sprintf('%s of %s %s at %s',measurename,CONN_x.Setup.l1covariates.names{ncovariates(icovariate)},measurename0,CONN_x.Setup.conditions.names{nconditions(icondition)});
                    else names{icondition,icovariate}=sprintf('%s of %s %s',measurename,CONN_x.Setup.l1covariates.names{ncovariates(icovariate)},measurename0);
                    end
                end
            end
            conn_importl2covariate(names(:),Z(:));
            if strcmpi(opt,'ok'), delete(thfig); end
        end
    end

    function conn_convertl12l2covariate_update(varargin)
        offon={'off','on'};
        nk={'n','k'};
        n1=get(ht0a,'value');
        n2=get(ht1a,'value');
        n3=get(ht1b,'value');
        order=get(ht1c,'value');
        conds=get(ht2a,'value');
        str=transforms{n1,4};
        K1='1';K2='K';K3='K';
        if ~isempty(regexp(transforms{n1,1},'scan-to-scan')), N1='2';N2='N';N3='N-1';
        else N1='1';N2='N';N3='N';
        end
        if order 
            str=['$$y = ',measures{n2,3}(measures{n3,3}(str,'k',K1,K2,K3),'n',N1,N2,N3),'$$'];
            set(hs1b,'position',[.1,.78,.8,.05]);
            set(ht1b,'position',[.1,.73,.8,.05]);
            set(hs1a,'position',[.1,.66,.8,.05]);
            set(ht1a,'position',[.1,.61,.8,.05]);
        else
            str=['$$y = ',measures{n3,3}(measures{n2,3}(str,'n',N1,N2,N3),'k',K1,K2,K3),'$$'];
            set(hs1a,'position',[.1,.78,.8,.05]);
            set(ht1a,'position',[.1,.73,.8,.05]);
            set(hs1b,'position',[.1,.66,.8,.05]);
            set(ht1b,'position',[.1,.61,.8,.05]);
        end
        set(ht2,'visible',offon{1+conds});
        set(ht3,'string',str);
    end
end

function y=userdefined(x,w,dim)
persistent f;
y=[];
if nargin==1, % clear
    f={};
else
    if isempty(f)&&~iscell(f), f={}; end
    if numel(f)<dim||isempty(f{dim})
        fcn={'h(x)','g(x)'};
        str={'timepoints','dimensions'};
        answ=conn_menu_inputdlg(sprintf('user-defined function %s aggregating across %s',fcn{dim},str{dim}),'',1,{sprintf('mean(x,%d)',dim)});
        if isempty(answ), return; end
        f{dim}=str2num(['@(x)(',answ{1},')']);
    end
    X=num2cell(x,dim);
    W=num2cell(w,dim);
    y=cellfun(@(y,w)f{dim}(y(w>0)),X,W);
end
end

        
