function conn_displaynetwork(option,varargin)
global CONN_x CONN_gui;
if ~nargin, option='init'; end
if ~ischar(option), option=varargin{2}; end;
cconditions=1;
switch(lower(option)),
    case 'init',
        if nargin>1, 
            data.source=varargin{1}; 
            data.side=3;
            h=conn_msgbox('computing adjacency matrices, please wait...','conn_displaynetwork');
            results=conn_process('results_roi');
            %data.cconditions=results(1).c2;%CONN_x.Results.xX.cconditions;
            close(h);
%             if 0,%numel(CONN_x.Results.xX.nconditions)>1
%                 h=conn_msgbox('computing adjacency matrices, please wait...','conn_displaynetwork');
%                 temp1=CONN_x.Results.xX.nconditions;
%                 temp2=CONN_x.Results.xX.cconditions;
%                 for n1=1:numel(CONN_x.Results.xX.nconditions)
%                     CONN_x.Results.xX.nconditions=temp1(n1);
%                     CONN_x.Results.xX.cconditions=1;
%                     results(n1,:)=conn_process('results_roi');
%                 end
%                 CONN_x.Results.xX.nconditions=temp1;
%                 CONN_x.Results.xX.cconditions=temp2;
%                 %data.cconditions=CONN_x.Results.xX.cconditions;
%                 close(h);
%             else
%             end
        else
            data.side=3;
            [filename,filepath]=uigetfile('*ROI*.mat');
            if ~ischar(filename), return; end
            results=conn_loadmatfile(fullfile(filepath,filename));results=results.ROI;
            if ~isfield(results(1),'c2'), results(1).c2=[]; end
            %data.cconditions=results(1).c2;
            data.source=0;%1:length(results);
        end
%         if nargin<2 || isempty(varargin{1}), 
%             answ=listdlg('Promptstring','Select ROIs','selectionmode','multiple','liststring',CONN_x.Analyses(CONN_x.Analysis).sources);
%             nsources=answ; 
%             data.source=0;
%         elseif iscell(varargin{1}), nsources=varargin{1}{1}; data.source=varargin{1}{2};
%         else, nsources=varargin{1}; data.source=0; end
        data.Z=permute(cat(4,results.y),[4,2,1,3]); % roi x roi x subject x conditions 
        %data.Z=permute(reshape(cat(3,results.y),[size(results(1).y),size(results)]),[4,2,1,3]); % roi x roi x subject x conditions
        data.names=results(1).names;
        data.names2=results(1).names2;
        %data.names2reduced=data.names2;[s1,s2]=regexp(data.names2,'\d+ \([LR]\)','start','end');for n1=1:numel(data.names2),if ~isempty(s1{n1}),data.names2reduced{n1}=[data.names2{n1}(s1{n1}:s2{n1}-4),data.names2{n1}(s2{n1}-1)];end;end
        data.namesreduced=regexprep(data.names2,{'^BA\.(\d+) \(([LR])\)\. .*','^\((-?\d+),(-?\d+),(-?\d+)\)$','^SLrois\.|^aal\.|^atlas\.|^networks\.','\s\(([LlRr])\)','([^\(\)]+)\(.+\)\s*$'},{'$1$2','($1 $2 $3)','',' ${lower($1)}','$1'});
        data.xyz=cat(1,results(1).xyz{:});
        data.xyz2=cat(1,results(1).xyz2{:});
        data.Zthr=.15;
        data.Zthrtype=3;
        data.Zthrside=1;
        data.thr=.05;
        data.thrtype=2;
        data.clusters=[];
        data.view=[];
        data.proj=[];
        data.x=[];
        data.y=[];
        data.z=[];
        data.bgz=0;
        data.display='connectivity';
        data.displayreduced=0;
        data.displaytheserois=1:length(data.names);
        data.displaylabels=0;
        data.displayeffectsize=0;
        data.display3d=0;
        data.visible='on';
        data.disptype=0;
        
        Zthr=data.Zthr; if data.Zthrtype==1,Zthr=atanh(Zthr);end
        if ~isfield(data,'Zthrside')||data.Zthrside==1,
            for n1=1:size(data.Z,4)
                data.Zresults(n1)=conn_network({data.Z(:,:,:,n1),data.names},data.displaytheserois,1,data.Zthrtype-1,Zthr);
            end
        elseif data.Zthrside==2,
            for n1=1:size(data.Z,4)
                data.Zresults(n1)=conn_network({-data.Z(:,:,:,n1),data.names},data.displaytheserois,1,data.Zthrtype-1,Zthr);
            end
        else
            for n1=1:size(data.Z,4)
                data.Zresults(n1)=conn_network({abs(data.Z(:,:,:,n1)),data.names},data.displaytheserois,1,data.Zthrtype-1,Zthr);
            end
        end
        data.Zmodel=struct('X',results(1).xX.X,'Xname',{results(1).xX.name},'C',{{results(1).c}},'Cname',{{results(1).cname}},'C2',results(1).c2);
        data.Zmodel=conn_network_results({data.Zresults,data.Zmodel});
        fnames=fieldnames(data.Zmodel.evaluate{1});
        data.hnames=fnames;
        data.hnamesroi=zeros(1,numel(fnames));
        data.h=[];data.F=[];data.p=[];data.dof=[];data.h_net=[];data.F_net=[];data.p_net=[];data.dof_net=[];
        for nfnames=1:numel(fnames),
            if numel(fnames{nfnames})>4&&strcmp(fnames{nfnames}(end-3:end),'_roi'),
                data.h=cat(1,data.h,data.Zmodel.evaluate{1}.(fnames{nfnames}).h);
                data.F=cat(1,data.F,data.Zmodel.evaluate{1}.(fnames{nfnames}).F);
                data.p=cat(1,data.p,data.Zmodel.evaluate{1}.(fnames{nfnames}).p);
                data.dof=cat(1,data.dof,data.Zmodel.evaluate{1}.(fnames{nfnames}).dof);
                data.hnames{nfnames}=['',fnames{nfnames}(1:end-4)];
                %data.hnames{nfnames}=['Analysis of ',fnames{nfnames}(1:end-4)];
                data.hnamesroi(nfnames)=1;
            else
                data.h_net=cat(1,data.h_net,data.Zmodel.evaluate{1}.(fnames{nfnames}).h);
                data.F_net=cat(1,data.F_net,data.Zmodel.evaluate{1}.(fnames{nfnames}).F);
                data.p_net=cat(1,data.p_net,data.Zmodel.evaluate{1}.(fnames{nfnames}).p);
                data.dof_net=cat(1,data.dof_net,data.Zmodel.evaluate{1}.(fnames{nfnames}).dof);
            end
            data.statsname=data.Zmodel.evaluate{1}.(fnames{nfnames}).statsname;
        end
        data.dof_net=permute(reshape(data.dof_net,[],size(data.p_net,1),size(data.p_net,2)),[2,3,1]);
        data.dof=permute(reshape(data.dof,[],size(data.p,1),size(data.p,2)),[2,3,1]);
        data.displaymeasure=1;
        
        if 0,%isfield(CONN_x,'Setup')&&isfield(CONN_x.Setup,'normalized')&&~CONN_x.Setup.normalized,
            filename=spm_select(1,'\.img$|\.nii$',['Select background anatomical image'],{},fileparts(CONN_x.Setup.structural{1}{1}{1}));
        elseif isfield(CONN_gui,'refs')&&isfield(CONN_gui.refs,'canonical')&&isfield(CONN_gui.refs.canonical,'filename')&&~isempty(CONN_gui.refs.canonical.filename)
            filename=CONN_gui.refs.canonical.filename;
        else
            filename=fullfile(fileparts(which('spm')),'canonical','avg152T1.nii');
        end
        data.ref=spm_vol(filename);
        %color1=[1,1,1];
        %color2=.85*[1,1,1];
        color1=[1,1,1];
        color2=[1,1,1];
        hfig=figure('visible','on','renderer','opengl');
        cmap=linspace(0,1,256)'; cmap=repmat(1-cmap.^10.*(1-cmap),1,3);
        set(hfig,'units','norm','position',[.05,.2,.9,.7],'numbertitle','off','name','Network theory second-level results','color',color1,'colormap',cmap,'menubar','none');
        %uicontrol('style','frame','units','norm','position',[.0,.85,.5,.15],'backgroundcolor',color2,'foregroundcolor',color2);
        uicontrol('style','frame','units','norm','position',[.5,0,.5,1],'backgroundcolor',color2,'foregroundcolor',color2);
        %uicontrol('style','frame','units','norm','position',[.53,.02,.44,.
        %96],'backgroundcolor',color2);
        data.handles=[...
            uicontrol('style','text','units','norm','position',[.70,.82,.25,.04],'string','Analysis threshold (p-value)','foregroundcolor','k','backgroundcolor',color2,'fontweight','bold','horizontalalignment','left'),...
            uicontrol('style','edit','units','norm','position',[.70,.78,.05,.04],'string',num2str(data.thr),'foregroundcolor','k','backgroundcolor',color2,'callback',{@conn_displaynetwork,'thr'},'tooltipstring','Second-level analysis results false-positive threshold value'),...
            uicontrol('style','popupmenu','units','norm','position',[.76,.78,.10,.04],'string',{'p-uncorrected','p-FDR corrected'},'foregroundcolor','k','backgroundcolor',color2,'callback',{@conn_displaynetwork,'thrtype'},'value',data.thrtype,'tooltipstring','Second-level analysis resultsfalse-positive control type'),...
            uicontrol('style','popupmenu','units','norm','position',[.87,.78,.11,.04],'string',{'one-sided (positive)','one-sided (negative)','two-sided'},'foregroundcolor','k','backgroundcolor',color2,'callback',{@conn_displaynetwork,'side'},'value',data.side,'tooltipstring','Second-level analysis resultsdirectionality'),...
            uicontrol('style','text','units','norm','position',[.53,.82,.16,.04],'string','Analysis measure','foregroundcolor','k','backgroundcolor',color2,'fontweight','bold','horizontalalignment','left'),...
            uicontrol('style','listbox','units','norm','position',[.53,.07,.45,.15],'string',' ','max',2,'value',1,'fontname','monospaced','foregroundcolor','k','backgroundcolor','w','callback',{@conn_displaynetwork,'list1'}),...
            uicontrol('style','text','units','norm','position',[.53,.71,.45,.04],'string',sprintf('%-20s %6s  %6s  %4s  %8s  %8s','ROI','beta','T','dof','p-unc','p-FDR'),'foregroundcolor',.5*[1 1 1],'backgroundcolor',color2,'fontname','monospaced','fontweight','bold','horizontalalignment','left'),...
            uicontrol('style','listbox','units','norm','position',[.53,.23,.45,.48],'string',' ','max',2,'value',1,'fontname','monospaced','foregroundcolor','k','backgroundcolor',color2,'backgroundcolor','w','callback',@(varargin)conn_displaynetwork('list2')),...
            0,...%uicontrol('style','popupmenu','units','norm','position',[.08,.0,.10,.05],'string',{'axial (x-y)','coronal (x-z)','sagittal (y-z)'},'foregroundcolor','k','backgroundcolor',color1,'callback',{@conn_displaynetwork,'view'},'value',1),...
            0,...%%%uicontrol('style','slider','units','norm','position',[.38,.0,.10,.05],'foregroundcolor','w','backgroundcolor',color1,'callback',{@conn_displaynetwork,'slider1'},'value',.5),...
            0,...%uicontrol('style','slider','units','norm','position',[.23,.0,.10,.05],'foregroundcolor','w','backgroundcolor',color1,'callback',{@conn_displaynetwork,'slider2'},'value',.5),...
            uicontrol('style','popupmenu','units','norm','position',[.53,.90,.16,.04],'string',{'Network of all ROIs','Network of selected ROIs'},'value',data.displayreduced+1,'foregroundcolor','k','backgroundcolor',color2,'callback',{@conn_displaynetwork,'displayreduced'},'tooltipstring','Defines which nodes (ROIs) are included in the network'),...
            uicontrol('style','popupmenu','units','norm','position',[.53,.78,.16,.04],'string',{data.hnames{find(data.hnamesroi>0)}},'value',data.displaymeasure,'foregroundcolor','k','backgroundcolor',color2,'callback',{@conn_displaynetwork,'selectmeasure'},'tooltipstring','Choose a network measure to analyze'),...
            uicontrol('style','text','units','norm','position',[.70,.94,.24,.04],'string','Network edges (adjacency matrix threshold)','foregroundcolor','k','backgroundcolor',color2,'fontweight','bold','horizontalalignment','left'),...
            uicontrol('style','edit','units','norm','position',[.70,.90,.05,.04],'string',num2str(data.Zthr),'foregroundcolor','k','backgroundcolor',color2,'callback',{@conn_displaynetwork,'Zthr'},'tooltipstring','Edge-defining threshold value (two nodes are connected by an edge if their connectivity strength is above this value)'),...
            uicontrol('style','popupmenu','units','norm','position',[.76,.90,.10,.04],'string',{'correlation coefficient','z-score','cost'},'foregroundcolor','k','backgroundcolor',color2,'callback',{@conn_displaynetwork,'Zthrtype'},'value',data.Zthrtype,'tooltipstring','Edge-defining threshold units (use ''cost'' to produce networks with the same number of edges across all subjects)'),...
            uicontrol('style','text','units','norm','position',[.53,.94,.16,.04],'string','Network nodes','foregroundcolor','k','backgroundcolor',color2,'fontweight','bold','horizontalalignment','left'),...
            uicontrol('style','pushbutton','units','norm','position',[.76,.02,.15,.04],'string','Export data','callback',{@conn_displaynetwork,'export'},'tooltipstring','Exports statistic results and network/ROI measures'),...
            uicontrol('style','popupmenu','units','norm','position',[.87,.905,.11,.04],'string',{'one-sided (positive)','one-sided (negative)','two-sided'},'foregroundcolor','k','backgroundcolor',color2,'callback',{@conn_displaynetwork,'Zthrside'},'value',data.Zthrside,'tooltipstring','Connectivity strength based on r values (one-sided) or abs(r) values (two-sided)'),...
            uicontrol('style','pushbutton','units','norm','position',[.60,.02,.15,.04],'string','Export Adjacency matrices','callback',{@conn_displaynetwork,'exportadj'},'tooltipstring','Exports adjacency matrix for each subject'),...
            uicontrol('style','text','units','norm','position',[0,0,.5,.03],'string','','foregroundcolor','k','backgroundcolor','w')];
        hc1=uicontextmenu('parent',hfig);
        uimenu(hc1,'Label','Export table','callback',@(varargin)conn_exportlist(data.handles(6)));
        set(data.handles(6),'uicontextmenu',hc1);
        hc1=uicontextmenu('parent',hfig);
        uimenu(hc1,'Label','Export table','callback',@(varargin)conn_exportlist(data.handles(8),'',get(data.handles(7),'string')));
        set(data.handles(8),'uicontextmenu',hc1);
        %uicontrol('style','text','units','norm','position',[.03,.01,.05,.04],'string','view:  ','foregroundcolor','k','backgroundcolor',color1,'fontweight','bold','horizontalalignment','right');
        %uicontrol('style','text','units','norm','position',[.35,.01,.02,.04],'string','q','fontname','symbol','horizontalalignment','right','foregroundcolor','k','backgroundcolor',color1,'fontweight','bold','horizontalalignment','right');
        %uicontrol('style','text','units','norm','position',[.20,.01,.02,.04],'string','z','horizontalalignment','right','foregroundcolor','k','backgroundcolor',color1,'fontweight','bold','horizontalalignment','right');
%         rcircle=[sin(linspace(0,2*pi,64)'),cos(linspace(0,2*pi,64))']*diag([5,5]);
%         h=subplot(121);set(h,'units','norm','position',[.13,.845,.25,.05],'xtick',[],'ytick',[],'box','on');axis equal;set(h,'xlim',[-10,200]);
%         h=patch(100+rcircle(:,1),0+rcircle(:,2),'w');set(h,'edgecolor','none','facecolor','r'); 
%         h=patch(150+rcircle(:,1),0+rcircle(:,2),'w');set(h,'edgecolor','none','facecolor','b'); 
%         h=text(0,0,'ROI-ROI connectivity:','horizontalalignment','left','fontsize',8);
%         h=text(100+10,0,'Positive','horizontalalignment','left','fontsize',8);
%         h=text(150+10,0,'Negative','horizontalalignment','left','fontsize',8);
    case 'thr',
        hfig=gcbf;
        data=get(hfig,'userdata');
        str=get(data.handles(2),'string');
        if ~isempty(str2num(str)),
            data.thr=max(0,str2num(str));
        end
        data.visible='on';
    case 'thrtype',
        hfig=gcbf;
        data=get(hfig,'userdata');
        value=get(data.handles(3),'value');
        data.thrtype=value;
        data.visible='on';
    case 'side',
        hfig=gcbf;
        data=get(hfig,'userdata');
        value=get(data.handles(4),'value');
        data.side=value;
        data.visible='on';
    case 'selectmeasure',
        hfig=gcbf;
        data=get(hfig,'userdata');
        value=get(data.handles(13),'value');
        data.displaymeasure=value;
        data.visible='on';
    case {'zthr','zthrtype','zthrside','displayreduced'},
        hfig=gcbf;
        data=get(hfig,'userdata');
        
        switch(lower(option)),
            case 'zthr',
                value=get(data.handles(15),'string');
                if ~isempty(value), 
                    value=str2num(value); 
                    if ~isempty(value), data.Zthr=value; end
                else
                    set(data.handles(15),'string',num2str(data.Zthr)); 
                    if ~isfield(data,'Zthrside')||data.Zthrside==1, 
                        conn_network({data.Z,data.names},data.displaytheserois,1,data.Zthrtype-1,[]);
                    elseif data.Zthrside==2, 
                        conn_network({-data.Z,data.names},data.displaytheserois,1,data.Zthrtype-1,[]);
                    else
                        conn_network({abs(data.Z),data.names},data.displaytheserois,1,data.Zthrtype-1,[]);
                    end
                    return;
                end
            case 'zthrtype',
                value=get(data.handles(16),'value');
                if ~isempty(value), data.Zthrtype=value; end
            case 'zthrside',
                value=get(data.handles(19),'value');
                if ~isempty(value), data.Zthrside=value; end
            case 'displayreduced',
                value=get(data.handles(12),'value');
                data.displayreduced=value-1;
                switch(data.displayreduced),
                    case 0,
                        data.displaytheserois=1:length(data.names);
                    case 1,
                        answ=listdlg('Promptstring','Select ROIs','selectionmode','multiple','liststring',data.names,'initialvalue',data.displaytheserois);
                        if ~isempty(answ), data.displaytheserois=answ; end
                end
                data.source=0;%data.source(data.source==0 || data.source<=length(data.displaytheserois));if isempty(data.source),data.soruce=0;end
                
        end
        
        data.visible='on';
        Zthr=data.Zthr; if data.Zthrtype==1,Zthr=atanh(Zthr);end
        if ~isfield(data,'Zthrside')||data.Zthrside==1,
            for n1=1:size(data.Z,4), data.Zresults(n1)=conn_network({data.Z(:,:,:,n1),data.names},data.displaytheserois,1,data.Zthrtype-1,Zthr); end
        elseif data.Zthrside==2,
            for n1=1:size(data.Z,4), data.Zresults(n1)=conn_network({-data.Z(:,:,:,n1),data.names},data.displaytheserois,1,data.Zthrtype-1,Zthr); end
        else
            for n1=1:size(data.Z,4), data.Zresults(n1)=conn_network({abs(data.Z(:,:,:,n1)),data.names},data.displaytheserois,1,data.Zthrtype-1,Zthr); end
        end
%         if ~isfield(data,'Zthrside')||data.Zthrside==1,
%             data.Zresults=conn_network({data.Z,data.names},data.displaytheserois,1,data.Zthrtype-1,Zthr);
%         else
%             data.Zresults=conn_network({abs(data.Z),data.names},data.displaytheserois,1,data.Zthrtype-1,Zthr);
%         end
        data.Zmodel=conn_network_results({data.Zresults,data.Zmodel});
        fnames=fieldnames(data.Zmodel.evaluate{1});
        data.hnames=fnames;
        data.hnamesroi=zeros(1,numel(fnames));
        data.h=[];data.F=[];data.p=[];data.dof=[];data.h_net=[];data.F_net=[];data.p_net=[];data.dof_net=[];
        for nfnames=1:numel(fnames),
            if numel(fnames{nfnames})>4&&strcmp(fnames{nfnames}(end-3:end),'_roi'),
                data.h=cat(1,data.h,data.Zmodel.evaluate{1}.(fnames{nfnames}).h);
                data.F=cat(1,data.F,data.Zmodel.evaluate{1}.(fnames{nfnames}).F);
                data.p=cat(1,data.p,data.Zmodel.evaluate{1}.(fnames{nfnames}).p);
                data.dof=cat(1,data.dof,data.Zmodel.evaluate{1}.(fnames{nfnames}).dof);
                data.hnames{nfnames}=['',fnames{nfnames}(1:end-4)];
                %data.hnames{nfnames}=['Analysis of ',fnames{nfnames}(1:end-4)];
                data.hnamesroi(nfnames)=1;
            else
                data.h_net=cat(1,data.h_net,data.Zmodel.evaluate{1}.(fnames{nfnames}).h);
                data.F_net=cat(1,data.F_net,data.Zmodel.evaluate{1}.(fnames{nfnames}).F);
                data.p_net=cat(1,data.p_net,data.Zmodel.evaluate{1}.(fnames{nfnames}).p);
                data.dof_net=cat(1,data.dof_net,data.Zmodel.evaluate{1}.(fnames{nfnames}).dof);
            end
            data.statsname=data.Zmodel.evaluate{1}.(fnames{nfnames}).statsname;
        end
        data.dof_net=permute(reshape(data.dof_net,[],size(data.p_net,1),size(data.p_net,2)),[2,3,1]);
        data.dof=permute(reshape(data.dof,[],size(data.p,1),size(data.p,2)),[2,3,1]);
        
    case 'list1',
        hfig=gcbf;
        data=get(hfig,'userdata');
        value=get(data.handles(6),'value');
        set(data.plotsadd2,'linewidth',1,'edgecolor','none');
        if all(value>0)&&all(value<=numel(data.plotsadd2)), 
            set(data.plotsadd2(value),'linewidth',2,'edgecolor','k'); 
            h=[];for n1=1:numel(value),h=[h,find(data.list2==value(n1))];end
            set(data.handles(8),'value',1+h);
        end
        return;
        if all(value>0&value<=length(data.displaytheserois)), 
            data.source=value;
            data.bgz=mean(data.z(data.displaytheserois(data.source)));
            %set(data.handles(11),'value',max(0,min(1,data.bgz/200+.5)));
            data.bgimage=spm_get_data(data.ref,pinv(data.ref.mat)*[data.proj(:,1:3)*[data.bgx(:),data.bgy(:),data.bgz+zeros(prod(size(data.bgx)),1)]';ones(1,prod(size(data.bgx)))]);
            data.bgimage(isnan(data.bgimage))=0;
            %set(data.refaxes,'cdata',convn(convn(reshape(data.bgimage,size(data.bgx)),conn_hanning(5),'same'),conn_hanning(5)','same'));
        else
            data.source=0; 
            data.bgz=0;
            %set(data.handles(11),'value',max(0,min(1,data.bgz/200+.5)));
            data.bgimage=spm_get_data(data.ref,pinv(data.ref.mat)*[data.proj(:,1:3)*[data.bgx(:),data.bgy(:),data.bgz+zeros(prod(size(data.bgx)),1)]';ones(1,prod(size(data.bgx)))]);
            data.bgimage(isnan(data.bgimage))=0;
            %set(data.refaxes,'cdata',convn(convn(reshape(data.bgimage,size(data.bgx)),conn_hanning(5),'same'),conn_hanning(5)','same'));
        end
        data.visible='on';
    case 'list2',
        if nargin>=2&&~isempty(varargin{1}), hfig=varargin{1}; 
        else hfig=gcbf;
        end
        data=get(hfig,'userdata');
        value=get(data.handles(8),'value');
%         if value==1||value-1>numel(data.list2), data.source=0;
%         else data.source=data.list2(value-1); end
        set(data.plotsadd2,'linewidth',1,'edgecolor','none');
        if all(value>1)&&all(value-1<=numel(data.plotsadd2)), 
            set(data.plotsadd2(data.list2(value-1)),'linewidth',2,'edgecolor','k'); 
            tstr=cellstr(get(data.handles(6),'string')); set(data.handles(21),'string',[tstr{data.list2(value-1)}]);
            set(data.handles(6),'value',data.list2(value-1));
        end
        return;
        data.visible='on';
    case {'view','view-axial','view-coronal','view-sagittal'}
        hfig=gcbf;
        data=get(hfig,'userdata');
        switch(lower(option))
            case 'view-axial', data.view=1;
            case 'view-coronal', data.view=2;
            case 'view-sagittal', data.view=3;
        end
        %data.view=get(data.handles(9),'value');
        data.proj=[];data.x=[];data.y=[];data.z=[];
        %set(data.handles(10),'value',.5);
        %set(data.handles(11),'value',.5);
        data.bgz=0;
    case {'displayefffectsize-on','displayefffectsize-off'}
        hfig=gcbf;
        data=get(hfig,'userdata');
        switch(lower(option))
            case 'displayefffectsize-on', data.displayeffectsize=1;
            case 'displayefffectsize-off', data.displayeffectsize=0;
        end
        data.proj=[];data.x=[];data.y=[];data.z=[];
        data.bgz=0;
    case 'display3d'
        hfig=gcbf;
        data=get(hfig,'userdata');
        data.view=1;
        data.proj=[];data.x=[];data.y=[];data.z=[];
        data.bgz=0;
        data.display3d=1;
    case 'changebackground',
        hfig=gcbf;
        data=get(hfig,'userdata');
        filename=spm_select(1,'\.img$|\.nii$',['Select background anatomical image'],{},fileparts(data.ref.fname));
        data.ref=spm_vol(filename);
        data.view=[];data.proj=[];data.x=[];data.y=[];data.z=[];
        data.bgz=0;
    case {'labelsoff','labelson'}
        hfig=gcbf;
        data=get(hfig,'userdata');
        data.displaylabels=strcmpi(option,'labelson');
    case 'slider1',
        hfig=gcbf;
        data=get(hfig,'userdata');
        value=get(data.handles(10),'value');
        ang=(value-.5)*pi;
        switch(data.view),
            case 1, data.proj=[1,0,0;0,cos(ang),-sin(ang);0,sin(ang),cos(ang)];
            case 2, data.proj=[1,0,0;0,sin(ang),-cos(ang);0,cos(ang),sin(ang)];
            case 3, data.proj=[0,sin(ang),cos(ang);1,0,0;0,cos(ang),-sin(ang)];
        end
        data.x=[];data.y=[];data.z=[];
    case 'slider2',
        hfig=gcbf;
        data=get(hfig,'userdata');
        value=get(data.handles(11),'value');
        data.bgz=(value-.5)*200;
        data.bgimage=spm_get_data(data.ref,pinv(data.ref.mat)*[data.proj(:,1:3)*[data.bgx(:),data.bgy(:),data.bgz+zeros(prod(size(data.bgx)),1)]';ones(1,prod(size(data.bgx)))]);
        data.bgimage(isnan(data.bgimage))=0;
        set(data.refaxes,'cdata',convn(convn(reshape(data.bgimage,size(data.bgx)),conn_hanning(5),'same'),conn_hanning(5)','same'));
        set(hfig,'userdata',data);
        return;
    case 'export',
        hfig=gcbf;
        data=get(hfig,'userdata');
        [filename,pathname]=uiputfile({'*.csv','comma-separated file (*.csv)'},'Save as');
        if isequal(filename,0), return; end
        Zthr=data.Zthr; if data.Zthrtype==1,Zthr=atanh(Zthr);end        
        if ~isfield(data,'Zthrside')||data.Zthrside==1,
            conn_network({data.Z,data.names,fullfile(pathname,filename)},data.displaytheserois,1,data.Zthrtype-1,Zthr);
        elseif data.Zthrside==2,
            conn_network({-data.Z,data.names,fullfile(pathname,filename)},data.displaytheserois,1,data.Zthrtype-1,Zthr);
        else
            conn_network({abs(data.Z),data.names,fullfile(pathname,filename)},data.displaytheserois,1,data.Zthrtype-1,Zthr);
        end            
        [nill,filename_name,filename_ext]=fileparts(filename);
        filename_out=fullfile(pathname,[filename_name,'.stats.mat']);
        stats=data.Zmodel.evaluate{1};
        stats.ROInames=data.names;
        conn_savematfile(filename_out,'stats');
        fprintf(1,'Statistics saved as: %s\n',filename_out);
        return;
    case 'exportadj',
        hfig=gcbf;
        data=get(hfig,'userdata');
        [filename,pathname]=uiputfile({'*.mat','Matlab file (*.mat)'; '*.dl','UCINET file (*.dl)'},'Save as');
        if isequal(filename,0), return; end
        Zthr=data.Zthr; if data.Zthrtype==1,Zthr=atanh(Zthr);end        
        if ~isfield(data,'Zthrside')||data.Zthrside==1,
            conn_network({data.Z,data.names,fullfile(pathname,filename)},data.displaytheserois,1,data.Zthrtype-1,Zthr);
        elseif data.Zthrside==2,
            conn_network({-data.Z,data.names,fullfile(pathname,filename)},data.displaytheserois,1,data.Zthrtype-1,Zthr);
        else
            conn_network({abs(data.Z),data.names,fullfile(pathname,filename)},data.displaytheserois,1,data.Zthrtype-1,Zthr);
        end            
        return;
    case 'disptype',
        hfig=gcf;
        data=get(hfig,'userdata');
        data.disptype=mod(1+data.disptype,2);
end

% selects optimal view
if isempty(data.view),
    a=std(data.xyz2(:,1:3),1,1);
    [nill,idx]=sort(a);
    if idx(end-1)==1, idx([end-1,end])=idx([end,end-1]); end
    if idx(end)==1&&idx(end-1)==2,data.view=1;
    elseif idx(end)==1&&idx(end-1)==3,data.view=2;
    else data.view=3; end
    %set(data.handles(9),'value',data.view);
end
% projector associated with view
if isempty(data.proj),
    switch(data.view),
        case 1, data.proj=[1,0,0;0,1,0;0,0,1];
        case 2, data.proj=[1,0,0;0,0,1;0,1,0];
        case 3, data.proj=[0,0,1;1,0,0;0,1,0];
    end
end
% projects coordinates and background image
scale=1;
scalecircle=1;
if isempty(data.x)||isempty(data.y),
    data.x=data.xyz2*data.proj(:,1);data.y=data.xyz2*data.proj(:,2);data.z=data.xyz2*data.proj(:,3);
    lim=[1,1,1;data.ref.dim];refminmax=sort([lim((dec2bin(0:7)-'0'+1)+repmat([0,2,4],[8,1])),ones(8,1)]*data.ref.mat(1:3,:)'*data.proj(:,1:2));
    [data.bgx,data.bgy]=meshgrid(refminmax(1,1):scale:refminmax(end,1),refminmax(1,2):scale:refminmax(end,2));
    data.bgimage=spm_get_data(data.ref,pinv(data.ref.mat)*[data.proj(:,1:3)*[data.bgx(:),data.bgy(:),data.bgz+zeros(prod(size(data.bgx)),1)]';ones(1,prod(size(data.bgx)))]);
    data.bgimage(isnan(data.bgimage))=0;
end
N=length(data.displaytheserois);%length(data.names);
N2=length(data.names2);
switch(data.display),
    case 'connectivity',
        cmap=jet(64);cmap=cmap(8:56,:);
        % computes stat threshold
        if ~isequal(data.statsname,'T'), data.side=1; end
        if isequal(data.statsname,'T'), set(data.handles(4),'enable','on'); else set(data.handles(4),'value',3,'enable','off'); end
        switch(data.side),
            case 1,p=data.p(data.displaymeasure,:);p_net=data.p_net(data.displaymeasure);
            case 2,p=1-data.p(data.displaymeasure,:);p_net=1-data.p_net(data.displaymeasure);
            case 3,p=2*min(data.p(data.displaymeasure,:),1-data.p(data.displaymeasure,:));p_net=2*min(data.p_net(data.displaymeasure),1-data.p_net(data.displaymeasure));
        end
        %if data.displayreduced, p=p(:,1:N); end
        %p=p(intersect(1:N,data.displaytheserois),data.displaytheserois);
        %p(setdiff(1:N,data.displaytheserois))=nan;
        P=reshape(conn_fdr(p(:)),size(p));
        data.P=P;
        switch(data.thrtype),
            case 1, z=(p<=data.thr).*(.1+.8*p/max(eps,data.thr));
            case 2, z=(P<=data.thr).*(.1+.8*P/max(eps,data.thr));
        end
        z(isnan(p))=nan;
        
        % text lists
        txt1={};for n1=1:N,%length(data.displaytheserois),
            if n1<=N&&z(n1)>0,
                %txt1{end+1}=(sprintf('%-5s= %-40s',['(',num2str(n1),')'],data.names{data.displaytheserois(n1)}));
                txt1{end+1}=[data.names{data.displaytheserois(n1)},...
                             '    ','x,y,z = (',num2str(data.xyz2(data.displaytheserois(n1),1),'%1.0f'),',',num2str(data.xyz2(data.displaytheserois(n1),2),'%1.0f'),',',num2str(data.xyz2(data.displaytheserois(n1),3),'%1.0f'),') mm'];
            end
        end;txt1{end+1}=' ';
        txt1=strvcat(txt1{:});
        %tp=[];sort2=[];txt2={};for n1=1:N,for n2=[1:n1-1,n1+1:N2+data.displayreduced*(N-N2)],
        tp=[];sort2=[];txt2={};
        if 1,
            txt2{end+1}=(sprintf('<HTML><pre><b>%-20s %6.2f  %6.2f  %4d  %8.6f</b></pre></HTML>',['network'],data.h_net(data.displaymeasure),data.F_net(data.displaymeasure),data.dof_net(data.displaymeasure,end),p_net));
        end
        for n1=1:N,
            if n1<=N&&z(n1)>0,
                temp=data.namesreduced{data.displaytheserois(n1)};
                if numel(temp)>20, temp=[temp(1:17),'...']; end
                txt2{end+1}=(sprintf('%-20s %6.2f  %6.2f  %4d  %8.6f  %8.6f',temp,data.h(data.displaymeasure,n1),data.F(data.displaymeasure,n1),data.dof(data.displaymeasure,n1,end),p(n1),P(n1))); 
                tp=cat(1,tp,P(n1)-1e-10*abs(data.F(data.displaymeasure,n1))); 
                sort2=cat(1,sort2,[n1]);
            end
        end
        if size(data.dof,3)>1
            set(data.handles(7),'string',sprintf('%-20s %6s  %6s  %4s  %8s  %8s','ROI','beta',[data.statsname,'(',num2str(data.dof(1)),')'],'dof','p-unc','p-FDR'));
        else
            set(data.handles(7),'string',sprintf('%-20s %6s  %6s  %4s  %8s  %8s','ROI','beta',data.statsname,'dof','p-unc','p-FDR'));
        end
        [nill,idxsort]=sort(tp); data.list2=idxsort;%sort2(idxsort,:); 
        txt2=strvcat(txt2{1},txt2{1+idxsort},' ');
        
        set(data.handles(6),'string',txt1,'value',max(1,min(size(txt1,1), get(data.handles(6),'value'))),'listboxtop',1); 
        set(data.handles(8),'string',txt2,'value',max(1,min(size(txt2,1), get(data.handles(8),'value'))),'listboxtop',1);
        
        % plots
        figure(hfig);set(hfig,'pointer','watch');drawnow;
        tcmap=linspace(0,1,256)'; tcmap=repmat(1-tcmap.*(1-tcmap),1,3);
        if ~data.disptype, set(hfig,'colormap',tcmap); else set(hfig,'colormap',gray(128)); end
        rcircle=[sin(linspace(0,2*pi,64)'),cos(linspace(0,2*pi,64))']*diag([5,5]*scalecircle);
        rtriang=[-1-j*.5;0;-1+j*.5];
        h=findobj(hfig,'tag','conn_displaynetwork_plot');
        if ~isempty(h),delete(h); end
        hax=axes('units','norm','position',[.03,.01,.45,.98],'parent',hfig);
        lim=[1,1,1;data.ref.dim];refminmax=sort([lim((dec2bin(0:7)-'0'+1)+repmat([0,2,4],[8,1])),ones(8,1)]*data.ref.mat(1:3,:)'*data.proj(:,1:2));
        if 0,%~data.disptype, data.refaxes=imagesc(refminmax(1,1):scale:refminmax(end,1),refminmax(1,2):scale:refminmax(end,2),convn(convn(reshape(data.bgimage,size(data.bgx)),conn_hanning(5),'same'),conn_hanning(5)','same')>3);hold on; 
        else data.refaxes=imagesc(refminmax(1,1):scale:refminmax(end,1),refminmax(1,2):scale:refminmax(end,2),convn(convn(reshape(data.bgimage,size(data.bgx)),conn_hanning(5),'same'),conn_hanning(5)','same'),'parent',hax);hold(hax,'on'); end
        hc1=uicontextmenu('parent',hfig);
        uimenu(hc1,'Label','Change background anatomical image','callback',{@conn_displaynetwork,'changebackground'});
        uimenu(hc1,'Label','Switch display type','callback',{@conn_displaynetwork,'disptype'});
        uimenu(hc1,'Label','Labels on','callback',{@conn_displaynetwork,'labelson'});
        uimenu(hc1,'Label','Labels off','callback',{@conn_displaynetwork,'labelsoff'});
        uimenu(hc1,'Label','View axial (x-y)','callback',{@conn_displaynetwork,'view-axial'});
        uimenu(hc1,'Label','View coronal (x-z)','callback',{@conn_displaynetwork,'view-coronal'});
        uimenu(hc1,'Label','View sagittal (y-z)','callback',{@conn_displaynetwork,'view-sagittal'});
        uimenu(hc1,'Label','Circle-sizes represent T-values','callback',{@conn_displaynetwork,'displayefffectsize-off'});
        uimenu(hc1,'Label','Circle-sizes represent beta-values','callback',{@conn_displaynetwork,'displayefffectsize-on'});
        uimenu(hc1,'Label','Menubar on','callback',['set(gcbf,''menubar'',''figure'');']);
        uimenu(hc1,'Label','Menubar off','callback',['set(gcbf,''menubar'',''none'');']);
        uimenu(hc1,'Label','Display 3d view','callback',{@conn_displaynetwork,'display3d'});
        set(data.refaxes,'uicontextmenu',hc1);
        set(hfig,'uicontextmenu',hc1);
        h=hax;data.buttondown=struct('h1',h);set(h,'tag','conn_displaynetwork_plot');
        
        EPS=1e-0;
        data.plotsadd2=[];
        idxtext=[];
        hold(hax,'on');
        if isfield(data,'displayeffectsize')&&data.displayeffectsize
            K=data.h(data.displaymeasure,:)/max(eps,max(abs(data.h(data.displaymeasure,:))));
        else
            K=data.F(data.displaymeasure,:)/max(eps,max(abs(data.F(data.displaymeasure,:))));
        end
        n1a=1;
        Zthrdisplay=mean(data.Zresults(1).Z_thr,3);
        Zthrdisplay_show=.25;
        for n1=1:N,%1:length(data.displaytheserois),%size(z,2),%N2,
%             n1=data.displaytheserois(na1);
            if (n1<=numel(z)&&z(n1)>0),
                k=K(n1);
                if isnan(k),k=0; end
                if ~data.disptype, h=patch(data.x(data.displaytheserois(n1))+rcircle(:,1)*abs(k),data.y(data.displaytheserois(n1))+rcircle(:,2)*abs(k),max(1,data.z(data.displaytheserois(n1)))*EPS+200+zeros(size(rcircle,1),1),'w','parent',hax);
                else h=patch(data.x(data.displaytheserois(n1))+rcircle(:,1)*max(.5,1+1e-3*data.z(data.displaytheserois(n1))),data.y(data.displaytheserois(n1))+rcircle(:,2)*max(.5,1+1e-3*data.z(data.displaytheserois(n1))),max(1,data.z(data.displaytheserois(n1)))*EPS+200+zeros(size(rcircle,1),1),'w','parent',hax); end
                if any(n1a==data.source), set(h,'edgecolor','k','linewidth',2); else set(h,'edgecolor','none'); end
                if ~data.disptype, set(h,'facecolor',cmap(round(1+48*(1+sign(k))/2),:));
                else set(h,'facecolor',cmap(round(1+48*(1+k)/2),:)); end
                data.plotsadd2(end+1)=h;
                %hold(hax,'on');
                set(h,'buttondownfcn',@conn_displaynetwork_menubuttondownfcn,'userdata',data);
                n1a=n1a+1;
            end
            if ~data.disptype, idx=find(Zthrdisplay(n1,:)>Zthrdisplay_show);
            else idx=find(Zthrdisplay(n1,:)>Zthrdisplay_show); end
            if ~isempty(idx),
                x=data.x([data.displaytheserois(n1)+zeros(1,numel(idx));data.displaytheserois(idx)]);
                y=data.y([data.displaytheserois(n1)+zeros(1,numel(idx));data.displaytheserois(idx)]);
                if ~data.disptype,h=plot3(x,y,max(0,repmat(Zthrdisplay(n1,idx),[2,1]))*EPS+.5,'parent',hax); %hold on;
                    for n2a=1:numel(idx),set(h(n2a),'color',.5+.5*(1-Zthrdisplay(n1,idx(n2a)))*[1,1,1],'visible',data.visible); end
                else h=plot3(x,y,max(0,repmat(Zthrdisplay(n1,idx),[2,1]))*EPS+.5,'parent',hax); %hold on;
                    for n2a=1:numel(idx),set(h(n2a),'color',(1-Zthrdisplay(n1,idx(n2a)))*[1,1,1],'visible',data.visible); end
                end
            end
            if (n1<=N&&z(n1)>0), 
                idxtext=[idxtext;n1];
            end
        end
        if isfield(data,'displaylabels')&&data.displaylabels
            if ~data.disptype,  
                h=text(data.x(data.displaytheserois(idxtext)),data.y(data.displaytheserois(idxtext))+scalecircle*5*abs(K(idxtext))'+3,max(1,data.z(data.displaytheserois(idxtext)))*EPS+202,{data.namesreduced{data.displaytheserois(idxtext)}},'parent',hax);
                set(h,'fontsize',8,'color','k','horizontalalignment','center','interpreter','none','fontweight','normal','backgroundcolor','none');
            else h=text(data.x(data.displaytheserois(idxtext)),data.y(data.displaytheserois(idxtext)),max(1,data.z(data.displaytheserois(idxtext)))*EPS+202,num2str(idxtext),'parent',hax);
                set(h,'fontsize',9,'color','w','horizontalalignment','center','interpreter','none','fontweight','bold','backgroundcolor','none'); 
            end
            set(h,'buttondownfcn',@conn_displaynetwork_menubuttondownfcn,'userdata',data);
        end
        %data.plotsadd2=data.plotsadd2(idxsort);
        hold(hax,'off');
        set(hax,'ydir','normal');
        axis(hax,'equal','off'); %axis equal; axis off;
        if data.display3d, data.display3d=0; datadisplay3d=1;else datadisplay3d=0; end
        set(hfig,'userdata',data);
        set(hfig,'pointer','arrow');
        if datadisplay3d
            tidxtext=1:N;
            c=mat2cell(cmap(round(1+48*(1+sign(K(tidxtext)))/2),:),ones(numel(tidxtext),1),3);
            for n1=1:numel(data.source),idxc=find(tidxtext==data.source(n1));c(idxc,:)=repmat({[.25,.25,.25]},[numel(idxc),1]); end
            conn_mesh_display('','',[],...
                struct('sph_names',{data.namesreduced(data.displaytheserois(tidxtext))},'sph_xyz',[data.x(data.displaytheserois(tidxtext)),data.y(data.displaytheserois(tidxtext)),data.z(data.displaytheserois(tidxtext))],...
                'sph_r',max(1,6*abs(K(tidxtext))).*(abs(K(tidxtext))>0).*ismember(tidxtext,idxtext),... %6*ones(numel(idxtext),1),...
                'sph_c',{c}), ...
                (1*Zthrdisplay(tidxtext,tidxtext)).*(Zthrdisplay(tidxtext,tidxtext)>Zthrdisplay_show), ...
                [], .25, [0,-.01,1],{[.5,.5,.5],[.5,.5,.5]});
        end
end
end

function conn_displaynetwork_menubuttondownfcn(varargin)
data=get(gcbo,'userdata');
xyz=get(data.buttondown.h1,'currentpoint'); 
x=xyz(1,1);y=xyz(1,2);z=0;
[nill,idx]=min(sqrt(abs(data.x(data.displaytheserois)-x).^2+abs(data.y(data.displaytheserois)-y).^2)+1e-10*abs(data.z(data.displaytheserois)-z).^2);
idx2a=find(data.list2(:,1)==data.displaytheserois(idx));
if ~isempty(idx2a),
    set(data.handles(8),'value',1+idx2a,'listboxtop',min(1+idx2a));
    conn_displaynetwork('list2',gcbf);
end
end
