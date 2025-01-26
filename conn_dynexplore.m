function ok=conn_dynexplore
global CONN_x CONN_gui CONN_h;

if isempty(CONN_gui)||~isfield(CONN_gui,'font_offset'), conn_font_init; end
filepathresults=fullfile(CONN_x.folders.firstlevel_dyn,CONN_x.dynAnalyses(CONN_x.dynAnalysis).name);
[B,H,B0,H0,IDX_subject,IDX_session]=deal([]); [ROInames,names]=deal({}); conn_loadmatfile(fullfile(filepathresults,'dyn_Base.mat'),'B','H','B0','H0','IDX_subject','IDX_session','ROInames','names');
if ~exist('ROInames','var'), conn_msgbox('Sorry, this option is not available until the first-level dynamic FC analyses have been re-run','',2); return; end
%names=regexprep(names,'Dynamic factor','Circuit');
names=arrayfun(@(n)sprintf('Circuit_%d',n),1:size(B,3),'uni',0);
B0=permute(B0,[3 1 2]);
B=permute(B,[3 1 2]);
IDX_scan=ones(size(IDX_session));
temp=find(IDX_subject(max(1,0:numel(IDX_subject)-1))~=IDX_subject|IDX_session(max(1,0:numel(IDX_session)-1))~=IDX_session);
IDX_scan(temp)=1-diff([1;temp]);
IDX_scan=cumsum(IDX_scan);
ROInamesDisplay=regexprep(ROInames,'^atlas\.','');
ROInamesDisplay2=regexprep(ROInames,'^atlas\.|\([^\(\)]+\)','');
TR=conn_get_rt;
KURx=zeros(1,size(B,1));
SKWx=zeros(1,size(B,1));
for n=1:size(B,1),
    b=shiftdim(B(n,:,:),1);
    %b(1:size(b,1)+1:end)=nan;
    mask=~isnan(b);
    vfact=sqrt(mean(b(mask).^2));
    b(mask)=b(mask)./vfact;
    B(n,:,:)=shiftdim(b,-1); H(:,n)=H(:,n)*vfact;
    mb=mean(b(mask));
    sb=std(b(mask),1);
    KURx(n)=mean(((b(mask)-mb)/sb).^4);
    SKWx(n)=mean(((b(mask)-mb)/sb).^3);
end
nsubs=max(IDX_subject);
VARt=zeros(1,size(B,1));
ZCt=zeros(1,size(B,1));
validc=0;
for n1=1:nsubs
    ts=unique(IDX_session(IDX_subject==n1));
    b=H(IDX_subject==n1,:);
    if ~isempty(ts)
        mb=mean(b,1);%./sum(b~=0,1);
        VARt(1,:)=VARt(1,:)+std(b,1,1);
        for n2=1:numel(ts)
            b=H(IDX_subject==n1&IDX_session==ts(n2),:);
            b=conn_bsxfun(@minus,b,mb);
            ZCt(1,:)=ZCt(1,:)+((mean((b(1:end-1,:)<=0&b(2:end,:)>0)|(b(1:end-1,:)>=0&b(2:end,:)<0),1))/TR(min(numel(TR),n1)))/numel(ts)/2;
        end
        validc=validc+1;
    end
end
VARt=VARt/validc;
ZCt=ZCt/validc;
varB=abs(B); %sqrt(max(0,conn_bsxfun(@times,mean(H.^2,1)',abs(B).^2)));
[Bhf_values,Bhf_order]=sort(varB,1,'descend');
[Bhfc_values,Bhfc_order]=sort(sqrt(mean(varB.^2,3)),1,'descend');
nfacshown=1:numel(names);
nfacselected=nfacshown;
try
    load(fullfile(fileparts(which('conn')),'connROIorder.mat'),'ROIconfiguration');
    [ok,idx1]=ismember(ROInames,ROIconfiguration.names2);
    if mean(ok)>.9
        ok=find(ok);
        z=angle(ROIconfiguration.xy2(idx1(ok),:)*[1;1i]).';
        [nill,idx2]=sort(z);
        idxorder=[ok(idx2), setdiff(1:numel(ROInames), ok)];
    elseif numel(ok)<=500 && numel(ok)>50
        %x=conn_bsxfun(@times,cat(1,.0*B0,B),sqrt(mean(cat(2,H0,H).^2,1))');
        x=tanh(B); %double(abs(B)>2);
        x=reshape(x,[],size(B,3));
        nx=sum(x.^2,1);
        y=sqrt(max(0,conn_bsxfun(@plus,nx,nx')-2*(x'*x)));
        y=y(tril(ones(size(y,1)),-1)>0)';
        z=conn_statslinkage(y,'av');
        [nill,nill,idxorder]=conn_statsdendrogram(z,0);
    else
        idxorder=1:numel(ROInames);
    end
catch
    idxorder=1:numel(ROInames);
end
Bhf_max={}; for n2=1:size(varB,1), 
    [idx1,idx2]=find(shiftdim(Bhf_order(1,idxorder,idxorder),1)==n2); 
    Bhf_max{n2}=idx1+size(Bhf_order,2)*(idx2-1); 
    Bhf_maxval{n2}=Bhf_values(1,idxorder(idx1)+size(Bhf_order,2)*(idxorder(idx2)-1)); 
    maskout=Bhf_maxval{n2}>=prctile(Bhf_maxval{n2},90); Bhf_max{n2}=Bhf_max{n2}(maskout); Bhf_maxval{n2}=Bhf_maxval{n2}(maskout); 
end

B0(:,1:size(B0,2)+1:end)=nan;
Hscale=max(abs(H(:)));
idxlock=[1,2,0]; 

boffset=[0 0 0 0];
conn_menu('frame',boffset+[.02,.42,.46,.47],'');%'Component loadings');
%conn_menu('frame',boffset+[.44 .46 .16 .46],'');
poslist=boffset+[.07 .42 .40 .36];
ht2=conn_menu('listbox',poslist,'',names,'<HTML>Select dynamic circuits for display</HTML>',@(varargin)conn_dynexplore_update([0 0 1 0]));
for n=6:20, set(ht2,'string',repmat(' ',1,5*n),'fontname','monospaced','fontsize',8+CONN_gui.font_offset); if get(ht2,'extent')*[0 0 1 0]'>poslist(3), break; end; end
ht3fieldsize=sprintf('%d',n-1);
temp=boffset+[.07 .78 .38/5 .04];
ht2a(1)=conn_menu('pushbutton',temp+[0*temp(3) 0 0 0],'','circuit','Select circuit(s) for display',@(varargin)conn_dynexplore_update([0 0 0 0 1]));
ht2a(2)=conn_menu('pushbutton',temp+[1*temp(3) 0 0 0],'','kurtosis','Spatial kurtosis (click to sort)',@(varargin)conn_dynexplore_update([0 0 0 0 2]));
ht2a(3)=conn_menu('pushbutton',temp+[2*temp(3) 0 0 0],'','skewness','Spatial skewness (click to sort)',@(varargin)conn_dynexplore_update([0 0 0 0 3]));
ht2a(4)=conn_menu('pushbutton',temp+[3*temp(3) 0 0 0],'','variability','Temporal component timeseries standard deviation averaged across all subjects (click to sort)',@(varargin)conn_dynexplore_update([0 0 0 0 4]));
ht2a(5)=conn_menu('pushbutton',temp+[4*temp(3) 0 0 0],'','frequency','Temporal component timeseries frequency (Hz) averaged across all subjects (click to sort)',@(varargin)conn_dynexplore_update([0 0 0 0 5]));
ht1=conn_menu('listbox2',boffset+[.80 .26 .10 .08],'Subjects',arrayfun(@(n)sprintf('Subject %d',n),1:CONN_x.Setup.nsubjects,'uni',0),'Select subject(s) for display',@(varargin)conn_dynexplore_update);
ht7=conn_menu('edit2',boffset+[.78 .38 .05 .04],'','2','display threshold',@(varargin)conn_dynexplore_update([0 0 1 0 0]));
ht7b=conn_menu('text2',boffset+[.73 .38 .05 .035],'','threshold');
[ht6a(1) ht6a(2)]=conn_menu('popup2',boffset+[.80 .18 .16 .04],'Seed ROI',ROInamesDisplay,'<HTML>Select seed ROI for dynamic connectivity timeseries display</HTML>',@(varargin)conn_dynexplore_selectroi);
[ht6b(1) ht6b(2)]=conn_menu('popup2',boffset+[.80 .11 .16 .04],'Target ROI',ROInamesDisplay,'<HTML>Select target ROI for dynamic connectivity timeseries display</HTML>',@(varargin)conn_dynexplore_selectroi);
ht3=conn_menu('slider',boffset+[.20 .10 .55 .04],'Time',1,'<HTML>Select time point</HTML>',@(varargin)conn_dynexplore_update);
try, addlistener(ht3, 'ContinuousValueChange',@(varargin)conn_dynexplore_update); end
set(ht3,'visible','off');
conn_menumanager('onregion',ht3,1,boffset+[.18 .08 .59 .08]);
ht4=conn_menu('popup2',boffset+[.60 .84 .18 .04],'',{'Spatial components','Estimated dynamic connectivity','Static vs. Dynamic connectivity plot'},'<HTML>Select connectivity display <br/> - <i>Spatial components</i> shows the spatial component scores for each selected circuit(s) (connections showing joint modulation of functional connectivity)<br/> - <i>Estimated dynamic connectivity</i> shows the estimated ROI-to-ROI connectivity at each timepoint/scan <br/> - <i>Static vs. Dynamic connectivity plot</i> shows a scatter plot comparing the static vs. dynamic connectivity (at each timepoint/scan) for each connection</HTML>',@(varargin)conn_dynexplore_update);
ht4b=conn_menu('popup2',boffset+[.60 .34 .20 .04],'',{'Display full connectivity','Display selected circuit(s) contribution only'},'<HTML>Select full connectivity to display the full dynamic ROI-to-ROI connectivity values or <br/> only the portion explained by the selected circuit(s)</HTML>',@(varargin)conn_dynexplore_update);
ht4c=conn_menu('pushbutton2',boffset+[.78 .84 .05 .04],'','Play','<HTML>Display dynamic connectivity over time</HTML>',@(varargin)conn_dynexplore_play);
%conn_menumanager('onregion',ht4c,1,boffset+[.65 .05 .3 .775]);
set(ht1,'max',2,'value',1);
set(ht2,'max',2,'value',1:numel(names));
set(ht4,'value',1);
set(ht4c,'interruptible','on','userdata',0);
set(ht6a(1),'value',idxlock(1));
set(ht6b(1),'value',idxlock(2));

%cmap=[[flipud(gray(128))*diag([0 0 1]);(gray(128))*diag([1 0 0])]; gray(256); CONN_gui.backgroundcolor];
cmap=[jet(256); gray(256); CONN_gui.backgroundcolor];
hfig=gcf;
handleplot=[];
handleimag=[];
handleload=[];
handlescor=[];
colororder=[];
handleloads1s2=[];
handleloaddoff=4;
c0=[]; c=[]; 
sortby=0;
showdisp=0;
conn_dynexplore_update([1 1 1 1 1]);
%conn_dynexplore_play;

    function conn_dynexplore_update(doreset)
        if nargin<1, doreset=[0 0 0 0 0]; end
        if numel(doreset)<5, doreset=[doreset zeros(1,5-numel(doreset))]; end
        try, if isfield(CONN_h,'menus')&&isfield(CONN_h.menus,'waiticonObj'), CONN_h.menus.waiticonObj.start; end; end
        nsub=get(ht1,'value');
        ifac=get(ht2,'value');
        nt=round(get(ht3,'value'));
        dtype=get(ht4,'value');
        ptype=get(ht4b,'value');
        thr=str2num(get(ht7,'string')); if isempty(thr), thr=2; set(ht7,'string',num2str(thr)); end
        nfacselected=nfacshown(ifac);
        itime=find(ismember(IDX_subject,nsub));
        Nt=numel(itime);
        Nf=size(B,1);
        ntold=nt; nt=max(1,min(Nt, nt));
        set(ht3,'min',1,'max',Nt,'sliderstep',1/max(1,Nt)*[1 2]);
        if nt~=ntold, set(ht3,'value',nt); end
        if isempty(itime), ntime=[]; 
        else ntime=itime(max(1,min(numel(itime),nt)));
        end

        %disp(idxlock);
        if isempty(ntime)
            c0=nan(size(B0,2),size(B0,3));
            c1=c0;
            c=c0;
        else
            c0=reshape(H0(ntime,:)*B0(:,:),size(B0,2),size(B0,3));
            c1=reshape(H(ntime,nfacselected)*B(nfacselected,:),size(B,2),size(B,3));
            if ptype==1, c=reshape(H(ntime,:)*B(:,:),size(B,2),size(B,3));
            else c=c1;
            end
        end
        imagepos=[.50 .42 .45 .41];
        %c0(1:size(c0,1)+1:end)=nan;
        
        if doreset(5)
            if doreset(5)==abs(sortby), sortby=-sortby;
            else sortby=doreset(5);
            end
            switch(abs(sortby))
                case 1, m=1:numel(names);
                case 2, m=KURx;
                case 3, m=SKWx;
                case 4, m=mean(VARt(1,:),1); 
                case 5, m=mean(ZCt(1,:),1);  
            end
            if sortby>0, [nill,idx]=sort(m);
            else         [nill,idx]=sort(m,'descend');
            end
            nfacshown=idx;
            invidx=idx; invidx(idx)=1:numel(idx);
            ifac=sort(invidx(nfacselected));
            nfacselected=nfacshown(ifac);
            doreset(3)=1;
        end
        str=cellfun(@(a,b,c,d,e)sprintf(regexprep('%-0s%-0s%-0s%-0s%-0s','0',ht3fieldsize),a,num2str(b),regexprep(num2str(c),'NaN','-'),num2str(d),num2str(e)),names,num2cell(KURx),num2cell(SKWx),num2cell(mean(VARt(1,:),1)),num2cell(mean(ZCt(1,:),1)),'uni',0);
        set(ht2,'string',str(nfacshown),'value',ifac,'fontname','monospaced','fontsize',8+CONN_gui.font_offset);
        if doreset(1)||isempty(handleplot)||any(~ishandle(handleplot))
            conn_menumanager('onregionremove',handleplot(ishandle(handleplot)));delete(handleplot(ishandle(handleplot)));
            h1=axes('units','norm','position',imagepos);
            colororder=get(h1,'colororder');
            %h2=plot(c0(idxorder,idxorder),c0(idxorder,idxorder)+c(idxorder,idxorder),'.',[-2 2 nan -2 2 nan 0 0],[-2 2 nan 0 0 nan -2 2],'k','markersize',6);
            temp=max(-2,min(2,c0(idxorder,idxorder)+c(idxorder,idxorder)));
            temp0=c0(idxorder,idxorder);
            h2=[]; hold on;
            for n2=1:Nf, h2=[h2 plot(temp0(Bhf_max{n2}),temp(Bhf_max{n2}),'.','markersize',6)]; end
            h2=[h2 plot([-2 2 nan -2 2 nan 0 0],[-2 2 nan 0 0 nan -2 2],'k','color',1-CONN_gui.backgroundcolor)]; hold off;
            hold on; h4=plot(0,0,'wo','visible','off','markersize',12); hold off;
            h3=uicontrol('units','norm','position',[.65 .1 .0001 .0001],'style','text','fontsize',7+CONN_gui.font_offset,'foregroundcolor','k','backgroundcolor','w','visible','off','horizontalalignment','left'); 
            axis equal tight;
            conn_menumanager('onregion',h3,1,imagepos,h1,@conn_dynexplore_mtnplot);
            set(h1,'xlim',[-1.5 1.5],'ylim',[-1.5 1.5],'color',CONN_gui.backgroundcolor,'xcolor',CONN_gui.backgroundcolor,'ycolor',CONN_gui.backgroundcolor,'xtick',[],'ytick',[]);
            set([h2(:)' h1 h4],'buttondownfcn',@conn_dynexplore_buttonpress);
            h1a=xlabel('Static connectivity'); set(h1a,'fontsize',CONN_gui.font_offset+8,'color',.1*[1 1 1]+.8*(.5>mean(CONN_gui.backgroundcolor)));
            h1b=ylabel('Dynamic connectivity'); set(h1b,'fontsize',CONN_gui.font_offset+8,'color',.1*[1 1 1]+.8*(.5>mean(CONN_gui.backgroundcolor)));
            handleplot=[h2(:)' h1a h1b h1 h4 h3];
        end
        if doreset(2)||isempty(handleimag)||any(~ishandle(handleimag))
            conn_menumanager('onregionremove',handleimag(ishandle(handleimag)));delete(handleimag(ishandle(handleimag)));
            h1=axes('units','norm','position',imagepos,'color',CONN_gui.backgroundcolor);
            if numel(idxorder)>100, 
                h2=image(conn_ind2rgb(128.5+127*(c0(idxorder,idxorder)+c(idxorder,idxorder)),cmap(1:256,:)));
            else
                [nill,h2]=conn_menu_plotmatrix(c0(idxorder,idxorder)+c(idxorder,idxorder),'scaletoarea',true);
            end
            hold on; h4=plot(0,0,'k+','visible','off','markersize',12); hold off;
            h3=uicontrol('units','norm','position',[.65 .1 .0001 .0001],'style','text','fontsize',7+CONN_gui.font_offset,'foregroundcolor','k','backgroundcolor','w','visible','off','horizontalalignment','left'); 
            %h3=text(0,0,1,''); set(h3,'visible','off','backgroundcolor','w','color','k','fontsize',CONN_gui.font_offset+7);
            axis equal tight;
            conn_menumanager('onregion',h3,1,imagepos,h1,@conn_dynexplore_mtnimag);
            set(h1,'color',CONN_gui.backgroundcolor,'xcolor',CONN_gui.backgroundcolor,'ycolor',CONN_gui.backgroundcolor,'xtick',[],'ytick',[],'ydir','reverse','xlim',[0,numel(idxorder)+1],'ylim',[0,numel(idxorder)+1]);
            set([h2 h1 h4],'buttondownfcn',@conn_dynexplore_buttonpress);
            %h=title({'Connectivity matrix','(total: static + dynamic)'}); set(h,'fontsize',CONN_gui.font_offset+8,'color',.1*[1 1 1]+.8*(.5>mean(CONN_gui.backgroundcolor)));
            handleimag=[h2 h1 h4 h3];
        end
        if doreset(3)||isempty(handleload)||any(~ishandle(handleload))
            conn_menumanager('onregionremove',handleload(ishandle(handleload)));delete(handleload(ishandle(handleload)));
            h1=axes('units','norm','position',imagepos,'color',CONN_gui.backgroundcolor);
            temp=permute(B(nfacselected,idxorder,idxorder),[2 3 1]);
            belowthr=abs(temp)<thr;
            temp=temp/max(abs(temp(:)));
            temp(belowthr)=1.5+.4*temp(belowthr);
            [temp2,handleloads1s2]=conn_menu_montage(h1,permute(temp,[1,2,4,3]));
            temp2=reshape(temp2,[numel(idxorder),handleloads1s2(2),numel(idxorder),handleloads1s2(1)]);
            temp=nan([numel(idxorder)+handleloaddoff handleloads1s2(2) numel(idxorder)+handleloaddoff handleloads1s2(1)]);
            temp(handleloaddoff/2+1:end-handleloaddoff/2,:,handleloaddoff/2+1:end-handleloaddoff/2,:)=temp2;
            temp=reshape(temp,[(numel(idxorder)+handleloaddoff)*handleloads1s2(2),(numel(idxorder)+handleloaddoff)*handleloads1s2(1)]);
            h2=image(conn_ind2rgb(round(128.5+127*temp),cmap));
            hold on; h4=plot(0,0,'k+','visible','off','markersize',12); hold off;
            h3=uicontrol('units','norm','position',[.65 .1 .0001 .0001],'style','text','fontsize',7+CONN_gui.font_offset,'foregroundcolor','k','backgroundcolor','w','visible','off','horizontalalignment','left'); 
            axis equal tight;
            conn_menumanager('onregion',h3,1,imagepos,h1,@conn_dynexplore_mtnload);
            set(h1,'color',CONN_gui.backgroundcolor,'xcolor',CONN_gui.backgroundcolor,'ycolor',CONN_gui.backgroundcolor,'xtick',[],'ytick',[],'ydir','reverse');
            set([h2 h1 h4],'buttondownfcn',@conn_dynexplore_buttonpress);
            %hold on; h2a=text(get(gca,'xlim')*[1.05;-.05],mean(get(gca,'ylim')),{'Spatial','components'},'color',.1*[1 1 1]+.8*(.5>mean(CONN_gui.backgroundcolor)),'horizontalalignment','right','fontsize',CONN_gui.font_offset+8); hold off;
            %h=ylabel('component scores'); set(h,'fontsize',CONN_gui.font_offset+8,'color',.1*[1 1 1]+.8*(.5>mean(CONN_gui.backgroundcolor)));
            handleload=[h2 h1 h4 h3];
        end
        if doreset(4)||isempty(handlescor)||any(~ishandle(handlescor))
            delete(handlescor(ishandle(handlescor)));
            h1=axes('units','norm','position',[.20 .12 .55 .24]);
            if isempty(itime), h2=plot(nan(2,Nf),nan(2,Nf),'-',nan,nan,'w-',[0 size(H,1)+1 nan 0 size(H,1)+1],[0 0 nan -2 -2],'k-',[0 0],[-4 2],'w:');
            else h2=plot(repmat((1:Nt)',1,Nf),H(itime,:),'-',1:Nt,zeros(1,Nt),'w-',[0 size(H,1)+1 nan 0 size(H,1)+1],[0 0 nan -2 -2],'k-',[0 0],[-4 2],'w:');
            end
            set(h2(end),'color',.1*[1 1 1]+.8*(.5>mean(CONN_gui.backgroundcolor)));
            set(h2(end-1),'color',.5*[1 1 1]);
            set(h2(Nf+1),'linewidth',2,'color',.1*[1 1 1]+.8*(.5>mean(CONN_gui.backgroundcolor)));
            axis tight;
            hold on; h3=text(0,1.6,'','horizontalalignment','left','color',.1*[1 1 1]+.8*(.5>mean(CONN_gui.backgroundcolor)),'fontsize',CONN_gui.font_offset+7); hold off;
            hold on; h4=text(0,1.6-2,'','interpreter','tex','horizontalalignment','right','color',.1*[1 1 1]+.8*(.5>mean(CONN_gui.backgroundcolor)),'fontsize',CONN_gui.font_offset+7); hold off;
            set(h1,'xlim',[0 Nt+1],'ylim',[-3.5 1.5],'color',CONN_gui.backgroundcolor,'xcolor',CONN_gui.backgroundcolor,'ycolor',CONN_gui.backgroundcolor,'xtick',[],'ytick',[]);
            h=xlabel('Time / scans'); set(h,'fontsize',CONN_gui.font_offset+8,'color',.1*[1 1 1]+.8*(.5>mean(CONN_gui.backgroundcolor)));
            hold on; h2a=text(get(gca,'xlim')*[1.025;-.025],0,{'Temporal','components'},'color',.1*[1 1 1]+.8*(.5>mean(CONN_gui.backgroundcolor)),'horizontalalignment','right','fontsize',CONN_gui.font_offset+8); 
            h2b=text(get(gca,'xlim')*[1.025;-.025],-2,{'ROI-to-ROI','connectivity','(dynamic)^(^*^)'},'interpreter','tex','color',.1*[1 1 1]+.8*(.5>mean(CONN_gui.backgroundcolor)),'horizontalalignment','right','fontsize',CONN_gui.font_offset+8); hold off;
            handlescor=[h2(:)' h3 h4 h1 h2b];
        end
        switch(dtype)
            case 1, set([handleplot handleimag],'visible','off'); set(handleload(1:end-1),'visible','on'); set([ ht7 ht7b],'visible','on'); set([ht3 ht4b ht4c],'visible','off');
            case 2, set([handleplot handleload],'visible','off'); set(handleimag(1:end-1),'visible','on'); set([ ht7 ht7b],'visible','off'); set([ht3 ht4b ht4c],'visible','on');
            case 3, set([handleimag handleload],'visible','off'); set(handleplot(1:end-1),'visible','on'); set([ ht7 ht7b],'visible','off'); set([ht3 ht4b ht4c],'visible','on');
        end
        temp=max(-2,min(2,c0(idxorder,idxorder)+c(idxorder,idxorder))); %temp(1:size(temp,1)+1:end)=nan;
        %temp=c(idxorder,idxorder); %temp(1:size(temp,1)+1:end)=nan;
        temp0=c0(idxorder,idxorder); 
        %temp0=max(-2,min(2, c0(idxorder,idxorder)+c1(idxorder,idxorder) ));
        
        %for n=1:size(temp,2), set(handleplot(n),'xdata',temp0(:,n),'ydata',temp(:,n)); end
        for n=1:Nf, set(handleplot(n),'xdata',temp0(Bhf_max{n}),'ydata',temp(Bhf_max{n})); end
        colororder=get(get(handleplot(1),'parent'),'colororder');
        %if idxlock(3), 
        %    for n=1:numel(idxorder), set(handleplot(n),'markersize',6,'zdata',zeros(1,numel(idxorder)),'color',.5*colororder(1+mod(n-1,size(colororder,1)),:)+.5*CONN_gui.backgroundcolor); end
        %    n=find(idxorder==idxlock(2),1); set(handleplot(n),'markersize',12,'zdata',ones(1,numel(idxorder)),'color','y');
        %else
        %    for n=1:numel(idxorder), set(handleplot(n),'markersize',6,'zdata',zeros(1,numel(idxorder)),'color',colororder(1+mod(n-1,size(colororder,1)),:)); end
        %end
        %if numel(nfacselected)<Nf, for n=1:Nf, set(handleplot(n),'zdata',ismember(n,nfacselected)+zeros(size(Bhf_max{n})), 'markersize',6,'color',colororder(1+mod(max([0,find(nfacselected==n,1)]),size(colororder,1)),:)); end
        if numel(nfacselected)<Nf, for n=1:Nf, ismembern=ismember(n,nfacselected); set(handleplot(n),'zdata',ismembern+zeros(size(Bhf_max{n})), 'markersize',6,'color',.25*colororder(1+mod(n-1,size(colororder,1)),:)+.75*(ismembern+(1-2*ismembern)*CONN_gui.backgroundcolor)); end
        else for n=1:Nf, set(handleplot(n),'zdata',Bhf_maxval{n},'markersize',6,'color',colororder(1+mod(n-1,size(colororder,1)),:)); end
        end
        set(handleplot(end-1),'xdata',c0(idxlock(1),idxlock(2)),'ydata',c0(idxlock(1),idxlock(2))+c(idxlock(1),idxlock(2)),'zdata',1);
        set(handleimag(end-1),'xdata',find(idxorder==idxlock(2),1),'ydata',find(idxorder==idxlock(1),1));
        set(handleload(end-1),'xdata',handleloaddoff/2+find(idxorder==idxlock(2),1)+(handleloaddoff+handleloads1s2(3))*mod(0:numel(nfacselected)-1,handleloads1s2(1)),'ydata',handleloaddoff/2+find(idxorder==idxlock(1),1)+(handleloaddoff+handleloads1s2(4))*floor((0:numel(nfacselected)-1)/handleloads1s2(1)));
        if numel(idxorder)>100, 
            temp=conn_ind2rgb(128.5+127*(c0(idxorder,idxorder)+c(idxorder,idxorder)),cmap(1:256,:));
            set(handleimag(1),'cdata',temp);
        else
            temp=conn_menu_plotmatrix(c0(idxorder,idxorder)+c(idxorder,idxorder),'scaletoarea',true);
            set(handleimag(1),'vertices',temp.vertices,'faces',temp.faces,'facevertexcdata',temp.facevertexcdata);
            set(handleimag(2),'xlim',[0,numel(idxorder)+1],'ylim',[0,numel(idxorder)+1]);
        end           
        for n=1:Nf, if any(nfacselected==n), set(handlescor(n),'xdata',1:Nt,'ydata',H(itime,n)/Hscale); else set(handlescor(n),'xdata',1:Nt,'ydata',nan(1,Nt)); end; end
        if ptype==1, tempr=H0(itime,:)*B0(:,idxlock(1),idxlock(2))+H(itime,:)*B(:,idxlock(1),idxlock(2));  % complete fit
        else tempr=H0(itime,:)*B0(:,idxlock(1),idxlock(2))+H(itime,nfacselected)*B(nfacselected,idxlock(1),idxlock(2)); % only portion explained by selected factor
        end
        set(handlescor(Nf+1),'xdata',1:Nt,'ydata',-2+tempr);
        set(handlescor(Nf+3),'xdata',nt+[0 0]); 
        set(handlescor(end-3),'position',[nt 1.6 0],'string',sprintf('Subject %d Session %d Scan %d',IDX_subject(ntime),IDX_session(ntime),IDX_scan(ntime))); 
        hext=get(handlescor(end-3),'extent'); hext=hext(end-1:end); set(handlescor(end-3),'position',[min(Nt,nt+hext(1))-hext(1) 1.6 0]); 
        set(handlescor(end-2),'position',[Nt -1.5-2 0],'string',sprintf('^(^*^) Connectivity between %s and %s',ROInamesDisplay2{idxlock(1:2)})); 
        set(handlescor(end-1),'xlim',[0 Nt+1]);
        if dtype==1, set(handlescor([Nf+1 end-2 end]),'visible','off'); set(handlescor(end-1),'ylim',[-1.5 1.5]);
        else set(handlescor([Nf+1 end-2 end]),'visible','on'); set(handlescor(end-1),'ylim',[-3.5 1.5]);
        end
        try, if isfield(CONN_h,'menus')&&isfield(CONN_h.menus,'waiticonObj'), CONN_h.menus.waiticonObj.stop; end; end
    end

    function conn_dynexplore_mtnimag(varargin)
        pos0=get(hfig,'currentpoint');
        if ~idxlock(3)
            pos=get(handleimag(end-2),'currentpoint');
            pos=pos(1,1:2);
        else pos=[find(idxorder==idxlock(2),1) find(idxorder==idxlock(1),1)];
        end
        if all(pos>0&pos<=numel(ROInames)),
            pos=max(1,min(numel(ROInames), round(pos)));
            hf=arrayfun(@(a,b)sprintf(' Circuit %d (%.3f)',a,b),Bhf_order(1:min(size(B,1),3),idxorder(pos(1)),idxorder(pos(2))),Bhf_values(1:min(size(B,1),3),idxorder(pos(1)),idxorder(pos(2))),'uni',0)';
            hf1=arrayfun(@(a,b)sprintf(' Circuit %d (%.3f)',a,b),Bhfc_order(1:min(size(B,1),3),idxorder(pos(1))),Bhfc_values(1:min(size(B,1),3),idxorder(pos(1))),'uni',0)';
            hf2=arrayfun(@(a,b)sprintf(' Circuit %d (%.3f)',a,b),Bhfc_order(1:min(size(B,1),3),idxorder(pos(2))),Bhfc_values(1:min(size(B,1),3),idxorder(pos(2))),'uni',0)';
            set(handleimag(end),'units','pixels','string',[{'Connectivity between:' ROInamesDisplay{idxorder(pos(2))} ['and ',ROInamesDisplay{idxorder(pos(1))}] ' ' 'Highest component scores for this connection:'} hf]); % {' ' ['Highest component score for ' ROInamesDisplay{idxorder(pos(2))}]} hf2 {' ' ['Highest component score for ' ROInamesDisplay{idxorder(pos(1))}]} hf1]);
            hext=get(handleimag(end),'extent'); hext=hext(end-1:end);
            set(handleimag(end),'position',[pos0(1:2)+[10 -hext(2)-10] hext]);
            set(handleimag(end-1),'xdata',pos(1),'ydata',pos(2),'visible','on');
            idxlock(1:2)=[idxorder(pos(2)) idxorder(pos(1))];
            set(ht6a(1),'value',idxlock(1));
            set(ht6b(1),'value',idxlock(2));
            conn_dynexplore_update;
        end
    end
    function conn_dynexplore_mtnload(varargin)
        pos0=get(hfig,'currentpoint');
        if ~idxlock(3)
            pos=get(handleload(end-2),'currentpoint');
            pos=pos(1,1:2);
        else pos=[handleloaddoff/2+find(idxorder==idxlock(2),1) handleloaddoff/2+find(idxorder==idxlock(1),1)];
        end
        if all(pos>0&pos<=(numel(ROInames)+handleloaddoff)*handleloads1s2(1:2)),
            pos1=max(1,min((numel(ROInames)+handleloaddoff)*handleloads1s2(1:2), round(pos)));
            idxf=max(1,min(numel(nfacselected), ceil(pos1(1)/(numel(ROInames)+handleloaddoff))+handleloads1s2(1)*(ceil(pos1(2)/(numel(ROInames)+handleloaddoff))-1) ));
            pos=max(1,min(numel(ROInames), 1+mod(pos1-1,numel(ROInames)+handleloaddoff)-handleloaddoff/2 ));
            hf=arrayfun(@(a,b)sprintf(' Circuit %d (%.3f)',a,b),Bhf_order(1:min(size(B,1),3),idxorder(pos(1)),idxorder(pos(2))),Bhf_values(1:min(size(B,1),3),idxorder(pos(1)),idxorder(pos(2))),'uni',0)';
            ht=sprintf('Circuit %d (%.3f)',nfacselected(idxf),B(nfacselected(idxf),idxorder(pos(1)),idxorder(pos(2))));
            hf1=arrayfun(@(a,b)sprintf(' Circuit %d (%.3f)',a,b),Bhfc_order(1:min(size(B,1),3),idxorder(pos(1))),Bhfc_values(1:min(size(B,1),3),idxorder(pos(1))),'uni',0)';
            hf2=arrayfun(@(a,b)sprintf(' Circuit %d (%.3f)',a,b),Bhfc_order(1:min(size(B,1),3),idxorder(pos(2))),Bhfc_values(1:min(size(B,1),3),idxorder(pos(2))),'uni',0)';
            set(handleload(end),'units','pixels','string',[{ht 'Connectivity between:' ROInamesDisplay{idxorder(pos(2))} ['and ',ROInamesDisplay{idxorder(pos(1))}] ' ' 'Highest component scores for this connection:'} hf]);% {' ' ['Highest component score for ' ROInamesDisplay{idxorder(pos(2))}]} hf2 {' ' ['Highest component score for ' ROInamesDisplay{idxorder(pos(1))}]} hf1]);
            %set(handleload(end),'units','pixels','string',[{ht 'Connectivity between:' ROInamesDisplay{idxorder(pos(2))} ['and ',ROInamesDisplay{idxorder(pos(1))}] ' ' 'Highest component scores for this connection:'} hf]);
            hext=get(handleload(end),'extent'); hext=hext(end-1:end);
            set(handleload(end),'position',[pos0(1:2)+[10 -hext(2)-10] hext]);
            %set(handleload(end-1),'xdata',pos1(1),'ydata',pos1(2),'visible','on');
            idxlock(1:2)=[idxorder(pos(2)) idxorder(pos(1))];
            set(ht6a(1),'value',idxlock(1));
            set(ht6b(1),'value',idxlock(2));
            conn_dynexplore_update;
        end
    end
    function conn_dynexplore_mtnplot(varargin)
        pos0=get(hfig,'currentpoint');
        if ~idxlock(3)
            pos=get(handleplot(end-2),'currentpoint');
            pos=pos(1,1:2);
        else pos=[c0(idxlock(1),idxlock(2)) c0(idxlock(1),idxlock(2))+c(idxlock(1),idxlock(2))];
        end
        [nill1,idx1]=min((c0-pos(1)).^2+(c+c0-pos(2)).^2,[],1);
        [nill2,idx2]=min(nill1,[],2);
        idx1=idx1(idx2);
        if nill2<.25^2,
            pos=max(1,min(numel(ROInames), [idx2 idx1]));
            hf=arrayfun(@(a,b)sprintf(' Circuit %d (%.3f)',a,b),Bhf_order(1:min(size(B,1),3),pos(1),pos(2)),Bhf_values(1:min(size(B,1),3),pos(1),pos(2)),'uni',0)';
            hf1=arrayfun(@(a,b)sprintf(' Circuit %d (%.3f)',a,b),Bhfc_order(1:min(size(B,1),3),pos(1)),Bhfc_values(1:min(size(B,1),3),pos(1)),'uni',0)';
            hf2=arrayfun(@(a,b)sprintf(' Circuit %d (%.3f)',a,b),Bhfc_order(1:min(size(B,1),3),pos(2)),Bhfc_values(1:min(size(B,1),3),pos(2)),'uni',0)';
            set(handleplot(end),'units','pixels','string',[{'Connectivity between:' ROInamesDisplay{pos(2)} ['and ',ROInamesDisplay{pos(1)}] ' ' 'Highest component scores for this connection:'} hf]);% {' ' ['Highest component score for ' ROInamesDisplay{pos(2)}]} hf2 {' ' ['Highest component score for ' ROInamesDisplay{pos(1)}]} hf1]);
            %set(handleplot(end),'units','pixels','string',[{'Connectivity between:' ROInamesDisplay{pos(2)} ['and ' ROInamesDisplay{pos(1)}] ' ','Highest component scores:'} hf]);
            hext=get(handleplot(end),'extent'); hext=hext(end-1:end);
            set(handleplot(end),'position',[pos0(1:2)+[10 -hext(2)-10] hext]);
            set(handleplot(end-1),'xdata',c0(idx1,idx2),'ydata',c0(idx1,idx2)+c(idx1,idx2),'visible','on');
            idxlock(1:2)=[idx1 idx2];
            set(ht6a(1),'value',idxlock(1));
            set(ht6b(1),'value',idxlock(2));
            conn_dynexplore_update;
        end
    end
    function conn_dynexplore_buttonpress(varargin)
        idxlock(3)=~idxlock(3);
    end
    function conn_dynexplore_selectroi(varargin)
        idxlock(1)=get(ht6a(1),'value');
        idxlock(2)=get(ht6b(1),'value');
        idxlock(3)=1;
        conn_dynexplore_update;
    end
    function conn_dynexplore_play(varargin)
        state=get(ht4c,'userdata');
        if state
            set(ht4c,'string','Play','userdata',0);
        else
            nsub=get(ht1,'value');
            nt=round(get(ht3,'value'));
            itime=find(ismember(IDX_subject,nsub));
            Nt=numel(itime);
            ntold=nt; nt=max(1,min(Nt, nt));
            set(ht3,'min',1,'max',Nt,'sliderstep',1/Nt*[1 2]);
            if nt~=ntold, set(ht3,'value',nt); end
            set(ht4c,'string','Stop','userdata',1);
            for n1=nt:Nt
                set(ht3,'value',n1);
                conn_dynexplore_update;
                drawnow;
                pause(.02);
                if ~ishandle(ht4c), return; end
                state=get(ht4c,'userdata');
                if ~state, break; end
            end
            if n1==Nt, set(ht3,'value',1); end
            set(ht4c,'string','Play','userdata',0);
            conn_dynexplore_update;
        end
    end
end

function rout=conn_ind2rgb(a,cm,dim)
a = max(1,min(size(cm,1),round(a)));
rout=reshape(cm(a,:),size(a,1),[],size(cm,2));
if nargin>2, 
    if dim==1,      rout=permute(rout, [3,1,2]);
    elseif dim==2,  rout=permute(rout, [1,3,2]);
    end
end
end


