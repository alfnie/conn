function fh=conn_montage_display(x,xlabels,style, xcovariate, xcovariate_name, xborders, xborders_name, xclusters, xclusters_name, xrois, xrois_name)

% CONN_MONTAGE_DISPLAY displays 3d data as montage
% fh = conn_montage_display(x [,labels]);
%   displays a montage of each matrix x(:,:,:,i), labeled by labels{i}
%   size(x,3)==3 : for rgb values
%   size(x,3)==1 : for colormap scaling
%

global CONN_x CONN_gui;
if isempty(CONN_gui)||~isfield(CONN_gui,'font_offset'), conn_font_init; end

if isstruct(x)
    state=x;
else
    if nargin<2, xlabels={}; end
    if ischar(x)||iscell(x)
        fname=cellstr(x);
        txyz=[];
        dispdata={};
        displabel={};
        donewarning='conn_montage_display';
        for n=1:numel(fname)
            temp=cellstr(fname{n});
            if numel(temp)==1,temp=cellstr(conn_expandframe(temp{1}));end
            nvols=unique([1 numel(temp)]);
            for nvol=nvols(:)'
                files=conn_fileutils('spm_vol',temp{nvol});
                if conn_surf_dimscheck(files), files=conn_surf_parent(files,donewarning); donewarning=''; end % surface
                if numel(txyz)<=1
                    dim=files(1).dim(1:2);
                    [tx,ty]=ndgrid(1:dim(1),1:dim(2));
                    if numel(txyz)==1, zslice=txyz;
                    else zslice=round(files(1).dim(3)/2);
                    end
                    txyz=files(1).mat*[tx(:) ty(:) zslice+zeros(numel(tx),1) ones(numel(tx),1)]';
                end
                dispdata{end+1}=fliplr(flipud(reshape(conn_fileutils('spm_get_data',files(1),pinv(files(1).mat)*txyz),dim(1:2))'));
                if isempty(xlabels)
                    displabel{end+1}=temp{nvol};
                else
                    if numel(nvols)>1, displabel{end+1}=sprintf('%s,%d',xlabels{n},nvol);
                    else displabel{end+1}=sprintf('%s',xlabels{n});
                    end
                end
            end
        end
        x=cat(4,dispdata{:});
        xlabels=displabel;
    end
    if nargin>=3&&~isempty(style)&&~isempty(regexp(style,'_scaled'))&&size(x,3)==1,
        tempx=abs(x);
        tempx(isnan(x))=0; 
        x=x./repmat(mean(mean(tempx,1),2),[size(x,1),size(x,2),1,1]); 
        style=regexprep(style,'_scaled','');
    end
    state.x=x;
    state.xlabels=xlabels;
    cmap2=repmat(1-linspace(1,0,128)'.^2,[1,3]).*hot(128)+repmat(linspace(1,0,128)'.^2,[1,3])*.1;
    cmap2=cmap2([1,21,35:end],:);
    state.colormap_default=cmap2;
    cmap2=[flipud(cmap2(:,[2,3,1]));cmap2];
    state.colormap=cmap2;
    if nargin>=3&&isequal(style,'matrix'), state.colormap=jet(2*96).^repmat(1+0*abs(linspace(1,-1,2*96))',1,3); end
    %state.colormap_default=[flipud(fliplr(hot))*diag([.5 .5 1]);gray];
    %state.cmap=state.colormap_default;
    state.bookmark_filename='';
    state.bookmark_descr='';
end
if nargin<3||isempty(style), style='montage'; end
if nargin<4||isempty(xcovariate), xcovariate=[]; end
if nargin<5||isempty(xcovariate_name), xcovariate_name={}; end
if nargin<6, xborders=[]; end
if nargin<7, xborders_name={}; end
if nargin<8, xclusters=[]; end
if nargin<9, xclusters_name={}; end
if nargin<10, xrois=[]; end
if nargin<11, xrois_name={}; end
if isfield(CONN_gui,'slice_display_skipbookmarkicons'), SKIPBI=CONN_gui.slice_display_skipbookmarkicons;
else SKIPBI=false;
end
fh=@(varargin)conn_montage_display_refresh([],[],varargin{:});
state.handles.fh=fh;
if ~isfield(CONN_x,'folders')||~isfield(CONN_x.folders,'bookmarks')||isempty(CONN_x.folders.bookmarks), state.dobookmarks=false;
else state.dobookmarks=true;
end
state.loop=0;
if ~isfield(state,'style'), state.style=style; end
if ~isfield(state,'style0'), state.style0=style; end
if ~isfield(state,'slide'), state.slide=1; end
if ~isfield(state,'xcov'), 
    state.xcov=xcovariate; 
    state.xcov_name=xcovariate_name;
    %if ~isempty(xcovariate), state.xcov=state.xcov./repmat(max(eps,max(abs(state.xcov),[],1)),size(state.xcov,1),1); end
end
if ~isfield(state,'x_border'), 
    state.x_border=xborders; 
    state.x_border_name=xborders_name; 
end
if ~isfield(state,'xclusters'), 
    state.xclusters=xclusters; 
    state.xclusters_name=xclusters_name; 
end
if ~isfield(state,'xrois'), 
    state.xrois=xrois; 
    state.xrois_name=xrois_name; 
end
if ~isfield(state,'x_orig'), state.x_orig=state.x; end
if ~isfield(state,'fontsize'), 
    if isempty(xborders_name), state.fontsize=6+CONN_gui.font_offset; 
    else state.fontsize=[6 4]+CONN_gui.font_offset; 
    end
end
if ~isfield(state,'skiplabels'), state.skiplabels=1; end
if ~isfield(state,'xnonzero'),%||~isfield(state,'xnonzeroorder'),
    state.xnonzero=repmat(any(any(state.x~=0,4),3),[1,1,size(state.x,3),size(state.x,4)]);
%     idx=find(any(any(state.x~=0,4),3));
%     if numel(idx)<=2e3
%         clX=reshape(detrend(reshape(state.x(state.xnonzero),[],size(state.x,4))','constant')',numel(idx),[]);
%         clk1=sum(clX.^2,2);
%         clk2=clX*clX';
%         clZ=conn_bsxfun(@plus,clk1,clk1')-2*clk2;
%         clI=tril(true(size(clX,1)),-1);
%         clL=conn_statslinkage(clZ(clI>0)', 'co');
%         [nill,nill,idx]=conn_statsdendrogram(clL,0);
%         state.xnonzeroorder=idx;
%     else
%         state.xnonzeroorder=1:numel(idx);
%     end
end

%x=x./(max(max(max(abs(x),[],1),[],2),[],3));
if strcmp(state.style0,'matrix'), pos=[.4 .25 .6 .7];
elseif any(strcmp(state.style,{'movie','moviereplay','timeseries'})), pos=[.35 .35 .5 .5];
else pos=[.4 .15 .6 .8];
end
if any(strcmp(state.style,{'matrix','timeseries'}))||strcmp(state.style0,'matrix'), fcolor=[1 1 1];
else fcolor=[0 0 0];
end
state.handles.hfig=conn_figure('units','norm','position',pos,'color',fcolor,'menubar','none','name',sprintf('%s display',state.style),'numbertitle','off','colormap',state.colormap);
hc=state.handles.hfig;
hc1=uimenu(hc,'Label','Effects');
if size(state.x,3)==1,
    hc2=uimenu(hc1,'Label','colormap');
    for n1={'normal','red','jet','hot','gray','bone','cool','hsv','spring','summer','autumn','winter','bluewhitered','random','brighter','darker','manual','color'}
        uimenu(hc2,'Label',n1{1},'callback',{@conn_montage_display_refresh,'colormap',n1{1}});
    end
end
hc2=uimenu(hc1,'Label','colorscale');
uimenu(hc2,'Label','colorbar limits','callback',{@conn_montage_display_refresh,'colorscale','rescale'});
uimenu(hc2,'Label','colorscale direct','callback',{@conn_montage_display_refresh,'colorscale','direct'});
uimenu(hc2,'Label','colorscale equalized','callback',{@conn_montage_display_refresh,'colorscale','equalize'});
hc2=uimenu(hc1,'Label','background');
uimenu(hc2,'Label','white background','callback',{@conn_montage_display_refresh,'background',[1 1 1]});
uimenu(hc2,'Label','light background','callback',{@conn_montage_display_refresh,'background',[.95 .95 .9]});
uimenu(hc2,'Label','dark background','callback',{@conn_montage_display_refresh,'background',[.11 .11 .11]});
uimenu(hc2,'Label','black background','callback',{@conn_montage_display_refresh,'background',[0 0 0]});
if ~isempty(state.xrois_name)||~isempty(state.x_border_name)||~isempty(state.xcov_name)
    hc2=uimenu(hc1,'Label','fontsize');
    uimenu(hc2,'Label','increase labels fontsize','callback',{@conn_montage_display_refresh,'fontsize','+'});
    uimenu(hc2,'Label','decrease labels fontsize','callback',{@conn_montage_display_refresh,'fontsize','-'});
    uimenu(hc2,'Label','set labels fontsize','callback',{@conn_montage_display_refresh,'fontsize','?'});
    uimenu(hc2,'Label','reduce number of labels shown','callback',{@conn_montage_display_refresh,'skiplabels','-'});
    uimenu(hc2,'Label','increase number of labels shown','callback',{@conn_montage_display_refresh,'skiplabels','+'});
end
hc2=uimenu(hc1,'Label','style');
if ~strcmp(style,'timeseries')
    uimenu(hc2,'Label','montage','callback',{@conn_montage_display_refresh,'style','montage'});
    uimenu(hc2,'Label','movie','callback',{@conn_montage_display_refresh,'style','movie'});
end
if size(state.x,3)==1,
    uimenu(hc2,'Label','timeseries','callback',{@conn_montage_display_refresh,'style','timeseries'});
    uimenu(hc2,'Label','matrix','callback',{@conn_montage_display_refresh,'style','matrix'});
end
hc1=uimenu(hc,'Label','Print');
uimenu(hc1,'Label','current view','callback',{@conn_montage_display_refresh,'print',1});
uimenu(hc1,'Label','as video','callback',{@conn_montage_display_refresh,'printvideo'});
if state.dobookmarks
    hc1=uimenu(hc,'Label','Bookmark');
    hc2=uimenu(hc1,'Label','Save','callback',{@conn_montage_display_refresh,'bookmark'});
    if ~isempty(state.bookmark_filename),
        hc2=uimenu(hc1,'Label','Save as copy','callback',{@conn_montage_display_refresh,'bookmarkcopy'});
    end
end
drawnow;
state.handles.hax=axes('units','norm','position',[0 0 1 1],'parent',state.handles.hfig);
[state.y,state.nX]=conn_menu_montage(state.handles.hax,state.x(:,:,:,end));
state.datalim=max([0;abs(state.x(:))]);
state.datalim_orig=state.datalim;
if strcmp(state.style,'matrix')
    tempy=state.y;
    tempy(isnan(state.y))=0;
    tempy=state.colormap(max(1,min(size(state.colormap,1), round((size(state.colormap,1)+1)/2+(size(state.colormap,1)-1)/2*tempy/state.datalim))),:);
    tempy(state.y==0,:)=.95;
    tempy(isnan(state.y),:)=1;
    %tempy(1:size(state.y,1)+1:size(state.y,1)*size(state.y,2),:)=1;
    state.y=reshape(tempy,size(state.y,1),size(state.y,2),3);
else    
    state.y(isnan(state.y))=0;
end
if size(state.y,3)==1, state.handles.him=image((size(state.colormap,1)+1)/2+(size(state.colormap,1)-1)/2*state.y/state.datalim,'parent',state.handles.hax); 
else state.handles.him=imagesc(state.y,'parent',state.handles.hax);
end
axis(state.handles.hax,'equal','off');
[state.handles.hyticks,state.handles.hylticks]=deal([]);
if strcmp(state.style,'matrix'), 
    hold(state.handles.hax,'on'); 
    if size(state.y,2)<=256, 
        plot3(repmat(.5:size(state.y,2)+.5,2,1),repmat([.5;size(state.y,1)+.5],1,size(state.y,2)+1),ones(2,size(state.y,2)+1),'w-','linewidth',1+(size(state.y,1)<32),'parent',state.handles.hax);
        plot3(repmat([.5;size(state.y,2)+.5],1,size(state.y,1)+1),repmat(.5:size(state.y,1)+.5,2,1),ones(2,size(state.y,1)+1),'w-','linewidth',1+(size(state.y,1)<32),'parent',state.handles.hax);
    end
    hold(state.handles.hax,'off'); 
end
if ~isempty(state.xrois_name) % rows/columns names
    hold(state.handles.hax,'on');
    state.handles.hyticks=text(0-.01*size(state.y,2)+zeros(1,numel(state.xrois_name)),1:numel(state.xrois_name),state.xrois_name(:)','horizontalalignment','right','fontsize',max(1,state.fontsize(end)-3),'color',.25*[1 1 1],'interpreter','none','parent',state.handles.hax);
    hold(state.handles.hax,'off'); 
end
if ~isempty(state.x_border) % clusters of rows/columns
    hold(state.handles.hax,'on');
    plot3(repmat(state.x_border(:)',2,1),repmat([.5;size(state.y,1)+.5],1,numel(state.x_border)),1+ones(2,numel(state.x_border)),'k-','color',.75*[1 1 1],'linewidth',1+(size(state.y,1)<32),'parent',state.handles.hax);
    plot3(repmat([.5;size(state.y,2)+.5],1,numel(state.x_border)),repmat(state.x_border(:)',2,1),1+ones(2,numel(state.x_border)),'k-','color',.75*[1 1 1],'linewidth',1+(size(state.y,1)<32),'parent',state.handles.hax);
    tb=[.5 state.x_border(:)' size(state.y,1)+.5];
    for n1=1:numel(tb)-1
        plot3([.4 0 0 .4],tb(n1+[0 0 1 1])+[0 .2 -.2 0],1+zeros(1,4),'k-','color',0*[1 1 1],'linewidth',1+(size(state.y,1)<32),'parent',state.handles.hax);
    end
    if ~isempty(state.x_border_name)
        state.handles.hylticks=text(.5-.30*size(state.y,2)+zeros(1,numel(state.x_border_name)),convn([.5,state.x_border(:)',size(state.y,1)+.5],[.5 .5],'valid'),state.x_border_name(:)','horizontalalignment','right','fontsize',max(1,state.fontsize(1)-1),'parent',state.handles.hax);
    end
    hold(state.handles.hax,'off'); 
end
if ~isempty(state.xclusters) % clusters of matrix elements
    hold(state.handles.hax,'on');
    for n1=reshape(unique(state.xclusters(state.xclusters>0)),1,[])
        ta=state.xclusters==n1;
        [ti,tj]=find(ta); % crop to save mem
        ti=[min(ti) max(ti)];
        tj=[min(tj) max(tj)];
        ta=ta(ti(1):ti(2),tj(1):tj(2));
        if numel(ta)<5e3, xfact=16; else xfact=4; end
        tb=double(ta(ceil((1:xfact*size(ta,1))/xfact),ceil((1:xfact*size(ta,2))/xfact)));
        tb=double(convn(tb,ones(3))>=4); % include diagonal neighb
        tc=contourc(tb,[.5 .5]);
        while ~isempty(tc)
            plot3((tc(1,1+(1:tc(2,1)))-1.5)/xfact+.5+tj(1)-1,(tc(2,1+(1:tc(2,1)))-1.5)/xfact+.5+ti(1)-1,3+zeros(1,tc(2,1)),'k-','color',.15*[1 1 1],'linewidth',1+(size(state.y,1)<128),'parent',state.handles.hax);
            tc=tc(:,1+tc(2,1)+1:end);
        end
    end
    hold(state.handles.hax,'off');
%set(state.handles.hax,'xcolor',[1 1 1],'ycolor',[1 1 1],'xtick',.5:size(state.y,2)+.5,'ytick',.5:size(state.y,1)+.5,'xticklabel',[],'yticklabel',[],'visible','on','box','on'); grid(state.handles.hax,'on');
end
clim=get(state.handles.hax,'clim');
if clim(1)/max(abs(clim))<-1e-2, clim=max(abs(clim))*[-1 1]; set(state.handles.hax,'clim',clim); set(state.handles.hfig,'colormap',state.colormap); end
set(state.handles.hfig,'resizefcn',{@conn_montage_display_refresh,'refresh'});
if strcmp(state.style,'matrix')
    state.handles.haxcolorbar=axes('units','norm','position',[.94 .4 .02 .2],'parent',state.handles.hfig);
    state.handles.hcolorbar=image(permute(flipud(state.colormap),[1,3,2]),'parent',state.handles.haxcolorbar);
    set(state.handles.haxcolorbar,'xlim',[.5 1.5],'ylim',[.5 size(state.colormap,1)+.5],'xtick',[],'ytick',[],'box','on');
    hold(state.handles.haxcolorbar,'on');
    state.handles.htxtcolorbar1=text(1,.5+size(state.colormap,1)*1.05,num2str(-state.datalim),'fontsize',max(1,state.fontsize(end)-3),'color',.25*[1 1 1],'horizontalalignment','center');
    state.handles.htxtcolorbar2=text(1,.5-size(state.colormap,1)*0.05,num2str(state.datalim),'fontsize',max(1,state.fontsize(end)-3),'color',.25*[1 1 1],'horizontalalignment','center');
    hold(state.handles.haxcolorbar,'off');
else
    [state.handles.haxcolorbar,state.handles.hcolorbar,state.handles.htxtcolorbar1,state.handles.htxtcolorbar2]=deal([]);
end 
if ~isempty(state.xlabels)||~isempty(state.xrois_name)||strcmp(state.style,'matrix'),
    state.handles.hlabel=uicontrol('style','text','horizontalalignment','left','visible','off','parent',state.handles.hfig);
    set(state.handles.hfig,'units','pixels','windowbuttonmotionfcn',@conn_menu_montage_figuremousemove);
    %drawnow;
end
state.handles.slider=uicontrol('style','slider','units','norm','position',[.4 .01 .2 .02],'foregroundcolor','w','backgroundcolor',[.95 .95 .9],'parent',state.handles.hfig,'callback',{@conn_montage_display_refresh,'slider'});
try, addlistener(state.handles.slider, 'ContinuousValueChange',@(varargin)conn_montage_display_refresh([],[],'slider')); end
set(state.handles.slider,'min',1,'max',size(state.x,4),'sliderstep',min(1,[1,10]/max(1,size(state.x,4)-1)),'value',1);
state.handles.startstop=uicontrol('style','togglebutton','units','norm','position',[.05 0 .1 .05],'string','Play','parent',state.handles.hfig,'callback',{@conn_montage_display_refresh,'startstop'});
state.handles.singleloop=uicontrol('style','checkbox','units','norm','position',[.85 0 .15 .05],'string','loop','parent',state.handles.hfig,'backgroundcolor',[.95 .95 .9],'foregroundcolor',.85*[1 1 1],'value',0);
state.handles.movietitle=uicontrol('style','text','units','norm','position',[.05 .95 .90 .05],'string','','parent',state.handles.hfig,'horizontalalignment','left','foregroundcolor',.25*[1 1 1],'backgroundcolor',[.95 .95 .9]);
if ~isempty(state.xcov)
    state.handles.haxcov=axes('units','norm','position',[.2 .075 .6 .20],'parent',state.handles.hfig);
    maxxcov=max(eps,max(abs(state.xcov),[],1));
    maxxcov=maxxcov.*(max(maxxcov.^.25)./(maxxcov.^.25));
    xcov=repmat(size(state.xcov,2)-1:-1:0, size(state.xcov,1),1) + .65*state.xcov./repmat(max(eps,maxxcov),size(state.xcov,1),1);
    state.handles.himcov=plot(xcov,'parent',state.handles.haxcov,'color',.5*[1 1 1]);%state.handles.himcov=image((size(state.colormap,1)+1)/2+(size(state.colormap,1)-1)/2*state.xcov','parent',state.handles.haxcov);
    hold(state.handles.haxcov,'on'); state.handles.refcov=plot([0 0],[min(xcov(:)) max(xcov(:))],'b','parent',state.handles.haxcov); hold(state.handles.haxcov,'off');%hold(state.handles.haxcov,'on'); state.handles.refcov=plot([1 1],[.5,size(state.xcov,2)+.5],'b','parent',state.handles.haxcov); hold(state.handles.haxcov,'off');
    hold(state.handles.haxcov,'on'); state.handles.scalecov=text(repmat(-5,[1,size(xcov,2)*2]),[min(xcov,[],1) max(xcov,[],1)],arrayfun(@num2str,[min(state.xcov,[],1) max(state.xcov,[],1)],'uni',0),'color',.75*[1 1 1],'horizontalalignment','right','fontsize',max(1,state.fontsize(end)-3),'parent',state.handles.haxcov); hold(state.handles.haxcov,'off');
    axis(state.handles.haxcov,'tight','off');
    if ~isempty(state.xcov_name),
        if numel(state.xcov_name)==1, 
            state.handles.titlecov=text(1.02*size(state.xcov,1),size(state.xcov,2)/2,state.xcov_name,'parent',state.handles.haxcov); 
            set(state.handles.titlecov,'rotation',90,'color',.5*[1 1 1],'horizontalalignment','center','fontsize',max(1,state.fontsize(1)-2),'interpreter','none');
        else
            state.handles.titlecov=text(1.02*size(state.xcov,1)+zeros(1,numel(state.xcov_name)),.25+size(state.xcov,2)-(1:numel(state.xcov_name)),state.xcov_name,'parent',state.handles.haxcov);
            set(state.handles.titlecov,'color',.5*[1 1 1],'horizontalalignment','left','fontsize',max(1,state.fontsize(1)-2),'interpreter','none');
        end
    end
else [state.handles.haxcov,state.handles.himcov,state.handles.scalecov,state.handles.refcov,state.handles.titlecov]=deal([]);
end

set(state.handles.hfig,'userdata',state);
fh('refresh');
if strcmp(state.style,'moviereplay')
    %drawnow;
    set(state.handles.startstop,'value',1);
    fh('startstop');
end

    
    function out=conn_montage_display_refresh(hObject,eventdata,option,varargin)
        out=[];
        doreset=false;
        if nargin<3||isempty(option), option='refresh'; end
        switch(option)
            case 'refresh', doreset=true;
            case 'slider',
            case 'close', state.loop=0; close(state.handles.hfig); return;
            case 'figurehandle', out=state.handles.hfig; return;
            case 'printvideo',
                if isempty(which('VideoWriter')), uiwait(errordlg('Sorry. VideoWriter functionality only supported on newer Matlab versions')); return; end
                videoformats={'*.avi','Motion JPEG AVI (*.avi)';'*.mj2','Motion JPEG 2000 (*.mj2)';'*.mp4;*.m4v','MPEG-4 (*.mp4;*.m4v)';'*.avi','Uncompressed AVI (*.avi)'; '*.avi','Indexed AVI (*.avi)'; '*.avi','Grayscale AVI (*.avi)'};
                [filename, pathname,filterindex]=uiputfile(videoformats,'Save video as','conn_video01.avi');
                if isequal(filename,0), return; end
                state.style='movie';
                defs_videowriteframerate=20; % fps
                objvideo = VideoWriter(fullfile(pathname,filename),regexprep(videoformats{filterindex,2},'\s*\(.*$',''));
                set(objvideo,'FrameRate',defs_videowriteframerate);
                open(objvideo);
                ss=1;
                set(state.handles.startstop,'string','Stop','value',1);
                for n=1:size(state.x,4),
                    state.slide=n;
                    try, 
                        set(state.handles.slider,'value',state.slide);
                        conn_montage_display_refresh([],[],'refresh');
                        set([state.handles.slider state.handles.startstop state.handles.singleloop],'visible','off');
                        ss=get(state.handles.startstop,'value');
                        currFrame=getframe(state.handles.hfig);
                        writeVideo(objvideo,currFrame);
                        drawnow;
                    catch, ss=0; break; 
                    end
                    if ~ss, break; end
                end
                close(objvideo);
                try, set([state.handles.slider state.handles.startstop state.handles.singleloop],'visible','on'); end
                if ~ss, return; end
                objvideo=[];
                objvideoname=get(objvideo,'Filename');
                try, set(state.handles.startstop,'string','Play','value',0); end
                conn_msgbox(sprintf('File %s created',fullfile(pathname,filename)),'');
                %try, if ispc, winopen(objvideoname); else system(sprintf('open %s',objvideoname)); end; end
                return;
            case 'print',
                conn_print(state.handles.hfig,varargin{:});
                return;
            case 'skiplabels'
                opt=varargin{1};
                if isequal(opt,'-'), opt=state.skiplabels+1;
                elseif isequal(opt,'+'), opt=state.skiplabels-1;
                end
                opt=max(1,min(numel(state.handles.hyticks)-1, opt));
                state.skiplabels=opt;
                if isfield(state.handles,'hyticks')&&all(ishandle(state.handles.hyticks)), set(state.handles.hyticks,'visible','off'); set(state.handles.hyticks(1:state.skiplabels:end),'visible','on'); end
            case 'fontsize'
                opt=varargin{1};
                if isequal(opt,'+'), opt=state.fontsize+1;
                elseif isequal(opt,'-'), opt=state.fontsize-1;
                elseif isequal(opt,'?'), opt=conn_menu_inputdlg('Enter fontsize','conn_montage_display',1,{num2str(state.fontsize)}); if ~isempty(opt), opt=str2num(opt{1}); end; if isempty(opt), return; end
                end
                state.fontsize=opt;
                if isfield(state.handles,'scalecov')&&all(ishandle(state.handles.scalecov)), set(state.handles.scalecov,'fontsize',max(1,state.fontsize(end)-3)); end
                if isfield(state.handles,'hyticks')&&all(ishandle(state.handles.hyticks)), set(state.handles.hyticks,'fontsize',max(1,state.fontsize(end)-3)); end
                if isfield(state.handles,'hylticks')&&all(ishandle(state.handles.hylticks)), set(state.handles.hylticks,'fontsize',max(1,state.fontsize(1)-1)); end
                if isfield(state.handles,'titlecov')&&all(ishandle(state.handles.titlecov)), set(state.handles.titlecov,'fontsize',max(1,state.fontsize(1)-2)); end
            case {'colorscale','colorbar'}
                opt=varargin{1};
                switch(opt)
                    case 'rescale'
                        if ~isequalwithequalnans(state.x,state.x_orig), return; end
                        if numel(varargin)>=2, val=varargin{2}; 
                        else
                            val=state.datalim;
                            val=conn_menu_inputdlg({'Enter new colorbar limit:'},'Rescale colorbar',1,{num2str(val)});
                            if ~isempty(val), val=str2num(val{1}); end
                        end
                        if isempty(val), return; end
                        if isequal(val,'symmetric'), val=max(abs(state.datalim)); end
                        state.datalim=max(abs(val));
                        set(state.handles.htxtcolorbar1,'string',num2str(-state.datalim));
                        set(state.handles.htxtcolorbar2,'string',num2str(state.datalim));
                    case 'equalize'
                        if ~isfield(state,'x_equalized')||~isfield(state,'colormap_equalized')
                            temp=state.x_orig;
                            temp(isnan(temp))=0;
                            [ut,nill,idx]=unique(abs(temp));
                            nidx=cumsum(accumarray(idx(:),1));
                            nidx=nidx-nidx(1);
                            temp=reshape(sign(temp(:)).*nidx(idx(:)),size(temp));
                            mat=max(abs(temp(:)));
                            state.x_equalized=temp/mat;
                            tidx=reshape(interp1(ut(:),nidx(:)/mat,linspace(0,max(abs(state.x_orig(:))),92)),[],1);
                            state.colormap_equalized=state.colormap(max(1,min(size(state.colormap,1), round((size(state.colormap,1)+1)/2+(size(state.colormap,1)-1)/2*cat(1,-flipud(tidx),tidx)))),:);
                        end
                        state.x=state.x_equalized;
                        state.colormap_plot=state.colormap_equalized;
                        state.datalim=max(abs(state.x(:)));
                        set(state.handles.htxtcolorbar1,'string',num2str(-state.datalim_orig));
                        set(state.handles.htxtcolorbar2,'string',num2str(state.datalim_orig));
                    case 'direct'
                        state.x=state.x_orig;
                        state.colormap_plot=state.colormap;
                        state.datalim=max(abs(state.x(:)));
                        set(state.handles.htxtcolorbar1,'string',num2str(-state.datalim_orig));
                        set(state.handles.htxtcolorbar2,'string',num2str(state.datalim_orig));
                end
            case {'start','stop','startstop'}
                if strcmp(option,'start'), state.loop=1; set(state.handles.startstop,'value',state.loop); state.style='moviereplay'; 
                elseif strcmp(option,'stop'), state.loop=0; set(state.handles.startstop,'value',state.loop);
                else state.loop=get(state.handles.startstop,'value');
                end
                if state.loop, % stop
                    set(state.handles.startstop,'string','Stop','value',1);
                    while 1,
                        state.slide=round(max(1,min(size(state.x,4), get(state.handles.slider,'value'))));
                        state.slide=1+mod(state.slide,size(state.x,4));
                        set(state.handles.slider,'value',state.slide);
                        conn_montage_display_refresh([],[],'refresh');
                        drawnow;
                        try, state.loop=state.loop & get(state.handles.startstop,'value'); infloop=get(state.handles.singleloop,'value');
                        catch, return;
                        end
                        if ~state.loop || (~infloop && state.slide==1), break; end % force single loop
                    end
                    state.loop=0;
                    try, set(state.handles.startstop,'string','Play','value',0); end
                else
                    set(state.handles.startstop,'string','Play','value',0);
                end
                return;
            case 'colormap'
                cmap=varargin{1};
                if ischar(cmap)
                    switch(cmap)
                        case 'normal', if strcmp(state.style0,'matrix'), cmap=jet(2*96).^repmat(1+0*abs(linspace(1,-1,2*96))',1,3); else cmap=state.colormap_default; end
                        case 'red', cmap=[linspace(0,1,96)',zeros(96,2)];
                        case 'hot', cmap=hot(96);
                        case 'jet', cmap=fixedge(jet(256));
                        case 'gray', cmap=gray(256);
                        case 'bone', cmap=bone(256);
                        case 'cool',cmap=fixedge(cool(256));
                        case 'hsv',cmap=fixedge(hsv(256));
                        case 'bluewhitered', cmap=[zeros(1,48) linspace(0,1,48) ones(1,48) linspace(1,.5,48); linspace(0,1,96) linspace(1,0,48) zeros(1,48); linspace(.5,1,48) ones(1,48) linspace(1,0,48) zeros(1,48)]'; cmap=repmat(abs(linspace(-1,1,192)'),1,3).*cmap+(1-repmat(abs(linspace(-1,1,192)'),1,3))*1;
                        case 'parulawhite', cmap=parula(192); cmap=repmat(abs(linspace(-1,1,192)'),1,3).*cmap+(1-repmat(abs(linspace(-1,1,192)'),1,3))*1;
                        case 'spring',cmap=spring(96);
                        case 'summer',cmap=summer(96);
                        case 'autumn',cmap=autumn(96);
                        case 'winter',cmap=winter(96);
                        case 'random',cmap=rand(96,3);
                        case 'brighter',cmap=min(1,1/sqrt(.95)*get(state.handles.hfig,'colormap').^(1/2)); cmap=cmap(round(size(cmap,1)/2)+1:end,:);
                        case 'darker',cmap=.95*get(state.handles.hfig,'colormap').^2; cmap=cmap(round(size(cmap,1)/2)+1:end,:);
                        case 'manual',answer=conn_menu_inputdlg({'colormap (96x3)'},'',1,{mat2str(state.colormap(round(size(state.colormap,1)/2)+1:end,:))});if ~isempty(answer), answer=str2num(answer{1}); end;if ~any(size(answer,1)==[96,2*96]), return; end;cmap=max(0,min(1,answer));
                        case 'color',cmap=uisetcolor([],'Select color'); if isempty(cmap)||isequal(cmap,0), return; end;
                        otherwise, disp('unknown value');
                    end
                end
                if ~isempty(cmap)
                    if size(cmap,2)<3, cmap=cmap(:,min(size(cmap,2),1:3)); end
                    if size(cmap,1)==1, cmap=linspace(0,1,96)'*cmap; end
                    if size(cmap,1)~=2*96, cmap=[flipud(cmap(:,[2,3,1]));cmap]; end
                    state.colormap=cmap;
                    state.colormap_plot=state.colormap;
                    state.x=state.x_orig;
                    state.datalim=max(abs(state.x(:)));
                    set(state.handles.hfig,'colormap',cmap);
                end
            case 'background'
                value=varargin{1};
                set([state.handles.slider state.handles.singleloop state.handles.movietitle],'backgroundcolor',value); 
                set(state.handles.hfig,'color',value);
                return;
            case 'style',
                state.loop=0; 
                state.style=varargin{1};
                doreset=true;
            case 'getstate',
                out=state;
                out=rmfield(out,'handles');
                return;
            case {'bookmark','bookmarkcopy'},
                tfilename=[];
                if numel(varargin)>0&&~isempty(varargin{1}), tfilename=varargin{1};
                elseif ~isempty(state.bookmark_filename)&&strcmp(option,'bookmark'), tfilename=state.bookmark_filename;
                end
                if numel(varargin)>1&&~isempty(varargin{2}), descr=cellstr(varargin{2});
                else descr=state.bookmark_descr;
                end
                fcn=regexprep(mfilename,'^conn_','');
                conn_args={fcn,conn_montage_display_refresh([],[],'getstate')};
                [fullfilename,tfilename,descr]=conn_bookmark('save',...
                    tfilename,...
                    descr,...
                    conn_args);
                if isempty(fullfilename), return; end
                if ~SKIPBI, 
                    tht=conn_msgbox('Printing bookmark icon. Please wait...','',-1);
                    conn_print(state.handles.hfig,conn_prepend('',fullfilename,'.jpg'),'-nogui','-r50','-nopersistent'); 
                    if ishandle(tht), delete(tht); end
                end
                state.bookmark_filename=tfilename;
                state.bookmark_descr=descr;
                conn_args={fcn,conn_montage_display_refresh([],[],'getstate')}; % re-save to include bookmark info
                conn_savematfile(conn_prepend('',fullfilename,'.mat'),'conn_args');
                if 0, conn_msgbox(sprintf('Bookmark %s saved',fullfilename),'',2);
                else out=fullfilename;
                end
                return;
        end
        switch(state.style)
            case 'montage'
                set(state.handles.hax,'position',[0 0 1 1]);
                if doreset, set(state.handles.hfig,'color',[0 0 0]); end
                [state.y,state.nX]=conn_menu_montage(state.handles.hax,state.x);
                set([state.handles.slider state.handles.startstop state.handles.singleloop state.handles.movietitle state.handles.haxcov state.handles.himcov(:)' state.handles.scalecov(:)' state.handles.refcov(:)'],'visible','off');
                axis(state.handles.hax,'equal'); 
                datalim=state.datalim;
            case 'matrix'
                set(state.handles.hax,'position',[.3 .1 .6 .8]);
                if doreset, set(state.handles.hfig,'color',[.95 .95 .9]); end
                state.slide=round(max(1,min(size(state.x,4), get(state.handles.slider,'value'))));
                [state.y,state.nX]=conn_menu_montage(state.handles.hax,state.x(:,:,state.slide));
                if size(state.x,4)>1, set([state.handles.startstop state.handles.singleloop state.handles.haxcov state.handles.himcov(:)' state.handles.scalecov(:)' state.handles.refcov(:)'],'visible','off');
                else set([state.handles.slider state.handles.startstop state.handles.singleloop state.handles.movietitle state.handles.haxcov state.handles.himcov(:)' state.handles.scalecov(:)' state.handles.refcov(:)'],'visible','off');
                end
                axis(state.handles.hax,'equal');
                datalim=state.datalim;
            case {'movie','moviereplay'}
                if strcmp(state.style0,'matrix'), set(state.handles.hax,'position',[.3 .1 .6 .8]); 
                elseif ~isempty(state.xcov), set(state.handles.hax,'position',[.05 .3 .90 .65]); 
                else set(state.handles.hax,'position',[.05 .05 .90 .90]); 
                end
                if doreset, set(state.handles.hfig,'color',[.95 .95 .9]); end
                if ~isempty(state.xcov), set([state.handles.himcov(:)' state.handles.scalecov(:)' state.handles.refcov(:)'],'visible','on'); end                
                state.slide=round(max(1,min(size(state.x,4), get(state.handles.slider,'value'))));
                [state.y,state.nX]=conn_menu_montage(state.handles.hax,state.x(:,:,:,state.slide));
                set([state.handles.slider state.handles.startstop state.handles.singleloop state.handles.movietitle],'visible','on');
                %axis(state.handles.hax,'equal');
                datalim=state.datalim;
            case {'timeseries'}
                if ~isempty(state.xcov), set([state.handles.himcov(:)' state.handles.scalecov(:)' state.handles.refcov(:)'],'visible','on'); set(state.handles.hax,'position',[.2 .3 .6 .55]);
                else set(state.handles.hax,'position',[.2 .05 .6 .80]);
                end
                if doreset, set(state.handles.hfig,'color',[.95 .95 .9]); end
                state.slide=1;%round(max(1,min(size(state.x,4), get(state.handles.slider,'value'))));
                temp=detrend(reshape(state.x(state.xnonzero),[],size(state.x,4))','constant');
                %temp=temp(:,state.xnonzeroorder);
                [state.y,state.nX]=conn_menu_montage(state.handles.hax,temp');
                set([state.handles.slider state.handles.startstop state.handles.singleloop],'visible','off'); %state.handles.refcov(:)'
                axis(state.handles.hax,'normal');
                state.y=.5+.5*state.y/max(abs(state.y(:)));
                datalim=1;
        end
        if strcmp(state.style0,'matrix')
            tempy=state.y;
            tempy(isnan(state.y))=0;
            tempy=state.colormap(max(1,min(size(state.colormap,1), round((size(state.colormap,1)+1)/2+(size(state.colormap,1)-1)/2*tempy/state.datalim))),:);
            tempy(state.y==0,:)=.95; %tempy(state.y==0,3)=.90;
            tempy(isnan(state.y),:)=1; 
            %tempy(1:size(state.y,1)+1:size(state.y,1)*size(state.y,2),:)=1;
            state.y=reshape(tempy,size(state.y,1),size(state.y,2),3);
        else
            state.y(isnan(state.y))=0;
        end
        if size(state.y,3)==1, set(state.handles.him,'cdata',(size(state.colormap,1)+1)/2+(size(state.colormap,1)-1)/2*state.y/datalim);
        else set(state.handles.him,'cdata',state.y);
        end
        set(state.handles.hax,'xlim',[-.05 size(state.y,2)+.55],'ylim',[.45 size(state.y,1)+.55],'clim',clim); 
        if any(strcmp(state.style0,{'matrix'})) % colorbar
            if size(state.x,4)>1&&state.slide<=numel(state.xlabels), 
                set(state.handles.movietitle,'string',state.xlabels{state.slide});
            end
            if isfield(state,'colormap_plot'), cmap=state.colormap_plot; 
            else cmap=state.colormap;
            end
            set(state.handles.hcolorbar,'cdata',permute(flipud(cmap),[1,3,2]));
            set(state.handles.haxcolorbar,'ylim',[.5,size(cmap,1)+.5]);
            try, ypos=get(state.handles.htxtcolorbar1,'position'); ypos(2)=.5+size(cmap,1)*1.05; set(state.handles.htxtcolorbar1,'position',ypos); end
            try, ypos=get(state.handles.htxtcolorbar2,'position'); ypos(2)=.5-size(cmap,1)*0.05; set(state.handles.htxtcolorbar2,'position',ypos); end
        end
        if any(strcmp(state.style,{'timeseries'})),
            if state.slide==numel(state.xlabels), 
                set(state.handles.movietitle,'string',state.xlabels{state.slide},'visible','on');
            else set(state.handles.movietitle,'visible','off');
            end
        end
        if any(strcmp(state.style,{'movie','moviereplay'})),
            if state.slide<=numel(state.xlabels), 
                set(state.handles.movietitle,'string',state.xlabels{state.slide});
            end
            if ~isempty(state.xcov)
                set(state.handles.refcov,'xdata',state.slide+[0 0]);
            end
        end
    end

    function conn_menu_montage_figuremousemove(varargin)
        if ~strcmp(state.style,'montage')&&~strcmp(state.style0,'matrix'), 
            set(state.handles.hlabel,'visible','off');
            return; 
        end
        p1=get(0,'pointerlocation');
        p2=get(state.handles.hfig,'position');
        p3=get(0,'screensize');
        p4=p2(1:2)+p3(1:2)-1; % note: fix issue when connecting to external monitor/projector
        pos0=(p1-p4);
        set(state.handles.hfig,'currentpoint',pos0);
        pos=(get(state.handles.hax,'currentpoint')); pos=pos(1,1:3);
        set(state.handles.hax,'units','pixels');posax=get(state.handles.hax,'position');set(state.handles.hax,'units','norm');
        if strcmp(state.style,'montage'), txyz=conn_menu_montage('coords2xyz',state.nX,pos(1:2)'); txyz=round(txyz(3));
        elseif strcmp(state.style0,'matrix'), txyz=round(pos(1:2));
        else txyz=state.slide;
        end
        tlabel={};
        if numel(txyz)>1&&size(state.x,1)==numel(state.xrois_name)&&size(state.x,2)==numel(state.xrois_name)&&all(txyz>=1&txyz<=numel(state.xrois_name))&&((~isnan(state.x(txyz(1),txyz(2)))&&state.x(txyz(1),txyz(2))~=0)||(~isnan(state.x(txyz(2),txyz(1)))&&state.x(txyz(2),txyz(1))~=0)), 
            tlabel=[state.xrois_name{txyz(2)}, ' - ', state.xrois_name{txyz(1)}];
            if ~isempty(state.xclusters)&&~isempty(state.xclusters_name)&&state.xclusters(txyz(2),txyz(1))>0, tlabel={tlabel, ['(in ',state.xclusters_name{state.xclusters(txyz(2),txyz(1))},')']}; end
        elseif numel(txyz)==2&&all(txyz>=1&txyz<=[size(state.x,2),size(state.x,1)])&&(~isnan(state.x(txyz(2),txyz(1)))&&state.x(txyz(2),txyz(1))~=0), 
            tlabel=mat2str(state.x(txyz(2),txyz(1)),4);
            if ~isempty(state.xclusters)&&~isempty(state.xclusters_name)&&state.xclusters(txyz(2),txyz(1))>0, tlabel={tlabel, ['(in ',state.xclusters_name{state.xclusters(txyz(2),txyz(1))},')']}; end
        elseif numel(txyz)==1&&txyz>=1&&txyz<=numel(state.xlabels)&&pos(1)>=1&&pos(1)<=state.nX(3)*state.nX(1)&&pos(2)>=1&&pos(2)<=state.nX(4)*state.nX(2)
            tlabel=state.xlabels{txyz};
        end
        if ~isempty(tlabel)
            set(state.handles.hlabel,'units','pixels','position',[pos0+[10 -10] 20 20],'visible','on','string',tlabel);
            hext=get(state.handles.hlabel,'extent');
            nlines=ceil(hext(3)/(p2(3)/2));
            ntlabel=numel(tlabel);
            newpos=[pos0+[-0*min(p2(3)/2,hext(3))/2 +20] min(p2(3)/2,hext(3)) nlines*hext(4)];
            newpos(1)=max(posax(1),newpos(1)-max(0,newpos(1)+newpos(3)-posax(1)-posax(3)));
            newpos(2)=max(posax(2),newpos(2)-max(0,newpos(2)+newpos(4)-posax(2)-posax(4)));
            %newpos(1)=max(0,newpos(1)-max(0,newpos(1)+newpos(3)-posax(3)));
            %newpos(2)=max(0,newpos(2)-max(0,newpos(2)+newpos(4)-posax(4)));
            set(state.handles.hlabel,'position',newpos,'string',reshape([tlabel,repmat(' ',1,nlines*ceil(ntlabel/nlines)-ntlabel)]',[],nlines)');
        else
            set(state.handles.hlabel,'visible','off');
        end
    end

end

function c=fixedge(c)
k=linspace(1,0,size(c,1))'.^4;
c=repmat(1-k,1,size(c,2)).*c+repmat(k.*mean(c,2),1,size(c,2));
end


