function fh = conn_polar_display(data, varargin)
% CONN_POLAR_DISPLAY displays polar data
%
% conn_polar_display(filename[, templatefile, ispercentage, manualscale, title])
% displays the overlap between ROIs or masks in 3D NIFTI volume "filename"
% and the 7-networks whole-brain parcellation (Yeo&Buckner 7-networks atlas)
%
%    filename     : 3D volume NIFTI file (ROI or mask file)
%    templatefile : 3D volume NIFTI template with ROIs/networks
%                   (defaults to conn/utils/surf/YeoBuckner2011.nii)
%    ispercentage : 0/1 0=raw numbers; 1=percentage overlap [0]
%    manualscale  : maximum value displayed
%
% conn_polar_display(data [, names, scaletoone, manualscale, title])
% displays the numbers in the data matrix data (one row per sample,
% one column per observation/variable) in polar display.
%    data         : Nx1 or NxM data matrix
%    names        : 1xN cell array with labels for each polar direction
%    ispercentage : 0/1 0=raw numbers; 1=percentage overlap [0]
%    manualscale  : maximum value displayed
% 
% fh=CONN_POLAR_DISPLAY(...) returns function handle implementing all GUI functionality
% note: if a conn_polar_display window is already open, its function handle can also be read using the syntax fh = gcf().UserData.handles.fh;
%
% VIEW OPTIONS
%  fh('style', option)                  : selects display style
%                                           'stacked'       : areas of individual ROIs are stacked on top of each other (like a stacked bar chart but in polar coordinates)
%                                           'layered'       : areas of individual ROIs are shown separately and superimposed
%  fh('select.roi', idx)                : displays only a subset of all ROIs (or columns of data) (idx: index to existing ROIs in input filename or columns of input data matrix)
%
% EFFECTS OPTIONS
%  fh('background',color)               : sets background color (color: [1x3] RGB values)
%  fh('colormap',type)                  : changes colormap (type: 'normal','red','jet','hot','gray','bone','cool','hsv','spring',
%                                           'summer','autumn','winter','bluewhitered','random','brighter','darker','manual','color')
%  fh('scaling', option)                : selects radial direction scaling. option may take the values:
%                                           'linear'        : intersection between bars and each polar direction scale linearly with K (the number of voxels in the intersection of each ROI and each network)
%                                           'equalized'     : intersection between bars and each polar direction scale with the square root of K (so that bar areas proportional to K)
%  fh('lines', option)                  : selects interpolation method:
%                                           'straight'      : interpolate values between two networks with straight lines
%                                           'curved'        : interpolate values between two networks with curved lines (pchip angular interpolation)
%  fh('edges', state)                   : displays bar edges (state: 'on' 'off')
%  fh('colormap',type)                  : changes colormap (type: 'normal','red','jet','hot','gray','bone','cool','hsv','spring',
%                                           'summer','autumn','winter','bluewhitered','random','brighter','darker','manual','color')
%  fh('fontsize', val)                  : changes text fontsize
%  fh('legend', state)                  : displays legend (state: 'on' 'off')
%
% PRINT OPTIONS
%  fh('togglegui',val);                  : shows/hides GUI (1/0)
%  fh('print',filename,'-nogui')         : prints display to high-resolution .jpg file
%
%
% e.g.
%   conn_polar_display('results.ROI.nii');
% e.g.
%   conn_polar_display(rand(7,1));
%

STYLE=2; % 1:linear; 2:curved
INTERP='pchip';
CUMSUM=false; 
DOLEGEND=true;
DOEDGES=true;
DOSQRT=true;
fh=@(varargin)conn_polar_display_refresh([],[],varargin{:});
if ~nargin, help(mfilename); return; end

if numel(varargin)>=1&&~isempty(varargin{1}), filetemplate=varargin{1}; else filetemplate=''; end
if numel(varargin)>=2&&~isempty(varargin{2}), ispercentage=varargin{2}; else ispercentage=[]; end
if numel(varargin)>=3&&~isempty(varargin{3}), manualscale=varargin{3}; else manualscale=[]; end
if numel(varargin)>=4&&~isempty(varargin{4}), datatitle=varargin{4}; else datatitle=''; end
datanames={};
names={};
templateorder=[];
if iscell(data)||ischar(data)
    if isempty(ispercentage), ispercentage=false; end
    filedata=data;
    if ischar(filedata), filedata={filedata}; end
    if isequal(filetemplate,'?'), 
        filetemplate='';
        answ=conn_questdlg('Reference atlas:','conn_polar_display','Yeo&Buckner 7 networks (2011)','Other','Yeo&Buckner 7 networks (2011)');
        if strcmp(answ,'Other'), 
            [tfilename,tpathname]=conn_fileutils('uigetfile','*.nii; *.img','Select reference atlas file');
            if ischar(tfilename), filetemplate=fullfile(tpathname,tfilename); end
        end
    end
    if isempty(filetemplate), 
        filetemplate=fullfile(fileparts(which(mfilename)),'utils','surf','YeoBuckner2011.nii'); 
        templateorder=[6,7,5,1,2,4,3]; % t(n) in what position should network#n should go
    end
    if conn_surf_dimscheck(filetemplate)
        assert(all(cellfun(@conn_surf_dimscheck,filedata)),'mismatch between template file (a surface NIFTI file) and data files (with one or serveral volume NIFTI files)');
    else
        for n=1:numel(filedata)
            if conn_surf_dimscheck(filedata{n}), filedata{n}=conn_surf_surf2vol(filedata{n},[],[],[], @mode); end
        end
    end
    labels=conn_vol_read(filetemplate,filedata{1});
    nlabels=max(labels(:));
    data=[];
    for nfiledata=1:numel(filedata)
        dataval=conn_vol_read(filedata{nfiledata});
        mask = dataval>0&labels>0;
        if size(dataval,4)>1
            if conn_existfile(conn_prepend('',filedata{nfiledata},'.txt'))
                datanames=regexp(conn_fileutils('fileread', conn_prepend('',filedata{nfiledata},'.txt')),'\n+','split');
                datanames=datanames(cellfun('length',datanames)>0);
            else
                datanames=arrayfun(@(n)sprintf('cluster%d',n),1:max(dataval(mask)),'uni',0);
            end
            if ispercentage, for n4=1:size(dataval,4), data=cat(2, data, accumarray(labels(mask(:,:,:,n4)),1,[nlabels,1])/nnz(dataval(:,:,:,:,n4)>0)/nnz(dataval(:,:,:,n4)>0)); end
            else for n4=1:size(dataval,4), data=cat(2, data, accumarray(labels(mask(:,:,:,n4)),1,[nlabels,1])); end
            end
        elseif all(rem(dataval(mask),1)==0), % ROI data file
            if conn_existfile(conn_prepend('',filedata{nfiledata},'.txt'))
                datanames=regexp(conn_fileutils('fileread',conn_prepend('',filedata{nfiledata},'.txt')),'\n+','split');
                datanames=datanames(cellfun('length',datanames)>0);
            else
                datanames=arrayfun(@(n)sprintf('cluster%d',n),1:max(dataval(mask)),'uni',0);
            end
            if ispercentage, data = cat(2, data, accumarray([labels(mask),dataval(mask)],1,[nlabels,numel(datanames)])./repmat(accumarray(dataval(dataval>0),1,[numel(datanames),1])',[nlabels,1]));
            else data = cat(2, data, accumarray([labels(mask),dataval(mask)],1,[nlabels,numel(datanames)]));
            end
        else
            if ispercentage, data = cat(2, data, accumarray(labels(mask),1,[nlabels,1])/nnz(dataval>0));
            else data = cat(2, data, accumarray(labels(mask),1,[nlabels,1]));
            end
        end
    end
    [names,ROIidx]=conn_roilabels(filetemplate);
    names(ROIidx)=names;
    if ~isempty(templateorder),
        data=data(templateorder,:);
        names=names(templateorder);
    end
else 
    if isempty(ispercentage), ispercentage=false; end
    names=filetemplate;
end
if size(data,1)==1, data=data.'; end

Nx=size(data,1);
Ny=size(data,2);
isint=~nnz(rem(data(:),1)~=0);
strdescrip={}; 
strdescrip_title='';
if ~isempty(names)
    strdescrip_title=sprintf('%-16s   ',datatitle);
    for nx=1:Nx, strdescrip_title=[strdescrip_title,sprintf(' %10s',upper(names{nx}(1:min(10,numel(names{nx})))))]; end
end
if isempty(datanames), datanames=arrayfun(@(n)sprintf('measure%d',n),1:Ny,'uni',0); end
if isint, numdescrip=' %10d'; else numdescrip=' %10f'; end
for ny=1:Ny
    strdescrip{end+1}=sprintf('%-16s   ',datanames{ny});
    for nx=1:Nx, strdescrip{end}=[strdescrip{end},sprintf(numdescrip,data(nx,ny))]; end %│
end
strdescrip{end+1}=sprintf('%-16s   ','(total)');
for nx=1:Nx, strdescrip{end}=[strdescrip{end},sprintf(numdescrip,sum(data(nx,:)))]; end
data=[data;data(1,:)];
data0=data;

state.colormap=jet(2*96).^repmat(1+0*abs(linspace(1,-1,2*96))',1,3); %state.colormap=state.colormap(end+1-(1:2*96),:);
state.background=[1 1 1];
state.fontsize=10;
state.handles.hfig=figure('units','norm','position',[.4 .25 .4 .7],'color',state.background,'menubar','none','name','polar display','numbertitle','off','colormap',state.colormap);
hc=state.handles.hfig;
hc1=uimenu(hc,'Label','View');
hc2=uimenu(hc1,'Label','style');
uimenu(hc2,'Label','stacked view','callback',{@conn_polar_display_refresh,'style','stacked'},'tag','style');
uimenu(hc2,'Label','layered view','checked',true,'callback',{@conn_polar_display_refresh,'style','layered'},'tag','style');
hc2=uimenu(hc1,'Label','Select ROIs','callback',{@conn_polar_display_refresh,'select.ROIs'});
hc1=uimenu(hc,'Label','Effects');
hc2=uimenu(hc1,'Label','scaling');
uimenu(hc2,'Label','linear','callback',{@conn_polar_display_refresh,'scaling','linear'},'tag','scaling');
uimenu(hc2,'Label','equalized','checked',true,'callback',{@conn_polar_display_refresh,'scaling','equalized'},'tag','scaling');
hc2=uimenu(hc1,'Label','lines');
uimenu(hc2,'Label','straight','callback',{@conn_polar_display_refresh,'lines','straight'},'tag','lines');
uimenu(hc2,'Label','curved','checked',true,'callback',{@conn_polar_display_refresh,'lines','curved'},'tag','lines');
hc2=uimenu(hc1,'Label','edges');
uimenu(hc2,'Label','on','checked',true,'callback',{@conn_polar_display_refresh,'edges','on'},'tag','edges');
uimenu(hc2,'Label','off','callback',{@conn_polar_display_refresh,'edges','off'},'tag','edges');
if Ny>1,
    hc2=uimenu(hc1,'Label','colormap');
    for n1={'normal','red','jet','hot','gray','bone','cool','hsv','spring','summer','autumn','winter','bluewhitered','random','brighter','darker','manual','color'}
        uimenu(hc2,'Label',n1{1},'callback',{@conn_polar_display_refresh,'colormap',n1{1}},'tag','colormap');
    end
end
hc2=uimenu(hc1,'Label','background');
uimenu(hc2,'Label','white background','callback',{@conn_polar_display_refresh,'background',[1 1 1]},'tag','background');
uimenu(hc2,'Label','light background','callback',{@conn_polar_display_refresh,'background',[.95 .95 .9]},'tag','background');
uimenu(hc2,'Label','dark background','callback',{@conn_polar_display_refresh,'background',[.11 .11 .11]},'tag','background');
uimenu(hc2,'Label','black background','callback',{@conn_polar_display_refresh,'background',[0 0 0]},'tag','background');
hc2=uimenu(hc1,'Label','fontsize');
uimenu(hc2,'Label','increase fontsize','callback',{@conn_polar_display_refresh,'fontsize','+'});
uimenu(hc2,'Label','decrease fontsize','callback',{@conn_polar_display_refresh,'fontsize','-'});
uimenu(hc2,'Label','set fontsize','callback',{@conn_polar_display_refresh,'fontsize','?'});
hc2=uimenu(hc1,'Label','legend');
uimenu(hc2,'Label','legend on','checked',true,'callback',{@conn_polar_display_refresh,'legend','on'},'tag','legend');
uimenu(hc2,'Label','legend off','callback',{@conn_polar_display_refresh,'legend','off'},'tag','legend');
hc1=uimenu(hc,'Label','Print');
uimenu(hc1,'Label','Print plot only','callback',{@conn_polar_display_refresh,'print'});
uimenu(hc1,'Label','Print full window','callback',{@conn_polar_display_refresh,'printall'});

if ~isempty(strdescrip_title), state.handles.table_title=uicontrol('units','norm','position',[.05 .20 .9 .03],'style','text','visible','off','parent',state.handles.hfig); end
state.handles.table=uicontrol('units','norm','position',[.05 .05 .9 .15],'style','listbox','max',2,'visible','off','parent',state.handles.hfig);
if ~isempty(which('conn_menu_search')),
    set(state.handles.table,'keypressfcn',@conn_menu_search);
    hc1=uicontextmenu('parent',state.handles.hfig);
    uimenu(hc1,'Label','Export table','callback',{@conn_polar_display_refresh,'export'});
    uimenu(hc1,'Label','Export overlap data','callback',{@conn_polar_display_refresh,'export_data'});
    set(state.handles.table,'uicontextmenu',hc1);
end
set(state.handles.table,'callback',{@conn_polar_display_refresh,'list'});
state.handles.fh=fh;
set(state.handles.hfig,'userdata',state);

selectedROIs=1:Ny;
fh('refresh');


    function out=conn_polar_display_refresh(hObject,eventdata,option,varargin)
        try
            if numel(hObject)==1&&ishandle(hObject)&&~isempty(get(hObject,'tag'))
                str=get(hObject,'tag');
                set(findobj(state.handles.hfig,'type','uimenu','tag',str),'checked','off');
                set(hObject,'checked','on');
            elseif isempty(hObject)&&isempty(eventdata)&&ischar(option)
                h=findobj(state.handles.hfig,'type','uimenu');
                s1=get(h,'callback');
                idx1=find(cellfun('length',s1)>1);
                idx2=find(cellfun(@iscell,s1(idx1)));
                idx3=find(cellfun(@(x)strcmp(x{2},option)&isequal(x(3:end),varargin),s1(idx1(idx2))));
                if numel(idx3)==1
                    h2=findobj(state.handles.hfig,'type','uimenu','tag',get(h(idx1(idx2(idx3))),'tag'));
                    set(h2,'checked','off');
                    set(h(idx1(idx2(idx3))),'checked','on');
                end
            end
        end
        out=[];
        doreset=false;
        if nargin<3||isempty(option), option='refresh'; end
        switch(lower(option))
            case 'refresh', 
                if isfield(state.handles,'axes')&&ishandle(state.handles.axes), delete(state.handles.axes); end
                if isfield(state.handles,'title')&&ishandle(state.handles.title), delete(state.handles.title); end
                data=data0(:,selectedROIs);
                Ny=numel(selectedROIs);
                if CUMSUM&&size(data,2)>1, data=cumsum(data,2); end
                ang=0*pi/2+(0:Nx)'*(2*pi)/Nx;
                state.facecolor=[.25 .25 .22; state.background+sign([.5 .5 .5]-state.background).*[.025 .025 .025]]; % in; out; 
                state.linecolor=[.55 .55 .55; .75 .75 .75]; % in; out
                if mean(state.background)<.5, 
                   state.linecolor=1-state.linecolor; 
                end
                state.handles.patch=[];
                state.handles.line={};
                if numel(datanames)>1
                    color=[state.colormap(round(linspace(1,size(state.colormap,1),numel(datanames))),:)];
                    color=color(selectedROIs,:);
                else color=state.facecolor(1,:);
                end
                state.handles.axes=axes('units','norm','position',[.2 .3 .6 .6],'color',state.background);
                hold(state.handles.axes,'on');
                if ~isempty(manualscale), R=manualscale; % manually scale
                elseif CUMSUM&&size(data0,2)>1, R=1.25*max(max(cumsum(abs(data0),2)));
                else R=1.25*max(abs(data0(:))); % scale set to 25% above maximum value
                end
                if ispercentage, R=min(1, R); end
                if DOSQRT, rtrans=@sqrt; 
                else rtrans=@(x)x; 
                end

                state.handles.backcircle=patch(rtrans(R)*cos(linspace(0,2*pi,1e3)),rtrans(R)*sin(linspace(0,2*pi,1e3)),'k-','facecolor',state.facecolor(2,:),'edgecolor','none','parent',state.handles.axes);
                %plot(R*cos(linspace(0,2*pi,1e3)),R*sin(linspace(0,2*pi,1e3)),'k:','color',state.linecolor(2,:),'linewidth',2,'parent',state.handles.axes);
                plot([zeros(1,Nx);rtrans(R)*cos(ang(1:end-1))'],[zeros(1,Nx);rtrans(R)*sin(ang(1:end-1))'],'k-','color',state.linecolor(2,:),'linewidth',1,'parent',state.handles.axes);
                for nx=1:Nx,
                    patch(.975*rtrans(R)*cos(ang(nx))+.025*rtrans(R)*cos(ang(nx)+2*pi/3*[-1 0 1 -1]),.975*rtrans(R)*sin(ang(nx))+.025*rtrans(R)*sin(ang(nx)+2*pi/3*[-1 0 1 -1]),'k','facecolor',state.linecolor(2,:),'edgecolor',state.linecolor(2,:),'parent',state.handles.axes);
                end
                if ispercentage
                    for nx=.1:.1:R
                        plot(rtrans(nx)*cos(linspace(0,2*pi,1e3)),rtrans(nx)*sin(linspace(0,2*pi,1e3)),'k:','color',state.linecolor(2,:),'linewidth',1,'parent',state.handles.axes);
                        text(0,+1*rtrans(nx+.01),sprintf('%d%%',round(nx*100)),'horizontalalignment','center','color',0*state.linecolor(2,:),'fontsize',round(state.fontsize*.8),'parent',state.handles.axes);
                    end
                else
                    for nx=10.^floor(log10(R)):10.^floor(log10(R)):R
                        plot(rtrans(nx)*cos(linspace(0,2*pi,1e3)),rtrans(nx)*sin(linspace(0,2*pi,1e3)),'k:','color',state.linecolor(2,:),'linewidth',1,'parent',state.handles.axes);
                        text(0,+1*rtrans(nx+.01),num2str(round(nx,1)),'horizontalalignment','center','color',0*state.linecolor(2,:),'fontsize',round(state.fontsize*.8),'parent',state.handles.axes);
                    end
                end
                state.handles.text=[];
                if ~isempty(names)
                    for nx=1:Nx,
                        state.handles.text(nx)=text(1.05*rtrans(R)*cos(ang(nx)),1.05*rtrans(R)*sin(ang(nx)),upper(names{nx}),'fontsize',round(state.fontsize*1.2),'color',state.linecolor(1,:),'parent',state.handles.axes);
                        if cos(ang(nx))>=0, set(state.handles.text(nx),'horizontalalignment','left','rotation',ang(nx)/pi*180);
                        else set(state.handles.text(nx),'horizontalalignment','right','rotation',ang(nx)/pi*180+180);
                        end
                    end
                end

                if CUMSUM, order=Ny:-1:1;
                else order=1:Ny;
                end
                for n=order
                    if STYLE==1
                        datax=cos(ang).*rtrans(data(:,n));
                        datay=sin(ang).*rtrans(data(:,n));
                        state.handles.patch(n)=patch(1.0*datax,1.0*datay,'k','parent',state.handles.axes);
                        %set(state.handles.patch(n),'facecolor',color(n,:),'edgecolor',state.linecolor(1,:),'facealpha',.5);
                        %state.handles.line{n}=plot([zeros(1,Nx);datax(1:end-1)'],[zeros(1,Nx);datay(1:end-1)'],'k-','color',state.linecolor(1,:),'linewidth',2,'parent',state.handles.axes);
                        if CUMSUM, set(state.handles.patch(n),'facecolor',.5*color(n,:)+.5*state.background,'edgecolor',state.linecolor(1,:));
                        else set(state.handles.patch(n),'facecolor',color(n,:),'edgecolor',state.linecolor(1,:),'facealpha',.5);
                        end
                    else
                        datax2=rtrans(interp1(1:Nx+1,data(:,n),1:1/32:Nx+1,INTERP)).*cos(interp1(1:Nx+1,ang,1:1/32:Nx+1));
                        datay2=rtrans(interp1(1:Nx+1,data(:,n),1:1/32:Nx+1,INTERP)).*sin(interp1(1:Nx+1,ang,1:1/32:Nx+1));
                        % datax2=reshape(((data(1:end-1,n)*cos((-8:8)*pi/16)).*cos(repmat(ang(1:end-1),1,17)+repmat((2*pi)/Nx*(-8:8)/16,Nx,1)))',[],1);
                        % datay2=reshape(((data(1:end-1,n)*cos((-8:8)*pi/16)).*sin(repmat(ang(1:end-1),1,17)+repmat((2*pi)/Nx*(-8:8)/16,Nx,1)))',[],1);
                        state.handles.patch(n)=patch(1.0*datax2,1.0*datay2,'k','parent',state.handles.axes);
                        if CUMSUM, set(state.handles.patch(n),'facecolor',.5*color(n,:)+.5*state.background,'edgecolor',state.linecolor(1,:));
                        else set(state.handles.patch(n),'facecolor',color(n,:),'edgecolor',state.linecolor(1,:),'facealpha',.5);
                        end
                    end
                    %plot(1.05*datax,1.05*datay,'k-','color',[.65 .65 .65],'linewidth',1,'parent',state.handles.axes);
                    % for nx=1:Nx,
                    %     if datax(nx)~=0||datay(nx)~=0,
                    %         tang=angle(datax(nx)+1i*datay(nx));
                    %         patch(1*datax(nx)+.025*cos(tang+2*pi/3*[-1 0 1 -1]),1*datay(nx)+.025*sin(tang+2*pi/3*[-1 0 1 -1]),'k','facecolor',state.linecolor(1,:),'edgecolor',state.linecolor(1,:),'parent',state.handles.axes);
                    %     end
                    % end
                end
                if DOEDGES
                    if STYLE==1
                        if CUMSUM
                            datax=cos(ang).*rtrans(data(:,Ny));
                            datay=sin(ang).*rtrans(data(:,Ny));
                            plot(1.0*datax,1.0*datay,'k','linewidth',2,'parent',state.handles.axes);
                        else
                            for n=1:Ny
                                datax=cos(ang).*rtrans(data(:,n));
                                datay=sin(ang).*rtrans(data(:,n));
                                plot(1.0*datax,1.0*datay,'k','linewidth',2,'parent',state.handles.axes);
                            end
                        end
                    else
                        if CUMSUM
                            datax2=rtrans(interp1(1:Nx+1,data(:,Ny),1:1/32:Nx+1,INTERP)).*cos(interp1(1:Nx+1,ang,1:1/32:Nx+1));
                            datay2=rtrans(interp1(1:Nx+1,data(:,Ny),1:1/32:Nx+1,INTERP)).*sin(interp1(1:Nx+1,ang,1:1/32:Nx+1));
                            plot(1.0*datax2,1.0*datay2,'k','linewidth',2,'parent',state.handles.axes);
                        else
                            for n=1:Ny
                                datax2=rtrans(interp1(1:Nx+1,data(:,n),1:1/32:Nx+1,INTERP)).*cos(interp1(1:Nx+1,ang,1:1/32:Nx+1));
                                datay2=rtrans(interp1(1:Nx+1,data(:,n),1:1/32:Nx+1,INTERP)).*sin(interp1(1:Nx+1,ang,1:1/32:Nx+1));
                                plot(1.0*datax2,1.0*datay2,'k','linewidth',2,'parent',state.handles.axes);
                            end
                        end
                    end
                end
                hold(state.handles.axes,'off');
                if DOLEGEND&&numel(datanames)>0
                    h=legend(state.handles.patch, datanames(selectedROIs), 'Location', 'SouthWest');
                    if mean(state.background)>.5, set(h,'box','off','color',[1 1 1]);
                    else set(h,'box','on','color',[.25 .25 .25]);
                    end
                    set(h,'fontsize',state.fontsize)
                end
                axis(state.handles.axes,'equal');
                set(state.handles.axes,'xlim',rtrans(R)*1.25*[-1 1],'ylim',rtrans(R)*1.25*[-1 1],'xtick',[],'ytick',[],'xcolor',state.background,'ycolor',state.background,'visible','off');
                set(state.handles.hfig,'color',state.background);

                if ispercentage, state.handles.title=uicontrol('units','norm','position',[0 .95 1 .05],'style','text','string','percentage of voxels in each network','horizontalalignment','center','fontsize',round(state.fontsize*1.2),'fontweight','bold','backgroundcolor',state.background,'parent',state.handles.hfig);
                else state.handles.title=uicontrol('units','norm','position',[0 .95 1 .05],'style','text','string','number of voxels in each network','horizontalalignment','center','fontsize',round(state.fontsize*1.2),'fontweight','bold','backgroundcolor',state.background,'parent',state.handles.hfig);
                end
                if ~isempty(strdescrip)
                    if ~isempty(strdescrip_title), set(state.handles.table_title,'visible','on','string',strdescrip_title,'fontname','monospaced','horizontalalignment','left','fontsize',round(state.fontsize),'backgroundcolor',state.background,'foregroundcolor',state.linecolor(1,:)); end
                    set(state.handles.table,'visible','on','string',strdescrip,'fontname','monospaced','fontsize',round(state.fontsize),'backgroundcolor',state.background,'foregroundcolor',state.linecolor(1,:));
                    if numel(selectedROIs)==size(data0,2), set(state.handles.table,'value',Ny+1); else set(state.handles.table,'value',selectedROIs); end
                end
                
            case 'close', close(state.handles.hfig);
            case 'figurehandle', out=state.handles.hfig;
            case 'print',
                set([state.handles.table, state.handles.table_title],'visible','off');
                conn_print(state.handles.hfig,varargin{:});
                set([state.handles.table, state.handles.table_title],'visible','on');
            case 'printall',
                conn_print(state.handles.hfig,varargin{:});
            case 'export'
                conn_exportlist(state.handles.table,'',get(state.handles.table_title,'string'));
            case 'export_data'
                [filename,filepath]=uiputfile({'*.mat','MAT-files (*.mat)'; '*.txt','text files (*.txt)'; '*.csv','CSV-files (*.csv)'; '*',  'All Files (*)'},'Save data as');
                if ~ischar(filename), return; end
                filename=fullfile(filepath,filename);
                if ~isempty(regexp(filename,'\.mat$')), conn_savematfile(filename,'-struct',struct('data',data0(:,selectedROIs),'names',{datanames(selectedROIs)},'networks',{names}));
                else conn_savetextfile(filename,data0(:,selectedROIs),datanames(selectedROIs));
                end
                fprintf('Data exported to %s\n',filename);
            case 'colormap'
                cmap=varargin{1};
                if ischar(cmap)
                    switch(cmap)
                        case 'normal', cmap=jet(2*96).^repmat(1+0*abs(linspace(1,-1,2*96))',1,3);
                        case 'red', cmap=[linspace(0,1,2*96)',zeros(2*96,2)];
                        case 'hot', cmap=hot(2*96);
                        case 'jet', cmap=fixedge(jet(2*96));
                        case 'gray', cmap=gray(2*96);
                        case 'bone', cmap=bone(2*96);
                        case 'cool',cmap=fixedge(cool(2*96));
                        case 'hsv',cmap=fixedge(hsv(2*96));
                        case 'bluewhitered', cmap=[zeros(1,48) linspace(0,1,48) ones(1,48) linspace(1,.5,48); linspace(0,1,96) linspace(1,0,48) zeros(1,48); linspace(.5,1,48) ones(1,48) linspace(1,0,48) zeros(1,48)]'; cmap=repmat(abs(linspace(-1,1,192)'),1,3).*cmap+(1-repmat(abs(linspace(-1,1,192)'),1,3))*1;
                        case 'parulawhite', cmap=parula(192); cmap=repmat(abs(linspace(-1,1,192)'),1,3).*cmap+(1-repmat(abs(linspace(-1,1,192)'),1,3))*1;
                        case 'spring',cmap=spring(2*96);
                        case 'summer',cmap=summer(2*96);
                        case 'autumn',cmap=autumn(2*96);
                        case 'winter',cmap=winter(2*96);
                        case 'random',cmap=rand(2*96,3);
                        case 'brighter',cmap=min(1,1/sqrt(.95)*state.colormap.^(1/2));
                        case 'darker',cmap=.95*state.colormap.^2; 
                        case 'manual',answer=conn_menu_inputdlg({'colormap (192x3)'},'',1,{mat2str(state.colormap(round(size(state.colormap,1)/2)+1:end,:))});if ~isempty(answer), answer=str2num(answer{1}); end;if size(answer,2)~=3, return; end;cmap=max(0,min(1,answer));
                        case 'color',cmap=uisetcolor([],'Select color'); if isempty(cmap)||isequal(cmap,0), return; end; cmap=repmat(cmap,2*96,1);
                        otherwise, disp('unknown value');
                    end
                end
                if ~isempty(cmap)
                    if size(cmap,2)<3, cmap=cmap(:,min(size(cmap,2),1:3)); end
                    if size(cmap,1)==1, cmap=linspace(0,1,96)'*cmap; end
                    if size(cmap,1)==1*96, cmap=[flipud(cmap(:,[2,3,1]));cmap]; end
                    if size(cmap,1)~=2*96, cmap=cmap(ceil((1:2*96)/(2*96)*size(cmap,1)),:); end
                    state.colormap=cmap;
                end
                fh('refresh');
            case 'background'
                value=varargin{1};
                state.background=value;
                fh('refresh');
            case 'legend'
                onoff=varargin{1};
                if ischar(onoff), DOLEGEND=strcmpi(onoff,'on');
                else DOLEGEND=onoff;
                end
                fh('refresh');
            case 'fontsize'
                opt=varargin{1};
                if isequal(opt,'+'), opt=state.fontsize+1;
                elseif isequal(opt,'-'), opt=state.fontsize-1;
                elseif isequal(opt,'?'), opt=conn_menu_inputdlg('Enter fontsize','conn_polar_display',1,{num2str(state.fontsize)}); if ~isempty(opt), opt=str2num(opt{1}); end; if isempty(opt), return; end
                end
                state.fontsize=opt;
                fh('refresh');
            case 'scaling'
                onoff=varargin{1};
                if ischar(onoff), DOSQRT=~strcmpi(onoff,'linear');
                else DOSQRT=onoff;
                end
                fh('refresh');
            case 'edges'
                onoff=varargin{1};
                if ischar(onoff), DOEDGES=strcmpi(onoff,'on');
                else DOEDGES=onoff;
                end
                fh('refresh');
            case 'select.rois'
                if numel(varargin)>0, s=varargin{1};
                else
                    [s,v] = listdlg('PromptString',['Select ROIs '],...
                        'SelectionMode','multiple',...
                        'ListString',datanames,...
                        'InitialValue',selectedROIs,...
                        'ListSize',[600,300]);
                    if isempty(s), return; end
                end
                selectedROIs=s;
                fh('refresh');
            case 'list'
                value=get(state.handles.table,'value');
                if any(value>numel(datanames)), selectedROIs=1:numel(datanames);
                else selectedROIs=value;
                end
                fh('refresh');
            case 'style',
                switch(lower(varargin{1}))
                    case {'style.stacked','stacked'}, CUMSUM=true;
                    case {'style.layered','layered'}, CUMSUM=false;
                    case 'lines.straight', STYLE=1;
                    case 'lines.curved', STYLE=2;
                    otherwise, error('unrecognized option %s',varargin{1});
                end
                fh('refresh');
            case 'lines',
                switch(lower(varargin{1}))
                    case 'straight', STYLE=1;
                    case 'curved', STYLE=2;
                    otherwise, error('unrecognized option %s',varargin{1});
                end
                fh('refresh');
            case 'getstate',
                out=state;
                out=rmfield(out,'handles');
                return;
        end
    end
end

function c=fixedge(c)
k=linspace(1,0,size(c,1))'.^4;
c=repmat(1-k,1,size(c,2)).*c+repmat(k.*mean(c,2),1,size(c,2));
end
