function fh = conn_polar_display(data, varargin)
% CONN_POLAR_DISPLAY displays polar data
%
% conn_polar_display(filename[, templatefile, ispercentage, manualscale])
% displays the overlap between non-zero values in 3D NIFTI volume "filename"
% and the 7-networks whole-brain parcellation (Yeo&Buckner 7-networks atlas)
%    filename     : 3D volume NIFTI file 
%    templatefile : 3D volume NIFTI template with ROIs/networks
%                   (defaults to conn/utils/surf/YeoBuckner2011.nii)
%    ispercentage : 0/1 0=raw numbers; 1=percentage overlap [0]
%    manualscale  : maximum value displayed
%
% conn_polar_display(data [, names, scaletoone, manualscale])
% displays the numbers in the data matrix data (one row per sample,
% one column per observation/variable) in polar display.
%    data         : Nx1 or NxM data matrix
%    names        : 1xN cell array with labels for each polar direction
%    ispercentage : 0/1 0=raw numbers; 1=percentage overlap [0]
%    manualscale  : maximum value displayed
% 
% e.g.
%   conn_polar_display('results.ROI.nii');
% e.g.
%   conn_polar_display(rand(7,1));
%

STYLE=2; % 1:surface; 2:petals
if ~nargin, help(mfilename); return; end

if numel(varargin)>=1&&~isempty(varargin{1}), filetemplate=varargin{2}; else filetemplate=''; end
if numel(varargin)>=2&&~isempty(varargin{2}), ispercentage=varargin{2}; else ispercentage=[]; end
if numel(varargin)>=3&&~isempty(varargin{3}), manualscale=varargin{3}; else manualscale=[]; end
datanames={};
templateorder=[];
if iscell(data)||ischar(data)
    if isempty(ispercentage), ispercentage=false; end
    filedata=data;
    if ischar(filedata), filedata={filedata}; end
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
        if all(rem(dataval(mask),1)==0), % ROI data file
            datanames=arrayfun(@(n)sprintf('cluster%d',n),1:max(dataval(mask)),'uni',0);
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

state.colormap=jet(2*96).^repmat(1+0*abs(linspace(1,-1,2*96))',1,3); 
state.handles.hfig=figure('units','norm','position',[.4 .25 .4 .5],'color',[1 1 1],'menubar','none','name','polar display','numbertitle','off','colormap',state.colormap);
state.facecolor=[.25 .25 .22; .975 .975 .975]; % in; out;
state.linecolor=[.65 .65 .65; .75 .75 .75]; % in; out

Nx=size(data,1);
Ny=size(data,2);
ang=0*pi/2+(0:Nx)'*(2*pi)/Nx;
data=[data;data(1,:)];
state.handles.patch=[];
state.handles.line={};
if Ny>1, color=[state.colormap(round(linspace(1,size(state.colormap,1),Ny)),:)];
else color=state.facecolor(1,:);
end
state.handles.axes=axes('units','norm','position',[.2 .2 .6 .6]);
hold(state.handles.axes,'on');
if ~isempty(manualscale), R=manualscale; % manually scale
else R=1.25*max(abs(data(:))); % scale set to 25% above maximum value
end
if ispercentage, R=min(1, R); end
patch(R*cos(linspace(0,2*pi,1e3)),R*sin(linspace(0,2*pi,1e3)),'k-','facecolor',state.facecolor(2,:),'edgecolor','none','parent',state.handles.axes);
%plot(R*cos(linspace(0,2*pi,1e3)),R*sin(linspace(0,2*pi,1e3)),'k:','color',state.linecolor(2,:),'linewidth',2,'parent',state.handles.axes);
plot([zeros(1,Nx);R*cos(ang(1:end-1))'],[zeros(1,Nx);R*sin(ang(1:end-1))'],'k-','color',state.linecolor(2,:),'linewidth',1,'parent',state.handles.axes);
for nx=1:Nx,
    patch(.975*R*cos(ang(nx))+.025*R*cos(ang(nx)+2*pi/3*[-1 0 1 -1]),.975*R*sin(ang(nx))+.025*R*sin(ang(nx)+2*pi/3*[-1 0 1 -1]),'k','facecolor',state.linecolor(2,:),'edgecolor',state.linecolor(2,:),'parent',state.handles.axes);
end
if ispercentage
    for nx=.1:.1:R
        plot(nx*cos(linspace(0,2*pi,1e3)),nx*sin(linspace(0,2*pi,1e3)),'k:','color',state.linecolor(2,:),'linewidth',1,'parent',state.handles.axes);
        text(0,-nx-.01,sprintf('%d%%',round(nx*100)),'horizontalalignment','center','color',0*state.linecolor(2,:),'fontsize',8,'parent',state.handles.axes);
    end
else
    for nx=10.^floor(log10(R)):10.^floor(log10(R)):R
        plot(nx*cos(linspace(0,2*pi,1e3)),nx*sin(linspace(0,2*pi,1e3)),'k:','color',state.linecolor(2,:),'linewidth',1,'parent',state.handles.axes);
        text(0,-nx-.01,num2str(round(nx,1)),'horizontalalignment','center','color',0*state.linecolor(2,:),'fontsize',8,'parent',state.handles.axes);
    end
end
state.handles.text=[];
if ~isempty(names)
    for nx=1:Nx,
        state.handles.text(nx)=text(1.05*R*cos(ang(nx)),1.05*R*sin(ang(nx)),upper(names{nx}),'fontsize',12,'color',min(state.linecolor(1,:),state.facecolor(1,:)),'parent',state.handles.axes);
        if cos(ang(nx))>=0, set(state.handles.text(nx),'horizontalalignment','left','rotation',ang(nx)/pi*180);
        else set(state.handles.text(nx),'horizontalalignment','right','rotation',ang(nx)/pi*180+180);
        end
    end
end

for n=1:Ny
    datax=cos(ang).*data(:,n);
    datay=sin(ang).*data(:,n);
    if STYLE==1
        state.handles.patch(n)=patch(1.0*datax,1.0*datay,'k','parent',state.handles.axes);
        set(state.handles.patch(n),'facecolor',color(n,:),'edgecolor',state.linecolor(1,:),'facealpha',.5);
        state.handles.line{n}=plot([zeros(1,Nx);datax(1:end-1)'],[zeros(1,Nx);datay(1:end-1)'],'k-','color',state.linecolor(1,:),'linewidth',2,'parent',state.handles.axes);
    else
        datax2=interp1(1:Nx+1,data(:,n),1:1/32:Nx+1,'pchip').*cos(interp1(1:Nx+1,ang,1:1/32:Nx+1));
        datay2=interp1(1:Nx+1,data(:,n),1:1/32:Nx+1,'pchip').*sin(interp1(1:Nx+1,ang,1:1/32:Nx+1));
        % datax2=reshape(((data(1:end-1,n)*cos((-8:8)*pi/16)).*cos(repmat(ang(1:end-1),1,17)+repmat((2*pi)/Nx*(-8:8)/16,Nx,1)))',[],1);
        % datay2=reshape(((data(1:end-1,n)*cos((-8:8)*pi/16)).*sin(repmat(ang(1:end-1),1,17)+repmat((2*pi)/Nx*(-8:8)/16,Nx,1)))',[],1);
        state.handles.patch(n)=patch(1.0*datax2,1.0*datay2,'k','parent',state.handles.axes);
        set(state.handles.patch(n),'facecolor',color(n,:),'edgecolor',state.linecolor(1,:),'facealpha',.5);
    end
    %plot(1.05*datax,1.05*datay,'k-','color',[.65 .65 .65],'linewidth',1,'parent',state.handles.axes);
    % for nx=1:Nx, 
    %     if datax(nx)~=0||datay(nx)~=0,
    %         tang=angle(datax(nx)+1i*datay(nx));
    %         patch(1*datax(nx)+.025*cos(tang+2*pi/3*[-1 0 1 -1]),1*datay(nx)+.025*sin(tang+2*pi/3*[-1 0 1 -1]),'k','facecolor',state.linecolor(1,:),'edgecolor',state.linecolor(1,:),'parent',state.handles.axes);
    %     end
    % end
end
if Ny>1&&numel(datanames)==Ny
    h=legend(state.handles.patch, datanames, 'Location', 'NorthEastOutside');
    set(h,'box','off');
end

hold(state.handles.axes,'off');
axis(state.handles.axes,'equal','tight');
set(state.handles.axes,'xtick',[],'ytick',[],'xcolor','w','ycolor','w');
if ispercentage, uicontrol('units','norm','position',[0 .95 1 .05],'style','text','string','Percentage of voxels in each network','horizontalalignment','center','fontsize',12,'fontweight','bold','parent',state.handles.hfig);
else uicontrol('units','norm','position',[0 .95 1 .05],'style','text','string','Number of voxels in each network','horizontalalignment','center','fontsize',12,'fontweight','bold','parent',state.handles.hfig);
end