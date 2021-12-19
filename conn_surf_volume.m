function p=conn_surf_volume(filename,DODISP,THR,SMOOTH,DISREGARDZEROS,HEMSEPARATE,SIGNSHIFT)
% CONN_SURF_VOLUME computes iso-surface from 3d volume data
% internal function

if nargin<1||isempty(filename)
    [file_name,file_path]=uigetfile('*.nii;*.img','Select volume',pwd);
    if isequal(file_name,0), return; end
    filename=fullfile(file_path,file_name);
end
DATA=[];
XYZ=[];
VOL=[];
if isstruct(filename), VOL=filename; fname=VOL.fname;
elseif ischar(filename), 
    [file_path,file_name,file_ext]=fileparts(filename);
    singlefile=isequal(file_name,'lh')|isequal(file_name,'rh')|strncmp(file_name,'lh.',3)|strncmp(file_name,'rh.',3); 
    if singlefile
        p=[conn_surf_readsurf(fullfile(file_path,[regexprep(file_name,'^(lh|rh)(\.)?','lh$2'),file_ext])),conn_surf_readsurf(fullfile(file_path,[regexprep(file_name,'^(lh|rh)(\.)?','rh$2'),file_ext]))];
        return
    else
        VOL=conn_fileutils('spm_vol',filename); fname=filename;
    end
elseif iscell(filename)
    DATA=filename{1};
    XYZ=filename{2};
    fname=[];
end
if isempty(DATA), DATA=conn_fileutils('spm_read_vols',VOL); end
if nargin<2||isempty(DODISP), DODISP=true; end
if nargin<3||isempty(THR), THR=0; end %median(DATA(~isnan(DATA)&DATA~=0)); end
if nargin<4||isempty(SMOOTH), SMOOTH=10; end
if nargin<5||isempty(DISREGARDZEROS), DISREGARDZEROS=true; end
if nargin<6||isempty(HEMSEPARATE), HEMSEPARATE=false; end
if nargin<7||isempty(SIGNSHIFT), SIGNSHIFT=false; end
if SIGNSHIFT, DATA=-DATA; end
%if size(DATA,4)>1, DATA=max(DATA,[],4); end
if numel(VOL)>1, VOL=VOL(1); end
refnames={}; refidx=[];
if ~isempty(fname)&&~any(rem(DATA(:),1))
    try
        [trefnames,trefidx]=conn_roilabels(fname);
        if ~isempty(trefnames)
            maxdata=max(DATA(:));
            if max(trefidx)==maxdata, refnames=trefnames; refidx=trefidx; reftype=1; if DODISP, THR=refidx(1); else THR=trefidx; end
            elseif numel(trefnames)==size(DATA,4), refnames=trefnames; refidx=trefidx; reftype=2; if DODISP, THR=refidx(1); else THR=trefidx; end
            end
        end
    end
end
DOLOGHIST=false;
VIEW=[1 1 1];
OLDTHR=nan;
OLDSMOOTH=nan;
OLDDISREGARDZEROS=nan;
OLDHEMSEPARATE=nan;
HANDLES=[];
p=[];
p1=[];
p2=[];
A={};
HISTa=[];
HISTb=[];

if isempty(XYZ)
    [x,y,z]=ndgrid(1:VOL.dim(1),1:VOL.dim(2),1:VOL.dim(3));
    XYZ=reshape([x(:),y(:),z(:),ones(numel(x),1)]*VOL.mat',size(x,1),size(x,2),size(x,3),[]);
end
hfig=figure('menubar','none','color','w','numbertitle','off','name','','units','norm','position',[.3,.4,.5,.4]);
conn_surf_volume_update;
if DODISP, uiwait(hfig); end
if ishandle(hfig)
    p=p2;
    delete(hfig);
    if DODISP, drawnow; end
end

    function h=conn_surf_volume_update(varargin)
        if numel(HANDLES)>=1, 
            if isempty(refnames), THR=str2num(get(HANDLES(1),'string')); 
            else THR=refidx(get(HANDLES(1),'value'));
            end
        end
        if numel(HANDLES)>=2, SMOOTH=str2num(get(HANDLES(2),'string')); end
        if numel(HANDLES)>=3, DISREGARDZEROS=get(HANDLES(3),'value'); end
        if numel(HANDLES)>=4, HEMSEPARATE=get(HANDLES(4),'value'); end
        clf(hfig);
        bg=.9*[1 1 1];
        uicontrol('style','frame','units','norm','position',[0,.83,1,.17],'foregroundcolor',bg,'backgroundcolor',bg,'parent',hfig);
        if isempty(refnames), uicontrol('style','text','units','norm','position',[.05,.90,.3,.05],'string','Intensity threshold','horizontalalignment','right','backgroundcolor',bg,'parent',hfig); end
        uicontrol('style','text','units','norm','position',[.05,.85,.3,.05],'string','Smoothing level','horizontalalignment','right','backgroundcolor',bg,'parent',hfig);
        if isempty(refnames), HANDLES(1)=uicontrol('style','edit','units','norm','position',[.36,.90,.1,.05],'string',mat2str(THR),'tooltipstring','<HTML>Select voxels based on intensity values <br/> - Entering v selects voxels with I>=v<br/> - Entering v w selects voxels with I>=v & I<=w<br/> - Entering multiple pairs v1 w1 ... vn wn selects voxels with (I>=v1 & I<=w1) | ... | (I>=vn & I<=wn)</HTML>','callback',@conn_surf_volume_update,'parent',hfig); 
        else HANDLES(1)=uicontrol('style','listbox','units','norm','position',[.05,.30,.4,.40],'string',refnames,'max',2,'value',find(ismember(refidx,THR)),'tooltipstring','<HTML>Select ROIs</HTML>','callback',@conn_surf_volume_update,'parent',hfig); 
            HANDLES(5)=uicontrol('style','pushbutton','units','norm','position',[.05,.25,.4,.05],'string','Select all','value',0,'callback','h=get(gcbo,''userdata''); set(h,''value'',1:numel(get(h,''string''))); feval(get(h,''callback''));','userdata',HANDLES(1),'parent',hfig);
        end
        HANDLES(2)=uicontrol('style','edit','units','norm','position',[.36,.85,.1,.05],'string',num2str(SMOOTH),'backgroundcolor',bg,'tooltipstring','Smooth suprathreshold surfaces (enter number of diffusion iterations)','callback',@conn_surf_volume_update,'parent',hfig);
        HANDLES(3)=uicontrol('style','checkbox','units','norm','position',[.55,.90,.3,.05],'string','Disregard Intensity=0 values','value',DISREGARDZEROS,'backgroundcolor',bg,'callback',@conn_surf_volume_update,'parent',hfig);
        if ~isempty(refnames), set(HANDLES(3),'visible','off'); end
        HANDLES(4)=uicontrol('style','checkbox','units','norm','position',[.55,.85,.3,.05],'string','Separate by hemisphere','value',HEMSEPARATE,'backgroundcolor',bg,'callback',@conn_surf_volume_update,'parent',hfig);
        %uicontrol('style','frame','units','norm','position',[0,.0,1,.15],'foregroundcolor','w','parent',hfig);
        uicontrol('style','pushbutton','units','norm','position',[.65,.03,.15,.09],'string','Cancel','callback','delete(gcbf)','parent',hfig);
        uicontrol('style','pushbutton','units','norm','position',[.80,.03,.15,.09],'string','Ok','callback','uiresume(gcbf)','parent',hfig);
        try
            if DODISP, set(hfig,'pointer','watch');pause(.001); hmsg=conn_msgbox({'Computing contours','Please wait...'}); end
            if ~isequal(DISREGARDZEROS,OLDDISREGARDZEROS)||~isequal(THR,OLDTHR)||~isequal(HEMSEPARATE,OLDHEMSEPARATE)
                p1=[];
                if isempty(refnames), validvolumes=1:size(DATA,4); thr=THR;
                elseif reftype==1, validvolumes=1:size(DATA,4); thr=reshape(repmat(THR(:)',2,1),1,[]);
                else validvolumes=THR; thr=0;
                end
                for m=validvolumes
                    dat=DATA(:,:,:,m);
                    if DISREGARDZEROS, dat(dat==0)=nan; end
                    if HEMSEPARATE
                        if numel(thr)==1
                            p1a=isosurface(XYZ(:,:,:,1),XYZ(:,:,:,2),XYZ(:,:,:,3),dat>=thr&XYZ(:,:,:,1)<=0,.5);
                            p1b=isosurface(XYZ(:,:,:,1),XYZ(:,:,:,2),XYZ(:,:,:,3),dat>=thr&XYZ(:,:,:,1)>=0,.5);
                            if isfield(p1a,'vertices')&&~isempty(p1a.vertices), p1=[p1 p1a]; end
                            if isfield(p1b,'vertices')&&~isempty(p1b.vertices), p1=[p1 p1b]; end
                        else
                            %tdat=0; for nthr=1:2:numel(thr), tdat=tdat|(dat>=thr(nthr)&dat<=thr(nthr+1)); end
                            for nthr=1:2:numel(thr)
                                tdat=dat>=thr(nthr)&dat<=thr(nthr+1);
                                p1a=isosurface(XYZ(:,:,:,1),XYZ(:,:,:,2),XYZ(:,:,:,3),tdat&XYZ(:,:,:,1)<=0,.5);
                                p1b=isosurface(XYZ(:,:,:,1),XYZ(:,:,:,2),XYZ(:,:,:,3),tdat&XYZ(:,:,:,1)>=0,.5);
                                if isfield(p1a,'vertices')&&~isempty(p1a.vertices), p1=[p1 p1a]; end
                                if isfield(p1b,'vertices')&&~isempty(p1b.vertices), p1=[p1 p1b]; end
                            end
                        end
                    else
                        if numel(thr)==1
                            p1a=isosurface(XYZ(:,:,:,1),XYZ(:,:,:,2),XYZ(:,:,:,3),dat>=thr,.5);
                            p1=[p1,p1a];
                        else
                            %tdat=0; for nthr=1:2:numel(thr), tdat=tdat|(dat>=thr(nthr)&dat<=thr(nthr+1)); end
                            for nthr=1:2:numel(thr), 
                                tdat=dat>=thr(nthr)&dat<=thr(nthr+1); 
                                p1a=isosurface(XYZ(:,:,:,1),XYZ(:,:,:,2),XYZ(:,:,:,3),tdat,.5);
                                p1=[p1,p1a];
                            end
                        end
                    end
                end
                for m=1:numel(p1),
                    A{m}=double(sparse(repmat(p1(m).faces,3,1),repmat(p1(m).faces,1,3), 1)>0);
                    A{m}=sparse(1:size(A{m},1),1:size(A{m},1),1./sum(A{m},2))*A{m};
                end
                dat=DATA;
                if DISREGARDZEROS, dat(dat==0)=nan; end
                minmaxdata=[min(dat(:)) max(dat(:))]*[1.1  -.1; -.1 1.1];
                [HISTa,HISTb]=hist(dat(~isnan(dat)),linspace(minmaxdata(1),minmaxdata(2),100));
            end
            if ~isequal(DISREGARDZEROS,OLDDISREGARDZEROS)||~isequal(THR,OLDTHR)||~isequal(SMOOTH,OLDSMOOTH)
                p2=p1;
                for m=1:numel(p2),for n=1:SMOOTH,p2(m).vertices=A{m}*p2(m).vertices; end;end
                OLDTHR=THR;
                OLDSMOOTH=SMOOTH;
            end
            
            if DODISP
                set(hfig,'name','conn_surf_volume: Compute isosurface','numbertitle','off','color','w');
                if isempty(refnames)
                    h1=subplot(1,2,1,'parent',hfig);
                    cla(h1);
                    patch([HISTb(1),HISTb(:)',HISTb(end)],[DOLOGHIST,HISTa(:)',DOLOGHIST],'b','parent',h1);
                    hold(h1,'on'); plot([1;1]*THR,repmat([DOLOGHIST;max(HISTa)],1,numel(THR)),'y-','parent',h1); hold(h1,'off');
                    xlabel(h1,'Intensity values');
                    set(h1,'units','norm','position',[.1,.4,.35,.3],'color','w','xcolor',.5*[1 1 1],'ycolor',.5*[1 1 1],'box','off','handlevisibility','off','xlim',[HISTb(1)-.5 HISTb(end)+.5],'ylim',[DOLOGHIST,max(DOLOGHIST+1e-4,max(HISTa))]);
                    if DOLOGHIST, set(h1,'yscale','log'); end
                end
                
                h2=subplot(1,2,2,'parent',hfig);
                for m=1:numel(p2),
                    patch(p2(m),'edgecolor','none','facevertexcdata',repmat([1 1 1],size(p2(m).vertices,1),1),...
                        'facecolor','interp','alphadatamapping','none','FaceLighting', 'gouraud','parent',h2);
                end
                axis(h2,'equal');
                view(h2,VIEW);
                %axis tight;
                light('position',[1 0 0],'parent',h2);
                light('position',[-1 0 0],'parent',h2);
                set(h2,'units','norm','position',[.55,.25,.4,.5],'color','w','xcolor',.5*[1 1 1],'ycolor',.5*[1 1 1],'zcolor',.5*[1 1 1]);
                xlabel(h2,'x'); ylabel(h2,'y'); zlabel(h2,'z');
                grid(h2,'on');
                title(h2,'Isosurface','color',.5*[1 1 1]);
                set(rotate3d(h2),'ActionPostCallback',@conn_surf_volume_rotate,'enable','on');

                drawnow;
                delete(hmsg(ishandle(hmsg)));
                set(hfig,'pointer','arrow');pause(.001);
            end
        catch
            conn_disp(lasterr);
        end
    end

    function conn_surf_volume_rotate(varargin)
        pos=get(gca,'cameraposition');
        VIEW=pos./max(eps,sqrt(sum(pos.^2)));
    end
end