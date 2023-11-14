function fh=conn_mesh_display(filenameSURF,filenameVOL,FSfolder,sphplots,connplots,thr,facealpha,position,defaultcolors,defaultfilepath,Vrange,domask,dosub)
% CONN_MESH_DISPLAY surface display in CONN
%
% CONN_MESH_DISPLAY(fileSURF) displays surface- or volume-level data in fileSURF projected to reference cortical surface
%      fileSURF         : surface- or volume- level 3D NIFTI file (values in this file will be projected and displayed on the 
%                       reference cortical surface -surfaces derived from the ICBM MNI 2009b NLIN asymmetric template-)
%
% CONN_MESH_DISPLAY(fileSURF [,fileVOL, dirFS, rois, connections, threshold]) 
%      fileVOL          : volume-level 3D NIFTI mask file (display clusters of non-zero values as custom masks) (default [])
%      dirFS            : directory containing reference freesurfer-generated surfaces (leave empty for default conn/utils/surf/)
%      rois             : structure defining ROIs to be displayed (default [])
%                           rois.sph_xyz    : (nx3) coordinates XYZ values
%                           rois.sph_c      : (nx3) color RGB values
%                           rois.sph_r      : (nx1) sphere radius
%                           rois.sph_shapes : (optional) patch structure with .vertices and .faces fields
%      connections      : (nxn) matrix of ROI-to-ROI connections to be displayed (values represent connection
%                         strength, with 0 or NaN for connections not to be displayed) (default [])
%      threshold        : only display voxels in 'fileSURF' above 'threshold' (default 0)
%                           set to NaN for automatic threshold (90% percentile absolute values)
%
% fh = CONN_MESH_DISPLAY(...) returns function handle for additional options 
% note: if a conn_mesh_display window is already open, its function handle can also be read using the syntax fh = gcf().UserData.handles.fh;
%
% VIEW OPTIONS
%  fh('view',viewdir [,viewpos,side])    : controls camera view
%                                           viewdir : camera direction (e.g. [1 0 0])
%                                           viewpos : camera-up vector []
%                                           side    : 0 both hemispheres; -1 left-hemisphere only; +1 right-hemisphere only
%  fh('copy')                            : copies current camera view
%  fh('paste')                           : paste camera view
%  fh('zoomin')                          : zoom-in display
%  fh('zoomout')                         : zoom-out display
%  fh('info')                            : shows details of fileSURF/fileVOL data being displayed 
%
% REFERENCE SURFACE OPTIONS
%   fh('brain',type);                    : chooses reference brain surface
%                                           type=1: White Matter
%                                           type=2: Grey Matter 
%                                           type=3: Semi-inflated WM
%                                           type=4: Inflated WM
%  fh('brain_transparency',val)          : set reference brain surface transparency level (val: 0-1)
%  fh('brain_color',color)               : set reference brain surface color (color: [1x3] RGB values) 
%  fh('repaint',fileSURF [,threshold])   : change reference brain surface activation file (fileSURF: NIFTI file)
%  fh('mask',state)                      : displays reference brain surface medial mask (state: 'on' 'off')
%  fh('act_transparency',val)            : set reference brain surface activation transparency level (val: 0-1)
%  fh('act_pos')                         : reference brain surface activation displays only positive acitivation values
%  fh('act_neg')                         : reference brain surface activation displays only negative acitivation values
%  fh('act_posneg')                      : reference brain surface activation displays both positive/negative acitivation values
%  fh('act_color',color)                 : reference brain surface activation colormap (color: [Mx3] RGB values)
%  fh('sub_transparency',val)            : set subcortical surface transparency level (val: 0-1)
%  fh('sub_color',color)                 : set subcortical surface color (color: [1x3] RGB values)
%  fh('brainmask_transparency',val)      : set brainmask surface transparency level (val: 0-1)
%  fh('brainmask_color',color)           : brainmask surface activation color (color: [1x3] RGB values)
%  fh('ref_sag',state)                   : displays reference sagittal slice (state: 'on' 'off')
%  fh('ref_cor',state)                   : displays reference coronal slice (state: 'on' 'off')
%  fh('ref_ax',state)                    : displays reference axial slice (state: 'on' 'off')
%  fh('ref_transparency',val)            : set reference slices transparency level (val: 0-1)
%  fh('ref_pos',pos)                     : set reference slices position (pos: [1x3] vector of x/y/z coordinates)
%  fh('ref_file',file)                   : set reference slices file (default: referenceT1_icbm)
%  fh('ref_method',val)                  : set reference slices display method 
%                                           val=1: sample values at each individual slice; 
%                                           val=2: Maximum Intensity Projection for each slide direction (default 1)
%  fh('ud_select',filename)              : adds custom reference surface (filename: NIFTI mask file)
%  fh('ud_delete')                       : removes last custom reference surface
%  fh('ud_transparency',val)             : set custom reference surface transparency level (val: 0-1)
%  fh('ud_color',color)                  : set custom reference surface color (color: [1x3] RGB values)
%
% ROI OPTIONS
%  fh('roi_transparency',val)            : set ROIs surface transparency level (val: 0-1)
%  fh('roi_color',color)                 : set ROIs surface color (color: [Mx3] RGB values)
%  fh('roi_shape',type)                  : set ROI shape (type: 'real' 'sphere' 'cube' 'star')
%  fh('spheres',scale)                   : increases/decreases all ROIs size (scale 0-inf)
%  fh('labels',state)                    : displays ROI labels (state: 'on' 'off')
%  fh('labels_font',scale)               : increases/decreases all ROI-labels font (scale 0-inf)
%  fh('labels_edit',labels)              : changes ROI labels (labels: cellstring array)
%  fh('con_transparency',val)            : set connectivity-lines transparency level (val: 0-1)
%  fh('con_color',color)                 : set connectivity-lines color (color: [1x3] RGB values)
%  fh('con_bundling',level)              : set connectivity-lines bundling level (level: integers 0-inf)
%  fh('con_width',scale)                 : increases/decreases all connectivity-lines width (scale 0-inf)
%
% EFFECTS OPTIONS
%  fh('material',option)                 : sets material reflectance mode (option: 'dull' 'shiny' 'metal' [1x5]; see "help material")
%  fh('light',level)                     : sets lighting strength (level: 0-1)
%  fh('background',color)                : sets background color (color: [1x3] RGB values)
%  fh('colorbar',state [, title])        : displays reference colorbar (state: 'on' 'off')
%  fh('colorbar','rescale',lim)          : changes colorbar limits (lim: [1x2] values in fileSURF NIFTI file)
%  fh('colorbar','rescale','symmetric')  : changes colorbar limits to force symmetric negative/positive limits 
%  fh('colormap',type)                   : changes colormap (type: 'normal','red','jet','hot','gray','bone','cool','hsv','spring',
%                                           'summer','autumn','winter','bluewhitered','random','brighter','darker','manual','color')
%  fh('smoother',state)                  : 'on' for smoother display (interpolates values in fileSURF NIFTI file when displaying on 
%                                           surface); 'off' for raw-data display
%  fh('axis',state)                      : display reference axis lines (state: 'on' 'off')
%  fh('menubar',state)                   : display menubar (state: 'on' 'off')
%
% PRINT OPTIONS
%  fh('print',type,filename)             : prints display to high-resolution .jpg file
%                                           type=1: prints current view
%                                           type=2: prints 2-view row display (left&right exterior views)
%                                           type=3: prints 3-view mosaic display (superior,right,posterior views; e.g. glass-displays)
%                                           type=4: prints 4-view mosaic display (left, right, medial-left, medial-right views)
%                                           type=5: prints 4-view column display (same as above all views in a single column)
%                                           type=6: prints 4-view row display (same as above all views in a single row)
%                                           type=7: prints 8-view mosaic display (left, right, medial-left, medial-right, superior, 
%                                                   inferior, anterior, posterior views)
%


global CONN_x CONN_gui;
persistent data;
if isempty(CONN_gui)||~isfield(CONN_gui,'font_offset'), conn_font_init; end

doinit=true;
if nargin==1&&isstruct(filenameSURF)&&isfield(filenameSURF,'structural'), % struct from getstate
    state=filenameSURF;
    filenameSURF={};
    if ~conn_fileutils('isdir',state.defaultfilepath), state.defaultfilepath=pwd; end
    if ~conn_fileutils('isdir',state.FSfolder), state.FSfolder=fullfile(fileparts(which('conn')),'utils','surf'); end
    if ~isfield(state,'fontclose'), state.fontclose=1; end
    if ~isfield(state,'Prange'), state.Prange=[]; end
    doinit=false;
elseif nargin>0&&isstruct(filenameSURF),data=filenameSURF;return
end
if isfield(CONN_gui,'slice_display_skipbookmarkicons'), SKIPBI=CONN_gui.slice_display_skipbookmarkicons; 
else SKIPBI=false;
end
fh=@(varargin)conn_mesh_display_refresh([],[],varargin{:});

%if nargin<1||(isempty(filenameSURF)&&~ischar(filenameSURF)&&~iscell(filenameSURF)), filenameSURF=spm_select(1,'image','Select a file'); end
%if nargin<2||(isempty(filenameVOL)&&~ischar(filenameVOL)), filenameVOL=spm_select(1,'image','Select a file'); end
if doinit
    if nargin<1, filenameSURF=''; end
    if nargin<2, filenameVOL=''; end
    if nargin<3||isempty(FSfolder), FSfolder=fullfile(fileparts(which('conn')),'utils','surf'); end
    if nargin<4, sphplots=[]; end
    if nargin<5, connplots=[]; end
    if nargin<6||isempty(thr), thr=0; end
    if nargin<7||isempty(facealpha), facealpha=1; end
    if nargin<8||isempty(position), position=[-1,0,0];end
    if nargin<9||isempty(defaultcolors), defaultcolors={[1,.1,.1],[.1,.1,1]}; end
    if nargin<10||isempty(defaultfilepath), defaultfilepath=pwd; end
    if nargin<11, Vrange=[]; end
    if nargin<12||isempty(domask), domask=true; end
    if nargin<13||isempty(dosub), dosub=isempty(filenameSURF); end
    state.actthr=thr;
    state.FSfolder=FSfolder;
    state.connplots=connplots;
    state.sphplots=sphplots;
    state.facealpha=facealpha;
    state.position=position;
    state.defaultcolors=defaultcolors;
    state.defaultfilepath=defaultfilepath;
    state.Vrange=Vrange;
    state.Prange=[];
    state.domask=domask;
    state.dosub=dosub;
end
if isfield(state,'actthr'), THR=state.actthr;
else THR=0;
end
emph=1;
if ~isempty(filenameSURF)&&ischar(filenameSURF), filenameSURF=conn_mesh_display_expand(filenameSURF); end
if ~isfield(CONN_x,'folders')||~isfield(CONN_x.folders,'bookmarks')||isempty(CONN_x.folders.bookmarks), state.dobookmarks=[];
else state.dobookmarks=CONN_x.folders.bookmarks;
end

hmsg=conn_msgbox('Initializing. Please wait...','',-1);
if ~isfield(CONN_gui,'refs'), conn initrefs; end
if isempty(data), 
    files={'white','cortex','pial.smoothed','inflated','subcortical'};
    for nfiles=1:numel(files)
        data.rend{nfiles}=conn_surf_readsurf(fullfile(state.FSfolder,[files{nfiles},'.surf']));
    end
    %data.rend{end+1}=
    data.curv{3}{1}=sign(conn_freesurfer_read_curv(fullfile(state.FSfolder,'lh.curv.paint')));
    data.curv{3}{2}=sign(conn_freesurfer_read_curv(fullfile(state.FSfolder,'rh.curv.paint')));
    data.curv{3}=cell2mat(data.curv{3});
    data.curv{1}=data.curv{3}; %zeros(size(data.curv{3}));
    data.curv{2}=zeros(size(data.curv{3}));
    data.curv{4}=2*data.curv{3};
    a=spm_vol(fullfile(fileparts(which(mfilename)),'utils','surf','mask.surface.brainmask.nii'));
    data.mask=spm_read_vols(a);
    data.mask=reshape(data.mask,[],2);
    filename=fullfile(fileparts(which(mfilename)),'utils','surf','mask.volume.brainmask.nii');
    %if CONN_x.Setup.analysismask==1&&isfield(CONN_x.Setup,'explicitmask'), filename=CONN_x.Setup.explicitmask{1}; end
    if conn_existfile(conn_prepend('',filename,'.patch'))
        load(conn_prepend('',filename,'.patch'),'patchs','-mat');
    else
        patchs=conn_surf_volume(filename,false,.5,10,true,true);
        try
            save(conn_prepend('',filename,'.patch'),'patchs');
        end
    end
    data.masksurface=patchs;
    %for n=1:numel(data.masksurface),data.masksurface(n).faces=fliplr(data.masksurface(n).faces); end
end

if doinit
    state.selectedsurface=3;
    if ~isempty(filenameSURF)
        if isstruct(filenameSURF{1})
            a=filenameSURF{1}.vol;
            state.info.surf=filenameSURF{1}.filename;
            filepath=fileparts(filenameSURF{1}.filename);
        else
            a=conn_fileutils('spm_vol',char(filenameSURF));
            state.info.surf=char(filenameSURF);
            filepath=fileparts(filenameSURF{1});
        end
        if ~isempty(filepath), state.defaultfilepath=filepath; end
        V=conn_mesh_display_getV(filenameSURF,a,state.FSfolder);
        state.Vrange=[min(V(V~=0)) max(V(V~=0))];
        state.reducedpatch=1;
    else
        state.info.surf='none';
        V=zeros(size(data.curv{1}));
        state.reducedpatch=1;
    end
    if ~isempty(filenameVOL)
        if iscell(filenameVOL),
            otherVOL=filenameVOL(2:end);
            filenameVOL=filenameVOL{1};
        else otherVOL={};
        end
        state.info.vol=char(filenameVOL);
        %if facealpha<=.5, state.selectedsurface=4;
        %else state.selectedsurface=2;
        %end
        state.selectedsurface=2;
        if ~isempty(filenameVOL)
            thr=THR;
            if isnan(thr), thr=0; end
            state.pVOL1=conn_surf_volume(filenameVOL,0,thr,[],1,1,0); 
            state.pVOL2=conn_surf_volume(filenameVOL,0,thr,[],1,1,1); 
            filepath=fileparts(filenameVOL);
            if ~isempty(filepath), state.defaultfilepath=filepath; end
        else
            state.pVOL1=[];
            state.pVOL2=[];
        end
        if facealpha>=.5, state.reducedpatch=1; end
    else
        state.pVOL1=[];
        state.pVOL2=[];
        otherVOL={};
        state.info.vol='none';
    end
    if ~isempty(otherVOL), state.pVOLother=conn_surf_volume(otherVOL{:});
    else state.pVOLother=[];
    end
    
    cmapstart=repmat(1-linspace(1,0,96)'.^8,[1,3]).*autumn(96);
    cmapstart=[flipud(cmapstart(:,[2,3,1]));cmapstart];
    %cmapstart=linspace(-1,1,128)'; 
    %cmapstart=cat(2,1-(cmapstart<0)-cmapstart.*max(0,cmapstart-.5),max(0,1-2*abs(cmapstart)),1-(cmapstart>0)-cmapstart.*min(0,cmapstart+.5));
    cmapdefault=repmat(1-linspace(1,0,128)'.^2,[1,3]).*hot(128)+repmat(linspace(1,0,128)'.^2,[1,3])*.1;
    cmapdefault=cmapdefault(33:end,:);
    %cmapdefault=cmapdefault([1,34:end],:);
    
    if isnan(THR), 
        sthr=sort(abs(V(~isnan(V))));
        V=V.*(abs(V)>sthr(round(.9*numel(sthr))));
    else
        V=V.*(abs(V)>THR);
    end
    %conn_surf_write('test.surf.nii',V);
    if size(V,3)==2, V=abs(V(:,:,1))-abs(V(:,:,2));
    elseif size(V,3)>2, [nill,V]=max(abs(conn_surf_smooth(V)),[],3); V(nill==0)=0; cmapstart=rand(size(cmapstart)); state.Vrange=[min(V(V~=0)) max(V(V~=0))];
    end
    
    if any(V(:)<0)&&any(V(:)>0),
        state.dotwosided=true;
    else
        state.dotwosided=false;
    end
    state.V0=V;
    state.colormap_default=cmapdefault;
    state.colormap=cmapstart;
    state.cdat0=cellfun(@(x)conn_bsxfun(@times,1-.05*x,shiftdim([.7,.65,.6],-1)),data.curv,'uni',0);%cdat=cellfun(@(x)conn_bsxfun(@times,.75-.04*x,shiftdim([1,.9,1],-1)),data.curv,'uni',0);
    state.smoother=1;
    state.facealpha=facealpha;
    state.facealphasub=min(facealpha,.5);
    if ~state.dosub, state.facealphasub=0; end
    state.facealphablob=1;
    state.facealphamask=0;
    state.facealpharef=1;
    state.facealphaud=1;
    state.facealpharoi=1;
    state.facealphatxt=0;
    state.facealphacon=.5;
    state.brain_color=[];
    state.act_color=[];
    state.sub_color=[];
    state.mask_color=[];
    state.ref_color=[];
    state.ud_color=[];
    state.roi_color=[];
    state.con_color=[];
    state.light_color=.5*[1 1 1];
    state.material='dull';
    state.showaxref=[0 0 0];
    state.fontsize=10;
    state.fontclose=1;
    state.bundling=[];
    state.roishape='sphere';
    state.pointer_mm=[40 -60 -10];
    state.up=[1 1 1]/sqrt(3);
    state.background=[.95 .95 .9];%.14*[1 1 1];%[.2,.6,.7];
    state.bookmark_filename='';
    state.bookmark_descr='';
    conn_mesh_display_refresh([],[],'remap');
end

hax0=findobj(0,'type','axes','tag','conn_mesh_display_axes');
if isempty(hax0)
    state.handles.hfig=figure('numbertitle','off','color',state.background,'units','norm','position',[.3 .4 .4 .4],'menubar','none','render','opengl','name','rendering, please wait...','colormap',state.colormap,'visible','off');
    figname='conn 3d display';
    axes('units','norm','position',[.95 .1 .04 .8]);
    state.handles.fullfigure=true;
else
    if numel(hax0)>1, hax0=hax0(end); end
    state.handles.hfig=hax0; while ~isequal(get(state.handles.hfig,'type'),'figure'), state.handles.hfig=get(state.handles.hfig,'parent'); end
    figure(state.handles.hfig);
    set(state.handles.hfig,'menubar','none'); delete(findobj(state.handles.hfig,'type','uimenu')); 
    figname=get(state.handles.hfig,'name');
    phax0=get(hax0,'position');
    axes('units','norm','position',[phax0(1)+phax0(3) phax0(2) .05*phax0(3) phax0(4)]);
    state.handles.fullfigure=false;
end
if state.dotwosided,
    temp=imagesc(max(0,min(1, ind2rgb(round((size(state.colormap,1)+1)/2+emph*(size(state.colormap,1)-1)/2*linspace(-1,1,128)'),state.colormap))));
    set(gca,'ydir','normal','ytick',[.5,64.5,128.5],'yticklabel',arrayfun(@(x)num2str(x,'%.2f'),state.Vrange,'uni',0),'xtick',[],'box','on','fontsize',8);
elseif any(state.Vrange>0)
    temp=imagesc(max(0,min(1, ind2rgb(round((size(state.colormap,1)+1)/2+emph*(size(state.colormap,1)-1)/2*linspace(0,1,128)'),state.colormap))));
    set(gca,'ydir','normal','ytick',[.5,128.5],'yticklabel',arrayfun(@(x)num2str(x,'%.2f'),state.Vrange,'uni',0),'xtick',[],'box','on','fontsize',8);
elseif any(state.Vrange<0)
    temp=imagesc(max(0,min(1, ind2rgb(round((size(state.colormap,1)+1)/2+emph*(size(state.colormap,1)-1)/2*linspace(-1,0,128)'),state.colormap))));
    set(gca,'ydir','normal','ytick',[.5,128.5],'yticklabel',arrayfun(@(x)num2str(x,'%.2f'),state.Vrange,'uni',0),'xtick',[],'box','on','fontsize',8);
else
    temp=imagesc(max(0,min(1, ind2rgb(round((size(state.colormap,1)+1)/2+emph*(size(state.colormap,1)-1)/2*linspace(-1,1,128)'),state.colormap))));
    set(gca,'ydir','normal','xtick',[],'ytick',[],'box','on','fontsize',8);
end
state.handles.hhide=[];

state.handles.colorbar=[gca temp];
set(state.handles.colorbar,'visible','off');
if isempty(hax0), state.handles.hax=axes('parent',state.handles.hfig);
else state.handles.hax=hax0; axes(state.handles.hax);
end
state.selected_vertices={1:size(data.rend{1}(1).vertices,1), CONN_gui.refs.surf.default2reduced};
state.selected_faces={data.rend{1}(1).faces, CONN_gui.refs.surf.spherereduced.faces};
state.handles.patch(1)=patch(struct('vertices',data.rend{state.selectedsurface}(1).vertices(state.selected_vertices{state.reducedpatch},:),'faces',state.selected_faces{state.reducedpatch}),'parent',state.handles.hax,...
    'facevertexcdata',permute(state.cdat{state.selectedsurface}(state.selected_vertices{state.reducedpatch},1,:),[1 3 2]),'facecolor','interp','edgecolor','none','alphadatamapping','none','FaceLighting', 'phong','facealpha',state.facealpha,'backfacelighting','unlit');
state.handles.patch(2)=patch(struct('vertices',data.rend{state.selectedsurface}(2).vertices(state.selected_vertices{state.reducedpatch},:),'faces',state.selected_faces{state.reducedpatch}),'parent',state.handles.hax,...
    'facevertexcdata',permute(state.cdat{state.selectedsurface}(state.selected_vertices{state.reducedpatch},2,:),[1 3 2]),'facecolor','interp','edgecolor','none','alphadatamapping','none','FaceLighting', 'phong','facealpha',state.facealpha,'backfacelighting','unlit');
if ~isempty(state.brain_color), set(state.handles.patch,'facecolor',state.brain_color); end
set(state.handles.patch,'tag','brain_surface');

state.handles.subpatch(1)=patch(data.rend{end}(1),'facecolor',.6*[1 1 1],'edgecolor','none','alphadatamapping','none','FaceLighting', 'phong','facealpha',state.facealphasub,'backfacelighting','unlit','parent',state.handles.hax);
state.handles.subpatch(2)=patch(data.rend{end}(2),'facecolor',.6*[1 1 1],'edgecolor','none','alphadatamapping','none','FaceLighting', 'phong','facealpha',state.facealphasub,'backfacelighting','unlit','parent',state.handles.hax);
%if ~state.dosub, set(state.handles.subpatch,'visible','off'); end
if ~isempty(state.sub_color), set(state.handles.subpatch,'facecolor',state.sub_color); end
set(state.handles.subpatch,'tag','subcortical_surface');

state.handles.maskpatch(1)=patch(data.masksurface(1),'facecolor',.9*[1 1 1],'edgecolor','none','alphadatamapping','none','FaceLighting', 'phong','facealpha',state.facealphamask,'backfacelighting','unlit','parent',state.handles.hax);
state.handles.maskpatch(2)=patch(data.masksurface(2),'facecolor',.9*[1 1 1],'edgecolor','none','alphadatamapping','none','FaceLighting', 'phong','facealpha',state.facealphamask,'backfacelighting','unlit','parent',state.handles.hax);
if ~isempty(state.mask_color), set(state.handles.maskpatch,'facecolor',state.mask_color); end
set(state.handles.maskpatch,'tag','mask_surface');

state.handles.patchblob1=[]; state.handles.patchblob2=[];
if ~isempty(state.pVOL1) % positive blob
    for npvol=1:2:numel(state.pVOL1), state.handles.patchblob1=[state.handles.patchblob1 patch(state.pVOL1(npvol),'facecolor','r','edgecolor','none','alphadatamapping','none','FaceLighting', 'phong','facealpha',state.facealphablob,'backfacelighting','reverselit','parent',state.handles.hax)]; end
    for npvol=2:2:numel(state.pVOL1), state.handles.patchblob2=[state.handles.patchblob2 patch(state.pVOL1(npvol),'facecolor','r','edgecolor','none','alphadatamapping','none','FaceLighting', 'phong','facealpha',state.facealphablob,'backfacelighting','reverselit','parent',state.handles.hax)]; end
end
state.handles.patchblob3=[]; state.handles.patchblob4=[];
if ~isempty(state.pVOL2) % negative blob
    for npvol=1:2:numel(state.pVOL2), state.handles.patchblob3=[state.handles.patchblob3 patch(state.pVOL2(npvol),'facecolor','b','edgecolor','none','alphadatamapping','none','FaceLighting', 'phong','facealpha',state.facealphablob,'backfacelighting','reverselit','parent',state.handles.hax)]; end
    for npvol=2:2:numel(state.pVOL2), state.handles.patchblob4=[state.handles.patchblob4 patch(state.pVOL2(npvol),'facecolor','b','edgecolor','none','alphadatamapping','none','FaceLighting', 'phong','facealpha',state.facealphablob,'backfacelighting','reverselit','parent',state.handles.hax)]; end
end
if ~isempty(state.act_color), 
    set([state.handles.patchblob1 state.handles.patchblob2],'facecolor',state.act_color); 
    set([state.handles.patchblob3 state.handles.patchblob4],'facecolor',state.act_color(:,[2 3 1])); 
end
set([state.handles.patchblob1 state.handles.patchblob2],'tag','activation_positive_surface');
set([state.handles.patchblob3 state.handles.patchblob4],'tag','activation_negative_surface');
hold(state.handles.hax,'on');
state.handles.leftrightlabel=text([-60 60],55+[0 0],70+[0 0],{'L','R'},'fontsize',12,'color',.5*[1 1 1],'parent',state.handles.hax,'visible','off');
hold(state.handles.hax,'off');
if ~isempty(state.pVOLother)
    state.handles.patchblobother=[]; 
    state.patchblobother_x=[];
    if numel(state.pVOLother)>1, cmap=parula(numel(state.pVOLother));
    else cmap=[0 .45 .74];
    end
    for n1=1:numel(state.pVOLother), 
        state.handles.patchblobother=[state.handles.patchblobother patch(state.pVOLother(n1),'facecolor',cmap(n1,:),'edgecolor','none','alphadatamapping','none','FaceLighting', 'phong','facealpha',state.facealphaud,'backfacelighting','reverselit','parent',state.handles.hax)]; 
        state.patchblobother_x=[state.patchblobother_x mean(state.pVOLother(n1).vertices(:,1))];
    end
    if ~isempty(state.ud_color),
        if size(state.ud_color,1)==1, set(state.handles.patchblobother,'facecolor',state.ud_color);
        else for n1=1:numel(state.handles.patchblobother), set(state.handles.patchblobother(n1),'facecolor',state.ud_color(1+rem(n1-1,size(state.ud_color,1)),:)); end
        end
    else state.ud_color=cmap;
    end
    set(state.handles.patchblobother,'tag','userdefined_surface');
else state.handles.patchblobother=[];state.patchblobother_x=[];
end
grid(state.handles.hax,'on');
th1=xlabel(state.handles.hax,{'X (mm)','left - right'});th2=ylabel(state.handles.hax,{'Y (mm)','posterior - anterior'});th3=zlabel(state.handles.hax,{'Z (mm)','inferior - superior'});set([th1,th2,th3],'fontsize',8);
axis(state.handles.hax,'off');
view(state.handles.hax,state.position);
if state.handles.fullfigure, set(state.handles.hax,'units','norm','position',[0,0,1,1],'color',get(state.handles.hfig,'color'),'xcolor',.5*[1 1 1],'ycolor',.5*[1 1 1],'zcolor',.5*[1 1 1],'box','on'); end
state.handles.light=[light light];set(state.handles.light,'position',state.position,'visible','on','color',state.light_color);set(state.handles.light(end),'position',-state.position);
tgca=state.handles.hax;
%axes(tgca);

state.handles.sphplots=[];
%state.sphplots_xyz=zeros(0,3);
%state.sphplots_r=[];
state.handles.sphplots_txt=[];
state.handles.sphplots_shapes=[];
state.handles.sphplots_shapesok=[];
troifiles=[];
if ~isempty(state.sphplots)
    if ~isstruct(state.sphplots), state.sphplots=struct('sph_xyz',state.sphplots); end
    if ~isfield(state.sphplots,'sph_c'), state.sphplots.sph_c=ones(size(state.sphplots.sph_xyz,1),3); end
    if ~isfield(state.sphplots,'sph_r'), state.sphplots.sph_r=5*ones(size(state.sphplots.sph_xyz,1),1); end
    if ~iscell(state.sphplots.sph_c), state.sphplots.sph_c=num2cell(state.sphplots.sph_c,2); end
    if ~isfield(state.sphplots,'sph_names'), state.sphplots.sph_names=repmat({''},1,size(state.sphplots.sph_xyz,1)); end
    state.sphplots.sph_xyz(isnan(state.sphplots.sph_xyz))=0;
    [x,y,z]=sphere(32);
    f=1;
    switch(state.roishape)
        case 'sphere', f=1;
        case 'cube', f=.25;
        case 'star', f=4;
    end
    x=sign(x).*abs(x).^f;y=sign(y).*abs(y).^f;z=sign(z).*abs(z).^f;
    bgc=1-round(mean(state.background))*[1 1 1];
    hold(state.handles.hax,'on');
    for n1=1:size(state.sphplots.sph_xyz,1),
        cdata=.6+.4*(state.sphplots.sph_xyz(n1,1)*x+state.sphplots.sph_xyz(n1,2)*y+state.sphplots.sph_xyz(n1,3)*z)/sqrt(sum(abs(state.sphplots.sph_xyz(n1,:)).^2));
        %cdata=max(0,min(1,bsxfun(@plus,cdata,shiftdim(state.sphplots.sph_c{n1}(:),-2)) ));
        cdata=conn_bsxfun(@times,cdata,shiftdim(state.sphplots.sph_c{n1}(:),-2));
        state.handles.sphplots(n1)=surf(state.handles.hax,state.sphplots.sph_xyz(n1,1)+state.sphplots.sph_r(n1)*x,state.sphplots.sph_xyz(n1,2)+state.sphplots.sph_r(n1)*y,state.sphplots.sph_xyz(n1,3)+state.sphplots.sph_r(n1)*z,zeros(size(x)),'cdata',cdata,'edgecolor','none','facealpha',state.facealpharoi);
        if isfield(state.sphplots,'sph_names')
            state.handles.sphplots_txt(n1)=text(state.sphplots.sph_xyz(n1,1)+state.up(1)*state.fontclose*state.sphplots.sph_r(n1),state.sphplots.sph_xyz(n1,2)+state.up(2)*state.fontclose*state.sphplots.sph_r(n1),state.sphplots.sph_xyz(n1,3)+state.up(3)*state.fontclose*state.sphplots.sph_r(n1),state.sphplots.sph_names{n1},'fontsize',state.fontsize,'horizontalalignment','center','color',bgc,'interpreter','none','visible','off','parent',state.handles.hax);
        end
        if isfield(state.sphplots,'sph_shapes')
            state.handles.sphplots_shapes(n1)=patch(struct('vertices',zeros(0,3),'faces',zeros(0,3)),'facecolor',state.sphplots.sph_c{n1},'edgecolor','none','facealpha',state.facealpharoi,'visible','off','parent',state.handles.hax);
            state.handles.sphplots_shapesok(n1)=false;
        end
    end
    hold(state.handles.hax,'off');
    if isempty(state.roi_color),
        try, state.roi_color=cell2mat(state.sphplots.sph_c(:)); end
    else
        if size(state.roi_color,1)==1, set([state.handles.sphplots state.handles.sphplots_shapes(ishandle(state.handles.sphplots_shapes))],'facecolor',state.roi_color); 
        else
            for n1=1:numel(state.handles.sphplots),
                set(state.handles.sphplots(n1),'facecolor',state.roi_color(1+rem(n1-1,size(state.roi_color,1)),:)); 
                if numel(state.handles.sphplots_shapes)>=n1&&ishandle(state.handles.sphplots_shapes(n1)), set(state.handles.sphplots_shapes(n1),'facecolor',state.roi_color(1+rem(n1-1,size(state.roi_color,1)),:)); end
            end
        end
    end
    set(state.handles.sphplots,'tag','rois');
end
if all(isnan(state.connplots(:))), state.connplots=[]; end
if ~isempty(state.connplots)
    hold(state.handles.hax,'on');
    [x,y,z]=cylinder([1,1],20);
    xyz=[x(:),y(:),z(:)]';
    state.handles.connplots=zeros(max(size(state.connplots)));
    state.connplotsLine=zeros(max(size(state.connplots)));
    if 1
        done=false(max(size(state.connplots)));
        lines=[]; distl=[]; state.ValueLines=[];
        wp=linspace(0,1,50)'; wp=[wp 1-wp];
        for n1=1:size(state.connplots,1)
            for n2=1:size(state.connplots,2)
                if n1~=n2&&~isnan(state.connplots(n1,n2))&&state.connplots(n1,n2)~=0,
                    if (~done(n2,n1)&&~done(n1,n2)),%&&(size(state.connplots,1)<n2||size(state.connplots,2)<n1||isnan(state.connplots(n2,n1))||abs(state.connplots(n1,n2))>abs(state.connplots(n2,n1))),
                        done(n1,n2)=true;
                        done(n2,n1)=true;
                        xx=state.sphplots.sph_xyz([n1,n2],:);
                        lines=cat(3,lines,wp*xx);
                        distl=cat(3,distl,sqrt(max(eps,sum(diff(xx,1,1).^2,2)))/(50-1));
                        state.ValueLines=cat(3,state.ValueLines,state.connplots(n1,n2));
                    end
                end
            end
        end
        if 0 % ring placeholder 1
            lines=conn_bsxfun(@times,1-.25*sin(linspace(0,pi,size(lines,1))'),lines);
        end
        if isempty(state.bundling), if size(lines,3)<1e3, state.bundling=5; else state.bundling=0; end; end
        if ~isfield(state,'Lines')||isempty(state.Lines), state.Lines=lines; end
        if ~isfield(state,'WidthLines'), state.WidthLines=[]; end
        if state.bundling==0||(1+max([0,state.bundling])>size(state.Lines,4)),
            [state.Lines,state.WidthLines] = conn_menu_bundle(state.Lines,state.WidthLines,state.ValueLines,max([0,state.bundling]));
        end
        nlines=0;
        done=false(max(size(state.connplots)));
        for n1=1:size(state.connplots,1)
            for n2=1:size(state.connplots,2)
                if n1~=n2&&~isnan(state.connplots(n1,n2))&&state.connplots(n1,n2)~=0,
                    if (~done(n2,n1)&&~done(n1,n2)),%&&(size(state.connplots,1)<n2||size(state.connplots,2)<n1||isnan(state.connplots(n2,n1))||abs(state.connplots(n1,n2))>abs(state.connplots(n2,n1))),
                        done(n1,n2)=true;
                        done(n2,n1)=true;
                        nlines=nlines+1;
                        %state.handles.connplots(n1,n2)=plot3(lines(:,1,nlines),lines(:,2,nlines),lines(:,3,nlines),'k-','linewidth',2,'parent',state.handles.hax);
                        [nill,state.handles.connplots(n1,n2)]=conn_mesh_display_patch(state.Lines(:,1,nlines,1+state.bundling),state.Lines(:,2,nlines,1+state.bundling),state.Lines(:,3,nlines,1+state.bundling),state.WidthLines(:,1,nlines,1+state.bundling),'facecolor','none','edgecolor','none','facealpha',state.facealphacon,'backFaceLighting','lit','parent',state.handles.hax);
                        if state.connplots(n1,n2)>0, set(state.handles.connplots(n1,n2),'facecolor',state.defaultcolors{1});
                        else set(state.handles.connplots(n1,n2),'facecolor',state.defaultcolors{2});
                        end
                        state.connplotsLine(n1,n2)=nlines;
                    end
                end
            end
        end
    else
        done=false(max(size(state.connplots)));
        for n1=1:size(state.connplots,1)
            for n2=1:size(state.connplots,2)
                if n1~=n2&&~isnan(state.connplots(n1,n2))&&state.connplots(n1,n2)~=0,
                    if (~done(n2,n1)&&~done(n1,n2))||abs(state.connplots(n1,n2))>abs(state.connplots(n2,n1)),
                        done(n1,n2)=true;
                        xx=state.sphplots.sph_xyz([n1,n2],:);
                        dx=xx(2,:)-xx(1,:);%dx=dx./repmat(max(eps,sqrt(sum(abs(dx).^2,2))),[1,3]);
                        ndx=null(dx(1,:));
                        txyz=reshape(xyz'*[(0+2*abs(state.connplots(n1,n2)))*ndx,dx(1,:)']',2,size(x,2),3);
                        %txyz=conn_bsxfun(@plus,txyz,lines(:,:,nlines));
                        state.handles.connplots(n1,n2)=mesh(state.handles.hax,xx(1,1)+txyz(:,:,1),xx(1,2)+txyz(:,:,2),xx(1,3)+txyz(:,:,3),'edgecolor','none');
                        %state.handles.connplots(n1,n2)=plot3(xx(:,1)+dx(:,1).*9,xx(:,2)+dx(:,2).*9,xx(:,3)+dx(:,3).*9,'r-','linewidth',round(1+1*abs(connplots(n1,n2))));
                        if state.connplots(n1,n2)>0, set(state.handles.connplots(n1,n2),'facecolor',state.defaultcolors{1});
                        else set(state.handles.connplots(n1,n2),'facecolor',state.defaultcolors{2});
                        end
                        set(state.handles.connplots(n1,n2),'userdata',[ndx dx(1,:)'/sqrt(sum(dx(1,:).^2))]);
                    end
                end
            end
        end
    end
    if ~isempty(state.con_color), set(state.handles.connplots(state.handles.connplots~=0),'facecolor',state.con_color); end
    set(state.handles.connplots(state.handles.connplots~=0),'tag','connections');
    hold off;
    if doinit
        t=conn_bsxfun(@plus,sum(state.connplots~=0,2),sum(state.connplots~=0,1));
        [ti,tj]=find(t==max(t(:)),1);
        state.pointer_mm=state.sphplots.sph_xyz(ti,:);
    end
else
    state.handles.connplots=[];%zeros(numel(state.handles.sphplots));
end
if 1 % reference axes
    state.structural=fullfile(fileparts(which('conn')),'utils','surf','referenceT1_icbm.nii');
    state.structural_method=1;
    %state.structural=fullfile(fileparts(which('conn')),'utils','surf','referenceT1_trans.nii');
    if isfield(CONN_gui,'refs')&&isfield(CONN_gui.refs,'canonical')&&isfield(CONN_gui.refs.canonical,'filename')&&~isempty(CONN_gui.refs.canonical.filename)
        if ~isequal(CONN_gui.refs.canonical.filename,fullfile(fileparts(which('conn')),'utils','surf','referenceT1_trans.nii')), %handles conn defaults (high-res for slice display)
            state.structural=CONN_gui.refs.canonical.filename;
        end
    end
    state.strvol=conn_fileutils('spm_vol',state.structural);
    state.structuraldata=conn_fileutils('spm_read_vols',state.strvol);
    %state.structuraldata(isnan(state.structuraldata))=0; %%
    [hpatch,state.Prange]=conn_mesh_display_axref(state.strvol,state.structuraldata,state.pointer_mm, state.structural_method,state.Prange);
    for nview=1:3
        state.handles.patchref(nview)=patch;
        if isfield(state,'showaxref')&&numel(state.showaxref)>=nview&&state.showaxref(nview), onoff='on'; else onoff='off'; end
        set(state.handles.patchref(nview),'faces',hpatch(nview).faces,'vertices',hpatch(nview).vertices,'facevertexcdata',hpatch(nview).facevertexcdata,'facevertexalpha',state.facealpharef*hpatch(nview).facevertexalpha,'userdata',hpatch(nview).facevertexalpha,'facecolor','flat','edgecolor','none','facealpha','flat','alphadatamapping','none','FaceLighting', 'gouraud','visible',onoff);
    end
    if ~isempty(state.ref_color), set(state.handles.patchref,'facecolor',state.ref_color); end
    set(state.handles.patchref,'tag','reference_surface');
end
axis(state.handles.hax,'equal','tight');
lighting(state.handles.hax,'phong');
if isempty(state.material), set(state.handles.light,'visible','off');
else
    set(state.handles.light,'visible','on');
    material(state.material);
end
if state.handles.fullfigure, state.handles.hhide=uicontrol('style','frame','unit','norm','position',[0 0 1 1],'backgroundcolor',state.background,'foregroundcolor',state.background); end

if state.handles.fullfigure, hc=state.handles.hfig;
else hc=uicontextmenu(state.handles.hfig);
end
hc1=uimenu(hc,'Label','View');
uimenu(hc1,'Label','Left view','callback',{@conn_mesh_display_refresh,'view',[-1,0,0],[],-1},'tag','view');
uimenu(hc1,'Label','Right view','callback',{@conn_mesh_display_refresh,'view',[1,0,0],[],1},'tag','view');
uimenu(hc1,'Label','Left medial view','callback',{@conn_mesh_display_refresh,'view',[1,0,0],[1,0,.5],-1},'tag','view');
uimenu(hc1,'Label','Right medial view','callback',{@conn_mesh_display_refresh,'view',[-1,0,0],[-1,0,.5],1},'tag','view');
uimenu(hc1,'Label','Anterior view','callback',{@conn_mesh_display_refresh,'view',[0,1,0],[],0},'tag','view');
uimenu(hc1,'Label','Posterior view','callback',{@conn_mesh_display_refresh,'view',[0,-1,0],[],0},'tag','view');
uimenu(hc1,'Label','Superior view','callback',{@conn_mesh_display_refresh,'view',[0,-1e-8,1],[],0},'tag','view');
uimenu(hc1,'Label','Superior view (flip)','callback',{@conn_mesh_display_refresh,'view',[+1e-8,0,1],[],0},'tag','view','visible','off');
uimenu(hc1,'Label','Inferior view','callback',{@conn_mesh_display_refresh,'view',[0,-1e-8,-1],[],0},'tag','view');
uimenu(hc1,'Label','Left view (both hem)','callback',{@conn_mesh_display_refresh,'view',[-1,0,0],[],0},'tag','view');
uimenu(hc1,'Label','Right view (both hem)','callback',{@conn_mesh_display_refresh,'view',[1,0,0],[],0},'tag','view');
uimenu(hc1,'Label','Camera-view copy','callback',{@conn_mesh_display_refresh,'copy'},'separator','on');
uimenu(hc1,'Label','Camera-view paste','callback',{@conn_mesh_display_refresh,'paste'});
uimenu(hc1,'Label','Zoom in','separator','on','callback',{@conn_mesh_display_refresh,'zoomin'});
uimenu(hc1,'Label','Zoom out','callback',{@conn_mesh_display_refresh,'zoomout'});
uimenu(hc1,'Label','info','separator','on','callback',{@conn_mesh_display_refresh,'info'});
hc0=uimenu(hc,'Label','Surfaces');
hc1=uimenu(hc0,'Label','Brain surface');
thdl=[uimenu(hc1,'Label','White Matter','callback', {@conn_mesh_display_refresh,'brain',1},'tag','brain');
    uimenu(hc1,'Label','Grey Matter','callback',{@conn_mesh_display_refresh,'brain',2},'tag','brain');
    uimenu(hc1,'Label','Semi-inflated WM','callback',{@conn_mesh_display_refresh,'brain',3},'tag','brain');
    uimenu(hc1,'Label','Inflated WM','callback',{@conn_mesh_display_refresh,'brain',4},'tag','brain')];
set(thdl,'checked','off');set(thdl(state.selectedsurface),'checked','on');
thdl=[uimenu(hc1,'Label','High-resolution brain surface','separator','on','callback',{@conn_mesh_display_refresh,'brain',[],1},'tag','res');
    uimenu(hc1,'Label','Low-resolution brain surface','callback',{@conn_mesh_display_refresh,'brain',[],2},'tag','res')];
set(thdl,'checked','off');set(thdl(state.reducedpatch),'checked','on');
tvalues=[1 .9:-.1:.1 .05 0];
thdl=uimenu(hc1,'Label','Brain surface on','separator','on','callback',{@conn_mesh_display_refresh,'brain_transparency',1},'tag','brain_transparency');
hc2=uimenu(hc1,'Label','Brain surface transparent');
for n1=1:numel(tvalues)-1, thdl=[thdl,uimenu(hc2,'Label',num2str(n1-1),'callback',{@conn_mesh_display_refresh,'brain_transparency',tvalues(n1)},'tag','brain_transparency')]; end
thdl=[thdl,uimenu(hc1,'Label','Brain surface off','callback',{@conn_mesh_display_refresh,'brain_transparency',0},'tag','brain_transparency')];
[nill,idx]=min(abs(state.facealpha-[1 tvalues]));
set(thdl,'checked','off');set(thdl(max(1,min(numel(thdl),idx))),'checked','on');
uimenu(hc1,'Label','Brain surface color','callback',{@conn_mesh_display_refresh,'brain_color'});
uimenu(hc1,'Label','Brain surface activation','callback',{@conn_mesh_display_refresh,'repaint'});
uimenu(hc1,'Label','show medial mask','callback',{@conn_mesh_display_refresh,'mask','on'},'checked','on','tag','brain_mask');
uimenu(hc1,'Label','hide medial mask','callback',{@conn_mesh_display_refresh,'mask','off'},'tag','brain_mask');
if ~isempty(state.pVOL1)||~isempty(state.pVOL2)
    hc1=uimenu(hc0,'Label','Activation surface');
    tvalues=[1 .9:-.1:.1 .05 0];
    thdl=uimenu(hc1,'Label','Activation surface on','callback',{@conn_mesh_display_refresh,'act_transparency',1},'tag','act_transparency');
    hc2=uimenu(hc1,'Label','Activation surface transparent');
    for n1=1:numel(tvalues)-1, thdl=[thdl,uimenu(hc2,'Label',num2str(n1-1),'callback',{@conn_mesh_display_refresh,'act_transparency',tvalues(n1)},'tag','act_transparency')]; end
    thdl=[thdl,uimenu(hc1,'Label','Activation surface off','callback',{@conn_mesh_display_refresh,'act_transparency',0},'tag','act_transparency')];
    [nill,idx]=min(abs(state.facealphablob-[1 tvalues]));
    set(thdl,'checked','off');set(thdl(max(1,min(numel(thdl),idx))),'checked','on');
    uimenu(hc1,'Label','Positive activation surface only','callback',{@conn_mesh_display_refresh,'act_pos'});
    uimenu(hc1,'Label','Negative activation surface only','callback',{@conn_mesh_display_refresh,'act_neg'});
    uimenu(hc1,'Label','Positive & Negative activation surfaces','callback',{@conn_mesh_display_refresh,'act_posneg'});
    uimenu(hc1,'Label','Activation surface color','callback',{@conn_mesh_display_refresh,'act_color'});
end

% thdl=uimenu(hc1,'Label','Brain surface on','separator','on','callback',{@conn_mesh_display_refresh,'brain_transparency',0},'tag','brain_transparency');
% hc2=uimenu(hc1,'Label','Brain surface transparency');
% thdl=[];
% for n1=[0 .05 .1:.1:.9],thdl=[thdl,uimenu(hc2,'Label',num2str(1-n1),'callback',{@conn_mesh_display_refresh,'brain_transparency',n1},'tag','brain_transparency')]; end
% thdl=[thdl,uimenu(hc1,'Label','Brain surface opaque','callback',{@conn_mesh_display_refresh,'brain_transparency',1},'tag','brain_transparency')];
% set(thdl,'checked','off');set(thdl(max(1,min(numel(thdl),1+round(state.facealpha*10)))),'checked','on');
% %     thdl=[uimenu(hc1,'Label','Brain surface subcortical mask on','callback',{@conn_mesh_display_refresh,'mask','on'});
% %           uimenu(hc1,'Label','Brain surface subcortical mask off','callback',{@conn_mesh_display_refresh,'mask','off'})];
% %     set(thdl,'checked','off');set(thdl(2-(domask>0)),'checked','on');

hc1=uimenu(hc0,'Label','Subcortical surface');
tvalues=[1 .9:-.1:.1 .05 0];
thdl=uimenu(hc1,'Label','Subcortical surface on','callback',{@conn_mesh_display_refresh,'sub_transparency',1},'tag','sub_transparency');
hc2=uimenu(hc1,'Label','Subcortical surface transparent');
for n1=1:numel(tvalues)-1, thdl=[thdl,uimenu(hc2,'Label',num2str(n1-1),'callback',{@conn_mesh_display_refresh,'sub_transparency',tvalues(n1)},'tag','sub_transparency')]; end
thdl=[thdl,uimenu(hc1,'Label','Subcortical surface off','callback',{@conn_mesh_display_refresh,'sub_transparency',0},'tag','sub_transparency')];
[nill,idx]=min(abs(state.facealphasub-[1 tvalues]));
set(thdl,'checked','off');set(thdl(max(1,min(numel(thdl),idx))),'checked','on');
uimenu(hc1,'Label','Subcortical surface color','callback',{@conn_mesh_display_refresh,'sub_color'});
% thdl=[uimenu(hc1,'Label','Subcortical reference on','separator','on','callback',{@conn_mesh_display_refresh,'sub','on'},'tag','sub');
%     uimenu(hc1,'Label','Subcortical reference off','callback',{@conn_mesh_display_refresh,'sub','off'},'tag','sub')];
% set(thdl,'checked','off');set(thdl(2-(dosub>0)),'checked','on');

hc1=uimenu(hc0,'Label','Brainmask surface');
tvalues=[1 .9:-.1:.1 .05 0];
thdl=uimenu(hc1,'Label','Brainmask surface on','callback',{@conn_mesh_display_refresh,'mask_transparency',1},'tag','mask_transparency');
hc2=uimenu(hc1,'Label','Brainmask surface transparent');
for n1=1:numel(tvalues)-1, thdl=[thdl,uimenu(hc2,'Label',num2str(n1-1),'callback',{@conn_mesh_display_refresh,'mask_transparency',tvalues(n1)},'tag','mask_transparency')]; end
thdl=[thdl,uimenu(hc1,'Label','Brainmask surface off','callback',{@conn_mesh_display_refresh,'mask_transparency',0},'tag','mask_transparency')];
[nill,idx]=min(abs(state.facealphamask-[1 tvalues]));
set(thdl,'checked','off');set(thdl(max(1,min(numel(thdl),idx))),'checked','on');
uimenu(hc1,'Label','Brainmask surface color','callback',{@conn_mesh_display_refresh,'mask_color'});

hc1=uimenu(hc0,'Label','Other surfaces');
tvalues=[1 .9:-.1:.1 .05 0];
uimenu(hc1,'Label','Other surfaces add','callback',{@conn_mesh_display_refresh,'ud_select'});
uimenu(hc1,'Label','Other surfaces delete','callback',{@conn_mesh_display_refresh,'ud_delete'});
thdl=uimenu(hc1,'Label','Other surfaces on','callback',{@conn_mesh_display_refresh,'ud_transparency',1},'tag','ud_transparency','separator','on');
hc2=uimenu(hc1,'Label','Other surfaces transparent');
for n1=1:numel(tvalues)-1, thdl=[thdl,uimenu(hc2,'Label',num2str(n1-1),'callback',{@conn_mesh_display_refresh,'ud_transparency',tvalues(n1)},'tag','ud_transparency')]; end
thdl=[thdl,uimenu(hc1,'Label','Other surfaces off','callback',{@conn_mesh_display_refresh,'ud_transparency',0},'tag','ud_transparency')];
[nill,idx]=min(abs(state.facealphaud-[1 tvalues]));
set(thdl,'checked','off');set(thdl(max(1,min(numel(thdl),idx))),'checked','on');
uimenu(hc1,'Label','Other surfaces color','callback',{@conn_mesh_display_refresh,'ud_color'});

hc1=uimenu(hc0,'Label','Reference slices');
thdl=[uimenu(hc1,'Label','Sagittal reference on','callback',{@conn_mesh_display_refresh,'sag_ref','on'},'tag','sag_ref');
    uimenu(hc1,'Label','Sagittal reference off','callback',{@conn_mesh_display_refresh,'sag_ref','off'},'tag','sag_ref')];
set(thdl,'checked','off');set(thdl(2-(state.showaxref(1)>0)),'checked','on');
thdl=[uimenu(hc1,'Label','Coronal reference on','callback',{@conn_mesh_display_refresh,'cor_ref','on'},'tag','cor_ref');
    uimenu(hc1,'Label','Coronal reference off','callback',{@conn_mesh_display_refresh,'cor_ref','off'},'tag','cor_ref')];
set(thdl,'checked','off');set(thdl(2-(state.showaxref(2)>0)),'checked','on');
thdl=[uimenu(hc1,'Label','Axial reference on','callback',{@conn_mesh_display_refresh,'ax_ref','on'},'tag','ax_ref');
    uimenu(hc1,'Label','Axial reference off','callback',{@conn_mesh_display_refresh,'ax_ref','off'},'tag','ax_ref')];
set(thdl,'checked','off');set(thdl(2-(state.showaxref(3)>0)),'checked','on');
tvalues=[1 .9:-.1:.1 .05 0];
hc2=uimenu(hc1,'Label','Reference slices transparent');
thdl=[];for n1=1:numel(tvalues)-1, thdl=[thdl,uimenu(hc2,'Label',num2str(n1-1),'callback',{@conn_mesh_display_refresh,'ref_transparency',tvalues(n1)},'tag','ref_transparency')]; end
[nill,idx]=min(abs(state.facealpharef-tvalues));
set(thdl,'checked','off');set(thdl(max(1,min(numel(thdl),idx))),'checked','on');
uimenu(hc1,'Label','Reference slices position','callback',{@conn_mesh_display_refresh,'pos_ref'});

if ~isempty(state.handles.sphplots)&&any(any(state.handles.sphplots))||~isempty(state.handles.sphplots_txt)&&any(any(state.handles.sphplots_txt))
    hc0=uimenu(hc,'Label','ROIs');
    tvalues=[1 .9:-.1:.1 .05 0];
    thdl=uimenu(hc0,'Label','ROIs on','callback',{@conn_mesh_display_refresh,'roi_transparency',1},'tag','roi_transparency');
    hc2=uimenu(hc0,'Label','ROIs transparent');
    for n1=1:numel(tvalues)-1, thdl=[thdl,uimenu(hc2,'Label',num2str(n1-1),'callback',{@conn_mesh_display_refresh,'roi_transparency',tvalues(n1)},'tag','roi_transparency')]; end
    thdl=[thdl,uimenu(hc0,'Label','ROIs off','callback',{@conn_mesh_display_refresh,'roi_transparency',0},'tag','roi_transparency')];
    [nill,idx]=min(abs(state.facealpharoi-[1 tvalues]));
    set(thdl,'checked','off');set(thdl(max(1,min(numel(thdl),idx))),'checked','on');
    uimenu(hc0,'Label','ROIs color','callback',{@conn_mesh_display_refresh,'roi_color'});
end
if ~isempty(state.handles.sphplots)&&any(any(ishandle(state.handles.sphplots)))
    hc1=uimenu(hc0,'Label','Shapes','separator','on');
    if isfield(state.sphplots,'sph_shapes'), uimenu(hc1,'Label','Real ROI shape','callback',{@conn_mesh_display_refresh,'roi_shape','real'},'tag','roi_shape'); end
    uimenu(hc1,'Label','Sphere @ ROI center','callback',{@conn_mesh_display_refresh,'roi_shape','sphere'},'tag','roi_shape');
    uimenu(hc1,'Label','Cube @ ROI center','callback',{@conn_mesh_display_refresh,'roi_shape','cube'},'tag','roi_shape');
    uimenu(hc1,'Label','Star @ ROI center','callback',{@conn_mesh_display_refresh,'roi_shape','star'},'tag','roi_shape');
    uimenu(hc1,'Label','Increase size','callback',{@conn_mesh_display_refresh,'spheres',1.25});
    uimenu(hc1,'Label','Decrease size','callback',{@conn_mesh_display_refresh,'spheres',1/1.25});
end
if ~isempty(state.handles.sphplots_txt)&&any(any(state.handles.sphplots_txt))
    hc1=uimenu(hc0,'Label','Labels');
    uimenu(hc1,'Label','Labels on','callback',{@conn_mesh_display_refresh,'labels',1},'tag','labels');
    uimenu(hc1,'Label','Labels off','callback',{@conn_mesh_display_refresh,'labels',0},'tag','labels','checked','on');
    uimenu(hc1,'Label','Increase size','callback',{@conn_mesh_display_refresh,'labelsfont',1.25});
    uimenu(hc1,'Label','Decrease size','callback',{@conn_mesh_display_refresh,'labelsfont',1/1.25});
    uimenu(hc1,'Label','Bring closer','callback',{@conn_mesh_display_refresh,'labelscloser',2});
    uimenu(hc1,'Label','Reset position','callback',{@conn_mesh_display_refresh,'labelscloser',[]});
    uimenu(hc1,'Label','Edit labels','callback',{@conn_mesh_display_refresh,'labelsedit'});
end
if ~isempty(state.handles.connplots)&&any(any(state.handles.connplots~=0))
    hc0=uimenu(hc,'Label','Connections');
    tvalues=[1 .9:-.1:.1 .05 0];
    thdl=uimenu(hc0,'Label','Connections on','callback',{@conn_mesh_display_refresh,'con_transparency',1},'tag','con_transparency');
    hc2=uimenu(hc0,'Label','Connections transparent');
    for n1=1:numel(tvalues)-1, thdl=[thdl,uimenu(hc2,'Label',num2str(n1-1),'callback',{@conn_mesh_display_refresh,'con_transparency',tvalues(n1)},'tag','con_transparency')]; end
    thdl=[thdl,uimenu(hc0,'Label','Connections off','callback',{@conn_mesh_display_refresh,'con_transparency',0},'tag','con_transparency')];
    [nill,idx]=min(abs(state.facealphacon-[1 tvalues]));
    set(thdl,'checked','off');set(thdl(max(1,min(numel(thdl),idx))),'checked','on');
    uimenu(hc0,'Label','Connections color','callback',{@conn_mesh_display_refresh,'con_color'});
    hc1=uimenu(hc0,'Label','Bundling','separator','on');
    thdl=[];
    for n1=0:max([10,state.bundling]),thdl=[thdl,uimenu(hc1,'Label',num2str(n1),'callback',{@conn_mesh_display_refresh,'bundling',n1},'tag','bundling')]; end
    thdl=[thdl,uimenu(hc1,'Label','more','callback',{@conn_mesh_display_refresh,'bundling'},'tag','bundling')]; 
    set(thdl,'checked','off');set(thdl(max(1,min(numel(thdl),1+round(max([0,state.bundling]))))),'checked','on');
    hc1=uimenu(hc0,'Label','Thickness');
    uimenu(hc1,'Label','Increase thickness','callback',{@conn_mesh_display_refresh,'lines',1.5});
    uimenu(hc1,'Label','Decrease thickness','callback',{@conn_mesh_display_refresh,'lines',1/1.5});
end
hc1=uimenu(hc,'Label','Effects');
uimenu(hc1,'Label','normal','callback',{@conn_mesh_display_refresh,'material','dull'},'tag','material','checked','on');
uimenu(hc1,'Label','emphasis','callback',{@conn_mesh_display_refresh,'material',[.1 .75 .5 1 .5]},'tag','material');
uimenu(hc1,'Label','sketch','callback',{@conn_mesh_display_refresh,'material',[.1 1 1 .25 0]},'tag','material');
uimenu(hc1,'Label','shiny','callback',{@conn_mesh_display_refresh,'material',[.3 .6 .9 20 1]},'tag','material');
uimenu(hc1,'Label','metal','callback',{@conn_mesh_display_refresh,'material',[.3 .3 1 25 .5]},'tag','material');
uimenu(hc1,'Label','flat','callback',{@conn_mesh_display_refresh,'material',[]},'tag','material');
uimenu(hc1,'Label','bright','callback',{@conn_mesh_display_refresh,'light',1},'separator','on','tag','light');
uimenu(hc1,'Label','medium','callback',{@conn_mesh_display_refresh,'light',.5},'tag','light','checked','on');
uimenu(hc1,'Label','dark','callback',{@conn_mesh_display_refresh,'light',.2},'tag','light');
uimenu(hc1,'Label','white background','callback',{@conn_mesh_display_refresh,'background',[1 1 1]},'separator','on','tag','background');
uimenu(hc1,'Label','light background','callback',{@conn_mesh_display_refresh,'background',[.95 .95 .9]},'tag','background');
uimenu(hc1,'Label','dark background','callback',{@conn_mesh_display_refresh,'background',[.11 .11 .11]},'tag','background');
uimenu(hc1,'Label','black background','callback',{@conn_mesh_display_refresh,'background',[0 0 0]},'tag','background');
uimenu(hc1,'Label','color background','callback',{@conn_mesh_display_refresh,'background',[]},'tag','background','checked','on');
if ~isempty(state.handles.colorbar)
    uimenu(hc1,'Label','show colorbar','callback',{@conn_mesh_display_refresh,'colorbar','on'},'separator','on','tag','colorbar');
    uimenu(hc1,'Label','hide colorbar','callback',{@conn_mesh_display_refresh,'colorbar','off'},'tag','colorbar','checked','on');
    uimenu(hc1,'Label','colorbar limits','callback',{@conn_mesh_display_refresh,'colorbar','rescale'});
    hc2=uimenu(hc1,'Label','activation colormap');
    for n1={'normal','red','jet','hot','gray','bone','cool','hsv','spring','summer','autumn','winter','bluewhitered','random','brighter','darker','manual','color'}
        uimenu(hc2,'Label',n1{1},'callback',{@conn_mesh_display_refresh,'colormap',n1{1}},'tag','overlaycolormap');
    end
    hc2=uimenu(hc1,'Label','smoother display','checked','on','callback',{@conn_mesh_display_refresh,'black_transparency','on'},'tag','black_transparency');
    hc2=uimenu(hc1,'Label','raw data display','callback',{@conn_mesh_display_refresh,'black_transparency','off'},'tag','black_transparency');
end
uimenu(hc1,'Label','show axis coordinates','callback',{@conn_mesh_display_refresh,'axis','on'},'separator','on','tag','axis');
uimenu(hc1,'Label','hide axis coordinates','callback',{@conn_mesh_display_refresh,'axis','off'},'tag','axis','checked','on');
uimenu(hc1,'Label','show menubar','callback',{@conn_mesh_display_refresh,'menubar','on'},'separator','on','tag','menubar');
uimenu(hc1,'Label','hide menubar','callback',{@conn_mesh_display_refresh,'menubar','off'},'tag','menubar','checked','on');

hc1=uimenu(hc,'Label','Print');
uimenu(hc1,'Label','current view','callback',{@conn_mesh_display_refresh,'print',1});
uimenu(hc1,'Label','2-view row','callback',{@conn_mesh_display_refresh,'print',2});
uimenu(hc1,'Label','3-view mosaic','callback',{@conn_mesh_display_refresh,'print',3});
uimenu(hc1,'Label','4-view mosaic','callback',{@conn_mesh_display_refresh,'print',4});
uimenu(hc1,'Label','4-view column','callback',{@conn_mesh_display_refresh,'print',5});
uimenu(hc1,'Label','4-view row','callback',{@conn_mesh_display_refresh,'print',6});
uimenu(hc1,'Label','8-view mosaic','callback',{@conn_mesh_display_refresh,'print',7});

if ~isempty(state.dobookmarks)
    hc1=uimenu(hc,'Label','Bookmark');
    hc2=uimenu(hc1,'Label','Save','callback',{@conn_mesh_display_refresh,'bookmark'});
    if ~isempty(state.bookmark_filename),
        hc2=uimenu(hc1,'Label','Save as copy','callback',{@conn_mesh_display_refresh,'bookmarkcopy'});
    end
end

state.handles.fh=fh;
set(state.handles.hfig,'userdata',state);%'uicontextmenu',hc,
try
    state.handles.rotate3d=rotate3d(state.handles.hfig);
    set(state.handles.rotate3d,'ActionPostCallback',{@conn_mesh_display_refresh,'position'});
    set(state.handles.rotate3d,'enable','on');
catch
    state.handles.rotate3d=[];
end
conn_mesh_display_refresh([],[],'view',[-1,0,0],[],-1);
conn_mesh_display_refresh([],[],'zoomout');
%drawnow;
cameraviewfields={'cameraposition','cameratarget','cameraviewangle','cameraupvector'};
for n=1:numel(cameraviewfields)
    set(state.handles.hax,[cameraviewfields{n},'mode'],'manual');
end
set(state.handles.hfig,'name',figname);
conn_mesh_display_refresh([],[],'view',[-1,0,0],[],0);
if state.handles.fullfigure, set(state.handles.hfig,'visible','on'); end %drawnow; end
set(state.handles.hhide,'visible','off');
if ishandle(hmsg), delete(hmsg); end
%rotate3d on;

    function out=conn_mesh_display_refresh(hObject,eventdata,option,varargin)
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
        redrawnowcolorbar=false;
        switch(option)
            case 'none', return;
            case 'close', close(state.handles.hfig); return;
            case 'figurehandle',out=state.handles.hfig; return;
            case 'getstate'
                out=state;
                out=rmfield(out,{'handles','structuraldata'});
            case {'bookmark','bookmarkcopy'}
                tfilename=[];
                if numel(varargin)>0&&~isempty(varargin{1}), tfilename=varargin{1}; 
                elseif ~isempty(state.bookmark_filename)&&strcmp(option,'bookmark'), tfilename=state.bookmark_filename;
                end
                if numel(varargin)>1&&~isempty(varargin{2}), descr=cellstr(varargin{2}); 
                else descr=state.bookmark_descr;
                end
                fcn=regexprep(mfilename,'^conn_','');
                conn_args={fcn,conn_mesh_display_refresh([],[],'getstate')};
                [fullfilename,tfilename,descr]=conn_bookmark('save',...
                    tfilename,...
                    descr,...
                    conn_args);
                if isempty(fullfilename), return; end
                if ~SKIPBI, 
                    ht=conn_msgbox('Printing bookmark icon. Please wait...','',-1);
                    conn_print(state.handles.hfig,conn_prepend('',fullfilename,'.jpg'),'-nogui','-r50','-nopersistent'); 
                    if ishandle(ht), delete(ht); end
                end
                state.bookmark_filename=tfilename;
                state.bookmark_descr=descr;
                   %conn_args={fcn,conn_mesh_display_refresh([],[],'getstate')}; % re-save to include bookmark info?
                   %save(conn_prepend('',fullfilename,'.mat'),'conn_args');
                if 0, conn_msgbox(sprintf('Bookmark %s saved',fullfilename),'',2);
                else out=fullfilename;
                end
                return;
            case 'info',
                conn_msgbox([{'Volume:'},cellstr(state.info.vol),{' ','Activation surface:'},cellstr(state.info.surf)],'3d display info');
                return;

            case 'visible',
                if numel(varargin)<1||isempty(varargin{1}), v=1;
                else v=varargin{1};
                end
                if ischar(v), v=strcmp(lower(v),'on'); end
                if v, 
                    set(state.handles.hfig,'name','conn 3d display');
                    set(state.handles.hhide,'visible','off'); 
                else
                    set(state.handles.hfig,'name','rendering, please wait...');
                    set(state.handles.hhide,'visible','on'); 
                end
                
            case {'mask','remap','remap&draw','colormap','black_transparency','smoother'}
                if ~strcmp(option,'remap')&&state.handles.fullfigure, set(state.handles.hfig,'name','rendering, please wait...');drawnow; end
                if strcmp(option,'colormap')
                    cmap=varargin{1};
                    if ischar(cmap)
                        switch(cmap)
                            case 'normal', cmap=state.colormap_default; 
                            case 'red', cmap=[linspace(0,1,96)',zeros(96,2)];
                            case 'hot', cmap=hot(96);
                            case 'jet', cmap=fixedge(jet(96));
                            case 'gray', cmap=gray(96);
                            case 'bone', cmap=bone(96);
                            case 'cool',cmap=fixedge(cool(96));
                            case 'hsv',cmap=fixedge(hsv(256));
                            %case 'jet', cmap=jet(2*96);
                            %case 'gray', cmap=gray(2*96);
                            %case 'cool',cmap=cool(96);
                            %case 'hsv',cmap=hsv(2*96);
                            case 'bluewhitered', cmap=[zeros(1,48) linspace(0,1,48) ones(1,48) linspace(1,.5,48); linspace(0,1,96) linspace(1,0,48) zeros(1,48); linspace(.5,1,48) ones(1,48) linspace(1,0,48) zeros(1,48)]'; cmap=repmat(abs(linspace(-1,1,192)'),1,3).*cmap+(1-repmat(abs(linspace(-1,1,192)'),1,3))*1;
                            case 'parulawhite', cmap=parula(192); cmap=repmat(abs(linspace(-1,1,192)'),1,3).*cmap+(1-repmat(abs(linspace(-1,1,192)'),1,3))*1;
                            case 'spring',cmap=repmat(1-linspace(1,0,96)'.^8,[1,3]).*spring(96);
                            case 'summer',cmap=repmat(1-linspace(1,0,96)'.^8,[1,3]).*summer(96);
                            case 'autumn',cmap=repmat(1-linspace(1,0,96)'.^8,[1,3]).*autumn(96);
                            case 'winter',cmap=repmat(1-linspace(1,0,96)'.^8,[1,3]).*winter(96);
                            case 'random',cmap=rand(96,3);
                            case 'brighter',cmap=min(1,1/sqrt(.95)*get(state.handles.hfig,'colormap').^(1/2)); cmap=cmap(round(size(cmap,1)/2)+1:end,:);
                            case 'darker',cmap=.95*get(state.handles.hfig,'colormap').^2; cmap=cmap(round(size(cmap,1)/2)+1:end,:);
                            case 'manual',answer=conn_menu_inputdlg({'colormap (96x3)'},'',1,{mat2str(state.colormap(round(size(state.colormap,1)/2)+1:end,:))});if ~isempty(answer), answer=str2num(answer{1}); end;if ~any(size(answer,1)==[96,2*96]), return; end;cmap=max(0,min(1,answer));
                            case 'color',cmap=uisetcolor([],'Select color'); if isempty(cmap)||isequal(cmap,0), return; end;
                            otherwise, disp('unknown value');
                        end
                    end
                    if size(cmap,2)<3, cmap=cmap(:,min(size(cmap,2),1:3)); end
                    if size(cmap,1)==1, cmap=repmat(cmap,96,1); end
                    if size(cmap,1)~=2*96, cmap=[flipud(cmap(:,[2,3,1]));cmap]; end
                    state.colormap=cmap;
                    set(state.handles.hfig,'colormap',cmap);
                elseif strcmp(option,'black_transparency')||strcmp(option,'smoother')
                    str=varargin{1};
                    if strcmp(str,'on'), state.smoother=true;
                    else state.smoother=false;
                    end
                elseif strcmp(option,'mask'),
                    state.domask=strcmp(varargin{1},'on');
                end
                V=state.V0;
                show=~isnan(V)&V~=0;
                V(~show)=0;
                if ~isempty(state.Vrange)
                    if sign(state.Vrange(1))*sign(state.Vrange(end))==-1, state.dotwosided=true; state.Vrange=[min(state.Vrange) 0 max(state.Vrange)];
                    else state.dotwosided=false; state.Vrange=sort(state.Vrange([1 end]));
                    end
                else
                    if state.dotwosided
                        state.Vrange=[min(V(show)) 0 max(V(show))];
                    elseif any(V(show)),
                        if nnz(V(show)>0), state.Vrange=[0 max(V(show))];
                        else              state.Vrange=[min(V(show)) 0];
                        end
                    else state.Vrange=[];
                    end
                end
                if numel(state.Vrange)>2
                    m1=show&V>0;V(m1)=V(m1)/abs(state.Vrange(3));
                    m1=show&V<0;V(m1)=V(m1)/abs(state.Vrange(1));
                    V(show)=max(-1,min(1, V(show) ));
                elseif numel(state.Vrange)==2&&state.Vrange(2)>0, V(show)=max(0,min(1, (V(show)-state.Vrange(1))/(state.Vrange(2)-state.Vrange(1)) ));
                elseif numel(state.Vrange)==2&&state.Vrange(1)<0, V(show)=max(-1,min(0, (V(show)-state.Vrange(2))/(state.Vrange(2)-state.Vrange(1)) ));
                end
                
                alpha=1;
                cdat2=max(0,min(1, ind2rgb(round((size(state.colormap,1)+1)/2+emph*(size(state.colormap,1)-1)/2*V),state.colormap)));
                if ~isempty(state.brain_color), cdat0=cellfun(@(x)conn_bsxfun(@times,1-0*x,shiftdim(state.brain_color,-1)),data.curv,'uni',0); % cellfun(@(x)conn_bsxfun(@times,1-.05*x,shiftdim([.7,.65,.6],-1)),data.curv,'uni',0);
                else cdat0=state.cdat0; 
                end
                cdat=cellfun(@(x)conn_bsxfun(@times,1-alpha*(show),x) + conn_bsxfun(@times,alpha*(show),cdat2),cdat0,'uni',0);
                if state.domask, cdat=cellfun(@(x)conn_bsxfun(@times,data.mask,x)+conn_bsxfun(@times,~data.mask,shiftdim(.6*[1 1 1],-1)),cdat,'uni',0); end
                if state.smoother
                    smoothing=1;
                    for n1=1:2
                        A=spm_mesh_adjacency(data.rend{1}(n1));
                        A=double(speye(size(A,1))|A);
                        Ax=sparse(1:size(A,1),1:size(A,1),1./max(eps,sum(A,2)))*A;
                        mask=sqrt(max(0,double(show(:,n1))));
                        for n=1:smoothing,mask=Ax*mask;end; mask=mask/max(eps,max(mask(:)));
                        mask=2*min(mask,1-mask);
                        Ax=sparse(1:size(A,1),1:size(A,1),1-mask)+sparse(1:size(A,1),1:size(A,1),mask./max(eps,sum(A,2)))*A;
                        for n=1:5*smoothing,for n2=1:numel(cdat), cdat{n2}(:,n1,:)=Ax*permute(cdat{n2}(:,n1,:),[1 3 2]);end;end
                    end
                end
                state.cdat=cdat;
                if ~strcmp(option,'remap')
                    for n2=1:2, 
                        set(state.handles.patch(n2),'vertices',data.rend{state.selectedsurface}(n2).vertices(state.selected_vertices{state.reducedpatch},:),...
                            'faces',state.selected_faces{state.reducedpatch},...
                            'facevertexcdata',permute(state.cdat{state.selectedsurface}(state.selected_vertices{state.reducedpatch},n2,:),[1 3 2])); 
                    end
                    set(state.handles.patch,'facecolor','interp');
                    redrawnowcolorbar=true;
%                     if state.dotwosided,
%                         set(state.handles.colorbar(1),'ytick',[.5,64.5,128.5],'yticklabel',arrayfun(@(x)num2str(x,'%.2f'),state.Vrange,'uni',0),'ycolor',max(0,min(1, 1-state.background)));
%                         set(state.handles.colorbar(2),'cdata',max(0,min(1, ind2rgb(round((size(state.colormap,1)+1)/2+emph*(size(state.colormap,1)-1)/2*linspace(-1,1,128)'),state.colormap))));
%                     elseif any(state.Vrange>0)
%                         set(state.handles.colorbar(1),'ytick',[.5,128.5],'yticklabel',arrayfun(@(x)num2str(x,'%.2f'),state.Vrange,'uni',0),'ycolor',max(0,min(1, 1-state.background)));
%                         set(state.handles.colorbar(2),'cdata',max(0,min(1, ind2rgb(round((size(state.colormap,1)+1)/2+emph*(size(state.colormap,1)-1)/2*linspace(0,1,128)'),state.colormap))));
%                     elseif any(state.Vrange<0)
%                         set(state.handles.colorbar(1),'ytick',[.5,128.5],'yticklabel',arrayfun(@(x)num2str(x,'%.2f'),state.Vrange,'uni',0),'ycolor',max(0,min(1, 1-state.background)));
%                         set(state.handles.colorbar(2),'cdata',max(0,min(1, ind2rgb(round((size(state.colormap,1)+1)/2+emph*(size(state.colormap,1)-1)/2*linspace(-1,0,128)'),state.colormap))));
%                     end
                end
                if ~strcmp(option,'remap')&&state.handles.fullfigure, set(state.handles.hfig,'name','conn 3d display'); end
                
            case 'view'
                if numel(varargin)<1||isempty(varargin{1}), v=get(state.handles.hax,'cameraposition');
                else v=varargin{1};
                end
                if numel(varargin)<2||isempty(varargin{2}), vl=v; else vl=varargin{2}; end
                if numel(varargin)<3||isempty(varargin{3}), side=0; else side=varargin{3}; end
                view(state.handles.hax,v);
                %up=null(v); up=sum(up,2)'; state.up=1.5*up/norm(up);
                state.up=v(:)'/norm(v);
                try, vup=get(state.handles.hax,'cameraupvector');
                catch, vup=state.up;
                end
                %up=get(state.handles.hax,'cameraup'); up=up(:)';
                if numel(state.handles.light)>2, set(state.handles.light(2:end-1),'position',vl); end
                set(state.handles.light(1),'position',1000*(vl(:)'/norm(vl)-.1*vup(:)'/norm(vup)));
                set(state.handles.light(end),'position',1000*(vl(:)'/norm(vl)+.1*vup(:)'/norm(vup)));
                set(state.handles.patch,'visible','off'); 
                set(state.handles.subpatch,'visible','off'); 
                set(state.handles.maskpatch,'visible','off');
                set(state.handles.sphplots,'visible','off');
                set(state.handles.sphplots_shapes(ishandle(state.handles.sphplots_shapes)),'visible','off');
                set(state.handles.sphplots_txt,'visible','off');
                set(state.handles.connplots(state.handles.connplots~=0),'visible','off');
                set([state.handles.patchblob1 state.handles.patchblob2 state.handles.patchblob3 state.handles.patchblob4],'visible','off');
                set(state.handles.patchblobother,'visible','off');
                if side<=0, %left
                    set(state.handles.patch(1),'visible','on'); 
                    set(state.handles.subpatch(1),'visible','on');
                    set(state.handles.maskpatch(1),'visible','on'); 
                    if ~isempty(state.sphplots)&&~isempty(state.sphplots.sph_xyz), 
                        thandles=state.handles.sphplots(state.sphplots.sph_xyz(:,1)<=5);
                        if isequal(state.roishape,'real'), 
                            thandles2=state.handles.sphplots_shapes(state.sphplots.sph_xyz(:,1)<=5); 
                            tok2=state.handles.sphplots_shapesok(state.sphplots.sph_xyz(:,1)<=5)&ishandle(thandles2); 
                            set(thandles2(tok2),'visible','on');
                            set(thandles,'visible','on'); %thandles(~tok2)
                        else
                            set(thandles(ishandle(thandles)),'visible','on');
                        end
                        if state.facealphatxt, set(state.handles.sphplots_txt(state.sphplots.sph_xyz(:,1)<=5),'visible','on'); end
                    end
                    if ~isempty(state.handles.connplots), temp=state.handles.connplots(state.sphplots.sph_xyz(:,1)<=5,state.sphplots.sph_xyz(:,1)<=5);set(temp(temp~=0),'visible','on'); end
                    set([state.handles.patchblob1 state.handles.patchblob3],'visible','on');
                    set(state.handles.patchblobother(state.patchblobother_x<=5),'visible','on');
                end
                if side>=0, %right
                    set(state.handles.patch(2),'visible','on'); 
                    set(state.handles.subpatch(2),'visible','on');
                    set(state.handles.maskpatch(2),'visible','on'); 
                    if ~isempty(state.sphplots)&&~isempty(state.sphplots.sph_xyz)
                        thandles=state.handles.sphplots(state.sphplots.sph_xyz(:,1)>=-5);
                        if isequal(state.roishape,'real'), 
                            thandles2=state.handles.sphplots_shapes(state.sphplots.sph_xyz(:,1)>=-5); 
                            tok2=state.handles.sphplots_shapesok(state.sphplots.sph_xyz(:,1)>=-5)&ishandle(thandles2); 
                            set(thandles2(tok2),'visible','on');
                            set(thandles,'visible','on'); %thandles(~tok2)
                        else
                            set(thandles(ishandle(thandles)),'visible','on');
                        end
                        if state.facealphatxt, set(state.handles.sphplots_txt(state.sphplots.sph_xyz(:,1)>=-5),'visible','on'); end
                    end
                    if ~isempty(state.handles.connplots), temp=state.handles.connplots(state.sphplots.sph_xyz(:,1)>=-5,state.sphplots.sph_xyz(:,1)>=-5); set(temp(temp~=0),'visible','on'); end
                    set([state.handles.patchblob2 state.handles.patchblob4],'visible','on');
                    set(state.handles.patchblobother(state.patchblobother_x>=-5),'visible','on');
                end
                if side==0, set(state.handles.connplots(state.handles.connplots~=0),'visible','on'); end
                for n=1:numel(state.handles.sphplots_txt), set(state.handles.sphplots_txt(n),'position',state.sphplots.sph_xyz(n,:)+state.up*state.fontclose*state.sphplots.sph_r(n)); end
                if ~isempty(state.handles.leftrightlabel)
                    if 1,%isequal(get(state.handles.hax,'visible'),'on')
                        if abs(v(1))/max(eps,norm(v))>.10, set(state.handles.leftrightlabel,'visible','off','color',.5*[1 1 1]);
                        else, set(state.handles.leftrightlabel,'visible','on','color',.5*[1 1 1]+(abs(v(1))/max(eps,norm(v))-.01)/(.10-.01)*(state.background-.5*[1 1 1]));
                        end
                    end
                end
            case 'copy'
                for n=1:numel(cameraviewfields)
                    tempv.(cameraviewfields{n})=get(state.handles.hax,cameraviewfields{n});
                end
                assignin('base','conn_mesh_display_view',tempv);
            case 'paste'
                try
                    tempv=evalin('base','conn_mesh_display_view');
                    for n=1:numel(cameraviewfields)
                        set(state.handles.hax,cameraviewfields{n},tempv.(cameraviewfields{n}));
                    end
                    conn_mesh_display_refresh([],[],'position');
                end
            case 'zoomin'
                x=get(gca,'cameraviewangle'); set(gca,'cameraviewangle',x*.75); 
                %x=get(gca,{'xlim','ylim','zlim'}); x=cellfun(@(x)1/1.25*x,x,'uni',0); set(gca,{'xlim','ylim','zlim'},x);
            case 'zoomout'
                x=get(gca,'cameraviewangle'); set(gca,'cameraviewangle',x/.75); 
                %x=get(gca,{'xlim','ylim','zlim'}); x=cellfun(@(x)1.25*x,x,'uni',0); set(gca,{'xlim','ylim','zlim'},x);
            case 'brain'
                N=1;
                if numel(varargin)>0&&~isempty(varargin{1}), N=1; newselected=varargin{1}; else newselected=state.selectedsurface; end
                if numel(varargin)>1&&~isempty(varargin{2}), state.reducedpatch=varargin{2}; end
                str=get(state.handles.hfig,'name'); 
                if N>1, 
                    set(state.handles.hfig,'name','press <spacebar> to stop'); 
                    set(state.handles.hfig,'currentcharacter','0');
                end
                for n1=1:N, 
                    if N>1&&isequal(get(state.handles.hfig,'CurrentCharacter'),' '), break; end; 
                    for n2=1:2, 
                        set(state.handles.patch(n2),'vertices',(N-n1)/N*data.rend{state.selectedsurface}(n2).vertices(state.selected_vertices{state.reducedpatch},:)+n1/N*data.rend{newselected}(n2).vertices(state.selected_vertices{state.reducedpatch},:),...
                            'faces',state.selected_faces{state.reducedpatch},...
                            'facevertexcdata',permute((N-n1)/N*state.cdat{state.selectedsurface}(state.selected_vertices{state.reducedpatch},n2,:)+n1/N*state.cdat{newselected}(state.selected_vertices{state.reducedpatch},n2,:),[1 3 2])); 
                    end
                    if N>1, drawnow; end
                end
                state.selectedsurface=newselected; 
                set(state.handles.hfig,'name',str);
            case {'spheres','sphereshape','roi_shape'}
                if strcmp(option,'spheres'),
                    scale=varargin{1};
                    state.sphplots.sph_r=state.sphplots.sph_r*scale;
                else
                    state.roishape=varargin{1};
                end
                if isequal(state.roishape,'real')
                    t0={}; troifiles=[];
                    for n1=1:numel(sphplots.sph_shapes)
                        if ischar(state.sphplots.sph_shapes{n1}),
                            t1=regexprep(state.sphplots.sph_shapes{n1},'\..*$','');
                            t2=regexprep(state.sphplots.sph_shapes{n1},'^[^\.]*\.','');
                            if isempty(troifiles), 
                                try, [troifiles,troinames]=conn_module('get','rois'); 
                                catch, fprintf('Unable to find ROIs in current CONN project\n'); return;
                                end
                            end
                            [tok,troiidx]=ismember(t1,troinames);
                            if tok, t0=[t0; {troifiles{troiidx}{1}{1}, t2, n1}];
                            else fprintf('Unable to find ROI %s in current CONN project\n',state.sphplots.sph_shapes{n1});
                            end
                        end
                    end
                    if ~isempty(t0), 
                        thmsg=conn_msgbox('Preparing ROI surfaces. Please wait...','',-1);
                        [ut1,nill,it1]=unique(t0(:,1));
                        for n1=1:numel(ut1)
                            try
                                p2=find(it1==n1);
                                p1=conn_surf_volume(ut1{n1},0,t0(it1==n1,2));
                                for n2=1:numel(p2), state.sphplots.sph_shapes{p2(n2)}=p1(n2); end
                            catch
                                fprintf('Unable to compute surface ROIs from %s\n',ut1{n1});
                            end
                        end
                        if ishandle(thmsg), delete(thmsg); end
                    end
                    for n1=1:numel(state.sphplots.sph_shapes)
                        if isstruct(state.sphplots.sph_shapes{n1})
                            set(state.handles.sphplots_shapes(n1),'vertices',state.sphplots.sph_shapes{n1}.vertices,'faces',state.sphplots.sph_shapes{n1}.faces);
                            state.handles.sphplots_shapesok(n1)=true;
                        end
                    end
                    for n1=1:numel(state.handles.sphplots)
                        if numel(state.handles.sphplots_shapes)>=n1&&ishandle(state.handles.sphplots_shapes(n1))&&state.handles.sphplots_shapesok(n1)
                            if isequal(get(state.handles.sphplots(n1),'visible'),'on')
                                set(state.handles.sphplots_shapes(n1),'visible','on');
                                %set(state.handles.sphplots(n1),'visible','off');
                            end
                        end
                    end
                    [x,y,z]=sphere(32);
                    f=1;
                    x=sign(x).*abs(x).^f;y=sign(y).*abs(y).^f;z=sign(z).*abs(z).^f;
                    for n=1:numel(state.handles.sphplots),
                        set(state.handles.sphplots(n),'xdata',state.sphplots.sph_xyz(n,1)+1*x,'ydata',state.sphplots.sph_xyz(n,2)+1*y,'zdata',state.sphplots.sph_xyz(n,3)+1*z,'facealpha',1,'facecolor',.5*[1 1 1]);
                    end
                    %set(state.handles.sphplots_shapes(ishandle(state.handles.sphplots_shapes)),'visible','on');
                    %set(state.handles.sphplots(ishandle(state.handles.sphplots)),'visible','off');
                    for n=1:numel(state.handles.sphplots_txt), set(state.handles.sphplots_txt(n),'position',state.sphplots.sph_xyz(n,:)+state.up*state.fontclose*state.sphplots.sph_r(n)); end
                else
                    [x,y,z]=sphere(32);
                    switch(state.roishape)
                        case 'sphere', f=1;
                        case 'cube', f=.25;
                        case 'star', f=4;
                    end
                    x=sign(x).*abs(x).^f;y=sign(y).*abs(y).^f;z=sign(z).*abs(z).^f;
                    for n=1:numel(state.handles.sphplots),
                        set(state.handles.sphplots(n),'xdata',state.sphplots.sph_xyz(n,1)+state.sphplots.sph_r(n)*x,'ydata',state.sphplots.sph_xyz(n,2)+state.sphplots.sph_r(n)*y,'zdata',state.sphplots.sph_xyz(n,3)+state.sphplots.sph_r(n)*z,'facealpha',state.facealpharoi,'facecolor',state.roi_color(1+rem(n-1,size(state.roi_color,1)),:));
                    end
                    for n=1:numel(state.handles.sphplots_txt), set(state.handles.sphplots_txt(n),'position',state.sphplots.sph_xyz(n,:)+state.up*state.fontclose*state.sphplots.sph_r(n)); end
                    for n1=1:numel(state.handles.sphplots)
                        if numel(state.handles.sphplots_shapes)>=n1&&ishandle(state.handles.sphplots_shapes(n1))&&state.handles.sphplots_shapesok(n1), 
                            if isequal(get(state.handles.sphplots_shapes(n1),'visible'),'on')
                                set(state.handles.sphplots(n1),'visible','on');
                                set(state.handles.sphplots_shapes(n1),'visible','off');
                            end
                        end
                    end
                    %set(state.handles.sphplots(ishandle(state.handles.sphplots)),'visible','on');
                    %set(state.handles.sphplots_shapes(ishandle(state.handles.sphplots_shapes)),'visible','off');
                end
            case 'labels'
                if ischar(varargin{1}), varargin{1}=strcmpi(varargin{1},'on'); end 
                if varargin{1}, str='on'; state.facealphatxt=1; else str='off'; state.facealphatxt=0; end
                set(state.handles.sphplots_txt,'visible',str);
            case {'labels_font','labelsfont'},
                state.fontsize=state.fontsize*varargin{1};
                set(state.handles.sphplots_txt,'fontsize',max(1,round(state.fontsize)));
            case {'labels_closer','labelscloser'},
                if isempty(varargin{1}), state.fontclose=1;
                else state.fontclose=state.fontclose*varargin{1};
                end
                for n=1:numel(state.handles.sphplots_txt), set(state.handles.sphplots_txt(n),'position',state.sphplots.sph_xyz(n,:)+state.up*state.fontclose*state.sphplots.sph_r(n)); end
            case {'labels_edit','labelsedit'},
                if numel(varargin)>0&&~isempty(varargin{1}), 
                    newname=varargin{1}; 
                    for n=1:numel(newname), set(state.handles.sphplots_txt(n),'string',newname{n}); end
                    state.sphplots.sph_names=newname;
                else
                    name=state.sphplots.sph_names;
                    %name=get(state.handles.sphplots_txt,'string');
                    ok=true;
                    thfig=dialog('units','norm','position',[.3,.4,.6,.4],'windowstyle','normal','name','ROI labels','color','w','resize','on');
                    uicontrol(thfig,'style','text','units','norm','position',[.1,.85,.8,.10],'string',sprintf('New ROI label names (%d)',numel(name)),'backgroundcolor','w','fontsize',9+CONN_gui.font_offset,'fontweight','bold');
                    ht1=uicontrol(thfig,'style','edit','units','norm','position',[.1,.30,.8,.55],'max',2,'string',name,'fontsize',8+CONN_gui.font_offset,'horizontalalignment','left','tooltipstring','manually edit the ROI labels');
                    ht2=uicontrol(thfig,'style','edit','units','norm','position',[.1,.20,.8,.1],'max',1,'string','','fontsize',8+CONN_gui.font_offset,'horizontalalignment','left','tooltipstring','enter Matlab command for fast-editing all ROIs simultaneously (str is input variable cell array; ouput is cell array; e.g. "lower(str)")','callback','ht1=get(gcbo,''userdata''); set(ht1,''string'',feval(inline(get(gcbo,''string''),''str''),get(ht1,''string'')))','userdata',ht1);
                    uicontrol(thfig,'style','pushbutton','string','Apply','units','norm','position',[.1,.01,.38,.10],'callback','uiresume','fontsize',8+CONN_gui.font_offset);
                    uicontrol(thfig,'style','pushbutton','string','Cancel','units','norm','position',[.51,.01,.38,.10],'callback','delete(gcbf)','fontsize',8+CONN_gui.font_offset);
                    while ok
                        uiwait(thfig);
                        ok=ishandle(thfig);
                        if ok,
                            newname=get(ht1,'string');
                            if numel(newname)~=numel(name), conn_msgbox(sprintf('Number of labels entered (%d) does not match expected value (%d)',numel(newname),numel(name)),'',2);
                            else
                                delete(thfig);
                                for n=1:numel(newname), set(state.handles.sphplots_txt(n),'string',newname{n}); end
                                state.sphplots.sph_names=newname;
                                ok=false;
                            end
                        end
                    end
                end
            case {'con_width','lines'}
                scale=varargin{1};
                state.WidthLines=state.WidthLines*scale;
                for n=find(state.handles.connplots(:)'~=0),
                    nlines=state.connplotsLine(n);
                    tlines=conn_mesh_display_patch(state.Lines(:,1,nlines,1+max([0,state.bundling])),state.Lines(:,2,nlines,1+max([0,state.bundling])),state.Lines(:,3,nlines,1+max([0,state.bundling])),state.WidthLines(:,1,nlines,1+max([0,state.bundling])));
                    set(state.handles.connplots(n),'vertices',tlines.vertices);
                end
%                 for n=find(state.handles.connplots(:)'~=0), 
%                     xyz=get(state.handles.connplots(n),{'xdata','ydata','zdata'}); 
%                     xyz2=reshape(cat(3,xyz{:}),[],3); 
%                     B=get(state.handles.connplots(n),'userdata'); 
%                     xyz2=conn_bsxfun(@plus,mean(xyz2,1),conn_bsxfun(@minus,xyz2,mean(xyz2,1))*B*diag([scale scale 1])*B'); 
%                     set(state.handles.connplots(n),{'xdata','ydata','zdata'},cellfun(@(x)reshape(x,size(xyz{1})),num2cell(xyz2,1),'uni',0)); 
%                 end
            case {'con_bundling','bundling'}
                if numel(varargin)>=1, str=varargin{1};
                else
                    answer=conn_menu_inputdlg({'Bundling level?'},'',1,{num2str(max([0,state.bundling]))});
                    if isempty(answer), return; end
                    str=max(0,round(str2num(answer{1})));
                    if ~(numel(str)==1&&~isnan(str)&&~isinf(str)), return; end
                end
                state.bundling=str;
                if 1+max([0,state.bundling])>size(state.Lines,4),
                    [state.Lines,state.WidthLines] = conn_menu_bundle(state.Lines,state.WidthLines,state.ValueLines,max([0,state.bundling]),true);
                end
                for n=find(state.handles.connplots(:)'~=0),
                    nlines=state.connplotsLine(n);
                    tlines=conn_mesh_display_patch(state.Lines(:,1,nlines,1+max([0,state.bundling])),state.Lines(:,2,nlines,1+max([0,state.bundling])),state.Lines(:,3,nlines,1+max([0,state.bundling])),state.WidthLines(:,1,nlines,1+max([0,state.bundling])));
                    set(state.handles.connplots(n),'vertices',tlines.vertices);
                end
            case 'material'
                str=varargin{1};
                if isempty(str), set(state.handles.light,'visible','off');
                else
                    set(state.handles.light,'visible','on');
                    material(str);
                end
                state.material=str;
            case 'light'
                color=varargin{1};
                if numel(color)==1, color=color*[1 1 1]; end
                set(state.handles.light,'color',color);
                state.light_color=color;
            case 'axis',
                str=varargin{1};
                set(state.handles.hax,'visible',str);
                %set(state.handles.leftrightlabel,'visible',str);
                if strcmp(str,'on'), set(state.handles.hax,'projection','orth'); %perspective'); 
                else set(state.handles.hax,'projection','orth'); 
                end
            case 'menubar',
                str=varargin{1};
                if strcmp(str,'on'), set(state.handles.hfig,'menubar','figure'); 
                else set(state.handles.hfig,'menubar','none'); 
                end
            case 'background'
                if numel(varargin)>0&&~isempty(varargin{1}), color=varargin{1};
                else color=uisetcolor(state.background,'Select color'); if isempty(color)||isequal(color,0), return; end; 
                end
                set([state.handles.hfig state.handles.hax],'color',color);
                set(state.handles.sphplots_txt,'color',1-round(mean(color))*[1 1 1]);
                state.background=color;
                redrawnowcolorbar=true;
            case 'brain_transparency'
                scale=varargin{1};
                oldscale=state.facealpha;
                state.facealpha=max(eps,scale);
                set([state.handles.patch],'facealpha',state.facealpha); 
                if state.facealpha<1e-10&&oldscale>1e-10&&state.reducedpatch~=2, 
                    reducedpatch=state.reducedpatch;
                    conn_mesh_display_refresh([],[],'brain',[],2);
                    state.reducedpatch=reducedpatch;
                elseif state.facealpha>1e-10&&oldscale<1e-10&&state.reducedpatch~=2, 
                    conn_mesh_display_refresh([],[],'brain',[],state.reducedpatch);
                end
            case 'brain_color'
                if numel(varargin)>0, color=varargin{1};
                else color=uisetcolor([],'Select color'); 
                end
                if isequal(color,0), color=[]; end
                state.brain_color=color;
                if isempty(color)||~isequal(state.info.surf,'none'), conn_mesh_display_refresh([],[],'remap&draw'); 
                else set(state.handles.patch,'facecolor',color); 
                end

            case 'act_transparency'
                scale=varargin{1};
                state.facealphablob=max(eps,scale);
                set([state.handles.patchblob1 state.handles.patchblob2 state.handles.patchblob3 state.handles.patchblob4],'facealpha',state.facealphablob);
            case 'act_pos'
                set([state.handles.patchblob1 state.handles.patchblob2],'visible','on');
                set([state.handles.patchblob3 state.handles.patchblob4],'visible','off');
            case 'act_neg'
                set([state.handles.patchblob1 state.handles.patchblob2],'visible','off');
                set([state.handles.patchblob3 state.handles.patchblob4],'visible','on');
            case 'act_posneg'
                set([state.handles.patchblob1 state.handles.patchblob2],'visible','on');
                set([state.handles.patchblob3 state.handles.patchblob4],'visible','on');
            case 'act_color'
                if numel(varargin)>0, color=varargin{1};
                else color=uisetcolor([],'Select color'); if isempty(color)||isequal(color,0), return; end; 
                end
                state.act_color=color;
                set([state.handles.patchblob1 state.handles.patchblob2],'facecolor',state.act_color);
                set([state.handles.patchblob3 state.handles.patchblob4],'facecolor',state.act_color(:,[2 3 1]));
            case 'sub_transparency'
                scale=varargin{1};
                state.facealphasub=max(eps,scale);
                set([state.handles.subpatch],'facealpha',state.facealphasub); 
            case 'sub_color'
                if numel(varargin)>0, color=varargin{1};
                else color=uisetcolor([],'Select color'); if isempty(color)||isequal(color,0), return; end; 
                end
                set(state.handles.subpatch,'facecolor',color); 
                state.sub_color=color;
            case {'brainmask_transparency','mask_transparency'}
                scale=varargin{1};
                state.facealphamask=max(eps,scale);
                set([state.handles.maskpatch],'facealpha',state.facealphamask); 
            case {'brainmask_color','mask_color'}
                if numel(varargin)>0, color=varargin{1};
                else color=uisetcolor([],'Select color'); if isempty(color)||isequal(color,0), return; end; 
                end
                set(state.handles.maskpatch,'facecolor',color); 
                state.mask_color=color;
            case 'ref_transparency'
                scale=varargin{1};
                state.facealpharef=max(eps,scale);
                for tnview=1:3
                    c=get(state.handles.patchref(tnview),'userdata');
                    set(state.handles.patchref(tnview),'facevertexalpha',state.facealpharef*c);
                end
            case 'ref_color'
                if numel(varargin)>0, color=varargin{1};
                else color=uisetcolor([],'Select color'); if isempty(color)||isequal(color,0), return; end; 
                end
                set(state.handles.patchref,'facecolor',color); 
                state.ref_color=color;
            case 'ref_material',
                str=varargin{1};
                material(state.handles.patchref,str);
            case 'ud_transparency'
                scale=varargin{1};
                state.facealphaud=max(eps,scale);
                if state.facealphaud>0&&isempty(state.handles.patchblobother), conn_mesh_display_refresh([],[],'ud_select'); end
                set(state.handles.patchblobother,'facealpha',state.facealphaud); 
            case 'ud_color'
                if numel(varargin)>0, color=varargin{1};
                elseif numel(state.handles.patchblobother)>1
                    color=[];
                    if isfield(state,'ud_color'), color=state.ud_color; end
                    if isempty(color), color=[0 .45 .74]; end
                    answer=conn_menu_inputdlg({sprintf('Surface color (%dxRGB)  e.g. rand(%d,3)',numel(state.handles.patchblobother),numel(state.handles.patchblobother))},'',1,{mat2str(color)});
                    if isempty(answer), return; end
                    color=str2num(answer{1});
                    if size(color,2)~=3, return; end
                else color=uisetcolor([],'Select color'); if isempty(color)||isequal(color,0), return; end; 
                end
                if isequal(color,'rand'), color=rand(numel(state.handles.patchblobother),3); end
                if size(color,1)==1, set(state.handles.patchblobother,'facecolor',color); 
                else for n1=1:numel(state.handles.patchblobother), set(state.handles.patchblobother(n1),'facecolor',color(1+rem(n1-1,size(color,1)),:)); end
                end
                state.ud_color=color;
            case 'roi_transparency'
                scale=varargin{1};
                state.facealpharoi=max(eps,scale);
                set([state.handles.sphplots state.handles.sphplots_shapes(ishandle(state.handles.sphplots_shapes))],'facealpha',state.facealpharoi); 
            case 'roi_color'
                if numel(varargin)>0, color=varargin{1};
                elseif numel(state.handles.sphplots)>1
                    color=[];
                    if isfield(state,'roi_color'), color=state.roi_color; end
                    if isempty(color), color=[0 .45 .74]; end
                    answer=conn_menu_inputdlg({sprintf('ROI color (%dxRGB)  e.g. rand(%d,3)',numel(state.handles.sphplots),numel(state.handles.sphplots))},'',1,{mat2str(color)});
                    if isempty(answer), return; end
                    color=str2num(answer{1});
                    if size(color,2)~=3, return; end
                else color=uisetcolor([],'Select color'); if isempty(color)||isequal(color,0), return; end; 
                end
                if isequal(color,'rand'), color=rand(numel(state.handles.sphplots),3); end
                %[nill,idx]=sort(state.sphplots.sph_r);color(idx,:)=color;
                state.roi_color=color;
                if size(state.roi_color,1)==1, set([state.handles.sphplots state.handles.sphplots_shapes(ishandle(state.handles.sphplots_shapes))],'facecolor',state.roi_color);
                else
                    for n1=1:numel(state.handles.sphplots),
                        set(state.handles.sphplots(n1),'facecolor',state.roi_color(1+rem(n1-1,size(state.roi_color,1)),:));
                        if numel(state.handles.sphplots_shapes)>=n1&&ishandle(state.handles.sphplots_shapes(n1)), set(state.handles.sphplots_shapes(n1),'facecolor',state.roi_color(1+rem(n1-1,size(state.roi_color,1)),:)); end
                    end
                end
                %if size(color,1)==1, set(state.handles.sphplots,'facecolor',color); 
                %else for n1=1:numel(state.handles.sphplots), set(state.handles.sphplots(n1),'facecolor',color(1+rem(n1-1,size(color,1)),:)); end
                %end
%                 if numel(varargin)>0, color=varargin{1};
%                 else color=uisetcolor([],'Select color'); if isempty(color)||isequal(color,0), return; end; 
%                 end
%                 set(state.handles.sphplots,'facecolor',color); 
            case 'con_transparency'
                scale=varargin{1};
                state.facealphacon=max(eps,scale);
                set(state.handles.connplots(state.handles.connplots~=0),'facealpha',state.facealphacon); 
            case 'con_color'
                if numel(varargin)>0, color=varargin{1};
                elseif nnz(state.handles.connplots~=0)>1
                    color=[];
                    if isfield(state,'con_color'), color=state.con_color; end
                    if isempty(color), color=[0 .45 .74]; end
                    answer=conn_menu_inputdlg({sprintf('Connection color (%dxRGB)  e.g. rand(%d,3)',numel(state.handles.sphplots),numel(state.handles.sphplots))},'',1,{mat2str(color)});
                    if isempty(answer), return; end
                    color=str2num(answer{1});
                    if size(color,2)~=3, return; end
                else color=uisetcolor([],'Select color'); if isempty(color)||isequal(color,0), return; end; 
                end
                if isequal(color,'rand'), color=rand(nnz(state.handles.connplots~=0),3); end
                state.con_color=color;
                if size(state.con_color,1)==1, set(state.handles.connplots(state.handles.connplots~=0),'facecolor',state.con_color); 
                else
                    tidx=find(state.handles.connplots~=0);
                    for n1=1:numel(tidx), 
                        set(state.handles.connplots(tidx(n1)),'facecolor',state.con_color(1+rem(n1-1,size(state.con_color,1)),:));
                    end
                end
            case 'ud_delete'
                delete(state.handles.patchblobother(ishandle(state.handles.patchblobother)));
                state.handles.patchblobother=[];
                state.patchblobother_x=[];
            case 'ud_select'
                %otherVOL=spm_select(inf,'image','Select a file');
                %if isempty(otherVOL), return; end
                if numel(varargin)>0, 
                    totherVOL=varargin{1};
                    tpVOLother=conn_surf_volume(totherVOL,varargin{2:end});
                    if isempty(tpVOLother), return; end
                    facecolor=[.5 .5 0];
                else
                    [tfilename,tpathname]=conn_fileutils('uigetfile','*.nii; *.img','Select file',fullfile(fileparts(which(mfilename)),'rois'));
                    if ischar(tfilename), totherVOL=fullfile(tpathname,tfilename);
                    else return;
                    end
                    tpVOLother=conn_surf_volume(totherVOL);
                    if isempty(tpVOLother), return; end
                    facecolor=uisetcolor([.5 .5 0],'Select surface color');
                end
                if isempty(facecolor)||isequal(facecolor,0), return; end
                %answer=conn_menu_inputdlg({'Surface color (RGB)'},'',1,{mat2str(facecolor)});
                %if isempty(answer), return; end
                %answer=str2num(answer{1}); 
                %if numel(answer)==3, facecolor=answer; end
                for np=1:numel(tpVOLother), 
                    if ~isempty(tpVOLother(np).vertices)&&~isempty(tpVOLother(np).faces)
                        state.handles.patchblobother=cat(2,state.handles.patchblobother,patch(tpVOLother(np),'facecolor',facecolor,'edgecolor','none','alphadatamapping','none','FaceLighting', 'phong','facealpha',state.facealphaud,'tag','userdefined_surface','parent',state.handles.hax));
                        state.patchblobother_x=[state.patchblobother_x mean(tpVOLother(np).vertices(:,1))];
                        state.pVOLother=[state.pVOLother tpVOLother(np)];
                    end
                end
            case {'ax_ref', 'sag_ref', 'cor_ref', 'ref_ax', 'ref_sag', 'ref_cor', 'ref_all'}
                switch(option)
                    case {'ref_sag','sag_ref'}, nref=1;
                    case {'ref_cor','cor_ref'}, nref=2;
                    case {'ref_ax','ax_ref'},  nref=3;
                    case 'ref_all', nref=1:3;
                end
                axtemp=varargin{1};
                if ischar(axtemp), axtemp=strcmp(axtemp,'on'); end
                if axtemp, dosubstr='on'; else dosubstr='off'; end
                state.showaxref(nref)=axtemp;
                set(state.handles.patchref(nref),'visible',dosubstr);
            case {'ref_file','file_ref'}
                if isempty(varargin)
                    [tfilename,tpathname]=conn_fileutils('uigetfile','*.nii; *.img','Select reference file');
                    if ischar(tfilename), state.structural=fullfile(tpathname,tfilename);
                    else return;
                    end
                else state.structural=varargin{1};
                end
                if numel(varargin)>1&&~isempty(varargin{2}), state.structural_method=varargin{2}; end
                state.strvol=conn_fileutils('spm_vol',state.structural);
                state.structuraldata=conn_fileutils('spm_read_vols',state.strvol);
                %state.structuraldata(isnan(state.structuraldata))=0; %%
                [thpatch,state.Prange]=conn_mesh_display_axref(state.strvol,state.structuraldata,state.pointer_mm, state.structural_method);
                for tnview=1:3
                    set(state.handles.patchref(tnview),'faces',thpatch(tnview).faces,'vertices',thpatch(tnview).vertices,'facevertexcdata',thpatch(tnview).facevertexcdata,'facevertexalpha',state.facealpharef*thpatch(tnview).facevertexalpha,'userdata',thpatch(tnview).facevertexalpha,'facecolor','flat');
                end
                redrawnowcolorbar=true;
            case {'ref_method'}
                if isempty(varargin)
                    answer=conn_menu_inputdlg({'Reference axes position (x/y/z coordinates in mm)'},'',1,{num2str(state.pointer_mm)});
                    if isempty(answer), return; end
                    state.structural_method=answer{1};
                else state.structural_method=varargin{1};
                end
                if ischar(state.structural_method), state.structural_method=str2num(state.structural_method); end
                [thpatch,state.Prange]=conn_mesh_display_axref(state.strvol,state.structuraldata,state.pointer_mm, state.structural_method,state.Prange);
                for tnview=1:3
                    set(state.handles.patchref(tnview),'faces',thpatch(tnview).faces,'vertices',thpatch(tnview).vertices,'facevertexcdata',thpatch(tnview).facevertexcdata,'facevertexalpha',state.facealpharef*thpatch(tnview).facevertexalpha,'userdata',thpatch(tnview).facevertexalpha,'facecolor','flat');
                end
                redrawnowcolorbar=true;
            case {'ref_pos','pos_ref'}
                if isempty(varargin)
                    answer={''};
                    while numel(str2num(answer{1}))~=3
                        answer=conn_menu_inputdlg({'Reference axes position (x/y/z coordinates in mm)'},'',1,{num2str(state.pointer_mm)});
                        if isempty(answer), return; end
                    end
                    state.pointer_mm=str2num(answer{1});
                elseif ischar(varargin{1})
                    switch(varargin{1})
                        case 'x+', state.pointer_mm(1)=state.pointer_mm(1)+5;
                        case 'x-', state.pointer_mm(1)=state.pointer_mm(1)-5;
                        case 'y+', state.pointer_mm(2)=state.pointer_mm(2)+5;
                        case 'y-', state.pointer_mm(2)=state.pointer_mm(2)-5;
                        case 'z+', state.pointer_mm(3)=state.pointer_mm(3)+5;
                        case 'z-', state.pointer_mm(3)=state.pointer_mm(3)-5;
                    end
                else state.pointer_mm=varargin{1};
                end
                [thpatch,state.Prange]=conn_mesh_display_axref(state.strvol,state.structuraldata,state.pointer_mm, state.structural_method,state.Prange);
                for tnview=1:3
                    set(state.handles.patchref(tnview),'faces',thpatch(tnview).faces,'vertices',thpatch(tnview).vertices,'facevertexcdata',thpatch(tnview).facevertexcdata,'facevertexalpha',state.facealpharef*thpatch(tnview).facevertexalpha,'userdata',thpatch(tnview).facevertexalpha,'facealpha','flat');
                end
            case 'colorbar',
                if isempty(state.Vrange)&&isempty(state.Prange), conn_msgbox({'No activation overlay present','Use menu Surface.BrainSurface.BrainSurfaceActivation to define a new activation overlay'},'',2); return; end
                if strcmp(varargin{1},'rescale')
                    if isempty(state.Prange), Vrange=state.Vrange; else Vrange=state.Prange; end
                    if numel(varargin)>1, answ=varargin{2}; 
                    else answ=conn_menu_inputdlg({'Enter new colorbar limits:'},'Rescale colorbar',1,{mat2str(Vrange([1 end]),6)});
                        if ~isempty(answ), answ=str2num(answ{1}); end
                    end
                    if isequal(answ,'symmetricmin'), answ=min(abs(state.Vrange(state.Vrange~=0)))*[-1 1];
                    elseif isequal(answ,'symmetric')||isequal(answ,'symmetricmax'), answ=max(abs(state.Vrange))*[-1 1]; 
                    end
                    if ~isempty(answ)&&numel(answ)==2
                        if isempty(state.Prange), 
                            state.Vrange([1 end])=answ;
                            conn_mesh_display_refresh([],[],'remap&draw');
                        else
                            state.Prange([1 end])=answ;
                            [thpatch,state.Prange]=conn_mesh_display_axref(state.strvol,state.structuraldata,state.pointer_mm, state.structural_method, state.Prange);
                            for tnview=1:3
                                set(state.handles.patchref(tnview),'faces',thpatch(tnview).faces,'vertices',thpatch(tnview).vertices,'facevertexcdata',thpatch(tnview).facevertexcdata,'facevertexalpha',state.facealpharef*thpatch(tnview).facevertexalpha,'userdata',thpatch(tnview).facevertexalpha,'facecolor','flat');
                            end
                        end
                    end
                %elseif isempty(state.Vrange)
                %    set(state.handles.colorbar,'visible',varargin{1});
                else
                    set(state.handles.colorbar,'visible',varargin{1});
                    if numel(varargin)>1, ylabel(varargin{2},'parent',state.handles.colorbar(1)); end
                end
                redrawnowcolorbar=true;
            case 'repaint'
                if numel(varargin)>=1&&~isempty(varargin{1})
                    filenameSURF=varargin{1};
                    if numel(varargin)>=2&&~isempty(varargin{2}), THR=varargin{2}; end
                else
                    [tfilename,tpathname]=conn_fileutils('uigetfile','*.nii; *.img','Select file');
                    if ischar(tfilename), filenameSURF=fullfile(tpathname,tfilename);
                    else return;
                    end
                end
                if ~isempty(filenameSURF)&&ischar(filenameSURF), filenameSURF=conn_mesh_display_expand(filenameSURF); end
                if ~isempty(filenameSURF)
                    if isstruct(filenameSURF{1})
                        a=filenameSURF{1}.vol;
                        state.info.surf=filenameSURF{1}.filename;
                        filepath=fileparts(filenameSURF{1}.filename);
                    else
                        a=conn_fileutils('spm_vol',char(filenameSURF));
                        state.info.surf=char(filenameSURF);
                        filepath=fileparts(filenameSURF{1});
                    end
                    if ~isempty(filepath), state.defaultfilepath=filepath; end
                    V=conn_mesh_display_getV(filenameSURF,a,state.FSfolder);
                    state.reducedpatch=1;
                else
                    state.info.surf='none';
                    V=zeros(size(data.curv{1}));
                end
                if isnan(THR),
                    sthr=sort(abs(V(~isnan(V))));
                    V=V.*(abs(V)>sthr(round(.9*numel(sthr))));
                else
                    V=V.*(abs(V)>THR);
                end
                %if ~isnan(THR), V=V.*(V>THR); end
                if size(V,3)==2, V=abs(V(:,:,1))-abs(V(:,:,2));
                elseif size(V,3)>2, [nill,V]=max(abs(conn_surf_smooth(V)),[],3); V(nill==0)=0; cmapstart=rand(size(cmapstart));
                end
                if any(V(:)<0)&&any(V(:)>0),
                    state.dotwosided=true;
                else
                    state.dotwosided=false;
                end
                state.V0=V;
                state.Vrange=[];
                conn_mesh_display_refresh([],[],'remap&draw');
                
            case 'repaintborder'
                %values=round(state.V0);
                values=double(state.V0~=0);
                faces=state.selected_faces{1};
                border=zeros(size(values));
                A=spm_mesh_adjacency(faces);
                for n=1:2
                    for n1=1:size(border,1), border(n1,n)=any(values(n1,n)~=values(A(:,n1)>0,n)); end
                    %for n2=2:4, for n1=find(border(:,n)==n2-1)'; border(~border(:,n)&A(:,n1)>0,n)=n2; end; end
                end
                border=1+border;
                %border=(5-border)/5; % plot background
                state.dotwosided=false;
                state.V0=border;
                state.Vrange=[];
                conn_mesh_display_refresh([],[],'remap&draw');
                
            case 'position'
                v=get(state.handles.hax,'cameraposition');
                try, h=state.handles.light;
                catch, h=findobj(gcbf,'type','light');
                end
                try, vup=get(state.handles.hax,'cameraupvector');
                catch, vup=state.up;
                end
                if numel(h)>2, set(h(2:end-1),'position',v); end
                state.up=reshape(v(1:3),1,3)/norm(v(1:3));
                set(h(1),'position',1000*(v(:)'/norm(v)-.1*vup(:)'/norm(vup)));
                set(h(end),'position',1000*(v(:)'/norm(v)+.1*vup(:)'/norm(vup)));
                %disp(v/norm(v));
                %disp(state.up/norm(state.up));
                %disp([1000*(v(:)'/norm(v)+[1;-1]*state.up(:)')])
                if ~isempty(state.handles.sphplots_txt)
                    for n=1:numel(state.handles.sphplots_txt), set(state.handles.sphplots_txt(n),'position',state.sphplots.sph_xyz(n,:)+state.up*state.fontclose*state.sphplots.sph_r(n)); end
                end
                if ~isempty(state.handles.leftrightlabel)
                    if 1,%isequal(get(state.handles.hax,'visible'),'on')
                        if abs(v(1))/max(eps,norm(v))>.10, set(state.handles.leftrightlabel,'visible','off');
                        else set(state.handles.leftrightlabel,'visible','on','color',.5*[1 1 1]+(abs(v(1))/max(eps,norm(v))-.01)/(.10-.01)*(state.background-.5*[1 1 1]));
                        end
                    end
                end

            case 'print_fcn1'
                fcn=varargin{1};
                feval(fcn{1},[],[],fcn{2:end});
                conn_mesh_display_refresh([],[],'colorbar','on');
                
            case 'print'
                str=varargin{1};
                if numel(varargin)>1, pfilename=varargin{2}; options=varargin(3:end);
                else pfilename=fullfile(state.defaultfilepath,'print01.jpg'); options={};
                end
                back=state.reducedpatch;
                conn_mesh_display_refresh([],[],'brain',[],1);
                if str~=1&&isequal(get(state.handles.colorbar(1),'visible'),'on'), colorbaron=true; conn_mesh_display_refresh([],[],'colorbar','off'); addcolorbar=@(x){@conn_mesh_display_refresh,'print_fcn1',x};
                else colorbaron=false; addcolorbar=@(x)x;
                end
                switch(str)
                    case 1, conn_print(state.handles.hfig,pfilename,options{:});
                    case 2, conn_print(state.handles.hfig,pfilename,options{:},'-row',get(findobj(state.handles.hfig,'label','Left view'),'callback'),addcolorbar(get(findobj(state.handles.hfig,'label','Right view'),'callback')));
                    case 3, conn_print(state.handles.hfig,pfilename,options{:},'-mosaic3',get(findobj(state.handles.hfig,'label','Right view (both hem)'),'callback'),get(findobj(state.handles.hfig,'label','Superior view (flip)'),'callback'),addcolorbar(get(findobj(state.handles.hfig,'label','Posterior view'),'callback')));
                    case 4, conn_print(state.handles.hfig,pfilename,options{:},'-mosaic',get(findobj(state.handles.hfig,'label','Left view'),'callback'),get(findobj(state.handles.hfig,'label','Left medial view'),'callback'),get(findobj(state.handles.hfig,'label','Right view'),'callback'),addcolorbar(get(findobj(state.handles.hfig,'label','Right medial view'),'callback')));
                    case 5, conn_print(state.handles.hfig,pfilename,options{:},'-column',get(findobj(state.handles.hfig,'label','Left view'),'callback'),get(findobj(state.handles.hfig,'label','Right view'),'callback'),get(findobj(state.handles.hfig,'label','Left medial view'),'callback'),addcolorbar(get(findobj(state.handles.hfig,'label','Right medial view'),'callback')));
                    case 6, conn_print(state.handles.hfig,pfilename,options{:},'-row',get(findobj(state.handles.hfig,'label','Left view'),'callback'),get(findobj(state.handles.hfig,'label','Left medial view'),'callback'),get(findobj(state.handles.hfig,'label','Right medial view'),'callback'),addcolorbar(get(findobj(state.handles.hfig,'label','Right view'),'callback')));
                    case 7, conn_print(state.handles.hfig,pfilename,options{:},'-mosaic8',get(findobj(state.handles.hfig,'label','Left view'),'callback'),get(findobj(state.handles.hfig,'label','Left medial view'),'callback'),get(findobj(state.handles.hfig,'label','Anterior view'),'callback'),get(findobj(state.handles.hfig,'label','Superior view'),'callback'),get(findobj(state.handles.hfig,'label','Inferior view'),'callback'),get(findobj(state.handles.hfig,'label','Right view'),'callback'),get(findobj(state.handles.hfig,'label','Right medial view'),'callback'),addcolorbar(get(findobj(state.handles.hfig,'label','Posterior view'),'callback')));
                end
                if back~=1, conn_mesh_display_refresh([],[],'brain',[],back); end
        end
        if redrawnowcolorbar&&~isempty(state.handles.colorbar)
            if ~isempty(state.Prange)
                c1=linspace(-1,1,127)'; 
                cmap=cat(2,1-(c1<0)-c1.*max(0,c1-.5),max(0,1-2*abs(c1)),1-(c1>0)-c1.*min(0,c1+.5));
                %cmap=cat(2,.5+.5*tanh(50*(c1+.05))-c1.*max(0,c1-.5),max(0,1-2*abs(c1)),.5-.5*tanh(50*(c1-.05))-c1.*min(0,c1+.5));
                if sign(state.Prange(1))*sign(state.Prange(end))==-1,
                    set(state.handles.colorbar(1),'ytick',[.5,64.5,128.5],'yticklabel',arrayfun(@(x)num2str(x,'%.2f'),state.Prange,'uni',0),'ycolor',.6-.2*round(mean(state.background))*[1 1 1]);
                    set(state.handles.colorbar(2),'cdata',max(0,min(1, ind2rgb(round((size(cmap,1)+1)/2+emph*(size(cmap,1)-1)/2*linspace(-1,1,128)'),cmap))));
                elseif any(state.Prange>0)
                    set(state.handles.colorbar(1),'ytick',[.5,128.5],'yticklabel',arrayfun(@(x)num2str(x,'%.2f'),[min(state.Prange) max(state.Prange)],'uni',0),'ycolor',.6-.2*round(mean(state.background))*[1 1 1]);
                    set(state.handles.colorbar(2),'cdata',max(0,min(1, ind2rgb(round((size(cmap,1)+1)/2+1+emph*(size(cmap,1)-3)/2*linspace(0,1,128)'),cmap))));
                elseif any(state.Prange<0)
                    set(state.handles.colorbar(1),'ytick',[.5,128.5],'yticklabel',arrayfun(@(x)num2str(x,'%.2f'),[min(state.Prange) max(state.Prange)],'uni',0),'ycolor',.6-.2*round(mean(state.background))*[1 1 1]);
                    set(state.handles.colorbar(2),'cdata',max(0,min(1, ind2rgb(round((size(cmap,1)+1)/2-1+emph*(size(cmap,1)-3)/2*linspace(-1,0,128)'),cmap))));
                end
            elseif state.dotwosided,
                set(state.handles.colorbar(1),'ytick',[.5,64.5,128.5],'yticklabel',arrayfun(@(x)num2str(x,'%.2f'),state.Vrange,'uni',0),'ycolor',.6-.2*round(mean(state.background))*[1 1 1]);
                set(state.handles.colorbar(2),'cdata',max(0,min(1, ind2rgb(round((size(state.colormap,1)+1)/2+emph*(size(state.colormap,1)-1)/2*linspace(-1,1,128)'),state.colormap))));
            elseif any(state.Vrange>0)
                set(state.handles.colorbar(1),'ytick',[.5,128.5],'yticklabel',arrayfun(@(x)num2str(x,'%.2f'),[min(state.Vrange) max(state.Vrange)],'uni',0),'ycolor',.6-.2*round(mean(state.background))*[1 1 1]);
                set(state.handles.colorbar(2),'cdata',max(0,min(1, ind2rgb(round((size(state.colormap,1)+1)/2+emph*(size(state.colormap,1)-1)/2*linspace(0,1,128)'),state.colormap))));
            elseif any(state.Vrange<0)
                set(state.handles.colorbar(1),'ytick',[.5,128.5],'yticklabel',arrayfun(@(x)num2str(x,'%.2f'),[min(state.Vrange) max(state.Vrange)],'uni',0),'ycolor',.6-.2*round(mean(state.background))*[1 1 1]);
                set(state.handles.colorbar(2),'cdata',max(0,min(1, ind2rgb(round((size(state.colormap,1)+1)/2+emph*(size(state.colormap,1)-1)/2*linspace(-1,0,128)'),state.colormap))));
            end
        end
    end
end

function tfilenameSURF=conn_mesh_display_expand(filenameSURF);
tfilenameSURF=cellstr(conn_expandframe(filenameSURF));
try
    [trefnames,trefidx]=conn_roilabels(filenameSURF);
    if ~isempty(trefnames)
        if numel(trefnames)>1&&numel(tfilenameSURF)==1, %3d-atlas
            tempvol=conn_fileutils('spm_vol',char(tfilenameSURF));
            tempdata=conn_fileutils('spm_read_vols',tempvol);
            maxdata=max(tempdata(:));
            if max(trefidx)==maxdata
                idata=listdlg('liststring',trefnames,'selectionmode','multiple','initialvalue',1:numel(trefnames),'promptstring','Select ROI:','ListSize',[300 200]);
                if isempty(idata), return; end
                tfilenameSURF=num2cell(struct('filename',char(filenameSURF),'data',arrayfun(@(n)tempdata==trefidx(n),idata,'uni',0),'vol',tempvol,'label',reshape(trefnames(idata),1,[])));
            end
        elseif numel(trefnames)>1 && numel(trefnames)==numel(tfilenameSURF), %4d-atlas
            idata=listdlg('liststring',trefnames,'selectionmode','multiple','initialvalue',1:numel(trefnames),'promptstring','Select ROI:','ListSize',[300 200]);
            if isempty(idata), return; end
            try
                tempvol=conn_fileutils('spm_vol',char(tfilenameSURF));
                tempdata=conn_fileutils('spm_read_vols',tempvol);
                temp={}; for n=1:numel(idata)
                    ndata=idata(n);
                    temp{n}=struct('filename',char(tfilenameSURF{ndata}),'data',tempdata(:,:,:,ndata),'vol',tempvol(ndata),'label',trefnames{ndata});
                end
                tfilenameSURF=temp;
            catch, 
                tfilenameSURF=tfilenameSURF(trefidx(idata));
            end
        end
    end
end
end

function V=conn_mesh_display_getV(filenameSURF,a,FSfolder)
    if conn_surf_dimscheck(a(1).dim),%if isequal(a.dim,conn_surf_dims(8).*[1 1 2])
        if isstruct(filenameSURF{1})
            V=zeros(size(filenameSURF{1}.data));
            for n=1:numel(filenameSURF), V(filenameSURF{n}.data>0)=n; end
            V=reshape(V,[],2);
        else,
            V=conn_fileutils('spm_read_vols',a);
            V=reshape(V,[],2,numel(a));
        end
    else
        V=[];
        if isequal(FSfolder,fullfile(fileparts(which('conn')),'utils','surf')), files={'lh.pial.surf','rh.pial.surf','lh.mid.surf','rh.mid.surf','lh.white.surf','rh.white.surf'}; 
            %files={fullfile(FSfolder,'lh.mid.surf'),fullfile(FSfolder,'rh.mid.surf')}; % use this for accurate values
        else files={fullfile(FSfolder,'lh.pial.surf'),fullfile(FSfolder,'rh.pial.surf'),fullfile(FSfolder,'lh.mid.surf'),fullfile(FSfolder,'rh.mid.surf'),fullfile(FSfolder,'lh.white.surf'),fullfile(FSfolder,'rh.white.surf')};
        end
        for nfilenameSURF=1:numel(filenameSURF)
            tV=conn_surf_extract(filenameSURF{nfilenameSURF},files,[],10); % use this for smoother display
            tV=cell2mat(tV);
            maskV=isnan(tV)|(tV==0);
            tV(maskV)=0;
            tV=sum(reshape(tV,[],2,3),3)./max(1,sum(reshape(~maskV,[],2,3),3));                                                                                                                                                                                % use this for smoother display
            %tV=sum(reshape(tV,[],2,1),3)./max(1,sum(reshape(~maskV,[],2,1),3));                                                                                                                                                                               % use this for accurate values
            V=cat(3,V,tV);
        end
    end
end
        

function [out,h]=conn_mesh_display_patch(x,y,z,w,varargin)
persistent tstruct
if isempty(tstruct)
    Np=50;
    temp=[reshape([1:Np-1; 2:Np],1,[])' reshape([2:Np; 2*Np-1:-1:Np+1],1,[])' reshape([2*Np:-1:Np+2;2*Np:-1:Np+2],1,[])'];
    i1=[2*Np:-1:Np+1 3*Np:-1:2*Np+1];
    i2=[2*Np+1:3*Np Np:-1:1];
    tstruct=struct('vertices',zeros(3*Np,3), 'faces',fliplr([temp;i1(temp);i2(temp)]));
end
if nargin<4||isempty(w), w=ones(size(x)); end
w(1)=0;w(end)=0;
xyz=[x(:) y(:) z(:)];
dx=xyz(end,:)-xyz(1,:);
ndx=null(dx)'; 
dxyz1=sqrt(2*w)*ndx(1,:);
dxyz2=sqrt(2*w)*ndx(2,:);
tstruct.vertices=[xyz+dxyz1*sqrt(3)/2-dxyz2/2; flipud(xyz-dxyz1*sqrt(3)/2-dxyz2/2); xyz+dxyz2];
out=tstruct;
if nargout>1, h=patch(tstruct,varargin{:}); end
end

function varargout=conn_mesh_display_surf2patch(x,y,z,refslice,varargin)
if isempty(refslice), refslice=1; end
x1=[1.5*x(:,1,refslice)-.5*x(:,2,refslice) .5*x(:,1:end-1,refslice)+.5*x(:,2:end,refslice) 1.5*x(:,end,refslice)-.5*x(:,end-1,refslice)];
x1=[1.5*x1(1,:)-.5*x1(2,:);.5*x1(1:end-1,:)+.5*x1(2:end,:);1.5*x1(end,:)-.5*x1(end-1,:)];
y1=[1.5*y(:,1,refslice)-.5*y(:,2,refslice) .5*y(:,1:end-1,refslice)+.5*y(:,2:end,refslice) 1.5*y(:,end,refslice)-.5*y(:,end-1,refslice)];
y1=[1.5*y1(1,:)-.5*y1(2,:);.5*y1(1:end-1,:)+.5*y1(2:end,:);1.5*y1(end,:)-.5*y1(end-1,:)];
z1=[1.5*z(:,1,refslice)-.5*z(:,2,refslice) .5*z(:,1:end-1,refslice)+.5*z(:,2:end,refslice) 1.5*z(:,end,refslice)-.5*z(:,end-1,refslice)];
z1=[1.5*z1(1,:)-.5*z1(2,:);.5*z1(1:end-1,:)+.5*z1(2:end,:);1.5*z1(end,:)-.5*z1(end-1,:)];
vertices=[x1(:) y1(:) z1(:)];
a=reshape(1:numel(x1),size(x1));
faces=reshape(cat(3, a(1:end-1,1:end-1), a(2:end,1:end-1),a(2:end,2:end),a(1:end-1,2:end)),[],4);

x0={x,y,z};
for n=1:3,
    d3(n)=mean(mean(x0{n}(:,:,min(size(x0{n},3),2))-x0{n}(:,:,1),1),2);
end

for n=1:numel(varargin),
    if n==1
        varargout{n}=struct('vertices',zeros(size(vertices,1)*size(x,3),3),'faces',zeros(size(faces,1)*size(x,3),4),'facevertexcdata',zeros(size(faces,1)*size(x,3),1));
        for m=1:size(x,3)
            newvertices=vertices+repmat(d3*(m-refslice),size(vertices,1),1);
            newfaces=(m-1)*size(vertices,1)+faces;
            newc=reshape(varargin{n}(:,:,m),[],1);
            varargout{n}.faces((m-1)*size(faces,1)+(1:size(faces,1)),:)=newfaces;
            varargout{n}.vertices((m-1)*size(vertices,1)+(1:size(vertices,1)),:)=newvertices;
            varargout{n}.facevertexcdata((m-1)*size(faces,1)+(1:size(faces,1)),:)=newc;
        end
    else
        varargout{n}=varargout{1};
        for m=1:size(x,3)
            newc=reshape(varargin{n}(:,:,m),[],1);
            varargout{n}.facevertexcdata((m-1)*size(faces,1)+(1:size(faces,1)),:)=newc;
        end
    end
end
end

function [hpatch,Vrange] = conn_mesh_display_axref(strvol,structural,pointer_mm, proj_type, Vrange)
if nargin<4||isempty(proj_type), proj_type=1; end % 1: value at slice; 2: max value in direction orth to slice
if nargin<5, Vrange=[]; end
nslices=1;
dslices=1;
time=1;
pointer_vox=round(pinv(strvol.mat)*[pointer_mm(:)' 1]')';
[x,y,z]=ndgrid(1:strvol.dim(1),1:strvol.dim(2),1:strvol.dim(3));
xyz=strvol.mat*[x(:) y(:) z(:) ones(numel(x),1)]';
xyz_x=reshape(xyz(1,:),strvol.dim(1:3));
xyz_y=reshape(xyz(2,:),strvol.dim(1:3));
xyz_z=reshape(xyz(3,:),strvol.dim(1:3));
kslices=1;
if proj_type>1, pointer_vox(:)=1; end
for nview=1:3
    tslices=dslices*((1:nslices)-round((nslices+1)/2));
    switch nview
        case 1,
            if proj_type==1, tslices(pointer_vox(1)+tslices<1|pointer_vox(1)+tslices>size(xyz_x,1))=[];
            else tslices=[0 size(xyz_x,1)-1];
            end
            x=permute(xyz_x(max(1,min(size(xyz_x,1),pointer_vox(1)+kslices*tslices)),:,:),[2 3 1]);
            y=permute(xyz_y(max(1,min(size(xyz_x,1),pointer_vox(1)+tslices)),:,:),[2 3 1]);
            z=permute(xyz_z(max(1,min(size(xyz_x,1),pointer_vox(1)+tslices)),:,:),[2 3 1]);
            if proj_type==1, 
                z1=permute(structural(max(1,min(size(xyz_x,1),pointer_vox(1)+tslices)),:,:,min(size(structural,4),time)),[2 3 1]);
            else
                z1=permute(max(structural(:,:,:,min(size(structural,4),time)),[],1),[2 3 1]); 
                z1min=permute(min(structural(:,:,:,min(size(structural,4),time)),[],1),[2 3 1]);
                z1(-z1min>z1)=z1min(-z1min>z1); %x(:)=pointer_mm(1);
                z1=repmat(z1,[1 1 2]);
            end
        case 2,
            if proj_type==1, tslices(pointer_vox(2)+tslices<1|pointer_vox(2)+tslices>size(xyz_x,2))=[];
            else tslices=[0 size(xyz_x,2)-1];
            end
            x=permute(xyz_x(:,max(1,min(size(xyz_x,2),pointer_vox(2)+tslices)),:),[1 3 2]);
            y=permute(xyz_y(:,max(1,min(size(xyz_x,2),pointer_vox(2)+kslices*tslices)),:),[1 3 2]);
            z=permute(xyz_z(:,max(1,min(size(xyz_x,2),pointer_vox(2)+tslices)),:),[1 3 2]);
            if proj_type==1, 
                z1=permute(structural(:,max(1,min(size(xyz_x,2),pointer_vox(2)+tslices)),:,min(size(structural,4),time)),[1 3 2]);
            else
                z1=permute(max(structural(:,:,:,min(size(structural,4),time)),[],2),[1 3 2]); 
                z1min=permute(min(structural(:,:,:,min(size(structural,4),time)),[],2),[1 3 2]); 
                z1(-z1min>z1)=z1min(-z1min>z1); %y(:)=pointer_mm(2);
                z1=repmat(z1,[1 1 2]);
            end
        case 3,
            if proj_type==1, tslices(pointer_vox(3)+tslices<1|pointer_vox(3)+tslices>size(xyz_x,3))=[];
            else tslices=[0 size(xyz_x,3)-1];
            end
            x=permute(xyz_x(:,:,max(1,min(size(xyz_x,3),pointer_vox(3)+tslices))),[1 2 3]);
            y=permute(xyz_y(:,:,max(1,min(size(xyz_x,3),pointer_vox(3)+tslices))),[1 2 3]);
            z=permute(xyz_z(:,:,max(1,min(size(xyz_x,3),pointer_vox(3)+kslices*tslices))),[1 2 3]);
            if proj_type==1, 
                z1=permute(structural(:,:,max(1,min(size(xyz_x,3),pointer_vox(3)+tslices)),min(size(structural,4),time)),[1 2 3]);
            else
                z1=permute(max(structural(:,:,:,min(size(structural,4),time)),[],3),[1 2 3]); 
                z1min=permute(min(structural(:,:,:,min(size(structural,4),time)),[],3),[1 2 3]); 
                z1(-z1min>z1)=z1min(-z1min>z1); %z(:)=pointer_mm(3);
                z1=repmat(z1,[1 1 2]);
            end
    end
    if proj_type==1, z2=convn(convn(isnan(z1),conn_hanning(7)/4,'same'),conn_hanning(7)'/4,'same');
    else
        if ~isempty(Vrange)
            if all(Vrange>=0), minc1=Vrange(1); maxc1=Vrange(end);
            elseif all(Vrange<=0), minc1=-Vrange(2); maxc1=-Vrange(1);
            else minc1=Vrange(1); maxc1=Vrange(end);
            end
        else
            if nnz(z1<0)&nnz(z1>0), minc1=min(z1(:)); maxc1=max(z1(:));
            else minc1=min(abs(z1(z1~=0))); maxc1=max(abs(z1(:)));
            end
        end
        if 1,%state.smoother
            x=cat(3,interp2(x(:,:,1),2),interp2(x(:,:,2),2));
            y=cat(3,interp2(y(:,:,1),2),interp2(y(:,:,2),2));
            z=cat(3,interp2(z(:,:,1),2),interp2(z(:,:,2),2));
            z1=interp2(z1(:,:,1),2,'nearest');
            mz1=z1(:,:,1)~=0;
            z1=convn(convn(z1,conn_hanning(9)/5,'same'),conn_hanning(9)'/5,'same');
            mz1=convn(convn(mz1,conn_hanning(9)/5,'same'),conn_hanning(9)'/5,'same');
            z1=(mz1>.5).*z1./max(eps,mz1);
            z1=repmat(z1,[1 1 2]);
        end
        z2=1-abs(convn(convn(sign(max(-1e16,z1)),conn_hanning(7)/4,'same'),conn_hanning(7)'/4,'same')); z1(z1==0)=nan;
    end
    [f1,f2]=conn_mesh_display_surf2patch(x,y,z,find(tslices==0),z1,z2);
    c1=f1.facevertexcdata;
    if proj_type==1, 
        c1=(c1-min(c1))/max(eps,max(c1)-min(c1)); % min:black, max:white
        c=repmat(c1,[1 3]);
        hpatch(nview)=struct('faces',f1.faces,'vertices',f1.vertices,'facevertexcdata',c,'facevertexalpha',max(0,1-2*f2.facevertexcdata));
        Vrange=[];
    else
        %c1=c1/max(abs(c1)); % positive:red, negative:blue
        if nnz(z1<0)&nnz(z1>0), c1(c1>0)=c1(c1>0)/abs(maxc1); c1(c1<0)=c1(c1<0)/abs(minc1);
        else c1=sign(c1).*max(0,min(1,(abs(c1)-minc1)/max(eps,maxc1-minc1)));
        end
        c=cat(2,1-(c1<0)-c1.*max(0,c1-.5),max(0,1-2*abs(c1)),1-(c1>0)-c1.*min(0,c1+.5));
        %c=cat(2,.5+.5*tanh(50*(c1+.05))-c1.*max(0,c1-.5),max(0,1-2*abs(c1)),.5-.5*tanh(50*(c1-.05))-c1.*min(0,c1+.5));
        %c=cat(2,(c1>0),abs(c1),(c1<0));
        hpatch(nview)=struct('faces',f1.faces,'vertices',f1.vertices,'facevertexcdata',c,'facevertexalpha',.25*(c1~=0)+.75*abs(c1).*(1-1*f2.facevertexcdata)); %max(0,min(1, 2*abs(c1).*(1-.1*f2.facevertexcdata))));
        %hpatch(nview).faces=[hpatch(nview).faces;size(hpatch(nview).vertices,1)+hpatch(nview).faces];
        %hpatch(nview).vertices=[hpatch(nview).vertices;hpatch(nview).vertices]; hpatch(nview).vertices(end/2+1,nview)=size(xyz_x,nview);
        %hpatch(nview).facevertexcdata=[hpatch(nview).facevertexcdata;hpatch(nview).facevertexcdata];
        %hpatch(nview).facevertexalpha=[hpatch(nview).facevertexalpha;hpatch(nview).facevertexalpha];
        if isempty(Vrange)
            if ~any(c1<0), Vrange=[minc1 maxc1];
            elseif ~any(c1>0), Vrange=[-maxc1 -minc1];
            else Vrange=[minc1 0 maxc1];
            end
        end
    end
end
end

function conn_mesh_display_keypress(h,ARG,varargin)
if nargin>1&&isfield(ARG,'Character')
    minusplus='-+';
    conn_mesh_display_refresh('pos_ref',[lower(ARG.Character), minusplus(1+(length(ARG.Modifier)>0))]);
end
end

function map = parula(n)
%   Copyright 2013 The MathWorks, Inc.
if nargin < 1
   n = size(get(gcf,'Colormap'),1);
end
values = [0.2081 0.1663 0.5292;0.2091 0.1721 0.5411;0.2101 0.1779 0.553;0.2109 0.1837 0.565;0.2116 0.1895 0.5771;0.2121 0.1954 0.5892;0.2124 0.2013 0.6013;0.2125 0.2072 0.6135;0.2123 0.2132 0.6258;0.2118 0.2192 0.6381;0.2111 0.2253 0.6505;0.2099 0.2315 0.6629;0.2084 0.2377 0.6753;0.2063 0.244 0.6878;0.2038 0.2503 0.7003;0.2006 0.2568 0.7129;0.1968 0.2632 0.7255;0.1921 0.2698 0.7381;0.1867 0.2764 0.7507;0.1802 0.2832 0.7634;0.1728 0.2902 0.7762;0.1641 0.2975 0.789;0.1541 0.3052 0.8017;0.1427 0.3132 0.8145;0.1295 0.3217 0.8269;0.1147 0.3306 0.8387;0.0986 0.3397 0.8495;0.0816 0.3486 0.8588;0.0646 0.3572 0.8664;0.0482 0.3651 0.8722;0.0329 0.3724 0.8765;0.0213 0.3792 0.8796;0.0136 0.3853 0.8815;0.0086 0.3911 0.8827;0.006 0.3965 0.8833;0.0051 0.4017 0.8834;0.0054 0.4066 0.8831;0.0067 0.4113 0.8825;0.0089 0.4159 0.8816;0.0116 0.4203 0.8805;0.0148 0.4246 0.8793;0.0184 0.4288 0.8779;0.0223 0.4329 0.8763;0.0264 0.437 0.8747;0.0306 0.441 0.8729;0.0349 0.4449 0.8711;0.0394 0.4488 0.8692;0.0437 0.4526 0.8672;0.0477 0.4564 0.8652;0.0514 0.4602 0.8632;0.0549 0.464 0.8611;0.0582 0.4677 0.8589;0.0612 0.4714 0.8568;0.064 0.4751 0.8546;0.0666 0.4788 0.8525;0.0689 0.4825 0.8503;0.071 0.4862 0.8481;0.0729 0.4899 0.846;0.0746 0.4937 0.8439;0.0761 0.4974 0.8418;0.0773 0.5012 0.8398;0.0782 0.5051 0.8378;0.0789 0.5089 0.8359;0.0794 0.5129 0.8341;0.0795 0.5169 0.8324;0.0793 0.521 0.8308;0.0788 0.5251 0.8293;0.0778 0.5295 0.828;0.0764 0.5339 0.827;0.0746 0.5384 0.8261;0.0724 0.5431 0.8253;0.0698 0.5479 0.8247;0.0668 0.5527 0.8243;0.0636 0.5577 0.8239;0.06 0.5627 0.8237;0.0562 0.5677 0.8234;0.0523 0.5727 0.8231;0.0484 0.5777 0.8228;0.0445 0.5826 0.8223;0.0408 0.5874 0.8217;0.0372 0.5922 0.8209;0.0342 0.5968 0.8198;0.0317 0.6012 0.8186;0.0296 0.6055 0.8171;0.0279 0.6097 0.8154;0.0265 0.6137 0.8135;0.0255 0.6176 0.8114;0.0248 0.6214 0.8091;0.0243 0.625 0.8066;0.0239 0.6285 0.8039;0.0237 0.6319 0.801;0.0235 0.6352 0.798;0.0233 0.6384 0.7948;0.0231 0.6415 0.7916;0.023 0.6445 0.7881;0.0229 0.6474 0.7846;0.0227 0.6503 0.781;0.0227 0.6531 0.7773;0.0232 0.6558 0.7735;0.0238 0.6585 0.7696;0.0246 0.6611 0.7656;0.0263 0.6637 0.7615;0.0282 0.6663 0.7574;0.0306 0.6688 0.7532;0.0338 0.6712 0.749;0.0373 0.6737 0.7446;0.0418 0.6761 0.7402;0.0467 0.6784 0.7358;0.0516 0.6808 0.7313;0.0574 0.6831 0.7267;0.0629 0.6854 0.7221;0.0692 0.6877 0.7173;0.0755 0.6899 0.7126;0.082 0.6921 0.7078;0.0889 0.6943 0.7029;0.0956 0.6965 0.6979;0.1031 0.6986 0.6929;0.1104 0.7007 0.6878;0.118 0.7028 0.6827;0.1258 0.7049 0.6775;0.1335 0.7069 0.6723;0.1418 0.7089 0.6669;0.1499 0.7109 0.6616;0.1585 0.7129 0.6561;0.1671 0.7148 0.6507;0.1758 0.7168 0.6451;0.1849 0.7186 0.6395;0.1938 0.7205 0.6338;0.2033 0.7223 0.6281;0.2128 0.7241 0.6223;0.2224 0.7259 0.6165;0.2324 0.7275 0.6107;0.2423 0.7292 0.6048;0.2527 0.7308 0.5988;0.2631 0.7324 0.5929;0.2735 0.7339 0.5869;0.2845 0.7354 0.5809;0.2953 0.7368 0.5749;0.3064 0.7381 0.5689;0.3177 0.7394 0.563;0.3289 0.7406 0.557;0.3405 0.7417 0.5512;0.352 0.7428 0.5453;0.3635 0.7438 0.5396;0.3753 0.7446 0.5339;0.3869 0.7454 0.5283;0.3986 0.7461 0.5229;0.4103 0.7467 0.5175;0.4218 0.7473 0.5123;0.4334 0.7477 0.5072;0.4447 0.7482 0.5021;0.4561 0.7485 0.4972;0.4672 0.7487 0.4924;0.4783 0.7489 0.4877;0.4892 0.7491 0.4831;0.5 0.7491 0.4786;0.5106 0.7492 0.4741;0.5212 0.7492 0.4698;0.5315 0.7491 0.4655;0.5418 0.749 0.4613;0.5519 0.7489 0.4571;0.5619 0.7487 0.4531;0.5718 0.7485 0.449;0.5816 0.7482 0.4451;0.5913 0.7479 0.4412;0.6009 0.7476 0.4374;0.6103 0.7473 0.4335;0.6197 0.7469 0.4298;0.629 0.7465 0.4261;0.6382 0.746 0.4224;0.6473 0.7456 0.4188;0.6564 0.7451 0.4152;0.6653 0.7446 0.4116;0.6742 0.7441 0.4081;0.683 0.7435 0.4046;0.6918 0.743 0.4011;0.7004 0.7424 0.3976;0.7091 0.7418 0.3942;0.7176 0.7412 0.3908;0.7261 0.7405 0.3874;0.7346 0.7399 0.384;0.743 0.7392 0.3806;0.7513 0.7385 0.3773;0.7596 0.7378 0.3739;0.7679 0.7372 0.3706;0.7761 0.7364 0.3673;0.7843 0.7357 0.3639;0.7924 0.735 0.3606;0.8005 0.7343 0.3573;0.8085 0.7336 0.3539;0.8166 0.7329 0.3506;0.8246 0.7322 0.3472;0.8325 0.7315 0.3438;0.8405 0.7308 0.3404;0.8484 0.7301 0.337;0.8563 0.7294 0.3336;0.8642 0.7288 0.33;0.872 0.7282 0.3265;0.8798 0.7276 0.3229;0.8877 0.7271 0.3193;0.8954 0.7266 0.3156;0.9032 0.7262 0.3117;0.911 0.7259 0.3078;0.9187 0.7256 0.3038;0.9264 0.7256 0.2996;0.9341 0.7256 0.2953;0.9417 0.7259 0.2907;0.9493 0.7264 0.2859;0.9567 0.7273 0.2808;0.9639 0.7285 0.2754;0.9708 0.7303 0.2696;0.9773 0.7326 0.2634;0.9831 0.7355 0.257;0.9882 0.739 0.2504;0.9922 0.7431 0.2437;0.9952 0.7476 0.2373;0.9973 0.7524 0.231;0.9986 0.7573 0.2251;0.9991 0.7624 0.2195;0.999 0.7675 0.2141;0.9985 0.7726 0.209;0.9976 0.7778 0.2042;0.9964 0.7829 0.1995;0.995 0.788 0.1949;0.9933 0.7931 0.1905;0.9914 0.7981 0.1863;0.9894 0.8032 0.1821;0.9873 0.8083 0.178;0.9851 0.8133 0.174;0.9828 0.8184 0.17;0.9805 0.8235 0.1661;0.9782 0.8286 0.1622;0.9759 0.8337 0.1583;0.9736 0.8389 0.1544;0.9713 0.8441 0.1505;0.9692 0.8494 0.1465;0.9672 0.8548 0.1425;0.9654 0.8603 0.1385;0.9638 0.8659 0.1343;0.9623 0.8716 0.1301;0.9611 0.8774 0.1258;0.96 0.8834 0.1215;0.9593 0.8895 0.1171;0.9588 0.8958 0.1126;0.9586 0.9022 0.1082;0.9587 0.9088 0.1036;0.9591 0.9155 0.099;0.9599 0.9225 0.0944;0.961 0.9296 0.0897;0.9624 0.9368 0.085;0.9641 0.9443 0.0802;0.9662 0.9518 0.0753;0.9685 0.9595 0.0703;0.971 0.9673 0.0651;0.9736 0.9752 0.0597;0.9763 0.9831 0.0538];
P = size(values,1);
map = interp1(1:size(values,1), values, linspace(1,P,n), 'linear');
end
function c=fixedge(c)
k=linspace(1,0,size(c,1))'.^4;
c=repmat(1-k,1,size(c,2)).*c+repmat(k.*mean(c,2),1,size(c,2));
end

