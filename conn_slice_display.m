function fh=conn_slice_display(data,structural,defaultfilepath,actthr,titlestr)
% CONN_SLICE_DISPLAY slice display in CONN
%
% CONN_SLICE_DISPLAY(fileDATA) displays volume-level data in fileDATA (overlaid on default reference structural image -ICBM MNI 2009b NLIN asymmetric template-)
% CONN_SLICE_DISPLAY(fileDATA,fileSTRUCT) displays volume-level data in fileDATA overlaid on structural image fileSTRUCT
% CONN_SLICE_DISPLAY('',fileSTRUCT) displays structural image fileSTRUCT
%
%  h=CONN_SLICE_DISPLAY(...) returns function handle implementing all GUI functionality
%

global CONN_x CONN_gui;
if isempty(CONN_gui)||~isfield(CONN_gui,'font_offset'), conn_font_init; end
if nargin<2||isempty(structural), structural=''; end
if nargin<3||isempty(defaultfilepath), defaultfilepath=pwd; end
if nargin<4||isempty(actthr), actthr=[0 0]; end
if nargin<5||isempty(titlestr), titlestr=''; end
if isfield(CONN_gui,'slice_display_dovol'), DOVOL=CONN_gui.slice_display_dovol; 
else DOVOL=false;
end
if isfield(CONN_gui,'slice_display_skipbookmarkicons'), SKIPBI=CONN_gui.slice_display_skipbookmarkicons; 
else SKIPBI=false;
end
if isfield(CONN_gui,'slice_display_showroitype2asroitype1'), EXPACT=CONN_gui.slice_display_showroitype2asroitype1; 
else EXPACT=true;
end
if numel(actthr)<2, actthr(2)=0; end

fh=@(varargin)conn_slice_display_refresh([],[],varargin{:});
hmsg=conn_msgbox('Initializing. Please wait...','',-1);drawnow;
state.structural=structural;
state.refnames={};
state.actthr=actthr(1);
state.volthr=actthr(2);
V=[];
doinit=true;
if ~isempty(data)&&ischar(data)&&~isempty(regexp(data,'\.mat$|\.jpg$')), data=conn_loadmatfile(conn_prepend('',data,'.mat'),'state'); data=data.state; end
if isstruct(data)&&isfield(data,'faces'), state.surf=data; data=[]; else state.surf=[]; end
if isstruct(data)&&isfield(data,'structural'), % struct from getstate
    state=data;
    if ~conn_fileutils('isdir',state.defaultfilepath), state.defaultfilepath=pwd; end
    doinit=false;
elseif isstruct(data)  % conn_slice_display( struct:T,p,stats,dof,mat ...
    state.isstat=true;
    state.isvol=true;
    state.T=data.T;
    state.p=data.p;
    state.stats=data.stats;
    state.dof=data.dof;
    state.mat=data.mat;
    if isfield(data,'clusters'), state.clusters=data.clusters;
    else state.clusters=[];
    end
    state.supra=data.supra.*sign(state.T);
    state.size=size(state.supra);
    state.info.structural='none';
    if isfield(data,'filename'), state.info.vol=data.filename;
    else state.info.vol='manually defined';
    end
    state.nslices=1;%9;
    state.dslices=8;
elseif ~isempty(data) % conn_slice_display(datafile [,structural_file])
    state.isstat=false;
    state.isvol=true;
    try, V=conn_fileutils('spm_vol',char(data));
    catch, V=conn_fileutils('spm_vol',conn_dir(data,'-R'));
    end
    if isempty(V), 
        if ishandle(hmsg), delete(hmsg); end
        error('file not found %s',char(data)); 
    end
    if conn_surf_dimscheck(V(1)),  V=conn_surf_parent(V,'conn_slice_display',true); end % data is in surface space
    try
        [ok,nill,fsfiles]=conn_checkFSfiles(char(data),false);
        if ok, state.surf=reshape(conn_surf_readsurf(fsfiles([2,5,1,4]),[],fsfiles{7}),[2,2]); tV=conn_fileutils('spm_vol',fsfiles{7}); state.freesurfertransparency=double(max(max(abs(tV(1).mat-V(1).mat)))<1e-4); end
    end
    %V=V(1);
    if isequal(state.structural,0) % resample to reference grid?
        state.structural='';
        temp=spm_vol(fullfile(fileparts(which('conn')),'utils','surf','referenceT1_icbm.nii'));
        state.mat=temp.mat;
        state.size=temp.dim;
        [x,y,z]=ndgrid(1:state.size(1),1:state.size(2),1:state.size(3));
        xyz=state.mat*[x(:) y(:) z(:) ones(numel(x),1)]';
        if numel(V)==1, 
            txyz=pinv(V(1).mat)*xyz;
            state.supra=reshape(conn_fileutils('spm_sample_vol',V,txyz(1,:),txyz(2,:),txyz(3,:),1),state.size(1),state.size(2),state.size(3));
        else
            state.supra=reshape(conn_fileutils('spm_get_data',V,pinv(V(1).mat)*xyz)',state.size(1),state.size(2),state.size(3),[]);
        end
    else
        state.supra=conn_fileutils('spm_read_vols',V);
        state.mat=V(1).mat;
        state.size=size(state.supra);
    end
    state.info.structural='none';
    state.info.vol=char(data);
    state.T=state.supra;
    state.nslices=1;%9;
    state.dslices=8;
else                   % conn_slice_display([],structural_file)
    state.isstat=false;
    state.isvol=false;
    try, V=conn_fileutils('spm_vol',char(state.structural));
    catch, V=conn_fileutils('spm_vol',conn_dir(state.structural,'-R'));%
    end
    if isempty(V), 
        if ishandle(hmsg), delete(hmsg); end
        error('file not found %s',char(state.structural));
    end
    if conn_surf_dimscheck(V(1)), V=conn_surf_parent(V,'conn_slice_display'); end % surface
    state.info.structural=char(state.structural);
    try
        [ok,nill,fsfiles]=conn_checkFSfiles(char(state.structural),false);
        if ok, state.surf=reshape(conn_surf_readsurf(fsfiles([2,5,1,4]),[],fsfiles{7}),[2,2]); tV=conn_fileutils('spm_vol',fsfiles{7}); state.freesurfertransparency=double(max(max(abs(tV(1).mat-V(1).mat)))<1e-4); end
    end
    %V=V(1);
    state.structural=conn_fileutils('spm_read_vols',V);
    state.mat=V(1).mat;
    state.size=size(state.structural);
    state.info.vol='none';
    state.T=state.structural;
    state.nslices=16;
    state.dslices=8;
    state.refnames={};
end
try % read .txt labels from ROI files
    if ~state.isstat&&~isempty(V)&&~any(rem(state.T(:),1))
        [trefnames,trefidx]=conn_roilabels(V(1).fname);
        if ~isempty(trefnames)
            maxdata=max(state.T(:));
            if max(trefidx)==maxdata,
                state.refnames=trefnames; state.refidx=trefidx; state.reftype=1;
            elseif numel(trefnames)==size(state.T,4),
                state.refnames=trefnames; state.refidx=trefidx; state.reftype=2;
                if EXPACT
                    if ~isfield(state,'supra'), state.supra=state.T; end
                    [nill,state.supra]=max(state.supra,[],4);
                    state.supra(nill==0)=0;
                    state.size=size(state.supra);
                    if isequal(state.structural,state.T), state.structural=state.supra; end
                    state.T=state.supra;
                    state.reftype=1;
                end
            end
        end
    end
end
state.handles.fh=fh;
if ~isfield(CONN_x,'folders')||~isfield(CONN_x.folders,'bookmarks')||isempty(CONN_x.folders.bookmarks), state.dobookmarks=false;
else state.dobookmarks=true;
end
if ~isfield(state,'cameraviews')
    tci=[];
    state.cameraviews=[-1 0 .001; 0 -1 .001; 0 -.001 1];
    state.cameraviewnames={'sagittal','coronal','axial'};
    tc=abs(state.cameraviews*state.mat(1:3,1:3)*diag(1./sqrt(sum(state.mat(1:3,1:3).^2,1))))';
    for tcn=1:3, [nill,tci(tcn)]=max(tc(tcn,:)); tc(:,tci)=-1; end
    state.cameraviews=state.cameraviews(tci,:);       % rows = world-space camera coordinates (when selecting slice across x/y/z in image-space)
    state.cameraviewnames=state.cameraviewnames(tci); % names of each view (when selecting slice across x/y/z in image-space)
    %[nill,state.cameraviewdirs]=sort(tci);            % which view is sagittal/coronal/axial (i.e. slice across x/y/z in world-space)
    [nill,state.cameraviewdirs]=max(abs(state.cameraviews),[],2);state.cameraviewdirs=state.cameraviewdirs'; % most influencial x/y/z dimension in world-space (for each x/y/z dimension in image-space)
end

if doinit
    if isempty(state.structural),
        state.structural=fullfile(fileparts(which('conn')),'utils','surf','referenceT1_icbm.nii');
        if isfield(CONN_gui,'refs')&&isfield(CONN_gui.refs,'canonical')&&isfield(CONN_gui.refs.canonical,'filename')&&~isempty(CONN_gui.refs.canonical.filename)
            if ~isequal(CONN_gui.refs.canonical.filename,fullfile(fileparts(which('conn')),'utils','surf','referenceT1_trans.nii')), %handles conn defaults (high-res for slice display)
                state.structural=CONN_gui.refs.canonical.filename;
            end
        end
    end
    state.time=1;
    if numel(state.size)>3, state.endtime=state.size(4); else state.endtime=1; end
    if 0,%state.isstat
        state.view=[1 1 1 1 0];
        %if ~state.isstat, state.view(4)=1; end
        state.cameraview=[1 1 1];
    else
        state.view=[0 0 1 0 0];
        state.cameraview=state.cameraviews(3,:);%[0 -.001 1];
    end
    state.transparency=1;
    state.slice_transparency=1;
    state.background=[.95 .95 .9];%.14*[1 1 1];%[.2,.6,.7];
    if ~isempty(state.refnames), state.cmap=[0.78 0.42 0.33;0.98 0.69 0.8;0.35 0.72 0.57;0.72 0.89 0.58;0.51 0.57 0.027;0.95 0.81 0.084;0.054 0.97 0.6;0.52 0.92 0.63;0.56 0.66 0.34;0.031 0.7 0.94;0.73 0.066 0.018;0.76 0.88 0.39;0.45 0.52 0.062;0.22 0.96 0.91;0.31 0.49 0.25;0.62 0.56 0.43;0.48 0.96 0.8;0.27 0.45 0.68;0.74 0.034 0.38;0.76 0.15 0.86;0.064 0.74 0.73;0.49 0.74 0.54;0.052 0.0057 0.43;0.67 0.0079 0.23;0.99 0.98 0.93;0.83 0.068 0.34;0.28 0.42 0.85;0.36 0.61 0.098;0.59 0.59 0.84;0.69 0.83 0.59;0.032 0.92 0.74;0.71 0.88 0.98;0.19 0.45 0.57;0.42 0.81 0.21;0.95 0.43 0.54;0.96 0.17 0.87;0.71 0.93 0.22;0.76 0.025 0.24;0.26 0.73 0.79;0.094 0.68 0.085;0.94 0.23 0.4;0.95 0.84 0.53;0.19 0.61 0.31;0.48 0.99 0.17;0.67 0.06 0.19;0.84 0.51 0.62;0.51 0.064 0.69;0.99 0.46 0.59;0.41 0.2 0.75;0.62 0.02 0.93;0.94 0.61 0.3;0.46 0.23 0.27;0.24 0.2 0.24;0.38 0.081 0.56;0.59 0.097 0.31;0.86 0.97 0.24;0.077 0.2 0.42;0.89 0.013 0.89;0.97 0.79 0.56;0.28 0.6 0.098;0.77 0.64 0.43;0.033 0.57 0.19;0.9 0.25 0.66;0.66 0.42 0.071;0.35 0.13 0.81;0.99 0.4 0.37;0.92 0.37 0.23;0.99 0.86 0.93;0.94 0.99 0.8;0.64 0.26 0.36;0.74 0.76 0.9;0.31 0.65 0.72;0.21 0.0045 0.72;0.14 0.99 0.29;0.46 0.25 0.56;0.35 0.27 0.71;0.77 0.13 0.2;0.36 0.6 0.87;0.33 0.15 0.7;0.71 0.4 0.82;0.37 0.98 0.66;0.95 0.13 0.16;0.12 0.77 0.28;0.33 0.75 0.84;0.61 0.59 0.5;0.16 0.61 0.37;0.52 0.22 0.42;0.35 0.96 0.41;0.77 0.32 0.013;0.39 0.29 0.43;0.77 0.81 0.69;0.071 0.49 0.85;0.36 0.064 0.3;0.059 0.76 0.028;0.38 0.31 0.82;0.85 0.011 0.099;0.02 0.88 0.48;0.53 0.16 0.66;0.5 0.34 0.85;0.16 0.99 0.6;0.61 0.77 0.72;0.65 0.39 0.58;0.41 0.42 0.25;0.14 0.22 0.61;0.04 0.24 0.29;0.63 0.89 0.17;0.95 0.57 0.00033;0.77 0.33 0.3;0.36 0.49 0.097;0.83 0.34 0.48;0.023 0.74 0.98;0.86 0.92 0.73;0.26 0.68 0.13;0.043 0.12 0.47;0.0035 0.65 0.9;0.38 0.85 0.56;0.75 0.69 0.58;0.86 0.56 0.48;0.96 0.76 0.43;0.21 0.41 0.93;0.21 0.58 0.077;0.2 0.015 0.52;0.17 0.68 0.94;0.75 0.11 0.29;0.53 0.68 0.45;0.88 0.78 0.95;0.65 0.86 0.32;0.72 0.85 0.65;0.79 0.45 0.028;0.41 0.74 0.2;0.27 0.47 0.61;0.18 0.07 0.89;0.96 0.97 0.59;0.27 0.25 0.76;0.54 0.13 0.29;0.49 0.57 0.27;0.28 0.26 0.2;0.72 0.041 1;0.094 0.19 0.96;0.49 0.61 0.99;0.63 0.26 0.88;0.31 0.57 0.33;0.11 0.78 0.11;0.35 0.3 0.8;0.58 0.034 0.97;0.67 0.54 0.39;0.43 0.56 0.52;0.14 0.82 0.48;0.08 0.33 0.75;0.46 0.41 0.2;0.41 0.1 0.12;0.29 0.35 0.39;0.051 0.61 0.56;0.94 0.74 0.013;0.38 0.81 0.24;0.94 0.77 0.89;0.23 0.5 0.36;0.4 0.1 0.95;0.27 0.53 0.72;0.6 0.017 0.35;0.81 0.23 0.91;0.86 0.21 0.97;0.11 0.28 0.84;0.14 0.89 0.17;0.38 0.93 0.87;0.69 0.64 0.43;0.1 0.76 0.49;0.53 0.4 0.86;0.68 0.98 0.91;0.5 0.32 0.64;0.75 0.53 0.29;1 0.53 0.67;0.79 0.98 0.67;0.49 0.15 0.17;0.82 0.35 0.34;0.12 0.87 0.75;0.011 0.49 0.66;0.12 0.95 0.4;0.57 0.6 0.048;0.87 0.69 0.23;0.1 0.53 0.9;0.75 0.46 0.32;0.16 0.5 0.76;0.96 0.46 0.76;0.91 0.75 0.82;0.12 0.32 0.14;0.5 0.49 0.57;0.46 0.31 0.74;0.18 0.1 0.94;0.54 0.23 0.52;0.32 0.86 0.74;0.39 0.12 0.34;0.8 0.0032 0.35;0.81 0.77 0.076;0.15 0.29 0.78;0.85 0.94 0.85;0.47 0.15 0.35;0.48 0.58 0.77;0.83 0.42 0.034;0.57 0.66 0.092;0.96 0.86 0.79;0.015 0.8 0.42;0.45 0.12 0.36;0.42 0.91 0.85;0.28 0.44 0.051;0.18 0.5 0.059;0.23 0.51 0.56;0.51 0.64 0.59;0.48 0.6 0.75;0.11 0.032 0.94;0.63 0.062 0.39;0.11 0.24 0.21;0.97 0.62 0.43;0.85 0.32 0.35;0.51 0.6 0.73;0.19 0.63 0.87;0.38 0.034 0.25;0.89 0.096 0.23;0.92 0.73 0.15;0.2 0.85 0.2;0.61 0.24 0.78;0.66 0.77 0.19;0.11 0.28 0.89;0.97 0.86 0.57;0.55 0.67 0.62;0.005 0.51 0.46;0.31 0.98 0.89;0.42 0.7 0.35;0.35 0.38 0.63;0.27 0.87 0.83;0.31 0.82 0.49;0.52 0.76 0.13;0.32 0.28 0.62;0.85 0.67 0.47;0.4 0.27 0.069;0.35 0.31 0.03;0.48 0.97 0.56;0.033 0.65 0.52;0.24 0.65 0.27;0.59 0.3 0.65;0.92 0.87 0.49;0.47 0.36 0.25;0.014 0.11 0.29;0.19 0.91 0.59;0.68 0.27 0.5;0.92 0.44 0.35;0.22 0.33 0.57;0.75 0.044 0.33;0.31 0.015 0.78;0.75 0.27 0.092;0.4 0.67 0.42;0.062 0.69 0.23;0.57 0.16 0.086;0.15 0.39 0.88;0.72 0.37 0.045;0.74 0.29 0.7]; 
    else state.cmap=autumn(256); %[linspace(0,1,256)',zeros(256,2)]; 
    end
    state.blackistransparent=true;
    state.contourtransparency=0;
    if ~isfield(state,'freesurfertransparency'), state.freesurfertransparency=0; end
    state.expandmultiview=true;
    state.colorbrightness=0;
    state.colorcontrast=1;
    state.viewoverlay=1;
    state.viewslicetitle=1;
    state.defaultfilepath=defaultfilepath;
    state.bookmark_filename='';
    state.bookmark_descr='';
    [x,y,z]=ndgrid(1:state.size(1),1:state.size(2),1:state.size(3));
    xyz=state.mat*[x(:) y(:) z(:) ones(numel(x),1)]';
    state.xyz_x=reshape(xyz(1,:),state.size(1:3));
    state.xyz_y=reshape(xyz(2,:),state.size(1:3));
    state.xyz_z=reshape(xyz(3,:),state.size(1:3));
    state.xyz_range=[min(state.xyz_x(:)) max(state.xyz_x(:)); min(state.xyz_y(:)) max(state.xyz_y(:)); min(state.xyz_z(:)) max(state.xyz_z(:))];
    state.resamplestructural=false;
    if ischar(state.structural)
        V=conn_fileutils('spm_vol',state.structural);
        if conn_surf_dimscheck(V(1)), V=conn_surf_parent(V,'conn_slice_display'); end % surface
        state.info.structural=state.structural;
        try
            [ok,nill,fsfiles]=conn_checkFSfiles(state.structural,false);
            if ok, state.surf=reshape(conn_surf_readsurf(fsfiles([2,5,1,4]),[],fsfiles{7}),[2,2]); tV=conn_fileutils('spm_vol',fsfiles{7}); state.freesurfertransparency=double(max(max(abs(tV(1).mat-V(1).mat)))<1e-4); end
        end
        fact=abs(det(state.mat)/det(V(1).mat));
        while fact>1.1
            state.size(1:3)=2*state.size(1:3)-1;
            state.mat=state.mat*[.5 0 0 .5;0 .5 0 .5;0 0 .5 .5;0 0 0 1];
            try, state.supra=state.supra(round(1:.5:end),round(1:.5:end),round(1:.5:end),:); end
            try, state.T=state.T(round(1:.5:end),round(1:.5:end),round(1:.5:end),:); end
            try, state.p=state.p(round(1:.5:end),round(1:.5:end),round(1:.5:end),:); end
            fact=fact/8;
        end
        [x,y,z]=ndgrid(1:state.size(1),1:state.size(2),1:state.size(3));
        xyz=state.mat*[x(:) y(:) z(:) ones(numel(x),1)]';
        state.xyz_x=reshape(xyz(1,:),state.size(1:3));
        state.xyz_y=reshape(xyz(2,:),state.size(1:3));
        state.xyz_z=reshape(xyz(3,:),state.size(1:3));
        state.xyz_range=[min(state.xyz_x(:)) max(state.xyz_x(:)); min(state.xyz_y(:)) max(state.xyz_y(:)); min(state.xyz_z(:)) max(state.xyz_z(:))];
        state.resamplestructural=true;
        if numel(V)==1, 
            txyz=pinv(V(1).mat)*xyz;
            state.structural=reshape(conn_fileutils('spm_sample_vol',V,txyz(1,:),txyz(2,:),txyz(3,:),1),state.size(1),state.size(2),state.size(3));
        else
            state.structural=reshape(conn_fileutils('spm_get_data',V,pinv(V(1).mat)*xyz)',state.size(1),state.size(2),state.size(3),[]);
        end
        %V=V(1);
    end
    state.Srange=[min(state.structural(:)) max(state.structural(:))];
    state.Srange=state.Srange+diff(state.Srange)*.01*[-1 1];
    if isnan(state.volthr), 
        sthr=sort(state.structural(~isnan(state.structural)));
        state.volthr=sthr(round(.9*numel(sthr)));
        %state.volthr=state.Srange*[2/3;1/3]; 
    end
    state.pointer_mm=[0 0 18];
    if any(state.pointer_mm'<state.xyz_range(:,1))||any(state.pointer_mm'>state.xyz_range(:,2))
        state.pointer_mm=mean(state.xyz_range,2)';
    end
    state.pointer_vox=round([state.pointer_mm 1]*pinv(state.mat)'*[1 0 0;0 1 0;0 0 1;0 0 0]+.001);
    if state.isvol, 
        if isnan(state.actthr), 
            sthr=sort(abs(state.supra(~isnan(state.supra))));
            state.actthr=sthr(round(.9*numel(sthr)));
            %state.actthr=max(abs(state.supra(:)))/2; 
        end; 
        state.Vrange=[min([0,min(state.T(state.supra>state.actthr|state.supra<-state.actthr))]) max([0,max(state.T(state.supra>state.actthr|state.supra<-state.actthr))])]; end
    %if ~isempty(state.structural), state.structural(isnan(state.structural))=0; end
else
    [x,y,z]=ndgrid(1:state.size(1),1:state.size(2),1:state.size(3));
    xyz=state.mat*[x(:) y(:) z(:) ones(numel(x),1)]';
    state.xyz_x=reshape(xyz(1,:),state.size(1:3));
    state.xyz_y=reshape(xyz(2,:),state.size(1:3));
    state.xyz_z=reshape(xyz(3,:),state.size(1:3));
    state.xyz_range=[min(state.xyz_x(:)) max(state.xyz_x(:)); min(state.xyz_y(:)) max(state.xyz_y(:)); min(state.xyz_z(:)) max(state.xyz_z(:))];
    if state.isvol&&~isfield(state,'Vrange'), state.Vrange=[min([0,min(state.T(state.supra>state.actthr|state.supra<-state.actthr))]) max([0,max(state.T(state.supra>state.actthr|state.supra<-state.actthr))])]; end
    if ~isfield(state,'bookmark_filename'), state.bookmark_filename=''; end % backwards compatibility
    if ~isfield(state,'bookmark_descr'), state.bookmark_descr=''; end
    if ~isfield(state,'slice_transparency'), state.slice_transparency=1; end
    if ~isfield(state,'volthr'), state.volthr=0; end
    if ~isfield(state,'Srange'), state.Srange=[min(state.structural(:)) max(state.structural(:))]; state.Srange=state.Srange+diff(state.Srange)*.01*[-1 1]; end
    if ~isfield(state,'viewoverlay'), state.viewoverlay=1; end
    if ~isfield(state,'viewslicetitle'), state.viewslicetitle=1; end
    if ~isfield(state,'freesurfertransparency'), state.freesurfertransparency=0; end
    if ~isfield(state,'surf'), state.surf=[]; end
end

if DOVOL&&state.isvol
    for n=1:size(state.supra,4)
        pVOL1{n}=conn_surf_volume({state.supra(:,:,:,n),cat(4,state.xyz_x,state.xyz_y,state.xyz_z)},0,0,[],1,0,0);
        pVOL2{n}=conn_surf_volume({state.supra(:,:,:,n),cat(4,state.xyz_x,state.xyz_y,state.xyz_z)},0,0,[],1,0,1);
    end
end
state.handles.hfig=figure('units','norm','position',[.1 .25 .8 .5],'name',['conn slice display ',titlestr],'numbertitle','off','menubar','none','color',state.background,'tag','conn_slice_display','interruptible','off','busyaction','cancel','renderer','opengl','colormap',state.cmap,'visible','off');
uicontrol('style','frame','units','norm','position',[.5 .65 .5 .35],'foregroundcolor',[.5 .5 .5]);
uicontrol('style','frame','units','norm','position',[.5 0 .5 .65],'foregroundcolor',[.5 .5 .5]);
uicontrol('style','text','units','norm','position',[.55 .55 .4 .05],'string','Reference point','horizontalalignment','center','fontweight','bold');
uicontrol('style','text','units','norm','position',[.55 .475 .10 .07],'string','Coordinates (mm):');
fontsize=get(0,'defaultuicontrolfontsize')+2;
for n=1:3,
    state.handles.pointer_mm(n)=uicontrol('style','edit','units','norm','position',[.7+.05*(n-1) .50 .05 .05],'string',num2str(state.pointer_mm(n)),'fontsize',fontsize,'callback',{@conn_slice_display_refresh,'pointer_mm'});
end
for n=1:3,
    state.handles.pointer_mm_delta(2*n-1)=uicontrol('style','pushbutton','units','norm','position',[.7+.025*(2*n-2) .470 .025 .03],'string','-','callback',{@conn_slice_display_refresh,'pointer_mm','-',n});
    state.handles.pointer_mm_delta(2*n-0)=uicontrol('style','pushbutton','units','norm','position',[.7+.025*(2*n-1) .470 .025 .03],'string','+','callback',{@conn_slice_display_refresh,'pointer_mm','+',n});
end
uicontrol('style','text','units','norm','position',[.55 .375 .10 .07],'string','Coordinates (voxels):');
for n=1:3,
    state.handles.pointer_vox(n)=uicontrol('style','edit','units','norm','position',[.7+.05*(n-1) .40 .05 .05],'string',num2str(state.pointer_vox(n)),'fontsize',fontsize,'callback',{@conn_slice_display_refresh,'pointer_vox'});
end
for n=1:3,
    state.handles.pointer_vox_delta(2*n-1)=uicontrol('style','pushbutton','units','norm','position',[.7+.025*(2*n-2) .37 .025 .03],'string','-','callback',{@conn_slice_display_refresh,'pointer_vox','-',n});
    state.handles.pointer_vox_delta(2*n-0)=uicontrol('style','pushbutton','units','norm','position',[.7+.025*(2*n-1) .37 .025 .03],'string','+','callback',{@conn_slice_display_refresh,'pointer_vox','+',n});
end
state.handles.view(1)=uicontrol('style','checkbox','units','norm','position',[.55 .90 .15 .05],'string',sprintf('View yz plane (%s)',state.cameraviewnames{1}),'value',state.view(1),'callback',{@conn_slice_display_refresh,'view'});
state.handles.view(2)=uicontrol('style','checkbox','units','norm','position',[.55 .85 .15 .05],'string',sprintf('View xz plane (%s)',state.cameraviewnames{2}),'value',state.view(2),'callback',{@conn_slice_display_refresh,'view'});
state.handles.view(3)=uicontrol('style','checkbox','units','norm','position',[.55 .80 .15 .05],'string',sprintf('View xy plane (%s)',state.cameraviewnames{3}),'value',state.view(3),'callback',{@conn_slice_display_refresh,'view'});
state.handles.view(4)=uicontrol('style','checkbox','units','norm','position',[.55 .70 .15 .05],'string','View axis','value',state.view(4),'callback',{@conn_slice_display_refresh,'view'});
% if state.isvol
%     if state.isstat, str='View activation volume'; else str='View volume'; end
%     state.handles.view(5)=uicontrol('style','checkbox','units','norm','position',[.55 .75 .15 .05],'string',str,'value',state.view(5),'callback',{@conn_slice_display_refresh,'view'},'visible','off');
% end
state.handles.multislice1=uicontrol('style','edit','units','norm','position',[.85 .875 .05 .075],'string',num2str(state.nslices),'callback',{@conn_slice_display_refresh,'multislice'},'tooltipstring','<HTML>Maximum number of slices shown in multi-slice display<br/> - Set to 1 for single-slice display</HTML>');
state.handles.multislice2=uicontrol('style','edit','units','norm','position',[.90 .875 .05 .075],'string',num2str(state.dslices),'callback',{@conn_slice_display_refresh,'multislice'},'tooltipstring','Distance between displayed slices (in voxels)');
state.handles.multislice1_delta(1)=uicontrol('style','pushbutton','units','norm','position',[.85+.025*0 .845 .025 .03],'string','-','callback',{@conn_slice_display_refresh,'multislice','-',1});
state.handles.multislice1_delta(2)=uicontrol('style','pushbutton','units','norm','position',[.85+.025*1 .845 .025 .03],'string','+','callback',{@conn_slice_display_refresh,'multislice','+',1});
state.handles.multislice2_delta(1)=uicontrol('style','pushbutton','units','norm','position',[.90+.025*0 .845 .025 .03],'string','-','callback',{@conn_slice_display_refresh,'multislice','-',2});
state.handles.multislice2_delta(2)=uicontrol('style','pushbutton','units','norm','position',[.90+.025*1 .845 .025 .03],'string','+','callback',{@conn_slice_display_refresh,'multislice','+',2});
state.handles.multislice0=uicontrol('style','checkbox','units','norm','position',[.85 .925 .15 .05],'string','Multi-slice display','horizontalalignment','center','fontweight','bold','callback',{@conn_slice_display_refresh,'multisliceset'});
if all(state.nslices==1), set(state.handles.multislice0,'value',0); set([state.handles.multislice1 state.handles.multislice2 state.handles.multislice1_delta state.handles.multislice2_delta],'visible','off'); else set(state.handles.multislice0,'value',1); end
state.handles.text1=uicontrol('style','edit','units','norm','position',[.55 .20 .4 .05],'string','','horizontalalignment','center');
if state.isstat
    uicontrol('style','text','units','norm','position',[.55 .25 .4 .05],'string','Statistics','horizontalalignment','center','fontweight','bold');
    state.handles.text2=uicontrol('style','edit','units','norm','position',[.55 .15 .4 .05],'string','','horizontalalignment','center');
else state.handles.text2=[];
end
state.handles.mode=uicontrol('style','togglebutton','units','norm','position',[0 0 .5 .05],'string','Click on image to select reference point','value',1,'callback',{@conn_slice_display_refresh,'togglepointer'},'tooltipstring','switch between click-to-rotate and click-to-select behavior');
state.handles.gui=uicontrol('style','togglebutton','units','norm','position',[.5 0 .5 .05],'string','Hide GUI','value',0,'callback',{@conn_slice_display_refresh,'togglegui'},'tooltipstring','shows/hides GUI elements');

axes_handle=axes('units','norm','position',[.05 .05 .4 .9]);
state.handles.patch=[patch(0,0,0,'w') patch(0,0,0,'w') patch(0,0,0,'w') patch(0,0,0,'w')];
state.handles.slidetext=[];
%hold on; state.handles.patchcontour=[plot3(0,0,0,'k','linewidth',2) plot3(0,0,0,'k','linewidth',2) plot3(0,0,0,'k','linewidth',2)]; hold off;
if state.isvol
    if DOVOL 
        for n=1:numel(pVOL1), state.handles.act1(n)=patch(pVOL1{n}); end
        for n=1:numel(pVOL1), state.handles.act2(n)=patch(pVOL2{n}); end
        set([state.handles.act1],'edgecolor','none','facecolor','r','facealpha',state.transparency,'visible','off');
        set([state.handles.act2],'edgecolor','none','facecolor','b','facealpha',state.transparency,'visible','off');
    else state.handles.act1=[]; state.handles.act2=[];
    end
    hold on; state.handles.patchcontour1=[plot3(0,0,0,'k-','linewidth',2) plot3(0,0,0,'k-','linewidth',2) plot3(0,0,0,'k-','linewidth',2)]; state.handles.patchcontour2=[plot3(0,0,0,'k-','linewidth',2) plot3(0,0,0,'k-','linewidth',2) plot3(0,0,0,'k-','linewidth',2)]; hold off;
else state.handles.act1=[]; state.handles.act2=[];state.handles.patchcontour1=[]; state.handles.patchcontour2=[];
end
if ~isempty(state.surf)
    hold on; state.handles.patchcontour3=[plot3(0,0,0,'k-','linewidth',2) plot3(0,0,0,'k-','linewidth',2) plot3(0,0,0,'k-','linewidth',2)]; state.handles.patchcontour4=[plot3(0,0,0,'k-','linewidth',2) plot3(0,0,0,'k-','linewidth',2) plot3(0,0,0,'k-','linewidth',2)]; hold off;
    if ~state.freesurfertransparency, set([state.handles.patchcontour3(:)' state.handles.patchcontour4(:)'],'visible','off'); end
else state.handles.patchcontour3=[]; state.handles.patchcontour4=[];
end
set(state.handles.patch,'edgecolor','none');
hold on;
state.handles.line1=plot3(state.pointer_mm(1)+[0 0],state.pointer_mm(2)+[0 0],state.xyz_range(3,:),'b-');
state.handles.line2=plot3(state.pointer_mm(1)+[0 0],state.xyz_range(2,:),state.pointer_mm(3)+[0 0],'b-');
state.handles.line3=plot3(state.xyz_range(1,:),state.pointer_mm(2)+[0 0],state.pointer_mm(3)+[0 0],'b-');
hold off;
axis equal tight off;
state.handles.axes=gca;
state.handles.light=[light light];set(state.handles.light,'position',[1 1 1],'visible','off','color',.5*[1 1 1]);

if state.isvol, 
    axes('units','norm','position',[.453 .15 .015 .75]);
    temp=imagesc(max(0,min(1, ind2rgb(round((size(state.cmap,1)+1)/2+(size(state.cmap,1)-1)/2*linspace(-1,1,128)'),state.cmap))));
    set(gca,'ydir','normal','ytick',[],'xtick',[],'box','off','yaxislocation','left');
    temp2=text([1,1,1],[1-128*.05,64.5,128+128*.05],{' ',' ',' '},'horizontalalignment','center');
    state.handles.colorbar=[gca temp temp2(:)'];
    set(state.handles.colorbar,'visible','off');
else state.handles.colorbar=[]; 
end

state.handles.slider=uicontrol('style','slider','units','norm','position',[.47 .1 .025 .8],'callback',{@conn_slice_display_refresh,'pointer_mm','x'},'tooltipstring','Select reference slice');
if state.endtime>1, 
    state.handles.time=uicontrol('style','slider','units','norm','position',[.55 .05 .40 .05],'min',0,'max',1,'sliderstep',1/(state.endtime-1)*[1 2],'callback',{@conn_slice_display_refresh,'time'},'tooltipstring',sprintf('Volume/scan %d/%d',state.time,state.endtime)); 
    state.handles.timestr=uicontrol('style','text','units','norm','position',[.95 .05 .05 .05],'string',sprintf('%d/%d',state.time,state.endtime));
end
if state.isvol
    uicontrol('style','text','units','norm','position',[.69 .775 .15 .05],'string','Overlay threshold','horizontalalignment','right');
    state.handles.actthr=uicontrol('style','slider','units','norm','position',[.85 .775 .10 .05],'callback',{@conn_slice_display_refresh,'actthr'},'tooltipstring','Select overlay threshold');
    set(state.handles.actthr,'value',max(0,min(1, abs(state.actthr)/max(eps,max(abs(state.Vrange))))));
    state.handles.actthr_txt=uicontrol('style','text','units','norm','position',[.95 .775 .05 .05],'string',mat2str(state.actthr,2));
    try, addlistener(state.handles.actthr, 'ContinuousValueChange',@(varargin)conn_slice_display_refresh(state.handles.actthr,[],'actthr')); end
    state.handles.viewoverlay=uicontrol('style','checkbox','units','norm','position',[.55 .75 .15 .05],'string','View overlay','value',state.viewoverlay,'callback',{@conn_slice_display_refresh,'viewoverlay'});    
else state.handles.actthr=[]; state.handles.viewoverlay=[];
end
uicontrol('style','text','units','norm','position',[.69 .725 .15 .05],'string','Background threshold','horizontalalignment','right');
state.handles.volthr=uicontrol('style','slider','units','norm','position',[.85 .725 .10 .05],'callback',{@conn_slice_display_refresh,'volthr'},'tooltipstring','Select background image threshold');
set(state.handles.volthr,'value',max(0,min(1, (state.volthr-state.Srange(1))/max(eps,state.Srange(2)-state.Srange(1)))));
state.handles.volthr_txt=uicontrol('style','text','units','norm','position',[.95 .725 .05 .05],'string',mat2str(state.volthr,2));
uicontrol('style','text','units','norm','position',[.69 .675 .15 .05],'string','Background brightness','horizontalalignment','right');
state.handles.bright=uicontrol('style','slider','units','norm','position',[.85 .675 .10 .05],'callback',{@conn_slice_display_refresh,'bright'},'tooltipstring','Select background image brightness');
set(state.handles.bright,'value',max(0,min(1, .5+.5*tanh(state.colorbrightness))));
try, addlistener(state.handles.slider, 'ContinuousValueChange',@(varargin)conn_slice_display_refresh(state.handles.slider,[],'pointer_mm_refresh','x')); end
try, addlistener(state.handles.time, 'ContinuousValueChange',@(varargin)conn_slice_display_refresh(state.handles.time,[],'time')); end
try, addlistener(state.handles.bright, 'ContinuousValueChange',@(varargin)conn_slice_display_refresh(state.handles.bright,[],'bright')); end
try, addlistener(state.handles.volthr, 'ContinuousValueChange',@(varargin)conn_slice_display_refresh(state.handles.volthr,[],'volthr')); end
hc=state.handles.hfig;
hc1=uimenu(hc,'Label','Effects');
uimenu(hc1,'Label','white background','callback',{@conn_slice_display_refresh,'background',[1 1 1]},'tag','background');
uimenu(hc1,'Label','light background','callback',{@conn_slice_display_refresh,'background',[.95 .95 .9]},'tag','background');
uimenu(hc1,'Label','dark background','callback',{@conn_slice_display_refresh,'background',[.11 .11 .11]},'tag','background');
uimenu(hc1,'Label','black background','callback',{@conn_slice_display_refresh,'background',[0 0 0]},'tag','background');
uimenu(hc1,'Label','color background','callback',{@conn_slice_display_refresh,'background',[]},'tag','background','checked','on');
uimenu(hc1,'Label','multi-slice montage','separator','on','checked','on','callback',{@conn_slice_display_refresh,'expandmultiview','on'},'tag','expandmultiview');
uimenu(hc1,'Label','multi-slice stack','checked','off','callback',{@conn_slice_display_refresh,'expandmultiview','off'},'tag','expandmultiview');
%uimenu(hc1,'Label','brighter','separator','on','callback',{@conn_slice_display_refresh,'brighter',.25});
%uimenu(hc1,'Label','darker','callback',{@conn_slice_display_refresh,'brighter',-.25});
tvalues=[1 0];
thdl=uimenu(hc1,'Label','slice title on','separator','on','callback',{@conn_slice_display_refresh,'slice_title',1},'tag','slice_title');
thdl=[thdl,uimenu(hc1,'Label','slice title off','callback',{@conn_slice_display_refresh,'slice_title',0},'tag','slice_title')];
[nill,idx]=min(abs(state.viewslicetitle-tvalues));
set(thdl,'checked','off');set(thdl(max(1,min(numel(thdl),idx))),'checked','on');
tvalues=[1 .9:-.1:.1 .05 0];
thdl=uimenu(hc1,'Label','background image on','separator','on','callback',{@conn_slice_display_refresh,'slice_transparency',1},'tag','slice_transparency');
hc2=uimenu(hc1,'Label','background image transparent');
for n1=1:numel(tvalues)-1, thdl=[thdl,uimenu(hc2,'Label',num2str(n1-1),'callback',{@conn_slice_display_refresh,'slice_transparency',tvalues(n1)},'tag','slice_transparency')]; end
thdl=[thdl,uimenu(hc1,'Label','background image off','callback',{@conn_slice_display_refresh,'slice_transparency',0},'tag','slice_transparency')];
[nill,idx]=min(abs(state.slice_transparency-[1 tvalues]));
set(thdl,'checked','off');set(thdl(max(1,min(numel(thdl),idx))),'checked','on');
if state.isvol
    tvalues=[1 .9:-.1:.1 .05 0];
    thdl=uimenu(hc1,'Label','overlay image on','separator','on','callback',{@conn_slice_display_refresh,'act_transparency',1},'tag','act_transparency');
    hc2=uimenu(hc1,'Label','overlay image transparent');
    for n1=1:numel(tvalues)-1, thdl=[thdl,uimenu(hc2,'Label',num2str(n1-1),'callback',{@conn_slice_display_refresh,'act_transparency',tvalues(n1)},'tag','act_transparency')]; end
    thdl=[thdl,uimenu(hc1,'Label','overlay image off','callback',{@conn_slice_display_refresh,'act_transparency',0},'tag','act_transparency')];
    [nill,idx]=min(abs(state.transparency-[1 tvalues]));
    set(thdl,'checked','off');set(thdl(max(1,min(numel(thdl),idx))),'checked','on');
%     hc2=uimenu(hc1,'Label','overlay surface transparency','separator','on');
%     thdl=[];
%     for n1=0:.1:.9,thdl=[thdl,uimenu(hc2,'Label',num2str(1-n1),'callback',{@conn_slice_display_refresh,'act_transparency',n1},'tag','act_transparency')]; end
%     thdl=[thdl,uimenu(hc1,'Label','overlay surface opaque','callback',{@conn_slice_display_refresh,'act_transparency',1},'tag','act_transparency')];
%     set(thdl,'checked','off');set(thdl(max(1,min(numel(thdl),1+round(state.transparency*10)))),'checked','on');
    hc2=uimenu(hc1,'Label','overlay colormap','separator','on');
    for n1={'red','jet','hot','gray','bone','cool','hsv','spring','summer','autumn','winter','random','brighter','darker','equalize','manual','color'}
        uimenu(hc2,'Label',n1{1},'callback',{@conn_slice_display_refresh,'colormap',n1{1}},'tag','overlaycolormap');
    end
    uimenu(hc1,'Label','colorbar on','callback',{@conn_slice_display_refresh,'colorbar','on'},'tag','colorbar');
    uimenu(hc1,'Label','colorbar off','callback',{@conn_slice_display_refresh,'colorbar','off'},'tag','colorbar','checked','on');
    uimenu(hc1,'Label','colorbar limits','callback',{@conn_slice_display_refresh,'colorbar','rescale'});
    if DOVOL&& (any(cellfun('length',pVOL1))||any(cellfun('length',pVOL2)))
        uimenu(hc1,'Label','overlay volume on','separator','on','callback',{@conn_slice_display_refresh,'vol_transparency',1},'tag','vol_transparency');
        uimenu(hc1,'Label','overlay volume off','checked','on','callback',{@conn_slice_display_refresh,'vol_transparency',0},'tag','vol_transparency');
        hc2=uimenu(hc1,'Label','overlay volume lighting');
        uimenu(hc2,'Label','normal','callback',{@conn_slice_display_refresh,'material','dull'},'tag','material');
        uimenu(hc2,'Label','emphasis','callback',{@conn_slice_display_refresh,'material',[.1 .75 .5 1 .5]},'tag','material');
        uimenu(hc2,'Label','sketch','callback',{@conn_slice_display_refresh,'material',[.1 1 1 .25 0]},'tag','material');
        uimenu(hc2,'Label','shiny','callback',{@conn_slice_display_refresh,'material',[.3 .6 .9 20 1]},'tag','material');
        uimenu(hc2,'Label','metal','callback',{@conn_slice_display_refresh,'material',[.3 .3 1 25 .5]},'tag','material');
        uimenu(hc2,'Label','flat','callback',{@conn_slice_display_refresh,'material',[]},'tag','material','checked','on');
        uimenu(hc2,'Label','bright','callback',{@conn_slice_display_refresh,'light',.8},'separator','on','tag','light');
        uimenu(hc2,'Label','medium','callback',{@conn_slice_display_refresh,'light',.5},'tag','light','checked','on');
        uimenu(hc2,'Label','dark','callback',{@conn_slice_display_refresh,'light',.2},'tag','light');
    end
    uimenu(hc1,'Label','overlay contour on','separator','on','callback',{@conn_slice_display_refresh,'contour_transparency',1},'tag','contour_transparency');
    uimenu(hc1,'Label','overlay contour off','checked','on','callback',{@conn_slice_display_refresh,'contour_transparency',0},'tag','contour_transparency');
end
if ~isempty(state.surf)
    ht1=uimenu(hc1,'Label','freesurfer contour on','separator','on','callback',{@conn_slice_display_refresh,'freesurfer_transparency',1},'tag','freesurfer_transparency');
    ht2=uimenu(hc1,'Label','freesurfer contour off','callback',{@conn_slice_display_refresh,'freesurfer_transparency',0},'tag','freesurfer_transparency');
    if state.freesurfertransparency, set(ht1,'checked','on'); else set(ht2,'checked','on'); end
end
hc2=uimenu(hc1,'Label','smoother display','separator','on','checked','on','callback',{@conn_slice_display_refresh,'black_transparency','on'},'tag','black_transparency');
hc2=uimenu(hc1,'Label','raw data display','callback',{@conn_slice_display_refresh,'black_transparency','off'},'tag','black_transparency');
hc2=uimenu(hc1,'Label','info','separator','on','callback',{@conn_slice_display_refresh,'info'});
hc1=uimenu(hc,'Label','Print');
uimenu(hc1,'Label','current view','callback',{@conn_slice_display_refresh,'print'});
if state.dobookmarks
    hc1=uimenu(hc,'Label','Bookmark');
    hc2=uimenu(hc1,'Label','Save','callback',{@conn_slice_display_refresh,'bookmark'});
    if ~isempty(state.bookmark_filename),
        hc2=uimenu(hc1,'Label','Save as copy','callback',{@conn_slice_display_refresh,'bookmarkcopy'});
    end
end
set(state.handles.hfig,'userdata',state);%'uicontextmenu',hc,
try, 
    state.handles.rotate3d=rotate3d(state.handles.hfig);
    set(state.handles.rotate3d,'ActionPostCallback',{@conn_slice_display_refresh,'position'}); 
    set(state.handles.rotate3d,'enable','on'); 
catch
    state.handles.rotate3d=[];
end
        
conn_slice_display_refresh([],[],'init');
if ishandle(hmsg), delete(hmsg); end
set(state.handles.hfig,'visible','on');
try, set(state.handles.hfig,'resizefcn',{@conn_slice_display_refresh,'init'}); end

    function out=conn_slice_display_refresh(hObject,eventdata,option,varargin)
        out=[];
        try
            if numel(hObject)==1&&ishandle(hObject)&&~isempty(get(hObject,'tag'))
                str=get(hObject,'tag');
                set(findobj(state.handles.hfig,'tag',str),'checked','off');
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
        redrawnow=false;
        redrawnowcolorbar=false;
        switch(option)
            case 'none', return;
            case 'close', close(state.handles.hfig); return;
            case 'figurehandle',out=state.handles.hfig; return;
            case 'init',
                redrawnow=true;
            case 'getstate'
                out=state;
                out=rmfield(out,{'handles','xyz_x','xyz_y','xyz_z'});
            case {'bookmark','bookmarkcopy'},
                tfilename=[];
                if numel(varargin)>0&&~isempty(varargin{1}), tfilename=varargin{1}; 
                elseif ~isempty(state.bookmark_filename)&&strcmp(option,'bookmark'), tfilename=state.bookmark_filename;
                end
                if numel(varargin)>1&&~isempty(varargin{2}), descr=cellstr(varargin{2}); 
                else descr=state.bookmark_descr;
                end
                fcn=regexprep(mfilename,'^conn_','');
                conn_args={fcn,conn_slice_display_refresh([],[],'getstate')};
                [fullfilename,tfilename,descr]=conn_bookmark('save',...
                    tfilename,...
                    descr,...
                    conn_args);
                if isempty(fullfilename), return; end
                if ~SKIPBI, 
                    ht=conn_msgbox('Printing bookmark icon. Please wait...','',-1);
                    conn_slice_display_refresh([],[],'print',conn_prepend('',fullfilename,'.jpg'),'-nogui','-r50','-nopersistent'); 
                    if ishandle(ht), delete(ht); end
                end
                state.bookmark_filename=tfilename;
                state.bookmark_descr=descr;
                   %conn_args={fcn,conn_slice_display_refresh([],[],'getstate')}; % re-save to include bookmark info?
                   %save(conn_prepend('',fullfilename,'.mat'),'conn_args');
                if 0, conn_msgbox(sprintf('Bookmark %s saved',fullfilename),'',2);
                else out=fullfilename;
                end
                return;
            case 'title',
                set(state.handles.hfig,'name',varargin{1});
            case 'info',
                voxelsize=round(1e3*sqrt(sum(state.mat(1:3,1:3).^2,1)))/1e3;
                if nargout>0, out=state.info;
                else conn_msgbox([{'Background image:'},reshape(cellstr(state.info.structural),1,[]),{' ','Overlay image/contour:'},reshape(cellstr(state.info.vol),1,[]),{' ','Voxel size used for display (mm):'},{mat2str(voxelsize)}],'Slice display info');
                end
                return;
            case 'togglepointer'
                if get(state.handles.mode,'value')==1, set(state.handles.mode,'string','Click on image to select reference point');
                else set(state.handles.mode,'string','Click on image to rotate');
                end
                redrawnow=true;
            case 'togglegui',
                if numel(varargin)>0, onoff=varargin{1};
                else onoff=get(state.handles.gui,'value')==1;
                end
                if onoff, 
                    set(state.handles.gui,'string','Show GUI','value',1);
                    h=findobj(state.handles.hfig,'type','uicontrol');
                    set(h(~strcmp(get(h,'style'),'togglebutton')),'visible','off');
                    set(state.handles.axes,'units','norm','position',[.05 .05 .90 .9]);
                    set(state.handles.slider,'units','norm','position',[.97 .1 .025 .8]);
                    if ~isempty(state.handles.colorbar), set(state.handles.colorbar(1),'unit','norm','position',[.95 .15 .015 .75]); end
                else
                    set(state.handles.gui,'string','Hide GUI','value',0);
                    h=findobj(state.handles.hfig,'type','uicontrol');
                    set(h(~strcmp(get(h,'style'),'togglebutton')),'visible','on');
                    set(state.handles.axes,'units','norm','position',[.05 .05 .40 .9]);
                    set(state.handles.slider,'units','norm','position',[.47 .1 .025 .8]);
                    if ~isempty(state.handles.colorbar), set(state.handles.colorbar(1),'unit','norm','position',[.453 .15 .015 .75]); end
                end
                redrawnow=true;
            case {'pointer_mm','pointer_mm_refresh'}
                value=[str2num(get(state.handles.pointer_mm(1),'string')) str2num(get(state.handles.pointer_mm(2),'string')) str2num(get(state.handles.pointer_mm(3),'string'))];
                if numel(value)==3
                    if nargin>4,
                        if strcmp(varargin{1},'+'), d=1; else d=-1; end
                        npointer=varargin{2};
                    elseif nargin>3&&isequal(varargin{1},'x'),
                        v=get(state.handles.slider,'value');
                        npointer=state.cameraviewdirs(find(state.view(1:3),1));
                        d=round(state.xyz_range(npointer,1)+v*(state.xyz_range(npointer,2)-state.xyz_range(npointer,1)))-value(npointer);
                    elseif nargin>3&&numel(varargin{1})==3,
                        value=varargin{1}(:)';
                        npointer=1; d=0;
                    else
                        npointer=1; d=0;
                    end
                    value(npointer)=value(npointer)+d;
                    state.pointer_mm=value;
                    state.pointer_vox=round([state.pointer_mm 1]*pinv(state.mat)'*[1 0 0;0 1 0;0 0 1;0 0 0]+.001);
                end
                for n=1:3, set(state.handles.pointer_mm(n),'string',num2str(state.pointer_mm(n)));end
                for n=1:3, set(state.handles.pointer_vox(n),'string',num2str(state.pointer_vox(n)));end
                if strcmp(option,'pointer_mm_refresh'), redrawnow=2;
                else redrawnow=true;
                end
            case {'pointer_vox','pointer_vox_refresh'}
                value=[str2num(get(state.handles.pointer_vox(1),'string')) str2num(get(state.handles.pointer_vox(2),'string')) str2num(get(state.handles.pointer_vox(3),'string'))];
                if numel(value)==3
                    if nargin>4,
                        if strcmp(varargin{1},'+'), d=1; else d=-1; end
                        npointer=varargin{2};
                    elseif nargin>3&&isequal(varargin{1},'x'),
                        v=get(state.handles.slider,'value');
                        npointer=state.cameraviewdirs(find(state.view(1:3),1));
                        d=round(1+v*(state.size(npointer)-1))-value(npointer);
                    elseif nargin>3&&numel(varargin{1})==3,
                        value=varargin{1}(:)'; 
                        npointer=1; d=0;
                    else
                        npointer=1; d=0;
                    end
                    value(npointer)=value(npointer)+d;
                    state.pointer_vox=value;
                    state.pointer_mm=round([state.pointer_vox 1]*(state.mat)'*[1 0 0;0 1 0;0 0 1;0 0 0]);
                end
                for n=1:3, set(state.handles.pointer_mm(n),'string',num2str(state.pointer_mm(n)));end
                for n=1:3, set(state.handles.pointer_vox(n),'string',num2str(state.pointer_vox(n)));end
                if strcmp(option,'pointer_vox_refresh'), redrawnow=2;
                else redrawnow=true;
                end
            case 'time'
                value=get(state.handles.time,'value');
                state.time=round(1+(state.endtime-1)*max(0,min(1,value)));
                set(state.handles.time,'tooltipstring',sprintf('Volume/scan %d/%d',state.time,state.endtime));
                set(state.handles.timestr,'string',sprintf('%d/%d',state.time,state.endtime));
                redrawnow=true;
            case 'buttondown',
                p=get(gca,'cameraposition'); 
                pos=get(state.handles.axes,'currentpoint');
                pos=pos(1,1:3);
                mp=-inf;mpos=[];
                for nview=1:3,
                    if state.view(nview)
                        if 0,
                            txyz=get(state.handles.patch(nview),'vertices');
                            txyz=conn_bsxfun(@minus,txyz,pos);
                            ftxyz=(txyz/p)*p;
                            [mind,idx]=min(sum(abs(txyz-ftxyz).^2,2));
                            tpos=txyz(idx,:)+pos;
                            tp=-mind;
                            if tp>mp&&tp>-4, mpos=tpos; mp=tp; end
                        else
                            npointer=state.cameraviewdirs(nview);
                            k=(state.pointer_mm(npointer)-pos(npointer))/p(npointer);
                            tpos=pos+p*k;
                            tp=p*tpos';
                            if all(tpos>=state.xyz_range(:,1)')&&all(tpos<=state.xyz_range(:,2)')&&tp>mp, mpos=tpos; mp=tp; end
                        end
                    end
                end
                if ~isempty(mpos)
                    state.pointer_mm=round(mpos);
                    state.pointer_vox=round([state.pointer_mm 1]*pinv(state.mat)'*[1 0 0;0 1 0;0 0 1;0 0 0]+.001);
                end
%                 nview=find(state.view(1:3));
%                 if numel(nview)==1
%                     pos=get(state.handles.axes,'currentpoint');
%                     switch nview
%                         case 1, state.pointer_mm([2 3])=round(pos(1,[2 3]));
%                         case 2, state.pointer_mm([1 3])=round(pos(1,[1 3]));
%                         case 3, state.pointer_mm([1 2])=round(pos(1,[1 2]));
%                     end
%                     state.pointer_vox=round([state.pointer_mm 1]*pinv(state.mat)'*[1 0 0;0 1 0;0 0 1;0 0 0]+.001);
%                     %state.pointer_mm=round([state.pointer_vox 1]*(state.mat)'*[1 0 0;0 1 0;0 0 1;0 0 0]);
%                 end
                for n=1:3, set(state.handles.pointer_mm(n),'string',num2str(state.pointer_mm(n)));end
                for n=1:3, set(state.handles.pointer_vox(n),'string',num2str(state.pointer_vox(n)));end
                redrawnow=true;
            case 'view',
                if numel(varargin)>0, for nview=1:min(length(state.handles.view),numel(varargin{1})), set(state.handles.view(nview),'value',varargin{1}(nview)); end; end
                oldview=state.view;
                for nview=1:length(state.handles.view), state.view(nview)=get(state.handles.view(nview),'value'); end
                if any(oldview(1:3)~=state.view(1:3))
                    nview=find(state.view(1:3));
                    if (all(state.nslices==1)||state.expandmultiview)&&isequal(nview,1), state.cameraview=state.cameraviews(1,:); set(state.handles.mode,'value',1,'string','Click on image to select reference point');
                    elseif (all(state.nslices==1)||state.expandmultiview)&&isequal(nview,2), state.cameraview=state.cameraviews(2,:); set(state.handles.mode,'value',1,'string','Click on image to select reference point');
                    elseif (all(state.nslices==1)||state.expandmultiview)&&isequal(nview,3), state.cameraview=state.cameraviews(3,:); set(state.handles.mode,'value',1,'string','Click on image to select reference point');
                    %elseif nnz(abs(view*[1 0 0;0 1 0;0 0 1;0 0 0])>.01)==3, state.cameraview=[1 1 1]; set(state.handles.mode,'value',0,'string','Click on image to rotate');
                    else state.cameraview=[1 1 1]; set(state.handles.mode,'value',0,'string','Click on image to rotate');
                    end
                end
                redrawnow=true;
            case 'multisliceset',
                if numel(varargin)>0&&~isempty(varargin{1}), set(state.handles.multislice0,'value',varargin{1}); end
                if get(state.handles.multislice0,'value')==1, 
                    set([state.handles.multislice1 state.handles.multislice2 state.handles.multislice1_delta state.handles.multislice2_delta],'visible','on');
                    state.nslices=16;
                    set(state.handles.multislice1,'string',num2str(state.nslices));
                else
                    set([state.handles.multislice1 state.handles.multislice2 state.handles.multislice1_delta state.handles.multislice2_delta],'visible','off');
                    state.nslices=1;
                    set(state.handles.multislice1,'string',num2str(state.nslices));
                end
                if numel(varargin)>1&&~isempty(varargin{2}), state.nslices=varargin{2}; set(state.handles.multislice1,'string',num2str(state.nslices)); end
                if numel(varargin)>2&&~isempty(varargin{3}), state.dslices=varargin{3}; set(state.handles.multislice2,'string',num2str(state.dslices)); end
                nview=find(state.view(1:3));
                if (all(state.nslices==1)||state.expandmultiview)&&isequal(nview,1), state.cameraview=state.cameraviews(1,:); set(state.handles.mode,'value',1,'string','Click on image to select reference point');
                elseif (all(state.nslices==1)||state.expandmultiview)&&isequal(nview,2), state.cameraview=state.cameraviews(2,:); set(state.handles.mode,'value',1,'string','Click on image to select reference point');
                elseif (all(state.nslices==1)||state.expandmultiview)&&isequal(nview,3), state.cameraview=state.cameraviews(3,:); set(state.handles.mode,'value',1,'string','Click on image to select reference point');
                else state.cameraview=[1 1 1]; set(state.handles.mode,'value',0,'string','Click on image to rotate');
                end
                redrawnow=true;
            case 'multislice'
                npointer=0;
                if nargin>4,
                    if strcmp(varargin{1},'+'), d=1; else d=-1; end
                    npointer=varargin{2};
                end
                val=str2num(get(state.handles.multislice1,'string'));
                if npointer==1, val=val+d; end
                if ~isempty(val), state.nslices=min(1e3,max(1,round(val))); end
                set(state.handles.multislice1,'string',num2str(state.nslices));
                val=str2num(get(state.handles.multislice2,'string'));
                if npointer==2, val=val+d; end
                if ~isempty(val), state.dslices=max(1,round(val)); end
                set(state.handles.multislice2,'string',num2str(state.dslices));
                redrawnow=true;
            case 'viewoverlay'
                if numel(varargin)>0, state.viewoverlay=varargin{1};
                else state.viewoverlay=get(state.handles.viewoverlay,'value');
                end
                set(state.handles.viewoverlay,'value',state.viewoverlay);
                if state.isvol&&state.viewoverlay, set([state.handles.actthr state.handles.actthr_txt],'enable','on'); 
                else set([state.handles.actthr state.handles.actthr_txt],'enable','off'); 
                end
                redrawnow=true;
            case 'actthr'
                if numel(varargin)>0, state.actthr=varargin{1};
                else state.actthr=max(eps,max(abs(state.Vrange)))*get(state.handles.actthr,'value');
                end
                set(state.handles.actthr,'value',max(0,min(1, abs(state.actthr)/max(eps,max(abs(state.Vrange))))));
                set(state.handles.actthr_txt,'string',mat2str(state.actthr,2));
                redrawnow=true;
            case 'volthr'
                if numel(varargin)>0, state.volthr=varargin{1};
                else state.volthr=state.Srange(1)+max(eps,state.Srange(2)-state.Srange(1))*get(state.handles.volthr,'value');
                end
                set(state.handles.volthr,'value',max(0,min(1, (state.volthr-state.Srange(1))/max(eps,state.Srange(2)-state.Srange(1)))));
                set(state.handles.volthr_txt,'string',mat2str(state.volthr,2));
                redrawnow=true;
                
            case 'colorbar',
                if strcmp(varargin{1},'rescale')
                    if numel(varargin)>1, answ=varargin{2}; 
                    else answ=conn_menu_inputdlg({'Enter new colorbar limits:'},'Rescale colorbar',1,{mat2str(state.Vrange([1 end]),6)});
                        if ~isempty(answ), answ=str2num(answ{1}); end
                    end
                    if ~isempty(answ)&&numel(answ)==2
                        state.Vrange([1 end])=answ;
                        redrawnow=true;
                    end
                    set(state.handles.actthr,'value',max(0,min(1, abs(state.actthr)/max(eps,max(abs(state.Vrange))))));
                    redrawnow=true;
                else
                    set(state.handles.colorbar(2:end),'visible',varargin{1});
                end
                redrawnowcolorbar=true;
            case 'material'
                str=varargin{1};
                if isempty(str), set(state.handles.light,'visible','off');
                else
                    set(state.handles.light,'visible','on');
                    material(str);
                end
            case 'light'
                scale=varargin{1};
                set(state.handles.light,'color',scale*[1 1 1]);
            case 'colormap'
                cmap=varargin{1};
                if ischar(cmap)
                    switch(cmap)
                        case 'red', cmap=[linspace(0,1,256)',zeros(256,2)];
                        case 'hot', cmap=hot(256);
                        case 'jet', cmap=fixedge(jet(256));
                        case 'gray', cmap=gray(256);
                        case 'bone', cmap=bone(256);
                        case 'cool',cmap=fixedge(cool(256));
                        case 'hsv',cmap=fixedge(hsv(256));
                        case 'spring',cmap=repmat(1-linspace(1,0,256)'.^8,[1,3]).*spring(256);
                        case 'summer',cmap=repmat(1-linspace(1,0,256)'.^8,[1,3]).*summer(256);
                        case 'autumn',cmap=repmat(1-linspace(1,0,256)'.^8,[1,3]).*autumn(256);
                        case 'winter',cmap=repmat(1-linspace(1,0,256)'.^8,[1,3]).*winter(256);
                        case 'random',cmap=rand(256,3);
                        case 'brighter',cmap=min(1,1/sqrt(.95)*get(state.handles.hfig,'colormap').^(1/2)); 
                        case 'darker',cmap=.95*get(state.handles.hfig,'colormap').^2; 
                        case 'equalize',
                            sp=sort(abs(state.supra(~isnan(state.supra)&state.supra~=0)));
                            sp=sp(1)+[0;cumsum(diff(sp)+1e-10)];
                            isp=interp1([min(sp(1)-1e-10,0); sp],linspace(1,size(state.cmap,1),numel(sp)+1),linspace(0,sp(end),256));
                            cmap=interp1((1:size(state.cmap,1))',state.cmap,isp);
                        case 'manual',answer=conn_menu_inputdlg({'colormap (256x3)'},'',1,{mat2str(state.cmap)});if ~isempty(answer), answer=str2num(answer{1}); end;if ~any(size(answer,1)==[256]), return; end;cmap=max(0,min(1,answer));
                        case 'color',cmap=uisetcolor([],'Select color'); if isempty(cmap)||isequal(cmap,0), return; end; 
                        otherwise, disp('unknown value');
                    end
                end
                if size(cmap,2)<3, cmap=cmap(:,min(size(cmap,2),1:3)); end
                if size(cmap,1)==1, cmap=repmat(cmap,2,1); end
                state.cmap=cmap;
                set(state.handles.hfig,'colormap',cmap);
                redrawnow=true;
                redrawnowcolorbar=true;
            case 'background'
                if numel(varargin)>0&&~isempty(varargin{1}), color=varargin{1};
                else color=uisetcolor(state.background,'Select color'); if isempty(color)||isequal(color,0), return; end; 
                end
                state.background=color;
                set(state.handles.hfig,'color',color);
                redrawnowcolorbar=true;
            case 'slice_title'
                scale=varargin{1};
                state.viewslicetitle=max(0,min(1, scale));
                redrawnow=true;
            case 'slice_transparency'
                scale=varargin{1};
                state.slice_transparency=max(eps,scale);
                redrawnow=true;
            case 'act_transparency'
                scale=varargin{1};
                state.transparency=max(eps,scale);
                set([state.handles.act1 state.handles.act2],'facealpha',state.transparency);
                redrawnow=true;
            case 'black_transparency'
                str=varargin{1};
                if strcmp(str,'on'), state.blackistransparent=true; 
                else state.blackistransparent=false; 
                end
                redrawnow=true;
            case 'vol_transparency'
                scale=varargin{1};
                state.view(5)=max(0,min(1,round(scale)));
                %set(state.handles.view(5),'value',state.view(5));
                redrawnow=true;
            case 'contour_transparency'
                scale=varargin{1};
                state.contourtransparency=scale;
                redrawnow=true;
            case 'freesurfer_transparency'
                scale=varargin{1};
                state.freesurfertransparency=scale;
                redrawnow=true;
            case 'bright'
                if numel(varargin)>0, scale=varargin{1};
                else scale=get(state.handles.bright,'value'); 
                end
                state.colorbrightness=atanh(max(-1+1e-4,min(1-1e-4,2*scale-1)));
                redrawnow=true;
            case 'brighter'
                scale=varargin{1};
                state.colorbrightness=state.colorbrightness+scale;
                set(state.handles.bright,'value',max(0,min(1, .5+.5*tanh(state.colorbrightness))));
                redrawnow=true;
            case 'contrast'
                scale=varargin{1};
                state.colorcontrast=state.colorcontrast+scale;
                redrawnow=true;
            case 'expandmultiview',
                str=varargin{1};
                if strcmp(str,'on'), state.expandmultiview=true; 
                else state.expandmultiview=false; 
                end
                nview=find(state.view(1:3));
                if (all(state.nslices==1)||state.expandmultiview)&&isequal(nview,1), state.cameraview=state.cameraviews(1,:); set(state.handles.mode,'value',1,'string','Click on image to select reference point');
                elseif (all(state.nslices==1)||state.expandmultiview)&&isequal(nview,2), state.cameraview=state.cameraviews(2,:); set(state.handles.mode,'value',1,'string','Click on image to select reference point');
                elseif (all(state.nslices==1)||state.expandmultiview)&&isequal(nview,3), state.cameraview=state.cameraviews(3,:); set(state.handles.mode,'value',1,'string','Click on image to select reference point');
                else state.cameraview=[1 1 1]; set(state.handles.mode,'value',0,'string','Click on image to rotate');
                end
                redrawnow=true;
            case 'position'
                p=get(state.handles.axes,'cameraposition'); 
                set(findobj(gcbf,'type','light'),'position',p);
                state.cameraview=[];
            case 'print'
                set([state.handles.gui state.handles.mode state.handles.text1 state.handles.text2],'visible','off');
                value=get(state.handles.gui,'value');
                set(state.handles.gui,'value',1);
                if numel(varargin)>0, conn_print(state.handles.hfig,varargin{:});
                else conn_print(state.handles.hfig,fullfile(state.defaultfilepath,'print01.jpg'));
                end
                set([state.handles.gui state.handles.mode],'visible','on');
                set(state.handles.gui,'value',value);
        end
        
        if redrawnow
            %tslices=state.dslices*((1:state.nslices)-1);
            kslices=1;
            minmax=[];
            cmap0=get(state.handles.hfig,'colormap');
            ccolor=get(state.handles.hfig,'color');
            dcamera=get(state.handles.axes,'cameraposition')-get(state.handles.axes,'cameratarget'); dcamera=dcamera(:)/max(eps,norm(dcamera));
            cmap=@(alpha,rgb)interp1(linspace(0,1,size(cmap0,1)),cmap0(:,rgb),max(0,min(1,alpha)));
            for nview=1:3
                if state.view(nview)
                    tslices=state.dslices(min(numel(state.dslices),nview))*((1:state.nslices(min(numel(state.nslices),nview)))-round((state.nslices(min(numel(state.nslices),nview))+1)/2));
                    switch nview
                        case 1,
                            tslices(state.pointer_vox(1)+tslices<1|state.pointer_vox(1)+tslices>state.size(1))=[];
                            x=permute(state.xyz_x(max(1,min(state.size(1),state.pointer_vox(1)+kslices*tslices)),:,:),[2 3 1]);
                            y=permute(state.xyz_y(max(1,min(state.size(1),state.pointer_vox(1)+tslices)),:,:),[2 3 1]);
                            z=permute(state.xyz_z(max(1,min(state.size(1),state.pointer_vox(1)+tslices)),:,:),[2 3 1]);
                            z1=permute(state.structural(max(1,min(state.size(1),state.pointer_vox(1)+tslices)),:,:,min(size(state.structural,4),state.time)),[2 3 1]);
                            if state.isvol, z2=permute(state.supra(max(1,min(state.size(1),state.pointer_vox(1)+tslices)),:,:,min(size(state.supra,4),state.time)),[2 3 1]); end
                        case 2,
                            tslices(state.pointer_vox(2)+tslices<1|state.pointer_vox(2)+tslices>state.size(2))=[];
                            x=permute(state.xyz_x(:,max(1,min(state.size(2),state.pointer_vox(2)+tslices)),:),[1 3 2]);
                            y=permute(state.xyz_y(:,max(1,min(state.size(2),state.pointer_vox(2)+kslices*tslices)),:),[1 3 2]);
                            z=permute(state.xyz_z(:,max(1,min(state.size(2),state.pointer_vox(2)+tslices)),:),[1 3 2]);
                            z1=permute(state.structural(:,max(1,min(state.size(2),state.pointer_vox(2)+tslices)),:,min(size(state.structural,4),state.time)),[1 3 2]);
                            if state.isvol, z2=permute(state.supra(:,max(1,min(state.size(2),state.pointer_vox(2)+tslices)),:,min(size(state.supra,4),state.time)),[1 3 2]); end
                        case 3,
                            tslices(state.pointer_vox(3)+tslices<1|state.pointer_vox(3)+tslices>state.size(3))=[];
                            x=permute(state.xyz_x(:,:,max(1,min(state.size(3),state.pointer_vox(3)+tslices))),[1 2 3]);
                            y=permute(state.xyz_y(:,:,max(1,min(state.size(3),state.pointer_vox(3)+tslices))),[1 2 3]);
                            z=permute(state.xyz_z(:,:,max(1,min(state.size(3),state.pointer_vox(3)+kslices*tslices))),[1 2 3]);
                            z1=permute(state.structural(:,:,max(1,min(state.size(3),state.pointer_vox(3)+tslices)),min(size(state.structural,4),state.time)),[1 2 3]);
                            if state.isvol, z2=permute(state.supra(:,:,max(1,min(state.size(3),state.pointer_vox(3)+tslices)),min(size(state.supra,4),state.time)),[1 2 3]); end
                    end
                    if isempty(tslices)
                        set([state.handles.patch(nview)],'visible','off');
                        %set([state.handles.patchcontour(nview)],'visible','off');
                        if state.isvol, set([state.handles.patchcontour1(nview),state.handles.patchcontour2(nview)],'visible','off'); end
                        if ~isempty(state.surf), set([state.handles.patchcontour3(nview),state.handles.patchcontour4(nview)],'visible','off'); end
                        f1=[];
                        continue
                    end
                    z1(z1<state.volthr)=nan;
                    if state.viewoverlay, tactthr=max(-eps,state.actthr);
                    else tactthr=max(abs(state.Vrange))+1;
                    end
                    if isfield(state,'Vrange')&&state.Vrange(1)==0,      tactthr_pos=tactthr; tactthr_neg=inf;
                    elseif isfield(state,'Vrange')&&state.Vrange(2)==0,  tactthr_pos=inf; tactthr_neg=tactthr;
                    else                        tactthr_pos=tactthr; tactthr_neg=tactthr; 
                    end
                    if state.isvol, 
                        z0=convn(convn(z2>tactthr_pos|z2<-tactthr_neg,conn_hanning(7)/4,'same'),conn_hanning(7)'/4,'same');
                        %if state.isstat, z0=convn(convn(z2>tactthr_pos|z2<-tactthr_neg,conn_hanning(7)/4,'same'),conn_hanning(7)'/4,'same');
                        %else z0=convn(convn((z2>tactthr_pos|z2<-tactthr_neg),conn_hanning(7)/4,'same'),conn_hanning(7)'/4,'same')&~convn(z2,[-1;0;1],'same')&~convn(z2,[-1 0 1],'same');
                        %end
                        if isempty(state.surf),
                            [f1,f2,x0,x1,x2,f3]=conn_slice_display_surf2patch(x,y,z,state.expandmultiview,state.blackistransparent,state.handles.axes,nview,find(tslices==0),[1 1 2 2 2 1],z1,z2,~isnan(z1),{(z2>tactthr_pos).*z2,tactthr},{-z2.*(z2<-tactthr_neg),tactthr},z0);
                        else
                            [f1,f2,x0,x1,x2,x3,x4,f3]=conn_slice_display_surf2patch(x,y,z,state.expandmultiview,state.blackistransparent,state.handles.axes,nview,find(tslices==0),[1 1 2 2 2 3 3 1],z1,z2,~isnan(z1),{(z2>tactthr_pos).*z2,tactthr},{-z2.*(z2<-tactthr_neg),tactthr},state.surf(:,1),state.surf(:,end),z0);
                        end
                    else
                        if isempty(state.surf),
                            [f1,x0]=conn_slice_display_surf2patch(x,y,z,state.expandmultiview,state.blackistransparent,state.handles.axes,nview,find(tslices==0),[1 2],z1,~isnan(z1));
                        else
                            [f1,x0,x3,x4]=conn_slice_display_surf2patch(x,y,z,state.expandmultiview,state.blackistransparent,state.handles.axes,nview,find(tslices==0),[1 2 3 3],z1,~isnan(z1),state.surf(:,1),state.surf(:,end));
                        end
                    end
                    if isempty(minmax), minmax=[min(f1.vertices,[],1);max(f1.vertices,[],1)];
                    else minmax=[min(minmax(1,:),min(f1.vertices,[],1));max(minmax(2,:),max(f1.vertices,[],1))];
                    end
                    c1=f1.facevertexcdata;
                    %h=conn_hanning(5);h=h/sum(h); c1=convn(convn(c1,h,'same'),h','same');
                    c0=c1;
                    %c1=(c1-min(c1))/max(eps,max(c1)-min(c1));
                    c1=(c1-state.Srange(1))/max(eps,state.Srange(2)-state.Srange(1));
                    if 0, 
                        disp([size(c1) nnz(c1)]); % placeholder for transparency threshold/masking
                        c1(c1<.25)=nan;
                    end
                    c=repmat(c1,[1 3]);
                    c=max(0,min(1,state.colorbrightness+c.^state.colorcontrast));
                    if state.isvol, 
                        %c0=c;
                        s2=f2.facevertexcdata;
                        maxs2=max(1e-4,max(abs(s2)));
                        mask1=s2>tactthr_pos;
                        mask2=s2<-tactthr_neg;
                        c2a=max(0,min(1,(s2-max(0,state.Vrange(1)))/max(eps,state.Vrange(2)-max(0,state.Vrange(1)))));
                        c2b=max(0,min(1,(-s2+min(0,state.Vrange(2)))/max(eps,min(0,state.Vrange(2))-state.Vrange(1))));
                        %c2=.1+.9*abs(s2)/maxs2;
                        if state.blackistransparent, fb=f3.facevertexcdata; %.^4;
                        else fb=ones(size(f3.facevertexcdata)); 
                        end
                        alphamix=1-state.transparency*fb;
                        %fb=zeros(size(s2));
                        %fb(mask1)=max(0,(s2(mask1)-tactthr)/max(eps,max(s2)-tactthr)).^.5;
                        %fb(mask2)=max(0,(-s2(mask2)-tactthr)/max(eps,max(-s2)-tactthr)).^.5;
                        %alphamix=max(0,1-2*state.transparency)+2*min(state.transparency,(1-state.transparency))*(1-fb); 
                        %%alphamix=max(0,1-state.transparency);
                        c(mask1,1)=alphamix(mask1).*max(0,c(mask1,1))+(1-alphamix(mask1)).*cmap(c2a(mask1,1),1);
                        c(mask1,2)=alphamix(mask1).*max(0,c(mask1,2))+(1-alphamix(mask1)).*cmap(c2a(mask1,1),2);
                        c(mask1,3)=alphamix(mask1).*max(0,c(mask1,3))+(1-alphamix(mask1)).*cmap(c2a(mask1,1),3);
                        c(mask2,3)=alphamix(mask2).*max(0,c(mask2,3))+(1-alphamix(mask2)).*cmap(c2b(mask2),1);
                        c(mask2,1)=alphamix(mask2).*max(0,c(mask2,1))+(1-alphamix(mask2)).*cmap(c2b(mask2),2);
                        c(mask2,2)=alphamix(mask2).*max(0,c(mask2,2))+(1-alphamix(mask2)).*cmap(c2b(mask2),3);
                        %c=conn_bsxfun(@times,abs(s2)<=tactthr,c0)+conn_bsxfun(@times,abs(s2)>tactthr,c);
                        if state.viewoverlay&&state.contourtransparency
                            set(state.handles.patchcontour1(nview),'xdata',x1(1,:),'ydata',x1(2,:),'zdata',x1(3,:),'color',cmap0(1,:),'visible','on');
                            set(state.handles.patchcontour2(nview),'xdata',x2(1,:),'ydata',x2(2,:),'zdata',x2(3,:),'color',cmap0(1,[2 3 1]),'visible','on');
                        else set([state.handles.patchcontour1(nview),state.handles.patchcontour2(nview)],'visible','off'); 
                        end
                    else fb=1;
                    end
                    if ~isempty(state.surf)
                        if state.viewoverlay&&state.freesurfertransparency
                            set(state.handles.patchcontour3(nview),'xdata',x3(1,:),'ydata',x3(2,:),'zdata',x3(3,:),'color',[1 0 0],'visible','on'); %cmap0(1,:)
                            set(state.handles.patchcontour4(nview),'xdata',x4(1,:),'ydata',x4(2,:),'zdata',x4(3,:),'color',[0 0 1],'visible','on'); %cmap0(1,[2 3 1])
                        else set([state.handles.patchcontour3(nview),state.handles.patchcontour4(nview)],'visible','off'); 
                        end
                    end
                    %c(~all(c>0,2),:)=nan;
                    set(state.handles.patch(nview),'faces',f1.faces,'vertices',f1.vertices,'facevertexcdata',c,'facecolor','flat','edgecolor','none','FaceLighting', 'gouraud','visible','on');
                    if state.blackistransparent, fa=min(1,max(0,10*(c0-state.volthr)/max(eps,state.Srange(2)-state.Srange(1)))).^2;
                    else fa=double(~isnan(c0));
                    end
                    if 1, set(state.handles.patch(nview),'facevertexalpha',state.slice_transparency*fa,'facealpha','flat','AlphaDataMapping','none'); 
                    else set(state.handles.patch(nview),'facevertexalpha',[],'facealpha',state.slice_transparency); 
                    end
                    %set(state.handles.patchcontour(nview),'xdata',x0(1,:),'ydata',x0(2,:),'zdata',x0(3,:),'color',ccolor,'visible','on');
                else
                    set([state.handles.patch(nview)],'visible','off');
                    %set([state.handles.patchcontour(nview)],'visible','off');
                    if state.isvol, set([state.handles.patchcontour1(nview),state.handles.patchcontour2(nview)],'visible','off'); end
                    if ~isempty(state.surf), set([state.handles.patchcontour3(nview),state.handles.patchcontour4(nview)],'visible','off'); end
                end
            end
            set([state.handles.act1 state.handles.act2],'visible','off');
            if state.view(4), set([state.handles.line1 state.handles.line2 state.handles.line3],'visible','on');
            else set([state.handles.line1 state.handles.line2 state.handles.line3],'visible','off');
            end
            if numel(state.view)>4&&state.view(5), 
                if ~isempty(state.handles.act1), set(state.handles.act1(min(state.time,numel(state.handles.act1))),'visible','on'); end
                if ~isempty(state.handles.act2), set(state.handles.act2(min(state.time,numel(state.handles.act2))),'visible','on'); end
            end
            
            if isempty(state.cameraview), state.cameraview=get(gca,'cameraposition'); state.cameraview=state.cameraview(:)'; 
            else view(state.handles.axes,state.cameraview); 
            end
            set(findobj(state.handles.hfig,'type','light'),'position',state.cameraview);
            set([state.handles.line1 state.handles.line2 state.handles.line3],'xdata',[],'ydata',[],'zdata',[]);
            try, set(state.handles.line1,'xdata',state.xyz_x(state.pointer_vox(1),state.pointer_vox(2),:),'ydata',state.xyz_y(state.pointer_vox(1),state.pointer_vox(2),:),'zdata',state.xyz_z(state.pointer_vox(1),state.pointer_vox(2),:)); end
            try, set(state.handles.line2,'xdata',state.xyz_x(state.pointer_vox(1),:,state.pointer_vox(3)),'ydata',state.xyz_y(state.pointer_vox(1),:,state.pointer_vox(3)),'zdata',state.xyz_z(state.pointer_vox(1),:,state.pointer_vox(3))); end
            try, set(state.handles.line3,'xdata',state.xyz_x(:,state.pointer_vox(2),state.pointer_vox(3)),'ydata',state.xyz_y(:,state.pointer_vox(2),state.pointer_vox(3)),'zdata',state.xyz_z(:,state.pointer_vox(2),state.pointer_vox(3))); end
            if state.isstat, 
                set([state.handles.text1 state.handles.text2],'string','','visible','off');
                try, 
                    if get(state.handles.gui,'value')~=1
                        if isequal(state.stats,'T'), set(state.handles.text1,'string',sprintf('Voxel-level: (%d,%d,%d)  %s(%s) = %.2f  p = %.6f (two-sided p = %.6f)',round(state.pointer_mm(1)),round(state.pointer_mm(2)),round(state.pointer_mm(3)), state.stats,mat2str(round(1e4*state.dof(end))/1e4),state.T(state.pointer_vox(1),state.pointer_vox(2),state.pointer_vox(3),min(state.time,size(state.T,4))),state.p(state.pointer_vox(1),state.pointer_vox(2),state.pointer_vox(3),min(state.time,size(state.p,4))),2*min(state.p(state.pointer_vox(1),state.pointer_vox(2),state.pointer_vox(3),min(state.time,size(state.p,4))),1-state.p(state.pointer_vox(1),state.pointer_vox(2),state.pointer_vox(3),min(state.time,size(state.p,4))))),'visible','on');
                        elseif ~isempty(state.p) set(state.handles.text1,'string',sprintf('Voxel-level: (%d,%d,%d)  %s(%s) = %.2f  p = %.6f',round(state.pointer_mm(1)),round(state.pointer_mm(2)),round(state.pointer_mm(3)), state.stats,mat2str(round(1e4*state.dof)/1e4),state.T(state.pointer_vox(1),state.pointer_vox(2),state.pointer_vox(3),min(state.time,size(state.T,4))),state.p(state.pointer_vox(1),state.pointer_vox(2),state.pointer_vox(3),min(state.time,size(state.p,4)))),'visible','on');
                        else set(state.handles.text1,'string',sprintf('Voxel-level: (%d,%d,%d)  %s(%s) = %.2f',round(state.pointer_mm(1)),round(state.pointer_mm(2)),round(state.pointer_mm(3)), state.stats,mat2str(round(1e4*state.dof)/1e4),state.T(state.pointer_vox(1),state.pointer_vox(2),state.pointer_vox(3),min(state.time,size(state.T,4)))),'visible','on');
                        end
                    end
                end
                if ~isempty(state.clusters)&&get(state.handles.gui,'value')~=1
                    supra=state.supra(:,:,:,min(state.time,size(state.supra,4)))~=0;
                    d=inf(state.size(1:3));
                    d(supra)=(state.xyz_x(supra)-state.pointer_mm(1)).^2+(state.xyz_y(supra)-state.pointer_mm(2)).^2+(state.xyz_z(supra)-state.pointer_mm(3)).^2;
                    [mind,idxd]=min(d(:));
                    for n=1:numel(state.clusters.idx)
                        if any(d(state.clusters.idx{n})==mind),
                            set(state.handles.text2,'string',sprintf('Closest-cluster: %s',state.clusters.stats{n}),'visible','on');
                            break;
                        end
                    end
                end
            else
                set([state.handles.text1],'string','','visible','off');
                if get(state.handles.gui,'value')~=1, 
                    try, 
                        val=state.T(state.pointer_vox(1),state.pointer_vox(2),state.pointer_vox(3),min(state.time,size(state.T,4)));
                        lname='';
                        if ~isempty(state.refnames), 
                            if state.reftype==1&&round(val)>0, lname=state.refnames{find(state.refidx==round(val),1)};
                            elseif state.reftype==2, lname=state.refnames{find(state.refidx==state.time,1)};
                            end
                        end
                        set(state.handles.text1,'string',sprintf('Value: (%d,%d,%d)  = %s %s',round(state.pointer_mm(1)),round(state.pointer_mm(2)),round(state.pointer_mm(3)), mat2str(1e-6*round(1e6*val)), lname),'visible','on'); 
                    end
                end
            end
            if get(state.handles.mode,'value')==1,
                set(state.handles.rotate3d,'enable','off');
                set(state.handles.patch,'buttondownfcn',{@conn_slice_display_refresh,'buttondown'});
                set([state.handles.line1 state.handles.line2 state.handles.line3 state.handles.act1 state.handles.act2],'buttondownfcn',{@conn_slice_display_refresh,'buttondown'});
            else
                set(state.handles.rotate3d,'enable','on');
                set(state.handles.patch,'buttondownfcn',[]);
                set([state.handles.line1 state.handles.line2 state.handles.line3 state.handles.act1 state.handles.act2],'buttondownfcn',[]);
            end
            if sum(state.view(1:3))==1&&all(state.nslices==1), 
                npointer=state.cameraviewdirs(find(state.view(1:3),1));
                d=max(0,min(1, (state.pointer_mm(npointer)-state.xyz_range(npointer,1))/(state.xyz_range(npointer,2)-state.xyz_range(npointer,1)) ));
                set(state.handles.slider,'value',d,'sliderstep',[1/max(1,state.xyz_range(npointer,2)-state.xyz_range(npointer,1)), 10/max(1,state.xyz_range(npointer,2)-state.xyz_range(npointer,1))]);
                if strcmp(get(state.handles.pointer_mm(1),'visible'),'on'), set(state.handles.slider,'visible','on'); else set(state.handles.slider,'visible','off'); end
            else set(state.handles.slider,'visible','off');
            end
            delete(state.handles.slidetext(ishandle(state.handles.slidetext)));
            state.handles.slidetext=[];
            if state.viewslicetitle&&sum(state.view(1:3))==1,%&&size(x,3)>1, 
                if isfield(f1,'titletext')
                    for nslice=1:numel(f1.titletext)
                        state.handles.slidetext(nslice)=text(f1.titlepos(1,nslice),f1.titlepos(2,nslice),f1.titlepos(3,nslice),f1.titletext{nslice},'color',mod(state.background+.5,1),'horizontalalignment','center','verticalalignment','middle','parent',state.handles.axes);
                    end
                end
            end
            if ~isempty(minmax), set(state.handles.axes,'xlim',minmax(:,1)'+[-.9 .9],'ylim',minmax(:,2)'+[-.9 .9],'zlim',minmax(:,3)'+[-.9 .9]); end
        end
        if ~isempty(state.handles.colorbar)&&redrawnowcolorbar
            if state.Vrange(1)>=0,
                set(state.handles.colorbar(2),'cdata',max(0,min(1, ind2rgb(round(linspace(1,size(state.cmap,1),128)'),state.cmap))));
            elseif state.Vrange(2)<=0,
                set(state.handles.colorbar(2),'cdata',max(0,min(1, ind2rgb(round(linspace(size(state.cmap,1),1,128)'),state.cmap(:,[2 3 1])))));
            else
                set(state.handles.colorbar(2),'cdata',cat(1,max(0,min(1, ind2rgb(round(linspace(size(state.cmap,1),1,64)'),state.cmap(:,[2 3 1])))),max(0,min(1, ind2rgb(round(linspace(1,size(state.cmap,1),64)'),state.cmap)))));
            end
            set(state.handles.colorbar(3),'string',num2str(state.Vrange(1),'%.2f'),'color',1-round(mean(state.background))*[1 1 1]);
            set(state.handles.colorbar(4),'string','');
            set(state.handles.colorbar(5),'string',num2str(state.Vrange(2),'%.2f'),'color',1-round(mean(state.background))*[1 1 1]);
        end
        if redrawnow>1, drawnow; end
    end
end

function varargout=conn_slice_display_surf2patch(x,y,z,expand,smoothc,axesh,nview,refslice,opts,varargin)
if isempty(refslice), refslice=1; end
minl=[0,8]; % remove contours with length (pixels) below minl
x1=[1.5*x(:,1,refslice)-.5*x(:,2,refslice) .5*x(:,1:end-1,refslice)+.5*x(:,2:end,refslice) 1.5*x(:,end,refslice)-.5*x(:,end-1,refslice)];
x1=[1.5*x1(1,:)-.5*x1(2,:);.5*x1(1:end-1,:)+.5*x1(2:end,:);1.5*x1(end,:)-.5*x1(end-1,:)];
y1=[1.5*y(:,1,refslice)-.5*y(:,2,refslice) .5*y(:,1:end-1,refslice)+.5*y(:,2:end,refslice) 1.5*y(:,end,refslice)-.5*y(:,end-1,refslice)];
y1=[1.5*y1(1,:)-.5*y1(2,:);.5*y1(1:end-1,:)+.5*y1(2:end,:);1.5*y1(end,:)-.5*y1(end-1,:)];
z1=[1.5*z(:,1,refslice)-.5*z(:,2,refslice) .5*z(:,1:end-1,refslice)+.5*z(:,2:end,refslice) 1.5*z(:,end,refslice)-.5*z(:,end-1,refslice)];
z1=[1.5*z1(1,:)-.5*z1(2,:);.5*z1(1:end-1,:)+.5*z1(2:end,:);1.5*z1(end,:)-.5*z1(end-1,:)];
vertices=[x1(:) y1(:) z1(:)];
mvertices=mean(vertices,1);
%mvertices=[mean(mean(x,1),2) mean(mean(y,1),2) mean(mean(z,1),2)];
%svertices=std(mvertices,0,1);
%Mvertices=max(vertices,[],1);
%[nill,Ivertices]=min(svertices);
a=reshape(1:numel(x1),size(x1));
faces=reshape(cat(3, a(1:end-1,1:end-1), a(2:end,1:end-1),a(2:end,2:end),a(1:end-1,2:end)),[],4);

[d1,d2,d3]=deal(zeros(1,3));
x0={x,y,z};
mx0=zeros(3,size(x,3));sx0=mx0;
if expand, zoomout=1.1;
else zoomout=1;
end
for n=1:3,
    d1(n)=zoomout*(size(x0{n},1)+1)*mean(mean(diff(x0{n}(:,:,1),1,1),1),2);
    d2(n)=zoomout*(size(x0{n},2)+1)*mean(mean(diff(x0{n}(:,:,1),1,2),1),2);
    d3(n)=mean(mean(x0{n}(:,:,min(size(x0{n},3),2))-x0{n}(:,:,1),1),2);
    tx=reshape(x0{n},[],size(x0{n},3));
    mx0(n,:)=mean(tx,1);
    sx0(n,:)=std(tx,0,1);
end
[nill,ix0]=min(sum(sx0,2));
switch(nview) % d1: left-to-right in plot; d2: top-to-bottom in plot
    case 1, d1=-abs(d1);d2=-abs(d2);
    case 2, d1=abs(d1); d2=-abs(d2); 
    case 3, d1=abs(d1); d2=-abs(d2); 
end
d=[norm(d1) norm(d2) norm(d3) zeros(1,min(2,numel(varargin)))];

units=get(axesh,'units');
set(axesh,'units','points');
emphx=get(axesh,'position');
set(axesh,'units',units);
emphx=emphx(4)/emphx(3);
[n1,n2]=ndgrid(1:size(x,3),1:size(x,3)); n1=n1(:); n2=n2(:);
[nill,idx]=max(min(d(2)./n1(:), emphx*d(1)./n2(:))-1e10*(n1(:).*n2(:)<size(x,3)));
n1=n1(idx); n2=n2(idx);
[i1,i2]=ind2sub([n1 n2],1:size(x,3));
if ~expand, i1(:)=1; i2(:)=1; end

done=false;
for n=find(opts==1),
    if ~done
        done=true;
        varargout{n}=struct('vertices',zeros(size(vertices,1)*numel(i1),3),'faces',zeros(size(faces,1)*numel(i1),4),'facevertexcdata',zeros(size(faces,1)*numel(i1),1),'titletext',{{}},'titlepos',zeros(3,0));
        for m=1:numel(i1)
            %newvertices=vertices+repmat(d1*(i1(m)-ceil(n1/2))+d2*(i2(m)-ceil(n2/2))+d3*(m-1),size(vertices,1),1);
            newvertices=vertices+repmat(d1*(i1(m)-i1(refslice))+d2*(i2(m)-i2(refslice))+d3*(m-refslice),size(vertices,1),1);
            newfaces=(m-1)*size(vertices,1)+faces;
            newc=reshape(varargin{n}(:,:,m),[],1);
            varargout{n}.faces((m-1)*size(faces,1)+(1:size(faces,1)),:)=newfaces;
            varargout{n}.vertices((m-1)*size(vertices,1)+(1:size(vertices,1)),:)=newvertices;
            varargout{n}.facevertexcdata((m-1)*size(faces,1)+(1:size(faces,1)),:)=newc;
            varargout{n}.titlepos(:,m)=(mvertices+d1*(i1(m)-i1(refslice)-0)+d2*(i2(m)-i2(refslice)-.49)+d3*(m-refslice))';
            varargout{n}.titletext{m}=sprintf('%c = %.0f',char('x'+ix0-1),mx0(ix0,m));
        end
    else
        varargout{n}=varargout{1};
        for m=1:numel(i1)
            newc=reshape(varargin{n}(:,:,m),[],1);
            varargout{n}.facevertexcdata((m-1)*size(faces,1)+(1:size(faces,1)),:)=newc;
        end
    end
end
if any(opts==2)
    if size(x,1)==1, x=cat(1,x,x); y=cat(1,y,y); z=cat(1,z,z); end
    if size(x,2)==1, x=cat(2,x,x); y=cat(2,y,y); z=cat(2,z,z); end
    if size(x,3)==1, x=cat(3,x,x); y=cat(3,y,y); z=cat(3,z,z); end
    for n=find(opts==2)
        c1=zeros(3,0); 
        %tmin=min(varargin{n}(:));tmax=max(varargin{n}(:));
        if 1,%~isnan(tmin)&&~isnan(tmax)&&tmin~=tmax
            if iscell(varargin{n}), data=varargin{n}{1}; thr=varargin{n}{2};
            else data=varargin{n}; thr=.5;
            end
            for m=1:size(data,3)
                tdata=zeros(size(data,1)+2,size(data,2)+2);
                tdata(2:end-1,2:end-1)=double(data(:,:,m));
                %tdata=double(data(:,:,m));
                tdata(isnan(tdata))=0;
                ct=contourc(0:size(data,2)+1,0:size(data,1)+1,tdata,thr*[1 1]); %tmin+(tmax-tmin)*[.10 .10]); %[.5 .5]);
                if ~isempty(ct), c1=[c1 [ct;m+zeros(1,size(ct,2))]]; end
            end
        end
        x1=zeros(3,0);
        try, x1=cat(1,interpn(x,c1(2,:),c1(1,:),refslice+0*c1(3,:)),interpn(y,c1(2,:),c1(1,:),refslice+0*c1(3,:)),interpn(z,c1(2,:),c1(1,:),refslice+0*c1(3,:))); end
        if ~isempty(x1), 
            m=c1(3,:);
            x1=x1+(d1'*(i1(m)-i1(refslice))+d2'*(i2(m)-i2(refslice))+d3'*(m-refslice));
            h=conn_hanning(5)'.^2/2.25;%h=conn_hanning(5)'/3;
            ci=1;while ci<size(c1,2), cj=ci+c1(2,ci)+1; if cj-ci<=minl(1+(smoothc~=0)), x1(:,ci:cj-1)=nan; else x1(:,ci)=nan; if smoothc, x1(:,ci+1:cj-2)=convn(x1(:,ci+1+mod(0:cj-ci-3+4,cj-ci-2)),h,'valid'); x1(:,cj-1)=x1(:,ci+1); end; end; ci=cj; end; 
        end
        varargout{n}=x1;
    end
end
if any(opts==3)
    dnull=null([d1;d2]);
    wb=dnull'*mx0;
    for n=find(opts==3)
        alldata=varargin{n}; 
        x1=zeros(3,0);
        for ndata=1:numel(alldata)
            data=alldata(ndata);
            wa=data.vertices*dnull;
            for m=1:numel(i1)
                wc=(wa-wb(m));
                [ws,wi]=sort(wc(data.faces),2);
                wd=reshape(find(ws(:,1).*ws(:,end)<-1e-4),[],1);
                if ~isempty(wd)
                    tx1=data.vertices(data.faces(wd+size(data.faces,1)*(wi(wd,1)-1)),:);
                    tx2=data.vertices(data.faces(wd+size(data.faces,1)*(wi(wd,2)-1)),:);
                    tx3=data.vertices(data.faces(wd+size(data.faces,1)*(wi(wd,end)-1)),:);
                    doswitch=find(ws(wd,2).*ws(wd,3)<0); tx0=tx1(doswitch,:); tx1(doswitch,:)=tx3(doswitch,:); tx3(doswitch,:)=tx0;
                    w1=tx1*dnull-wb(m);w2=tx2*dnull-wb(m);w3=tx3*dnull-wb(m);
                    w=abs(w2)./max(eps,abs(w2-w1)); txa=tx1.*repmat(w,1,3) + tx2.*repmat(1-w,1,3);
                    w=abs(w3)./max(eps,abs(w3-w1)); txb=tx1.*repmat(w,1,3) + tx3.*repmat(1-w,1,3);
                    txa=txa-(txa*dnull-wb(m))*dnull';
                    txb=txb-(txb*dnull-wb(m))*dnull';
                    txc=reshape([txa txb nan(size(txa))]',3,[]);
                    x1=[x1 txc+repmat(d1'*(i1(m)-i1(refslice))+d2'*(i2(m)-i2(refslice))+0*d3'*(m-refslice),1,size(txc,2))]; % txc+repmat(d1'*(i1(m)-i1(refslice))+d2'*(i2(m)-i2(refslice))-.01*d3'*(m-refslice),1,size(txc,2))];
                end
            end
        end
        varargout{n}=x1;
    end
end
end

function c=spring(m)
r = (0:m-1)'/max(m-1,1); 
c = [ones(m,1) r 1-r];
end
function c=summer(m)
r = (0:m-1)'/max(m-1,1); 
c = [r .5+r/2 .4*ones(m,1)];
end
function c=winter(m)
r = (0:m-1)'/max(m-1,1); 
c = [zeros(m,1) r .5+(1-r)/2];
end
function c=autumn(m)
r = (0:m-1)'/max(m-1,1);
c = [ones(m,1) r zeros(m,1)];
end

function c=fixedge(c)
k=linspace(1,0,size(c,1))'.^4;
c=repmat(1-k,1,size(c,2)).*c+repmat(k.*mean(c,2),1,size(c,2));
end





    