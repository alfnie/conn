function hfigure=conn_display(SPMfilename,varargin)
% CONN_DISPLAY displays second-level analysis results
%
% BASIC SYNTAX:
%
% hf = conn_display(SPMfilename)
%   SPMfilename : input file name (typically SPM.mat) containing second-level analysis results
%                 either from voxel-based analyses, surface-based analyses, or ROI-to-ROI analyses
%   hf          : output handle to results explorer figure
%
% e.g. conn_display('/results/SPM.mat');
%
% ADVANCED OPTIONS FOR VOXEL- AND SURFACE- LEVEL RESULTS FIGURES (applied when creating a new results figure)
%
% hf = conn_display(SPMfilename, ncon, style)
%   SPMfilename : input file name (typically SPM.mat) containing second-level analysis results
%   ncon        : input contrast number (if multiple contrasts are defined in 2nd-level results SPM.mat file)
%   style       : inpput voxel/cluster threshold settings (see Nieto-Castanon, 2020 for details about these methods; www.conn-toolbox.org/fmri-methods)
%                   [style=1] for default settings #1 : parametric statistics (Gaussian Random Field theory, Worsley et al. 1996) (default for volume-based results) 
%                   [style=2] for default settings #2 : nonparametric statistics (permutation/randomization analyses, Bullmore et al. 1999) (default for surface-based results)
%                   [style=3] for default settings #3 : nonparametric statistics (Threshold Free Cluster Enhancement, Smith and Nichols, 2007)
%   hf          : output results explorer window handle (which can be used in advanced options described below)
%
% e.g. conn_display('/results/SPM.mat',1,3);
%
% ADVANCED OPTIONS FOR VOXEL- AND SURFACE- LEVEL RESULTS FIGURES (applied to an existing results figure)
% 
% conn_display(hf, 'surface_print', filename)       prints image with current results projected to the cortical-surface
% conn_display(hf, 'surface_view', filename)        displays image with current results projected to the cortical-surface
% conn_display(hf, 'volume_print', filename)        prints image with current results on 3d brain
% conn_display(hf, 'volume_view', filename)         displays image with current results on 3d brain
% conn_display(hf, 'slice_print', filename)         prints image with current results on individual slices
% conn_display(hf, 'slice_view', filename)          displays image with current results on individual slices
% conn_display(hf, 'glass_print', filename)         prints image with current results on 3d glass-brain
% conn_display(hf, 'glass_view', filename)          displays image with current results on 3d glass-brain
% conn_display(hf, 'export_mask', filename)         exports NIFTI mask file with current results
% conn_display(hf, 'ref_atlas', filename)           uses alternative reference atlas for anatomical descriptions of each cluster
% conn_display(hf, 'close')                         closes results explorer window
% conn_display(hf, 'fwec.option', style)            modifies default threshold settings (see 'style' valid values above)
% conn_display(hf, 'fwec.voxellevel.value', value)  modifies voxel-level threshold value
% conn_display(hf, 'fwec.voxellevel.type', type)    modifies voxel-level threshold type
%                                                   [type=1] p-uncorrected
%                                                   [type=2] p-FDR corrected
%                                                   [type=3] p-FDR corrected (TFCE)
%                                                   [type=4] p-FWE corrected (TFCE)
%                                                   [type=5] F/T/X stat
% conn_display(hf, 'fwec.voxellevel.side', type)    modifies voxel-level threshold directionality
%                                                   [side=1] one-sided (positive)
%                                                   [side=2] one-sided (negative)
%                                                   [side=3] two-sided
% conn_display(hf, 'fwec.clusterlevel.value', value) modifies cluster-level threshold value
% conn_display(hf, 'fwec.clusterlevel.type', type)    modifies cluster-level threshold type 
%                                                   [type=1] cluster-size
%                                                   [type=2] cluster-size p-FWE corrected
%                                                   [type=3] cluster-size p-FDR corrected
%                                                   [type=4] cluster-size p-uncorrected
%                                                   [type=5] peak-voxel p-FWE corrected
%                                                   [type=6] peak-voxel p-uncorrected
%                                                   [type=7] cluster-mass
%                                                   [type=8] cluster-mass p-FWE corrected
%                                                   [type=9] cluster-mass p-FDR corrected
%                                                   [type=10] cluster-mass p-uncorrected
% conn_display(hf, 'fwec.clusterlevel.parametric', type) modifies cluster-level parametric statistics option
%                                                   [type=1] parametric statistics (Random Field Theory)
%                                                   [type=2] non-parametric statistics (Randomization/Permutation tests)
%
% ADVANCED OPTIONS FOR ROI-to-ROI RESULTS FIGURES (applied when creating a new results figure)
%
% hf = conn_display(SPMfilename, ncon, style)
%   SPMfilename : input file name (either SPM.mat or ROI.mat file) containing second-level analysis results
%   ncon        : contrast number (if multiple contrasts are defined in 2nd-level results SPM.mat file)
%   style       : voxel/cluster threshold settings (see Nieto-Castanon, 2020 for details about these methods; www.conn-toolbox.org/fmri-methods)
%                   [style=1] for default settings #1 : parametric multivariate statistics (Functional Network Connectivity, Jafri et al., 2008)
%                   [style=2] for default settings #2 : non-parametric statistics (Spatial Pairwise Clustering, Zalesky et al., 2012)
%                   [style=3] for default settings #3 : non-parametric statistics (Threshold Free Cluster Enhancement, Smith and Nichols, 2007)
%                   [style=4] for alternative settings for connection-based inferences : parametric univariate statistics (connection-level inferences, FDR corrected, Benjamini & Hochberg, 1995)
%                   [style=5] for alternative settings for ROI-level inferences: parametric multivariate statistics (ROI-level inferences, FDR corrected, Benjamini & Hochberg, 1995)
%                   [style=6] for alternative settings for network-level inferences: non-parametric statistics (network-level inferences, Network Based Statistics, Zalesky et al., 2010)
%   hf          : output results explorer window handle (which can be used in advanced options described below)
%
% e.g. conn_display('/results/SPM.mat',1,3);
%
% ADVANCED OPTIONS FOR ROI-to-ROI RESULTS FIGURES (applied to an existing results figure)
%
% conn_display(hf, 'ring_print', filename)          prints image with current results on ring/circle view
% conn_display(hf, 'ring_view', filename)           display image with current results on ring/circle view 
% conn_display(hf, 'glass_print', filename)         prints image with current results on 3d glass-brain
% conn_display(hf, 'glass_view', filename)          displays image with current results on 3d glass-brain
% conn_display(hf, 'matrix_print', filename)        prints image with current results on ROI-to-ROI matrix view
% conn_display(hf, 'matrix_view', filename)         displays image with current results on ROI-to-ROI matrix view
% conn_display(hf, 'export_mask', filename)         exports NIFTI mask file with current results
% conn_display(hf, 'menubar')                       displays/hides menubar with advanced display options
% conn_display(hf, 'close')                         closes results explorer window
% conn_display(hf, 'roi.select' [,ROInames])        modifies selection of ROIs
% conn_display(hf, 'roi.order',option [,filename])  modifies ROI sorting procedure
%                                                   [option=1] Use hierarchical clustering method (default)
%                                                   [option=2] Use CONN atlas apriori order/groups (atlas.nii ROIs only)
%                                                   [option=3] Use CONN networks apriori order/groups (networks.nii ROIs only)
%                                                   [option=4] Load ROI order/groups from input file filename (see 'doc conn_roiclusters' for details)
%                                                   [option=5] Load ROI order/groups from clipboard
%                                                   [option=6] Save ROI order/groups to output file filename
%                                                   [option=7] Save ROI order/groups to clipboard
%                                                   [option=8] Manually define ROI order/groups
% conn_display(hf, 'fwec.option', style)            modifies default threshold settings (see 'style' valid values above)
% conn_display(hf, 'fwec.connectionlevel.value', value) modifies voxel-level threshold value 
% conn_display(hf, 'fwec.connectionlevel.type', type) modifies voxel-level threshold type
%                                                   [type=1] p-uncorrected
%                                                   [type=2] p-FDR corrected
%                                                   [type=3] p-FDR corrected (TFCE)
%                                                   [type=4] p-FWE corrected (TFCE)
%                                                   [type=5] F/T/X stat
% conn_display(hf, 'fwec.connectionlevel.side', type) modifies voxel-level threshold directionality
%                                                   [side=1] one-sided (positive)
%                                                   [side=2] one-sided (negative)
%                                                   [side=3] two-sided
% conn_display(hf, 'fwec.clusterlevel.value', value) modifies cluster-level threshold value
% conn_display(hf, 'fwec.clusterlevel.type', type)    modifies cluster-level threshold type 
%                                                   [type=1] network-level p-uncorrected (NBS mass/intensity)
%                                                   [type=2] network-level p-FDR corrected (NBS mass/intensity)
%                                                   [type=3] network-level p-FWE corrected (NBS mass/intensity)
%                                                   [type=4] cluster-level p-uncorrected (MVPA omnibus test)
%                                                   [type=5] cluster-level p-FDR corrected (MVPA omnibus test)
%                                                   [type=6] cluster-level p-uncorrected (SPC mass/intensity)
%                                                   [type=7] cluster-level p-FDR corrected (SPC mass/intensity)
%                                                   [type=8] cluster-level p-FWE corrected (SPC mass/intensity)
%                                                   [type=9] ROI-level p-uncorrected (MVPA omnibus test)
%                                                   [type=10] ROI-level p-FDR corrected (MVPA omnibus test)
%                                                   [type=11] ROI-level p-uncorrected (ROI mass/intensity)
%                                                   [type=12] ROI-level p-FDR corrected (ROI mass/intensity)
%                                                   [type=13] ROI-level p-FWE corrected (ROI mass/intensity)
%                                                   [type=14] none
%
%
% see also CONN_MODULE
%

% 03/09 alfnie@gmail.com
%

% note: for voxel-based results
% conn_display(SPMfilename, ncon, {vox_height,vox_type,clu_height,clu_type},side,parametric) uses user-defined thresholding options
% note: for ROI-to-ROI results
% conn_display(SPMfilename, ncon, {con_height,con_type,clu_height,clu_type}) uses user-defined thresholding options
%

global CONN_x CONN_gui;
if isempty(CONN_gui)||~isfield(CONN_gui,'font_offset'), conn_font_init; end

hfigure=[];
fulld=true;
if nargin<1 || isempty(SPMfilename),
    [filename,filepath]=conn_fileutils('uigetfile','SPM.mat');
    if ~ischar(filename), return; end
    SPMfilename=fullfile(filepath,filename);
elseif ishandle(SPMfilename) % conn_display(hfig,opts)
    str=get(SPMfilename,'tag');
    if strcmp(str,'conn_displayroi'), 
        conn_displayroi(SPMfilename,[],varargin{:});
    else
        conn_vproject(SPMfilename,[],varargin{:});
    end
    return;
else
    if conn_fileutils('isdir',SPMfilename), 
        if conn_existfile(fullfile(SPMfilename,'SPM.mat')), SPMfilename=fullfile(SPMfilename,'SPM.mat'); 
        elseif conn_existfile(fullfile(SPMfilename,'ROI.mat')), SPMfilename=fullfile(SPMfilename,'ROI.mat'); 
        else SPMfilename=fullfile(SPMfilename,'SPM.mat'); 
        end
    end
    [filepath,filename,fileext]=fileparts(SPMfilename);
    if isempty(filepath), filepath=pwd;end
    if isempty(filename), SPMfilename=conn_fullfile(filepath,'SPM.mat');
    else SPMfilename=conn_fullfile(filepath,[filename,fileext]);
    end
end
ncon=1;
DOCROP=true;
THR=1;%{.001,1,.05,3};
side=1;
parametric=[];
if numel(varargin)>=1&&~isempty(varargin{1}), ncon=varargin{1}; end
if numel(varargin)>=2&&~isempty(varargin{2}), THR=varargin{2}; end
if numel(varargin)>=3&&~isempty(varargin{3}), side=varargin{3}; end
if numel(varargin)>=4&&~isempty(varargin{4}), parametric=varargin{4}; end

forceusespmresults=false;
try, if nargin>2&&any(strcmp(varargin,'forceusespmresults')), forceusespmresults=true; end; end
[filepath,filename,fileext]=fileparts(SPMfilename);

cwd=pwd;

if ~isempty(filepath), conn_fileutils('cd',filepath); end
if strcmp([filename,fileext],'ROI.mat'), % ROI-to-ROI analyses
    hfigure=conn_displayroi('initfile',fullfile(filepath,[filename,fileext]),varargin{:});
    cd(cwd);
    return;
end
hm=conn_msgbox({sprintf('Loading %s',SPMfilename),'please wait...'},''); 
if fulld, SPM=struct; conn_loadmatfile(fullfile(filepath,[filename,fileext]),'SPM'); end
if ~isfield(SPM.xX,'isSurface'), SPM.xX.isSurface=false; end
issurface=SPM.xX.isSurface;
if issurface&~iscell(THR), THR=max(2,THR); end % skips parametric options
if isfield(SPM.xX,'isMtx')&&SPM.xX.isMtx % ROI-to-ROI analyses
    close(hm);
    hfigure=conn_displayroi('initspm',fullfile(filepath,[filename,fileext]),varargin{:});
    cd(cwd);
    return    
% elseif isfield(SPM.xX,'isSurface')&&SPM.xX.isSurface
%     close(hm);
%     conn_surf_results(fullfile(filepath,[filename,fileext]));
else % voxel-based
    voxeltovoxel=0;
    vol=SPM.xY.VY(1);
    % if isfield(CONN_x,'Setup')&&isfield(CONN_x.Setup,'steps')&&numel(CONN_x.Setup.steps)>2, voxeltovoxel=CONN_x.Setup.steps([3]);
    % else voxeltovoxel=0;
    % end
    [x,y,z]=ndgrid(1:vol.dim(1),1:vol.dim(2),1:vol.dim(3));
    xyz=vol.mat*[x(:),y(:),z(:),ones(numel(x),1)]';
    if 0,
        filename=CONN_x.Setup.explicitmask{1};
        strfile=spm_vol(filename);
        a=mean(reshape(spm_get_data(strfile,pinv(strfile.mat)*xyz),vol.dim(1),vol.dim(2),vol.dim(3),[]),4);
    elseif ~issurface&&conn_existfile(fullfile(fileparts(which(mfilename)),'utils','surf','mask.volume.brainmask.nii')),
        filename=fullfile(fileparts(which(mfilename)),'utils','surf','mask.volume.brainmask.nii');
        strfile=spm_vol(filename);
        a=reshape(spm_get_data(strfile,pinv(strfile.mat)*xyz),vol.dim(1:3));
    elseif issurface&&conn_existfile(fullfile(fileparts(which(mfilename)),'utils','surf','mask.surface.brainmask.nii')),
        filename=fullfile(fileparts(which(mfilename)),'utils','surf','mask.surface.brainmask.nii');
        strfile=spm_vol(filename);
        a=reshape(spm_get_data(strfile,pinv(strfile.mat)*xyz),vol.dim(1:3));
    else a=ones(vol.dim(1:3));
    end
    ab0=(a>0).*a;
%     if isfield(CONN_gui,'refs')&&isfield(CONN_gui.refs,'canonical')&&isfield(CONN_gui.refs.canonical,'filename')&&~isempty(CONN_gui.refs.canonical.filename)
%         strfile=spm_vol(CONN_gui.refs.canonical.filename);
%     else
%         strfile=spm_vol(fullfile(fileparts(which('spm')),'canonical','avg305T1.nii'));
%     end
%     b=reshape(spm_get_data(strfile,pinv(strfile.mat)*xyz),vol.dim(1:3));
%     ab0=(a>0).*(a.*b);
    isparam=~isempty(SPM)&&isfield(SPM,'xVol');
    isnonparam=~isempty(SPM)&&(forceusespmresults||(isfield(SPM,'xX_multivariate')&&isfield(SPM.xX_multivariate,'F')));
    [param,nonparam]=deal(struct('p',[],'logp',[],'F',[],'stats',[],'backg',[]));
    for doparam=find([isparam isnonparam])
        if doparam==1 % SPM-based analyses, for parametric stats
            if isfield(SPM,'xCon')&&length(SPM.xCon)>1,
                if nargin<2||isempty(ncon)||isequal(ncon,'?'), ncon=spm_conman(SPM); end
            else ncon=1;
            end
            if ischar(ncon), 
                tcons={SPM.xCon.name};
                tncon=find(ismember(tcons,ncon)); 
                assert(~isempty(tncon), 'unable to find contrast %s; available contrasts are %s',ncon,sprintf('%s, ',tcons{:}));
                ncon=tncon;
            end
            statsname=SPM.xCon(ncon).STAT;
            [nill,tfname,tfext]=fileparts(SPM.xCon(ncon).Vspm.fname);
            tvol=conn_fileutils('spm_vol',fullfile(filepath,[tfname,tfext]));
            %tvol=SPM.xCon(ncon).Vspm;
            if length(tvol.dim)>3,tvol.dim=tvol.dim(1:3); end;
            if ~isfield(tvol,'dt'), tvol.dt=[spm_type('float32') spm_platform('bigend')]; end
            try
                T=conn_fileutils('spm_read_vols',tvol);
                df=[SPM.xCon(ncon).eidf,SPM.xX.erdf];
                R=SPM.xVol.R;
                S=SPM.xVol.S;
                stat_thresh_fwhm=1;
                try, v2r=1/prod(SPM.xVol.FWHM(1:3)); catch, v2r=[]; end
            catch
                isparam=false;
                continue;
            end
            %             dof=regexp(a(1).descrip,'SPM{T_\[([\d\.])*\]','tokens');
        else % non SPM-based analyses, for nonparametric stats
            if ~forceusespmresults
                T=reshape(SPM.xX_multivariate.F,vol.dim);
                statsname=SPM.xX_multivariate.statsname;
                df=SPM.xX_multivariate.dof;
                %ncon=1;
            else
                T=conn_fileutils('spm_read_vols',tvol);
                SPM.xX_multivariate.X=SPM.xX.X;
                SPM.xX_multivariate.C=SPM.xCon(ncon).c';
                SPM.xX_multivariate.M=1;
                SPM.xX_multivariate.statsname=statsname;
                SPM.xX_multivariate.dof=df;
                SPM.xX_multivariate.F=permute(T,[4,5,1,2,3]);
                [nill,tfname,tfext]=fileparts(SPM.xCon(ncon).Vcon.fname);
                tvol=conn_fileutils('spm_vol',fullfile(filepath,[tfname,tfext]));
                SPM.xX_multivariate.h=permute(conn_fileutils('spm_read_vols',tvol),[4,5,1,2,3]);
                %SPM.xX_multivariate.h=permute(conn_fileutils('spm_read_vols',SPM.xCon(ncon).Vcon),[4,5,1,2,3]);
                SPM.xX_multivariate.derivedfromspm=true;
                conn_savematfile(SPMfilename,'SPM','-v7.3');
            end
            R=[];S=[];v2r=[];stat_thresh_fwhm=[];
        end
        p=nan+zeros(size(T));idxvalid=find(~isnan(T));
        if ~isempty(idxvalid)
            if size(df(:,:,:),3)==numel(T)
                switch(statsname),
                    case 'T', p(idxvalid)=1-spm_Tcdf(T(idxvalid),squeeze(df(1,end,idxvalid))); side=3;
                    case 'F', p(idxvalid)=1-spm_Fcdf(T(idxvalid),squeeze(df(1,1,idxvalid)),squeeze(df(1,2,idxvalid)));
                    case 'X', p(idxvalid)=1-spm_Xcdf(T(idxvalid),squeeze(df(1,end,idxvalid)));
                end
            else
                switch(statsname),
                    case 'T', p(idxvalid)=1-spm_Tcdf(T(idxvalid),df(end)); side=3;
                    case 'F', p(idxvalid)=1-spm_Fcdf(T(idxvalid),df(1),df(2));
                    case 'X', p(idxvalid)=1-spm_Xcdf(T(idxvalid),df);
                end
            end
        end
        if ~DOCROP||issurface
            x=[1,size(T,1)];
            y=[1,size(T,2)];
            z=[1,size(T,3)];
        else
            x=find(any(any(ab0,2),3)|any(any(T,2),3)); x=max(1,min(size(T,1), [x(1)-1 x(end)+1])); %crop
            y=find(any(any(ab0,1),3)|any(any(T,1),3)); y=max(1,min(size(T,2), [y(1)-1 y(end)+1]));
            z=find(any(any(ab0,1),2)|any(any(T,1),2)); z=max(1,min(size(T,3), [z(1)-1 z(end)+1]));
        end
        p=p(x(1):x(end),y(1):y(end),z(1):z(end));
        ab=ab0(x(1):x(end),y(1):y(end),z(1):z(end));
        ab(:,:,[1,end],1)=0;ab(:,[1,end],:,1)=0;ab([1,end],:,:,1)=0;
        T=T(x(1):x(end),y(1):y(end),z(1):z(end));
        logp=-log(max(eps,p));%logp=p;logp(logp==0)=nan; logp=-log(logp);
        logp(isnan(T))=nan;
        stats={vol.mat*(eye(4)+[zeros(4,3),[x(1)-1;y(1)-1;z(1)-1;0]]),R,df,S,v2r,statsname,stat_thresh_fwhm};
        if doparam==1, param=struct('p',p,'logp',logp,'F',T,'stats',{stats},'backg',ab);
        else nonparam=struct('p',p,'logp',logp,'F',T,'stats',{stats},'backg',ab);
        end
    end
    h0=get(0,'screensize');
    hfigure=figure('menubar','none','numbertitle','off','color','w','units','pixels','position',[h0(3)-.75*h0(3)+2,h0(4)-.9*h0(4)-48,.75*h0(3)-2,.9*h0(4)]);
    close(hm);
    conn_vproject(param,nonparam,[],[],THR,side,parametric,[],[],[],.50,[],SPMfilename,voxeltovoxel,issurface);
end

cd(cwd);
