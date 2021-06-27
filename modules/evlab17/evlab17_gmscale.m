function varargout=evlab17_gmscale(spmfile,contrastname,gmtype,overwrite)
% EVLAB17_GMSCALE grand-mean normalization of contrast images
%
% evlab17_gmscale(spmfile, contrastname, gmtype)
%   spmfile      : location of first-level SPM.mat file (or alternatively SPM structure)
%   contrastname : existing contrast name(s) (or alternatively index to SPM.xCon array) 
%   gmtype       : type of grand-mean scaling 
%      'SPMgs'      : SPM grand-mean definition based on global signal mask 
%                     SPMgs = average BOLD signal in voxels where BOLD signal is above 1/8th of the average BOLD signal across all voxels (note: this procedure is sensitive to volume size)
%      'MNIbb1'     : SPM grand-mean definition based on global signal mask if the data had been resampled to bounding box x=[-78 78] mm, y=[-112 76] mm and z=[-50 86] mm (79x95x69 voxels)
%                     MNIbb1 = average BOLD signal in voxels where BOLD signal is above 1/8th of the average BOLD signal across all voxels within small bounding box
%      'MNIbb2'     : SPM grand-mean definition based on global signal mask if the data had been resampled to bounding box x=[-90 90] mm, y=[-126 90] mm, and z=[-72 108] (91x109x91 voxels)
%                     MNIbb2 = average BOLD signal in voxels where BOLD signal is above 1/8th of the average BOLD signal across all voxels within large bounding box
%

if nargin<4||isempty(overwrite), overwrite=true; end
opts={'SPMgs','MNIbb1','MNIbb2','SPMam0','SPMam'};
pwd0=pwd;

if nargin<1||isempty(spmfile), 
    [spmfile,spmpath]=uigetfile('SPM.mat','Select SPM analysis file');
    if ~ischar(spmfile)||isempty(spmfile), return; end
    spmfile=fullfile(spmpath,spmfile);
elseif isstruct(spmfile)
    spmpath=spmfile.swd;
else
    [spmpath,spmname]=fileparts(spmfile);
    if isempty(spmpath), spmpath=pwd; end
end
cd(spmpath);
if isstruct(spmfile), SPM=spmfile;
else load(spmpfile,'SPM');
end
if nargin<2||isempty(contrastname)
    ncon=listdlg('name',mfilename,'PromptString','Select contrast to grand-mean normalize','ListString',{SPM.xCon.name},'SelectionMode','multiple','ListSize',[500 200]);
    if isempty(ncon), return; end
    contrastname={SPM.xCon(ncon).name};
end
if ischar(contrastname), contrastname=cellstr(contrastname); 
elseif isnumeric(contrastname), contrastname={SPM.xCon(contrastname).name};
end
[ok,icon]=ismember(contrastname,{SPM.xCon.name});
assert(all(ok),'unable to match contrast name(s) %s',sprintf('%s ',contrastname{~ok}));
if nargin<3||isempty(gmtype)
    nopt=listdlg('name',mfilename,'PromptString','Select grand-mean scaling method','ListString',opts(1:3),'SelectionMode','single','ListSize',[500 200]);
    if isempty(nopt), return; end
    gmtype=opts{nopt};
end
assert(ismember(gmtype,opts),'unable to match grand-mean scaling method %s (valid names are %s)',gmtype,sprintf('%s ',opts{:}));

if ~overwrite
    for ncon=1:numel(contrastname)
        fname=conn_prepend('',SPM.xCon(icon(ncon)).Vcon.fname,['.',gmtype,'.nii']);
        done(ncon)=conn_existfile(fname);
    end
    contrastname=contrastname(~done);
    if isempty(contrastname), fprintf('All output files already exist. Skipping grand-mean scaling\n'); return; end
    [ok,icon]=ismember(contrastname,{SPM.xCon.name});
end

switch(gmtype)
    case {'SPMgs','SPMam0','SPMam'}
        dim=SPM.Vbeta(1).dim;
        mat=SPM.Vbeta(1).mat;
        [x,y,z]=ndgrid(1:dim(1),1:dim(2),1:dim(3));
        xyz=[x(:) y(:) z(:) ones(numel(x),1)]';
    case 'MNIbb1'
        dim=[79 95 69];
        mat=[-2 0 0 80;0 2 0 -114;0 0 2 -52;0 0 0 1];
        [x,y,z]=ndgrid(1:dim(1),1:dim(2),1:dim(3));
        xyz=pinv(SPM.Vbeta(1).mat)*mat*[x(:) y(:) z(:) ones(numel(x),1)]';
    case 'MNIbb2'
        dim=[91 109 91];
        mat=[-2 0 0 92;0 2 0 -128;0 0 2 -74;0 0 0 1];
        [x,y,z]=ndgrid(1:dim(1),1:dim(2),1:dim(3));
        xyz=pinv(SPM.Vbeta(1).mat)*mat*[x(:) y(:) z(:) ones(numel(x),1)]';
end

% computes average BOLD signal
%note: 100/mean(SPM.xGX.rg(SPM.Sess(n).row)) == SPM.xGX.gSF(SPM.Sess(n).row(k)) == SPM.xY.VY(SPM.Sess(n).row(k)).pinfo(1)

fprintf('computing average BOLD signal...');
if strcmp(gmtype,'SPMam')
    % computes BOLD average from design matrix and regressor coefficients
    M=0; % average BOLD signal from (possibly already grand-mean scaled) regressor images / contrast images
    X=SPM.xX.X;
    c0=mean(X,1);
    for n=1:numel(c0)
        if c0(n)~=0
            try, b=spm_get_data(SPM.Vbeta(n),xyz);
            catch, b=spm_get_data(spm_vol(SPM.Vbeta(n).fname),xyz);
            end
            M=M+c0(n)*reshape(b,dim);
        end
    end
else
    % computes BOLD average from original functional data
    try, SPM.xY.P=char(SPM.xY.Y); end
    VY=spm_vol(SPM.xY.P);
    gs=zeros(1,numel(VY));
    factor=zeros(1,numel(VY));
    for nrun=1:numel(SPM.Sess),
        M=0; % average BOLD signal from raw data (before any SPM-based grand-mean scaling)
        for n=reshape(SPM.Sess(nrun).row,1,[])
            b=spm_get_data(VY(n),xyz);
            mask=isfinite(b);
            gs(n)=mean(b(mask&b>mean(b(mask))/8));
            factor(n)=SPM.xY.VY(n).pinfo(1)/VY(n).pinfo(1);
            M=M+reshape(b,dim);
        end
        gm(nrun)=mean(gs(SPM.Sess(nrun).row));  % global mean
        fm(nrun)=mean(factor(SPM.Sess(nrun).row));  % factor between original raw data and original grand-mean-scaled data (typically equal to 100/gm)
    end
    M=M/numel(VY);
end
fprintf(' Done\n');

% compute grand mean scaling factors
if ismember(gmtype,{'SPMam','SPMam0'})
    mask=isfinite(M);
    maskfile=fullfile(spmpath,'mask.nii');
    if ~conn_existfile(maskfile), maskfile=fullfile(spmpath,'mask.img'); end
    if ~conn_existfile(maskfile), analysismask=true;
    else analysismask=reshape(spm_get_data(spm_vol(maskfile),xyz),dim)>0;
    end
    am=mean(M(mask&analysismask)); % mean within analysis mask
    if strcmp(gmtype,'SPMam0'), am=am*mean(fm); end % analysis-mean (in units of current regressors/contrasts)
    krun=100./am; 
else
    gm=gm.*fm; % grand-mean (in units of current regressors/contrasts)
    krun=100./gm;
end

% apply scaling factor to contrast images
for ncon=1:numel(contrastname)
    if numel(krun)==1, 
        try, b=spm_read_vols(SPM.xCon(icon(ncon)).Vcon);
        catch, b=spm_read_vols(spm_vol(SPM.xCon(icon(ncon)).Vcon.fname));
        end
        kb=krun*b;
    else
        c=SPM.xCon(icon(ncon)).c;
        assert(size(c,2)==1,'scaling in multivariate contrasts not yet implemented')
        b=0;
        for nrun=1:numel(SPM.Sess),
            for n=reshape(SPM.Sess(nrun).col,1,[])
                if c(n)~=0
                    try, tb=spm_get_data(SPM.Vbeta(n),xyz);
                    catch, tb=spm_get_data(spm_vol(SPM.Vbeta(n).fname),xyz);
                    end
                    b=b+c(n)*krun(nrun)*reshape(tb,dim);
                end
            end
        end
        kb=b;
    end
    if nargout
        varargout{ncon}=kb;
    else
        fname=conn_prepend('',SPM.xCon(icon(ncon)).Vcon.fname,['.',gmtype,'.nii']);
        V=struct('mat',mat,'dim',dim,'pinfo',[1;0;0],'dt',[spm_type('float32') spm_platform('bigend')],'fname',fname,'descrip',sprintf('evlab17_gmscale %s grand-mean scaled contrast k(run)=%s',gmtype,mat2str(krun,6)));
        spm_write_vol(V,kb);
        fprintf('output file %s saved; k(run)=%s\n',fname,mat2str(krun,6));
    end
end

cd(pwd0);

