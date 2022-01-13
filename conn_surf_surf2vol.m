function fileout=conn_surf_surf2vol(filein,fileout,FSfolder,interp, aggreg)
% CONN_SURF_SURF2VOL converts surface nifti file to volume nifti file
%
% conn_surf_surf2vol(filename [, fileout, folderREF, interp, aggreg])
%     filename  : input surface nifti file (containing fsaverage 2*nvertices voxels; see help conn_surf_write)
%     fileout   : output volume nifti file (in MNI space)
%     folderREF : folder with reference surfaces (default conn/utils/surf)
%     interp    : interpolation sample(s) along cortical surface perpendicular (0:white 1:pial; default .05:.10:.95 taking the average across 10 samples along the normal between the white and pial surfaces)
%     aggreg    : aggregation function (default @mean; switch to @mode for categorical data to avoid the resampling operation to interpolate across different categories)
%
% e.g. conn_surf_curv2nii('curv.surf.nii')
%      creates curv.vol.nii volume nifti file
%

if nargin<2, fileout=[]; end
if nargin<3, FSfolder=[]; end
if nargin<4, interp=[]; end
if nargin<5, aggreg=[]; end
if any(conn_server('util_isremotefile',filein)), fileout=conn_server('util_remotefile',conn_server('run',mfilename,conn_server('util_localfile',filein),conn_server('util_localfile',fileout),FSfolder,interp,aggreg)); return; end

if size(filein,1)>1, 
    fileout=char(cellfun(@conn_surf_surf2vol,cellstr(filein),'uni',0));
    return
end
if iscell(filein), 
    fileout=cellfun(@conn_surf_surf2vol,filein,'uni',0);
    return
end

if nargin<2||isempty(fileout), 
    if ~isempty(regexp(filein,'\.surf.nii$')), fileout=regexprep(filein,'\.surf.nii$','.vol.nii');
    else
        [filepath,filename,fileext]=fileparts(filein);
        fileout=fullfile(filepath,[filename,'.vol',fileext]);
    end
end        
isfsaverage=false;
if nargin<3||isempty(FSfolder), FSfolder=fullfile(fileparts(which(mfilename)),'utils','surf'); isfsaverage=true; end
if nargin<4||isempty(interp), interp=.05:.10:.95; end % interpolation samples (0=white 1=pial)
if nargin<5||isempty(aggreg), aggreg=[]; end
    
filein=conn_server('util_localfile',filein);
fileout=conn_server('util_localfile',fileout);
a1=spm_vol(filein);
b1=spm_read_vols(a1);
mat=[-2 0 0 92;0 2 0 -128;0 0 2 -74;0 0 0 1];
imat=pinv(mat);
dim=[91 109 91];
if ~isfsaverage
    files2=cellfun(@(x)fullfile(fileparts(FSfolder),'mri',x),{'T1.nii','T1.mgh','T1.mgz','brain.nii','brain.mgh','brain.mgz'},'uni',0);
    existfiles2=conn_existfile([files2]);
    if any(existfiles2) % entering alternative mri/T1.nii file (for non-fsaverage surfaces)
        out=conn_file(files2{find(existfiles2,1)});
        vol=spm_vol(out{1});
        a.vox2ras1=vol.mat;
        a.volsize=vol.dim([2 1 3]);
        a.volres = sqrt(sum(vol.mat(:,1:3).^2,1));
        a.vox2ras0=conn_freesurfer_vox2ras_1to0(vol.mat);
        a.tkrvox2ras=conn_freesurfer_vox2ras_tkreg(a.volsize,a.volres);
        T=a.vox2ras0*pinv(a.tkrvox2ras); %*[xyz_data(:,:);ones(1,size(xyz_data,2)*size(xyz_data,3))];
        %imat=imat*T;
        % note: create volume in subject-specific space (same space as T1.nii file)
        mat=vol.mat;
        imat=pinv(mat)*T;
        dim=vol.dim;
        %FSfolder=fullfile(fileparts(fileparts(FSfolder)),'surf');
    else isfsaverage=true;
    end
end
    
try, surfparams={conn_surf_readsurf(fullfile(FSfolder,'lh.white.surf')),conn_surf_readsurf(fullfile(FSfolder,'lh.pial.surf')),conn_surf_readsurf(fullfile(FSfolder,'rh.white.surf')),conn_surf_readsurf(fullfile(FSfolder,'rh.pial.surf'))};
catch, surfparams={conn_surf_readsurf(fullfile(FSfolder,'lh.white')),conn_surf_readsurf(fullfile(FSfolder,'lh.pial')),conn_surf_readsurf(fullfile(FSfolder,'rh.white')),conn_surf_readsurf(fullfile(FSfolder,'rh.pial'))};
end

resolution=8;
nvertices2=2+10*2^(2*resolution-2);
matchdims=cellfun(@(a)size(a.vertices,1)==nvertices2,surfparams);
if conn_surf_dimscheck(a1)&&~all(matchdims) % fsaverage input surface nifti file but non-fsaverage folderREF surfaces
    if ~conn_existfile(fullfile(FSfolder,'lh.sphere.reg'))||~conn_existfile(fullfile(FSfolder,'rh.sphere.reg')), error('unable to find [lr]h.sphere.reg files to match fsaverage file to non-fsaverage surfaces'); end
    xyz_ref1=conn_freesurfer_read_surf(fullfile(FSfolder,'lh.sphere.reg'));
    xyz_ref2=conn_freesurfer_read_surf(fullfile(FSfolder,'rh.sphere.reg'));
    [xyz_sphere1,sphere2ref1,ref2sphere1]=conn_surf_sphere(resolution,xyz_ref1);
    [xyz_sphere2,sphere2ref2,ref2sphere2]=conn_surf_sphere(resolution,xyz_ref2);
    b1=reshape(b1,[],2);
    b1=[b1(sphere2ref1,1);b1(sphere2ref2,2)];
end
M=0;
N=0;
Recoords=[]; B1=[];
for alpha=interp(:)'
    surfcoords=[alpha*surfparams{1}.vertices+(1-alpha)*surfparams{2}.vertices;alpha*surfparams{3}.vertices+(1-alpha)*surfparams{4}.vertices];
    recoords=round([surfcoords,ones(size(surfcoords,1),1)]*imat(1:3,:)');
    Recoords=[Recoords; recoords];
    B1=[B1;b1(:)];
end
M=accumarray(Recoords,B1(:),dim,aggreg);
%    M=M+accumarray(recoords,b1(:),dim);
%    N=N+accumarray(recoords,1,dim);
%end
%M=M./max(eps,N);
if numel(alpha)==1, dt=a1.dt;
else dt=[spm_type('float32') spm_platform('bigend')];
end
a2=struct('mat',mat,'dim',dim,'fname',fileout,'pinfo',[1;0;0],'n',[1,1],'dt',dt);
spm_write_vol(a2,M);

%exts=[exts {'.txt','.csv','.xls','.info','.icon.jpg','.json'}];

fprintf('created file %s\n',fileout);
