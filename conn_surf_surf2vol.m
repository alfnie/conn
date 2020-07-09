function fileout=conn_surf_surf2vol(filein,fileout,FSfolder)
% conn_surf_surf2vol converts surface nifti file to volume nifti file
%
% conn_surf_surf2vol(filename, fileout)
%     filename : input surface nifti file (containing fsaverage 2*nvertices voxels; see help conn_surf_write)
%     fileout  : output volume nifti file (in MNI space)
%
% e.g. conn_surf_curv2nii('curv.surf.nii')
%      creates curv.vol.nii volume nifti file
%

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
if nargin<3||isempty(FSfolder), FSfolder=fullfile(fileparts(which(mfilename)),'utils','surf'); end
    
a1=spm_vol(filein);
b1=spm_read_vols(a1);
mat=[-2 0 0 92;0 2 0 -128;0 0 2 -74;0 0 0 1];
imat=pinv(mat);
dim=[91 109 91];
if ~isdir(FSfolder) % entering alterenative mri/T1.nii file (for non-fsaverage surfaces)
    vol=spm_vol(FSfolder);
    a.vox2ras1=vol.mat;
    a.volsize=vol.dim([2 1 3]);
    a.volres = sqrt(sum(vol.mat(:,1:3).^2,1));
    a.vox2ras0=conn_freesurfer_vox2ras_1to0(vol.mat);
    a.tkrvox2ras=conn_freesurfer_vox2ras_tkreg(a.volsize,a.volres);
    T=a.vox2ras0*pinv(a.tkrvox2ras); %*[xyz_data(:,:);ones(1,size(xyz_data,2)*size(xyz_data,3))];
    imat=imat*T;
    FSfolder=fullfile(fileparts(fileparts(FSfolder)),'surf');
end
try, surfparams={conn_surf_readsurf(fullfile(FSfolder,'lh.white.surf')),conn_surf_readsurf(fullfile(FSfolder,'lh.pial.surf')),conn_surf_readsurf(fullfile(FSfolder,'rh.white.surf')),conn_surf_readsurf(fullfile(FSfolder,'rh.pial.surf'))};
catch, surfparams={conn_surf_readsurf(fullfile(FSfolder,'lh.white')),conn_surf_readsurf(fullfile(FSfolder,'lh.pial')),conn_surf_readsurf(fullfile(FSfolder,'rh.white')),conn_surf_readsurf(fullfile(FSfolder,'rh.pial'))};
end

M=0;
N=0;
for alpha=.05:.10:.95
    surfcoords=[alpha*surfparams{1}.vertices+(1-alpha)*surfparams{2}.vertices;alpha*surfparams{3}.vertices+(1-alpha)*surfparams{4}.vertices];
    recoords=round([surfcoords,ones(size(surfcoords,1),1)]*imat(1:3,:)');
    M=M+accumarray(recoords,b1(:),dim);
    N=N+accumarray(recoords,1,dim);
end
M=M./max(eps,N);
a2=struct('mat',mat,'dim',dim,'fname',fileout,'pinfo',[1;0;0],'n',[1,1],'dt',[spm_type('float32') spm_platform('bigend')]);
spm_write_vol(a2,M);
fprintf('created file %s\n',fileout);
