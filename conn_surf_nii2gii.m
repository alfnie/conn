function fileout=conn_surf_nii2gii(filein,fileout)
% conn_surf_nii2gii converts surface nifti files to freesurfer gifti files
%
% conn_surf_nii2gii(filename)
%     filename : input surface nifti file (containing fsaverage nvertices voxels from both hemispheres)
%
% e.g. conn_surf_nii2gii('bold.surf.nii')
%      creates lh.bold.surf.gii and rh.bold.surf.gii freesurfer gifti files
%

if nargin<2||isempty(fileout), fileout=[]; end
if any(conn_server('util_isremotefile',filein)), fileout=conn_server('util_remotefile',conn_server('run',mfilename,conn_server('util_localfile',filein),conn_server('util_localfile',fileout))); return; end

if size(filein,1)>1, 
    fileout=char(cellfun(@conn_surf_nii2gii,cellstr(filein),'uni',0));
    return
end
if iscell(filein), 
    fileout=cellfun(@conn_surf_nii2gii,filein,'uni',0);
    return
end
if nargin<2||isempty(fileout), fileout=filein; end

data=conn_surf_read(filein);
save(gifti(permute(data(:,1,:),[1,3,2])),conn_prepend('lh.',fileout,'.gii'));
save(gifti(permute(data(:,2,:),[1,3,2])),conn_prepend('rh.',fileout,'.gii'));
fprintf('created files %s & %s\n',conn_prepend('lh.',fileout,'.gii'),conn_prepend('rh.',fileout,'.gii'));
