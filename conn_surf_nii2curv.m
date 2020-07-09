function fileout=conn_surf_nii2curv(filein)
% conn_surf_nii2curv converts surface nifti file to freesurfer paint files
%
% conn_surf_nii2curv(filename)
%     filename : input surface nifti/analyze file (containing fsaverage nvertices*2 voxels; see help conn_surf_sphere)
%
% e.g. conn_surf_nii2curv('con_0001.nii')
%      creates lh.con_0001 and rh.con_0001 freesurfer "curvature" files
%

fileout={};
a=spm_vol(filein);
b=spm_read_vols(a);
b=reshape(b,[],2); % note: assumming both-hemispheres data
nvertices=size(b,1);
resolution=round(log2((nvertices-2)*2/5)/2);
nfaces=5*4^resolution;
fileout={conn_prepend('lh.',filein,''),conn_prepend('rh.',filein,'')};
conn_freesurfer_write_curv(fileout{1},b(:,1),nfaces);
conn_freesurfer_write_curv(fileout{2},b(:,2),nfaces);
