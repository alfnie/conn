function [fileout,Nrois]=conn_createvoxelroi(filein, fileout)
% CONN_CREATEVOXELROI creates nifti ROI file/atlas identifying individual voxels within a mask
% 
% conn_createvoxelroi(filein, fileout)
%    filein     : input mask filename
%    fileout    : output ROI filename (by default [filein].ROI.nii)

if nargin<2||isempty(fileout), fileout=conn_prepend('',filein,'.ROI.nii'); end
if any(conn_server('util_isremotefile',{fileout, filein})), 
    [fileout,Nrois]=conn_server('run',mfilename,conn_server('util_localfile',fileout),conn_server('util_localfile',filein)); 
    fileout=conn_server('util_remotefile',fileout);
    return
end

vol=spm_vol(filein);
mask=spm_read_vols(vol);
idx=reshape(find(mask>0),[],1);
Nrois=numel(idx);
assert(Nrois<=65535,'unable to create atlas file with more than 65535 ROIs');
b=zeros(size(mask));
b(idx)=1:numel(idx);
[x,y,z]=ndgrid(1:vol(1).dim(1),1:vol(1).dim(2),1:vol(1).dim(3));
xyz=vol.mat*[x(idx) y(idx) z(idx) ones(size(idx))]';
vol=struct('fname',fileout, 'mat',vol.mat, 'dim',vol.dim, 'dt', [spm_type('uint16') spm_platform('bigend')], 'pinfo',[1;0;0],'descrip',''); 
try, spm_unlink(fileout); end
vol=spm_write_vol(vol,b);
if Nrois>1
    fh={};
    for n=1:Nrois, fh{end+1}=sprintf('(%d,%d,%d)\n',round(xyz(1,n)),round(xyz(2,n)),round(xyz(3,n))); end
    conn_fileutils('filewrite_raw',conn_prepend('',fileout,'.txt'), fh);
end
end
