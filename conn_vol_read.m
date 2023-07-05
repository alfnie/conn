function [data, vol] = conn_vol_read(filename, fileref, hold)
% CONN_VOL_READ reads volume nifti file
%
% [data, vol] = conn_vol_read(filename) reads 3D/4D volume data from filename
%  filename : input filename *.nii
%  data     : output data [Ni, Nj, Nk, nobservations]
%  vol      : header information (see "help spm_vol")
%              vol.mat      : [4,4] affine voxel-to-world transformation matrix
%              vol.dt(1)    : datatype (see "help spm_type")
%

if ~nargin, help(mfilename); return; end
isremotefile=conn_server('util_isremotefile',filename);
if isremotefile, remotefilename=filename; filename=conn_cache('pull',remotefilename); end

vol = spm_vol(filename);
if nargin>1&&~isempty(fileref)
    if nargin<=2||isempty(hold), hold=0; end
    volref = spm_vol(fileref);
    [x,y,z]=ndgrid(1:volref(1).dim(1), 1:volref(1).dim(2), 1:volref(1).dim(3));
    xyz=volref(1).mat*[x(:) y(:) z(:) ones(numel(x),1)]';
    data = [];
    for nvol=1:numel(vol), 
        txyz = pinv(vol(nvol).mat)*xyz;
        data = cat(4,data, reshape(spm_sample_vol(vol(nvol),txyz(1,:),txyz(2,:),txyz(3,:),hold),volref(1).dim)); 
    end
else
    try
        data = spm_read_vols(vol);
    catch
        data=[];
        for nvol=1:numel(vol)
            tdata=spm_read_vols(vol(nvol));
            data=cat(4,data,tdata);
        end
    end
end

