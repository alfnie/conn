function [data, vol] = conn_vol_read(filename)
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
data = spm_read_vols(vol);

