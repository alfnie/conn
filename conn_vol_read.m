function [data, vol] = conn_vol_read(filename)
% CONN_VOL_READ reads volume nifti file
%
% [data, vol] = conn_vol_read(filename) reads 3D/4D volume data from filename
%  filename : input filename *.nii
%  data     : output data [Ni, Nj, Nk, nobservations]
%  vol      : header information (see "help spm_vol")
%

if ~nargin, help(mfilename); return; end
isremotefile=conn_server('util_isremotefile',filename);
if isremotefile, remotefilename=filename; filename=conn_cache('pull',remotefilename); end

vol = spm_vol(filename);
data = spm_read_vols(vol);

