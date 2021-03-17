function data = conn_surf_read(filename)
% CONN_SURF_READ reads surface nifti file
%
% data = conn_surf_read(filename) reads surface numeric data from filename
%  filename : input filename *.surf.nii
%  data     : output data [163842 vertices, 1|2 hemispheres, nsamples] 3d matrix
%

if ~nargin, help(mfilename); return; end
isremotefile=conn_server('util_isremotefile',filename);
if isremotefile, remotefilename=filename; filename=conn_cache('pull',remotefilename); end

vol = spm_vol(filename);
data = spm_read_vols(vol);
data = reshape(data, 163842, [], size(data,4));

