function data = conn_surf_read(filename)
% CONN_SURF_READ reads surface nifti file
%
% data = conn_surf_read(filename) reads surface numeric data from filename
%  filename : input filename *.surf.nii
%  data     : output data [nvertices, nobservations]
%             (note: full density fsaverage space, nvertices=163842*2 for both hemispheres data)
%

if ~nargin, help(mfilename); return; end
if conn_server('util_isremotefile',filename), filename=conn_cache('pull',filename); end

vol = spm_vol(filename);
data = spm_read_vols(vol);
data = reshape(data, [], size(data,4));

