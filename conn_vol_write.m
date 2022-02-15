function a = conn_vol_write(filename,data,M,ftype)
% CONN_VOL_WRITE writes volume data to nifti file
% conn_vol_write(filename, data [, mat, ftype])
%  filename : *.nii
%  data     : [Ni, Nj, Nk, nobservations] data 3D/4D matrix
%  mat      : [4,4] affine voxel-to-world transformation matrix (e.g. vol.mat; see "help spm_vol")
%  ftype    : datatype (e.g. vol.dt(1); see "help spm_type")


if ~nargin, help(mfilename); return; end

if nargin<3||isempty(M), M=eye(4); end
if nargin<4||isempty(ftype), if isstruct(M), ftype=M.dt(1); else ftype=spm_type('float32'); end; end
if isstruct(M), M=M.mat; end
isremotefile=conn_server('util_isremotefile',filename);
if isremotefile, remotefilename=filename; filename=conn_cache('new',remotefilename); end

N=size(data,4);
dims=[size(data,1) size(data,2) size(data,3)];
a=struct('fname',filename,'mat',M,'dim',dims,'n',[1,1],'pinfo',[1;0;0],'dt',[ftype spm_platform('bigend')]);
spm_unlink(filename);
a=repmat(a,1,N); for n=1:N, a(n).n=[n,1]; end
a=spm_create_vol(a);
for n=1:N, a(n)=spm_write_vol(a(n),data(:,:,:,n)); end

if isremotefile, conn_cache('push',remotefilename); a=[]; end

