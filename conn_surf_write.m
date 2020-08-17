function a = conn_surf_write(filename,data)
% CONN_SURF_WRITE writes surface data to surface nifti file
% conn_surf_write(filename,data)
%  filename : *.surf.nii
%  data     : [nvertices, nobservations]
%             (note: full density fsaverage space, nvertices=163842*2 for both hemispheres data)

if ~nargin, help(mfilename); return; end

dims=conn_surf_dims(8);
N=size(data,2);
switch(size(data,1))
    case prod(dims),    
    case 2*prod(dims),  dims=dims.*[1 1 2];
    otherwise,          error('incompatible number of vertices %d (expected %d or %d vertices for full density fsaverage space)',size(data,1),prod(dims),2*prod(dims));
end
data=reshape(data,[dims N]);
a=struct('fname',filename,'mat',eye(4),'dim',dims,'n',[1,1],'pinfo',[1;0;0],'dt',[spm_type('float32') spm_platform('bigend')]);
spm_unlink(filename);
a=repmat(a,1,N); for n=1:N, a(n).n=[n,1]; end
a=spm_create_vol(a);
for n=1:N, a(n)=spm_write_vol(a(n),data(:,:,:,n)); end


