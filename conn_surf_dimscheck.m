function [issurf,ismtx] = conn_surf_dimscheck(dim)
if isstruct(dim), dim=dim(1).dim; 
elseif iscell(dim)||ischar(dim), dim=conn_fileutils('nifti',dim); dim=size(dim(1).dat); 
end
if numel(dim)<3, dim=[dim(:)',ones(1,3)]; end
dim=reshape(dim(1:3),1,3);
issurf = isequal(sort(dim),sort(conn_surf_dims(8).*[1 1 2]))||isequal(sort(dim),sort(conn_surf_dims(8).*[1 2 1]));
ismtx = dim(3)==1&dim(1)==dim(2);
