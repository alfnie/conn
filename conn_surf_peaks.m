function varargout=conn_surf_peaks(filename,ref)
% CONN_SURF_PEAKS finds local maxima in .surf.nii file
%
% xyz_surf=conn_surf_peaks(FILENAME)
% returns list of local maxima indexes (xyz_surf)
%
% note: 
%  additional parameters "conn_surf_peaks(...,REF)" reference file defining surface faces (default 'conn/utils/surf/lh.pial.surf')
%
% see also: conn_surf_coords
%

if nargin<2, ref=[]; end
if any(conn_server('util_isremotefile',filename)), filename=conn_server('util_remotefile',conn_server('run',mfilename,conn_server('util_localfile',filename),nsmooth,ref)); return; end

if isempty(ref), ref=fullfile(fileparts(which('conn')),'utils','surf','lh.pial.surf'); end
[xyz,faces]=conn_freesurfer_read_surf(ref);
faces=faces+1;
iscellfilename=iscell(filename);
if iscellfilename, filename=char(filename); end
if ~ischar(filename), 
    a=[];
    b=filename;
    filename='manual data'; 
else 
    filename=conn_server('util_localfile',filename);
    a=spm_vol(filename);
    b=spm_read_vols(a);
end
if rem(numel(b),2*163842), error('Incorrect dimensions in file %s (%d voxels, expected 163842*2)',filename,numel(b)); end
sB=size(b);
b=reshape(b,163842,[]);
b(isnan(b))=0;
A=double(sparse(repmat(faces,3,1),repmat(faces,1,3), 1)>0);
A=double(A|speye(size(A,1)));
[i,j]=find(A);
c=false(size(b));
for n=1:size(b,2)
    tc=accumarray(i,b(j,n),[size(A,1),1],@max);
    c(:,n)=tc==b(:,n)&b(:,n)~=0;
end
b=reshape(b,2*163842,[]);
c=reshape(c,2*163842,[]);
for n=1:size(c,2)
    idx1=find(c(:,n));
    [nill,idx2]=sort(b(idx1),'descend');
    varargout{n}=idx1(idx2);
end


