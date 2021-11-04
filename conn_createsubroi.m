function [fileout,Nrois]=conn_createsubroi(filein, fileout, Nrois, D, fileinSelect)
% CONN_CREATESUBROI creates nifti ROI file/atlas dividing an existing mask into smaller parcels along one axis (e.g. x/y/z)
% 
% conn_createsubroi(filein, fileout, N, xyz [,val])
%    filein     : input mask filename
%    fileout    : output ROI filename (by default [filein].ROI.nii)
%    N          : number of parcels/ROIs in the output file
%    xyz        : [1x3] direction vector defining spatial axis, output ROIs will divide the original mask using planes 
%                perpendicular to the defined xyz direction vector (e.g. xyz = [1 0 0] to divide an ROI along the x- axis)
%                and placed such that there is an approximate equal number of voxels within each parcel
%                note: set xyz = [] to select the principal/longest axis of the original mask/ROI
%    val        : (optional, for input ROI files) only select from the input file voxels with values = val
%                (e.g. to select a single ROI within an atlas file). The default behavior is to select all 
%                voxels with values > 0 in the input mask file
%

if nargin<2||isempty(fileout), fileout=conn_prepend('',filein,'.ROI.nii'); end
if nargin<3||isempty(Nrois), Nrois=2; end
if nargin<4||isempty(D), D=[]; end
if nargin<5||isempty(fileinSelect), fileinSelect=[]; end
if any(conn_server('util_isremotefile',{fileout, filein})), 
    [fileout,Nrois]=conn_server('run',mfilename,conn_server('util_localfile',fileout),conn_server('util_localfile',filein), Nrois, D); 
    fileout=conn_server('util_remotefile',fileout);
    return
end

assert(Nrois<=65535,'unable to create atlas file with more than 65535 ROIs');
vol=spm_vol(filein);
mask=spm_read_vols(vol);
if isempty(fileinSelect), idx=reshape(find(mask>0),[],1);
else idx=reshape(find(mask==fileinSelect),[],1);
end
if numel(idx)<Nrois, fprintf('warning: Nrois reduced to %d (total number of voxels in mask)\n',numel(idx)); Nrois=numel(idx); end

[x,y,z]=ndgrid(1:vol(1).dim(1),1:vol(1).dim(2),1:vol(1).dim(3));
xyz=vol.mat*[x(idx) y(idx) z(idx) ones(size(idx))]';
if isempty(D) % principal axis (note: consider adding graph-based distance partitioning)
    [nill,nill,D]=svd(xyz'-repmat(mean(xyz,2)',size(xyz,2),1),0);
    D=D(1:3,1);
end
[nill,I]=sort(D(:)'*xyz(1:3,:));
I(I)=ceil((1:numel(idx))/numel(idx)*Nrois);

b=zeros(size(mask));
b(idx)=I;
vol=struct('fname',fileout, 'mat',vol.mat, 'dim',vol.dim, 'dt', [spm_type('uint16') spm_platform('bigend')], 'pinfo',[1;0;0],'descrip',''); 
try, spm_unlink(fileout); end
vol=spm_write_vol(vol,b);
if Nrois>1
    fh={};
    for n=1:Nrois, mxyz=mean(xyz(:,I==n),2); fh{end+1}=sprintf('(%d,%d,%d)\n',round(mxyz(1)),round(mxyz(2)),round(mxyz(3))); end
    conn_fileutils('filewrite_raw',conn_prepend('',fileout,'.txt'), fh);
end
end
