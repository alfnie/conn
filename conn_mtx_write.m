function a = conn_mtx_write(filename,data,names,coords,samples)
% CONN_MTX_WRITE writes matrix data to nifti file
%
% conn_mtx_write(filename,data [,names, coords, labels])
%  filename : output filename *.mtx.nii
%  data     : input data [nrois, nrois, nsamples] matrix
%  names    : (optional) list of ROI names {1,nrois} cell array (default ROI#1,ROI#2,...)
%  coords   : (optional) list of ROI xyz coordinates {1,nrois} matrix (default [0 0 0],[0 0 0],...)
%  samples  : (optional) list of sample labels {1,nsamples} cell array (default Sample#1,Sample#2,...)
%
% creates *.mtx.nii and *.mtx.json file with ROI-to-ROI matrix data
%

if ~nargin, help(mfilename); return; end

M=size(data,1);
N=size(data,3);
assert(M==size(data,2),'input matrix must be square'); 
dims=[M,size(data,2),1];
data=reshape(data,[dims,N]);
a=struct('fname',filename,'mat',eye(4),'dim',dims,'n',[1,1],'pinfo',[1;0;0],'dt',[spm_type('float32') spm_platform('bigend')]);
spm_unlink(filename);
a=repmat(a,1,N); for n=1:N, a(n).n=[n,1]; end
a=spm_create_vol(a);
for n=1:N, a(n)=spm_write_vol(a(n),data(:,:,:,n)); end

if nargin<3||isempty(names), names=arrayfun(@(n)sprintf('ROI#%04d',n),1:M,'uni',0); end
if nargin<4||isempty(coords), coords=repmat({[0 0 0]},1,M); end
if nargin<5||isempty(samples), samples=arrayfun(@(n)sprintf('Sample#%04d',n),1:N,'uni',0); end
spm_jsonwrite(conn_prepend('',filename,'.json'),struct('names',{names},'coords',{coords},'samples',{samples}),struct('indent',' '));

