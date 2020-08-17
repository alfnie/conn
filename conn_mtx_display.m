function fh = conn_mtx_display(filename, varargin)
% CONN_MTX_DISPLAY displays matrix data from nifti file
%
% conn_mtx_display(filename) displays matrix numeric data from filename
%  filename : input filename *.mtx.nii
%
% alternative syntax:
% conn_mtx_display(data [, names, coords, samples])
% SEE ALSO: conn_mtx_read, conn_mtx_write
% 

if ~nargin, help(mfilename); return; end

if ischar(filename)
    [data,names,coords,samples] = conn_mtx_read(filename);
    if numel(varargin)>=1, clusters=varargin{1}; else clusters=[]; end
    if numel(varargin)>=2, clusters_names=varargin{2}; else clusters_names={}; end
else
    data=filename;
    if numel(varargin)>=1, names=varargin{1}; else names={}; end
    if numel(varargin)>=2, coords=varargin{2}; else coords={}; end
    if numel(varargin)>=3, samples=varargin{3}; else samples={}; end
    if numel(varargin)>=4, clusters=varargin{4}; else clusters=[]; end
    if numel(varargin)>=5, clusters_names=varargin{5}; else clusters_names={}; end
end
fh=conn_montage_display(...
    permute(data,[1,2,4,3]),...
    samples,...
    'matrix',...
    [],{},... % cov
    [],{},... % borders
    clusters,clusters_names,... % clusters
    coords, names ... % rois
    );

