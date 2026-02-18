function fh = conn_mtx_display(filename, varargin)
% CONN_MTX_DISPLAY displays matrix data from nifti file
%
% conn_mtx_display(filename) displays matrix numeric data from filename
%  filename : input filename *.mtx.nii
%
% alternative syntax:
% conn_mtx_display(data [, names, coords, samples, borders, borders_names])
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
if isempty(coords)&&isempty(samples)&&isempty(clusters) % tries to get this info from available CONN files
    info=[];
    if all(cellfun('length',regexp(names,'^schaefer\.'))>0) % schaefer atlas files
        info=conn_roiclusters_load(fullfile(fileparts(which('conn')),'rois','schaefer.groups.mat'),names);
    elseif all(cellfun('length',regexp(names,'^atlas\.'))>0) % schaefer atlas files
        info=conn_roiclusters_load(fullfile(fileparts(which('conn')),'rois','atlas.groups.mat'),names);
    elseif all(cellfun('length',regexp(names,'^networks\.'))>0) % schaefer atlas files
        info=conn_roiclusters_load(fullfile(fileparts(which('conn')),'rois','networks.groups.mat'),names);
    end
    if ~isempty(info)
        idx1=info.displaytheserois;
        [nill,idx2]=sort(info.clusters(idx1));
        idx=idx1(idx2);
        data=data(idx,idx,:,:);
        clusters=find(diff(info.clusters(idx)))+.5;
        clusters_names=info.names_clusters;
        coords=info.xyz2(idx,:);
        names=info.names(idx);
    end
end
fh=conn_montage_display(...
    permute(data,[1,2,4,3]),...
    samples,...
    'matrix',...
    [],{},... % cov
    clusters,clusters_names,... % borders
    [],{},... % clusters
    coords, names ... % rois
    );

