
function fh = conn_mtx_braindisplay(filename, threshold)
% conn_mtx_braindisplay displays matrix data from nifti file
%
% conn_mtx_braindisplay(filename) displays matrix numeric data from filename
%  filename : input filename *.mtx.nii
%
% alternative syntax:
% conn_mtx_braindisplay(data, threshold) % displays only connections with abs(data)>threshold
% conn_mtx_braindisplay(data, 'top100')  % displays only top 100 strongest connections
% conn_mtx_braindisplay(data, 'perc95')  % displays only connections with abs(data) above its 95% percentile
%
% SEE ALSO: conn_mtx_read, conn_mtx_write
%

if ~nargin, help(mfilename); return; end
[data,names,coords,samples] = conn_mtx_read(filename);
if nargin<2, threshold=prctile(abs(data(:)),99); end

n=size(data,1);
if size(data,3)>1, data=mean(data,3); end
if n~=3&&size(coords,1)==3, coords=coords'; end
if ischar(threshold)&&~isempty(regexp(threshold,'^top'))
    [nill,idx]=sort(abs(data(:)),'descend');
    idx(idx)=1:numel(idx);
    mask=reshape(idx,size(data))<=str2double(regexprep(threshold,'^top',''));
elseif ischar(threshold)&&~isempty(regexp(threshold,'^perc'))
    mask=abs(data)>prctile(abs(data(:)),str2double(regexprep(threshold,'^perc','')));
else
    mask=abs(data)>threshold;
end
[i,j]=find(mask); 
idx=unique(union(i,j)); 
fh=conn_mesh_display('','',[],...
    struct('sph_names',{names(idx)},'sph_xyz',coords(idx,:),...
    'sph_r',3*ones(numel(idx),1)),...
    data(idx,idx).*mask(idx,idx),... 
    [], .2, [0,-1e-8,1]);
