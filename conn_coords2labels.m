function labels = conn_coords2label(coords, filename)
% labels = conn_coords2label(coords[, atlasfile])
%   coords : 3xM matrix of spatial coordinates (in MNI space, mm units)
%   labels : 1xM cell array of labels (returns an empty label if a voxel is not labeled)
%   atlasfile : (optional) alternative atlas file 
%               default CONN atlas file (conn/rois/atlas.nii)
%

% alfnie@gmail.com 21/01

global CONN_gui
persistent default_atlas;

atlas=[];
isdefault=false;
if size(coords,2)==3&&size(coords,1)~=3, coords=coords.'; end % note: transpose if entering Mx3 coords
if nargin<2||isempty(filename)
    if isfield(CONN_gui,'refs')&&isfield(CONN_gui.refs,'rois')&&isfield(CONN_gui.refs.rois,'filename')&&~isempty(CONN_gui.refs.rois.filename), atlas=CONN_gui.refs.rois; % uses current CONN reference atlas
    else filename=fullfile(fileparts(which('conn')),'rois','atlas.nii'); isdefault=true; % uses default CONN reference atlas
    end
end
if isempty(atlas)
    if isdefault&&~isempty(default_atlas), atlas=default_atlas;
    else 
        [filename_path,filename_name,filename_ext]=fileparts(filename);
        V=spm_vol(filename);
        [idxlabels,strlabels]=rex(filename,filename,'level','clusters','disregard_zeros',false); strlabels=regexprep(strlabels,['^',filename_name,'\.'],'');
        tempdata=spm_read_vols(V); if numel(V)>1, [nill,tempdata]=max(tempdata,[],4); tempdata(~nill)=0; idxlabels=1:numel(strlabels); end % note: for 4D atlases with multiple labels per voxel, return only first
        atlas=struct('filename',filename,'filenameshort',filename_name,'V',V,'data',tempdata,'labels',{strlabels},'labelsidx',full(sparse(1,round(idxlabels(:)'),1:numel(idxlabels))));
        if isdefault, default_atlas=atlas; end
    end
end

coords = round(pinv(atlas.V(1).mat)*[coords; ones(1,size(coords,2))]); % coordinates in voxels
ok = all(coords(1:3,:)>=1,1)&coords(1,:)<=size(atlas.data,1)&coords(2,:)<=size(atlas.data,2)&coords(3,:)<=size(atlas.data,3);
idx = sub2ind(size(atlas.data),coords(1,ok),coords(2,ok),coords(3,ok));
val = atlas.data(idx); % ROI ids
labels=repmat({''},1,numel(val));
labels(ok)=atlas.labels(atlas.labelsidx(val));


