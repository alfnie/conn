function xyz_vol = conn_surf_coords(xyz_surf,ctype)
% conn_surf_coords converts surface- to volume- coordinates
%
% xyz_vol = conn_surf_coords(xyz_surf,ctype)
%   xyz_surf : [Nx1] 1d index to node in two-hemisphere surface (indexes from 1 to 163842*2)
%              OR [Nx3] 3d index to node in two-hemisphere surface volume (indexes from [1 1 1] to [42 83 94])
%   ctype    : reference surface: 'volume','sphere','inflated','semiinflated'
%   xyz_vol  : [Nx3] 3d world-coordinates of corresponding node in reference surface
%

persistent xyz_sphere xyz_volume xyz_inflated;

if nargin<2||isempty(ctype), ctype='inflated'; end
if isempty(xyz_sphere),
    xyz_sphere=conn_surf_sphere;
    xyz_sphere=60*[xyz_sphere.vertices',xyz_sphere.vertices'];
    xyz_volume_left=conn_surf_readsurf(fullfile(fileparts(which('conn')),'utils','surf','lh.white.surf'));
    xyz_volume_right=conn_surf_readsurf(fullfile(fileparts(which('conn')),'utils','surf','rh.white.surf'));
    xyz_volume=[xyz_volume_left.vertices',xyz_volume_right.vertices'];
    xyz_inflated_left=conn_surf_readsurf(fullfile(fileparts(which('conn')),'utils','surf','lh.inflated.surf'));
    xyz_inflated_right=conn_surf_readsurf(fullfile(fileparts(which('conn')),'utils','surf','rh.inflated.surf'));
    xyz_inflated=[xyz_inflated_left.vertices',xyz_inflated_right.vertices'];
    xyz_seminflated_left=conn_surf_readsurf(fullfile(fileparts(which('conn')),'utils','surf','lh.pial.smoothed.surf'));
    xyz_seminflated_right=conn_surf_readsurf(fullfile(fileparts(which('conn')),'utils','surf','rh.pial.smoothed.surf'));
    xyz_seminflated=[xyz_seminflated_left.vertices',xyz_seminflated_right.vertices'];
end

if isempty(xyz_surf), idx=1:size(xyz_volume,2); 
elseif size(xyz_surf,1)==3, idx=sub2ind(conn_surf_dims(8).*[1,1,2],xyz_surf(1,:),xyz_surf(2,:),xyz_surf(3,:));
else idx=xyz_surf;
end


switch(lower(ctype))
    case 'volume'
        xyz_vol=xyz_volume(:,idx);
    case 'sphere'
        xyz_vol=xyz_sphere(:,idx);
    case 'inflated',
        xyz_vol=xyz_inflated(:,idx);
    case {'semiinflated','semi-inflated','seminflated'}
        xyz_vol=xyz_seminflated(:,idx);
    otherwise
        error('unrecognized option %s (valid options: volume/sphere/inflated)',lower(ctype));
end
end