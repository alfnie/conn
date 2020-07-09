function surf=conn_surf_readsurf(filename,refvol,singlefile)
% conn_surf_readsurf reads freesurface surface files
%
% patch = conn_surf_readsurf(filename)
%   reads individual or lh|rh pair of surface files and returns patch structure (with fields faces and vertices)
%   e.g. patch=conn_surf_readsurf('/data/surf/lh.white') % reads individual file
%   e.g. patch=conn_surf_readsurf('/data/surf/white')    % reads lh|rh pair of files
%
% patch = conn_surf_readsurf(filename,vol)
%   returns surface coordinates in world-space using 'vol' reference-volume header information
%   e.g. patch=conn_surf_readsurf('/data/surf/white','/data/mri/T1.mgz')
%

if nargin<1||isempty(filename), filename=fullfile(fileparts(which(mfilename)),'utils','surf','pial.smoothed.surf'); end
if nargin<2, refvol=[]; end
if nargin<3, singlefile=[]; end
if ~isempty(refvol)&&~isstruct(refvol), 
    refvol=conn_file(refvol);
    refvol=refvol{3}(1);
end
if iscell(filename)
    for n=1:numel(filename),
        tsurf=conn_surf_readsurf(filename{n},singlefile,refvol);
        if n==1, surf=tsurf; else surf=[surf, tsurf]; end
    end
else
    if ~isempty(refvol),% convert coordiantes to world coordinates using reference volume header info
        a.vox2ras1=refvol.mat;
        a.volsize=refvol.dim([2 1 3]);
        a.volres = sqrt(sum(refvol.mat(:,1:3).^2,1));
        a.tkrvox2ras=conn_freesurfer_vox2ras_tkreg(a.volsize,a.volres);
        a.vox2ras0=conn_freesurfer_vox2ras_1to0(refvol.mat);
    end
    [file_path,file_name,file_ext]=fileparts(filename);
    
    if isempty(singlefile), singlefile=isequal(file_name,'lh')|isequal(file_name,'rh')|strncmp(file_name,'lh.',3)|strncmp(file_name,'rh.',3); end
    if singlefile
        [vertices,faces]=conn_freesurfer_read_surf(filename);
        surf=struct('vertices',vertices,'faces',faces+1);
    else
        [vertices,faces]=conn_freesurfer_read_surf(fullfile(file_path,['lh.',file_name,file_ext]));
        surf=struct('vertices',vertices,'faces',faces+1);
        [vertices,faces]=conn_freesurfer_read_surf(fullfile(file_path,['rh.',file_name,file_ext]));
        surf(2)=struct('vertices',vertices,'faces',faces+1);
    end
    if ~isempty(refvol)
        for n=1:numel(surf),
            xyz=surf(n).vertices';
            xyz=a.vox2ras0(1:3,:)*pinv(a.tkrvox2ras)*[xyz;ones(1,size(xyz,2))];
            surf(n).vertices=xyz';
        end
    end
end


