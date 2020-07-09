function surf=conn_surf_readsurfresampled(filenames,dowrite,varargin)
if nargin<2||isempty(dowrite), dowrite=false; end
if ~iscell(filenames), filenames=cellstr(filenames); end

for n=1:numel(filenames)
    filename=filenames{n};
    [file_path,file_name,file_ext]=fileparts(filename);
    file_name=[file_name,file_ext];
    if strncmp(file_name,'lh.',3), hem='lh';
    elseif strncmp(file_name,'rh.',3), hem='rh';
    else error(['unable to determine hemisphere of file ',filename]);
    end
    fileref=fullfile(file_path,[hem,'.sphere.reg']);
    if ~conn_existfile(fileref), error(['unable to find file ',fileref]); end
    
    % resample at sphere reference grid
    data_ref=conn_surf_readsurf(filename,varargin{:});
    resolution=8;
    ref=conn_surf_readsurf(fileref);
    [xyz_sphere,sphere2ref,ref2sphere]=conn_surf_sphere(resolution,ref.vertices);
    Mref2sphere = conn_surf_interpmtx(ref,100*xyz_sphere.vertices,ref2sphere);
    xyz=Mref2sphere*data_ref.vertices;
    %xyz=data_ref(ref2sphere,:);
    faces=xyz_sphere.faces;
    surf(n)=struct('vertices',xyz,'faces',faces);
    if dowrite
        conn_freesurfer_write_surf([filename,'.fsavg'],xyz,faces-1);
    end
end
    
