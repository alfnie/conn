function surf=conn_surf_readcurvresampled(filenames,dowrite)
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
    data_ref=conn_freesurfer_read_curv(filename);
    resolution=8;
    ref=conn_surf_readsurf(fileref);
    [xyz_sphere,sphere2ref,ref2sphere]=conn_surf_sphere(resolution,ref.vertices);
    Mref2sphere = conn_surf_interpmtx(ref,xyz_sphere.vertices,ref2sphere);
    data_new=Mref2sphere*data_ref;
    %xyz_ref=conn_freesurfer_read_surf(fileref);
    %[xyz_sphere,sphere2ref,ref2sphere]=conn_surf_sphere(resolution,xyz_ref);
    %data_new=data_ref(ref2sphere);
    if numel(filenames)==1, surf=data_new;
    else surf{n}=data_new;
    end
    if dowrite
        nvertices=numel(data_new);
        resolution=round(log2((nvertices-2)*2/5)/2);
        nfaces=5*4^resolution;
        conn_freesurfer_write_curv([filename,'.fsavg'],data_new,nfaces);
    end
end
    
