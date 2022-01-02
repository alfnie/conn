function fileout=conn_surf_gii2nii(filein,filein2,fileout)
% conn_surf_gii2nii converts freesurfer gifti files to surface nifti file
%
% conn_surf_gii2nii(filename)
%     filename : input lh gifti file (containing fsaverage nvertices voxels; see help conn_freesurfer_read_curv)
%
% e.g. conn_surf_gii2nii('lh.curv.gii')
%      creates curv.surf.nii surface nifti file
%      from data in lh.curv.gii and rh.curv.gii gifti files
%

if nargin<2||isempty(filein2), filein2=[]; end
if nargin<3||isempty(fileout), fileout=[]; end
if any(conn_server('util_isremotefile',filein)), fileout=conn_server('util_remotefile',conn_server('run',mfilename,conn_server('util_localfile',filein),conn_server('util_localfile',filein2),conn_server('util_localfile',fileout))); return; end

if size(filein,1)>1, 
    fileout=char(cellfun(@conn_surf_gii2nii,cellstr(filein),'uni',0));
    return
end
if iscell(filein), 
    fileout=cellfun(@conn_surf_gii2nii,filein,'uni',0);
    return
end

filein=conn_server('util_localfile',filein);
[filepath,filename,fileext]=fileparts(filein);
if nargin<2||isempty(filein2)
    if ~isempty(regexp(filename,'^lh\.'))
        filein2=fullfile(filepath,[regexprep(filename,'^lh\.','rh.'),fileext]);
    elseif ~isempty(regexp(filename,'_hemi-L'))
        filein2=fullfile(filepath,[regexprep(filename,'_hemi-L','_hemi-R'),fileext]);
    else error('input filename must be of the form lh.* or *_hemi-L*');
    end
end
if nargin<3||isempty(fileout)
    if ~isempty(regexp(filename,'^lh\.'))
        fileout=fullfile(filepath,[regexprep(filename,'^lh\.',''),fileext,'.surf.nii']);
    elseif ~isempty(regexp(filename,'_hemi-L'))
        fileout=fullfile(filepath,[regexprep(filename,'_hemi-L',''),'.surf.nii']);
    else error('input filename must be of the form lh.* or *_hemi-L*');
    end
end
a1=gifti(filein);
a2=gifti(filein2);
conn_surf_write(fileout,[a1.cdata;a2.cdata]);
fprintf('created file %s\n',fileout);
