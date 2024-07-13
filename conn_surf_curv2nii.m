function fileout=conn_surf_curv2nii(filein,filein2,fileout)
% conn_surf_curv2nii converts freesurfer paint files to surface nifti file
%
% conn_surf_curv2nii(filename)
%     filename : input lh curvature file (containing fsaverage nvertices voxels; see help conn_freesurfer_read_curv)
%
% e.g. conn_surf_curv2nii('lh.curv')
%      creates curv.surf.nii surface nifti file
%

if size(filein,1)>1, 
    fileout=char(cellfun(@conn_surf_curv2nii,cellstr(filein),'uni',0));
    return
end
if iscell(filein), 
    fileout=cellfun(@conn_surf_curv2nii,filein,'uni',0);
    return
end

fileout={};
[filepath,filename,fileext]=fileparts(filein);
filename=[filename,fileext];
if nargin<2||isempty(filein2)
    assert(~isempty(regexp(filename,'^lh\.')),'input filename must be of the form lh.*');
    filein2=fullfile(filepath,[regexprep(filename,'^lh\.','rh.')]);
end
if nargin<3||isempty(fileout)
    assert(~isempty(regexp(filename,'^lh\.')),'input filename must be of the form lh.*');
    fileout=fullfile(filepath,[regexprep(filename,'^lh\.',''),'.surf.nii']);
end
a1=conn_freesurfer_read_curv(filein);
a2=conn_freesurfer_read_curv(filein2);
conn_surf_write(fileout,[a1(:);a2(:)]);
fprintf('created file %s\n',fileout);