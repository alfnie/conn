function fileout=conn_surf_mgh2nii(filein,filein2,fileout)
% conn_surf_mgh2nii converts freesurfer .mgh surface paint files to .nii surface nifti file
%
% conn_surf_mgh2nii(filename)
%     filename : input lh curvature file (containing fsaverage nvertices voxels; see help conn_freesurfer_read_curv)
%
% e.g. conn_surf_mgh2nii('lh.data.mgh')
%      creates data.surf.nii surface nifti file
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
if nargin<2||isempty(filein2)
    assert(~isempty(regexp(filename,'^lh\.')),'input filename must be of the form lh.*');
    filein2=fullfile(filepath,[regexprep(filename,'^lh\.','rh.'),fileext]);
end
if nargin<3||isempty(fileout)
    assert(~isempty(regexp(filename,'^lh\.')),'input filename must be of the form lh.*');
    fileout=fullfile(filepath,[regexprep(filename,'^lh\.',''),fileext,'.surf.nii']);
end
a1=conn_freesurfer_load_mgh(filein);
a2=conn_freesurfer_load_mgh(filein2);
conn_surf_write(fileout,[a1(:);a2(:)]);
fprintf('created file %s\n',fileout);