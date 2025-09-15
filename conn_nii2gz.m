function filename=conn_nii2gz(filename)
% CONN_GZ2NII converts NIFTI files to .nii.gz format
%
% filename = conn_nii2gz(filename)
%
% note: conn_nii2gz will not overwrite output file if it exists
%

if any(conn_server('util_isremotefile',filename)), filename=conn_server('util_remotefile',conn_server('run',mfilename,conn_server('util_localfile',filename))); return; end

filename=conn_server('util_localfile',filename);
ischarfilename=ischar(filename);
filename=cellfun(@strtrim,cellstr(filename),'uni',0);
filenameout=regexprep(regexprep(filename,'\.gz\s*$|\.zip\s*$',''),'.*','$0.gz');
redo=~cellfun(@conn_existfile,filenameout);
if ~any(redo), filename=filenameout; if ischarfilename, filename=char(filename); end; return; end

fprintf('compressing files...');
[pathname,name,ext]=spm_fileparts(filename{1});
filename(redo)=gzip(filename(redo));
filename=filenameout; 
if ischarfilename, filename=char(filename); end; 
fprintf('done\n');

