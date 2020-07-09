function filename=conn_gz2nii(filename)

ischarfilename=ischar(filename);
filename=cellfun(@strtrim,cellstr(filename),'uni',0);
filenameout=regexprep(filename,'\.gz\s*$|\.zip\s*$','');
redo=~cellfun(@conn_existfile,filenameout);
if ~any(redo), filename=filenameout; if ischarfilename, filename=char(filename); end; return; end

fprintf('unzipping gz files...');
[pathname,name,ext]=spm_fileparts(filename{1});
if strcmp(ext,'.gz')
    filename(redo)=gunzip(filename(redo));
elseif strcmp(ext,'.zip')
    filename(redo)=unzip(filename(redo));
end
filename=filenameout; 
if ischarfilename, filename=char(filename); end; 
fprintf('done\n');

