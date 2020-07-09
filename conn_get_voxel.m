function x=conn_get_voxel(V,nv)

if nargin<2, error('insufficient arguments'); end
if isfield(V,'softlink')&&~isempty(V.softlink), 
    str1=regexp(V.fname,'Subject\d+','match'); if ~isempty(str1), V.softlink=regexprep(V.softlink,'Subject\d+',str1{end}); end
    [file_path,file_name,file_ext]=fileparts(V.fname);
    matcfilename=fullfile(file_path,V.softlink); 
else
    matcfilename=[V.fname,'c'];
end
handle=fopen(matcfilename,'r+b');
ok=fseek(handle,4*((nv-1)*V.size.Nt),-1);
if ok<0 error('error while writing to file'); end 
x=fread(handle,V.size.Nt,'float');
fclose(handle);
