function conn_write_voxel(V,x,nv)

if nargin<3, error('insufficient arguments'); end
if isempty(x), return; end
if numel(x)~=V.size.Nt, error('mismatch dimensions'); end;
if isfield(V,'softlink')&&~isempty(V.softlink), 
    str1=regexp(V.fname,'Subject\d+','match'); if ~isempty(str1), V.softlink=regexprep(V.softlink,'Subject\d+',str1{end}); end
    [file_path,file_name,file_ext]=fileparts(V.fname);
    matcfilename=fullfile(file_path,V.softlink); 
    overwrite=V.overwritesoftlink;
else
    matcfilename=[V.fname,'c'];
    overwrite=true;
end
if overwrite
    handle=fopen(matcfilename,'r+b');
    ok=fseek(handle,4*((nv-1)*V.size.Nt),-1);
    if ok<0 error('error while writing to file'); end
    fwrite(handle,x(:),'float');
    fclose(handle);
end
