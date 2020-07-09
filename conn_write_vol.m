function [V]=conn_write_vol(V,x)

if ~isequal(size(x),[V.size.Nt sum(V.size.Nv)]), error('mismatch dimensions'); end;
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
    handle=fopen(matcfilename,'wb');
    fwrite(handle,x,'float');
    fclose(handle);
end
