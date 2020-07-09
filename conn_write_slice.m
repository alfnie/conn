function conn_write_slice(V,x,slice,time0)

if nargin<4, time0=1; end
if isempty(x), return; end
if time0-1+size(x,1)>V.size.Nt || size(x,2)~=V.size.Nv(slice), error('mismatch dimensions'); end;
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
    ok=fseek(handle,4*(sum(V.size.Nv(1:slice-1))*V.size.Nt+time0-1),-1);
    if ok<0 error('error while writing to file'); end
    if size(x,1)~=V.size.Nt,
        for n1=1:size(x,2), fwrite(handle,x(:,n1),'float'); fseek(handle,4*(V.size.Nt-size(x,1)),0); end
    else
        fwrite(handle,x(:),'float');
    end
    fclose(handle);
end
