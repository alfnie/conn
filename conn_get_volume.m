function [x,idx]=conn_get_volume(V,allowlinks)
if nargin<2||isempty(allowlinks), allowlinks=false; end

if any(conn_server('util_isremotefile',V.fname)), 
    V.fname=conn_server('util_localfile',V.fname);
    if ischar(allowlinks)&&~isempty(regexp(allowlinks,'^.*run_keepas:')), [x,idx]=conn_server('run_keepas',regexprep(allowlinks,'^.*run_keepas:',''),mfilename,V);
    elseif allowlinks, [x,idx]=conn_server('run_keep',mfilename,V);
    else [x,idx]=conn_server('run',mfilename,V);
    end
    return
end
V.fname=conn_server('util_localfile',V.fname);

if isfield(V,'softlink')&&~isempty(V.softlink), 
    str1=regexp(V.fname,'Subject\d+','match'); if ~isempty(str1), V.softlink=regexprep(V.softlink,'Subject\d+',str1{end}); end
    [file_path,file_name,file_ext]=fileparts(V.fname);
    matcfilename=fullfile(file_path,V.softlink); 
else
    matcfilename=[V.fname,'c'];
end
handle=fopen(matcfilename,'rb');
x=fread(handle,V.size.Nt*sum(V.size.Nv),'float');
x=V.scale*reshape(x,[V.size.Nt,sum(V.size.Nv)]);
fclose(handle);
idx=V.voxels;



