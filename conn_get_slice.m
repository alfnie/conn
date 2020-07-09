function [x,idx]=conn_get_slice(V,slice,ntime)
if nargin<3, ntime=V.size.Nt; 
else ntime=max(1,min(V.size.Nt,ntime)); 
end
if nargin<2, slice=1; end
if isfield(V,'softlink')&&~isempty(V.softlink), 
    str1=regexp(V.fname,'Subject\d+','match'); if ~isempty(str1), V.softlink=regexprep(V.softlink,'Subject\d+',str1{end}); end
    [file_path,file_name,file_ext]=fileparts(V.fname);
    matcfilename=fullfile(file_path,V.softlink); 
else
    matcfilename=[V.fname,'c'];
end

for nrepeat=1:4
    handle=fopen(matcfilename,'rb');
    fseek(handle,4*sum(V.size.Nv(1:slice-1))*V.size.Nt,-1);
    if ntime<V.size.Nt
        x=fread(handle,V.size.Nv(slice)*ntime,[num2str(ntime),'*float'],4*(V.size.Nt-ntime));
    else
        x=fread(handle,V.size.Nv(slice)*V.size.Nt,'float');
    end
    if numel(x)==ntime*V.size.Nv(slice), break; 
    else disp('I/O error. retyring...'); 
    end
end
x=V.scale*reshape(x,[ntime,V.size.Nv(slice)]);
fclose(handle);
idx=1+rem(V.voxels(sum(V.size.Nv(1:slice-1))+(1:V.size.Nv(slice)))-1,prod(V.matdim.dim(1:2)));



