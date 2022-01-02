function [x,idx]=conn_get_time(V,time,slice,reshapeoutput)

if nargin<2 || isempty(time), time=1; end
if nargin<3, slice=[]; end
if nargin<4, reshapeoutput=true; end
if any(conn_server('util_isremotefile',V.fname)), V.fname=conn_server('util_localfile',V.fname); [x,idx]=conn_server('run',mfilename,V,time,slice,reshapeoutput); return; end
V.fname=conn_server('util_localfile',V.fname);

if isfield(V,'softlink')&&~isempty(V.softlink), 
    str1=regexp(V.fname,'Subject\d+','match'); if ~isempty(str1), V.softlink=regexprep(V.softlink,'Subject\d+',str1{end}); end
    [file_path,file_name,file_ext]=fileparts(V.fname);
    matcfilename=fullfile(file_path,V.softlink); 
else
    matcfilename=[V.fname,'c'];
end

handle=fopen(matcfilename,'rb');
if isempty(slice), 
	fseek(handle,4*(time-1),-1);
	x=fread(handle,sum(V.size.Nv),'float',4*(V.size.Nt-1));
	idx=V.voxels;
else
	fseek(handle,4*(sum(V.size.Nv(1:slice-1))*V.size.Nt + time-1),-1);
	x=fread(handle,V.size.Nv(slice),'float',4*(V.size.Nt-1));
	idx=1+rem(V.voxels(sum(V.size.Nv(1:slice-1))+(1:V.size.Nv(slice)))-1,prod(V.matdim.dim(1:2)));
end
x=V.scale*x;
if reshapeoutput
    if isempty(slice), y=zeros(V.matdim.dim); 
    else y=zeros(V.matdim.dim(1:2)); 
    end
    y(idx)=x;
    x=y;
end
fclose(handle);


