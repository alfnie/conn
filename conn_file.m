function [out,nV,filename]=conn_file(filename,doconvert)
if nargin<2||isempty(doconvert), doconvert=true; end
if isempty(filename)
    nV=[];
    out={[],[],[]};
elseif ~ischar(filename)&&~iscell(filename),
    nV=[];
    out={'[raw values]',[],filename};
else 
    filename=char(filename);
    [nV,str,icon,filename]=conn_getinfo(filename,doconvert);
    out={fliplr(deblank(fliplr(deblank(conn_fullfile(filename))))),str,icon};
end

