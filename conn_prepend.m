function fileout=conn_prepend(str1,file,varargin)
if iscell(file),
    fileout=cell(size(file));
    for n1=1:numel(file),
        fileout{n1}=conn_prepend(str1,file{n1},varargin{:});
    end
elseif ischar(file),
    [fpath,ffile,fext]=fileparts(file);
    if nargin<3, str2=fext; else, str2=varargin{1}; end
    if numel(str1)==1&&~ischar(str1), 
        fileout=fullfile(fpath,[ffile(1-str1:end),str2]);
    else
        fileout=fullfile(fpath,[str1,ffile,str2]);
    end
else,
    fileout=file;
end
