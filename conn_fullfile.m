function filename=conn_fullfile(varargin)

if ~nargin, filename=''; 
elseif nargin==1, filename=varargin{1}; 
else filename=fullfile(varargin{:}); 
end
if isempty(filename), return; end
if iscell(filename), 
    filename=cellfun(@conn_fullfile,filename,'uni',0); 
    return; 
elseif size(filename,1)>1, 
    filename=char(cellfun(@conn_fullfile,cellstr(filename),'uni',0)); 
    return; 
end

[filename_path,filename_name,filename_ext]=fileparts(filename);
if isempty(filename_path),
    filename_path=pwd;
else
    ok=false;
    try, ok=java.io.File(filename_path).isAbsolute; end
    if ~ok
        cwd=pwd;
        try
            cd(filename_path);
            filename_path=pwd;
            cd(cwd);
        end
    end
end
filename=fullfile(filename_path,[filename_name,filename_ext]);
