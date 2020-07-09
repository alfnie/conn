function [ok,files]=conn_checkdistributionfiles(filename,varargin)
if isdeployed, ok=true; files={}; return; end
if ~nargin||isempty(filename), 
    conn_checkdistributionfiles('spm',varargin{:});
    conn_checkdistributionfiles('conn',varargin{:});
    return
end
thispath=which(filename);
if any(strcmpi(varargin(cellfun(@ischar,varargin)),'nolog')), conndisp=@(varargin)fprintf(varargin{:}); else conndisp=@(varargin)conn_disp('fprintf',varargin{:}); end
if isempty(thispath), 
    conndisp('Error: %s not found!\n',filename); 
    switch(filename)
        case 'spm',
            conndisp(...
                ['To install SPM follow these instructions:\n',...
                    '  SPM installation:\n',...
                    '    1) Download the SPM installation spm*.zip file from http://www.fil.ion.ucl.ac.uk/spm and uncompress the contents of this file in your desired target installation folder\n',...
                    '    2) Start Matlab, \n',...
                    '        Click on the "Set path" menu-item or button in the main Matlab command window (in older Matlab versions this item is in the ''File'' menu-list, in newer Matlab versions this button is the ''Home'' tab)\n',...
                    '        Click on "Add folder" (not "Add with subfolders") and select the target installation folder (make sure the selected folder is the one containing the file spm.m)\n',...
                    '        Click on "Save" (this will keep these changes for future Matlab sessions)\n']);
        case 'conn'
            conndisp(...
                ['To install CONN follow these instructions:\n',...
                    '  CONN installation:\n',...
                    '    1) Download the CONN installation conn*.zip file from http://www.nitrc.org/projects/conn and uncompress the contents of this file in your desired target installation folder\n',...
                    '    2) Start Matlab, \n',...
                    '        Click on the "Set path" menu-item or button in the main Matlab command window (in older Matlab versions this item is in the ''File'' menu-list, in newer Matlab versions this button is the ''Home'' tab)\n',...
                    '        Click on "Add folder" (not "Add with subfolders") and select the target installation folder (make sure the selected folder is the one containing the file conn.m)\n',...
                    '        Click on "Save" (this will keep these changes for future Matlab sessions)\n']);
    end
    ok=false;files={};return; 
end
thispath=fileparts(thispath);
if nargin<=1||~any(strcmpi(varargin(cellfun(@ischar,varargin)),'silent')), conndisp('%s @ %s\n',filename,thispath); end
names=dir(fullfile(thispath,'*.m'));
files={names.name};
okpath=cellfun(@(x)strcmpi(fileparts(which(x)),thispath)|strncmp(x,'.',1)|strncmpi(x,'contents.m',1),files);
if ~nargout
    for n=find(~okpath(:)')
        foldername=fileparts(which(files{n}));
        if isempty(foldername)||strcmp(foldername,pwd), conndisp('Warning: %s overloaded by version in current folder (%s). If experiencing unexpected errors please try cd''ing to a different folder before starting CONN, or removing the duplicated file %s from the folder %s\n',files{n},pwd,files{n},pwd);
        else conndisp('Warning: %s overloaded by version in folder %s. If experiencing unexpected errors please try deleting the folder %s from the Matlab path (see ''pathtool'')\n',files{n},foldername,foldername);
        end
    end
end
files=files(~okpath);
ok=all(okpath);
