function fh=evlab17_resultsplots(filename,varargin)
% EVLAB17_RESULTSPLOTS displays second-level analysis results
%   evlab17_resultsplots('/myfolder/SPM.mat');
%      displays second-level results stored in /myfolder
%   evlab17_resultsplots('secondlevel.cfg');
%      displays second-level results defined in secondlevel.cfg file (cfg file needs to contain #folder field specifying target folder)
%

evlab17_module init silent;

% interprets info
cwd=pwd;
if conn_existfile(filename,2), filename=fullfile(filename,'SPM.mat'); end
assert(conn_existfile(filename),'file %s not found',filename);
[nill,nill,ext]=spm_fileparts(filename);
if strcmp(ext,'.cfg')
    options=conn_loadcfgfile(filename);
    assert(isfield(options,'folder'),'missing #folder information in %s',filename);
    filename=fullfile(char(options.folder),'SPM.mat');
end
fh=conn_display(filename,varargin{:});
cd(cwd);


end



