function evlab17_qaplots(varargin)
% EVLAB17_QAPLOTS displays Quality Assurance plots
%   evlab17_qaplots('/myfolder/evlab17_foo.mat'); 
%      displays QA plots associated with dataset evlab17_foo
%   evlab17_qaplots('/myfolder/nii'); 
%      displays all QA plots in /myfolder/nii folder
%

evlab17_module init silent;

% interprets info
if numel(varargin)<1, pathname=pwd;
elseif isequal(varargin{1},'dataset')&&numel(varargin)>1, pathname=varargin{2};
else pathname=varargin{1};
end

cwd=pwd;
if isdir(pathname)
    cd(pathname)
    conn_qaplotsexplore thesefolders;
elseif conn_existfile(pathname)&&~isempty(regexp(pathname,'\.mat$'))
    cd(fileparts(pathname));
    evlab17_module('load',pathname);
    conn_qaplotsexplore thesefolders isconn;
else error('evlab17_qaplots argument must be a directory or a evlab17_*.mat file');
end
cd(cwd);


end

