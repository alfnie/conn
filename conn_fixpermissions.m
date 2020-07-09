function conn_fixpermissions(dataset,perms,opt)
% CONN_FIXPERMISSIONS changes file permissions of project files/folders
%
% conn_fixpermissions(folder)
%   sets file permissions to grant group-level write permissions to all files/folders within "folder" (i.e. perm = g+rw)
%
% conn_fixpermissions(folder,perms)
%   lets user specify permissions (see "man chmod")
%   default perms='g+rw'
%
% conn_fixpermissions([],perms)
%   changes default permissions (in following uses of conn_fixpermissions within current Matlab session)
%

persistent bak_perms;

if isempty(bak_perms), bak_perms='g+rw'; end
if nargin<3||isempty(opt), opt=false; end
if nargin<2||isempty(perms), perms=bak_perms; end
if nargin<1, return; end
if isempty(dataset), bak_perms=perms; return; end

if opt
    disp(['changing permissions in ',conn_prepend('',dataset,'')]);
    [ok,msg]=system(sprintf('chmod -R %s ''%s''',perms,conn_prepend('',dataset,'')));
    [ok,msg]=system(sprintf('chmod -R %s ''%s''',perms,conn_prepend('',dataset,'.mat')));
    [ok,msg]=system(sprintf('chmod -R %s ''%s''',perms,conn_prepend('',dataset,'.dmat')));
    [ok,msg]=system(sprintf('chmod -R %s ''%s''',perms,conn_prepend('',dataset,'.emat')));
    [ok,msg]=system(sprintf('chmod -R %s ''%s''',perms,conn_prepend('',dataset,'.qlog')));
else
    disp(['changing permissions in ',dataset]);
    [ok,msg]=system(sprintf('chmod -R %s ''%s''',perms,dataset));
end
