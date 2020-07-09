% This script creates a stand-alone version of CONN
% 
% This requires Matlab compiler properly installed and configured in your system
%

spmpath=fileparts(which('spm'));
connpath=fileparts(which('conn'));
if ispc, [ok,nill]=system('del conn.exe');
elseif ismac, [ok,nill]=system('rm -rf conn.app');
else     [ok,nill]=system('rm conn'); 
end

% miscelaenous SPM standalone issues
spm_jobman('initcfg');
fid = fopen(fullfile(spm('dir'),'config','spm_cfg_static_tools.m'),'wt');
fprintf(fid,'function values = spm_cfg_static_tools\n');
fprintf(fid,...
    '%% Static listing of all batch configuration files in the SPM toolbox folder\n');
% create code to insert toolbox config
%-Toolbox autodetection
%-Get the list of toolbox directories
tbxdir = fullfile(spm('Dir'),'toolbox');
d  = dir(tbxdir); d = {d([d.isdir]).name};
dd = regexp(d,'^\.');
%(Beware, regexp returns an array if input cell array is of dim 0 or 1)
if ~iscell(dd), dd = {dd}; end
d  = {'' d{cellfun('isempty',dd)}};
ft = {};
%-Look for '*_cfg_*.m' files in these directories
for i=1:length(d)
    d2 = fullfile(tbxdir,d{i});
    di = dir(d2); di = {di(~[di.isdir]).name};
    f2 = regexp(di,'.*_cfg_.*\.m$');
    if ~iscell(f2), f2 = {f2}; end
    fi = di(~cellfun('isempty',f2));
    if ~isempty(fi)
        ft = [ft(:); fi(:)];
    end
end
if ~isempty(ft)
    if isempty(ft)
        ftstr = '';
    else
        ft = cellfun(@(cft)strtok(cft,'.'),ft,'UniformOutput',false);
        ftstr  = sprintf('%s ', ft{:});
    end
    fprintf(fid,'values = {%s};\n', ftstr);
end
fclose(fid);
cfg_util('dumpcfg');
copyfile(fullfile(spmpath,'Contents.m'),fullfile(spmpath,'Contents.txt'));

% miscelaenous CONN standalone issues
conn_premakestandalone;

% compilation
mcc('-mv','-o','conn','conn.m','-N','-p',fullfile(matlabroot,'toolbox','signal'),'-R','-singleCompThread','-R','-startmsg,Loading MCR. Please wait...','-a',spmpath,'-a',connpath);
if ~ispc, [ok,nill]=system('chmod 755 conn run_conn.sh'); end

