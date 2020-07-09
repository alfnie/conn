% creates files for standalone release in PC
path_current=pwd;
path_source=fileparts(which('conn'));
path_target='\Program Files\conn_standalone'; 
[ok,msg]=mkdir(path_target); % NOTE: probably need write permissions for this (e.g. icacls "C:\Program Files\conn_standalone" /grant YOURUSERNAME:(OI)(CI)F /T)
ver=conn('ver');
cd(fullfile(path_source,'standalone'));
conn_makestandalone
[ok,msg]=system(sprintf('move conn.exe "%s\"',path_target));
[ok,msg]=system(sprintf('move readme.txt "%s\"',path_target));
[ok,msg]=system(sprintf('copy installation_pc.txt "%s\"',path_target));
[ok,msg]=system('del mccExcludedFiles.log');
[ok,msg]=system('del requiredMCRProducts.txt');
cd(path_target)
[ok,msg]=system('del -i conn*.zip');
fileout=sprintf('conn%s_win64.zip',regexprep(ver,'\.',''));
[ok,msg]=system(sprintf('zip %s *',fileout));
fprintf('Finished compilation. Created file %s\n',fullfile(path_target,fileout));
cd(path_current);
