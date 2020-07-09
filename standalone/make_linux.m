% creates files for standalone release in linux
path_current=pwd;
path_source=fileparts(which('conn'));
path_target=fullfile(fileparts(path_source),'conn_standalone');
[ok,msg]=mkdir(path_target); % NOTE: may need write permissions for this
ver=conn('ver');
cd(fullfile(path_source,'standalone'));
conn_makestandalone
[ok,msg]=system(sprintf('mv -f conn ''%s/''',path_target));
[ok,msg]=system(sprintf('mv -f run_conn.sh ''%s/''',path_target));
[ok,msg]=system(sprintf('mv -f readme.txt ''%s/''',path_target));
[ok,msg]=system(sprintf('cp -f modulefile.txt ''%s/''',path_target));
[ok,msg]=system(sprintf('cp -f installation_linux.txt ''%s/''',path_target));
[ok,msg]=system('rm -f mccExcludedFiles.log');
[ok,msg]=system('rm -f requiredMCRProducts.txt');
cd(path_target)
[ok,msg]=system('rm -f conn*.zip');
fileout=sprintf('conn%s_glnxa64.zip',regexprep(ver,'\.',''));
[ok,msg]=system(sprintf('zip -r %s *',fileout));
fprintf('Finished compilation. Created file %s\n',fullfile(path_target,fileout));
cd(path_current);
