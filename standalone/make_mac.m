% creates files for standalone release in Mac
path_current=pwd;
path_source=fileparts(which('conn'));
path_target='/Applications/conn_standalone';
[ok,msg]=mkdir(path_target); % NOTE: may need write permissions for this
ver=conn('ver');
cd(fullfile(path_source,'standalone'));
conn_makestandalone
[ok,msg]=system(sprintf('rm -fr ''%s/conn.app''',path_target));
[ok,msg]=system(sprintf('mv -f conn.app ''%s/''',path_target));
[ok,msg]=system(sprintf('mv -f run_conn.sh ''%s/''',path_target));
[ok,msg]=system(sprintf('mv -f readme.txt ''%s/''',path_target));
[ok,msg]=system(sprintf('cp -f installation_mac.txt ''%s/''',path_target));
[ok,msg]=system('rm -f mccExcludedFiles.log');
[ok,msg]=system('rm -f requiredMCRProducts.txt');
cd(path_target)
[ok,msg]=system('rm -f conn*.zip');
fileout=sprintf('conn%s_maci64.zip',regexprep(ver,'\.',''));
[ok,msg]=system(sprintf('zip -r %s *',fileout));
fprintf('Finished compilation. Created file %s\n',fullfile(path_target,fileout));
cd(path_current);
