function varargout=conn_fileutils(option,varargin)
% conn_fileutils(option,...) performs basic fileio operations on files/directories (system-independent, client/server-independent)
%
% GENERAL FILEIO FUNCTIONS:
% conn_fileutils('fileread', filename)              : reads text file
% conn_fileutils('filewrite', filename, str)        : creates text file with cell array str (one line per element)
% conn_fileutils('fileappend', filename, str)       : appends text file with cell array str (one line per element)
% conn_fileutils('textread', filename)              : see "help textread"
% conn_fileutils('copyfile', source, target)        : copies file
% conn_fileutils('movefile', source, target)        : moves file
% conn_fileutils('renamefile', source, target)      : renames file
% conn_fileutils('deletefile', filename)            : deletes file
% conn_fileutils('emptyfile', filename)             : creates empty file
% conn_fileutils('mkdir', dirname)                  : creates new folder
% conn_fileutils('rmdir', dirname)                  : removes folder
% conn_fileutils('linkdir', dirname, link)          : soft-link folders
% conn_fileutils('dir', name)                       : lists files in folder
% conn_fileutils('isdir', dirname)                  : returns true if dirname points to an existing directory
% conn_fileutils('cd', dirname)                     : changes current working directory
% conn_fileutils('homedir')                         : returns user-specific home directory
% conn_fileutils('imread',...)                      : see "help imread"
% conn_fileutils('imwrite',...)                     : see "help imwrite"
% conn_fileutils('imopen', filename)                : open image in system viewer
% conn_fileutils('uigetfile', ...)                  : see "help uigetfile"
% conn_fileutils('uiputfile', ...)                  : see "help uiputfile"
%
% SPM-SPECIFIC functions: overloaded SPM functions to support /CONNSERVER/[filepath] nomenclature (see below)
% conn_fileutils('nifti',...)
% conn_fileutils('spm_vol',...) 
% conn_fileutils('spm_localvol',...)                : like spm_vol but it creates a cache copy of the file first and returns a handle to the local cache copy
% conn_fileutils('spm_read_vols',...)
% conn_fileutils('spm_get_data',...)
% conn_fileutils('spm_sample_vol',...)
% conn_fileutils('spm_write_vol',...)
% conn_fileutils('spm_unlink',...)                      
% conn_fileutils('spm_file_merge',...)
% conn_fileutils('spm_jsonwrite',...)
%
% note: CONN_FILEUTILS accepts /CONNSERVER/[filepath] nomenclature for remote files or folders in conn_server 
% machine (but when using multiple files/directories (e.g. filecopy) all files must be in the same machine; 
% see CONN_REMOTELY CONN_SERVER and CONN_CACHE for options for transfering files among different machines)
%

varargout=cell(1,nargout);
switch(lower(option))
    case 'file',
        if any(conn_server('util_isremotefile',varargin{1})), varargout={conn_cache('pull',varargin{1})};
        else varargout={conn_server('util_localfile',varargin{1})};
        end
        
    case {'fileread','readfile'}
        if any(conn_server('util_isremotefile',varargin{1})), varargout={conn_server('run',mfilename,option,conn_server('util_localfile',varargin{1}),varargin{2:end})};
        else varargout={fileread(conn_server('util_localfile',varargin{1}),varargin{2:end})};
        end
        
    case {'textread','readtext'}
        if any(conn_server('util_isremotefile',varargin{1})), varargout={conn_server('run',mfilename,option,conn_server('util_localfile',varargin{1}),varargin{2:end})};
        else varargout={textread(conn_server('util_localfile',varargin{1}),varargin{2:end})};
        end
        
    case {'filewrite','writefile','filewrite_raw','writefile_raw','fileappend','appendfile','fileappend_raw','appendfile_raw'}
        if any(conn_server('util_isremotefile',varargin{1})), conn_server('run',mfilename,option,conn_server('util_localfile',varargin{1}),varargin{2:end});
        else
            if isempty(regexp(lower(option),'_raw$')), strcmd='%s\n'; else strcmd='%s'; end
            if isempty(regexp(lower(option),'append')), filecmd='wt'; else filecmd='at'; end
            filename=conn_server('util_localfile',varargin{1});
            str=varargin{2};
            fh=fopen(filename,filecmd);
            if isequal(fh,-1), error('unable to open file %s',filename); end
            if ischar(str), 
                fprintf(fh,'%s',str);
            else 
                for n=1:numel(str),
                    fprintf(fh,strcmd,str{n});
                end
            end
            fclose(fh);
        end
        
    case {'filecopy','copyfile'}
        if any(conn_server('util_isremotefile',varargin{1})), 
            conn_server('run',mfilename,option,conn_server('util_localfile',varargin{1}),conn_server('util_localfile',varargin{2}));
        elseif iscell(varargin{1}), cellfun(@(a,b)conn_fileutils(option,a,b),varargin{1},varargin{2},'uni',0);
        else
            varargin(1:2)=conn_server('util_localfile',varargin(1:2));
            if ispc, [ok,nill]=system(['copy "',varargin{1},'" "',varargin{2},'"']);
            else, [ok,nill]=system(['''cp'' -f ''',varargin{1},''' ''',varargin{2},'''']);
            end
            if ~isequal(ok,0), error('Error copying file %s to %s, missing file or invalid file permissions',varargin{1},varargin{2}); end
        end
        
    case {'filelink','linkfile'}
        if any(conn_server('util_isremotefile',varargin{1})), 
            conn_server('run',mfilename,option,conn_server('util_localfile',varargin{1}),conn_server('util_localfile',varargin{2}));
        elseif iscell(varargin{1}), cellfun(@(a,b)conn_fileutils(option,a,b),varargin{1},varargin{2},'uni',0);
        else
            varargin(1:2)=conn_server('util_localfile',varargin(1:2));
            [ok,nill]=system(['ln -fs ''',varargin{1},''' ''',varargin{2},'''']);
            if ~isequal(ok,0), error('Error linking file %s to %s, missing file or invalid file permissions',varargin{1},varargin{2}); end
        end
        
    case {'dirlink','linkdir'}
        if any(conn_server('util_isremotefile',varargin{1})), 
            conn_server('run',mfilename,option,conn_server('util_localfile',varargin{1}),conn_server('util_localfile',varargin{2}));
        elseif iscell(varargin{1}), cellfun(@(a,b)conn_fileutils(option,a,b),varargin{1},varargin{2},'uni',0);
        else
            varargin(1:2)=conn_server('util_localfile',varargin(1:2));
            if ispc, [ok,nill]=system(['mklink /d "',varargin{2},'" "',varargin{1},'"']);
            else [ok,nill]=system(['ln -fs ''',varargin{1},''' ''',varargin{2},'''']);
            end
            if ~isequal(ok,0), error('Error linking file %s to %s, missing file or invalid file permissions',varargin{1},varargin{2}); end
        end
        
    case {'filemove','movefile'}
        if any(conn_server('util_isremotefile',varargin{1})), 
            conn_server('run',mfilename,option,conn_server('util_localfile',varargin{1}),conn_server('util_localfile',varargin{2}));
        elseif iscell(varargin{1}), cellfun(@(a,b)conn_fileutils(option,a,b),varargin{1},varargin{2},'uni',0);
        else
            varargin(1:2)=conn_server('util_localfile',varargin(1:2));
            if ispc, [ok,nill]=system(['move "',varargin{1},'" "',varargin{2},'"']);
            else, [ok,nill]=system(['mv -f ''',varargin{1},''' ''',varargin{2},'''']);
            end
            if ~isequal(ok,0), error('Error moving file %s to %s, missing file or invalid file permissions',varargin{1},varargin{2}); end
        end
        
    case {'filerename','renamefile'}
        if any(conn_server('util_isremotefile',varargin{1})), 
            conn_server('run',mfilename,option,conn_server('util_localfile',varargin{1}),conn_server('util_localfile',varargin{2}));
        elseif iscell(varargin{1}), cellfun(@(a,b)conn_fileutils(option,a,b),varargin{1},varargin{2},'uni',0);
        else
            varargin(1:2)=conn_server('util_localfile',varargin(1:2));
            if ispc, [ok,nill]=system(['ren "',varargin{1},'" "',varargin{2},'"']);
            else, [ok,nill]=system(['mv -f ''',varargin{1},''' ''',varargin{2},'''']);
            end
            if ~isequal(ok,0), error('Error renaming file %s to %s, missing file or invalid file permissions',varargin{1},varargin{2}); end
        end
        
    case {'filedelete','deletefile','filedelete_multiple','deletefile_multiple'}
        if any(conn_server('util_isremotefile',varargin{1})), 
            conn_server('run',mfilename,option,conn_server('util_localfile',varargin{1}));
        elseif iscell(varargin{1}), cellfun(@(x)conn_fileutils(option,x),varargin{1},'uni',0);
        else
            varargin{1}=conn_server('util_localfile',varargin{1});
            if ~isempty(regexp(option,'_multiple$'))
                if ispc, [ok,nill]=system(sprintf('del "%s"*',varargin{1}));
                else [ok,nill]=system(sprintf('rm -f ''%s''*',varargin{1}));
                end
            else
                if ispc, [ok,nill]=system(sprintf('del "%s"',varargin{1}));
                else [ok,nill]=system(sprintf('rm -f ''%s''',varargin{1}));
                end
            end
            if ~isequal(ok,0), error('Error deleting file %s, missing file or invalid file permissions',varargin{1}); end
        end
        
    case {'fileempty','emptyfile'}
        if any(conn_server('util_isremotefile',varargin{1})), 
            conn_server('run',mfilename,option,conn_server('util_localfile',varargin{1}));
        elseif iscell(varargin{1}), cellfun(@(x)conn_fileutils(option,x),varargin{1},'uni',0);
        else
            varargin{1}=conn_server('util_localfile',varargin{1});
            fh=fopen(varargin{1},'wb');
            if isequal(fh,-1), error('Error creating file %s, check file permissions',varargin{1}); end
            fclose(fh);
        end
        
    case 'mkdir' % note: no error if fails
        if any(conn_server('util_isremotefile',varargin{1})), 
            [ok,msg]=conn_server('run',mfilename,option,conn_server('util_localfile',varargin{1}),varargin{2:end});
        else
            varargin{1}=conn_server('util_localfile',varargin{1});
            [ok,msg]=mkdir(varargin{:});
        end
        if nargout>=1, varargout{1}=ok; end
        if nargout>=2, varargout{2}=msg; end
        
    case {'rmdir','rmdir_dironly','rmdir_recursive'} % note: no error if fails
        if any(conn_server('util_isremotefile',varargin{1})), 
            [ok,msg]=conn_server('run',mfilename,option,conn_server('util_localfile',varargin{1}));
        elseif iscell(varargin{1}), [ok,msg]=cellfun(@(x)conn_fileutils(option,x),varargin{1},'uni',0);
        else
            varargin{1}=conn_server('util_localfile',varargin{1});
            if isdir(varargin{1})
                if strcmpi(option,'rmdir_recursive'), [ok,msg]=rmdir(varargin{1},'s');
                elseif strcmpi(option,'rmdir_dironly'), [ok,msg]=rmdir(varargin{1});
                elseif ispc,
                    [ok,msg]=system(sprintf('del /Q "%s"',fullfile(varargin{1},'*')));
                    [ok,msg]=system(sprintf('rmdir "%s"',varargin{1})); ok=isequal(ok,0); % note: rmdir success output ~= system output
                else
                    [ok,msg]=system(sprintf('rm -f ''%s''/*',varargin{1}));
                    [ok,msg]=system(sprintf('rmdir ''%s''',varargin{1})); ok=isequal(ok,0);
                end
            else, ok=false; msg='directory not found'; 
            end
        end
        if nargout>=1, varargout{1}=ok; end 
        if nargout>=2, varargout{2}=msg; end
        
    case 'isdir' 
        if any(conn_server('util_isremotefile',varargin{1})), 
            out=conn_server('run',mfilename,option,conn_server('util_localfile',varargin{1}));
        else
            out=isdir(conn_server('util_localfile',varargin{1}));
        end
        varargout{1}=out;
        
    case 'dir'
        if numel(varargin)<1||isempty(varargin{1}), filename='.'; else filename=varargin{1}; end
        if numel(varargin)<2||isempty(varargin{2}), option2=''; else option2=varargin{2}; end
        
        if any(conn_server('util_isremotefile',filename)), out=conn_server('run',mfilename,option, conn_server('util_localfile',filename),option2);
        elseif isequal(option2,'-ls'), out=ls('-al',conn_server('util_localfile',filename));
        else out=dir(conn_server('util_localfile',filename));
        end
        if ~nargout,
            if isequal(option2,'-ls'), disp(out);
            elseif isstruct(out) disp(char({out.name}));
            end
            varargout={};
        else varargout{1}=out;
        end
        
    case 'cd'
        if any(conn_server('util_isremotefile',varargin{1})), 
            conn_server('run',mfilename,option,conn_server('util_localfile',varargin{1}));
        else cd(conn_server('util_localfile',varargin{1}));
        end
                
    case 'homedir'
        if ispc, varargout{1}=conn_fullfile(getenv('USERPROFILE'));
        else varargout{1}=conn_fullfile('~/');
        end
        
    case 'java.io.file'
        if any(conn_server('util_isremotefile',varargin{1})), [varargout{1:nargout}]=conn_server('run',mfilename,option,conn_server('util_localfile',varargin{1}),varargin{2:end});
        else [varargout{1:nargout}]=java.io.File(conn_server('util_localfile',varargin{1}),varargin{2:end});
        end
        
    case 'getdiskspace'
        if any(conn_server('util_isremotefile',varargin{1})), [varargout{1:nargout}]=conn_server('run',mfilename,option,conn_server('util_localfile',varargin{1}),varargin{2:end});
        else
            a=java.io.File(conn_server('util_localfile',varargin{1}),varargin{2:end});
            varargout={struct('getUsableSpace',a.getUsableSpace,'getTotalSpace',a.getTotalSpace,'canWrite',a.canWrite)};
        end

    case 'imread',
        if any(conn_server('util_isremotefile',varargin{1})), [varargout{1:nargout}]=imread(conn_cache('pull',varargin{1}),varargin{2:end});
        else [varargout{1:nargout}]=imread(conn_server('util_localfile',varargin{1}),varargin{2:end});
        end
        
    case 'imwrite',
        if ischar(varargin{2})
            if any(conn_server('util_isremotefile',varargin{2})), [varargout{1:nargout}]=conn_server('run',mfilename,option,varargin{1},conn_server('util_localfile',varargin{2}),varargin{3:end});
            else [varargout{1:nargout}]=imwrite(varargin{1},conn_server('util_localfile',varargin{2}),varargin{3:end});
            end
        elseif ischar(varargin{3})
            if any(conn_server('util_isremotefile',varargin{3})), [varargout{1:nargout}]=conn_server('run',mfilename,option,varargin{1:2},conn_server('util_localfile',varargin{3}),varargin{4:end});
            else [varargout{1:nargout}]=imwrite(varargin{1:2},conn_server('util_localfile',varargin{3}),varargin{4:end});
            end
        else error('unsupported imwrite syntax');
        end
        
    case 'imopen',
        if any(conn_server('util_isremotefile',varargin{1})), [ok,msg]=system(sprintf('open %s',conn_cache('pull',varargin{1}),varargin{2:end}));
        else [ok,msg]=system(sprintf('open %s',conn_server('util_localfile',varargin{1})));
        end
        
    case 'uigetfile'
        if conn_projectmanager('inserver')
            if 1 % type-in remote or local file
                if numel(varargin)<1||isempty(varargin{1}), varargin{1}=''; end
                if numel(varargin)<2||isempty(varargin{2}), varargin{2}=''; end
                if numel(varargin)<3||isempty(varargin{3}), varargin{3}='/CONNSERVER/';
                else varargin{3}=conn_server('util_remotefile',varargin{3});
                end           
                filename=conn_fileutils_uifile('get',varargin{:});
                if ~isempty(filename),
                    [tfilepath,t1,t2]=fileparts(filename);
                    varargout={[t1,t2],tfilepath};
                    %new_tpathname=conn_projectmanager('homedir');
                    %if any(conn_server('util_isremotefile',filename))
                    %    varargout={[t1,t2],tfilepath};
                    %else % if selected local file push to remote
                    %    new_tfilename=[t1,'_',char(conn_tcpip('hash',[t1, mat2str(now)])),t2];
                    %    conn_server('push',conn_server('util_localfile',filename),fullfile(new_tpathname,new_tfilename));
                    %    varargout={new_tfilename,new_tpathname};
                    %end
                else varargout={0,0};
                end                    
            else % select local file and push to remote
                [tfilename,tpathname]=uigetfile(varargin{1},varargin{2},conn_server('util_localfile',varargin{3}),varargin{4:end});
                varargout={tfilename,tpathname};
                if ischar(tfilename),
                    [nill,t1,t2]=fileparts(tfilename);
                    new_tfilename=[t1,'_',char(conn_tcpip('hash',[t1, mat2str(now)])),t2];
                    new_tpathname=CONN_x.folders.data;
                    conn_server('push',fullfile(tpathname,tfilename),fullfile(new_tpathname,new_tfilename));
                    varargout={new_tfilename,new_tpathname};
                end
            end
        else [varargout{1:nargout}]=uigetfile(varargin{:});
        end        
        
    case 'uiputfile'
        if conn_projectmanager('inserver') % type-in remote or local file
            if numel(varargin)<1||isempty(varargin{1}), varargin{1}=''; end
            if numel(varargin)<2||isempty(varargin{2}), varargin{2}=''; end
            if numel(varargin)<3||isempty(varargin{3}), varargin{3}='/CONNSERVER/';
            else varargin{3}=conn_server('util_remotefile',varargin{3});
            end
            [filename,checked]=conn_fileutils_uifile('put',varargin{:});
            if ~isempty(filename),
                [tfilepath,t1,t2]=fileparts(filename);
                if ~checked&&conn_existfile(filename), 
                    tansw=conn_questdlg({'Output file already exist','Overwrite existing file?'},'','Yes','No','Yes');
                    if ~strcmp(tansw,'Yes'), varargout={0,0}; return; end
                end
                varargout={[t1,t2],tfilepath};
            else varargout={0,0};
            end
        else [varargout{1:nargout}]=uiputfile(varargin{:});
        end
        
    case 'spm_unlink'
        files=conn_server('util_localfile',varargin);
        if any(conn_server('util_isremotefile',varargin)), 
            conn_server('run',mfilename,option,files{:});
        else spm_unlink(files{:});
        end
        
    case 'spm_jsonwrite'
        if any(conn_server('util_isremotefile',varargin{1})), [varargout{1:nargout}]=conn_server('run',mfilename,option,conn_server('util_localfile',varargin{1}),varargin{2:end});
        else [varargout{1:nargout}]=spm_jsonwrite(conn_server('util_localfile',varargin{1}),varargin{2:end});
        end
                
    case 'nifti'
        if any(conn_server('util_isremotefile',varargin{1})), varargout{1}=conn_server('run',mfilename,option,conn_server('util_localfile',varargin{1}),varargin{2:end});
        else varargout{1}=struct(nifti(conn_server('util_localfile',varargin{1}),varargin{2:end}));
        end
        
    case 'nifti_nvol'
        if any(conn_server('util_isremotefile',varargin{1})), varargout{1}=conn_server('run',mfilename,option,conn_server('util_localfile',varargin{1}),varargin{2:end});
        else
            nfilename=nifti(conn_server('util_localfile',varargin{1}),varargin{2:end});
            nV=0; for n=1:numel(nfilename), tV=size(nfilename(n).dat,4); nV=nV+tV; end
            varargout={nV};
        end
        
    case 'spm_vol'
        if any(conn_server('util_isremotefile',varargin{1})), 
            V=conn_server('run',mfilename,option,conn_server('util_localfile',varargin{1}),varargin{2:end});
            if iscell(V), for n=1:numel(V), V{n}.fname=conn_server('util_remotefile',V{n}.fname); end
            else for n=1:numel(V), V(n).fname=conn_server('util_remotefile',V(n).fname); end
            end
        else V=spm_vol(conn_server('util_localfile',varargin{1}),varargin{2:end});
        end
        varargout{1}=V;
        
    case 'spm_localvol' % pull file to local cache and returns spm_vol of local file
        if isstruct(varargin{1}), % convert struct to local
            if isstruct(varargin{1})&&isfield(varargin{1},'fname')&&any(conn_server('util_isremotefile',varargin{1}(1).fname)), vol=spm_vol(char(conn_cache('pull',conn_expandframe(varargin{1}))));
            else vol=varargin{1};
            end
        elseif any(conn_server('util_isremotefile',varargin{1})), vol=spm_vol(char(conn_cache('pull',varargin{1})));
        else vol=spm_vol(char(conn_server('util_localfile',varargin{1})));
        end
        varargout{1}=vol;
        
    case 'spm_read_vols'
        DOCACHE=true;
        if DOCACHE
            vol=conn_fileutils('spm_localvol',varargin{1});
            data=spm_read_vols(vol,varargin{2:end});
        else
            if ischar(varargin{1})
                if any(conn_server('util_isremotefile',varargin{1})), data=conn_server('run',mfilename,option,conn_server('util_localfile',varargin{1}),varargin{2:end});
                else data=spm_read_vols(spm_vol(conn_server('util_localfile',varargin{1})),varargin{2:end});
                end
            elseif isstruct(varargin{1})&&isfield(varargin{1},'fname')
                V=varargin{1};
                if iscell(V), for n=1:numel(V), V{n}.fname=conn_server('util_localfile',V{n}.fname); end
                else for n=1:numel(V), V(n).fname=conn_server('util_localfile',V(n).fname); end
                end
                if any(conn_server('util_isremotefile',varargin{1}(1).fname)),
                    data=conn_server('run',mfilename,option,V,varargin{2:end});
                else 
                    data=spm_read_vols(V,varargin{2:end});
                end
            else data=spm_read_vols(varargin{:});
            end
        end
        varargout{1}=data;
        
    case 'spm_get_data'
        DOCACHE=false;
        try, DOCACHE=size(varargin{2},2)>.10*sum(arrayfun(@(n)prod(size(varargin{1}(n).private.dat)),1:numel(varargin{1}))); end
        if DOCACHE
            vol=conn_fileutils('spm_localvol',varargin{1});
            data=spm_get_data(vol,varargin{2:end});
        else
            if isstruct(varargin{1})&&isfield(varargin{1},'fname')
                V=varargin{1};
                if iscell(V), for n=1:numel(V), V{n}.fname=conn_server('util_localfile',V{n}.fname); end
                else for n=1:numel(V), V(n).fname=conn_server('util_localfile',V(n).fname); end
                end
                if any(conn_server('util_isremotefile',varargin{1}(1).fname)), data=conn_server('run',mfilename,option,V,varargin{2:end});
                else data=spm_get_data(V,varargin{2:end});
                end
            else data=spm_get_data(varargin{:});
            end
        end
        varargout{1}=data;
        
    case 'spm_sample_vol'
        DOCACHE=false;
        try, DOCACHE=numel(varargin{2})>.10*sum(arrayfun(@(n)prod(size(varargin{1}(n).private.dat)),1:numel(varargin{1}))); end
        if DOCACHE
            vol=conn_fileutils('spm_localvol',varargin{1});
            [varargout{1:nargout}]=spm_sample_vol(vol,varargin{2:end});
        else
            if isstruct(varargin{1})&&isfield(varargin{1},'fname')
                V=varargin{1};
                if iscell(V), for n=1:numel(V), V{n}.fname=conn_server('util_localfile',V{n}.fname); end
                else for n=1:numel(V), V(n).fname=conn_server('util_localfile',V(n).fname); end
                end
                if any(conn_server('util_isremotefile',varargin{1}(1).fname)), [varargout{1:nargout}]=conn_server('run',mfilename,option,V,varargin{2:end});
                else [varargout{1:nargout}]=spm_sample_vol(V,varargin{2:end});
                end
            else [varargout{1:nargout}]=spm_sample_vol(varargin{:});
            end
        end
        
    case 'spm_write_vol'
        if isstruct(varargin{1})&&isfield(varargin{1},'fname')
            V=varargin{1};
            if iscell(V), for n=1:numel(V), V{n}.fname=conn_server('util_localfile',V{n}.fname); end
            else for n=1:numel(V), V(n).fname=conn_server('util_localfile',V(n).fname); end
            end
            if any(conn_server('util_isremotefile',varargin{1}(1).fname)), 
                V=conn_server('run',mfilename,option,V,varargin{2:end});
                if iscell(V), for n=1:numel(V), V{n}.fname=conn_server('util_remotefile',V{n}.fname); end
                else for n=1:numel(V), V(n).fname=conn_server('util_remotefile',V(n).fname); end
                end
            else V=spm_write_vol(V,varargin{2:end});
            end
        else V=spm_write_vol(varargin{:});
        end
        if nargout>=1, varargout{1}=V; end
        
    case 'spm_file_merge'
        if isstruct(varargin{1})&&isfield(varargin{1},'fname')
            V=varargin{1};
            if iscell(V), for n=1:numel(V), V{n}.fname=conn_server('util_localfile',V{n}.fname); end
            else for n=1:numel(V), V(n).fname=conn_server('util_localfile',V(n).fname); end
            end
            if any(conn_server('util_isremotefile',varargin{1}(1).fname)), conn_server('run',mfilename,option,V,conn_server('util_localfile',varargin{2}),varargin{3:end});
            else spm_file_merge(V,conn_server('util_localfile',varargin{2}),varargin{3:end});
            end
        elseif ischar(varargin{1})
            if any(conn_server('util_isremotefile',varargin{1})), conn_server('run',mfilename,option,conn_server('util_localfile',varargin{1}),conn_server('util_localfile',varargin{2}),varargin{3:end});
            else spm_file_merge(conn_server('util_localfile',varargin{1}),conn_server('util_localfile',varargin{2}),varargin{3:end});
            end
        else spm_file_merge(varargin{:});
        end
            
    otherwise
        error('unknown option %s',option);
end
if ~nargout, varargout={}; end
end


function [filename,skipoverwritequestion]=conn_fileutils_uifile(style,varargin)
if nargin<1||isempty(style), style='get'; end
if numel(varargin)<1||isempty(varargin{1}), varargin{1}='*'; end
if numel(varargin)<2||isempty(varargin{2}), if strcmp(style,'put'), varargin{2}='Select a file to write'; else varargin{2}='Select a file'; end; end
if numel(varargin)<3||isempty(varargin{3}), varargin{3}=''; end
options=varargin;

filename='';
skipoverwritequestion=strcmp(style,'get');
thfig=figure('units','norm','position',[.4,.5,.35,.15],'color',1*[1 1 1],'name','file dialog','numbertitle','off','menubar','none');
ht1a=uicontrol('style','text','units','norm','position',[.1,.8,.8,.15],'string',sprintf('%s (%s)',varargin{2},varargin{1}),'horizontalalignment','left','backgroundcolor',1*[1 1 1],'fontweight','bold');
ht1=uicontrol('style','edit','units','norm','position',[.1,.65,.7,.15],'string',varargin{3},'tooltipstring','<HTML>enter full path to target file<br/> - click the button at the right to do this using the system dialog box (only for local files) <br/> - use /CONNSERVER/[remote_filepath] syntax to refer to files in remote server</HTML>','userdata',skipoverwritequestion,'callback',['set(gcbo,''userdata'',',num2str(skipoverwritequestion),');']);
ht1b=uicontrol('style','pushbutton','units','norm','position',[.8,.65,.1,.15],'string','...','backgroundcolor',1*[1 1 1],'tooltipstring','select target file using system dialog box (only for local files)','callback',@conn_uifile_callback);
%ht2a=uicontrol('style','text','units','norm','position',[.1,.45,.8,.15],'string','(use /CONNSERVER/* for remote files)','horizontalalignment','center','backgroundcolor',1*[1 1 1]);
uicontrol('style','pushbutton','string','OK','units','norm','position',[.1,.01,.38,.25],'callback','uiresume');
uicontrol('style','pushbutton','string','Cancel','units','norm','position',[.51,.01,.38,.25],'callback','delete(gcbf)');
uicontrol(ht1);
while isempty(filename),
    uiwait(thfig);
    if ~ishandle(thfig), return; end
    filename=get(ht1,'string');
    skipoverwritequestion=isequal(get(ht1,'userdata'),1);
    if ~skipoverwritequestion&&~isempty(filename)&&conn_existfile(filename)
        tansw=conn_questdlg({'Output file already exist','Overwrite existing file?'},'','Yes','No','Yes');
        if ~strcmp(tansw,'Yes'), filename=''; 
        else skipoverwritequestion=true;
        end
    end
end
delete(thfig);
drawnow;

    function conn_uifile_callback(varargin)
        filename=conn_server('util_localfile',get(ht1,'string'));
        if strcmp(style,'put'), [tfilename,tpathname]=uiputfile(options{1:2},filename);
        else [tfilename,tpathname]=uigetfile(options{1:2},filename);
        end
        if ~isequal(tpathname,0), 
            filename=fullfile(tpathname,tfilename); 
            set(ht1,'string',filename,'userdata',1); 
        end
    end
end
