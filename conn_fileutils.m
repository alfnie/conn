function varargout=conn_fileutils(option,varargin)
% conn_fileutils(option,...) performs basic fileio operations on files/directories (system-independent)
%
% GENERAL FILEIO FUNCTIONS:
% conn_fileutils('fileread', filename)              : reads text file
% conn_fileutils('filewrite', filename, str)        : creates text file with cell array str (one line per element)
% conn_fileutils('fileappend', filename, str)       : appends text file with cell array str (one line per element)
% conn_fileutils('copyfile', source, target)        : copies file
% conn_fileutils('movefile', source, target)        : moves file
% conn_fileutils('renamefile', source, target)      : renames file
% conn_fileutils('deletefile', filename)            : deletes file
% conn_fileutils('emptyfile', filename)             : creates empty file
% conn_fileutils('mkdir', dirname)                  : creates new folder
% conn_fileutils('rmdir', dirname)                  : removes folder
% conn_fileutils('dir', name)                       : lists files in folder
% conn_fileutils('isdir', dirname)                  : returns true if dirname points to an existing directory
% conn_fileutils('cd', dirname)                     : changes current working directory
% conn_fileutils('homedir')                         : returns user-specific home directory
% conn_fileutils('imread')                          : see "help imread"
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
    case {'fileread','readfile'}
        if any(conn_server('util_isremotefile',varargin{1})), varargout={conn_server('run',mfilename,option,conn_server('util_localfile',varargin{1}),varargin{2:end})};
        else varargout={fileread(varargin{:})};
        end
        
    case {'filewrite','writefile','filewrite_raw','writefile_raw','fileappend','appendfile','fileappend_raw','appendfile_raw'}
        if any(conn_server('util_isremotefile',varargin{1})), conn_server('run',mfilename,option,conn_server('util_localfile',varargin{1}),varargin{2:end});
        else
            if isempty(regexp(lower(option),'_raw$')), strcmd='%s\n'; else strcmd='%s'; end
            if isempty(regexp(lower(option),'append')), filecmd='wt'; else filecmd='at'; end
            filename=varargin{1};
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
            [ok,nill]=system(['ln -fs ''',varargin{1},''' ''',varargin{2},'''']);
            if ~isequal(ok,0), error('Error linking file %s to %s, missing file or invalid file permissions',varargin{1},varargin{2}); end
        end
        
    case {'filemove','movefile'}
        if any(conn_server('util_isremotefile',varargin{1})), 
            conn_server('run',mfilename,option,conn_server('util_localfile',varargin{1}),conn_server('util_localfile',varargin{2}));
        elseif iscell(varargin{1}), cellfun(@(a,b)conn_fileutils(option,a,b),varargin{1},varargin{2},'uni',0);
        else
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
            fh=fopen(varargin{1},'wb');
            if isequal(fh,-1), error('Error creating file %s, check file permissions',varargin{1}); end
            fclose(fh);
        end
        
    case 'mkdir' % note: no error if fails
        if any(conn_server('util_isremotefile',varargin{1})), 
            [ok,msg]=conn_server('run',mfilename,option,conn_server('util_localfile',varargin{1}),varargin{2:end});
        else
            [ok,msg]=mkdir(varargin{:});
        end
        if nargout>=1, varargout{1}=ok; end
        if nargout>=2, varargout{2}=msg; end
        
    case {'rmdir','rmdir_dironly','rmdir_recursive'} % note: no error if fails
        if any(conn_server('util_isremotefile',varargin{1})), 
            [ok,msg]=conn_server('run',mfilename,option,conn_server('util_localfile',varargin{1}));
        elseif iscell(varargin{1}), [ok,msg]=cellfun(@(x)conn_fileutils(option,x),varargin{1},'uni',0);
        else
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
            out=isdir(varargin{1});
        end
        varargout{1}=out;
        
    case 'dir'
        if numel(varargin)<1||isempty(varargin{1}), filename='.'; else filename=varargin{1}; end
        if numel(varargin)<2||isempty(varargin{2}), option2=''; else option2=varargin{2}; end
        
        if any(conn_server('util_isremotefile',filename)), out=conn_server('run',mfilename,option, conn_server('util_localfile',filename),option2);
        elseif isequal(option2,'-ls'), out=ls('-al',filename);
        else out=dir(filename);
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
        else cd(varargin{1});
        end
                
    case 'homedir'
        if ispc, varargout{1}=conn_fullfile(getenv('USERPROFILE'));
        else varargout{1}=conn_fullfile('~/');
        end
        
    case 'java.io.file'
        if any(conn_server('util_isremotefile',varargin{1})), [varargout{1:nargout}]=conn_server('run',mfilename,option,conn_server('util_localfile',varargin{1}),varargin{2:end});
        else [varargout{1:nargout}]=java.io.File(varargin{:});
        end
        
    case 'getdiskspace'
        if any(conn_server('util_isremotefile',varargin{1})), [varargout{1:nargout}]=conn_server('run',mfilename,option,conn_server('util_localfile',varargin{1}),varargin{2:end});
        else
            a=java.io.File(varargin{:});
            varargout={struct('getUsableSpace',a.getUsableSpace,'getTotalSpace',a.getTotalSpace,'canWrite',a.canWrite)};
        end

    case 'imread',
        if any(conn_server('util_isremotefile',varargin{1})), [varargout{1:nargout}]=imread(conn_cache('pull',varargin{1}),varargin{2:end});
        else [varargout{1:nargout}]=imread(varargin{:});
        end
        
    case 'spm_unlink'
        if any(conn_server('util_isremotefile',varargin)), 
            files=conn_server('util_localfile',varargin); 
            conn_server('run',mfilename,option,files{:});
        else spm_unlink(varargin{:});
        end
        
    case 'spm_jsonwrite'
        if any(conn_server('util_isremotefile',varargin{1})), [varargout{1:nargout}]=conn_server('run',mfilename,option,conn_server('util_localfile',varargin{1}),varargin{2:end});
        else [varargout{1:nargout}]=spm_jsonwrite(varargin{:});
        end
                
    case 'nifti'
        if any(conn_server('util_isremotefile',varargin{1})), varargout{1}=conn_server('run',mfilename,option,conn_server('util_localfile',varargin{1}),varargin{2:end});
        else varargout{1}=struct(nifti(varargin{:}));
        end
        
    case 'spm_vol'
        if any(conn_server('util_isremotefile',varargin{1})), 
            V=conn_server('run',mfilename,option,conn_server('util_localfile',varargin{1}),varargin{2:end});
            if iscell(V), for n=1:numel(V), V{n}.fname=conn_server('util_remotefile',V{n}.fname); end
            else for n=1:numel(V), V(n).fname=conn_server('util_remotefile',V(n).fname); end
            end
        else V=spm_vol(varargin{:});
        end
        varargout{1}=V;
        
    case 'spm_localvol' % pull file to local cache and returns spm_vol of local file
        if isstruct(varargin{1}), % convert struct to local
            if isstruct(varargin{1})&&isfield(varargin{1},'fname')&&any(conn_server('util_isremotefile',varargin{1}(1).fname)), vol=spm_vol(char(conn_cache('pull',conn_expandframe(varargin{1}))));
            else vol=varargin{1};
            end
        elseif any(conn_server('util_isremotefile',varargin{1})), vol=spm_vol(char(conn_cache('pull',varargin{1})));
        else vol=spm_vol(char(varargin{1}));
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
                else data=spm_read_vols(spm_vol(varargin{1}),varargin{2:end});
                end
            elseif isstruct(varargin{1})&&isfield(varargin{1},'fname')&&any(conn_server('util_isremotefile',varargin{1}(1).fname)),
                V=varargin{1};
                if iscell(V), for n=1:numel(V), V{n}.fname=conn_server('util_localfile',V{n}.fname); end
                else for n=1:numel(V), V(n).fname=conn_server('util_localfile',V(n).fname); end
                end
                data=conn_server('run',mfilename,option,V,varargin{2:end});
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
            if isstruct(varargin{1})&&isfield(varargin{1},'fname')&&any(conn_server('util_isremotefile',varargin{1}(1).fname)),
                V=varargin{1};
                if iscell(V), for n=1:numel(V), V{n}.fname=conn_server('util_localfile',V{n}.fname); end
                else for n=1:numel(V), V(n).fname=conn_server('util_localfile',V(n).fname); end
                end
                data=conn_server('run',mfilename,option,V,varargin{2:end});
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
            if isstruct(varargin{1})&&isfield(varargin{1},'fname')&&any(conn_server('util_isremotefile',varargin{1}(1).fname)),
                V=varargin{1};
                if iscell(V), for n=1:numel(V), V{n}.fname=conn_server('util_localfile',V{n}.fname); end
                else for n=1:numel(V), V(n).fname=conn_server('util_localfile',V(n).fname); end
                end
                [varargout{1:nargout}]=conn_server('run',mfilename,option,V,varargin{2:end});
            else [varargout{1:nargout}]=spm_sample_vol(varargin{:});
            end
        end
        
    case 'spm_write_vol'
        if isstruct(varargin{1})&&isfield(varargin{1},'fname')&&any(conn_server('util_isremotefile',varargin{1}(1).fname)),
            V=varargin{1};
            if iscell(V), for n=1:numel(V), V{n}.fname=conn_server('util_localfile',V{n}.fname); end
            else for n=1:numel(V), V(n).fname=conn_server('util_localfile',V(n).fname); end
            end
            V=conn_server('run',mfilename,option,V,varargin{2:end});
            if iscell(V), for n=1:numel(V), V{n}.fname=conn_server('util_remotefile',V{n}.fname); end
            else for n=1:numel(V), V(n).fname=conn_server('util_remotefile',V(n).fname); end
            end
        else V=spm_write_vol(varargin{:});
        end
        if nargout>=1, varargout{1}=V; end
        
    case 'spm_file_merge'
        if isstruct(varargin{1})&&isfield(varargin{1},'fname')&&any(conn_server('util_isremotefile',varargin{1}(1).fname)),
            V=varargin{1};
            if iscell(V), for n=1:numel(V), V{n}.fname=conn_server('util_localfile',V{n}.fname); end
            else for n=1:numel(V), V(n).fname=conn_server('util_localfile',V(n).fname); end
            end
            conn_server('run',mfilename,option,V,conn_server('util_localfile',varargin{2}),varargin{3:end});
        elseif ischar(varargin{1})&&any(conn_server('util_isremotefile',varargin{1})),
            conn_server('run',mfilename,option,conn_server('util_localfile',varargin{1}),conn_server('util_localfile',varargin{2}),varargin{3:end});
        else spm_file_merge(varargin{:});
        end
            
    otherwise
        error('unknown option %s',option);
end
if ~nargout, varargout={}; end
end
