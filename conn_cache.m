function varargout=conn_cache(option, varargin)
% CONN_CACHE manages local cache & server drives
%
% commands:
%   filename_local=conn_cache('pull',filename_remote)           : copies file from remote storage to local cache (note: 
%                                                                 filename_remote should include the FULL path to the file)
%   filename_local=conn_cache('new',filename_remote)            : creates local cache for a yet-to-be-created file in remote 
%                                                                 storage (note: file in remote storage will only be created 
%                                                                 after a push)
%   conn_cache('push',filename_remote)                          : copies file from local cache to remote storage
%   conn_cache('pushall')                                       : copies all files in local cache to remote storage
%   
%   conn_cache                                                  : initializes CONN drive (clear all memory and cached files)
%   conn_cache('setlocal', folder_local)                        : defines local cache folder (default: ~/.conn_cache)
%   conn_cache('rename',filename1_remote,filename2_remote)      : reassigns remote storage target of file "filename1_remote" 
%                                                                 to "filename2_remote"
%   conn_cache('sethash',method)                                : defines algorithm used to identify file-changes ('md5' 
%                                                                 or 'timestamp'; see "help conn_tcpip")
%
%   note: CONN_CACHE accepts /CONNSERVER/[filepath] nomenclature for remote files in conn_server machine, e.g.
%     conn_cache('pull','/Volumes/usb-disk/data/myfile.nii')    pulls file from the /data folder within the Volumes/usb-disk drive
%     conn_cache('pull','/CONNSERVER/data/myfile.nii')          pulls file from the /data folder within the machine running conn_server
%

% note: conn_server related functions
%   conn_server('util_isremotefile',filename)       % true if filename is file in server (i.e. filepath is /CONNSERVER/...)
%   conn_server('util_remotefile',filename)         % adds /CONNSERVER/ to filepath
%   conn_server('util_localfile',filename)          % removes /CONNSERVER/ from filepath


persistent params

varargout={};
if isempty(params)
    params=struct(...
        'cachedserver',true,...
        'hash','timestamp',...  % see conn_tcpip hash field
        'local_folder',[],...
        'local_files',{{}},...
        'local_hashes',{{}},...
        'remote_files',{{}},...
        'remote_hashes',{{}});        
    if ispc, params.local_folder=conn_fullfile(getenv('USERPROFILE'),'.conn_cache');
    else params.local_folder=conn_fullfile('~/.conn_cache');
    end
    try, if ~conn_existfile(params.local_folder,2), conn_fileutils('mkdir',params.local_folder); end; end
end
if ~nargin||isempty(option), option='init'; end

switch(lower(option))
    case 'init'
        params.local_files={};
        params.remote_files={};
        params.local_hashes={};
        params.remote_hashes={};
        if ispc, 
            [ok,nill]=system(sprintf('del /Q "%s"',fullfile(params.local_folder,'conncache_*')));
            [ok,nill]=system(sprintf('del /Q "%s"',fullfile(params.local_folder,'cachetmp_*')));
            [ok,nill]=system(sprintf('del /Q "%s"',fullfile(params.local_folder,'conntcpipwrite_*')));
        else
            [ok,nill]=system(sprintf('rm -f ''%s''/conncache_*',params.local_folder));
            [ok,nill]=system(sprintf('rm -f ''%s''/cachetmp_*',params.local_folder));
            [ok,nill]=system(sprintf('rm -f ''%s''/conntcpipwrite_*',params.local_folder));
        end
        fprintf('CONN drive initialized\nLocal/cache folder: %s\n', params.local_folder);
        
    case {'pull','new'}   % conn drive pull <remotefile>
        filename_remote=varargin{1};
        if iscell(filename_remote),  
            fnum=regexp(filename_remote,',\d+\s*$','match','once');
            filename_remote=regexprep(filename_remote,',\d+\s*$','');
            [ufilename_remote,nill,idx]=unique(filename_remote);
            filename_local=cellfun(@(x)conn_cache(option,x),ufilename_remote,'uni',0); 
            filename_local=cellfun(@(a,b)[a,b],reshape(filename_local(idx),size(filename_remote)),fnum,'uni',0);
        else
            fnum=regexp(filename_remote,',\d+\s*$','match','once');
            filename_remote=regexprep(filename_remote,',\d+\s*$','');
            remote_inserver=conn_server('util_isremotefile',filename_remote);
            filename_remote=conn_fullfile(filename_remote);
            [fext,fexts]=conn_cache_exts(filename_remote);
            idx=find(strcmp(filename_remote,params.remote_files),1,'last');
            if isempty(idx), filename_local=conn_fullfile(params.local_folder, ['conncache_', char(mlreportgen.utils.hash([filename_remote, mat2str(now)])), fext]);
            else filename_local=params.local_files{idx};
            end
            if isempty(idx), idx=numel(params.remote_files)+1; end
            if strcmpi(option,'new'),
                local_hash=[];
                remote_hash=[];
            else
                f1=cellfun(@(x)conn_prepend('',filename_remote,x),fexts,'uni',0);
                f2=cellfun(@(x)conn_prepend('',filename_local,x),fexts,'uni',0);
                in=conn_existfile(f1);
                assert(in(1),'unable to find file %s',filename_remote);
                if remote_inserver
                    local_hash=[];
                    remote_hash=conn_server('run','conn_cache','hash',conn_server('util_localfile',f1(in)),params.hash);
                    changed=in;
                    if params.cachedserver&&idx<=numel(params.remote_files)&&~isempty(params.local_hashes{idx})&&~isempty(params.remote_hashes{idx}) % tries to see if already matched
                        local_hash=conn_cache_hash(f2(in),params.hash);
                        if isequal(params.local_hashes{idx},local_hash)&&isequal(params.remote_hashes{idx},remote_hash)&&(isequal(params.hash,'timestamp')||isequal(local_hash,remote_hash))
                            varargout={[filename_local fnum]};
                            return
                        end
                        if all([size(local_hash,2), size(remote_hash,2), size(params.local_hashes{idx},2), size(params.remote_hashes{idx},2)]==nnz(in))
                            ok1=cellfun(@(a,b)isequal(a,b),num2cell(params.local_hashes{idx},1),num2cell(local_hash,1));
                            ok2=cellfun(@(a,b)isequal(a,b),num2cell(params.remote_hashes{idx},1),num2cell(remote_hash,1));
                            if isequal(params.hash,'timestamp'), ok3=true; else ok3=cellfun(@(a,b)isequal(a,b),num2cell(local_hash,1),num2cell(remote_hash,1)); end
                            changed(in)=~(ok1&ok2&ok3);
                        end
                    end
                    [hash{1:nnz(in&changed)}]=conn_server('pull',f1(in&changed),f2(in&changed));
                    assert(all(cellfun('length',hash)),'Failed connection');
                    local_hash(:,changed(in))=cat(2,hash{cellfun('length',hash)>0});
                else
                    for n=reshape(find(in),1,[])
                        conn_cache_copyfile(f1{n}, f2{n});
                    end
                    local_hash=[];
                    remote_hash=[];
                end
            end
            params.local_files{idx}=filename_local;
            params.local_hashes{idx}=local_hash;
            params.remote_files{idx}=filename_remote;
            params.remote_hashes{idx}=remote_hash;
            filename_local=[filename_local fnum];
        end
        if nargout, varargout={filename_local}; end
        
    case 'pushall'
        for n=1:numel(params.remote_files)
            conn_cache('push',params.remote_files{n});
        end
        if nargout, varargout={params.local_files}; end
        
    case {'push'}   % conn drive push <remotefile> [<localfile>]
        filename_remote=varargin{1};
        if iscell(filename_remote),  
            fnum=regexp(filename_remote,',\d+\s*$','match','once');
            filename_remote=regexprep(filename_remote,',\d+\s*$','');
            [ufilename_remote,nill,idx]=unique(filename_remote);
            filename_local=cellfun(@(x)conn_cache(option,x),ufilename_remote,'uni',0); 
            filename_local=cellfun(@(a,b)[a,b],reshape(filename_local(idx),size(filename_remote)),fnum,'uni',0);
        else
            fnum=regexp(filename_remote,',\d+\s*$','match','once');
            filename_remote=regexprep(filename_remote,',\d+\s*$','');
            remote_inserver=conn_server('util_isremotefile',filename_remote);
            filename_remote=conn_fullfile(filename_remote);
            if numel(varargin)>=2&&~isempty(varargin{2}),
                filename_local=varargin{2};
                idx=find(strcmp(filename_remote, params.remote_files),1,'last');
                if isempty(idx), % init new entry
                    idx=numel(params.remote_files)+1;
                    params.local_files{idx}=filename_local;
                    params.local_hashes{idx}=[];
                    params.remote_files{idx}=filename_remote;
                    params.remote_hashes{idx}=[];
                end
            else
                idx=find(strcmp(filename_remote, params.remote_files),1,'last');
                if isempty(idx), error('Could not find a match for file %s. Please enter explicit filename_local argument',filename_remote);
                else filename_local=params.local_files{idx};
                end
            end
            [fext,fexts]=conn_cache_exts(filename_remote);
            f1=cellfun(@(x)conn_prepend('',filename_local,x),fexts,'uni',0);
            f2=cellfun(@(x)conn_prepend('',filename_remote,x),fexts,'uni',0);
            in=conn_existfile(f1);
            assert(in(1),'unable to find file %s',filename_local);
            if remote_inserver
                local_hash=conn_cache_hash(f1(in),params.hash);
                remote_hash=[];
                changed=in;
                if params.cachedserver&&idx<=numel(params.remote_files)&&~isempty(params.local_hashes{idx})&&~isempty(params.remote_hashes{idx}) % tries to see if already matched
                    remote_hash=conn_server('run','conn_cache','hash',conn_server('util_localfile',f2(in)),params.hash);
                    if isequal(params.local_hashes{idx},local_hash)&&isequal(params.remote_hashes{idx},remote_hash)&&(isequal(params.hash,'timestamp')||isequal(local_hash,remote_hash))
                        varargout={[filename_local fnum]};
                        return
                    end
                    if all([size(local_hash,2), size(remote_hash,2), size(params.local_hashes{idx},2), size(params.remote_hashes{idx},2)]==nnz(in))
                        ok1=cellfun(@(a,b)isequal(a,b),num2cell(params.local_hashes{idx},1),num2cell(local_hash,1));
                        ok2=cellfun(@(a,b)isequal(a,b),num2cell(params.remote_hashes{idx},1),num2cell(remote_hash,1));
                        if isequal(params.hash,'timestamp'), ok3=true; else ok3=cellfun(@(a,b)isequal(a,b),num2cell(local_hash,1),num2cell(remote_hash,1)); end
                        changed(in)=~(ok1&ok2&ok3);
                    end
                end
                conn_server('push',f1(in&changed),f2(in&changed));
                remote_hash(:,changed(in))=conn_server('run','conn_cache','hash',conn_server('util_localfile',f2(in&changed)),params.hash);
            else
                for n=reshape(find(in),1,[])
                    conn_cache_copyfile(f1{n}, f2{n});
                end
                local_hash=[];
                remote_hash=[];
            end
            params.local_hashes{idx}=local_hash;
            params.remote_hashes{idx}=remote_hash;
            filename_local=[filename_local fnum];
        end
        if nargout, varargout={filename_local}; end
        
    case 'clear'   % conn cache clear filename_remote
        if isempty(varargin), conn_cache('init');
        else
            filename_remote=varargin{1};
            remote_inserver=conn_server('util_isremotefile',filename_remote);
            filename_remote=conn_fullfile(filename_remote);
            idx=find(strcmp(filename_remote, params.remote_files),1,'last');
            if isempty(idx), error('Could not find a match for file %s. Please enter explicit filename_local argument',filename_remote);
            else
                %conn_fileutils('deletefile',params.local_files{idx});
                params.local_files(idx)=[];
                params.remote_files(idx)=[];
                params.local_hashes(idx)=[];
                params.remote_hashes(idx)=[];
            end
        end
        
    case 'rename'   % conn drive rename <remotefile> <remotefile_newname>
        filename_remote=varargin{1};
        remote_inserver=conn_server('util_isremotefile',filename_remote);
        filename_remote=conn_fullfile(filename_remote);
        idx=find(strcmp(filename_remote, params.remote_files),1,'last');
        if isempty(idx), error('Could not find a match for file %s. Please enter explicit filename_local argument',filename_remote);
        else
            params.remote_files{idx}=conn_fullfile(varargin{2});
            params.remote_hashes{idx}=[];
        end
        
    case 'hash' % conn cache has filename
        varargout={conn_cache_hash(varargin{:},params.hash)};
        
    case 'clearhash'  % conn drive clearhash
        params.local_hashes=cell(size(params.local_hashes));
        params.remote_hashes=cell(size(params.remote_hashes));
        
    case 'setlocal'  % conn drive setlocal <path>
        params.local_folder=varargin{1};
        
    case 'cachedserver',
        params.cachedserver=varargin{1};
        if ischar(params.cachedserver), params.cachedserver=str2double(params.cachedserver); end
        
    case 'sethash',
        params.hash=varargin{1};
        conn_tcpip('sethash',varargin{1});
        
    case 'private.local_folder'
        varargout={params.local_folder};
    case 'private',
        varargout={params};
        
    otherwise
        error('unrecognized option %s',option);
        
end
        
end

function conn_cache_copyfile(a,b,varargin)
if ispc, [ok,nill]=system(['copy "',a,'" "',b,'"']);
else, [ok,nill]=system(['''cp'' -f ''',a,''' ''',b,'''']);
end
if ~isequal(ok,0), error('Error copying file %s to %s, check target permissions',a,b); end
end

function hash=conn_cache_hash(filename,htype,varargin)
if iscell(filename)
    hash=cellfun(@(x)conn_cache_hash(x,htype),filename,'uni',0);
    hash=cat(2,hash{:});
elseif strcmp(htype,'timestamp')
    fh=dir(filename);
    hash=fh.datenum;
else
    maxsize=65536;
    try
        fh=fopen(filename,'rb');
        hh=java.security.MessageDigest.getInstance(htype);
        while 1
            tdata=fread(fh, maxsize, 'uint8')';
            if isempty(tdata), break; end
            hh.update(tdata);
        end
        fclose(fh);
        hash=typecast(hh.digest,'uint8');
    catch
        hash=[];
    end
end
end

function [fext,fexts]=conn_cache_exts(filename)
[nill,nill,fext]=spm_fileparts(filename);
fexts={fext};
if ismember(fext,{'.nii','.img'}), fexts=[fexts, {'.mat','.json','.txt','.csv','.xls','.info','.icon.jpg'}]; end
if strcmp(fext,'.img'), fexts=[fexts, {'.hdr'}];
elseif strcmp(fext,'.mat'), fexts=[fexts, {'.txt','.jpg'}];
elseif strcmp(fext,'.txt'), fexts=[fexts, {'.mat','.jpg'}];
elseif strcmp(fext,'.jpg'), fexts=[fexts, {'.mat','.txt'}];
elseif strcmp(fext,'.matc'), fexts=[fexts, {'.mat'}];
end
end
