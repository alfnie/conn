function varargout = conn_server(option, varargin)
% CONN_SERVER server over TCP/IP
% 
% commands from server side:
%   conn_server('start' [,PORT])                    : starts CONN server (e.g. run this command from a cluster environment)
%
% commands from client side:
%   conn_server('connect', IP , SERVER_ID)          : connects to CONN server (e.g. run this command from your laptop)
%   conn_server('run', fcn, arg1, arg2, ...)        : runs on server the command fcn(arg1,arg2,...); note: fcn must be a CONN function name, e.g. conn_server('run','conn','load','myfile.mat')
%   [var1, var2,...]=conn_server('run',...)         : as 'run' but also collecting the output(s) of the fcn call
%   conn_server('run_immediatereturn',...)          : as 'run' but without waiting for the command to finish on the server (or check for possible error messages)
%   conn_server('push',file_local,file_remote)      : copies file from client to server
%   conn_server('pull', file_remote,file_local)     : copies file from server to client
%   conn_server('ping')                             : pings the server and shows round-trip times
%   conn_server('disconnect')                       : disconnects from CONN server (a new connection with the same server can be stablished using "conn_server('connect',...)" again from this other machine)
%   conn_server('isconnected')                      : tests bidirectional transfer and returns true if connection is alive
%   conn_server('exit')                             : stops the CONN server and disconnects (no new connections are allowed until the server is restarted using "conn_server('start',...)" again from the server side)
%

%
% development reference notes on use of conn_server and conn_cache in CONN internal functions:
%
%   usage model #1: (a function will run remotely if it needs to work with remote files) (optimal when the fcn call requires minimal transfer of information)
%   function fileout = fcn(filein, varargin)
%      if any(conn_server('util_isremotefile',filin)), fileout=conn_server('util_remotefile',conn_server('run',mfilename,conn_server('util_localfile',filein),varargin{:})); return; end
%
%   usage model #2: (a function will work on a local cache/copy of remote files) (optimal when working with mixture of local and remote files)
%   function fileout = fcn(filein, varargin)
%      isremotefile=conn_server('util_isremotefile',filein);
%      if isremotefile, remotefilename=filename; filename=conn_cache('new',filename); end
%      ...
%      if isremotefile, conn_cache('push',remotefilename); end
%   end
%
%   usage model #3: (a function that accesses files only through calls to functions that follow model #1 or model #2)
%

persistent params
if ~nargin||isempty(option), option='start'; end

varargout={};
if isempty(params)
    params=struct(...
        'isserver',false,...
        'state','off');
end

switch(lower(option))
    case 'start' % init server
        params.isserver=true;
        if numel(varargin)>=1&&~isempty(varargin{1}), port=varargin{1}; else port=[]; end
        if numel(varargin)>=2&&~isempty(varargin{2}), id=varargin{2}; else id=char(mlreportgen.utils.hash(mat2str(now))); id=id(1:8); end
        ok=false;
        while ~ok
            try
                conn_tcpip('open','server',port,id,2);
                params.state='on';
                ok=true;
            catch me
                disp(me.message)
            end
            pause(rand);
        end
        conn_cache clear;
        conn_jobmanager clear;
        conn_server('continue');
        
    case 'connect' % init client
        if numel(varargin)>=1&&~isempty(varargin{1}), ip=varargin{1}; else ip=[]; end
        if numel(varargin)>=2&&~isempty(varargin{2}), sid=varargin{2}; else sid=[]; end
        params.isserver=false;
        port=str2double(regexprep(sid,'CONN.*',''));
        id=regexprep(sid,'^\d*CONN','');
        ok=false;
        while ~ok
            try
                params.state='off';
                conn_tcpip('open','client',ip,port,id);
                ok=true;
                remote_ver=conn_server('run','conn','ver');
                local_ver=conn('ver');
                if ~isequal(local_ver, remote_ver)
                    fprintf('Unable to continue: this machine is running a different CONN version from the server (local = %s, server = %s)\n',local_ver,remote_ver);
                    conn_server('disconnect');
                    return
                end
                params.state='on';
            end
            pause(rand);
        end
        conn_cache clear;
        conn_jobmanager clear;
        
    case {'run','run_immediatereturn'}
        data.cmd=varargin;
        data.nargout=nargout;
        data.immediatereturn=strcmp(lower(option),'run_immediatereturn');
        conn_tcpip('write',data);
        if ~data.immediatereturn;
            ok=false;
            while ~ok
                try
                    var=conn_tcpip('read');
                    ok=true;
                catch me,
                    if ~isempty(regexp(me.message,'SocketTimeoutException')), fprintf('.');  % timeout
                    elseif ~isempty(regexp(me.message,'EOFException|IOException|SocketException')) % restart
                        fprintf('\n Connection to server lost.'); 
                        conn_server restart;
                        conn_tcpip('write',data);
                    else error('ERROR: connection problem: %s',me.message); 
                    end
                end
            end
            assert(iscell(var)&~isempty(var), 'ERROR: unexpected data from server');
            if isequal(var{1},'ok')
                varargout=var(2:end);
            else
                error('ERROR MESSAGE RETURNED BY SERVER: %s\n',var{2});
            end
        end
        
    case {'continue','test'}
        if params.isserver % server
            % listening
            dispstr=sprintf(' CONN SERVER ACTIVE %s',datestr(now));
            fprintf('%s',dispstr);
            while 1
                data=[];
                try, 
                    data=conn_tcpip('read'); 
                catch me,
                    if ~isempty(regexp(me.message,'EOFException|IOException|SocketException')) % restart
                        dispstr='';
                        fprintf('\n Connection to client lost.');
                        conn_server restart;
                    elseif ~isempty(regexp(me.message,'SocketTimeoutException')) % timeout
                        fprintf(repmat('\b',[1,length(dispstr)]));
                        dispstr=sprintf(' CONN SERVER ACTIVE %s',datestr(now));
                        fprintf(dispstr);
                    else
                        disp(me.message);
                        pause(rand);
                    end
                end
                if isempty(data)
                elseif isequal(data,'restart')
                    conn_server restart;
                elseif isequal(data,'exit')
                    fprintf('\n Server closed by client\n'); dispstr='';
                    conn_tcpip('close');
                    return
                elseif isequal(data,'ping')
                    conn_tcpip('write','ok');
                elseif isstruct(data)&&numel(data)==1&&isfield(data,'cmd')
                    if ischar(data.cmd), data.cmd={data.cmd}; end
                    data.cmd{1}=regexprep(data.cmd{1},'\.m$','');
                    if isequal(data.cmd{1},'conn')||~isempty(regexp(data.cmd{1},'^conn_'))||isequal(data.cmd{1},'spm')||~isempty(regexp(data.cmd{1},'^spm_')) % run only conn or spm commands
                        fprintf('\n -'); dispstr='';
                        fh=eval(sprintf('@%s',data.cmd{1}));
%disp([data.cmd]);
                        if ~isfield(data,'immediatereturn')||~data.immediatereturn
                            if isfield(data,'nargout')&&data.nargout>0,
                                argout=cell(1,data.nargout+1);
                                try
                                    [argout{2:data.nargout+1}]=feval(fh,data.cmd{2:end});
                                    argout{1}='ok';
                                catch me
                                    argout={'ko', char(getReport(me,'basic','hyperlinks','off'))};
                                end
                            else
                                try
                                    feval(fh,data.cmd{2:end});
                                    argout={'ok'};
                                catch me
                                    argout={'ko', char(getReport(me,'basic','hyperlinks','off'))};
                                end
                            end
                            disp([data.cmd, argout]);
                            conn_tcpip('write',argout);
                        else
                            try
                                feval(fh,data.cmd{2:end});
                            end
                            disp(data.cmd);
                        end
                    elseif ~isfield(data,'immediatereturn')||~data.immediatereturn
                        argout={'ko', sprintf('unrecognized option %s ',data.cmd{1})};
                        conn_tcpip('write',argout);
                    end
                else % testing
                    fprintf('\n received test data:\n');  dispstr='';
                    disp(data)
                end
            end
        else % client testing
            % sending rand
            for niter=1:1
                if niter>1, pause(1); end
                data={niter, rand(1,1e6)};
                fprintf('\n sending '); disp(data)
                conn_tcpip('write',data);
                fprintf(' sent\n'); dispstr='';
            end
        end
        
    case 'push'
        filename1=varargin{1};
        filename2=varargin{2};
        if ischar(filename1), filename1=cellstr(filename1); end
        if ischar(filename2), filename2=cellstr(filename2); end
        filename2=regexprep(filename2,'[\\\/]CONNSERVER','');
        assert(numel(filename1)==numel(filename2), 'mismatched number of source/destination files');
        ok=false(1,numel(filename1));
        for n=1:numel(filename1)
            if conn_existfile(filename1{n}), 
                conn_server('run_immediatereturn','conn_tcpip','readtofile',filename2{n});
                ok(n)=conn_tcpip('writefromfile',filename1{n});
            else error('file %s not found', filename1{n}); 
            end
        end
        assert(all(ok), 'communication with server interrupted during file push');
        
    case 'pull'
        filename1=varargin{1};
        filename2=varargin{2};
        if ischar(filename1), filename1=cellstr(filename1); end
        if ischar(filename2), filename2=cellstr(filename2); end
        filename1=regexprep(filename1,'[\\\/]CONNSERVER','');
        assert(numel(filename1)==numel(filename2), 'mismatched number of source/destination files');
        hash=cell(1,numel(filename1));
        ok=false(1,numel(filename1));
        for n=1:numel(filename1)
            conn_server('run_immediatereturn','conn_tcpip','writefromfile',filename1{n});
        end
        for n=1:numel(filename1)
            hash{n}=conn_tcpip('readtofile',filename2{n});
            ok(n)=~isempty(hash{n});
        end
        assert(ok(1), 'communication with server interrupted during file pull');
        varargout=hash;
        
    case 'disconnect' % disconnect
        fprintf('Disconnecting from server\n');
        conn_tcpip('close');
        params.state='off';
        conn_cache clear;
        conn_jobmanager clear;
        
    case 'exit' % disconnect & close server
        fprintf('Terminating server and disconnecting\n');
        conn_tcpip('write','exit');
        conn_tcpip('close');
        params.state='off';
        conn_cache clear;
        conn_jobmanager clear;

    case 'restart'
        if params.isserver
            fprintf(' Restarting connection...\n'); 
            connection=conn_tcpip('private');
            ok=false;
            while ~ok
                try
                    conn_tcpip('open','server',connection.port,connection.id,0);
                    ok=true;
                end
                pause(10+rand);
            end
        else
            fprintf(' Restarting connection...\n');
            connection=conn_tcpip('private');
            ok=false;
            while ~ok
                try
                    conn_tcpip('open','client',connection.ip,connection.port,connection.id,0);
                    ok=true;
                end
                pause(10+rand);
            end
        end
        
    case 'isserver'
        if numel(varargin)>=1, params.isserver=varargin{1}; end
        varargout={params.isserver};
        
    case 'state'
        if numel(varargin)>=1, params.state=varargin{1}; end
        varargout={params.state};
        
    case 'ping'
        for niter=1:10
            if niter>1, pause(rand); end
            t1=clock;  
            conn_tcpip('write','ping');
            try, pong=conn_tcpip('read');
            catch, pong='-';
            end
            t2=etime(clock,t1);
            fprintf('received %s from %s: round-trip time = %d ms\n',pong,conn_tcpip('private.ip'),ceil(1000*t2));
        end
        
    case 'isconnected'
        try, 
            remote_ver=conn_server('run','conn','ver');
            local_ver=conn('ver');
            assert(isequal(local_ver, remote_ver));
            params.state='on';
            varargout={true};
        catch
            params.state='off';
            varargout={false};
        end
        
    case 'util_isremotefile'
        varargout={false};
        try, varargout={cellfun(@(x)~isempty(regexp(x,'^\s*[\\\/]CONNSERVER')),cellstr(varargin{1}))};
        end
        
    case 'util_localfile'
        varargout=varargin;
        try
            ischarfilename=ischar(varargin{1});
            filename=regexprep(cellstr(varargin{1}),'^\s*[\\\/]CONNSERVER','');
            if ischarfilename, filename=char(filename); end
            varargout={filename};
        end

    case {'util_remotefile','util_remotefile_keeprelative'}
        varargout=varargin;
        try
            ischarfilename=ischar(varargin{1});
            filename=cellstr(varargin{1});
            if strcmpi(option,'util_remotefile_keeprelative'), change=cellfun('length',regexp(filename,'[\\\/]'))>0;
            else change=true(size(filename));
            end            
            filename(change)=regexprep(filename(change),'^\s*[\\\/]CONNSERVER','');
            filename(change)=regexprep(filename(change),'^\s*[\\\/]*(.*)',[filesep,'CONNSERVER',filesep,'$1']);
            if ischarfilename, filename=char(filename); end
            varargout={filename};
        end

    otherwise,
        error('unrecognized option %s',option);
end
end