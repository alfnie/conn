function varargout = conn_server(option, varargin)
% CONN_SERVER server over TCP/IP
%
% commands from server side: (e.g. from a computer at work/university where the data lives)
%   conn_server('start' [,PORT])                    : starts CONN server and prints out instructions for connecting to this server
%
% commands from client side: (e.g. from a computer at home/office where you live)
%   conn_server('connect', IP , SERVER_ID)          : stablishes connection to CONN server
%   conn_server('run', fcn, arg1, arg2, ...)        : runs on server the command fcn(arg1,arg2,...); note: fcn must be a CONN function name, e.g. conn_server('run','conn','load','myfile.mat')
%   [var1, var2,...]=conn_server('run',...)         : as 'run' but also collecting the output(s) of the fcn call
%   conn_server('run_immediatereturn',...)          : as 'run' but without waiting for the command to finish on the server (or check for possible error messages)
%   conn_server('push',file_local,file_remote)      : copies file from client to server
%   conn_server('pull', file_remote,file_local)     : copies file from server to client
%   conn_server('ping')                             : pings the server and shows round-trip times
%   conn_server('isconnected')                      : tests bidirectional transfer and returns true if connection is alive
%   conn_server('disconnect')                       : disconnects from CONN server (a new connection with the same server can be stablished using "conn_server('connect',...)" again from this or other machine)
%   conn_server('exit')                             : stops the remote CONN server and disconnects (no new connections are allowed until the server is restarted, e.g. using "conn_server('start',...)" again from the server side)
%
% alternative procedure for SSH-accessible cluster/HPC enviornment only:
% command from server side: (from a computer within the cluster/HPC environment)
%   conn_server('HPC_save' [, fileout])             : HPC initialization, creates a .json file with basic information necessary to connect to this server (machine IP, and where to find Matlab/SPM/CONN) [~/connserverinfo.json]
% commands from client side: (from a computer outside of the cluster/HPC environment)
%   conn_server('HPC_start' [, serverip])           : SSH to HPC login node, submits a job to start a CONN server, and connects this computer to that remote server (note: expects to find serverip:~/connserverinfo.json file; run "conn_server HPC_save" in HPC once)
%   conn_server('HPC_start', filein)                : same as above but loading HPC server information explicitly from the provided .json file
%   conn_server('HPC_submitstart')                  : submits a job to start a new CONN server (run internally by HPC_start)
%   conn_server('HPC_restart')                      : restarts a dropped connection
%   conn_server('HPC_exit')                         : terminates CONN server and disconnects
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
%

persistent params
if ~nargin||isempty(option), option='start'; end

varargout={};
if isempty(params)
    params=struct(...
        'isserver',false,...
        'info',struct(),...
        'state','off');
end

switch(lower(option))
    case 'start' % init server
        params.isserver=true;
        if numel(varargin)>=1&&~isempty(varargin{1}), port=varargin{1}; else port=0; end
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
        data.hpc=false;
        if isfield(params.info,'host')&&~isempty(params.info.host)&&isfield(params.info,'scp')&&params.info.scp>0, data.hpc=true; end
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
                        fprintf('\n Idle connection to server. ');
                        conn_server restart;
                        conn_tcpip('write',data);
                    else error('ERROR: connection problem: %s',me.message);
                    end
                end
            end
            assert(iscell(var)&~isempty(var), 'ERROR: unexpected data from server');
            if isequal(var{1},'ok_hasattachment') % server is waiting for us to scp its response
                tmpfile1=fullfile(conn_cache('private.local_folder'),['cachetmp_', char(mlreportgen.utils.hash(mat2str(now))),'.mat']);
                tmpfile2=var{2};
                if conn_server('HPC_pull',tmpfile2,tmpfile1), var=load(tmpfile1,'arg'); var=var.arg;
                else var={'ko','unable to run HPC_pull to read response from server'};
                end
                conn_fileutils('deletefile',tmpfile1);
                conn_tcpip('write','ok'); % signal server to continue
            end
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
            ntimedout=0; % minutes
            while 1
                data=[];
                try,
                    data=conn_tcpip('read');
                    ntimedout=0;
                catch me,
                    if ~isempty(regexp(me.message,'EOFException|IOException|SocketException')) %||ntimedout>15 % restart
                        dispstr='';
                        fprintf('\n Idle connection to client.');
                        conn_server restart;
                    elseif ~isempty(regexp(me.message,'SocketTimeoutException')) % timeout
                        fprintf(repmat('\b',[1,length(dispstr)]));
                        dispstr=sprintf(' CONN SERVER ACTIVE %s',datestr(now));
                        fprintf(dispstr);
                        conn_tcpip('writekeepalive');
                        ntimedout=ntimedout+1;
                    else
                        disp(me.message);
                        pause(rand);
                    end
                end
                if isempty(data)
                elseif isequal(data,'handshake')||isequal(data,'restart')
                    conn_server restart;
                elseif isequal(data,'exit')
                    fprintf('\n Server closed by client\n'); dispstr='';
                    conn_tcpip('close');
                    return
                elseif isequal(data,'ping')
                    conn_tcpip('write','ok');
                elseif isstruct(data)&&numel(data)==1&&isfield(data,'cmd') % run command
                    if ischar(data.cmd), data.cmd={data.cmd}; end
                    data.cmd{1}=regexprep(data.cmd{1},'\.m$','');
                    if isequal(data.cmd{1},'conn')||~isempty(regexp(data.cmd{1},'^conn_'))||isequal(data.cmd{1},'spm')||~isempty(regexp(data.cmd{1},'^spm_')) % run only conn or spm commands
                        fprintf('\n -'); dispstr='';
                        try
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
                                if numel(argout)>1&&isequal(argout{1},'ok')&&isfield(data,'hpc')&&data.hpc>0&&getfield(whos('argout'),'bytes')>1e6, % send response using ssh/scp?
                                    tmpfile=fullfile(conn_cache('private.local_folder'),['cachetmp_', char(mlreportgen.utils.hash(mat2str(now))),'.mat']);
                                    arg=argout; save(tmpfile,'arg');
                                    argout={'ok_hasattachment',tmpfile};
                                    conn_tcpip('write',argout);
                                    conn_tcpip('read'); % wait for client to finish
                                    conn_fileutils('deletefile',tmpfile);
                                else % send response through communication socket
                                    conn_tcpip('write',argout);
                                end
                            else
                                try
                                    feval(fh,data.cmd{2:end});
                                end
                                disp(data.cmd);
                            end
                        catch me
                            argout={'ko', char(getReport(me,'basic','hyperlinks','off'))};
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
                if isfield(params.info,'host')&&~isempty(params.info.host)&&isfield(params.info,'scp')&&params.info.scp>0
                    ok(n)=conn_server('HPC_push',filename1{n},filename2{n});
                else
                    conn_server('run_immediatereturn','conn_tcpip','readtofile',filename2{n});
                    ok(n)=conn_tcpip('writefromfile',filename1{n});
                end
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
        if isfield(params.info,'host')&&~isempty(params.info.host)&&isfield(params.info,'scp')&&params.info.scp>0
            for n=1:numel(filename1)
                ok(n)=conn_server('HPC_pull',filename1{n},filename2{n});
                if ok(n), hash{n}=conn_cache('hash',filename2{n}); end
            end
        else        
            for n=1:numel(filename1)
                conn_server('run_immediatereturn','conn_tcpip','writefromfile',filename1{n});
            end
            for n=1:numel(filename1)
                hash{n}=conn_tcpip('readtofile',filename2{n});
                ok(n)=~isempty(hash{n});
            end
        end
        assert(ok(1), 'communication with server interrupted during file pull');
        varargout=hash;
        
    case 'clear'
        conn_cache clear;
        conn_tcpip clear;
        
    case 'hpc_save' % saves .json info
        if numel(varargin)>=1, filename=varargin{1};
        else filename=fullfile(conn_fileutils('homedir'),'connserverinfo.json');
        end
        if ispc, [nill,str1]=system('hostname');
        else [nill,str1]=system('hostname -f');
        end
        osquotes=char('"'*ispc+''''*~ispc);
        [nill,str2]=system('whoami');
        matlabpath=fullfile(matlabroot,'bin');
        whichfolders=cellfun(@(x)fileparts(which(x)),{'spm','conn','evlab17'},'uni',0);
        isdep=false;
        try, isdep=isdeployed; end
        if conn_jobmanager('options','cmd_rundeployed'), isdep=true; end
        if isdep&&~isempty(conn_jobmanager('options','cmd_deployedfile')), fun_callback=conn_jobmanager('options','cmd_deployedfile');
        elseif isdep,                             fun_callback=conn_jobmanager('checkdeployedname');
        else                                      fun_callback=[osquotes fullfile(matlabpath,'matlab') osquotes];
        end
        if isdep,   cmd=sprintf('%s %s',fun_callback,'%s');
        elseif ~isempty(whichfolders{3}), cmd=sprintf('%s -nodesktop -noFigureWindows -nosplash -r "addpath ''%s''; addpath ''%s''; addpath ''%s''; conn %s; exit"',...
                fun_callback, whichfolders{3}, whichfolders{1}, whichfolders{2},'%s');
        else cmd=sprintf('%s -nodesktop -noFigureWindows -nosplash -r "addpath ''%s''; addpath ''%s''; conn %s; exit"',...
                fun_callback, whichfolders{1}, whichfolders{2},'%s');
        end
        info = struct(...
            'host', regexprep(str1,'\n',''),...
            'CONNcmd',cmd);
        if ~isempty(filename), spm_jsonwrite(filename,info); fprintf('Host information saved in %s\n',filename); end
        if nargout>0, varargout={info}; end
        if ~nargout&&isempty(filename), disp(info); end
        
    case {'hpc_start','hpc_restart'} % init server remotely and connect to it
        % development reference notes on ssh tunneling
        %    $ ssh -fN -o ServerAliveInterval=60 -o ServerAliveCountMax=10 -o ControlMaster=yes -o ControlPath=<local_filename> <login_node>        % authenticate first
        %    $ ssh -o ControlPath=<local_filename> -L<local_port>:<server_ip>:<server_port> <login_node>                                           % port forwarding on shared connection
        %    $ ssh -o ControlPath=<local_filename> -O forward -L<local_port>:<server_ip>:<server_port> <login_node>
        %    $ ssh -o ControlPath=<local_filename> -O check <login_node>                                                                             % exit/close connection
        %    $ ssh -o ControlPath=<local_filename> -O exit <login_node>                                                                             % exit/close connection
        if numel(varargin)>=1,
            filename=varargin{1};
            if isempty(regexp(filename,'.json$')), params.info.host=filename;
            else params.info=conn_jsonload(filename);
            end
        else params.info.CONNcmd='';
        end
        if ~isfield(params.info,'host')||isempty(params.info.host), params.info.host=input('Server address [none]: ','s');
        else
            temp=input(sprintf('Server address [%s]: ',params.info.host),'s');
            if ~isempty(temp), params.info.host=regexprep(temp,'\s+',''); end
        end
        if strcmp(params.info.host,'none'), params.info.host=''; end
        if numel(varargin)>=2&&~isempty(varargin{2}), params.info.local_port=varargin{2}; if ischar(params.info.local_port), params.info.local_port=str2double(params.info.local_port); end
        else params.info.local_port=[];
        end
        if ~isfield(params.info,'user')||isempty(params.info.user), [nill,str2]=system('whoami'); params.info.user=regexprep(str2,'\n',''); end
        if isempty(params.info.host),
            params.info.user='';
            params.info.login_ip='';
            params.info.CONNcmd='';
            params.info.filename_ctrl='';
        else
            temp=input(sprintf('Username [%s]: ',params.info.user),'s');
            if ~isempty(temp), params.info.user=temp; end
            params.info.login_ip=sprintf('%s@%s',params.info.user,params.info.host);
            localcachefolder=conn_cache('private.local_folder');
        end
        params.info.scp=false;
        if ispc&&~isempty(params.info.host),
            error('in development');
        else
            if ~isempty(params.info.host)
                params.info.filename_ctrl=fullfile(localcachefolder,sprintf('connserver_ctrl_%s_%s',params.info.host,params.info.user));
                try, conn_fileutils('deletefile',params.info.filename_ctrl); end
                fprintf('Connecting to %s... ',params.info.login_ip);
                % starts a shared SSH connection
                system(sprintf('ssh -f -N -o ControlMaster=yes -o ControlPath=''%s'' %s', params.info.filename_ctrl,params.info.login_ip)); % starts a shared connection (note: use "sleep 600" instead of -N?)
                [ok,msg]=system(sprintf('ssh -o ControlPath=''%s'' -O check %s', params.info.filename_ctrl,params.info.login_ip)); if ok~=0, error(msg); end
                if ~isfield(params.info,'CONNcmd')||isempty(params.info.CONNcmd) % attempts to load server info from remote ~/connserverinfo.json file
                    filename=fullfile(conn_cache('private.local_folder'),['conncache_', char(mlreportgen.utils.hash(mat2str(now)))]);
                    conn_fileutils('deletefile',filename);
                    fprintf('\nDownloading configuration information from %s:%s to %s\n',params.info.login_ip,'~/connserverinfo.json',filename);
                    [ok,msg]=system(sprintf('scp -q -o ControlPath=''%s'' %s:~/connserverinfo.json %s',...
                        params.info.filename_ctrl,params.info.login_ip,filename));
                    assert(conn_existfile(filename),'unable to find ~/connserverinfo.json file in %s',params.info.login_ip);
                    params.info.CONNcmd=conn_jsonread(filename,'CONNcmd');
                end
            end
            if strcmpi(option,'hpc_restart')
                if ~isfield(params.info,'remote_ip')||isempty(params.info.remote_ip), params.info.remote_ip=input('Remote session host address: ','s'); end
                if ~isfield(params.info,'remote_port')||isempty(params.info.remote_port), params.info.remote_port=str2double(input('Remote session access port: ','s')); end
                if ~isfield(params.info,'remote_id')||isempty(params.info.remote_id), params.info.remote_id=input('Remote session id: ','s'); end
                if ~isfield(params.info,'remote_log')||isempty(params.info.remote_log), params.info.remote_log=input('Remote session log folder: ','s'); end
                params.info.local_port=[];
            elseif isempty(params.info.host)
                if ~isfield(params.info,'remote_ip'), params.info.remote_ip=''; end; if isempty(params.info.remote_ip), temp=input('Remote session host address: ','s'); else temp=input(sprintf('Remote session host address [%s]: ',params.info.remote_ip),'s'); end; if ~isempty(temp), params.info.remote_ip=temp; end
                if ~isfield(params.info,'remote_port'), params.info.remote_port=[]; end; if isempty(params.info.remote_port), temp=input('Remote session access port: ','s'); else temp=input(sprintf('Remote session access port [%d]: ',params.info.remote_port),'s'); end; if ~isempty(temp), params.info.remote_port=str2double(temp); end
                if ~isfield(params.info,'remote_id'), params.info.remote_id=''; end; if isempty(params.info.remote_id), temp=input('Remote session id: ','s'); else temp=input(sprintf('Remote session id [%s]: ',params.info.remote_id),'s'); end; if ~isempty(temp), params.info.remote_id=temp; end
                if ~isfield(params.info,'remote_log'), params.info.remote_log=''; end; if isempty(params.info.remote_log), temp=input('Remote session log folder: ','s'); else temp=input(sprintf('Remote session log folder [%s]: ',params.info.remote_log),'s'); end; if ~isempty(temp), params.info.remote_log=temp; end
                params.info.local_port=[];
            else
                fprintf('Requesting a new Matlab session in %s. This may take a few minutes, please be patient as your job currently sits in a queue. CONN will resume automatically when the new Matlab session becomes available\n',params.info.login_ip);
                % submit jobs to start server#2 in arbitrary remote node using HPC scheduler
                [ok,msg]=system(sprintf('ssh -o ControlPath=''%s'' %s "%s"', params.info.filename_ctrl,params.info.login_ip, regexprep(sprintf(params.info.CONNcmd,sprintf('server HPC_submitstart')),'"','\\"')));
                if ~isempty(regexp(msg,'HPC_SUBMITSTART error')), error('Error initiating server job\n %s',msg);
                else
                    keys=regexp(msg,'HOST:([^\n]+)\nPORT:([^\n]+)\nID:([^\n]+)\nLOG:([^\n]+)\n','tokens','once');
                    params.info.remote_ip=keys{1};
                    params.info.remote_port=str2double(keys{2});
                    params.info.remote_id=keys{3};
                    params.info.remote_log=keys{4};
                    fprintf('Remote session started:\n  Host address = %s\n  Access port = %d\n  ID = %s\n  Log folder = %s\n',params.info.remote_ip,params.info.remote_port,params.info.remote_id,params.info.remote_log);
                end
            end
            if isempty(params.info.host)
                if strcmpi(option,'hpc_restart'), conn_tcpip('open','client',params.info.remote_ip,params.info.remote_port,params.info.remote_id,0); params.state='on';
                else conn_server('connect',params.info.remote_ip,sprintf('%dCONN%s',params.info.remote_port,params.info.remote_id));
                end
            else
                % stablishes port-forward link between this computer and server#2
                if isempty(params.info.local_port) % finds first available local port
                    tsocket=java.net.ServerSocket(0);
                    params.info.local_port=tsocket.getLocalPort;
                    tsocket.close;
                    clear tsocket
                    pause(1);
                end
                
                fprintf('Establishing secure communication path to remote session (%d:%s:%d)\n',params.info.local_port,params.info.remote_ip,params.info.remote_port);
                system(sprintf('ssh -o ControlPath=''%s'' -O forward -L%d:%s:%d %s', params.info.filename_ctrl,params.info.local_port,params.info.remote_ip,params.info.remote_port,params.info.login_ip)); 
                %system(sprintf('ssh -f -N -T -o ExitOnForwardFailure=yes -o ControlPath=''%s'' -L%d:%s:%d %s', params.info.filename_ctrl,params.info.local_port,params.info.remote_ip,params.info.remote_port,params.info.login_ip)); 
                 %fprintf('Connecting to server\n');
                if strcmpi(option,'hpc_restart'), conn_tcpip('open','client','localhost',params.info.local_port,params.info.remote_id,0); params.state='on';
                else conn_server('connect','localhost',sprintf('%dCONN%s',params.info.local_port,params.info.remote_id));
                end
            end
        end
        
    case 'hpc_exit' % send exit signal to server to stop running (this will also cause the remote Matlab session to exit)
        if conn_server('isconnected'), 
            fprintf('Exiting remote CONN session and disconnecting\n');
            conn_tcpip('write','exit');
            try, [ok,msg]=system(sprintf('ssh -o ControlPath=''%s'' -O exit %s', params.info.filename_ctrl,params.info.login_ip)); end
            conn_tcpip('close');
            params.state='off';
            conn_cache clear;
            conn_jobmanager clear;
        else
            fprintf('unable to connect to server, please terminate the server manually or use "conn_server HPC_restart" to restart the connection with the server and try "conn_server HPC_exit" again\n');
            system(sprintf('ssh -o ControlPath=''%s'' -O exit %s', params.info.filename_ctrl,params.info.login_ip));
        end
        
    case 'hpc_exitforce' % run remotely a command to forcibly delete the server's job
        if isfield(params.info,'host')&&~isempty(params.info.host)&&isfield(params.info,'login_ip')&&~isempty(params.info.login_ip)&&isfield(params.info,'remote_log')&&~isempty(params.info.remote_log)
            localcachefolder=conn_cache('private.local_folder');
            params.info.filename_ctrl=fullfile(localcachefolder,sprintf('connserver_ctrl_%s_%s',params.info.host,params.info.user));
            fprintf('Connecting to %s... ',params.info.login_ip);
            % starts a shared SSH connection
            [ok,msg]=system(sprintf('ssh -o ControlPath=''%s'' -O check %s', params.info.filename_ctrl,params.info.login_ip)); 
            if ok~=0, 
                try, conn_fileutils('deletefile',params.info.filename_ctrl); end
                system(sprintf('ssh -f -N -o ControlMaster=yes -o ControlPath=''%s'' %s', params.info.filename_ctrl,params.info.login_ip));
            end
            if ~isfield(params.info,'CONNcmd')||isempty(params.info.CONNcmd) % attempts to load server info from remote ~/connserverinfo.json file
                filename=fullfile(conn_cache('private.local_folder'),['conncache_', char(mlreportgen.utils.hash(mat2str(now)))]);
                conn_fileutils('deletefile',filename);
                fprintf('\nDownloading configuration information from %s:%s to %s\n',params.info.login_ip,'~/connserverinfo.json',filename);
                [ok,msg]=system(sprintf('scp -q -o ControlPath=''%s'' %s:~/connserverinfo.json %s',...
                    params.info.filename_ctrl,params.info.login_ip,filename));
                assert(conn_existfile(filename),'unable to find ~/connserverinfo.json file in %s',params.info.login_ip);
                params.info.CONNcmd=conn_jsonread(filename,'CONNcmd');
            end            
            [ok,msg]=system(sprintf('ssh -o ControlPath=''%s'' %s "%s"', params.info.filename_ctrl,params.info.login_ip, regexprep(sprintf(params.info.CONNcmd,sprintf('server HPC_submitexit ''%s''',params.info.remote_log)),'"','\\"')));
            if ok~=0, disp(msg); 
            else fprintf('Remote CONN session deletion requested successfully\n');
            end
        end
        
    case 'hpc_submitstart'
        info=conn_jobmanager('submit','orphan_conn',[],1,[],'server','start',varargin{:});
        fprintf('HPC_log in %s\n',info.pathname);
        fprintf('waiting for job to start...\n')
        for n=1:600
            pause(5.5+rand); % check log files every 6s
            info=conn_jobmanager('statusjob',info,[],~rem(n,10),true); % check queue every minute
            if all(ismember(info.tagmsg,{'running','error','finished'})), break; end
        end
        if all(ismember(info.tagmsg,{'running'}))
            ok=false;
            while ~ok
                str=conn_fileutils('fileread',info.stdout{1});
                ok=~isempty(regexp(str,'Waiting for client connection'));
                if ~ok, pause(1+rand); end
            end
            if ok,
                fprintf('\nHPC_SUBMITSTART finished succesfully\n');
                match1=regexp(str,'\<ssh -L 6111:([^:]*):(\d+)','tokens');
                match2=regexp(str,'\<conn_server connect localhost \d+CONN(\w+)','tokens');
                fprintf('HOST:%s\nPORT:%s\nID:%s\nLOG:%s\n',match1{1}{1},match1{1}{2},match2{1}{1},info.pathname);
            else fprintf('HPC_SUBMITSTART error\n');
            end
        else fprintf('HPC_SUBMITSTART error\n');
        end
        
    case 'hpc_submitexit' % internal use only
        try
            pathname=varargin{1};
            info=struct; conn_loadmatfile(fullfile(pathname,'info.mat'),'info');
            info=conn_jobmanager('canceljob',info);
            fprintf('HPC_EXIT finished\n');
        catch
            fprintf('HPC_EXIT error\n');
        end
        
    case 'hpc_info'
        if numel(varargin)>=1, params.info=varargin{1}; end
        varargout={params.info};
        
    case 'hpc_details'
        clear h;
        info=struct; conn_loadmatfile(fullfile(conn_server('util_remotefile',params.info.remote_log),'info.mat'),'info');
        tfiles={info.stdout, info.stderr, info.stdlog, info.scripts};
        [nill,names]=cellfun(@fileparts,tfiles{1},'uni',0);
        names=regexp(names,'^.{9}','match','once');
        h.hfig=figure('units','norm','position',[.7 .3 .3 .6],'name','log details','numbertitle','off','menubar','none','color','w');
        h.files=uicontrol(h.hfig,'style','popupmenu','units','norm','position',[.1 .95 .9 .05],'string',names,'value',1);
        h.types=uicontrol(h.hfig,'style','popupmenu','units','norm','position',[.1 .90 .9 .05],'string',{'console output (stdout)','error output (stderr)','Matlab log','submission script','submission command','submission command output','status command output'},'value',1);
        h.str=uicontrol(h.hfig,'style','listbox','units','norm','position',[.05 .1 .9 .75],'string','','max',2,'horizontalalignment','left','fontname','monospaced');
        h.refresh=uicontrol(h.hfig,'style','pushbutton','units','norm','position',[.25 .025 .5 .05],'string','refresh');
        h.info=info; h.tfiles=tfiles; set(h.hfig,'userdata',h);
        set([h.files, h.types, h.refresh],'callback',@(varargin)conn_projectmanager_update_details(h.hfig));
        if ishandle(h.hfig), conn_projectmanager_update_details(h.hfig); end
        
    case 'hpc_push'
%         clear h;
%         h.hfig=figure('units','norm','position',[.3 .4 .5 .2],'name','file transfer (scp)','numbertitle','off','menubar','none','color','w');
%         uicontrol('style','text','units','norm','position',[.1 .70 .19 .15],'string','copy FROM:','backgroundcolor','w','horizontalalignment','right','parent',h.hfig);
%         h.filesFROM=uicontrol('style','edit','max',1,'units','norm','position',[.3 .70 .6 .15],'string','','backgroundcolor','w','horizontalalignment','left','parent',h.hfig);
%         uicontrol('style','text','units','norm','position',[.1 .50 .19 .15],'string','TO:','backgroundcolor','w','horizontalalignment','right','parent',h.hfig);
%         h.filesTO=uicontrol('style','edit','max',1,'units','norm','position',[.3 .50 .6 .15],'string','','backgroundcolor','w','horizontalalignment','left','parent',h.hfig);
%         uicontrol('style','pushbutton','string','OK','units','norm','position',[.1,.01,.38,.25],'callback','uiresume');
%         uicontrol('style','pushbutton','string','Cancel','units','norm','position',[.51,.01,.38,.25],'callback','delete(gcbf)');
        filelocal=varargin{1};
        fileremote=varargin{2};
        ok=system(sprintf('scp -C -q -o ControlPath=''%s'' ''%s'' %s:''%s''', params.info.filename_ctrl,filelocal,params.info.login_ip,fileremote));
        varargout={isequal(ok,0)}; 
        
    case 'hpc_pull'
        fileremote=varargin{1};
        filelocal=varargin{2};
        [ok,msg]=system(sprintf('scp -C -q -o ControlPath=''%s'' %s:''%s'' ''%s''', params.info.filename_ctrl,params.info.login_ip,fileremote,filelocal));
        if ~isequal(ok,0), disp(msg); end
        varargout={isequal(ok,0)}; 
        
    case 'disconnect' % disconnect
        fprintf('Disconnecting from server\n');
        conn_tcpip('close');
        params.state='off';
        conn_cache clear;
        conn_jobmanager clear;
        
    case 'exit' % disconnect & close server
        fprintf('Exiting remote CONN session and disconnecting\n');
        if conn_server('isconnected')
            conn_tcpip('write','exit');
        end
        conn_tcpip('close');
        params.state='off';
        conn_cache clear;
        conn_jobmanager clear;
        
    case 'restart'
        if params.isserver
            fprintf(' Restarting connection...\n');
            connection=conn_tcpip('private');
            ok=false;
            for ntry=1:3
                try
                    conn_tcpip('open','server',connection.port,connection.id,0);
                    ok=true;
                    break;
                end
                pause(10+rand);
            end
            if ~ok, fprintf(' Restart failed. Please restart connection manually\n'); end
        else
            fprintf(' Restarting connection...\n');
            if isfield(params.info,'filename_ctrl')&&~isempty(params.info.filename_ctrl)
                try
                    [tok,tmsg]=system(sprintf('ssh -o ControlPath=''%s'' -O check %s', params.info.filename_ctrl,params.info.login_ip));
                    if tok~=0, system(sprintf('ssh -f -N -o ControlMaster=yes -o ControlPath=''%s'' %s', params.info.filename_ctrl,params.info.login_ip)); end
                    system(sprintf('ssh -o ControlPath=''%s'' -O forward -L%d:%s:%d %s', params.info.filename_ctrl,params.info.local_port,params.info.remote_ip,params.info.remote_port,params.info.login_ip));
                end
            end
            connection=conn_tcpip('private');
            ok=false;
            for ntry=1:3
                try
                    conn_tcpip('open','client',connection.ip,connection.port,connection.id,0);
                    ok=true;
                    break;
                end
                pause(10+rand);
            end
            if ~ok, fprintf(' Restart failed. Please restart connection manually\n'); end
        end
        
    case 'isserver'
        if numel(varargin)>=1, params.isserver=varargin{1}; end
        varargout={params.isserver};
        
    case 'state'
        if numel(varargin)>=1, params.state=varargin{1}; end
        varargout={params.state};
        
    case 'ping'
        if nargout>0
            t1=clock;
            conn_tcpip('write','ping');
            try, pong=conn_tcpip('read');
            catch, pong='-';
            end
            t2=etime(clock,t1);
            str=sprintf('sent ping received %s from %s:%d: round-trip time = %d ms\n',pong,conn_tcpip('private.ip'),conn_tcpip('private.port'),ceil(1000*t2));
            varargout={str};
        else
            for niter=1:10
                if niter>1, pause(rand); end
                t1=clock;
                conn_tcpip('write','ping');
                try, pong=conn_tcpip('read');
                catch, pong='-';
                end
                t2=etime(clock,t1);
                fprintf('sent ping received %s from %s:%d: round-trip time = %d ms\n',pong,conn_tcpip('private.ip'),conn_tcpip('private.port'),ceil(1000*t2));
            end
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

function conn_projectmanager_update_details(hfig,varargin) 
if nargin<1, hfig=gcbf; end
h=get(hfig,'userdata');
i=get(h.files,'value');
j=get(h.types,'value');
switch(j)
    case {1,2,3,4},
        try,
            str=conn_fileutils('fileread',h.tfiles{j}{i});
            str=char(str(:)');
            b=find(diff([0 str==8 0]));
            for n=1:2:numel(b)-1,
                str(max(1,b(n)-(b(n+1)-b(n))):b(n+1)-1)=0;
            end
            str=str(str~=0);
            str=regexp(str,'[\r\n]','split');
        catch
            str={' '};
        end
    case 5,
        if isfield(h.info,'submitcmd')&&numel(h.info.submitcmd)>=i, str=regexp(h.info.submitcmd{i},'[\r\n]','split');
        else str={' '};
        end
    case 6,
        if isfield(h.info,'submitmsg')&&numel(h.info.submitmsg)>=i, str=regexp(h.info.submitmsg{i},'[\r\n]','split');
        else str={' '};
        end
    case 7,
        if isfield(h.info,'statemsg')&&numel(h.info.statemsg)>=i, str=regexp(h.info.statemsg{i},'[\r\n]','split');
        else str={' '};
        end
end
set(h.str,'string',str,'value',numel(str),'listboxtop',numel(str));
%uiwait(h.hfig);
end