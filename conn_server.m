function varargout = conn_server(option, varargin)
% CONN_SERVER server over TCP/IP
%
% COMMANDS FROM SERVER SIDE: (e.g. from a computer at work/university where the data lives)
%   conn_server('start' [,PORT, PUBLICKEY])          : starts CONN server and prints out instructions for connecting to this server (only a single client connection will be permitted)
%                                                     PORT : port number (default []). Leaving PORT empty will bind to first available port
%                                                     PUBLICKEY : hashed password (default []). If a non-empty PUBLICKEY string is provided the server will reject any client connection that
%                                                       cannot provide the associated PRIVATEKEY (note: use [PUBLICKEY,PRIVATEKEY]=conn_tcpip('keypair') to generate a one-time-use KEY pair)
%                                                     e.g.
%                                                       >> [PUBLICKEY,PRIVATEKEY]=conn_tcpip('keypair')
%                                                             PUBLICKEY =
%                                                                   '137227a467c8e62c85e8f3d91767af1a'
%                                                             PRIVATEKEY =
%                                                                   'e5823002ac891dacbf3c48ae54d8f438'
%                                                       >> conn_server('start', [], '137227a467c8e62c85e8f3d91767af1a');
%                                                             Opening port 60039...
%                                                             ******************************************************************************
%                                                             To connect to this server, use the Matlab syntax:
%                                                                conn_server connect mbp-2018.lan 60039CONNprivatekey
%
%                                                             To connect to this server using ssh-tunneling, use the following syntax instead:
%                                                                $(OS-command)     :   ssh -L 6111:alfonsosmbp2018.lan:60039 alfnie@YOUR-INSTITUTION-SSH-LOGIN-NODE
%                                                                >(Matlab-command) :   conn_server connect localhost 6111CONNprivatekey
%                                                             ******************************************************************************%
% COMMANDS FROM CLIENT SIDE: (e.g. from a computer at home/office where you work)
%   conn_server('connect', IP , SERVER_ID)          : stablishes connection to CONN server
%                                                     IP : address of computer running conn_server
%                                                     SERVER_ID: keyword of the form <PORT>CONN<PRIVATEKEY>
%                                                     e.g.
%                                                       >> conn_server('connect', 'mbp-2018.lan', '60039CONNe5823002ac891dacbf3c48ae54d8f438');
%                                                            Connecting to 127.0.0.1:60039...
%                                                            Succesfully established connection to server
%
%   conn_server('run', fcn, arg1, arg2, ...)        : runs on server the function fcn(arg1,arg2,...) and wait for its termination, e.g.
%                                                        >> tic; conn_server('run','pause',10); toc
%                                                        Elapsed time is 10.019476 seconds.
%   [var1, var2,...]=conn_server('run',...)         : as 'run' but also collecting back the output(s) of the fcn call, e.g.
%                                                        >> a = conn_server('run','fft',1:10)
%                                                        a =
%                                                          55.0000 + 0.0000i  -5.0000 +15.3884i  -5.0000 + 6.8819i  -5.0000 + 3.6327i  -5.0000 + 1.6246i  -5.0000 + 0.0000i  -5.0000 - 1.6246i  -5.0000 - 3.6327i  -5.0000 - 6.8819i  -5.0000 -15.3884i
%                                                        >> b = conn_server('run','sum',a)
%                                                        b =
%                                                          10.0000 - 0.0000i
%   conn_server('run_immediatereturn',...)          : as 'run' but without waiting for the command to finish on the server (or check for possible error messages), e.g.
%                                                        >> tic; conn_server('run_immediatereturn','pause',10); toc
%                                                        Elapsed time is 0.006219 seconds.
%   varLink = conn_server('run_keep',...)           : as 'run' but keeping the result in server, returning only a link/label to the result.
%                                                     note: variable links can be used in subsequent conn_server 'run' or 'run_keep' commands, e.g.
%                                                        >> a = conn_server('run_keep','fft',1:10)
%                                                        a =
%                                                         struct with fields: conn_server_variable: 'var_a5f0baeb5e8d3d4fc18dae3553a61c8b_1'
%                                                        >> b = conn_server('run','sum',a)
%                                                        b =
%                                                          10.0000 - 0.0000i
%                                                        >> conn_server('clear',a);
%                                                     note: use varLink = conn_server('run_keepas',label,...) to fix the remote variable link/label (e.g. to save memory by re-using the same remote variable-space), e.g.
%                                                        >> a = conn_server('run_keepas','a','fft',1:10)
%                                                        a =
%                                                         struct with fields: conn_server_variable: 'labeled_a_1'
%                                                        >> a = conn_server('run_keepas','a','abs',a)
%                                                        a =
%                                                         struct with fields: conn_server_variable: 'labeled_a_1'
%                                                        >> b = conn_server('run','sum',a)
%                                                        b =
%                                                          132.2490
%                                                        >> conn_server('clear',a);
%   conn_server cmd                                 : interactive command-line execution in remote CONN server
%   [...]=conn_server('cmd',...)                    : single-command execution in remote CONN server, e.g.
%                                                   	 >> y=conn_server('cmd','fft(1:10)');
%                                                        y =
%                                                          55.0000 + 0.0000i  -5.0000 +15.3884i  -5.0000 + 6.8819i  -5.0000 + 3.6327i  -5.0000 + 1.6246i  -5.0000 + 0.0000i  -5.0000 - 1.6246i  -5.0000 - 3.6327i  -5.0000 - 6.8819i  -5.0000 -15.3884i
%
%   conn_server('push',file_local,file_remote)      : copies file from client to server
%   conn_server('pull', file_remote,file_local)     : copies file from server to client
%   conn_server('ping')                             : pings the server and shows round-trip times
%   conn_server('isconnected')                      : tests bidirectional transfer and returns true if connection is alive
%   conn_server('disconnect')                       : disconnects from CONN server (a new connection with the same server can be stablished using "conn_server('connect',...)" again from this or other machine)
%   conn_server('exit')                             : stops the remote CONN server and disconnects (no new connections are allowed until the server is restarted, e.g. using "conn_server('start',...)" again from the server side)
%



persistent params
persistent local_vars
if ~nargin||isempty(option), option='start'; end

varargout={};
if isempty(params)
    params=struct(...
        'isserver',false,...
        'state','off');
end
if isempty(local_vars), local_vars=struct; end

switch(lower(option))
    case {'start','startpersistent'} % init server
        params.isserver=true;
        disphelp=2;
        if numel(varargin)>=1&&~isempty(varargin{1}), port=varargin{1}; else port=0; end
        if numel(varargin)>=2&&~isempty(varargin{2}), id=char(varargin{2}); else id=[]; end
        if numel(varargin)>=3&&isequal(varargin{3},'silent'), varargin{3}=[]; disphelp=0; end
        if numel(varargin)>=3&&~isempty(varargin{3}), disphdl=varargin{3}; else disphdl=[]; end
        ok=false;
        while ~ok
            try
                conn_tcpip('open','server',port,id,disphelp);
                params.state='on';
                ok=true;
            catch me
                disp(me.message)
            end
            pause(rand);
            if ~isempty(disphdl)&&(any(~ishandle(disphdl))||(numel(disphdl)>1&&get(disphdl(2),'value')>0)), return; end
        end
        conn_cache clear;
        conn_jobmanager clear;
        if strcmpi(option,'startpersistent'), conn_server('continuepersistent',varargin{3:end});
        else conn_server('continue',varargin{3:end});
        end
        
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
        
    case {'run','run_immediatereturn','run_withwaitbar','run_keep','run_keepas'}
        data.type='run';
        data.id=char(conn_tcpip('hash',mat2str(now)));
        if strcmpi(option,'run_withwaitbar'), statushandle=varargin{1}; data.cmd=varargin(2:end); data.withwaitbar=true; data.keeplocal=false;
        elseif strcmpi(option,'run_keepas'), statushandle=[]; data.cmd=varargin(2:end); data.withwaitbar=true; data.keeplocal=varargin{1};
        else statushandle=[]; data.cmd=varargin; data.withwaitbar=false; data.keeplocal=strcmpi(option,'run_keep');
        end
        data.nargout=nargout;
        data.immediatereturn=strcmpi(option,'run_immediatereturn');
        %data.hpc=false;
        %if isfield(params.info,'host')&&~isempty(params.info.host)&&isfield(params.info,'scp')&&params.info.scp>0, data.hpc=true; end
        conn_tcpip('flush');
        conn_tcpip('write',data);
        if ~data.immediatereturn;
            ok=false;
            while ~ok
                try
                    var=conn_tcpip('read');
                catch me,
                    if ~isempty(regexp(me.message,'SocketTimeoutException')), fprintf('.');  % timeout
                    elseif ~isempty(regexp(me.message,'EOFException|IOException|SocketException'))
                        error('ERROR: connection may be down\n');
                        %% restart
                        %fprintf('\n Idle connection to server. ');
                        %conn_server restart;
                        %conn_tcpip('write',data);
                    else error('ERROR: connection problem %s',me.message);
                    end
                    var=struct('type','null','id',data.id);
                end
                if isfield(var,'type')&&isequal(var.type,'status')
                    str=regexprep(var.msg,'[\r\n\s]*',' ');
                    if data.withwaitbar,
                        conn_disp('__nolog',str);
                        if ~isempty(statushandle),
                            try
                                set(statushandle(1),'string',str); drawnow;
                                if numel(statushandle)>1&&get(statushandle(2),'value')>0
                                    conn_tcpip('poke','STOP');
                                end
                            end
                        end
                    else
                        if iscell(str), for n=1:numel(str), fprintf('%s\n',str{n}); end
                        else fprintf('%s',str);
                        end
                    end
                elseif isempty(var)||~isstruct(var)||~isfield(var,'type')||~isfield(var,'id')||~isequal(var.id, data.id),
                    disp(var);
                    fprintf('WARNING: unexpected response\n');
                else
                    %                     if isequal(var.type,'ok_hasattachment') % server is waiting for us to scp its response
                    %                         tmpfile1=fullfile(conn_cache('private.local_folder'),['cachetmp_', char(conn_tcpip('hash',mat2str(now))),'.mat']);
                    %                         tmpfile2=var.msg;
                    %                         if conn_server_ssh('pull',tmpfile2,tmpfile1), var=load(tmpfile1,'arg'); var=var.arg;
                    %                         else var=struct('type','ko','id',var.id,'msg','unable to run SSH_pull to read response from server');
                    %                         end
                    %                         conn_fileutils('deletefile',tmpfile1);
                    %                         conn_tcpip('write',struct('type','ack_hasattachment','id',var.id)); % signal server to continue
                    %                     end
                    if isequal(var.type,'ok')
                        if isfield(var,'msg'), varargout=var.msg;
                        else varargout={};
                        end
                        ok=true;
                    elseif isequal(var.type,'ko')
                        if isfield(var,'msg'), error('ERROR at server: %s\n',var.msg);
                        else error('ERROR at server\n');
                        end
                        %ok=true;
                    end
                end
            end
        end
        
    case {'continue','continuepersistent','test'}
        if params.isserver % server
            conn_disp('__portcomm',true);
            try
                % listening
                if numel(varargin)>=1&&~isempty(varargin{1}), disphdl=varargin{1}; else disphdl=[]; end
                continueonexit=strcmpi(option,'continuepersistent');
                dispstr=sprintf(' CONN SERVER ACTIVE %s',datestr(now));
                conn_server_fprintf(disphdl,'fprintf','%s',dispstr);
                ntimedout=0; % minutes
                while 1
                    data=[];
                    drawnow;
                    try,
                        data=conn_tcpip('read');
                        ntimedout=0;
                    catch me,
                        if ~isempty(regexp(me.message,'EOFException|IOException|SocketException')) %||ntimedout>15 % restart
                            dispstr='';
                            conn_server_fprintf(disphdl,'fprintf','\n Idle connection to client.');
                            if (~isempty(disphdl)&&(any(~ishandle(disphdl))||(numel(disphdl)>1&&get(disphdl(2),'value')>0)))||~continueonexit
                                conn_server_fprintf(disphdl,'fprintf',' Closing server.\n');
                                conn_tcpip('close');
                                conn_disp('__portcomm',false);
                                return
                            else
                                conn_server('restart',disphdl);
                            end
                        elseif ~isempty(regexp(me.message,'SocketTimeoutException')) % timeout
                            conn_server_fprintf(disphdl,'fprintf',repmat('\b',[1,length(dispstr)]));
                            dispstr=sprintf(' CONN SERVER ACTIVE %s',datestr(now));
                            conn_server_fprintf(disphdl,'fprintf',dispstr);
                            conn_tcpip('writekeepalive');
                            ntimedout=ntimedout+1;
                        else
                            conn_server_fprintf(disphdl,me.message);
                            pause(rand);
                        end
                    end
                    if isempty(data)
                    elseif isequal(data,'handshake')||isequal(data,'restart')||(isequal(data,'exit')&&continueonexit)
                        conn_server('restart',disphdl);
                    elseif isequal(data,'exit')||isequal(data,'exitforce')
                        fprintf('\n Server closed by client\n'); dispstr='';
                        conn_tcpip('close');
                        conn_disp('__portcomm',false);
                        return
                    elseif ~isempty(disphdl)&&(any(~ishandle(disphdl))||(numel(disphdl)>1&&get(disphdl(2),'value')>0))
                        fprintf('\n Server closed\n'); dispstr='';
                        conn_tcpip('close');
                        conn_disp('__portcomm',false);
                        return
                    elseif isequal(data,'ping')
                        conn_tcpip('write','ok');
                    elseif isstruct(data)&&numel(data)==1&&isfield(data,'cmd') % run command
                        if ischar(data.cmd), data.cmd={data.cmd}; end
                        try, data.cmd{1}=regexprep(data.cmd{1},'\.m$',''); end
                        if ischar(data.cmd{1})||conn_server('util_isremotevar',data.cmd{1}), %isequal(data.cmd{1},'conn')||~isempty(regexp(data.cmd{1},'^conn_'))||isequal(data.cmd{1},'spm')||~isempty(regexp(data.cmd{1},'^spm_')) % run only conn or spm commands
                            conn_server_fprintf(disphdl,'fprintf','\n -'); dispstr=''; timer=[];
                            try
                                conn_server_fprintf(disphdl,data.cmd);
                                if conn_server('util_isremotevar',data.cmd{1}), fh=@(x)x; data.cmd=[{''},data.cmd];
                                elseif isempty(data.cmd{1}), fh=@(x)x;
                                else fh=eval(sprintf('@%s',data.cmd{1}));
                                end
                                is_server_variable=find(conn_server('util_isremotevar',data.cmd(2:end)));
                                for nvar=reshape(is_server_variable,1,[]) % use variables defined using run_keep
                                    if numel(data.cmd)>2&&isequal(data.cmd{1},'conn_server')&&isequal(data.cmd{2},'clear'), data.cmd{1+nvar}=data.cmd{1+nvar}.conn_server_variable;
                                    elseif isfield(local_vars,data.cmd{1+nvar}.conn_server_variable), data.cmd{1+nvar}=local_vars.(data.cmd{1+nvar}.conn_server_variable);
                                    else data.cmd{1+nvar}=[]; % error, variable does not exist
                                    end
                                end
                                %                             if isfield(data,'displaystatus')&&data.displaystatus,
                                %                                 htimer=timer('name','jobmanager','startdelay',1,'period',10,'executionmode','fixedspacing','taskstoexecute',inf,'busymode','drop','timerfcn',@conn_server_update);
                                %                                 start(htimer);
                                %                             end
                                if ~isfield(data,'immediatereturn')||~data.immediatereturn
                                    if isfield(data,'nargout')&&data.nargout>0,
                                        varout=cell(1,data.nargout);
                                        try
                                            [varout{1:data.nargout}]=feval(fh,data.cmd{2:end});
                                            argout=struct('type','ok','id',data.id,'msg',{varout});
                                        catch me
                                            str=conn_errormessage(me);
                                            str=sprintf('%s\n',str{:});
                                            argout=struct('type','ko','id',data.id,'msg',str);
                                        end
                                    else
                                        try
                                            feval(fh,data.cmd{2:end});
                                            argout=struct('type','ok','id',data.id);
                                        catch me
                                            str=conn_errormessage(me);
                                            str=sprintf('%s\n',str{:});
                                            argout=struct('type','ko','id',data.id,'msg',str);
                                        end
                                    end
                                    conn_server_fprintf(disphdl,argout);
                                    if isequal(argout.type,'ok')&&isfield(argout,'msg')&&isfield(data,'keeplocal')&&all(data.keeplocal>0)
                                        for nvar=1:numel(argout.msg)
                                            if ischar(data.keeplocal), varname=['labeled_',data.keeplocal,'_',num2str(nvar)];
                                            else varname=['var_',char(conn_tcpip('hash',mat2str(now))),'_',num2str(nvar)];
                                            end
                                            local_vars.(varname)=argout.msg{nvar};
                                            argout.msg{nvar}=struct('conn_server_variable',varname);
                                        end
                                        conn_tcpip('write',argout);
                                        %                                 elseif isequal(argout.type,'ok')&&isfield(argout,'msg')&&isfield(data,'hpc')&&data.hpc>0&&getfield(whos('argout'),'bytes')>1e6, % send response using ssh/scp?
                                        %                                     tmpfile=fullfile(conn_cache('private.local_folder'),['cachetmp_', char(conn_tcpip('hash',mat2str(now))),'.mat']);
                                        %                                     arg=argout; save(tmpfile,'arg');
                                        %                                     argout=struct('type','ok_hasattachment','id',data.id,'msg',tmpfile);
                                        %                                     conn_tcpip('write',argout);
                                        %                                     while 1,
                                        %                                         rsp=[];
                                        %                                         try, rsp=conn_tcpip('read'); end % wait for client to finish
                                        %                                         if isfield(rsp,'id')&&isequal(rsp.id,data.id)&&isfield(rsp,'type')&&isequal(rsp.type,'ack_hasattachment'), break; end
                                        %                                         fprintf('.'); pause(rand);
                                        %                                     end
                                        %                                     conn_fileutils('deletefile',tmpfile);
                                    else % send response through communication socket
                                        conn_tcpip('write',argout);
                                    end
                                else
                                    try
                                        feval(fh,data.cmd{2:end});
                                    end
                                end
                            catch me
                                str=conn_errormessage(me);
                                str=sprintf('%s\n',str{:});
                                argout=struct('type','ko','id',data.id,'msg',str);
                            end
                            %                         if ~isempty(htimer), try, stop(htimer); delete(htimer); end; end
                            
                        elseif ~isfield(data,'immediatereturn')||~data.immediatereturn
                            argout=struct('type','ko', 'id', data.id, 'msg', sprintf('unrecognized option %s ',data.cmd{1}));
                            conn_tcpip('write',argout);
                        end
                    else % testing (obsolete)
                        conn_server_fprintf(disphdl,'fprintf','\n received test data:\n');  dispstr='';
                        disp(data)
                    end
                end
            catch me
                conn_disp('__portcomm',false);
                str=regexprep(char(getReport(me,'extended','hyperlinks','off')),'[\t ]+',' ');
                conn_server_fprintf(disphdl,'fprintf','\nERROR (terminating server)\n%s',str);
            end
            conn_disp('__portcomm',false);
            
        else % client testing (obsolete)
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
        filename1=conn_server('util_localfile',filename1);
        ok=false(1,numel(filename1));
        for n=1:numel(filename1)
            if conn_existfile(filename1{n}),
                %if isfield(params.info,'host')&&~isempty(params.info.host)&&isfield(params.info,'scp')&&params.info.scp>0
                %    ok(n)=conn_server_ssh('push',filename1{n},filename2{n});
                %else
                conn_server('run_immediatereturn','conn_tcpip','readtofile',filename2{n});
                ok(n)=conn_tcpip('writefromfile',filename1{n});
                %end
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
        %if isfield(params.info,'host')&&~isempty(params.info.host)&&isfield(params.info,'scp')&&params.info.scp>0
        %    for n=1:numel(filename1)
        %        ok(n)=conn_server_ssh('pull',filename1{n},filename2{n});
        %        if ok(n), hash{n}=conn_cache('hash',filename2{n}); end
        %    end
        %else
        conn_tcpip('flush');
        for n=1:numel(filename1)
            conn_server('run_immediatereturn','conn_tcpip','writefromfile',filename1{n});
        end
        for n=1:numel(filename1)
            hash{n}=conn_tcpip('readtofile',filename2{n});
            ok(n)=~isempty(hash{n});
        end
        %end
        assert(ok(1), 'communication with server interrupted during file pull');
        varargout=hash;
        
    case {'cmd','command'}
        singlecommand=nargout>0|~isempty(varargin);
        tnameserver=conn_tcpip('private.ip');
        while 1
            if ~isempty(varargin), cmd=varargin{1}; opts=varargin(2:end);
            else
                fprintf('%s',[tnameserver,' >> ']);
                cmd=input('','s'); 
                opts={};
            end
            switch(lower(cmd))
                case {'quit','exit'}, break;
                otherwise
                    try
                        if singlecommand,
                            [varargout{1:nargout}]=conn_server('run','conn_server_cmd',cmd,opts{:});
                            %[varargout{1:nargout}]=conn_server('run','conn','cmd',cmd,opts{:});
                        else
                            str=conn_server('run','conn_server_cmd_capture',cmd);
                            %str=conn_server('run','conn','cmd_capture',cmd);
                            disp(str);
                        end
                    catch me
                        str=regexprep(char(getReport(me,'extended','hyperlinks','off')),'[\t ]+',' ');
                        disp(str);
                    end
            end
            if singlecommand, break; end
        end
        
    case 'clear_cache'
        conn_cache clearall;
        %conn_tcpip clear;
        
    case {'clear','clear_keep'}
        if params.isserver,
            if strcmpi(option,'clear_keep')||isempty(varargin)||isequal(varargin,{'all'})
                local_vars=struct;
            else
                local_vars=rmfield(local_vars,varargin);
            end
        else
            try, conn_server('run_immediatereturn','conn_server','clear',varargin{:}); end
        end
        
    case 'disconnect' % disconnect
        fprintf('Disconnecting from server\n');
        conn_tcpip('close');
        params.state='off';
        conn_cache clear;
        conn_jobmanager clear;
        
    case 'exit' % disconnect & close server
        if params.isserver
            conn_tcpip('close');
            conn_disp('__portcomm',false);
        else
            fprintf('Exiting remote CONN session and disconnecting\n');
            if conn_server('isconnected')
                conn_tcpip('write','exit');
            end
            conn_tcpip('close');
            params.state='off';
            conn_cache clear;
            conn_jobmanager clear;
        end
        
    case 'restart'
        if params.isserver
            if numel(varargin)>=1&&~isempty(varargin{1}), disphdl=varargin{1}; else disphdl=[]; end
            fprintf(' Restarting connection...\n');
            connection=conn_tcpip('private');
            ok=false;
            for ntry=1:3
                try
                    conn_tcpip('open','server',connection.port,connection.id,0);
                    ok=true;
                    break;
                end
                if ~isempty(disphdl)&&(any(~ishandle(disphdl))||(numel(disphdl)>1&&get(disphdl(2),'value')>0)), return; end
                pause(10+rand);
            end
            if ~ok, fprintf(' Restart failed. Please restart connection manually\n'); end
        else
            fprintf(' Restarting connection...\n');
            %if ~isempty(params.info.host)&&isfield(params.info,'filename_ctrl')&&~isempty(params.info.filename_ctrl)
            %    try
            %        [tok,tmsg]=system(sprintf('%s -o ControlPath=''%s'' -O check %s', params.options.cmd_ssh, params.info.filename_ctrl,params.info.login_ip));
            %        if tok~=0, system(sprintf('%s -f -N -o ControlMaster=yes -o ControlPath=''%s'' %s', params.options.cmd_ssh, params.info.filename_ctrl,params.info.login_ip)); end
            %        system(sprintf('%s -o ControlPath=''%s'' -O forward -L%d:%s:%d %s', params.options.cmd_ssh, params.info.filename_ctrl,params.info.local_port,params.info.remote_ip,params.info.remote_port,params.info.login_ip));
            %    end
            %end
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
            if isequal(pong,'ok'), failed=''; else failed='FAILED: '; end
            str=sprintf('%ssent ping received %s from %s:%d: round-trip time = %d ms\n',failed,pong,conn_tcpip('private.ip'),conn_tcpip('private.port'),ceil(1000*t2));
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
            conn_tcpip('write','ping');
            assert(isequal(conn_tcpip('read'),'ok'));
            params.state='on';
            varargout={true};
        catch
            params.state='off';
            varargout={false};
        end
        
    case 'util_who'
        if params.isserver
            names=fieldnames(local_vars);
        else
            names=conn_server('run','conn_server','util_who');
        end
        if nargout>0, varargout={names};
        else disp(char(names));
        end
        
    case 'util_isremotevar',
        var=varargin{1};
        if ~iscell(var), var={var}; end
        out=false(size(var));
        for nvar=1:numel(var), if isstruct(var{nvar})&&isfield(var{nvar},'conn_server_variable'), out(nvar)=true; end; end
        varargout={out};
        
    case 'util_cleanremotevar'
        var=varargin{1};
        if iscell(var), var=cellfun(@(var)conn_server(option,var),var,'uni',0);
        elseif isstruct(var)&&isfield(var,'conn_server_variable'), var=struct('conn_server_variable',var.conn_server_variable);
        end
        varargout={var};
        
    case 'util_isremotefile'
        varargout={false};
        try, varargout={cellfun(@(x)~isempty(regexp(x,'^\s*[\\\/]CONNSERVER')),cellstr(varargin{1}))};
        end
        
    case 'util_localfile'
        varargout=varargin;
        try
            ischarfilename=ischar(varargin{1});
            filename=regexprep(cellstr(varargin{1}),'^\s*[\\\/]CONNSERVER','');
            filename=regexprep(filename,'\\|\/',['\',filesep]); % localfiles use convention of this computer
            if ischarfilename, filename=char(filename); end
            varargout={filename};
        end
        
    case 'util_localformat'
        varargout=varargin;
        try
            ischarfilename=ischar(varargin{1});
            filename=regexprep(cellstr(varargin{1}),'\\|\/',['\',filesep]); % localfiles use convention of this computer
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
            filename(change)=regexprep(filename(change),'^\s*[\\\/]*(.*)',['\',filesep,'CONNSERVER','\',filesep,'$1']);
            filename(change)=regexprep(filename(change),'\\|\/','/'); % remotefiles use common / convention
            if ischarfilename, filename=char(filename); end
            varargout={filename};
        end
        
    otherwise,
        error('unrecognized option %s',option);
end
end

function varargout=conn_server_cmd(varargin)
if numel(varargin)>1, [varargout{1:nargout}]=feval(varargin{:});
else [str,varargout{1:nargout}]=evalc(sprintf('evalin(''base'',''%s'')',regexprep(varargin{1},'''',''''''))); disp(str);
end
end

function varargout=conn_server_cmd_capture(varargin)
[varargout{1:nargout}]=evalc(sprintf('evalin(''base'',''%s'')',regexprep(varargin{1},'''','''''')));
disp(varargout{1});
end

function endnow=conn_server_fprintf(disphdl, varargin)
endnow=false;
try
    if ~isequal(varargin{1},'fprintf'), newstr=regexp(evalc('disp(varargin{1})'),'[\r\n]+','split');
    else newstr=regexp(sprintf(varargin{2:end}),'[\r\n]+','split');
    end
    newstr=newstr(cellfun('length',newstr)>0);
    if isempty(disphdl)
        disp(char(newstr));
    else
        if any(~ishandle(disphdl)), endnow=true;
        else
            str=cellstr(get(disphdl(1),'string'));
            str=[str(:);newstr(:)];
            str=str(max(1,numel(str)-10):end);
            set(disphdl(1),'string',str);
            if numel(disphdl)>1, endnow=get(disphdl(2),'value'); end
        end
    end
catch
end
end