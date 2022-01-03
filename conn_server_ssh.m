function varargout = conn_server_ssh(option, varargin)
% internal function

%
% alternative procedure for SSH-accessible networks only:
% command from server side: (from a computer within the SSH-accessible network)
%   conn_server_ssh('save' [, fileout])             : server configuration (run only once in order to configure this computer), creates a .json file with basic information necessary to connect to this server (machine IP, and where to find Matlab/SPM/CONN) [~/connserverinfo.json]
%
% commands from client side: (from a computer with a SSH client)
%   conn_server_ssh('start' [, IP])                 : SSH to network login node, submits a job to start a CONN server, and connects this computer to that remote server (note: expects to find IP:~/connserverinfo.json file; run "conn_server_ssh setup" in server once)
%   conn_server_ssh('start' [, filein])             : same as above but loading CONN server information explicitly from the provided .json file
%   conn_server_ssh('submitstart', P, K)            : submits a job to start a new CONN server (this function is run internally by SSH_start)
%                                                     P: profile name used to submit jobs (e.g. 'background'; see "conn_jobmanager profiles")
%                                                     K: PUBLICKEY string
%   conn_server_ssh('restart')                      : restarts a dropped connection
%   conn_server_ssh('exit')                         : terminates CONN server and disconnects
%


persistent params
if ~nargin||isempty(option), option='start'; end

varargout={};
if isempty(params)
    %'isserver',false,...
    params=struct(...
        'info',struct(),...
        'options',struct(),...
        'state','off');
    filename=fullfile(conn_fileutils('homedir'),'connclientinfo.json');
    if conn_existfile(filename), params.options=conn_jsonread(filename); else params.options=struct; end
    if ~isfield(params.options,'use_ssh'), params.options.use_ssh=true; end
    if ~isfield(params.options,'cmd_ssh'), params.options.cmd_ssh='ssh'; end
    if ~isfield(params.options,'cmd_scp'), params.options.cmd_scp='scp'; end
end

switch(lower(option))
    
    case {'start','restart'} % init server remotely and connect to it
        % development reference notes on ssh tunneling
        %    $ ssh -fN -o ServerAliveInterval=60 -o ServerAliveCountMax=10 -o ControlMaster=yes -o ControlPath=<local_filename> <login_node>        % authenticate first
        %    $ ssh -o ControlPath=<local_filename> -L<local_port>:<server_ip>:<server_port> <login_node>                                           % port forwarding on shared connection
        %    $ ssh -o ControlPath=<local_filename> -O forward -L<local_port>:<server_ip>:<server_port> <login_node>
        %    $ ssh -o ControlPath=<local_filename> -O check <login_node>                                                                             % exit/close connection
        %    $ ssh -o ControlPath=<local_filename> -O exit <login_node>                                                                             % exit/close connection
        if numel(varargin)>=1,
            filename=varargin{1};
            if isequal(filename,'recent'), filename=fullfile(conn_fileutils('homedir'),'conn_recentservers.json'); end
            if isempty(regexp(filename,'.json$')), params.info.host=filename;
            else params.info=conn_jsonread(filename);
            end
        else params.info.CONNcmd='';
        end
        if params.options.use_ssh, 
            allthesame=true;
            if ~isfield(params.info,'host')||isempty(params.info.host), params.info.host=conn_server_ssh_input('Server address [local]: ','s'); allthesame=false;
            else
                temp=conn_server_ssh_input(sprintf('Server address [%s]: ',params.info.host),'s');
                if ~isempty(temp),
                    temp=regexprep(temp,'\s+','');
                    if ~isequal(params.info.host,temp), allthesame=false; end
                    params.info.host=temp;
                end
            end
            if strcmp(params.info.host,'local')||strcmp(params.info.host,'none'), params.info.host=''; end
        else
            params.info.host=''; 
        end
        if numel(varargin)>=2&&~isempty(varargin{2}), params.info.local_port=varargin{2}; if ischar(params.info.local_port), params.info.local_port=str2double(params.info.local_port); end
        else params.info.local_port=[];
        end
        if ~isfield(params.info,'user')||isempty(params.info.user), [nill,str2]=system('whoami'); params.info.user=regexprep(str2,'\n',''); allthesame=false; end
        if isempty(params.info.host),
            params.info.user='';
            params.info.login_ip='';
            params.info.CONNcmd='';
            params.info.filename_ctrl='';
            allthesame=false;
        else
            temp=conn_server_ssh_input(sprintf('Username [%s]: ',params.info.user),'s');
            if ~isempty(temp),
                if ~isequal(params.info.user,temp), allthesame=false; end
                params.info.user=temp;
            end
            params.info.login_ip=sprintf('%s@%s',params.info.user,params.info.host);
            localcachefolder=conn_cache('private.local_folder');
        end
        params.info.scp=false;
        if isempty(params.options.cmd_ssh), 
            error('No SSH client found (see Tools.RemoteOptions.Configuration)');
        else
            if ~isempty(params.info.host)
                params.info.filename_ctrl=fullfile(localcachefolder,sprintf('connserver_ctrl_%s_%s',params.info.host,params.info.user));
                try, conn_fileutils('deletefile',params.info.filename_ctrl); end
                % starts a shared SSH connection
                if ispc
                    fprintf('Initializing connection to %s... ',params.info.login_ip);
                    system(sprintf('%s -f -N -o ControlMaster=yes -o ControlPath=''%s'' %s &', params.options.cmd_ssh,params.info.filename_ctrl,params.info.login_ip),'-echo'); % starts a shared connection
                    ok=1; while ok~=0, pause(1); [ok,msg]=system(sprintf('%s -o ControlPath=''%s'' -O check %s', params.options.cmd_ssh, params.info.filename_ctrl,params.info.login_ip)); end
                elseif 1, 
                    fprintf('Connecting to %s... ',params.info.login_ip);
                    system(sprintf('%s -f -N -o ControlMaster=yes -o ControlPath=''%s'' %s', params.options.cmd_ssh,params.info.filename_ctrl,params.info.login_ip),'-echo'); % starts a shared connection
                    [ok,msg]=system(sprintf('%s -o ControlPath=''%s'' -O check %s', params.options.cmd_ssh, params.info.filename_ctrl,params.info.login_ip)); 
                    if ok~=0, error(msg); end
                else
                    %[msg]=evalc('!ssh -o ControlPath=''/Users/alfnie/.conn_cache/connserver_ctrl_scc1.bu.edu_alfnie'' -O check alfnie@scc1.bu.edu') % skips authentication if ControlMaster already exists
                    %if isempty(regexp(msg,'Master running'))
                    [ok,msg]=system(sprintf('%s -o ControlPath=''%s'' -O check %s', params.options.cmd_ssh, params.info.filename_ctrl,params.info.login_ip)); 
                    fprintf('Connecting to %s... ',params.info.login_ip);
                    if ok~=0,
                        system(sprintf('%s -f -N -o ControlMaster=yes -o ControlPath=''%s'' %s', params.options.cmd_ssh,params.info.filename_ctrl,params.info.login_ip),'-echo'); % starts a shared connection
                        [ok,msg]=system(sprintf('%s -o ControlPath=''%s'' -O check %s', params.options.cmd_ssh, params.info.filename_ctrl,params.info.login_ip));
                        if ok~=0, error(msg); end
                    end
                end
                if 1,%~isfield(params.info,'CONNcmd')||isempty(params.info.CONNcmd) % attempts to load server info from remote ~/connserverinfo.json file
                    filename=fullfile(conn_cache('private.local_folder'),['conncache_', char(conn_tcpip('hash',mat2str(now)))]);
                    conn_fileutils('deletefile',filename);
                    fprintf('\nDownloading configuration information from %s:%s to %s\n',params.info.login_ip,'~/connserverinfo.json',filename);
                    [ok,msg]=system(sprintf('%s -q -o ControlPath=''%s'' %s:~/connserverinfo.json %s',...
                        params.options.cmd_scp, params.info.filename_ctrl,params.info.login_ip,filename));
                    assert(conn_existfile(filename),'unable to find ~/connserverinfo.json file in %s',params.info.login_ip);
                    tjson=conn_jsonread(filename);
                    params.info.CONNcmd=tjson.CONNcmd;
                    if isfield(tjson,'SERVERcmd'), params.info.SERVERcmd=tjson.SERVERcmd; end
                    if isfield(tjson,'SERVERpersistent'), params.info.SERVERpersistent=tjson.SERVERpersistent; end
                end
            end
            ntries=1;
            startnewserver=false;
            if allthesame ...
                    &&isfield(params.info,'remote_ip')&&~isempty(params.info.remote_ip) ...
                    &&isfield(params.info,'remote_port')&&~isempty(params.info.remote_port) ...
                    &&isfield(params.info,'remote_id')&&~isempty(params.info.remote_id) ...
                    &&isfield(params.info,'remote_log')&&~isempty(params.info.remote_log)
                params.info.local_port=[];
                ntries=2;
                if ~isfield(params.info,'start_time')||isempty(params.info.start_time), params.info.start_time=datestr(now); end
                fprintf('Attempting to reconnect to last remote CONN session\n');
            elseif strcmpi(option,'restart')
                if ~isfield(params.info,'remote_ip')||isempty(params.info.remote_ip), params.info.remote_ip=conn_server_ssh_input('Remote session host address: ','s'); end
                if ~isfield(params.info,'remote_port')||isempty(params.info.remote_port), params.info.remote_port=str2double(conn_server_ssh_input('Remote session access port: ','s')); end
                if ~isfield(params.info,'remote_id')||isempty(params.info.remote_id), params.info.remote_id=conn_server_ssh_input('Remote session id: ','s'); end
                if ~isfield(params.info,'remote_log')||isempty(params.info.remote_log), params.info.remote_log=conn_server_ssh_input('Remote session log folder: ','s'); end
                if ~isfield(params.info,'start_time')||isempty(params.info.start_time), params.info.start_time=datestr(now); end
                params.info.local_port=[];
            elseif isempty(params.info.host)
                if ~isfield(params.info,'remote_ip'), params.info.remote_ip=''; end; if isempty(params.info.remote_ip), temp=conn_server_ssh_input('Remote session host address: ','s'); else temp=conn_server_ssh_input(sprintf('Remote session host address [%s]: ',params.info.remote_ip),'s'); end; if ~isempty(temp), params.info.remote_ip=temp; end
                if ~isfield(params.info,'remote_port'), params.info.remote_port=[]; end; if isempty(params.info.remote_port), temp=conn_server_ssh_input('Remote session access port: ','s'); else temp=conn_server_ssh_input(sprintf('Remote session access port [%d]: ',params.info.remote_port),'s'); end; if ~isempty(temp), params.info.remote_port=str2double(temp); end
                if ~isfield(params.info,'remote_id'), params.info.remote_id=''; end; if isempty(params.info.remote_id), temp=conn_server_ssh_input('Remote session id: ','s'); else temp=conn_server_ssh_input(sprintf('Remote session id [%s]: ',params.info.remote_id),'s'); end; if ~isempty(temp), params.info.remote_id=deblank(temp); end
                params.info.remote_log=''; %if ~isfield(params.info,'remote_log'), params.info.remote_log=''; end; if isempty(params.info.remote_log), temp=conn_server_ssh_input('Remote session log folder: ','s'); else temp=conn_server_ssh_input(sprintf('Remote session log folder [%s]: ',params.info.remote_log),'s'); end; if ~isempty(temp), params.info.remote_log=deblank(temp); end
                if ~isfield(params.info,'start_time')||isempty(params.info.start_time), params.info.start_time=datestr(now); end
                params.info.local_port=[];
            else startnewserver=true;
            end
            while ntries>0
                if startnewserver
                    fprintf('Requesting a new Matlab session in %s. This may take a few minutes, please be patient as your job currently sits in a queue. CONN will resume automatically when the new Matlab session becomes available\n',params.info.login_ip);
                    [keys_public,keys_private]=conn_tcpip('keypair');
                    % submit jobs to start server#2 in arbitrary remote node using HPC scheduler
                    if isfield(params.info,'SERVERpersistent')&&all(params.info.SERVERpersistent>0), sbc='submitstartpersistent'; else sbc='submitstart'; end
                    if isfield(params.info,'SERVERcmd')&&~isempty(params.info.SERVERcmd), tstr=sprintf('server_ssh %s ''%s'' ''%s''',sbc,params.info.SERVERcmd, keys_public);
                    else tstr=sprintf('server_ssh %s '''' ''%s''', sbc,keys_public);
                    end
                    [ok,msg]=system(sprintf('%s -o ControlPath=''%s'' %s "%s"', params.options.cmd_ssh, params.info.filename_ctrl,params.info.login_ip, regexprep(sprintf(params.info.CONNcmd,tstr),'"','\\"')));
                    if ~isempty(regexp(msg,'SSH_SUBMITSTART error')), error('Error initiating server job\n %s',msg);
                    else
                        try
                            keys=regexp(msg,'HOST:([^\n]+)\nPORT:([^\n]+)\nID:([^\n]+)\nLOG:([^\n]+)\n','tokens','once');
                            params.info.remote_ip=keys{1};
                            params.info.remote_port=str2double(keys{2});
                            params.info.remote_id=keys{3};
                            params.info.remote_log=keys{4};
                            params.info.start_time=datestr(now);
                            fprintf('Remote session started:\n  Host address = %s\n  Access port = %d\n  ID = %s\n  Log folder = %s\n',params.info.remote_ip,params.info.remote_port,params.info.remote_id,params.info.remote_log);
                            params.info.remote_id=keys_private;
                        catch
                            disp(msg)
                            error('unable to start remote CONN server')
                        end
                    end
                end
                if isempty(params.info.host)
                    if strcmpi(option,'restart'), conn_tcpip('open','client',params.info.remote_ip,params.info.remote_port,params.info.remote_id,0); params.state='on';
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
                    [ok,msg]=system(sprintf('%s -o ControlPath=''%s'' -O forward -L%d:%s:%d %s', params.options.cmd_ssh, params.info.filename_ctrl,params.info.local_port,params.info.remote_ip,params.info.remote_port,params.info.login_ip));
                    if ok~=0, 
                        params.info.local_port=[]; 
                    else
                        %system(sprintf('ssh -f -N -T -o ExitOnForwardFailure=yes -o ControlPath=''%s'' -L%d:%s:%d %s', params.info.filename_ctrl,params.info.local_port,params.info.remote_ip,params.info.remote_port,params.info.login_ip));
                        %fprintf('Connecting to server\n');
                        if strcmpi(option,'restart'), conn_tcpip('open','client','localhost',params.info.local_port,params.info.remote_id,0); params.state='on';
                        else conn_server('connect','localhost',sprintf('%dCONN%s',params.info.local_port,params.info.remote_id));
                        end
                    end
                end
                ntries=max(0,ntries-1);
                filename=fullfile(conn_fileutils('homedir'),'conn_recentservers.json');
                spm_jsonwrite(filename,params.info,struct('indent',' '));
                if isunix, try, system(sprintf('chmod 600 ''%s''',filename)); end; end
                fprintf('Connection information saved in %s\n',filename);
                if conn_server('isconnected'), 
                    ntries=0; 
                elseif ntries>0, 
                    startnewserver=true; params.info.local_port=[]; 
                    fprintf('No remote session found, starting a new session\n');
                else
                    fprintf('Unable to connect to remote CONN session\n',filename);
                end
            end
        end
        
     case 'setup' % saves .json info
        if numel(varargin)>=1&&~isempty(varargin{1}), filename=varargin{1};
        else filename=fullfile(conn_fileutils('homedir'),'connserverinfo.json');
        end
        if numel(varargin)>=2&&~isempty(varargin{2}), profilename=varargin{2};
        else profilename=conn_jobmanager('getdefault');
        end
        if numel(varargin)>=3&&~isempty(varargin{3}), ispersistent=varargin{3};
        else ispersistent=false;
        end
        if ispc, [nill,str1]=system('hostname');
        else [nill,str1]=system('hostname -f');
        end
        osquotes=char('"'*ispc+''''*~ispc);
        [nill,str2]=system('whoami');
        matlabpath=fullfile(matlabroot,'bin');
        whichfolders=cellfun(@(x)fileparts(which(x)),{'spm','conn'},'uni',0);
        isdep=false;
        try, isdep=isdeployed; end
        if conn_jobmanager('options','cmd_rundeployed'), isdep=true; end
        if isdep&&~isempty(conn_jobmanager('options','cmd_deployedfile')), fun_callback=conn_jobmanager('options','cmd_deployedfile');
        elseif isdep,                             fun_callback=conn_jobmanager('checkdeployedname');
        else                                      fun_callback=[osquotes fullfile(matlabpath,'matlab') osquotes];
        end
        if isdep,   cmd=sprintf('%s %s',fun_callback,'%s');
        else
            addpaths=cellfun(@(x)sprintf('addpath ''%s'';',x),whichfolders(cellfun('length',whichfolders)>0),'uni',0);
            cmd=sprintf('%s -nodesktop -noFigureWindows -nosplash -r "%s conn %s; exit"',...
                fun_callback, sprintf('%s ',addpaths{:}),'%s');
        end
        info = struct(...
            'host', regexprep(str1,'\n',''),...
            'CONNcmd',cmd,...
            'SERVERcmd',profilename,...
            'SERVERpersistent',ispersistent);
        if ~isempty(filename), spm_jsonwrite(filename,info); fprintf('Host information saved in %s\n',filename); end
        if nargout>0, varargout={info}; end
        if ~nargout&&isempty(filename), disp(info); end
        
    case {'exit','softexit'} % send exit signal to server to stop running (this will also cause the remote Matlab session to exit)
        if strcmpi(option,'softexit')||conn_server('isconnected'), 
            if strcmpi(option,'exit')
                fprintf('Exiting remote CONN session and disconnecting\n');
                conn_tcpip('write','exit');
            else
                fprintf('Disconnecting from remote CONN session\n');
            end
            if ~isempty(params.info.host), try, [ok,msg]=system(sprintf('%s -o ControlPath=''%s'' -O exit %s', params.options.cmd_ssh, params.info.filename_ctrl,params.info.login_ip)); end; end
            conn_tcpip('close');
            conn_cache clear;
            conn_jobmanager clear;
            %params.state='off';
            params=[];
        else
            fprintf('unable to connect to server, please terminate the server manually or use "conn_server_ssh restart" to restart the connection with the server and try "conn_server_ssh exit" again\n');
            if ~isempty(params.info.host), system(sprintf('%s -o ControlPath=''%s'' -O exit %s', params.options.cmd_ssh, params.info.filename_ctrl,params.info.login_ip)); end
        end
        
    case 'forceexit' % run remotely a command to forcibly delete the server's job
        if isfield(params.info,'host')&&~isempty(params.info.host)&&isfield(params.info,'login_ip')&&~isempty(params.info.login_ip)&&isfield(params.info,'remote_log')&&~isempty(params.info.remote_log)
            localcachefolder=conn_cache('private.local_folder');
            params.info.filename_ctrl=fullfile(localcachefolder,sprintf('connserver_ctrl_%s_%s',params.info.host,params.info.user));
            fprintf('Connecting to %s... ',params.info.login_ip);
            % starts a shared SSH connection
            [ok,msg]=system(sprintf('%s -o ControlPath=''%s'' -O check %s', params.options.cmd_ssh, params.info.filename_ctrl,params.info.login_ip)); 
            if ok~=0, 
                try, conn_fileutils('deletefile',params.info.filename_ctrl); end
                system(sprintf('%s -f -N -o ControlMaster=yes -o ControlPath=''%s'' %s', params.options.cmd_ssh, params.info.filename_ctrl,params.info.login_ip));
            end
            if ~isfield(params.info,'CONNcmd')||isempty(params.info.CONNcmd) % attempts to load server info from remote ~/connserverinfo.json file
                filename=fullfile(conn_cache('private.local_folder'),['conncache_', char(conn_tcpip('hash',mat2str(now)))]);
                conn_fileutils('deletefile',filename);
                fprintf('\nDownloading configuration information from %s:%s to %s\n',params.info.login_ip,'~/connserverinfo.json',filename);
                [ok,msg]=system(sprintf('%s -q -o ControlPath=''%s'' %s:~/connserverinfo.json %s',...
                    params.options.cmd_scp, params.info.filename_ctrl,params.info.login_ip,filename));
                assert(conn_existfile(filename),'unable to find ~/connserverinfo.json file in %s',params.info.login_ip);
                params.info.CONNcmd=conn_jsonread(filename,'CONNcmd');
            end            
            [ok,msg]=system(sprintf('%s -o ControlPath=''%s'' %s "%s"', params.options.cmd_ssh, params.info.filename_ctrl,params.info.login_ip, regexprep(sprintf(params.info.CONNcmd,sprintf('server_ssh submitexit ''%s''',params.info.remote_log)),'"','\\"')));
            if ok~=0, disp(msg); 
            else fprintf('Remote CONN session deletion requested successfully\n');
            end
        end
        
    case {'submitstart','submitstartpersistent'}
        if ~isempty(varargin)&&~isempty(varargin{1}), conn_jobmanager('setprofile',varargin{1}); end
        if ~isempty(varargin)&&~isempty(varargin{2}), id=char(varargin{2}); else id=[]; end
        info=conn_jobmanager('submit','orphan_conn',[],1,[],'server',regexprep(lower(option),'^submit',''),[],id); 
        fprintf('SSH_log in %s\n',info.pathname);
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
                fprintf('\nSSH_SUBMITSTART finished succesfully\n');
                match1=regexp(str,'\<ssh -L 6111:([^:]*):(\d+)','tokens');
                match2=regexp(str,'\<conn_server connect localhost \d+CONN(\w+)','tokens');
                fprintf('HOST:%s\nPORT:%s\nID:%s\nLOG:%s\n',match1{1}{1},match1{1}{2},match2{1}{1},info.pathname);
            else fprintf('SSH_SUBMITSTART error\n');
            end
        else fprintf('SSH_SUBMITSTART error\n');
        end
        
    case 'submitexit' % internal use only
        try
            pathname=varargin{1};
            info=struct; conn_loadmatfile(fullfile(pathname,'info.mat'),'info');
            info=conn_jobmanager('canceljob',info);
            fprintf('SSH_EXIT finished\n');
        catch
            fprintf('SSH_EXIT error\n');
        end
        
    case 'info'
        if numel(varargin)>=1, params.info=varargin{1}; end
        varargout={params.info};
        
    case 'options'
        if numel(varargin)>=1, 
            params.options=varargin{1}; 
            filename=fullfile(conn_fileutils('homedir'),'connclientinfo.json');
            spm_jsonwrite(filename,params.options,struct('indent',' '));
        end
        varargout={params.options};
        
    case 'details'
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

    case {'push','folderpush'}
%         clear h;
%         h.hfig=figure('units','norm','position',[.3 .4 .5 .2],'name','file transfer (scp)','numbertitle','off','menubar','none','color','w');
%         uicontrol('style','text','units','norm','position',[.1 .70 .19 .15],'string','copy FROM:','backgroundcolor','w','horizontalalignment','right','parent',h.hfig);
%         h.filesFROM=uicontrol('style','edit','max',1,'units','norm','position',[.3 .70 .6 .15],'string','','backgroundcolor','w','horizontalalignment','left','parent',h.hfig);
%         uicontrol('style','text','units','norm','position',[.1 .50 .19 .15],'string','TO:','backgroundcolor','w','horizontalalignment','right','parent',h.hfig);
%         h.filesTO=uicontrol('style','edit','max',1,'units','norm','position',[.3 .50 .6 .15],'string','','backgroundcolor','w','horizontalalignment','left','parent',h.hfig);
%         uicontrol('style','pushbutton','string','OK','units','norm','position',[.1,.01,.38,.25],'callback','uiresume');
%         uicontrol('style','pushbutton','string','Cancel','units','norm','position',[.51,.01,.38,.25],'callback','delete(gcbf)');
        filelocal=conn_server('util_localfile',varargin{1});
        fileremote=conn_server('util_localfile',varargin{2});
        if strcmpi(option,'folderpush'), ok=system(sprintf('%s -C -r -o ControlPath=''%s'' ''%s'' %s:''%s''', params.options.cmd_scp, params.info.filename_ctrl,regexprep(filelocal,'[\\\/]+$',''),params.info.login_ip,regexprep(fileremote,'[^\\\/]$','$0/')));
        else ok=system(sprintf('%s -C -q -o ControlPath=''%s'' ''%s'' %s:''%s''', params.options.cmd_scp, params.info.filename_ctrl,filelocal,params.info.login_ip,fileremote));
        end
        varargout={isequal(ok,0)}; 
        
    case {'pull','folderpull'}
        fileremote=conn_server('util_localfile',varargin{1});
        filelocal=conn_server('util_localfile',varargin{2});
        if strcmpi(option,'folderpull'), [ok,msg]=system(sprintf('%s -C -r -o ControlPath=''%s'' %s:''%s'' ''%s''', params.options.cmd_scp, params.info.filename_ctrl,params.info.login_ip,regexprep(fileremote,'[\\\/]+$',''),regexprep(filelocal,'[^\\\/]$',['$0','\',filesep])));
        else [ok,msg]=system(sprintf('%s -C -q -o ControlPath=''%s'' %s:''%s'' ''%s''', params.options.cmd_scp, params.info.filename_ctrl,params.info.login_ip,fileremote,filelocal));
        end
        if ~isequal(ok,0), disp(msg); end
        varargout={isequal(ok,0)};

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

function y = conn_server_ssh_input(x,varargin)
try
    fprintf('%s',char(x));
    y = input('', varargin{:});
catch
    if isempty(regexp(x,'\s+\[[^\]]+\]\:\s*$'))
        y = inputdlg(x);
    else
        y = inputdlg(regexprep(x,'\s+\[[^\]]+\]\:\s*$',' :'),'',1,regexp(x,'\s+\[([^\]]+)\]\:\s*$','tokens','once'));
    end
    if isempty(y), error('user-input canceled');
    elseif numel(varargin)>0&&isequal(varargin{1},'s'), y = y{1}; 
    else y = str2num(y{1});
    end
end
end

