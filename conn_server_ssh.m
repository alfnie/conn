function varargout = conn_server_ssh(option, varargin)
% internal function

%
% alternative procedure for SSH-accessible networks only:
% command from server side: (from a computer within the SSH-accessible network)
%   conn_server_ssh('save' [, fileout])             : server configuration (run only once in order to configure this computer), creates a .json file with basic information necessary to connect to this server (machine IP, and where to find Matlab/SPM/CONN) [~/connserverinfo.json]
%
% commands from client side: (from a computer with a SSH client)
%   conn_server_ssh('start' [, IP])                 : SSH to network login node, submits a job to start a CONN server, and connects this computer to that remote server (note: expects to find IP:~/connserverinfo.json file; run "conn_server_ssh install" in server once)
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
        'options',struct('use_ssh',true,'cmd_ssh','ssh','cmd_scp','scp','use_key',false,'file_key',''),...
        'state','off');
    filename=fullfile(conn_fileutils('homedir'),'connclientinfo.json');
    if conn_existfile(filename), params.options=conn_jsonread(filename); end
    if ~isfield(params.options,'use_ssh'), params.options.use_ssh=true; end
    if ~isfield(params.options,'cmd_ssh'), params.options.cmd_ssh='ssh'; end
    if ~isfield(params.options,'cmd_scp'), params.options.cmd_scp='scp'; end
    if ~isfield(params.options,'use_key'), params.options.use_key=false; end
    if ~isfield(params.options,'file_key'), params.options.file_key=''; end
    params=conn_server_ssh_updatefilekey(params);
end

switch(lower(option))
    
    case {'start','startwithgui','restart'} % init server remotely and connect to it
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
        startwithgui=strcmpi(option,'startwithgui')&params.options.use_ssh;
        startwithgui_hmsg=[];
        if startwithgui, 
            allthesame=true;
            if ~isfield(params.info,'user')||isempty(params.info.user), [nill,str2]=system('whoami'); params.info.user=regexprep(str2,'\n',''); allthesame=false; end
            if ~isfield(params.info,'host')||isempty(params.info.host), params.info.host=''; allthesame=false; end
            clear h;
            h.hfig=figure('units','norm','position',[.3 .6 .3 .2],'name','Start SSH connection','numbertitle','off','menubar','none','color','w');
            uicontrol('style','text','units','norm','position',[0 .85 1 .15],'string','Connect using SSH to the network hosting your CONN projects','backgroundcolor','w','horizontalalignment','center','parent',h.hfig);
            uicontrol('style','text','units','norm','position',[.1 .62 .37 .14],'string','Remote server address:','backgroundcolor','w','horizontalalignment','right','parent',h.hfig);
            h.answer_host=uicontrol('style','edit','max',1,'units','norm','position',[.5 .62 .4 .15],'string',params.info.host,'backgroundcolor','w','horizontalalignment','left','parent',h.hfig,'tooltipstring','Enter the address of your SSH-accessible remote server or your HPC cluster login node');
            uicontrol('style','text','units','norm','position',[.1 .46 .37 .14],'string','Username:','backgroundcolor','w','horizontalalignment','right','parent',h.hfig);
            h.answer_user=uicontrol('style','edit','max',1,'units','norm','position',[.5 .46 .4 .15],'string',params.info.user,'backgroundcolor','w','horizontalalignment','left','parent',h.hfig,'tooltipstring','Enter your username (in the remote server or HPC cluster)');
            uicontrol('style','text','units','norm','position',[.1 .30 .37 .14],'string','Authentication method:','backgroundcolor','w','horizontalalignment','right','parent',h.hfig);
            h.answer_password=uicontrol('style','popupmenu','max',1,'units','norm','position',[.5 .30 .4 .15],'string',{'password','public key'},'backgroundcolor','w','horizontalalignment','left','parent',h.hfig,'tooltipstring',conn_menu_formathtml('<HTML>Select your remote server authentication method<br/>When selecting <i>password</i> you will be prompted to enter your password when the connection to the server is initiated (CONN does not store your password).<br/>When selecting <i>public key</i> you will be asked to select your SSH identity file (e.g. an id_rsa file containing a private key provided by your institution or one generated using ssh-keygen)</HTML>'),'callback','if get(gcbo,''value'')==2, conn_msgbox({''Click ''''Continue'''' to select the identity file containing'',''your private key(s) for RSA or DSA authentication.'',''CONN will use the syntax ''''ssh -i identity_file ...'''' when connecting to your server'','' '',''e.g. id_rsa file created by ssh-keygen'',''e.g. key-pair-name.pem file from AWS''},''conn'',2); [tfilename1,tfilepath1]=conn_fileutils(''uigetfile'',''*'',''Select identity file'',get(gcbo,''userdata'')); if ischar(tfilename1), tfilename1=conn_server(''util_localfile_filesep'',[],fullfile(tfilepath1,tfilename1)); set(gcbo,''userdata'',fullfile(tpathname,tfilename)); end; end', 'userdata',params.options.file_key,'value',1+(params.options.use_key&~isempty(params.options.file_key)));
            uicontrol('style','pushbutton','string','OK','units','norm','position',[.1,.01,.38,.15],'callback','uiresume','parent',h.hfig);
            uicontrol('style','pushbutton','string','Cancel','units','norm','position',[.51,.01,.38,.15],'callback','delete(gcbf)','parent',h.hfig);
            uiwait(h.hfig);
            if ~ishandle(h.hfig), error('Connection procedure canceled by user'); end
            params_info_host=get(h.answer_host,'string'); if ~isequal(params.info.host,params_info_host), allthesame=false; end; params.info.host=params_info_host;
            params_info_user=get(h.answer_user,'string'); if ~isequal(params.info.user,params_info_user), allthesame=false; end; params.info.user=params_info_user;
            params_options_use_key=get(h.answer_password,'value')==2; if ~isequal(params.options.use_key,params_options_use_key), allthesame=false; end; params.options.use_key=params_options_use_key;
            params_options_file_key=get(h.answer_password,'userdata'); if ~isequal(params.options.file_key,params_options_file_key), allthesame=false; end; params.options.file_key=params_options_file_key;
            params=conn_server_ssh_updatefilekey(params);
            delete(h.hfig);
            if params.options.use_key&&~isempty(params.options.file_key), startwithgui_hmsg=conn_msgbox('Connecting to remote server. Please wait','');
            elseif ispc, startwithgui_hmsg=conn_msgbox({'Connecting to remote server','Please enter your access credentials in a new OS command-line windows when prompted      '},'');
            else startwithgui_hmsg=conn_msgbox({'Connecting to remote server','Please enter your access credentials in Matlab''s command-line window when prompted      '},'');
            end
        end
        if params.options.use_ssh, 
            allthesame=true;
            if ~isfield(params.info,'host')||isempty(params.info.host), params.info.host=conn_server_ssh_input('Server address [local]: ','s'); allthesame=false;
            elseif startwithgui
                fprintf('Server address [%s]\n',params.info.host);
            else
                temp=conn_server_ssh_input(sprintf('Server address [%s]: ',params.info.host),'s');
                if ~isempty(temp),
                    temp=regexprep(temp,'\s+','');
                    %if ~isempty(temp,'ssh-?only$'), params.options.use_scp=false; temp=regexprep(temp,'ssh-?only$',''); end
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
            if startwithgui
                fprintf('Username [%s]\n',params.info.user);
            else
                temp=conn_server_ssh_input(sprintf('Username [%s]: ',params.info.user),'s');
                if ~isempty(temp),
                    if ~isequal(params.info.user,temp), allthesame=false; end
                    params.info.user=temp;
                end
            end
            params.info.login_ip=sprintf('%s@%s',params.info.user,params.info.host);
            localcachefolder=conn_cache('private.local_folder');
        end
        params.info.scp=false;
        params.info.windowscmbugfixed=false; % allows use of SSH ControlMaster option in Windows; set to true when Microsoft's OpenSSH for Windows fully supports -o ControlMaster
        if isempty(params.options.cmd_ssh), 
            error('No SSH client found (see Tools.RemoteOptions.Configuration)');
        else
            if ~isempty(params.info.host)
                params.info.filename_ctrl=fullfile(localcachefolder,conn_tcpip('hash',sprintf('connserver_ctrl_%s_%s',params.info.host,params.info.user)));
                try, conn_fileutils('deletefile',params.info.filename_ctrl); end
                % starts a shared SSH connection
                if ~ispc||~isfield(params.info,'windowscmbugfixed')||params.info.windowscmbugfixed
                    fprintf('Connecting to %s... ',params.info.login_ip);
                    system(sprintf('%s -f -N -o ControlMaster=yes -o ControlPath=''%s'' %s', params.options.cmd_ssh,params.info.filename_ctrl,params.info.login_ip),'-echo'); % starts a shared connection
                    [ok,msg]=system(sprintf('%s -o ControlPath=''%s'' -O check %s', params.options.cmd_ssh, params.info.filename_ctrl,params.info.login_ip)); 
                    if ok~=0, error(msg); end
%                     [ok,msg]=system(sprintf('%s -o ControlPath=''%s'' -O check %s', params.options.cmd_ssh, params.info.filename_ctrl,params.info.login_ip)); 
%                     fprintf('Connecting to %s... ',params.info.login_ip);
%                     if ok~=0,
%                         system(sprintf('%s -f -N -o ControlMaster=yes -o ControlPath=''%s'' %s', params.options.cmd_ssh,params.info.filename_ctrl,params.info.login_ip),'-echo'); % starts a shared connection
%                         [ok,msg]=system(sprintf('%s -o ControlPath=''%s'' -O check %s', params.options.cmd_ssh, params.info.filename_ctrl,params.info.login_ip));
%                         if ok~=0, error(msg); end
%                     end
                end
                if 1,%~isfield(params.info,'CONNcmd')||isempty(params.info.CONNcmd) % attempts to load server info from remote ~/connserverinfo.json file
                    filename=fullfile(conn_cache('private.local_folder'),['conncache_', char(conn_tcpip('hash',mat2str(now)))]);
                    conn_fileutils('deletefile',filename);
                    optionrepeat=true;
                    while optionrepeat
                        optionrepeat=false;
                        if ~ispc||~isfield(params.info,'windowscmbugfixed')||params.info.windowscmbugfixed
                            fprintf('\nReading configuration information from %s:%s to %s\n',params.info.login_ip,'~/connserverinfo.json',filename);
                            if 0,%params.options.use_scp % note: use ssh in case scp not enabled at server
                                [ok,msg]=system(sprintf('%s -q -o ControlPath=''%s'' %s:~/connserverinfo.json %s',...
                                    params.options.cmd_scp, params.info.filename_ctrl,params.info.login_ip,filename));
                            else
                                [ok,msg]=system(sprintf('%s -q -o ControlPath=''%s'' %s ''cat ~/connserverinfo.json'' > %s',...
                                    params.options.cmd_ssh, params.info.filename_ctrl,params.info.login_ip,filename));
                            end
                        else
                            tstr=sprintf('SSH downloading configuration information from %s:%s',params.info.login_ip,'~/connserverinfo.json');
                            if 0,%params.options.use_scp
                                [ok,msg]=system(sprintf('start "%s" /WAIT cmd /c "%s -q %s:~/connserverinfo.json %s"',...
                                    tstr, params.options.cmd_scp, params.info.login_ip,filename));
                            else
                                [ok,msg]=system(sprintf('start "%s" /WAIT cmd /c %s %s "cat ~/connserverinfo.json" ^> %s 2^>^&1', ...
                                    tstr, params.options.cmd_ssh, params.info.login_ip, filename));
                            end
                        end
                        tjson=[]; 
                        if conn_existfile(filename), 
                            tjson=conn_jsonread(filename); 
                            if ~isstruct(tjson)||~isfield(tjson,'CONNcmd') % if non-empty but still non-interpretable json file, removes non-json header (anything preceeding the "{" symbol in file matching the last "}" symbol in file) 
                                tjsonstr=fliplr(conn_fileutils('fileread',filename));
                                tjsonlvl1=cumsum(tjsonstr=='{');
                                tjsonlvl2=cumsum(tjsonstr=='}');
                                conn_fileutils('filewrite_raw',filename,fliplr(tjsonstr(1:find((tjsonlvl1>0)&(tjsonlvl1==tjsonlvl2),1,'first'))));
                                tjson=conn_jsonread(filename);
                            end
                        end
                        if ~conn_existfile(filename)||isempty(tjson)||~isstruct(tjson)||~isfield(tjson,'CONNcmd')
                            fprintf('Unable to find CONN distribution in %s.\nIf this is your first time connecting to %s, please use the following steps to confirm that CONN is available there and then try connecting again:\n   1. Log in to %s as user %s\n   2. Launch Matlab\n   3. From Matlab command-window type "conn remotely setup" (without quotes) and confirm that the file ~/connserverinfo.json is created\n   4. Log out from %s\n(see "in the server computer" section at https://web.conn-toolbox.org/resources/remote-configuration for details)\n\n',params.info.login_ip,params.info.host,params.info.host,params.info.user,params.info.host);
                            optionrepeat=isequal(input(sprintf('Try connecting again to %s? (yes|no) : ',params.info.login_ip),'s'),'yes');
                            assert(optionrepeat,'Unable to proceed. Missing ~/connserverinfo.json file in %s',params.info.login_ip);
                        end
                    end
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
                    CONNcmd=regexprep(params.info.CONNcmd,'-singleCompThread','-singleCompThread -nojvm'); % note: adds -nojvm to reduce memory load of Matlab startup in login server
                    if ~ispc||~isfield(params.info,'windowscmbugfixed')||params.info.windowscmbugfixed
                        [ok,msg]=system(sprintf('%s -o ControlPath=''%s'' %s "%s"', params.options.cmd_ssh, params.info.filename_ctrl,params.info.login_ip, regexprep(sprintf(CONNcmd,tstr),'"','\\"')));
                    else
                        tfilename=fullfile(conn_cache('private.local_folder'),['conncache_', char(conn_tcpip('hash',mat2str(now)))]);
                        [ok,msg]=system(sprintf('start "SSH requesting a new Matlab session in %s. This may take a few minutes. Please wait" /WAIT cmd /c %s %s "%s" ^> %s 2^>^&1', params.info.login_ip, params.options.cmd_ssh, params.info.login_ip, regexprep(sprintf(CONNcmd,tstr),'"','\\"'),tfilename));
                        msg=conn_fileutils('fileread',tfilename);
                        conn_fileutils('deletefile',tfilename);
                    end
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
                    else conn_server('connect',params.info.remote_ip,sprintf('%d:%s',params.info.remote_port,params.info.remote_id));
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
                    
                    if ~ispc||~isfield(params.info,'windowscmbugfixed')||params.info.windowscmbugfixed
                        fprintf('Establishing secure communication path to remote session (%d:%s:%d)\n',params.info.local_port,params.info.remote_ip,params.info.remote_port);
                        [ok,msg]=system(sprintf('%s -o ControlPath=''%s'' -O forward -L%d:%s:%d %s', params.options.cmd_ssh, params.info.filename_ctrl,params.info.local_port,params.info.remote_ip,params.info.remote_port,params.info.login_ip));
                        try, disp(char(msg)); end
                    else
                        tstr=sprintf('SSH secure communication channel %d:%s:%d',params.info.local_port,params.info.remote_ip,params.info.remote_port);
                        [ok,msg]=system(sprintf('start "%s" cmd /c "%s -f -N -L%d:%s:%d %s"', tstr, params.options.cmd_ssh, params.info.local_port,params.info.remote_ip,params.info.remote_port,params.info.login_ip));
                        try, disp(char(msg)); end
                    end
                    if ok~=0, 
                        params.info.local_port=[]; 
                    else
                        %system(sprintf('ssh -f -N -T -o ExitOnForwardFailure=yes -o ControlPath=''%s'' -L%d:%s:%d %s', params.info.filename_ctrl,params.info.local_port,params.info.remote_ip,params.info.remote_port,params.info.login_ip));
                        %fprintf('Connecting to server\n');
                        if strcmpi(option,'restart'), conn_tcpip('open','client','localhost',params.info.local_port,params.info.remote_id,0); params.state='on';
                        else conn_server('connect','localhost',sprintf('%d:%s',params.info.local_port,params.info.remote_id));
                        end
                    end
                end
                ntries=max(0,ntries-1);
                if conn_server('isconnected'), 
                    ntries=0; 
                    filename=fullfile(conn_fileutils('homedir'),'conn_recentservers.json');
                    spm_jsonwrite(filename,params.info,struct('indent',' '));
                    if isunix, try, system(sprintf('chmod 600 ''%s''',filename)); end; end
                    fprintf('Connection information saved in %s\n',filename);
                elseif ntries>0, 
                    startnewserver=true; params.info.local_port=[]; 
                    fprintf('No remote session found, starting a new session\n');
                else
                    fprintf('Unable to connect to remote CONN session\n');
                end
                if ~isempty(params.info.host)&&ispc&&isfield(params.info,'windowscmbugfixed')&&~params.info.windowscmbugfixed, try, for n1=1:2, pause(1); [ok,msg]=system(sprintf('taskkill /FI "WindowTitle eq %s" /F',tstr)); end; end; end
            end
        end
        if ~isempty(startwithgui_hmsg), delete(startwithgui_hmsg); end

     case {'setup','install'} % saves .json info
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
            cmd=sprintf('%s -singleCompThread -nodesktop -noFigureWindows -nosplash -r "%s conn %s; exit"',...
                fun_callback, sprintf('%s ',addpaths{:}),'%s');
        end

        %if ispersistent, OScmd=sprintf('%s -singleCompThread -nodesktop -noFigureWindows -nosplash -r "%s conn server start '''' ''%s''; exit"', fun_callback, sprintf('%s ',addpaths{:}), '%s');
        %else             OScmd=sprintf('%s -singleCompThread -nodesktop -noFigureWindows -nosplash -r "%s conn server startpersistent '''' ''%s''; exit"', fun_callback, sprintf('%s ',addpaths{:}), '%s');
        %end
        %info=conn_jobmanager('submit','orphan_conn',[],1,[],'server','$START$',[],'$ID$'); 
        % job=conn_jobmanager('job','process_orphan','conn',[],'server',cmdstart,[],'$ID$');
        % info=conn_jobmanager('createjob',job,{1});
        % info=conn_jobmanager('submitjob',info,[],false);
        % conn_savematfile(fullfile(info.pathname,'info.mat'),'info');

        info = struct(...
            'host', regexprep(str1,'\n',''),...
            'CONNcmd',cmd,...
            'SERVERcmd',profilename,...
            'SERVERpersistent',ispersistent);
        if ~isempty(filename), spm_jsonwrite(filename,info,struct('indent',' ')); fprintf('CONN information saved in %s\n',filename); end
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
            try, if ~isempty(params.info.host)&& (~ispc||~isfield(params.info,'windowscmbugfixed')||params.info.windowscmbugfixed), [ok,msg]=system(sprintf('%s -o ControlPath=''%s'' -O exit %s', params.options.cmd_ssh, params.info.filename_ctrl,params.info.login_ip)); end; end
            conn_tcpip('close');
            conn_cache clear;
            conn_jobmanager clear;
            %params.state='off';
            params=[];
        else
            fprintf('unable to connect to server, please terminate the server manually or use "conn_server_ssh restart" to restart the connection with the server and try "conn_server_ssh exit" again\n');
            try, if ~isempty(params.info.host)&& (~ispc||~isfield(params.info,'windowscmbugfixed')||params.info.windowscmbugfixed), system(sprintf('%s -o ControlPath=''%s'' -O exit %s', params.options.cmd_ssh, params.info.filename_ctrl,params.info.login_ip)); end; end
        end
        
    case 'forceexit' % run remotely a command to forcibly delete the server's job
        if isfield(params.info,'host')&&~isempty(params.info.host)&&isfield(params.info,'login_ip')&&~isempty(params.info.login_ip)&&isfield(params.info,'remote_log')&&~isempty(params.info.remote_log)
            localcachefolder=conn_cache('private.local_folder');
            params.info.filename_ctrl=fullfile(localcachefolder,conn_tcpip('hash',sprintf('connserver_ctrl_%s_%s',params.info.host,params.info.user)));
            if ~ispc||~isfield(params.info,'windowscmbugfixed')||params.info.windowscmbugfixed
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
                    fprintf('\nReading configuration information from %s:%s to %s\n',params.info.login_ip,'~/connserverinfo.json',filename);
                    %[ok,msg]=system(sprintf('%s -q -o ControlPath=''%s'' %s:~/connserverinfo.json %s',...
                    %    params.options.cmd_scp, params.info.filename_ctrl,params.info.login_ip,filename));
                    [ok,msg]=system(sprintf('%s -q -o ControlPath=''%s'' %s ''cat ~/connserverinfo.json'' > %s',...
                        params.options.cmd_ssh, params.info.filename_ctrl,params.info.login_ip,filename));
                    assert(conn_existfile(filename),'unable to find ~/connserverinfo.json file in %s',params.info.login_ip);
                    params.info.CONNcmd=conn_jsonread(filename,'CONNcmd');
                end
                [ok,msg]=system(sprintf('%s -o ControlPath=''%s'' %s "%s"', params.options.cmd_ssh, params.info.filename_ctrl,params.info.login_ip, regexprep(sprintf(params.info.CONNcmd,sprintf('server_ssh submitexit ''%s''',params.info.remote_log)),'"','\\"')));
            else
                [ok,msg]=system(sprintf('start "SSH exiting remote CONN server" /WAIT cmd /c %s %s "%s"', params.options.cmd_ssh, params.info.login_ip, regexprep(sprintf(params.info.CONNcmd,sprintf('server_ssh submitexit ''%s''',params.info.remote_log)),'"','\\"')));
            end
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
            if ok, % reads from conn_tcpip printout the target HOST&PORT
                fprintf('\nSSH_SUBMITSTART finished succesfully\n');
                match1=regexp(str,'\<ssh -L 6111:([^:]*):(\d+)','tokens');
                match2=regexp(str,'\<conn_server connect localhost \d+:(\w+)','tokens');
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

    case 'checkcontrolmaster'
        if ~(isempty(params.info.filename_ctrl)||exist(params.info.filename_ctrl,'file')), % if ControlMaster is missing or timed-out re-authenticate to create a new SSH ControlMaster socket
            try, conn_fileutils('deletefile',params.info.filename_ctrl); end
            fprintf('re-authenticating with %s... ',params.info.login_ip);
            system(sprintf('%s -f -N -o ControlMaster=yes -o ControlPath=''%s'' %s', params.options.cmd_ssh,params.info.filename_ctrl,params.info.login_ip),'-echo'); % starts a shared connection
            [ok,msg]=system(sprintf('%s -o ControlPath=''%s'' -O check %s', params.options.cmd_ssh, params.info.filename_ctrl,params.info.login_ip));
            if ok~=0, error(msg); end
            varargout={true};
        else
            varargout={false};
        end

    case {'push','folderpush'}
%         clear h;
%         h.hfig=figure('units','norm','position',[.3 .4 .5 .2],'name','file transfer (scp)','numbertitle','off','menubar','none','color','w');
%         uicontrol('style','text','units','norm','position',[.1 .70 .19 .15],'string','copy FROM:','backgroundcolor','w','horizontalalignment','right','parent',h.hfig);
%         h.filesFROM=uicontrol('style','edit','max',1,'units','norm','position',[.3 .70 .6 .15],'string','','backgroundcolor','w','horizontalalignment','left','parent',h.hfig);
%         uicontrol('style','text','units','norm','position',[.1 .50 .19 .15],'string','TO:','backgroundcolor','w','horizontalalignment','right','parent',h.hfig);
%         h.filesTO=uicontrol('style','edit','max',1,'units','norm','position',[.3 .50 .6 .15],'string','','backgroundcolor','w','horizontalalignment','left','parent',h.hfig);
%         uicontrol('style','pushbutton','string','OK','units','norm','position',[.1,.01,.38,.25],'callback','uiresume');
%         uicontrol('style','pushbutton','string','Cancel','units','norm','position',[.51,.01,.38,.25],'callback','delete(gcbf)');
        filelocal=conn_server('util_localfile_filesep','/',regexprep(varargin{1},'^(\w*):',''));
        fileremote=conn_server('util_localfile_filesep','/',regexprep(varargin{2},'^(\w*):',''));
        if ~ispc||~isfield(params.info,'windowscmbugfixed')||params.info.windowscmbugfixed
            conn_server_ssh('checkcontrolmaster');
            if strcmpi(option,'folderpush'), ok=system(sprintf('%s -C -r -o ControlPath=''%s'' ''%s'' %s:''%s''', params.options.cmd_scp, params.info.filename_ctrl,regexprep(filelocal,'[\\\/]+$',''),params.info.login_ip,regexprep(fileremote,'[^\\\/]$','$0/')));
            else ok=system(sprintf('%s -C -q -o ControlPath=''%s'' ''%s'' %s:''%s''', params.options.cmd_scp, params.info.filename_ctrl,filelocal,params.info.login_ip,fileremote));
            end
        else 
            tstr=sprintf('SSH copying to %s',params.info.login_ip);
            if strcmpi(option,'folderpush'), ok=system(sprintf('start "%s" /WAIT cmd /c %s -C -r "%s" %s:"%s"', tstr, params.options.cmd_scp, regexprep(filelocal,'[\\\/]+$',''),params.info.login_ip,regexprep(fileremote,'[^\\\/]$','$0/')));
            else ok=system(sprintf('start "%s" /WAIT cmd /c %s -C -q "%s" %s:"%s"', tstr, params.options.cmd_scp, filelocal,params.info.login_ip,fileremote));
            end
        end
        varargout={isequal(ok,0)}; 
        
    case {'pull','folderpull'}
        fileremote=conn_server('util_localfile_filesep','/',regexprep(varargin{1},'^(\w*):',''));
        filelocal=conn_server('util_localfile_filesep','/',regexprep(varargin{2},'^(\w*):',''));
        if ~ispc||~isfield(params.info,'windowscmbugfixed')||params.info.windowscmbugfixed
            conn_server_ssh('checkcontrolmaster');
            if strcmpi(option,'folderpull'), [ok,msg]=system(sprintf('%s -C -r -o ControlPath=''%s'' %s:''%s'' ''%s''', params.options.cmd_scp, params.info.filename_ctrl,params.info.login_ip,regexprep(fileremote,'[\\\/]+$',''),regexprep(filelocal,'[^\\\/]$',['$0','\',filesep])));
            else [ok,msg]=system(sprintf('%s -C -q -o ControlPath=''%s'' %s:''%s'' ''%s''', params.options.cmd_scp, params.info.filename_ctrl,params.info.login_ip,fileremote,filelocal));
            end
        else
            tstr=sprintf('SSH copying to %s',params.info.login_ip);
            if strcmpi(option,'folderpull'), [ok,msg]=system(sprintf('start "%s" /WAIT cmd /c %s -C -r %s:"%s" "%s"', tstr, params.options.cmd_scp, params.info.login_ip,regexprep(fileremote,'[\\\/]+$',''),regexprep(filelocal,'[^\\\/]$',['$0','\',filesep])));
            else [ok,msg]=system(sprintf('start "%s" /WAIT cmd /c %s -C -q %s:"%s" "%s"', tstr, params.options.cmd_scp, params.info.login_ip,fileremote,filelocal));
            end
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

function params=conn_server_ssh_updatefilekey(params)
if ~isempty(regexp(params.options.cmd_ssh, 'ssh\s\-i\s')),
    params.options.use_key=true;
    params.options.file_key=regexprep(regexprep(params.options.cmd_ssh,'.*ssh\s\-i\s+',''),'[\''\"]','');
    params.options.cmd_ssh=regexprep(params.options.cmd_ssh,'^(.*ssh)\s\-i\s.*$','$1');
end
if ~isempty(regexp(params.options.cmd_scp, 'scp\s\-i\s')),
    params.options.cmd_scp=regexprep(params.options.cmd_scp,'^(.*scp)\s\-i\s.*$','$1');
end
if params.options.use_ssh&&params.options.use_key&&~isempty(params.options.file_key)
    if ispc
        params.options.cmd_ssh=sprintf('%s -i "%s"', params.options.cmd_ssh, params.options.file_key);
        params.options.cmd_scp=sprintf('%s -i "%s"', params.options.cmd_scp, params.options.file_key);
    else
        params.options.cmd_ssh=sprintf('%s -i ''%s''', params.options.cmd_ssh, params.options.file_key);
        params.options.cmd_scp=sprintf('%s -i ''%s''', params.options.cmd_scp, params.options.file_key);
    end
end
end
