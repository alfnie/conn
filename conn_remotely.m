function varargout=conn_remotely(option,varargin)
%
% Using CONN remotely
%
% CONN allows you to work and interact with CONN projects stored in any SSH-accessible computer
% (e.g. projects stored in your own server, or in your institution's cluster-computing environment)
% using the standard SSH protocol for tunneling all communications between your local computer and
% the computer where the CONN project data is stored and analyzed
%
% To use it simply select in the CONN gui the menu 'Project. Connect to remote projects' and log-in to
% your institution cluster-computing environment (or your personal server). Once connected CONN will be 
% able to work with projects and datasets stored in the server computer/network, and, if that server computer
% is part of a larger cluster-computing environment, CONN will also be able to use that larger environment to
% run preprocessing and analysis steps in parallel
%
%
% Advanced Configuration/Installation: %!
%
% To use CONN remotely you only need a Linux computer (referred to as the 'host') that 
%   a) it is located in a network-accessible location (e.g. fixed IP address) 
%   b) it has an OpenSSH-compatible SSH server/daemon running
%   c) it has Matlab/SPM/CONN installed
%
% You also need a second Linux/Mac/Windows computer (referred to as the 'client') that
%   a) it has regular network/internet access
%   b) it has an OpenSSH-compatible SSH client installed (e.g. this is already installed by default as 
%      part of the operating system in most Linux/Mac/Windows computers)
%
% Installation in the 'host' computer requires the standard installation of SPM/CONN, and also running
% just once as part of the installation procedure the Matlab command "conn_remotely install" (without quotes),
% which will create a file in your home directory [~/connserverinfo.json] describing the location of
% Matlab, CONN, and SPM packages (note: this command needs to be run for any user that wishes to access this
% computer remotely; or alternatively a fixed ~/connserverinfo.json file may be manually copied to the home
% directory of each user).
%
% Installation in the 'client' computer simply requires the standard installation of SPM/CONN.
%
% If the 'host' computer is part of a HPC or cluster environment, see
%    www.conn-toolbox.org/resources/cluster-configuration
% for additional HPC/cluster installation instructions. While optional, setting up CONN's cluster/HPC
% configuration in your 'host' computer allows CONN_REMOTELY to also submit jobs to your cluster remotely
% (jobs submitted from the 'client' computer which will be run in the local network of the 'host' computer).
%
% Advanced Usage: %!
%
% Manually start CONN from the 'client' computer using the Matlab command:
%
%    conn remotely
%
% When prompted enter:
%
%   "Server address" prompt  :      Enter the IP-address of the 'host' computer (this may be the target computer
%                                   you want to connect to, or simply a login-node within a target network you
%                                   want to connect to (e.g. when connecting to a computer cluster environment)
%   "Username" prompt        :      Enter your username in the 'host' computer (used to establish an SSH-connection
%                                   with the host)
%   "Password" prompt        :      Enter your password in the 'host' computer (note: CONN does not read nor store
%                                   your password, this will be read by your system's ssh program directly from tty)
%
% After this CONN will launch the standard CONN GUI and allow you to load and work with any CONN projects that may be
% accessible from the 'host' computer.
%
% Additional options: %!
%
%    conn remotely install     : (run in 'host' computer just once to configure this host to accept incoming requests)
%                                This step simply creates a file ~/connserverinfo.json with the following information:
%                                 CONNcmd     : OS-level command used to start/run CONN on 'host' computer
%                                 host        : name of 'host' computer (or login-node for HPC/cluster environments)
%                                 SERVERcmd   : name of HPC profile to use when starting a CONN 'server'
%                                             see "conn_jobmanager profiles" to display a list of available profiles
%                                             see "help conn_jobmanager" for additional options
%                                    (e.g. SERVERcmd='Background' will run CONN server process in the same computer
%                                    as the 'server' computer that we are connecting to, while SERVERcmd='Slurm' will
%                                    run CONN server process in a different computer which will be requested using
%                                    the options defined in the HPC.configuration 'Slurm' profile)
%                                 SERVERpersistent : 1/0 to restart server after connection with client ends
%                                When conn_remotely connects to a host it will try to first download this file from
%                                ~/connserverinfo.json for information on how to proceed
%
%    conn remotely start <IP>  : Same as 'conn remotely' but specfying the host IP address.
%                                This step starts CONN from the 'client' computer (this will attempt first to create a 
%                                new connection to the same CONN server used in our last connection to this host, and 
%                                if that fails it will ask the host computer to start a new CONN server and then it 
%                                will connect to it)
%
%    conn remotely startwithgui: Same as 'conn remotely' but with a graphical user interface
%
%    conn remotely restart     : re-starts CONN GUI after a dropped connection (this will not start a new CONN server
%                                but it will only attempt to create a new connection to the CONN server used in our
%                                last connection to this host)
%
%    conn isremotely           : re-starts CONN GUI (this will not start a new CONN server and it will not create a
%                                new connection)
%
%    conn remotely start local : Manually start conn client without SSH communications framework (see information
%                                below in "Alternative (non-SSH) structure)" section)
%
%    conn remotely start server : Manually start conn server without SSH communications framework (see information 
%                                below in "Alternative (non-SSH) structure)" section)
%
%    conn remotely settings    : Launches GUI for defining default client/server communication options
%
% Additional information: %!
%
% When conn_remotely starts, the 'client' computer will attempt to connect using ssh to the 'host' computer and it will
% establish a secure tunnel for communications. The host computer will then use its SERVERcmd profile in HPC/cluster
% configuration settings to launch a new Matlab session in a 'server' computer, which may be the same as the 'host'
% computer or a different computer within the same network (e.g. if using default HPC settings this session will be
% simply run as a new background Matlab process in the server computer; if using non-default HPC setting this session
% will be submitted to your cluster scheduler and run in a work-node as appropriate different from the login-node)
% which will in turn wait for the client's connection. Once the connection is established, the 'host' computer will
% only act as a gateway for communications between the 'client' and the 'server' computers. The 'server' computer
% will wait for instructions from the 'client', while the 'client' runs CONN normally (typically for anything that
% relates to GUI-interaction), querying the 'server' when needed (e.g. when needing to load data from remote files),
% or requesting the 'server' to run longer or more complex steps when appropriate (e.g. when preprocessing the data
% or running any analysis steps that requires significant interaction with the data). When the CONN GUI is closed in
% the 'client' computer, the 'server' will exit and the SSH tunnel will be closed.
%
%   Standard (SSH-based) structure of CONN REMOTELY components:
%                                 ----------------------------------------------------
%      --------                  |   --------                      ---------          |
%     | Client | INIT: -----------> | Host   | INIT: -----------> | Server  |         |
%     |________| <-- TCP/SSH -----> |________| <-- TCP/SSH -----> |_________|         |
%  (local CONN session)          |  (gateway)                (remote CONN session)    |
%                                |____________________________________________________|
%                                             (remote network)
%
% When working with a remote project, all of the elements in CONN's gui will be run locally by the 'client' computer
% and therefore will be as responsive as usual, but some steps in the GUI do require largers amount of data being
% transferred from the 'server' (e.g. when switching to a new tab, or when launching a 3D-view renderer) so these
% may take longer than usual to get started while that data is being transferred, please be patient. Also when
% working with a remote project and running preprocessing or analysis steps, you will be offered the option to have
% those steps run on the server Matlab session directly, or have the remote session use instead its own HPC/Cluster
% configuration options to submit a job to its HPC/Cluster environment (to clarify: there is no option to run those
% processes locally on your 'client' computer since that would require large amounts of data being transfered from
% the 'server' where the data is stored to the 'client' where the analysis would be run, which would make that option
% impractical in most scenarios).
%
% If the 'host' computer is NOT SSH-accessible (e.g. another computer in your home/office network), the
% following workaround allows you to still use CONN_REMOTELY functionality while setting up the server and the
% connection manually. Note that this will skip SSH entirely, so in this case your communication data will not be
% tunneled nor encrypted. This procedure is only recommended for communications WITHIN a local/secured network:
%
%   step 1) in the 'host/server' computer, run the Matlab command below to launch a CONN server process:
%
%               conn_remotely start server;
%
%           and then enter when prompted the desired one-time-use password and TCP-port number (leave empty to 
%           use the automatically search for an available port)
%
%   step 2) in the 'client' computer, close any active CONN windows and run the command:
%
%               conn_remotely start local;
%
%           and then enter when prompted the host/server IP address, port number, and password as defined above
%
%   Alternative (non-SSH) structure of CONN REMOTELY components:
%      -------------------------------------------------------
%     |      --------                      ---------------    |
%     |     | Client | <------TCP------>  | Host/Server   |   |
%     |     |________|                    |_______________|   |
%     |(local CONN session)             (remote CONN session) |
%     |_______________________________________________________|
%                         (local network)
%
%
% Advanced functions for debugging and command-line use: %!
%
% CONNECTION
%
% conn_remotely on [IP]       : starts remote CONN server and connect to it (equivalent to running "conn remotely" but
%                               without launching CONN GUI)
% conn_remotely off           : disconnects from remote CONN server (equivalent to closing GUI when using "conn remotely"
%                               syntax to start CONN); the server may remain listening for new connections or it may be 
%                               closed after disconnection depending on the server "conn_remotely settings" options
% conn_remotely offandon      : restarts communication with remote CONN server after dropped connection (equivalent
%                               to running "conn remotely restart" but without launching CONN GUI)
% conn_remotely forceoff      : disconnects from remote CONN server and closes server (note: this forces the server to close
%                               irrespective of the server "conn_remotely settings" options)
%
% REMOTE EXECUTION
%
% conn_remotely cmd                : interactive command-line execution in remote CONN server
% [...]=conn_remotely('cmd',...)   : single-command execution in remote CONN server
%                                    e.g. y=conn_remotely('cmd','fft',1:4);
%                                    e.g. conn_remotely('cmd','run','/usr/me/mybatch01.m')
%
% REMOTE FILE COPY
%
% conn_remotely('push',localfile,remotefolder)                : copies local file 'localfile' into remote folder 'remotefolder'
%                                                               e.g. conn_remotely push /data/file1.nii /usr/me/data/
% conn_remotely('folderpush',localfolder,remotefolder)        : copies the local folder 'localfolder' (and all of its contents)
%                                                               into remote folder 'remotefolder'
%                                                               e.g. conn_remotely push /data /usr/me
% conn_remotely('pull',remotefile,localfolder)                : copies remote file 'localfile' into local folder 'localfolder'
%                                                               e.g. conn_remotely pull /usr/me/data/file1.nii /data/
% conn_remotely('folderpull',remotefolder,localfolder)        : copies the remote folder 'remotefolder' (and all of its contents)
%                                                               into local folder 'localfolder'
%                                                               e.g. conn_remotely pull /usr/me/data /
% conn_remotely('push',localfile,remotefile)                  : copies local file 'localfile' to remote storage and rename it
%                                                               as 'remotefile'
%                                                               e.g. conn_remotely push /data/file1.nii /usr/me/data/file1.nii
% conn_remotely('pull',remotefile,localfile)                  : copies remote file 'remotefile' to local storage and rename it
%                                                               as 'localfile'
%                                                               e.g. conn_remotely pull /usr/me/data/file1.nii /data/file1.nii
%
% LOCAL EXECUTION
%
% When connected to a remote CONN server and running locally any conn_* function, most CONN functions that take
% files or directories as input arguments will accept filepaths of the form: '/CONNSERVER/[filepath]' to refer
% to a file/folder [filepath] in the remote server storage. For example:
%
% conn_loadmatfile /CONNSERVER/usr/me/file1.mat             : loads the variables in the remote file1.mat file into
%                                                             your (local) workspace
% conn_mesh_display /CONNSERVER/usr/me/data/file1.nii       : displays (locally) the remote file file1.nii
% conn_cache pull /CONNSERVER/usr/me/data/file1.nii         : creates a local-cache copy of the remote file file1.nii
% conn_display /CONNSERVER/usr/me/results/SPM.mat           : displays (locally) the remote 2nd-level results
% conn_fileutils mkdir /CONNSERVER/usr/me/newfolder         : creates a new remote directory (equivalent to: conn_remotely
%                                                             cmd mkdir /usr/me/newfolder)
% data=conn_surf_read('/CONNSERVER/usr/me/file1.nii');      : reads the 2D-surface nifti remote file file1.nii
% data=conn_mtx_read('/CONNSERVER/usr/me/file1.nii');       : reads the 2D-matrix nifti remote file file1.nii
% data=conn_fileutils('spm_read_vols','/CONNSERVER/usr/me/file1.nii'); : will read the 3D/4D nifti remote file file1.nii
%
% When connected to a remote CONN server and running locally CONN_HPC functions, CONN will always use the remote-server
% HPC options to submit jobs instead of the local HPC options. For example:
%
% conn('submit','disp','hello world')                       : will have the remote-server use its own HPC options
%                                                             (e.g. slurm profile) to submit a job to the remote network
%                                                             where the function disp('hello world') will be executed
% conn('submit','run','/usr/me/mybatch01.m')                : will have the remote-server use its own HPC options
%                                                             (e.g. slurm profile) to submit a job to the remote network where
%                                                             the script mybatch01 will be executed
% conn_jobmanager('getprofile')                             : lists the remote-server default HPC profile (just the
%                                                               same as str=conn_remotely('cmd','conn_jobmanager','getprofile'))
%

% notes: development reference notes on use of conn_server and conn_cache in CONN internal functions:
%
%
%   usage model #1: (a function will run remotely if it needs to work with remote files) (optimal when the fcn call requires minimal transfer of information)
%
%   function fileout = fcn(filein, varargin)
%      if any(conn_server('util_isremotefile',filein)), fileout=conn_server('util_remotefile',conn_server('run',mfilename,conn_server('util_localfile',filein),varargin{:})); return; end
%      % the code below this point works with normal/local files
%      ...
%
%
%   usage model #2: (a function will work on a local cache/copy of remote files) (optimal when working with mixture of local and remote files)
%
%   function fileout = fcn(filein, varargin)
%      isremotefile=conn_server('util_isremotefile',filein);
%      OPTION1 (when remote file already exists):   if isremotefile, remotefilename=filename; filename=conn_cache('pull',filename); end
%      OPTION2 (when remote file is to be created): if isremotefile, remotefilename=filename; filename=conn_cache('new',filename); end
%      ...
%      ... % the code here works with filename locally/normally
%      ...
%      if isremotefile, conn_cache('push',remotefilename); end
%   end
%
%
%   usage model #3: (a function will run remotely if it needs to work with remote variables) (optimal when working with large intermediate variables)
%
%   function output = fcn(input, varargin)
%      if conn_server('util_isremotevar',input), output=conn_server('run_keep',mfilename,input,varargin{:}); return; end
%      % the code below this point works with normal/local variable "input"
%      ...
%
%
%   usage model #4: (a function that accesses files only through calls to functions that follow model #1 or model #2)
%      [no changes needed]
%

global CONN_gui;
if ~nargin||isempty(option), option='start'; end

switch(option)
    case {'start','startwithgui','restart','end','softstart','softrestart','softend'}         % GUI interaction
        if isequal(varargin,{'server'}), conn_remotely('startserver',varargin{2:end}); return; end
        if isfield(CONN_gui,'isremote'), isremote=CONN_gui.isremote;
        else isremote=false; 
        end
        dosoft=~isempty(regexp(option,'^soft'));
        option=regexprep(option,'^soft','');
        dorestart=strcmpi(option,'restart');
        doend=strcmpi(option,'end');
        dogui=~isempty(regexp(option,'withgui$'));
        keepgui=conn('findgui');
        connversion={'CONN functional connectivity toolbox',' (',conn('ver'),') '};
        hfig=findobj('tag',connversion{1});
        if isempty(hfig), keepgui=false; end
        if keepgui
            if ~conn('gui_isready')||(~doend&&dosoft), Answ='Proceed';
            elseif doend, Answ=conn_questdlg({'Proceeding will close the current project and lose any unsaved progress','Do you want to proceed with starting CONN_locally?'},'CONN locally','Proceed','Cancel','Proceed');
            else Answ=conn_questdlg({'Proceeding will close the current project and lose any unsaved progress','Do you want to proceed with starting CONN_remotely?'},'CONN remotely','Proceed','Cancel','Proceed');
            end
        else Answ='Proceed';
        end
        if strcmp(Answ,'Proceed')
            if isequal(varargin,{'serverwithgui'}),
                conn initfromgui;
                conn importrois;
                conn gui_recent_init;
                conn gui_setup
                conn_remotely('startserverwithgui',varargin{2:end});
            elseif doend
                conn init;
                if dosoft, conn_server_ssh softexit; % note: softend does not stop server
                else conn_server_ssh exit;
                end
                CONN_gui.isremote=false;
                if keepgui
                    conn initfromgui;
                    conn importrois;
                    conn gui_recent_init;
                    conn gui_setup;
                else
                    conn;
                end
            else
                if dorestart, conn_server_ssh restart recent;
                elseif isequal(varargin,{'local'}), conn_server_ssh start local;
                elseif dogui, 
                    try, conn_server_ssh startwithgui recent; 
                    catch, return;
                    end
                else conn_server_ssh start recent; % note: change to "conn_server_ssh('start',varargin{:})" to allow user-defined .json files?
                end
                CONN_gui.isremote=true;
                conn gui_isconnected; % checks if connected
                if keepgui
                    if ~dosoft % note: softstart/softrestart does not close project
                        conn initfromgui;
                        conn importrois;
                    end
                    conn gui_recent_init;
                    conn gui_setup;
                else
                    conn isremotely;
                end
            end
        end
    case {'command','cmd'}
        if nargout>0, [varargout{1:nargout}]=conn_server('cmd',varargin{:});
        else conn_server('cmd',varargin{:});
        end
    case 'settings'
        tjson1=conn_server_ssh('options');
        filename2=fullfile(conn_projectmanager('homedir'),'connserverinfo.json');
        if conn_existfile(filename2), use_ssh=true; tjson2=conn_jsonread(filename2); else use_ssh=false; tjson2=struct; end
        if ~isfield(tjson2,'SERVERcmd'), tjson2.SERVERcmd=conn_jobmanager('getdefault'); end
        if ~isfield(tjson2,'SERVERpersistent'), tjson2.SERVERpersistent=false; end
        [tstr,tidx]=conn_jobmanager('profiles');
        profiles=tstr;
        tnull=find(strcmp('Null profile',tstr));
        tlocal=find(strcmp('Background process (Unix,Mac)',tstr),1);
        tvalid=setdiff(1:numel(tstr),tnull);
        tstr=cellfun(@(x)sprintf('distributed processing (run on %s)',x),tstr,'uni',0);
        if 1, tvalid=tidx; if isunix&&~isempty(tlocal)&&~ismember(tlocal,tvalid), tvalid=[tvalid(:)' tlocal]; end % show default+local profiles
        elseif 1, tvalid=tidx; % show only default profile
        else tstr{tidx}=sprintf('<HTML><b>%s</b></HTML>',tstr{tidx}); % show all profiles
        end
        toptions=tstr(tvalid);
        dprofile=strmatch(tjson2.SERVERcmd,profiles(tvalid),'exact');
        if isempty(dprofile), dprofile=strmatch(tjson2.SERVERcmd,profiles(tvalid)); end
        if isempty(dprofile), dprofile=1; else dprofile=dprofile(1); end
        
        handles.hfig=figure('units','norm','position',[.2 .35 .6 .35],'name','Remote connection settings','numbertitle','off','menubar','none','color','w','userdata',false);
        handles.save=uicontrol(handles.hfig,'style','pushbutton','units','norm','position',[.55,.01,.2,.10],'string','Save','callback','uiresume(gcbf)','tooltipstring','Save all profile changes for future Matlab sessions');
        handles.exit=uicontrol(handles.hfig,'style','pushbutton','units','norm','position',[.75,.01,.2,.10],'string','Exit','callback','close(gcbf)');
        
        strclient={{'With the currently-selected options, you may access CONN projects stored in remote computers directly over (unencrypted) TCP/IP (note: without SSH','security and encryption) by clicking the button labeled ''Manually connect to CONN server now''',' ','note: for this option to work, the remote computer needs to have the associated option ''Connect to this computer using SSH'' also unchecked, and',' a CONN server needs to be manually started there before you can connect to that server (see help conn_remotely for additional details)'} ...
                    {'With the currently-selected options, you may directly access CONN projects stored in remote computers','simply by selecting in the main CONN GUI the menu ''Project -> Connect to remote projects''',' ','note: for this option to work, the remote computer needs to be accessible using SSH and have in CONN the associated option','''Connect to this computer using SSH'' checked (see help conn_remotely for additional details)'}};
        uicontrol(handles.hfig,'style','text','units','norm','position',[.05 .90 .9 .075],'string','When this computer is client: (e.g. working from this computer on a CONN project stored remotely)','fontweight','normal','backgroundcolor','w','horizontalalignment','left');
        handles.client=uicontrol(handles.hfig,'style','checkbox','units','norm','position',[.15,.80,.75,.075],'string','Connect from this computer using SSH (secure/encrypted communications, servers are automatically started)','value',tjson1.use_ssh,'backgroundcolor','w','horizontalalignment','left','tooltipstring','<HTML>When this option is checked this computer uses secure/encrypted communications and servers are automatically started<br/>When this option is unchecked this computer uses unsecure/unencrypted communications and servers are manually started</HTML>');
        uicontrol(handles.hfig,'style','pushbutton','units','norm','position',[.10,.80,.025,.075],'string','?','backgroundcolor','w','callback',@(varargin)conn_msgbox(strclient{1+get(handles.client,'value')},'',-1),'tooltipstring','<HTML>help / how-to connect to remote servers</HTML>','userdata',handles.client);
        %handles.cmd_ssh_str=uicontrol(handles.hfig,'style','text','units','norm','position',[.20 .70 .45 .075],'string','Local command for logging into remote machine :','backgroundcolor','w','horizontalalignment','left');
        handles.cmd_ssh_str=uicontrol(handles.hfig,'style','popupmenu','units','norm','position',[.20,.70,.45,.075],'string',{'SSH using password authentication','SSH using public key authentication','manually defined connection method'},'backgroundcolor','w','fontname','monospaced','horizontalalignment','left','tooltipstring','<HTML>Method used to connect with server (SSH or SSH-compatible client application)<br/></HTML>','userdata',false);
        handles.cmd_ssh=uicontrol(handles.hfig,'style','edit','units','norm','position',[.65,.70,.25,.075],'string',tjson1.cmd_ssh,'backgroundcolor','w','fontname','monospaced','horizontalalignment','left','tooltipstring','<HTML>System command used to call SSH-client application (remote login)<br/>e.g. /usr/bin/ssh -i identity_file <br/> - note: this application command syntax must be compatible with OpenSSH SSH clients, including options for remote execution, control sockets for connection sharing, and TCP forwarding</HTML>');
        handles.cmd_scp_str=uicontrol(handles.hfig,'style','text','units','norm','position',[.20 .60 .45 .075],'string','Local command for remote file transfer :','backgroundcolor','w','horizontalalignment','left');
        handles.cmd_scp=uicontrol(handles.hfig,'style','edit','units','norm','position',[.65,.60,.25,.075],'string',tjson1.cmd_scp,'backgroundcolor','w','fontname','monospaced','horizontalalignment','left','tooltipstring','<HTML>System command used to call SCP application (secure copy)<br/>e.g. /usr/bin/scp -i identity_file <br/> - note: this application command syntax must be compatible with OpenSSH SCP</HTML>');
        handles.now_client=uicontrol(handles.hfig,'style','pushbutton','units','norm','position',[.20,.70,.7,.075],'string','Manually connect to CONN server now','callback','set(gcbf,''userdata'',true); uiresume(gcbf)','tooltipstring','<HTML>Connect to an existing CONN server over TCP/IP (note: without SSH secured/encrypted communications)</HTML>');
        hdls={[handles.cmd_ssh_str,handles.cmd_ssh,handles.cmd_scp_str,handles.cmd_scp], handles.now_client};
        if tjson1.use_ssh, set(hdls{1},'visible','on'); set(hdls{2},'visible','off'); else set(hdls{1},'visible','off'); set(hdls{2},'visible','on'); end
        set(handles.client,'callback','hdl=get(gcbo,''userdata''); if get(gcbo,''value''), set(hdl{1},''visible'',''on''); set(hdl{2},''visible'',''off''); else set(hdl{1},''visible'',''off''); set(hdl{2},''visible'',''on''); end','userdata',hdls);
        if isequal(tjson1.cmd_ssh,'ssh')&&isequal(tjson1.cmd_scp,'scp'), set([handles.cmd_ssh,handles.cmd_scp_str,handles.cmd_scp],'visible','off'); set(handles.cmd_ssh_str,'value',1);
        elseif ~isempty(regexp(tjson1.cmd_ssh,'^ssh -i'))&&~isempty(regexp(tjson1.cmd_scp,'^scp -i'))&&strcmp(regexprep(tjson1.cmd_ssh,'^ssh -i',''),regexprep(tjson1.cmd_scp,'^scp -i','')), set([handles.cmd_ssh,handles.cmd_scp_str,handles.cmd_scp],'visible','off'); set(handles.cmd_ssh_str,'value',2);
        else set([handles.cmd_ssh,handles.cmd_scp_str,handles.cmd_scp],'visible','on'); set(handles.cmd_ssh_str,'value',3);
        end
        set(handles.cmd_ssh_str,'callback',@(varargin)[set(handles.cmd_ssh_str,'userdata',true) set([handles.cmd_ssh,handles.cmd_scp_str,handles.cmd_scp],'visible',char(getfield({'off','off','on'},{get(gcbo,'value')})))]);
        
        if isfield(CONN_gui,'isremote')&&CONN_gui.isremote, % note: when connected to remote session options displayed are the options at the remote/server computer, not the local/client machine
            info=conn_remotely('info');
            if isfield(info,'host')&&~isempty(info.host), tnameserver=info.host;
            elseif isfield(info,'remote_ip')&&~isempty(info.remote_ip), tnameserver=info.remote_ip;
            else tnameserver='CONN server';
            end
            set(handles.hfig,'name',sprintf('Remote connection settings (at %s)',tnameserver));
        end
        if 0,%isfield(CONN_gui,'isremote')&&CONN_gui.isremote, 
            handles.server=[];
            handles.cmd_server=[];
            handles.now_server=[];
        else
            strserver={{'With the currently-selected options, you can only remotely access projects stored in this computer','when manually clicking below on the button labeled ''Manually start CONN server now''',' ','note: for this option to work, the client/other computer needs to have the associated option ''Connect from this computer using SSH'' unchecked (see help conn_remotely for additional details)'} ...
                {'With the currently-selected options, you can remotely access projects stored in this computer', '(as long as you have remote SSH access to this computer)',' ','note: for this option to work, this computer needs to be accessible using SSH and the client/other computer needs','to have in CONN the associated option ''Connect from this computer using SSH'' checked (see help conn_remotely for additional details)'}};
            uicontrol(handles.hfig,'style','text','units','norm','position',[.05 .40 .9 .075],'string','When this computer is server: (e.g. working from another computer on CONN projects stored on this machine)','fontweight','normal','backgroundcolor','w','horizontalalignment','left');
            handles.server=uicontrol(handles.hfig,'style','checkbox','units','norm','position',[.15,.30,.75,.05],'string','Connect to this computer using SSH (secure/encrypted communications, servers are automatically started)','value',use_ssh,'backgroundcolor','w','horizontalalignment','left','tooltipstring','<HTML>When this option is checked this computer uses secure/encrypted communications and servers are automatically started<br/>When this option is unchecked this computer uses unsecure/unencrypted communications and servers are manually started</HTML>');
            uicontrol(handles.hfig,'style','pushbutton','units','norm','position',[.10,.30,.025,.075],'string','?','backgroundcolor','w','callback',@(varargin)conn_msgbox(strserver{1+get(handles.server,'value')},'',-1),'tooltipstring','<HTML>help / how-to connect to remote servers</HTML>','userdata',handles.server);
            handles.cmd_server=uicontrol(handles.hfig,'style','popupmenu','units','norm','position',[.20,.20,.7,.075],'string',toptions,'value',dprofile,'backgroundcolor','w','horizontalalignment','left','tooltipstring','<HTML>HPC profile to use when starting a CONN server</HTML>');
            handles.now_server=uicontrol(handles.hfig,'style','pushbutton','units','norm','position',[.20,.20,.7,.075],'string','Manually start CONN server now','callback','close(gcbf);conn_remotely(''start'',''serverwithgui'');','tooltipstring','<HTML>Starts a CONN server on this computer and waits for client to connect over TCP/IP (note: without SSH secured/encrypted communications)</HTML>');
            handles.persistent_server=uicontrol(handles.hfig,'style','checkbox','units','norm','position',[.20,.125,.7,.075],'string','Close server when connection to client ends','value',~tjson2.SERVERpersistent,'backgroundcolor','w','horizontalalignment','left','tooltipstring','<HTML>Behavior of CONN server when client disconnects (when unchecked, CONN server will wait for a new client after a connection ends)</HTML>');
            hdls={handles.cmd_server handles.now_server handles.persistent_server};
            if use_ssh, set([hdls{1} hdls{3}],'visible','on'); set(hdls{2},'visible','off'); else set([hdls{1} hdls{3}],'visible','off'); set(hdls{2},'visible','on'); end
            set(handles.server,'callback','hdl=get(gcbo,''userdata''); if get(gcbo,''value''), set([hdl{1} hdl{3}],''visible'',''on''); set(hdl{2},''visible'',''off''); else set([hdl{1} hdl{3}],''visible'',''off''); set(hdl{2},''visible'',''on''); end','userdata',hdls);
        end
        while 1
            uiwait(handles.hfig);
            if ishandle(handles.hfig) % save
                tjson1.use_ssh=get(handles.client,'value');
                if tjson1.use_ssh
                    val=get(handles.cmd_ssh_str,'value');
                    switch(val)
                        case 1,
                            tjson1.cmd_ssh='ssh';
                            tjson1.cmd_scp='scp';
                            tjson1.use_key=false;
                        case 2,
                            tfilename1=char(regexp(tjson1.cmd_ssh,'ssh -i\s*(.*)','once','tokens'));
                            tfilename1=regexprep(tfilename1,'[''"]','');
                            if isempty(tfilename1)||get(handles.cmd_ssh_str,'userdata')>0
                                try; if isempty(tfilename1)&&exist('~/.ssh','dir'), tfilename1='~/.ssh/'; end; end
                                conn_msgbox({'Click ''Continue'' to select the identity file containing','your private key(s) for RSA or DSA authentication.','CONN will use the syntax ''ssh -i identity_file ...'' when connecting to your server',' ','e.g. id_rsa file created by ssh-keygen','e.g. key-pair-name.pem file from AWS'},'conn',2);
                                [tfilename1,tfilepath1]=conn_fileutils('uigetfile','*','Select identity file',tfilename1);
                                if ~ischar(tfilename1), continue; end
                                tfilename1=conn_server('util_localfile_filesep',[],fullfile(tfilepath1,tfilename1));
                            end
                            tjson1.use_key=true;
                            tjson1.file_key=tfilename1;
                            if ispc, 
                                tjson1.cmd_ssh=sprintf('ssh -i "%s"',tfilename1);
                                tjson1.cmd_scp=sprintf('scp -i "%s"',tfilename1);
                            else
                                tjson1.cmd_ssh=sprintf('ssh -i ''%s''',tfilename1);
                                tjson1.cmd_scp=sprintf('scp -i ''%s''',tfilename1);
                            end
                        case 3,
                            tjson1.cmd_ssh=get(handles.cmd_ssh,'string');
                            tjson1.cmd_scp=get(handles.cmd_scp,'string');
                            tjson1.use_key=false;
                    end
                end
                if isfield(CONN_gui,'isremote')&&CONN_gui.isremote, conn_server('run','conn_server_ssh','options',tjson1);
                else conn_server_ssh('options',tjson1); % client info
                end
                if ~isempty(handles.server)
                    if get(handles.server,'value')>0
                        dprofile=get(handles.cmd_server,'value');
                        tjson2.SERVERcmd=profiles{tvalid(dprofile)};
                        tjson2.SERVERpersistent=~get(handles.persistent_server,'value');
                        if isfield(CONN_gui,'isremote')&&CONN_gui.isremote, conn_server('run','conn_remotely','setup',conn_server('util_localfile',filename2),tjson2.SERVERcmd,tjson2.SERVERpersistent);
                        else conn_remotely('setup',filename2,tjson2.SERVERcmd,tjson2.SERVERpersistent); % server info
                        end
                    else
                        try, conn_fileutils('deletefile',filename2); end
                    end
                end
                if get(handles.hfig,'userdata'), % run CONN client now
                    conn_remotely;
                end
                delete(handles.hfig);
            end
            break;
        end
        
    case {'settings_runserver','startserver','startserverwithgui'}
        str=char(conn_tcpip('hash',mat2str(now))); str=str(1:8);
        isgui=usejava('desktop')&strcmp(lower(option),'startserverwithgui');
        port='';
        if isgui, 
            answ=conn_menu_inputdlg({'Create a one-time-use password to access this server:','Specify local TCP port: (leave empty for automatic selection)'},'',1,{str,port},struct('Resize','on'));
            if numel(answ)~=2||isempty(answ{1}),return; end
            str=answ{1};
            port=answ{2};
        else 
            str2=input(sprintf('Create a one-time-use password to access this server [%s]: ',str),'s');
            if ~isempty(str2), str=str2; end
            port=input('Specify local TCP port [auto]: ','s');
        end
        if ispc, [nill,str1]=system('hostname');
        else [nill,str1]=system('hostname -f');
        end
        hotsname=regexprep(str1,'\n','');
        if isempty(port)||isequal(port,'auto')
            tsocket=java.net.ServerSocket(0);
            portnumber=tsocket.getLocalPort;
            tsocket.close;
            clear tsocket
        else
            portnumber=str2double(port);
            assert(numel(portnumber)==1&isfinite(portnumber), 'unrecognized port number');
        end
        pause(1);
        strdisp={'Perform the following steps to connect to this server from a client/remote machine',' ','Step 1: start CONN using the syntax "conn remotely start local" (without quotes)',['Step 2: in prompt ''Remote session host address'' enter ',hotsname,' (or the IP address of this machine)'],['Step 3: in prompt ''Remote session port number'' enter ',num2str(portnumber)],['Step 4: in prompt ''Remote session id'' enter ',str]};
        if isgui
            [hmsg,hmsg2]=conn_msgbox([{'                                  CONN server running                                  '},repmat({' '},1,10)],'',-1, 2);
            set(hmsg2(1),'style','edit','max',2,'string',strdisp);
            set(hmsg2(2),'callback','if get(gcbo,''value''), conn_tcpip close; delete(gcbf); end');
            set(hmsg,'closerequestfcn','conn_tcpip close; delete(gcbf)');
            drawnow;
            conn_server('startpersistent',portnumber,conn_tcpip('keypair',str),hmsg2);
            if ishandle(hmsg), delete(hmsg); end
        else
            disp(char(strdisp));
            conn_server('startpersistent',portnumber,conn_tcpip('keypair',str),'silent');
        end
        
    case 'on'
        conn_server_ssh('start','recent',varargin{:});
        CONN_gui.isremote=true;
    case {'offandon','off&on','offon'}
        conn_server_ssh('restart','recent',varargin{:});
        CONN_gui.isremote=true;
    case {'off','softoff','forceoff'}
        conn_server_ssh(regexprep(option,'off$','exit'),varargin{:});
        CONN_gui.isremote=false;
    case {'push','pull','folderpush','folderpull','info','details','setup','install'}
        if nargout>0, [varargout{1:nargout}]=conn_server_ssh(option,varargin{:});
        else conn_server_ssh(option,varargin{:});
        end
    case 'speedtest'
        figure;
        ntest=60;
        n=1e3; 
        list=nan(ntest,2); 
        for itest=1:ntest, 
            fprintf('%d/%d\n',itest,ntest);
            t1=clock; 
            a=conn_server('run','randn',1,n); 
            t=etime(clock,t1); 
            list(itest,:)=[n*8*8/1024^2 t]; 
            if t>1, n=ceil(n/(1+rand)); 
            else n=ceil(n*(1+rand)); 
            end
            plot(list(:,1),list(:,2),'o'); xlabel('Mbits'); ylabel('seconds'); 
            if itest>10, 
                b=([1+0*list(1:itest,2) list(1:itest,2)]\list(1:itest,1));
                title(sprintf('%.2f Mbps',b(2))); 
            end
            drawnow;
        end
        
    otherwise
        error('unrecognized option %s',option);
end

end

