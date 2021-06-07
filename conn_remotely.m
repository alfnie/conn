function varargout=conn_remotely(option,varargin)
% 
% Using CONN remotely
% 
% CONN_REMOTELY allows you to work and interact with CONN projects stored in any SSH-accessible 
% computer (e.g. projects stored in your own server, or in your institution's cluster-computing environment)
% using the standard SSH protocol for tunneling all communications between your local computer and 
% the computer where the CONN project data is stored and analyzed
% 
% Configuration/Installation: %!
%
% To use this functionality you need to have access to a Linux computer ('host') that is in a network-
% accessible location (e.g. fixed IP address) and that has a OpenSSH-compatible SSH server/daemon running, 
% as well as access to a second Linux/Mac/Windows computer ('client') with regular network/internet access 
% and a OpenSSH-compatible SSH client installed (e.g. this is already installed by default as part of the
% operating system in most Linux & Mac computers). 
%
% Installation in the 'host' computer requires the standard installation of SPM/CONN, and also running
% just once as part of the installation procedure the Matlab command "conn_remotely setup" (without quotes),
% which will create a file in your home directory [~/connserverinfo.json] describing the location of
% Matlab, CONN, and SPM packages. This command needs to be run for any user that wishes to access this
% computer remotely (alternatively another ~/connserverinfo.json file may be manually copied to a user's
% home directory).
%
% If the 'host' computer is part of a HPC or cluster environment, see 
%    www.conn-toolbox.org/resources/cluster-configuration 
% for additional HPC/cluster installation instructions. While optional, setting up CONN's cluster/HPC 
% configuration in your 'host' computer allows CONN_REMOTELY to also submit jobs to your cluster remotely 
% (jobs submitted from the 'client' computer which will be run in the local network of the 'host' computer). 
%
% Installation in the 'client' computer simply requires the standard installation of SPM/CONN. 
%
% Usage: %!
%
% Start CONN from the 'client' computer using the Matlab command:
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
%                                   your password, this will be read by your system's ssh command directly from tty)
%
% After this CONN will launch the standard CONN GUI and allow you to load and work with any CONN projects that may be 
% accessible from the 'host' computer. 
%
% Options: %!
%
%    conn remotely setup       : (run in 'host' computer just once) creates file ~/connserverinfo.json with the 
%                                following information:
%                                 CONNcmd     : OS-level command used to start/run CONN on 'host' computer
%                                 host        : name of 'host' computer (or login-node for HPC/cluster environments)
%                                 SERVERcmd   : name of HPC profile to use when starting a CONN 'server' 
%                                             see "conn_jobmanager profiles" to display a list of available profiles
%                                             see "help conn_jobmanager" for additional options
%                                    (e.g. SERVERcmd='Background' will run CONN server process in the same computer
%                                    as the 'server' computer that we are connecting to, while SERVERcmd='Slurm' will
%                                    run CONN server process in a different computer which will be requested using
%                                    the options defined in the HPC.configuration 'Slurm' profile)
%
%    conn remotely             : starts CONN from the 'client' computer (this will start a new CONN server and connect
%                                to it)
%
%    conn remotely restart     : re-starts CONN GUI after a dropped connection (this will not start a new CONN server 
%                                but instead will create a new connection to and existing server)
%
%    conn isremotely           : re-starts CONN GUI (this will not start a new CONN server and it will not create a 
%                                new connection)
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
%                                          ------------------------------------------------------------------------
%      --------                           |     --------                                       ---------           |     
%     | Client | INIT: ------(ssh)-----------> | Host   | INIT: ---(bash/sbatch/qsub/etc.)--> | Server  |          |
%     |________| COMM: <--TCP / SSH tunnel---> |________| COMM: <----TCP / SSH tunnel-------> |_________|          |
% (local CONN session)                    |    (gateway)                                  (remote CONN session)    |
%                                         |________________________________________________________________________|
%                                                   (remote network)
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
% tunneled nor encrypted:
%
%   step 1) in the 'client' computer, run the following Matlab command and write-down/store the contents of the 
%           resulting public and private keys:
%               [PUBLICKEY, PRIVATEKEY] = conn_tcpip('keypair')
%
%   step 2) in the 'host/server' computer, run the Matlab command below to launch a CONN server process:
%               conn_server('start',[],'enter-your-public-key-here');
%           From the output printed by this command, note the first occurrence of a keyphrase of the form: 
%               conn_server <IP_address> <numeric_code>CONNprivatekey
%
%   step 3) in the 'client' computer, run "conn_remotely" and then enter:
%           2.a) in the "Server address" prompt:                 enter 'none' (without quotes)
%           2.b) in the "Remote session host address" prompt:    enter the <IP_address> portion of the keyphrase above
%           2.c) in the "Remote session access port" prompt:     enter the <numeric_code> portion of the keyphrase above
%           2.d) in the "Remote session id" prompt:              enter your private key here
%           2.e) in the "Remote session log folder" prompt:      leave empty (press Return)
% 
%   Alternative (non-SSH) structure of CONN REMOTELY components:
%      -------------------------------------------------------
%     |      --------                      ---------------    |
%     |     | Client | COMM: <---TCP--->  | Host/Server   |   |
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
% conn_remotely on            : starts remote CONN server and connect to it (equivalent to running "conn remotely" but 
%                               without launching CONN GUI)
% conn_remotely off           : disconnects from remote CONN server (equivalent to closing GUI when using "conn remotely" 
%                               syntax to start CONN)
% conn_remotely offandon      : restarts communication with remote CONN server after dropped connection (equivalent 
%                               to running "conn remotely restart" but without launching CONN GUI)
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

global CONN_gui;
if ~nargin||isempty(option), option='init'; end

switch(option)
    case {'init','start','restart'}         % GUI interaction
        dorestart=strcmpi(option,'restart');
        dostart=strcmpi(option,'start');
        if dostart||dorestart, keepgui=conn('findgui');
        else keepgui=false;
        end
        connversion={'CONN functional connectivity toolbox',' (',conn('ver'),') '};
        hfig=findobj('tag',connversion{1});
        if keepgui||isempty(hfig), Answ='Proceed';
        else Answ=conn_questdlg({'Proceeding will close the current project and loose any unsaved progress','Do you want to proceed with starting CONN remotely?'},'CONN remotely','Proceed','Cancel','Proceed');
        end
        if strcmp(Answ,'Proceed')
            if ~keepgui&&~isempty(hfig), conn forceclose; end
            doinitconnection=true;
            if doinitconnection
                try
                    if dorestart, conn_server SSH_restart recent;
                    else conn_server SSH_start % note: change to "conn_server('SSH_start',varargin{:})" to allow user-defined .json files
                    end
                catch me
                    fprintf('ERROR: Unable to start remote session: %s',me.message);
                end
            end
            CONN_gui.isremote=true;
            if keepgui 
                conn gui_setup
            else
                conn gui_isconnected;
                conn isremotely;
            end
        end
    case 'setup'
        conn_server SSH_save;
    case 'on'
        conn_server('SSH_start',varargin{:});
        CONN_gui.isremote=true;
    case {'offandon','off&on','offon'}
        conn_server SSH_restart;
        CONN_gui.isremote=true;
    case 'off'
        conn_server SSH_exit;
        CONN_gui.isremote=false;
    case {'push','pull','folderpush','folderpull'}
        conn_server(['SSH_',option],varargin{:});
    case {'command','cmd'}
        singlecommand=nargout>0|~isempty(varargin);
        info=conn_server('SSH_info');
        if isfield(info,'host')&&~isempty(info.host), tnameserver=info.host;
        else tnameserver='none';
        end
        while 1
            if ~isempty(varargin), cmd=varargin{1}; opts=varargin(2:end);
            else cmd=input([tnameserver,' >> '],'s'); opts={};
            end
            switch(lower(cmd))
                case {'quit','exit'}, break;
                otherwise
                    try
                        if singlecommand,
                            [varargout{1:nargout}]=conn_server('run','conn','cmd',cmd,opts{:});
                        else
                            str=conn_server('run','conn','cmd_capture',cmd);
                            disp(str);
                        end
                    catch me
                        str=regexprep(char(getReport(me,'extended','hyperlinks','off')),'[\t ]+',' ');
                        disp(str);
                    end
            end
            if singlecommand, break; end
        end
    otherwise
        error('unrecognized option %s',option);
end

end