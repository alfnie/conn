function varargout=conn_remotely(varargin)
% 
% Using CONN remotely
% 
% CONN_REMOTELY allows you to work and interact with CONN projects stored in any SSH-accessible 
% computer (e.g. projects stored in your own server, or in your institution's cluster-computing environment)
% 
% Configuration/Installation: %!
%
% To use this functionality you need to have access to a Linux computer ("server") that is in a network-
% accessible location and which has a ssh server/daemon running, as well as access to a second Linux/Mac 
% computer ("client") with regular network/internet access and a ssh client installed (e.g. this is installed 
% by default in Linux/Mac). 
%
% Installation in the "server" computer requires the standard installation of SPM/CONN, and then simply
% running (just once as part of the installation procedure) the Matlab command: 
%
%    conn_remotely setup
% 
% which will create a file in your home directory [~/connserverinfo.json] describing the location of
% Matlab, CONN, and SPM installations. This command needs to be run for any user that wishes to access
% this server remotely (alternatively the ~/connserverinfo.json file may be copied across different users).
% If the "server" computer is part of a HPC or cluster environment, CONN_REMOTELY will also be able to initiate
% jobs in your cluster remotely (see https://web.conn-toolbox.org/resources/cluster-configuration for HPC/cluster 
% installation instructions). In this case, the "server" computer may simply be any login-node of your cluster, 
% and CONN will use this login-node only as a proxy while all CPU-intensive processes will be run in automatically-
% spawned work-nodes as required. 
%
% Installation in the "client" computer simply requires the standard installation of SPM/CONN. 
%
% Usage: %!
%
% Start CONN from the "client" computer using the Matlab command:
%
%    conn_remotely
% 
% When prompted enter:
%
%   "Server address" prompt  :      Enter the IP-address of the "server" computer (this may be a login-node in a 
%                                   computing cluster environment) 
%   "Username" prompt        :      Enter your username in the "server" computer (used to establish an SSH-connection with 
%                                   the "server")
%   "Password" prompt        :      Enter your password in the "server" computer (note: CONN does not read or record your 
%                                   password, this is read by your system's ssh command directly from tty)
%
% After this CONN will launch the standard CONN gui and allow you to load and work with remote CONN projects
% stored in your "server" computer. 
%
%
% note: if the "server" computer is NOT SSH-accessible (e.g. another computer in your home/office network), the following 
% workaround allows you to still use CONN_REMOTELY while setting up the server and the connection manually:
%
%   step 1) in the "server" computer, run the Matlab command "conn server" to manually launch a CONN server process
%           From the output printed by this command, note the first occurrence of a keyphrase of the form 
%           conn_server <IP_address> <numeric_code>CONN<alphanumeric_code> (to be used below)
%
%   step 2) in the "client" computer, run "conn_remotely" and then enter:
%           2.a) in the "Server address" prompt, enter 'none' (without quotes)
%           2.b) in the "Remote session host address" prompt, enter the <IP_address> portion of the keyphrase above
%           2.c) in the "Remote session access port" prompt, enter the <numeric_code> portion of the keyphrase above
%           2.d) in the "Remote session id" prompt, enter the <alphanumeric_code> portion of the keyphrase above
%           2.e) in the "Remote session log folder" prompt, leave empty (press Return)
% 


if nargin>=1&&isequal(varargin{1},'setup'),conn_server('HPC_save'); return; end

keepgui=false;
if nargin>=1, keepgui=true; end
if nargin>=1&&isequal(varargin{1},'restart'),dorestart=true;
else dorestart=false;
end
connversion={'CONN functional connectivity toolbox',' (',conn('ver'),') '};
hfig=findobj('tag',connversion{1});
if keepgui||isempty(hfig), Answ='Proceed';
else Answ=conn_questdlg({'Proceeding will close the current project and loose any unsaved progress','Do you want to proceed with starting CONN remotely?'},'CONN remotely','Proceed','Cancel','Proceed');
end
if strcmp(Answ,'Proceed')
    if ~keepgui&&~isempty(hfig), conn forceclose; end
    doinitconnection=true;
    if 0, %isequal(conn_server('state'),'on')
        fprintf('Testing current connection...\n');
        if conn_server('isconnected'), doinitconnection=false;
        else fprintf('Starting new connection\n');
        end
    end
    if doinitconnection
        try
            if dorestart, conn_server HPC_restart;
            else conn_server HPC_start;
            end
        catch me
            fprintf('ERROR: Unable to start remote session: %s',me.message);
        end
    end
    if keepgui conn gui_setup
    else conn guiisremote;
    end
end


end