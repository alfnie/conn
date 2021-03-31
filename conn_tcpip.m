
function varargout = conn_tcpip(option, varargin)
% CONN_TCPIP manages TCP/IP client/server communication
%
% commands:
%   conn_tcpip('open','server' [,PORT, ID])         : open TCP/IP port in local machine and waits for client to connect (default PORT=6111 ID='')
%                                                    If no PORT is entered, conn_tcpip will first attempt port 6111 and if that fails it will
%                                                    bind instead to the first available port
%
%   conn_tcpip('open','client' [,IP, PORT, ID])     : connect to server TCP/IP port at IP:PORT (default IP='127.0.0.1' PORT=6111 ID='')
%                                                    If an ID string was provided when opening the server, the server will reject any 
%                                                    client connection that is not open using the exact same ID (case sensitive)
%
%   conn_tcpip('write',VAR)                         : sends data (from server/client to client/server)
%                                                    VAR may be any Matlab variable (numeric/string/cell/struct/etc.)
%
%   VAR = conn_tcpip('read')                        : receives data (from server/client to client/server)
%
%   conn_tcpip('writefromfile',FILENAME)            : sends raw data from file (continuous stream / no memory-limitation)
%
%   TIMESTAMP = conn_tcpip('readtofile',FILENAME)    : receives raw data and save it to file (continuous stream / no memory-limitation)
%                                                    Returns TIMESTAMP (conn_tcpip sethash timestamp) or MD5-hash (conn_tcpip sethash md5) of received data 
%
%   conn_tcpip('writekeepalive')                    : sends an empty/keepalive packet
%
%   conn_tcpip('close')                             : stops server/client
%   conn_tcpip('clear')                             : clears connection buffers
%
%   conn_tcpip('sethash','md5')                     : uses MD5-hash to identify changes in files (slower)
%   conn_tcpip('sethash','timestamp')               : uses filesystem last-modified timestamps to identify changes in files (faster)
%
% e.g.
%
%   ON MACHINE #1 (server)
%
%     >> conn tcpip open server
%     Opening port 6111...
%     Waiting for client connection...
%     *************************************************************************************************************
% 
%     To connect to this server from a different machine, use the Matlab syntax:
%        conn_tcpip open client alfonsosmbp2018.lan
% 
%     To connect to this server using ssh-tunneling, use the following syntax instead:
%        $(OS-command)     :   ssh -L 6111:localhost:6111 alfnie@alfonsosmbp2018.lan
%        >(Matlab-command) :   conn_tcpip open client 127.0.0.1 6111
% 
%     *************************************************************************************************************
%
%   ON MACHINE #2 (client)
%
%     >> conn_tcpip open client alfonsosmbp2018.lan
%     Connecting to alfonsosmbp2018.lan:6111...
%     Succesfully established connection to server
%
%   ON MACHINE #1 or #2
%   
%     >> conn_tcpip('write', rand(1,10));
%
%   ON MACHINE #2 or #1
%   
%     >> disp(conn_tcpip('read'));
%     0.3902    0.8448    0.3720    0.7076    0.3457    0.8772    0.7505    0.4344    0.4007    0.3729
%

persistent connection base64

varargout={};
if isempty(connection)
    connection=struct(...
        'isserver',false,...
        'ip','127.0.0.1',...
        'port',6111,...
        'id',[],...
        'hash','timestamp',... % timestamp|md5|sha-1|sha-256|sha-384|sha-512
        'buffer',[],...
        'length',nan,...
        'header1',[],...
        'header2',[],...
        'maxlength',65532); % note: want divisible by 12 and <=65535
    try, base64=org.apache.commons.codec.binary.Base64; end
    if ispc, connection.cache=conn_fullfile(getenv('USERPROFILE'),'.conn_cache');
    else connection.cache=conn_fullfile('~/.conn_cache');
    end
    try, [nill,nill]=mkdir(connection.cache); end
    try, if ~conn_existfile(connection.cache,true), connection.cache=pwd; end; end
end

switch(lower(option))
    case 'open' % open client/server socket
        type=varargin{1};        
        connection.isserver=strcmpi(type,'server');
        connection.buffer=[];
        connection.length=nan;
        connection.header1=[];
        connection.header2=[];
        forceport=false;
        if connection.isserver          % open server PORT
            if numel(varargin)>=2 && ~isempty(varargin{2}), connection.port=varargin{2}; forceport=true; end
            if numel(varargin)>=3 && ~isempty(varargin{3}), connection.id=varargin{3}; else connection.id=''; end
            if numel(varargin)>=4 && ~isempty(varargin{4}), disphelp=varargin{4}; else disphelp=true; end
        else                            % open client IP PORT ID
            if numel(varargin)>=2 && ~isempty(varargin{2}), connection.ip=varargin{2}; else connection.ip='127.0.0.1'; end
            if numel(varargin)>=3 && ~isempty(varargin{3}), connection.port=varargin{3}; else connection.port=6111; end
            if numel(varargin)>=4 && ~isempty(varargin{4}), connection.id=varargin{4}; else connection.id=''; end
        end
        if ischar(connection.port), connection.port=str2double(connection.port); end
        try, connection.socket.close; pause(rand); end
        try, connection.channel.close; if connection.isserver, pause(rand); else pause(5+rand); end; end
        if connection.isserver
            try
                if isempty(connection.port)||isequal(connection.port,0)
                    connection.socket=java.net.ServerSocket(0);
                    connection.port=connection.socket.getLocalPort;
                    fprintf('Opening port %d...\n',connection.port);
                else
                    fprintf('Opening port %d...\n',connection.port);
                    connection.socket=java.net.ServerSocket(connection.port);
                end
            catch
                if forceport
                    error('Port %d unavailable. Try a different port, or leave the port entry empty to have conn_tcpip bind to the first available port',connection.port);
                else
                    connection.socket=java.net.ServerSocket(0);
                    connection.port=connection.socket.getLocalPort;
                    fprintf('Port %d unavailable. Switching to port %d...\n',6111,connection.port);
                end
            end
            connection.socket.setSoTimeout(1000*10);
            %if ~isempty(connection.id), fprintf('Server ID = %s\n',connection.id); end
            if disphelp
                if ispc, [nill,str1]=system('hostname');
                else [nill,str1]=system('hostname -f');
                end
                [nill,str2]=system('whoami');
                if disphelp>1 % help for conn_server use
                    fullid=[num2str(connection.port), 'CONN', connection.id];
                    fprintf('******************************************************************************\n\n');
                    fprintf('To connect to this server, use the Matlab syntax:\n');
                    fprintf('   conn_server connect %s %s\n\n',regexprep(str1,'\n',''),fullid);
                    fprintf('To connect to this server using ssh-tunneling, use the following syntax instead:\n');
                    fprintf('   $(OS-command)     :   ssh -L 6111:%s:%d %s@%s\n',regexprep(str1,'\n',''),connection.port,regexprep(str2,'\n',''),'YOUR-INSTITUTION-SSH-LOGIN-NODE');
                    fprintf('   >(Matlab-command) :   conn_server connect localhost %s\n\n',regexprep(fullid,'^\d+','6111'));
                    fprintf('******************************************************************************\n');
                else % help for conn_tcpip use
                    fprintf('******************************************************************************\n\n');
                    fprintf('To connect to this server, use the Matlab syntax:\n');
                    if ~isempty(connection.id), fprintf('   conn_tcpip open client %s %d ''%s''\n\n',regexprep(str1,'\n',''),connection.port,connection.id);
                    elseif connection.port==6111, fprintf('   conn_tcpip open client %s\n\n',regexprep(str1,'\n',''));
                    else fprintf('   conn_tcpip open client %s %d\n\n',regexprep(str1,'\n',''),connection.port);
                    end
                    fprintf('******************************************************************************\n');
                end
            end
            fprintf('Waiting for client connection...');
            ok=false;
            while ~ok
                try
                    connection.channel=connection.socket.accept;
                    ok=true;
                    connection.channel.setSoTimeout(1000*60);
                catch
                    fprintf('.');
                end
            end
            fprintf('\nClient connected...\n');
        else
            connection.socket=[];
            fprintf('Connecting to %s:%d...\n',connection.ip,connection.port);
            connection.channel=java.net.Socket();
            connection.channel.setSoTimeout(1000*60);
            addr=java.net.InetSocketAddress(connection.ip,connection.port);
            connection.channel.connect(addr);
            %connection.channel=java.net.Socket(connection.ip,connection.port);
        end
        ostream   = connection.channel.getOutputStream;        
        connection.output.stream=java.io.DataOutputStream(ostream);
        istream   = connection.channel.getInputStream;  
        connection.input.stream = java.io.DataInputStream(istream);
        
        if isempty(connection.id), handshake='-'; 
        else handshake=connection.id;
        end
        try,
            if connection.isserver % handshake
                varcheck=conn_tcpip('read');
                if isequal(varcheck,'handshake'), varcheck=conn_tcpip('read'); end
                if isequal(varcheck,handshake)
                    conn_tcpip('write','ok');
                    if true||isempty(connection.id), fprintf('Succesfully established connection to client\n');
                    else fprintf('Succesfully established connection to client. Connection ID = %s\n',connection.id);
                    end
                else
                    conn_tcpip('close');
                    fprintf('Client was unable to match CONN server ID. Try starting a server at a different port\n');
                end
            else
                for ntries=1:3
                    try
                        pause(2^(ntries-1)+rand);
                        conn_tcpip('write','handshake');
                        conn_tcpip('write',handshake);
                        ok=conn_tcpip('read');
                        if isequal(ok,'ok')
                            if true||isempty(connection.id), fprintf('Succesfully established connection to server\n');
                            else fprintf('Succesfully established connection to server. Connection ID = %s\n',connection.id);
                            end
                        else
                            conn_tcpip('close');
                            fprintf('Client was unable to match CONN server ID\n');
                        end
                        break;
                    end
                end
            end
        catch
            conn_tcpip('close');
            fprintf('Communication failure. Try again at a later time\n');
        end
                
    case 'clear' % clear buffers
        connection.buffer=[];
        connection.length=nan;
        connection.header1=[];
        connection.header2=[];
        try
            while connection.input.stream.available>0
                data=char(connection.input.stream.readUTF);
            end
        end
        
    case 'close' % close socket
        try, connection.channel.close; end
        try, connection.socket.close; end
        connection.buffer=[];
        connection.length=nan;
        connection.header1=[];
        connection.header2=[];
        
    case {'read','readtofile'} % reads Matlab variable from remote
        hash=[];
        readtofile=strcmpi(option,'readtofile');
        if readtofile
            filename=varargin{1}; % reads raw data to file (data does not need to fit in memory)
            if ~isempty(connection.hash)&&nargout>0&&~isequal(connection.hash,'timestamp'), hash=java.security.MessageDigest.getInstance(connection.hash); end
        elseif 1 % read to file then load var
            filename=fullfile(connection.cache,['conntcpipread_',char(mlreportgen.utils.hash(mat2str(now))),'.mat']);
        else % read to var directly
            filename=[];
        end
        filehandle=[];
        bored=false;
        
        while 1
            data=[];
            bytes_available = connection.input.stream.available;
            if bytes_available>0 || bored
                %fprintf('.');
                data=char(connection.input.stream.readUTF);
            end
            bored=true;
            if ~isempty(data)
                connection.header1=[connection.header1,numel(connection.buffer)+reshape(find(data=='<'),1,[])];
                connection.header2=[connection.header2,numel(connection.buffer)+reshape(find(data=='>'),1,[])];
                connection.buffer=[connection.buffer, data];
                bored=false;
            end
            %disp(connection.buffer);
            if ~isempty(connection.header1)&&connection.header1(1)>1&&connection.header1(1)<=connection.length % found a new header, yet incomplete data
                i2=connection.header1(1)-1;
                connection.buffer=connection.buffer(i2+1:end);
                connection.header1=connection.header1-i2; connection.header1(connection.header1<=0)=[];
                connection.header2=connection.header2-i2; connection.header2(connection.header2<=0)=[];
                connection.length=nan;
                fprintf('Incomplete TCP packet. Disregarding\n');
                break;
            end
            if ~isempty(filename)&&~isnan(connection.length), % writing to file
                while connection.length>0 && numel(connection.buffer)>=min(connection.maxlength,connection.length)
                    rlength=min(connection.maxlength,connection.length);
                    if isempty(base64), tdata=matlab.net.base64decode(connection.buffer(1:rlength)); %(slow)
                    else tdata=typecast(base64.decode(uint8(connection.buffer(1:rlength))')','uint8');
                    end
                    if isempty(filehandle), 
                        [tpath,tname,text]=fileparts(filename);
                        cachefilename=fullfile(tpath,['cachetmp_',tname,text]);
                        filehandle=fopen(cachefilename,'wb'); 
                    end
                    fwrite(filehandle,tdata,'uint8');
                    if ~isempty(hash), hash.update(tdata); end
                    i2=rlength;
                    connection.buffer=connection.buffer(i2+1:end);
                    connection.header1=connection.header1-i2; connection.header1(connection.header1<=0)=[];
                    connection.header2=connection.header2-i2; connection.header2(connection.header2<=0)=[];
                    connection.length=connection.length-rlength;
                end
                if connection.length<=0,
                    varargout={[]};
                    if connection.length==0 % finished 
                        if isempty(filehandle),
                            [tpath,tname,text]=fileparts(filename);
                            cachefilename=fullfile(tpath,['cachetmp_',tname,text]);
                            filehandle=fopen(cachefilename,'wb');
                        end
                        fclose(filehandle);
                        if readtofile
                            if ispc, [ok,nill]=system(['move "',cachefilename,'" "',filename,'"']);
                            else, [ok,nill]=system(['mv -f ''',cachefilename,''' ''',filename,'''']);
                            end
                            if ~isequal(ok,0), error('Error moving file from %s to %s, check target permissions',cachefilename,filename); end
                            if ~isempty(hash), varargout={typecast(hash.digest,'uint8')};
                            elseif isequal(connection.hash,'timestamp'), fd=dir(filename); varargout={fd.datenum};
                            end
                        else
                            loadok=false;
                            try
                                load(cachefilename,'msg');
                                loadok=true;
                            end
                            if ispc, [ok,nill]=system(sprintf('del "%s"',cachefilename));
                            else [ok,nill]=system(sprintf('rm -f ''%s''',cachefilename));
                            end
                            if loadok, varargout={msg};
                            else fprintf('Error reading variable (possibly incomplete data). Disregarding');
                            end
                        end
                    end
                    connection.length=nan;
                    return;
                end
                bored=false;
            elseif isempty(filename)&&numel(connection.buffer)>=connection.length % finished reading var
                if connection.length>0
                    if isempty(base64), tdata=matlab.net.base64decode(connection.buffer(1:connection.length)); %(slow)
                    else tdata=typecast(base64.decode(uint8(connection.buffer(1:connection.length))')','uint8');
                    end
                    try
                        msg=getArrayFromByteStream(tdata);
                        varargout={msg};
                        connection.buffer=connection.buffer(connection.length+1:end);
                        connection.header1=connection.header1-connection.length; connection.header1(connection.header1<=0)=[];
                        connection.header2=connection.header2-connection.length; connection.header2(connection.header2<=0)=[];
                        connection.length=nan;
                        return
                    end
                    fprintf('Unable to decode TCP packet (possibly incomplete data). Disregarding\n');
                end
                if ~isempty(connection.header1)
                    i2=connection.header1(1)-1;
                    connection.buffer=connection.buffer(i2+1:end);
                    connection.header1=connection.header1-i2; connection.header1(connection.header1<=0)=[];
                    connection.header2=connection.header2-i2; connection.header2(connection.header2<=0)=[];
                    connection.length=nan;
                else
                    connection.buffer=[]; % note: any data in buffer between endofmessage and next header is disregarded
                    connection.header1=[];
                    connection.header2=[];
                    connection.length=nan;
                end
                bored=false;
            end
            while isnan(connection.length)&&~isempty(connection.header1)&&~isempty(connection.header2)&&any(connection.header2>connection.header1(1)) % read header
                i1=connection.header1(1);
                i2=connection.header2(find(connection.header2>i1,1));
                i1=connection.header1(find(connection.header1<i2,1,'last'));
                connection.length=str2double(connection.buffer(i1+1:i2-1)); % note: empty or non-numeric are keepalive packets (length=NaN)
                connection.buffer=connection.buffer(i2+1:end);
                connection.header1=connection.header1-i2; connection.header1(connection.header1<=0)=[];
                connection.header2=connection.header2-i2; connection.header2(connection.header2<=0)=[];
                bored=false;
            end
            %if bored, pause(rand*bored); end
        end
        if ~isempty(filename), 
            if ~isempty(filehandle), 
                fclose(filehandle); 
                if ispc, [ok,nill]=system(sprintf('del "%s"',cachefilename));
                else [ok,nill]=system(sprintf('rm -f ''%s''',cachefilename));
                end
            end
            varargout={[]};
        else varargout={};
        end
        
    case 'writekeepalive'
        try, connection.output.stream.writeUTF('<>'); end
            
    case 'write' % writes Matlab variable to remote
        if numel(varargin)>=1, msg=varargin{1};
        else msg='hello world';
        end
        if 1 % save to file then send file
            filename=fullfile(connection.cache,['conntcpipwrite_',char(mlreportgen.utils.hash(mat2str(now))),'.mat']);
            info=whos('msg');
            if info.bytes>2e9, save(filename,'msg','-v7.3'); 
            else save(filename,'msg','-v7'); 
            end
            ok=conn_tcpip('writefromfile',filename);
            if ispc, [nill,nill]=system(sprintf('del "%s"',filename));
            else [nill,nill]=system(sprintf('rm -f ''%s''',filename));
            end
        else % send data directly (unused failback)
            if isempty(msg), data='';
            elseif isempty(base64), data=matlab.net.base64encode(mygetByteStreamFromArray(msg)); %(slow)
            else, data=char(base64.encode(uint8(mygetByteStreamFromArray(msg))')');
            end
            header=sprintf('<%d>',numel(data));
            ok=false;
            try
                connection.output.stream.writeUTF(header);
                i1=0;i2=numel(data);
                while i1<i2
                    tdata=data(i1+1:min(i2,i1+connection.maxlength));
                    connection.output.stream.writeUTF(tdata);
                    i1=i1+connection.maxlength;
                end
                ok=true;
            catch
                fprintf('Unable to send TCP packet (possibly unresponsive server). Disregarding\n');
            end
        end
        varargout={ok};
        
    case 'writefromfile' % writes raw data from file to remote (data does not need to fit in memory)
        filename=varargin{1}; 
        assert(isempty(regexp(filename,'<>')),'invalid characters in filename');
        fh=fopen(filename,'rb');
        if isequal(fh,-1)
            %fprintf('unable to find file %s. Disregarding\n',filename);
            fsize=-1;
        else
            fseek(fh,0,'eof');
            fsize=4*ceil(ftell(fh)/3); % size in base64
            fseek(fh,0,'bof');
        end
        ok=false;
        data=char([]);
        header=sprintf('<%d>',fsize);
        try
            connection.output.stream.writeUTF(header);
        catch
            fprintf('Unable to send TCP packet (possibly unresponsive server). Disregarding\n');
        end
        while fsize>0&&~ok
            if numel(data)<connection.maxlength, 
                msg=fread(fh, connection.maxlength, 'uint8')'; 
                if isempty(msg), 
                elseif isempty(base64), data=[data matlab.net.base64encode(msg)]; %(slow)
                else, data=[data char(base64.encode(msg')')];
                end
            end        
            if isempty(data), ok=true; break; end
            idx=1:min(numel(data),connection.maxlength);
            try
                connection.output.stream.writeUTF(data(idx));
                data(idx)=[];
            catch
                fprintf('Unable to send TCP packet (possibly unresponsive server). Disregarding\n');
                break;
            end
        end
        varargout={ok};
        
    case 'sethash',
        if numel(varargin)>=1, connection.hash=varargin{1}; end
        varargout={connection.hash};
        
    case 'setmaxlength',
        connection.maxlength=varargin{1};
        if ischar(connection.maxlength), connection.maxlength=str2double(connection.maxlength); end
        
    case 'private.ip',
        %if numel(varargin)>=1, connection.ip=varargin{1}; end
        varargout={connection.ip};
        
    case 'private.port',
        %if numel(varargin)>=1, connection.port=varargin{1}; end
        varargout={connection.port};
        
    case 'private'
        varargout={connection};
        
    otherwise,
        error('unrecognized option %s',option);
end
end


