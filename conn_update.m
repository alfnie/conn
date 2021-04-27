function [updated,newver]=conn_update(thisver,force,silent)
global CONN_x

if nargin<2||isempty(force), force=false; end % 0: checks if new version available; 1: forces download of latest version any way; 2: forces download (and no questions asked)
if nargin<3||isempty(silent), silent=false; end % 1: silent mode (minimizes user display); 2: really-silent mode (no user interaction at all)
if nargin<1||isempty(thisver), try, thisver=conn('ver'); catch, thisver=[]; force=max(force,1); end; end 
if isdeployed, path0=matlabroot;
else path0=fileparts(which('conn'));
end
url0='https://www.nitrc.org';
url1='/projects/conn/';

updated=false;
newver=false;
if ~isempty(thisver),
    version_conn=regexp(thisver,'(\d.*?)\.(\w+)','tokens','once');
    if iscell(version_conn)&&length(version_conn)==2,txt=['CONN functional connectivity toolbox (',version_conn{1},'.',version_conn{2},')'];
    else txt='';
    end
else
    version_conn={};
    txt='';
end
if ~silent, conn_update_disp({txt,'Checking for updates. Please wait...'}); end
if ~isempty(which('webread')), try, s=webread([url0,url1]); ok=true; catch, ok=false; end
else [s,ok]=urlread([url0,url1]);
end

if ~ok, if ~silent, conn_update_disp(['Cannot access ',[url0,url1],' to check for updates. Please try again later'],2); end
else
    existprev=isfield(CONN_x,'ver')&&~isempty(txt);
    if isdeployed, version_url=regexp(s,['conn standalone_precompiled v\.(\d.*?)\.(\w+) ',lower(computer)],'tokens','once');
    else version_url=regexp(s,'conn v\.(\d.*?)\.(\w+)','tokens','once');
    end
    if isempty(version_url), 
        if ~silent, conn_update_disp('Sorry, automatic updates not available. Please check back later',2); end
        return;
    end
    if force||existprev, 
        if existprev, 
            if iscell(version_url) && length(version_url)==2 && (str2num(version_url{1})<str2num(version_conn{1}) || (str2num(version_url{1})==str2num(version_conn{1}) && base2dec(version_url{2},36)<=base2dec(version_conn{2},36))), samever=1; 
            else samever=0;
            end
        else samever=0;
        end
        newver=samever==0;
        if force,samever=0; end
        if samever,
            if ~silent, conn_update_disp({['Your version is up to date']},1); end
        else
            if silent==2&&force~=2,% only check / do not offer to download 
                %conn_update_disp({txt,['Update (',version_url{1},'.',version_url{2},') available at ',url0,url1]});
            else,
                %group_id=regexp(s,'frs/\?group_id=(\d+)">Files','tokens','once');
                if ~silent&&force==2, conn_update_disp({txt,['Update (',version_url{1},'.',version_url{2},') available at ',url0,url1]}); end
                if force==2, answ='Yes'; 
                elseif ~isempty(which('conn_questdlg')), answ=conn_questdlg({['Update (',version_url{1},'.',version_url{2},') available at ',url0,url1],'Download and install this update now?'},'conn update','Yes', 'No', 'release notes','Yes');
                else answ=questdlg({['Update (',version_url{1},'.',version_url{2},') available at ',url0,url1],'Download and install this update now?'},'conn update','Yes', 'No', 'release notes','Yes');
                end
                if strcmp(answ,'release notes'),
                    try
                        group_id=regexp(s,'frs/\?group_id=(\d+)">','tokens','once');
                        release_id=regexp(s,['frs/\?group_id=',group_id{1},'\&amp\;release_id=(\d+)'],'tokens','once');
                        url2=['/frs/shownotes.php?release_id=',release_id{1}];
                        conn('gui_help','url',[url0,url2]);
                    catch
                        conn('gui_help','url',[url0,url1]);
                    end
                    conn_update_disp;
                elseif strcmp(answ,'Yes'),
                    connversion={'CONN functional connectivity toolbox'};
                    hfig=findobj('tag',connversion{1});
                    if ~isempty(hfig)&&ishandle(hfig)&&~silent,
                        answ=conn_questdlg({'Warning. Your current CONN session needs to be closed to proceed with this update.','Any unsaved data will be lost. Please save your project before proceeding.','Do you want to continue?'},'conn update','Yes', 'No', 'Yes');
                    else, answ='Yes';end
                    if strcmp(answ,'Yes'),
                        group_id=regexp(s,'frs/\?group_id=(\d+)">','tokens','once');
                        url2=['/frs/?group_id=',group_id{1}];
                        conn_update_disp('Accessing download page... ');
                        if ~isempty(which('webread')), try, s2=webread([url0,url2]); ok=true; catch, ok=false; end
                        else [s2,ok]=urlread([url0,url2]);
                        end
                        if isdeployed, 
                            file_id=regexp(s2,['/frs/download\.php/(\d+)/conn',version_url{1},version_url{2},'_',lower(computer),'\.zip'],'tokens','once');
                            url3=['/frs/download.php/',file_id{1},'/conn',version_url{1},version_url{2},'_',lower(computer),'.zip'];
                        else
                            file_id=regexp(s2,['/frs/download\.php/(\d+)/conn',version_url{1},version_url{2},'\.zip'],'tokens','once');
                            url3=['/frs/download.php/',file_id{1},'/conn',version_url{1},version_url{2},'.zip'];
                        end
                        
                        try
                            if ~isempty(strfind(path0,['conn',version_conn{1},version_conn{2}])),
                                path1=regexprep(path0,['conn',version_conn{1},version_conn{2},'.*'],['conn',version_url{1},version_url{2}]);
                                [ok,ko]=mkdir(path1);
                            else
                                [path1,name1]=fileparts(path0);
                            end
                        catch
                            path1=pwd;
                        end
                        ok=silent|isequal(path0,fullfile(path1,'conn'));
%                         [nill,temp]=fileparts(path1); if length(version_conn)==2&&strcmp(temp,['conn',version_conn{1},version_conn{2}]), path1=fileparts(path1); end
%                         %[nill,temp]=fileparts(fileparts(path0));if length(version_conn)==2&&strcmp(temp,['conn',version_conn{1},version_conn{2}]),path1=fileparts(path0);else,path1=path0;end
%                         path1=fullfile(fileparts(path1),['conn',version_url{1},version_url{2}]);
%                         
%                         if existprev,
%                             [ok,ko]=mkdir(fileparts(path1),['conn',version_url{1},version_url{2}]);
%                         end
%                         ok=0;
                        while ~ok,
                            if force==2, answ=input(['Installation path for conn directory ' path1 '? (Y/N) : '],'s'); if ~isequal(answ,'Y'), path1=input(['Installation path for conn directory: '],'s'); end
                            else disp('Select a directory to install the new version'); path1=uigetdir(path1,'Select a directory to download/install CONN');
                            end
                            if isequal(path1,0),break;end
                            [ok,ko]=mkdir(path1,'conn');
                            %[ok,ko]=mkdir(fileparts(path1),['conn',version_url{1},version_url{2}]);
                            if ~ok,
                                conn_update_disp(['Unable to create directory ',path1,' - check permissions'],2);
                            end
                        end
                        if ~isequal(path1,0),
                            [ok,ko]=mkdir(path1,'conn');
                            path2=fullfile(path1,'conn');
                            if existprev,
                                askagain=conn_register('nill');
                            else
                                askagain=true;
                            end
                            if existprev&&isequal(path0,path2)
                                try,
                                    if ispc, [nill,nill]=system(sprintf('copy "%s" "%s"',fullfile(path0,'conn_font_default.dat'),fullfile(path0,'bakconn_font_default.dat')));
                                             [nill,nill]=system(sprintf('copy "%s" "%s"',fullfile(path0,'conn_jobmanager.mat'),fullfile(path0,'bakconn_jobmanager.mat')));
                                    else     [nill,nill]=system(sprintf('''cp'' ''%s'' ''%s''',fullfile(path0,'conn_font_default.dat'),fullfile(path0,'bakconn_font_default.dat')));
                                             [nill,nill]=system(sprintf('''cp'' ''%s'' ''%s''',fullfile(path0,'conn_jobmanager.mat'),fullfile(path0,'bakconn_jobmanager.mat')));
                                    end
                                end
                            end
                            
                            fileattrib(path2,'+w','a','s');
                            conn_update_disp({'CONN is downloading updates','This process may take around 5 minutes depending on your connection speed','Please wait...'});
                            unzip([url0,url3],path1);
                            
                            conn_update_disp({['CONN successfully updated (release ',version_url{1},'.',version_url{2},')'],['Installation folder ',path2]},0);
                            try
                                path2=fileparts(conn_dir(fullfile(path2,'conn.m')));
                            end
                            updated=true;
                            
                            if existprev,
                                conn forceclose;
                            end
                            if ~isdeployed&&~isequal(path0,path2)
                                if ~isempty(path0), rmpath(path0); end
                                addpath(path2);
                                rehash path;
                            end
                            if ~askagain, conn_register('donotaskagain'); end
                            if existprev
                                try,
                                    if isequal(path0,path2)
                                        if ispc, [nill,nill]=system(sprintf('copy "%s" "%s"',fullfile(path0,'bakconn_font_default.dat'),fullfile(path2,'conn_font_default.dat')));
                                                 [nill,nill]=system(sprintf('copy "%s" "%s"',fullfile(path0,'bakconn_jobmanager.mat'),fullfile(path2,'conn_jobmanager.mat')));
                                        else     [nill,nill]=system(sprintf('''cp'' ''%s'' ''%s''',fullfile(path0,'bakconn_font_default.dat'),fullfile(path2,'conn_font_default.dat')));
                                                 [nill,nill]=system(sprintf('''cp'' ''%s'' ''%s''',fullfile(path0,'bakconn_jobmanager.mat'),fullfile(path2,'conn_jobmanager.mat')));
                                        end
                                    else
                                        if ispc, [nill,nill]=system(sprintf('copy "%s" "%s"',fullfile(path0,'conn_font_default.dat'),fullfile(path2,'conn_font_default.dat')));
                                                 [nill,nill]=system(sprintf('copy "%s" "%s"',fullfile(path0,'conn_jobmanager.mat'),fullfile(path2,'conn_jobmanager.mat')));
                                        else     [nill,nill]=system(sprintf('''cp'' ''%s'' ''%s''',fullfile(path0,'conn_font_default.dat'),fullfile(path2,'conn_font_default.dat')));
                                                 [nill,nill]=system(sprintf('''cp'' ''%s'' ''%s''',fullfile(path0,'conn_jobmanager.mat'),fullfile(path2,'conn_jobmanager.mat')));
                                        end
                                    end
                                    if isequal(path0,pwd)||isequal(path2,pwd)
                                        cd('..');
                                    end
                                end
                            end
                            if isdeployed,
                                [nill,nill]=system('conn');
                            else
                                if ~isequal(path0,path2)
                                    outputpathfile=which('pathdef.m');
                                    if force==2, conn_update_disp({['Update (',version_url{1},'.',version_url{2},') successfully installed'],['Warning: MATLAB path not yet saved for future sessions.'],['To use this new release in future MATLAB sessions type:'],['   addpath ',path2]});
                                    else
                                        [outputpathfilename, outputpathpathname]=uiputfile('*.m','Save updated MATLAB path?',outputpathfile);
                                        if ~isequal(outputpathfilename,0)
                                            ko=savepath(fullfile(outputpathpathname,outputpathfilename));
                                            if ko, conn_update_disp({['Update (',version_url{1},'.',version_url{2},') successfully installed'],['Warning: Not possible to save MATLAB path. Check file permissions.'],['To use this new release in future MATLAB sessions type:'],' ',['   addpath ',path2]},1);
                                            elseif ~isempty(path0), conn_update_disp({['Update (',version_url{1},'.',version_url{2},') successfully installed'],['Saved new MATLAB path: ',path2],['You can delete if you wish the previous version at ',path0]},1); 
                                            else conn_update_disp({['Update (',version_url{1},'.',version_url{2},') successfully installed'],['Saved new MATLAB path: ',path2]},1); 
                                            end
                                        else
                                            if ~isempty(path0), conn_update_disp({['Update (',version_url{1},'.',version_url{2},') successfully installed'],['You may delete if you wish the previous version at ',path0],'Warning: MATLAB path not saved',['To use this new release in future MATLAB sessions type:'],' ',['   addpath ',path2]},1);
                                            else conn_update_disp({['Update (',version_url{1},'.',version_url{2},') successfully installed'],'Warning: MATLAB path not saved',['To use this new release in future MATLAB sessions type:'],' ',['   addpath ',path2]},1);
                                            end
                                        end
                                    end
                                end
                                try, rehash; end
                                conn;                                
                            end
                        else conn_update_disp(['Installation cancelled'],2);
                        end
                    else conn_update_disp; 
                    end
                else conn_update_disp;
                end
            end
        end
    else
        newver=true;
        if silent~=2, conn_update_disp({['CONN version ',version_url{1},'.',version_url{2},' available at ',url0,url1]},1); end
    end
end

    function conn_update_disp(varargin)
        persistent h;
        if isfield(CONN_x,'gui')&&~isstruct(CONN_x.gui)&&CONN_x.gui,
            if ~isempty(h)&&ishandle(h),close(h);end
            if ~isempty(varargin), h=conn_msgbox(varargin{1},'conn_update',varargin{2:end}); end
        else,
            if ~isempty(varargin), disp(char(varargin{1})); end
        end
        
    end

end


