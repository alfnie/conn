function varargout=conn_disp(varargin)
% internal function: logged "disp"

global CONN_x CONN_h CONN_gui;
persistent MAXHISTORY SAVELOG CWCOPY TBSPACE PORTCOMM lastdate lastlogfile;

if isempty(lastdate), lastdate=''; end
if isempty(lastlogfile), lastlogfile=''; end
if isempty(MAXHISTORY), MAXHISTORY=1e4; end
if isempty(SAVELOG), SAVELOG=true; end
if isempty(CWCOPY), CWCOPY=true; end
if isempty(TBSPACE), TBSPACE=0; end
if isempty(PORTCOMM), PORTCOMM=false; end
silent=false;
varargout={};
%try
mirrorscreen=false;
if nargin>=1&&ischar(varargin{1})&&size(varargin{1},1)==1&&~isempty(regexp(varargin{1},'^__'))
    switch(lower(varargin{1}))
        case '__maxhistory'
            if nargin>1, MAXHISTORY=varargin{2}; end
            if nargout>0, varargout={MAXHISTORY}; end
            return
        case '__savelog'
            if nargin>1, SAVELOG=varargin{2}; end
            if nargout>0, varargout={SAVELOG}; end
            return
        case '__cwcopy'
            if nargin>1, CWCOPY=varargin{2}; end
            if nargout>0, varargout={CWCOPY}; end
            return
        case '__tbspace'
            if nargin>1, TBSPACE=varargin{2}; end
            if nargout>0, varargout={TBSPACE}; end
            return
        case '__portcomm',
            if nargin>1, PORTCOMM=varargin{2}; end
            if nargout>0, varargout={PORTCOMM}; end
            return
        case '__lastdate'
            if nargin>1, lastdate=varargin{2}; end
            if nargout>0, varargout={lastdate}; end
            return
        case '__exit'
            if isfield(CONN_h,'screen')&&isfield(CONN_h.screen,'hlog')&&ishandle(CONN_h.screen.hlog),delete(CONN_h.screen.hlog); end
            return
        case '__clear'
            try
                if isfield(CONN_x,'gui')&&(isnumeric(CONN_x.gui)&&CONN_x.gui || isfield(CONN_x.gui,'display')&&CONN_x.gui.display)&&isfield(CONN_h,'screen')&&isfield(CONN_h.screen,'hfig')&&ishandle(CONN_h.screen.hfig)
                    str={' '};
                    set(CONN_h.screen.hlogstr,'string',str,'value',numel(str),'listboxtop', numel(str),'backgroundcolor',CONN_gui.backgroundcolor,'foregroundcolor',CONN_gui.fontcolorA);
                    drawnow;
                    set(CONN_h.screen.hlogstr,'value',[]);
                    lastdate='';
                end
            end
            return
        case '__copy'
            str=get(CONN_h.screen.hlogstr,'string');
            if ~iscell(str), str=cellstr(str); end
            val=get(CONN_h.screen.hlogstr,'value');
            str=str(val);
            clipboard('copy',sprintf('%s\n',str{:}));
            return
        case '__show'
            conn_disp;
            figure(CONN_h.screen.hlog);
            return
        case '__init'
            conn_disp;
            conn_disp('__restart');
            return
        case '__nolog'
            savelog=SAVELOG;
            SAVELOG=false;
            conn_disp(varargin{2:end});
            SAVELOG=savelog;
            return
        case '__restart'
            try,
                if isfield(CONN_x,'gui')&&(isnumeric(CONN_x.gui)&&CONN_x.gui || isfield(CONN_x.gui,'display')&&CONN_x.gui.display)&&isfield(CONN_h,'screen')&&isfield(CONN_h.screen,'hfig')&&ishandle(CONN_h.screen.hfig)
                    filename=fullfile(conn_prepend('',conn_projectmanager('projectfile'),''),'logfile.txt');
                    if conn_existfile(filename)
                        str=conn_fileutils('fileread',filename);
                        str=regexp(str,'[\r\n]+','split');
                        str=[{' '} str];
                        if numel(str)>MAXHISTORY, iscropped=numel(str); str=str(end-MAXHISTORY+1:end); 
                        else iscropped=0;
                        end
                        %pivot=cellfun('length',regexp(str,'^\d+-\w+-\d+\s+\d+:\d+'));
                        set(CONN_h.screen.hlogstr,'string',str,'value',numel(str),'listboxtop', numel(str),'backgroundcolor',CONN_gui.backgroundcolor,'foregroundcolor',CONN_gui.fontcolorA);
                        drawnow;
                        set(CONN_h.screen.hlogstr,'value',[]);
                        lastdate='';
                        if iscropped>0, conn_msgbox({sprintf('Log loaded correctly from %s',filename),sprintf('note: log display only showing last %d lines (out of %d)',MAXHISTORY,iscropped)},'Log display',true);
                        else conn_msgbox(sprintf('Log loaded correctly from %s',filename),'Log display',true);
                        end
                        if ~isequal(lastlogfile,filename), set(CONN_h.screen.hlog,'name',sprintf('CONN log history (%s)',filename)); end
                        lastlogfile=filename;
                    end
                end
            end
            return
        otherwise, error('unrecognized option %s',varargin{1});
    end
end
savelog=SAVELOG;
%if isfield(CONN_gui,'isremote')&&~isempty(CONN_gui.isremote)&&CONN_gui.isremote>0, savelog=false; end
if ~(isfield(CONN_x,'filename')&&~isempty(CONN_x.filename)&&ischar(CONN_x.filename)&&isfield(CONN_x,'pobj')&&isfield(CONN_x.pobj,'isextended')&&~CONN_x.pobj.isextended), savelog=false; end
if isfield(CONN_x,'gui')&&(isnumeric(CONN_x.gui)&&CONN_x.gui || isfield(CONN_x.gui,'display')&&CONN_x.gui.display)&&isfield(CONN_h,'screen')&&isfield(CONN_h.screen,'hfig')&&ishandle(CONN_h.screen.hfig)
    mirrorscreen=true;
    if ~isfield(CONN_h,'screen')||~isfield(CONN_h.screen,'hlog')||~ishandle(CONN_h.screen.hlog),
        pos=get(CONN_h.screen.hfig,'position');
        try, fntname=CONN_gui.fontname; catch, fntname=get(0,'FixedWidthFontName'); end
        CONN_h.screen.hlog=figure('units','pixels','position',[pos(1),1,pos(3),max(200,pos(2))],'color',CONN_gui.backgroundcolor,'doublebuffer','on','tag','conn_logwindow','name','CONN log history','numbertitle','off','menubar','none','resize','on','interruptible','off');
        CONN_h.screen.hlogstr=uicontrol('units','norm','position',[0 0 1 1],'style','listbox','string',{' '},'horizontalalignment','left','max',2,'keypressfcn',@conn_menu_search,'fontname',fntname,'backgroundcolor',CONN_gui.backgroundcolor,'foregroundcolor',CONN_gui.fontcolorA,'fontsize',8+CONN_gui.font_offset,'tooltipstring',['<HTML>Log of CONN''s processing and analysis steps<br/> - ',CONN_gui.rightclick,'-click for additional options<br> - note: keyboard shortcuts: ''',CONN_gui.keymodifier,'-F'' finds match to keyword; ''right arrow'' next match; ''left arrow'' previous match; ''',CONN_gui.keymodifier,'-A'' select all</HTML>']);
        set(CONN_h.screen.hlogstr,'units','pixels');
        tpos=get(CONN_h.screen.hlogstr,'position');
        tpos2=get(CONN_h.screen.hlogstr,'extent');
        set(CONN_h.screen.hlogstr,'units','norm');
        htb=uicontrol('style','frame','units','pixels','position',tpos+[tpos(3)-1*CONN_gui.uicontrol_border-0*18,0,1*CONN_gui.uicontrol_border+0*18-tpos(3),0],'foregroundcolor',CONN_gui.backgroundcolor,'backgroundcolor',CONN_gui.backgroundcolor,'units','norm');
        %conn_menumanager('onregion',htb,-1,get(CONN_h.screen.hlogstr,'position'),CONN_h.screen.hlogstr);
        uicontrol('style','frame','units','pixels','position',tpos+[0,0,CONN_gui.uicontrol_border-tpos(3),0],'foregroundcolor',CONN_gui.backgroundcolor,'backgroundcolor',CONN_gui.backgroundcolor,'units','norm');
        uicontrol('style','frame','units','pixels','position',tpos+[0,0,0,1*CONN_gui.uicontrol_border+0*tpos2(4)-tpos(4)],'foregroundcolor',CONN_gui.backgroundcolor,'backgroundcolor',CONN_gui.backgroundcolor,'units','norm');
        uicontrol('style','frame','units','pixels','position',tpos+[0,tpos(4)-CONN_gui.uicontrol_border,0,CONN_gui.uicontrol_border-tpos(4)],'foregroundcolor',CONN_gui.backgroundcolor,'backgroundcolor',CONN_gui.backgroundcolor,'units','norm');
        hc1=uicontextmenu;
        uimenu(hc1,'Label','Clear display','callback',@(varargin)conn_disp('__clear'));
        uimenu(hc1,'Label','Reload log history of currently open project','callback',@(varargin)conn_disp('__restart'));
        uimenu(hc1,'Label','Copy selected log history lines to clipboard','callback',@(varargin)conn_disp('__copy'));
        uimenu(hc1,'Label','Export selected log history lines to file','callback',@(varargin)conn_exportlist(CONN_h.screen.hlogstr,[],[],true));
        set(CONN_h.screen.hlogstr,'uicontextmenu',hc1)
        figure(CONN_h.screen.hfig);
    end
elseif isfield(CONN_h,'screen')&&isfield(CONN_h.screen,'hlog')&&ishandle(CONN_h.screen.hlog)&&isfield(CONN_h.screen,'hlogstr')&&ishandle(CONN_h.screen.hlogstr), mirrorscreen=true;
end
if nargin>=1
    if mirrorscreen||savelog||PORTCOMM||nargout>0
        if iscell(varargin{1})&&all(cellfun(@ischar,varargin{1}))
            newstr=varargin{1}; newstr=newstr(cellfun('length',newstr)>0);
        elseif ~ischar(varargin{1})
            newstr=regexp(evalc('disp(varargin{1})'),'[\r\n]+','split'); newstr=newstr(cellfun('length',newstr)>0);
        elseif strcmp(varargin{1},'fprintf'),
            if nargin>1&&ischar(varargin{2}), newstr=regexp(sprintf(varargin{2:end}),'[\r\n]+','split'); newstr=newstr(cellfun('length',newstr)>0);
            else newstr='';
            end
        elseif strcmp(varargin{1},'struct')
            newstr=conn_disp_struct(varargin{2:end});
        else newstr=varargin{1};
        end
        if ~isempty(newstr)
            if ischar(newstr), newstr=cellstr(newstr); end
            newdate=datestr(now,'yyyy-mmm-dd HH:MM');
            newstr=regexprep(newstr,{'(initializing)?[\s\.]*please[\s\.]*wait[\s\.]*','[^\r\n]*\b[^\r\n]*'},{'',''},'ignorecase');
            str={};
            for n=1:numel(newstr),
                if ~isempty(newstr{n})
                    if isempty(lastdate)||~strcmp(lastdate,newdate), str{end+1}=[newdate '  : ' repmat(' ',1,TBSPACE) newstr{n}]; lastdate=newdate;
                    else str{end+1}=[repmat(' ',1,numel(newdate)) '    ' repmat(' ',1,TBSPACE) newstr{n}];
                    end
                end
            end
            if mirrorscreen&&~isempty(str)
                try
                    oldstr=get(CONN_h.screen.hlogstr,'string');
                    if isempty(oldstr), oldstr={}; end
                    if ischar(oldstr), oldstr=cellstr(oldstr); end
                    oldstr=[reshape(oldstr,[],1);reshape(str,[],1)];
                    if numel(oldstr)>MAXHISTORY, oldstr=oldstr(end-MAXHISTORY+1:end); end
                    set(CONN_h.screen.hlogstr,'string',oldstr,'value',numel(oldstr),'listboxtop', numel(oldstr));
                    drawnow;
                    set(CONN_h.screen.hlogstr,'value',[]);
                end
            end
            if savelog&&~isempty(str)
                try
                    if isfield(CONN_x,'pobj')&&isfield(CONN_x.pobj,'holdsdata')&&~CONN_x.pobj.holdsdata&&isfield(CONN_x.pobj,'isextended')&&~CONN_x.pobj.isextended
                        filename=conn_prepend('',conn_projectmanager('projectfile'),'.log');
                        filepath=fileparts(filename);
                    else
                        filepath=conn_prepend('',conn_projectmanager('projectfile'),'');
                        filename=fullfile(filepath,'logfile.txt');
                    end
                    try, 
                        conn_fileutils('fileappend',filename, str);
                    catch
                        conn_fileutils('mkdir',filepath);
                        conn_fileutils('fileappend',filename, str);
                    end
                    %try
                    %    fh=fopen(filename,'at');
                    %catch
                    %    [ok,nill]=mkdir(filepath);
                    %    fh=fopen(filename,'at');
                    %end
                    %for n=1:numel(str), fprintf(fh,'%s\n',str{n}); end
                    %fclose(fh);
                    if mirrorscreen&&~isequal(lastlogfile,filename), try, set(CONN_h.screen.hlog,'name',sprintf('CONN log history (%s)',filename)); end; end
                    if ~isequal(lastlogfile,filename),
                        tname=conn_dirn(filename);
                        if tname.bytes>104857600 % rename if above 100Mb limit
                            newfilename=conn_prepend('',filename,[datestr(now,'_yyyy_mm_dd_HHMMSSFFF'),'.txt']);
                            try, conn_fileutils('renamefile',filename,newfilename); end
                            %if ispc, [ok,msg]=system(sprintf('ren "%s" "%s"',filename,newfilename));
                            %else [ok,msg]=system(sprintf('mv ''%s'' ''%s''',filename,newfilename));
                            %end
                        end
                    end
                    lastlogfile=filename;
                end
            elseif ~savelog&&mirrorscreen, try, set(CONN_h.screen.hlog,'name','CONN log history'); end; lastlogfile=''; 
            end
            if PORTCOMM&&~isempty(str) % from server
                conn_tcpip('write',struct('type','status','id','unknown','msg',{str}));
                try, tcodes=conn_tcpip('peek'); 
                catch, fprintf('warning: communications failure; disconnecting conn_disp tunneling\n'); PORTCOMM=false; tcodes=[];
                end
                if ~isempty(tcodes)&&ismember('STOP',tcodes), error('<DisregardMessage>Process stopped by user'); end
            end
            if nargout>0&&~isempty(newstr), varargout={newstr}; end
        end
    end
    if ~mirrorscreen||CWCOPY
        if iscell(varargin{1})&&all(cellfun(@ischar,varargin{1})), newstr=varargin{1}; newstr=newstr(cellfun('length',newstr)>0); disp(sprintf('%s\n',newstr{:})); 
        elseif ~ischar(varargin{1}), disp(varargin{:});
        elseif strcmp(varargin{1},'fprintf'), fprintf(varargin{2:end});
        elseif strcmp(varargin{1},'struct'), disp(varargin{2});
        else disp(varargin{:});
        end
    end
elseif mirrorscreen
    try, set(CONN_h.screen.hlogstr,'backgroundcolor',CONN_gui.backgroundcolor,'foregroundcolor',CONN_gui.fontcolorA); end
end
%end
end

function out=conn_disp_struct(matlabbatch,str)
out={};
if nargin<2, str=''; end
try
    if isempty(matlabbatch),
    elseif iscell(matlabbatch)&&ischar(matlabbatch{1})
        if numel(matlabbatch)>4
            out{end+1}=sprintf('%s(%d) = %s',str,1,matlabbatch{1});
            out{end+1}=sprintf('%s(%d) = %s',str,2,matlabbatch{2});
            out{end+1}=sprintf('%s(%d) = %s',str,numel(matlabbatch)-1,matlabbatch{end-1});
            out{end+1}=sprintf('%s(%d) = %s',str,numel(matlabbatch),matlabbatch{end});
        elseif numel(matlabbatch)>1
            for n=1:numel(matlabbatch)
                out{end+1}=sprintf('%s(%d) = %s',str,n,matlabbatch{n});
            end
        else
            out{end+1}=sprintf('%s = %s',str,matlabbatch{1});
        end
    elseif iscell(matlabbatch)
        if numel(matlabbatch)==1
            out=cat(2,out,conn_disp_struct(matlabbatch{1},str));
        else
            for n=1:numel(matlabbatch)
                out=cat(2,out,conn_disp_struct(matlabbatch{n},sprintf('%s(%d)',str,n)));
            end
        end
    elseif ~isempty(matlabbatch)&&isnumeric(matlabbatch)
        if numel(matlabbatch)>100
            out{end+1}=sprintf('%s = [%s %s %s %s %s ... %s %s %s %s %s]',str,mat2str(matlabbatch(1)),mat2str(matlabbatch(2)),mat2str(matlabbatch(3)),mat2str(matlabbatch(4)),mat2str(matlabbatch(5)),mat2str(matlabbatch(end-4)),mat2str(matlabbatch(end-3)),mat2str(matlabbatch(end-2)),mat2str(matlabbatch(end-1)),mat2str(matlabbatch(end)));
        else
            out{end+1}=sprintf('%s = %s',str,mat2str(matlabbatch));
        end
    elseif ~isempty(matlabbatch)&&ischar(matlabbatch)
        if numel(matlabbatch)>200, out{end+1}=sprintf('%s = %s...%s',str,1,matlabbatch(1:50),matlabbatch(end-50+1:end));
        else out{end+1}=sprintf('%s = %s',str,matlabbatch);
        end
    elseif numel(matlabbatch)>1
        for n=1:numel(matlabbatch)
            out=cat(2,out,conn_disp_struct(matlabbatch(n),sprintf('%s(%d)',str,n)));
        end
    elseif isstruct(matlabbatch)
        names=fieldnames(matlabbatch);
        for n=1:numel(names)
            out=cat(2,out,conn_disp_struct(matlabbatch.(names{n}),sprintf('%s.%s',str,names{n})));
        end
    end
end
end


