function conn_msghelp(option,varargin)
persistent dates dates_num titles msgs ids dlg titles_fmt selected keys sortby;
global CONN_gui;
if isempty(CONN_gui)||~isfield(CONN_gui,'font_offset'), CONN_gui.font_offset=0; end

if isempty(titles)
    if conn_existfile(conn_prepend('',which(mfilename),'.mat')), load(conn_prepend('',which(mfilename),'.mat'),'strall','titles','msgs','ids','dates','dates_num','titles_fmt');
    else conn_msghelp compile;
    end
    selected=1:numel(titles);
    keys={};
    sortby=1;
end
if ~nargin, option='init'; end

switch(option)

    case 'init'
        conn_msghelp('show',1);
        conn_msghelp('key');
        conn_msghelp('show',ceil(numel(msgs)*rand^2));
        
    case 'compile'
        disp('Compiling message database. Please wait...');
        %strall=urlread('https://www.nitrc.org/forum/forum.php?set=custom&forum_id=1144&style=flat&max_rows=10000');
        if ~isempty(which('webread')), 
            if 0, % stopped working after too many posts...
                strall=webread('https://www.nitrc.org/forum/forum.php?set=custom&forum_id=1144&style=flat&max_rows=10000',weboptions('Timeout',60));
            else
                strall=''; ntotal=[];
                for nblock=1:100
                    if isempty(ntotal), fprintf('reading block %d\n',nblock);
                    else fprintf('reading block %d/%d\n',nblock,ceil(ntotal/1000));
                    end
                    temp='';
                    try, temp=webread(['https://www.nitrc.org/forum/forum.php?set=custom&forum_id=1144&style=flat&max_rows=1000&offset=',num2str(1000*(nblock-1))],weboptions('Timeout',60)); end
                    if isempty(temp), break; end
                    strall=[strall,temp];
                    if isempty(ntotal), try, ntotal=str2num(char(regexp(temp,'Showing \d+-\d+ of (\d+)','tokens','once'))); end; end
                    if ~isempty(ntotal)&nblock*1000>=ntotal, break; end
                end
            end
        else strall=urlread('https://www.nitrc.org/forum/forum.php?set=custom&forum_id=1144&style=flat&max_rows=10000');
        end
        str=regexprep(strall,{'<div class="quote">(.*?)</div>|<div class="attachment">(.*?)</div>|<a href(.*?)</a>'},{'$1'});
        str=regexprep(str,'</?td.*?>|</?br>|<!.*?>|</?span.*?>|</?strong>|</?tr.*?>|</?table.*?>','');
        str=regexprep(str,'<xml.*?>.*?</xml>|<!--[if.*?<![endif','');
        msg=regexp(str,'<div class="forum-post[^"]*"><div class="header">((Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)\s+\d+,\s*\d+).*?</div><div class="subject">(.*?)</div><div class="body">(.*?)</div>.*?<div class="footer">.*?msg_id=(\d+)','tokens');
        if numel(msg)<1e3, error('There was a problem reading from the NITRC forum site. Please try again later'); end
        [dates,titles,msgs,ids]=cellfun(@(x)deal(x{:}),msg,'uni',0);
        %str=regexprep(strall,{'<div class="quote">(.*?)</div>','<div class="bbcode.*?">(.*?)</div>'},{'$1','<DIV CLASS="STARTMESSAGE">$1</DIV>'});
        %str=regexprep(str,'</?td.*?>|</?br>|</?div.*?>|<!.*?>|</?span.*?>|</?strong>|</?tr.*?>|</?table.*?>','');
        %str=regexprep(str,'<xml.*?>.*?</xml>|<!--[if.*?<![endif|','');
        %msg=regexp(str,'>([^<]*)<DIV CLASS="STARTMESSAGE">(.*?)</DIV>.*?msg_id=(\d+).*?((Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)\s+\d+,\s*\d+)','tokens');
        %if numel(msg)<1e3, error('There was a problem reading from the NITRC forum site. Please try again later'); end
        %[titles,msgs,ids,dates]=cellfun(@(x)deal(x{:}),msg,'uni',0);
        titles=regexprep(titles,{'&quot;','&nbsp;','&gt;','&lt;','&amp;'},{'"',' ','>','<','&'});
        msgs=regexprep(msgs,{'&quot;','&nbsp;','&gt;','&lt;','&amp;'},{'"',' ','>','<','&'});
        titles=regexprep(titles,'^\s+|\s+$|<.*?>|\n','');
        titles_fmt=cellfun(@(a,b)sprintf('%s (%s)',a,regexprep(b,'\s*\d+,','')),titles,dates,'uni',0);
        dates_num=reshape(datenum(dates),size(dates));
        save(conn_prepend('',which(mfilename),'.mat'),'strall','titles','msgs','ids','dates','dates_num','titles_fmt');
        fprintf('Done (%d entries)\n',numel(msgs));

    case 'show',
        if nargin>1&&~isempty(varargin{1}), kmsg=varargin{1};
        else kmsg=ceil(numel(msgs)*rand); 
        end
        kmsg=max(1,min(numel(selected),kmsg));
        dlg.fig=findobj(0,'tag','conn_msghelp');
        if isempty(dlg.fig), 
            dlg.fig=figure('units','norm','position',[.2,.1,.6,.8],'menubar','none','numbertitle','off','name','Support questions search','color',[1 1 1],'tag','conn_msghelp'); 
            bg=.9*[1 1 1];
            uicontrol(dlg.fig,'style','frame','units','norm','position',[0,.85,1,.15],'backgroundcolor',bg,'foregroundcolor',bg);
            %uicontrol(dlg.fig,'units','norm','position',[.1 .95 .8 .025],'style','text','string','Search :','backgroundcolor',bg,'fontweight','bold','horizontalalignment','left','fontsize',CONN_gui.font_offset+10);
            dlg.key=uicontrol(dlg.fig,'units','norm','position',[.1 .90 .7 .05],'style','edit','max',1,'backgroundcolor',bg,'horizontalalignment','left','fontsize',CONN_gui.font_offset+10,'tooltipstring','<HTML>Enter search keywords <br/> - Enter words or partial words to match (e.g. <i>analys</i>)<br/> - Enter multiple keywords (separated by spaces) to match only posts containing <i>all</i> keywords (e.g. <i>artifact motion</i>)<br/> - Use single quotes to search for exact word matches (no partial-word matches) (e.g. <i>''art''</i>) <br/> - Use double-quotes to search for a multi-word keyword (e.g. <i>"motion artifact"</i>) <br/> - Use regexp strings for more complex search commands (e.g. <i>t.?test</i>)</HTML>','callback','conn_msghelp(''key'')');
            uicontrol(dlg.fig,'units','norm','position',[.8 .90 .1 .05],'style','pushbutton','string','Search','fontsize',CONN_gui.font_offset+8,'tooltipstring','Search database of support questions/answers','callback','conn_msghelp(''key'')');
            dlg.titlelist=uicontrol(dlg.fig,'units','norm','position',[.1 .78 .8 .025],'style','text','string','Posts:','backgroundcolor','w','fontsize',CONN_gui.font_offset+10,'fontweight','bold','horizontalalignment','left');
            dlg.list=uicontrol(dlg.fig,'units','norm','position',[.1 .50 .8 .27],'style','listbox','max',1,'fontname','monospaced','fontsize',CONN_gui.font_offset+8,'tooltipstring','select post','callback','conn_msghelp(''show'',get(gcbo,''value''))');
            dlg.sort=uicontrol(dlg.fig,'units','norm','position',[.7 .46 .2 .04],'style','pushbutton','string','sorted by relevance','fontsize',CONN_gui.font_offset+8,'tooltipstring','Switch between ''sorted by relevance'' and ''sorted by date''','callback','conn_msghelp(''sort'')');
            dlg.title=uicontrol(dlg.fig,'units','norm','position',[.1 .40 .8 .025],'style','text','string','','backgroundcolor','w','horizontalalignment','left','fontsize',CONN_gui.font_offset+10,'fontweight','bold');
            dlg.box=uicontrol(dlg.fig,'units','norm','position',[.1 .1 .8 .29],'style','listbox','max',2,'string','','backgroundcolor','w','horizontalalignment','left','fontsize',CONN_gui.font_offset+8);
            dlg.goto=uicontrol(dlg.fig,'units','norm','position',[.7 .06 .2 .04],'style','pushbutton','string','original post','fontsize',CONN_gui.font_offset+8,'tooltipstring','See this post in the NITRC CONN Forum website');
            set(dlg.list,'string',titles_fmt(selected));
            uicontrol(dlg.key);
        else
            dlg.fig=dlg.fig(1);
        end
        if ~isempty(selected)
            set(dlg.fig,'pointer','watch'); 
            set(dlg.box,'string',''); 
            drawnow;
            imsg=selected(kmsg);
            strdate=dates{imsg};
            numdate=dates_num(imsg);
            strid=ids{imsg};
            strtitle=titles_fmt{imsg};
            str=msgs{imsg};
            str=regexprep(str,'Originally posted by(.*?:)','\n\nOriginally posted by$1\n');
            str=regexp(str,'[\r\n]','split');
            str=regexprep(str,{'(CONN|conn|v\.|\s)(\d\d[a-z])(\W)'},{'<b>$1$2</b>$3'});
            for n=1:numel(keys), if numel(keys{n})>3, str=regexprep(str,keys{n},'<b><FONT color=rgb(0,0,255)>$1</FONT></b>','ignorecase'); end; end
            str=regexprep(str,{'(.*)'},{'<HTML>$1</HTML>'});
            idx=strmatch('<HTML>Originally posted by',str);
            if ~isempty(idx), str(idx(1):end)=regexprep(str(idx(1):end),'<HTML>(.*)</HTML>','<HTML><FONT color=rgb(100,100,100)>$1</FONT></HTML>'); end
            if ~ishandle(dlg.fig), return; end
            set(dlg.title,'string',strtitle);
            set(dlg.box,'string',str,'value',[]);
            set(dlg.goto,'callback',sprintf('conn gui_help url http://www.nitrc.org/forum/message.php?msg_id=%s',strid));
            set(dlg.list,'value',kmsg);
            set([dlg.title dlg.box dlg.goto],'visible','on');
            set(dlg.fig,'pointer','arrow');
        else
            set([dlg.title dlg.box dlg.goto],'visible','off');
        end
              
    case 'sort'
        sortby=1+(sortby==1);
        conn_msghelp('key');
        
    case 'key'
        if ~ishandle(dlg.fig), return; end
        str=get(dlg.key,'string');
        
        str=regexprep(str,'"(.*?)"|([^\s]+)','<separator>$1<separator>');
        str=regexprep(str,'''(.*?)''','\\<$1\\>');
        keys=regexp(str,'<separator>','split');
        keys=keys(cellfun('length',keys)>0&~cellfun(@(x)all(x==' '),keys));
        keys=cellfun(@(x)['(' x ')'],keys,'uni',0);
        ok=zeros(numel(keys),numel(msgs));
        for n=1:numel(keys)
            ok(n,:)=2*cellfun('length',regexpi(titles_fmt,keys{n}))+cellfun('length',regexpi(msgs,keys{n}));
        end
        selected=find(all(ok,1));
        if isempty(keys)||isempty(selected), thissortby=2; set(dlg.sort,'visible','off');
        else thissortby=sortby; set(dlg.sort,'visible','on'); 
        end
        switch(thissortby)
            case 1, % sort by relevance
                [nill,idx]=sort(prod(ok(:,selected),1)+1e-10*selected,'descend');
                set(dlg.sort,'string','sorted by relevance');
            case 2, % sort by date
                [nill,idx]=sort(dates_num(selected)+1e-10*selected,'descend');
                set(dlg.sort,'string','sorted by date');
        end
        selected=selected(idx);
        set(dlg.list,'string',titles_fmt(selected),'value',max(1,min(numel(selected), get(dlg.list,'value'))));
        set(dlg.titlelist,'string',sprintf('Posts: (%d matching records)',numel(selected)));
        conn_msghelp('show',1);
end
end
                        