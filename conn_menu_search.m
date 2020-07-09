function conn_menu_search(hdl,event)
% internal function for gui search in listbox uicontrol objects
% adds key actions to uicontrol object:
%    control-F or command-F: find search string
%    right arrow: next match
%    left arrow: previous match
% e.g. h=uicontrol('style','listbox','keypressfcn',@conn_menu_search,...);

persistent lasthdl searchstring lastposition;

%try
    DEBUG=false;
    %if DEBUG&&~isequal(get(hdl,'style'),'listbox'), return; end
    option='';
    if isequal(event.Key,'rightarrow'),     option='next';
    elseif isequal(event.Key,'leftarrow'),  option='prev';
    elseif isequal(event.Key,'f')&&(isequal(event.Modifier,{'control'})||isequal(event.Modifier,{'command'})), option='find';
    elseif DEBUG
            disp('event.Modifier:');
            disp(event.Modifier)
            disp('event.Key:');
            disp(event.Key)
    end
    if ~isempty(option)
        if ~isequal(lasthdl,hdl),
            lasthdl=hdl;
            searchstring='';
            lastposition=[];
        end
        match1=[];
        selectall=false;
        str=get(hdl,'string');
        val=get(hdl,'value');
        if ischar(str), str=cellstr(str); end
%         if strcmp(option,'down'), 
%             if get(hdl,'max')>1&&isequal(event.Modifier,{'shift'}), match1=unique([val(:);min(numel(str),max(val)+1)]); % extend down
%             else match1=unique(min(numel(str),val+1)); % move down
%             end
%             set(hdl,'value',match1);
%         elseif strcmp(option,'up'), 
%             if get(hdl,'max')>1&&isequal(event.Modifier,{'shift'}), match1=unique([val(:);max(1,min(val)-1)]); % extend up
%             else match1=unique(max(1,val-1)); % move up
%             end
%             set(hdl,'value',match1);
%         elseif strcmp(option,'all'), 
%             if get(hdl,'max')>1, match1=1:numel(str); % select all
%             else return;
%             end
%             set(hdl,'value',match1);
        if strcmp(option,'find')||isempty(searchstring),
            searchstring='';
            %answ=inputdlg('Enter search string','',1,{''});
            %if ~isempty(answ), searchstring=answ{1}; end
            thfig=figure('units','norm','position',[.4,.4,.15,.15],'color',1*[1 1 1],'name','Search keyword matches in list','numbertitle','off','menubar','none');
            ht1a=uicontrol('style','text','units','norm','position',[.1,.75,.8,.15],'string','Enter search string:','horizontalalignment','left','backgroundcolor',1*[1 1 1],'fontweight','bold');
            ht1=uicontrol('style','edit','units','norm','position',[.1,.55,.8,.2],'string','','tooltipstring','enter case-insensitive keyword string to search (or enter regexp pattern for advance searches)');%,'callback','uiresume');
            ht2=uicontrol('style','checkbox','units','norm','position',[.1,.3,.8,.15],'string','Select all matches','value',0,'backgroundcolor',1*[1 1 1],'tooltipstring','<HTML>Checking this option will select in the original list all entries that match the above pattern<br/>Unchecking this option will jump in the original list to the first/closest entry that matches the above pattern (and you may then use left/right arrows to jump to other matches)</HTML>');
            uicontrol('style','pushbutton','string','OK','units','norm','position',[.1,.01,.38,.25],'callback','uiresume');
            uicontrol('style','pushbutton','string','Cancel','units','norm','position',[.51,.01,.38,.25],'callback','delete(gcbf)');
            uicontrol(ht1);
            uiwait(thfig);
            if ~ishandle(thfig), return; end
            searchstring=get(ht1,'string');
            selectall=get(ht2,'value');
            delete(thfig);
            drawnow;
        end
        if isempty(match1)&&~isempty(searchstring),
            if isempty(val), val=lastposition; end
            if isempty(val), val=get(hdl,'listboxtop'); end
            if isempty(val), val=numel(str); end
            match=find(cellfun('length',regexpi(str,searchstring))>0);
            if isempty(match), disp('no match found');
            else
                if strcmp(option,'next'),    [nill,idx]=min(abs(match-max(val))+1e10*(match<=max(val))); % note: remove abs if you prefer 'next' to wrap after reaching end
                elseif strcmp(option,'prev'),[nill,idx]=min(abs(min(val)-match)+1e10*(match>=min(val))); % note: remove abs if you prefer 'prev' to wrap after reaching start
                elseif selectall,            idx=1:numel(match);
                else                         [nill,idx]=min(abs(match-mean(val)));
                end
                match1=match(idx);
                set(hdl,'value',max(1,min(numel(str),match1)),'listboxtop', max([1,min(numel(str),min(match1))-1]));
            end
        end
        if ~isempty(match1)
            h=get(hdl,'callback');
            if ~isempty(h)
                drawnow;
                if isa(h,'function_handle'), feval(h,hdl,event);
                elseif iscell(h)&&isa(h{1},'function_handle'), feval(h{1},h{2:end}); 
                else eval(h);
                end
            end
        end
    end
%end
end