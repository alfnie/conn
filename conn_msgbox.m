function [h,h2] = conn_msgbox(txt,title,ok,lowerquarter)
global CONN_gui;
if isempty(CONN_gui)||~isfield(CONN_gui,'font_offset'), CONN_gui.font_offset=0; end

if nargin<2||isempty(title), title=''; end
if nargin<3||isempty(ok), ok=false; end
if nargin<4||isempty(lowerquarter), lowerquarter=false; end
if ok==2, bg=.925*[1 1 1]; fg=[.25 0 0]; ws='modal'; else bg=.925*[1 1 1]; fg=[0 0 0]; ws='normal'; end
h2=[];
% if ok==2, h=dialog('units','norm','position',[.5 .7 .2 .2],'color',bg,'menubar','none','numbertitle','off','name',title,'resize','off','units','pixels');
% else h=figure('units','norm','position',[.5 .7 .2 .2],'color',bg,'menubar','none','numbertitle','off','name',title,'resize','off','units','pixels');
% end
h=figure('units','norm','position',[.5 .75 .2 .2],'color',bg,'menubar','none','numbertitle','off','name',title,'resize','off','units','pixels');
if ok>0 % wait for user confirmation
    ha=uicontrol('style','text','units','norm','position',[0 .35 1 .5],'backgroundcolor',bg,'horizontalalignment','center','string',txt,'units','pixels','fontsize',9+CONN_gui.font_offset,'foregroundcolor',fg,'parent',h);
    hb=uicontrol('style','pushbutton','units','norm','position',[.25 .05 .5 .20],'string','Continue','callback','uiresume(gcbf)','parent',h,'visible','off');
    hext=get(ha,'extent');
    hext2=max([150 60],hext(end-1:end)+[60 90]);
    hpos=get(h,'position');
    set(h,'position',[hpos(1)-hext2(1)/2,hpos(2)-hext2(2)/2,hext2(1),hext2(2)]);
    set(ha,'position',[30 60 hext(end-1:end)]);
    %uicontrol(hb);
    %if ok==2, conn_menu_plotmatrix('',h,10,[.25 .25 .5 .1]); end
    %set(h,'windowstyle',ws);
    set(hb,'visible','on');
    if ok==2, conn('modalfig',h);end % highlighted message
    %drawnow;
    checkdesktop=true; try, checkdesktop=checkdesktop&usejava('awt'); end 
    conn_disp(char(txt));
    if ~checkdesktop, if ok==2, fprintf(2,'ERROR: CONN requires user interaction to continue. Please re-run locally and with a user-display available\n'); end
    else uiwait(h);
    end
    %if ok==2, [nill,hc]=conn_menu_plotmatrix('',h,10,[.25 .25 .5 .1]); delete(hc(ishandle(hc))); end
    delete(h(ishandle(h)));
else % do not wait for user confirmation
    ha=uicontrol('style','text','units','norm','position',[0 0 1 1],'backgroundcolor',bg,'horizontalalignment','center','string',txt,'units','pixels','fontsize',9+CONN_gui.font_offset,'foregroundcolor',fg,'parent',h);
    hext=get(ha,'extent');
    hext2=hext(end-1:end)+[60 60];
    hpos=get(h,'position');
    set(h,'position',[hpos(1)-hext2(1)/2,hpos(2)-hext2(2)/2,hext2(1),hext2(2)]);
    set(ha,'position',[30 30 hext(end-1:end)]);
    if lowerquarter, h2=uicontrol('style','text','units','norm','position',[.1 .10 .8 .15],'backgroundcolor',.9*bg,'horizontalalignment','left','string','','fontsize',6+CONN_gui.font_offset,'foregroundcolor',.25*[1 1 1],'parent',h); end
    if ok==0, conn_disp(char(txt)); end
    drawnow;
end

