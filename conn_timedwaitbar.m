function h=conn_timedwaitbar(a,b,str)
% CONN_TIMEDWAITBAR
%
% See also WAITBAR
%

persistent t0a t0b a0b n0 atick hsteps;
if isstr(b),
    if nargin<3||isempty(str), str=true; end
    t0a=clock;
    n0=0;
    atick=0;
    t0b=[];
    [h hsteps]=waitbar(a,b);
    if str, conn_disp(b); end
    %    hAxes = findobj(h,'type','axes');
    %    hTitle = get(hAxes,'title');
    %    set(hTitle,'color','w');
    set(findobj(h,'type','patch'),'edgecolor','none');
    set(findobj(h,'type','line'),'color',.8*get(0,'defaultuicontrolbackgroundcolor'));
    set(findobj(h,'type','axes'),'color',.8*get(0,'defaultuicontrolbackgroundcolor'));
    drawnow;
else
    t2=clock;
    t1a=etime(t2,t0a);
    n0=n0+1;
    a=min(1,max(0,a));
    if ~isempty(t0b)&&a>.2,
        t1b=etime(t2,t0b);
        t1=(1-a)*t1a.*(1-a)./(a+1e-4) + a*t1b.*(1-a)./(a-a0b+1e-4);
    else
        t1=t1a.*(1-a)./(a+1e-4);
    end
    if isempty(t0b)&&a>.1
        t0b=clock;
        a0b=a;
    end
    t1h=floor(t1/60/60);
    t1m=floor((t1-t1h*60*60)/60);
    t1s=floor((t1-t1h*60*60-t1m*60));
    waitbar(a,b);
    ht3=findobj(b,'style','togglebutton');
    if isequal(get(ht3,'value'),1),
        delete(b);
        error('<DisregardMessage>Process stopped by user');
    end
    if n0<5&&a<.2, tstr=''; 
    else
        if t1h>=1      tstr=['ETA ~ ',num2str(t1h+1,'%d'),' hours'];
        elseif t1m>=30, tstr=['ETA ~ 1 hour'];
        elseif t1m>=20, tstr=['ETA ~ 30 minutes'];
        elseif t1m>=10, tstr=['ETA ~ 20 minutes'];
        elseif t1m>=5, tstr=['ETA ~ 10 minutes'];
        elseif t1m>=1, tstr=['ETA ~ ',num2str(t1m+1,'%d'),' minutes'];
        elseif t1s>=10,tstr=['ETA ~ 1 minute'];
        else           tstr=['ETA soon'];
        end
    end
    if ~isequal(hsteps,[1 1]), tstr=['Step ',num2str(hsteps(1)),' ',tstr]; end
    %tstr=['Time left ',num2str(t1h,'%02d'),':',num2str(t1m,'%02d'),':',num2str(t1s,'%02d')];
    if nargin>2, tstr=[tstr,'  [',str,']']; end
    set(b,'name',tstr,'windowstyle','normal');
    set(findobj(b,'tag','patch'),'facecolor',(1-a)*[.2 .2 .8]+a*[.8 .8 .2]); %(1-a)*[5/6,2/6,1.5/6]+a*[1.5/6,5/6,2/6]);
    if a-atick>=.01, 
        atick=a;
        hax=findobj(b,'tag','conn_timedwaitbar_plotmatrix');
        if ~isempty(hax), conn_menu_plotmatrix('',hax(1),[1 1 10 1+9*a]); end %,[],'colormap',[.95 .95 .95; (1-a)*[.8 .2 .2]+a*[.8 .8 .2]]); %(1-a)*[5/6,2/6,1.5/6]+a*[1.5/6,5/6,2/6]]);
    end
    h=b;
end
%close(hw);
end

function [ht,steps] = waitbar(a,b)
global CONN_gui;
if isempty(CONN_gui)||~isfield(CONN_gui,'font_offset'), CONN_gui.font_offset=0; end
if isstr(b)
    if 0,%~isempty(regexp(b,'^Step \d+/\d+:\s*'))
        b1=regexp(b,'^Step (\d+)/(\d+):\s*','tokens','once');
        steps=str2double(b1);
        %b=regexprep(b,'^Step \d+/\d+:\s*','');
    else steps=[1 1];
    end
    color=[1 1 1];%get(0,'defaultuicontrolbackgroundcolor');
    ht=dialog('units','norm','position',[.4,.5,.3,.15],'windowstyle','modal','name','','handlevisibility','on','color',color,'tag','conn_timedwaitbar'); 
    htcancel=uicontrol('units','norm','position',[.4 .05 .2 .2],'style','togglebutton','string','Stop','callback','if get(gcbo,''value''), set(gcbo,''string'', ''Stopping...''); else set(gcbo,''string'',''Stop''); end; drawnow;');
    ha=[]; for n=1:steps(2), 
        ha(n)=axes('units','norm','position',[.3+.4*(n-1+.0)/steps(2) .445 .4*1/steps(2) .01]); 
        ht2=patch([0 1 1 0],[0 0 1 1],'k','facecolor',color*1,'edgecolor','none'); 
        if n==steps(1),    ht2=patch([0 0 0 0],[0 0 1 1],'k','facecolor',[5 2 1.5]/6,'edgecolor','none'); set(gca,'xlim',[0 1],'ylim',[0 1]); axis off; set(ht2,'tag','patch'); 
        elseif n<steps(1), ht2=patch([0 1 1 0],[0 0 1 1],'k','facecolor',[1.5 5 2]/6,'edgecolor','none'); set(gca,'xlim',[0 1],'ylim',[0 1]); axis off; 
        else               ht2=patch([0 1 1 0],[0 0 1 1],'k','facecolor',color,'edgecolor','none'); set(gca,'xlim',[0 1],'ylim',[0 1]); axis off;
        end
    end
    ht3=uicontrol('units','norm','position',[0 .6 1 .3],'style','text','string',b,'fontsize',8+CONN_gui.font_offset,'backgroundcolor',color);
    [nill,hax2]=conn_menu_plotmatrix('',ht,[-3,1,10,1],[.3 .4 .4 .1]);
    set(hax2,'tag','conn_timedwaitbar_plotmatrix');
else
    ht2=findobj(b,'tag','patch');
    set(ht2,'xdata',[0 a a 0]); drawnow;
end
end
