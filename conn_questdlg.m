function opt = conn_questdlg(txt,title,varargin)
global CONN_gui;
if isempty(CONN_gui)||~isfield(CONN_gui,'font_offset'), conn_font_init; end

opt='';
if nargin<2, title=''; end
if nargin<3, varargin={'Continue'}; end
tooltip={};
if ~isempty(varargin), 
    n=find(strcmp(lower(varargin),'tooltipstring'),1);
    if ~isempty(n), 
        tooltip=varargin(n+1:end);
        varargin=varargin(1:n-1);
    end
end
bg=.925*[1 1 1]; fg=[.25 0 0];
Nc=max(1,numel(varargin)-1);
checkdesktop=true; try, checkdesktop=checkdesktop&usejava('awt'); end
if ~checkdesktop, fprintf(2,'ERROR: CONN requires user interaction to continue. Please re-run locally and with a display available\n'); end
h=conn_figure('units','norm','position',[.5 .7 .4 .2],'color',bg,'menubar','none','numbertitle','off','name',title,'resize','off');
set(h,'units','pixels');
% sh=get(h,'position');
% ht=axes('units','pixels','position',[0 2 max(1,sh(3)-0) max(1,sh(4)-4)]);
% set(ht,'color',bg,'xtick',[],'ytick',[],'xcolor',bg,'ycolor',bg);
ttxt=txt;
if isempty(ttxt)
    dovert=true;
    ttxt=varargin(1:Nc);
else
    dovert=false;
    if ~iscell(ttxt),ttxt={ttxt}; end
    [nill,iv]=max(cellfun('length',varargin(1:Nc))); ttxt2=repmat(varargin{iv},1,2*Nc);
    %ttxt2=repmat([varargin{1:Nc}],1,2);
    if numel(ttxt2)>numel(ttxt{end}), ttxt{end}=ttxt2; end
end
selected=false(1,Nc);
if dovert
    pos=get(h,'position');
    set(h,'position',pos+[0 -20*Nc 0 40*(Nc+1)-pos(4)]);
    for nc=1:Nc,
        hb(nc)=uicontrol('style','pushbutton','units','norm','position',[.05 ((Nc-nc+1)-.5+.01)/(Nc+1) .9 .99/(Nc+1)],'string',varargin{nc},'callback',['set(gcbf,''userdata'',',num2str(nc),');uiresume(gcbf)'],'parent',h,'visible','off');
        if strcmp(varargin{nc},varargin{end}),set(hb(nc),'fontweight','bold'); selected(nc)=true; end
        if ~isempty(tooltip), set(hb(nc),'tooltipstring',tooltip{min(numel(tooltip),nc)}); end
    end
    set(hb,'units','pixels');
    hext=get(hb,'extent');
    if iscell(hext), hext=cat(1,hext{:}); end
    hext=max(hext(:,end-1:end),[],1);
    set(hb,'units','norm');
    set(h,'units','pixels');
    hpos=get(h,'position');
    set(h,'position',[hpos(1)-1.5*hext(1)/2,hpos(2),1.5*hext(1),hpos(4)]);
    set(h,'units','norm');
    %set(h,'position',[hpos(1)-hext2(1)/2,hpos(2),hext2(1),hpos(4)]);
    %set(ha,'position',[30 60 hext(end-1:end)],'string',txt);
    %set(hb,'units','pixels'); for nc=1:numel(hb), ipos=get(hb(nc),'position'); set(hb(nc),'position',[ipos(1:3) ipos0(4)]); end
else
    ha=uicontrol('style','text','units','norm','position',[0 .35 1 .5],'backgroundcolor',bg,'horizontalalignment','center','string',ttxt,'units','pixels','fontsize',9+CONN_gui.font_offset,'foregroundcolor',fg,'parent',h);
    for nc=1:Nc,
        hb(nc)=uicontrol('style','pushbutton','units','norm','position',[(nc-.5+.05)/(Nc+1) .05 .9/(Nc+1) .15],'string',varargin{nc},'callback',['set(gcbf,''userdata'',',num2str(nc),');uiresume(gcbf)'],'parent',h,'visible','off');
        if strcmp(varargin{nc},varargin{end}),set(hb(nc),'fontweight','bold'); selected(nc)=true; end
        if ~isempty(tooltip), set(hb(nc),'tooltipstring',tooltip{min(numel(tooltip),nc)}); end
    end
    set(hb,'units','pixels');
    ipos0=get(hb(1),'position');
    set(hb,'units','norm');
    hext=get(ha,'extent')+[0 0 0 30];
    hext2=max([100*Nc 60],hext(end-1:end)+[60 90]);
    hpos=get(h,'position');
    set(h,'position',[hpos(1)-hext2(1)/2,hpos(2)-hext2(2)/2,hext2(1),hext2(2)]);
    set(ha,'position',[15 60 hext(end-1)+15 hext(end)],'string',txt);
    set(hb,'units','pixels'); for nc=1:numel(hb), ipos=get(hb(nc),'position'); set(hb(nc),'position',[ipos(1:3) ipos0(4)]); end
end
conn('modalfig',h);
set(h,'userdata','');
set(hb,'visible','on');
if any(selected), uicontrol(hb(find(selected,1))); end
uiwait(h);
if ishandle(h), nopt=get(h,'userdata'); if ~isempty(nopt), opt=varargin{nopt}; end; delete(h); end
