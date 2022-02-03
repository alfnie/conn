function [ANSW]=conn_menu_inputdlg(varargin) %PROMPT,NAME,NUMLINES,DEFAULTANSWER,OPTIONS
% internal fcn
if numel(varargin)<1||isempty(varargin{1}), PROMPT={'Input:'}; else PROMPT=varargin{1}; end
if ~iscell(PROMPT), PROMPT={PROMPT}; end
    N=numel(PROMPT);

if N<=16&&(numel(varargin)<3||isequal(varargin{3},1))
    if numel(varargin)<2||isempty(varargin{2}), NAME=''; else NAME=varargin{2}; end
    if numel(varargin)<4||isempty(varargin{4}), DEFAULTANSWER={''}; else DEFAULTANSWER=varargin{4}; end
    if ~iscell(DEFAULTANSWER), DEFAULTANSWER={DEFAULTANSWER}; end
    if numel(DEFAULTANSWER)<N, DEFAULTANSWER=DEFAULTANSWER(min(numel(DEFAULTANSWER), 1:N)); end
    scale=.5+.5*N;
    dy=(1-.15/scale-.45/scale)/N/2;
    ANSW={};
    thfig=figure('units','norm','position',[.4,.6-.1*scale,.35,.15*scale],'color',1*[1 1 1],'name',NAME,'numbertitle','off','menubar','none');
    for n1=1:N
        ht1a(n1)=uicontrol('style','text','units','norm','position',[.1,.45/scale+(2*(N+1-n1)-1)*dy,.8,.7*dy],'string',char(PROMPT{n1}),'horizontalalignment','left','backgroundcolor',1*[1 1 1],'fontweight','bold');
        ht1(n1)=uicontrol('style','edit','units','norm','position',[.11,.45/scale+(2*(N+1-n1)-2)*dy,.78,dy],'string',char(DEFAULTANSWER{n1}),'tooltipstring',char(PROMPT{n1}));
    end
    if N==1, set(ht1,'horizontalalignment','center'); else set(ht1,'horizontalalignment','left'); end
    uicontrol('style','pushbutton','string','OK','units','norm','position',[.1,.01,.38,.25/scale],'callback','uiresume');
    uicontrol('style','pushbutton','string','Cancel','units','norm','position',[.51,.01,.38,.25/scale],'callback','delete(gcbf)');
    uicontrol(ht1);
    while isempty(ANSW),
        uiwait(thfig);
        if ~ishandle(thfig), ANSW={}; return; end
        ANSW=get(ht1,'string');
        if ~iscell(ANSW), ANSW={ANSW}; end
    end
    delete(thfig);
    drawnow;
else
    ANSW=inputdlg(varargin{:});
end


2*(.45+.2*N)/(1+N)