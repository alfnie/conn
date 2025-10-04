
function handle=conn_menumanager(handle,varargin)

% CONN_MENUMANAGER User interface manager
%   conn_menumanager on; Initializes the figure callback manager. 
%   conn_menumanager off; Terminates the figure callback manager. 
%   h=conn_menumanager([],'argname1',argvalue1,'argname2',argvalue2,...) creates
%     a user interface control and returns a handle to it;
%   conn_menumanager(h,'argname1',argvalue1,'argname2',argvalue2,...) sets
%     user interface control properties.
%   argvalue=conn_menumanager(h,'argname') returns user interface control property 'argname'.
%
% See the code for details of valid properties.
%
% Example:
%     set(gcf,'color','k');
%     conn_menumanager;
%   	h1=conn_menumanager([],              'n',2,...
%   									'string',{'Apples','Oranges'},...
%   									'position',[.20,.70,.1,.1]);
%   	h2=conn_menumanager([],              'n',3,...
%   									'string',{'Red','Blue','Orange'},...
%   									'position',[.30,.65,.1,.15]);
%   	h3=conn_menumanager([],              'n',3,...
%   									'string',{'Green','Yellow','White'},...
%   									'position',[.40,.65,.1,.15]);
%   	h=conn_menumanager([],               'n',3,...
%   									'string',{'Fruits','Colors1','Colors2'},...
%   									'help',{'Select fruits','Select colors','Select other colors'},...
%   									'position',[.20,.81,.3,.05],...
%                                     'order','horizontal',...
%   									'callback',{h1,h2,h3} );
%     conn_menumanager(h,'on',1);
%

% alfnie@gmail.com 04/07

global CONN_x CONN_MM CONN_gui;
if ~isfield(CONN_gui,'font_offset'), conn_font_init; end
if ~isfield(CONN_gui,'waiticon'), CONN_gui.waiticon='watch'; end
if ~isfield(CONN_gui,'fontname'), CONN_gui.fontname=get(0,'FixedWidthFontName'); ; end
if nargin<1, handle='on'; end
%AVOIDTEXTBUG=CONN_gui.dounixGUIbugfix; % avoids issue with remote X display fonts not resizing correctly (note: disregards rotate field)
UIXCONTROLTOOLTIP=false; 

if ischar(handle),
    handle=lower(handle);
	switch(handle),
        case 'cursordown'
            CONN_gui.iscursordown=true;
		case {'cursormove','cursorup','cursorout'}
            if strcmp(handle,'cursorup'), CONN_gui.iscursordown=false; end
            if ~isequal(gcf,CONN_MM.gcf), return; end
            active=find(CONN_MM.ACTIVE); if isempty(active), return; end
            if CONN_MM.private.busy&&strcmp(handle,'cursormove'), CONN_MM.private.busy=mod(CONN_MM.private.busy+1,10); return; end
            CONN_MM.private.busy=1;
			p1=get(0,'pointerlocation');
            p2=get(CONN_MM.gcf,'position'); 
            p3=get(0,'screensize');
            p2(1:2)=p2(1:2)+p3(1:2)-1; % note: fix issue when connecting to external monitor/projector
            pos=(p1-p2(1:2))./p2(3:4);
            if (strcmp(handle,'cursormove')||strcmp(handle,'cursorup'))&&(sum(abs(pos-CONN_MM.restpos))<.1||now-CONN_MM.resttime<5e-5), return; end
            CONN_MM.restpos=[nan nan];
            CONN_MM.resttime=nan;
            xpos=repmat(pos,[length(active),1]);
            [md,n0]=max(sign(min(min(xpos-CONN_MM.POSITION(active,1:2),CONN_MM.POSITION(active,1:2)+CONN_MM.POSITION(active,3:4)-xpos),[],2)),[],1);
            if ~strcmp(handle,'cursorout')&&md>-.01, % move inside active area 
                n0=active(n0);
                if CONN_MM.MENU{n0}.order(1)=='h',
                    k1=max(eps,min(1, (pos(1)-CONN_MM.MENU{n0}.position(1))./CONN_MM.MENU{n0}.position(3) ));
                    n1=CONN_MM.MENU{n0}.BINDEX(ceil(k1*numel(CONN_MM.MENU{n0}.BINDEX)));
                    %n1=max(1,min(CONN_MM.MENU{n0}.n, ceil((pos(1)-CONN_MM.MENU{n0}.position(1))./CONN_MM.MENU{n0}.position(3)*CONN_MM.MENU{n0}.n) ));
                else
                    k1=max(eps,min(1, 1-(pos(2)-CONN_MM.MENU{n0}.position(2))./CONN_MM.MENU{n0}.position(4) ));
                    n1=CONN_MM.MENU{n0}.BINDEX(ceil(k1*numel(CONN_MM.MENU{n0}.BINDEX)));
                    %n1=max(1,min(CONN_MM.MENU{n0}.n, CONN_MM.MENU{n0}.n+1-ceil((pos(2)-CONN_MM.MENU{n0}.position(2))./CONN_MM.MENU{n0}.position(4)*CONN_MM.MENU{n0}.n) ));
                end
                if ~CONN_MM.MENU{n0}.enable(n1), return; end
                if any(CONN_MM.INSIDE~=[n0,n1]),
                    if CONN_MM.INSIDE(1),
                        if ~iscell(CONN_MM.MENU{CONN_MM.INSIDE(1)}.callback{CONN_MM.INSIDE(2)}) && CONN_MM.MENU{CONN_MM.INSIDE(1)}.callback{CONN_MM.INSIDE(2)}~=n0, %leave an open menu
                            CONN_MM.MENU{CONN_MM.MENU{CONN_MM.INSIDE(1)}.callback{CONN_MM.INSIDE(2)}}.value=0;
                            conn_menumanager(CONN_MM.MENU{CONN_MM.INSIDE(1)}.callback{CONN_MM.INSIDE(2)},'off',1);
                        end
                        if CONN_MM.INSIDE(1)~=n0, %change menus
                            if ~iscell(CONN_MM.MENU{CONN_MM.INSIDE(1)}.callback{CONN_MM.INSIDE(2)}) && CONN_MM.MENU{CONN_MM.INSIDE(1)}.callback{CONN_MM.INSIDE(2)}==n0, % change to a forward-linked menu
                                if isempty(CONN_MM.GINSIDE), CONN_MM.GINSIDE=[CONN_MM.INSIDE(1),n0]; 
                                else CONN_MM.GINSIDE=[CONN_MM.GINSIDE,n0]; end
                            elseif ~isempty(CONN_MM.GINSIDE),
                                idx=find(CONN_MM.GINSIDE==n0); %change to a back-linked menu
                                if isempty(idx), idx=1; end; % change to a not-linked menu
                                for n2=length(CONN_MM.GINSIDE):-1:idx+1, conn_menumanager(CONN_MM.GINSIDE(n2),'off',1); CONN_MM.MENU{CONN_MM.GINSIDE(n2)}.value=0; end
                                if idx==1, CONN_MM.GINSIDE=[]; else CONN_MM.GINSIDE=CONN_MM.GINSIDE(1:idx); end
                            elseif CONN_MM.INSIDE(1), % change to a different menu
                                CONN_MM.MENU{CONN_MM.INSIDE(1)}.value=0;
                                conn_menumanager(CONN_MM.INSIDE(1),'off',1,'on',1);
                            end
                        end
                    end
                    
                    CONN_MM.HELP.string=CONN_MM.MENU{n0}.help{n1};
                    CONN_MM.MENU{n0}.value=n1;
                    conn_menumanager(n0,'off',1,'on',1);
                    if ~iscell(CONN_MM.MENU{n0}.callback{n1}),
                        conn_menumanager(CONN_MM.MENU{n0}.callback{n1},'on',1,'linkon',1);
                    end
                end
                CONN_MM.INSIDE=[n0,n1];
                if strcmp(handle,'cursorup'),
                    if CONN_MM.MENU{n0}.toggle==1, CONN_MM.MENU{n0}.state(:)=0; CONN_MM.MENU{n0}.state(n1)=1;
					elseif CONN_MM.MENU{n0}.toggle, CONN_MM.MENU{n0}.state(n1)=mod(CONN_MM.MENU{n0}.state(n1)+1,2); 
					else, CONN_MM.MENU{n0}.state(n1)=1; end
                    conn_menumanager(n0,'off',1,'on',1);
                    donecallback=false;
                    if  ~isempty(CONN_MM.MENU{n0}.callback{n1}),
                        if iscell(CONN_MM.MENU{n0}.callback{n1}) && CONN_MM.MENU{n0}.state(n1), 
                            conn_menumanager('cursorout');
                            CONN_MM.restpos=pos;
                            CONN_MM.resttime=now;
                            set(CONN_MM.gcf,'pointer',CONN_gui.waiticon);pause(.001);
                            if ~CONN_MM.MENU{n0}.toggle
                                CONN_MM.MENU{n0}.state(n1)=mod(CONN_MM.MENU{n0}.state(n1)+1,2);
                                if CONN_MM.ACTIVE(n0), conn_menumanager(n0,'off',1,'on',1); end
                            end
							feval(CONN_MM.MENU{n0}.callback{n1}{1},CONN_MM.MENU{n0}.callback{n1}{2:end}); 
                            set(CONN_MM.gcf,'pointer','arrow');
                            CONN_MM.restpos=[nan nan];
                            CONN_MM.resttime=nan;
                            donecallback=true;
                        elseif iscell(CONN_MM.MENU{n0}.callback2{n1})&&numel(CONN_MM.MENU{n0}.callback2)>=n1&&~isempty(CONN_MM.MENU{n0}.callback2{n1}) && CONN_MM.MENU{n0}.state(n1), 
                            conn_menumanager('cursorout');
                            CONN_MM.restpos=pos;
                            CONN_MM.resttime=now;
                            set(CONN_MM.gcf,'pointer',CONN_gui.waiticon);pause(.001);
                            if ~CONN_MM.MENU{n0}.toggle
                                CONN_MM.MENU{n0}.state(n1)=mod(CONN_MM.MENU{n0}.state(n1)+1,2);
                                if CONN_MM.ACTIVE(n0), conn_menumanager(n0,'off',1,'on',1); end
                            end
							feval(CONN_MM.MENU{n0}.callback2{n1}{1},CONN_MM.MENU{n0}.callback2{n1}{2:end}); 
                            set(CONN_MM.gcf,'pointer','arrow');
                            CONN_MM.restpos=[nan nan];
                            CONN_MM.resttime=nan;
                            donecallback=true;
						end
                    end
                    if ~CONN_MM.MENU{n0}.toggle&&~donecallback
                        CONN_MM.MENU{n0}.state(n1)=mod(CONN_MM.MENU{n0}.state(n1)+1,2);
                        if CONN_MM.ACTIVE(n0), conn_menumanager(n0,'off',1,'on',1); end
                    end
                    if donecallback, return; end
                end
            else % move outside active area
                if any(CONN_MM.INSIDE~=[0,0]),
                    CONN_MM.HELP.string='';
                    if ~iscell(CONN_MM.MENU{CONN_MM.INSIDE(1)}.callback{CONN_MM.INSIDE(2)}),
                        CONN_MM.MENU{CONN_MM.MENU{CONN_MM.INSIDE(1)}.callback{CONN_MM.INSIDE(2)}}.value=0;
                        conn_menumanager(CONN_MM.MENU{CONN_MM.INSIDE(1)}.callback{CONN_MM.INSIDE(2)},'off',1,'linkon',0);
                    end
                    if ~isempty(CONN_MM.GINSIDE),
                        for n2=length(CONN_MM.GINSIDE):-1:2, CONN_MM.MENU{CONN_MM.GINSIDE(n2)}.value=0; conn_menumanager(CONN_MM.GINSIDE(n2),'off',1); end
                        CONN_MM.MENU{CONN_MM.GINSIDE(1)}.value=0; conn_menumanager(CONN_MM.GINSIDE(1),'off',1,'on',1);
                        CONN_MM.GINSIDE=[];
                    else % change to a different menu
                        CONN_MM.MENU{CONN_MM.INSIDE(1)}.value=0;
                        conn_menumanager(CONN_MM.INSIDE(1),'off',1,'on',1);
                    end
                    CONN_MM.INSIDE=[0,0];
                end
            end
            if ~any(CONN_MM.INSIDE)&&~isempty(CONN_MM.onregionarea)
                xpos=repmat(pos,[size(CONN_MM.onregionarea,1),1]);
                md=find(min(min(xpos-CONN_MM.onregionarea(:,1:2),CONN_MM.onregionarea(:,1:2)+CONN_MM.onregionarea(:,3:4)-xpos),[],2)>-.001)';
                onregions=[];
                offregions=[];
                colorB=[];
                colorC=[];
                if ~strcmp(handle,'cursormove'),
                    for n0=reshape(md,1,[])
                        if CONN_MM.onregionvisible(n0)==2&&~isempty(CONN_MM.onregioncallback{n0})
                            if ischar(CONN_MM.onregioncallback{n0}), eval(CONN_MM.onregioncallback{n0});
                            elseif isa(CONN_MM.onregioncallback{n0},'function_handle'), feval(CONN_MM.onregioncallback{n0});
                            end
                        end
                    end
                    for n0=reshape(setdiff(1:numel(CONN_MM.onregioncallback),[md,CONN_MM.onregioninside]),1,[]),
                        if all(ishandle(CONN_MM.onregioncallback{n0}))
                            if CONN_MM.onregionvisible(n0)==1&&all(strcmp(get(CONN_MM.onregioncallback{n0},'visible'),'on')), offregions=[offregions CONN_MM.onregionhandle{n0}];
                            elseif CONN_MM.onregionvisible(n0)==-1&&~all(strcmp(get(CONN_MM.onregioncallback{n0},'visible'),'on')), onregions=[onregions CONN_MM.onregionhandle{n0}];
                            end
                        end
                    end
                end
                %set(CONN_MM.onregionhandle(CONN_MM.onregionvisible>0),'visible','off');
                %set(CONN_MM.onregionhandle(CONN_MM.onregionvisible<0),'visible','on');

                if ~isequal(md,CONN_MM.onregioninside)
                    CONN_MM.private.time=clock;
                    match=conn_bsxfun(@eq,CONN_MM.onregioninside',md);
                    for n0=CONN_MM.onregioninside(~any(match,2)'),%reshape(setdiff(CONN_MM.onregioninside,md),1,[]) % going out of these regions
                        if n0>0
                            if ~isempty(CONN_MM.onregioncallback{n0})&&CONN_MM.onregionvisible(n0)~=2,
                                if ishandle(CONN_MM.onregioncallback{n0}), cond=all(strcmp(get(CONN_MM.onregioncallback{n0},'visible'),'on')); %if ~cond, offregions=[offregions CONN_MM.onregionhandle{n0}]; end
                                else cond=all(feval(CONN_MM.onregioncallback{n0},CONN_MM.onregionhandle{n0}));
                                end
                            else cond=1; end
                            if cond&&all(ishandle(CONN_MM.onregionhandle{n0})),
                                if CONN_MM.onregionvisible(n0)==1, offregions=[offregions CONN_MM.onregionhandle{n0}]; %set(CONN_MM.onregionhandle{n0},'visible','off');
                                elseif CONN_MM.onregionvisible(n0)==-1, onregions=[onregions CONN_MM.onregionhandle{n0}]; %set(CONN_MM.onregionhandle{n0},'visible','on');
                                elseif CONN_MM.onregionvisible(n0)==0 colorB=[colorB CONN_MM.onregionhandle{n0}]; 
                                end
                            end
                        end
                    end
%                 elseif etime(clock,CONN_MM.private.time)>.025
%                     for n0=md,%setdiff(md,CONN_MM.onregioninside) %in
                    for n0=md(~any(match,1)),%reshape(setdiff(md,CONN_MM.onregioninside),1,[]) % going into these regions
                        if ~isempty(CONN_MM.onregioncallback{n0})&&CONN_MM.onregionvisible(n0)~=2,
                            if ishandle(CONN_MM.onregioncallback{n0}), cond=all(strcmp(get(CONN_MM.onregioncallback{n0},'visible'),'on')); %if ~cond, offregions=[offregions CONN_MM.onregionhandle{n0}]; end
                            else cond=all(feval(CONN_MM.onregioncallback{n0},CONN_MM.onregionhandle{n0}));
                            end
                        else cond=1; end
                        if cond&&all(ishandle(CONN_MM.onregionhandle{n0})),
                            if CONN_MM.onregionvisible(n0)==1, onregions=[onregions CONN_MM.onregionhandle{n0}]; %set(CONN_MM.onregionhandle{n0},'visible','on'); % make visible when over this region
                            elseif CONN_MM.onregionvisible(n0)==-1, offregions=[offregions CONN_MM.onregionhandle{n0}]; %set(CONN_MM.onregionhandle{n0},'visible','off'); % make invisible when over this region
                            elseif CONN_MM.onregionvisible(n0)==0, colorC=[colorC CONN_MM.onregionhandle{n0}]; % change color when over this region
                            end
                        end
                    end
                end
                if ~isempty(onregions), set(onregions,'visible','on'); end
                if ~isempty(offregions), set(offregions,'visible','off'); end
                if CONN_gui.domacGUIbugfix==2 
                    if ~isempty(colorB), colorB=colorB(arrayfun(@(x)~strcmp(get(x,'style'),'popupmenu'),colorB)); end
                    if ~isempty(colorC), colorC=colorC(arrayfun(@(x)~strcmp(get(x,'style'),'popupmenu'),colorC)); end
                end
                if 0,%~CONN_gui.isjava
                    if ~isempty(colorB), colorB=colorB(arrayfun(@(x)~strcmp(get(x,'style'),'pushbutton'),colorB)); end
                    if ~isempty(colorC), colorC=colorC(arrayfun(@(x)~strcmp(get(x,'style'),'pushbutton'),colorC)); end
                end
                %if ~isempty(colorB), set(colorB,'foregroundcolor',[.4 .4 .4]+.2*(mean(CONN_gui.backgroundcolor)<.5)); end
                %if ~isempty(colorC), set(colorC,'foregroundcolor',[0 0 0]+1*(mean(CONN_gui.backgroundcolor)<.5)); end
                %if ~isempty(colorB), for n1=1:numel(colorB), set(colorB(n1),'foregroundcolor',[.4 .4 .4]+.2*(mean(get(colorB(n1),'backgroundcolor'))<.5)); % color when going out
                %end; end
                %if ~isempty(colorC), for n1=1:numel(colorC), set(colorC(n1),'foregroundcolor',[0 0 0]+1*(mean(get(colorC(n1),'backgroundcolor'))<.5)); % color when going in
                %end; end
                deltacolor=.025;%*(1-2*(CONN_gui.backgroundcolor>.5));
                if ~isempty(colorB), % color when going out
                    for n1=1:numel(colorB),
                        tbg=get(colorB(n1),'backgroundcolor');
                        tal=0;%~CONN_gui.isjava&isequal(get(colorB(n1),'style'),'popupmenu');
                        if all(tbg==max(0,min(1,CONN_gui.backgroundcolor+deltacolor))), set(colorB(n1),'backgroundcolor',CONN_gui.backgroundcolor,'foregroundcolor',CONN_gui.fontcolor);
                        elseif all(tbg==max(0,min(1,CONN_gui.backgroundcolorA+deltacolor))), set(colorB(n1),'backgroundcolor',CONN_gui.backgroundcolorA,'foregroundcolor',CONN_gui.fontcolorA);
                        elseif all(tbg==max(0,min(1,CONN_gui.backgroundcolorE+deltacolor))), set(colorB(n1),'backgroundcolor',CONN_gui.backgroundcolorE,'foregroundcolor',CONN_gui.fontcolorA);
                        end 
                    end
                end
                if ~isempty(colorC), % color when going in
                    for n1=1:numel(colorC),
                        tbg=get(colorC(n1),'backgroundcolor');
                        tal=0;%~CONN_gui.isjava&isequal(get(colorC(n1),'style'),'popupmenu');
                        if all(tbg==CONN_gui.backgroundcolor), set(colorC(n1),'backgroundcolor',max(0,min(1,CONN_gui.backgroundcolor+deltacolor)),'foregroundcolor',tal*[.5 .5 .5]+(1-tal)*round(CONN_gui.fontcolor));
                        elseif all(tbg==CONN_gui.backgroundcolorA), set(colorC(n1),'backgroundcolor',max(0,min(1,CONN_gui.backgroundcolorA+deltacolor)),'foregroundcolor',tal*[.5 .5 .5]+(1-tal)*round(CONN_gui.fontcolorA));
                        elseif all(tbg==CONN_gui.backgroundcolorE), set(colorC(n1),'backgroundcolor',max(0,min(1,CONN_gui.backgroundcolorE+deltacolor)),'foregroundcolor',tal*[.5 .5 .5]+(1-tal)*round(CONN_gui.fontcolorA));
                        end 
                    end
                end
                    %if ~isempty(colorB), for n1=1:numel(colorB), set(colorB(n1),'foregroundcolor',[.0 .0 .0]+1*(mean(get(colorB(n1),'backgroundcolor'))>=.5)); % color when going out
                    %end; end
                    %if ~isempty(colorC), for n1=1:numel(colorC), set(colorC(n1),'foregroundcolor',[0 0 0]+1*(mean(get(colorC(n1),'backgroundcolor'))<.5)); % color when going in
                    %end; end
                if CONN_gui.doemphasis3
                    if ~isempty(colorB), colorB=colorB(arrayfun(@(x)isequal(get(x,'backgroundcolor'),max(0,min(1,.9*CONN_gui.backgroundcolorA+.1*[.5 .5 .5]))),colorB)); end
                    if ~isempty(colorC), colorC=colorC(arrayfun(@(x)isequal(get(x,'backgroundcolor'),CONN_gui.backgroundcolorA),colorC)); end
                    if ~isempty(colorB), set(colorB,'backgroundcolor',CONN_gui.backgroundcolorA); end
                    if ~isempty(colorC), set(colorC,'backgroundcolor',max(0,min(1,.9*CONN_gui.backgroundcolorA+.1*[.5 .5 .5]))); end
                end
                CONN_MM.onregioninside=md;
                if isempty(CONN_MM.onregioninside), CONN_MM.onregioninside=0; end
                anymtn=false;
                for n0=md(:)' % mouse-movement inside region
                    if numel(CONN_MM.onregionmotioncallback)>=n0&&~isempty(CONN_MM.onregionmotioncallback{n0}), 
                        if ~isempty(CONN_MM.onregioncallback{n0}),
                            if ishandle(CONN_MM.onregioncallback{n0}), cond=all(strcmp(get(CONN_MM.onregioncallback{n0},'visible'),'on')); %if ~cond, offregions=[offregions CONN_MM.onregionhandle{n0}]; end
                            else
                                %disp(CONN_MM.onregioncallback{n0})
                                cond=all(feval(CONN_MM.onregioncallback{n0},CONN_MM.onregionhandle{n0}));
                            end
                        else cond=1; end
                        if cond, anymtn=true; conn_menumanager_mousemove(CONN_MM.onregionmotioncallback{n0},handle); end
                    end
                end
%                 if anymtn, set(CONN_MM.gcf,'pointer','crosshair');
%                 else       set(CONN_MM.gcf,'pointer','arrow');
%                 end
            end
            CONN_MM.private.busy=0;
        case 'ison',
            handle=false;
            try, handle=ishandle(CONN_MM.gcf); end
		case 'on',
            CONN_MM=struct('MENU',[],'ACTIVE',zeros(1,0),'INSIDE',[0,0],'GINSIDE',[],'POSITION',zeros(0,4),'BINDEX',1);
            CONN_MM.HELP=struct('handle',[],'string','');
            CONN_MM.gcf=gcf;
            CONN_MM.map=get(gcf,'colormap');
            set([0,CONN_MM.gcf],'units','pixels');
            CONN_MM.private=struct('busy',0,'time',clock);
            CONN_MM.onregionarea=[];
            CONN_MM.onregionhandle={};
            CONN_MM.onregioncallback={};
            CONN_MM.onregionmotioncallback={};
            CONN_MM.onregionvisible=[];
            CONN_MM.onregioninside=0;
            CONN_MM.restpos=[nan nan];
            CONN_MM.resttime=nan;
            set(CONN_MM.gcf,'windowbuttonmotionfcn','conn_menumanager((''cursormove''));','windowbuttonupfcn','conn_menumanager(''cursorup'');','windowbuttondownfcn','conn_menumanager(''cursordown'');');
		case 'off',
            for n0=find(CONN_MM.ACTIVE), conn_menumanager(n0,'off',1); end
        	set(CONN_MM.gcf,'windowbuttonmotionfcn','','windowbuttonupfcn','');
 		case 'clf',
            CONN_MM.private.busy=1;
            %figure(CONN_MM.gcf);
            if ~isfield(CONN_gui,'isresizing')||~CONN_gui.isresizing
                if 0
                    hax=axes('units','norm','position',[.91,.92,.07,.015],'parent',CONN_MM.gcf);
                    text(0,0,'please wait...','horizontalalignment','center','color',CONN_gui.fontcolorA,'parent',hax); set(hax,'xlim',[-1 1],'ylim',[-1 1],'visible','off');
                    [nill,hc]=conn_menu_plotmatrix('',CONN_MM.gcf,[2 1 8],[.91 .93 .07 .02]);
                    delete(hc(ishandle(hc)));
                else
                    hax=axes('units','norm','position',[.78,.93,.20,.02],'parent',CONN_MM.gcf);
                    text(1,0,sprintf('GUI busy\\fontsize{%d} (perhaps%s)\\fontsize{%d} please wait',5+CONN_gui.font_offset,regexprep(conn_msg(1),', please wait.*$',''),6+CONN_gui.font_offset),'horizontalalignment','right','interpreter','tex','color',CONN_gui.fontcolorA, 'fontsize',6+CONN_gui.font_offset,'fontangle','normal','parent',hax); set(hax,'xlim',[-1 1],'ylim',[-1 1],'visible','off');
                    %hax=axes('units','norm','position',[.4,.05,.20,.025],'parent',CONN_MM.gcf);
                    %text(0,0,['GUI busy',conn_msg(1)],'horizontalalignment','center','color',CONN_gui.fontcolorA, 'fontsize',8+CONN_gui.font_offset,'fontangle','normal','parent',hax); set(hax,'xlim',[-1 1],'ylim',[-1 1],'visible','off');
                    drawnow;
                end
            end
 			clf(CONN_MM.gcf);
            for n1=1:length(CONN_MM.MENU), CONN_MM.MENU{n1}.handle.axes={}; CONN_MM.MENU{n1}.handle.prevstate=nan; CONN_MM.MENU{n1}.value=0; end
 			if ~isempty(CONN_MM.ACTIVE), CONN_MM.ACTIVE(:)=0; end
            CONN_MM.INSIDE=[0,0]; CONN_MM.GINSIDE=[];
            CONN_MM.onregionarea=[];
            CONN_MM.onregionhandle={};
            CONN_MM.onregioncallback={};
            CONN_MM.onregionmotioncallback={};
            CONN_MM.onregionvisible=[];
            CONN_MM.onregioninside=0;
            CONN_MM.private.busy=0;
        case 'updatebackgroundcolor'
            color=get(CONN_MM.gcf,'color');
            for handle=1:numel(CONN_MM.MENU)
                CONN_MM.MENU{handle}.backgroundcolor=color;
                CONN_MM.CDATA{handle}=conn_menumanager_buttonshape(CONN_MM.MENU{handle},CONN_MM.gcf);
            end
        case 'onregion',
            if ~isempty(varargin{1}),%&&~any(ismember(varargin{1},[CONN_MM.onregionhandle{:}]))
                treg=varargin{1}; if ~iscell(treg), treg={treg}; end
                if numel(varargin)>1&&~isempty(varargin{2}), tvis=varargin{2}; else tvis=1; end
                if numel(varargin)>2&&~isempty(varargin{3}), tpos=varargin{3}; else tpos=get(varargin{1},'position'); end
                if numel(varargin)>3&&~isempty(varargin{4}), tcall=varargin{4}; else tcall=[]; end
                if numel(varargin)>4&&~isempty(varargin{5}), tmcall=varargin{5}; else tmcall=[]; end
                CONN_MM.onregionarea=cat(1,CONN_MM.onregionarea,tpos);
                CONN_MM.onregionhandle=[CONN_MM.onregionhandle,treg];
                CONN_MM.onregioncallback=[CONN_MM.onregioncallback,{tcall}];
                CONN_MM.onregionmotioncallback=[CONN_MM.onregionmotioncallback,{tmcall}];
                CONN_MM.onregionvisible=[CONN_MM.onregionvisible,tvis];
            end
        case 'onregionremove',
            if numel(varargin)<1||isempty(varargin{1}), return; end
            for n1=1:numel(varargin{1})
                idx=find(cellfun(@(x)any(x==varargin{1}(n1)),CONN_MM.onregionhandle,'ErrorHandler',@(varargin)false));
                if ~isempty(idx)
                    CONN_MM.onregionarea(idx,:)=[];
                    CONN_MM.onregionhandle(idx)=[];
                    CONN_MM.onregioncallback(idx)=[];
                    CONN_MM.onregionmotioncallback(idx)=[];
                    CONN_MM.onregionvisible(idx)=[];
                    CONN_MM.onregioninside(ismember(CONN_MM.onregioninside,idx))=[];
                    CONN_MM.onregioninside=CONN_MM.onregioninside(:)'-sum(bsxfun(@gt,CONN_MM.onregioninside(:)',idx(:)),1);
                end
            end
        case 'helpstring'
            bg=get(CONN_MM.gcf,'color');
            if ishandle(CONN_MM.HELP.handle), 
                set(CONN_MM.HELP.handle,'string',[varargin{:}],'visible','on','color',mod(bg-.3,1));
                %set(CONN_MM.HELP.handle,'string',[varargin{:}],'visible','on','foregroundcolor',1-bg);
            else
                ha=axes('units','norm','position',[.15,.02,.7,.03],'visible','off','parent',CONN_MM.gcf);
                CONN_MM.HELP.handle=text(0,1,[varargin{:}],'color',mod(bg-.3,1),'fontname','default','fontsize',8+CONN_gui.font_offset,'horizontalalignment','center','verticalalignment','middle','interpreter','none','parent',ha);
                set(ha,'xlim',[-1 1],'ylim',[0 2],'visible','off');
                %CONN_MM.HELP.handle=uicontrol('style','text','units','norm','position',[.15,.02,.7,.03],'backgroundcolor',bg,'foregroundcolor',1-bg,'fontname','default','string',[varargin{:}],'fontsize',8+CONN_gui.font_offset);
            end
            
        case 'private',
            handle=CONN_MM;
        otherwise
            error(['Incorrect parameter ',handle]);
    end
else
    if isempty(handle), % creates new menu entry
        if isempty(CONN_MM), CONN_MM=struct('MENU',[],'ACTIVE',zeros(1,0),'INSIDE',[0,0],'GINSIDE',[],'POSITION',zeros(0,4),'BINDEX',1); end
        handle=length(CONN_MM.MENU)+1; 
        CONN_MM.POSITION(handle,:)=0;
        CONN_MM.ACTIVE(handle)=0;
        fntname=CONN_gui.fontname;
        % note1: possible status are: [Mouse outside of item, Mouse inside item, Item selected and mouse inside item, Item selected and mouse outside item]
        params=struct(...
            'n',0,...                                           % number of items
            'position',[0,0,.1,.1],...						    % 1 x 4 normalized-unit position vector
            'dposition',[0,0],...                               % 1 x 2 normalized-unit position change vector (item changes positions when selected)
            'dfont',[0],...                                     % font size change vector (item changes positions when selected)
            'string',{{}},...                                   % cell array with one text string per item
            'order','vertical',...                              % 'vertical'/'horizontal' order of items
            'callback',{{}},...									% cell array specifying the callback functions (one per item), callback functions can be defined as cell arrays (first element is callback function, next elements are optional callback arguments; note: conn_menumanager automatically will add the "state" status as a last argument to this call), or as other conn_menumanager handles
            'callback2',{{}},...                                % cell array specifying the alternative callback functions (one per item), callback functions can be defined as cell arrays (first element is callback function, next elements are optional callback arguments; note: conn_menumanager automatically will add the "state" status as a last argument to this call), or as other conn_menumanager handles
            'help',{{}},...										% cell array with one help text string per item
            'rotate',0,...										% 1/0 to indicate whether to rotate 90 degrees the button
            'backgroundcolor',[],...                            % 1x3 vector of color components
            'color',[0.2083 0.2125 0.2167;0.2083 0.2125 0.2167;0.08333 0.3333 0.5417],...       % 3x3 matrix of color components. Rows are status (cursor-out,cursor-in,selected) and columns are rgb components
            'colorb',[0.9167 0.6667 0.2917;0.9167 0.6667 0.2917;0.9167 0.6667 0.2917],...
            'horizontalalignment','center',...
            'transparent',1,...
            'fontcolor',[.75,.75,.75; 1,1,1;1,1,.9],...                             % 
            'fontname',fntname,...                              % 
            'fontsize',8,...									% 
            'fontweight','normal',...                             % 
            'fontangle','normal',...                             % 
            'roll',0,...
            'bordertype','round',...                          % 
            'state',[],...                                     %  array of 1/0 status values (one value per item)
            'linkon',0,...
            'enable',[],...
            'displayed',1,...
            'prevstate',[],...
            'toggle',2,...										% 0/1 to indicate whether toggle button (2 indicates possible multiple buttons selected at a time)
            'on',0,...                                          % Set to 1 to display the object
            'off',0,...                                         % Set to 1 to hide the object
            'value',0,...                                       %  private (what item is mouse on)
            'handle',struct('axes',[],'image',[],'text',[]));   %  private (handles to draw objects)
        for n1=1:2:nargin-2, params.(lower(varargin{n1}))=varargin{n1+1}; end
        if isempty(params.string), params.string=cell(1,params.n); end
        if isempty(params.callback), for n1=1:params.n, params.callback{n1}={}; end; end
        if isempty(params.callback2), for n1=1:params.n, params.callback2{n1}={}; end; end
        if isempty(params.help), params.help=cell(1,params.n); end
        if isempty(params.state), params.state=zeros(1,params.n); end
        if isempty(params.enable), params.enable=ones(1,params.n); end
        if isempty(params.prevstate), params.prevstate=nan(1,params.n); end
        if strcmp(params.bordertype,'round')|strcmp(params.bordertype,'xround'), 
            if isfield(params,'colorb'), 
                params.color=params.colorb;
            else params.color=max(0,min(1, params.color(:,[3 3 1])*diag([2 2 2]) ));
            end
            %%params.fontcolor(2,:)=params.fontcolor(3,:);%(:,[3 2 1]); 
        else
            %params.color=params.color(:,[3 2 1]); 
        end
        if params.displayed
            if isempty(params.backgroundcolor), params.backgroundcolor=get(CONN_MM.gcf,'color'); end
            CONN_MM.CDATA{handle}=conn_menumanager_buttonshape(params,CONN_MM.gcf);
        end
        if isempty(CONN_MM.MENU), CONN_MM.MENU={params}; handle=1;
        else CONN_MM.MENU{handle}=params; end
    elseif all(handle>0)
        for nhandle=1:length(handle),
            for n1=1:2:nargin-2, CONN_MM.MENU{handle(nhandle)}.(lower(varargin{n1}))=varargin{n1+1}; end
        end
        if ~mod(nargin,2),
            handle=CONN_MM.MENU{handle(1)}.(lower(varargin{end}));
            return;
        end
    else handle=1;
    end 

    GCF=[];
    KEEPHANDLES=true;
    redrawnow=false;
    turnon=[];
    for nhandle=1:length(handle),
        thishandle=handle(nhandle);
        if ~CONN_MM.MENU{thishandle}.displayed, continue; end
        if CONN_MM.MENU{thishandle}.off || CONN_MM.MENU{thishandle}.on
            %disp([CONN_MM.MENU{thishandle}.off,CONN_MM.MENU{thishandle}.on,thishandle]);
            %disp(['off ',num2str(thishandle)]);
            if ~CONN_MM.MENU{thishandle}.off||~CONN_MM.MENU{thishandle}.on, changed=1:CONN_MM.MENU{thishandle}.n;
            else
                changed=CONN_MM.MENU{thishandle}.state~=CONN_MM.MENU{thishandle}.prevstate;
                CONN_MM.MENU{thishandle}.prevstate=CONN_MM.MENU{thishandle}.state;
                if CONN_MM.MENU{thishandle}.value>0&&CONN_MM.MENU{thishandle}.value<=CONN_MM.MENU{thishandle}.n, changed(CONN_MM.MENU{thishandle}.value)=1; CONN_MM.MENU{thishandle}.prevstate(CONN_MM.MENU{thishandle}.value)=nan; end
                changed=find(changed);
            end
            for n1=changed(:)',%1:length(CONN_MM.MENU{thishandle}.handle.axes),
				if numel(CONN_MM.MENU{thishandle}.handle.axes)>=n1&&any(ishandle(CONN_MM.MENU{thishandle}.handle.axes{n1})),
					%delete(CONN_MM.MENU{thishandle}.handle.text{n1}(ishandle(CONN_MM.MENU{thishandle}.handle.text{n1})));
					%if ~CONN_gui.dounixGUIbugfix, delete(CONN_MM.MENU{thishandle}.handle.image(n1)); end
                    if KEEPHANDLES, set(findobj(CONN_MM.MENU{thishandle}.handle.axes{n1}(ishandle(CONN_MM.MENU{thishandle}.handle.axes{n1}))),'visible','off'); 
                    else delete(findobj(CONN_MM.MENU{thishandle}.handle.axes{n1}(ishandle(CONN_MM.MENU{thishandle}.handle.axes{n1})))); 
                    end
				end
            end
            CONN_MM.ACTIVE(thishandle)=0;
            CONN_MM.MENU{thishandle}.off=0;
        end
        if CONN_MM.MENU{thishandle}.on,
            %disp(['on ',num2str(thishandle)]);
            %if ~isequal(CONN_MM.gcf,gcf), GCF=gcf; figure(CONN_MM.gcf); end
            t1=cellfun(@iscell,CONN_MM.MENU{thishandle}.string);t2=cellfun('length',CONN_MM.MENU{thishandle}.string);t0=1+t1.*(t2-1); t3=(0:CONN_MM.MENU{thishandle}.n)+cumsum([0,t0-1]);t4=t3/t3(end);
            CONN_MM.MENU{thishandle}.BINDEX=zeros(1,t3(end));
            for n1=1:CONN_MM.MENU{thishandle}.n,
                CONN_MM.MENU{thishandle}.BINDEX(t3(n1)+1:t3(n1+1))=n1;
            end
            for n1=changed(:)',%1:CONN_MM.MENU{thishandle}.n,
                x=1+(CONN_MM.MENU{thishandle}.value==n1 & ~CONN_MM.MENU{thishandle}.state(n1))+2*CONN_MM.MENU{thishandle}.state(n1);
                ximage=CONN_MM.CDATA{thishandle}{n1}{x};
                if ~CONN_gui.doemphasis3&&CONN_MM.MENU{thishandle}.linkon, 
                    ximage=max(0,min(1,ximage*1)); ximage=ximage.*repmat(.75+.20*tanh((size(ximage,1):-1:1)/2)'*tanh(min(0:size(ximage,2)-1,size(ximage,2)-1:-1:0)/4),[1,1,size(ximage,3)]); 
                    if 0,%CONN_MM.MENU{thishandle}.order(1)~='h'
                        ximage0=.5+0*round(ximage);
                        if t4(n1)==0, ximagei1=ceil(min(size(ximage,1)*.05,8));
                        else ximagei1=1;
                        end
                        if t4(n1)==1, ximagei3=size(ximage,1)+1-ceil(min(size(ximage,1)*.15,8));
                        else ximagei3=size(ximage,1);
                        end
                        ximagei2=1;
                        ximagei4=size(ximage,2)+1-ximagei2;
                        ximage(ximagei1:ximagei3,ximagei2,:)=ximage0(ximagei1:ximagei3,ximagei2,:); % left
                        ximage(ximagei1:ximagei3,ximagei4,:)=ximage0(ximagei1:ximagei3,ximagei4,:); % right
                        if t4(n1)==0, % top
                            ximage(ximagei1,ximagei2:ximagei4,:)=ximage0(ximagei1,ximagei2:ximagei4,:);
                        end
                        if t4(n1+1)==1 % bottom
                            ximage(ximagei3,ximagei2:ximagei4,:)=ximage0(ximagei3,ximagei2:ximagei4,:);
                        end
                    end
                end
                if ~iscell(CONN_MM.MENU{thishandle}.callback{n1})&&CONN_MM.MENU{thishandle}.order(1)~='h'
                    try, ximage(round(size(ximage,1)/2)+(1:5),end-16:end-8,:)=.5; end
                end
                %CONN_MM.MENU{thishandle}.BINDEX(t3(n1)+1:t3(n1+1))=n1;
                if CONN_MM.MENU{thishandle}.order(1)=='h', pos=[CONN_MM.MENU{thishandle}.position(1)+t4(n1)*CONN_MM.MENU{thishandle}.position(3), CONN_MM.MENU{thishandle}.position(2), CONN_MM.MENU{thishandle}.position(3)*(t4(n1+1)-t4(n1)), CONN_MM.MENU{thishandle}.position(4)];
                else pos=[CONN_MM.MENU{thishandle}.position(1), CONN_MM.MENU{thishandle}.position(2)+CONN_MM.MENU{thishandle}.position(4)-t4(n1+1)*CONN_MM.MENU{thishandle}.position(4), CONN_MM.MENU{thishandle}.position(3), CONN_MM.MENU{thishandle}.position(4)*(t4(n1+1)-t4(n1))]; end
                if CONN_MM.MENU{thishandle}.state(n1), pos(1:2)=pos(1:2)+CONN_MM.MENU{thishandle}.dposition; end
                if length(CONN_MM.MENU{thishandle}.handle.axes)>=n1 && any(ishandle(CONN_MM.MENU{thishandle}.handle.axes{n1})), 
                    if KEEPHANDLES, 
                        set(findobj(CONN_MM.MENU{thishandle}.handle.axes{n1}(ishandle(CONN_MM.MENU{thishandle}.handle.axes{n1}))),'visible','off'); 
                    else delete(findobj(CONN_MM.MENU{thishandle}.handle.axes{n1}(ishandle(CONN_MM.MENU{thishandle}.handle.axes{n1})))); 
                    end
                end
                fontsize=CONN_MM.MENU{thishandle}.fontsize; 
                %if px1>px0&&px1>px2, fontcolor=fontcolor1; 
                %elseif px2>px0, fontcolor=fontcolor2; 
                %end
%                 if (mx1<.4&&mx2<.4)||(mx1>.6&&mx2>.6), fontcolor=1-fontcolor; 
%                 elseif mx1>.3&&mx1<.7&&mx2>.3&&mx2<.7, fontcolor=max(0,min(1,fontcolor+.5*(2*((mx2>mx1)-1)))); 
%                 end
                fontweight=CONN_MM.MENU{thishandle}.fontweight;
                if x==3, 
                    fontsize=fontsize+CONN_MM.MENU{thishandle}.dfont; 
                    %fontweight='bold';
                %elseif x==2, fontsize=fontsize+round(CONN_MM.MENU{thishandle}.dfont/2);
                end
                
                if 1,%CONN_gui.dounixGUIbugfix,
                    %CONN_MM.MENU{thishandle}.handle.text{n1}=[];
                    sstrings=cellstr(CONN_MM.MENU{thishandle}.string{n1});
                    nstrings=numel(sstrings);
                    for n2=1:nstrings
                        if CONN_MM.MENU{thishandle}.order(1)=='h', tpos=[pos(1)+(n2-1)/nstrings*pos(3)-1e-3 pos(2) pos(3)/nstrings+2e-3 pos(4)];
                        else tpos=[pos(1) pos(2)+pos(4)-n2/nstrings*pos(4)-1e-3 pos(3) pos(4)/nstrings+2e-3];
                        end
                        %htemp=uicontrol('style','text','units','norm','position',tpos,'string',sstrings{n2},'buttondownfcn','conn_menumanager(''cursorup'');','enable','inactive',...
                        if CONN_gui.isjava % draws menus as uicontrol objects in R2024 or below
                            if numel(CONN_MM.MENU{thishandle}.handle.axes)>=n1&&numel(CONN_MM.MENU{thishandle}.handle.axes{n1})>=n2&&ishandle(CONN_MM.MENU{thishandle}.handle.axes{n1}(n2)), 
                                htemp=CONN_MM.MENU{thishandle}.handle.axes{n1}(n2); 
                            else
                                htemp=uicontrol('style','togglebutton','visible','off','units','norm','position',tpos,'string',sstrings{n2},...
                                    'backgroundcolor',CONN_gui.backgroundcolor,...
                                    'fontname',CONN_MM.MENU{thishandle}.fontname,...
                                    'fontsize',fontsize+CONN_gui.font_offset,...
                                    'fontunits','points',...
                                    'horizontalalignment',CONN_MM.MENU{thishandle}.horizontalalignment,...
                                    'fontweight',fontweight,...
                                    'fontangle',CONN_MM.MENU{thishandle}.fontangle,...
                                    'parent',CONN_MM.gcf,'callback','conn_menumanager(''cursorup'');');
                            end
                            if UIXCONTROLTOOLTIP&&CONN_gui.tooltips,set(htemp,'tooltipstring',CONN_MM.MENU{thishandle}.help{n1}); end
                        else 
                            if numel(CONN_MM.MENU{thishandle}.handle.axes)>=n1&&numel(CONN_MM.MENU{thishandle}.handle.axes{n1})>=n2&&ishandle(CONN_MM.MENU{thishandle}.handle.axes{n1}(n2)), 
                                htemp=CONN_MM.MENU{thishandle}.handle.axes{n1}(n2); 
                            else
                                if CONN_MM.MENU{thishandle}.linkon % draws floating menus as uipanel objects in R2025a and above (avoids flickering effect when using uicontrols, and "behind-uicontrols" effect when using axes)
                                    htemp=uipanel('units','norm','position',tpos,'backgroundcolor',CONN_gui.backgroundcolor,'foregroundcolor',CONN_gui.backgroundcolor,'bordercolor',CONN_gui.backgroundcolor,'visible','off','parent',CONN_MM.gcf);
                                    htemp1=axes('units','norm','position',[0 0 1 1],'xtick',[],'ytick',[],'color',CONN_gui.backgroundcolor,'xcolor',CONN_gui.backgroundcolor,'ycolor',CONN_gui.backgroundcolor,'visible','off','parent',htemp);
                                else % draws fixed menus as axes objects in R2025a and above (avoids flickering effect when using uicontrols, and off/on effect when using uipanels)
                                    htemp1=axes('units','norm','position',tpos,'xtick',[],'ytick',[],'color',CONN_gui.backgroundcolor,'xcolor',CONN_gui.backgroundcolor,'ycolor',CONN_gui.backgroundcolor,'visible','off','parent',CONN_MM.gcf);
                                    htemp=htemp1;
                                end
                                htemp2=image(0,'parent',htemp1,'visible','off');
                                htemp3=text(0,0,sstrings{n2},'clipping','on','parent',htemp1,'visible','off',...
                                    'interpreter','none',...
                                    'fontname',CONN_MM.MENU{thishandle}.fontname,...
                                    'fontsize',fontsize+CONN_gui.font_offset,...
                                    'fontunits','points',...
                                    'horizontalalignment','center',...
                                    'fontweight',CONN_MM.MENU{thishandle}.fontweight,...
                                    'fontangle',CONN_MM.MENU{thishandle}.fontangle);
                                if CONN_MM.MENU{thishandle}.linkon, set(htemp,'userdata',[htemp1 htemp2 htemp3]);
                                else set(htemp,'userdata',[htemp2 htemp3]);
                                end
                            end
                        end
                        xextent=round((get(0,'screensize')*[0 0 0 0;0 0 0 0;1 0 1 0;0 1 0 1]).*tpos + 4);
                        %set(htemp(1),'units','pixels');
                        %xextent=round(get(htemp(1),'position'))+4;
%                         if CONN_MM.MENU{thishandle}.order(1)=='h', tximage=ximage(round(linspace(1,size(ximage,1),xextent(4)*3)),round(linspace(1+(n2-1)/nstrings*size(ximage,2),n2/nstrings*size(ximage,2),xextent(3)*3)),:);
%                         else tximage=ximage(round(linspace(1+(n2-1)/nstrings*size(ximage,1),n2/nstrings*size(ximage,1),xextent(4)*3)),round(linspace(1,size(ximage,2),xextent(3)*3)),:);
%                         end
%                         tximage=permute(mean(mean(reshape(tximage,3,size(tximage,1)/3,3,size(tximage,2)/3,size(tximage,3)),1),3),[2,4,5,1,3]);
                        if CONN_MM.MENU{thishandle}.order(1)=='h', tximage=ximage(round(linspace(1,size(ximage,1),xextent(4)*1)),round(linspace(1+(n2-1)/nstrings*size(ximage,2),n2/nstrings*size(ximage,2),xextent(3)*1)),:);
                        else tximage=ximage(round(linspace(1+(n2-1)/nstrings*size(ximage,1),n2/nstrings*size(ximage,1),xextent(4)*1)),round(linspace(1,size(ximage,2),xextent(3)*1)),:);
                        end
                        fontcolor=CONN_MM.MENU{thishandle}.fontcolor(min(size(CONN_MM.MENU{thishandle}.fontcolor,1),x),:);
                        mx1=shiftdim(mean(mean(tximage,1),2),1);mx2=fontcolor;
                        px0=abs(mx1-mx2);
                        fontcolor1=fontcolor.*px0+min(1,mx1+1).*(1-px0); px1=mean(abs(mx1-fontcolor1));
                        fontcolor2=fontcolor.*px0+max(0,mx1-1).*(1-px0); px2=mean(abs(mx1-fontcolor2));
                        fontcolor3=fontcolor.*px0+(1-round(mx1)).*(1-px0); px3=mean(abs(mx1-fontcolor3));
                        px=[mean(px0) px1 px2 0*px3].^4; [nill,ipx]=max(px); px(ipx)=px(ipx)+.1; px=px/max(eps,sum(px));
                        fontcolor=max(0,min(1,px*[fontcolor;fontcolor1;fontcolor2;fontcolor3]));
                        CONN_MM.MENU{thishandle}.handle.axes{n1}(n2)=htemp;
                        %if n2==1, CONN_MM.MENU{thishandle}.handle.axes{n1}=htemp;
                        %else CONN_MM.MENU{thishandle}.handle.axes{n1}=[CONN_MM.MENU{thishandle}.handle.axes{n1}(:)',htemp(:)'];
                        %end
                        if CONN_gui.isjava   % draws menus as uicontrol objects
                            set(htemp,'cdata',tximage);%,'backgroundcolor',mx1);
                            if ~CONN_MM.MENU{thishandle}.enable(n1), set(htemp,'fontsize',fontsize+CONN_gui.font_offset,'fontangle','italic','fontweight','normal','foregroundcolor',.5*fontcolor+.5*mean(ximage(:)));
                            else set(htemp,'fontsize',fontsize+CONN_gui.font_offset,'foregroundcolor',fontcolor,'fontweight',fontweight);
                            end
                            turnon=[turnon htemp(:)']; %set(htemp,'visible','on');
                        else
                            if CONN_MM.MENU{thishandle}.linkon % draws menus as uipanel objects
                                htempn=get(htemp,'userdata');
                                htemp1=htempn(1); htemp2=htempn(2); htemp3=htempn(3);
                                turnon=[turnon htemp htemp2 htemp3]; %set(htemp,'visible','on');
                            else % draws menus as axes objects
                                htempn=get(htemp,'userdata');
                                htemp1=htemp; htemp2=htempn(1); htemp3=htempn(2);
                                turnon=[turnon htemp2 htemp3]; %set(htemp,'visible','on');
                            end
                            set(htemp2,'cdata',tximage);
                            set(htemp1,'xlim',[.5 size(tximage,2)+.5],'ylim',[.5 size(tximage,1)+.5],'ydir','normal','visible','off','xtick',[],'ytick',[],'color',CONN_gui.backgroundcolor,'xcolor',CONN_gui.backgroundcolor,'ycolor',CONN_gui.backgroundcolor);
                            if ~CONN_MM.MENU{thishandle}.enable(n1), set(htemp3,'position',[size(tximage,2)/2 size(tximage,1)/2 0],'fontangle','italic','color',.5*fontcolor+.5*mean(ximage(:)));
                            else set(htemp3,'position',[size(tximage,2)/2 size(tximage,1)/2 0],'color',fontcolor);
                            end
                        end
                        %if CONN_MM.MENU{thishandle}.roll&&numel(changed)==CONN_MM.MENU{thishandle}.n, drawnow; end
                    end
                end
            end
            CONN_MM.MENU{thishandle}.on=0;
            CONN_MM.ACTIVE(thishandle)=1;
            CONN_MM.POSITION(thishandle,:)=CONN_MM.MENU{thishandle}.position;

            bg=get(CONN_MM.gcf,'color');
            %if ~(UIXCONTROLTOOLTIP&&CONN_gui.dounixGUIbugfix)&&CONN_gui.tooltips
            if ~UIXCONTROLTOOLTIP&&CONN_gui.tooltips
                if ishandle(CONN_MM.HELP.handle),
                    if isempty(CONN_MM.HELP.string),
                        if isempty(CONN_x.filename), set(CONN_MM.HELP.handle,'string','Project: undefined','visible','on');
                        else set(CONN_MM.HELP.handle,'string',['Project: ',CONN_x.filename],'visible','on');
                        end
                    else set(CONN_MM.HELP.handle,'string',CONN_MM.HELP.string,'visible','on');
                    end
                else
                    ha=axes('units','norm','position',[.3,.015,.4,.03],'visible','off','parent',CONN_MM.gcf);
                    CONN_MM.HELP.handle=text(0,1,CONN_MM.HELP.string,'color',mod(bg-.3,1),'fontname','default','fontsize',8+CONN_gui.font_offset,'horizontalalignment','center','verticalalignment','middle','interpreter','none','parent',ha);
                    set(ha,'xlim',[-1 1],'ylim',[0 2],'visible','off');
                    %CONN_MM.HELP.handle=uicontrol('style','text','units','norm','position',[.15,.02,.7,.03],'backgroundcolor',bg,'foregroundcolor',1-bg,'fontname','default','string',CONN_MM.HELP.string,'fontsize',8+CONN_gui.font_offset); 
                end
            else
                if ishandle(CONN_MM.HELP.handle),
                    if isempty(CONN_x.filename), set(CONN_MM.HELP.handle,'string','Project: undefined','visible','on');
                    else set(CONN_MM.HELP.handle,'string',['Project: ',CONN_x.filename],'visible','on');
                    end
                else
                    ha=axes('units','norm','position',[.3,.015,.4,.03],'visible','off','parent',CONN_MM.gcf);
                    CONN_MM.HELP.handle=text(0,1,'','color',mod(bg-.3,1),'fontname','default','fontsize',8+CONN_gui.font_offset,'horizontalalignment','center','verticalalignment','middle','interpreter','none','parent',ha);
                    set(ha,'xlim',[-1 1],'ylim',[0 2],'visible','off');
                    %CONN_MM.HELP.handle=uicontrol('style','text','units','norm','position',[.15,.02,.7,.03],'backgroundcolor',bg,'foregroundcolor',1-bg,'fontname','default','string','','fontsize',8+CONN_gui.font_offset); 
                end
            end
            redrawnow=true;
            %if nhandle==length(handle), drawnow; end
        end
    end
    set(turnon,'visible','on');
    %if ~isempty(GCF)&&ishandle(GCF), figure(GCF); end
    %if redrawnow, drawnow; end
end
end

function y=conn_menumanager_buttonshape(params,hfig)
global CONN_gui;
p0=get(hfig,'position'); 
p1=p0(3:4).*params.position(1:2);
p2=p0(3:4).*params.position(3:4); 
nstrings=cellfun(@(x)numel(cellstr(x)),params.string);
if params.order(1)=='h', p2(1)=p2(1)/sum(nstrings); 
else p2(2)=p2(2)/sum(nstrings); 
end
kx=2*ceil(.75*p2(1))+1;ky=2*ceil(.75*p2(2))+1; 
[tx,ty]=meshgrid(1:kx,1:ky); c1=(ky+1)/2; c2=kx-(ky-1)/2;
kz=c1:-1:.0*c1; t=0;
%try, if isequal(datestr(now,'mmdd'),'0401'), params.color=params.color*min(1,sparse(1:3,randperm(3),1,3,3)); end; end
rns=8;
for n1=kz, t=t+double((abs(tx-c1).^rns+abs(ty-c1).^rns)<abs(n1).^rns | (abs(tx-c2).^rns+abs(ty-c1).^rns)<abs(n1).^rns | (tx>=c1&tx<=c2&abs(ty-c1)<n1))/length(kz); end;
t=max(0,min(1,t));
%t=max(0,t-.20).^.125;%.5+.5*tanh(10*(t-.1));%t.^.05;%.25;
t=(t+19*max(0,t-.20).^.0125)/20;%.5+.5*tanh(10*(t-.1));%t.^.05;%.25;
t(tx>=kx-2)=0;
switch(lower(params.bordertype))
    case 'round'
        t2=repmat(.5+0*.5*linspace(1,0,size(t,1))'.^2,1,size(t,2)).*t;
        t1=repmat(.75+0*.25*linspace(1,0,size(t,1))'.^2,1,size(t,2)).*t;
        %bg=[0 0 0];%min(1,4*params.backgroundcolor);
        bg=params.backgroundcolor;
        t={t2,t1,t2};
        if CONN_gui.doemphasis3, transparent=max(0,min(1,min(.0,params.transparent)*[1 .5 0]));
        else transparent=max(0,min(1,min(.0,params.transparent)*[1 .5 0]));
        end
    case 'xround'
        t2=repmat(.5+0*.5*linspace(1,0,size(t,1))'.^2,1,size(t,2)).*t;
        t1=repmat(.75+0*.25*linspace(1,0,size(t,1))'.^2,1,size(t,2)).*t;
        %bg=[0 0 0];%min(1,4*params.backgroundcolor);
        bg=params.backgroundcolor;
        t={t2,t1,t2};
        if CONN_gui.doemphasis3, transparent=max(0,min(1,min(.25,params.transparent)*[1 .5 0]));
        else transparent=max(0,min(1,min(.25,params.transparent)*[1 .5 0]));
        end
    case 'square'
        if 0, t=1+zeros(size(tx)); end
        if 0, t2=1-.5*repmat(linspace(1,-1,size(t,1))'.^8,1,size(t,2)); t1=.75-.5*repmat(linspace(1,-1,size(t,1))'.^8,1,size(t,2));
        elseif strcmp(params.order,'horizontal'), t2=t-0*.5*repmat(linspace(1,-1,size(t,1))'.^10,1,size(t,2)); t1=.75*t-0*.5*repmat(linspace(1,-1,size(t,1))'.^10,1,size(t,2));
        else t2=t-0*.5*repmat(linspace(1,-1,size(t,2)).^50,size(t,1),1); t1=.75*t-0*.5*repmat(linspace(1,-1,size(t,2)).^50,size(t,1),1);
        end
        %if strcmp(params.order,'horizontal'), t2(1:end-10,:)=0;
        %else  t2(:,11:end)=0;
        %end
        bg=params.backgroundcolor;
        t={t1,.25*t1+.75*t2,t2};
        if CONN_gui.doemphasis3, transparent=max(0,min(1, params.transparent*[1 .25 0] ));
        else transparent=max(0,min(1, params.transparent*[1 .25 0] ));
        end
end
y=repmat({cell(1,3)},1,params.n);
for n=1:params.n
for x=1:3, 
    if any(params.color(x,:)<0)
        y{n}{x}=max(0,min(1,cat(3,...
            reshape(interp1([0,1],[bg(1),min(1,abs(params.color(x,1))*bg(1))],t{x}(:)*((x==3)+(x~=3)*(1-params.transparent)),'linear'),size(t{x})),...
            reshape(interp1([0,1],[bg(2),min(1,abs(params.color(x,2))*bg(2))],t{x}(:)*((x==3)+(x~=3)*(1-params.transparent)),'linear'),size(t{x})),...
            reshape(interp1([0,1],[bg(3),min(1,abs(params.color(x,3))*bg(3))],t{x}(:)*((x==3)+(x~=3)*(1-params.transparent)),'linear'),size(t{x})))));
    else
        pos=[p1 p2];
        if params.order(1)=='h', pos(1)=pos(1)+pos(3)*sum(nstrings(1:n-1)); pos(3)=pos(3)*nstrings(n);
        else pos(2)=pos(2)+pos(4)*(sum(nstrings)-sum(nstrings(1:n))); pos(4)=pos(4)*nstrings(n);
        end
        p=min(1,transparent(x)+1-t{x}.^.5); %max(params.transparent,1-t{x}.^.5);
        bg=conn_guibackground('get',pos,size(p));
        y{n}{x}=conn_bsxfun(@plus,conn_bsxfun(@times,p,bg), conn_bsxfun(@times,(1-p).*t{x},shiftdim(params.color(x,:),-1)));
        y{n}{x}=max(0,min(1,y{n}{x}));
%         p=t{x}*((x==0)+(x~=0)*(1-params.transparent));
%         bg=conn_guibackground('get',pos,size(p));
%         y{n}{x}=conn_bsxfun(@plus,conn_bsxfun(@times,1-p,bg), conn_bsxfun(@times,p,shiftdim(params.color(x,:),-1)));
%         y{x}=max(0,min(1,cat(3,...
%             (1-p)*bg(1) + p*params.color(x,1)
%             reshape(interp1([0,1],[bg(1),params.color(x,1)],t{x}(:)*((x==0)+(x~=0)*(1-params.transparent)),'linear'),size(t{x})),...
%             reshape(interp1([0,1],[bg(2),params.color(x,2)],t{x}(:)*((x==0)+(x~=0)*(1-params.transparent)),'linear'),size(t{x})),...
%             reshape(interp1([0,1],[bg(3),params.color(x,3)],t{x}(:)*((x==0)+(x~=0)*(1-params.transparent)),'linear'),size(t{x})))));
    end
end
end
% for x=1:3, y{x}=max(0,min(1,cat(3,reshape(interp1([0,1],[1.15*bg(1),params.color(x,1)],t(:)*((x~=1)+(x==1)*(1-params.transparent)),'linear'),size(t)),...
%         reshape(interp1([0,1],[1.15*bg(2),params.color(x,2)],t(:)*((x~=1)+(x==1)*(1-params.transparent)),'linear'),size(t)),...
%         reshape(interp1([0,1],[1.15*bg(3),params.color(x,3)],t(:)*((x~=1)+(x==1)*(1-params.transparent)),'linear'),size(t)))));
% end
end

function conn_menumanager_mousemove(varargin)
global CONN_MM;
p1=get(0,'pointerlocation');
p2=get(CONN_MM.gcf,'position');
p3=get(0,'screensize');
p2(1:2)=p2(1:2)+p3(1:2)-1; % note: fix issue when connecting to external monitor/projector
pos=(p1-p2(1:2));
set(CONN_MM.gcf,'currentpoint',pos);
if nargin>0, feval(varargin{:}); end
end


