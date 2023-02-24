function fh=conn_referenceexplore(varargin)
% CONN_REFERENCEEXPLORE: GUI for conn_reference
%

global CONN_x CONN_gui;

if isfield(CONN_gui,'font_offset'),font_offset=CONN_gui.font_offset; else font_offset=0; end
hfig=[];
fh=@(varargin)conn_referenceexplore_update([],[],varargin{:});

if ~isfield(CONN_x,'SetupPreproc')||~isfield(CONN_x.SetupPreproc,'log')||isempty(CONN_x.SetupPreproc.log), 
    options.preproclog_text={};
    options.preproclog={};
else
    str1={}; idx1=[];
    for n=1:numel(CONN_x.SetupPreproc.log),
        if isequal(CONN_x.SetupPreproc.log{n}{1},'timestamp'), str1{end+1}=CONN_x.SetupPreproc.log{n}{2}; idx1(end+1)=n; end
    end
    options.preproclog_text=str1;
    options.preproclog=CONN_x.SetupPreproc.log(idx1);
end

options.analyses_listnames={};options.analyses_listtype=[];options.analyses_listidx=[];
txt={CONN_x.Analyses(:).name};
if 1, options.shownanalyses=find(cellfun(@(x)isempty(regexp(x,'^(.*\/|.*\\)?Dynamic factor .*\d+$')),txt));
else options.shownanalyses=1:numel(txt);
end
for n1=reshape(options.shownanalyses,1,[]),
    switch(CONN_x.Analyses(n1).type)
        case 1,
            options.analyses_listnames=[options.analyses_listnames,{CONN_x.Analyses(n1).name}];
            options.analyses_listtype=[options.analyses_listtype, 1];
            options.analyses_listidx=[options.analyses_listidx, n1];
        case 2,
            options.analyses_listnames=[options.analyses_listnames,{CONN_x.Analyses(n1).name}];
            options.analyses_listtype=[options.analyses_listtype, 2];
            options.analyses_listidx=[options.analyses_listidx, n1];
        case 3,
            options.analyses_listnames=[options.analyses_listnames,{CONN_x.Analyses(n1).name}];
            options.analyses_listtype=[options.analyses_listtype, 0];
            options.analyses_listidx=[options.analyses_listidx, n1];
    end
end
for n1=1:numel(CONN_x.vvAnalyses)
    options.analyses_listnames=[options.analyses_listnames,{CONN_x.vvAnalyses(n1).name}];
    options.analyses_listtype=[options.analyses_listtype, 3+1*all(ismember(conn_v2v('fieldtext',CONN_x.vvAnalyses(n1).measures,1),{'3','4'}))+3*all(ismember(conn_v2v('fieldtext',CONN_x.vvAnalyses(n1).measures,1),{'2'}))];
    options.analyses_listidx=[options.analyses_listidx, n1];
end
temp={CONN_x.dynAnalyses.name};
options.analyses_listnames=[options.analyses_listnames,temp(:)'];
options.analyses_listtype=[options.analyses_listtype, 5+zeros(1,numel(temp))];
options.analyses_listidx=[options.analyses_listidx, 1:numel(temp)];
%                 idx=[];
%                 if ~isempty(state)&&state(1)==1, idx=find(options.analyses_listtype==1&options.analyses_listidx==CONN_x.Analysis,1); end     % RRC
%                 if ~isempty(state)&&state(1)==2, idx=find(options.analyses_listtype==2&options.analyses_listidx==CONN_x.Analysis,1); end     % SBC
%                 if ~isempty(state)&&state(1)==3, idx=find(options.analyses_listtype==3&options.analyses_listidx==CONN_x.vvAnalysis,1); end   % V2V
%                 if ~isempty(state)&&state(1)==4, idx=find(options.analyses_listtype==4&options.analyses_listidx==CONN_x.vvAnalysis,1); end   % ICA/PCA
%                 if ~isempty(state)&&state(1)==5, idx=find(options.analyses_listtype==5&options.analyses_listidx==CONN_x.dynAnalysis,1); end  % dyn
%                 if ~isempty(state)&&state(1)==6, idx=find(options.analyses_listtype==6&options.analyses_listidx==CONN_x.vvAnalysis,1); end   % MVPA
%                 if isempty(idx)&&~isempty(state), 
%                     idx=find(options.analyses_listtype==state(1),1); 
%                     if ~isempty(idx)&&(state(1)==1||state(1)==2), CONN_x.Analysis=options.analyses_listidx(idx);
%                     elseif ~isempty(idx)&&(state(1)==3||state(1)==4||state(1)==6), CONN_x.vvAnalysis=options.analyses_listidx(idx);
%                     elseif ~isempty(idx)&&state(1)==5, CONN_x.dynAnalysis=options.analyses_listidx(idx);
%                     end
%                 end 
%options.firstlevelinfo

bgc=.9*[1 1 1];
figcolor=[1 1 1];%[.95 .95 .9];
dlg.handles.hfig=hfig;
if isempty(dlg.handles.hfig)||~ishandle(dlg.handles.hfig), dlg.handles.hfig=figure('units','norm','position',[.1,.3,.8,.6],'menubar','none','numbertitle','off','name','Write methods','color',figcolor,'colormap',gray(256),'interruptible','off','busyaction','cancel','tag','conn_referenceexplore','userdata',fh);
else figure(dlg.handles.hfig); clf(dlg.handles.hfig);
end
uicontrol('style','frame','units','norm','position',[0,.71,1,.29],'backgroundcolor',bgc,'foregroundcolor',bgc,'fontsize',9+font_offset);
dlg.handles.cb1=uicontrol('style','checkbox','units','norm','position',[.025,.925,.20,.05],'backgroundcolor',bgc,'foregroundcolor','k','horizontalalignment','left','value',1,'string','<HTML>Include <b>Preprocessing</b> details</HTML>','fontweight','normal','fontsize',9+font_offset,'callback',{@conn_referenceexplore_update,'preprocessing'},'interruptible','off');
dlg.handles.cb2=uicontrol('style','checkbox','units','norm','position',[.275,.925,.20,.05],'backgroundcolor',bgc,'foregroundcolor','k','horizontalalignment','left','value',1,'string','<HTML>Include <b>Denoising (1st-level)</b> details</HTML>','fontweight','normal','fontsize',9+font_offset,'callback',{@conn_referenceexplore_update,'denoising'},'interruptible','off');
dlg.handles.cb3=uicontrol('style','checkbox','units','norm','position',[.525,.925,.20,.05],'backgroundcolor',bgc,'foregroundcolor','k','horizontalalignment','left','value',1,'string','<HTML>Include <b>Analyses (1st-level)</b> details</HTML>','fontweight','normal','fontsize',9+font_offset,'callback',{@conn_referenceexplore_update,'firstlevel'},'interruptible','off');
dlg.handles.cb4=uicontrol('style','checkbox','units','norm','position',[.775,.925,.20,.05],'backgroundcolor',bgc,'foregroundcolor','k','horizontalalignment','left','value',0,'string','<HTML>Include <b>Results (2nd-level)</b> details</HTML>','fontweight','normal','fontsize',9+font_offset,'callback',{@conn_referenceexplore_update,'secondlevel'},'interruptible','off');
if ~isfield(CONN_x,'SetupPreproc')||~isfield(CONN_x.SetupPreproc,'log')||isempty(CONN_x.SetupPreproc.log), set(dlg.handles.cb1,'value',0,'enable','off'); end
if ~CONN_x.isready(2), set(dlg.handles.cb2,'value',0,'enable','off'); end
if ~CONN_x.isready(3), set(dlg.handles.cb3,'value',0,'enable','off'); end
if ~CONN_x.isready(4), set(dlg.handles.cb4,'value',0,'enable','off'); end
set(dlg.handles.cb4,'visible','off');

dlg.handles.preproclog=uicontrol('style','listbox','units','norm','position',[.025,.75,.20,.17],'max',1,'backgroundcolor',bgc,'foregroundcolor','k','string',options.preproclog_text,'tooltipstring','Select a preprocessing pipeline to include in the methods description','fontsize',9+font_offset,'callback',{@conn_referenceexplore_update,'preprocessing'},'interruptible','off');
dlg.handles.firstlevelinfo=uicontrol('style','listbox','units','norm','position',[.525,.75,.20,.17],'max',2,'backgroundcolor',bgc,'foregroundcolor','k','string',options.analyses_listnames,'tooltipstring','Select a first-level analysis to include in the methods description','fontsize',9+font_offset,'callback',{@conn_referenceexplore_update,'firstlevel'},'interruptible','off');

jBrowserPanel = javaObjectEDT(com.mathworks.mlwidgets.help.LightweightHelpPanel);
[dlg.handles.html, dlg.handles.text] = javacomponent(jBrowserPanel, [], gcf);
set(dlg.handles.text, 'Units','norm','Position',[.025,.05,.95,.625]);
%dlg.handles.text=uicontrol('style','listbox','units','norm','position',[.025,.05,.95,.625],'string','','max',2,'backgroundcolor',figcolor,'foregroundcolor','k','HorizontalAlignment','left','fontsize',32+font_offset);

conn_referenceexplore_update([],[],'init');
if ~ishandle(dlg.handles.hfig), return; end

    function conn_referenceexplore_update(hObject,eventdata,option,varargin)
        if isfield(dlg,'handles')&&isfield(dlg.handles,'hfig')&&~ishandle(dlg.handles.hfig), return; end
        switch(lower(option))
            case {'init','preprocessing','denoising','firstlevel','secondlevel'}
                steps={'init'};
                if get(dlg.handles.cb1,'value')>0, steps{end+1}='preprocessing'; set(dlg.handles.preproclog,'enable','on'); 
                else set(dlg.handles.preproclog,'enable','off'); 
                end                
                if get(dlg.handles.cb2,'value')>0, steps{end+1}='denoising'; end
                if get(dlg.handles.cb3,'value')>0, steps{end+1}='firstlevel'; set(dlg.handles.firstlevelinfo,'enable','on'); 
                else set(dlg.handles.firstlevelinfo,'enable','off');
                end
                if get(dlg.handles.cb4,'value')>0, steps{end+1}='secondlevel'; end
                idx1=get(dlg.handles.preproclog,'value');
                idx2=get(dlg.handles.firstlevelinfo,'value');
                if isempty(idx1)||isempty(options.preproclog), preproclog={};
                else preproclog=options.preproclog{idx1};
                end
                if isempty(idx2)||isempty(options.analyses_listnames), firstlevelinfo={};
                else firstlevelinfo=regexprep(options.analyses_listnames(idx2),' \(SBC\)$| \(RRC\)$','');
                end
                %options.analyses_listtype=[options.analyses_listtype, 1];
                %options.analyses_listidx=[options.analyses_listidx, n1];
                str=conn_reference(steps,'preproclog',preproclog,'firstlevelinfo',firstlevelinfo);
                dlg.handles.html.getLightweightBrowser.load(['file://',str]);
                %set(dlg.handles.text,'string',str);
                %set(dlg.handles.text,'string',cellfun(@(x)['<HTML>',x,'</HTML>'],regexp(str,'(\n|<br>)+','split'),'uni',0));
                %set(dlg.handles.text,'string',['<HTML>',regexprep(str,'\n',''),'</HTML>']);
        end
    end
end



