function fh=conn_qaplotsexplore(varargin)
% CONN_QAPLOTSEXPLORE: QA plots display
%

global CONN_x CONN_gui;

if isfield(CONN_gui,'font_offset'),font_offset=CONN_gui.font_offset; else font_offset=0; end
if ~isfield(CONN_x,'folders')||~isfield(CONN_x.folders,'qa')||isempty(CONN_x.folders.qa),
    qafolder=pwd;
    isCONN=false;
else qafolder=CONN_x.folders.qa;
    isCONN=true;
end
fh=@(varargin)conn_qaplotsexplore_update([],[],varargin{:});
if nargin&&any(strcmp(varargin(cellfun(@ischar,varargin)),'initdenoise')), dlg.forceinitdenoise=true;
else dlg.forceinitdenoise=false;
end
if nargin&&any(strcmp(varargin(cellfun(@ischar,varargin)),'createdenoise')), dlg.forceinitdenoise=true; dlg.createdenoise=true;
else dlg.createdenoise=false;
end
if nargin&&any(strcmp(varargin(cellfun(@ischar,varargin)),'createreport')), dlg.createreport=true;
else dlg.createreport=false;
end
if nargin&&any(strcmp(varargin(cellfun(@ischar,varargin)),'keepcurrentfigure')), dlg.keepcurrentfigure=true;
else dlg.keepcurrentfigure=false;
end
if nargin&&any(strcmp(varargin(cellfun(@ischar,varargin)),'overwritecurrentfigure')), dlg.overwritecurrentfigure=true;
else dlg.overwritecurrentfigure=false;
end
hfig=[];
if dlg.overwritecurrentfigure||dlg.keepcurrentfigure
    hfig=findobj('tag','conn_qaplotsexplore');
    if ~isempty(hfig),
        hfig=hfig(end);
        if dlg.keepcurrentfigure,
            figure(hfig);
            fh=get(hfig,'userdata');
            fh(varargin{~strcmp(varargin(cellfun(@ischar,varargin)),'keepcurrentfigure')});
            return;
        end
    end
end
if nargin&&any(strcmp(varargin(cellfun(@ischar,varargin)),'thesefolders')),   % each subfolder within current folder contains a different QA report
    qafolder=pwd;
    qagroups=conn_sortfilenames(conn_dir(fullfile(qafolder,'QA_*'),'-dir','-cell'));
    dlg.sets=regexprep(cellfun(@(x)x(numel(qafolder)+1:end),qagroups,'uni',0),'^[\\\/]','');
    isCONN=false;
elseif nargin&&any(strcmp(varargin(cellfun(@ischar,varargin)),'thisfolder')), % single QA report in current folder
    [qafolder,qagroups]=fileparts(pwd);
    dlg.sets={qagroups};
    isCONN=false;
elseif nargin&&any(strcmp(varargin(cellfun(@ischar,varargin)),'flfolders')),  % each conn folder contains a different subject
    if ischar(varargin{1})&&conn_fileutils('isdir',varargin{1}), qafolder=varargin{1};
    else qafolder=pwd;
    end
    qagroups=conn_sortfilenames(conn_dir(fullfile(qafolder,'QA_*'),'-dir','-cell'));
    if isempty(qagroups), conn_disp('no QA plots found'); return; end
    dlg.sets=regexprep(cellfun(@(x)x(numel(qafolder)+1:end),qagroups,'uni',0),'^[\\\/]','');
    nfilesep=sum(char(dlg.sets)==filesep,2);
    mnfilesep=mode(nfilesep);
    switch(mnfilesep)
        case {0,1} % unknown
            isCONN=false;
        case 2,   % single subject
            conn_disp('fl single subject display');
        case 3,   % single dataset
            conn_disp('fl single dataset display');
            [qafolder,qagroups]=fileparts(qafolder);
            if isempty(qagroups), [qafolder,qagroups]=fileparts(qafolder); end
            if isempty(qafolder), qafolder=pwd; end
            dlg.sets={qagroups};
            isCONN=false;
        otherwise % multiple datasets
            conn_disp('fl multiple dataset display');
            for nk=1:numel(dlg.sets), ik=find(dlg.sets{nk}==filesep); dlg.sets{nk}=dlg.sets{nk}(1:ik(end-3)-1); end
            dlg.sets=unique(dlg.sets);
            isCONN=false;
    end
elseif nargin&&ischar(varargin{1})&&conn_fileutils('isdir',varargin{1})                        % single QA report in specified folder
    [qafolder,qagroups]=fileparts(varargin{1});
    if isempty(qagroups), [qafolder,qagroups]=fileparts(qafolder); end
    if isempty(qafolder), qafolder=pwd; end
    dlg.sets={qagroups};
    isCONN=false;
else
    qagroups=conn_sortfilenames(conn_dir(fullfile(qafolder,'QA_*'),'-dir','-cell'));
    dlg.sets=regexprep(cellfun(@(x)x(numel(qafolder)+1:end),qagroups,'uni',0),'^[\\\/]','');
end
if nargin&&any(strcmp(varargin(cellfun(@ischar,varargin)),'isconn')), isCONN=true; end
if dlg.forceinitdenoise&&(dlg.createdenoise||isempty(dlg.sets)),
    conn_qaplotsexplore_update([],[],'newsetinit');
end
dlg.iset=numel(dlg.sets);
dlg.dispsize=[];
dlg.showavg=1;
bgc=.9*[1 1 1];
dlg.handles.fh=fh;
figcolor=[.95 .95 .9];
dlg.handles.hfig=hfig;
if isempty(dlg.handles.hfig)||~ishandle(dlg.handles.hfig), dlg.handles.hfig=figure('units','norm','position',[.1,.3,.8,.6],'menubar','none','numbertitle','off','name','Quality Assurance reports','color',figcolor,'colormap',gray(256),'interruptible','off','busyaction','cancel','tag','conn_qaplotsexplore','userdata',fh);
else figure(dlg.handles.hfig); clf(dlg.handles.hfig);
end
%dlg.handles.menuprint=uimenu(dlg.handles.hfig,'Label','Print');
%uimenu(dlg.handles.menuprint,'Label','Print entire report (pdf/html)','callback',{@conn_qaplotsexplore_update,'createreport'});
%uimenu(dlg.handles.menuprint,'Label','Print individual plot (jpeg/tiff)','callback',{@conn_qaplotsexplore_update,'print'});
uicontrol('style','frame','units','norm','position',[0,.91,1,.09],'backgroundcolor',bgc,'foregroundcolor',bgc,'fontsize',9+font_offset);
uicontrol('style','text','units','norm','position',[.025,.925,.1,.05],'backgroundcolor',bgc,'foregroundcolor','k','horizontalalignment','left','string','QA reports:','fontweight','bold','fontsize',9+font_offset);
uicontrol('style','text','units','norm','position',[.025,.835,.1,.05],'backgroundcolor',figcolor,'foregroundcolor','k','horizontalalignment','left','string','Plots:','fontweight','bold','fontsize',9+font_offset);
dlg.handles.set=uicontrol('style','popupmenu','units','norm','position',[.125,.925,.5,.05],'string',dlg.sets,'value',dlg.iset,'backgroundcolor',bgc,'foregroundcolor','k','tooltipstring','<HTML>Select a Quality Assurance report<br/> - each report contains one or multiple plots created to visually assess the quality of the structural/functional data<br/> and/or easily identify potential outlier subjects or failed preprocessing steps<br/> - choose one existing report from this list, or select <i>new report</i> to create a new report instead</HTML>','callback',{@conn_qaplotsexplore_update,'set'},'fontsize',9+font_offset,'interruptible','off');
dlg.handles.settxt=uicontrol('style','text','units','norm','position',[.125,.925,.5,.05],'string','No QA sets found in this CONN project. Select ''New report'' to get started','backgroundcolor',bgc,'foregroundcolor','k','fontsize',9+font_offset,'visible','off');
dlg.handles.analysis=uicontrol('style','popupmenu','units','norm','position',[.125,.835,.5,.05],'string',' ','backgroundcolor',figcolor,'foregroundcolor','k','tooltipstring','<HTML>Select a Quality Assurance plot within this report<br/> - each report may contain one or multiple plots<br/> - choose one existing plot from this report, or select <i>new plot</i> to create a new plot and add it to this report</HTML>','callback',{@conn_qaplotsexplore_update,'plot'},'fontsize',9+font_offset,'interruptible','off');
dlg.handles.analysistxt=uicontrol('style','text','units','norm','position',[.125,.835,.5,.05],'string','No plots found in this QA report. Select ''New plot'' to get started','backgroundcolor',figcolor,'foregroundcolor','k','fontsize',9+font_offset,'visible','off');
dlg.handles.recomputeset=uicontrol('style','pushbutton','units','norm','position',[.625 .925 .12 .05],'string','Recreate report','tooltipstring','Recomputes all plots in current QA report','fontsize',9+font_offset,'callback',{@conn_qaplotsexplore_update,'recreateset'},'interruptible','off');
dlg.handles.recomputeplot=uicontrol('style','pushbutton','units','norm','position',[.625 .835 .12 .05],'string','Recreate plot','tooltipstring','Recomputes the current plot','fontsize',9+font_offset,'callback',{@conn_qaplotsexplore_update,'recreateplot'},'interruptible','off');
dlg.handles.addnewset=uicontrol('style','pushbutton','units','norm','position',[.745 .925 .12 .05],'string','Create new report','tooltipstring','Starts a new QA report','fontsize',9+font_offset,'callback',{@conn_qaplotsexplore_update,'newset'},'interruptible','off');
dlg.handles.addnewplot=uicontrol('style','pushbutton','units','norm','position',[.745 .835 .12 .05],'string','Create new plot','tooltipstring','Adds a new plot to the current report','fontsize',9+font_offset,'callback',{@conn_qaplotsexplore_update,'newplot'},'interruptible','off');
dlg.handles.deleteset=uicontrol('style','pushbutton','units','norm','position',[.875 .925 .05 .05],'string','Delete','tooltipstring','Deletes this report (and all plots in it)','fontsize',9+font_offset,'callback',{@conn_qaplotsexplore_update,'deleteset'},'interruptible','off');
dlg.handles.deleteplot=uicontrol('style','pushbutton','units','norm','position',[.875 .835 .05 .05],'string','Delete','tooltipstring','Removes this plot from the current report','fontsize',9+font_offset,'callback',{@conn_qaplotsexplore_update,'deleteplot'},'interruptible','off');
dlg.handles.printset=uicontrol('style','pushbutton','units','norm','position',[.925 .925 .05 .05],'string','Export','tooltipstring','Exports entire report to html/pdf file','fontsize',9+font_offset,'callback',{@conn_qaplotsexplore_update,'printset'},'interruptible','off');
dlg.handles.printplot=uicontrol('style','pushbutton','units','norm','position',[.925 .835 .05 .05],'string','Print','tooltipstring','Prints current plot (to jpeg/tiff file)','fontsize',9+font_offset,'callback',{@conn_qaplotsexplore_update,'printplot'},'interruptible','off');

%uicontrol('style','frame','units','norm','position',[0,0,.17,.86],'backgroundcolor',figcolor,'foregroundcolor',figcolor,'fontsize',9+font_offset);
dlg.handles.subjects=uicontrol('style','listbox','units','norm','position',[.035,.57,.1,.16],'max',2,'backgroundcolor',figcolor,'foregroundcolor','k','string','','tooltipstring','Select one or multiple subjects for display','fontsize',9+font_offset,'callback',{@conn_qaplotsexplore_update,'subjects'},'interruptible','off');
dlg.handles.sessions=uicontrol('style','listbox','units','norm','position',[.035,.41,.1,.10],'max',2,'backgroundcolor',figcolor,'foregroundcolor','k','string','','tooltipstring','Select one or multiple sessions for display','fontsize',9+font_offset,'callback',{@conn_qaplotsexplore_update,'sessions'},'interruptible','off');
dlg.handles.measures=uicontrol('style','listbox','units','norm','position',[.035,.25,.1,.10],'max',2,'backgroundcolor',figcolor,'foregroundcolor','k','string','','tooltipstring','Select one or multiple measures for display','fontsize',9+font_offset,'callback',{@conn_qaplotsexplore_update,'measures'},'interruptible','off');
dlg.handles.selectall1=uicontrol('style','pushbutton','units','norm','position',[.035 .52 .1 .05],'string','Select all','tooltipstring','Selects all subjects','fontsize',9+font_offset,'callback',{@conn_qaplotsexplore_update,'allsubjects'},'interruptible','off');
dlg.handles.selectall2=uicontrol('style','pushbutton','units','norm','position',[.035 .36 .1 .05],'string','Select all','tooltipstring','Selects all sessions','fontsize',9+font_offset,'callback',{@conn_qaplotsexplore_update,'allsessions'},'interruptible','off');
dlg.handles.selectall3=uicontrol('style','pushbutton','units','norm','position',[.035 .20 .1 .05],'string','Select all','tooltipstring','Selects all measures','fontsize',9+font_offset,'callback',{@conn_qaplotsexplore_update,'allmeasures'},'interruptible','off');
dlg.handles.showannot=uicontrol('style','checkbox','units','norm','position',[.02 .12 .15 .05],'string','show annotations','backgroundcolor',figcolor,'foregroundcolor','k','tooltipstring','show/hide plot annotations','fontsize',9+font_offset,'callback',{@conn_qaplotsexplore_update,'subjects'},'interruptible','off','visible','on');
dlg.handles.showdiff=uicontrol('style','checkbox','units','norm','position',[.02 .07 .15 .05],'string','show diff','backgroundcolor',figcolor,'foregroundcolor','k','fontsize',9+font_offset,'tooltipstring','<HTML>Show z score maps: normalized differences between each image and the average of all images in this plot<br/>this can be useful to highlight differences between a selected subject and the group</HTML>','value',0,'callback',{@conn_qaplotsexplore_update,'subjects'},'interruptible','off');
dlg.handles.invertim=uicontrol('style','checkbox','units','norm','position',[.02 .02 .15 .05],'string','transparent display','value',1,'backgroundcolor',figcolor,'foregroundcolor','k','tooltipstring','invert image colors if necessary to keep light background and dark foreground','fontsize',9+font_offset,'callback',{@conn_qaplotsexplore_update,'refresh'},'interruptible','off','visible','on');
dlg.handles.hax=[];%axes('units','norm','position',[.20 .10 .75 .70],'visible','off');
dlg.handles.han=[];
dlg.handles.text1=uicontrol('style','text','units','norm','position',[.30,.035,.52,.04],'backgroundcolor',figcolor,'foregroundcolor','k','horizontalalignment','center','string','','fontsize',9+font_offset);
dlg.handles.text2=uicontrol('style','text','units','norm','position',[.30,.005,.52,.03],'backgroundcolor',figcolor,'foregroundcolor','k','horizontalalignment','center','string','','fontsize',6+font_offset);
dlg.handles.textoptions=uicontrol('style','popupmenu','units','norm','position',[.45,.025,.2,.05],'backgroundcolor',figcolor,'foregroundcolor','k','horizontalalignment','center','string','','fontsize',9+font_offset,'tooltipstring','<HTML>Choose what to display when selecting multiple subjects/sessions<br/> -<i>average/variability</i> computes the average/variability of the selected plots (e.g. across multiple subjects or sessions)</HTML>','visible','off','callback',{@conn_qaplotsexplore_update,'textoptions'},'interruptible','off');
dlg.handles.details=uicontrol('style','pushbutton','units','norm','position',[.825 .015 .15 .07],'string','Details','tooltipstring','go to interactive or high-resolution version of this plot','fontsize',10+font_offset,'fontweight','bold','callback',{@conn_qaplotsexplore_update,'details'},'visible','off');
conn_qaplotsexplore_update([],[],'set',dlg.iset);
conn_qaplotsexplore_update([],[],'selectall');
%conn_qaplotsexplore_update([],[],'allsubjects');
%conn_qaplotsexplore_update([],[],'allsessions');
%conn_qaplotsexplore_update([],[],'allmeasures');
dlg.handles.hlabel=uicontrol('style','text','horizontalalignment','left','visible','off','fontsize',9+font_offset);
if ~ishandle(dlg.handles.hfig), return; end
set(dlg.handles.hfig,'units','pixels','windowbuttonmotionfcn',@conn_qaplotsexplore_figuremousemove,'windowbuttonupfcn',{@conn_qaplotsexplore_figuremousemove,'buttonup'},'resizefcn',{@conn_qaplotsexplore_update,'resize'});
if dlg.createreport, conn_qaplotsexplore_update([],[],'printset','nogui'); conn_qaplotsexplore_update([],[],'close'); end

    function conn_qaplotsexplore_update(hObject,eventdata,option,varargin)
        if isfield(dlg,'handles')&&isfield(dlg.handles,'hfig')&&~ishandle(dlg.handles.hfig), return; end
        switch(lower(option))
            case 'resize',
                try
                    if dlg.dispsize(end)>1, conn_qaplotsexplore_update([],[],'refresh'); end
                end
            case 'close'
                delete(dlg.handles.hfig);
            case 'printplot'
                conn_qaplotsexplore_update([],[],'togglegui');
                conn_print;
                conn_qaplotsexplore_update([],[],'togglegui');
            case 'printset'
                if isdeployed, conn_msgbox({'Sorry, exporting QA reports to html/pdf is still unavailable in standalone applications','We are work to resolve this issue'},'exporting report',2); return; end
                cwd=conn_projectmanager('pwd');
                if numel(varargin)>0&&any(strcmp(varargin(cellfun(@ischar,varargin)),'nogui')), nogui=true;
                else nogui=false;
                end
                outputDir=fullfile(qafolder,dlg.sets{dlg.iset});
                if ispc, options={'*.html','HTML report (*.html)'; '*.pdf','PDF report (*.pdf)'; '*.ppt','PowerPoint report (*.ppt)'; '*.doc','Word report (*.doc)'; '*.xml','XML report (*.xml)'; '*.tex','LaTeX report (*.tex)'};
                else options={'*.html','HTML report (*.html)'; '*.pdf','PDF report (*.pdf)'; '*.xml','XML report (*.xml)'; '*.tex','LaTeX report (*.tex)'};
                end
%                 if nogui, answ=options{1,2};
%                 else answ=conn_questdlg('','report type',options{:,2},options{1,2});
%                 end
%                 if isempty(answ), return; end
%                 tfilename=conn_prepend('',fullfile(outputDir,'report.html'),regexprep(options{strmatch(answ,options(:,2),'exact'),1},'*',''));
                qaDir=outputDir;
                if nogui, tfilename=fullfile(outputDir,'report.html');
                else [tfilename,outputDir]=uiputfile(options(:,1:2),'Output report to file:',fullfile(outputDir,'report.html'));
                end
                if ~ischar(tfilename)||isempty(tfilename), return; end
                plots_summary=find(dlg.uanalysestype~=1); 
                plots_subject=1:numel(dlg.uanalysestype);
                if ~nogui, 
                    plots_summary=listdlg('liststring',dlg.uanalyses_long,'selectionmode','multiple','initialvalue',plots_summary,'promptstring',{'Select plot(s) to include in SUMMARY section','(in this section each plot will display all subjects/measures combined)'},'ListSize',[600 250]);
                    plots_subject=listdlg('liststring',dlg.uanalyses_long,'selectionmode','multiple','initialvalue',plots_subject,'promptstring',{'Select plot(s) to include in SUBJECTS section','(in this section a separate plot per individual subject will be displayed)'},'ListSize',[600 250]);
                end
                if isempty(plots_summary)&&isempty(plots_subject), return; end
                [nill,nill,tfileext]=fileparts(tfilename);
                format=regexprep(tfileext,{'^\.','^tex$'},{'','latex'});
                invertim=get(dlg.handles.invertim,'value')==0;
                hmsg=conn_msgbox({sprintf('Generating %s report',format),'This may take several minutes. Please wait...'});
                conn_fileutils('cd',outputDir);
                if strcmp(format,'html')&&~isempty(which('private/mxdom2simplehtml.xsl')),
                    stylesheet=conn_prepend('',fullfile(outputDir,tfilename),'_style.xsl');
                    stylestr=textread(which('private/mxdom2simplehtml.xsl'),'%s');
                    stylestr=regexprep(stylestr,'background:#fff',['background:#',sprintf('%x',ceil(figcolor*255))]);
                    %tfh=fopen(stylesheet,'wt');for n1=1:numel(stylestr), fprintf(tfh,'%s\n',stylestr{n1}); end; fclose(tfh);
                    conn_fileutils('filewrite',stylesheet,stylestr);
                    stylesheet={'stylesheet',stylesheet};
                else stylesheet={};
                end
                %
                filename=conn_prepend('',fullfile(outputDir,tfilename),'.m');
                tfh={};
                %tfh=fopen(filename,'wt');
                tfh{end+1}=sprintf('%%%% Quality Control report\n');
                tfh{end+1}=sprintf('%% %s\n',qaDir);
                tfh{end+1}=sprintf('%%\n%% auto-generated by <http://www.conn-toolbox.org CONN> @ %s\n',datestr(now));
                tfh{end+1}=sprintf('%%\n%% <matlab:conn_qaplotsexplore(''%s'') [open report in Matlab]>\n',regexprep(qaDir,' ','%20'));
                tfh{end+1}=sprintf('\n%%%% *SUMMARY PLOTS*\n');
                maxns=0;maxnm=0;
                for np=1:numel(dlg.uanalyses)
                    in=find(ismember(dlg.ianalyses,np));
                    usubjects_shown=unique(dlg.isubjects(in));
                    usessions_shown=unique(dlg.isessions(in));
                    umeasures_shown=unique(dlg.imeasures(in));
                    maxns=max(maxns,numel(usubjects_shown));
                    maxnm=max(maxnm,numel(umeasures_shown));
                    if ismember(np,plots_summary)
                        tfh{end+1}=sprintf('\n%%%% %s\n',dlg.uanalyses_long{np});
                        tfh{end+1}=sprintf('%% <matlab:options=conn_qaplotsexplore(''%s'');options(''plot'',%d);options(''selectall''); [open plot in Matlab]>\n',regexprep(qaDir,' ','%20'),np);
                        if numel(usubjects_shown)>1, tfh{end+1}=sprintf('%% %d subjects\n',numel(usubjects_shown)); end
                        if numel(usessions_shown)>1, tfh{end+1}=sprintf('%% %d sessions\n',numel(usessions_shown)); end
                        if numel(umeasures_shown)>1, tfh{end+1}=sprintf('%% %d measures\n',numel(umeasures_shown)); end
                        tfh{end+1}=sprintf('options=conn_qaplotsexplore(''%s'');\n',regexprep(qaDir,' ','%20'));
                        if invertim, tfh{end+1}=sprintf('options(''invertimage'',0);\n'); end
                        tfh{end+1}=sprintf('options(''plot'',%d);\n',np);
                        tfh{end+1}=sprintf('options(''selectall'');\n');
                        tfh{end+1}=sprintf('options(''displayannotation'');\n');
                        tfh{end+1}=sprintf('options(''togglegui'');\n');
                        tfh{end+1}=sprintf('snapnow;\n');
                        tfh{end+1}=sprintf('options(''close'');\n');
                    end
                end
                if ~isempty(plots_subject)
                    tfh{end+1}=sprintf('\n%%%% *INDIVIDUAL SUBJECT PLOTS*\n');
                    for np=plots_subject(:)'
                        in=find(ismember(dlg.ianalyses,np));
                        usubjects_shown=unique(dlg.isubjects(in));
                        usessions_shown=unique(dlg.isessions(in));
                        umeasures_shown=unique(dlg.imeasures(in));
                        if numel(usubjects_shown)>1||(numel(usubjects_shown)==1&&numel(umeasures_shown)==1)
                            tfh{end+1}=sprintf('\n%%%% %s\n',regexprep(dlg.uanalyses_long{np},'\(\d+\)\s*$',''));
                            for ns=1:numel(usubjects_shown)
                                tfh{end+1}=sprintf('\n%%%%\n%% <matlab:options=conn_qaplotsexplore(''%s'');options(''plot'',%d);options(''preselectall'');options(''selectsubjects'',%d); [%s]>\n',regexprep(qaDir,' ','%20'),np,ns,dlg.usubjects{usubjects_shown(ns)});
                                if ns==1,
                                    tfh{end+1}=sprintf('options=conn_qaplotsexplore(''%s'');\n',regexprep(qaDir,' ','%20'));
                                    if invertim, tfh{end+1}=sprintf('options(''invertimage'',0);\n'); end
                                    tfh{end+1}=sprintf('options(''plot'',%d);\n',np);
                                else tfh{end+1}=sprintf('options(''togglegui'');\n');
                                end
                                tfh{end+1}=sprintf('options(''preselectall'');\n');
                                tfh{end+1}=sprintf('options(''selectsubjects'',%d);\n',ns);
                                tfh{end+1}=sprintf('options(''displayannotation'');\n');
                                tfh{end+1}=sprintf('options(''togglegui'');\n');
                                tfh{end+1}=sprintf('snapnow;\n');
                            end
                            tfh{end+1}=sprintf('options(''close'');\n');
                        elseif numel(umeasures_shown)>1
                            tfh{end+1}=sprintf('\n%%%% %s\n',regexprep(dlg.uanalyses_long{np},'\(\d+\)\s*$',''));
                            for ns=1:numel(umeasures_shown)
                                tfh{end+1}=sprintf('\n%%%%\n%% <matlab:options=conn_qaplotsexplore(''%s'');options(''plot'',%d);options(''preselectall'');options(''selectmeasures'',%d); [measure %s]>\n',regexprep(qaDir,' ','%20'),np,ns,dlg.umeasures{umeasures_shown(ns)});
                                if ns==1,
                                    tfh{end+1}=sprintf('options=conn_qaplotsexplore(''%s'');\n',regexprep(qaDir,' ','%20'));
                                    if invertim, tfh{end+1}=sprintf('options(''invertimage'',0);\n'); end
                                    tfh{end+1}=sprintf('options(''plot'',%d);\n',np);
                                else tfh{end+1}=sprintf('options(''togglegui'');\n');
                                end
                                tfh{end+1}=sprintf('options(''preselectall'');\n');
                                tfh{end+1}=sprintf('options(''selectmeasures'',%d);\n',ns);
                                tfh{end+1}=sprintf('options(''displayannotation'');\n');
                                tfh{end+1}=sprintf('options(''togglegui'');\n');
                                tfh{end+1}=sprintf('snapnow;\n');
                            end
                            tfh{end+1}=sprintf('options(''close'');\n');
                        end
                    end
                end
                %fclose(tfh);
                conn_fileutils('filewrite_raw',filename,tfh);
                if conn_server('util_isremotefile',filename)
                    conn_server('run','publish',conn_server('util_localfile',filename),struct('format',format,stylesheet{:},'showCode',false,'useNewFigure',false,'figureSnapMethod','getframe','outputDir',conn_server('util_localfile',outputDir)));
                else
                    publish(filename,struct('format',format,stylesheet{:},'showCode',false,'useNewFigure',false,'figureSnapMethod','getframe','outputDir',outputDir));
                end
                conn_disp('fprintf','Report created\nOpen %s to view\n',fullfile(outputDir,tfilename));
                if ishandle(hmsg), delete(hmsg); end
                tfile=fullfile(outputDir,tfilename);
                if conn_server('util_isremotefile',tfile), tfile=conn_cache('pull',tfile); end
                switch(format)
                    case 'html', web(tfile);
                    case 'pdf',  if ispc, system(sprintf('open "%s"',tfile));
                        else system(sprintf('open ''%s''',tfile));
                        end
                    otherwise, try, system(sprintf('open %s',tfile)); end
                end
                conn_fileutils('cd',cwd);
            case 'togglegui'
                if ~isfield(dlg,'togglegui')||isempty(dlg.togglegui)
                    h=findobj(dlg.handles.hfig,'type','uicontrol','visible','on');
                    if ishandle(dlg.handles.hax), pos=get(dlg.handles.hax,'position'); else pos=[]; end
                    dlg.togglegui=struct('handles',h,'position',pos);
                    set(dlg.togglegui.handles,'visible','off');
                    set(dlg.handles.hlabel,'visible','off');
                    if ishandle(dlg.handles.hax), set(dlg.handles.hax,'position',[.1 .2 .8 .7]); end
                else
                    set(dlg.togglegui.handles,'visible','on');
                    if ishandle(dlg.handles.hax)&&~isempty(dlg.togglegui.position), set(dlg.handles.hax,'position',dlg.togglegui.position); end
                    dlg.togglegui=[];
                end
            case 'annotate'
                descr=get(dlg.handles.han(end),'string');
                fname=dlg.files_txt{dlg.dataIDXplots(dlg.dataIDXsubjects)};
                dlg.txtA{dlg.dataIDXsubjects}=descr;
                conn_fileutils('filewrite',fname,regexprep(descr,'\n',''));
                %fh=fopen(fname,'wt');
                %for n=1:numel(descr), fprintf(fh,'%s\n',regexprep(descr{n},'\n','')); end
                %fclose(fh);
            case 'invertimage',
                set(dlg.handles.invertim,'value',varargin{1});
            case {'preselectall','selectall'}
                set(dlg.handles.subjects,'value',1:numel(cellstr(get(dlg.handles.subjects,'string'))));
                set(dlg.handles.sessions,'value',1:numel(cellstr(get(dlg.handles.sessions,'string'))));
                set(dlg.handles.measures,'value',1:numel(cellstr(get(dlg.handles.measures,'string'))));
                if strcmpi(option,'selectall'), conn_qaplotsexplore_update([],[],'refresh'); end
            case 'selectsubjects'
                set(dlg.handles.subjects,'value',varargin{1});
                conn_qaplotsexplore_update([],[],'subjects');
            case 'selectsessions'
                set(dlg.handles.sessions,'value',varargin{1});
                conn_qaplotsexplore_update([],[],'sessions');
            case 'selectmeasures'
                set(dlg.handles.measures,'value',varargin{1});
                conn_qaplotsexplore_update([],[],'measures');
            case 'allsubjects'
                set(dlg.handles.subjects,'value',1:numel(cellstr(get(dlg.handles.subjects,'string'))));
                conn_qaplotsexplore_update([],[],'subjects');
            case 'allsessions'
                set(dlg.handles.sessions,'value',1:numel(cellstr(get(dlg.handles.sessions,'string'))));
                conn_qaplotsexplore_update([],[],'sessions');
            case 'allmeasures'
                set(dlg.handles.measures,'value',1:numel(cellstr(get(dlg.handles.measures,'string'))));
                conn_qaplotsexplore_update([],[],'measures');
            case 'textoptions',
                dlg.showavg=get(dlg.handles.textoptions,'value');
                conn_qaplotsexplore_update([],[],'subjects');
            case {'deleteset','newset','newsetinit'}
                if strcmp(lower(option),'deleteset')
                    answ=conn_questdlg(sprintf('Are you sure you want to delete report %s?',dlg.sets{dlg.iset}),'','Delete','Cancel','Delete');
                    if ~isequal(answ,'Delete'), return; end
                    f=conn_dir(fullfile(qafolder,dlg.sets{dlg.iset},'QA_*'),'-R');
                    if isempty(f),
                        f=conn_dir(fullfile(qafolder,dlg.sets{dlg.iset},'QA_*'),'\.mat$|\.txt|\.jpg');
                        if ~isempty(f),
                            answ=listdlg('liststring',cellstr(f),'selectionmode','multiple','initialvalue',[],'promptstring','Select file(s) to delete:','ListSize',[600 250]);
                            if isempty(answ), f='';
                            else f=f(answ,:);
                            end
                        end
                    end
                    if ~isempty(f),
                        f=cellstr(f);
                        conn_fileutils('spm_unlink',f{:});
                    end
                    if strncmp(fullfile(qafolder,dlg.sets{dlg.iset}),pwd,numel(fullfile(qafolder,dlg.sets{dlg.iset}))), cd(qafolder); end
                    [ok,nill]=conn_fileutils('rmdir_dironly',fullfile(qafolder,dlg.sets{dlg.iset}));
                    try, if ~ok, conn_disp(['warning: unable to remove report directory; ',nill]); end; end
                    tag='';
                else
                    if numel(varargin)>=1, answ={varargin{1}};
                    elseif strcmp(lower(option),'newsetinit'), answer={[]};
                    else answer=inputdlg({'Name of new QA report: (must be valid folder name)'},'',1,{datestr(now,'yyyy_mm_dd_HHMMSSFFF')});
                    end
                    if isempty(answer), return; end
                    if isempty(answer{1}), answer={datestr(now,'yyyy_mm_dd_HHMMSSFFF')}; end
                    tag=['QA_',answer{1}];
                    conn_fileutils('mkdir',qafolder,tag);
                end
                qagroups=conn_sortfilenames(conn_dir(fullfile(qafolder,'QA_*'),'-dir','-cell'));
                dlg.sets=regexprep(cellfun(@(x)x(numel(qafolder)+1:end),qagroups,'uni',0),'^[\\\/]','');
                [nill,qanames]=cellfun(@fileparts,dlg.sets,'uni',0);
                dlg.iset=find(strcmp(qanames,tag),1);
                if isempty(dlg.iset), dlg.iset=1; end
                if ~strcmp(lower(option),'newsetinit')
                    set(dlg.handles.set,'string',dlg.sets);
                    conn_qaplotsexplore_update([],[],'set',dlg.iset);
                end
            case 'deleteplot'
                answ=conn_questdlg(sprintf('Are you sure you want to delete plot %s?',dlg.uanalyses_long{dlg.ianalysis}),'','Delete','Cancel','Delete');
                if ~isequal(answ,'Delete'), return; end
                in=find(ismember(dlg.ianalyses,dlg.ianalysis));
                tfiles={};
                for n=1:numel(in),
                    tfiles{end+1}=dlg.files_jpg{in(n)};
                    tfiles{end+1}=dlg.files{in(n)};
                end
                conn_fileutils('spm_unlink',tfiles{:});
                conn_qaplotsexplore_update([],[],'set');
            case {'newplot','recreateplot','recreateset'}
                if ~isCONN, conn_msgbox('Load existing CONN project before proceeding','Error',2); return; end
                tag=dlg.sets{dlg.iset};
                analyses={'QA_COV','QA_NORM_structural','QA_NORM_functional','QA_NORM_ROI','QA_REG_functional','QA_REG__structural','QA_REG__functional','QA_REG__mni','QA_COREG_functional','QA_TIME_functional','QA_TIMEART_functional','QA_DENOISE_timeseries','QA_DENOISE_QC-FC','QA_DENOISE','QA_DENOISE_scatterplot','QA_DENOISE_QC-FC_scatterplot','QA_SPM_design','QA_SPM_contrasts','QA_SPM_results'};
                %analyses_numbers=[31,1,2,3,10,4,5,6,7,8,9,12,13,11,14,15,21,22,23];
                defaultset=[1,2,4,5,9,11,12,13,31];
                if ~isfield(CONN_x,'isready')||~CONN_x.isready(2),
                    dokeep=~cellfun('length',regexp(lower(analyses),'denoise'));
                    analyses=analyses(dokeep);%analyses_numbers=analyses_numbers(dokeep); % disable analyses that require having run Setup step
                end
                if isempty(CONN_x.Setup.spm)||isempty(CONN_x.Setup.spm{1})||isempty(CONN_x.Setup.spm{1}{1})
                    dokeep=~cellfun('length',regexp(lower(analyses),'spm'));
                    analyses=analyses(dokeep);%analyses_numbers=analyses_numbers(dokeep); % disable analyses that require SPM.mat files
                end
                [uanalyses_long,analyses_numbers]=conn_qaplotsexplore_translate(analyses);
                %uanalyses_long = regexprep(analyses,...
                %    {'^QA_NORM_(.*)','^QA_REG_functional','^QA_REG_(.*?)_?functional','^QA_REG_(.*?)_?structural','^QA_REG_(.*?)_?mni','^QA_COREG_(.*)','^QA_TIME_(.*)','^QA_TIMEART_(.*)','^QA_DENOISE_timeseries','^QA_DENOISE_QC-FC','^QA_DENOISE_scatterplot','^QA_DENOISE','^QA_SPM_design','^QA_SPM_contrasts'},...
                %    {'QA normalization: $1 data + outline of MNI TPM template','QA registration: functional data + structural overlay','QA registration: functional data + outline of ROI $1','QA registration: structural data + outline of ROI $1','QA registration: mni reference template + outline of ROI $1','QA realignment: $1 center-slice across multiple sessions/datasets','QA artifacts: $1 movie across all timepoints/acquisitions','QA artifacts: BOLD GS changes & subject motion timeseries with $1 movie','QA denoising: BOLD signal traces (carpetplot) before and after denoising + ART timeseries','QA denoising: distribution of QC-FC associations before and after denoising','QA denoising: scatterplot of functional correlations (FC) vs. distance (mm) before and after denoising','QA denoising: distribution of functional correlations (FC) before and after denoising','QA SPM design: review SPM first-level design matrix','QA SPM contrasts: review SPM first-level contrasts'});
                if strcmpi(option,'recreateset'), answ1=find(ismember(analyses_numbers,dlg.uanalysesnumber));
                elseif strcmpi(option,'recreateplot'), answ1=find(analyses_numbers==dlg.uanalysesnumber(dlg.ianalysis));
                elseif isempty(dlg.uanalyses), answ1=find(ismember(analyses_numbers,defaultset));
                else answ1=[];
                end
                
                [tstr,tidx]=conn_jobmanager('profiles');
                tnull=find(strcmp('Null profile',tstr));
                tlocal=find(strcmp('Background process (Unix,Mac)',tstr),1);
                tvalid=setdiff(1:numel(tstr),tnull);
                tstr=cellfun(@(x)sprintf('distributed processing (run on %s)',x),tstr,'uni',0);
                if 1, tvalid=tidx; if isunix&&~isempty(tlocal)&&~ismember(tlocal,tvalid), tvalid=[tvalid(:)' tlocal]; end % show default+local profiles
                elseif 1, tvalid=tidx; % show only default profile
                else tstr{tidx}=sprintf('<HTML><b>%s</b></HTML>',tstr{tidx}); % show all profiles
                end
                
                nl2covariates=[];
                [x,nl2covariates]=conn_module('get','l2covariates','^[^_]');
                il2covariates=find(cellfun('length',regexp(nl2covariates,'^QC_.*(QCOR_|MeanMotion|MeanGlobal|MeanGSchange|ValidScans)')));
                
                fh=figure('units','norm','position',[.4,.4,.4,.5],'menubar','none','numbertitle','off','name','compute QA plots','color','w');
                h1=uicontrol('units','norm','position',[.1,.90,.4,.05],'style','text','string','Plot type: ','horizontalalignment','left','fontweight','bold','backgroundcolor',get(fh,'color'));
                h2=uicontrol('units','norm','position',[.1,.65,.8,.25],'style','listbox','max',2,'string',uanalyses_long,'value',answ1,'tooltipstring','<HTML>Select type of plot(s) to create</HTML>');
                h4D=uicontrol('units','norm','position',[.1,.55,.35,.05],'style','checkbox','string','All subjects','value',1,'horizontalalignment','left','fontweight','bold','backgroundcolor',get(fh,'color'));
                h4d=uicontrol('units','norm','position',[.1,.25,.25,.30],'style','listbox','max',2,'string',arrayfun(@(n)sprintf('subject %d',n),1:CONN_x.Setup.nsubjects,'uni',0),'value',1:CONN_x.Setup.nsubjects,'tooltipstring','<HTML>Select subjet(s) to include in these plots<br/> - right-click for additional options</HTML>');
                h4A=uicontrol('units','norm','position',[.45,.55,.45,.05],'style','text','string','Plot options: ','horizontalalignment','left','fontweight','bold','backgroundcolor',get(fh,'color'));
                h4a=uicontrol('units','norm','position',[.45,.45,.45,.10],'style','listbox','max',1,'string',[{'primary dataset'},arrayfun(@(n)sprintf('secondary dataset #%d %s',n,regexprep(CONN_x.Setup.secondarydataset(n).label,'(.+)','($1)')),1:numel(CONN_x.Setup.secondarydataset),'uni',0)],'value',min(2,numel(CONN_x.Setup.secondarydataset)+1),'tooltipstring','<HTML>Select functional dataset(s) to include in functional data plots</HTML>');
                h4b=uicontrol('units','norm','position',[.45,.35,.45,.10],'style','listbox','max',2,'string',[CONN_x.Setup.rois.names(1:end-1), regexprep(CONN_x.Setup.rois.names(1:3),'^(.*)$','eroded $1')],'value',1,'tooltipstring','<HTML>Select ROI(s) to include in ROI data plots</HTML>');
                h4c=uicontrol('units','norm','position',[.45,.25,.45,.10],'style','listbox','max',2,'string',nl2covariates,'value',il2covariates,'tooltipstring','<HTML>Select QC variable(s) to include in QC/FC plots</HTML>');
                toptions=[{'local processing (run on this computer)'} tstr(tvalid)];
                if CONN_gui.isremote
                    info=conn_server('SSH_info');
                    if isfield(info,'host')&&~isempty(info.host), tnameserver=info.host;
                    else tnameserver='CONN server';
                    end
                    toptions=regexprep(toptions,'\<run on (this computer)?',['run on ',tnameserver,' ']);
                end
                h3=uicontrol('style','popupmenu','units','norm','position',[.1,.125,.8,.05],'string',toptions,'value',1);
                h6=uicontrol('units','norm','position',[.1,.025,.4,.08],'style','pushbutton','string','Start','callback','uiresume(gcbf)');
                h7=uicontrol('units','norm','position',[.5,.025,.4,.08],'style','pushbutton','string','Cancel','callback','close(gcbf)');
                hc1=uicontextmenu;uimenu(hc1,'Label','select group (2nd-level covariate)','callback',@(varargin)conn_qaplotsexplore_callbackgui('group',h4d)); set(h4d,'uicontextmenu',hc1);
                set([h2,h4D],'callback',@(varargin)conn_qaplotsexplore_callbackgui('update',h2,h4a,h4b,h4c,h4A,analyses_numbers,h4D,h4d));
                conn_qaplotsexplore_callbackgui('update',h2,h4a,h4b,h4c,h4A,analyses_numbers,h4D,h4d);
                uiwait(fh);
                if ishandle(fh),
                    procedures=analyses_numbers(get(h2,'value'));
                    if get(h4D,'value'), validsubjects=1:CONN_x.Setup.nsubjects;
                    else validsubjects=get(h4d,'value');
                    end
                    validsets=get(h4a,'value')-1;
                    validrois=get(h4b,'value');
                    nl2covariates=nl2covariates(get(h4c,'value'));
                    nalt=get(h3,'value');
                    nl1contrasts='?';
                    close(fh);
                    if nalt==1,
                        %conn_batch('subjects',validsubjects,'QA.foldername',fullfile(qafolder,tag),'QA.plots',procedures,'QA.rois',validrois,'QA.sets',validsets,'QA.l2covariates',nl2covariates,'QA.l1contrasts',nl1contrasts);
                        %conn_qaplots(conn_server('util_localfile',fullfile(qafolder,tag)),procedures,validsubjects,validrois,validsets,nl2covariates,nl1contrasts);
                        conn_process('qaplots',conn_server('util_localfile',fullfile(qafolder,tag)),procedures,validsubjects,validrois,validsets,nl2covariates,nl1contrasts);
                    else
                        if numel(validsubjects)>1
                            answer=inputdlg('Number of parallel jobs?','',1,{num2str(min(50,numel(validsubjects)))});
                            if isempty(answer), return; end
                            N=str2num(answer{1});
                        else N=1;
                        end
                        %conn_batch('subjects',validsubjects,'parallel.N',N,'parallel.profile',tvalid(nalt-1),'QA.foldername',fullfile(qafolder,tag),'QA.plots',procedures,'QA.rois',validrois,'QA.sets',validsets,'QA.l2covariates',nl2covariates,'QA.l1contrasts',nl1contrasts);
                        conn_jobmanager('setprofile',tvalid(nalt-1));
                        conn_jobmanager('submit','qaplots',validsubjects,N,[],conn_server('util_localfile',fullfile(qafolder,tag)),procedures,validsubjects,validrois,validsets,nl2covariates,nl1contrasts);
                    end
                    conn_qaplotsexplore_update([],[],'set');
                    conn_msgbox('Finished creating new plot','',2);
                end
                
                %                 if 1||isempty(answ), answ=listdlg('liststring',uanalyses_long,'selectionmode','multiple','initialvalue',answ,'promptstring','Select type of plot(s) to create:','ListSize',[600 250]); end
                %                 if isempty(answ), return; end
                %                 procedures=analyses_numbers(answ);
                %                 validsets=[];
                %                 if any(ismember(procedures,[2,7,8,9,10]))&&numel(CONN_x.Setup.secondarydataset)>0,
                %                     if numel(CONN_x.Setup.secondarydataset)==0, nalt=1;
                %                     else nalt=listdlg('liststring',arrayfun(@(n)sprintf('dataset %d',n),0:numel(CONN_x.Setup.secondarydataset),'uni',0),'selectionmode','multiple','initialvalue',min(2,numel(CONN_x.Setup.secondarydataset)+1),'promptstring',{'Select functional dataset(s)','to include in functional data plots:'},'ListSize',[300 200]);
                %                     end
                %                     if isempty(nalt), return; end
                %                     validsets=nalt-1;
                %                 end
                %                 validrois=[];
                %                 if any(ismember(procedures,[3:6])),
                %                     %nalt=listdlg('liststring',CONN_x.Setup.rois.names(1:end-1),'selectionmode','multiple','initialvalue',2,'promptstring',{'Select ROI(s)','to include in ROI data plots:'},'ListSize',[300 200]);
                %                     nalt=listdlg('liststring',[CONN_x.Setup.rois.names(1:end-1), regexprep(CONN_x.Setup.rois.names(1:3),'^(.*)$','eroded $1')],'selectionmode','multiple','initialvalue',1,'promptstring',{'Select ROI(s)','to include in ROI data plots:'},'ListSize',[300 200]);
                %                     if isempty(nalt), return; end
                %                     temp=numel(CONN_x.Setup.rois.names)-1;
                %                     nalt(nalt>temp)=-(nalt(nalt>temp)-temp);
                %                     validrois=nalt;
                %                 end
                %                 nl2covariates=[];
                %                 if any(ismember(procedures,[13,31]))
                %                     [x,nl2covariates]=conn_module('get','l2covariates','^[^_]');
                %                     if any(procedures==13), il2covariates=find(cellfun('length',regexp(nl2covariates,'^QC_.*(Motion|Global|GSchange|ValidScans)')));
                %                     else il2covariates=find(cellfun('length',regexp(nl2covariates,'^QC_')));
                %                     end
                %                     il2covariates=listdlg('liststring',nl2covariates,'selectionmode','multiple','initialvalue',il2covariates,'promptstring',{'Select QC variable(s)','to include in QC-FC analyses:'},'ListSize',[300 200]);
                %                     if isempty(il2covariates), return; end
                %                     nl2covariates=nl2covariates(il2covariates);
                %                 end
                %                 nl1contrasts='?';
                %                 if CONN_x.Setup.nsubjects==1, nalt=1;
                %                 else  nalt=listdlg('liststring',arrayfun(@(n)sprintf('subject %d',n),1:CONN_x.Setup.nsubjects,'uni',0),'selectionmode','multiple','initialvalue',1:CONN_x.Setup.nsubjects,'promptstring',{'Select subject(s)','to include in these plots:'},'ListSize',[300 200]);
                %                 end
                %                 if isempty(nalt), return; end
                %                 validsubjects=nalt;
                %
                %                 [tstr,tidx]=conn_jobmanager('profiles');
                %                 tnull=find(strcmp('Null profile',tstr));
                %                 tlocal=find(strcmp('Background process (Unix,Mac)',tstr),1);
                %                 tvalid=setdiff(1:numel(tstr),tnull);
                %                 tstr=cellfun(@(x)sprintf('distributed processing (run on %s)',x),tstr,'uni',0);
                %                 if 1, tvalid=tidx; if isunix&&~isempty(tlocal)&&~ismember(tlocal,tvalid), tvalid=[tvalid(:)' tlocal]; end % show default+local profiles
                %                 elseif 1, tvalid=tidx; % show only default profile
                %                 else tstr{tidx}=sprintf('<HTML><b>%s</b></HTML>',tstr{tidx}); % show all profiles
                %                 end
                %                 nalt=listdlg('liststring',[{'local processing (run on this computer)'} tstr(tvalid)],'selectionmode','single','initialvalue',1,'promptstring','','ListSize',[300 200]);
                %                 if isempty(nalt), return; end
                %conn_qaplots(fullfile(qafolder,tag),procedures,validsubjects,validrois,validsets,nl2covariates,nl1contrasts);
                figure(dlg.handles.hfig);
            case 'displayannotation'
                if isfield(dlg,'dataIDXsubjects')
                    in=dlg.dataIDXsubjects;
                    if numel(in)==1
                        txt=dlg.txtA(in);
                    elseif numel(in)>1
                        txt=arrayfun(@(n)sprintf('%s %s %s: %s',dlg.usubjects{dlg.isubjects(dlg.dataIDXplots(n))},dlg.usessions{dlg.isessions(dlg.dataIDXplots(n))},dlg.umeasures{dlg.imeasures(dlg.dataIDXplots(n))},sprintf('%s ',dlg.txtA{n}{:})),in,'uni',0);
                    end
                    for n1=1:numel(txt),
                        if ischar(txt{n1}), txt{n1}=cellstr(txt{n1}); end
                        for n2=1:numel(txt{n1})
                            if ~isempty(txt{n1}{n2})&&isempty(regexp(txt{n1}{n2},'\[empty\]$')), conn_disp('fprintf','%s\n',txt{n1}{n2}); end;
                        end
                    end
                end
            case 'details'
                filename=dlg.filethis;
                if isempty(filename)||~conn_existfile(filename), conn_msgbox(sprintf('Data file %s not found',filename),'Details not available',2);
                else
                    conn_bookmark('open',filename);
                    %load(filename,'state');
                    %conn_slice_display(state);
                end
                return;
                
            case 'set'
                if numel(varargin)>=1&&~isempty(varargin{1}), dlg.iset=varargin{1}; set(dlg.handles.set,'value',dlg.iset);
                else dlg.iset=get(dlg.handles.set,'value');
                end
                if isempty(dlg.sets),
                    set(dlg.handles.settxt,'visible','on');
                    set([dlg.handles.set dlg.handles.deleteset dlg.handles.recomputeset dlg.handles.printset],'visible','off');
                    set(dlg.handles.addnewset,'fontweight','bold');
                    set([dlg.handles.selectall1 dlg.handles.selectall2 dlg.handles.selectall3 dlg.handles.subjects dlg.handles.sessions dlg.handles.measures dlg.handles.showdiff dlg.handles.showannot dlg.handles.invertim dlg.handles.analysis dlg.handles.analysistxt dlg.handles.addnewplot dlg.handles.deleteplot dlg.handles.recomputeplot dlg.handles.printplot dlg.handles.text1 dlg.handles.text2 dlg.handles.textoptions dlg.handles.details],'visible','off');
                    delete(dlg.handles.hax(ishandle(dlg.handles.hax)));
                    delete(dlg.handles.han(ishandle(dlg.handles.han)));
                    dlg.dispsize=[];
                    return
                else
                    set(dlg.handles.settxt,'visible','off');
                    set([dlg.handles.set dlg.handles.deleteset dlg.handles.recomputeset dlg.handles.printset],'visible','on');
                    set(dlg.handles.addnewset,'fontweight','normal');
                end
                qafiles=conn_dir(fullfile(qafolder,dlg.sets{dlg.iset},'QA_*.mat'));%,'-R'); % note: remove -R to allow search of recursive subfolders in "plots"
                if isempty(qafiles), qanames={};
                else qanames=cellstr(qafiles);
                end
                jpgok=conn_existfile(conn_prepend('',qanames,'.jpg'))|cellfun('length',regexp(qanames,'QA_DENOISE\.|QA_DENOISE_QC-FC\.|QA_DENOISE_scatterplot\.|QA_DENOISE_QC-FC_scatterplot\.|QA_COV\.'))>0;
                txtok=conn_existfile(conn_prepend('',qanames,'.txt'));
                qanames=qanames(jpgok);
                txtok=txtok(jpgok);
                dlg.files=qanames;
                dlg.files_jpg=conn_prepend('',qanames,'.jpg');
                dlg.files_txt=conn_prepend('',qanames,'.txt');
                try
                    if ~all(txtok),conn_fileutils('emptyfile',dlg.files_txt(~txtok)); end
                    %if ~all(txtok),cellfun(@(s)conn_fileutils('emptyfile',s),dlg.files_txt(~txtok),'uni',0); end
                    %if ~all(txtok),cellfun(@(s)fclose(fopen(s,'wt')),dlg.files_txt(~txtok),'uni',0); end
                catch
                    conn_disp('warning: unable to create annotation files; please check QA folder write-permissions');
                end
                if isempty(qanames)
                    qanames_parts1={};qanames_parts2={};qanames_parts3={};qanames_parts4={};
                    qanames_parts2isnumber=true;
                    qanames_parts3isnumber=true;
                else
                    [qafolders,qanames]=cellfun(@fileparts,qanames,'uni',0);
                    qanames=regexp(qanames,'\.','split');
                    
                    qanames_parts=repmat({''},numel(qanames),3);
                    for n=1:numel(qanames), i=1:min(3,numel(qanames{n})); qanames_parts(n,i)=qanames{n}(i); end % analyses / subject /session
                    qanames_parts1=qanames_parts(:,1);
                    %qanames_parts2isnumber=true;
                    %qanames_parts2=str2double(regexprep(qanames_parts(:,2),'^subject',''));
                    qanames_parts2=regexprep(qanames_parts(:,2),'^subject','');
                    temp=str2double(qanames_parts2);
                    qanames_parts2isnumber=any(~isnan(temp));
                    if qanames_parts2isnumber, qanames_parts2=temp; end
                    qanames_parts3isnumber=true;
                    qanames_parts3=str2double(regexprep(qanames_parts(:,3),'^session',''));
                    qanames_parts4valid=cellfun('length',regexp(qanames_parts(:,2),'^measure'))>0;
                    qanames_parts4=regexprep(qanames_parts(:,2),'^measure','');
                    if numel(qafolders)>1&any(any(diff(char(qafolders),1,1))),
                        tidx=find(any(diff(char(qafolders),1,1)),1);
                        tidx=max([1,find(ismember(qafolders{1}(1:tidx),'\/'),1,'last')+1]);
                        if max(qanames_parts2)==1 % folders are subjects
                            these=find(qanames_parts2==1);
                            tqanames_parts2=repmat({'---'},size(qanames_parts2));
                            tqanames_parts2(these)=regexprep(cellfun(@(b)b(min(numel(b),tidx):end),qafolders(these),'uni',0),'[\/\\].*','');
                            qanames_parts2isnumber=false;
                            qafolders_str=regexp(tqanames_parts2(these),'\d+','match');
                            if all(cellfun('length',qafolders_str)==1)
                                qafolders_n=str2double([qafolders_str{:}]);
                                if all(~isnan(qafolders_n))
                                    qanames_parts2(these)=qafolders_n;
                                    qanames_parts2isnumber=true;
                                end
                            end
                            if ~qanames_parts2isnumber, qanames_parts2=tqanames_parts2; end
                            %                             qafolders_str=regexp(qafolders(these),'\d+','match','once');
                            %                             qafolders_n=str2double(qafolders_str);
                            %                             if all(~isnan(qafolders_n)), eqafolders=cellfun(@(a,b)[a,b],qafolders_str,qanames_parts1(these),'uni',0);
                            %                             else eqafolders=cellfun(@(a,b)[a,b],qafolders(these),qanames_parts1(these),'uni',0);
                            %                             end
                            %                             try, [i1,i2,i3]=unique(eqafolders,'last');
                            %                             catch, [i1,i2,i3]=unique(eqafolders);
                            %                             end
                            %                             qanames_parts2(these)=nan;
                            %                             if ~all(~isnan(qafolders_n)), qanames_parts2(these(i2))=1:numel(i1); else qanames_parts2(these(i2))=qafolders_n(i2); end % keeps last
                            %                             %if ~all(~isnan(qafolders_n)), qanames_parts2(these)=i3; else qanames_parts2(these)=qafolders_n; end % keeps all (aggregate)
                            %                             %qanames_parts2(these)=1:numel(these); % keeps all (expand)
                            %                             for n1=reshape(unique(qanames_parts2(these(~isnan(qanames_parts2(these))))),1,[]),temp=unique(qafolders(these(qanames_parts2(these)==n1))); conn_disp('fprintf','note: subject %d = %s\n',n1,sprintf('%s ',temp{:})); end
                            %                             if any(isnan(qanames_parts2(these))), temp=unique(qafolders(these(isnan(qanames_parts2(these))))); conn_disp('fprintf','note: subject --- = %s\n',sprintf('%s ',temp{:})); end
                        else                      % folders are plots
                            qanames_parts1=cellfun(@(a,b)sprintf('%s (%s)',a,regexprep(b(min(numel(b),tidx):end),'[\/\\].*','')),qanames_parts1,qafolders,'uni',0);
                        end
                    end
                    if qanames_parts2isnumber, qanames_parts2(isnan(qanames_parts2))=0; end
                    if qanames_parts3isnumber, qanames_parts3(isnan(qanames_parts3))=0; end
                    [qanames_parts4{~qanames_parts4valid}]=deal('');
                end
                [dlg.uanalyses,nill,dlg.ianalyses]=unique(qanames_parts1);
                [dlg.usubjects,nill,dlg.isubjects]=unique(qanames_parts2);
                [dlg.usessions,nill,dlg.isessions]=unique(qanames_parts3);
                [dlg.umeasures,nill,dlg.imeasures]=unique(qanames_parts4);
                allids=[dlg.ianalyses, dlg.isubjects, dlg.isessions, dlg.imeasures]; % show only last plot if multiple exist for each subject/session/plot/measure
                try, [nill,valid]=unique(allids,'rows','last');
                catch, [nill,valid]=unique(allids,'rows');
                end
                repeated=setdiff(1:size(allids,1),valid);
                dlg.ianalyses(repeated)=nan;
                %dlg.uanalysestype=ones(size(dlg.uanalyses)); % QA_NORM/QA_REG/QA_COREG/QA_TIME/QA_TIMEART
                %dlg.uanalysestype(cellfun('length',regexp(dlg.uanalyses,'^QA_DENOISE$'))>0)=2; %QA_DENOISE_FC
                %dlg.uanalysestype(cellfun('length',regexp(dlg.uanalyses,'^QA_DENOISE_QC-FC$'))>0)=3; %QA_DENOISE_QC-FC
                %dlg.uanalysestype(cellfun('length',regexp(dlg.uanalyses,'^QA_DENOISE_scatterplot$'))>0)=4; %QA_DENOISE_scatterplot
                %dlg.uanalysestype(cellfun('length',regexp(dlg.uanalyses,'^QA_COV$'))>0)=5; %QA_COV
                if qanames_parts2isnumber, dlg.usubjects=regexprep(arrayfun(@(n)sprintf('subject %d',n),dlg.usubjects,'uni',0),'^subject 0$','---'); end
                if qanames_parts3isnumber, dlg.usessions=regexprep(arrayfun(@(n)sprintf('session %d',n),dlg.usessions,'uni',0),'^session 0$','---'); end
                %dlg.umeasures=regexprep(arrayfun(@(n)sprintf('measure %d',n),dlg.umeasures,'uni',0),'^measure 0$','---');
                dlg.ianalysis=0;
                [dlg.uanalyses_long,dlg.uanalysesnumber,dlg.uanalysestype]=conn_qaplotsexplore_translate(dlg.uanalyses);
                if dlg.forceinitdenoise,
                    dlg.ianalysis=find(dlg.uanalysestype==2,1);
                    %if ~dlg.createdenoise&&~isempty(dlg.ianalysis)
                    %    answ=conn_questdlg({'Overwrite existing denoising plot?'},'','Yes','No','Yes');
                    %    if strcmp(answ,'Yes'), dlg.createdenoise=true; end
                    %end
                    if dlg.createdenoise||isempty(dlg.ianalysis)
                        conn_qaplots(fullfile(qafolder,dlg.sets{dlg.iset}),11);
                        conn_qaplotsexplore_update([],[],'set');
                        return;
                    end
                    if isempty(dlg.ianalysis), dlg.ianalysis=find(dlg.uanalysestype>1,1); end
                    dlg.forceinitdenoise=false;
                    %else
                    %    temp=find(dlg.uanalysestype>1,1); % note: tries loading denoising by default (faster)
                    %    if ~isempty(temp), dlg.ianalysis=temp; end
                end
                if numel(dlg.uanalyses)==1, dlg.ianalysis=1; end
                %if ~isfield(dlg,'ianalysis')||isempty(dlg.ianalysis)||dlg.ianalysis<1||dlg.ianalysis>numel(dlg.uanalyses), dlg.ianalysis=1; end
                %dlg.uanalyses_long = regexprep(dlg.uanalyses,...
                %    {'^QA_NORM_(.*)','^QA_REG_functional','^QA_REG_(.*?)_?functional','^QA_REG_(.*?)_?structural','^QA_REG_(.*?)_?mni','^QA_COREG_(.*)','^QA_TIME_(.*)','^QA_TIMEART_(.*)','^QA_DENOISE_timeseries','^QA_DENOISE_QC-FC','^QA_DENOISE_scatterplot','^QA_DENOISE','^QA_SPM_design','^QA_SPM_contrasts'},...
                %    {'QA normalization: $1 data + outline of MNI TPM template','QA registration: functional data + structural overlay','QA registration: functional data + outline of ROI $1','QA registration: structural data + outline of ROI $1','QA registration: mni reference template + outline of ROI $1','QA realignment: $1 center-slice across multiple sessions/datasets','QA artifacts: $1 movie across all timepoints/acquisitions','QA artifacts: BOLD GS changes & subject motion timeseries with $1 movie','QA denoising: BOLD signal traces (carpetplot) before and after denoising + ART timeseries','QA denoising: distribution of QC-FC associations before and after denoising','QA denoising: scatterplot of functional correlations (FC) vs. distance (mm) before and after denoising','QA denoising: distribution of functional correlations (FC) before and after denoising','QA SPM design: review SPM first-level design matrix','QA SPM contrasts: review SPM first-level contrasts'});
                dlg.uanalyses_long=arrayfun(@(n,m)sprintf('%s (%d)',dlg.uanalyses_long{n},m),1:numel(dlg.uanalyses_long),accumarray(dlg.ianalyses(dlg.ianalyses>0),1)','uni',0);
                set(dlg.handles.analysis,'string',[{'<HTML><i>choose an existing QA plot to display</i></HTML>'},dlg.uanalyses_long],'value',dlg.ianalysis+1);
                conn_qaplotsexplore_update([],[],'plot');
                
            case 'plot'
                if numel(varargin)>=1, dlg.ianalysis=varargin{1}; set(dlg.handles.analysis,'value',dlg.ianalysis+1);
                else dlg.ianalysis=get(dlg.handles.analysis,'value')-1;
                end
                if isempty(dlg.sets)||isempty(dlg.uanalyses),
                    set([dlg.handles.selectall1 dlg.handles.selectall2 dlg.handles.selectall3 dlg.handles.subjects dlg.handles.sessions dlg.handles.measures dlg.handles.showdiff dlg.handles.showannot dlg.handles.invertim dlg.handles.analysis dlg.handles.text1 dlg.handles.text2 dlg.handles.textoptions dlg.handles.details],'visible','off');
                    delete(dlg.handles.hax(ishandle(dlg.handles.hax)));
                    delete(dlg.handles.han(ishandle(dlg.handles.han)));
                    set([dlg.handles.analysis dlg.handles.deleteplot dlg.handles.recomputeplot dlg.handles.printplot dlg.handles.recomputeset dlg.handles.printset],'visible','off')
                    set(dlg.handles.analysistxt,'visible','on');
                    set(dlg.handles.addnewplot,'fontweight','bold','visible','on');
                    dlg.dispsize=[];
                    return;
                else
                    set([dlg.handles.selectall1 dlg.handles.selectall2 dlg.handles.selectall3 dlg.handles.subjects dlg.handles.sessions dlg.handles.measures dlg.handles.analysis dlg.handles.deleteplot dlg.handles.recomputeplot dlg.handles.printplot dlg.handles.recomputeset dlg.handles.printset],'visible','on')
                    set(dlg.handles.analysistxt,'visible','off');
                    set(dlg.handles.addnewplot,'fontweight','normal','visible','on');
                end
                if isequal(dlg.ianalysis,0)
                    set([dlg.handles.selectall1 dlg.handles.selectall2 dlg.handles.selectall3 dlg.handles.subjects dlg.handles.sessions dlg.handles.measures dlg.handles.showdiff dlg.handles.showannot dlg.handles.invertim dlg.handles.text1 dlg.handles.text2 dlg.handles.textoptions,dlg.handles.details,dlg.handles.deleteplot dlg.handles.recomputeplot dlg.handles.printplot],'visible','off');
                    delete(dlg.handles.hax(ishandle(dlg.handles.hax)));
                    delete(dlg.handles.han(ishandle(dlg.handles.han)));
                elseif dlg.uanalysestype(dlg.ianalysis)==1 % QA_NORM/QA_REG/QA_COREG/QA_TIME/QA_TIMEART
                    set([dlg.handles.subjects dlg.handles.sessions dlg.handles.measures dlg.handles.selectall1 dlg.handles.selectall2 dlg.handles.selectall3],'visible','off');
                    set([dlg.handles.showdiff dlg.handles.showannot dlg.handles.invertim],'visible','off');
                    set(dlg.handles.hfig,'pointer','watch');
                    in=find(ismember(dlg.ianalyses,dlg.ianalysis));% & ismember(dlg.isubjects, dlg.isubject) & ismember(dlg.isessions, dlg.isession);
                    %ht=conn_msgbox(sprintf('Loading %d plots. Please wait...',numel(in)),'');
                    dlg.dataIDXplots=in;
                    [dlg.dataM,dlg.dataS,dlg.dataA,dlg.txtA,dlg.dataD,dataDmax,dataDidx]=conn_qaplotsexplore_readfiles(dlg.files_jpg,dlg.files_txt,in,false,[],true);
                    %if ishandle(ht),delete(ht); end
                    if ~ishandle(dlg.handles.hfig), return; end
                    dlg.usubjects_shown=unique(dlg.isubjects(in));
                    set(dlg.handles.subjects,'string',dlg.usubjects(dlg.usubjects_shown),'value',unique(max(1,min(numel(dlg.usubjects_shown),get(dlg.handles.subjects,'value')))));
                    dlg.usessions_shown=unique(dlg.isessions(in));
                    set(dlg.handles.sessions,'string',dlg.usessions(dlg.usessions_shown),'value',unique(max(1,min(numel(dlg.usessions_shown),get(dlg.handles.sessions,'value')))));
                    dlg.umeasures_shown=unique(dlg.imeasures(in));
                    set(dlg.handles.measures,'string',dlg.umeasures(dlg.umeasures_shown),'value',unique(max(1,min(numel(dlg.umeasures_shown),get(dlg.handles.measures,'value')))));
                    set([dlg.handles.showdiff dlg.handles.showannot dlg.handles.invertim dlg.handles.analysis],'visible','on');
                    conn_qaplotsexplore_update([],[],'subjects');
                    set(dlg.handles.hfig,'pointer','arrow');
                elseif dlg.uanalysestype(dlg.ianalysis)>1 %QA_DENOISE/QA_COV
                    set([dlg.handles.subjects dlg.handles.sessions dlg.handles.measures dlg.handles.selectall1 dlg.handles.selectall2 dlg.handles.selectall3],'visible','off');
                    set([dlg.handles.showdiff dlg.handles.showannot dlg.handles.invertim dlg.handles.text1 dlg.handles.text2 dlg.handles.textoptions dlg.handles.details],'visible','off');
                    set(dlg.handles.hfig,'pointer','watch');
                    in=find(ismember(dlg.ianalyses,dlg.ianalysis));% & ismember(dlg.isubjects, dlg.isubject) & ismember(dlg.isessions, dlg.isession);
                    %ht=conn_msgbox(sprintf('Loading %d plots. Please wait...',numel(in)),'');
                    dlg.dataIDXplots=in;
                    [dlg.dataA,dlg.labelA,dlg.infoA,dlg.lineA,dlg.txtA,miny,maxy]=conn_qaplotsexplore_readplots(dlg.files,dlg.files_txt,in, dlg.usubjects(dlg.isubjects(in)));
                    
                    %miny1=inf;miny2=inf;miny3=inf;maxy1=0;maxy2=0;maxy3=0;
                    %for n=1:numel(dlg.dataA), miny1=min(miny1,min(dlg.dataA{n}{1})); miny2=min(miny2,min(dlg.dataA{n}{2})); miny3=min(miny3,min(dlg.dataA{n}{3})); maxy1=max(maxy1,max(dlg.dataA{n}{1})); maxy2=max(maxy2,max(dlg.dataA{n}{2})); maxy3=max(maxy3,max(dlg.dataA{n}{3})); end
                    dlg.plotminmax=[miny; maxy]; %[miny1 miny2 miny3; maxy1 maxy2 maxy3];
                    dlg.plothistinfo=[0 maxy(2)*1.1 (maxy(2)+maxy(min(numel(maxy),3)))*1.1 max(maxy(2),maxy(min(numel(maxy),3))) maxy(1)*1.1];
                    %if ishandle(ht),delete(ht); end
                    if ~ishandle(dlg.handles.hfig), return; end
                    dlg.usubjects_shown=unique(dlg.isubjects(in));
                    set(dlg.handles.subjects,'string',dlg.usubjects(dlg.usubjects_shown),'value',unique(max(1,min(numel(dlg.usubjects_shown),get(dlg.handles.subjects,'value')))));
                    dlg.usessions_shown=unique(dlg.isessions(in));
                    set(dlg.handles.sessions,'string',dlg.usessions(dlg.usessions_shown),'value',unique(max(1,min(numel(dlg.usessions_shown),get(dlg.handles.sessions,'value')))));
                    dlg.umeasures_shown=unique(dlg.imeasures(in));
                    set(dlg.handles.measures,'string',dlg.umeasures(dlg.umeasures_shown),'value',unique(max(1,min(numel(dlg.umeasures_shown),get(dlg.handles.measures,'value')))));
                    set(dlg.handles.showannot,'visible','on');
                    conn_qaplotsexplore_update([],[],'subjects'); 
                    set(dlg.handles.hfig,'pointer','arrow');
                end
                
            case {'subjects','sessions','measures','selectannotation','refresh'}
                if isempty(dlg.sets)||isempty(dlg.uanalyses), return; end
                if isequal(dlg.ianalysis,0)
                    delete(dlg.handles.hax(ishandle(dlg.handles.hax)));
                    delete(dlg.handles.han(ishandle(dlg.handles.han)));
                    return;
                elseif isempty(dlg.usubjects_shown)||isempty(dlg.usessions_shown)||isempty(dlg.umeasures_shown),
                    delete(dlg.handles.hax(ishandle(dlg.handles.hax)));
                    delete(dlg.handles.han(ishandle(dlg.handles.han)));
                    set([dlg.handles.text1 dlg.handles.text2 dlg.handles.textoptions dlg.handles.details],'visible','off');
                    return;
                end
                isvisible=[numel(dlg.usubjects_shown)>=2 numel(dlg.usessions_shown)>=2 numel(dlg.umeasures_shown)>=2];
                %.73 .57 - .51 .41 - .35 .25
                set([dlg.handles.subjects dlg.handles.sessions dlg.handles.measures dlg.handles.selectall1 dlg.handles.selectall2 dlg.handles.selectall3],'visible','off');
                if isvisible(1), set(dlg.handles.subjects,'position',[.035,.57*isvisible(2)+.41*~isvisible(2)*isvisible(3)+.25*~isvisible(2)*~isvisible(3),.1,.16*isvisible(2)+.32*~isvisible(2)*isvisible(3)+.48*~isvisible(2)*~isvisible(3)]); end
                if isvisible(2), set(dlg.handles.sessions,'position',[.035,.41*isvisible(3)+.25*~isvisible(3),.1,.10*isvisible(3)+.26*~isvisible(3)+~.22*isvisible(1)]); end
                if isvisible(3), set(dlg.handles.measures,'position',[.035,.25,.1,.10*isvisible(2)+.26*~isvisible(2)*isvisible(1)+.48*~isvisible(2)*~isvisible(1)]); end
                set(dlg.handles.selectall1,'position',get(dlg.handles.subjects,'position').*[1 1 1 0]+[0 -.05 0 .05]);
                set(dlg.handles.selectall2,'position',get(dlg.handles.sessions,'position').*[1 1 1 0]+[0 -.05 0 .05]);
                set(dlg.handles.selectall3,'position',get(dlg.handles.measures,'position').*[1 1 1 0]+[0 -.05 0 .05]);
                if isvisible(1), set([dlg.handles.subjects dlg.handles.selectall1],'visible','on'); end
                if isvisible(2), set([dlg.handles.sessions dlg.handles.selectall2],'visible','on'); end
                if isvisible(3), set([dlg.handles.measures dlg.handles.selectall3],'visible','on'); end
                if strcmp(lower(option),'selectannotation')
                    n=get(dlg.handles.han(end),'value');
                    subjects=dlg.isubjects(dlg.dataIDXplots(dlg.dataIDXsubjects(n)));
                    sessions=dlg.isessions(dlg.dataIDXplots(dlg.dataIDXsubjects(n)));
                    measures=dlg.imeasures(dlg.dataIDXplots(dlg.dataIDXsubjects(n)));
                    set(dlg.handles.subjects,'value',find(dlg.usubjects_shown==subjects));
                    set(dlg.handles.sessions,'value',find(dlg.usessions_shown==sessions));
                    set(dlg.handles.measures,'value',find(dlg.umeasures_shown==measures));
                else
                    subjects=get(dlg.handles.subjects,'value');
                    sessions=get(dlg.handles.sessions,'value');
                    measures=get(dlg.handles.measures,'value');
                    if isempty(subjects)||any(subjects>numel(dlg.usubjects_shown)), subjects=1:numel(dlg.usubjects_shown); set(dlg.handles.subjects,'value',subjects); end
                    if isempty(sessions)||any(sessions>numel(dlg.usessions_shown)), sessions=1:numel(dlg.usessions_shown); set(dlg.handles.sessions,'value',sessions); end
                    if isempty(measures)||any(measures>numel(dlg.umeasures_shown)), measures=1:numel(dlg.umeasures_shown); set(dlg.handles.measures,'value',measures); end
                    subjects=dlg.usubjects_shown(subjects);
                    sessions=dlg.usessions_shown(sessions);
                    measures=dlg.umeasures_shown(measures);
                end
                in=find(ismember(dlg.isubjects(dlg.dataIDXplots),subjects)&ismember(dlg.isessions(dlg.dataIDXplots),sessions)&ismember(dlg.imeasures(dlg.dataIDXplots),measures));
                dlg.dataIDXsubjects=in;
                switch(dlg.uanalysestype(dlg.ianalysis))
                    case 1,
                        if numel(in)>1, set(dlg.handles.text1,'string','computing. please wait...','visible','on');set(dlg.handles.textoptions,'visible','off');set(dlg.handles.text2,'visible','off'); drawnow; end
                        if numel(in)==1,%if size(dlg.dataA,4)>1,
                            set(dlg.handles.showdiff,'visible','on');
                            showdiff=get(dlg.handles.showdiff,'value');
                        else
                            set(dlg.handles.showdiff,'visible','off');
                            showdiff=false;
                        end
                        val=(numel(in)>1)*dlg.showavg + (numel(in)<=1);
                        delete(dlg.handles.hax(ishandle(dlg.handles.hax)));
                        pos=[.20 .10 .75 .65];
                        if get(dlg.handles.showannot,'value'), pos=[pos(1)+.225 pos(2) pos(3)-.20 pos(4)]; end
                        dlg.handles.hax=axes('units','norm','position',pos,'color',figcolor,'visible','off','parent',dlg.handles.hfig);
                        switch(val)
                            case {1,2,3,4},
                                if showdiff&&~isempty(dlg.dataD), data=dlg.dataD;
                                else data=dlg.dataA;
                                end
                                if iscell(data), sd4=numel(data);
                                else sd4=size(data,4);
                                end
                                if val==1,
                                    if ~showdiff&&isequal(in(:)',1:sd4), dlg=conn_qaplotsexplore_updatemean(dlg); data=dlg.dataM;
                                    elseif showdiff&&isequal(in(:)',1:sd4),
                                        dlg=conn_qaplotsexplore_updatemean(dlg); 
                                        dlg.dlgD=abs(conn_qaplotsexplore_getdata(data,in)-repmat(dlg.dataM,[1,1,1,numel(in)]));
                                        data=mean(dlg.dlgD,4);
                                    elseif showdiff, dlg=conn_qaplotsexplore_updatemean(dlg); data=mean(abs(conn_qaplotsexplore_getdata(data,in)-repmat(dlg.dataM,[1,1,1,numel(in)])),4);
                                    else data=mean(conn_qaplotsexplore_getdata(data,in),4);
                                    end
                                    dlg.dispsize=[size(data,2) size(data,1)];
                                elseif val==2,
                                    if ~showdiff&&isequal(in(:)',1:sd4), dlg=conn_qaplotsexplore_updatemean(dlg); data=dlg.dataS;
                                    elseif showdiff&&isequal(in(:)',1:sd4),
                                        dlg=conn_qaplotsexplore_updatemean(dlg); 
                                        dlg.dlgD=abs(conn_qaplotsexplore_getdata(data,in)-repmat(dlg.dataM,[1,1,1,numel(in)]));
                                        data=std(dlg.dlgD,1,4);
                                    elseif showdiff, dlg=conn_qaplotsexplore_updatemean(dlg); data=std(abs(conn_qaplotsexplore_getdata(data,in)-repmat(dlg.dataM,[1,1,1,numel(in)])),1,4);
                                    else data=std(conn_qaplotsexplore_getdata(data,in),1,4);
                                    end
                                    data=sqrt(sum(data.^2,3));
                                    data=data/max(data(:));
                                    dlg.dispsize=[size(data,2) size(data,1)];
                                elseif val==3,
                                    if showdiff&&isequal(in(:)',1:sd4),
                                        dlg=conn_qaplotsexplore_updatemean(dlg); 
                                        dlg.dlgD=abs(conn_qaplotsexplore_getdata(data,in)-repmat(dlg.dataM,[1,1,1,numel(in)]));
                                        data=dlg.dlgD;
                                    elseif showdiff, dlg=conn_qaplotsexplore_updatemean(dlg); data=abs(conn_qaplotsexplore_getdata(data,in)-repmat(dlg.dataM,[1,1,1,numel(in)]));
                                    else data=conn_qaplotsexplore_getdata(data,in);
                                    end
                                elseif val==4,
                                    if showdiff, dlg=conn_qaplotsexplore_updatemean(dlg); data=abs(conn_qaplotsexplore_getdata(data,in(end))-dlg.dataM);
                                    else data=conn_qaplotsexplore_getdata(data,in(end));
                                    end
                                end
                                if get(dlg.handles.invertim,'value')>0,
                                    if size(data,3)==3&&min(data(:))>=0&&max(data(:))<=1,
                                        maxdata=mode(round(reshape(permute(data,[1,2,4,3]),[],3)*100))/100;
                                        if maxdata<.5, data=1-data; maxdata=1-maxdata; end
                                        data=max(0,min(1, data.*repmat(shiftdim(figcolor./maxdata,-1),[size(data,1),size(data,2),1,size(data,4)]) ));
                                    elseif size(data,3)==1, data=-data;
                                    end
                                end
                                [data,dlg.dispsize]=conn_menu_montage(dlg.handles.hax,data);
                                masknan=find(any(isnan(data),3)); if ~isempty(masknan), mostnan=data(1,1,:); for n=1:size(data,3), data(masknan+(n-1)*size(data,1)*size(data,2))=mostnan(n); end; end
                                cla(dlg.handles.hax);
                                him=imagesc(data,'parent',dlg.handles.hax);
                                axis(dlg.handles.hax,'equal');
                                set(dlg.handles.hax,'ydir','reverse','visible','off');
                                %if numel(in)==1, set(him,'buttondownfcn',{@conn_qaplotsexplore_update,'details'} );
                                %else set(him,'buttondownfcn','conn_disp(''select individual subject/session first'')');
                                %end
                            case 5, % placeholder
                                data=reshape(dlg.dataD(:,:,:,in),[],numel(in));
                                cla(dlg.handles.hax);
                                for n=1:size(data,2),
                                    [b,a]=hist(log10(data(data(:,n)>0,n)),linspace(-3,1,100));
                                    plot(a,b,'parent',dlg.handles.hax);
                                    hold(dlg.handles.hax,'on');
                                end
                                hold(dlg.handles.hax,'off');
                        end
                        if numel(in)==1,
                            if showdiff, str='diff-map of '; else str=''; end
                            [tpath,tstr]=fileparts(dlg.files_jpg{dlg.dataIDXplots(in)});
                            set(dlg.handles.text1,'string',[str,tstr],'visible','on');
                            %set(dlg.handles.text2,'string',tpath,'visible','on');
                            set(dlg.handles.text2,'string',tpath(max(numel(tpath)-100,numel(qafolder)+1):end),'visible','on');
                            set(dlg.handles.textoptions,'visible','off');
                            dlg.filethis=dlg.files{dlg.dataIDXplots(in)};
                            set([dlg.handles.details],'visible','on');
                            %if conn_existfile(fullfile(qafolder,dlg.sets{dlg.iset},conn_prepend('',dlg.filethis,'.mat'))), set(dlg.handles.details,'visible','on'); else set(dlg.handles.details,'visible','on'); end
                        else
                            if showdiff, str='z-maps'; else str='images'; end
                            if dlg.showavg==4
                                [tpath,tstr]=fileparts(dlg.files_jpg{dlg.dataIDXplots(in(end))});
                                set(dlg.handles.text1,'string',[str,tstr],'visible','on');
                                set(dlg.handles.text2,'string',tpath(max(numel(tpath)-100,numel(qafolder)+1):end),'visible','on');
                            else
                                set([dlg.handles.text1 dlg.handles.text2],'visible','off');
                            end
                            set(dlg.handles.textoptions,'visible','on','value',dlg.showavg,'string',{sprintf('average of %d %s',numel(in),str),sprintf('variability of %d %s',numel(in),str),sprintf('Montage of %d %s',numel(in),str)}); %,sprintf('Last of %d %s',numel(in),str)});
                            set([dlg.handles.details],'visible','off');
                            dlg.filethis='';
                        end
                        dlg.followmouse=get(dlg.handles.showdiff,'value');
                        
                    case {2,3,4} % QA_DENOISE
                        delete(dlg.handles.hax(ishandle(dlg.handles.hax)));
                        pos=[.30 .175 .55 .575];
                        if get(dlg.handles.showannot,'value'), pos=[pos(1)+.225 pos(2) pos(3)-.20 pos(4)]; end
                        dlg.handles.hax=axes('units','norm','position',pos,'color',figcolor,'visible','off','parent',dlg.handles.hfig);
                        dlg.results_patch=dlg.dataA(in); %%%
                        dlg.results_label=dlg.labelA(in);
                        dlg.results_info=dlg.infoA(in);
                        dlg.results_line=dlg.lineA(in);
                        dlg.handles.resultspatch=[];
                        dlg.handles.resultsline=[];
                        hold(dlg.handles.hax,'on');
                        if dlg.uanalysestype(dlg.ianalysis)==4, dlg.plothistinfo2=dlg.plothistinfo(5); dlg.plothistinfo3=2*dlg.plothistinfo2; dlg.plothistinfo4=dlg.plothistinfo2;
                        else dlg.plothistinfo2=dlg.plothistinfo(2); dlg.plothistinfo3=dlg.plothistinfo(3); dlg.plothistinfo4=dlg.plothistinfo(4);
                        end
                        refshow1=0;refshow1n=0;refshow2=0;refshow2n=0;
                        for n=1:numel(dlg.results_patch),
                            if dlg.uanalysestype(dlg.ianalysis)==4,
                                dlg.handles.resultspatch(n,1)=patch(dlg.results_patch{n}{3},dlg.results_patch{n}{1}+dlg.plothistinfo2,'k','edgecolor','none','linestyle','-','facecolor',.9*[.8 .8 1],'facealpha',.25/sqrt(numel(dlg.results_patch)),'parent',dlg.handles.hax); % title('Connectivity histogram before denoising'); xlabel('Correlation (r)');
                                dlg.handles.resultspatch(n,2)=patch(dlg.results_patch{n}{2},dlg.results_patch{n}{1},'k','edgecolor','none','linestyle','-','facecolor',.9*[.8 .8 1],'facealpha',.25/sqrt(numel(dlg.results_patch)),'parent',dlg.handles.hax); %title('Connectivity histogram after denoising'); xlabel('Correlation (r)');
                            else
                                dlg.handles.resultspatch(n,1)=patch(dlg.results_patch{n}{1},dlg.results_patch{n}{3}+dlg.plothistinfo2,'k','edgecolor','k','linestyle',':','facecolor',.9*[.8 .8 1],'facealpha',.25/sqrt(numel(dlg.results_patch)),'parent',dlg.handles.hax); % title('Connectivity histogram before denoising'); xlabel('Correlation (r)');
                                dlg.handles.resultspatch(n,2)=patch(dlg.results_patch{n}{1},dlg.results_patch{n}{2},'k','edgecolor','k','linestyle',':','facecolor',.9*[.8 .8 1],'facealpha',.25/sqrt(numel(dlg.results_patch)),'parent',dlg.handles.hax); %title('Connectivity histogram after denoising'); xlabel('Correlation (r)');
                            end
                            if numel(dlg.results_patch{n})>4
                                refshow2=refshow2+dlg.results_patch{n}{5}; refshow2n=refshow2n+1;
                                refshow1=refshow1+dlg.results_patch{n}{4}; refshow1n=refshow1n+1;
                            end
                        end
                        for n=1:numel(dlg.results_patch),
                            if numel(dlg.results_line)>=n&&~isempty(dlg.results_line{n})
                                dlg.handles.resultsline(n,1)=plot(dlg.results_line{n}{2},dlg.results_line{n}{1},'k.','linewidth',1,'color',.75*[.8 .8 1],'parent',dlg.handles.hax);
                                dlg.handles.resultsline(n,2)=plot(dlg.results_line{n}{3},dlg.plothistinfo2+dlg.results_line{n}{1},'k.','linewidth',1,'color',.75*[.8 .8 1],'parent',dlg.handles.hax);
                            end
                        end
                        plot([-1 1 nan -1 1],[0 0 nan dlg.plothistinfo2 dlg.plothistinfo2],'k-','linewidth',1,'parent',dlg.handles.hax);
                        if numel(dlg.results_line)>0
                            dlg.handles.resultsline_add=[plot(0,0,'k-','visible','off','parent',dlg.handles.hax), plot(0,0,'k-','visible','off','parent',dlg.handles.hax)];
                        else dlg.handles.resultsline_add=[];
                        end
                        htlegend=[];
                        if numel(dlg.results_patch)>0
                            if dlg.uanalysestype(dlg.ianalysis)==4,
                                dlg.handles.resultspatch_add=[patch(dlg.results_patch{n}{2},dlg.results_patch{n}{1},'k','edgecolor','none','linestyle','-','facecolor','k','facealpha',.25','visible','off','parent',dlg.handles.hax),...
                                    patch(dlg.results_patch{n}{2},dlg.results_patch{n}{1},'k','edgecolor','none','linestyle','-','facecolor','k','facealpha',.25,'visible','off','parent',dlg.handles.hax)];
                            else
                                dlg.handles.resultspatch_add=[patch(dlg.results_patch{n}{1},dlg.results_patch{n}{2},'k','edgecolor','k','linestyle','-','facecolor','k','facealpha',.25','visible','off','parent',dlg.handles.hax),...
                                    patch(dlg.results_patch{n}{1},dlg.results_patch{n}{2},'k','edgecolor','k','linestyle','-','facecolor','k','facealpha',.25,'visible','off','parent',dlg.handles.hax)];
                            end
                            if refshow1n
                                if dlg.uanalysestype(dlg.ianalysis)==4,
                                    plot(refshow2/refshow2n,dlg.results_patch{n}{1}+dlg.plothistinfo2,'r.','linewidth',2,'parent',dlg.handles.hax);
                                    plot(refshow1/refshow1n,dlg.results_patch{n}{1},'r.','linewidth',2,'parent',dlg.handles.hax);
                                    htlegend=plot(refshow1/refshow1n,dlg.results_patch{n}{1},'r-','linewidth',2,'parent',dlg.handles.hax,'visible','off');
                                    plot([0 0],[min(dlg.results_patch{n}{1}) max(dlg.results_patch{n}{1})+eps]+dlg.plothistinfo2,'r-','linewidth',1);
                                    plot([0 0],[min(dlg.results_patch{n}{1}) max(dlg.results_patch{n}{1})+eps],'r-','linewidth',1);
                                else
                                    plot(dlg.results_patch{n}{1},refshow2/refshow2n+dlg.plothistinfo2,'r--','linewidth',2,'parent',dlg.handles.hax);
                                    htlegend=plot(dlg.results_patch{n}{1},refshow1/refshow1n,'r--','linewidth',2,'parent',dlg.handles.hax);
                                end
                            else
                                plot([0 0],[0 dlg.plothistinfo3],'k-','linewidth',1);
                            end
                        else dlg.handles.resultspatch_add=[];
                        end
                        %plot([-1 1;-1 1],[ylim;ylim]','k-',[-1 -1;1 1],[ylim;ylim],'k-');
                        
                        if ~isempty(dlg.results_info)&&isstruct(dlg.results_info{1})&&isfield(dlg.results_info{1},'units'), text(0,-dlg.plothistinfo3*.1,dlg.results_info{1}.units,'horizontalalignment','center','fontsize',11+font_offset,'parent',dlg.handles.hax);
                        else text(0,-dlg.plothistinfo3*.1,'Correlation coefficients (r)','horizontalalignment','center','fontsize',11+font_offset,'parent',dlg.handles.hax);
                        end
                        text(-.95,dlg.plothistinfo2*.25,'After denoising','horizontalalignment','left','fontsize',10+font_offset,'fontweight','bold','parent',dlg.handles.hax);
                        text(-.95,dlg.plothistinfo2+(dlg.plothistinfo3-dlg.plothistinfo2)*.25,'Before denoising','horizontalalignment','left','fontsize',10+font_offset,'fontweight','bold','parent',dlg.handles.hax);
                        if numel(in)==1,
                            ttitle=dlg.results_label{1};
                            tlabel={};
                            if iscell(ttitle), tlabel=ttitle(2:end); ttitle=ttitle{1}; end
                            text(0,-dlg.plothistinfo3*.175,tlabel,'horizontalalignment','center','fontsize',5+font_offset,'interpreter','none','parent',dlg.handles.hax);
                            if ~isempty(dlg.results_info)&&isstruct(dlg.results_info{1}),
                                if isfield(dlg.results_info{1},'IntersectionBefore')
                                    if isfield(dlg.results_info{1},'PercentSignificantBefore')
                                        text(-.95,dlg.plothistinfo2*.25-min(dlg.plothistinfo2,dlg.plothistinfo3-dlg.plothistinfo2)*.15,{sprintf('%.2f%c%.2f (%.1f%% match with NH)',dlg.results_info{1}.MeanAfter,177,dlg.results_info{1}.StdAfter,100*dlg.results_info{1}.IntersectionAfter),sprintf('%.1f%% edges with p<.05, %.1f%% edges with q<.05',100*dlg.results_info{1}.PercentSignificantAfter(1),100*dlg.results_info{1}.PercentSignificantAfter(2))},'horizontalalignment','left','fontsize',5+font_offset,'parent',dlg.handles.hax);
                                        text(-.95,dlg.plothistinfo2+(dlg.plothistinfo3-dlg.plothistinfo2)*.25-min(dlg.plothistinfo2,dlg.plothistinfo3-dlg.plothistinfo2)*.15,{sprintf('%.2f%c%.2f (%.1f%% match with NH)',dlg.results_info{1}.MeanBefore,177,dlg.results_info{1}.StdBefore,100*dlg.results_info{1}.IntersectionBefore),sprintf('%.1f%% edges with p<.05, %.1f%% edges with q<.05',100*dlg.results_info{1}.PercentSignificantBefore(1),100*dlg.results_info{1}.PercentSignificantBefore(2))},'horizontalalignment','left','fontsize',5+font_offset,'parent',dlg.handles.hax);
                                    else
                                        text(-.95,dlg.plothistinfo2*.25-min(dlg.plothistinfo2,dlg.plothistinfo3-dlg.plothistinfo2)*.15,{sprintf('%.2f%c%.2f (%.1f%% match with NH)',dlg.results_info{1}.MeanAfter,177,dlg.results_info{1}.StdAfter,100*dlg.results_info{1}.IntersectionAfter)},'horizontalalignment','left','fontsize',5+font_offset,'parent',dlg.handles.hax);
                                        text(-.95,dlg.plothistinfo2+(dlg.plothistinfo3-dlg.plothistinfo2)*.25-min(dlg.plothistinfo2,dlg.plothistinfo3-dlg.plothistinfo2)*.15,{sprintf('%.2f%c%.2f (%.1f%% match with NH)',dlg.results_info{1}.MeanBefore,177,dlg.results_info{1}.StdBefore,100*dlg.results_info{1}.IntersectionBefore)},'horizontalalignment','left','fontsize',5+font_offset,'parent',dlg.handles.hax);
                                    end
                                elseif isfield(dlg.results_info{1},'DofBefore')
                                    text(-.95,dlg.plothistinfo2*.25-min(dlg.plothistinfo2,dlg.plothistinfo3-dlg.plothistinfo2)*.15,{sprintf('%.2f%c%.2f (df=%.1f)',dlg.results_info{1}.MeanAfter,177,dlg.results_info{1}.StdAfter,dlg.results_info{1}.DofAfter)},'horizontalalignment','left','fontsize',5+font_offset,'parent',dlg.handles.hax);
                                    text(-.95,dlg.plothistinfo2+(dlg.plothistinfo3-dlg.plothistinfo2)*.25-min(dlg.plothistinfo2,dlg.plothistinfo3-dlg.plothistinfo2)*.15,{sprintf('%.2f%c%.2f (df=%.1f)',dlg.results_info{1}.MeanBefore,177,dlg.results_info{1}.StdBefore,dlg.results_info{1}.DofBefore)},'horizontalalignment','left','fontsize',5+font_offset,'parent',dlg.handles.hax);
                                elseif isfield(dlg.results_info{1},'CorrBefore')
                                    text(-.95,dlg.plothistinfo2*.25-min(dlg.plothistinfo2,dlg.plothistinfo3-dlg.plothistinfo2)*.15,{sprintf('R^2 = %.3f',dlg.results_info{1}.CorrAfter.^2)},'horizontalalignment','left','fontsize',5+font_offset,'parent',dlg.handles.hax);
                                    text(-.95,dlg.plothistinfo2+(dlg.plothistinfo3-dlg.plothistinfo2)*.25-min(dlg.plothistinfo2,dlg.plothistinfo3-dlg.plothistinfo2)*.15,{sprintf('R^2 = %.3f',dlg.results_info{1}.CorrBefore.^2)},'horizontalalignment','left','fontsize',5+font_offset,'parent',dlg.handles.hax);
                                else
                                    text(-.95,dlg.plothistinfo2*.25-min(dlg.plothistinfo2,dlg.plothistinfo3-dlg.plothistinfo2)*.15,{sprintf('%.2f%c%.2f',dlg.results_info{1}.MeanAfter,177,dlg.results_info{1}.StdAfter)},'horizontalalignment','left','fontsize',5+font_offset,'parent',dlg.handles.hax);
                                    text(-.95,dlg.plothistinfo2+(dlg.plothistinfo3-dlg.plothistinfo2)*.25-min(dlg.plothistinfo2,dlg.plothistinfo3-dlg.plothistinfo2)*.15,{sprintf('%.2f%c%.2f',dlg.results_info{1}.MeanBefore,177,dlg.results_info{1}.StdBefore)},'horizontalalignment','left','fontsize',5+font_offset,'parent',dlg.handles.hax);
                                end
                            end
                            if ~isempty(ttitle), text(0,dlg.plothistinfo3*1.05,ttitle,'horizontalalignment','center','fontsize',10+font_offset,'fontweight','bold','interpreter','none','parent',dlg.handles.hax); end
                        elseif ~isempty(dlg.results_info)&&isfield(dlg.results_info{1},'DofBefore')
                            text(-.95,dlg.plothistinfo2*.25-min(dlg.plothistinfo2,dlg.plothistinfo3-dlg.plothistinfo2)*.15,{sprintf('mean r=%.2f%c%.2f',mean(cellfun(@(x)x.MeanAfter,dlg.results_info)),177,std(cellfun(@(x)x.MeanAfter,dlg.results_info)))},'horizontalalignment','left','fontsize',5+font_offset,'parent',dlg.handles.hax);
                            text(-.95,dlg.plothistinfo2+(dlg.plothistinfo3-dlg.plothistinfo2)*.25-min(dlg.plothistinfo2,dlg.plothistinfo3-dlg.plothistinfo2)*.15,{sprintf('mean r=%.2f%c%.2f',mean(cellfun(@(x)x.MeanBefore,dlg.results_info)),177,std(cellfun(@(x)x.MeanBefore,dlg.results_info)))},'horizontalalignment','left','fontsize',5+font_offset,'parent',dlg.handles.hax);
                        elseif ~isempty(dlg.results_info)&&isfield(dlg.results_info{1},'MeanBefore')
                            text(-.95,dlg.plothistinfo2*.25-min(dlg.plothistinfo2,dlg.plothistinfo3-dlg.plothistinfo2)*.15,{sprintf('%.2f%c%.2f',mean(cellfun(@(x)x.MeanAfter,dlg.results_info)),177,mean(cellfun(@(x)x.StdAfter,dlg.results_info)))},'horizontalalignment','left','fontsize',5+font_offset,'parent',dlg.handles.hax);
                            text(-.95,dlg.plothistinfo2+(dlg.plothistinfo3-dlg.plothistinfo2)*.25-min(dlg.plothistinfo2,dlg.plothistinfo3-dlg.plothistinfo2)*.15,{sprintf('%.2f%c%.2f',mean(cellfun(@(x)x.MeanBefore,dlg.results_info)),177,mean(cellfun(@(x)x.StdBefore,dlg.results_info)))},'horizontalalignment','left','fontsize',5+font_offset,'parent',dlg.handles.hax);
                        end
                        if dlg.uanalysestype(dlg.ianalysis)==4,
                            plot([.99 1 1 .99 nan .99 1 1 .99],[dlg.plotminmax([1 1 2 2],1)' nan dlg.plothistinfo2+dlg.plotminmax([1 1 2 2],1)'],'k-','linewidth',1,'parent',dlg.handles.hax);
                            text(1.05*[1 1 1 1],[dlg.plotminmax(1:2,1)' dlg.plothistinfo2+dlg.plotminmax(1:2,1)'],arrayfun(@(n)sprintf('%d mm',round(n)),[dlg.plotminmax(1:2,1)' dlg.plotminmax(1:2,1)'],'uni',0),'parent',dlg.handles.hax);
                            text(1.10,dlg.plothistinfo2/2,'Distance (mm)','rotation',90,'horizontalalignment','center','fontsize',11+font_offset,'parent',dlg.handles.hax);
                            text(1.10,dlg.plothistinfo2*1.5,'Distance (mm)','rotation',90,'horizontalalignment','center','fontsize',11+font_offset,'parent',dlg.handles.hax);
                        end
                        hold(dlg.handles.hax,'off');
                        if ~isempty(htlegend)
                            if dlg.uanalysestype(dlg.ianalysis)==4,
                                try, ht=legend(htlegend,{sprintf('expected mean%cstd under Null Hypothesis (if no QC-FC associations exist at any distance level; random permutations)',177)},'location','northwest'); set(ht,'box','off','fontsize',5+font_offset); catch, legend(ht,'randomised reference'); end
                            else
                                try, ht=legend(htlegend,{'expected shape of distribution under Null Hypothesis (if no QC-FC associations exist; random permutations)'},'location','northwest'); set(ht,'box','off','fontsize',5+font_offset); catch, legend(ht,'randomised reference'); end
                            end
                        end
                        if dlg.uanalysestype(dlg.ianalysis)==4, set(dlg.handles.hax,'ylim',[0 2*dlg.plothistinfo2]);
                        else set(dlg.handles.hax,'ylim',dlg.plothistinfo([1 3]))
                        end
                        set(dlg.handles.hax,'xlim',[-1,1],'ytick',[],'ycolor',figcolor,'ydir','normal','visible','on');
                        conn_qaplotsexplore_figuremousemove([],[],'updatemousetrack');
                    case 5 % QA_COV
                        delete(dlg.handles.hax(ishandle(dlg.handles.hax)));
                        pos=[.30 .175 .55 .575];
                        if get(dlg.handles.showannot,'value'), pos=[pos(1)+.225 pos(2) pos(3)-.20 pos(4)]; end
                        dlg.plothistinfo4=1;
                        dlg.handles.hax=axes('units','norm','position',pos,'color',figcolor,'visible','off','parent',dlg.handles.hfig);
                        dlg.results_patch=dlg.dataA(in); %%%
                        dlg.results_label=dlg.labelA(in);
                        dlg.results_info=dlg.infoA(in);
                        dlg.results_line=dlg.lineA(in);
                        dlg.handles.resultspatch=[];
                        dlg.handles.resultsline=[];
                        hold(dlg.handles.hax,'on');
                        refshow1=0;refshow1n=0;refshow2=0;refshow2n=0;
                        dlg.handles.resultspatch=[];
                        tx=[];ty=[];
                        for n=1:numel(dlg.results_patch),
                            if numel(in)>1
                                mask=dlg.results_line{n}{2}>dlg.results_info{n}.InterquartilesDisplay(5,:) | dlg.results_line{n}{2}<dlg.results_info{n}.InterquartilesDisplay(1,:);
                                if nnz(mask),
                                    ttx=reshape(dlg.results_line{n}{1}(mask),[],1);
                                    tty=reshape(dlg.results_line{n}{2}(mask),[],1);
                                    plot([ttx ttx-.25]',[tty tty+sign(tty-.5)*.05]','-','color',.5*[.8 .8 1],'parent',dlg.handles.hax);
                                    text(ttx-.25,tty+sign(tty-.5)*.05,regexprep(dlg.results_label{n}{1}{1},'^[sS]ubject ','S'),'color',.5*[.8 .8 1],'horizontalalignment','right','fontsize',5+font_offset,'parent',dlg.handles.hax);
                                end
                            end
                            if ~isequal(tx,dlg.results_patch{n}{1})||~isequal(ty,dlg.results_patch{n}{2});
                                tx=dlg.results_patch{n}{1};
                                ty=dlg.results_patch{n}{2};
                                dlg.handles.resultspatch(n,1)=patch(tx(:),ty(:),-ones(numel(tx),1),'k','edgecolor','none','linestyle','-','facecolor',.9*[.8 .8 1],'facealpha',.75,'parent',dlg.handles.hax);
                                tx2=[tx;nan(1,size(tx,2))];
                                ty2=[ty;nan(1,size(ty,2))];
                                plot(tx2(:),ty2(:),'k-','parent',dlg.handles.hax);
                                %for n2=1:size(dlg.results_patch{n}{1},2)
                                %    dlg.handles.resultspatch(n,n2)=patch(dlg.results_patch{n}{1}(:,n2),dlg.results_patch{n}{2}(:,n2),'k','edgecolor',.8*[.8 .8 1],'linestyle','-','facecolor',.9*[.8 .8 1],'facealpha',.25,'parent',dlg.handles.hax);
                                %end
                            end
                        end
                        dlg.handles.resultsline=[];
                        for n=1:numel(dlg.results_patch),
                            dlg.handles.resultsline(n,1)=plot(dlg.results_line{n}{1},dlg.results_line{n}{2},'ko','markerfacecolor',.25*[.8 .8 1],'markeredgecolor',.25*[.8 .8 1],'linewidth',1,'color',.75*[.8 .8 1],'parent',dlg.handles.hax);
                        end
                        if numel(in)==1,
                            ttitle=dlg.results_label{n}{1}{1};
                            text(numel(dlg.results_line{1}{1})/2+.5,1.1*1.05,ttitle,'horizontalalignment','center','fontsize',10+font_offset,'fontweight','bold','interpreter','none','parent',dlg.handles.hax);
                            set(dlg.handles.resultsline(:,1),'linestyle','-');
                        end
                        dlg.handles.resultsline_add=plot(0,0,'ko-','markeredgecolor',[1 1 1],'visible','off');
                        dlg.handles.resultspatch_add=[];
                        ht=[];tx=[];
                        for n=1:numel(dlg.results_info),
                            if ~isequal(tx,dlg.results_info{n}.InterquartilesDisplay)
                                tx=dlg.results_info{n}.InterquartilesDisplay;
                                ht(1)=plot([.5 1:size(tx,2) size(tx,2)+.5],tx(1,[1 1:end end]),'r--','linewidth',2,'parent',dlg.handles.hax);
                                ht(2)=plot([.5 1:size(tx,2) size(tx,2)+.5],tx(2,[1 1:end end]),'k:','linewidth',1,'parent',dlg.handles.hax);
                                %ht(3)=plot([.5 1:size(tx,2) size(tx,2)+.5],tx(3,[1 1:end end]),'k:','linewidth',1,'parent',dlg.handles.hax);
                                ht(3)=plot([.5 1:size(tx,2) size(tx,2)+.5],tx(4,[1 1:end end]),'k:','linewidth',1,'parent',dlg.handles.hax);
                                ht(4)=plot([.5 1:size(tx,2) size(tx,2)+.5],tx(5,[1 1:end end]),'r--','linewidth',2,'parent',dlg.handles.hax);
                                text(size(tx,2)+.55+zeros(1,4),tx([5,4,2,1],end)',{'3rd Q + 1.5 IQR','3rd Quartile','1st Quartile','1st Q - 1.5 IQR'},'horizontalalignment','left','fontsize',5+font_offset,'parent',dlg.handles.hax);
                                if numel(dlg.results_info{n}.Variables)>1,
                                    ty=dlg.results_info{n}.Interquartiles;
                                    text((1:size(tx,2))+.2,tx(1,:)-.02,arrayfun(@(x)mat2str(x,max(ceil(log10(abs(x))),2)),ty(1,:),'uni',0),'horizontalalignment','right','color',.5*[1 1 1],'fontsize',5+font_offset,'rotation',90,'parent',dlg.handles.hax);
                                    text((1:size(tx,2))+.2,tx(5,:)+.02,arrayfun(@(x)mat2str(x,max(ceil(log10(abs(x))),2)),ty(5,:),'uni',0),'horizontalalignment','left','color',.5*[1 1 1],'fontsize',5+font_offset,'rotation',90,'parent',dlg.handles.hax);
                                    text((1:size(tx,2)),-.15+zeros(1,size(tx,2)),regexprep(dlg.results_info{n}.Variables,'^QC_',''),'horizontalalignment','right','color','k','fontsize',5+font_offset,'rotation',90,'interpreter','none','parent',dlg.handles.hax);
                                end
                            end
                            %ht(1)=plot(1:size(dlg.results_info{n}.InterquartilesDisplay,2),dlg.results_info{n}.InterquartilesDisplay(2,:),'k--','linewidth',2,'parent',dlg.handles.hax);
                            %ht(2)=plot(1:size(dlg.results_info{n}.InterquartilesDisplay,2),dlg.results_info{n}.InterquartilesDisplay(4,:),'k--','linewidth',2,'parent',dlg.handles.hax);
                        end
                        if numel(dlg.results_info{1}.Variables)==1,
                            tx=dlg.results_info{n}.InterquartilesDisplay;
                            ty=dlg.results_info{n}.Interquartiles;
                            for n1=[1,2,4,5], text(.45,tx(n1,:),arrayfun(@(x)mat2str(x,max(ceil(log10(abs(x))),2)),ty(n1,:),'uni',0),'horizontalalignment','right','color',.5*[1 1 1],'fontsize',5+font_offset,'parent',dlg.handles.hax); end
                            tx=dlg.results_info{1}.Variables;
                            text(1,-.15,[regexprep(dlg.results_info{1}.Variables,'^QC_',''),dlg.results_info{1}.Variables_descr],'horizontalalignment','center','color','k','fontsize',8+font_offset,'interpreter','none','parent',dlg.handles.hax);
                        end
                        hold(dlg.handles.hax,'off');
                        set(dlg.handles.hax,'xlim',[.5,numel(dlg.results_line{1}{1})+.5],'ytick',[],'ycolor',figcolor,'ydir','normal','visible','off');
                        conn_qaplotsexplore_figuremousemove([],[],'updatemousetrack');

                end
                if get(dlg.handles.showannot,'value')
                    delete(dlg.handles.han(ishandle(dlg.handles.han)));
                    %idx=find(cellfun('length',dlg.txtA(in))>0);
                    if numel(in)==1
                        dlg.handles.han=[uicontrol('style','text','units','norm','position',[.20 .75 .20 .05],'string','annotations','horizontalalignment','center','backgroundcolor',figcolor,'foregroundcolor','k','parent',dlg.handles.hfig),...
                            uicontrol('style','edit','units','norm','position',[.20 .10 .20 .65],'max',2,'string',dlg.txtA{in},'horizontalalignment','left','backgroundcolor',figcolor,'foregroundcolor','k','callback',{@conn_qaplotsexplore_update,'annotate'},'parent',dlg.handles.hfig)];
                    elseif numel(in)>1
                        txt=arrayfun(@(n)sprintf('%s %s %s: %s',dlg.usubjects{dlg.isubjects(dlg.dataIDXplots(n))},dlg.usessions{dlg.isessions(dlg.dataIDXplots(n))},dlg.umeasures{dlg.imeasures(dlg.dataIDXplots(n))},sprintf('%s ',dlg.txtA{n}{:})),in,'uni',0);
                        dlg.handles.han=[uicontrol('style','text','units','norm','position',[.20 .75 .20 .05],'string','annotations','horizontalalignment','center','backgroundcolor',figcolor,'foregroundcolor','k','parent',dlg.handles.hfig),...
                            uicontrol('style','listbox','units','norm','position',[.20 .10 .20 .65],'max',1,'string',txt,'horizontalalignment','left','backgroundcolor',figcolor,'foregroundcolor','k','callback',{@conn_qaplotsexplore_update,'selectannotation'},'interruptible','off','parent',dlg.handles.hfig)];
                    end
                else
                    delete(dlg.handles.han(ishandle(dlg.handles.han)));
                end
                drawnow;
                
        end
    end

    function conn_qaplotsexplore_figuremousemove(hObject,eventdata,option,varargin)
        persistent lastpos_counter newpos refpos nlines mousetrack
        if nargin>=3&&isequal(option,'updatemousetrack')
            xlim=get(dlg.handles.hax,'xlim');ylim=get(dlg.handles.hax,'ylim');
            if dlg.uanalysestype(dlg.ianalysis)==5,     ix=1; iy=2; ii=dlg.results_line; bpos=0;
            elseif dlg.uanalysestype(dlg.ianalysis)==4, ix=[2,3]; iy=[1,1]; ii=dlg.results_line; bpos=[0 dlg.plothistinfo2];
            else,                                       ix=[1,1]; iy=[2,3]; ii=dlg.results_patch; bpos=[0 dlg.plothistinfo2];
            end
            zlim=256;
            mousetrack.idx=zeros(zlim,zlim);
            mousetrack.n1=zeros(zlim,zlim);
            mousetrack.xlim=xlim;
            mousetrack.ylim=ylim;
            mousetrack.zlim=zlim;
            for n1=1:numel(dlg.results_patch)
                for n2=1:numel(ix)
                    x1=1+round((zlim-1)*max(0,min(1, (ii{n1}{ix(n2)}-xlim(1))/(xlim(2)-xlim(1)) )));
                    y1=1+round((zlim-1)*max(0,min(1, (ii{n1}{iy(n2)}+bpos(n2)-ylim(1))/(ylim(2)-ylim(1)) )));
                    mousetrack.idx((x1-1)*zlim+y1)=1:numel(ii{n1}{1});
                    mousetrack.n1((x1-1)*zlim+y1)=n1;
                    mousetrack.labels((x1-1)*zlim+y1)=n2;
                end
            end
            [i1,i2,v]=find(mousetrack.n1);
            [n1,n2]=ndgrid(-5:5,-5:5);
            [nill,i3]=sort(sqrt(abs(n1(:)).^2+abs(n2(:)).^2));
            for n3=reshape(i3,1,[])
                idx=max(1,min(zlim, i1+n1(n3)))+zlim*(max(1,min(zlim, i2+n2(n3)))-1);
                val=mousetrack.n1(idx)==0;
                mousetrack.n1(idx(val))=v(val);
            end
            return
        end
        try
            p1=get(0,'pointerlocation');
            p2=get(dlg.handles.hfig,'position');
            p3=get(0,'screensize');
            p4=p2(1:2)+p3(1:2)-1; % note: fix issue when connecting to external monitor/projector
            pos0=(p1-p4);
            if pos0(1)/p2(3)<0.25, return; end
            set(dlg.handles.hfig,'currentpoint',pos0);
            pos=(get(dlg.handles.hax,'currentpoint'));
            pos=pos(1,1:3);
            switch(dlg.uanalysestype(dlg.ianalysis))
                case 1, % QA_NORM/QA_REG
                    pos=round(pos);
                    set(dlg.handles.hax,'units','pixels');posax=get(dlg.handles.hax,'position');set(dlg.handles.hax,'units','norm');
                    nX=dlg.dispsize;
                    if numel(nX)<5, return; end
                    txyz=conn_menu_montage('coords2xyz',nX,pos(1:2)');
                    if txyz(3)>=1&&txyz(3)<=nX(end)&&txyz(1)>=1&&pos(1)<=nX(3)*nX(1)&&pos(2)>=1&&pos(2)<=nX(4)*nX(2)&&isfield(dlg,'followmouse')&&dlg.followmouse>0
                        f1=dlg.dataDidx(txyz(2),txyz(1));
                        f2=dlg.dataDmax(txyz(2),txyz(1));
                        if f2>0||nX(end)>1,
                            tlabel={};
                            if nX(end)>1, tlabel=[{[dlg.usubjects{dlg.isubjects(dlg.dataIDXplots(dlg.dataIDXsubjects(txyz(3))))},' ',dlg.usessions{dlg.isessions(dlg.dataIDXplots(dlg.dataIDXsubjects(txyz(3))))}],' '},tlabel];
                            elseif f2>0, tlabel=[tlabel {'Most different from average at this location:',sprintf('%s %s (diff z=%.2f)',dlg.usubjects{dlg.isubjects(dlg.dataIDXplots(f1))},dlg.usessions{dlg.isessions(dlg.dataIDXplots(f1))},f2)}];
                            end
                            if isempty(lastpos_counter), lastpos_counter=inf; end
                            if lastpos_counter>100,
                                lastpos_counter=0;
                                set(dlg.handles.hlabel,'units','pixels','position',[pos0+[10 -10] 20 20],'visible','on','string',tlabel);%,'fontsize',8+4*f2);
                                hext=get(dlg.handles.hlabel,'extent');
                                nlines=ceil(hext(3)/(p2(3)/2));
                                refpos=pos0;
                                newpos=[pos0+[-min(p2(3)/2,hext(3))/2 +10] min(p2(3)/2,hext(3)) nlines*hext(4)]; % text position figure coordinates
                                newpos(1)=max(posax(1),newpos(1)-max(0,newpos(1)+newpos(3)-posax(1)-posax(3)));
                                newpos(2)=max(posax(2),newpos(2)-max(0,newpos(2)+newpos(4)-posax(2)-posax(4)));
                                ntlabel=numel(tlabel);
                                set(dlg.handles.hlabel,'position',newpos,'string',reshape([tlabel,repmat(' ',1,nlines*ceil(ntlabel/nlines)-ntlabel)]',[],nlines)');
                            else
                                lastpos_counter=lastpos_counter+1;
                                ntlabel=numel(tlabel);
                                set(dlg.handles.hlabel,'visible','on','position',[newpos(1:2)+(pos0(1:2)-refpos(1:2)) newpos(3:4)],'string',reshape([tlabel,repmat(' ',1,nlines*ceil(ntlabel/nlines)-ntlabel)]',[],nlines)');
                            end
                        else
                            set(dlg.handles.hlabel,'visible','off');
                        end
                    else
                        set(dlg.handles.hlabel,'visible','off');
                    end
                case {2,3,4,5}, %QA_DENOISE
                    posb=pos;
                    if dlg.uanalysestype(dlg.ianalysis)==5, labels=2;
                    elseif pos(2)>=dlg.plothistinfo2&&pos(2)<=dlg.plothistinfo3, labels=3; %posb(2)=posb(2)-dlg.plothistinfo2; 
                    elseif pos(2)>=dlg.plothistinfo(1)&&pos(2)<=dlg.plothistinfo2, labels=2; 
                    else pos=[]; 
                    end
                    if ~isempty(pos)
                        if ~isempty(mousetrack)&&isfield(mousetrack,'n1')
                            dwin=[];
                            x1=1+round((mousetrack.zlim-1)*max(0,min(1, (posb(1)-mousetrack.xlim(1))/(mousetrack.xlim(2)-mousetrack.xlim(1)) )));
                            y1=1+round((mousetrack.zlim-1)*max(0,min(1, (posb(2)-mousetrack.ylim(1))/(mousetrack.ylim(2)-mousetrack.ylim(1)) )));
                            n1=mousetrack.n1(y1,x1);
                            if n1>0
                                idx=mousetrack.idx(y1,x1);
                                dwin=n1;
                                dpos=posb(1:2);
                                %if dlg.uanalysestype(dlg.ianalysis)==4, dpos=[dlg.results_line{n1}{labels}(idx) dlg.results_line{n1}{1}(idx)]; 
                                %else dpos=[dlg.results_patch{n1}{1}(idx) dlg.results_patch{n1}{labels}(idx)]; 
                                %end
                                %posb(1:2)=dpos;
                            end
                        elseif 0
                            dwin=[];dmin=inf;
                            if dlg.uanalysestype(dlg.ianalysis)==5,
                                for n1=1:numel(dlg.results_patch)
                                    [d,idx]=min(abs(dlg.results_line{n1}{1}-posb(1))+abs((dlg.results_line{n1}{2}-posb(2))));
                                    if d<dmin, dmin=d; dwin=n1; dpos=[dlg.results_line{n1}{1}(idx) dlg.results_line{n1}{2}(idx)]; end
                                end
                            elseif dlg.uanalysestype(dlg.ianalysis)==4,
                                for n1=1:numel(dlg.results_patch)
                                    [d,idx]=min(abs(dlg.results_line{n1}{labels}-posb(1))+abs((dlg.results_line{n1}{1}-posb(2))/dlg.plothistinfo4));
                                    if d<dmin, dmin=d; dwin=n1; dpos=[dlg.results_line{n1}{labels}(idx) dlg.results_line{n1}{1}(idx)]; end
                                end
                            elseif dlg.uanalysestype(dlg.ianalysis)==2,
                                %if nargin<=2&&rand<.9, return; end
                                [nill,idx]=min(abs(dlg.results_patch{1}{1}-posb(1)));
                                d=zeros(1,numel(dlg.results_patch));
                                for n1=1:numel(dlg.results_patch), d(n1)=dlg.results_patch{n1}{labels}(idx); end
                                [dmin,dwin]=min(abs(d-posb(2)));
                                dpos=[dlg.results_patch{dwin}{1}(idx) dlg.results_patch{dwin}{labels}(idx)];
                            else
                                %if nargin<=2&&rand<.9, return; end
                                for n1=1:numel(dlg.results_patch)
                                    [d,idx]=min(abs(dlg.results_patch{n1}{1}-posb(1))+abs((dlg.results_patch{n1}{labels}-posb(2))/dlg.plothistinfo4));
                                    if d<dmin, dmin=d; dwin=n1; dpos=[dlg.results_patch{n1}{1}(idx) dlg.results_patch{n1}{labels}(idx)]; end
                                end
                            end
                        end
                        if ~isempty(dwin)&&max(abs(dpos-posb(1:2))./[1 dlg.plothistinfo4])<.10
                            if numel(dlg.results_patch)>1&&nargin>2&&isequal(option,'buttonup')
                                n=dwin;
                                subjects=dlg.isubjects(dlg.dataIDXplots(dlg.dataIDXsubjects(n)));
                                sessions=dlg.isessions(dlg.dataIDXplots(dlg.dataIDXsubjects(n)));
                                measures=dlg.imeasures(dlg.dataIDXplots(dlg.dataIDXsubjects(n)));
                                set(dlg.handles.subjects,'value',find(dlg.usubjects_shown==subjects));
                                set(dlg.handles.sessions,'value',find(dlg.usessions_shown==sessions));
                                set(dlg.handles.measures,'value',find(dlg.umeasures_shown==measures));
                                conn_qaplotsexplore_update([],[],'refresh');
                            else
                                %tlabel=[dlg.usubjects{dlg.isubjects(dlg.dataIDXplots(dlg.dataIDXsubjects(dwin)))},' ',dlg.usessions{dlg.isessions(dlg.dataIDXplots(dlg.dataIDXsubjects(dwin)))}];
                                tlabel=dlg.results_label{dwin};
                                if iscell(tlabel)&&~isempty(tlabel), tlabel=tlabel{1}; end
                                if dlg.uanalysestype(dlg.ianalysis)==5, tlabel=tlabel([1,max(2,min(numel(tlabel), 1+round(dpos(1))))]); end
                                set(dlg.handles.hlabel,'units','pixels','position',[pos0+[10 -10] 20 20],'visible','on','string',tlabel);
                                hext=get(dlg.handles.hlabel,'extent');
                                nlines=ceil(hext(3)/(p2(3)/2));
                                ntlabel=numel(tlabel);
                                set(dlg.handles.hlabel,'position',[pos0+[-min(p2(3)/2,hext(3))-10 -10] min(p2(3)/2,hext(3)) nlines*hext(4)],'string',reshape([tlabel,repmat(' ',1,nlines*ceil(ntlabel/nlines)-ntlabel)]',[],nlines)');
                                if dlg.uanalysestype(dlg.ianalysis)==5
                                    set(dlg.handles.resultsline_add,'xdata',get(dlg.handles.resultsline(dwin,1),'xdata'),'ydata',get(dlg.handles.resultsline(dwin,1),'ydata'),'zdata',get(dlg.handles.resultsline(dwin,1),'zdata'),'visible','on');
                                else
                                    set(dlg.handles.resultspatch_add(1),'xdata',get(dlg.handles.resultspatch(dwin,1),'xdata'),'ydata',get(dlg.handles.resultspatch(dwin,1),'ydata'),'zdata',get(dlg.handles.resultspatch(dwin,1),'zdata'),'visible','on');
                                    set(dlg.handles.resultspatch_add(2),'xdata',get(dlg.handles.resultspatch(dwin,2),'xdata'),'ydata',get(dlg.handles.resultspatch(dwin,2),'ydata'),'zdata',get(dlg.handles.resultspatch(dwin,2),'zdata'),'visible','on');
                                    if size(dlg.handles.resultsline,1)>=dwin,
                                        set(dlg.handles.resultsline_add(1),'xdata',get(dlg.handles.resultsline(dwin,1),'xdata'),'ydata',get(dlg.handles.resultsline(dwin,1),'ydata'),'zdata',get(dlg.handles.resultsline(dwin,1),'zdata'),'visible','on');
                                        set(dlg.handles.resultsline_add(2),'xdata',get(dlg.handles.resultsline(dwin,2),'xdata'),'ydata',get(dlg.handles.resultsline(dwin,2),'ydata'),'zdata',get(dlg.handles.resultsline(dwin,2),'zdata'),'visible','on');
                                    end
                                end
                            end
                        else
                            set(dlg.handles.hlabel,'visible','off','string','');
                            set(dlg.handles.resultsline_add,'visible','of');
                            set(dlg.handles.resultspatch_add,'visible','of');
                        end
                    else
                        set(dlg.handles.hlabel,'visible','off','string','');
                        set(dlg.handles.resultspatch_add,'visible','of');
                        set(dlg.handles.resultsline_add,'visible','of');
                    end
            end
        end
    end
end

function [descrip, procedure, proceduretype]=conn_qaplotsexplore_translate(root)
% root_list={'^QA_NORM_(.*)','^QA_REG_functional','^QA_REG_(.*?)_?functional','^QA_REG_(.*?)_?structural','^QA_REG_(.*?)_?mni','^QA_COREG_(.*)','^QA_TIME_(.*)','^QA_TIMEART_(.*)','^QA_DENOISE_timeseries','^QA_DENOISE_QC-FC','^QA_DENOISE_scatterplot','^QA_DENOISE','^QA_SPM_design','^QA_SPM_contrasts'};
% root_descrip={'QA normalization: $1 data + outline of MNI TPM template','QA registration: functional data + structural overlay','QA registration: functional data + outline of ROI $1','QA registration: structural data + outline of ROI $1','QA registration: mni reference template + outline of ROI $1','QA realignment: $1 center-slice across multiple sessions/datasets','QA artifacts: $1 movie across all timepoints/acquisitions','QA artifacts: BOLD GS changes & subject motion timeseries with $1 movie','QA denoising: BOLD signal traces (carpetplot) before and after denoising + ART timeseries','QA denoising: distribution of QC-FC associations before and after denoising','QA denoising: scatterplot of functional correlations (FC) vs. distance (mm) before and after denoising','QA denoising: distribution of functional correlations (FC) before and after denoising','QA SPM design: review SPM first-level design matrix','QA SPM contrasts: review SPM first-level contrasts'};
%
root_list={...
    '^QA_NORM_structural.*?(\(.*?\)\s*)?$','QA normalization: structural data + outline of MNI TPM template $1','1','1';
    '^QA_NORM_functional.*?(\(.*?\)\s*)?$','QA normalization: functional data + outline of MNI TPM template $1','2','1';
    '^QA_NORM_(.*?)(\(.*?\)\s*)?$','QA normalization: $1 data + outline of MNI TPM template $2','3','1';
    '^QA_REG_functional.*?(\(.*?\)\s*)?$','QA registration: functional data + structural overlay $1','10','1';
    '^QA_REG_(.*?)_?functional.*?(\(.*?\)\s*)?$','QA registration: functional data + outline of ROI $1 $2','5','1';
    '^QA_REG_(.*?)_?structural.*?(\(.*?\)\s*)?$','QA registration: structural data + outline of ROI $1 $2','4','1';
    '^QA_REG_(.*?)_?mni.*?(\(.*?\)\s*)?$','QA registration: mni reference template + outline of ROI $1 $2','6','1';
    '^QA_COREG_(.*?)(\(.*?\)\s*)?$','QA realignment: $1 center-slice across multiple sessions/datasets $2','7','1';
    '^QA_TIME_(.*?)(\(.*?\)\s*)?$','QA artifacts: $1 movie across all timepoints/acquisitions $2','8','1';
    '^QA_TIMEART_(.*?)(\(.*?\)\s*)?$','QA artifacts: BOLD GS changes & subject motion timeseries with $1 movie $2','9','1';
    '^QA_DENOISE_timeseries.*?(\(.*?\)\s*)?$','QA denoising: BOLD signal traces (carpetplot) + ART timeseries $1','12','1';
    '^QA_DENOISE_QC-FC_scatterplot.*?(\(.*?\)\s*)?$','QA denoising: scatterplot of QC-FC associations (QC-FC) vs. distance (mm) $1','15','4';
    '^QA_DENOISE_QC-FC.*?(\(.*?\)\s*)?$','QA denoising: distribution of QC-FC associations (QC-FC) $1','13','3';
    '^QA_DENOISE_scatterplot.*?(\(.*?\)\s*)?$','QA denoising: scatterplot of functional connectivity values (FC) vs. distance (mm) $1','14','4';
    '^QA_DENOISE.*?(\(.*?\)\s*)?$','QA denoising: distribution of functional connectivity values (FC) $1','11','2';
    '^QA_COV.*?(\(.*?\)\s*)?$','QA variables: distribution of subject-level QC measures $1','31','5';
    '^QA_SPM_design.*?(\(.*?\)\s*)?$','QA SPM design: review SPM first-level design matrix $1','21','1';
    '^QA_SPM_contrasts.*?(\(.*?\)\s*)?$','QA SPM contrasts: review SPM first-level contrasts $1','22','1';
    '^QA_SPM_results.*?(\(.*?\)\s*)?$','QA SPM results: review SPM contrast effect-sizes $1','23','1';
    '^(\D.*)','$1','0','1'};

descrip=regexprep(root,...
    root_list(:,1)',root_list(:,2)');
procedure=cellfun(@str2num,regexprep(root,...
    root_list(:,1)',root_list(:,3)'));
proceduretype=cellfun(@str2num,regexprep(root,...
    root_list(:,1)',root_list(:,4)'));
end

function conn_qaplotsexplore_callbackgui(option,h1,h2,h3,h4,h0,procedures,h5,h6)
global CONN_x;
switch(option)
    case 'update'
        v=get(h1,'value');
        h=[h2,h3,h4];
        if isempty(v)
            set([h h0],'visible','off');
        else
            v=procedures(v);
            vv=[any(ismember(v,[2,7,8,9,10])),any(ismember(v,[3:6])),any(ismember(v,[13,15,31]))];
            set(h(vv==0),'visible','off');
            set(h(vv>0),'visible','on');
            b=0;
            for n=reshape(find(vv),1,[])
                set(h(n),'position',[.45,.25+b,.45,.30/(1.2*sum(vv)-0.2)]);
                b=b+.30/sum(vv);
            end
            if any(vv>0), set(h0,'visible','on');
            else set(h0,'visible','off');
            end
        end
        if get(h5,'value'), set(h6,'visible','off');
        else set(h6,'visible','on');
        end
    case 'group'
        v=listdlg('liststring',CONN_x.Setup.l2covariates.names(1:end-1),'selectionmode','single','initialvalue',1,'promptstring',{'Select group-defining covariate','(0/1 values defining subjects to include in QA plots)'},'ListSize',[400 250]);
        if ~isempty(v),
            values=conn_module('get','l2covariates',CONN_x.Setup.l2covariates.names{v});
            valid=find(~isnan(values)&values~=0);
            set(h1,'value',valid);
        end
end
end

function m=conn_qaplotsexplore_getdata(data,in);
if iscell(data)
    m=[];
    for n=1:numel(in),
        a=conn_fileutils('imread',data{in(n)});
        if isa(a,'uint8'), a=double(a)/255; end
        if isempty(m), m=zeros([size(a,1),size(a,2),size(a,3),numel(in)]); end
        if size(a,1)>size(m,1), m(size(a,1),1,1,1)=0; end
        if size(a,2)>size(m,2), m(1,size(a,2),1,1)=0; end
        if size(a,3)>size(m,3), m(1,1,size(a,3),1)=0; end
        if size(m,1)>size(a,1), a(size(m,1),1,1,1)=0; end
        if size(m,2)>size(a,2), a(1,size(m,2),1,1)=0; end
        if size(m,3)>size(a,3), a(1,1,size(m,3),1)=0; end
        %if mean(data(:))<.5, data=1-data; end
        %if size(data,3)==3, data=data.*repmat(shiftdim(figcolor,-1),[size(data,1),size(data,2)]); end
        m(:,:,:,n)=a;
    end
else
    m=data(:,:,:,in);
end
end

function dlg=conn_qaplotsexplore_updatemean(dlg)
if isequal(dlg.dataM,0), [dlg.dataM,dlg.dataS,dlg.dataA,dlg.txtA,dlg.dataD,dataDmax,dataDidx]=conn_qaplotsexplore_readfiles(dlg.files_jpg,dlg.files_txt,dlg.dataIDXplots,false,[],false); end
end

