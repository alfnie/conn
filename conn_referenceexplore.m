function fh=conn_referenceexplore(varargin)
% CONN_REFERENCEEXPLORE: GUI for conn_reference
%

global CONN_x CONN_gui;
guifields=struct(...
    'include',[],...        % [0/1,0/1,0/1,0/1] describe fields to display initially
    'preproclog',[],...     % number of SetupPreproc.log entry to display initially
    'firstlevelinfo',[],... % number/name of Analysis.name entry to display initially
    'secondlevelinfo',[],...% number/name of thresholding option to display initially
    'fileout',[]);

for n=1:2:numel(varargin)-1, assert(isfield(guifields,varargin{n}),'unrecognized option %s',varargin{n}); guifields.(varargin{n})=varargin{n+1}; end

if isfield(CONN_gui,'font_offset'),font_offset=CONN_gui.font_offset; else font_offset=0; end
hfig=[];
fh=@(varargin)conn_referenceexplore_update([],[],varargin{:});

if isempty(guifields.fileout), 
    guifields.fileout=conn_fullfile('referencesfile.html');
    try, if ~CONN_gui.isremote, guifields.fileout=conn_fullfile(CONN_x.folders.methods,['referencesfile_',datestr(now,'yyyy_mm_dd'),'.html']); end; end
end
if ~isfield(CONN_x,'SetupPreproc')||~isfield(CONN_x.SetupPreproc,'log')||isempty(CONN_x.SetupPreproc.log), 
    options.preproclog_text={};
    options.preproclog={};
else
    str1={}; idx1=[];
    for n=1:numel(CONN_x.SetupPreproc.log),
        if isequal(CONN_x.SetupPreproc.log{n}{1},'timestamp'), 
            try
                tidx=find(strcmp(CONN_x.SetupPreproc.log{n}(1:2:end-1),'subjects'),1);
                str1{end+1}=[CONN_x.SetupPreproc.log{n}{2}, sprintf(' (%d subjects)',numel(CONN_x.SetupPreproc.log{n}{2*tidx}))]; idx1(end+1)=n;
            catch
                str1{end+1}=CONN_x.SetupPreproc.log{n}{2}; idx1(end+1)=n;
            end
        end
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
options.secondlevel_txt{1}={'standard settings #1: Random Field Theory parametric statistics','standard settings #2: Permutation/randomization analysis','standard settings #3: Threshold Free Cluster Enhancement'};
options.secondlevel_names{1}={'RFT','RANDOM','TFCE'};
options.secondlevel_txt{2}={'standard settings #1: Functional Network Connectivity','standard settings #2: Spatial Pairwise Clustering','standard settings #3: Threshold Free Cluster Enhancement','alternative settings for connection-based inferences: parametric univariate statistics ','alternative settings for ROI-based inferences: parametric multivariate statistics','alternative settings for network-based inferences: Network Based Statistics'};
options.secondlevel_names{2}={'RRC_FNC','RRC_SPC','RRC_TFCE','RRC_CON','RRC_ROI','RRC_NET'};

bgc=.9*[1 1 1];
figcolor=[1 1 1];%[.95 .95 .9];
dlg.handles.hfig=hfig;
if isempty(dlg.handles.hfig)||~ishandle(dlg.handles.hfig), dlg.handles.hfig=figure('units','norm','position',[.1,.3,.8,.6],'menubar','none','numbertitle','off','name','Write methods','color',figcolor,'colormap',gray(256),'interruptible','off','busyaction','cancel','tag','conn_referenceexplore','userdata',fh);
else figure(dlg.handles.hfig); clf(dlg.handles.hfig);
end
uicontrol('style','frame','units','norm','position',[0,.71,1,.29],'backgroundcolor',bgc,'foregroundcolor',bgc,'fontsize',9+font_offset);
if isempty(guifields.include), guifields.include=true(1,4); end
if ~isfield(CONN_x,'SetupPreproc')||~isfield(CONN_x.SetupPreproc,'log')||isempty(CONN_x.SetupPreproc.log), guifields.include(1)=false; end
if ~CONN_x.isready(2), guifields.include(2)=false; end
if ~CONN_x.isready(3), guifields.include(3)=false; end
if ~CONN_x.isready(4), guifields.include(4)=false; end
dlg.handles.cb1=uicontrol('style','checkbox','units','norm','position',[.025,.925,.20,.05],'backgroundcolor',bgc,'foregroundcolor','k','horizontalalignment','left','value',guifields.include(1),'string','Describe Preprocessing steps','fontweight','normal','fontsize',9+font_offset,'callback',{@conn_referenceexplore_update,'preprocessing'},'interruptible','off');
dlg.handles.cb2=uicontrol('style','checkbox','units','norm','position',[.275,.925,.20,.05],'backgroundcolor',bgc,'foregroundcolor','k','horizontalalignment','left','value',guifields.include(2),'string','Describe Denoising (1st-level) steps','fontweight','normal','fontsize',9+font_offset,'callback',{@conn_referenceexplore_update,'denoising'},'interruptible','off');
dlg.handles.cb3=uicontrol('style','checkbox','units','norm','position',[.525,.925,.20,.05],'backgroundcolor',bgc,'foregroundcolor','k','horizontalalignment','left','value',guifields.include(3),'string','Describe Analyses (1st-level) steps','fontweight','normal','fontsize',9+font_offset,'callback',{@conn_referenceexplore_update,'firstlevel'},'interruptible','off');
dlg.handles.cb4=uicontrol('style','checkbox','units','norm','position',[.775,.925,.20,.05],'backgroundcolor',bgc,'foregroundcolor','k','horizontalalignment','left','value',guifields.include(4),'string','Describe Results (2nd-level) steps','fontweight','normal','fontsize',9+font_offset,'callback',{@conn_referenceexplore_update,'secondlevel'},'interruptible','off');
if ~guifields.include(1), set(dlg.handles.cb1,'enable','off'); end
if ~guifields.include(2), set(dlg.handles.cb2,'enable','off'); end
if ~guifields.include(3), set(dlg.handles.cb3,'enable','off'); end
if ~guifields.include(4), set(dlg.handles.cb4,'enable','off'); end

if ischar(guifields.preproclog), guifields.preproclog=find(strcmp(guifields.preproclog,options.preproclog_text),1); end
if isempty(guifields.preproclog), guifields.preproclog=~isempty(options.preproclog_text); end
if ischar(guifields.firstlevelinfo), guifields.firstlevelinfo=find(strcmp(guifields.firstlevelinfo,options.analyses_listnames),1); end
if isempty(guifields.firstlevelinfo), guifields.firstlevelinfo=1; end
if ismember(options.analyses_listtype(guifields.firstlevelinfo),[1,5]), guifields.secondleveltype=2; %connection-level
else guifields.secondleveltype=1; %voxel-level
end
if ischar(guifields.secondlevelinfo), guifields.secondlevelinfo=find(strcmp(guifields.secondlevelinfo,options.secondlevel_names{guifields.secondleveltype}),1); end
if isempty(guifields.secondlevelinfo), guifields.secondlevelinfo=1; end
dlg.handles.txt1=uicontrol('style','text','units','norm','position',[.025,.88,.20,.04],'backgroundcolor',bgc,'foregroundcolor','k','string','Select preprocessing log:','fontsize',9+font_offset,'horizontalalignment','left');
dlg.handles.preproclog=uicontrol('style','popupmenu','units','norm','position',[.025,.75,.20,.13],'max',1,'backgroundcolor',bgc,'foregroundcolor','k','string',[{'select preprocessing log date'},options.preproclog_text],'value',1+guifields.preproclog,'tooltipstring','Select a preprocessing pipeline to include in the methods description','fontsize',9+font_offset,'callback',{@conn_referenceexplore_update,'preprocessing'},'interruptible','off');
dlg.handles.txt2=uicontrol('style','text','units','norm','position',[.525,.88,.20,.04],'backgroundcolor',bgc,'foregroundcolor','k','string','Select 1st-level analysis:','fontsize',9+font_offset,'horizontalalignment','left');
dlg.handles.firstlevelinfo=uicontrol('style','listbox','units','norm','position',[.525,.75,.20,.13],'max',2,'backgroundcolor',bgc,'foregroundcolor','k','string',options.analyses_listnames,'value',guifields.firstlevelinfo,'tooltipstring','Select one or several first-level analyses to include in the methods description','fontsize',9+font_offset,'callback',{@conn_referenceexplore_update,'firstlevel'},'interruptible','off');
dlg.handles.txt3=uicontrol('style','text','units','norm','position',[.775,.88,.20,.04],'backgroundcolor',bgc,'foregroundcolor','k','string','Select inferential method used:','fontsize',9+font_offset,'horizontalalignment','left');
dlg.handles.secondlevelinfo=uicontrol('style','popupmenu','units','norm','position',[.775,.75,.20,.13],'max',1,'backgroundcolor',bgc,'foregroundcolor','k','string',options.secondlevel_txt{guifields.secondleveltype},'value',guifields.secondlevelinfo,'tooltipstring','Select inferential method / false-positive control option used in second-level analyses','fontsize',9+font_offset,'callback',{@conn_referenceexplore_update,'secondlevel'},'interruptible','off');

warning('off','MATLAB:ui:javacomponent:FunctionToBeRemoved');
jBrowserPanel = javaObjectEDT(com.mathworks.mlwidgets.help.LightweightHelpPanel);
[dlg.handles.html, dlg.handles.text] = javacomponent(jBrowserPanel, [], gcf);
set(dlg.handles.text, 'Units','norm','Position',[.025,.05,.95,.625]);
uicontrol('style','pushbutton','units','norm','position',[.575,.0,.20,.04],'string','Export to Word','fontsize',9+font_offset,'horizontalalignment','center','callback',{@conn_referenceexplore_export,'word'});
uicontrol('style','pushbutton','units','norm','position',[.775,.0,.20,.04],'string','Export to html','fontsize',9+font_offset,'horizontalalignment','center','callback',{@conn_referenceexplore_export,'html'});

%dlg.handles.text=uicontrol('style','listbox','units','norm','position',[.025,.05,.95,.625],'string','','max',2,'backgroundcolor',figcolor,'foregroundcolor','k','HorizontalAlignment','left','fontsize',32+font_offset);

conn_referenceexplore_update([],[],'init');
if ~ishandle(dlg.handles.hfig), return; end

    function conn_referenceexplore_update(hObject,eventdata,option,varargin)
        if isfield(dlg,'handles')&&isfield(dlg.handles,'hfig')&&~ishandle(dlg.handles.hfig), return; end
        switch(lower(option))
            case {'init','preprocessing','denoising','firstlevel','secondlevel'}
                steps={'init'};
                if get(dlg.handles.cb1,'value')>0, steps{end+1}='preprocessing'; set([dlg.handles.txt1,dlg.handles.preproclog],'enable','on'); 
                else set([dlg.handles.txt1,dlg.handles.preproclog],'enable','off'); 
                end                
                if get(dlg.handles.cb2,'value')>0, steps{end+1}='denoising'; end
                if get(dlg.handles.cb3,'value')>0, steps{end+1}='firstlevel'; set([dlg.handles.txt2,dlg.handles.firstlevelinfo],'enable','on'); 
                else set([dlg.handles.txt2,dlg.handles.firstlevelinfo],'enable','off');
                end
                if get(dlg.handles.cb4,'value')>0, steps{end+1}='secondlevel'; set([dlg.handles.txt3,dlg.handles.secondlevelinfo],'enable','on'); 
                else set([dlg.handles.txt3,dlg.handles.secondlevelinfo],'enable','off'); 
                end
                idx1=get(dlg.handles.preproclog,'value')-1;
                if idx1==0&&~isempty(options.preproclog_text), idx1=1; set(dlg.handles.preproclog,'value',idx1+1); end
                idx2=get(dlg.handles.firstlevelinfo,'value');
                idx3=get(dlg.handles.secondlevelinfo,'value');
                if isempty(idx1)||isempty(options.preproclog), preproclog={};
                else preproclog=options.preproclog{idx1};
                end
                if isempty(idx2)||isempty(options.analyses_listnames), firstlevelinfo={};
                else firstlevelinfo=regexprep(options.analyses_listnames(idx2),' \(SBC\)$| \(RRC\)$','');
                end
                if ~isempty(idx2)&&all(ismember(options.analyses_listtype(idx2),[1,5])), secondleveltype=2; %connection-level
                else secondleveltype=1; %voxel-level
                end
                if isempty(idx3), secondlevelinfo={};
                else
                    idx3=max(1,min(numel(options.secondlevel_names{secondleveltype}),idx3));
                    secondlevelinfo=options.secondlevel_names{secondleveltype}(idx3);
                end
                set(dlg.handles.secondlevelinfo,'string',options.secondlevel_txt{secondleveltype},'value',idx3);
                %options.analyses_listtype=[options.analyses_listtype, 1];
                %options.analyses_listidx=[options.analyses_listidx, n1];
                str=conn_reference(steps,'preproclog',preproclog,'firstlevelinfo',firstlevelinfo,'secondlevelinfo',secondlevelinfo,'fileout',guifields.fileout);
                dlg.handles.html.getLightweightBrowser.load(['file://',guifields.fileout]);
                %set(dlg.handles.text,'string',str);
                %set(dlg.handles.text,'string',cellfun(@(x)['<HTML>',x,'</HTML>'],regexp(str,'(\n|<br>)+','split'),'uni',0));
                %set(dlg.handles.text,'string',['<HTML>',regexprep(str,'\n',''),'</HTML>']);
        end
    end

    function conn_referenceexplore_export(hObject,eventdata,option,varargin)
        switch(option)
            case 'html'
                [tfilename,tfilepath]=uiputfile('*.html','Save description as');
                if isequal(tfilename,0), return; end
                conn_fileutils('copyfile',guifields.fileout,fullfile(tfilepath,tfilename));
                fprintf('document exported to %s\n',fullfile(tfilepath,tfilename));
            case 'word'
                [tfilename,tfilepath]=uiputfile('*.docx','Save description as');
                if isequal(tfilename,0), return; end
                rpt = mlreportgen.dom.Document(fullfile(tfilepath,tfilename),'docx');
                htmlFileObj = mlreportgen.dom.HTMLFile(guifields.fileout);
                append(rpt,htmlFileObj);
                close(rpt);
                fprintf('document exported to %s\n',fullfile(tfilepath,tfilename));
        end
    end
end



