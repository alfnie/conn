
function h=conn_filesearchtool(varargin)
persistent cwd;
global CONN_gui CONN_h;
if ~isfield(CONN_gui,'font_offset'), conn_font_init; end
if ~isfield(CONN_gui,'parse_html'), CONN_gui.parse_html={'<HTML><FONT color=rgb(100,100,100)>','</FONT></HTML>'}; end
if ~isfield(CONN_gui,'rightclick'), CONN_gui.rightclick='right'; end
if ~isfield(CONN_gui,'backgroundcolor'), CONN_gui.backgroundcolor=[0.1400 0.1400 0.1400]; end
if ~isfield(CONN_gui,'fontname'), CONN_gui.fontname=get(0,'FixedWidthFontName'); end
if ~isfield(CONN_h,'screen')||~isfield(CONN_h.screen,'hfig'), CONN_h.screen.hfig=gcf; end
if isempty(cwd), cwd=conn_projectmanager('pwd'); end

if nargin<1 || ischar(varargin{1}),
    fntname=CONN_gui.fontname;
    fields={'position',[.5,.1,.4,.8],...
        'backgroundcolor',.9*[1,1,1],...
        'foregroundcolor',[0,0,0],...
        'fontname',fntname,...
        'title','',...
        'button','Import',...
        'buttonhelp','',...
        'type','files',...
        'folder',cwd,...
        'filter','*',...
        'regexp','.',...
        'callback','',...
        'reduced',false,...
        'filename','',...
        'localcopy',[],...
        'inserver',[],...
        'max',2};
    params=[];
    for n1=1:2:nargin, params=setfield(params,lower(varargin{n1}),varargin{n1+1}); end
    for n1=1:2:length(fields), if ~isfield(params,fields{n1}) | isempty(getfield(params,fields{n1})), params=setfield(params,fields{n1},fields{n1+1}); end; end;
    M=[params.position(1),params.position(3),0,0,0;params.position(2),0,params.position(4),0,0;0,0,0,params.position(3),0;0,0,0,0,params.position(4)]';
    if params.reduced, h.frame=conn_menu('frame2noborder',[1,0,0,1,1]*M);
    else h.frame=conn_menu('frame2',[1,0,0,1,1]*M);
    end
    h.strbutton=params.button;
    h.reduced=params.reduced;
    %uicontrol('style','frame','units','norm','position',[1,0,0,1,1]*M,'backgroundcolor',params.backgroundcolor);
    %axes('units','norm','position',[1,0,0,1,1]*M,'color',params.backgroundcolor,'xcolor',min(1,1*params.backgroundcolor),'ycolor',min(1,1*params.backgroundcolor),'xtick',[],'ytick',[],'box','on');
    if h.reduced
        h.filename=conn_menu('edit2',[1,.05,.90,.9,.10]*M,params.title,params.filename);
    else
        if ~isempty(params.title), uicontrol('style','text','units','norm','position',[1,.05,.8,.9,.175]*M,'foregroundcolor',params.foregroundcolor,'backgroundcolor',params.backgroundcolor,'string',params.title,'fontname',params.fontname,'fontsize',9+CONN_gui.font_offset,'fontweight','bold','horizontalalignment','left','parent',CONN_h.screen.hfig); end
        %if ~isempty(params.title), uicontrol('style','text','units','norm','position',[1,.05,.8,.9,.175]*M,'foregroundcolor',params.foregroundcolor,'backgroundcolor',params.backgroundcolor,'string',params.title,'fontname',params.fontname,'fontsize',9+CONN_gui.font_offset,'fontweight','bold','horizontalalignment','left','parent',CONN_h.screen.hfig); end
        h.filename=[];
    end
    %uicontrol('style','text','units','norm','position',[1,.05,.85,.2,.05]*M,'foregroundcolor',params.foregroundcolor,'backgroundcolor',params.backgroundcolor,'string','Folder','fontangle','normal','fontname','default','fontsize',8+CONN_gui.font_offset,'horizontalalignment','center');
    %uicontrol('style','text','units','norm','position',[1,.05,.15,.2,.05]*M,'foregroundcolor',params.foregroundcolor,'backgroundcolor',params.backgroundcolor,'string','Filter','fontangle','normal','fontname','default','fontsize',8+CONN_gui.font_offset,'horizontalalignment','center');
    %uicontrol('style','text','units','norm','position',[1,.05,.10,.2,.05]*M,'foregroundcolor',params.foregroundcolor,'backgroundcolor',params.backgroundcolor,'string','Regexp','fontangle','normal','fontname','default','fontsize',8+CONN_gui.font_offset,'horizontalalignment','center');
    %h.filter=uicontrol('style','edit','units','norm','position',[1,.25,.8,.7,.05]*M,'foregroundcolor',params.foregroundcolor,'backgroundcolor',params.backgroundcolor,'string',params.filter,'tooltipstring','Select a file name filter (wildcards may be used)','fontsize',8+CONN_gui.font_offset);
    %h.folder=uicontrol('style','edit','units','norm','position',[1,.25,.85,.7,.05]*M,'foregroundcolor',params.foregroundcolor,'backgroundcolor',params.backgroundcolor,'string',params.folder,'fontname','default','fontsize',8+CONN_gui.font_offset,'tooltipstring','Select the root folder');
    h.selectfile=uicontrol('style','edit','position',[1,.25,.85,.7,.05]*M,'string','','max',2,'visible','off','parent',CONN_h.screen.hfig);
    %h.find=uicontrol('style','togglebutton','units','norm','position',[1,.7,.925,.25,.05]*M,'value',0,'string','Find','fontname','default','fontsize',8+CONN_gui.font_offset,'horizontalalignment','center','tooltipstring','Recursively searchs file names matching the filter starting from the current folder');
    %h.files=uicontrol('style','listbox','units','norm','position',[1,.05,.2,.9,.55]*M,'foregroundcolor',params.foregroundcolor,'backgroundcolor',params.backgroundcolor,'string','','max',params.max,'fontname','default','fontsize',8+CONN_gui.font_offset,'tooltipstring','Displays file matches. Double-click a folder for browsing to a different location. Double-click a file to import it to the toolbox');
    if h.reduced, h.files=conn_menu('listbox2',[1,.05,.43,.9,.38]*M,'',''); 
    else h.files=conn_menu('listbox2',[1,.05,.38,.9,.48]*M,'',''); %,'<HTML>Displays file matches<br/> - Double-click a folder for browsing to a different location<br/> - Double-click a file to import it to the toolbox</HTML>');
    end
    set(h.files,'max',params.max);
    h.selected=uicontrol('style','text','units','norm','position',[1,.04,.20,.9,.10]*M,'foregroundcolor',params.foregroundcolor,'backgroundcolor',params.backgroundcolor,'string','','fontname',params.fontname,'fontsize',8+CONN_gui.font_offset,'horizontalalignment','center','parent',CONN_h.screen.hfig);
    %h.select=uicontrol('style','pushbutton','units','norm','position',[1,.7,.14,.25,.05]*M,'string','Select','fontname','default','fontsize',8+CONN_gui.font_offset,'horizontalalignment','center','tooltipstring','Enter selected file(s) or open selected folder','callback',{@conn_filesearchtool,'files',true});
    if ~isempty(params.inserver), h.inserver=params.inserver; else h.inserver=conn_projectmanager('inserver'); end
    if h.inserver, params.folder=conn_server('util_remotefile',params.folder); % switch to proper local/remote folder when initializing
    else params.folder=conn_server('util_localfile',params.folder);
    end
    cwd=params.folder;
    h.isfiles=strcmp(params.type,'files');
    if h.isfiles
        if ~isempty(params.buttonhelp), h.strbuttonhelp=params.buttonhelp; 
        else h.strbuttonhelp=[h.strbutton,'s selected file(s) or open selected folder'];
        end
        h.select=conn_menu('pushbuttonblue2',[1,.05,.31,.45,.07]*M,'',h.strbutton,h.strbuttonhelp,{@conn_filesearchtool,'files',true});
        set(h.select,'fontsize',10+CONN_gui.font_offset);
        if h.reduced
            h.find=conn_menu('pushbutton2',[1,.50,.31,.45,.07]*M,'','Cancel','');
            set(h.find,'value',0,'fontsize',10+CONN_gui.font_offset);
            h.selectspm=[];
        else
            h.find=conn_menu('pushbutton2',[1,.50,.31,.45,.07]*M,'','Find','Recursively searchs file names matching the files/filter options below starting from the current folder');
            set(h.find,'value',0,'fontsize',10+CONN_gui.font_offset);
            h.selectspm=conn_menu('pushbutton2',[1,.7,.92,.29,.07]*M,'','ALTselect','<HTML>Alternative GUIs for file selection (OS-specific, spm_select GUI, select data from this or other CONN projects, etc.)</HTML>');
        end
    else
        if ~isempty(params.buttonhelp), h.strbuttonhelp=params.buttonhelp; 
        else h.strbuttonhelp=[h.strbutton,' highlighted folder'];
        end
        h.select=conn_menu('pushbuttonblue2',[1,.05,.31,.45,.07]*M,'',h.strbutton,h.strbuttonhelp,{@conn_filesearchtool,'selectfolder',true});
        set(h.select,'fontsize',10+CONN_gui.font_offset);
        if h.reduced
            h.find=conn_menu('pushbutton',[1,.50,.31,.45,.07]*M,'','Cancel','');
            set(h.find,'value',0,'fontsize',10+CONN_gui.font_offset);
        else
            h.find=[];
        end
        h.selectspm=[];
    end
    h.find_state=0;
    h.titles_folder=conn_menu('text2',[1,.10,.170,.20,.04]*M,'','path : '); set(h.titles_folder,'horizontalalignment','left');
    h.folder=conn_menu('edit2',[1,.3,.170,.6,.04]*M,'',params.folder,'<HTML>Current folder<br/> - Type a full path to change to a different folder</HTML>');
    str=conn_filesearch_breakfolder(params.folder);
    if h.reduced, h.folderpopup=conn_menu('popup2',[1,.05,.83,.9,.04]*M,'',regexprep(str,'(.+)[\/\\]$','$1'),'<HTML>Current folder<br/> - Select to change to a different location in current directory tree</HTML>');
    else h.folderpopup=conn_menu('popup2',[1,.05,.88,.9,.04]*M,'',regexprep(str,'(.+)[\/\\]$','$1'),'<HTML>Current folder<br/> - Select to change to a different location in current directory tree</HTML>');
    end
    set(h.folderpopup,'value',numel(str));
    h.titles_filter=conn_menu('text2',[1,.10,.120,.20,.04]*M,'','files : '); set(h.titles_filter,'horizontalalignment','left');
    h.filter=conn_menu('edit2',[1,.30,.120,.6,.04]*M,'',params.filter,'<HTML>File name filter (standard format)<br/> - Show only files with names matching any of these patterns<br/> - Use ";" to define multiple filters</HTML>');
    h.titles_regexp=conn_menu('text2',[1,.10,.070,.20,.04]*M,'','filter : '); set(h.titles_regexp,'horizontalalignment','left');
    h.regexp=conn_menu('edit2',[1,.30,.070,.6,.04]*M,'',params.regexp,'<HTML>Additional full-path filename filter (regexp format)<br/> - Show only files with full-path filename matching this pattern<br/> - See <i>help regexp</i> for additional information');
    if ~isempty(params.localcopy)&&~h.reduced,  
        h.titles_localcopy=conn_menu('text2',[1,.10,.020,.20,.04]*M,'','options : '); set(h.titles_localcopy,'horizontalalignment','left');
        h.localcopy=conn_menu('popup2',[1,.30,.020,.65,.04]*M,'',{'import selected files','import copy of selected files'},['<HTML>Controls behavior of ''',h.strbutton,''' button:<br/> - <i>import selected files</i> : (default) selected files will be imported into your CONN project directly from their original locations/folders<br/> - <i>import copy of selected files</i> : selected files will be first copied to your project conn_*/data/BIDS folder and then imported into your CONN project <br/>(e.g. use this when importing data from read-only folders if the files need to be further modified or processed)</HTML>']); set(h.localcopy,'value',params.localcopy);
    else 
        h.titles_localcopy=[];
        h.localcopy=[];
    end
%     hc1=uicontextmenu;
%     h.selectspmalt0=uimenu(hc1,'Label','<HTML>Use OS gui to select individual file(s)</HTML>');%,'callback',{@conn_filesearchtool,'selectspmalt2'});
%     h.selectspmalt1=uimenu(hc1,'Label','<HTML>Use SPM gui to select individual file(s) (disregards filter field)</HTML>');%,'callback',{@conn_filesearchtool,'selectspmalt2'});
%     h.selectspmalt2=uimenu(hc1,'Label','<HTML>Use SPM gui to select individual volume(s) from 4d nifti files</HTML>');%,'callback',{@conn_filesearchtool,'selectspmalt2'});
%     h.selectspmalt3=uimenu(hc1,'Label','<HTML>Select file(s) already entered in this or a different CONN project</HTML>');%,'callback',{@conn_filesearchtool,'selectspmalt2'});
%     set(h.selectspm,'uicontextmenu',hc1);
    h.callback=params.callback;
    set([h.files,h.find,h.folder,h.folderpopup,h.filter,h.regexp,h.localcopy,h.select,h.selectspm, h.filename],'userdata',h);
    names={'files','find','selectspm'}; for n1=1:length(names), set(h.(names{n1}),'callback',{@conn_filesearchtool,names{n1}}); end
    names={'folder','folderpopup','filter','regexp','filename'}; for n1=1:length(names), set(h.(names{n1}),'callback',{@conn_filesearchtool,names{n1},true}); end
    set([h.find h.select h.selectspm h.filter h.regexp h.localcopy h.folder h.titles_folder h.titles_filter h.titles_regexp h.titles_localcopy],'visible','off');
    if h.reduced, conn_menumanager('onregion',[h.find h.select h.selectspm h.localcopy],1,params.position+[-.05 -.05 +.10 +.10],h.files);
    else conn_menumanager('onregion',[h.find h.select h.selectspm h.filter h.regexp h.localcopy h.folder h.titles_folder h.titles_filter h.titles_regexp h.titles_localcopy],1,params.position+[-.05 -.05 +.10 +.10],h.files);
    end
    conn_filesearchtool(h.folder,[],'folder',true);
else,
    h=get(varargin{1},'userdata');
    set(h.selected,'string','');
    doubleclick=strcmp(get(gcbf,'SelectionType'),'open');
    selection=nargin>3|doubleclick;
    if strcmp(varargin{3},'selectspm'),
        opts={'Use SPM gui to select individual file(s)',...
            'Use SPM gui to select individual volume(s) from 4d nifti files',...
            'Use OS gui to select individual file(s)',...
            'Select file(s) entered in this or a different CONN project'};
        optsnames={'selectspmalt','selectspmalt2','selectspmalt0','selectspmalt3'};
        answ=conn_questdlg('','select files',opts{:},opts{1});
        if isempty(answ), return; end
        [nill,idx]=ismember(answ,opts);
        varargin{3}=optsnames{idx};
    end
    switch(varargin{3}),
        case {'selectspmalt','selectspmalt1','selectspmalt2','selectspmalt3','selectspmalt0'}
            if strcmp(varargin{3},'selectspmalt'), 
                regfilter=regexprep(get(h.filter,'string'),{'\s*',';([^;]+)',';','\.','*','([^\$])$'},{'','\$|$1','','\\.','.*','$1\$'});
                names=spm_select(inf,regfilter,[],[],get(h.folder,'string'));
            elseif strcmp(varargin{3},'selectspmalt1'), 
                names=spm_select(inf,'any',[],[],get(h.folder,'string'));
            elseif strcmp(varargin{3},'selectspmalt2'), 
                names=spm_select(inf,'image',[],[],get(h.folder,'string'));
            elseif strcmp(varargin{3},'selectspmalt3'), 
                names=conn_filesearch_selectconn(get(h.filter,'string'));
            elseif strcmp(varargin{3},'selectspmalt0'), 
                filter=regexp(get(h.filter,'string'),';','split');
                [tname,tpath]=uigetfile(filter','Select file(s)',get(h.folder,'string'),'multiselect','on');
                if isequal(tname,0), return; end
                tname=cellstr(tname);
                names=char(cellfun(@(x)fullfile(tpath,x),tname,'uni',0));
            end
            if ~isempty(names)
                if iscell(h.callback),
                    if length(h.callback)>1, feval(h.callback{1},h.callback{2:end},names); else, feval(h.callback{1},names); end
                else, feval(h.callback,names); end
            end
        case {'folder','folderpopup','filter','regexp','files','selectfolder','filename'},
            %parse={[regexprep(CONN_gui.parse_html{1},'<FONT color=rgb\(\d+,\d+,\d+\)>','<FONT color=rgb(150,100,100)><i>'),'-'],regexprep(CONN_gui.parse_html{2},'<\/FONT>','</FONT></i>')};
            parse={[regexprep(CONN_gui.parse_html{1},'<FONT color=rgb\(\d+,\d+,\d+\)>','<FONT color=rgb(180,180,180)>'),'-'],regexprep(CONN_gui.parse_html{2},'<\/FONT>','</FONT>')};
            pathname=fliplr(deblank(fliplr(deblank(get(h.folder,'string')))));
            if strcmp(varargin{3},'filename')
                filename=get(h.filename,'string');
                if conn_fileutils('isdir',filename)
                    pathname=filename; filename=''; set(h.filename,'string',filename);
                else
                    [tpath,tname,text]=fileparts(filename);
                    if ~isempty(tpath),
                        pathname=tpath; filename=[tname,text]; set(h.filename,'string',filename);
                        if ~conn_fileutils('isdir',pathname), return; end
                    end
                end
            end
            if strcmp(pathname,'.')||strncmp(pathname,'.\',2)||strncmp(pathname,'./',2), pathname=conn_fullfile(conn_projectmanager('pwd'),pathname(3:end)); end
            if ~conn_fileutils('isdir',pathname), pathname=cwd; end
            if strcmp(varargin{3},'folderpopup')
                str=conn_filesearch_breakfolder(pathname);
                set(h.folder,'string',[str{1:get(h.folderpopup,'value')}]);
                pathname=fliplr(deblank(fliplr(deblank(get(h.folder,'string')))));
            end
            cwd=pathname;
            selectfolder=strcmp(varargin{3},'selectfolder');
            if ismember(varargin{3},{'files','selectfolder'}),
                if h.reduced&&selection&&~doubleclick, % disregards if list points to a different folder
                    filename=fliplr(deblank(fliplr(deblank(get(h.filename,'string')))));
                    pathname=fullfile(pathname,filename);
                else
                    %disp(get(h.files,'value'))
                    filename=get(h.files,'string');
                    filename=filename(get(h.files,'value'),:);
                    if isempty(filename), return; end
                    filename=fliplr(deblank(fliplr(deblank(filename(1,:)))));
                    if strncmp(filename,parse{1},numel(parse{1})), filename=fliplr(deblank(fliplr(deblank(filename(numel(parse{1})+1:end-numel(parse{2})))))); end
                    if strcmp(filename,'..'),
                        selectfolder=false;
                        idx=find(pathname=='/'|pathname=='\'); idx(idx==length(pathname))=[];
                        if ~isempty(idx), pathname=pathname(1:idx(end)); else return; end
                    elseif ~selectfolder
                        pathname=fullfile(pathname,filename);
                    end
                end
            end
            isdirectory=(conn_fileutils('isdir',pathname) || (~h.reduced&&~conn_existfile(pathname)));
            if ismember(varargin{3},{'files','selectfolder'})&&h.reduced&&~isdirectory, set(h.filename,'string',filename); end
            if ~selectfolder&&isdirectory&&selection, % cd to a new directory
                if h.inserver, pathname=conn_server('util_remotefile',pathname); end % note: comment this line to allow access to local drives
                str=conn_filesearch_breakfolder(pathname);
                results={[parse{1},'   ..',parse{2}]}; 
                names=conn_dirn(fullfile(pathname,'*'));
                if ~isempty(names)&&numel(str)==1&&h.inserver, names=[names(1) names(:)']; names(1).name='CONNSERVER'; names(1).isdir=true; end
                for n1=1:length(names), if names(n1).isdir&&~strcmp(names(n1).name,'.')&&~strcmp(names(n1).name,'..'), results{end+1}=[parse{1},'   ',names(n1).name,parse{2}]; end; end
                n0results=numel(results);
                if n0results>0, results=conn_sortfilenames(results); end
                if h.isfiles
                    filter=get(h.filter,'string');
                    filter2=get(h.regexp,'string');
                    if isempty(filter), filter='*'; end
                    if isempty(filter2), filter2='.'; end
                    [filternow,filter]=strtok(filter,';');
                    allnames={names.name}; allnames=allnames([names.isdir]==0);
                    while ~isempty(filternow),
                        if 1 % faster
                            namematch=cellfun('length',regexp(allnames,['^',regexprep(fliplr(deblank(fliplr(deblank(filternow)))),{'[\[\]\(\)\{\}\\\^\$\.\|\?\+]','\*'},{'\\$0','.*'}),'$']))>0;
                            if any(namematch)&&~isequal(filter2,'.'), namematch(namematch)=cellfun('length',regexp(allnames(namematch),filter2))>0; end
                            if any(namematch), results=[results allnames(namematch)]; end
                        else
                            filename=fullfile(pathname,fliplr(deblank(fliplr(deblank(filternow)))));
                            names=conn_dirn(filename);
                            for n1=1:length(names),
                                if ~names(n1).isdir&&~isempty(regexp(names(n1).name,filter2)), results{end+1}=names(n1).name; end;
                            end
                        end
                        [filternow,filter]=strtok(filter,';');
                    end
                end
                if numel(results)>n0results, results(n0results+1:end)=conn_sortfilenames(results(n0results+1:end)); end
                idx=[];
                selectfile=get(h.selectfile,'string');
                if ~isempty(selectfile), 
                    idx=find(ismember(results,selectfile)); 
                    try, if numel(results)==n0results, idx=find(ismember(results,[parse{1},'   ',char(selectfile),parse{2}])); end; end
                end
                if isempty(idx), idx=1; end
                idx=unique(max(1,min(numel(results),idx)));
                set(h.files,'string',char(results),'value',idx,'listboxtop',1);
                set(h.folder,'string',pathname); %fullfile(pathname,filesep));
                set(h.folderpopup,'string',regexprep(str,'(.+)[\/\\]$','$1'),'value',numel(str));
                set(h.selectfile,'string','');
                set(h.select,'string',h.strbutton,'fontweight','bold','enable','on');
                cwd=pathname; %fullfile(pathname,filesep);
%             elseif selectfolder&&~isempty(h.callback)&&selection, % select folder
%                 tcolor=get(h.select,'backgroundcolor');
%                 set(h.select,'backgroundcolor',CONN_gui.backgroundcolor,'string','working...','fontweight','normal','enable','off');drawnow;
%                 if iscell(h.callback),
%                     if length(h.callback)>1, feval(h.callback{1},h.callback{2:end},pathname); else, feval(h.callback{1},pathname); end
%                 else, feval(h.callback,pathname); 
%                 end
%                 set(h.select,'backgroundcolor',tcolor,'string',h.strbutton,'fontweight','bold','enable','on');
            elseif (selectfolder||~isdirectory)&&~isempty(h.callback)&&selection, % select folder or select file
                if h.reduced&&~isdirectory
                    idx=1;
                    names=fliplr(deblank(fliplr(deblank(get(h.filename,'string')))));
                else
                    idx=get(h.files,'value');
                    names=get(h.files,'string');
                end
                if ~isempty(idx) & size(names,1)>=max(idx),
                    names=names(idx,:);
                    pathname=fliplr(deblank(fliplr(deblank(get(h.folder,'string')))));
                    if isempty(pathname)||(pathname(end)~='/'&&pathname(end)~='\'), pathname=[pathname,filesep]; end
                    names=regexprep(cellstr(names),'^\s+|\s+$','');
                    for n2=1:numel(names),
                        if strncmp(names{n2},parse{1},numel(parse{1})), names{n2}=names{n2}(numel(parse{1})+1:end-numel(parse{2})); end
                    end
                    names=char(regexprep(names,'^\s+|\s+$',''));
                    names=[repmat(pathname,[size(names,1),1]),names];
                    tcolor=get(h.select,'backgroundcolor');
                    set(h.select,'backgroundcolor',CONN_gui.backgroundcolor,'string','working...','fontweight','normal','enable','off');drawnow;
                    if iscell(h.callback),
                        if length(h.callback)>1, feval(h.callback{1},h.callback{2:end},names); else, feval(h.callback{1},names); end
                    else, feval(h.callback,names); end
                    try, set(h.select,'backgroundcolor',tcolor,'string',h.strbutton,'fontweight','bold','enable','on'); end
                end
            elseif ~selection&&~isdirectory, % show info
                idx=get(h.files,'value');
                names=get(h.files,'string');
                strselected=sprintf('%d files selected',numel(idx)); 
                if ~isempty(idx) & size(names,1)>=max(idx),
                    names=names(idx,:);
                    pathname=fliplr(deblank(fliplr(deblank(get(h.folder,'string')))));
                    if isempty(pathname)||(pathname(end)~='/'&&pathname(end)~='\'), pathname=[pathname,filesep]; end
                    names=[repmat(pathname,[size(names,1),1]),names];
                    try
                        if h.reduced, strselected='';
                        elseif size(names,1)>4,
                            strselected={sprintf('[%d files]',size(names,1))};
                            strselected{end+1}=['First: ',deblank(names(1,:))]; strselected{end+1}=['Last : ',deblank(names(end,:))];
                            for n1=1:length(strselected), if length(strselected{n1})>25+9, strselected{n1}=[strselected{n1}(1:4),' ... ',strselected{n1}(end-25+1:end)]; end; end; 
                        else
                            temp=conn_file(names,false);
                            if ~isempty(temp{3}), strselected=temp{2};
                            elseif ~isempty(temp{2}), strselected=temp{2};
                            else, strselected='unrecognized format';
                            end
                        end
                    catch
                        strselected='unrecognized format';
                    end
                end
                set(h.selected,'string',strselected);
            end
        case {'find'}
            if h.reduced
                if iscell(h.callback),
                    if length(h.callback)>1, feval(h.callback{1},h.callback{2:end},''); else, feval(h.callback{1},''); end
                else, feval(h.callback,''); end
            else
                state=xor(1,h.find_state);%get(h.find,'value');
                h.find_state=state;
                set(h.find,'userdata',h);
                if state,
                    results=get(h.files,'string');
                    results=results(find(results(:,1)=='<'),:);
                    resultsnew=[];
                    set(h.find,'string','Cancel');
                    pathname=fliplr(deblank(fliplr(deblank(get(h.folder,'string')))));
                    filter=get(h.filter,'string');
                    filter2=get(h.regexp,'string');
                    if strcmp(filter2,'.'), filter2=''; end
                    set(h.files,'string',resultsnew,'value',1);
                    ok=dirtree(pathname,filter,filter2,h,length(pathname));
                    resultsnew=get(h.files,'string');
                    resultsnew=conn_sortfilenames(resultsnew);
                    set(h.files,'string',strvcat(results,resultsnew));
                    if ok, h.find_state=0; end  % clears the results when clicking 'cancel' button
                    set(h.find,'value',0,'string','Find','userdata',h);
                else
                    set(h.find,'value',0,'string','Find');
                    resultsnew=get(h.files,'string');
                end
            end
    end
end
end

function ok=dirtree(pathname,filter,filter2,h,L)
persistent dcharcount
if isempty(dcharcount), dcharcount=0; end

ok=true;
if ~get(h.find,'value'), ok=false; return; end
dchar=' ...    ';
filterrest=filter;
[filternow,filterrest]=strtok(filterrest,';');
txt1=get(h.files,'string');
txt={};
while ~isempty(filternow),
    if size(txt1,1)>1e5, % Change this value to increase the maximum number of files displayed
        txt=strvcat(txt1,txt{:});
        set(h.files,'string',txt,'value',1);
        set(h.selected,'string',sprintf('(%d files found) %s',size(txt,1),dchar(ones(1,8))));
        return;
    end
    filename=fullfile(pathname,fliplr(deblank(fliplr(deblank(filternow)))));
    dir0=conn_dirn(filename);
    [names,idx]=sortrows(strvcat(dir0(:).name));
    for n1=1:length(dir0),
        if ~dir0(idx(n1)).isdir
            tfilename=fullfile(pathname(L+1:end),dir0(idx(n1)).name);
            if isempty(filter2)||~isempty(regexp(tfilename,filter2)), txt{end+1}=tfilename; end
        end
    end
    [filternow,filterrest]=strtok(filterrest,';');
end
txt=strvcat(txt1,txt{:});
set(h.files,'string',txt);
set(h.selected,'string',sprintf('(%d files found) %s',size(txt,1),dchar(mod(dcharcount+(1:8),length(dchar))+1)));
dcharcount=rem(dcharcount-1,8);
drawnow;
set(h.selected,'string',sprintf('(%d files found) %s',size(txt,1),dchar(ones(1,8))));
dir0=conn_dirn(pathname);
[names,idx]=sortrows(strvcat(dir0(:).name));
for n1=1:length(dir0),
    if dir0(idx(n1)).isdir && ~strcmp(dir0(idx(n1)).name,'.') && ~strcmp(dir0(idx(n1)).name,'..'),
        ok=dirtree(fullfile(pathname,dir0(idx(n1)).name),filter,filter2,h,L);
    end
end
end

function str=conn_filesearch_breakfolder(pathname)
idx=find(pathname=='/'|pathname=='\');
str=mat2cell(pathname,1,diff([0 idx(:)' numel(pathname)]));
% str={pathname};
% pbak='';
% while ~isempty(pathname)
%     [pathname,temp]=fileparts(pathname);
%     if isequal(pbak,pathname), str{end+1}=pathname; break; end
%     str{end+1}=pathname;
%     pbak=pathname;
% end
str=str(cellfun('length',str)>0);
end

function names=conn_filesearch_selectconn(varargin)
global CONN_x CONN_gui;
if isfield(CONN_gui,'font_offset'),font_offset=CONN_gui.font_offset; else font_offset=0; end
opt0_vals=[];
if nargin>=1&&ischar(varargin{1}), 
    [ok,iok]=ismember(varargin{1},{'*.img; *.nii; *.mgh; *.mgz; *.gz','*.img; *.nii; *.gz','*.img; *.nii; *.tal; *.mgh; *.mgz; *.annot; *.gz','*.mat; *.txt; *.par; *.1d; *.csv; *.tsv'});
    if ok, opt0_vals=iok; end
end
names={};
a.CONN_x=CONN_x;
dlg.fig=figure('units','norm','position',[.4,.3,.2,.4],'menubar','none','numbertitle','off','name','Select data from CONN project','color','w');
dlg.m0=uicontrol('style','popupmenu','units','norm','position',[.1,.90,.8,.05],'string',{'From current CONN project','From other CONN project'},'callback',@conn_filesearch_selectconn_select,'fontsize',9+font_offset,'parent',dlg.fig);
dlg.m1=uicontrol('style','popupmenu','units','norm','position',[.1,.825,.8,.05],'string',{'Structural','Functional','ROIs','Covariates'},'callback',@conn_filesearch_selectconn_update,'fontsize',9+font_offset,'parent',dlg.fig);
dlg.m2=uicontrol('style','popupmenu','units','norm','position',[.1,.75,.8,.05],'string',a.CONN_x.Setup.rois.names(1:end-1),'callback',@conn_filesearch_selectconn_update,'fontsize',9+font_offset,'visible','off','parent',dlg.fig);
dlg.m3=uicontrol('style','popupmenu','units','norm','position',[.1,.75,.8,.05],'string',a.CONN_x.Setup.l1covariates.names(1:end-1),'callback',@conn_filesearch_selectconn_update,'fontsize',9+font_offset,'visible','off','parent',dlg.fig);
dlg.m9=uicontrol('style','popupmenu','units','norm','position',[.1,.75,.8,.05],'string',arrayfun(@(n)sprintf('dataset %d',n),0:numel(a.CONN_x.Setup.secondarydataset),'uni',0),'callback',@conn_filesearch_selectconn_update,'fontsize',9+font_offset,'visible','off','parent',dlg.fig);
dlg.m4=uicontrol('style','checkbox','units','norm','position',[.1,.65,.4,.05],'value',1,'string','All subjects','backgroundcolor','w','callback',@conn_filesearch_selectconn_update,'fontsize',9+font_offset,'parent',dlg.fig);
dlg.m5=uicontrol('style','listbox','units','norm','position',[.1,.3,.4,.35],'max',2,'string',arrayfun(@(n)sprintf('Subject%d',n),1:a.CONN_x.Setup.nsubjects,'uni',0),'backgroundcolor','w','tooltipstring','Select subjects','visible','off','callback',@conn_filesearch_selectconn_update,'fontsize',9+font_offset,'parent',dlg.fig);
dlg.m6=uicontrol('style','checkbox','units','norm','position',[.5,.65,.8,.05],'value',1,'string','All sessions','backgroundcolor','w','callback',@conn_filesearch_selectconn_update,'fontsize',9+font_offset,'parent',dlg.fig);
dlg.m7=uicontrol('style','listbox','units','norm','position',[.5,.3,.4,.35],'max',2,'string',arrayfun(@(n)sprintf('Session%d',n),1:max(a.CONN_x.Setup.nsessions),'uni',0),'backgroundcolor','w','tooltipstring','Select subjects','visible','off','callback',@conn_filesearch_selectconn_update,'fontsize',9+font_offset,'parent',dlg.fig);
dlg.m8=uicontrol('style','text','units','norm','position',[.1,.20,.8,.05],'string','','horizontalalignment','center','fontsize',9+font_offset,'parent',dlg.fig);
dlg.m11=uicontrol('style','pushbutton','units','norm','position',[.55,.04,.2,.1],'string','Import','callback','uiresume(gcbf)','fontsize',9+font_offset,'parent',dlg.fig);
dlg.m12=uicontrol('style','pushbutton','units','norm','position',[.75,.04,.2,.1],'string','Cancel','callback','delete(gcbf)','fontsize',9+font_offset,'parent',dlg.fig);
if ~isempty(opt0_vals), set(dlg.m1,'value',opt0_vals); end
conn_filesearch_selectconn_update;
uiwait(dlg.fig);
if ~ishandle(dlg.fig), return; end

opt0=get(dlg.m1,'value');
if get(dlg.m4,'value'), nsubs=1:a.CONN_x.Setup.nsubjects; set(dlg.m5,'visible','off'); else set(dlg.m5,'visible','on'); nsubs=get(dlg.m5,'value'); end
if get(dlg.m6,'value'), nsessall=1:max(a.CONN_x.Setup.nsessions); set(dlg.m7,'visible','off'); else set(dlg.m7,'visible','on'); nsessall=get(dlg.m7,'value'); end
nsessmax=a.CONN_x.Setup.nsessions(min(length(a.CONN_x.Setup.nsessions),nsubs));
nfields=sum(sum(conn_bsxfun(@le,nsessall(:),nsessmax(:)')));
nroi=get(dlg.m2,'value');
ncov=get(dlg.m3,'value');
nset=get(dlg.m9,'value')-1;
delete(dlg.fig);
switch(opt0)
    case 1,
        n0=0;
        if ~a.CONN_x.Setup.structural_sessionspecific, nsess=1;
        else nsess=nsessall;
        end
        for n1=1:length(nsubs),
            nsub=nsubs(n1);
            for nses=intersect(nsess,1:nsessmax(n1))
                n0=n0+1;
                names{n0}=a.CONN_x.Setup.structural{nsub}{nses}{1};
            end
        end
    case 2,
        n0=0;
        for n1=1:length(nsubs),
            nsub=nsubs(n1);
            for nses=intersect(nsessall,1:nsessmax(n1))
                n0=n0+1;
                Vsource=a.CONN_x.Setup.functional{nsub}{nses}{1};
                if nset
                    try
                        if a.CONN_x.Setup.secondarydataset(nset).functionals_type==4
                            VsourceUnsmoothed=cellstr(a.CONN_x.Setup.secondarydataset(nset).functionals_explicit{nsub}{nses}{1});
                        else
                            Vsource1=cellstr(Vsource);
                            VsourceUnsmoothed=conn_rulebasedfilename(Vsource1,a.CONN_x.Setup.secondarydataset(nset).functionals_type,a.CONN_x.Setup.secondarydataset(nset).functionals_rule);
                        end
                        existunsmoothed=conn_existfile(VsourceUnsmoothed); %existunsmoothed=cellfun(@conn_existfile,VsourceUnsmoothed);
                        if ~all(existunsmoothed),
                            fprintf('warning: set-%d data for subject %d session %d not found. Using set-0 functional data instead for ROI extraction\n',nset,nsub,nses);
                        else
                            Vsource=char(VsourceUnsmoothed);
                        end
                    catch
                        fprintf('warning: error in CONN_x.Setup.secondarydataset for subject %d session %d not found. Using dataset-0 functional data instead for ROI extraction\n',nsub,nses);
                    end
                end
                names{n0}=Vsource; %a.CONN_x.Setup.functional{nsub}{nses}{1};
            end
        end
    case 3,
        n0=0;
        if nroi<=3, subjectspecific=1; sessionspecific=a.CONN_x.Setup.structural_sessionspecific;
        else subjectspecific=a.CONN_x.Setup.rois.subjectspecific(nroi); sessionspecific=a.CONN_x.Setup.rois.sessionspecific(nroi);
        end
        if ~subjectspecific, nsubs=1; end
        if ~sessionspecific, nsess=1;
        else nsess=nsessall;
        end
        for n1=1:length(nsubs),
            nsub=nsubs(n1);
            for nses=intersect(nsess,1:nsessmax(n1))
                n0=n0+1;
                names{n0}=a.CONN_x.Setup.rois.files{nsub}{nroi}{nses}{1};
            end
        end
    case 4,
        n0=0;
        for n1=1:length(nsubs),
            nsub=nsubs(n1);
            for nses=intersect(nsessall,1:nsessmax(n1))
                n0=n0+1;
                names{n0}=a.CONN_x.Setup.l1covariates.files{nsub}{ncov}{nses}{1};
            end
        end
end

    function conn_filesearch_selectconn_select(varargin)
        if get(dlg.m0,'value')==1, a.CONN_x=CONN_x;
        else
            [filename,pathname]=conn_fileutils('uigetfile',{'*.mat','conn-project files (conn_*.mat)';'*','All Files (*)'},'Select CONN project','conn_*.mat');
            if ischar(filename), a=conn_loadmatfile(fullfile(pathname,filename),'CONN_x'); end
        end
        set(dlg.m2,'string',a.CONN_x.Setup.rois.names(1:end-1),'value',min(numel(a.CONN_x.Setup.rois.names)-1,get(dlg.m2,'value')));
        set(dlg.m3,'string',a.CONN_x.Setup.l1covariates.names(1:end-1),'value',min(numel(a.CONN_x.Setup.l1covariates.names)-1,get(dlg.m3,'value')));
        set(dlg.m9,'string',arrayfun(@(n)sprintf('dataset %d',n),0:numel(a.CONN_x.Setup.secondarydataset),'uni',0),'value',min(numel(numel(a.CONN_x.Setup.secondarydataset))+1,get(dlg.m9,'value')));
        set(dlg.m5,'string',arrayfun(@(n)sprintf('Subject%d',n),1:a.CONN_x.Setup.nsubjects,'uni',0),'value',unique(min(a.CONN_x.Setup.nsubjects,get(dlg.m5,'value'))));
        set(dlg.m7,'string',arrayfun(@(n)sprintf('Session%d',n),1:max(a.CONN_x.Setup.nsessions),'uni',0),'value',unique(min(max(a.CONN_x.Setup.nsessions),get(dlg.m7,'value'))));
        conn_filesearch_selectconn_update;
    end

    function names=conn_filesearch_selectconn_update(varargin)
        opt0=get(dlg.m1,'value');
        if opt0==2, set([dlg.m9],'visible','on');
        else set([dlg.m9],'visible','off');
        end
        if opt0<=2, set([dlg.m2 dlg.m3],'visible','off');
        elseif opt0==3, set(dlg.m2,'visible','on');set(dlg.m3,'visible','off');
        elseif opt0==4, set(dlg.m2,'visible','off');set(dlg.m3,'visible','on');
        end
        if get(dlg.m4,'value'), nsubs=1:a.CONN_x.Setup.nsubjects; set(dlg.m5,'visible','off'); else set(dlg.m5,'visible','on'); nsubs=get(dlg.m5,'value'); end
        if get(dlg.m6,'value'), nsessall=1:max(a.CONN_x.Setup.nsessions); set(dlg.m7,'visible','off'); else set(dlg.m7,'visible','on'); nsessall=get(dlg.m7,'value'); end
        
        nsessmax=a.CONN_x.Setup.nsessions(min(length(a.CONN_x.Setup.nsessions),nsubs));
        nfields=sum(sum(conn_bsxfun(@le,nsessall(:),nsessmax(:)')));
        nroi=get(dlg.m2,'value');
        ncov=get(dlg.m3,'value');
        nset=get(dlg.m9,'value')-1;
        if opt0==1&&~a.CONN_x.Setup.structural_sessionspecific, nfields=numel(nsubs);
        elseif opt0==3
            if nroi<=3, subjectspecific=1; sessionspecific=a.CONN_x.Setup.structural_sessionspecific;
            else subjectspecific=a.CONN_x.Setup.rois.subjectspecific(nroi); sessionspecific=a.CONN_x.Setup.rois.sessionspecific(nroi);
            end
            if subjectspecific&&~sessionspecific, nfields=numel(nsubs);
            elseif ~subjectspecific&&~sessionspecific, nfields=1;
            elseif ~subjectspecific&&sessionspecific, nfields=numel(nsessall);
            end
        end
        set(dlg.m8,'string',sprintf('%d files',nfields));
    end
end


