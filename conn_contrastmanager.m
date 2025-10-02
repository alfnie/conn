function [ok,varargout]=conn_contrastmanager(str,varargin)
% CONN_CONTRASTMANAGER manages labeled 2nd-level analyses
%
% conn_contrastmanager; 
%   displays GUI to define/edit/delete labels for 2nd-level analyses
%
% 
% names = conn_contrastmanager('names')                     : returns list of existing 2nd-level analysis labels
% names = conn_contrastmanager('namesextended')             : returns list of existing 2nd-level analysis extended labels
% conn_contrastmanager('delete',idx)                        : deletes existing 2nd-level analysis label
% conn_contrastmanager('edit',idx)                          : edits existing 2nd-level analysis label
% conn_contrastmanager('add')                               : creates new 2nd-level analysis label with information from the main CONN gui (2nd-level results tab) 
% [doexist,doexisti]=conn_contrastmanager('check')          : checks if model defined in main CONN gui (2nd-level results tab) is already labeled
% [doexist,doexisti]=conn_contrastmanager('check',subjecteffects,contrastsubjecteffects,conditions,contrastconditions) : checks if model defined by input parameters is already labeled
%

global CONN_x CONN_gui;
if ~isfield(CONN_gui,'font_offset'), conn_font_init; end

ht1=[];
forcesave=false;
if nargin
    ok=false;
    varargout={[],[]};
    switch(str)
        case 'names', 
            ok=CONN_x.Results.saved.names; 
            return
        case 'namesextended', 
            ok=conn_strexpand(CONN_x.Results.saved.names,CONN_x.Results.saved.descrip); 
            return
        case 'get'

        case {'add','delete','edit','sort'}, 
            if isempty(varargin), conn_contrastmanager_update(str,0);
            else conn_contrastmanager_update(str,varargin{:});
            end
            return
        case 'check',
            if nargin>1, [nsubjecteffects,csubjecteffects,nconditions,cconditions]=deal(varargin{:});
            else [nsubjecteffects,csubjecteffects,nconditions,cconditions]=deal(CONN_x.Results.xX.nsubjecteffects,CONN_x.Results.xX.csubjecteffects,CONN_x.Results.xX.nconditions,CONN_x.Results.xX.cconditions);
            end
            if iscell(nsubjecteffects), [ok,nsubjecteffects]=ismember(nsubjecteffects,CONN_x.Setup.l2covariates.names(1:end-1)); assert(all(ok),'unrecognized 2nd-level covariate names'); end
            if iscell(nconditions), [ok,nconditions]=ismember(nconditions,CONN_x.Setup.conditions.names(1:end-1)); assert(all(ok),'unrecognized condition names'); end
            for ncontrast=1:numel(CONN_x.Results.saved.names)
                if isequal(CONN_x.Results.saved.nsubjecteffects{ncontrast},CONN_x.Setup.l2covariates.names(nsubjecteffects))&&...
                        isequal(CONN_x.Results.saved.csubjecteffects{ncontrast},csubjecteffects)&&...
                        isequal(CONN_x.Results.saved.nconditions{ncontrast},CONN_x.Setup.conditions.names(nconditions))&&...
                        isequal(CONN_x.Results.saved.cconditions{ncontrast},cconditions),
                    ok=true;
                    varargout={ncontrast, CONN_x.Results.saved.names{ncontrast}};
                    return;
                end
            end
            return
        case 'guiadd', forcesave=true;
        otherwise, error('invalid option %s',str);
    end
end
thfig=dialog('units','norm','position',[.3,.4,.6,.3],'windowstyle','normal','name','2nd-level design manager','color','w','resize','on');
uicontrol(thfig,'style','text','units','norm','position',[.05,.85,.7,.08],'string','2nd-level designs:','backgroundcolor','w','horizontalalignment','left','fontsize',9+CONN_gui.font_offset,'fontweight','bold');
ht1=uicontrol(thfig,'style','listbox','units','norm','position',[.05,.2,.7,.65],'max',1,'string','','value',[],'fontsize',8+CONN_gui.font_offset,'tooltipstring',conn_menu_formathtml('<HTML>List of 2nd-level designs of interest <br/> - This list can be used for quick access to commonly used 2nd-level designs (combinations of between-subjects and between-conditions effects and contrasts)<br/> - To add a new design to this list first define it in the main CONN gui <i>second-level results</i> tab and then click the <i>save new contrast</i> button</HTML>'));
ht_add=uicontrol(thfig,'style','pushbutton','string','New','units','norm','position',[.8,.75,.15,.10],'callback',@(varargin)conn_contrastmanager_update('add'),'fontsize',8+CONN_gui.font_offset,'tooltipstring','Adds current 2nd-level design definition (defined in the main CONN gui) to this list');
ht_del=uicontrol(thfig,'style','pushbutton','string','Delete','units','norm','position',[.8,.65,.15,.10],'callback',@(varargin)conn_contrastmanager_update('delete'),'fontsize',8+CONN_gui.font_offset,'tooltipstring','Deletes selected 2nd-level design from this list');
ht_ren=uicontrol(thfig,'style','pushbutton','string','Edit','units','norm','position',[.8,.55,.15,.10],'callback',@(varargin)conn_contrastmanager_update('edit'),'fontsize',8+CONN_gui.font_offset,'tooltipstring','Edits selected 2nd-level design');
ht_sort=uicontrol(thfig,'style','pushbutton','string','Sort','units','norm','position',[.8,.45,.15,.10],'callback',@(varargin)conn_contrastmanager_update('sort'),'fontsize',8+CONN_gui.font_offset,'tooltipstring','Resorts alphabetically all 2nd-level designs');
%uicontrol(thfig,'style','text','string','note: changes to this list are temporary until your project is saved','units','norm','position',[.1,.15,.8,.08],'backgroundcolor','w','fontsize',8+CONN_gui.font_offset);
ht_sel=[]; %uicontrol(thfig,'style','pushbutton','string','Select','units','norm','position',[.1,.01,.38,.10],'callback','uiresume','fontsize',8+CONN_gui.font_offset,'tooltipstring','Enters selected contrast definition into CONN-gui');
ht_ext=uicontrol(thfig,'style','pushbutton','string','Close','units','norm','position',[.81,.01,.18,.10],'callback','uiresume','fontsize',8+CONN_gui.font_offset);
conn_contrastmanager_update('');
if forcesave||isempty(CONN_x.Results.saved.names), conn_contrastmanager_update('add'); end
% set(ht1 ht2],'callback',@conn_orthogonalizemenuupdate);
% conn_orthogonalizemenuupdate;
uiwait(thfig);
ok=ishandle(thfig);
if ok,
    ok=get(ht1,'value');
    delete(thfig);
end

    function conn_contrastmanager_update(option,ncontrast,name)
        if nargin<2||isempty(ncontrast), ncontrast=get(ht1,'value'); dogui=true; 
        else dogui=false; 
        end
        switch(option)
            case 'sort'
                [nill,tidx]=sort(CONN_x.Results.saved.names);
                CONN_x.Results.saved.names=CONN_x.Results.saved.names(tidx);
                CONN_x.Results.saved.labels=CONN_x.Results.saved.labels(tidx);
                CONN_x.Results.saved.descrip=CONN_x.Results.saved.descrip(tidx);
                CONN_x.Results.saved.nsubjecteffects=CONN_x.Results.saved.nsubjecteffects(tidx);
                CONN_x.Results.saved.csubjecteffects=CONN_x.Results.saved.csubjecteffects(tidx);
                CONN_x.Results.saved.nconditions=CONN_x.Results.saved.nconditions(tidx);
                CONN_x.Results.saved.cconditions=CONN_x.Results.saved.cconditions(tidx);
            case 'add'
                if nargin<3||isempty(name), 
%                     newname=conn_resultsfolder('subjectsconditions',0,CONN_x.Results.xX.nsubjecteffects,CONN_x.Results.xX.csubjecteffects,CONN_x.Results.xX.nconditions,CONN_x.Results.xX.cconditions);
                    newname='NEW_MODEL';
                    answ=conn_menu_inputdlg({'2nd-level design name (alphanumeric, case sensitive, short)','Description (arbitrarily long)'},'New 2nd-level design',1,{newname,''});
                    if isempty(answ), return; end
                    name=regexprep(answ{1},'[^\w\d ]','');
                    descrip=answ{2};
                else
                    name=regexprep(char(name),'[^\w\d ]','');
                    descrip='';
                end
                idx=1:numel(CONN_x.Results.saved.names);
                if any(strcmp(CONN_x.Results.saved.names(idx),name))
                    if nargin<3, conn_msgbox('Duplicated design name. Unable to proceed','',2); end
                    return;
                end
                if isempty(descrip), nameext=sprintf('<b>%s</b>',name);
                else nameext=sprintf('<b>%s</b> <i>(%s)</i>',name,descrip);
                end
                label=['<HTML>',nameext,' <small>Subject effects: <i>',strjoinstr(CONN_x.Setup.l2covariates.names(CONN_x.Results.xX.nsubjecteffects),' & '),'{',mat2str(CONN_x.Results.xX.csubjecteffects),'}</i> ; Conditions: <i>',strjoinstr(CONN_x.Setup.conditions.names(CONN_x.Results.xX.nconditions),' & '),'{',mat2str(CONN_x.Results.xX.cconditions),'}</i></small></HTML>'];
                label=conn_menu_formathtml(label);
                ncontrast=numel(CONN_x.Results.saved.names)+1;
                CONN_x.Results.saved.names{ncontrast}=name;
                CONN_x.Results.saved.labels{ncontrast}=label;
                CONN_x.Results.saved.descrip{ncontrast}=descrip;
                CONN_x.Results.saved.nsubjecteffects{ncontrast}=CONN_x.Setup.l2covariates.names(CONN_x.Results.xX.nsubjecteffects);
                CONN_x.Results.saved.csubjecteffects{ncontrast}=CONN_x.Results.xX.csubjecteffects;
                CONN_x.Results.saved.nconditions{ncontrast}=CONN_x.Setup.conditions.names(CONN_x.Results.xX.nconditions);
                CONN_x.Results.saved.cconditions{ncontrast}=CONN_x.Results.xX.cconditions;
%                 [ok,nill]=mkdir(CONN_x.folders.secondlevel,dirname);
            case 'delete'
                if ncontrast
                    name=CONN_x.Results.saved.names{ncontrast};
                    answ=conn_questdlg(sprintf('Delete contrast %s?',name),'','Delete','Cancel','Delete');
                    if ~isequal(answ,'Delete'), return; end
                    idx=setdiff(1:numel(CONN_x.Results.saved.names),ncontrast);
                    CONN_x.Results.saved.names=CONN_x.Results.saved.names(idx);
                    CONN_x.Results.saved.labels=CONN_x.Results.saved.labels(idx);
                    CONN_x.Results.saved.descrip=CONN_x.Results.saved.descrip(idx);
                    CONN_x.Results.saved.nsubjecteffects=CONN_x.Results.saved.nsubjecteffects(idx);
                    CONN_x.Results.saved.csubjecteffects=CONN_x.Results.saved.csubjecteffects(idx);
                    CONN_x.Results.saved.nconditions=CONN_x.Results.saved.nconditions(idx);
                    CONN_x.Results.saved.cconditions=CONN_x.Results.saved.cconditions(idx);
                end
            case 'edit'
                if ncontrast==0
                    ncontrast=numel(CONN_x.Results.saved.names)+1;
                    CONN_x.Results.saved.names{ncontrast}='NEW_MODEL';
                    CONN_x.Results.saved.labels{ncontrast}='';
                    CONN_x.Results.saved.descrip{ncontrast}='';
                    CONN_x.Results.saved.nsubjecteffects{ncontrast}=CONN_x.Setup.l2covariates.names(1);
                    CONN_x.Results.saved.csubjecteffects{ncontrast}=1;
                    CONN_x.Results.saved.nconditions{ncontrast}=CONN_x.Setup.conditions.names(1);
                    CONN_x.Results.saved.cconditions{ncontrast}=1;
                end
                if ncontrast
                    answ={};
                    if nargin<3||isempty(name),
                        answ=conn_menu_inputdlg({'2nd-level design name (alphanumeric, case sensitive, short)','Description (e.g. Does the average connectivity during rest differ from zero?)','Subject effects','Between-subjects contrast','Conditions','Between-conditions contrast'},'Edit 2nd-level design',1,{CONN_x.Results.saved.names{ncontrast},CONN_x.Results.saved.descrip{ncontrast},strjoinstr(CONN_x.Results.saved.nsubjecteffects{ncontrast},';'),mat2str(CONN_x.Results.saved.csubjecteffects{ncontrast}),strjoinstr(CONN_x.Results.saved.nconditions{ncontrast},';'),mat2str(CONN_x.Results.saved.cconditions{ncontrast})});
                        if isempty(answ), return; end
                        name=regexprep(answ{1},'[^\w\d ]','');
                    end
                    idx=setdiff(1:numel(CONN_x.Results.saved.names),ncontrast);
                    if any(strcmp(CONN_x.Results.saved.names(idx),name))
                        conn_msgbox(sprintf('Duplicated design name %s. Unable to proceed',name),'',2);
                        return;
                    end
                    if numel(answ)>1
                        CONN_x.Results.saved.descrip{ncontrast}=answ{2};
                        CONN_x.Results.saved.nsubjecteffects{ncontrast}=regexp(answ{3},'\s*;\s*','split');
                        CONN_x.Results.saved.csubjecteffects{ncontrast}=str2num(answ{4});
                        CONN_x.Results.saved.nconditions{ncontrast}=regexp(answ{5},'\s*;\s*','split');
                        CONN_x.Results.saved.cconditions{ncontrast}=str2num(answ{6});
                    end
                    if isempty(CONN_x.Results.saved.descrip{ncontrast}), nameext=sprintf('<b>%s</b>',name);
                    else nameext=sprintf('<b>%s</b> <i>(%s)</i>',name,CONN_x.Results.saved.descrip{ncontrast});
                    end
                    label=['<HTML>',nameext,' <small>Subject effects: ',strjoinstr(CONN_x.Results.saved.nsubjecteffects{ncontrast},' & '),'{',mat2str(CONN_x.Results.saved.csubjecteffects{ncontrast}),'} ; Conditions: ',strjoinstr(CONN_x.Results.saved.nconditions{ncontrast},' & '),'{',mat2str(CONN_x.Results.saved.cconditions{ncontrast}),'}</small></HTML>'];
                    label=conn_menu_formathtml(label);
                    %label=['<HTML><b>',name,'</b>',regexprep(CONN_x.Results.saved.labels{ncontrast},'<HTML><b>.*?</b>','')];
                    CONN_x.Results.saved.names{ncontrast}=name;
                    CONN_x.Results.saved.labels{ncontrast}=label;
%                 cwd=pwd;
%                 cd(CONN_x.folders.secondlevel);
%                 if ispc, [ok,nill]=system(sprintf('ren "%s" "%s"',olddirname,dirname));
%                 else 	 [ok,nill]=system(sprintf('mv ''%s'' ''%s''',olddirname,dirname));
%                 end
%                 cd(cwd);
                end
        end
        ncontrast=min([numel(CONN_x.Results.saved.names),ncontrast]);
        if ~ncontrast, ncontrast=[]; end
        labels=CONN_x.Results.saved.labels;
        try, set(ht1,'string',labels,'value',ncontrast); end
        fileresultsnames=fullfile(CONN_x.folders.secondlevel,'_list_results.mat');
        results=CONN_x.Results.saved;
        conn_savematfile(fileresultsnames,'results');
        try
            if numel(labels)>0, set([ht_del ht_ren ht_sort ht_sel],'enable','on');
            else set([ht_del ht_ren ht_sort ht_sel],'enable','off');
            end
        end
    end
end


function str=strjoinstr(str1,str2)
str=[str1(:)';repmat({str2},1,length(str1))];
str=reshape(str(1:end-1),1,numel(str)-1);
str=[str{:}];
end

function str=conn_strexpand(varargin)
if nargin<1, str={}; return; end
str=varargin{1};
changed=false(size(str));
for n1=2:nargin
    for n2=1:min(numel(str),numel(varargin{n1}))
        if ~isempty(varargin{n1}{n2}), changed(n2)=true; str{n2}=[str{n2} ' <i>(' varargin{n1}{n2} ')</i>']; end
    end
end
for n2=find(changed(:))'
    str{n2}=['<HTML>' str{n2} '</HTML>'];
end
str=conn_menu_formathtml(str);
end


