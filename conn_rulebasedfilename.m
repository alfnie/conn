function [filename,ok]=conn_rulebasedfilename(filename,option,rule,filenamestotest,namestring,asktochange)
% CONN_RULEDBASEDFILENAME
% rule-based definition of functional data filenames
%

global CONN_gui CONN_x;

default_rule={{1,'^(s(?!ub-))+',''},{1,'^s+w',''},{1,'.*','s$0'},{1,'.*','sw$0'},{2,'func_smoothed','func_normalized'}}; % same as functional without leading 's' character (SPM convention for unsmoothed volume filenames)
if nargin<3||isempty(rule), rule=default_rule{1}; 
elseif ~iscell(rule), rule=default_rule{rule}; 
end
if nargin<4, filenamestotest=[]; end
if nargin<5||isempty(namestring), namestring='functional'; end
if nargin<6||isempty(asktochange), asktochange=false; end
if ~isfield(CONN_gui,'font_offset'), conn_font_init; end

ok=true;
if nargin<2
    if nargin<1||isempty(filename), filename='functional'; end
    switch(filename)
        case 'functional'
            filename=CONN_x.Setup.functional{1}{1}{1};
            [rule,ok]=conn_rulebasedfilename(filename,0,[],CONN_x.Setup.functional,'functional',true);
            filename={};
            if isequal(ok,true)
                ht=conn_msgbox('Changing Primary functional dataset file references. Please wait...','');
                for nsub=1:CONN_x.Setup.nsubjects
                    nsess=CONN_x.Setup.nsessions(min(length(CONN_x.Setup.nsessions),nsub));
                    for nses=1:nsess
                        Vsource=CONN_x.Setup.functional{nsub}{nses}{1};
                        Vsource1=cellstr(Vsource);
                        Vsource2=conn_rulebasedfilename(Vsource1,3,rule);
                        CONN_x.Setup.functional{nsub}{nses}=conn_file(char(Vsource2));
                        filename{nsub}{nses}=char(Vsource2);
                    end
                end
                if ishandle(ht),delete(ht); end
            end
        case 'structural'
            filename=CONN_x.Setup.structural{1}{1}{1};
            tfilename=CONN_x.Setup.structural;
            if ~CONN_x.Setup.structural_sessionspecific, tfilename=cellfun(@(x)x(1),tfilename,'uni',0); end
            [rule,ok]=conn_rulebasedfilename(filename,0,[],tfilename,'structural',true);
            filename={};
            if isequal(ok,true)
                ht=conn_msgbox('Changing structural file references. Please wait...','');
                for nsub=1:CONN_x.Setup.nsubjects
                    nsess=CONN_x.Setup.nsessions(min(length(CONN_x.Setup.nsessions),nsub));
                    if CONN_x.Setup.structural_sessionspecific, 
                        for nses=1:nsess
                            Vsource=CONN_x.Setup.structural{nsub}{nses}{1};
                            Vsource1=cellstr(Vsource);
                            Vsource2=conn_rulebasedfilename(Vsource1,3,rule);
                            CONN_x.Setup.structural{nsub}{nses}=conn_file(char(Vsource2));
                            filename{nsub}{nses}=char(Vsource2);
                        end
                    else 
                        Vsource=CONN_x.Setup.structural{nsub}{1}{1};
                        Vsource1=cellstr(Vsource);
                        Vsource2=conn_rulebasedfilename(Vsource1,3,rule);
                        for nses=1:nsess
                            CONN_x.Setup.structural{nsub}{nses}=conn_file(char(Vsource2));
                            filename{nsub}{nses}=char(Vsource2);
                        end
                    end
                end
                if ishandle(ht),delete(ht); end
            end
        case arrayfun(@(n)sprintf('dataset%d',n),0:numel(CONN_x.Setup.secondarydataset),'uni',0)
            nset=str2double(regexp(filename,'\d+$','match','once'));
            Vsource={};
            for nsub=1:CONN_x.Setup.nsubjects
                nsess=CONN_x.Setup.nsessions(min(length(CONN_x.Setup.nsessions),nsub));
                for nses=1:nsess
                    Vsource{nsub}{nses}{1}=conn_get_functional(nsub,nses,nset);
                end
            end
            filename=conn_get_functional(1,1,nset);
            [rule,ok]=conn_rulebasedfilename(filename,0,[],Vsource,'functional',true);
            filename={};
            if isequal(ok,true)
                ht=conn_msgbox(sprintf('Changing dataset #%d file references. Please wait...',nset),'');
                for nsub=1:CONN_x.Setup.nsubjects
                    nsess=CONN_x.Setup.nsessions(min(length(CONN_x.Setup.nsessions),nsub));
                    for nses=1:nsess
                        Vsource1=cellstr(Vsource{nsub}{nses}{1});
                        filename{nsub}{nses}=char(conn_rulebasedfilename(Vsource1,3,rule));
                    end
                end
                for nsub=1:CONN_x.Setup.nsubjects
                    nsess=CONN_x.Setup.nsessions(min(length(CONN_x.Setup.nsessions),nsub));
                    for nses=1:nsess
                        conn_set_functional(nsub,nses,nset,filename{nsub}{nses}); 
                    end
                end
                if ishandle(ht),delete(ht); end
            end
        case arrayfun(@(n)sprintf('l1covariate%d',n),1:numel(CONN_x.Setup.l1covariates.names)-1,'uni',0)
            ncov=str2double(regexp(filename,'\d+$','match','once'));
            filename=CONN_x.Setup.l1covariates.files{1}{ncov}{1}{1};
            files={};for n=1:numel(CONN_x.Setup.l1covariates.files), files{n}=CONN_x.Setup.l1covariates.files{n}{ncov}; end
            [rule,ok]=conn_rulebasedfilename(filename,0,[],files,'covariate',true);
            filename={};
            if isequal(ok,true)
                ht=conn_msgbox('Changing first-level covariate file references. Please wait...','');
                for nsub=1:CONN_x.Setup.nsubjects
                    nsess=CONN_x.Setup.nsessions(min(length(CONN_x.Setup.nsessions),nsub));
                    for nses=1:nsess
                        Vsource=CONN_x.Setup.l1covariates.files{nsub}{ncov}{nses}{1};
                        Vsource1=cellstr(Vsource);
                        Vsource2=conn_rulebasedfilename(Vsource1,3,rule);
                        CONN_x.Setup.l1covariates.files{nsub}{ncov}{nses}=conn_file(char(Vsource2));
                        filename{nsub}{nses}=char(Vsource2);
                    end
                end
                if ishandle(ht),delete(ht); end
            end
        case arrayfun(@(n)sprintf('roi%d',n),1:numel(CONN_x.Setup.rois.names)-1,'uni',0)
            nroi=str2double(regexp(filename,'\d+$','match','once'));
            filename=CONN_x.Setup.rois.files{1}{nroi}{1}{1};
            files={};for n=1:numel(CONN_x.Setup.rois.files), files{n}=CONN_x.Setup.rois.files{n}{nroi}; end
            [rule,ok]=conn_rulebasedfilename(filename,0,[],files,'ROI',true);
            filename={};
            if isequal(ok,true)
                ht=conn_msgbox('Changing ROI file references. Please wait...','');
                for nsub=1:CONN_x.Setup.nsubjects
                    nsess=CONN_x.Setup.nsessions(min(length(CONN_x.Setup.nsessions),nsub));
                    for nses=1:nsess
                        Vsource=CONN_x.Setup.rois.files{nsub}{nroi}{nses}{1};
                        Vsource1=cellstr(Vsource);
                        Vsource2=conn_rulebasedfilename(Vsource1,3,rule);
                        CONN_x.Setup.rois.files{nsub}{nroi}{nses}=conn_file(char(Vsource2));
                        filename{nsub}{nses}=char(Vsource2);
                    end
                end
                if ishandle(ht),delete(ht); end
            end
    end
    return
end

switch option
    case 0, % gui edit rules
        if isempty(filename), filename='functional'; end
        if ischar(filename), filename=cellstr(filename); end
        filename=filename{1};
        hfig=figure('units','norm','position',[.1,.3,.5,.4],'numbertitle','off','name','Find related files/datasets','menubar','none','color','w');
        uicontrol('style','frame','units','norm','position',[.0,.6,1,.4],'backgroundcolor',.9*[1 1 1],'foregroundcolor',.9*[1 1 1],'fontsize',9+CONN_gui.font_offset);
        
        uicontrol('units','norm','position',[.05,.85,.25,.07],'style','text','string','Find files by changing:','backgroundcolor',.9*[1 1 1],'horizontalalignment','left','fontweight','bold','fontsize',8+CONN_gui.font_offset);
        %uicontrol('units','norm','position',[.05,.75,.05,.07],'style','text','string','In','backgroundcolor',.9*[1 1 1],'horizontalalignment','left','fontweight','normal','fontsize',8+CONN_gui.font_offset);
        hm1=uicontrol('units','norm','position',[.05,.75,.25,.07],'style','popupmenu','string',{'In file name','In absolute path name'},'value',rule{1},'fontsize',8+CONN_gui.font_offset,'backgroundcolor',.9*[1 1 1]);
        hm2a=uicontrol('units','norm','position',[.30,.75,.20,.07],'style','popupmenu','string',{'remove','add','find','regexprep'},'backgroundcolor',.9*[1 1 1],'horizontalalignment','left','fontweight','normal','fontsize',8+CONN_gui.font_offset);
        hm2=uicontrol('units','norm','position',[.50,.75,.15,.07],'style','edit','string',rule{2},'backgroundcolor',.9*[1 1 1],'fontsize',8+CONN_gui.font_offset,'max',2,'tooltipstring',conn_menu_formathtml('<HTML>Search-string pattern (using regular expressions)<br/> - use * for 0 or more occurrences<br/> - use + for 1 or more occurrences<br/> - use . for any character<br/> - use ^ for beginning of string<br/> - use $ for end of string<br/> - see HELP REGEXP for advanced options</HTML>'));
        hm3a=uicontrol('units','norm','position',[.67,.75,.25,.07],'style','popupmenu','string',{'from beginning','from end'},'backgroundcolor',.9*[1 1 1],'horizontalalignment','left','fontweight','normal','fontsize',8+CONN_gui.font_offset);
        hm3b=uicontrol('units','norm','position',[.67,.75,.25,.07],'style','popupmenu','string',{'to beginning','to end'},'backgroundcolor',.9*[1 1 1],'horizontalalignment','left','fontweight','normal','fontsize',8+CONN_gui.font_offset);
        hm3c=uicontrol('units','norm','position',[.67,.75,.12,.07],'style','text','string','replace with','backgroundcolor',.9*[1 1 1],'horizontalalignment','left','fontweight','normal','fontsize',8+CONN_gui.font_offset);
        hm3=uicontrol('units','norm','position',[.79,.75,.15,.07],'style','edit','string',rule{3},'backgroundcolor',.9*[1 1 1],'fontsize',8+CONN_gui.font_offset,'max',2,'tooltipstring','Replacement-string patterns (see HELP REGEXPREP for advanced options)');
        hm6=uicontrol('units','norm','position',[.05,.60,.9,.07],'style','popupmenu','string',{'<HTML><i>example suggestions</i></HTML>','to remove leading ''s'' (SPM convention for unsmoothed volumes)','to remove leading ''sw'' (SPM convention for subject-space unsmoothed volumes)','to prepend ''s'' character (SPM convention for additionally-smoothed volumes)','to prepend ''sw'' characters (SPM convention for additionally-normalized&smoothed volumes)','to change folder'},'value',1,'fontsize',8+CONN_gui.font_offset,'backgroundcolor',.9*[1 1 1]);
        uicontrol('units','norm','position',[.05,.50,.9,.05],'style','text','string','If original data points to this file:','backgroundcolor','w','horizontalalignment','left','fontweight','bold','fontsize',8+CONN_gui.font_offset);
        hm4=uicontrol('units','norm','position',[.05,.40,.9,.1],'style','edit','string',filename,'backgroundcolor','w','horizontalalignment','left','fontsize',8+CONN_gui.font_offset,'tooltipstring','example of original filename');
        uicontrol('units','norm','position',[.05,.30,.9,.05],'style','text','string','I would like it instead to point to this file:','backgroundcolor','w','horizontalalignment','left','fontweight','bold','fontsize',8+CONN_gui.font_offset);
        hm5=uicontrol('units','norm','position',[.05,.20,.9,.1],'style','text','string','','backgroundcolor','w','horizontalalignment','left','fontsize',8+CONN_gui.font_offset,'tooltipstring','example of new filename (after string search/replacement)');
        uicontrol('style','pushbutton','string','OK','units','norm','position',[.26,.01,.34,.12],'callback','uiresume');
        uicontrol('style','pushbutton','string','Cancel','units','norm','position',[.63,.01,.34,.12],'callback','delete(gcbf)');
        set([hm6,hm1,hm2,hm3,hm4,hm2a,hm3a,hm3b,hm3c],'callback',@conn_rulebasedfilename_refresh);
        set(hm6,'visible','off'); %%% removed suggestions
        conn_rulebasedfilename_interpret(rule);
        conn_rulebasedfilename_refresh;
        uiwait(hfig);
        
        if ishandle(hfig)
            delete(hfig);
            if ~isempty(filenamestotest)&&~isempty(filenamestotest{1})&&~isempty(filenamestotest{1}{1})&&~isempty(filenamestotest{1}{1}{1})
                answ=conn_questdlg('Before applying this name-change, would you like first to double-check that the new filenames point to existing files?','','Yes, check first','No, proceed','No, proceed');
                if isequal(answ,'Yes, check first')
                    ko=0;
                    for nsub=1:numel(filenamestotest)
                        for nses=1:numel(filenamestotest{nsub})
                            temp1=cellstr(filenamestotest{nsub}{nses}{1});
                            temp1=temp1{1};
                            temp2=conn_rulebasedfilename(temp1,3,rule);
                            ttxt=sprintf('%s --> %s',temp1,temp2);
                            if conn_existfile(temp2)
                                ttxt=[ttxt sprintf('   ... Ok (file exists)\n')];
                            else
                                ttxt=[ttxt sprintf('\nWARNING!!!. File %s not found\n',temp2)];
                                ko=ko+1;
                            end
                            conn_disp('__nolog','fprintf',ttxt);
                        end
                    end
                    if ko, ok=false; conn_msgbox({sprintf('Filename convention results in several incorrect/non-existing filenames (%d)',ko),' See the command window for a full list, and modify the filename-rule if appropriate'},'ERROR!',true);
                    else
                        if asktochange
                            answ=conn_questdlg({'Filename convention results in correct/existing filenames. See the command window for a full list',sprintf('Change all %s files to these new filenames?',namestring)},'','Yes','No','Yes');
                            if isequal(answ,'Yes'), ok=true;
                            else ok=false;
                            end
                        else
                            ok=true;
                            conn_msgbox('Filename convention results in correct/existing filenames. See the command window for a full list','',true);
                        end
                    end
                end
            end
            filename=rule;
        else
            ok=false;
            filename=0;
        end
        
    case 1, % same as functional data files (do nothing)
        
    case 2, % use default programmatic rule
        filename=conn_rulebasedfilename(filename,3,{});
        
    case 3, % use user-defined programmatic rule
        waschar=false;
        if ischar(filename), waschar=true; filename=cellstr(filename); end
        if rule{1}==1,[file_path,filename,file_ext,file_num]=cellfun(@spm_fileparts,filename,'uni',0); end
        filename=regexprep(strtrim(filename),cellstr(rule{2}),cellstr(rule{3}));
        if rule{1}==1,filename=cellfun(@(a,b,c,d)fullfile(a,[b,c,d]),file_path,filename,file_ext,file_num,'uni',0); end
        if waschar, filename=char(filename); end
end

    function conn_rulebasedfilename_refresh(varargin)
        if get(hm6,'value')>1
            conn_rulebasedfilename_interpret(default_rule{get(hm6,'value')-1});
            set(hm6,'value',1);
        end
        switch(get(hm2a,'value'))
            case 1, % remove
                set([hm3a],'visible','on'); set([hm3b hm3c hm3],'visible','off'); set([hm2 hm3],'max',1);
                str1=get(hm2,'string');
                if get(hm3a,'value')==1, str1=['^',str1];
                else str1=[str1,'$'];
                end
                str2='';
            case 2, % add
                set([hm3b],'visible','on'); set([hm3a hm3c hm3],'visible','off'); set([hm2 hm3],'max',1);
                str1='.*';
                str2=get(hm2,'string');
                if get(hm3b,'value')==1, str2=[str2,'$0'];
                else str2=['$0',str2];
                end
            case 3, % find
                set([hm3c hm3],'visible','on'); set([hm3a hm3b],'visible','off'); set([hm2 hm3],'max',1);
                str1=get(hm2,'string');
                str2=get(hm3,'string');
            case 4, % find multiple
                set([hm3c hm3],'visible','on'); set([hm3a hm3b],'visible','off'); set([hm2 hm3],'max',2);
                str1=get(hm2,'string');
                str2=get(hm3,'string');
        end
        rule={get(hm1,'value'),str1,str2};
        try
            str0=get(hm4,'string');
            str=conn_rulebasedfilename(str0,3,rule);
            %if rule{1}==1, [nill,str]=fileparts(str); end
            %set(hm4,'string',str0);
            set(hm5,'string',str,'foregroundcolor','k');
        catch me
            set(hm5,'string',['-- error evaluating find/replace expression --- ' me.message],'foregroundcolor','r');
        end
    end

    function conn_rulebasedfilename_interpret(value)
        set(hm1,'value',value{1});
        if isempty(value{3}), % remove
            if ~isempty(value{2})&&value{2}(1)=='^',
                set(hm2a,'value',1);
                set(hm3a,'value',1);
                set(hm2,'string',value{2}(2:end));
            elseif ~isempty(value{2})&&value{2}(end)=='$',
                set(hm2a,'value',1);
                set(hm3a,'value',2);
                set(hm2,'string',value{2}(1:end-1));
            else
                set(hm2a,'value',3+iscell(value{2}));
                set(hm2,'string',value{2});
                set(hm3,'string',value{3});
            end
        elseif strcmp(value{2},'.*') % add
            if numel(value{3})>1&&strcmp(value{3}(end-1:end),'$0')
                set(hm2a,'value',2);
                set(hm3a,'value',1);
                set(hm2,'string',value{3}(1:end-2));
            elseif numel(value{3})>1&&strcmp(value{3}(1:2),'$0')
                set(hm2a,'value',2);
                set(hm3a,'value',2);
                set(hm2,'string',value{3}(3:end));
            else
                set(hm2a,'value',3+iscell(value{2}));
                set(hm2,'string',value{2});
                set(hm3,'string',value{3});
            end
        elseif iscell(value{2})||iscell(value{3}),
            set(hm2a,'value',4);
            set(hm2,'string',value{2});
            set(hm3,'string',value{3});
        else
            set(hm2a,'value',3);
            set(hm2,'string',value{2});
            set(hm3,'string',value{3});
        end
    end
end
