function Ok=conn_importl2covariate(name,y,dogui,validsubjects,descrip)
global CONN_x CONN_gui;
Ok=false;
if ~isfield(CONN_gui,'font_offset'), conn_font_init; end
%if ~isfield(CONN_x,'filename')||isempty(CONN_x.filename), conn_msgbox('No CONN toolbox project loaded','',2); end
if numel(name)~=numel(y), error('Inconsistent number of names/values'); end
if nargin<3||isempty(dogui), dogui=true; end
if nargin<4||isempty(validsubjects), validsubjects=1:CONN_x.Setup.nsubjects; end
if nargin<5||isempty(descrip), descrip=repmat({''},1,numel(name)); end

ny=cellfun(@numel,y);
if any(ny~=CONN_x.Setup.nsubjects)
    if any(rem(ny,CONN_x.Setup.nsubjects))
        error('Second-level covariates have incorrect number of subjects. Please re-run second-level analyses and import these values again');
    else
        y2={};
        name2={};
        for n=1:numel(y)
            if ny(n)~=CONN_x.Setup.nsubjects
                if isequal(size(y{n}),[1 CONN_x.Setup.nsubjects]), y{n}=y{n}.'; 
                else y{n}=reshape(y{n},CONN_x.Setup.nsubjects,[]);
                end
                for n1=1:size(y{n},2)
                    y2{end+1}=y{n}(:,n1);
                    name2{end+1}=[name{n},' measure',num2str(n1)];
                end
            else
                y2{end+1}=y{n};
                name2{end+1}=name{n};
            end
        end
        y=y2;
        name=name2;
    end
end
if dogui
    Ok=false;
    ok=true;
    ASKSELECT=true;
    if ASKSELECT
        thfig=dialog('units','norm','position',[.3,.4,.4,.4],'windowstyle','normal','name','Import 2nd-level covariate','color','w','resize','on');
        uicontrol(thfig,'style','text','units','norm','position',[.1,.85,.8,.10],'string','Select 2nd-level covariates to import','backgroundcolor','w','fontsize',9+CONN_gui.font_offset,'fontweight','bold');
        ht1=uicontrol(thfig,'style','listbox','units','norm','position',[.1,.25,.8,.6],'max',2,'string',name,'value',1:numel(name),'fontsize',8+CONN_gui.font_offset,'horizontalalignment','left');
        uicontrol(thfig,'style','pushbutton','string','Import','units','norm','position',[.1,.01,.38,.10],'callback','uiresume','fontsize',8+CONN_gui.font_offset);
        uicontrol(thfig,'style','pushbutton','string','Cancel','units','norm','position',[.51,.01,.38,.10],'callback','delete(gcbf)','fontsize',8+CONN_gui.font_offset);
        hc1=uicontextmenu('parent',thfig);
        uimenu(hc1,'label','edit taget variable names','callback',@conn_importl2covariate_changenames);
        set(ht1,'uicontextmenu',hc1);        
    else
        thfig=dialog('units','norm','position',[.3,.4,.4,.4],'windowstyle','normal','name','Import 2nd-level covariate','color','w','resize','on');
        uicontrol(thfig,'style','text','units','norm','position',[.1,.85,.8,.10],'string',sprintf('New 2nd-level covariate names (%d)',numel(name)),'backgroundcolor','w','fontsize',9+CONN_gui.font_offset,'fontweight','bold');
        ht1=uicontrol(thfig,'style','edit','units','norm','position',[.1,.25,.8,.6],'max',2,'string',name,'fontsize',8+CONN_gui.font_offset,'horizontalalignment','left');
        uicontrol(thfig,'style','text','string','note: changes to 2nd-level covariates are temporary until your project is saved','units','norm','position',[.1,.12,.8,.10],'backgroundcolor','w','fontsize',8+CONN_gui.font_offset);
        uicontrol(thfig,'style','pushbutton','string','Import','units','norm','position',[.1,.01,.38,.10],'callback','uiresume','fontsize',8+CONN_gui.font_offset);
        uicontrol(thfig,'style','pushbutton','string','Cancel','units','norm','position',[.51,.01,.38,.10],'callback','delete(gcbf)','fontsize',8+CONN_gui.font_offset);
    end
    while ok
        uiwait(thfig);
        ok=ishandle(thfig);
        if ok,
            if ASKSELECT
                newname=name;
                selected=get(ht1,'value');
            else
                newname=get(ht1,'string');
                selected=1:numel(name);
            end
            if numel(newname)~=numel(name), conn_msgbox(sprintf('Number of variable names entered (%d) does not match expected value (%d)',numel(newname),numel(name)),'',2);
            elseif numel(unique(newname(selected)))~=numel(newname(selected))||any(ismember(newname(selected),CONN_x.Setup.l2covariates.names(1:end-1))), conn_msgbox('All covariate names must be unique and different than existing second-level covariate names','',2);
            else
                for n1=1:numel(selected)
                    n=selected(n1);
                    icov=numel(CONN_x.Setup.l2covariates.names);
                    CONN_x.Setup.l2covariates.names{icov}=newname{n};
                    CONN_x.Setup.l2covariates.descrip{icov}=['CONN imported values ',newname{n}];
                    CONN_x.Setup.l2covariates.names{icov+1}=' ';
                    for n1=1:CONN_x.Setup.nsubjects, CONN_x.Setup.l2covariates.values{n1}{icov}=nan; end
                    for n1=validsubjects(:)', CONN_x.Setup.l2covariates.values{n1}{icov}=y{n}(n1); end
                    %for n1=1:CONN_x.Setup.nsubjects, CONN_x.Setup.l2covariates.values{n1}{icov}=y{n}(n1); end
                end
                delete(thfig);
                conn_msgbox({sprintf('%d new second-level covariates added',numel(selected)),'(see Tools->Calculator to explore their values)'},'',1);
                Ok=true;
                ok=false;
            end
        end
    end
else
    Ok=true;
    for n=1:numel(name)
        icov=find(strcmp(name{n},CONN_x.Setup.l2covariates.names(1:end-1)),1);
        if isempty(icov), 
            icov=numel(CONN_x.Setup.l2covariates.names); 
            CONN_x.Setup.l2covariates.names{icov}=name{n};
            CONN_x.Setup.l2covariates.names{icov+1}=' ';
            CONN_x.Setup.l2covariates.descrip{icov}=descrip{n};
            for n1=1:CONN_x.Setup.nsubjects, CONN_x.Setup.l2covariates.values{n1}{icov}=nan; end
            %validsubjects=1:CONN_x.Setup.nsubjects;
        elseif ~isempty(descrip{n})
            CONN_x.Setup.l2covariates.descrip{icov}=descrip{n};
        end
        for n1=validsubjects(:)', CONN_x.Setup.l2covariates.values{n1}{icov}=y{n}(n1); end
        conn_disp('fprintf','Updated 2nd-level covariate %s\n',name{n});
    end
end

    function conn_importl2covariate_changenames(varargin)
        int_thfig=dialog('units','norm','position',[.7,.4,.3,.4],'windowstyle','normal','name','Import 2nd-level covariate','color','w','resize','on');
        uicontrol(int_thfig,'style','text','units','norm','position',[.1,.85,.8,.10],'string',sprintf('New 2nd-level covariate names (%d)',numel(name)),'backgroundcolor','w','fontsize',9+CONN_gui.font_offset,'fontweight','bold');
        int_ht1=uicontrol(int_thfig,'style','edit','units','norm','position',[.1,.25,.8,.6],'max',2,'string',name,'fontsize',8+CONN_gui.font_offset,'horizontalalignment','left');
        uicontrol(int_thfig,'style','pushbutton','string','Ok','units','norm','position',[.1,.01,.38,.10],'callback','uiresume','fontsize',8+CONN_gui.font_offset);
        uicontrol(int_thfig,'style','pushbutton','string','Cancel','units','norm','position',[.51,.01,.38,.10],'callback','delete(gcbf)','fontsize',8+CONN_gui.font_offset);
        int_ok=true;
        while int_ok
            uiwait(int_thfig);
            int_ok=ishandle(int_thfig);
            if int_ok,
                int_name=get(int_ht1,'string');
                int_name=int_name(cellfun('length',int_name)>0);
                if numel(name)==numel(int_name),
                    name=int_name;
                    set(ht1,'string',name);
                    delete(int_thfig)
                    int_ok=false;
                end
            end
        end
    end
end
