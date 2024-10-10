function ok = conn_datasetcopy(dataset_from,dataset_to,subjects,doexchange,doallowundo,dogui)
global CONN_x CONN_gui;
if nargin<6||isempty(dogui), dogui=false; end
if nargin<5||isempty(doallowundo), doallowundo=false; end
if nargin<4||isempty(doexchange), doexchange=false; end
if nargin<3||isempty(subjects), subjects=1:CONN_x.Setup.nsubjects; end
if nargin<2, dataset_to=[]; end % note: -1 is structural, 0 is functional, 1,2,... are secondary datasets
if nargin<1, dataset_from=[]; end % note: -1 is structural, 0 is functional, 1,2,... are secondary datasets

if ischar(dataset_from), dataset_from=conn_datasetlabel(dataset_from); end
if ischar(dataset_to), dataset_to=conn_datasetlabel(dataset_to,'add'); end

if isempty(dataset_to)||isempty(dataset_from)
    if isequal(dataset_to,-1)||isequal(dataset_from,-1), dnamestype=[-1,1:numel(CONN_x.Setup.secondarydataset)]; dnames=[{'structural data'},arrayfun(@(n)sprintf('secondary dataset #%d %s',n,regexprep(CONN_x.Setup.secondarydataset(n).label,'(.+)','($1)')),1:numel(CONN_x.Setup.secondarydataset),'uni',0)];
    elseif isequal(dataset_to,0)||isequal(dataset_from,0), dnamestype=[0,1:numel(CONN_x.Setup.secondarydataset)]; dnames=[{'functional data'},arrayfun(@(n)sprintf('secondary dataset #%d %s',n,regexprep(CONN_x.Setup.secondarydataset(n).label,'(.+)','($1)')),1:numel(CONN_x.Setup.secondarydataset),'uni',0)];
    else dnamestype=[-1,0,1:numel(CONN_x.Setup.secondarydataset)]; dnames=[{'structural data','functional data'},arrayfun(@(n)sprintf('secondary dataset #%d %s',n,regexprep(CONN_x.Setup.secondarydataset(n).label,'(.+)','($1)')),1:numel(CONN_x.Setup.secondarydataset),'uni',0)];
    end
    if isempty(dataset_to), dataset_to=1; end
    if isempty(dataset_from), dataset_from=1; end
    hfig=figure('units','norm','position',[.1,.3,.4,.3],'numbertitle','off','name','Move/Exchange datasets','menubar','none','color','w');
    uicontrol('style','frame','units','norm','position',[.0,.30,1,.70],'backgroundcolor',.9*[1 1 1],'foregroundcolor',.9*[1 1 1],'fontsize',9+CONN_gui.font_offset);    
    uicontrol('units','norm','position',[.05,.85,.9,.10],'style','text','string','Reassign all imaging files between different datasets','backgroundcolor',.9*[1 1 1],'horizontalalignment','left','fontweight','bold','fontsize',8+CONN_gui.font_offset);
    hm2a=uicontrol('units','norm','position',[.05,.7,.2,.10],'style','text','string','From:','backgroundcolor',.9*[1 1 1],'horizontalalignment','left','fontweight','normal','fontsize',8+CONN_gui.font_offset);
    hm2=uicontrol('units','norm','position',[.25,.7,.7,.10],'style','popupmenu','string',dnames,'value',find(dnamestype==dataset_from),'backgroundcolor',.9*[1 1 1],'horizontalalignment','left','fontweight','normal','fontsize',8+CONN_gui.font_offset);
    hm3a=uicontrol('units','norm','position',[.05,.55,.2,.10],'style','text','string','To:','backgroundcolor',.9*[1 1 1],'horizontalalignment','left','fontweight','normal','fontsize',8+CONN_gui.font_offset);
    hm3=uicontrol('units','norm','position',[.25,.55,.7,.10],'style','popupmenu','string',dnames,'value',find(dnamestype==dataset_to),'backgroundcolor',.9*[1 1 1],'horizontalalignment','left','fontweight','normal','fontsize',8+CONN_gui.font_offset);
    hm1=uicontrol('units','norm','position',[.05,.4,.9,.10],'style','checkbox','string','Exchange datasets','value',doexchange,'fontsize',8+CONN_gui.font_offset,'backgroundcolor',.9*[1 1 1],'tooltipstring','<HTML>Check this option to have imaging files exchanged between the two datasets<br/>Uncheck this option to have imaging files copied from the origin to the target dataset</HTML>');
    uicontrol('style','pushbutton','string','OK','units','norm','position',[.26,.01,.34,.12],'callback','uiresume');
    uicontrol('style','pushbutton','string','Cancel','units','norm','position',[.63,.01,.34,.12],'callback','delete(gcbf)');
    uiwait(hfig);
    
    if ~ishandle(hfig), ok=false; return; end
    doexchange=get(hm1,'value');      
    dataset_from=dnamestype(get(hm2,'value'));
    dataset_to=dnamestype(get(hm3,'value'));
    delete(hfig);
    dogui=true;
    %if dataset_from==dataset_to, conn_msgbox({'Source and target datasets are already the same',' ','Dataset operation canceled'},'',2); return; end
end
ok=true;
F1={};F2={};
if dogui, hmsg=conn_msgbox('Loading files... please wait','',-1); end
for isub=1:numel(subjects)
    nsub=subjects(isub);
    for nses=1:CONN_x.Setup.nsessions(min(length(CONN_x.Setup.nsessions),nsub))
        f1=conn_get_functional(nsub,nses,dataset_from);
        if isempty(f1), 
            fprintf('warning: set-%d data for subject %d session %d not defined.\n',dataset_from,nsub,nses);
            ok=false;
            break;
        end
        existf=conn_existfile(cellstr(f1)); 
        if ~all(existf),
            fprintf('warning: set-%d data for subject %d session %d not found.\n',dataset_from,nsub,nses);
            ok=false;
            break;
        elseif doexchange
            f2=conn_get_functional(nsub,nses,dataset_to);
            if isempty(f2),
                fprintf('warning: set-%d data for subject %d session %d not defined.\n',dataset_to,nsub,nses);
                ok=false;
                break;
            end
            existf=conn_existfile(cellstr(f2)); 
            if ~all(existf),
                fprintf('warning: set-%d data for subject %d session %d not found.\n',dataset_to,nsub,nses);
                ok=false;
                break;
            else
                F1{nsub}{nses}=f1;
                F2{nsub}{nses}=f2;
            end
        else
            F1{nsub}{nses}=f1;
        end
    end
    if ~ok, break; end
end
if dataset_to<0, strset='structural';
elseif dataset_to==0, strset='functional';
else strset=sprintf('Set-%d',dataset_to);
end
if ~ok||(dogui&&doallowundo&&strcmp(conn_questdlg(sprintf('%s volumes successfully reassigned',strset),'','Ok','Undo','Ok'),'Undo')),
    fprintf('dataset operation canceled\n');
else
    for isub=1:numel(subjects)
        nsub=subjects(isub);
        for nses=1:CONN_x.Setup.nsessions(min(length(CONN_x.Setup.nsessions),nsub))
            conn_set_functional(nsub,nses,dataset_to,F1{nsub}{nses});
            if doexchange
                conn_set_functional(nsub,nses,dataset_from,F2{nsub}{nses});
            end
        end
    end
end
if dogui&&ishandle(hmsg), delete(hmsg); end



