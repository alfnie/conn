function dataplot=conn_displaysubject(datafiles,consubjects)
% CONN_DISPLAYSUBJECT
% conn_displaysubject SPM.mat
%   loads second-level model and display individual-subject volumes
% conn_displaysubject(datafiles [,subjectindexes])
%   display individual-subject volumes

global CONN_gui CONN_x;
if isempty(CONN_gui)||~isfield(CONN_gui,'font_offset'), CONN_gui.font_offset=0; end

if nargin<2, consubjects=[]; end
datafiles_names={};
if nargin<1||isempty(datafiles), 
    datafiles={};
    files=conn_module('get','spm');
    if isempty(consubjects), showsubjects=1:numel(files); 
    else showsubjects=consubjects;
    end
    consubjects=[];
    for n=showsubjects(:)',
        if ~isempty(files{n}), 
            try
                SPM=struct; conn_loadmatfile(files{n},'SPM');
                if isfield(SPM,'xCon'),
                    for n1=1:numel(SPM.xCon),
                        datafiles=[datafiles {fullfile(fileparts(files{n}),SPM.xCon(n1).Vspm.fname)}];
                        consubjects=[consubjects n];
                        datafiles_names=[datafiles_names {sprintf('subject %d %s',n,SPM.xCon(n1).name)}];
                    end
                end
            catch
                temp=conn_dir(fullfile(fileparts(files{n}),'spmT*.nii'),'-cell','-R');
                if isempty(temp), temp=conn_dir(fullfile(temp{n},'spmT*.img'),'-cell','-R'); end
                if ~isempty(temp),
                    datafiles=[datafiles temp(:)'];
                    consubjects=[consubjects n+zeros(1,numel(temp))];
                    [nill,temp,nill]=cellfun(@fileparts,temp,'uni',0);
                    datafiles_names=[datafiles_names cellfun(@(x)sprintf('subject %d contrast %s',n,x),temp,'uni',0)];
                end
            end
        end
    end
    if isempty(datafiles), return; end
end
hmsg=conn_msgbox('Initializing. Please wait...','',-1);
if ischar(datafiles), datafiles=cellstr(datafiles); end
if size(datafiles,1)==1, datafiles=datafiles'; end
[file_path,file_name,file_ext]=cellfun(@fileparts,datafiles,'uni',0);
cwd=[];
if numel(file_name)==1&&strcmp([file_name{1},file_ext{1}],'SPM.mat');
    cwd=pwd;
    spmfile_path=file_path{1};
    if isempty(spmfile_path), spmfile_path='.'; end
    conn_fileutils('cd',spmfile_path);
    SPM=struct; conn_loadmatfile(datafiles{1},'SPM');
    idxsubjects=1:size(SPM.xX.X,1);
    if isfield(SPM.xX,'SelectedSubjects'), consubjects=find(SPM.xX.SelectedSubjects); consubjects=consubjects(1+rem(idxsubjects-1,numel(consubjects)));
    else consubjects=idxsubjects;
    end
    datafiles=reshape({SPM.xY.VY.fname},size(SPM.xY.VY));
    if ~conn_existfile(datafiles{1}), conn_msgbox({sprintf('Unable to find file %s',datafiles{1}),'Please re-compute second-level model and try again'},'Outdated file references',2); return; end
    if isfield(SPM.xX,'isSurface')&&SPM.xX.isSurface, dispopts={'Surface display'};
    elseif isfield(SPM.xX,'SelectedSubjects'), dispopts={'Slice display (reference anatomical)','Slice display (own subject anatomical)','Surface display','Volume display','Glass display'};
    else dispopts={'Slice display','Surface display','Volume display','Glass display'};
    end
else
    idxsubjects=1:numel(datafiles);
    a=conn_fileutils('spm_vol',datafiles{1});
    if conn_surf_dimscheck(a), dispopts={'Surface display'};
    elseif ~isempty(consubjects), dispopts={'Slice display (reference anatomical)','Slice display (own subject anatomical)','Surface display','Volume display','Glass display'};
    else dispopts={'Slice display','Surface display','Volume display','Glass display'};
    end
    spmfile_path=[];
end
if isempty(consubjects), consubjects=nan(size(idxsubjects)); end
if isempty(datafiles_names), datafiles_names=arrayfun(@(n)sprintf('subject %d',n),consubjects,'uni',0); end
if ishandle(hmsg), delete(hmsg); end
ok=true;
thfig=dialog('units','norm','position',[.3,.4,.3,.3],'windowstyle','normal','name','Plot individual subject','color','w','resize','on');
uicontrol(thfig,'style','text','units','norm','position',[.1,.85,.8,.10],'string','Display type:','horizontalalignment','left','backgroundcolor','w','fontsize',9+CONN_gui.font_offset,'fontweight','bold');
ht1=uicontrol(thfig,'style','popup','units','norm','position',[.1,.75,.8,.10],'max',2,'string',dispopts,'fontsize',8+CONN_gui.font_offset,'horizontalalignment','left','tooltipstring','select type of display');
if ~all(isnan(consubjects)), uicontrol(thfig,'style','text','units','norm','position',[.1,.6,.8,.10],'string','Subject(s):','horizontalalignment','left','backgroundcolor','w','fontsize',9+CONN_gui.font_offset,'fontweight','bold'); end
ht2=uicontrol(thfig,'style','listbox','units','norm','position',[.1,.4,.8,.20],'max',2,'string',datafiles_names,'fontsize',8+CONN_gui.font_offset,'horizontalalignment','left','tooltipstring','select individual subject(s) to display');
uicontrol(thfig,'style','text','units','norm','position',[.1,.25,.40,.10],'string','Threshold:','horizontalalignment','left','backgroundcolor','w','fontsize',9+CONN_gui.font_offset,'fontweight','bold');
ht3=uicontrol(thfig,'style','edit','units','norm','position',[.5,.25,.4,.10],'max',1,'string','.25','fontsize',8+CONN_gui.font_offset,'horizontalalignment','left','tooltipstring','select voxel-level threshold (voxels with absolute connectivity values below this threshold are not displayed / color coded; e.g. 0.25)');
uicontrol(thfig,'style','text','units','norm','position',[.1,.15,.40,.10],'string','Colormap range:','horizontalalignment','left','backgroundcolor','w','fontsize',9+CONN_gui.font_offset,'fontweight','bold');
ht4=uicontrol(thfig,'style','edit','units','norm','position',[.5,.15,.4,.10],'max',1,'string','[-1 1]','fontsize',8+CONN_gui.font_offset,'horizontalalignment','left','tooltipstring','select colorbar range (range of connectivity values color coded; e.g. [-1 1])');
uicontrol(thfig,'style','pushbutton','string','Ok','units','norm','position',[.1,.01,.38,.10],'callback','uiresume','fontsize',8+CONN_gui.font_offset);
uicontrol(thfig,'style','pushbutton','string','Cancel','units','norm','position',[.51,.01,.38,.10],'callback','delete(gcbf)','fontsize',8+CONN_gui.font_offset);
if all(isnan(consubjects)), set(ht2,'value',1:numel(consubjects),'visible','off'); end
while ok
    uiwait(thfig);
    ok=ishandle(thfig);
    if ok,
        dispopt=dispopts{get(ht1,'value')};
        thr_str=get(ht3,'string');
        thr=str2num(thr_str);
        vrange_str=get(ht4,'string');
        vrange=str2num(vrange_str);
        nidx=get(ht2,'value');
        if (isempty(thr_str)||~isempty(thr))&&(isempty(vrange_str)||~isempty(vrange))
            ok=false;
            delete(thfig);
            if isempty(thr), thr=nan; end
        end
    else dispopt=[];
    end
end
if isempty(dispopt),dataplot=[];if ~isempty(cwd),cd(cwd);end;return; end
idxsubjects=idxsubjects(nidx);
consubjects=consubjects(nidx);
switch(dispopt)
    case 'Volume display'
        fhall={};
        for n=1:numel(idxsubjects)
            fh=conn_mesh_display([],datafiles(idxsubjects(n),:),[],[],[],thr);
            set(fh('figurehandle'),'name',['Subject ',mat2str(consubjects(n))]);
            fhall{end+1}=fh;
        end
        fh=fhall;
    case 'Surface display'
        fhall={};
        for n=1:numel(idxsubjects)
            fh=conn_mesh_display(datafiles(idxsubjects(n),:),[],[],[],[],thr);
            %fh('colormap','hot');
            if ~isempty(vrange), fh('colorbar','rescale',vrange); end
            fh('colorbar','on');
            set(fh('figurehandle'),'name',['Subject ',mat2str(consubjects(n))]);
            fhall{end+1}=fh;
        end
        fh=fhall;
    case {'Slice display','Slice display (reference anatomical)'},
        fh=conn_slice_display(datafiles(idxsubjects,:),'',...
            spmfile_path,thr);
        fh('contour_transparency',1);
        fh('colormap','hot');
        fh('colorbar','rescale',vrange);
        fh('colorbar','on');
        set(fh('figurehandle'),'name',['Subject(s) ',mat2str(consubjects(:)')]);
    case 'Slice display (own subject anatomical)'
        fhall={};
        for n=1:numel(idxsubjects)
            fh=conn_slice_display(datafiles(idxsubjects(n),:),CONN_x.Setup.structural{consubjects(n)}{1}{1},...
                spmfile_path,thr);
            fh('contour_transparency',1);
            fh('colormap','hot');
            fh('colorbar','rescale',vrange);
            fh('colorbar','on');
            set(fh('figurehandle'),'name',['Subject ',mat2str(consubjects(n))]);
            fhall{end+1}=fh;
        end
        fh=fhall;
    case 'Glass display'
        fhall={};
        for n=1:numel(idxsubjects)
            fh=conn_mesh_display([],datafiles(idxsubjects(n),:),[],[],[],thr);
            set(fh('figurehandle'),'name',['Subject ',mat2str(consubjects(n))]);
            fh('brain',4);
            fh('background',[1 1 1]);
            fh('brain_transparency',0);
            fh('sub_transparency',0);
            fh('mask_transparency',.2);
            fh('material',[.1 1 1 .25 0]);
            fhall{end+1}=fh;
        end
        fh=fhall;
end
dataplot=fh;
if ~isempty(cwd), cd(cwd); end
