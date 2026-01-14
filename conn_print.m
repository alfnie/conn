function conn_print(varargin)
% CONN_PRINT high-resolution print
% 
% conn_print                                    prints current figure, launches a GUI to specify filename and print options
% conn_print filename;                          specifies output filename
% conn_print(fighandle,...)                     specifies a figure handle
% conn_print(...,print_options)                 additional options to print command (e.g. -r600; see help print)
% conn_print(...,'-spm',...)                    selects SPM figure instead of current figure
% conn_print(...,'-nogui',...)                  skips GUI (uses default or last-used options)
% conn_print(...,'-nopersistent',...)           does not store selected options in "last-used options" list
% conn_print(...,'-mosaic',mosaic_commands,...) prints mosaic (multiple images; -mosaic|-column|-row|-mosaic3|-mosaic8;  before printing each sub-image the function handle mosaic_commands#n is first executed) 
%

persistent printoptions;

if isempty(printoptions), printoptions={'-djpeg90','-r600','-opengl'}; end
answ=varargin;
if ~isempty(answ)&&numel(answ{1})==1&&ishandle(answ{1})&&strcmp(get(answ{1},'type'),'figure'), hfig=answ{1}; answ=answ(2:end); 
elseif any(strcmp(answ,'-spm')), hfig=spm_figure('FindWin','Graphics');if isempty(hfig)||~ishandle(hfig), return; end; answ=answ(~strcmp(answ,'-spm'));
else hfig=gcf;
end
warning('off','MATLAB:print:CustomResizeFcnInPrint');
warning('off','MATLAB:print:ExcludesUIInFutureRelease');
%try, if ~isequal(get(hfig,'windowstate'),'normal'), set(hfig,'windowstate','normal','position',get(hfig,'position')); drawnow; end; end % note: temporal fix to avoid issues with maximized windows in latest Matlab versions
set(hfig,'inverthardcopy','off');
units=get(hfig,{'units','paperunits'});
set(hfig,'units','points');
set(hfig,'paperunits','points','paperpositionmode','manual','paperposition',get(hfig,'position'));
set(hfig,{'units','paperunits'},units);
if any(strcmp(answ,'-nogui')), dogui=false; answ=answ(~strcmp(answ,'-nogui'));
else dogui=true;
end
if any(strcmp(answ,'-noerror')), doerror=false; answ=answ(~strcmp(answ,'-noerror'));
else doerror=true;
end
if any(strcmp(answ,'-nopersistent')), persistentoptions=false; answ=answ(~strcmp(answ,'-nopersistent'));
else persistentoptions=true;
end
PRINT_OPTIONS=printoptions;
charansw=find(cellfun(@ischar,answ));
if any(cellfun('length',regexp(answ(charansw),'^-r\d+$'))), PRINT_OPTIONS{2}=answ{charansw(cellfun('length',regexp(answ(charansw),'^-r\d+$'))>0)}; answ(charansw(cellfun('length',regexp(answ(charansw),'^-r\d+$'))>0))=[]; charansw=find(cellfun(@ischar,answ)); end
if any(cellfun('length',regexp(answ(charansw),'^-ADD.+$'))), n1=charansw(cellfun('length',regexp(answ(charansw),'^-ADD.+$'))>0); PRINT_OPTIONS(end+(1:numel(n1)))=regexprep(answ(n1),'^-ADD','-'); answ(n1)=[]; end
if numel(answ)<1, answ{1}='print01.jpg'; end
charansw=find(cellfun(@ischar,answ));
if any(cellfun(@(x)any(strcmp(x,{'-mosaic','-column','-row','-mosaic3','-mosaic8'})),answ(charansw)))
    idx=charansw(find(cellfun(@(x)any(strcmp(x,{'-mosaic','-column','-row','-mosaic3','-mosaic8'})),answ(charansw)),1));
    mosaic_commands=answ(idx+1:end);
    if any(strcmp(answ,'-mosaic')), mosaic_type=1; if numel(mosaic_commands)~=4, error('incorrect number of arguments for print -mosaic option'); end
    elseif any(strcmp(answ,'-column')), mosaic_type=2; %if numel(mosaic_commands)~=4, error('incorrect number of arguments for print -column option'); end
    elseif any(strcmp(answ,'-row')), mosaic_type=3; %if numel(mosaic_commands)~=4, error('incorrect number of arguments for print -row option'); end
    elseif any(strcmp(answ,'-mosaic8')), mosaic_type=4; if numel(mosaic_commands)~=8, error('incorrect number of arguments for print -mosaic8 option'); end
    elseif any(strcmp(answ,'-mosaic3')), mosaic_type=5; if numel(mosaic_commands)~=3, error('incorrect number of arguments for print -mosaic3 option'); end
    end
    answ=answ(1:idx-1);
else
    mosaic_commands={};
end
if numel(answ)<2, answ(1+(1:numel(PRINT_OPTIONS)))=PRINT_OPTIONS; end

if dogui
%     answ=inputdlg({'Output file','Print options (see ''help print'')'},...
%         'print options',1,...
%         {answ{1},strtrim(sprintf('%s ',answ{2:end}))});
    thfig=figure('units','norm','position',[.4,.5,.25,.15],'color',1*[1 1 1],'name','CONN print','numbertitle','off','menubar','none');
    ht1a=uicontrol('style','text','units','norm','position',[.1,.8,.8,.15],'string','Save image as:','horizontalalignment','left','backgroundcolor',1*[1 1 1],'fontweight','bold');
    ht1=uicontrol('style','edit','units','norm','position',[.1,.65,.7,.15],'string',answ{1},'tooltipstring','enter full path to target image file, or click the button at the right to do this using the system dialog box','userdata',0,'callback','set(gcbo,''userdata'',0);');
    ht1b=uicontrol('style','pushbutton','units','norm','position',[.8,.65,.1,.15],'string','...','backgroundcolor',1*[1 1 1],'tooltipstring','select target image file/directory using system dialog box','callback','h=get(gcbo,''userdata'');filename=get(h,''string'');[tfilename,tpathname]=uiputfile({''*.jpg'',''JPG files (*.jpg)'';''*'',''All Files (*)''},''Save image as:'',filename);if ~isequal(tpathname,0), filename=fullfile(tpathname,tfilename); set(h,''string'',filename,''userdata'',1); end','userdata',ht1);
    ht2a=uicontrol('style','text','units','norm','position',[.1,.45,.8,.15],'string','Print options:','horizontalalignment','left','backgroundcolor',1*[1 1 1]);
    ht2=uicontrol('style','edit','units','norm','position',[.1,.3,.8,.15],'string',strtrim(sprintf('%s ',answ{2:end})),'tooltipstring','print command options (see ''help conn_print'' and ''help print'' for additional information)');%,'callback','uiresume');
    uicontrol('style','pushbutton','string','OK','units','norm','position',[.1,.01,.38,.25],'callback','uiresume');
    uicontrol('style','pushbutton','string','Cancel','units','norm','position',[.51,.01,.38,.25],'callback','delete(gcbf)');
    uicontrol(ht1);
    filename='';
    while isempty(filename),
        uiwait(thfig);
        if ~ishandle(thfig), return; end
        answ={get(ht1,'string'),get(ht2,'string')};
        skipoverwritequestion=isequal(get(ht1,'userdata'),1);
        filename=answ{1};
        if ~isempty(filename)&&conn_existfile(filename)&&~skipoverwritequestion
            tansw=conn_questdlg({'Output file already exist','Overwrite existing file?'},'','Yes','No','Yes');
            if ~strcmp(tansw,'Yes'), filename=''; end
        end
        if ~isempty(filename)
            try, conn_fileutils('emptyfile',filename); %fclose(fopen(filename,'wb'));
            catch,
                conn_msgbox('Unable to create target image file. Please check file/folder permissions and try again','',2);
                filename='';
            end
        end
    end
    delete(thfig);
    drawnow;
else
    answ={answ{1},strtrim(sprintf('%s ',answ{2:end}))};
    skipoverwritequestion=false;
end
if ~isempty(answ)
    filename=answ{1};
    isremotefile=conn_server('util_isremotefile',filename);
    if isremotefile, filename=conn_cache('new',filename); 
    else filename=conn_server('util_localfile',filename);
    end
    PRINT_OPTIONS=regexp(strtrim(answ{2}),'\s+','split');
    %if persistentoptions, printoptions=PRINT_OPTIONS; end
    oldpointer=get(hfig,'pointer');
    set(hfig,'pointer','watch');
    if ~isempty(mosaic_commands)
        hw=waitbar(0,'Printing. Please wait...');
        set(hw,'handlevisibility','off','hittest','off','color','w');
        domosaiccrop=true;
        a={};
        for n1=1:numel(mosaic_commands)
            if iscell(mosaic_commands{n1}), feval(mosaic_commands{n1}{1},[],[],mosaic_commands{n1}{2:end});
            elseif isa(mosaic_commands{n1},'function_handle'), feval(mosaic_commands{n1});
            else eval(mosaic_commands{n1});
            end
            drawnow;
            conn_print_internal(hfig,doerror,PRINT_OPTIONS{:},filename);
            b=imread(filename);
            if isa(b,'uint8'), b=double(b)/255; end
            if max(b(:))>1, b=double(b)/double(max(b(:))); end
            a{n1}=double(b);
            waitbar(n1/numel(mosaic_commands),hw);
            set(hw,'handlevisibility','off');
        end
        if domosaiccrop
            cropt_idx={};
            for n=1:numel(a)
                cropt=any(any(diff(a{n},1,2),2),3);
                cropt_idx{n,1}=max(1,sum(~cumsum(cropt))-64):size(a{n},1)-max(0,sum(~cumsum(flipud(cropt)))-64);
                cropt=any(any(diff(a{n},1,1),1),3);
                cropt_idx{n,2}=max(1,sum(~cumsum(cropt))-64):size(a{n},2)-max(0,sum(~cumsum(fliplr(cropt)))-64);
            end
        end
        switch mosaic_type
            case 1, %mosaic
                if domosaiccrop
                    cropt_idx13=union(cropt_idx{1,1},cropt_idx{3,1});
                    cropt_idx12=union(cropt_idx{1,2},cropt_idx{2,2});
                    cropt_idx24=union(cropt_idx{2,1},cropt_idx{4,1});
                    cropt_idx34=union(cropt_idx{3,2},cropt_idx{4,2});
                    a=[a{1}(cropt_idx13,cropt_idx12,:),a{3}(cropt_idx13,cropt_idx34,:);a{2}(cropt_idx24,cropt_idx12,:),a{4}(cropt_idx24,cropt_idx34,:)];
                else
                    a=[a{1},a{3};a{2},a{4}];
                end
            case 2, % col
                if domosaiccrop
                    cropt_idx1234=cropt_idx{1,2}; for n=2:numel(a), cropt_idx1234=union(cropt_idx1234,cropt_idx{n,2}); end
                    ta=[]; for n=1:numel(a), ta=cat(1,ta,a{n}(cropt_idx{n,1},cropt_idx1234,:)); end; a=ta;
                    %cropt_idx1234=union(union(union(cropt_idx{1,2},cropt_idx{2,2}),cropt_idx{3,2}),cropt_idx{4,2});
                    %a=[a{1}(cropt_idx{1,1},cropt_idx1234,:);a{2}(cropt_idx{2,1},cropt_idx1234,:);a{3}(cropt_idx{3,1},cropt_idx1234,:);a{4}(cropt_idx{4,1},cropt_idx1234,:)];
                else
                    a=cat(1,a{:});
                    %a=[a{1};a{2};a{3};a{4}];
                end
            case 3, % row
                if domosaiccrop
                    cropt_idx1234=cropt_idx{1,1}; for n=2:numel(a), cropt_idx1234=union(cropt_idx1234,cropt_idx{n,1}); end
                    ta=[]; for n=1:numel(a), ta=cat(2,ta,a{n}(cropt_idx1234,cropt_idx{n,2},:)); end; a=ta;
                    %cropt_idx1234=union(union(union(cropt_idx{1,1},cropt_idx{2,1}),cropt_idx{3,1}),cropt_idx{4,1});
                    %a=[a{1}(cropt_idx1234,cropt_idx{1,2},:),a{2}(cropt_idx1234,cropt_idx{2,2},:),a{3}(cropt_idx1234,cropt_idx{3,2},:),a{4}(cropt_idx1234,cropt_idx{4,2},:)];
                else
                    a=cat(2,a{:});
                    %a=[a{1},a{2},a{3},a{4}];
                end
            case 4,%mosaic8
                if domosaiccrop
                    a=a([1,2,3,6,7,8,4,5]); cropt_idx=cropt_idx([1,2,3,6,7,8,4,5],:); % note: resorts inputs
                    cropt_idx14=union(cropt_idx{1,1},cropt_idx{3,1});
                    cropt_idx25=union(cropt_idx{2,1},cropt_idx{4,1});
                    cropt_idx36=union(cropt_idx{3,1},cropt_idx{6,1});
                    cropt_idx123=union(union(cropt_idx{1,2},cropt_idx{2,2}),cropt_idx{3,2});
                    cropt_idx456=union(union(cropt_idx{4,2},cropt_idx{5,2}),cropt_idx{6,2});
                    cropt_idx78=union(cropt_idx{7,2},cropt_idx{8,2});
                    a={[a{1}(cropt_idx14,cropt_idx123,:);a{2}(cropt_idx25,cropt_idx123,:);a{3}(cropt_idx36,cropt_idx123,:)],[a{7}(cropt_idx{7,1},cropt_idx78,:);a{8}(cropt_idx{8,1},cropt_idx78,:)],[a{4}(cropt_idx14,cropt_idx456,:);a{5}(cropt_idx25,cropt_idx456,:);a{6}(cropt_idx36,cropt_idx456,:)]};
                else
                    a={[a{1};a{1};a{3}],[a{7};a{8}],[a{4};a{5};a{6}]};
                end
                a=[a{1} [a{2} ;repmat(a{1}(1,1,:),[size(a{1},1)-size(a{2},1),size(a{2},2)])] a{3}];
            case 5, %mosaic3
                if domosaiccrop
                    cropt_idx13=union(cropt_idx{1,1},cropt_idx{3,1});
                    cropt_idx12=union(cropt_idx{1,2},cropt_idx{2,2});
                    cropt_idx24=cropt_idx{2,1};
                    cropt_idx34=cropt_idx{3,2};
                    a=[a{1}(cropt_idx13,cropt_idx12,:),a{3}(cropt_idx13,cropt_idx34,:);a{2}(cropt_idx24,cropt_idx12,:),repmat(a{1}(cropt_idx13(1),cropt_idx12(1),:),numel(cropt_idx24),numel(cropt_idx34))];
                else
                    a=[a{1},a{3};a{2},repmat(a{1}(cropt_idx13(1),cropt_idx12(1),:),numel(cropt_idx24),numel(cropt_idx34))];
                end
        end
        imwrite(a,filename);
        delete(hw);
    else
        drawnow;
        conn_print_internal(hfig,doerror,PRINT_OPTIONS{:},filename);
    end
    if dogui
        try
            if 1
                if ispc, winopen(filename);
                else     system(['open "',filename,'"']);
                end
            else
                a=imread(filename);
                hf=figure('name',['printed file ',filename],'numbertitle','off','color','w');
                imagesc(a); ht=title(filename); set(ht,'interpreter','none'); axis equal tight; set(gca,'box','on','xtick',[],'ytick',[]); set(hf,'handlevisibility','off','hittest','off');
            end
        end
    end
    set(hfig,'pointer',oldpointer);
    if persistentoptions, printoptions=PRINT_OPTIONS; end
    if isremotefile, filename=answ{1}; conn_cache('push',filename); end
    conn_disp(['Saved file: ',filename]);
end
end

function conn_print_internal(hfig,doerror,varargin)
try,
    print(hfig,varargin{:});
catch
    warning('off','MATLAB:prnRenderer:opengl');
    warning('off','MATLAB:print:InvertHardcopyIgnoredDefaultColorUsed');
    if doerror
        print(hfig,'-noui',varargin{:});
    else
        try,
            print(hfig,'-noui',varargin{:});
        catch
            filename=varargin{end};
            conn_disp('fprintf','Unable to create figure %s. Using default icon\n',filename);
            try, conn_fileutils('filecopy',fullfile(fileparts(which(mfilename)),'conn_print_wrn.jpg'),filename); end
            %if ispc, [nill,nill]=system(sprintf('copy "%s" "%s"',fullfile(fileparts(which(mfilename)),'conn_print_wrn.jpg'),filename));
            %else     [nill,nill]=system(sprintf('cp -f ''%s'' ''%s''',fullfile(fileparts(which(mfilename)),'conn_print_wrn.jpg'),filename));
            %end
        end
    end
    warning('on','MATLAB:print:InvertHardcopyIgnoredDefaultColorUsed');
    warning('on','MATLAB:prnRenderer:opengl');
%     filename=varargin{end};
%     pos=get(hfig,'position');
%     p3=get(0,'screensize');
%     %pos(1:2)=pos(1:2)+p3(1:2)-1; % note: fix issue when connecting to external monitor/projector
%     pos(2)=p3(4)-pos(2)-pos(4);
%     rect = java.awt.Rectangle(pos(1), pos(2), pos(3), pos(4));
%     robot = java.awt.Robot;
%     jImage = robot.createScreenCapture(rect);
%     h = jImage.getHeight;
%     w = jImage.getWidth;
%     pixelsData = reshape(typecast(jImage.getData.getDataStorage, 'uint8'), 4, w, h);
%     img = cat(3, reshape(pixelsData(3, :, :), w, h)', reshape(pixelsData(2, :, :), w, h)', reshape(pixelsData(1, :, :), w, h)');
%     img = max(0,min(1, double(img)/255));
%     imwrite(img,filename);
end
end

