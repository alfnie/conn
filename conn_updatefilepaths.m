function var=conn_updatefilepaths(varargin)
% CONN_UPDATEFILEPATHS checks and updates filepaths if data has moved to a different folder
%
% conn_updatefilepaths
%   updates all files in current CONN structure
%   checks if files are present, and if not, tries to:
%     1) apply any existing search/replace string rule, if applicable (sucessful only if match is found AND the replacement change produces a filename that is present/exists)
%     or 2) ask the user to locate the missing file (which also implicitly defines a new search/replace rule)
%
% Search/replace rules are kept between sucessive calls to conn_updatefilepaths
% The following commands allow more direct access to search/replace rules:
%
% conn_updatefilepaths('init', root_searchstring, root_replacestring)
%   explicitly defines potential search/replace patterns
%   root_searchstring:   search string (char array or cell array)
%   root_replacestring:  replacement string strings (char array or cell array)
%   defines filename change as: [root_searchstring][body_filename] -> [root_replacestring][body_filename]
%     e.g. [/disk1/mydata][/subject01/session01/myfile.img] -> [/disk2/users/me/data][/subject01/session01/myfile.img]
%
% conn_updatefilepaths('add', root_searchstring, root_replacestring)
%    adds new search/replace pattern to existing ones
%
% conn_updatefilepaths('hold','off');
%   "forget" any potential search/replace patterns that may have been defined before
%
% conn_updatefilepaths('hold','on');
%   "remembers" from this point on any search/replace patterns defined through the GUI
%
% conn_updatefilepaths('silent', root_searchstring, root_replacestring)
%    like "conn_updatefilepaths" (no arguments) but without checking if files exist or not
%    and without user-interaction (only applying specified search/replace pattern)
%



global CONN_x CONN_gui;
persistent ht htcancel changed changes holdon clearon silent remote local hcount;

connfolder=fileparts(which(mfilename));
if isempty(CONN_gui)||~isfield(CONN_gui,'font_offset'), conn_font_init; end
if isempty(changed),changed=0;end
if isempty(holdon),holdon=0;end
if isempty(clearon),clearon=0;end
if isempty(silent),silent=0;end
if isempty(remote),remote=0;end
if isempty(local),local=0;end
if nargin==1,
    if ~isempty(ht),
        if ishandle(htcancel)&&get(htcancel,'value'), return; end
        %if ~all(ishandle(ht)), return; end
        if isempty(hcount)||hcount>=100
            hcount=0;
            conn_menu_plotmatrix('',ht(2),[1 1 10]);
            %set(ht(2),'cdata',circshift(get(ht(2),'cdata'),[0 1]));
            drawnow;
        else hcount=hcount+1;
        end
    end
    var=varargin{1};
    if ~iscell(var), return; end
    if length(var)==3 && ischar(var{1}) && iscell(var{2}) && (isstruct(var{3})||isnumeric(var{3})),
        if silent||remote||local
            if (changed|remote|local)&&~isempty(var{1})&&~strcmp(deblank(var{1}(1,:)),'[raw values]'),
                filename=var{1};
                switch(filesep),case '\',idx=find(filename=='/');case '/',idx=find(filename=='\');end; filename(idx)=filesep;
                filename=regexprep(cellstr(filename),'^\s+|\s+$','');
                if remote, [filename,changedfilename]=conn_updatefilepaths_change(filename,'remotefile');
                elseif local, [filename,changedfilename]=conn_updatefilepaths_change(filename,'localfile');
                else [filename,changedfilename]=conn_updatefilepaths_change(filename,changes);
                end
                if changedfilename
                    var{1}=char(filename);
                    if isstruct(var{3})&&isfield(var{3},'fname'),
                        for nfname=1:numel(var{3})
                            if remote, var{3}(nfname).fname=conn_updatefilepaths_change(var{3}(nfname).fname,'remotefile');
                            elseif local, var{3}(nfname).fname=conn_updatefilepaths_change(var{3}(nfname).fname,'localfile');
                            else var{3}(nfname).fname=conn_updatefilepaths_change(var{3}(nfname).fname,changes);
                            end
                        end
                    end
                    try,
                        if isempty(var{3})&&~isempty(varargin{1}{3}), var{3}=varargin{1}{3}; end
                    end
                end
            end
            if isstruct(var{3})&&isfield(var{3},'fname')&&isfield(var{3},'private'), [var{3}.private]=deal([]); end % space-saving measures
            if isnumeric(var{3})&&~isempty(regexp(var{1}(1,:),'\.tal\s*$|\.mat\s*$|\.txt\s*$|\.par\s*$|\.1d\s*$|\.csv\s*$|\.tsv\s*$')), var{3}=zeros(size(var{3},1),0); end % space-saving measures
        else
            if  ~strcmp(deblank(var{1}(1,:)),'[raw values]')&&~conn_existfile(var{1}(1,:),true),
                filename=fliplr(deblank(fliplr(deblank(var{1}))));
                iscleared=false;
                switch(filesep),case '\',idx=find(filename=='/');case '/',idx=find(filename=='\');end; filename(idx)=filesep;
                n=0; while n<size(filename,1),
                    ok=conn_existfile(filename(n+1,:),true);
                    if ok, n=n+1;
                    else
                        askthis=1;
                        if changed,
                            for nch=numel(changes):-1:1
                                if strncmp(changes(nch).key,filename(n+1,:),changes(nch).m1)
                                    filenamet=[changes(nch).fullname2(1:changes(nch).m2),filename(n+1,changes(nch).m1+1:end)];
                                    if conn_existfile(filenamet,true), ok=1; break; end
                                end
                            end
                        end
                        if ok
                            try,
                                conn_disp(['conn_updatefilepaths: updating reference ',deblank(filenamet)]);
                                fileinfo=conn_file(deblank(filenamet));
                                filename=strvcat(filename(1:n,:),filenamet,filename(n+2:end,:));
                                askthis=0;
                            end
                        end
                        if askthis,
                            fullname1=deblank(filename(n+1,:));
                            [pathname1,name1,ext1,num1]=spm_fileparts(fullname1);
                            [pathname1c,pathname1b]=fileparts(pathname1);
                            [pathname1d,pathname1c]=fileparts(pathname1c);
                            [nill,pathname1d]=fileparts(pathname1d);
                            if strcmp(pathname1c,'conn')&&strcmp(pathname1b,'rois')&&conn_existfile(fullfile(connfolder,'rois',[name1,ext1])) % automatic fixes: fix changes to conn directory
                                filename=strvcat(filename(1:n,:),fullfile(connfolder,'rois',[name1,ext1,num1]),filename(n+2:end,:));
                            elseif strcmp(pathname1d,'conn')&&strcmp(pathname1c,'utils')&&strcmp(pathname1b,'otherrois')&&conn_existfile(fullfile(connfolder,'utils','otherrois',[name1,ext1]))
                                filename=strvcat(filename(1:n,:),fullfile(connfolder,'utils','otherrois',[name1,ext1,num1]),filename(n+2:end,:));
                            elseif strcmp(pathname1d,'conn')&&strcmp(pathname1c,'utils')&&strcmp(pathname1b,'surf')&&conn_existfile(fullfile(connfolder,'utils','surf',[name1,ext1]))
                                filename=strvcat(filename(1:n,:),fullfile(connfolder,'utils','surf',[name1,ext1,num1]),filename(n+2:end,:));
                            elseif strcmp(pathname1,fullfile(connfolder,'rois'))&&conn_existfile(fullfile(connfolder,'utils','otherrois',[name1,ext1])) % automatic fixes: fix ROIs that moved from conn/rois to conn/utils/otherrois
                                filename=strvcat(filename(1:n,:),fullfile(connfolder,'utils','otherrois',[name1,ext1,num1]),filename(n+2:end,:));
                            elseif clearon
                                filename=[];
                                iscleared=true;
                            else
                                conn_disp(['conn_updatefilepaths: file ',fullname1,' not found']);
                                [name2,pathname2]=uigetfile(['*',ext1],['File not found: ',name1,ext1],['*',name1,ext1]);
                                if all(name2==0), filename=[]; break; end
                                fullname2=fullfile(pathname2,[name2,num1]);
                                fullnamematch=strvcat(fliplr(fullname1),fliplr(fullname2));
                                m=sum(cumsum(fullnamematch(1,:)~=fullnamematch(2,:))==0);
                                m1=max(0,length(fullname1)-m); m2=max(0,length(fullname2)-m);
                                changed=1;
                                changes=cat(2,changes, struct('key',filename(n+1,1:m1),'fullname2',fullname2,'m1',m1,'m2',m2));
                                filename=strvcat(filename(1:n,:),[fullname2(1:m2),filename(n+1,m1+1:end)],filename(n+2:end,:));
                                %filename=strvcat(filename(1:n,:),[repmat(fullname2(1:m2),[size(filename,1)-n,1]),filename(n+1:end,m1+1:end)]);
                            end
                        end
                    end
                end
                if ~isempty(filename),
                    %[V,str,icon]=conn_getinfo(filename);
                    var=conn_file(filename);
                    try,
                        if isempty(var{3})&&~isempty(varargin{1}{3}), var{3}=varargin{1}{3}; end
                    end
                elseif iscleared,
                    var={[],[],[]};
                else 
                    var=[]; 
                end
            else
                if ~isstruct(var{3})||(isfield(var{3},'fname')&&numel(var{3})==1)
                    try
                        filename=fliplr(deblank(fliplr(deblank(var{1}))));
                        if ~strcmp(filename,'[raw values]'),
                            %[V,str,icon]=conn_getinfo(filename);
                            var=conn_file(filename); %{filename,str,icon};
                        end
                    end
                end
                if isstruct(var{3})&&isfield(var{3},'fname')&&isfield(var{3},'private'), [var{3}.private]=deal([]); end % space-saving measures
                if isnumeric(var{3})&&~isempty(regexp(var{1}(1,:),'\.tal\s*$|\.mat\s*$|\.txt\s*$|\.par\s*$|\.1d\s*$|\.csv\s*$|\.tsv\s*$')), var{3}=zeros(size(var{3},1),0); end % space-saving measures
            end
        end
    else
        for nvar=1:length(var),
            temp=conn_updatefilepaths(var{nvar});
            if ~isempty(temp), var{nvar}=temp; elseif ~isempty(var{nvar}), var=[]; break; end
        end
    end
elseif nargin>1&&ischar(varargin{1}),
    switch(lower(varargin{1}))
        case {'init','add','silent','remotefile','localfile'}
            if nargin<3, error('Insufficient arguments. Usage conn_updatefilepaths(''init'',root_searchstring, root_replacestring)'); end
            match1=cellstr(varargin{2});
            match2=cellstr(varargin{3});
            if numel(match1)~=numel(match2), error('mismatch number of search/replace pairs. Usage conn_updatefilepaths(''init'',root_searchstring, root_replacestring)'); end
            if strcmp(lower(varargin{1}),'init')||strcmp(lower(varargin{1}),'silent')||strcmp(lower(varargin{1}),'remotefile')||strcmp(lower(varargin{1}),'localfile')||~holdon
                changed=0;
                changes=[];
            end
            holdon=false;
            for n1=1:numel(match1)
                holdon=true;
                changed=1;
                changes=cat(2,changes, struct('key',match1{n1},'fullname2',match2{n1},'m1',numel(match1{n1}),'m2',numel(match2{n1})));
            end
            if strcmp(lower(varargin{1}),'remotefile')
                silent=true;
                remote=true;
                local=false;
                conn_updatefilepaths;
            elseif strcmp(lower(varargin{1}),'localfile')
                silent=true;
                remote=false;
                local=true;
                conn_updatefilepaths;
            elseif changed&&strcmp(lower(varargin{1}),'silent')
                silent=true;
                remote=false;
                local=false;
                conn_updatefilepaths;
            end
            silent=false;
            remote=false;
            local=false;
        case 'hold'
            if ~isempty(varargin{2}), holdon=strcmp(lower(varargin{2}),'on'); end
            silent=false;
            remote=false;
            local=false;
            var=holdon;
        case 'clear'
            if ~isempty(varargin{2}), clearon=strcmp(lower(varargin{2}),'on'); end
            var=clearon;
        otherwise
            error('unknown option %s',varargin{1});
    end
else
    update={'Setup.rois.files',...
        'Setup.l1covariates.files',...
        'Setup.structural',...
        'Setup.functional',...
        'Setup.spm',...
        'Setup.dicom',...
        'Setup.bids',...
        'Setup.unwarp_functional',...       % note: obsolete (use 'vdm' secondary dataset)
        'Setup.coregsource_functional',...  % note: obsolete (use 'ref' secondary dataset)
        'Setup.explicitmask'};
    if ~holdon
        changed=0;
        changes=[];
    end
    ht=[];
    if ~silent, conn_disp('Checking if data files have been edited or moved. Please wait...'); end
    try
        if ~silent&&isfield(CONN_x,'gui')&&(isnumeric(CONN_x.gui)&&CONN_x.gui || isfield(CONN_x.gui,'display')&&CONN_x.gui.display),
            ht=dialog('units','norm','position',[.4,.5,.3,.15],'windowstyle','normal','name','','handlevisibility','on','color','w','colormap',conn_bsxfun(@min,[1 1 1],(flipud(gray(100)))));
            htcancel=uicontrol('units','norm','position',[.3 .15 .4 .2],'style','togglebutton','string','Cancel');%,'callback',@conn_updatefilepaths_stop);
            [nill,ht(2)]=conn_menu_plotmatrix('',ht,[1 1 10],[.3 .4 .4 .1]);
            %axes('units','norm','position',[.3 .4 .4 .1]);
            %ht(2)=image(max(0,.5*conn_hanning(16)*(0+50*(sin(16*pi*(0:199)/200))))); axis tight off;
            uicontrol('units','norm','position',[0 .6 1 .3],'style','text','backgroundcolor','w','string',{'Checking if data files have been edited or moved','Press ''Cancel'' to skip this step'},'fontsize',8+CONN_gui.font_offset);
            drawnow;
        end
        if ~ischar(CONN_x.filename) || ~conn_existfile(CONN_x.filename), CONN_x.filename=''; end
        for nupdate=1:length(update),
            str=regexp(update{nupdate},'\.','split');
            temp=getfield(CONN_x,str{:});
            if ~isempty(temp),
                temp=conn_updatefilepaths(temp);
                if ~isempty(temp), CONN_x=setfield(CONN_x,str{:},temp);
                else break;
                end
            end
        end
        for nalt=1:numel(CONN_x.Setup.secondarydataset)
            temp=CONN_x.Setup.secondarydataset(nalt).functionals_explicit;
            if ~isempty(temp)
                temp=conn_updatefilepaths(temp);
                if ~isempty(temp), CONN_x.Setup.secondarydataset(nalt).functionals_explicit=temp;
                else break;
                end
            end
        end
        if isfield(CONN_x,'vvAnalyses')&&isfield(CONN_x.vvAnalyses,'mask')
            for nalt=1:numel(CONN_x.vvAnalyses)
                temp=CONN_x.vvAnalyses(nalt).mask;
                if ~isempty(temp)
                    temp=conn_updatefilepaths(temp);
                    if ~isempty(temp), CONN_x.vvAnalyses(nalt).mask=temp;
                    else break;
                    end
                end
            end
        end
    catch
       conn_disp('warning: conn_updatefilepaths did not finish');
    end
    if any(ishandle(ht)), delete(ht(ishandle(ht))); drawnow; end
    ht=[];
    htcancel=[];
    var=[];
end
%     function conn_updatefilepaths_stop(varargin)
%         ht=[ht nan];
%     end
end

function [filename,changed]=conn_updatefilepaths_change(filename,changes)
waschar=false;
if ischar(filename), waschar=true; filename=cellstr(filename); end
notmatched=1:numel(filename);
if isequal(changes,'remotefile')
    notmatched=find(conn_server('util_isremotefile',filename));
    filename=conn_server('util_remotefile',filename);
elseif isequal(changes,'localfile')
    notmatched=find(~conn_server('util_isremotefile',filename));
    filename=conn_server('util_localfile',filename);
else 
    for nch=numel(changes):-1:1
        if isempty(notmatched), break; end
        ok=strncmp(changes(nch).key,filename(notmatched),changes(nch).m1);
        if any(ok),
            filename(notmatched(ok))=cellfun(@(x)[changes(nch).fullname2(1:changes(nch).m2),x(changes(nch).m1+1:end)],filename(notmatched(ok)),'uni',0);
            notmatched(ok)=[];
        end
    end
end
changed=numel(notmatched)<numel(filename);
if waschar, filename=char(filename); end
end 
    
    
