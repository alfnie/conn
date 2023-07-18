function varargout=conn_bookmark(option,varargin)
% CONN_BOOKMARK manages CONN bookmarks
%   conn_boorkmark('open',filename)
%     Opens previously saved bookmark filename
%   conn_boorkmark('edit',filename)
%     Edits properties of previously saved bookmark filename
%   conn_boorkmark('save',filename [,description,conn_args,conn_opts])
%     Saves new bookmark filename
%     Optional fields
%           description
%           conn_args   : conn(conn_args{:}) will be used to open this bookmark
%           conn_opts   : 'forcecd' cds to bookmark folder before running conn(conn_args{:}) command 
%

global CONN_x CONN_gui;
varargout={[]};

switch(option)
    case 'open'
        filename=varargin{1};
        filename=conn_prepend('',filename,'.mat');
        cwd1=pwd;
        hmsg=[]; try, hmsg=conn_msgbox({'Loading bookmark. Please wait...',regexprep(char(filename),'^.*[\\\/]','')},'',-1); end
        if ~conn_existfile(filename), error('file %s not found',filename); end
        conn_args={};opts={};
        try,
            warning('off','MATLAB:load:variableNotFound');
            conn_loadmatfile(filename,'conn_args'); 
            %warning('on','MATLAB:load:variableNotFound');
        end
        if isempty(conn_args), %back-compatibility check
            state={};
            try, conn_loadmatfile(filename,'state'); end 
            if ~isempty(state), conn_slice_display(state); end
            return;
        end
        try, 
            warning('off','MATLAB:load:variableNotFound');
            conn_loadmatfile(filename,'opts'); 
            %warning('on','MATLAB:load:variableNotFound');
        end
        if numel(conn_args)>=2&&isequal(conn_args{1},'image_display')&&isempty(fileparts(conn_args{2})),
            conn_args{2}=fullfile(fileparts(filename), conn_args{2});
        end
        if numel(conn_args)>=2&&isstruct(conn_args{2})
            if isfield(conn_args{2},'bookmark_filename'), conn_args{2}.bookmark_filename=filename; end
            if isfield(conn_args{2},'bookmark_descr'), 
                try
                    descr = conn_fileutils('fileread',conn_prepend('',filename,'.txt'));
                    descr=regexp(descr,'\n','split');
                    conn_args{2}.bookmark_descr=descr;
                end
            end
        end
        if any(strcmp(opts,'forcecd')),
            cwd2=fileparts(filename);
            if ~isempty(cwd2), try, cd(cwd2); end; end
        end
        if ishandle(hmsg), delete(hmsg); end
        conn(conn_args{:});
        if any(strcmp(opts,'forcecd')), cd(cwd1); end
        
    case 'edit'
        tfilename=varargin{1};
        [pathfilename,pfilename_name,pfilename_ext]=fileparts(tfilename);
        pfilename=[pfilename_name,pfilename_ext];          % bookmark filename
        [rootpath,pathfilename_name,pathfilename_ext]=fileparts(pathfilename);
        pathfilename=[pathfilename_name,pathfilename_ext]; % bookmark folder
        if strcmp(pathfilename,'bookmarks')&&isempty(regexp(rootpath,'bookmarks$')), rootpath=fullfile(rootpath,pathfilename); pathfilename=''; end
        if ~isempty(rootpath), filepath=rootpath; end      % root bookmarks folder
        
        tdirs=conn_dirn(fullfile(filepath,'*'));
        tdirs=tdirs([tdirs.isdir]&~ismember({tdirs.name},{'.','..'}));
        tdirs={tdirs.name};
        descr = fileread(conn_prepend('',tfilename,'.txt'));
        descr=regexp(descr,'\n','split');
        thfig=figure('units','norm','position',[.35,.5,.3,.3],'color','w','name','edit bookmark','numbertitle','off','menubar','none');
        uicontrol('style','text','units','norm','position',[.1,.8,.8,.1],'string','Description:','horizontalalignment','left','fontweight','bold','fontsize',8+CONN_gui.font_offset,'backgroundcolor','w');
        ht1=uicontrol('style','edit','units','norm','position',[.1,.5,.8,.3],'string',char(descr),'horizontalalignment','left','max',2,'fontsize',8+CONN_gui.font_offset,'backgroundcolor','w','tooltipstring','description of this bookmarked plots/results (optional)');
        uicontrol('style','text','units','norm','position',[.1,.35,.8,.1],'string','Folder:','horizontalalignment','left','fontweight','bold','fontsize',8+CONN_gui.font_offset,'backgroundcolor','w');
        uicontrol('style','text','units','norm','position',[.1,.25,.6,.1],'string','root bookmarks folder','horizontalalignment','left','backgroundcolor','w','fontsize',9+CONN_gui.font_offset);
        ht2=uicontrol('style','popupmenu','units','norm','position',[.1,.25,.6,.1],'string',[{'<HTML><i>root folder</i></HTML>'},tdirs],'backgroundcolor','w','fontsize',9+CONN_gui.font_offset,'fontweight','bold','tooltipstring','folder where this bookmark will be stored');
        ht2a=uicontrol('style','pushbutton','units','norm','position',[.75,.25,.15,.1],'string','new','fontsize',8+CONN_gui.font_offset,'tooltipstring','creates a new bookmarks folder','callback',{@conn_bookmark_update,'newfolder'});
        if ismember(pathfilename,tdirs), set(ht2,'value',1+find(strcmp(pathfilename,tdirs),1)); end
        if isempty(tdirs), set(ht2,'visible','off'); end
        uicontrol('style','pushbutton','string','Ok','units','norm','position',[.1,.01,.38,.13],'callback','uiresume','fontsize',9+CONN_gui.font_offset);
        uicontrol('style','pushbutton','string','Cancel','units','norm','position',[.51,.01,.38,.13],'callback','delete(gcbf)','fontsize',9+CONN_gui.font_offset);
        uiwait(thfig);
        if ~ishandle(thfig), varargout={false}; return; end
        descr=cellstr(get(ht1,'string'));
        if ~isempty(tdirs)&&get(ht2,'value')>1, pfilename=fullfile(tdirs{get(ht2,'value')-1},pfilename); end
        newtfilename=fullfile(filepath,pfilename);
        delete(thfig);
        if ~strcmp(tfilename,newtfilename)
            for ext={'.jpg','.mat','.txt'}
                try, conn_fileutils('movefile',conn_prepend('',tfilename,ext{1}),conn_prepend('',newtfilename,ext{1})); end
                %if ispc, [ok,nill]=system(['move "',conn_prepend('',tfilename,ext{1}),'" "',conn_prepend('',newtfilename,ext{1}),'"']);
                %else, [ok,nill]=system(['mv -f ''',conn_prepend('',tfilename,ext{1}),''' ''',conn_prepend('',newtfilename,ext{1}),'''']);
                %end
            end
        end
        descr=descr(cellfun('length',descr)>0);
        conn_fileutils('filewrite',newtfilename,descr);
        %cachefile=conn_cache('pull',conn_prepend('',newtfilename,'.txt'));
        %tfh=fopen(cachefile,'wt');
        %for n1=1:numel(descr), fprintf(tfh,'%s\n',regexprep(descr{n1},'\n','')); end
        %fclose(tfh);
        %conn_cache('push',conn_prepend('',newtfilename,'.txt'));
        varargout={true};
    case 'save'
        if numel(varargin)>=1&&~isempty(varargin{1}), pfilename=varargin{1};
        else pfilename=sprintf('%s.bookmark.jpg',datestr(now,'yyyy_mm_dd_HHMMSSFFF'));
        end
        if ~isfield(CONN_x,'folders')||~isfield(CONN_x.folders,'bookmarks')||isempty(CONN_x.folders.bookmarks), filepath=pwd;
        else filepath=CONN_x.folders.bookmarks;
        end
        descr={};
        conn_args={};
        opts={};
        if numel(varargin)>=2&&~isempty(varargin{2}), descr=varargin{2}; end
        if numel(varargin)>=3&&~isempty(varargin{3}), conn_args=varargin{3}; end
        if numel(varargin)>=4&&~isempty(varargin{4}), opts=varargin{4}; end
        [pathfilename,pfilename_name,pfilename_ext]=fileparts(pfilename);
        pfilename=[pfilename_name,pfilename_ext];          % bookmark filename
        [rootpath,pathfilename_name,pathfilename_ext]=fileparts(pathfilename);
        pathfilename=[pathfilename_name,pathfilename_ext]; % bookmark folder
        if strcmp(pathfilename,'bookmarks')&&isempty(regexp(rootpath,'bookmarks$')), rootpath=fullfile(rootpath,pathfilename); pathfilename=''; end
        if ~isempty(rootpath), filepath=rootpath; end      % root bookmarks folder
        tdirs=conn_dirn(fullfile(filepath,'*'));
        tdirs=tdirs([tdirs.isdir]&~ismember({tdirs.name},{'.','..'}));
        tdirs={tdirs.name};
        thfig=figure('units','norm','position',[.35,.5,.3,.3],'color','w','name','new bookmark','numbertitle','off','menubar','none');
        uicontrol('style','text','units','norm','position',[.1,.8,.8,.1],'string','Description:','horizontalalignment','left','fontweight','bold','fontsize',8+CONN_gui.font_offset,'backgroundcolor','w');
        ht1=uicontrol('style','edit','units','norm','position',[.1,.5,.8,.3],'string',char(descr),'horizontalalignment','left','max',2,'fontsize',8+CONN_gui.font_offset,'backgroundcolor','w','tooltipstring','description of this bookmarked plots/results (optional)');
        uicontrol('style','text','units','norm','position',[.1,.35,.8,.1],'string','Folder:','horizontalalignment','left','fontweight','bold','fontsize',8+CONN_gui.font_offset,'backgroundcolor','w');
        uicontrol('style','text','units','norm','position',[.1,.25,.6,.1],'string','root bookmarks folder','horizontalalignment','left','backgroundcolor','w','fontsize',9+CONN_gui.font_offset);
        ht2=uicontrol('style','popupmenu','units','norm','position',[.1,.25,.6,.1],'string',[{'<HTML><i>root folder</i></HTML>'},tdirs],'backgroundcolor','w','fontsize',9+CONN_gui.font_offset,'fontweight','bold','tooltipstring','folder where this bookmark will be stored');
        ht2a=uicontrol('style','pushbutton','units','norm','position',[.75,.25,.15,.1],'string','new','fontsize',8+CONN_gui.font_offset,'tooltipstring','creates a new bookmarks folder','callback',{@conn_bookmark_update,'newfolder'});
        if ismember(pathfilename,tdirs), set(ht2,'value',1+find(strcmp(pathfilename,tdirs),1)); end
        if isempty(tdirs), set(ht2,'visible','off'); end
        uicontrol('style','pushbutton','string','Ok','units','norm','position',[.1,.01,.38,.13],'callback','uiresume','fontsize',9+CONN_gui.font_offset);
        uicontrol('style','pushbutton','string','Cancel','units','norm','position',[.51,.01,.38,.13],'callback','delete(gcbf)','fontsize',9+CONN_gui.font_offset);
        uiwait(thfig);
        if ~ishandle(thfig), varargout=cell(1,3); return; end
        descr=cellstr(get(ht1,'string'));
        if ~isempty(tdirs)&&get(ht2,'value')>1, pfilename=fullfile(tdirs{get(ht2,'value')-1},pfilename); end
        tfilename=fullfile(filepath,pfilename);
        delete(thfig);
        %answer=conn_menu_inputdlg({'Plot description (optional)'},'',4,{char(descr)});
        %if isempty(answer), return; end
        %descr=cellstr(answer{1});
        descr=descr(cellfun('length',descr)>0);
        conn_fileutils('filewrite',conn_prepend('',tfilename,'.txt'),descr);
        %tfh=fopen(conn_prepend('',tfilename,'.txt'),'wt');
        %for n1=1:numel(descr), fprintf(tfh,'%s\n',descr{n1}); end
        %fclose(tfh);
        if isempty(opts), conn_savematfile(conn_prepend('',tfilename,'.mat'),'conn_args','-v7.3');
        else conn_savematfile(conn_prepend('',tfilename,'.mat'),'conn_args','opts','-v7.3');
        end
        try, conn_fileutils('copyfile',fullfile(fileparts(which('conn')),'conn_icon.jpg'),conn_prepend('',tfilename,'.jpg')); end
        %try, [nill,nill]=copyfile(fullfile(fileparts(which('conn')),'conn_icon.jpg'),conn_prepend('',tfilename,'.jpg')); end
        varargout={tfilename,pfilename,descr};
        conn_msgbox({'Bookmark saved',tfilename,'Manage your bookmarks @ Results->AllAnalyses->AllBookmarks'},'',0);
        %fprintf('Bookmark saved @ %s\nManage your bookmarks @ Results->AllAnalyses->AllBookmarks\n',tfilename);

end

    function conn_bookmark_update(hObject,eventdata,option,varargin)
        switch(option)
            case 'newfolder'
                answer=conn_menu_inputdlg({'Folder name (alphanumeric, case sensitive)'},'',1,{''});
                if isempty(answer), return; end
                tfname=answer{1};
                conn_fileutils('mkdir',filepath,tfname);
                tdirs{end+1}=tfname;
                tdirs=unique(tdirs);
                set(ht2,'string',[{'<HTML><i>root folder</i></HTML>'},tdirs],'value',1,'visible','on');
                if ismember(tfname,tdirs), set(ht2,'value',1+find(strcmp(tfname,tdirs),1)); end
        end
    end
end
