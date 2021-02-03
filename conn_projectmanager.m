function varargout = conn_projectmanager(option,varargin)
% CONN_PROJECTMANAGER 
%
% internal function: manages project access, delayed processing steps and parallelization/synchronization options
%

global CONN_x;
switch(lower(option))
    case 'null'
        pobj.isextended=false;    % is this project an extension of a different base project?
        pobj.id='';
        pobj.holdsdata=true;      % has this project an independent data folder?
        pobj.importedfiles={};    % has this project just imported associated extended projects?
        varargout={pobj};
        
    case 'extendedname', % checks if input filename indicates extended project
        filename=varargin{1};
        [jid,jidx]=regexp(filename,'\?id=(\d+)\s*(,\s*version\s*=[^,\?]+)?(,\s*subjects\s*=[^,\?]+)?(,\s*partition\s*=\d+-\d+)?\s*$','tokens','start','once');
        if ~isempty(jid),
            pobj.isextended=true;
            pobj.id=jid{1};
            if numel(jid)<2||isempty(jid{2}), pobj.ver=conn('ver'); else pobj.ver=regexprep(jid{2},'^,\s*version\s*=\s*',''); end
            if numel(jid)>=3&&~isempty(jid{3}), pobj.subjects=str2num(regexprep(jid{3},'^,\s*subjects\s*=\s*','')); pobj.holdsdata=false;
            else pobj.holdsdata=true;
            end
            if numel(jid)>=4&&~isempty(jid{4}), pobj.partition=str2double(regexp(regexprep(jid{4},'^,\s*partition\s*=\s*',''),'-','split')); end
            if ~isequal(pobj.ver,conn('ver')), error('Incorrect CONN version. Expected %s, found %s. Parallel processing only supported when all nodes are running the same version of CONN',pobj.ver,conn('ver')); end
            filename=filename(1:jidx-1);
            %if ~isfield(pobj,'subjects'), filename=conn_projectmanager('projectfile',filename,pobj); pobj.isextended=0; end % treat as normal project file
        else
            pobj=conn_projectmanager('null');
        end
        varargout={filename, pobj};

    case 'projectfile' % returns local filename associated with project
        if nargin>1, basefilename=varargin{1};
        else basefilename=CONN_x.filename;
        end
        if nargin>2, pobj=varargin{2};
        else pobj=CONN_x.pobj;
        end
        if nargin>3, fext=varargin{3};
        else fext='.dmat';
        end
        if ~pobj.isextended, localfilename=basefilename; 
        else localfilename=conn_prepend('',basefilename,['.' pobj.id fext]); 
        end
        varargout={localfilename,basefilename};
        
    case 'parentfile' % returns parent project filename associated with project
        if nargin>1, localfilename=varargin{1};
        else localfilename=CONN_x.filename;
        end
        if nargin>2, pobj=varargin{2};
        else pobj=CONN_x.pobj;
        end
        if nargin>3, fext=varargin{3};
        else fext='.dmat';
        end
        if ~pobj.isextended, basefilename=localfilename; 
        else basefilename=regexprep(localfilename,['\.' pobj.id '\' fext '$'],'.mat'); 
        end
        varargout={basefilename,localfilename};
        
    case 'ispending'
        if CONN_x.pobj.holdsdata
            localfilename=conn_projectmanager('projectfile',CONN_x.filename,struct('id','*','isextended',true));
            allfiles=conn_dir(localfilename,'-R'); % check .dmat
            varargout={~isempty(allfiles)};
        else
            varargout={false};
        end
            
    case 'updateproject' % merges any extended projects (run after loading base project)
        if nargin>1, dogui=varargin{1};
        else dogui=true;
        end
        if CONN_x.pobj.holdsdata
            localfilename=conn_projectmanager('projectfile',CONN_x.filename,struct('id','*','isextended',true));
            allfiles=conn_dir(localfilename,'-R'); % check .dmat
            allfiles=cellstr(allfiles);allfiles=char(allfiles(cellfun('length',regexp(cellstr(allfiles),'\d{4}(\d+)\.dmat$'))>0));
            alllogs={};
            if ~isempty(allfiles),
                tag=regexp(cellstr(allfiles),'\d{4}(\d+)\.dmat$','tokens','once');
                [utag,nill,itag]=unique([tag{:}]);
                vtag=true(size(utag));
                for n=1:numel(utag)
                    pathname=fullfile(conn_prepend('',CONN_x.filename,'.qlog'),utag{n});
                    if exist(pathname,'dir')&&conn_existfile(fullfile(pathname,'info.mat'))
                        load(fullfile(pathname,'info.mat'),'info'); % look at associated .qlog folders
                        info=conn_jobmanager('statusjob',info,[],true);
                        validlabels={'finished','canceled'}; %{'finished','stopped'};
                        vtag(n)=all(ismember(info.tagmsg,validlabels));
                        alllogs{n}=info.stdlog;
                    end
                end
                ivtag=find(vtag,1);
                if isempty(ivtag)
                    if dogui&&(isequal(CONN_x.gui,1)||(isstruct(CONN_x.gui)&&isfield(CONN_x.gui,'display')&&CONN_x.gui.display))
                        answ=conn_questdlg({'Warning: There are pending jobs submitted but not yet finished','Changes to this project will be disregarded until all pending jobs are finished or canceled', 'Do you want to see these pending jobs now?'},'Warning!','Yes','No','Yes');
                        if isequal(answ,'Yes'), conn_jobmanager(info); end
                        return;
                    else
                        conn_disp('fprintf','Warning: pending jobs in %s not finished yet. Until then, any modifications to this project may be overwritten once the pending jobs finish and they are merged back into this project\n',localfilename);
                        return;
                    end
                elseif numel(vtag)>1
                    allfiles=cellstr(allfiles);
                    allfiles=char(allfiles(itag==ivtag));
                end
                alllogs=alllogs{ivtag};
            end
            if ~isempty(allfiles)
                conn_disp(allfiles);
                conn_disp('fprintf','Merging finished jobs. Please wait...');
                filename=CONN_x.filename;
                pobj=CONN_x.pobj;
                temp=load(deblank(allfiles(1,:)),'CONN_x','-mat');
                if ~isfield(temp,'CONN_x')||~isfield(temp.CONN_x,'pobj')||~isfield(temp.CONN_x.pobj,'holdsdata')||temp.CONN_x.pobj.holdsdata
                    conn_merge(allfiles);
                else
                    CONN_x=temp.CONN_x;
                    CONN_x.filename=filename;
                    CONN_x.pobj=pobj;
                    conn_merge(allfiles);
                    %if size(allfiles,1)>1, conn_merge(allfiles(2:end,:)); end
                end
                CONN_x.pobj.importedfiles=[CONN_x.pobj.importedfiles;reshape(cellstr(allfiles),[],1)];
                if isfield(temp,'CONN_x')&&isfield(temp.CONN_x,'pobj')&&isfield(temp.CONN_x.pobj,'holdsdata')&&~temp.CONN_x.pobj.holdsdata
                    id=regexp(cellstr(allfiles),'(\d+)\.dmat$','tokens','once');
                    id=unique([id{:}]);
                    addfiles={};
                    for n=1:numel(id)
                        for ianalysis=1:numel(CONN_x.Analyses)
                            if isfield(CONN_x.Analyses(ianalysis),'name')&&isfield(CONN_x.Analyses(ianalysis),'sourcenames')
                                filesourcenames=fullfile(CONN_x.folders.firstlevel,CONN_x.Analyses(ianalysis).name,'_list_sources.mat');
                                filesourcenames=conn_projectmanager('projectfile',filesourcenames,struct('id',id{n},'isextended',true),'.mat');
                                addfiles{end+1}=filesourcenames;
                            end
                        end
                        if isfield(CONN_x,'vvAnalyses')
                            for ianalysis=1:numel(CONN_x.vvAnalyses)
                                if isfield(CONN_x.vvAnalyses(ianalysis),'name')&&isfield(CONN_x.vvAnalyses(ianalysis),'measurenames')
                                    filemeasurenames=fullfile(CONN_x.folders.firstlevel_vv,CONN_x.vvAnalyses(ianalysis).name,'_list_measures.mat');
                                    filemeasurenames=conn_projectmanager('projectfile',filemeasurenames,struct('id',id{n},'isextended',true),'.mat');
                                    addfiles{end+1}=filemeasurenames;
                                end
                            end
                        end
                        if isfield(CONN_x.Setup.conditions,'allnames')
                            fileconditionnames=fullfile(CONN_x.folders.preprocessing,'_list_conditions.mat');
                            fileconditionnames=conn_projectmanager('projectfile',fileconditionnames,struct('id',id{n},'isextended',true),'.mat');
                            addfiles{end+1}=fileconditionnames;
                        end
                    end
                    if ~isempty(addfiles)
                        CONN_x.pobj.importedfiles=[CONN_x.pobj.importedfiles;reshape(addfiles,[],1)];
                    end
                end
                if ~isempty(alllogs)&&iscell(alllogs)
                    for n=1:numel(alllogs)
                        if ischar(alllogs{n})
                            flog=alllogs{n};
                            if ~conn_existfile(flog), flog=regexprep(flog,'\.stdlog$','.stdout'); end % fix for PC/Mac background Matlab-based jobs
                            if conn_existfile(flog)
                                conn_disp('fprintf','Importing log descriptions from %s :\n',flog);
                                cwcopy=conn_disp('__cwcopy');
                                tbspace=conn_disp('__tbspace');
                                conn_disp('__cwcopy',false);
                                conn_disp('__tbspace',tbspace+6);
                                try
                                    str=regexp(fileread(flog),'[\r\n]+','split');
                                    str=sprintf('%s\n',str{:});
                                    conn_disp(char(str));
                                end
                                conn_disp('__tbspace',tbspace);
                                conn_disp('__cwcopy',cwcopy);
                            end
                        end
                    end
                end
                conn_disp('fprintf','Done\n');
            end
            localfilename=conn_projectmanager('projectfile',CONN_x.filename,struct('id','*','isextended',true),'.emat');
            allfiles=conn_dir(localfilename,'-R');
            if ~isempty(allfiles)
                conn_disp('fprintf','Performing delayed processing steps. Please wait...');
                conn_disp(unique(cellstr(allfiles)));
                psteps={};
                for n=1:size(allfiles,1)
                    temp=load(deblank(allfiles(n,:)),'-mat');
                    if isfield(temp,'process')&&~isempty(temp.process)
                        if ischar(temp.process), temp.process={temp.process};
                        elseif ~iscell(temp.process),temp.process=num2cell(temp.process); 
                        end
                        for n1=1:numel(temp.process)
                            if ~iscell(temp.process{n1}), temp.process{n1}={temp.process{n1}}; end
                            if ~any(cellfun(@(x)isequal(temp.process{n1},x),psteps))
                                psteps=[psteps, temp.process(n1)];
                            end
                        end
                    end
                end
                for n=1:numel(psteps)
                    conn_process(psteps{n}{:}); 
                end
                CONN_x.pobj.importedfiles=[CONN_x.pobj.importedfiles;reshape(cellstr(allfiles),[],1)];
                conn_disp('fprintf','Done\n');
            end
        end
        
    case 'cleanproject'  % removes any delayed-writing files
        if CONN_x.pobj.holdsdata
            for n=1:numel(CONN_x.pobj.importedfiles)
                if ispc, [ok,nill]=system(['del "',deblank(CONN_x.pobj.importedfiles{n}),'"']);
                else     [ok,nill]=system(['rm -f ''',deblank(CONN_x.pobj.importedfiles{n}),'''']);
                end
                if conn_existfile(deblank(CONN_x.pobj.importedfiles{n})), conn_disp('fprintf','Unable to delete file %s. Check file/folder permissions and try again\n',CONN_x.pobj.importedfiles{n}); end
            end
            CONN_x.pobj.importedfiles={};
        end
        
    case 'addstep'  % adds delayed functional step to base project (conn_process step) 
        if CONN_x.pobj.holdsdata,
            localfilename=CONN_x.filename;
        else
            localfilename=conn_projectmanager('projectfile');
        end
        filename=conn_prepend('',localfilename,'.emat');
        process={};
        if conn_existfile(filename), load(filename,'process','-mat'); end
        process=[process {varargin}];
        save(filename,'process');
        
    otherwise,
        error('unrecognized option',option);
            
end
end

