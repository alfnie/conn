function [icond,isnew]=conn_conditionnames(name,flag,varargin)
global CONN_x;
if ~isfield(CONN_x.Setup.conditions,'allnames'), CONN_x.Setup.conditions.allnames=CONN_x.Setup.conditions.names(1:end-1); end
if nargin<2, flag='-'; end
idx=strmatch(name,CONN_x.Setup.conditions.allnames,'exact');
if ~isempty(idx),
    icond=idx(1);
    isnew=0;
    if strcmp(flag,'delete')
        CONN_x.Setup.conditions.allnames{icond}=[CONN_x.Setup.conditions.allnames{icond}, '_deleted', char(floor(10*rand(1,16))+'0')];
        fileconditionnames=fullfile(CONN_x.folders.preprocessing,'_list_conditions.mat');
        fileconditionnames=conn_projectmanager('projectfile',fileconditionnames,CONN_x.pobj,'.mat');
        allnames=CONN_x.Setup.conditions.allnames;
        save(fileconditionnames,'allnames');
    elseif strcmp(flag,'renamebutnotsave')
        CONN_x.Setup.conditions.allnames{icond}=varargin{1};
    elseif strcmp(flag,'rename')
        CONN_x.Setup.conditions.allnames{icond}=varargin{1};
        fileconditionnames=fullfile(CONN_x.folders.preprocessing,'_list_conditions.mat');
        fileconditionnames=conn_projectmanager('projectfile',fileconditionnames,CONN_x.pobj,'.mat');
        allnames=CONN_x.Setup.conditions.allnames;
        save(fileconditionnames,'allnames');
    end
else,
    icond=length(CONN_x.Setup.conditions.allnames)+1;
    if strcmp(flag,'recover')
        match1=strncmp(CONN_x.Setup.conditions.allnames,name,numel(name));
        match2=cellfun('length',regexp(CONN_x.Setup.conditions.allnames,'_deleted\d{16}$'));
        if any(match1&match2)
            icond=find(match1&match2,1);
            CONN_x.Setup.conditions.allnames{icond}=name;
            fileconditionnames=fullfile(CONN_x.folders.preprocessing,'_list_conditions.mat');
            fileconditionnames=conn_projectmanager('projectfile',fileconditionnames,CONN_x.pobj,'.mat');
            allnames=CONN_x.Setup.conditions.allnames;
            save(fileconditionnames,'allnames');
        end
    elseif strcmp(flag,'+'),
        CONN_x.Setup.conditions.allnames{icond}=name;
        fileconditionnames=fullfile(CONN_x.folders.preprocessing,'_list_conditions.mat');
        fileconditionnames=conn_projectmanager('projectfile',fileconditionnames,CONN_x.pobj,'.mat');
        allnames=CONN_x.Setup.conditions.allnames;
        save(fileconditionnames,'allnames');
    end
    isnew=1;
end
end
