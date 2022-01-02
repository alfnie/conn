function data=conn_loadmatfile(filename,varargin)
% loads .mat file into structure
% data = conn_loadmatfile(filename)

DOKEEP=false;
DOCACHE=false;
idx=find(cellfun(@(x)isequal(x,'-run_keep'),varargin));
if ~isempty(idx),
    varargin(idx)=[];
    DOKEEP=true;
end
idx=find(cellfun(@(x)isequal(x,'-run_keepas'),varargin));
if ~isempty(idx),
    DOKEEP=varargin{idx(1)+1};
    varargin([idx idx+1])=[];
end
idx=find(cellfun(@(x)isequal(x,'-cache'),varargin));
if ~isempty(idx),
    varargin(idx)=[];
    DOCACHE=true;
end
if any(conn_server('util_isremotefile',filename)), 
    if ~DOCACHE&&isempty(varargin)&&isequal(DOKEEP,false)
        try
            info=conn_fileutils('dir',filename);
            DOCACHE=info.bytes>1e6;
        end
    end
    if DOCACHE
        filelocal=conn_cache('pull',filename);
        data=load(filelocal,varargin{:});
    elseif ischar(DOKEEP)&&~isempty(DOKEEP)
        data=conn_server('run_keepas',DOKEEP, mfilename,conn_server('util_localfile',filename),varargin{:});
    elseif DOKEEP
        data=conn_server('run_keep',mfilename,conn_server('util_localfile',filename),varargin{:});
    else
        data=conn_server('run',mfilename,conn_server('util_localfile',filename),varargin{:});
    end
    if ~isempty(regexp(filename,'\<SPM\.mat$'))&&isfield(data,'SPM')
        try, if isfield(data.SPM,'swd'), data.SPM.swd=conn_server('util_remotefile',data.SPM.swd); end; end
        try, if isfield(data.SPM.xY,'VY'), for n=1:numel(data.SPM.xY.VY), data.SPM.xY.VY(n).fname=conn_loadmatfile_utilremote(data.SPM.xY.VY(n).fname,data.SPM.swd); end; end; end
        try, if isfield(data.SPM.xY,'P'), data.SPM.xY.P=conn_loadmatfile_utilremote(data.SPM.xY.P,data.SPM.swd); end; end
        try, if isfield(data.SPM,'xCon'), for n=1:numel(data.SPM.xCon), data.SPM.xCon(n).Vspm.fname=conn_loadmatfile_utilremote(data.SPM.xCon(n).Vspm.fname,data.SPM.swd); data.SPM.xCon(n).Vcon.fname=conn_loadmatfile_utilremote(data.SPM.xCon(n).Vcon.fname,data.SPM.swd); end; end; end
        try, if isfield(data.SPM,'xX_multivariate'), data.SPM.xX_multivariate.Zfiles=conn_server('util_remotefile',data.SPM.xX_multivariate.Zfiles); end; end
    end
    if ~isempty(regexp(filename,'\<info\.mat$'))&&isfield(data,'info')
        for fnames={'stdout','stderr','stdlog','scripts','pathname'}
            try, if isfield(data.info,fnames{1}), data.info.(fnames{1})=conn_server('util_remotefile',data.info.(fnames{1})); end; end
        end
        for n=1:numel(data.info.private),
            for fnames={'project','pathname'}
                try, data.info.private{n}.(fnames{1})=conn_server('util_remotefile',data.info.private{n}.(fnames{1})); end
            end
        end
    end
    if ~isempty(regexp(filename,'\<node\.\d+\.mat$'))&&isfield(data,'job')
        for n=1:numel(data.job),
            for fnames={'project','tag','pathname'}
                try, data.job(n).(fnames{1})=conn_server('util_remotefile',data.job(n).(fnames{1})); end
            end
        end
    end
else
    data=load(conn_server('util_localfile',filename),varargin{cellfun(@(x)~isequal(x,'-cache'),varargin)});
end

if ~nargout, for fname=reshape(fieldnames(data),1,[]), assignin('caller',fname{1},data.(fname{1})); end; end

end

function filename=conn_loadmatfile_utilremote(filename,filepath)
ischarfilename=ischar(filename);
filename=cellstr(filename);
change=cellfun('length',regexp(filename,'[\\\/]'))==0;
for n=reshape(find(change),1,[]), filename{n}=fullfile(filepath,filename{n});  end
filename(~change)=conn_server('util_remotefile',filename(~change));
if ischarfilename, filename=char(filename); end
end
