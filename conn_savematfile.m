function conn_savematfile(filename,varargin)
% saves .mat file
% conn_savematfile(FILENAME,VARIABLE_NAMES,...,FORMAT)
% conn_savematfile(FILENAME,'-struct',STRUCT,FORMAT)
% SEE ALSO save


if isequal(varargin{1},'-struct')
    data=varargin{2};
    if ischar(data), data=evalin('caller',data); end
    opts=[false false true(1,numel(varargin)-2)];
else
    data=struct;
    opts=cellfun('length',regexp(varargin,'^-'))>0;
    for idx=reshape(find(~opts),1,[])
        data.(varargin{idx})=evalin('caller',varargin{idx});
    end
end
if any(conn_server('util_isremotefile',filename)), 
    if ~isempty(regexp(filename,'\<SPM\.mat$'))&&isfield(data,'SPM')
        try, if isfield(data.SPM,'swd'), data.SPM.swd=conn_server('util_localfile',data.SPM.swd); end; end
        try, if isfield(data.SPM.xY,'VY'), for n=1:numel(data.SPM.xY.VY), data.SPM.xY.VY(n).fname=conn_server('util_localfile',data.SPM.xY.VY(n).fname); end; end; end
        try, if isfield(data.SPM.xY,'P'), data.SPM.xY.P=conn_server('util_localfile',data.SPM.xY.P); end; end
        try, if isfield(data.SPM,'xCon'), for n=1:numel(data.SPM.xCon), data.SPM.xCon(n).Vspm.fname=conn_server('util_localfile',data.SPM.xCon(n).Vspm.fname); data.SPM.xCon(n).Vcon.fname=conn_server('util_localfile',data.SPM.xCon(n).Vcon.fname); end; end; end
        try, if isfield(data.SPM,'xX_multivariate'), data.SPM.xX_multivariate.Zfiles=conn_server('util_localfile',data.SPM.xX_multivariate.Zfiles); end; end
    end
    if ~isempty(regexp(filename,'\<info\.mat$'))&&isfield(data,'info')
        for fnames={'stdout','stderr','stdlog','scripts','pathname'}
            try, if isfield(data.info,fnames{1}), data.info.(fnames{1})=conn_server('util_localfile',data.info.(fnames{1})); end; end
        end
        for n=1:numel(data.info.private),
            for fnames={'project','pathname'}
                try, data.info.private{n}.(fnames{1})=conn_server('util_localfile',data.info.private{n}.(fnames{1})); end
            end
        end
    end
    if ~isempty(regexp(filename,'\<node\.\d+\.mat$'))&&isfield(data,'job')
        for n=1:numel(data.job),
            for fnames={'project','tag','pathname'}
                try, data.job(n).(fnames{1})=conn_server('util_localfile',data.job(n).(fnames{1})); end
            end
        end
    end
    conn_server('run',mfilename,conn_server('util_localfile',filename),'-struct',data,varargin{opts}); 
else        
    save(filename,'-struct','data',varargin{opts});
end
end
