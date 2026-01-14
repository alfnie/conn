function varargout=conn_menu_selectsubjects(hdl,varargin)
global CONN_x;
validcovariates=find(cellfun(@(x)isempty(regexp(x,'^_')),CONN_x.Setup.l2covariates.names(1:end-1)));
v=listdlg('liststring',CONN_x.Setup.l2covariates.names(validcovariates),'selectionmode','multiple','initialvalue',1,'promptstring',{'Select group-defining covariate(s)','(0/1 values defining subjects to include in these analyses)'},'ListSize',[400 250]);
if ~isempty(v),
    v=validcovariates(v);
    valid=[];
    for n=1:numel(v)
        values=conn_module('get','l2covariates',CONN_x.Setup.l2covariates.names{v(n)});
        valid=union(valid,find(~isnan(values)&values~=0));
    end
    set(hdl,'value',valid);
end
end

