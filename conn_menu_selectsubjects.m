function varargout=conn_menu_selectsubjects(hdl,varargin)
global CONN_x;
v=listdlg('liststring',CONN_x.Setup.l2covariates.names(1:end-1),'selectionmode','single','initialvalue',1,'promptstring',{'Select group-defining covariate','(0/1 values defining subjects to include in these analyses)'},'ListSize',[400 250]);
if ~isempty(v),
    values=conn_module('get','l2covariates',CONN_x.Setup.l2covariates.names{v});
    valid=find(~isnan(values)&values~=0);
    set(hdl,'value',valid);
end
end

