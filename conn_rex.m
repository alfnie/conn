function varargout = conn_rex(varargin)
persistent check;
if isempty(check)||~check,
    str1=fileparts(which('rex'));str2=fileparts(which('conn'));
    check=isequal(str1,str2);
end
if check,
    [varargout{1:nargout}]=rex(varargin{:});
else
    error(['Another REX function exists besides conn/rex. Please move the conn toolbox path up in Matlab''s ''File->Set Path'' to solve possible compatibility issues']);
    check=[];
end