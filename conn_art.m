function varargout = conn_art(varargin)
persistent check;
if isempty(check)||~check,
    str1=fileparts(which('art'));str2=fileparts(which('conn'));
    check=isequal(str1,str2);
end
if check,
    [varargout{1:nargout}]=art(varargin{:});
else
    error(['Another ART function exists besides conn/art. Please move the conn toolbox path up in Matlab''s ''File->Set Path'' to solve possible compatibility issues']);
    check=[];
end