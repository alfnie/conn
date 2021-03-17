function varargout=conn_dirn(filename,option,varargin)
% SEE ALSO dir

if nargin<1||isempty(filename), filename='.'; end
if nargin<2||isempty(option), option=''; end

if any(conn_server('util_isremotefile',filename)), out=conn_server('run',mfilename,conn_server('util_localfile',filename),option); 
elseif isequal(option,'-ls'), out=ls('-al',filename);
else out=dir(filename);
end

if ~nargout,
    if isequal(option,'-ls'), disp(out);
    elseif isstruct(out), disp(char({out.name})); 
    end
    varargout={}; 
else varargout={out}; 
end
end

