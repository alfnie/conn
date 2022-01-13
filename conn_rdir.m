function varargout=conn_rdir(filename,option,varargin)
% SEE ALSO dir

if nargin<1||isempty(filename), filename='.'; end
if nargin<2||isempty(option), option=''; end

out=conn_server('run','conn_dirn',conn_server('util_localfile',filename),option,varargin{:}); 

if ~nargout,
    if isequal(option,'-ls'), disp(out);
    elseif isstruct(out), disp(char({out.name})); 
    end
    varargout={}; 
else varargout={out}; 
end
end

