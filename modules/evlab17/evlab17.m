function varargout=evlab17(varargin)
% evlab17 internal function

if isempty(which(sprintf('evlab17_%s',varargin{1}))), varargin=[{'module'},varargin]; end
if ischar(varargin{1})&&~isempty(which(sprintf('evlab17_%s',varargin{1}))),
    fh=eval(sprintf('@evlab17_%s',varargin{1}));
    if ~nargout, feval(fh,varargin{2:end});
    else [varargout{1:nargout}]=feval(fh,varargin{2:end});
    end
else
    disp(sprintf('unrecognized evlab17 option %s or evlab17_%s function',varargin{1},varargin{1}));
end
