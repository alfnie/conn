function varargout = fl_basic(varargin)
% internal function : wrapper for backwards-compatibility
if nargout>0, [varargout{1:nargout}]=fl(varargin{:});
else fl(varargin{:});
end
end