function varargout = conn_link(varargin)
% internal function
% SEE CONN_REMOTELY
[varargout{1:nargout}]=conn_remotely(varargin{:});
