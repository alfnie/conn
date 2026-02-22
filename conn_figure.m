function varargout=conn_figure(varargin)
global CONN_gui
if ~isfield(CONN_gui,'themeoptsPopup'), CONN_gui.themeoptsPopup={}; end
[varargout{1:nargout}]=figure(varargin{:},CONN_gui.themeoptsPopup{:});
end