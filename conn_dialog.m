function varargout=conn_dialog(varargin)
global CONN_gui
if ~isfield(CONN_gui,'themeoptsPopup'), CONN_gui.themeoptsPopup={}; end
[varargout{1:nargout}]=dialog(varargin{:},CONN_gui.themeoptsPopup{:});
end