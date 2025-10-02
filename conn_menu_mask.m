function h=conn_menu_mask(varargin)

global CONN_gui
if isfield(CONN_gui,'isjava')&&~CONN_gui.isjava
    i=find(strcmpi(varargin(1:2:end-1),'backgroundcolor'),1);
    if ~isempty(i), h=uipanel(varargin{:},'bordercolor',varargin{2*i});
    else h=uipanel(varargin{:});
    end
else
    h=uicontrol('style','frame',varargin{:});
end
end
