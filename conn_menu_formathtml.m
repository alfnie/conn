function str=conn_menu_formathtml(str)
global CONN_gui
if isfield(CONN_gui,'isjava')&&~CONN_gui.isjava, str=regexprep(str,{'<\/?HTML>|<\/?b>|<\/?i>|<\/?p>|<\/?em>|<\/?small>','<br/>'},{'','\n'},'ignorecase'); end
end
