function str=conn_cell2html(str)
global CONN_gui

if ischar(str), str=cellstr(str); end
if ~iscell(str), str={str}; end
str=cellfun(@char,str,'uni',0);
str=regexprep(str,'\\|\/','\\');
if numel(str)>1, 
    if isfield(CONN_gui,'isjava')&&CONN_gui.isjava, str=['<HTML>',strjoinstr(str,'<br/>'),'</HTML>']; end
end
str=char(str);
end

function str=strjoinstr(str1,str2)
str=[str1(:)';repmat({str2},1,length(str1))];
str=reshape(str(1:end-1),1,numel(str)-1);
str=[str{:}];
end

