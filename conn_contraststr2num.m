function value = conn_contraststr2num(txt)

value=str2num(txt); 
if isempty(value), % kronecker rules
    value=str2num([regexprep(txt,{'e(\d+)','d(\d+)','a(\d+)','(.*?) x '},{'eye($1)','diff(eye($1))','ones(1,$1)/$1','kron([$1],['},'ignorecase'),repmat('])',1,numel(regexp(txt,' x ')))]); 
end 
if isempty(value), % eval in base
    try value=evalin('base',txt); 
    catch, value=[]; end; 
end