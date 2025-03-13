function [A,OK]=conn_replacestructs(A,B,FIELDS,CHANGE)
% for internal use (merging structures with conflicts)
%
% Anew=CONN_REPLACESTRUCTS(A,B,FIELDS,CHANGE)
% modifies A by changing the specific fields FIELDS in the way indicated in CHANGE
% (and using the information from B at the same fields when necessary)
% note: order of steps indicated in CHANGE matters (some of the steps are non-commutative)
% see ALSO CONN_COMPARESTRUCTS
%
% e.g.
% [FIELDS,CHANGE,str]=conn_comparestructs(a,b); 
% disp(char(str)); 
% anew=conn_replacestructs(a,b,FIELDS,CHANGE); 
% disp(isequaln(anew,b))

OK=true;
if ischar(A), A=conn_loadmatfile(A); end
if ischar(B), B=conn_loadmatfile(B); end
for n=1:numel(FIELDS)
    fields=regexp(FIELDS{n},'[\.\(\{]+','split');
    try
        switch(CHANGE{n})
            case {'change value','add element'}
                A=conn_replacestructs_change(A,B,fields);
            case 'add field'
                A=conn_replacestructs_addfield(A,fields);
            case 'remove field'
                A=conn_replacestructs_rmfield(A,fields);
            case 'remove element'
                A=conn_replacestructs_delete(A,fields);
        end
    catch
       fprintf('Warning: Project update conflict resolution is unable to incorporate the change %s (skipping)\n',FIELDS{n});
       OK=false;
    end
end
end

function a=conn_replacestructs_change(a,b,fields)
if isempty(fields)
    a=b;
elseif numel(fields)==1
    if fields{1}(end)==')'||fields{1}(end)=='}',
        idx=str2double(fields{1}(1:end-1));
        a(idx)=b(idx);
    else
        a.(fields{1})=b.(fields{1});
    end
else
    if fields{1}(end)==')',
        idx=str2double(fields{1}(1:end-1));
        a(idx)=conn_replacestructs_change(a(idx), b(idx), fields(2:end));
    elseif fields{1}(end)=='}',
        idx=str2double(fields{1}(1:end-1));
        a{idx}=conn_replacestructs_change(a{idx}, b{idx}, fields(2:end));
    else
        a.(fields{1})=conn_replacestructs_change(a.(fields{1}), b.(fields{1}), fields(2:end));
    end
end
end


function a=conn_replacestructs_addfield(a,fields)
if isempty(fields)
    a=[];
elseif numel(fields)==1
    [a.(fields{1})]=deal([]);
else
    if fields{1}(end)==')',
        idx=str2double(fields{1}(1:end-1));
        a(idx)=conn_replacestructs_addfield(a(idx), fields(2:end));
    elseif fields{1}(end)=='}',
        idx=str2double(fields{1}(1:end-1));
        a{idx}=conn_replacestructs_addfield(a{idx}, fields(2:end));
    else
        a.(fields{1})=conn_replacestructs_addfield(a.(fields{1}), fields(2:end));
    end
end
end

function a=conn_replacestructs_rmfield(a,fields)
if isempty(fields)
elseif numel(fields)==1
    a=rmfield(a,fields{1});
else
    if fields{1}(end)==')',
        idx=str2double(fields{1}(1:end-1));
        a(idx)=conn_replacestructs_rmfield(a(idx), fields(2:end));
    elseif fields{1}(end)=='}',
        idx=str2double(fields{1}(1:end-1));
        a{idx}=conn_replacestructs_rmfield(a{idx}, fields(2:end));
    else
        a.(fields{1})=conn_replacestructs_rmfield(a.(fields{1}), fields(2:end));
    end
end
end

function a=conn_replacestructs_delete(a,fields)
if isempty(fields)
elseif numel(fields)==1
    if fields{1}(end)==')'||fields{1}(end)=='}',
        idx=str2double(fields{1}(1:end-1));
        a(idx)=[];
    else
        disp('debug')
        disp(a)
        disp(fields{1})
    end
else
    if fields{1}(end)==')',
        idx=str2double(fields{1}(1:end-1));
        a(idx)=conn_replacestructs_delete(a(idx), fields(2:end));
    elseif fields{1}(end)=='}',
        idx=str2double(fields{1}(1:end-1));
        a{idx}=conn_replacestructs_delete(a{idx}, fields(2:end));
    else
        a.(fields{1})=conn_replacestructs_delete(a.(fields{1}), fields(2:end));
    end
end
end
