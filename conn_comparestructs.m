function [FIELDS,STATUS,DESCRIPTION]=conn_comparestructs(A,B,varargin)
% for internal use (identifying/merging structures with conflicts)
%
% CONN_COMPARESTRUCTS(A,B)
% compares recursively A to B and list the fields FIELDS and changes CHANGE
% that would need to be applied to make A equal to B
% see ALSO CONN_REPLACESTRUCTS
%
% e.g.
% [FIELDS,CHANGE,str]=conn_comparestructs(a,b); 
% disp(char(str)); 
% anew=conn_replacestructs(a,b,FIELDS,CHANGE); 
% disp(isequaln(anew,b))

OPTIONS=struct(...
    'seed','',...
    'CheckStructs',true,...         % goes into structure arrays
    'CheckCells',true,...           % goes into cell arrays
    'MaxLevel',inf);              % maximum number of levels

for nvarargin=1:2:numel(varargin)-1,
    assert(isfield(OPTIONS,varargin{nvarargin}),'unrecognized option %s',varargin{nvarargin})
    OPTIONS.(varargin{nvarargin})=varargin{nvarargin+1};
end

DODISP=~nargout;
if ischar(A), A=conn_loadmatfile(A); end
if ischar(B), B=conn_loadmatfile(B); end
conn_comparestructs_output('init');
conn_comparestructs_internal(A,B,OPTIONS.seed,0);

    function conn_comparestructs_internal(a,b,basestr,level)
        if level>OPTIONS.MaxLevel, return; end
        %if isempty(a), conn_comparestructs_output('disp',basestr,'change value', [basestr,' empty input A']); return; end
        %if isempty(b), conn_comparestructs_output('disp',basestr,'change value', [basestr,' empty input B']); return; end
        na=numel(a);
        nb=numel(b);
        %if na~=nb, conn_comparestructs_output('disp',basestr, 'change value', sprintf(' mismatch in %s number of elements (%d in A; %d in B)',[basestr],na,nb)); return; end
        if xor(iscell(a),iscell(b)),                                conn_comparestructs_output('disp',basestr, 'change value', sprintf('non-cell input %c',char('A'+iscell(a)))); return; end
        if xor(isstruct(a),isstruct(b)),                            conn_comparestructs_output('disp',basestr, 'change value', sprintf('non-structure input %c',char('A'+isstruct(a)))); return; end
        if ~iscell(a)&&~isstruct(a), % not a cell nor a struct -> compare directly
            if ~testequal(a,b),                                     conn_comparestructs_output('disp',basestr, 'change value', 'unequal elements'); end
        elseif iscell(a)&&~OPTIONS.CheckCells % two cell arrays -> compare directly
            if ~testequal(a,b),                                     conn_comparestructs_output('disp',basestr, 'change value', 'unequal cells'); end
        elseif isstruct(a)&&~OPTIONS.CheckStructs % two struct arrays -> compare directly
            if ~testequal(a,b),                                     conn_comparestructs_output('disp',basestr, 'change value', 'unequal structs'); end
        elseif isstruct(a),
            if ~testequal(a,b), % two struct arrays -> compare elements and each field
                f=union(fieldnames(a),fieldnames(b));
                for n=1:numel(f)
                    if ~isfield(a,f{n}),
                        if isempty(basestr),                conn_comparestructs_output('disp', sprintf('%s',f{n}), 'add field', 'field missing in A');
                        else,                               conn_comparestructs_output('disp', sprintf('%s.%s',basestr,f{n}), 'add field', 'field missing in A');
                        end
                    elseif ~isfield(b,f{n}),
                        if isempty(basestr),                conn_comparestructs_output('disp', sprintf('%s',f{n}), 'remove field', 'field missing in B');
                        else,                               conn_comparestructs_output('disp', sprintf('%s.%s',basestr,f{n}), 'remove field', 'field missing in B');
                        end
                    end
                end
                for nc=max(na,nb):-1:1
                    if nc>numel(a),                         conn_comparestructs_output('disp', sprintf('%s(%d)',basestr,nc), 'add element', 'element missing in A');
                    elseif nc>numel(b),                     conn_comparestructs_output('disp', sprintf('%s(%d)',basestr,nc), 'remove element', 'element missing in B');
                    else
                        for n=1:numel(f),
                            a1=a(nc);
                            b1=b(nc);
                            if ~isfield(a,f{n}),
                                if isempty(basestr)&&max(na,nb)==1, conn_comparestructs_output('disp', f{n}, 'change value', 'initialize this new-field value');
                                elseif max(na,nb)==1,               conn_comparestructs_output('disp', sprintf('%s.%s',basestr,f{n}), 'change value', 'initialize this new-field value');
                                else                                conn_comparestructs_output('disp', sprintf('%s(%d).%s',basestr,nc,f{n}), 'change value', 'initialize this new-field value');
                                end
                            elseif ~isfield(b,f{n}),
                            elseif ~testequal(a1.(f{n}),b1.(f{n})),
                                if OPTIONS.CheckStructs&&(isstruct(a1.(f{n}))||iscell(a1.(f{n}))),
                                    if isempty(basestr)&&max(na,nb)==1, conn_comparestructs_internal(a1.(f{n}),b1.(f{n}),sprintf('%s',f{n}),level+1);
                                    elseif max(na,nb)==1,           conn_comparestructs_internal(a1.(f{n}),b1.(f{n}),sprintf('%s.%s',basestr,f{n}),level+1);
                                    else                            conn_comparestructs_internal(a1.(f{n}),b1.(f{n}),sprintf('%s(%d).%s',basestr,nc,f{n}),level+1);
                                    end
                                else
                                    if isempty(basestr)&&max(na,nb)==1, conn_comparestructs_output('disp', f{n}, 'change value', 'unequal field values');
                                    elseif max(na,nb)==1,           conn_comparestructs_output('disp', sprintf('%s.%s',basestr,f{n}), 'change value', 'unequal field values');
                                    else                            conn_comparestructs_output('disp', sprintf('%s(%d).%s',basestr,nc,f{n}), 'change value', 'unequal field values');
                                    end
                                end
                            end
                        end
                    end
                end
            end
        elseif iscell(a) % two cell arrays -> look inside
            if ~testequal(a,b), % two struct arrays -> compare each field
                for nc=max(na,nb):-1:1
                    if nc>numel(a),                                 conn_comparestructs_output('disp', sprintf('%s{%d}',basestr,nc), 'add element', 'cell array element missing in A');
                    elseif nc>numel(b),                             conn_comparestructs_output('disp', sprintf('%s{%d}',basestr,nc), 'remove element', 'cell array element missing in B');
                    else
                        a1=a{nc};
                        b1=b{nc};
                        if ~testequal(a1,b1),
                            if OPTIONS.CheckCells&&(isstruct(a1)||iscell(a1)), conn_comparestructs_internal(a1,b1,sprintf('%s{%d}',basestr,nc),level+1);
                            else                                    conn_comparestructs_output('disp', sprintf('%s{%d}',basestr,nc), 'change value', 'unequal cell array element values');
                            end
                        end
                    end
                end
            end
        end
    end

    function conn_comparestructs_output(option,field,status,str)
        switch(option)
            case 'init'
                if DODISP, disp('Differences between A and B:');
                else 
                    FIELDS={};
                    DESCRIPTION={};
                    STATUS={};
                end
            case 'disp'
                str=sprintf('%s %s (%s)',status,field,str);
                if DODISP, disp(str);
                else
                    FIELDS{end+1}=field;
                    DESCRIPTION{end+1}=str;
                    STATUS{end+1}=status;
                end
        end
    end
end

function ok=testequal(a,b)
try, if isa(a,'function_handle'), a=func2str(a); end; end
try, if isa(b,'function_handle'), b=func2str(b); end; end
ok=isequalwithequalnans(a,b);
end