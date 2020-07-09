function conn_comparestructs(a,b,basestr)
CHECKSTRUCTS=true;
CHECKCELLS=true;
if nargin<3, 
    basestr=''; 
    disp('Differences between A and B:');
    if ischar(a), a=load(a); end
    if ischar(b), b=load(b); end
end
if isempty(a), disp([basestr,' empty input A']); return; end
if isempty(b), disp([basestr,' empty input B']); return; end
na=numel(a);
nb=numel(b);
if na~=nb, disp(sprintf(' mismatch in %s number of elements (%d in A; %d in B)',[basestr],na,nb)); end
if xor(iscell(a),iscell(b)), disp(sprintf(' mismatch in %s non-cell input %c',basestr,char('A'+iscell(a)))); return; end
if xor(isstruct(a),isstruct(b)), disp(sprintf(' mismatch in %s non-structure input %c',basestr,char('A'+isstruct(a)))); return; end
if ~iscell(a)&&~isstruct(a), 
    if ~isequalwithequalnans(a,b),disp(sprintf(' mismatch in %s',basestr)); end
    return
end
if iscell(a)||min(na,nb)>1
    for nc=1:min(na,nb)
        if iscell(a)&&~isequalwithequalnans(a{nc},b{nc}),
            if CHECKCELLS&&(isstruct(a{nc})||iscell(a{nc})), conn_comparestructs(a{nc},b{nc},sprintf('%s{%d}',basestr,nc));
            else disp(sprintf(' mismatch in %s{%d}',basestr,nc));
            end
        elseif isstruct(a)&&~isequalwithequalnans(a(nc),b(nc)),
            if CHECKSTRUCTS&&(isstruct(a(nc))||iscell(a(nc))), conn_comparestructs(a(nc),b(nc),sprintf('%s(%d)',basestr,nc));
            else disp(sprintf(' mismatch in %s(%d)',basestr,nc));
            end
        end
    end
    return
else
    a=a(1);
    b=b(1);
    f=union(fieldnames(a),fieldnames(b));
    for n=1:numel(f),
        ok=true;
        if ~isfield(a,f{n}), disp([' field ',basestr,f{n},' non-existing in A']); ok=false; end
        if ~isfield(b,f{n}), disp([' field ',basestr,f{n},' non-existing in B']); ok=false; end
        if ok&&~isequalwithequalnans(a.(f{n}),b.(f{n})),
            if isstruct(a.(f{n}))||iscell(a.(f{n})),
                if isempty(basestr), conn_comparestructs(a.(f{n}),b.(f{n}),f{n});
                else conn_comparestructs(a.(f{n}),b.(f{n}),[basestr,'.',f{n}]);
                end
            else
                disp([' mismatch in ',basestr,'.',f{n}]);
            end
        end
    end
end
