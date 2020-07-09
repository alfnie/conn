function y=conn_bsxfun(fun,a,b)
persistent check;
if isempty(check),
    check=exist('bsxfun','builtin');
end
if check
    y=bsxfun(fun,a,b);
else
    % ensures compatibility with older matlab versions
    sa=size(a);
    sb=size(b);
    sa(numel(sa)+1:numel(sb))=1;
    sb(numel(sb)+1:numel(sa))=1;
    nd=max(sa,sb);
    k=nd;
    k(nd==sa)=1;
    if any(k>1), a=repmat(a,k); end
    k=nd;
    k(nd==sb)=1;
    if any(k>1), b=repmat(b,k); end
    y=feval(fun,a,b);
end

    