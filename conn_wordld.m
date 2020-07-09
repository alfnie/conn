function d = conn_wordld(s1,s2)
% Levenshtein distance between two words
if ischar(s1)&&ischar(s2)
    w=0:numel(s1);
    for j=1:numel(s2)
        prev = w;
        w(1)=+j;
        for i=1:numel(s1)
            w(i+1)=min(min(w(i),prev(i+1))+1,prev(i)+(s1(i)~=s2(j)));
        end
    end
    d=w(end);
else
    s1=cellstr(s1);
    s2=cellstr(s2);
    d=cellfun(@conn_wordld,repmat(s1(:),1,numel(s2)),repmat(s2(:)',numel(s1),1));
end
end

