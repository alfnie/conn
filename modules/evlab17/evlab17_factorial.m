function STR=evlab17_factorial(cells)
% internal use
%

% e.g. evlab17_factorial({'A_1','A_2','A_3','B_1','B_2','B_3'})
% returns lists of all main effects and interactions
%

names=regexp(cells,'_','split');
names=cat(1,names{:});
nfactors=size(names,2);
nconditions=size(names,1);
for n=1:nfactors
    [u{n},nill,i(:,n)]=unique(names(:,n));
    j{n}=0:numel(u{n});
end
[out{1:nfactors}]=ndgrid(j{:});
out=permute(cat(nfactors+1,out{:}),[nfactors+1,1:nfactors]);
STR={};N=[];
for nc=1:numel(out)/nfactors,
    f=out(:,nc);
    nf=nnz(f);
    if nf>0&nf<nfactors
        i1=find(f);
        str=''; 
        idx=1:nconditions;
        for n1=1:numel(i1), 
            nfactor=i1(n1);
            str=[str,u{nfactor}{f(nfactor)}]; 
            if n1<numel(i1), str=[str,'_']; end
            idx(i(idx,nfactor)~=f(nfactor))=[];
        end
        for n1=1:numel(idx)
            str=[str,' ',cells{idx(n1)},' ',num2str(1/numel(idx))];
        end
        STR{end+1}=str;
        N(end+1)=numel(idx);
    end
end
[nill,idx]=sort(-N);
STR=STR(idx);
if ~nargout, disp(char(STR)); end


