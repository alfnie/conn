function [sources,idxsources,C]=conn_peaks(filename,DisMin,SizeMin)
% internal function

METHOD=1;
if nargin<2||isempty(DisMin), DisMin=0; end %12; end
if nargin<3||isempty(SizeMin), SizeMin=0; end %.05; end
if iscell(filename)
    b=filename{1};
    mat=filename{2};
else
    a=spm_vol(filename);
    b=spm_read_vols(a);
    mat=a.mat;
end
switch(METHOD)
    case 1
        C=[];
        c=b;
        for n1=1:3
            c=cat(1,max(c(1,:,:),c(2,:,:)),max(c(3:end,:,:),max(c(2:end-1,:,:),c(1:end-2,:,:))),max(c(end-1,:,:),c(end,:,:)));
            c=shiftdim(c,1);
        end
        c=(c==b).*b;
        idx0=find(c>0);
        [sc,idx1]=sort(-c(idx0));
        [x,y,z]=ind2sub(size(c),idx0(idx1));
        xyz=[x,y,z,ones(size(x))]*mat';
        sources=zeros(numel(idx1),3);
        idxsources=zeros(numel(idx1),1);
        nsources=0;
        for n1=1:numel(idx1)
            if DisMin>0, d=sqrt((sources(1:nsources,1)-xyz(n1,1)).^2+(sources(1:nsources,2)-xyz(n1,2)).^2+(sources(1:nsources,3)-xyz(n1,3)).^2); end
            if DisMin==0||all(d>DisMin)
                nsources=nsources+1;
                sources(nsources,:)=xyz(n1,1:3);
                idxsources(nsources)=idx0(idx1(n1));
            else
                c(idx0(n1))=0;
            end
        end
        sources=sources(1:nsources,:);
        idxsources=idxsources(1:nsources);
        C=c>0;
        
    case 2
        idxvalid=find(~isnan(b)&b~=0);
        [sb,idx]=sort(b(idxvalid),'descend');
        idxvalid=idxvalid(idx);
        [i1,i2,i3]=ind2sub(size(b),idxvalid);
        [j1,j2,j3]=ind2sub(size(b),(1:numel(b))');
        xyz=[j1,j2,j3];
        xyz2=[xyz,ones(size(j1))]*mat(1:3,:)';
        center=sub2ind(size(b),ceil(size(b,1)/2),ceil(size(b,2)/2),ceil(size(b,3)/2));
        d2=sum(conn_bsxfun(@minus,xyz2,xyz2(center,:)).^2,2);
        DI=find(d2<DisMin^2);
        DI1=xyz(DI,1)-xyz(center,1);
        DI2=xyz(DI,2)-xyz(center,2);
        DI3=xyz(DI,3)-xyz(center,3);
        DI=DI-center;
        NMin=SizeMin*numel(DI);
        
        C=zeros(size(b));
        L=zeros(numel(idxvalid),1);
        nL=0;
        for n1=1:numel(idxvalid)
            if ~C(idxvalid(n1))
                di1=DI1+i1(n1);
                di2=DI2+i2(n1);
                di3=DI3+i3(n1);
                di=DI+idxvalid(n1);
                iv=find(di1>0&di1<=size(b,1)&di2>0&di2<=size(b,2)&di3>0&di3<=size(b,3));
                iv(C(di(iv))>0)=[];
                if numel(iv)>=NMin,
                    nL=nL+1;
                    C(di(iv))=nL;
                    L(nL)=idxvalid(n1);
                end
            end
        end
        L=L(1:nL);
        idxsources=L;
        sources=xyz2(L,:);
end
end
