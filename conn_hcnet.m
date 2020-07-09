function [order,Labels]=conn_hcnet(C)
% CONN_HCNET hierarchical clustering based on connectivity networks
% p=CONN_HCNET(C)
%  returns node-permutation p such that for all values thr the connected 
%  subnetworks in the graph C(p,p)<thr are formed by contiguous nodes

[N1,N2]=size(C);
N=max(N1,N2);
[sC,idx]=sort(C(:),'ascend');
[i,j]=ind2sub(size(C),idx);
Labels=zeros(N,N);
order=1:N;
Labels(:,1)=(1:N)';
Index=ones(N,1);
Start=(1:N)';
End=(1:N)';
n0=1;
for n=find(~isnan(sC))'
    l1=Labels(i(n),n0);
    l2=Labels(j(n),n0);
    if l1~=l2
        n0=n0+1;
        c1=find(Labels(order,n0-1)==l1);
        c2=find(Labels(order,n0-1)==l2);
        i1=find(order(c1)==i(n));
        i2=find(order(c2)==j(n));
        l=l1;
        if c1(1)>c2(1)
            l=l2;
            ct=c1;c1=c2;c2=ct;
            it=i1;i1=i2;i2=it;
        end
        if i1<=numel(c1)/2, c1=flipud(c1); end
        if i2>numel(c2)/2, c2=flipud(c2); end
        minmax=[min(c1),max(c1),min(c2),max(c2)];
        Labels(:,n0)=Labels(:,n0-1);
        Labels(order(c1),n0)=l;
        Labels(order(c2),n0)=l;
        order=[order(1:minmax(1)-1),order(c1),order(c2),order(minmax(2)+1:minmax(3)-1),order(minmax(4)+1:end)];
    end
end
%disp('Label(i,thr) matrix: labels for connected subnetworks in graph C(p,p)<thr')
%disp(Labels(order,:));
end
