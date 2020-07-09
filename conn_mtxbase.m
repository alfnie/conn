function [ReducedC,MixingMtx]=conn_mtxbase(C)
MixingMtx=eye(size(C,1));
rankC=rank(C);
while rankC<size(MixingMtx,1) % reduce redundant dimensions
    for n2=size(MixingMtx,1):-1:1,
        rankcconditions2=rank(MixingMtx([1:n2-1,n2+1:size(MixingMtx,1)],:)*C);
        if rankcconditions2==rankC
            MixingMtx(n2,:)=[];
            break;
        end
    end
end
ReducedC=MixingMtx*C;
