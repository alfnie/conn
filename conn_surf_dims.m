function dim=conn_surf_dims(resolution)
persistent alldims

if nargin<1||isempty(resolution), resolution=8; end
if numel(resolution)==1
    if numel(alldims)>=resolution && ~isempty(alldims{resolution}), 
        dim=alldims{resolution};
        return
    end
    nvertices=2+10*2^(2*resolution-2);
    n=[1,1,factor(nvertices)];
    minx=inf;
    for n1=1:numel(n)-2,
        for n2=n1+1:numel(n)-1,
            for n3=n2+1:numel(n)
                x=[prod(n(1:n1)),prod(n(n2:n3-1)),prod(n(n3:end))];
                tminx=std(x);
                if tminx<minx, dim=x; tminx=minx; end
            end
        end
    end
    if resolution==8, dim=dim([1 3 2]); end % note: fix so that resolution=5 fits within one slice of resolution=8
    alldims{resolution}=dim;
else
    dim=resolution;
    nvertices=prod(dim);
    resolution=round((log((nvertices-2)/10)/log(2)+2)/2);
    nvertices2=2+10*2^(2*resolution-2);
    if nvertices~=nvertices2, if ~nargout, conn_disp('not recognized data dimensions'); end; dim=[]; return; end
    for n1=1:12
        dim2=conn_surf_dims(n1);
        if isequal(dim,dim2),
            dim=n1;
            return;
        end
    end
    dim=[];
end



