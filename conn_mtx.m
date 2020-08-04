function varargout=conn_mtx(option,varargin)
%  internal function
% conn_mtx: manages large matrices that do not fit in memory
%           supports eigenvector-computation functionality
%
% C = conn_mtx('init',[a, b1, b2, b3, b4, ...],filename);
%     creates virtual [a,b1*b2*...] matrix C
%     The matrix C is stored as a [b2,b3,...] array of [a,b1] matrices (blocks)
%     where each individual block is expected to fit in memory
%
% x = conn_mtx('getblock',C,idx); 
%      reads [a,b1] matrix x from C idx-th block
%      i.e. ~ x=C(idx)
%
% conn_mtx('addtoblock',C,idx,x); 
%      adds [a,b1] matrix x to C idx-th block
%      i.e. ~ C(idx)+=x
%
% blockcolumns = conn_mtx('getblockcolumns',C);
%      returns columns of C associated with each [a,b1] block
%
% fC=conn_mtx('multiplicationhandle',C);
%      returns function handle fC, such that Y = fC(X) returns Y=C*X for any input matrix X
%      additional options: Y = fC(X,'transp') returns Y=C'*X
%                          Y = fC(X,'notransp') returns Y=C*X
%
% [Q,D] = conn_mtx('svd',C,K,...);
%      returns [a,K] matrix Q containing the first K left singular vectors of virtual [a,b*n1*n2*...] matrix C
%      optional additional inputs for svds(FUN,N,K,...) function 
%
% [Q,D] = conn_mtx('eig',C,K,...);
%      returns [a,K] matrix Q containing the first K left eigenvectors of virtual [a,b*n1*n2*...] matrix C
%      optional additional inputs for eigs(FUN,N,K,...) function 
%


switch(lower(option))
    case 'init' 
        [dims,filename]=deal(varargin{:});
        dims=[dims ones(1,max(0,4-length(dims)))];
        nblocks=prod(dims(3:end));
        Msize=[dims(1) prod(dims(2:end))];
        filenames=arrayfun(@(n)[filename,sprintf('.%d.dat',n)],1:nblocks,'uni',0);
        C=struct('fname',filename,...
                'filenames',{filenames},...
                'size',Msize,...
                'blocksize_in',dims(1:2),...
                'blocksize_out',dims(3:end),...
                'numberofblocks',nblocks);
        M=zeros(dims(1:2));
        for nblock=1:nblocks
            conn_mtx_save(C,nblock,M);
        end
        varargout={C};
        
    case 'getblockcolumns'
        C=varargin{1};
        blockcolumns=mat2cell(1:C.blocksize_in(2)*C.numberofblocks,1,C.blocksize_in(2)+zeros(1,C.numberofblocks));
        varargout={blockcolumns};

    case 'addtoblock'
        [C,nblock,c]=deal(varargin{:});
        M=conn_mtx_load(C,nblock);
        M=M+c;
        conn_mtx_save(C,nblock,M);
        varargout={M};
        
    case 'getblock'
        [C,i]=deal(varargin{:});
        nblock=i(1);
        if numel(i)>1, nblock=nblock+(i(2:end)-1)*cumprod(C.blocksize_out(1:numel(i)-1))'; end
        M=conn_mtx_load(C,nblock);
        varargout={M};
        
    case 'multiplicationhandle'
        C=varargin{1};
        fh=@(x,varargin)conn_mtx_mult(x,C,varargin{:});
        varargout={fh};
        
    case 'eig'
        C=varargin{1};
        if numel(varargin)>=2, NdimsOut=varargin{2}; else NdimsOut=1; end
        if numel(varargin)>=3, opts=varargin(3:end); else opts={}; end
        fC=conn_mtx('multiplicationhandle',C);
        [Q0,D]=eigs(fC,C.size(1),NdimsOut,opts{:});
        varargout={Q0,D};
        
    case 'svd'
        C=varargin{1};
        if numel(varargin)>=2, NdimsOut=varargin{2}; else NdimsOut=1; end
        if numel(varargin)>=3, opts=varargin(3:end); else opts={}; end
        fC=conn_mtx('multiplicationhandle',C);
        [Q0,D]=svds(fC,C.size,NdimsOut,opts{:});
        varargout={Q0,D};
end
end

function y=conn_mtx_mult(x,C,option)
if nargin>2&&isequal(option,'transp')
    y=zeros(C.size(2),size(x,2));
    for nblock=1:C.numberofblocks,
        M=conn_mtx_load(C,nblock);
        idx=C.blocksize_in(2)*(nblock-1)+(1:C.blocksize_in(2));
        y(idx,:)=M'*x;
    end
else
    y=zeros(C.size(1),size(x,2));
    for nblock=1:C.numberofblocks,
        M=conn_mtx_load(C,nblock);
        idx=C.blocksize_in(2)*(nblock-1)+(1:C.blocksize_in(2));
        y=y+M*x(idx,:);
    end
end
end

function M=conn_mtx_load(C,nblock)
fname=C.filenames{nblock};
fh=fopen(fname,'rb');
M=reshape(fread(fh,inf,'float64'),C.blocksize_in);
fclose(fh);
end

function conn_mtx_save(C,nblock,M)
fname=C.filenames{nblock};
fh=fopen(fname,'wb');
fwrite(fh,M,'float64');
fclose(fh);
end



        