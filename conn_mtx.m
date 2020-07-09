function varargout=conn_mtx(option,varargin)
% conn_mtx: manages large matrices that do not fit in memory
%           supports eigenvector-computation functionality
%  internal function

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
        C=deal(varargin{:});
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
        C=deal(varargin{:});
        fh=@(x)conn_mtx_mult(x,C);
        varargout={fh};
end
end

function y=conn_mtx_mult(x,C)
    y=zeros(C.size(1),size(x,2));
    for nblock=1:C.numberofblocks,
        M=conn_mtx_load(C,nblock);
        idx=C.blocksize_in(2)*(nblock-1)+(1:C.blocksize_in(2));
        y=y+M*x(idx,:);
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



        