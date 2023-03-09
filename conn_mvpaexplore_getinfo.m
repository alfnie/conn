function varargout=conn_mvpaexplore_getinfo(option,varargin);

switch(option)
    case 1, % [X,Y,IDX,V1]=conn_mvpaexplore_getinfo(1,filename_B1,nslice,nv);
        filename_B1=varargin{1};
        nslice=varargin{2};
        nv=varargin{3};
        if ~isempty(filename_B1)&&any(conn_server('util_isremotefile',filename_B1{1}))
            filename_B1=conn_server('util_localfile',filename_B1);
            [varargout{1:nargout}]=conn_server('run',mfilename,option,filename_B1,nslice,nv); 
        else
            V1=conn_vol(filename_B1{1});
            X={};
            Y={};
            IDX=[];
            for n1=1:numel(filename_B1),
                y1=V1;
                y1.fname=filename_B1{n1};
                [Y{n1},IDX]=conn_get_slice(y1,nslice);
                X{n1}=conn_get_voxel(y1,nv(1:3));
            end
            varargout={X,Y,IDX,V1};
        end
    case 2, % [X]=conn_mvpaexplore_getinfo(2,filename_B1,V1,nv);
        filename_B1=varargin{1};
        V1=varargin{2};
        nv=varargin{3};
        if ~isempty(filename_B1)&&any(conn_server('util_isremotefile',filename_B1{1}))
            filename_B1=conn_server('util_localfile',filename_B1);
            [varargout{1:nargout}]=conn_server('run',mfilename,option,filename_B1,V1,nv); 
        else
            X={};
            for n1=1:numel(filename_B1),
                y1=V1;
                y1.fname=filename_B1{n1};
                X{n1}=conn_get_voxel(y1,nv(1:3));
            end
            varargout={X};
        end
    case 3, %[Y,IDX]=conn_mvpaexplore_getinfo(3,filename_B1,V1,nslice);
        filename_B1=varargin{1};
        V1=varargin{2};
        nslice=varargin{3};
        if ~isempty(filename_B1)&&any(conn_server('util_isremotefile',filename_B1{1}))
            filename_B1=conn_server('util_localfile',filename_B1);
            [varargout{1:nargout}]=conn_server('run',mfilename,option,filename_B1,V1,nslice); 
        else
            Y={};
            IDX=[];
            for n1=1:numel(filename_B1),
                y1=V1;
                y1.fname=filename_B1{n1};
                [Y{n1},IDX]=conn_get_slice(y1,nslice);
            end
            varargout={Y,IDX};
        end
end
end
