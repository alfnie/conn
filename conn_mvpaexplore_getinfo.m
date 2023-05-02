function varargout=conn_mvpaexplore_getinfo(option,varargin);

switch(option)
    case 0 % [filename_B1,V1,Nt]=conn_mvpaexplore_getinfo(0,filepath,icondition,validsubjects);
        filepath=varargin{1};
        icondition=varargin{2};
        validsubjects=varargin{3};
        if any(conn_server('util_isremotefile',filepath))
            filepath=conn_server('util_localfile',filepath);
            [varargout{1:nargout}]=conn_server('run_keepas','filename_B1_V1',mfilename,option,filepath,icondition,validsubjects);
        else
            tstr=num2str(icondition,'%03d');
            filename_B1=arrayfun(@(nsub)fullfile(filepath,['vvPC_Subject',num2str(nsub,'%03d'),'_Condition',tstr,'.mat']), validsubjects,'uni',0);
            if nargout>2
                V1=conn_vol(filename_B1{1});
                Nt=cellfun(@(x)conn_fileutils('dir',[x,'c']).bytes,filename_B1);
                Nt=round(Nt/Nt(1)*V1.size.Nt);
                varargout={filename_B1,V1,Nt};
            elseif nargout>1
                V1=conn_vol(filename_B1{1});
                varargout={filename_B1,V1};
            else varargout={filename_B1};
            end
        end

    case 1, % [X,Y,IDX]=conn_mvpaexplore_getinfo(1,filename_B1,V1,Nt,nslice,nv);
        filename_B1=varargin{1};
        V1=varargin{2};
        Nt=varargin{3};
        nslice=varargin{4};
        nv=varargin{5};
        if ~isempty(filename_B1)&&any(conn_server('util_isremotevar',filename_B1))
            [varargout{1:nargout}]=conn_server('run_keepas','X_Y_IDX',mfilename,option,filename_B1,V1,Nt,nslice,nv); 
            if nargout>2, varargout{3}=conn_server('run',varargout{3}); end % note: passes IDX directly instead of as a link
        else
            X={};
            Y={};
            IDX=[];
            for n1=1:numel(filename_B1),
                y1=V1;
                y1.fname=filename_B1{n1};
                y1.size.Nt=Nt(n1);
                [Y{n1},IDX]=conn_get_slice(y1,nslice);
                X{n1}=conn_get_voxel(y1,nv);
            end
            varargout={X,Y,IDX};
        end
    case 2, % [X]=conn_mvpaexplore_getinfo(2,filename_B1,V1,Nt,nv);
        filename_B1=varargin{1};
        V1=varargin{2};
        Nt=varargin{3};
        nv=varargin{4};
        if ~isempty(filename_B1)&&any(conn_server('util_isremotevar',filename_B1))
            [varargout{1:nargout}]=conn_server('run_keepas','X',mfilename,option,filename_B1,Nt,V1,nv); 
        else
            X={};
            for n1=1:numel(filename_B1),
                y1=V1;
                y1.fname=filename_B1{n1};
                y1.size.Nt=Nt(n1);
                X{n1}=conn_get_voxel(y1,nv);
            end
            varargout={X};
        end
    case 3, %[Y,IDX]=conn_mvpaexplore_getinfo(3,filename_B1,V1,Nt,nslice);
        filename_B1=varargin{1};
        V1=varargin{2};
        Nt=varargin{3};
        nslice=varargin{4};
        if ~isempty(filename_B1)&&any(conn_server('util_isremotevar',filename_B1))
            [varargout{1:nargout}]=conn_server('run_keepas','Y_IDX',mfilename,option,filename_B1,V1,Nt,nslice); 
            if nargout>1, varargout{2}=conn_server('run',varargout{2}); end % note: passes IDX directly instead of as a link
        else
            Y={};
            IDX=[];
            for n1=1:numel(filename_B1),
                y1=V1;
                y1.fname=filename_B1{n1};
                y1.size.Nt=Nt(n1);
                [Y{n1},IDX]=conn_get_slice(y1,nslice);
            end
            varargout={Y,IDX};
        end

    case 4, % xy=conn_mvpaexplore_getinfo(4,X,Y,W);
        X=varargin{1};
        Y=varargin{2};
        W=varargin{3};
        if ~isempty(X)&&any(conn_server('util_isremotevar',X))
            [varargout{1:nargout}]=conn_server('run',mfilename,option,X,Y,W); 
        else
            xy=0;
            for n1=1:numel(X)
                if any(W(:,n1)~=0)&&~isempty(X{n1}), xy=xy+W(:,n1)*X{n1}'*Y{n1}; end
            end
            varargout={xy};
        end

end
end
