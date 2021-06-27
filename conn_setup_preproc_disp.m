function conn_setup_preproc_disp(matlabbatch,str)
if nargin<2, str=''; end
try
    if isempty(matlabbatch),
    elseif iscell(matlabbatch)&&ischar(matlabbatch{1})
        if numel(matlabbatch)>4
            conn_disp('fprintf','%s(%d) = %s\n',str,1,matlabbatch{1});
            conn_disp('fprintf','%s(%d) = %s\n',str,2,matlabbatch{2});
            conn_disp('fprintf','%s(%d) = %s\n',str,numel(matlabbatch)-1,matlabbatch{end-1});
            conn_disp('fprintf','%s(%d) = %s\n',str,numel(matlabbatch),matlabbatch{end});
        elseif numel(matlabbatch)>1
            for n=1:numel(matlabbatch)
                conn_disp('fprintf','%s(%d) = %s\n',str,n,matlabbatch{n});
            end
        else
            conn_disp('fprintf','%s = %s\n',str,matlabbatch{1});
        end
    elseif iscell(matlabbatch)
        if numel(matlabbatch)==1
            conn_setup_preproc_disp(matlabbatch{1},str);
        else
            for n=1:numel(matlabbatch)
                conn_setup_preproc_disp(matlabbatch{n},sprintf('%s(%d)',str,n));
            end
        end
    elseif ~isempty(matlabbatch)&&isnumeric(matlabbatch)
        if numel(matlabbatch)>100
            conn_disp('fprintf','%s = [%s %s %s %s %s ... %s %s %s %s %s]\n',str,mat2str(matlabbatch(1)),mat2str(matlabbatch(2)),mat2str(matlabbatch(3)),mat2str(matlabbatch(4)),mat2str(matlabbatch(5)),mat2str(matlabbatch(end-4)),mat2str(matlabbatch(end-3)),mat2str(matlabbatch(end-2)),mat2str(matlabbatch(end-1)),mat2str(matlabbatch(end)));
        else
            conn_disp('fprintf','%s = %s\n',str,mat2str(matlabbatch));
        end
    elseif ~isempty(matlabbatch)&&ischar(matlabbatch)
        if numel(matlabbatch)>200, conn_disp('fprintf','%s = %s...%s\n',str,1,matlabbatch(1:50),matlabbatch(end-50+1:end));
        else conn_disp('fprintf','%s = %s\n',str,matlabbatch);
        end
    elseif numel(matlabbatch)>1
        for n=1:numel(matlabbatch)
            conn_setup_preproc_disp(matlabbatch(n),sprintf('%s(%d)',str,n));
        end
    elseif isstruct(matlabbatch)
        names=fieldnames(matlabbatch);
        for n=1:numel(names)
            conn_setup_preproc_disp(matlabbatch.(names{n}),sprintf('%s.%s',str,names{n}));
        end
    end
end
end