function conn_savetextfile(tfilename,data,names,descrip)
% conn_savetextfile saves numeric data data to text file
% conn_savetextfile(tfilename,data [,names,descrip])
%

if nargin<4||isempty(descrip), descrip={}; end
if nargin<3||isempty(names), names={}; end
if any(conn_server('util_isremotefile',tfilename)), conn_server('run',mfilename,conn_server('util_localfile',tfilename),data,names,descrip); return; end
tfilename=conn_server('util_localfile',tfilename);

[nill,nill,tfileext]=fileparts(tfilename);
switch(tfileext)
    case '.mat'
        if ~isempty(names)&&~isempty(descrip), save(tfilename,'data','names','descrip');
        elseif ~isempty(names), save(tfilename,'data','names');
        else save(tfilename,'data');
        end
    otherwise,
        if strcmp(tfileext,'.txt'), names=regexprep(names,'\s','');
        else                        names=regexprep(names,'\,','');
        end
        fh=fopen(tfilename,'wt');
        if isstruct(data)
            names=fieldnames(data);
            for n1=1:numel(names),
                fprintf(fh,'%s',names{n1});
                if n1<numel(names)&&strcmp(tfileext,'.csv'), fprintf(fh,',');
                elseif n1<numel(names), fprintf(fh,' ');
                else fprintf(fh,'\n');
                end
            end
            for n2=1:size(data.(names{1}),1),
                for n1=1:numel(names),
                    if iscell(data.(names{n1})), fprintf(fh,'%s',data.(names{n1}){n2});
                    else fprintf(fh,'%s',mat2str(data.(names{n1})(n2)));
                    end
                    if n1<numel(names)&&strcmp(tfileext,'.csv'), fprintf(fh,',');
                    elseif n1<size(data,2), fprintf(fh,' ');
                    else fprintf(fh,'\n');
                    end
                end
            end
        else
            for n1=1:numel(names),
                if isempty(names{n1}), names{n1}='-'; end
                fprintf(fh,'%s',names{n1});
                if n1<numel(names)&&strcmp(tfileext,'.csv'), fprintf(fh,',');
                elseif n1<numel(names), fprintf(fh,' ');
                else fprintf(fh,'\n');
                end
            end
            for n3=1:size(data,3)
                for n2=1:size(data,1),
                    for n1=1:size(data,2),
                        if iscell(data(n2,n1,n3))&&ischar(data{n2,n1,n3}), fprintf(fh,'%s',data{n2,n1,n3});
                        else fprintf(fh,'%s',mat2str(data(n2,n1,n3)));
                        end
                        if n1<size(data,2)&&strcmp(tfileext,'.csv'), fprintf(fh,',');
                        elseif n1<size(data,2), fprintf(fh,' ');
                        else fprintf(fh,'\n');
                        end
                    end
                end
                if n3<size(data,3), fprintf(fh,'#\n'); end
            end
        end
        fclose(fh);
end