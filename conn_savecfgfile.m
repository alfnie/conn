function conn_savecfgfile(filename,out);
% writes .cfg file from matlab structure or cell
% conn_savecfgfile(filename, out);
%   filename : .cfg file name
%   out      : matlab structure (or cell with format {fieldname1,fieldvalue2,...})
%
% e.g.
%  conn_savecfgfile('example.cfg',struct('a',[1 2 3;4 5 6],'b',{{'one line';'another line'}}))
%  will output file example.cfg
%    #a
%    1 2 3
%    4 5 6
%    #b
%    one line
%    another line
%

if any(conn_server('util_isremotefile',filename)), conn_server('run',mfilename,conn_server('util_localfile',filename),out); return; end
filename=conn_server('util_localfile',filename);

fh=fopen(filename,'wt');
if iscell(out),
    names=out(1:2:end-1);
    values=out(2:2:end);
else
    [names,values]=structlist(out);
end
for n=1:numel(names)
    fprintf(fh,'\n#%s\n',names{n});
    var=values{n};
    if iscell(var)
        for r=1:numel(var), fprintf(fh,'%s\n',var{r}); end
    elseif ischar(var)
        fprintf(fh,'%s\n',var); 
    else
        for m=1:size(var,1), fprintf(fh,'%s\n',num2str(var(m,:))); end
    end
end
fclose(fh);
end

function [names,values]=structlist(out)
names={};
values={};
fnames=fieldnames(out);
for n=1:numel(fnames),
    if isstruct(out.(fnames{n})), 
        [tn,tv]=structlist(out.(fnames{n}));
        tn=cellfun(@(x)[fnames{n} '.' x],tn,'uni',0);
    else
        tn={fnames{n}};
        tv={out.(fnames{n})};
    end
    names=[names;tn(:)];
    values=[values;tv(:)];
end
end