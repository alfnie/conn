function out=conn_loadcfgfile(filename,out);

% reads .cfg file into matlab structure
% out=conn_loadcfgfile(filename [,out] );
%   filename : .cfg file name
%   out      : output matlab structure
%
% e.g.
%  out=conn_loadcfgfile('example.cfg')
%  with example.cfg
%    #a
%    1 2 3
%    4 5 6
%    #b
%    one line
%    another line
%  will output
%   out.a=[1 2 3;4 5 6]
%   out.b={'one line';'another line'}
%

if nargin<2, out=struct; end
if isstruct(filename), for f=reshape(fieldnames(filename),1,[]), out=setfield(out,f{1},getfield(filename,f{1})); end; return; end
if any(conn_server('util_isremotefile',filename)), out=conn_server('run',mfilename,conn_server('util_localfile',filename),out); return; end
filename=char(filename);

fieldname={'arg'};
s=fileread(filename);
s=strtrim(regexp(s,'\n','split'));
s=s(cellfun('length',s)>0);
for n1=1:length(s),
    txt=strtrim(regexprep(s{n1},'\%.*',''));
    if numel(txt)>=1&&txt(1)=='%', txt=[]; % comment
    elseif numel(txt)>=1&&txt(1)=='#', % field name
        idx=regexp(txt,'\s+');
        if ~isempty(idx),
            fieldname=txt(2:idx(1)-1);
            txt=txt(idx(1)+1:end);
        else
            fieldname=txt(2:end);
            txt=[];
        end
        fieldname=regexp(fieldname,'\.','split');
        out=setfield(out,fieldname{:},[]);
    end
    if ~isempty(txt), % field value
        if ~any(ismember(regexprep(lower(txt),'nan|inf',''),'abcdefghijklmnopqrstuvwxyz')), [n,ok]=str2num(txt);
        else ok=false;
        end
        if ok, newvalue=n; else newvalue={txt}; end
        temp=[];
        try, temp=getfield(out,fieldname{:}); end
        try
            if ~isempty(temp),
                out=setfield(out,fieldname{:},cat(1,temp,newvalue));
            else
                out=setfield(out,fieldname{:},newvalue);
            end
        catch
            if isempty(temp), error('%s\nfile %s line %d\n problem entering %s in %s field (line %d)\n',which(mfilename),filename,n1,txt,sprintf('%s',fieldname{:}));
            else error('%s\nfile %s line %d\n problem concatenating entry "%s" (size %s) into existing field #%s (size %s)\n',which(mfilename),filename,n1,txt,mat2str(size(newvalue)),sprintf('%s.',fieldname{:}),mat2str(size(temp)));
            end
        end
    end
end
end

