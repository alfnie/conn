function data = conn_loadcsvfile(tfilename,doheader, showlines)
if nargin<2||isempty(doheader), doheader=true; end
if nargin<3||isempty(showlines), showlines=false; end
if any(conn_server('util_isremotefile',tfilename)), data=conn_server('run',mfilename,conn_server('util_localfile',tfilename),doheader); return; end
tfilename=conn_server('util_localfile',tfilename);

tdata=regexp(fileread(tfilename),'[\n\r]+','split');
tdata=tdata(cellfun('length',tdata)>0);
tdata=regexp(tdata,'[,\t]','split');
nfields=cellfun('length',tdata);
Nfields=mode(nfields);
if showlines&&~all(nfields==Nfields), fprintf('warning: not all lines contain the same number of fields (lines %s)\n',mat2str(find(nfields~=Nfields))); end
tdata=tdata(nfields<=Nfields);
tdata=cellfun(@(x)[x cell(1,Nfields-numel(x))],tdata,'uni',0);
tdata=cat(1,tdata{:});
tdata=regexprep(tdata,'^\s+|\s+$','');
ndata=cellfun('length',tdata)>0;
tdata=tdata(any(ndata,2),any(ndata,1)); % eliminates empty rows&cols
isndata=cellfun('length',regexpi(tdata,'^\s*([\-\+\d\.e]+|nan)\s*$'))|~cellfun('length',tdata);

fieldnames=tdata(1,:);
for n=1:numel(fieldnames),
    if doheader
        name=regexprep(fieldnames{n},{'\([^\(\)]*\)','[^\w\d]'},'');
        val=tdata(2:end,n);
        if all(isndata(2:end,n)), val=str2double(val); end
    else
        name=sprintf('field%d',n);
        val=tdata(1:end,n);
        if all(isndata(1:end,n)), val=str2double(val); end
    end
    try, validname=isvarname(name); catch, validname=true; end
    %if ~validname, fprintf('warning: column header %s is not a valid variable name (must begin with letter and be followed by letters/digits/underscores only)\n',name); name=sprintf('var%d',n); end
    assert(validname, 'column header %s is not a valid variable name (must begin with letter and be followed by letters/digits/underscores only)',name);
    data.(name)=val;
end
end
