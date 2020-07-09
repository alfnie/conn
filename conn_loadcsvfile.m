function data = conn_loadcsvfile(tfilename,doheader)
if nargin<2||isempty(doheader), doheader=true; end

tdata=regexp(fileread(tfilename),'[\n\r]+','split');
tdata=tdata(cellfun('length',tdata)>0);
tdata=regexp(tdata,',','split');
nfields=cellfun('length',tdata);
Nfields=mode(nfields);
tdata=tdata(nfields==Nfields);
tdata=cat(1,tdata{:});
tdata=regexprep(tdata,'^\s+|\s+$','');
ndata=cellfun('length',tdata)>0;
tdata=tdata(any(ndata,2),any(ndata,1)); % eliminates empty rows&cols
isndata=cellfun('length',regexpi(tdata,'^\s*([\-\+\d\.]+|nan)\s*$'))|~cellfun('length',tdata);

fieldnames=tdata(1,:);
for n=1:numel(fieldnames),
    if doheader
        name=regexprep(fieldnames{n},'[^\w\d]','');
        val=tdata(2:end,n);
        if all(isndata(2:end,n)), val=str2double(val); end
    else
        name=sprintf('field%d',n);
        val=tdata(1:end,n);
        if all(isndata(1:end,n)), val=str2double(val); end
    end
    try, validname=isvarname(name); catch, validname=true; end
    assert(validname, 'column header %s is not a valid variable name (must begin with letter and be followed by letters/digits/underscores only)',name);
    data.(name)=val;
end
end
