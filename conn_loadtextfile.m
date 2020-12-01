function tdata = conn_loadtextfile(tfilename,okstruct)
% loads text/mat file (.mat .txt .csv .tsv .json)
% data = conn_loadtextfile(filename)
%
% note: to read one specific field in a file
%  a) use '<filename>,<fieldname>' format in filename                   e.g. conn_loadtextfile('demographics.tsv,age')
%  or b) use data = conn_loadtextfile(filename,fieldname) syntax        e.g. conn_loadtextfile('demographics.tsv','age')
%
% see spm_load

if nargin<2||isempty(okstruct), okstruct=true; end
if ischar(okstruct), v=okstruct; okstruct=true; 
elseif ischar(tfilename)&&any(tfilename==',')
    tfields=regexp(tfilename,',','split');
    assert(numel(tfields)==2,'unable to interpret filename %s',tfields);
    tfilename=tfields{1};
    v=tfields{2};
else v='';
end
try, 
    if ischar(tfilename)&&~isempty(tfilename)&&~isempty(regexp(tfilename,'\.mat$')) % bugfix older spm
        tdata=load(tfilename,'-mat');
    else
        tdata=spm_load(tfilename,v);
    end
catch,
    tdata=regexp(fileread(tfilename),'[\r\n]+','split');
    tdata=regexprep(tdata,'\<n/a\>','nan');
    tdata=tdata(cellfun('length',tdata)>0);
    vdata=cellfun(@str2num,tdata,'uni',0);
    if numel(tdata)>1&&isequal(find(cellfun('length',vdata)==0),1), tnames=tdata{1}; tdata=cat(1,vdata{2:end}); 
    else tnames={}; tdata=cat(1,vdata{:}); 
    end
    if ~isempty(v)&&isempty(tnames), error('unable to find header line in %s',tfilename);
    elseif ~isempty(v),
        tnames=regexp(tnames,'[\s,\t]+','split');
        assert(numel(tnames)==size(tdata,2),'unexpected number of entries in header line in %s',tfilename);
        idx=find(strcmp(tnames,v));
        assert(numel(idx)==1, 'unable to find field %s in %s',v,tfilename);
        vdata=vdata(:,idx);
    elseif okstruct&&~isempty(tnames)
        try
            tnames=regexp(tnames,'[\s,\t]+','split');
            tdata=cell2struct(num2cell(tdata,1),tnames,2);
        end
    end        
end
if isstruct(tdata)&&~okstruct,
    tempnames=fieldnames(tdata);
    try, tdata=cell2mat(cellfun(@(n)tdata.(n),tempnames(:)','uni',0));
    catch, tdata=tdata.(tempnames{1});
    end
end
