function [files,idx,iidx]=conn_sortfilenames(files)
% sorts strings with numbers

if isempty(files), return; end
ischarfilename=ischar(files);
if ischarfilename, files=cellstr(files); end

n=max(cellfun('length',regexprep(files,'[^\d]+','')));
tfiles=regexprep(files,'\d+',['${[repmat(''0'',1,max(0,',num2str(n),'-numel($0))) $0]}']);
[nill,idx]=sort(tfiles);
files=files(idx);
if ischarfilename, files=char(files); end
if nargout>2, iidx=idx; iidx(idx)=1:numel(idx); end
