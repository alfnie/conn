function V=conn_vol(filename,docache)

if nargin>1&&~isempty(docache)&&docache, opts={'-cache'}; else opts={}; end
V=struct; 
conn_loadmatfile(filename,'V',opts{:});
V.fname=filename;
V.overwritesoftlink=true;
if isfield(V,'softlink')&&~isempty(V.softlink), 
    V.softlink=regexprep(V.softlink,'^cachetmp_',''); 
    str1=regexp(V.fname,'Subject\d+','match'); if ~isempty(str1), V.softlink=regexprep(V.softlink,'Subject\d+',str1{end}); end
end


