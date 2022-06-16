function conn_exportlist(cbo,filename,header,onlyselected,dodisp,dotabs)
% conn_exportlist
% exports uicontrol string

if nargin<1||isempty(cbo),cbo=gcbo; end
if nargin<2||isempty(filename)
    [filename,filepath]=uiputfile({'*.txt','*.txt (text file)'},'Save table as');
    if ~ischar(filename), return; end
    filename=fullfile(filepath,filename);
end
[filepath,filename,fileext]=fileparts(filename);
if nargin<6||isempty(dodisp),dodisp=true; end
if nargin<5||isempty(dotabs),dotabs=false; end

str=cellstr(get(cbo,'string'));
if nargin>=4&&~isempty(onlyselected), str=str(get(cbo,'value')); end
str=regexprep(str,'<[^<>]*>','');
if dotabs,
    header=regexprep(header,'\s(\s)+','\t'); 
    str=regexprep(str,'\s(\s)+','\t'); 
end
if nargin>=3&&~isempty(header), str=[{header} reshape(str,1,[])]; end
conn_fileutils('filewrite',fullfile(filepath,[filename,fileext]), str);
%fh=fopen(fullfile(filepath,[filename,fileext]),'wt');
%if nargin>=3&&~isempty(header), fprintf(fh,'%s\n',header); end
%for n=1:numel(str),
%    fprintf(fh,'%s\n',str{n});
%end
%fclose(fh);
if dodisp, fprintf('Table exported to %s\n',fullfile(filepath,[filename,fileext])); end
end