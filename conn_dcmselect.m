function i=conn_dcmselect(filenameout,filedescripout,folder,str,i,DOGUI)
% internal function

if nargin<6||isempty(DOGUI), DOGUI=false; end
if nargin<5||isempty(i), i=[]; end
if nargin<4||isempty(str), str='FUNCTIONAL'; end
if nargin<3||isempty(folder), folder='func'; end

if DOGUI, i=listdlg('PromptString',{sprintf('select %s runs if you wish to copy these runs to ./%s/ folder, or prease Cancel to skip this step',str,folder)},'SelectionMode','multiple','ListSize',[800 400],'ListString',filedescripout,'InitialValue',i); end
for n1=1:numel(i),
    tfilenameout=cellstr(filenameout{i(n1)});
    [tpath,tname,text]=fileparts(tfilenameout{1});
    [ok,nill]=mkdir(tpath,folder);
    for n2=1:numel(tfilenameout)
        [nill,tname,text]=fileparts(tfilenameout{n2});
        for text2={text,'.json','.mat'}
            if ispc, [nill,nill]=system(sprintf('copy "%s" "%s"',  conn_prepend('',tfilenameout{n2},text2{1}),fullfile(tpath,folder,[tname,text2{1}])));
            else     [nill,nill]=system(sprintf('cp ''%s'' ''%s''',conn_prepend('',tfilenameout{n2},text2{1}),fullfile(tpath,folder,[tname,text2{1}])));
            end
        end
    end
end
end
