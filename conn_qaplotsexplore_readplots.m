function [dataA,labelA,infoA,lineA,txtA,miny,maxy]=conn_qaplotsexplore_readplots(files,files_txt,in,usubjects, dogui, skipall);
if nargin<5||isempty(dogui), dogui=true; end
if nargin<6||isempty(skipall), skipall=false; end
if any(conn_server('util_isremotefile',files)), 
    hmsg=conn_msgbox(sprintf('Loading %d images. Please wait...',numel(in)),'',-1);
    [dataA,labelA,infoA,lineA,txtA,miny,maxy]=conn_server('run',mfilename,conn_server('util_localfile',files),conn_server('util_localfile',files_txt),in,usubjects, false, skipall); 
    if ishandle(hmsg), delete(hmsg); end
    return; 
end
files=conn_server('util_localfile',files);
files_txt=conn_server('util_localfile',files_txt);

if dogui, ht=conn_waitbar(0,sprintf('Loading %d plots. Please wait...',numel(in)),false); end
dataA={};
labelA={};
infoA={};
lineA={};
%dataB={};
dopull=conn_server('util_isremotefile',files(in));
if any(dopull), % deprecated
    files(in(dopull))=conn_cache('pull',files(in(dopull)));
    files_txt(in(dopull))=conn_prepend('',files(in(dopull)),'.txt');
end
for n=1:numel(in),
    data=conn_loadmatfile(files{in(n)});
    if isfield(data,'results_patch'), dataA{n}=data.results_patch; end
    %dataB{n}=data.results_label;
    descr=''; try, descr = conn_fileutils('fileread',files_txt{in(n)}); end
    if isfield(data,'results_label'), labelA{n}=data.results_label;
    else labelA{n}='';
    end
    if isfield(data,'results_info'), infoA{n}=data.results_info;
    else infoA{n}=[];
    end
    if isfield(data,'results_line'), lineA{n}=data.results_line;
    else lineA{n}=[];
    end
    if isempty(descr), txtA{n}={'[empty]'};
    else txtA{n}=regexp(descr,'\n+','split');
    end
    if dogui, conn_waitbar(n/numel(in),ht); end
end
if isempty(dataA)&&~isempty(infoA) % refresh plot info
    X=[];Xnames={};Xdescr={};Xsub=[];Xdir=[];
    inok=false(numel(in),1);
    for n=1:numel(infoA)
        if ~isempty(infoA{n}.Variables)
            if isempty(Xnames), Xnames=infoA{n}.Variables; Xdescr=infoA{n}.Variables_descr; inok(n)=true;
            elseif ~isequal(infoA{n}.Variables,Xnames), conn_disp('fprintf','warning %s mismatch variable names %s (expected %s)\n',files{in(n)},sprintf('%s ',Xnames{:}),sprintf('%s ',infoA{n}.Variables{:}));
            else inok(n)=true;
            end
        end
        if inok(n), 
            X=cat(1,X,infoA{n}.Values); 
            Xsub=cat(1,Xsub,infoA{n}.Subjects); 
            if isempty(Xdir)&&isfield(infoA{n},'Variables_dir'), Xdir=infoA{n}.Variables_dir; end
        end
    end
    tX=X;X=nan(numel(in),size(tX,2));X(inok,:)=tX;
    [dataA,labelA,infoA,lineA]=conn_qaplots_covupdate(X,Xnames,Xdescr,Xsub,Xdir); %usubjects);
end
assert(numel(dataA)==numel(in) && all(cellfun('length',dataA)>0),'missing information in plot files');
miny=inf(1,numel(dataA{1})); maxy=-inf(1,numel(dataA{1}));
for n=1:numel(dataA), miny=min(miny,cellfun(@(x)min(x(:)),dataA{n})); maxy=max(maxy,cellfun(@(x)max(x(:)),dataA{n})); end
if dogui, conn_waitbar('close',ht); end
end
