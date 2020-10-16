function [dataset,filter]=conn_bidsdir(filenames,varargin)
% CONN_BIDSDIR lists contents of BIDS folder 
%
% dataset=conn_bidsdir(rootfolder)
% 
% e.g. conn_bidsdir /data/myconnectome
%
%

if nargin<1||isempty(filenames), filenames=pwd; end

% bids filename spec 2018/7
%spec={'sub','ses','task','acq','ce','rec','dir','run','mod','echo','recording','proc'};
spec={'sub','ses','task','acq','ce','rec','dir','run','mod','echo','recording','proc','space','res','desc'};
filter={};

if isstruct(filenames)
    dataset=filenames;
    if ~isempty(varargin)&&isequal(varargin{1},'noempty'), noempty=true; varargin=varargin(2:end); else noempty=false; end
    if ~isempty(varargin)&&isequal(varargin{1},'forcefirst'), forcefirst=true; varargin=varargin(2:end); else forcefirst=false; end
    filter=varargin;
    filtervalid=true(size(filter));
    select=true;
    for n=1:2:numel(varargin)-1
        if isempty(varargin{n+1}), temp=true(size(dataset.data.(varargin{n}))); 
        else temp=ismember(dataset.data.(varargin{n}),varargin{n+1});
        end
        if ~isempty(temp)&&~any(temp)&&~isempty(varargin{n+1}) % special case regexp wildcards or #EMPTY keyword
            if isequal(varargin{n+1},'#EMPTY'), temp=select&cellfun('length',dataset.data.(varargin{n}))==0;
            elseif ischar(varargin{n+1})&&any(varargin{n+1}=='*'), temp=select&cellfun('length',regexp(dataset.data.(varargin{n}),varargin{n+1}));
            elseif iscell(varargin{n+1})&&numel(varargin{n+1})==1&&ischar(varargin{n+1}{1})&&any(char(varargin{n+1}{1})=='*'), temp=select&cellfun('length',regexp(dataset.data.(varargin{n}),varargin{n+1}{1}));
            end
        end
        temp=select&temp;
        if forcefirst&&n>1,
            if any(temp~=select), filtervalid([n,n+1])=false; end
        elseif noempty, % skip filter conditions that lead to empty selection
            if any(temp), select=temp;
            else filtervalid([n,n+1])=false;
            end
        else
            select=temp;
        end
    end
    filter=filter(filtervalid);
    if any(~select)
        fnames=fieldnames(dataset.data);
        for n=1:numel(fnames)
            dataset.data.(fnames{n})=dataset.data.(fnames{n})(select);
        end
    end
else
    if iscell(filenames), files_all=conn_sortfilenames(filenames);
    else
        dirs_all=conn_dir(fullfile(filenames,'sub-*'),'-cell','-R','-sort','-dir');
        files_all={};
        for n=1:numel(dirs_all)
            tfiles_all=conn_dir(fullfile(dirs_all{n},'sub-*'),'-cell','-inf','-sort');
            if ~isempty(tfiles_all), files_all=[files_all, reshape(tfiles_all,1,[])]; end
        end
        if isempty(files_all),fprintf(sprintf('warning: no sub-* files found in %s',filenames));dataset=[];return;end
    end
    isgz=cellfun('length',regexp(files_all,'\.gz$'))>0;
    files_all=regexprep(files_all,'\.gz$','');
    [files_all,nill,idx]=unique(files_all,'legacy'); isgz=accumarray(idx(:),isgz(:),[],@min)>0; % note: disregards [filename.ext].gz if [filename.ext] exists
    [files_all_path,files_all_name,files_all_ext]=cellfun(@fileparts,files_all,'uni',0);
    files_all_ext(isgz)=cellfun(@(x)[x,'.gz'],files_all_ext(isgz),'uni',0);
    files_all(isgz)=cellfun(@(x)[x,'.gz'],files_all(isgz),'uni',0);
    
    str='^';for n=1:numel(spec), if n>1, str=[str,'(_']; else str=[str,'(']; end; str=[str,spec{n},'-[^_\.]+)?']; end; str=[str,'(_[^_]+)*?$'];
    files_parts=regexp(files_all_name,str,'tokens','once');
    files_parts=cat(1,files_parts{:});
    dataset.data.file=files_all(:);
    dataset.data.description=regexprep(arrayfun(@(n)[files_parts{n,[size(files_parts,2),2:size(files_parts,2)-1]}],1:size(files_parts,1),'uni',0),'_+',' ');
    dataset.data.contents=regexprep(files_parts(:,end),'^_+','');
    [nill,dataset.data.folder]=cellfun(@fileparts,files_all_path(:),'uni',0);
    dataset.data.format=regexprep(files_all_ext(:),'^\.+','');
    dataset.data.series=regexprep(cellfun(@(a,b,c)[a,b,c],...
        files_parts(:,find(ismember(spec,'ses'),1)),...
        files_parts(:,find(ismember(spec,'task'),1)),...
        files_parts(:,find(ismember(spec,'run'),1)),...
        'uni',0),'^_+','');
    for n=1:numel(spec),
        dataset.data.(spec{n})=regexprep(files_parts(:,n),['^_?',spec{n},'-'],'');
    end
end
for n=[{'file','description','contents','folder','format','series'},spec] 
    [dataset.dict.(n{1}),nill,idx1]=unique(dataset.data.(n{1})(cellfun('length',dataset.data.(n{1}))>0)); 
    if isempty(dataset.dict.(n{1})), dataset.dict.(n{1})={}; 
    else
        [nill,idx2]=sort(accumarray(idx1(:),1),'descend');
        dataset.dict.(n{1})=reshape(dataset.dict.(n{1})(idx2),1,[]);
    end
end
