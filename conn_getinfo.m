function [nV,str,icon,filename]=conn_getinfo(filename,doconvert)

if nargin<2||isempty(doconvert), doconvert=true; end
str=[]; icon=[]; V=[]; nV=[];
if isempty(filename), return; end
[pathname,name,ext]=spm_fileparts(filename(1,:));
if doconvert
    if any(strcmp(ext(1:min(4,length(ext))),{'.mgh','.mgz'}))
        filename=conn_mgh2nii(filename);
        [pathname,name,ext]=spm_fileparts(filename(1,:));
    end
    if any(strcmp(ext(1:min(3,length(ext))),{'.gz'}))
        filename=conn_gz2nii(filename);
        [pathname,name,ext]=spm_fileparts(filename(1,:));
    end
    if any(strcmp(ext(1:min(6,length(ext))),{'.annot'}))
        filename=conn_annot2nii(filename);
        [pathname,name,ext]=spm_fileparts(filename(1,:));
    end
    if any(strcmp(ext(1:min(4,length(ext))),{'.dcm'}))
        filename=conn_dcm2nii(filename);
        [pathname,name,ext]=spm_fileparts(filename(1,:));
    end
    if any(strcmp(ext(1:min(4,length(ext))),{'.gii'}))
        filename=conn_surf_gii2nii(filename);
        [pathname,name,ext]=spm_fileparts(filename(1,:));
    end
end
if size(filename,1)>=1&&isdir(deblank(filename(1,:)))
    nV=size(filename,1);
    str=[{'Directory'},reshape(cellstr(filename),1,[])];
%     try
%         nfiles=numel(conn_dir(filename,'-R','-cell'));
%         nsubdirs=numel(conn_dir(filename,'-R','-cell','-dir','-skipdot'));
%         str{end+1}=sprintf('%d files; %d subdirectories',nfiles,nsubdirs);
%     end
else
    switch(ext(1:min(4,length(ext)))),
        case {'.mgh','.mgz'}
            nV=1;
            str={'MGH file'};
        case {'.gz'}
            nV=1;
            str={'compressed file'};
        case {'.annot'}
            nV=1;
            str={'Freesurfer annotation file'};
        case {'.gii'}
            nV=1;
            str={'Freesurfer gifti file'};
        case {'.dcm'}
            nV=1;
            str={'DICOM file'};
        case {'.img','.hdr','.nii'},
            %maxvolsdisplayed=2;
            issinglevolume=cellfun('length',regexp(cellstr(filename),',\d+$'))>0;
            nV=sum(issinglevolume);
            ok=false;
            if all(issinglevolume)
                try
                    if nV>1
                        V1=spm_vol([deblank(filename(1,:))]);
                        V2=spm_vol([deblank(filename(end,:))]);
                        V=[V1 V2];
                        ok=true;
                    end
                end
            elseif size(filename,1)==1
                nfilename=nifti(filename);
                nV=0; for n=1:numel(nfilename), tV=size(nfilename(n).dat,4); nV=nV+tV; end
                try
                    if nV>10
                        V1=spm_vol([deblank(filename),',1']);
                        V2=spm_vol([deblank(filename),',',num2str(nV)]);
                        V=[V1 V2];
                        ok=true;
                    end
                end
            else
                try
                    V1=spm_vol([deblank(filename(1,:))]);
                    V2=spm_vol([deblank(filename(end,:))]);
                    V=[V1(1) V2(end)];
                    nfilename=nifti(cellstr(filename(~issinglevolume,:)));
                    for n=1:numel(nfilename), tV=size(nfilename(n).dat,4); nV=nV+tV; end
                    ok=true;
                end
            end
            if ~ok, V=spm_vol(filename); nV=numel(V); end
            %V=spm_vol(filename);
            if length(V)==1, icon=V;
            else icon=V([1,end]);
                %else icon=V(reshape(unique(round(linspace(1,numel(V),maxvolsdisplayed))),1,[]));
            end
        case {'.tal','.mat','.txt','.par','.1d','.csv','.tsv'},
            nV=size(filename,1);
            for n1=1:nV,
                x=conn_loadtextfile(deblank(filename(n1,:)));
                if isstruct(x), 
                    names=fieldnames(x); 
                    try, x=cell2mat(cellfun(@(n)x.(n),names(:)','uni',0)); names=deblank(sprintf('%s ',names{:})); 
                    catch, names=names{1}; x=x.(names); 
                    end
                else names=''; 
                end
                tok=false;
                try
                    if strcmp(names,'CONN_x')&&isstruct(x)
                        V(n1).dim=x.Setup.nsubjects;
                        V(n1).fname=deblank(filename(n1,:));
                        temp=spm_vol(x.Setup.structural{1}{1}{1});
                        if x.Setup.nsubjects>1, temp=[temp spm_vol(x.Setup.structural{x.Setup.nsubjects}{1}{1})]; end
                        icon=cat(2,icon,temp);
                        tok=true;
                    elseif strcmp(names,'SPM')&&isstruct(x)
                        V(n1).dim=size(x.xX.X);
                        V(n1).fname=deblank(filename(n1,:));
                        %temp=x.xX.X;
                        %icon=cat(2,icon,temp(round(linspace(1,size(temp,1),128)),:));
                        thisicon=struct('X',x.xX.X,'conditions',[],'sessions',1,'functional',0,'VY',[]);
                        if isfield(x,'xY')&&isfield(x.xY,'P'), thisicon.functional=size(x.xY.P,1);
                        elseif isfield(x,'xY')&&isfield(x.xY,'VY'), thisicon.functional=numel(x.xY.VY);
                        end
                        if isfield(x,'xY')&&isfield(x.xY,'VY'), thisicon.VY=x.xY.VY([1 end]); end
                        if isfield(x,'Sess')
                            thisicon.sessions=numel(x.Sess);
                            conditions={};
                            for n2=1:numel(x.Sess)
                                conditions=[conditions,[x.Sess(n2).U.name]];
                            end
                            thisicon.conditions=conditions;
                        end
                        if isempty(icon), icon=thisicon;
                        else icon=cat(2,icon,thisicon);
                        end
                        tok=true;
                    end
                end
                if ~tok
                    V(n1).dim=size(x);
                    V(n1).fname=deblank(filename(n1,:));
                    if isnumeric(x) && (n1==1 || n1==size(filename,1)),
                        icon=cat(2,icon,x(:,:));
                    end
                end
            end
        otherwise,
            error(['File type ',ext,' not implemented']);
    end
end
if ~isempty(V), 
	str=['[',num2str(nV), ' file']; if nV>1, str=[str,'s']; end; str={[str,']']};
	if nV>1 && any(any(detrend(cat(1,V(:).dim),0))),str{end}=[str{end} ' Dimensions: NOT MATCHED'];
	else str{end}=[str{end} ' x [size',sprintf(' %1.0f ',V(1).dim),']']; end
	if nV==1, str{end+1}=V(1).fname;
    else
        str{end+1}=['First: ',V(1).fname]; str{end+1}=['Last : ',V(end).fname];
        if isfield(V,'n')&&(V(1).n(1)>1||V(end).n(1)>1), str{end-1}=[str{end-1},',',num2str(V(1).n(1))]; str{end}=[str{end},',',num2str(V(end).n(1))]; end
    end
	for n1=1:length(str), if length(str{n1})>30+9, str{n1}=[str{n1}(1:4),' ... ',str{n1}(end-30+1:end)]; end; end; %str{n1}=[str{n1}(1:17),' ... ',str{n1}(end-17+1:end)]; end; end
end
end


