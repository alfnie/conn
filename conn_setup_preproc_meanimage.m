function [fileout,filetested,fileother]=conn_setup_preproc_meanimage(filename,filetype)
% conn_setup_preproc_meanimage(filename)
% finds "mean", "norm_spm12", "norm_spm8", and "json" files in same location as functional data volume

if nargin<2||isempty(filetype), filetype='mean'; end
fileother=[];
switch(filetype)
    case 'mean'
        [file_path,file_name,file_ext,file_num]=spm_fileparts(filename);
        % [PREFIX r BASENAME] -> [PREFIX(minus 'a' or 's') mean BASENAME]
        idx1=find(file_name(1:end-1)=='r');
        ok1=false(size(idx1));
        str1=cell(size(idx1));
        for n=1:numel(idx1),
            str1{n}=fullfile(file_path,[regexprep(file_name(1:idx1(n)-1),'[as]','') 'mean' file_name(idx1(n)+1:end) file_ext]);
            if conn_existfile(str1{n}), ok1(n)=true; end
        end
        % [PREFIX u BASENAME] -> [PREFIX(minus 'a' or 's') meanu BASENAME]
        idx2=find(file_name(1:end-1)=='u');
        ok2=false(size(idx2));
        str2=cell(size(idx2));
        for n=1:numel(idx2),
            str2{n}=fullfile(file_path,[regexprep(file_name(1:idx2(n)-1),'[as]','') 'mean' file_name(idx2(n):end) file_ext]);
            if conn_existfile(str2{n}), ok2(n)=true; end
        end
        % [PREFIX BASENAME] -> [PREFIX(minus 'a' or 's') mean|art_mean_ BASENAME]
        idx3=1:numel(file_name);
        ok3=false(size(idx3));
        str3=cell(size(idx3));
        for n=1:numel(idx3),
            str3{n}=fullfile(file_path,[regexprep(file_name(1:idx3(n)-1),'[as]','') 'art_mean_' file_name(idx3(n):end) file_ext]);
            if conn_existfile(str3{n}), ok3(n)=true;
            else
                str3{n}=fullfile(file_path,[regexprep(file_name(1:idx3(n)-1),'[as]','') 'mean' file_name(idx3(n):end) file_ext]);
                if conn_existfile(str3{n}), ok3(n)=true;
                end
            end
        end
        i1=find(ok1,1);
        i2=find(ok2,1);
        i3=find(ok3,1);
        min1=min([inf idx1(i1)]);
        min2=min([inf idx2(i2)]);
        min3=min([inf idx3(i3)]);
        [nill,i]=min([min3 min1+1 min2]);
        if isinf(nill),   fileout=[];
        elseif i==1,      fileout=str3{i3};
        elseif i==2,      fileout=str1{i1};
        elseif i==3,      fileout=str2{i2};
        end
        str=[str1 str2 str3];
        if isempty(str)||isempty(fileout), filetested=filename;
        else filetested=str{1};
        end
        
    case 'norm_spm12'
        [file_path,file_name,file_ext,file_num]=spm_fileparts(filename);
        % [PREFIX wc0 BASENAME] -> [y_ BASENAME]
        idx1=strfind(file_name(1:end-3),'wc0');
        ok1=false(size(idx1));
        str1=cell(size(idx1));
        for n=1:numel(idx1),
            str1{n}=fullfile(file_path,['y_' file_name(idx1(n)+3:end) '.nii']);
            if conn_existfile(str1{n}), ok1(n)=true; end
        end
        % [PREFIX w BASENAME] -> [y_ BASENAME]
        idx2=find(file_name(1:end-1)=='w');
        ok2=false(size(idx2));
        str2=cell(size(idx2));
        for n=1:numel(idx2),
            str2{n}=fullfile(file_path,['y_' file_name(idx2(n)+1:end) '.nii']);
            if conn_existfile(str2{n}), ok2(n)=true; end
        end
        i1=find(ok1,1);
        i2=find(ok2,1);
        min1=min([inf idx1(i1)]);
        min2=min([inf idx2(i2)]);
        [nill,i]=min([min1 min2+1]);
        if isinf(nill),   fileout=[]; fileother=[];
        elseif i==1,      fileout=str1{i1}; fileother=fullfile(file_path,[file_name(idx1(i1)+3:end) '.nii']);
        elseif i==2,      fileout=str2{i2}; fileother=fullfile(file_path,[file_name(idx2(i2)+1:end) '.nii']);
        end
        str=[str1 str2];
        if isempty(str)||isempty(fileout), filetested=filename;
        else filetested=str{1};
        end

    case 'norm_spm8'
        [file_path,file_name,file_ext,file_num]=spm_fileparts(filename);
        % [PREFIX wc0 BASENAME] -> [BASENAME _seg_sn.mat]
        idx1=strfind(file_name(1:end-3),'wc0');
        ok1=false(size(idx1));
        str1=cell(size(idx1));
        for n=1:numel(idx1),
            str1{n}=fullfile(file_path,[file_name(idx1(n)+3:end) '_seg_sn.mat']);
            if conn_existfile(str1{n}), ok1(n)=true; end
        end
        % [PREFIX w BASENAME] -> [BASENAME _seg_sn.mat]
        idx2=find(file_name(1:end-1)=='w');
        ok2=false(size(idx2));
        str2=cell(size(idx2));
        for n=1:numel(idx2),
            str2{n}=fullfile(file_path,[file_name(idx2(n)+1:end) '_seg_sn.mat']);
            if conn_existfile(str2{n}), ok2(n)=true; end
        end
        i1=find(ok1,1);
        i2=find(ok2,1);
        min1=min([inf idx1(i1)]);
        min2=min([inf idx2(i2)]);
        [nill,i]=min([min1 min2+1]);
        if isinf(nill),   fileout=[]; fileother=[];
        elseif i==1,      fileout=str1{i1}; fileother=fullfile(file_path,[file_name(idx1(i1)+3:end) file_ext]);
        elseif i==2,      fileout=str2{i2}; fileother=fullfile(file_path,[file_name(idx2(i2)+1:end) file_ext]);
        end
        str=[str1 str2];
        if isempty(str)||isempty(fileout), filetested=filename;
        else filetested=str{1};
        end
    case 'json'
        if ~isempty(regexp(filename,'\.surf\.nii$')), checksurf=true; else checksurf=false; end
        [file_path,file_name,file_ext,file_num]=spm_fileparts(filename);
        file_ext=regexprep(file_ext,'\.surfnii$','.surf.nii');
        % [PREFIX BASENAME] -> [BASENAME.json]
        idx1=1:numel(file_name);
        ok1=zeros(size(idx1));
        str1=cell(size(idx1));
        for n=1:numel(idx1),
            str1{n}=fullfile(file_path,[file_name(idx1(n):end) '.json']);
            if conn_existfile(str1{n}), ok1(n)=1; 
            elseif checksurf&&conn_existfile(regexprep(str1{n},'\.surf.json$','.json')), ok1(n)=2; 
            end
        end
        i1=find(ok1>0,1);
        min1=min([inf idx1(i1)]);
        [nill,i]=min([min1]);
        if isinf(nill),   fileout=[]; fileother=[];
        elseif i==1,      fileout=str1{i1}; fileother=fullfile(file_path,[file_name(idx1(i1):end) file_ext]);
                          if checksurf&&ok1(i1)>1, fileout=regexprep(fileout,'\.surf.json$','.json'); end
        end
        str=[str1];
        if isempty(str)||isempty(fileout), filetested=filename;
        else filetested=str{1};
        end
end
end

