function Series=conn_dcmdir(filenames,DOGUI)
% CONN_DCMDIR lists contents of DICOM folder 
%
% series=conn_dcmdir(filepattern)
% 
% e.g. conn_dcmdir /data/MR*.dcm
%
%

if nargin<2||isempty(DOGUI), DOGUI=[]; end
if nargin<1||isempty(filenames), filenames=fullfile(pwd,'*.dcm'); end
if ischar(filenames), filenames=conn_dir(filenames,'-ls'); end
if ischar(filenames), filenames=cellstr(filenames); end

ESSENTIALS=true; 
OTHERESSENTIALS={'MosaicRefAcqTimes','Private_0019_1029'};
SeriesNumber=nan(1,numel(filenames));
SeriesInfo={};
HDR={};
ContentTime=nan(1,numel(filenames));
ImagePositionPatient=nan(3,numel(filenames));
dict=load(fullfile(fileparts(which('spm')),'spm_dicom_dict.mat'));
filesout_path='';
count=[];
if ~isempty(DOGUI), ht=conn_waitbar(0,sprintf('%s Reading DICOM directory information. Please wait',DOGUI)); end
for n=1:numel(filenames)
    [filesout_path,nill,nill]=fileparts(filenames{n});
    fprintf('.');
    warning('off','spm:dicom');
    a=spm_dicom_header(filenames{n},dict);
    if ESSENTIALS&&~isempty(a), 
        areduced=spm_dicom_essentials(a); 
        for back=OTHERESSENTIALS
            if isfield(a,back{1}), areduced.(back{1})=a.(back{1}); end
        end
    else areduced=a;
    end
    if ~isempty(a)&&isfield(a,'SeriesNumber'), 
        i=a.SeriesNumber;
        SeriesNumber(n)=i;
        try, if isfield(a,'ContentTime'), ContentTime(n)=a.ContentTime; end; end
        try, if isfield(a,'ImagePositionPatient'), ImagePositionPatient(:,n)=a.ImagePositionPatient; end; end
        if numel(SeriesInfo)<i||isempty(SeriesInfo{i})
            SeriesInfo{i}=areduced;
            FieldNames{i}=fieldnames(areduced);
            HDR{n}=a;
        elseif ~isempty(FieldNames{i})
            valid=cellfun(@(s)isfield(areduced,s),FieldNames{i});
            valid(valid)=cellfun(@(s)isequal(areduced.(s),SeriesInfo{i}.(s)),FieldNames{i}(valid));
            valid(ismember(FieldNames{i},OTHERESSENTIALS))=true;
            SeriesInfo{i}=rmfield(SeriesInfo{i},FieldNames{i}(~valid));
            FieldNames{i}=FieldNames{i}(valid);
            HDR{n}=areduced;
        end
        if ~ismember(i,count), 
            if ~isempty(DOGUI), 
                if isfield(a,'SeriesDescription'), tstr=a.SeriesDescription; else tstr=''; end
                conn_waitbar(n/numel(filenames),ht,sprintf('#%d: %s',i,tstr)); 
            end
            count=[count,i]; fprintf('\n'); 
        end
    else
        fprintf('\nSkipping file %s (not a DICOM file)\n',filenames{n});
        if ~isempty(DOGUI), conn_waitbar(n/numel(filenames),ht,sprintf('%s is not a DICOM file',filenames{n})); end
    end
end
fprintf('\n');
if ~isempty(DOGUI), conn_waitbar('close',ht); end
[useries,nill,iseries]=unique(SeriesNumber(~isnan(SeriesNumber)));
if isempty(useries), Series=struct([]); return; end 
Nseries=accumarray(iseries,1);
SeriesNfiles=nan(size(SeriesInfo));
SeriesNfiles(useries)=Nseries;
idxNseries=1:numel(Nseries);%[nill,idxNseries]=sort(Nseries,'descend');
fprintf('Found %d series and %d DICOM files\n',numel(useries),sum(~isnan(SeriesNumber)));

try, fhlog=fopen(fullfile(filesout_path,sprintf('conn_dcmdir_%s.log',datestr(now,'yyyy_mm_dd_HHMMSSFFF'))),'wt'); 
catch, fhlog=[]; 
end
clear Series;
fprintf('\n________________________________________________________\n')
fprintf2('SeriesNumber, SeriesDescription, check, size, volumes, files, Repetition Time (ms)\n'); %, Echo Times (ms)\n');
for n=1:numel(useries),
    i=useries(idxNseries(n));
    idx=find(SeriesNumber==i);
    Series(i)=struct('SeriesNumber',i,'SeriesDescription','','SeriesInfo',SeriesInfo{i},'Size',[],'check',false,'files',{filenames(idx)},'headers',{HDR(idx)});

    fprintf2('%3d',i);
    if isfield(SeriesInfo{i},'SeriesDescription'), fprintf2(',%64s',regexprep(SeriesInfo{i}.SeriesDescription,',','')); Series(i).SeriesDescription=SeriesInfo{i}.SeriesDescription;
    else fprintf2(',%64s','');
    end
    
    s=[];
    if isfield(SeriesInfo{i},'Columns'), s(end+1)=SeriesInfo{i}.Columns; end
    if isfield(SeriesInfo{i},'Rows'), s(end+1)=SeriesInfo{i}.Rows; end
    try
        tmp=ImagePositionPatient(:,idx)';
        tmp=size(unique(tmp,'rows'),1);
        s(end+1)=tmp;
    end
    try
        tmp=ImagePositionPatient(:,idx);
        tmp=sum(all(tmp==repmat(tmp(:,1),[1,size(tmp,2)]),1));
        s(end+1)=tmp;
    end
    Series(i).Size=s;
    if numel(s)==4&&prod(s(3:4))==numel(idx), 
        fprintf2(',ok,%16s',sprintf('%dx%dx%d (%d)',s(1),s(2),s(3),s(4))); Series(i).check=true;
    else
        fprintf2(',err,%16s',mat2str(s)); Series(i).check=false;
    end
    fprintf2(',%4d',numel(idx));
    if isfield(SeriesInfo{i},'RepetitionTime'), fprintf2(',%.2f',SeriesInfo{i}.RepetitionTime); %Series(i).SeriesInfo.dcmdir_RepetitionTime=SeriesInfo{i}.RepetitionTime/1000;
    else fprintf2(',');
    end
    fprintf2('\n');
%     try
%         tmp=ContentTime(idx);
%         Series(i).SeriesInfo.dcmdir_ContentTime=tmp;
%     end
%     if isfield(SeriesInfo{i},'ContentTime'), fprintf2(',%f',SeriesInfo{i}.ContentTime);
%     else
%         try
%             tmp=ContentTime(idx);
%             fprintf2('%.02f ',unique(tmp));
%         end
%         fprintf2(',');
%     end
end
for n=1:numel(useries)
    i=useries(idxNseries(n));
    idx=find(SeriesNumber==i);
    fprintf2('\n  DICOM series #%d (%d dcm files)\n',i,numel(idx));
    str=evalc('disp(HDR{idx(1)});');
    fprintf2('%s',str);
end
try, fclose(fhlog); end

    function fprintf2(varargin)
        fprintf(varargin{:});
        try, fprintf(fhlog,varargin{:}); end
    end
end


