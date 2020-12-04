function [outstruct,filenameout,filedescripout,filetypeout,fileRTout,Series]=conn_dcmconvert(filenames,varargin)
% CONN_DCMCONVERT converts contents of DICOM folder 
%
% series=conn_dcmconvert(filepattern)
%   converts all dicom files matching filepattern to run-#.nii 
%   nifti files (separated into individual dicom runs)
%
% conn_dcmconvert(...,optionname,optionvalue,...) sets CONN_DCM2NII options
%
% e.g. conn_dcmconvert /data/MR*.dcm folderout /data/nii
%
% see also CONN_DCMDIR, CONN_DCM2NII
%

% regexp ImageType DICOM fields to identify anat/func/etc.
ISANAT='PRIMARY.OTHER|ORIGINAL.PRIMARY.M.ND.NORM';
ISFUNC='ORIGINAL.PRIMARY.M.ND.MOSAIC';
ISDWI='ORIGINAL.PRIMARY.DIFFUSION.NONE.ND.MOSAIC';
ISFMAP='ORIGINAL.PRIMARY.(M|P).ND[^\w]*$';

outstruct.filename={};
outstruct.filedescrip={};
outstruct.filetype=[];
outstruct.isanat=[];
outstruct.isfunc=[];
outstruct.isdwi=[];
outstruct.isfmap=[];
outstruct.isselected=[];
DOGUI=false; %usejava('awt');
opts={'overwrite',true};
if nargin<1||isempty(filenames), 
    DOGUI=true;
    if 1
        conn_disp('Select root folder of DICOM data');
        pathname = uigetdir(pwd, 'Select root folder of DICOM data');
        conn_disp('Select target folder');
        folderout = uigetdir(pathname, 'Select target folder');
        if isequal(folderout,0), return; end
        filenames=conn_dir(fullfile(pathname,'*'),'-cell','-sort');
    else
        conn_disp('Select all DICOM files');
        [filename, pathname] = uigetfile( {'*.*', 'All Files (*.*)';'*.dcm', 'DCM Files (*.dcm)'}, 'Select all DICOM files','MultiSelect', 'on');
        if isequal(filename,0), return; end
        filename=conn_sortfilenames(cellstr(filename));
        filenames=cellfun(@(x)fullfile(pathname,x),filename,'uni',0);
        conn_disp('Select target folder');
        folderout = uigetdir(pathname, 'Select target folder');
        if isequal(folderout,0), return; end
    end
elseif isequal(filenames,'-init')
    return
else
    folderout='./nii';
end
opts=[{'folderout',folderout},opts];
if ischar(filenames), filenames=conn_dir(filenames,'-cell','-s'); end
if ischar(filenames), filenames=cellstr(filenames); end

if isstruct(filenames), Series=filenames;
else Series=conn_dcmdir(filenames);
end
[filenameout,filetypeout,fileRTout,filedescripout,Series]=conn_dcm2nii(Series,opts{:},'renamefiles',true,varargin{:});
valid=filetypeout>0;
i=find(valid);
try, 
    filesout_path=fileparts(filenameout{i(1)});
    logfilename=fullfile(filesout_path,sprintf('conn_dcmconvert_%s.log',datestr(now,'yyyy_mm_dd_HHMMSSFFF')));
    fhlog=fopen(logfilename,'wt'); 
    for n1=1:numel(i), fprintf(fhlog,'%s\n',filedescripout{i(n1)}); end
    fclose(fhlog);
    save(conn_prepend('',logfilename,'.mat'),'filenameout','filetypeout','fileRTout','filedescripout');
    conn_disp('fprintf','DCM2NII log information stored in %s\n',logfilename);
    %note: ImageType=regexp(filedescripout,'(ORIGINAL|DERIVED)\S+','match','once')
    %note: SeriesDescription=regexprep(filedescripout,'\(\d+x.*|^Series #\d+ \: ','')
end

isanat=find(filetypeout==1&arrayfun(@(n)~isempty(regexp(Series(n).SeriesInfo.ImageType,ISANAT)),1:numel(Series),'ErrorHandler',@(varargin)false));
isfunc=find(filetypeout>1&arrayfun(@(n)~isempty(regexp(Series(n).SeriesInfo.ImageType,ISFUNC)),1:numel(Series),'ErrorHandler',@(varargin)false));
isdwi=find(arrayfun(@(n)~isempty(regexp(Series(n).SeriesInfo.ImageType,ISDWI)),1:numel(Series),'ErrorHandler',@(varargin)false));
isfmap=find(arrayfun(@(n)~isempty(regexp(Series(n).SeriesInfo.ImageType,ISFMAP)),1:numel(Series),'ErrorHandler',@(varargin)false));

ORG={'FUNCTIONAL','func',isfunc; 'STRUCTURAL','anat',isanat; 'FIELDMAP','fmap',isfmap; 'DTI','dwi',isdwi};
valid=find(filetypeout>0);
for n=1:size(ORG,1)
    if isempty(valid), break; end
    i=find(ismember(valid,ORG{n,3}));
    j=conn_dcmselect(filenameout(valid),filedescripout(valid),ORG{n,2},ORG{n,1},i,DOGUI);
    valid(j)=[];
end

outstruct.filename=filenameout;
outstruct.filedescrip=filedescripout;
outstruct.filetype=filetypeout;
outstruct.isanat=isanat;
outstruct.isfunc=isfunc;
outstruct.isdwi=isdwi;
outstruct.isfmap=isfmap;
outstruct.isselected=isfunc;

