function conn_ROI2dl(filename,rois,thr)
% conn_ROI2dl creates .dl UCINET format adjacency matrix information from
% first-level ROI results. 
%

if nargin<1||isempty(filename),
    global CONN_x;
    if ~isempty(CONN_x)&&isfield(CONN_x,'folders')&&isfield(CONN_x.folders,'firstlevel'),pathname=fullfile(CONN_x.folders.firstlevel,CONN_x.Analyses(1).name);
    else, pathname=pwd; end
    [filename, pathname] = uigetfile('.mat', 'Select a condition file',fullfile(pathname,'resultsROI_Condition*.mat'));
    filename=fullfile(pathname,filename);
end

% loads connectivity values
load(filename,'Z','names')
[nROI,nROI2,nSubjects]=size(Z);
if nargin<2||isempty(rois),
    [rois,ok]=listdlg('PromptString','Select ROIs:',...
        'SelectionMode','multiple',...
        'ListString',char(names),...
        'initialvalue',1:nROI);
    if ~ok,return;end
end
names={names{rois}};Z=Z(rois,rois,:);nROI=length(rois);
%averages across all subjects
z=nanmean(Z(1:nROI,1:nROI,:),3);

% plots connectivity matrix
figure;subplot(121);imagesc(z);axis equal;axis tight;
title('\fontsize{14}raw connectivity \fontsize{12}(z- values)');
subplot(122);hist(z(:));axis square;
title('\fontsize{14}histogram');
xlabel('\fontsize{12}z- values');

% threshold connectivity values
if nargin<3||isempty(thr),
    answer=inputdlg({'Enter the z-value threshold'},'',1,{num2str(nanmean(z(:)))});
    if isempty(answer),return;end
    thr=str2double(answer{1});
end
z_thr=z>thr;

% plots thresholded connectivity matrix
subplot(121);imagesc(z_thr);axis equal;axis tight;
title('\fontsize{14}adjacency matrix');
subplot(122);hold on;plot(thr*[1,1],get(gca,'ylim'),'k-');hold off;

% creates .dl UCINET format file
[pathname,filename,fileext]=fileparts(filename);
filename_out=fullfile(pathname,[filename,'.dl']);
fh=fopen(filename_out,'wt');
while fh<0, [filename_out,filename_outpath]=uiputfile('*.dl','select an output file name'); filename_out=fullfile(filename_outpath,filename_out); fh=fopen(filename_out,'wt'); end
fprintf(fh,'dl format=fullmatrix,n=%d\n',nROI);
fprintf(fh,'labels:\n');
for n1=1:nROI,name=names{n1};name(name==' ')=[];name=name(1:min(numel(name),18));fprintf(fh,'%s',name);if n1==nROI,fprintf(fh,'\n');else,fprintf(fh,',');end;end
fprintf(fh,'data:\n');
for n1=1:nROI,fprintf(fh,'%d ',z_thr(n1,:));fprintf(fh,'\n');end
fclose(fh);

fprintf(1,'Output .dl file created :\n %s\n',filename_out);

