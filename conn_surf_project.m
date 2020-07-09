function [mxyz,idx]=conn_surf_project(xyz,filenames,matched,file_patch,file_offset,file_scale)
% internal function

if ~iscell(filenames), filenames={filenames}; end
for n1=1:numel(filenames)
    if nargin>=4&&numel(file_patch)>=n1
        xyz0{n1}=file_patch{n1}.vertices;
        if nargin>=6
            mp=mean(xyz0{n1},1);
            xyz0{n1}=conn_bsxfun(@plus,mp+file_offset{n1},conn_bsxfun(@minus,xyz0{n1},mp)*file_scale{n1});
        end
    else
        xyz0{n1}=conn_freesurfer_read_surf(filenames{n1});
    end
end
if matched
    filenamesmatched=cell(1,numel(filenames));
    for n1=1:numel(filenames), 
        [tfile_path,tfile_name,tfile_ext]=fileparts(filenames{n1}); 
        tfilename=regexprep(tfile_name,'^(lh\.|rh\.).+','$1');
        if any(strcmp(tfilename,{'lh.','rh.'})), filenamesmatched{n1}=fullfile(fileparts(which(mfilename)),'surf',[tfilename,'pial.surf']); end
    end
    if any(~cellfun('length',filenamesmatched)), error('Unable to identify hemisphere from filename (file names should start with lh. or rh.)'); end
    [filenames1,matchedfilesall]=unique(filenamesmatched);
    xyz1=xyz0;
    xyz0=cell(1,numel(filenames1));
    for n1=1:numel(filenames1)
        xyz0{n1}=conn_freesurfer_read_surf(filenames1{n1});
    end
else
end
d=zeros(size(xyz,1),numel(xyz0));
idx=d;
for n1=1:numel(xyz0),
    [d(:,n1),idx(:,n1)]=min(permute(sum(abs(conn_bsxfun(@minus,xyz0{n1},permute(xyz,[3,2,1]))).^2,2),[3,2,1]),[],3);
end
[mind,minidx]=min(d,[],2);
mind=sqrt(mind);

mxyz=xyz;
for n1=1:size(d,1),
    if matched,
        matchedfile=matchedfilesall(minidx(n1));
        matchedvertex=idx(n1,matchedfile);
        mxyz(n1,:)=xyz1{matchedfile}(matchedvertex,:);
        if mind(n1)>10, fprintf('Warning! Surface coordinates are far from original coordinates. '); end
        fprintf('Matched coordinates (%d,%d,%d) to surface %s at (%d,%d,%d) paired with surface %s at (%d,%d,%d)\n',round(xyz(n1,1)),round(xyz(n1,2)),round(xyz(n1,3)),filenamesmatched{minidx(n1)},round(xyz0{minidx(n1)}(matchedvertex,1)),round(xyz0{minidx(n1)}(matchedvertex,2)),round(xyz0{minidx(n1)}(matchedvertex,3)),filenames{matchedfile},round(mxyz(n1,1)),round(mxyz(n1,2)),round(mxyz(n1,3)));
    else
        matchedfile=minidx(n1);
        matchedvertex=idx(n1,matchedfile);
        mxyz(n1,:)=xyz0{matchedfile}(matchedvertex,:);
        if mind(n1)>10, fprintf('Warning! Surface coordinates are far from original coordinates. '); end
        fprintf('Matched coordinates (%d,%d,%d) to surface %s at (%d,%d,%d)\n',round(xyz(n1,1)),round(xyz(n1,2)),round(xyz(n1,3)),filenames{matchedfile},round(mxyz(n1,1)),round(mxyz(n1,2)),round(mxyz(n1,3)));
    end
end
if matched
    idx=matchedfilesall(minidx);
else
    idx=minidx;
end
