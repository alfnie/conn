function conn_imarrange(filesin,pos,fileout,crop)
% IMARRANGE arrange multiple images in matrix configuration
% conn_imarrange(filesin, pos, fileout)
%   filesin : string or cell array with input-image file names
%   pos     : two-dimensional matrix of indexes to filesin
%   fileout : name of output image file
%   crop    : (optional) 1/0 crop white-space borders from each image

if nargin<4||isempty(crop), crop=false; end
if ischar(filesin), filesin=cellstr(filesin); end
if isempty(pos), pos=reshape(1:numel(filesin),size(filesin)); end
for n1=1:numel(filesin),
    [a,b]=ind2sub(size(pos),find(pos==n1));
    A{a,b}=imread(deblank(filesin{n1}));
    if crop
        for nr=1:2,
            ok=find(any(min(A{a,b},[],1)~=max(A{a,b},[],1),3));
            if numel(ok)<2, ok=1:size(A{a,b},2); end
            A{a,b}=A{a,b}(:,ok(1):ok(end),:);
            A{a,b}=permute(A{a,b},[2,1,3]);
        end
    end
    sX=size(A{a,b});
    bg=A{a,b}(1);
end
for n1=1:prod(size(A)), 
    if isempty(A{n1}), A{n1}=repmat(bg,zeros(sX)); end 
    if n1>1&&size(A{n1},1)>size(A{1},1), A{n1}=A{n1}(1:size(A{1},1),:,:); end
    if n1>1&&size(A{n1},1)<size(A{1},1), A{n1}=cat(1,A{n1},repmat(A{n1}(end,:,:),size(A{1},1)-size(A{n1},1),1)); end
    if n1>1&&size(A{n1},2)>size(A{1},2), A{n1}=A{n1}(:,1:size(A{1},2),:); end
    if n1>1&&size(A{n1},2)<size(A{1},2), A{n1}=cat(2,A{n1},repmat(A{n1}(:,end,:),1,size(A{1},2)-size(A{n1},2))); end
end
A=cell2mat(A);
imwrite(A,fileout);

