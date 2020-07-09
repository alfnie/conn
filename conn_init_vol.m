function [V]=conn_init_vol(V,voxels,values,softlink,overwritesoftlink)

if nargin<5, softlink=[]; end
if nargin<4, overwritesoftlink=false; end
if nargin<3, values=[]; end
if nargin<2, voxels=[]; end
if ~isempty(voxels),
    nslice=1+floor((voxels-1)/prod(V.matdim.dim(1:2)));
    V.size.Nv=hist(nslice,1:V.matdim.dim(3));
    V.voxels=voxels(:);
    V.voxelsinv=zeros(V.matdim.dim(1:3));
    V.voxelsinv(V.voxels)=reshape(1:length(V.voxels),size(V.voxels));
end
if ~isempty(values)&&numel(values)~=sum(V.size.Nv)*V.size.Nt, error('mismatch dimensions'); end
if ~isempty(softlink), 
    softlink=[softlink,'c'];
    str1=regexp(V.fname,'Subject\d+','match'); if ~isempty(str1), softlink=regexprep(softlink,'Subject\d+',str1{end}); end
    [file_path,file_name,file_ext]=fileparts(V.fname);
    [softlink_path,softlink_name,softlink_ext]=fileparts(softlink);
    softlink=[softlink_name,softlink_ext];
    matcfilename=fullfile(file_path,softlink); 
else matcfilename=[V.fname,'c'];
end
if isempty(softlink)||overwritesoftlink
    handle=fopen(matcfilename,'wb');
    if isempty(values)
        fwrite(handle,[0],'float',4*(sum(V.size.Nv)*V.size.Nt));
    else
        fwrite(handle,values,'float');
    end
    fclose(handle);
end
V.scale=1;
V.GM=[];
V.gm=[];
V.softlink=softlink;
V.overwritesoftlink=overwritesoftlink;
save(V.fname,'V'); 
% if str2num(version('-release'))>=14, save(filename,'-V6','V');
% else, save(filename,'V'); end
