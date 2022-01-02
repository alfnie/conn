function [filename,Nrois]=conn_createmniroi(fname,xyz,rad,res,mtype,ftype);
% CONN_CREATEMNIROI creates nifti ROI file/atlas containing spherical ROIs center
% at user-defined coordinates
% 
% conn_createmniroi(fname,xyz,rad,res,mtype)
%    fname : output filename                            : char array
%    xyz   : MNI coordinates (in mm)                    : Nx3 matrix with x/y/z coordinates (N = number of ROIs)
%    rad   : sphere radius (in mm)                      : single value or Nx1 vector
%    res   : output file voxel resolution (in mm)       : single value (default 2mm)
%    mtype : type of atlas file                         : '3d' or '4d' (default '3d')

if nargin<4||isempty(res), res=2; end        % voxel resolution (mm)
if nargin<5||isempty(mtype), mtype='3d'; end % 3d, 4d
if nargin<6||isempty(ftype), ftype=0; end    % 1: default bounding box; 2: tight bounding box
if any(conn_server('util_isremotefile',fname)), 
    [filename,Nrois]=conn_server('run',mfilename,conn_server('util_localfile',fname),xyz,rad,res,mtype,ftype); 
    filename=conn_server('util_remotefile',filename);
    return
else fname=conn_server('util_localfile',fname);
end

if ftype==0, ftype=1+(res<1); end
if ischar(xyz), 
    try,
        xyz=spm_load(xyz);
    catch
        try
            xyz=regexp(fileread(xyz),'[\r\n]','split');
            xyz=regexp(xyz(cellfun('length',xyz)>0),'[-+]?\d*\.?\d*','match');
            xyz=str2double(cat(1,xyz{:}));
        catch
            error('unable to interpret file %s. Each row should contain only three numbers separated by spaces or commas',xyz);
        end
    end
end
Nrois=size(xyz,1);
if isempty(rad)&&size(xyz,2)==4, rad=xyz(:,4); xyz=xyz(:,1:3); end
if size(xyz,2)~=3, error('incorrect xyz dimensions (expected matrix with 3 columns)'); end
if Nrois>1||numel(rad)>1, ftype=1; end
if numel(rad)>1&&numel(rad)~=size(xyz,1), error('incorrect rad dimensions (expected single element, or vector with %d elemens)',size(xyz,1)); end
if ischar(mtype),
    switch(mtype)
        case '3d', mtype=1;
        case '4d', mtype=2;
        otherwise, error('unrecognized file type %s (valid types are ''3d'' or ''4d'' for 3d/4d nifti atlas files)',mtype);
    end
elseif mtype==0, mtype=1+(Nrois>1); 
end
rad3=rad;
if size(rad3,1)==1, rad3=rad3'; end
if size(rad3,1)==1, rad3=repmat(rad3,size(xyz,1),1); end
if size(rad3,2)==1, rad3=repmat(rad3,1,3); end
switch(ftype)
    case 1, % default bounding box
        bbox=[-90,-126,-72;90,90,108];
        bbox=[min(bbox(1,:),min(xyz-rad3,[],1)); max(bbox(2,:),max(xyz+rad3,[],1))];
        [x,y,z]=ndgrid(bbox(1,1):res:bbox(2,1), bbox(1,2):res:bbox(2,2), bbox(1,3):res:bbox(2,3));
        if mtype==1, b=zeros(size(x));
        else         b=zeros([size(x,1) size(x,2) size(x,3) Nrois]);
        end
        mind2=inf(size(x));
        for n=1:Nrois, 
            td2=(x-xyz(n,1)).^2+(y-xyz(n,2)).^2+(z-xyz(n,3)).^2;
            mask=(td2<rad(min(numel(rad),n))^2);
            if mtype==1, tmask=td2<mind2 & mask; b(tmask)=n; mind2(tmask)=td2(tmask); % note: if overlapping, voxel assigned to closest ROI centroid
            else         b(find(mask)+size(b,1)*size(b,2)*size(b,3)*(n-1))=1;
            end
        end
        mat=[diag([1 1 1])*res, bbox(1,:)'-res*[1;1;1]; zeros(1,3) 1];
    case 2, % tight bounding box
        n=ceil(rad/res);
        [x,y,z]=ndgrid((-n:n)*res);
        b=x.^2+y.^2+z.^2<rad^2;
        mat=[diag([-1 1 1])*res, xyz(:)-((n+1)*res)*[-1;1;1]; zeros(1,3) 1];
    otherwise
        error('unknown type %d (1:default bounding box; 2:tight bounding box)',ftype);
end
filename=fname;
vol=struct('fname',filename, 'mat',mat, 'dim',[size(b,1),size(b,2),size(b,3)], 'dt', [spm_type('uint16') spm_platform('bigend')], 'pinfo',[1;0;0],'descrip',''); 
if Nrois==1, vol.descrip=sprintf('Seed %s mm radius %d',mat2str(xyz),rad); end
try, spm_unlink(filename); end
if mtype==1, vol=spm_write_vol(vol,b);
else vol=repmat(vol,[Nrois,1]);
    for nh=1:Nrois,vol(nh).n=[nh,1];end
	vol=spm_create_vol(vol);
    for nh=1:Nrois,spm_write_vol(vol(nh),b(:,:,:,nh)); end
end
if Nrois>1
    %fh=fopen(conn_prepend('',filename,'.txt'),'wt');
    fh={};
    for n=1:Nrois, fh{end+1}=sprintf('(%d,%d,%d)\n',round(xyz(n,1)),round(xyz(n,2)),round(xyz(n,3))); end
    %fclose(fh);
    conn_fileutils('filewrite_raw',conn_prepend('',filename,'.txt'), fh);
end
end
