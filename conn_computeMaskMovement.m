function fileout = conn_computeMaskMovement(filein)
if ~ischar(filein), filein=char(filein); end
[filename_path,filename_name,filename_ext]=spm_fileparts(filein(1,:));
a=spm_vol(filein);
[tempx,tempy,tempz]=ind2sub(a(1).dim(1:3),1:prod(a(1).dim(1:3)));
xyz_voxel=[tempx(:),tempy(:),tempz(:),ones(numel(tempx),1)]';
xyz=a(1).mat*xyz_voxel;

MaskGradient=zeros(prod(a(1).dim(1:3)),3);
for n=1:numel(a)
    b=reshape(spm_get_data(a(n),pinv(a(n).mat)*xyz),a(1).dim(1:3));
    mask=b>mean(b(~isnan(b)&b~=0))/8;
    mask=b>0.80*mean(b(mask))/4;
    b(~mask)=nan;
    [dy,dx,dz]=gradient(b);
    dX=[dx(:) dy(:) dz(:)]*a(1).mat(1:3,1:3)';
    ok=~isnan(dX);
    MaskGradient(ok)=MaskGradient(ok)+dX(ok);
    if ~rem(n,10), fprintf('.'); end
end
if n>=10,fprintf('Done\n');end
MaskGradient=MaskGradient/numel(a);

% write output files
DIMS=6;
MaskGradient=[MaskGradient, ...
    MaskGradient(:,2).*xyz(3,:)'-MaskGradient(:,3).*xyz(2,:)',...
    MaskGradient(:,1).*xyz(3,:)'-MaskGradient(:,3).*xyz(1,:)',...
    MaskGradient(:,1).*xyz(2,:)'-MaskGradient(:,2).*xyz(1,:)'];
MaskGradient=MaskGradient./repmat(max(abs(MaskGradient),[],1),size(MaskGradient,1),1);
fileout=fullfile(filename_path,['MaskMotion_',filename_name,filename_ext]);
V=struct('mat',a(1).mat,'dim',a(1).dim(1:3),'pinfo',[1;0;0],'dt',[spm_type('float32') spm_platform('bigend')],'fname',fileout);
V=repmat(V,[DIMS,1]);for nh=1:DIMS,V(nh).n=[nh,1];end
V=spm_create_vol(V);
for nh=1:DIMS, V(nh)=spm_write_vol(V(nh),reshape(MaskGradient(:,nh),a(1).dim(1:3))); end

