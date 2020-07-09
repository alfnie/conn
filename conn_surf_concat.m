function filename=conn_surf_concat(filename_lh, filename_rh, filename)
if ~iscell(filename_lh), filename_lh=cellstr(filename_lh); end
if ~iscell(filename_rh), filename_rh=cellstr(filename_rh); end
if numel(filename_lh)~=numel(filename_rh), error('input lengths must match in size'); end
if nargin<3||isempty(filename), filename=cell(size(filename_lh)); end
ps=conn_surf_dims(8);
for n=1:numel(filename_lh)
    a1=spm_vol(filename_lh{n});
    a2=spm_vol(filename_rh{n});
    if prod(a1(1).dim(1:3))~=prod(ps), error('Incorrect volume dimensions in %s (prod(%s)~=%d)',filename_lh{n},num2str(a1(1).dim(1:3)), prod(ps)); end
    if prod(a2(1).dim(1:3))~=prod(ps), error('Incorrect volume dimensions in %s (prod(%s)~=%d)',filename_rh{n},num2str(a2(1).dim(1:3)), prod(ps)); end
    if isempty(filename{n})
        [tfilepath,tfilename1,tfileext]=spm_fileparts(filename_lh{n});
        [tfilepath,tfilename2,tfileext]=spm_fileparts(filename_lh{n});
        i1=regexp(tfilename1,'lh');
        i2=regexp(tfilename1,'lh');
        i=intersect(i1,i2);
        if isempty(i)
            tfilename=fullfile(tfilepath,[tfilename1, '_', tfilename2, tfileext]);
        else
            i=i(1);
            tfilename=fullfile(tfilepath,[tfilename2(1:i-1),'lhrh',tfilename2(i+2:end), tfileext]);
        end
        filename{n}=tfilename;
    end
    b1=spm_read_vols(a1);
    b2=spm_read_vols(a2);
    V=a1(1);
    V.fname=filename{n};
    V.pinfo=[1;0;0];
    V.dim=conn_surf_dims(8).*[1 1 2];
    V.mat=eye(4);
    b=reshape(cat(3,b1,b2), [V.dim size(b1,4)]);
    if size(b,4)>1
        V=repmat(V,[size(b,4),1]);for nh=1:size(b,4),V(nh).n=[nh,1];end
        V=spm_create_vol(V);
        for nh=1:size(b,4), V(nh)=spm_write_vol(V(nh),b(:,:,:,nh)); end
    else
        spm_write_vol(V,b);
    end
end

