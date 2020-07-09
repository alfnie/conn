function V=conn_create_vol(filename,V1,voxels,V0,sfilename,mask,doPSC,issurface,surfacesmooth)

% CONN_CREATE_VOL creates CONN volume data from SPM volume and list of voxel indices
% conn_create_vol(fileoutput,spmVsource,voxels,spmVreference,sfilename_reference_to_source,mask);
%

if ischar(V1), V1=spm_vol(V1); end
if nargin<4 || isempty(V0), V0=V1; end
if nargin<5, sfilename=[]; end
if nargin<6, mask=[]; end
if nargin<7||isempty(doPSC), doPSC=1; end
if nargin<8||isempty(issurface), issurface=0; end
if nargin<9||isempty(surfacesmooth), surfacesmooth=0; end
if ischar(V0), V0=spm_vol(V0); end
if ischar(mask), mask=spm_vol(mask); end
%if ~isempty(mask), maskX=spm_read_vols(mask); end

issourcesurface=conn_surf_dimscheck(V1);
issurface=issurface|issourcesurface;
Nt=numel(V1);
V.fname=filename;
V.issurface=issurface;
if issurface, V.surfacesmooth=surfacesmooth; end
if issurface
    dim=conn_surf_dims(8).*[1 1 2];
    V.matdim=struct('mat',eye(4),'dim',dim);
    V.size=struct('Nt',Nt,'Ns',dim(3),'Nv',zeros(1,dim(3)));
else
    V.matdim=struct('mat',V0(1).mat,'dim',V0(1).dim);
    V.size=struct('Nt',Nt,'Ns',V0(1).dim(3),'Nv',zeros(1,V0(1).dim(3)));
end
if isempty(voxels), voxels=1:prod(V.matdim.dim(1:3)); 
else voxels=sort(voxels); end
voxels=voxels(:);
V.voxels=[];
nslice=1+floor((voxels-1)/prod(V.matdim.dim(1:2)));
Nv=hist(nslice,1:V.matdim.dim(3));

handle=fopen([V.fname,'c'],'wb');
fclose(handle);

V.GM=0;
N=0;
V.gm=zeros(Nt,1);
refinfo=[];
if issurface % surface-based extraction
    idx=1:numel(voxels);
    xyz=conn_convertcoordinates('idx2tal',voxels(idx),V.matdim.mat,V.matdim.dim);
    if ~isempty(sfilename),
        xyzs=conn_obj_coords2(sfilename, xyz')';
        if size(xyzs,2)<4, xyzs=cat(2,xyzs,ones(size(xyzs,1),1)); end
    else xyzs=xyz;
    end
    if ~isempty(mask)
        if ~isequal(V.matdim.dim,mask(1).dim)||~isequal(V.matdim.mat,mask(1).mat), % invalidates masks not in surface space
            mask=spm_vol(fullfile(fileparts(which(mfilename)),'utils','surf','mask.surface.brainmask.nii'));
        end
        z=spm_get_data(mask,pinv(mask(1).mat)*xyz');
        %[z,nill,nill,refinfo]=conn_surf_extract(mask,[],V0,0,false,false,refinfo);
        idx0=find(z>.5);
    else idx0=1:length(idx);
    end
    nslice=1+floor((voxels(idx(idx0))-1)/prod(V.matdim.dim(1:2)));
    V.size.Nv=hist(nslice,1:V.matdim.dim(3));
    N=numel(idx0);
    V.voxels=voxels(idx(idx0));
%     try   % faster / more memory
    x=zeros(numel(idx0),numel(V1));
    for nscan=1:numel(V1),
        if issourcesurface
            xt=spm_get_data(V1(nscan),pinv(V1(nscan).mat)*xyzs(idx0,:)');
        else
            conn_disp('fprintf','Warning: volume-level functional data found. Extracting functional data at surface for surface-level analyses. Please apply the "functional surface resampling" (and optionally "functional surface smoothing") preprocessing steps to your functional data to avoid seeing this warning in the future');
            [data,nill,nill,refinfo]=conn_surf_extract(V1(nscan),[],V0,surfacesmooth,false,false,refinfo);
            xt=data(idx0);
        end
        xt(isnan(xt))=0;
        x(:,nscan)=xt;
        V.gm(nscan)=sum(xt);
        V.GM=V.GM+V.gm(nscan)/length(V1);
    end
%         mx=mean(x,2);
%         for n1=1:size(x,1), x(n1,:)=x(n1,:)-mx(n1); end
    conn_write_vol(V,x');
%     catch % slower / less memory
%         for nscan=1:numel(V1),
%             [data,nill,nill,refinfo]=conn_surf_extract(V1(nscan),[],V0,surfacesmooth,false,false,refinfo);
%             x=data(idx0);
%             x(isnan(x))=0;
%             V.gm(nscan)=sum(x);
%             V.GM=V.GM+V.gm(nscan)/length(V1);
%             conn_write_time(V,x,nscan)
%         end
%     end
elseif issourcesurface, error('Mismatch analysis-space: surface-level functional data (%s) cannot be used for volume-level analyses. Please change analysis-space (e.g. in Setup.Options tab) to "surface", or change functional data to 3d/4d nifti volumes',V1(1).fname);
else % volume-based extraction
    for slice=1:V.matdim.dim(3),
        idx=sum(Nv(1:slice-1))+(1:Nv(slice));
        if ~isempty(idx),
            xyz=conn_convertcoordinates('idx2tal',voxels(idx),V.matdim.mat,V.matdim.dim);
            if ~isempty(sfilename),
                xyzs=conn_obj_coords2(sfilename, xyz')';
                if size(xyzs,2)<4, xyzs=cat(2,xyzs,ones(size(xyzs,1),1)); end
            else xyzs=xyz;
            end
            
            if ~isempty(mask)
                z=spm_get_data(mask,pinv(mask(1).mat)*xyz');
                idx0=find(z>.5);
            else idx0=1:length(idx);
            end
            if ~isempty(idx0),
                x=zeros(length(V1),length(idx0));
                for nscan=1:length(V1),
                    x(nscan,:)=spm_get_data(V1(nscan),pinv(V1(nscan).mat)*xyzs(idx0,:)');
                end
                mx=mean(x,1);
                V.GM=V.GM+sum(mx(~isnan(mx)));
                N=N+size(x,2);
                V.gm=V.gm+sum(x(:,~isnan(mx)),2);
%                 for n1=1:size(x,2), x(:,n1)=x(:,n1)-mx(n1); end
                V.size.Nv(slice)=length(idx0);
                V.voxels=cat(1,V.voxels,voxels(idx(idx0)));
                conn_write_slice(V,x,slice);
            end
        end
    end
end
    
V.GM=V.GM/N;
V.gm=V.gm/N;
if doPSC, 
    V.scale=100/abs(V.GM);
else
    V.scale=1;
end
V.gm=V.gm*V.scale;
V.voxelsinv=zeros(V.matdim.dim(1:3)); 
V.voxelsinv(V.voxels)=reshape(1:length(V.voxels),size(V.voxels));

save(filename,'V'); 
% if str2num(version('-release'))>=14, save(filename,'-V6','V');
% else, save(filename,'V'); end
