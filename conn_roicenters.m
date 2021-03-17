function [ROIxyz,ROInames,ROIidx,ROIstruct]=conn_roicenters(filename,varargin)
% CONN_ROICENTERS returns voxel/vertex within ROI which is closest to center-of-mass of ROI
% [ROIxyz,ROInames,ROIidx,ROIstruct] = conn_roicenters(filename)
%   filename : ROI file
%   ROIxyz   : 3xN ROI center coordinates
%   ROInames : N ROI names
%   ROIidx   : N ROI numeric indexes
%   ROIstruct : ROI struct - e.g. for conn_mesh_display ROI plots; conn_mesh_display(filename,[],[],ROIstruct)
% 
% note: for surfaces enter conn_roicenters(filename,ctype) with ctype='volumes'|'sphere'|'inflated' (default inflated)
%       

if any(conn_server('util_isremotefile',filename)), [ROIxyz,ROInames,ROIidx,ROIstruct]=conn_server('run',mfilename,conn_server('util_localfile',filename),varargin{:}); return; end
ROIxyz=[];
ROIdata=[];
ROIstruct=[];
[ROInames,ROIidx]=conn_roilabels(filename);
if isempty(ROIidx), fprintf('warning: unable to interpret roi file %s\n',filename);
else
    tfilename=cellstr(conn_expandframe(filename));
    a=spm_vol(filename);
    b=spm_read_vols(a);
    if numel(ROInames)>1&&numel(tfilename)==1, %3d-atlas
        maxdata=max(b(:));
        if max(ROIidx)==maxdata
            ROIdata=arrayfun(@(n)b==n,ROIidx,'uni',0);
        end
    elseif numel(ROInames)>1 && numel(ROInames)==numel(tfilename), %4d-atlas
        ROIdata=num2cell(b(:,:,:,ROIidx)>0,4);
    else ROIdata={b>0};
    end
    if isempty(ROIdata), fprintf('warning: unable to interpret data in roi file %s\n',filename);
    else
        if conn_surf_dimscheck(a(1))
            xyz=conn_surf_coords([],varargin{:}); % note: for surfaces enter conn_roicenters(filename,ctype) with ctype='volumes'|'sphere'|'inflated' (default inflated)
        else
            [x,y,z]=ndgrid(1:a(1).dim(1),1:a(1).dim(2),1:a(1).dim(3));
            xyz=a(1).mat*[x(:) y(:) z(:) ones(numel(x),1)]';
        end
        ROIxyz=nan(3,numel(ROIdata));
        for n=1:numel(ROIdata)
            idx=find(ROIdata{n});
            if ~isempty(idx)
                txyz=xyz(:,idx);
                mxyz=mean(txyz,2);
                [nill,i]=min(sum(abs(txyz-repmat(mxyz,1,size(txyz,2))).^2,1));
                ROIxyz(:,n)=txyz(1:3,i); % return point within ROI closest to ROI-centroid
            end
        end
        ROIstruct=struct('sph_xyz',ROIxyz','sph_names',{ROInames},'sph_r',ones(size(ROInames)));
    end
end