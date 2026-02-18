function data=conn_roiclusters_load(filename,roinames,DOPLOT)
% conn_roiclusters_load
% loads ROI order&clusters from ROIorder.mat file
%
% data=conn_roiclusters_load(filename); % loads ROI order&clusters information
% data=conn_roiclusters_load(filename,roinames); % only loads information from subset of ROIs matching roinames
%
% data.roinames             : names of ROIs
% data.displaytheserois     : order of ROIs for display
% data.clusters             : cluster membership of each ROI
% data.names_clusters       : cluster names
% data.xy2                  : ring coordinates of ROIs
% data.xy2_clusters         : ring coordinates of clusters
% data.xyz2                 : 3D coordinates of ROIs
%
if nargin<3||isempty(DOPLOT), DOPLOT=true; end

conn_loadmatfile(filename,'ROIconfiguration');
extnames2=ROIconfiguration.names2(ROIconfiguration.displaytheserois);
if isempty(roinames), roinames=extnames2; end
data.names=roinames;
ok=ismember(extnames2,roinames);
lnames=lower(roinames);
for n1=reshape(find(~ok),1,[])
    if ismember(lower(extnames2{n1}),lnames), extnames2{n1}=roinames{idx};
    else
        idx=strmatch(extnames2{n1},roinames); % allows partial-name matches
        if numel(idx)==1, extnames2{n1}=roinames{idx}; end
    end
end
[ok,idx]=ismember(roinames,extnames2);
if ~nnz(ok),
    if DOPLOT, conn_msgbox('Unable to import ROI configuration information. No matching ROIs','',2); end
    return
end
data.displaytheserois=find(ok);
if ~isfield(data,'xy2'), data.xy2=nan(numel(roinames),2); end
data.xy2(:)=nan;
if isfield(ROIconfiguration,'xy2'),
    data.xy2(data.displaytheserois,:)=ROIconfiguration.xy2(ROIconfiguration.displaytheserois(idx(ok)),:);
else
    [idxok,tidx]=sort(idx(ok));
    clusters=1+cumsum(diff(ROIconfiguration.clusters(ROIconfiguration.displaytheserois(idxok)))~=0);
    a=-pi+2*pi*((1:numel(data.displaytheserois))'+[1;clusters(:)]-1)/(numel(data.displaytheserois)+max(clusters));
    data.xy2(data.displaytheserois(tidx),:)=200*[cos(a) sin(a)];
    %             [clusters,tidx]=sort(ROIconfiguration.clusters(ROIconfiguration.displaytheserois(idx(ok))));
    %             a=-pi+2*pi*((1:numel(data.displaytheserois))'+clusters-1)/(numel(data.displaytheserois)+max(clusters));
    %             data.xy2(data.displaytheserois(tidx),:)=200*[cos(a) sin(a)];
end
if ~isfield(data,'clusters'), data.clusters=nan(numel(roinames),1); end
data.clusters(:)=0;
data.clusters(data.displaytheserois)=ROIconfiguration.clusters(ROIconfiguration.displaytheserois(idx(ok)));
if isfield(ROIconfiguration,'xy2_clusters'),
    data.xy2_clusters=data.xy2;
    data.xy2_clusters(data.displaytheserois,:)=ROIconfiguration.xy2_clusters(ROIconfiguration.displaytheserois(idx(ok)),:);
    mxy2_clusters=[]; for n1=1:max(data.clusters), if any(data.clusters==n1), mxy2_clusters(n1,:)=mean(data.xy2_clusters(data.clusters==n1,:),1); end; end
    data.xy2_clusters(data.clusters>0,:)=mxy2_clusters(data.clusters(data.clusters>0),:);
elseif ~isempty(data.clusters) % note: automatic fill-in if missing
    data.xy2_clusters=data.xy2;
    mxy2_clusters=[]; for n1=1:max(data.clusters), if any(data.clusters==n1), mxy2_clusters(n1,:)=mean(data.xy2(data.clusters==n1,:),1); end; end
    data.xy2_clusters(data.clusters>0,:)=mxy2_clusters(data.clusters(data.clusters>0),:);
else data.xy2_clusters=[];
end
if isfield(ROIconfiguration,'names_clusters')&&~isempty(ROIconfiguration.names_clusters), data.names_clusters=ROIconfiguration.names_clusters;
else data.names_clusters={};
end
if ~isfield(data,'xyz2'), data.xyz2=nan(numel(roinames),3); end
if isfield(ROIconfiguration,'xyz2')&&~isempty(ROIconfiguration.xyz2), % note: automatic fill-in if missing
    data.xyz2(:)=nan;
    data.xyz2(data.displaytheserois,:)=ROIconfiguration.xyz2(ROIconfiguration.displaytheserois(idx(ok)),:);
end
%data.displaytheserois=ROIconfiguration.displaytheserois;
%data.xy2=ROIconfiguration.xy2;
%data.clusters=ROIconfiguration.clusters;
if isfield(ROIconfiguration,'x2'), % order of ROIs (takes precedence over order listed in displaytheserois)
    [nill,tidx]=sort(ROIconfiguration.x2(ROIconfiguration.displaytheserois(idx(ok))));
else
    [nill,tidx]=sort(mod(pi+angle(data.xy2(data.displaytheserois,:)*[1;1i]),2*pi));
end
data.displaytheserois=data.displaytheserois(tidx);