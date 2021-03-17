function filename=conn_roiclusters(rois,xy2,xyz2,groups,filename)
% CONN_ROICLUSTERS creates ROIorder.mat file defining list of ROIs as well as ROI order&clusters 
% that can be used by conn_display in ROI-to-ROI analyses
%
% fileout = conn_roiclusters(rois [,coords2, coords3, groups, fileout])
%   rois        : list of ROI names; (1xN) cell array describing N groups of ROIs (each group with an arbitrary number of ROIs) listed as rois{ngroup}{nroi}
%   coords2     : (optional) wheel coordinates for each ROI; (1xM) cell array or (2xM) matrix (default: distributed in unit circle based on group membership)
%   coords3     : (optional) spatial coordinates for each ROI; (1xM) cell array or (3xM) matrix (default: {})
%                 alternatively, filename of atlas file defining these ROIs (coordinates will be extracted using conn_roicenters)
%   groups      : (optional) list of group names; (1xN) cell array (default {})
%   fileout     : (optional) output filename (default 'ROIorder.mat')
%
% e.g. conn_roiclusters( {{'roi1','roi2'},{'roi3','roi4','roi5'}}, {}, {}, {'group1','group2'} );
%

global CONN_gui;

if nargin<2||isempty(xy2), xy2={}; end
if nargin<3||isempty(xyz2), xyz2={}; end
if nargin<4||isempty(groups),groups={}; end
if nargin<5||isempty(filename), filename=fullfile(pwd,'ROIorder.mat'); end

if isstruct(rois)  % accepts alternative format ROIconfiguration struct
    roinames=rois.names2(rois.displaytheserois);            % roi_names: list of valid roi-names (other names in roi_groups will be interpreted as group-names)
    groups={};
    lastcluster=nan;
    numcluster=1;
    for n1=1:numel(rois.displaytheserois)
        if rois.clusters(rois.displaytheserois(n1))~=lastcluster
            if isfield(rois,'names_clusters')&&numel(rois.names_clusters)>=rois.clusters(rois.displaytheserois(n1)), groups{end+1}=rois.names_clusters{rois.clusters(rois.displaytheserois(n1))}; 
            else groups{end+1}=sprintf('G%d',numcluster);
            end
            numcluster=numcluster+1;
        end
        groups{end+1}=rois.names2{rois.displaytheserois(n1)};
        lastcluster=rois.clusters(rois.displaytheserois(n1));
    end
    if 1
        if isempty(CONN_gui)||~isfield(CONN_gui,'font_offset'), conn_font_init; end
        thfig=dialog('units','norm','position',[.2,.3,.6,.5],'windowstyle','normal','name','ROI order/groups','color','w','resize','on');
        uicontrol(thfig,'style','text','units','norm','position',[.1,.72,.8,.25],'string',{'Manually resort the following ROIs in the desired order (e.g. use copy/paste), adding, where appropriate, header lines defining group-names for each separate group of ROIs',' ','(e.g. enter the name of the first group of ROIs in the first line followed by the names of the ROIs within this group, each on a separate line, then enter the name of the second group again followed by the names of the ROIs within the second group, etc.)'},'backgroundcolor','w','fontsize',8+CONN_gui.font_offset);
        ht1=uicontrol(thfig,'style','edit','units','norm','position',[.1,.15,.8,.55],'max',2,'string',groups,'fontsize',8+CONN_gui.font_offset,'horizontalalignment','left');
        uicontrol(thfig,'style','pushbutton','string','Apply','units','norm','position',[.1,.01,.38,.10],'callback','uiresume','fontsize',8+CONN_gui.font_offset);
        uicontrol(thfig,'style','pushbutton','string','Cancel','units','norm','position',[.51,.01,.38,.10],'callback','delete(gcbf)','fontsize',8+CONN_gui.font_offset);
        ok=true;
        while ok
            uiwait(thfig);
            ok=ishandle(thfig);
            if ok,
                groups=get(ht1,'string');
                ok=false;
                delete(thfig);
            else
                filename='';
                return
            end
        end
    end
    groups=groups(cellfun('length',groups)>0);
    groupname={'grouproi'};
    grouprois={[]};
    fprintf('available ROI labels: %s\n',sprintf('%s ',roinames{:}));
    for n1=1:numel(groups)
        idxc=find(strcmp(groups{n1},roinames));
        if numel(idxc)~=1, idxc=find(strncmp(groups{n1},roinames,numel(groups{n1}))); end
        if numel(idxc)~=1, idxc=find(strncmp(groups{n1},regexprep(roinames,'^.*\.',''),numel(groups{n1}))); end
        if numel(idxc)==1, grouprois{end}=[grouprois{end} idxc];
        else %if isempty(idxc),
            fprintf('ROI-label %s not found. Assuming label of new group\n',groups{n1});
            groupname{end+1}=groups{n1};
            grouprois{end+1}=[];
        %else fprintf('warning: multiple (%d) possible ROI-label matches found for %s in %s. Selecting first\n',numel(idxc),groups{n1},sprintf('%s ',roinames{idxc})); grouprois{end}=[grouprois{end} idxc(1)];
        end
    end
    newroinames=roinames(setdiff(1:numel(roinames),[grouprois{:}]));
    for n1=1:numel(groupname)
        if ~isempty(grouprois{n1})
            newroinames=[newroinames groupname(n1)];
        end
    end
    rois=grouprois(cellfun('length',grouprois)>0);
    rois=cellfun(@(x)roinames(x),rois,'uni',0);
    groups=newroinames;
end
M=numel(rois);
for n=1:M, if ischar(rois{n}), rois{n}=cellstr(rois{n}); end; end
names2=cellfun(@(x)reshape(x,1,[]),rois(:)','uni',0);
names2=cat(2,names2{:});
clusters=arrayfun(@(n)n+zeros(1,numel(rois{n})),1:M,'uni',0);
clusters=cat(2,clusters{:})';
N=numel(names2);
displaytheserois=1:N;
if isempty(xy2), a=-pi+2*pi*((1:N)'+clusters-1)/(N+M); xy2=[cos(a) sin(a)]';
elseif iscell(xy2), xy2=reshape(cat(2,xy2{:}),2,N);
end
if isempty(xyz2),
elseif iscell(xyz2), xyz2=reshape(cat(2,xyz2{:}),3,N);
end

ROIconfiguration=struct('names2',{names2},'xy2',200*xy2.','displaytheserois',displaytheserois,'clusters',clusters);
if ~isempty(groups), ROIconfiguration.names_clusters=groups; end
if ~isempty(xyz2), ROIconfiguration.xyz2=xyz2.'; end
assert(size(xy2,2)==N,'incorrect number of elements in coords2 input (found %d expected %d)',size(xy2,2),N); 
assert(isempty(xyz2)|size(xyz2,2)==N,'incorrect number of elements in coords3 input (found %d expected %d)',size(xyz2,2),N); 
assert(isempty(groups)|numel(groups)==M,'incorrect number of elements in groups input (found %d expected %d)',numel(groups),M); 
conn_savematfile(filename,'ROIconfiguration');


