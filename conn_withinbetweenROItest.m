function varargout=conn_withinbetweenROItest(gui)
% statistical tests within- and between- networks of ROIs
%

if nargin<1||isempty(gui), gui=1; end

%% gets options from current state of the GUI
global CONN_x CONN_h;
filepathresults1=fullfile(CONN_x.folders.firstlevel,CONN_x.Analyses(CONN_x.Analysis).name); % first-level data folder
nconditions=CONN_x.Results.xX.nconditions;            % selected conditions
cconditions=CONN_x.Results.xX.cconditions;            % between-conditions contrast
nsources=CONN_x.Results.xX.nsources;                  % selected source ROIs (indexes to names below)
csources=CONN_x.Results.xX.csources;                  % between-sources contrast (not used)
nsubjecteffects=CONN_x.Results.xX.nsubjecteffects;    % selected between-subject effects
csubjecteffects=CONN_x.Results.xX.csubjecteffects;    % selected between-subjects contrast
ntargets=CONN_x.Results.xX.roiselected2;              % selected target ROIs (indexes to names below)
X=cell2mat(cellfun(@(x)[x{1:length(CONN_x.Setup.l2covariates.names)-1}],CONN_x.Setup.l2covariates.values(1:CONN_x.Setup.nsubjects),'uni',0)'); % design matrix
allsources=CONN_h.menus.m_results.roiresults.names;   % names of all source ROIs
alltargets=CONN_h.menus.m_results.roiresults.names2;  % names of all target ROIs
if gui
    [sel,ok]=listdlg('PromptString','Select Set-1 ROIs','SelectionMode','multiple','ListString',allsources,'InitialValue',nsources);
    if ok
        nsources=sel;
        [ok,idx]=ismember(allsources(nsources),alltargets);
        ntargets=setdiff(ntargets,idx);
        [sel,ok]=listdlg('PromptString','Select Set-2 ROIs','SelectionMode','multiple','ListString',alltargets,'InitialValue',ntargets);
        if ok ntargets=sel;
        end
    end
end

%% loads the ROI-to-ROI data matrices
icondition=[];isnewcondition=[];for ncondition=nconditions(:)',[icondition(ncondition),isnewcondition(ncondition)]=conn_conditionnames(CONN_x.Setup.conditions.names{ncondition}); end
if any(isnewcondition), error(['Some conditions have not been processed yet. Re-run previous step']); end
Y=zeros([size(X,1),numel(allsources),numel(alltargets),length(nconditions)]);
for n0=1:length(nconditions),
    filename=fullfile(filepathresults1,['resultsROI_Condition',num2str(icondition(nconditions(n0)),'%03d'),'.mat']);
    load(filename,'Z','names','names2','xyz');
    [ok1,iroi1]=ismember(allsources,names);
    [ok2,iroi2]=ismember(alltargets,names2);
    if ~all(ok1)|~all(ok2), error('Error loading data'); return; end
    Y(:,:,:,n0)=permute(Z(iroi1,iroi2,:),[3,1,2]);
end
validsubjects=find(any(X(:,nsubjecteffects)~=0,2)&~any(isnan(X(:,nsubjecteffects)),2));     % remove subjects not included in the design
X=X(validsubjects,nsubjecteffects);
Y=Y(validsubjects,:,:,:);

%% breaks down correlation matrix into within/between matrices
[ok1,idx]=ismember(allsources(nsources),alltargets);
Zrois_within1=Y(:,nsources,idx,:);                       % Within-rois matrix  (subjects x source ROIs x source ROIs x condition)  
Zrois_between=Y(:,nsources,setdiff(ntargets,idx),:);     % Between-rois matrix (subjects x source ROIs x setdiff(target ROIs,source ROIs) x conditions) 
[ok2,idx]=ismember(alltargets(ntargets),allsources);
if any(~ok2), conn_disp('WARNING: some of the selected Set-2 ROIs have not been defined as ''sources'' in the CONN toolbox first-level analyses. Within-network computations can only be computed for ROIs defined as ''sources''. Within-Set-2 computations will now skip these ROIs'); end
Zrois_within2=Y(:,idx(ok2),ntargets(ok2),:);                       % Within-rois matrix  (subjects x source ROIs x source ROIs x condition)  
conn_disp('fprintf','\nSet-1 ROIs: %s\n',sprintf('%s ',allsources{nsources}));
conn_disp('fprintf','Set-2 ROIs: %s\n',sprintf('%s ',alltargets{ntargets}));
conn_disp('fprintf','Within matrix Set-1 (%d x %d ROIs; %d valid ROI pairs; %d subjects)\n',size(Zrois_within1,2),size(Zrois_within1,3),nnz(~any(any(isnan(Zrois_within1),1),4)),size(Zrois_within1,1));
conn_disp('fprintf','Within matrix Set-2 (%d x %d ROIs; %d valid ROI pairs; %d subjects)\n',size(Zrois_within2,2),size(Zrois_within2,3),nnz(~any(any(isnan(Zrois_within2),1),4)),size(Zrois_within2,1));
conn_disp('fprintf','Between matrix (%d x %d ROIs; %d valid ROI pairs; %d subjects)\n',size(Zrois_between,2),size(Zrois_between,3),nnz(~any(any(isnan(Zrois_between),1),4)),size(Zrois_between,1));

%% computes the within/between averages
mask_within1=isnan(Zrois_within1);
mask_within2=isnan(Zrois_within2);
mask_between=isnan(Zrois_between);
Zrois_within1(mask_within1)=0;
Zrois_within2(mask_within2)=0;
Zrois_between(mask_between)=0;
Z_within1=permute(sum(sum(Zrois_within1,2),3)./sum(sum(~mask_within1,2),3),[1,4,2,3]);
Z_within2=permute(sum(sum(Zrois_within2,2),3)./sum(sum(~mask_within2,2),3),[1,4,2,3]);
Z_between=permute(sum(sum(Zrois_between,2),3)./sum(sum(~mask_between,2),3),[1,4,2,3]);

%% performs stats test
conn_disp('fprintf','\nSet-1 Within-network test (average effect between Set-1 ROIs):\n');
[h1,F1,p1,dof1,statsname1]=conn_glm(X,Z_within1,csubjecteffects,cconditions);
conn_disp('fprintf',' %s(%s)=%.2f  (p=%.6f)\n',statsname1,sprintf('%d ',dof1),F1,p1);
conn_disp('fprintf','Set-2 Within-network test (average effect between Set-2 ROIs):\n');
[h2,F2,p2,dof2,statsname2]=conn_glm(X,Z_within2,csubjecteffects,cconditions);
conn_disp('fprintf',' %s(%s)=%.2f  (p=%.6f)\n',statsname2,sprintf('%d ',dof2),F2,p2);
conn_disp('fprintf','Between-networks test (average effect between Set-1 & Set-2 ROIs):\n');
[h0,F0,p0,dof0,statsname0]=conn_glm(X,Z_between,csubjecteffects,cconditions);
conn_disp('fprintf',' %s(%s)=%.2f  (p=%.6f)\n',statsname0,sprintf('%d ',dof0),F0,p0);

if ~nargout
    [filename,filepath]=uiputfile('*.csv','Save within/between data to file');
    if ~isequal(filename,0)
        fh=fopen(fullfile(filepath,filename),'wt');
        fprintf(fh,'Statistical test (within/between ROI network tests)\n');
        fprintf(fh,'Set-1 ROIs:,"%s"\n',sprintf('%s ',allsources{nsources}));
        fprintf(fh,'Set-2 ROIs:,"%s"\n',sprintf('%s ',alltargets{ntargets}));
        fprintf(fh,'Conditions:,"%s"\n',sprintf('%s ',CONN_x.Setup.conditions.names{nconditions}));
        fprintf(fh,'Between-conditions contrast:,%s\n',mat2str(cconditions));
        fprintf(fh,'Between-subject effects:,"%s"\n',sprintf('%s ',CONN_x.Setup.l2covariates.names{nsubjecteffects}));
        fprintf(fh,'Between-subjects contrast:,%s\n',mat2str(csubjecteffects));
        fprintf(fh,'\n,Within_network1,Within_network2,Between_networks\n');
        for n=1:size(Z_within1,1)
            fprintf(fh,'Subject %d,%s,%s,%s\n',n,mat2str(Z_within1(n,:)),mat2str(Z_within2(n,:)),mat2str(Z_between(n,:)));
        end
        fprintf(fh,'\nEffect size,%s,%s,%s\n',mat2str(h1),mat2str(h2),mat2str(h0));
        fprintf(fh,'Stat name,%s,%s,%s\n',statsname1,statsname2,statsname0);
        fprintf(fh,'Stat dof (degrees of freedom),%s,%s,%s\n',mat2str(dof1),mat2str(dof2),mat2str(dof0));
        fprintf(fh,'Stat value,%s,%s,%s\n',mat2str(F1),mat2str(F2),mat2str(F0));
        fprintf(fh,'p (p-value),%s,%s,%s\n',mat2str(p1),mat2str(p2),mat2str(p0));
        fclose(fh);
    end
    varargout={};
else
    varargout={Z_within1,Z_within2,Z_between};
end
