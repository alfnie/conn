function data=conn_randomise(X,Y,c1,c2,aopt,Pthr,Pthr_type,Pthr_side,niters,filename,overwrite,datatype,analysismask,groupingsamples,isdisplay)
% conn_randomise(X,Y,C1,C2,opt,Pthr,Pthr_type,Pthr_side,niters,
% X: [n,m] (# of samples x # of effects) design matrix
% Y: [n,v,r1,r2] (# of samples x # of variables x # of seed ROIs x # of target ROIs) data matrix
% Y: [n,v,r...]  (# of samples x # of variables x # of voxels) data matrix
% C1: [a,m] left-hand contrast matrix
% C2: [b,v] right-hand contrast matrix
% opt: conn_glm analysis option ([], 'collapse_none', 'collapse_predictors', 'collapse_outcomes', 'collapse_all', 'collapse_all_rao', 'collapse_all_bartlett', 'collapse_all_satterthwaite', 'collapse_all_sphericity')
%    Y(:,:,k) = X*B+e
%    C1*B*C2' = 0
% Pthr: height threshold [.05 .01 .005 .001 .0005 .0001]
% Pthr_type: threshold type 1 (p-unc) 2 (p-fdr rows -obsolete-) 3 (p-fdr all) 4 (T/F/X stat) 5 (p-TFCE)
% Pthr_side: threshold direction 1 (one-sided positive), 2 (one-sided negative), 3 (two-sided)


if nargin<5||isempty(aopt), aopt=[]; end
if nargin<6||isempty(Pthr), Pthr=[.05 .01 .005 .001 .0005 .0001]; end
if nargin<7||isempty(Pthr_type), Pthr_type=ones(size(Pthr)); end
if nargin<8||isempty(Pthr_side), Pthr_side=ones(size(Pthr)); end
if nargin<9||isempty(niters), niters=1000; end
if nargin<10||(ischar(filename)&&isempty(filename)), filename='temporal_conn_randomise.mat'; end
if nargin<11||isempty(overwrite), overwrite=false; end
if nargin<12||isempty(datatype), datatype='none'; end % 'none'|'matrix'|'surface'|volumecoords or [3xN] voxel indices
if nargin<13, analysismask=[]; end 
if nargin<14, groupingsamples=[]; end
if nargin<15, isdisplay=[]; end
if numel(niters)==1, BlockDistributed=0; Nseed0=0; 
else BlockDistributed=niters(2); niters=niters(1); Nseed0=(BlockDistributed-1)*niters; 
end
VERSION=2;
maxT=10; % stat-resolution
maxTmax=maxT*100; % stat-max
[Nn,Nv,Nr1,Nr2,Nr3]=size(Y);
Nr=Nr1*Nr2*Nr3; 

DOPREV=false; % back-compatibility (now obsolete after release 9.x)
DOCUMUL=true; % allow cumulative/distributed computations
DOSCALEF=true; % back-compatibility
doperminsteadofrand=false;
doflipsigninsteadofrand=false;
randdatainsteadofdesign=false;

% baseline results
model=struct('X',X,'c1',c1,'c2',c2,'dims',size(Y),'opts',aopt);
results=[]; % [results.h,results.F,results.p,results.dof,results.statsname]=conn_glm(X,Y,c1,c2);

% remove orthogonal contrasts
nc=null(c1);
xnc=X*nc;
iX=pinv(X'*X);
ix=pinv(xnc'*xnc);
r=zeros(Nn,size(c2,1),Nr);
%s=r;
for nr=1:Nr
    r(:,:,nr)=Y(:,:,nr)*c2'-X*(iX*(X'*Y(:,:,nr)*c2'));              % full-model residual
    %r(:,:,nr)=Y(:,:,nr)*c2'-xnc*(ix*(xnc'*Y(:,:,nr)*c2'));          % partial-model residual
    %s(:,:,nr)=xnc*(ix*(xnc'*Y(:,:,nr)*c2'));
end
x=X*c1'-xnc*(ix*(xnc'*X*c1'));        % removes from X
univariate=size(c1,1)==1&size(c2,1)==1&isempty(aopt);
if univariate&&numel(Pthr)==1, Pthr=repmat(Pthr,[1 3]); Pthr_type=repmat(Pthr_type,[1 3]); Pthr_side=1:3; end

% initialize
newPthr=Pthr;
newPthr_type=Pthr_type;
newPthr_side=Pthr_side;
if ~overwrite&&isstruct(filename)
    Hist_Voxel_stat=filename.Hist_Voxel_stat;
    Dist_Voxel_statmax=filename.Dist_Voxel_statmax;
    Hist_Cluster_size=filename.Hist_Cluster_size;
    Hist_Cluster_mass=filename.Hist_Cluster_mass;
    Hist_Cluster_score=filename.Hist_Cluster_score;
    Dist_Cluster_sizemax=filename.Dist_Cluster_sizemax;
    Dist_Cluster_massmax=filename.Dist_Cluster_massmax;
    Dist_Cluster_scoremax=filename.Dist_Cluster_scoremax;
    Hist_Seed_size=filename.Hist_Seed_size;
    Hist_Seed_mass=filename.Hist_Seed_mass;
    Hist_Seed_score=filename.Hist_Seed_score;
    Dist_Seed_sizemax=filename.Dist_Seed_sizemax;
    Dist_Seed_massmax=filename.Dist_Seed_massmax;
    Dist_Seed_scoremax=filename.Dist_Seed_scoremax;
    Hist_Network_size=filename.Hist_Network_size;
    Hist_Network_mass=filename.Hist_Network_mass;
    Hist_Network_score=filename.Hist_Network_score;
    Dist_Network_sizemax=filename.Dist_Network_sizemax;
    Dist_Network_massmax=filename.Dist_Network_massmax;
    Dist_Network_scoremax=filename.Dist_Network_scoremax;
    Pthr=filename.Pthr;
    Pthr_type=filename.Pthr_type;
    Pthr_side=filename.Pthr_side;
    if isfield(filename,'datatype'), datatype=filename.datatype; end
elseif ~overwrite&&ischar(filename)&&conn_existfile(filename)
    thisversion=VERSION;
    VERSION=[];
    conn_loadmatfile(filename);
    assert(isequal(VERSION,thisversion),'Incompatible conn_randomise data files. Please re-compute second-level analyses or delete existing nonparametric_*.mat files');
else
    Hist_Voxel_stat={};
    Dist_Voxel_statmax={};
    Hist_Cluster_size={};
    Hist_Cluster_mass={};
    Hist_Cluster_score={};
    Dist_Cluster_sizemax={};
    Dist_Cluster_massmax={};
    Dist_Cluster_scoremax={};
    Hist_Seed_size={};
    Hist_Seed_mass={};
    Hist_Seed_score={};
    Dist_Seed_sizemax={};
    Dist_Seed_massmax={};
    Dist_Seed_scoremax={};
    Hist_Network_size={};
    Hist_Network_mass={};
    Hist_Network_score={};
    Dist_Network_sizemax={};
    Dist_Network_massmax={};
    Dist_Network_scoremax={};
    Pthr=[];
    Pthr_type=[];
    Pthr_side=[];
end
idx=find(~ismember([newPthr(:),newPthr_type(:),newPthr_side(:)],[Pthr(:),Pthr_type(:),Pthr_side(:)],'rows'));
if DOCUMUL, 
    try
        [tok,tidx]=ismember([Pthr(:),Pthr_type(:),Pthr_side(:)],[newPthr(:),newPthr_type(:),newPthr_side(:)],'rows');
        Nseed0offset=0;
        for n1=1:numel(newPthr),
            tidx2=find(tok>0 & tidx==n1);
            if ~isempty(tidx2), 
                ndone=sum(arrayfun(@(n)size(Dist_Cluster_sizemax{n},1),tidx2));
                Nseed0offset=max(Nseed0offset,ndone);
                conn_disp('fprintf','%d runs (%d randomization samples) already matching threshold values %s. Adding %d new randomization samples\n',numel(tidx2),ndone,mat2str([newPthr(n1) newPthr_type(n1) newPthr_side(n1)]),niters); 
            end
        end
        Nseed0=Nseed0+Nseed0offset;
    end
    idx=1:numel(newPthr); 
    if BlockDistributed>0
        filename=conn_prepend(['parallel_',num2str(BlockDistributed),'_'],filename);
        Hist_Voxel_stat={};
        Dist_Voxel_statmax={};
        Hist_Cluster_size={};
        Hist_Cluster_mass={};
        Hist_Cluster_score={};
        Dist_Cluster_sizemax={};
        Dist_Cluster_massmax={};
        Dist_Cluster_scoremax={};
        Hist_Seed_size={};
        Hist_Seed_mass={};
        Hist_Seed_score={};
        Dist_Seed_sizemax={};
        Dist_Seed_massmax={};
        Dist_Seed_scoremax={};
        Hist_Network_size={};
        Hist_Network_mass={};
        Hist_Network_score={};
        Dist_Network_sizemax={};
        Dist_Network_massmax={};
        Dist_Network_scoremax={};
        Pthr=[];
        Pthr_type=[];
        Pthr_side=[];
    end
end
nthr=numel(idx);
if ~nthr,     
    data=struct('VERSION',VERSION,'Pthr',Pthr,'Pthr_type',Pthr_type,'Pthr_side',Pthr_side,'maxT',maxT,'Hist_Voxel_stat',{Hist_Voxel_stat},'Dist_Voxel_statmax',{Dist_Voxel_statmax},'Hist_Cluster_size',{Hist_Cluster_size},'Hist_Cluster_mass',{Hist_Cluster_mass},'Hist_Cluster_score',{Hist_Cluster_score},'Dist_Cluster_sizemax',{Dist_Cluster_sizemax},'Dist_Cluster_massmax',{Dist_Cluster_massmax},'Dist_Cluster_scoremax',{Dist_Cluster_scoremax},'Hist_Seed_size',{Hist_Seed_size},'Hist_Seed_mass',{Hist_Seed_mass},'Hist_Seed_score',{Hist_Seed_score},'Dist_Seed_sizemax',{Dist_Seed_sizemax},'Dist_Seed_massmax',{Dist_Seed_massmax},'Dist_Seed_scoremax',{Dist_Seed_scoremax},'Hist_Network_size',{Hist_Network_size},'Hist_Network_mass',{Hist_Network_mass},'Hist_Network_score',{Hist_Network_score},'Dist_Network_sizemax',{Dist_Network_sizemax},'Dist_Network_massmax',{Dist_Network_massmax},'Dist_Network_scoremax',{Dist_Network_scoremax});
    return;
end
thridx=numel(Pthr)+(1:nthr);
Hist_Voxel_stat=[Hist_Voxel_stat,repmat({sparse(1+maxTmax*100,1)},1,nthr)];
Dist_Voxel_statmax=[Dist_Voxel_statmax,repmat({zeros(niters,1)},1,nthr)];
Hist_Cluster_size=[Hist_Cluster_size,repmat({sparse(1+Nr,1)},1,nthr)];
Hist_Cluster_mass=[Hist_Cluster_mass,repmat({sparse(1+maxTmax*Nr,1)},1,nthr)];
Hist_Cluster_score=[Hist_Cluster_score,repmat({sparse(1+maxTmax*Nr,1)},1,nthr)];
Dist_Cluster_sizemax=[Dist_Cluster_sizemax,repmat({zeros(niters,1)},1,nthr)];
Dist_Cluster_massmax=[Dist_Cluster_massmax,repmat({zeros(niters,1)},1,nthr)];
Dist_Cluster_scoremax=[Dist_Cluster_scoremax,repmat({zeros(niters,1)},1,nthr)];
Hist_Seed_size=[Hist_Seed_size,repmat({sparse(1+Nr,1)},1,nthr)];
Hist_Seed_mass=[Hist_Seed_mass,repmat({sparse(1+maxTmax*Nr,1)},1,nthr)];
Hist_Seed_score=[Hist_Seed_score,repmat({sparse(1+maxTmax*Nr,1)},1,nthr)];
Dist_Seed_sizemax=[Dist_Seed_sizemax,repmat({zeros(niters,1)},1,nthr)];
Dist_Seed_massmax=[Dist_Seed_massmax,repmat({zeros(niters,1)},1,nthr)];
Dist_Seed_scoremax=[Dist_Seed_scoremax,repmat({zeros(niters,1)},1,nthr)];
Hist_Network_size=[Hist_Network_size,repmat({sparse(1+Nr,1)},1,nthr)];
Hist_Network_mass=[Hist_Network_mass,repmat({sparse(1+maxTmax*Nr,1)},1,nthr)];
Hist_Network_score=[Hist_Network_score,repmat({sparse(1+maxTmax*Nr,1)},1,nthr)];
Dist_Network_sizemax=[Dist_Network_sizemax,repmat({zeros(niters,1)},1,nthr)];
Dist_Network_massmax=[Dist_Network_massmax,repmat({zeros(niters,1)},1,nthr)];
Dist_Network_scoremax=[Dist_Network_scoremax,repmat({zeros(niters,1)},1,nthr)];
Pthr=[Pthr,newPthr(idx(:)')];
Pthr_type=[Pthr_type,newPthr_type(idx(:)')];
Pthr_side=[Pthr_side,newPthr_side(idx(:)')];
Nx=rank(X);
Nc0=rank(x);
Ns=rank(c2);
dof=size(Y,1)-Nx;
[nill,nill,nill,statsdof,statsname]=conn_glm(X,r(:,:,find(~any(any(isnan(r),1),2),1)),c1,[],aopt);
switch(statsname)
    case 'T'
        Xthr(thridx)=spm_invTcdf(max(0,min(1,1-Pthr(thridx))),statsdof);
        Xthrb(thridx)=spm_invTcdf(max(0,min(1,Pthr(thridx))),statsdof);
        Xthrc(thridx)=spm_invTcdf(max(0,min(1,1-Pthr(thridx)/2)),statsdof);
        spmXcdf=linspace(spm_invTcdf(1e-10,statsdof),spm_invTcdf(1-1e-10,statsdof),1e4);
        spmXcdf=cat(1,spmXcdf, spm_Tcdf(spmXcdf,statsdof));
    case 'F'
        Xthr(thridx)=spm_invFcdf(max(0,min(1,1-Pthr(thridx))),statsdof(1),statsdof(2));
        spmXcdf=linspace(spm_invFcdf(1e-10,statsdof(1),statsdof(2)),spm_invFcdf(1-1e-10,statsdof(1),statsdof(2)),1e4);
        spmXcdf=cat(1,spmXcdf, spm_Fcdf(spmXcdf,statsdof(1),statsdof(2)));
    case 'X'
        Xthr(thridx)=spm_invXcdf(max(0,min(1,1-Pthr(thridx))),statsdof);
        spmXcdf=linspace(spm_invXcdf(1e-10,statsdof),spm_invXcdf(1-1e-10,statsdof),1e4);
        spmXcdf=cat(1,spmXcdf, spm_Xcdf(spmXcdf,statsdof));
end
% if size(c1,1)==1&size(c2,1)==1,
%     Xthr(thridx)=spm_invTcdf(max(0,min(1,1-Pthr(thridx))),dof);
%     Xthrb(thridx)=spm_invTcdf(max(0,min(1,Pthr(thridx))),dof);
%     Xthrc(thridx)=spm_invTcdf(max(0,min(1,1-Pthr(thridx)/2)),dof);
%     statsdof=dof;
%     statsname='T';
% elseif size(c1,1)==1
%     Xthr(thridx)=spm_invFcdf(max(0,min(1,1-Pthr(thridx))),Ns,dof-Ns+1);
%     statsdof=[Ns,dof-Ns+1];
%     statsname='F';
% elseif size(c2,1)==1||Ns==1
%     Xthr(thridx)=spm_invFcdf(max(0,min(1,1-Pthr(thridx))),Nc0,dof);
%     statsdof=[Nc0,dof];
%     statsname='F';
% else
%     Xthr(thridx)=spm_invXcdf(max(0,min(1,1-Pthr(thridx))),Ns*Nc0);
%     statsdof=Ns*Nc0;
%     statsname='X';
% end
tmask=thridx(Pthr_type(thridx)==4); Xthr(tmask)=Pthr(tmask);
if DOPREV
    if univariate
        ix=pinv(x'*x);
        Cyy=reshape(sum(abs(r).^2,1),[1,Nr]);
        k_dof=sqrt(dof*x'*x);
        r=r(:,:);
    else
        conn_glm_steps(1,x,r);
    end
% else
%     conn_glm_steps(1,X,r,c1);
end
if any(Pthr_type(thridx)==5)
    if ischar(datatype)
        if strcmp(datatype,'volume')||strcmp(datatype,'matrix'), adj=[];
        else adj=load(fullfile(fileparts(which(mfilename)),'utils','surf','surf_top.mat'),'A'); adj=adj.A;
        end
    elseif (isstruct(datatype)&&isfield(datatype,'xyz'))||(~issparse(datatype)&&size(datatype,1)==3), % volume-case (entering xyz voxel coordinates)
        if isstruct(datatype), datatype=datatype.xyz; end
        tsize=2+max(datatype,[],2);
        adj=struct('size',tsize','index',sub2ind(tsize,1+datatype(1,:),1+datatype(2,:),1+datatype(3,:)));
    elseif (isstruct(datatype)&&isfield(datatype,'faces'))||(~issparse(datatype)&&size(datatype,2)==3), % surface-case (entering triangular faces info)
        adj=spm_mesh_adjacency(datatype); 
    else adj=datatype; % explicit sparse adjacency
    end
end
% if 1,%any(Pthr_type(thridx)>1)
%     switch(statsname)
%         case 'T',
%             spmXcdf=linspace(spm_invTcdf(1e-10,statsdof),spm_invTcdf(1-1e-10,statsdof),1e4);
%             spmXcdf=cat(1,spmXcdf, spm_Tcdf(spmXcdf,statsdof));
%         case 'F',
%             spmXcdf=linspace(spm_invFcdf(1e-10,statsdof(1),statsdof(2)),spm_invFcdf(1-1e-10,statsdof(1),statsdof(2)),1e4);
%             spmXcdf=cat(1,spmXcdf, spm_Fcdf(spmXcdf,statsdof(1),statsdof(2)));
%         case 'X',
%             spmXcdf=linspace(spm_invXcdf(1e-10,statsdof),spm_invXcdf(1-1e-10,statsdof),1e4);
%             spmXcdf=cat(1,spmXcdf, spm_Xcdf(spmXcdf,statsdof));
%     end
% end
if isempty(isdisplay), isdisplay=~all(get(0,'screensize')==1); end
if isdisplay>0, ht=conn_waitbar(0,'Updating non-parametric statistics. Please wait'); 
elseif isdisplay==0, conn_disp('Updating non-parametric statistics. Please wait'); 
end

% run simulations: randomise residual model
try, warning('off','MATLAB:RandStream:ActivatingLegacyGenerators'); warning('off','MATLAB:RandStream:ReadingInactiveLegacyGeneratorState'); warning('off','MATLAB:nearlySingularMatrix'); warning('off','MATLAB:singularMatrix'); end
try, randstate=rand('state'); end

for niter=1:niters
    rand('seed',Nseed0+niter-1);
    if DOPREV
        xnew=x;
        if doperminsteadofrand % permutation of residuals
            shiftsign=randperm(size(r,1))';
            xnew=x(shiftsign,:);
        else % randomisation of residuals
            if isempty(groupingsamples), shiftsign=rand(size(r,1),1)<.5;
            else shiftsign=full((groupingsamples*rand(size(groupingsamples,2),1))<.5);
            end
            xnew(shiftsign,:)=-xnew(shiftsign,:);
        end
        if univariate % speed-up for univariate tests (T stats)
            Cxy=xnew'*r;
            b=ix*Cxy;
            e2=Cyy-b.*Cxy;
            f=b./max(eps,sqrt(abs(e2)))*k_dof;
            f=f(:);
        else % general but much slower cases (F&Wilks lambda stats)
            [h,f]=conn_glm_steps(2,xnew,r,analysismask);
            f=f(:);
        end
    elseif randdatainsteadofdesign,
        rnew=r;
        if doperminsteadofrand % permutation of residuals
            shiftsign=randperm(size(r,1))';
            rnew(:,:)=rnew(shiftsign,:);
        elseif doflipsigninsteadofrand % randomisation of residuals (sign flip)
            if isempty(groupingsamples), shiftsign=rand(size(r,1),1)<.5;
            else shiftsign=full((groupingsamples*rand(size(groupingsamples,2),1))<.5);
            end
            rnew(shiftsign,:)=-rnew(shiftsign,:);
        else % randomisation of residuals (orth proj)
            if isempty(groupingsamples), rnew(:,:)=orth(rand(Nn)-.5)*rnew(:,:);
            else rnew(:,:)=full(groupingsamples)*orth(rand(size(groupingsamples,2))-.5)*pinv(full(groupingsamples))*rnew(:,:);
            end
        end
        %rnew=rnew*sqrt(size(r,1)/(size(r,1)-1))+s;
        %[h,f]=conn_glm_steps(2,X,rnew,analysismask);
        if ~isempty(analysismask), [h,f]=conn_glm(X,rnew(:,:,analysismask),c1,[],aopt);
        else                       [h,f]=conn_glm(X,rnew,c1,[],aopt);
        end
        f=f(:);
    else
        xnew=X;
        if doperminsteadofrand % permutation of residuals
            shiftsign=randperm(size(r,1))';
            xnew=xnew(shiftsign,:);
        elseif doflipsigninsteadofrand % randomisation of residuals (sign flip)
            if isempty(groupingsamples), shiftsign=rand(size(r,1),1)<.5;
            else shiftsign=full((groupingsamples*rand(size(groupingsamples,2),1))<.5);
            end
            xnew(shiftsign,:)=-xnew(shiftsign,:);
        else % randomisation of residuals (orth proj)
            if isempty(groupingsamples), xnew=orth(rand(Nn)-.5)'*xnew;
            else xnew=full(groupingsamples)*orth(rand(size(groupingsamples,2))-.5)'*pinv(full(groupingsamples))*xnew;
            end
        end
        if ~isempty(analysismask), [h,f]=conn_glm(xnew,r(:,:,analysismask),c1,[],aopt);
        else                       [h,f]=conn_glm(xnew,r,c1,[],aopt);
        end
        f=f(:);
    end
    if ~univariate||any(Pthr_type(thridx)==2|Pthr_type(thridx)==3),
        p=nan(size(f));
        p(~isnan(f))=1-max(0,min(1,interp1(spmXcdf(1,:),spmXcdf(2,:),f(~isnan(f)),'linear','extrap')));
        pb=1-p;
        pc=2*min(p,1-p);
        if univariate
            if any(Pthr_type(thridx)==2&Pthr_side(thridx)==1), p2=reshape(conn_fdr(reshape(p,[Nr1,Nr2]),2),size(p)); end % note: obsolete Pthr_type==2
            if any(Pthr_type(thridx)==2&Pthr_side(thridx)==2), pb2=reshape(conn_fdr(reshape(pb,[Nr1,Nr2]),2),size(p)); end
            if any(Pthr_type(thridx)==2&Pthr_side(thridx)==3), pc2=reshape(conn_fdr(reshape(pc,[Nr1,Nr2]),2),size(p)); end
            if any(Pthr_type(thridx)==3&Pthr_side(thridx)==1), p3=conn_fdr(p(:)); end
            if any(Pthr_type(thridx)==3&Pthr_side(thridx)==2), pb3=conn_fdr(pb(:)); end
            if any(Pthr_type(thridx)==3&Pthr_side(thridx)==3), pc3=conn_fdr(pc(:)); end
        else
            if any(Pthr_type(thridx)==2), p2=reshape(conn_fdr(reshape(p,[Nr1,Nr2]),2),size(p)); end
            if any(Pthr_type(thridx)==3), p3=conn_fdr(p(:)); end
        end
    end
    if any(Pthr_type(thridx)==5),
        if isequal(statsname,'T')
            if isempty(analysismask), tf=reshape(f,[Nr1,Nr2,Nr3]);
            else tf=zeros([Nr1,Nr2,Nr3]); tf(analysismask)=f;
            end
            if isequal(datatype,'matrix'), tf(1:size(tf,1),1:size(tf,1))=(tf(1:size(tf,1),1:size(tf,1))+tf(1:size(tf,1),1:size(tf,1))')/2; tf(~triu(ones(size(tf)),1))=0; end
            [TFCEpos,TFCEpospeaks,TFCEposd]=conn_tfce(abs(tf).*(tf>0),'Hmin',1,'adjacency',adj); % note: minH=1
            [TFCEneg,TFCEnegpeaks,TFCEnegd]=conn_tfce(abs(tf).*(tf<0),'Hmin',1,'adjacency',adj); 
            TFCE=TFCEpos.*(tf>0)-TFCEneg.*(tf<0);
            TFCEpeaks=(TFCEpospeaks&(tf>0))|(TFCEnegpeaks&(tf<0));
            TFCEd=TFCEposd.*(tf>0)-TFCEnegd.*(tf<0);
            if ~isempty(analysismask), TFCE=TFCE(analysismask); TFCEpeaks=TFCEpeaks(analysismask); TFCEd=TFCEd(analysismask); end
        else
            if isempty(analysismask), tf=reshape(f,[Nr1,Nr2,Nr3]);
            else tf=zeros([Nr1,Nr2,Nr3]); tf(analysismask)=f;
            end
            if isequal(datatype,'matrix'), tf(1:size(tf,1),1:size(tf,1))=(tf(1:size(tf,1),1:size(tf,1))+tf(1:size(tf,1),1:size(tf,1))')/2; tf(~triu(ones(size(tf)),1))=0; end
            [TFCE,TFCEpeaks,TFCEd]=conn_tfce(sqrt(abs(tf)),'Hmin',1,'adjacency',adj); 
            if ~isempty(analysismask), TFCE=TFCE(analysismask); TFCEpeaks=TFCEpeaks(analysismask); TFCEd=TFCEd(analysismask); end
        end
    end
    for i=thridx % all new threshold values
        if Pthr_type(i)==5, 
            if ~isequal(statsname,'T'), f=TFCE; fpeaks=TFCEpeaks; fd=TFCEd;
            elseif Pthr_side(i)==1,     f=TFCE; fpeaks=TFCEpeaks&(TFCE>0); fd=TFCEd;
            elseif Pthr_side(i)==2,     f=-TFCE; fpeaks=TFCEpeaks&(TFCE<0); fd=-TFCEd;
            elseif Pthr_side(i)==3,     f=abs(TFCE); fpeaks=TFCEpeaks; fd=abs(TFCEd);
            else error('incorrect option Pthr_side');
            end
            % maximum voxel-stat & histogram of voxel-stat
            Dist_Voxel_statmax{i}(niter)=max(f(:));  % maximum voxel-stat
            Hist_Voxel_stat{i}=Hist_Voxel_stat{i}+sparse(1+round(min(maxTmax*100,maxT*100*max(0,fd(fpeaks)/sqrt(Nr)))),1,1/max(eps,nnz(fpeaks)),maxTmax*100+1,1);   % histogram of peak stats (note: Hist_(1) corresponds to stats<=0)
        else
            % maximum voxel-stat & histogram of voxel-stat
            Dist_Voxel_statmax{i}(niter)=max(abs(f(:)));  % maximum voxel-stat
            Hist_Voxel_stat{i}=Hist_Voxel_stat{i}+sparse(1+round(min(maxTmax*100,maxT*100*max(0,f(~isnan(f))))),1,1/max(eps,nnz(~isnan(f))),maxTmax*100+1,1);   % histogram of voxel stats (note: Hist_(1) corresponds to stats<=0)

            if univariate
                if Pthr_type(i)==1
                    if Pthr_side(i)==1,     show=f(:)>=Xthr(i); %ithres=Xthr(i);
                    elseif Pthr_side(i)==2, show=f(:)<=Xthrb(i); %ithres=abs(Xthrb(i));
                    elseif Pthr_side(i)==3, show=abs(f(:))>=Xthrc(i); %ithres=Xthrc(i);
                    else error('incorrect option Pthr_side');
                    end
                elseif Pthr_type(i)==2
                    if Pthr_side(i)==1,     show=p2(:)<=Pthr(i); %ithres=Xthr(i);
                    elseif Pthr_side(i)==2, show=pb2(:)<=Pthr(i); %ithres=abs(Xthrb(i));
                    elseif Pthr_side(i)==3, show=pc2(:)<=Pthr(i); %ithres=Xthrc(i);
                    else error('incorrect option Pthr_side');
                    end
                elseif Pthr_type(i)==3
                    if Pthr_side(i)==1,     show=p3(:)<=Pthr(i); %ithres=Xthr(i);
                    elseif Pthr_side(i)==2, show=pb3(:)<=Pthr(i); %ithres=abs(Xthrb(i));
                    elseif Pthr_side(i)==3, show=pc3(:)<=Pthr(i); %ithres=Xthrc(i);
                    else error('incorrect option Pthr_side');
                    end
                elseif Pthr_type(i)==4
                    if Pthr_side(i)==1,     show=f(:)>=Xthr(i); %ithres=Xthr(i);
                    elseif Pthr_side(i)==2, show=f(:)<=-Xthr(i); %ithres=Xthr(i);
                    elseif Pthr_side(i)==3, show=abs(f(:))>=Xthr(i); %ithres=Xthr(i);
                    else error('incorrect option Pthr_side');
                    end
                else error('incorrect option Pthr_type');
                end
            elseif Pthr_type(i)==1, show=p(:)<=Pthr(i); %ithres=Xthr(i);
            elseif Pthr_type(i)==2, show=p2(:)<=Pthr(i);%ithres=Xthr(i);
            elseif Pthr_type(i)==3, show=p3(:)<=Pthr(i);%ithres=Xthr(i);
            elseif Pthr_type(i)==4, show=f(:)>=Xthr(i); %ithres=Xthr(i);
            else error('incorrect option Pthr_type');
            end
            
            if ischar(datatype)&&strcmp(datatype,'matrix'),nclstats=2; % additional stats for ROI-to-ROI analyses
            else nclstats=1;
            end
            for nclstat=1:nclstats
                if nclstat==1, tdatatype=datatype;
                else tdatatype='network';
                end
                
                % cluster-size & cluster-mass for each cluster
                if isequal(statsname,'T')&max(f(show))>0&min(f(show))<0&nclstat==1 % note: remove '&nclstat==1' to have network-level stats separate for positive/negative effects
                    if isempty(analysismask), tf=reshape(show&(f(:)>0),[Nr1,Nr2,Nr3]);
                    else tf=false([Nr1,Nr2,Nr3]); tf(analysismask)=show&(f(:)>0);
                    end
                    [nclL1,CLUSTER_labels1]=conn_clusters(tf,tdatatype);
                    if isempty(analysismask), tf=reshape(show&(f(:)<0),[Nr1,Nr2,Nr3]);
                    else tf=false([Nr1,Nr2,Nr3]); tf(analysismask)=show&(f(:)<0);
                    end
                    [nclL2,CLUSTER_labels2]=conn_clusters(tf,tdatatype);
                    CLUSTER_labels=CLUSTER_labels1;
                    CLUSTER_labels(CLUSTER_labels1==0&CLUSTER_labels2>0)=CLUSTER_labels2(CLUSTER_labels1==0&CLUSTER_labels2>0)+numel(nclL1);
                    if ~isempty(analysismask), CLUSTER_labels=CLUSTER_labels(analysismask); end
                    nclL=[reshape(nclL1,[],1);reshape(nclL2,[],1)];
                else
                    if isempty(analysismask), tf=reshape(show,[Nr1,Nr2,Nr3]);
                    else tf=false([Nr1,Nr2,Nr3]); tf(analysismask)=show;
                    end
                    [nclL,CLUSTER_labels]=conn_clusters(tf,tdatatype);
                    if ~isempty(analysismask), CLUSTER_labels=CLUSTER_labels(analysismask); end
                end
                mask=CLUSTER_labels>0;
                V=double(show);
                if ~DOSCALEF||~isequal(statsname,'T'), V(show)=abs(f(show));
                else V(show)=abs(f(show)).^2;
                end
                if nnz(mask), mclL=accumarray(CLUSTER_labels(mask),V(mask),[max([0,numel(nclL),max(CLUSTER_labels(mask))]),1]);
                else mclL=0;
                end
                sclL=mclL.^(3/2)./max(eps,nclL)/3; % cluster TFCE-score
                %mclL=accumarray(CLUSTER_labels(mask),V(mask)-ithres,[max([0,max(CLUSTER_labels(mask))]),1]);
                
                if nclstat==1
                    % maximum cluster-size & histogram of cluster-size
                    if isempty(nclL), nclL=0; end
                    Dist_Cluster_sizemax{i}(niter)=max(nclL); % maximum cluster-size
                    Hist_Cluster_size{i}=Hist_Cluster_size{i}+sparse(1+nclL,1,1/numel(nclL),Nr+1,1);   % histogram of cluster size (note: Hist_(1) corresponds to cluster size=0)
                    
                    % maximum cluster-mass & histogram of cluster-mass
                    if isempty(mclL), mclL=0; end
                    Dist_Cluster_massmax{i}(niter)=max(mclL);  % maximum cluster-mass
                    Hist_Cluster_mass{i}=Hist_Cluster_mass{i}+sparse(1+round(min(maxTmax*Nr,maxT*mclL)),1,1/numel(mclL),maxTmax*Nr+1,1);   % histogram of cluster mass (note: Hist_(1) corresponds to cluster mass=0)
                    
                    % maximum cluster-score & histogram of cluster-score
                    if isempty(sclL), sclL=0; end
                    Dist_Cluster_scoremax{i}(niter)=max(sclL);  % maximum cluster-score
                    Hist_Cluster_score{i}=Hist_Cluster_score{i}+sparse(1+round(min(maxTmax*Nr,maxT*sclL)),1,1/numel(sclL),maxTmax*Nr+1,1);   % histogram of cluster score (note: Hist_(1) corresponds to cluster score=0)
                else % additional stats for ROI-to-ROI analyses
                    % maximum seed-size & histogram of seed-size
                    if isempty(analysismask), tf=reshape(show,[Nr1,Nr2]);
                    else tf=zeros([Nr1,Nr2]); tf(analysismask)=show;
                    end
                    nsdL=sum(tf,2);
                    Dist_Seed_sizemax{i}(niter)=max(nsdL); % maximum seed-size
                    Hist_Seed_size{i}=Hist_Seed_size{i}+sparse(1+nsdL,1,1/numel(nsdL),Nr+1,1);   % histogram of seed size (note: Hist_(1) corresponds to seed size=0)
                    
                    % maximum seed-mass & histogram of seed-mass
                    if isempty(analysismask), tf=reshape(V,[Nr1,Nr2]);
                    else tf=zeros([Nr1,Nr2]); tf(analysismask)=V;
                    end
                    msdL=sum(tf,2);
                    Dist_Seed_massmax{i}(niter)=max(msdL);  % maximum seed-mass
                    Hist_Seed_mass{i}=Hist_Seed_mass{i}+sparse(1+round(min(maxTmax*Nr,maxT*msdL)),1,1/numel(msdL),maxTmax*Nr+1,1);   % histogram of seed mass (note: Hist_(1) corresponds to seed mass=0)
                    
                    % maximum seed-score & histogram of seed-score
                    ssdL=msdL.^(3/2)./max(eps,nsdL)/3; % seed TFCE-score
                    Dist_Seed_scoremax{i}(niter)=max(ssdL);  % maximum seed-score
                    Hist_Seed_score{i}=Hist_Seed_score{i}+sparse(1+round(min(maxTmax*Nr,maxT*ssdL)),1,1/numel(ssdL),maxTmax*Nr+1,1);   % histogram of seed score (note: Hist_(1) corresponds to seed score=0)
                    
                    % maximum network-size & histogram of network-size
                    if isempty(nclL), nclL=0; end
                    Dist_Network_sizemax{i}(niter)=max(nclL); % maximum cluster-size
                    Hist_Network_size{i}=Hist_Network_size{i}+sparse(1+nclL,1,1/numel(nclL),Nr+1,1);   % histogram of cluster size (note: Hist_(1) corresponds to cluster size=0)
                    
                    % maximum network-mass & histogram of network-mass
                    if isempty(mclL), mclL=0; end
                    Dist_Network_massmax{i}(niter)=max(mclL);  % maximum cluster-mass
                    Hist_Network_mass{i}=Hist_Network_mass{i}+sparse(1+round(min(maxTmax*Nr,maxT*mclL)),1,1/numel(mclL),maxTmax*Nr+1,1);   % histogram of cluster mass (note: Hist_(1) corresponds to cluster mass=0)
                    
                    % maximum network-score & histogram of network-score
                    if isempty(sclL), sclL=0; end
                    Dist_Network_scoremax{i}(niter)=max(sclL);  % maximum cluster-score
                    Hist_Network_score{i}=Hist_Network_score{i}+sparse(1+round(min(maxTmax*Nr,maxT*sclL)),1,1/numel(sclL),maxTmax*Nr+1,1);   % histogram of cluster score (note: Hist_(1) corresponds to cluster score=0)
                end
            end
        end
    end
    
    if ~rem(niter,max(1,ceil(niters/100))),
        if isdisplay>0, conn_waitbar(niter/niters,ht); 
        elseif isdisplay==0, fprintf('.'); 
        end
    end
end
for i=thridx 
    Hist_Voxel_stat{i}=Hist_Voxel_stat{i}/max(eps,sum(Hist_Voxel_stat{i}));
    Hist_Cluster_mass{i}=Hist_Cluster_mass{i}/max(eps,sum(Hist_Cluster_mass{i}));
    Hist_Cluster_size{i}=Hist_Cluster_size{i}/max(eps,sum(Hist_Cluster_size{i}));
    Hist_Cluster_score{i}=Hist_Cluster_score{i}/max(eps,sum(Hist_Cluster_score{i}));
    if isequal(datatype,'matrix')
        Hist_Seed_mass{i}=Hist_Seed_mass{i}/max(eps,sum(Hist_Seed_mass{i}));
        Hist_Seed_size{i}=Hist_Seed_size{i}/max(eps,sum(Hist_Seed_size{i}));
        Hist_Seed_score{i}=Hist_Seed_score{i}/max(eps,sum(Hist_Seed_score{i}));
        Hist_Network_mass{i}=Hist_Network_mass{i}/max(eps,sum(Hist_Network_mass{i}));
        Hist_Network_size{i}=Hist_Network_size{i}/max(eps,sum(Hist_Network_size{i}));
        Hist_Network_score{i}=Hist_Network_score{i}/max(eps,sum(Hist_Network_score{i}));
    end
end
if isdisplay>0, conn_waitbar('close',ht);
elseif isdisplay==0 fprintf('\n');
end
if ~nargout
    conn_savematfile(filename,'VERSION','model','results','Pthr','Pthr_type','Pthr_side','maxT','Hist_Voxel_stat','Dist_Voxel_statmax','Hist_Cluster_size','Hist_Cluster_mass','Hist_Cluster_score','Dist_Cluster_sizemax','Dist_Cluster_massmax','Dist_Cluster_scoremax','Hist_Seed_size','Hist_Seed_mass','Hist_Seed_score','Dist_Seed_sizemax','Dist_Seed_massmax','Dist_Seed_scoremax','Hist_Network_size','Hist_Network_mass','Hist_Network_score','Dist_Network_sizemax','Dist_Network_massmax','Dist_Network_scoremax','-v7.3');
else
    data=struct('VERSION',VERSION,'model',model,'results',results,'Pthr',Pthr,'Pthr_type',Pthr_type,'Pthr_side',Pthr_side,'maxT',maxT,'Hist_Voxel_stat',Hist_Voxel_stat,'Dist_Voxel_statmax',Dist_Voxel_statmax,'Hist_Cluster_size',{Hist_Cluster_size},'Hist_Cluster_mass',{Hist_Cluster_mass},'Hist_Cluster_score',{Hist_Cluster_score},'Dist_Cluster_sizemax',{Dist_Cluster_sizemax},'Dist_Cluster_massmax',{Dist_Cluster_massmax},'Dist_Cluster_scoremax',{Dist_Cluster_scoremax},'Hist_Seed_size',{Hist_Seed_size},'Hist_Seed_mass',{Hist_Seed_mass},'Hist_Seed_score',{Hist_Seed_score},'Dist_Seed_sizemax',{Dist_Seed_sizemax},'Dist_Seed_massmax',{Dist_Seed_massmax},'Dist_Seed_scoremax',{Dist_Seed_scoremax},'Hist_Network_size',{Hist_Network_size},'Hist_Network_mass',{Hist_Network_mass},'Hist_Network_score',{Hist_Network_score},'Dist_Network_sizemax',{Dist_Network_sizemax},'Dist_Network_massmax',{Dist_Network_massmax},'Dist_Network_scoremax',{Dist_Network_scoremax});
end
try, rand('state',randstate); end
try, warning('on','MATLAB:RandStream:ActivatingLegacyGenerators'); warning('on','MATLAB:RandStream:ReadingInactiveLegacyGeneratorState'); warning('on','MATLAB:nearlySingularMatrix'); warning('on','MATLAB:singularMatrix'); end
end

% note: code below no longer used (after release 9.x) and will be removed in future releases
function [h,F,p,dof,statsname,BB,EE,opt]=conn_glm_steps(steps,X,Y,mask)%,C,M,opt)
% CONN_GLM General Linear Model estimation and hypothesis testing.
%
%   [h,F,p,dof]=CONN_GLM(X,Y,C,M) estimates a linear model of the form Y = X*B + E
%   where Y is an observed matrix of response or output variables (rows are observations, columns are output variables)
%         X is an observed design or regressor matrix (rows are observations, columns are predictor variables)
%         B is a matrix of unknown regression parameters (rows are predictor variables, columns are output variables)
%         E is an matrix of unobserved multivariate normally distributed disturbances with zero mean and unknown covariance.
%   and tests a general null hypothesis of the form h = C*B*M' = 0
%   where C is matrix or vector of "predictor" contrasts (rows are contrasts, columns are predictor variables, defaults to C=eye(size(X,2)) )
% 		  M is matrix or vector of "outcome" contrasts (rows are contrasts, columns are output variables, defaults to M=eye(size(Y,2)) )
%
%   CONN_GLM returns the following information:
%		  h:   a matrix of estimated contrast effect sizes (h = C*B*M')
%		  F:   the test statistic(s) (T,F,or Chi2 value, depending on whether h is a scalar, a vector, or a matrix. See below)
%		  p:   p-value of the test(s)
%		  dof: degrees of freedom
%
%   Additional information:
%   By default CONN_GLM will use a T, F, or a Chi2 statistic for hypothesis testing depending on the size of h=C*B*M'. The default options are:
%                  when size(h)=[1,1]      -> T statistic (note: one-sided t-test)
%	                      			          Examples of use: one-sided two-sample t-test, linear regression
%                  when size(h)=[1,Ns]     -> F statistic (note: equivalent to two-sided t-test when Ns=1)
%   							  	 		  Examples of use: Hotelling's two sample t-square test, two-sided t-test, multivariate regression
%                  when size(h)=[Nc,1]     -> F statistic (note: equivalent to two-sided t-test when Nc=1)
%					  			     		  Examples of use: ANOVA, ANCOVA, linear regression omnibus test
%                  when size(h)=[Nc,Ns]    -> Wilks' Lambda statistic, Bartlett's Chi2 approximation
%								     		  Examples of use: MANOVA, MANCOVA, multivariate regression omnibus test, likelihood ratio test
%   The default option can be changed using the syntax CONN_GLM(X,Y,C,M,opt) where opt is one of the following character strings:
%               CONN_GLM(X,Y,C,M,'collapse_none') will perform a separate univariate-test on each of the elements of the matrix h=C*B*M'
%  			    CONN_GLM(X,Y,C,M,'collapse_outcomes') will perform a separate multivariate-test on each of the rows of the matrix h (collapsing across multiple outcome variables or outcome contrasts).
% 			    CONN_GLM(X,Y,C,M,'collapse_predictors') will perform a separate multivariate-test on each of the columns of the matrix h (collapsing across multiple predictor variables or predictor contrasts).
% 			    CONN_GLM(X,Y,C,M,'collapse_all') will perform a single multivariate test on the matrix h.
%
% Example of use:
%   % MANOVA (three groups, two outcome variables)
%   % Data preparation
%    N1=10;N2=20;N3=30;
%    Y1=randn(N1,2)+repmat([0,0],[N1,1]); % data for group 1 (N1 samples, population mean = [0,0])
%    Y2=randn(N2,2)+repmat([0,1],[N2,1]); % data for group 2 (N2 samples, population mean = [0,1])
%    Y3=randn(N3,2)+repmat([1,0],[N3,1]); % data for group 2 (N3 samples, population mean = [1,0])
%    Y=cat(1,Y1,Y2,Y3);
%    X=[ones(N1,1),zeros(N1,2); zeros(N2,1),ones(N2,1),zeros(N2,1); zeros(N3,2),ones(N3,1)];
%   % Sample data analyses
%    [h,F,p,dof]=conn_glm(X,Y,[1,-1,0;0,1,-1]); disp(['Multivariate omnibus test of non-equality of means across the three groups:']);disp([' Chi2(',num2str(dof),') = ',num2str(F),'   p = ',num2str(p)]);
%    [h,F,p,dof]=conn_glm(X,Y,[1,-1,0]); disp(['Multivariate test of non-equality of means between groups 1 and 2:']);disp([' F(',num2str(dof(1)),',',num2str(dof(2)),') = ',num2str(F),'   p = ',num2str(p)]);
%    [h,F,p,dof]=conn_glm(X,Y,[-1,1,0],eye(2),'collapse_none'); disp(['Univariate one-sided test of non-equality of means between groups 1 and 2 on each outcome variable:']);disp([' T(',num2str(dof),') = ',num2str(F(:)'),'   p = ',num2str(p(:)')]);
%

% alfnie@gmail.edu
% 04/03

persistent Nx Ns Na dofe Nc0 iX r ir kF1 C1; 
if nargin<4, mask=[]; end
opt=[];
if ~isempty(opt),switch(lower(opt)),case {'collapse_none','aa'},opt='AA';case {'collapse_outcomes','ab'},opt='AB';case {'collapse_predictors','ba'},opt='BA';case {'collapse_all','bb'},opt='BB';otherwise,error(['Unknown option ',opt]); end; end
if any(steps==1) %(1,X,Y,C)
    if nargin<4||isempty(mask), mask=1; end
    C=mask;
    [N1,Nx]=size(X);
    [N2,Ns,Na]=size(Y);
    if N1~=N2, error('wrong dimensions'); end
    %if nargin<4 || isempty(C), C=eye(Nx); end
    %if nargin<5 || isempty(M), M=speye(Ns,Ns); else, Ns=rank(M); end
    
    Nx=rank(X);
    dofe=N1-Nx;
    Nc0=rank(X*C');

    iX=pinv(X'*X);
    r=C*iX*C';
    ir=pinv(r);
    kF1=1./max(eps,r)*(dofe-Ns+1)/Ns;
    C1=C;
end

if any(steps==2) %(2,X,Y,mask)
    %if nargin<4 || isempty(C), C=eye(Nx); end
    %if nargin<5 || isempty(M), M=speye(Ns,Ns); end
    opt=[char(65+(size(C1,1)>1)),char(65+(Ns>1))]; 
    %if nargin<6 || isempty(opt), opt=[char(65+(size(C,1)>1)),char(65+(size(M,1)>1))]; else, opt=upper(opt); end
    if Na>1, Na_h=[]; nas=find(~any(any(isnan(Y)|isinf(Y),1),2)); else nas=1; end
    if ~isempty(mask), nas(~mask(nas))=[]; end
    univariate=(strcmp(opt,'AA')|strcmp(opt,'BA'));%&size(M,1)==Ns&size(M,2)==Ns&all(all(M==eye(Ns)));
    switch(opt),
        case 'AA',                          % h: [1,1]  T-test
            dof=              dofe;
            statsname=        'T';
        case 'AB',                          % h: [1,Ns] F-test
            dof=              [Ns,dofe-Ns+1];
            statsname=        'F';
        case 'BA',                          % h: [Nc,1] F-test
            dof=              [Nc0,dofe];
            statsname=        'F';
        case 'BB',                          % h: [Nc,Ns] Wilk's Lambda-test (Ns,dof,Nc0)
            if Ns==1 % rank-defficient EE case
                dof=              [Nc0,dofe];
                statsname=        'F';
            else
                dof=              [Ns*Nc0];
                statsname=        'X';
            end
    end
    Ball=reshape(iX*(X'*Y(:,:)),[size(iX,1),Ns,Na]);
    Eall=reshape(Y(:,:)-X*Ball(:,:),[size(Y,1),Ns,Na]);
    Ball(abs(Ball)<1e-10)=0;
    for na=nas(:)'
        h=C1*Ball(:,:,na);
        E=Eall(:,:,na);
        %B=iX*(X'*Y(:,:,na));
        %E=Y(:,:,na)-X*B;
        if univariate, EE=sum(abs(E).^2,1); % univariate case within matrix
        else, EE=E'*E; end          	% within matrix
        %else, if size(E,2)<size(E,1), EE=full(M*(E'*E)*M'); else, EE=E*M'; EE=full(EE'*EE); end; end          	% within matrix
        %h=B;%full(C*B*M');
        %h(abs(h)<1e-10)=0;
        
        switch(opt),
            case 'AA',                          % h: [1,1]  T-test
                k=                sqrt(diag(r)*diag(EE).');
                F=			      real(h./max(eps,k))*sqrt(dofe);
                if nargout>2, p=nan+zeros(size(F));idxvalid=find(~isnan(F));if ~isempty(idxvalid), p(idxvalid)=1-spm_Tcdf(F(idxvalid),dofe); end; end
            case 'AB',                          % h: [1,Ns] F-test
                F=			      (h*inv(EE)*h.')*kF1;
                %F=			      real(sum((h*pinv(EE)).*conj(h),2)./max(eps,diag(r)))*(dofe-Ns+1)/Ns;
                if nargout>2, p=nan+zeros(size(F));idxvalid=find(~isnan(F));if ~isempty(idxvalid), p(idxvalid)=1-spm_Fcdf(F(idxvalid),Ns,dofe-Ns+1); end; end
            case 'BA',                          % h: [Nc,1] F-test
                BB=h'*ir*h;                    % between matrix
                F=                (BB./max(eps,EE))*dofe/Nc0;
                %F=                real(diag(BB)./max(eps,diag(EE)))*dofe/Nc0;
                if nargout>2, p=nan+zeros(size(F));idxvalid=find(~isnan(F));if ~isempty(idxvalid), p(idxvalid)=1-spm_Fcdf(F(idxvalid),Nc0,dofe); end; end
            case 'BB',                          % h: [Nc,Ns] Wilk's Lambda-test (Ns,dof,Nc0)
                BB=h'*ir*h;                    % between matrix
                if Ns==1 % rank-defficient EE case
                    F=                (BB./max(eps,EE))*dofe/Nc0;
                    %F=                real(diag(BB)./max(eps,diag(EE)))*dofe/Nc0; F=F(1);
                    if nargout>2, p=nan+zeros(size(F));idxvalid=find(~isnan(F));if ~isempty(idxvalid), p(idxvalid)=1-spm_Fcdf(F(idxvalid),Nc0,dofe); end; end
                else
                    F=                -(dofe-1/2*(Ns-Nc0+1))*real(log(real(det(EE)./max(eps,det(EE+BB)))));
                    if nargout>2, p=nan+zeros(size(F));idxvalid=find(~isnan(F));if ~isempty(idxvalid), p(idxvalid)=1-spm_Xcdf(F(idxvalid),Ns*Nc0); end; end
                end
        end
        if Na>1
            if isempty(Na_h), Na_h=nan(size(h,1),size(h,2),Na);Na_F=nan(size(F,1),size(F,2),Na);Na_p=nan(size(F,1),size(F,2),Na); end
            Na_h(:,:,na)=h;
            Na_F(:,:,na)=F;
            if nargout>2, Na_p(:,:,na)=p; end
        end
    end
    if Na>1
        h=Na_h;
        F=Na_F;
        if nargout>2, p=Na_p; end
    end
end
end

