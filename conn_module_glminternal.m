function SPM = conn_module_glminternal(X,filenames,C1,C2,swd,effectnames,datanames,secondlevelanalyses,maskfile,subjectids)
% CONN_MODULE_GLMINTERNAL second-level model estimation internal function

%
% conn_module_glminternal(X,Y,c1,c2,folder) 
%   X     : design matrix (Nsubjects x Neffects)
%   Y     : data files (cell array Nsubjects x Nmeasures)
%   c1    : contrast between-subjects (Nc1 x Neffects) (default eye(size(X,2)))
%   c2    : contrast between-measures (Nc2 x Nmeasures) (default eye(size(Y,2)))
%   folder: folder where analysis are stored (default current folder)
%
% eg: conn_module_glminternal( ...
%    [1; 1; 1; 1] ,...
%    {'subject1.img'; 'subject2.img'; 'subject3.img'; 'subject4.img'} );
%     performs a one-sample t-test and stores the analysis results in the current folder
%
% eg: conn_module_glminternal( ...
%    [1 0; 1 0; 0 1; 0 1; 0 1],...
%    {'subject1_group1.img'; 'subject2_group1.img'; 'subject1_group2.img'; 'subject2_group2.img'; 'subject3_group2.img'},...
%    [1 -1]);
%     performs a two-sample t-test and stores the analysis results in the current folder
%
% eg: conn_module_glminternal( ...
%    [1; 1; 1; 1],...
%    {'subject1_time1.img', subject1_time2.img'; 'subject2_time1.img', subject2_time2.img'; 'subject3_time1.img', subject3_time2.img'; 'subject4_time1.img', subject4_time2.img'},...
%    1,...
%    [1 -1]);
%     performs a paired t-test and stores the analysis results in the current folder
%
% additional inputs: conn_module_glminternal(X,Y,c1,c2,folder,effectnames,datanames,secondlevelanalyses,maskfile)
%   effectnames         : names of subject effects (one name per column of X)
%   datanames           : names of measures (one name per column of Y)
%   secondlevelanalyses : type of second-level analyses  (1:all; 2:parametric only; 3:nonparametric only)
%   maskfile            : analysis mask
%

askall=~nargin;
if nargin<2||isempty(filenames), 
    filenames=cellstr(spm_select(inf,'any','Select source data (one or multiple images per subject)',[],pwd,'.*\.nii$')); 
    if isempty(filenames)||isempty(filenames{1}), return; end
end
teffectnames={};
if nargin<1||isempty(X), 
    str='Select design matrix (.csv, .txt, .mat file with one row per subject and one column per model regressor)';
    conn_disp(str); [filename,pathname]=uigetfile({'*',  'All Files (*)'},str);
    if isequal(filename,0), return; end
    X=fullfile(pathname,filename);
end
if ischar(X)
    conn_disp('fprintf','loading design information from %s\n',X);
    x=conn_loadtextfile(X);
    if isstruct(x)
        enames=fieldnames(x);
        X=[]; teffectnames={};
        for n=1:numel(enames),
            t=x.(enames{n});
            X=[X t]; 
            if size(t,2)>1, teffectnames=[teffectnames arrayfun(@(m)[enames{n},num2str(m)],1:size(t,2),'uni',0)];
            else teffectnames=[teffectnames enames(n)];
            end
        end
    else X=x; 
    end    
end
[Ns,Nx]=size(X);
if ischar(filenames), filenames=cellstr(filenames); end
if size(filenames,1)==1&&Ns>1, filenames=filenames'; end
[Ns2,Ny]=size(filenames);
if Ns~=Ns2&&Ns>0&&Ny==1&&~mod(Ns2,Ns)
    opts={'First all subjects for measure 1, followed by all subjects for measure 2, etc.','First all measures for subject 1, followed by all measures for subject 2, etc.'};
    answ=conn_questdlg('',sprintf('Order of data files (%d files, %d subjects, %d measures)',Ns2,Ns,Ns2/Ns),opts{[1,2,1]});
    if isempty(answ), return; end
    if strcmp(answ,opts{1}), filenames=reshape(filenames,Ns,[]);
    else filenames=reshape(filenames,[],Ns)';
    end
    [Ns2,Ny]=size(filenames);
end
if Ns~=Ns2, error('Incorrect dimensions (design matrix X and input data filenames should have the same number of rows)'); end
for n=1:numel(filenames), if isempty(fileparts(filenames{n})), filenames{n}=fullfile(pwd,filenames{n}); end; end
if nargin<3||isempty(C1), 
    if askall&&Nx>1
        answ=inputdlg({sprintf('Enter between-subjects contrast (vector/matrix with at least one row and %d columns)',Nx)},'',1,{mat2str(eye(Nx))},struct('Resize','on'));
        if numel(answ)~=1||isempty(answ{1}),return; end
        C1=str2num(answ{1});
    else C1=eye(Nx); 
    end
end
if nargin<4||isempty(C2), 
    if askall&&Ny>1
        answ=inputdlg({sprintf('Enter between-measures contrast (vector/matrix with at least one row and %d columns)',Ny)},'',1,{mat2str(eye(Ny))},struct('Resize','on'));
        if numel(answ)~=1||isempty(answ{1}),return; end
        C2=str2num(answ{1});
    else C2=eye(Ny); 
    end
end
if nargin<5||isempty(swd), 
    if askall, 
        str='Select output folder for SPM second-level analysis files'; 
        conn_disp(str); swd = uigetdir(pwd, str);
        if isequal(swd,0), return; end
    else swd=pwd; 
    end
end
if nargin<8||isempty(secondlevelanalyses), secondlevelanalyses=1; end % type of second-level analysis: 1:all; 2:param only; 3:nonparam only
if nargin<9, maskfile=''; end                                         % analysis mask file (restricts analyses to only within-mask voxels)
if nargin<10, subjectids=''; end                                      % subject id's (keeps track in SPM.xX.SubjectIDs field of subject IDs used in this analysis)
if size(C1,2)~=Nx, error('Incorrect dimensions (X and C1 should have the same number of columns)'); end
if size(C2,2)~=Ny, error('Incorrect dimensions (filenames and C2 should have the same number of columns)'); end
NANTO0=true; % skips SPM's explicit masking of NaN values
FORCEORTH=true; % enforces orthogonal rows in within-subject contrasts (C2)

% checks for missing data
MD=any(X~=0,2)&~any(isnan(X),2);
for nsub=1:Ns
    if MD(nsub)
        filename1=filenames(nsub,:);
        filename2=regexprep(filenames(nsub,:),',\d+$','');
        %for n1=1:numel(filename1),if isempty(fileparts(filename1{n1})), filename1{n1}=fullfile(pwd,filename1{n1}); filename2{n1}=fullfile(pwd,filename2{n1}); end; end
        SPM.xY.VY(nsub,:)=spm_vol(char(filename1))';
        for n1=1:numel(filename1),
            SPM.xY.VY(nsub,n1).fname=filename2{n1};
            if isfield(SPM.xY.VY(nsub,n1),'descrip')&&ischar(SPM.xY.VY(nsub,n1).descrip)&&strcmp(SPM.xY.VY(nsub,n1).descrip,'CONNlabel:MissingData'), MD(nsub)=false; end
        end
    end
end
if any(~MD)
    fprintf('warning: %d samples (%s) removed due to missing data\n',nnz(~MD),mat2str(find(~MD)')); 
    X=X(MD,:);
    SPM.xY.VY=SPM.xY.VY(MD,:);
    filenames=filenames(MD,:);
    if ~isempty(subjectids), subjectids=subjectids(MD); end
    Ns=nnz(MD);
end
if ~any(MD)
    fprintf('warning: no samples with valid data. Skipping 2nd-level analysis\n'); 
    return
end

pwd0=pwd;
swd=conn_fullfile(swd);
[ok,nill]=mkdir(swd);
cd(swd);

if nargin<6||isempty(effectnames), if numel(teffectnames)==Nx, effectnames=teffectnames; else effectnames=arrayfun(@(n)sprintf('effect %d',n),1:Nx,'uni',0); end; end % names of design matrix columns
if nargin<7||isempty(datanames), datanames=arrayfun(@(n)sprintf('measure %d',n),1:Ny,'uni',0); end % names of data matrix columns
SPM.xX_multivariate.Zfiles=filenames;
SPM.xX_multivariate.Znames=datanames;
SPM.xX_multivariate.Zcontr=C2;
[SPM.xX.isSurface,SPM.xX.isMtx]=conn_surf_dimscheck(SPM.xY.VY(1).dim); 
issurface=isfield(SPM.xX,'isSurface')&&SPM.xX.isSurface;
ismatrix=isfield(SPM.xX,'isMtx')&&SPM.xX.isMtx;

if ~ismatrix&&(NANTO0||~isequal(C2,eye(Ny))) % two-step procedure for SPM
    if FORCEORTH&&~isequal(C2,eye(Ny)), C2=spm_orth(C2','norm')'; end
    clear VY;
    nw=0;
    for nr=reshape(find(any(C2,2)),1,[]),
        nw=nw+1;
        datanames2{nw}='';
        for nc=reshape(find(C2(nr,:)),1,[])
            datanames2{nw}=[datanames2{nw},regexprep(sprintf(' %+g*%s',C2(nr,nc),datanames{nc}),{'\+1\*','\-1\*'},{'+','-'})];
        end
        for nsub=1:Ns
            b=0;
            for nc=reshape(find(C2(nr,:)),1,[])
                temp=spm_read_vols(SPM.xY.VY(nsub,nc));
                if NANTO0, temp(isnan(temp))=0; end
                b=b+C2(nr,nc)*temp;
            end
            a=SPM.xY.VY(nsub,1);
            a.fname=fullfile(swd,sprintf('data_%04d_%04d.nii',nsub,nw));
            a.n=[1 1];
            a.pinfo=[1;0;0];
            a.descrip=sprintf('within-subject contrast %s',mfilename);
            a=spm_write_vol(a,b);
            VY(nsub,nw)=a;
        end
    end
    if ~isequal(C2,eye(Ny)), datanames=regexprep(datanames2,'^ \+',''); end
    SPM.xY.VY=VY;
    [Ns2,Ny]=size(VY);
    C2=eye(Ny);
end
SPM.altestsmooth=0;

nrepeated=size(SPM.xY.VY,2);
if nrepeated>1
    SPM.xX.name={};
    for n01=1:numel(datanames),for n00=1:Nx,SPM.xX.name{n00,n01}=[effectnames{n00},'_',datanames{n01}]; end; end
    SPM.xX.name=SPM.xX.name(:)';
else SPM.xX.name=effectnames(:)';
end
SPM.xX_multivariate.X=X;
SPM.xX_multivariate.C=C1;
SPM.xX_multivariate.M=C2;
SPM.xX_multivariate.Xnames=effectnames;
SPM.xX_multivariate.Ynames=datanames;
SPM.xX.SelectedSubjects=logical(full(sparse(1:Ns,1,1,Ns,1)));
if ~isempty(subjectids), SPM.xX.SubjectIDs=subjectids; end
SPM.xX.X=kron(eye(nrepeated),X);
SPM.xX.iH     = [];
SPM.xX.iC     = 1:size(SPM.xX.X,2);
SPM.xX.iB     = [];
SPM.xX.iG     = [];
SPM.xGX       = [];
if nrepeated>1
    xVi=struct('I',[repmat((1:size(SPM.xY.VY,1))',[nrepeated,1]),reshape(repmat(1:nrepeated,[size(SPM.xY.VY,1),1]),[],1)],'var',[0,1],'dep',[1,0]);
    SPM.xVi=spm_non_sphericity(xVi);
end

save('SPM.mat','SPM');
spm_unlink('mask.img','mask.hdr','mask.nii');
files=cat(1,dir('spmT_*.cluster.mat'),dir('nonparametric_p*.mat'));
if ~isempty(files)
    files={files.name};
    spm_unlink(files{:});
end
if issurface&&isempty(maskfile), maskfile=fullfile(fileparts(which(mfilename)),'utils','surf','mask.surface.brainmask.nii'); end % note: surface anayses mask-out medial/subcortical segment
if ~isempty(maskfile), vmaskfile=spm_vol(char(maskfile)); end
if issurface||ismatrix||ismember(secondlevelanalyses,[1 3]) % nonparametric stats
    mask=ones(SPM.xY.VY(1).dim(1:3));
    [gridx,gridy]=ndgrid(1:SPM.xY.VY(1).dim(2),1:SPM.xY.VY(1).dim(3));
    xyz0=[gridx(:),gridy(:)]';
    donefirst=false;
    for n2=1:SPM.xY.VY(1).dim(1)
        xyz=[n2+zeros(1,size(xyz0,2)); xyz0; ones(1,size(xyz0,2))];
        y=spm_get_data(SPM.xY.VY(:)',xyz);
        maskthis=~any(isnan(y),1)&any(diff(y,1,1)~=0,1);
        if ~isempty(maskfile), maskthis=maskthis&spm_get_data(vmaskfile,pinv(vmaskfile.mat)*SPM.xY.VY(1).mat*xyz)>0; end
        mask(n2,:,:)=reshape(maskthis,[1 SPM.xY.VY(1).dim(2:3)]);
        if any(maskthis)
            y=reshape(y,size(SPM.xY.VY,1),size(SPM.xY.VY,2),SPM.xY.VY(1).dim(2),SPM.xY.VY(1).dim(3));
            if ~donefirst
                donefirst=true;
                [results_h,results_F,nill,SPM.xX_multivariate.dof,SPM.xX_multivariate.statsname]=conn_glm(SPM.xX_multivariate.X,y(:,:,maskthis),SPM.xX_multivariate.C,SPM.xX_multivariate.M);
                SPM.xX_multivariate.h=zeros([size(results_h,1),size(results_h,2),SPM.xY.VY(1).dim(1:3)]);
                SPM.xX_multivariate.F=zeros([size(results_F,1),size(results_F,2),SPM.xY.VY(1).dim(1:3)]);
            else
                [results_h,results_F]=conn_glm(SPM.xX_multivariate.X,y(:,:,maskthis),SPM.xX_multivariate.C,SPM.xX_multivariate.M);
            end
            SPM.xX_multivariate.h(:,:,n2,maskthis)=results_h;
            SPM.xX_multivariate.F(:,:,n2,maskthis)=results_F;
        end
    end
    if ~donefirst, error('Please check your data: There are no inmask voxels'); end
end
if issurface||ismatrix % surface- or matrix- based analyses
    V=struct('mat',SPM.xY.VY(1).mat,'dim',SPM.xY.VY(1).dim,'fname','mask.nii','pinfo',[1;0;0],'n',[1,1],'dt',[spm_type('uint8') spm_platform('bigend')]);
    spm_write_vol(V,double(mask));
    save('SPM.mat','SPM');
    conn_disp('fprintf','\nSecond-level results saved in folder %s\n',pwd);
    if ~nargout, conn_display('SPM.mat'); end
elseif ismember(secondlevelanalyses,[1 2]) % volume-based parametric stats
    save('SPM.mat','SPM');
    spm('Defaults','fmri');
    spm_unlink('mask.img','mask.hdr','mask.nii');
    if ~isempty(maskfile), 
        if isempty(mask)
            mask=ones(SPM.xY.VY(1).dim(1:3));
            [gridx,gridy]=ndgrid(1:SPM.xY.VY(1).dim(2),1:SPM.xY.VY(1).dim(3));
            xyz0=[gridx(:),gridy(:)]';
            for n2=1:SPM.xY.VY(1).dim(1)
                xyz=[n2+zeros(1,size(xyz0,2)); xyz0; ones(1,size(xyz0,2))];
                y=spm_get_data(SPM.xY.VY(:)',xyz);
                maskthis=~any(isnan(y),1)&any(diff(y,1,1)~=0,1);
                if ~isempty(maskfile), maskthis=maskthis&spm_get_data(vmaskfile,pinv(vmaskfile.mat)*SPM.xY.VY(1).mat*xyz)>0; end
                mask(n2,:,:)=reshape(maskthis,[1 SPM.xY.VY(1).dim(2:3)]);
            end
        end
        V=struct('mat',SPM.xY.VY(1).mat,'dim',SPM.xY.VY(1).dim,'fname','originalmask.nii','pinfo',[1;0;0],'n',[1,1],'dt',[spm_type('uint8') spm_platform('bigend')]);
        V=spm_write_vol(V,double(mask));
        SPM.xM=struct('T',[],'TH',-Inf(size(SPM.xX.X,1),1),'I',0,'VM',{V},'xs',struct('Masking','analysis threshold')); 
    end
    SPM=spm_spm(SPM);
    c=kron(C2,C1); 
    cname='connectivity result';
    if size(c,1)==1, Statname='T'; else Statname='F'; end
    if ~isfield(SPM.xX,'xKXs'), error('SPM analyses did not finish correctly'); end
    SPM.xCon = spm_FcUtil('Set',cname,Statname,'c',c',SPM.xX.xKXs);
    if isfield(SPM,'altestsmooth')&&SPM.altestsmooth, % modified smoothness estimation
        SPM=conn_est_smoothness(SPM);
        save('SPM.mat','SPM');
    end
    SPM=spm_contrasts(SPM,1:length(SPM.xCon));
    SPM.xY.VY=SPM.xY.VY(:);
    SPM.xsDes='';
    save('SPM.mat','SPM');
    conn_disp('fprintf','Second-level results saved in folder %s\n',pwd);
    if ~nargout,
        conn_display('SPM.mat',1);
    end
elseif ismember(secondlevelanalyses,[1 3]) % volume-based nonparametric stats
    V=struct('mat',SPM.xY.VY(1).mat,'dim',SPM.xY.VY(1).dim,'fname','mask.nii','pinfo',[1;0;0],'n',[1,1],'dt',[spm_type('uint8') spm_platform('bigend')]);
    spm_write_vol(V,double(mask));
    save('SPM.mat','SPM');
    conn_disp('fprintf','Second-level results saved in folder %s\n',pwd);
    if ~nargout,
        conn_display('SPM.mat',1);
    end
end
cd(pwd0);


end

