function out=conn_importspm(spmfiles,varargin)
% Imports experiment info from SPM.mat files
%
% CONN_IMPORTSPM(spmfiles);
%
% CONN_IMPORTSPM(spmfiles,parameter1_name,parameter1_value,...);
%  Valid parameter names are:
%    'addfunctional'             : 1/0 add functinal data (from SPM.xY) to Setup.Functional [true]
%    'addconditions'             : 1/0 add condition information (from SPM.Sess(nses).U) to Setup.Conditions [true]
%    'breakconditionsbysession'  : 1/0 creates session-specific conditions [false] 
%    'addrestcondition'          : 1/0 add default 'rest' condition (entire functional data) to Setup.Conditions [true]
%    'keeppreviousconditions'    : 1/0 keeps previously-defined conditions [true]
%    'addcovariates'             : 1/0 add additional covariates (from SPM.Sess(nses).C) to Setup.CovariatesFirstLevel [true]
%    'addrealignment'            : 1/0 add realignment covariates (from rp_<functional_filename>.txt) to Setup.CovariatesFirstLevel [false]
%    'addartfiles'               : 1/0 add ART covariates (from art_regression_outliers_<functional_filename>.mat) to Setup.CovariatesFirstLevel [true]
%    

global CONN_x;
options=struct('addfunctional',true,...
               'addconditions',true,...
               'addrestcondition',true,...
               'addcovariates',true,...
               'addrealignment',false,...
               'addartfiles',true,...
               'breakconditionsbysession',false,...
               'keeppreviousconditions',true,...
               'localcopy',false,...
               'nset',0,...
               'bidsname','',...
               'subjects',[]);
for n1=1:2:nargin-2, if ~isfield(options,lower(varargin{n1})), error('unknown option %s',lower(varargin{n1})); else options.(lower(varargin{n1}))=varargin{n1+1}; end; end
if ~isempty(options.subjects), SUBJECTS=options.subjects; else SUBJECTS=1:CONN_x.Setup.nsubjects; end
if nargin>0&&~isempty(spmfiles),
    if ~iscell(spmfiles),spmfiles=cellstr(spmfiles);end
    nsubjects=length(spmfiles);
    for nsub=1:nsubjects,
        CONN_x.Setup.spm{SUBJECTS(nsub)}=conn_file(char(spmfiles{nsub}));
    end
end
if isempty(options.bidsname)
    if ~options.nset, options.bidsname='func';
    else options.bidsname=sprintf('dataset%dfunc',options.nset);
    end
end

%CONN_x.Setup.nsessions=zeros(1,CONN_x.Setup.nsubjects);
changed=0;err=0;
importfunctional=options.addfunctional;
if ~options.keeppreviousconditions, CONN_x.Setup.conditions.names={' '}; end
for nsub=SUBJECTS,
    if ~isempty(CONN_x.Setup.spm{nsub}{1}),
        session_count=0;
        files=cellstr(CONN_x.Setup.spm{nsub}{1});
        for ifile=1:numel(files)
            %try,
            spmfile=conn_loadmatfile(files{ifile});
            if options.addconditions
                if isfield(spmfile.SPM,'xBF')&&isfield(spmfile.SPM.xBF,'UNITS'),
                    units=spmfile.SPM.xBF.UNITS;
                    if strcmp(units,'scans'), units=1;
                    elseif strcmp(units,'secs'), units=2;
                    else, conn_disp(['ERROR subject ',num2str(nsub),': invalid SPM.xBF.UNITS value (assuming scans)',' in file ',files{ifile}]); err=err+1; end
                else, conn_disp(['ERROR subject ',num2str(nsub),': SPM.xBF.UNITS not found (assuming scans)',' in file ',files{ifile}]); units=1; err=err+1; end
            end
            if importfunctional||(options.addconditions&&units~=2)
                if isfield(spmfile.SPM,'xY')&&isfield(spmfile.SPM.xY,'RT'),
                    if session_count>0&&CONN_x.Setup.RT(nsub)~=spmfile.SPM.xY.RT,
                        conn_disp(['warning subject ',num2str(nsub),': SPM.xY.RT different from previous session(s) (leaving unchanged)',' in file ',files{ifile}]);
                    else
                        CONN_x.Setup.RT(nsub)=spmfile.SPM.xY.RT;
                    end
                else,
                    if isempty(CONN_x.Setup.RT), CONN_x.Setup.RT=nan; end
                    CONN_x.Setup.RT(nsub)=CONN_x.Setup.RT(end);
                    conn_disp(['warning subject ',num2str(nsub),': SPM.xY.RT not found (leaving unchanged)',' in file ',files{ifile}]);
                end
            end
            if ~isfield(spmfile.SPM,'nscan')&&~isfield(spmfile.SPM,'Sess'), conn_disp(['warning subject ',num2str(nsub),': SPM.nscan / SPM.Sess not found (treating as single-session design)',' in file ',files{ifile}]); nsess=1;
            elseif isfield(spmfile.SPM,'nscan'), nsess=length(spmfile.SPM.nscan);
            else nsess=length(spmfile.SPM.Sess);
            end
            for nses=1:nsess,
                displayedtext=false;
                session_count=session_count+1;
                if nsess>1, idxscans=spmfile.SPM.Sess(nses).row;
                else idxscans=1:size(spmfile.SPM.xX.X,1); 
                end
                nscans=numel(idxscans);
                if ~isfield(spmfile.SPM.xY,'P')&&isfield(spmfile.SPM.xY,'Y'), 
                    temp=SPM.xY.Y;
                    if numel(temp)==1, temp=cellstr(conn_expandframe(temp{1})); end
                    spmfile.SPM.xY.P=char(temp);
                    spmfile.SPM.xY.VY=spm_data_hdr_read(SPM.xY.P);
                end
                if ~isfield(spmfile.SPM.xY,'P')&&isfield(spmfile.SPM.xY,'VY'), 
                    temp={spmfile.SPM.xY.VY.fname};
                    if numel(temp)==1, temp=cellstr(conn_expandframe(temp{1})); end
                    spmfile.SPM.xY.P=char(temp);
                end
                if isfield(spmfile.SPM.xY,'P')
                    temp=spmfile.SPM.xY.P;
                    if size(spmfile.SPM.xY.P,1)~=size(spmfile.SPM.xX.X,1), temp=char(cellfun(@conn_expandframe,cellstr(temp),'uni',0)); end
                    assert(size(temp,1)==size(spmfile.SPM.xX.X,1),'unexpected number of input volumes in SPM.xY structure (xY.P=%d, xX.X=%d)',size(temp,1),size(spmfile.SPM.xX.X,1));
                    filename=fliplr(deblank(fliplr(deblank(temp(idxscans,:)))));
                    switch(filesep),case '\',idx=find(filename=='/');case '/',idx=find(filename=='\');end; filename(idx)=filesep;
                else
                    if nses==1, conn_disp(['warning subject ',num2str(nsub),': SPM.xY.P not found ',' in file ',files{ifile}]); end
                    filename=[];
                end
                if importfunctional
                    lastfile='';
                    n=0; while n<size(filename,1),
                        ok=conn_existfile(filename(n+1,:));
                        if ok, n=n+1;
                        else
                            if changed,
                                fullnamematch=strvcat(fliplr(fullname1),fliplr(fullname2));
                                m=sum(cumsum(fullnamematch(1,:)~=fullnamematch(2,:))==0);
                                m1=max(0,length(fullname1)-m); m2=max(0,length(fullname2)-m);
                                %filename=strvcat(filename(1:n,:),[repmat(fullname2(1:m2),[size(filename,1)-n,1]),filename(n+1:end,m1+1:end)]);
                                newfile=[fullname2(1:m2),filename(n+1,m1+1:end)];
                                filenamet=strvcat(filename(1:n,:),newfile,filename(n+2:end,:));
                                if isequal(lastfile,regexprep(newfile,',\d+$','')), askthis=0;filename=filenamet;
                                elseif ~conn_existfile(filenamet(n+1,:)), askthis=1;
                                else
                                    try
                                        if ~displayedtext
                                            conn_disp(['conn_importspm: updating reference from ',deblank(filename(n+1,:)),' to ',deblank(filenamet(n+1,:))]);
                                            displayedtext=true;
                                        end                                          
                                        fileinfo=conn_file(deblank(filenamet(n+1,:)));
                                        %[V,str,icon]=conn_getinfo(deblank(filenamet(n+1,:)));
                                        askthis=0;filename=filenamet;
                                    catch
                                        askthis=1;
                                    end
                                end
                                lastfile=regexprep(newfile,',\d+$','');
                            else askthis=1;
                            end
                            if askthis,
                                conn_disp(['conn_importspm: file ',deblank(filename(n+1,:)),' not found']);
                                fullname1=deblank(filename(n+1,:));
                                [pathname1,name1,ext1,num1]=spm_fileparts(fullname1);
                                name2='';
                                while ~strcmp(name2,[name1,ext1])&&~isequal(name2,0)
                                    conn_disp(['File not found: ',name1,ext1]);
                                    [name2,pathname2]=conn_fileutils('uigetfile',['*',ext1],['File not found: ',name1,ext1],['*',name1,ext1]);
                                end
                                if isequal(name2,0), importfunctional=false; break; end
                                fullname2=fullfile(pathname2,[name2,num1]);
                                changed=1;
                                fullnamematch=strvcat(fliplr(fullname1),fliplr(fullname2));
                                m=sum(cumsum(fullnamematch(1,:)~=fullnamematch(2,:))==0);
                                m1=max(0,length(fullname1)-m); m2=max(0,length(fullname2)-m);
                                filename=strvcat(filename(1:n,:),[fullname2(1:m2),filename(n+1,m1+1:end)],filename(n+2:end,:));
                            end
                        end
                    end
                end
                if options.addfunctional||options.addcovariates||options.addrealignment||options.addartfiles
                    [filename1_path,filename1_name,filename1_ext,filename1_num]=spm_fileparts(filename(1,:));
                    filename1=fullfile(filename1_path,[filename1_name,filename1_ext]);
                end
                if importfunctional
                    if options.localcopy,
                        [nill,nill,nscans]=conn_importvol2bids(filename,nsub,session_count,options.bidsname);
                        if ~options.nset, CONN_x.Setup.nscans{nsub}{session_count}=nscans; end
                    else
                        nscans=conn_set_functional(nsub,session_count,options.nset,filename);
                    end
                    %[fileinfo,V]=conn_file(filename);
                    %CONN_x.Setup.functional{nsub}{session_count}=fileinfo;
                    %CONN_x.Setup.nscans{nsub}{session_count}=numel(V);
                else
                    if options.addfunctional, conn_disp(['ERROR subject ',num2str(nsub),' no functional data imported',' in file ',files{ifile}]); err=err+1; end
                    if options.addfunctional||options.addconditions||options.addcovariates||options.addrealignment||options.addartfiles, CONN_x.Setup.nscans{nsub}{session_count}=nscans; end
                end
                % adds rest condition
                if options.addconditions&&options.addrestcondition
                    name='rest';
                    idx=strmatch(name,CONN_x.Setup.conditions.names(1:end-1),'exact');
                    if isempty(idx), idx=length(CONN_x.Setup.conditions.names); CONN_x.Setup.conditions.names{end+1}=' '; end
                    CONN_x.Setup.conditions.model{idx}=[];
                    CONN_x.Setup.conditions.param(idx)=0;
                    CONN_x.Setup.conditions.filter{idx}=[];
                    CONN_x.Setup.conditions.names{idx}=name;
                    CONN_x.Setup.conditions.values{nsub}{idx}{session_count}{1}=0;
                    CONN_x.Setup.conditions.values{nsub}{idx}{session_count}{2}=inf;
                end
                % adds other conditions/covariates
                if 0, % spmascovariate option is no longer supported   isfield(CONN_x.Setup,'spmascovariate')&&CONN_x.Setup.spmascovariate
                    if isfield(spmfile.SPM.Sess(nses),'col')&&~isempty(spmfile.SPM.Sess(nses).col)
                        name='SPM effects';
                        idx=strmatch(name,CONN_x.Setup.l1covariates.names,'exact');
                        if isempty(idx), idx=length(CONN_x.Setup.l1covariates.names); CONN_x.Setup.l1covariates.names{end+1}=' '; end
                        CONN_x.Setup.l1covariates.names{idx}=name;
                        CONN_x.Setup.l1covariates.files{nsub}{idx}{session_count}={'[raw values]',[],spmfile.SPM.xX.X(spmfile.SPM.Sess(nses).row,spmfile.SPM.Sess(nses).col)};
                    end
                    if isfield(spmfile.SPM.Sess(nses),'C')&&~isempty(spmfile.SPM.Sess(nses).C.C)
                        name='SPM covariates';
                        idx=strmatch(name,CONN_x.Setup.l1covariates.names,'exact');
                        if isempty(idx), idx=length(CONN_x.Setup.l1covariates.names); CONN_x.Setup.l1covariates.names{end+1}=' '; end
                        CONN_x.Setup.l1covariates.names{idx}=name;
                        CONN_x.Setup.l1covariates.files{nsub}{idx}{session_count}={'[raw values]',[],spmfile.SPM.Sess(nses).C.C};
                    elseif conn_existfile(conn_prepend('rp_',filename1,'.txt')),
                        name='realignment';
                        idx=strmatch(name,CONN_x.Setup.l1covariates.names,'exact');
                        if isempty(idx), idx=length(CONN_x.Setup.l1covariates.names); CONN_x.Setup.l1covariates.names{end+1}=' '; end
                        CONN_x.Setup.l1covariates.names{idx}=name;
                        CONN_x.Setup.l1covariates.files{nsub}{idx}{session_count}=conn_file(conn_prepend('rp_',filename1,'.txt'));
                    end
                else
                    if options.addconditions,
                        if ~isfield(spmfile.SPM,'Sess')||~isfield(spmfile.SPM.Sess(nses),'U'), 
                            conn_disp(['ERROR subject ',num2str(nsub),': SPM.Sess.U not found (no condition information found/imported)',' in file ',files{ifile}]); err=err+1;
                        else
                            nconditions=length(spmfile.SPM.Sess(nses).U);
                            for ncondition=1:nconditions,
                                name=spmfile.SPM.Sess(nses).U(ncondition).name{1};
                                if isempty(name)||~ischar(name), name=sprintf('SPMcondition%d',ncondition); end
                                if options.breakconditionsbysession, name=sprintf('%s_Session%d',name,session_count); end
                                idx=strmatch(name,CONN_x.Setup.conditions.names,'exact');
                                if isempty(idx),
                                    idx=length(CONN_x.Setup.conditions.names);
                                    CONN_x.Setup.conditions.names{end+1}=' ';
                                end
                                CONN_x.Setup.conditions.model{idx}=[];
                                CONN_x.Setup.conditions.param(idx)=0;
                                CONN_x.Setup.conditions.filter{idx}=[];
                                CONN_x.Setup.conditions.names{idx}=name;
                                if units==2,
                                    CONN_x.Setup.conditions.values{nsub}{idx}{session_count}{1}=spmfile.SPM.Sess(nses).U(ncondition).ons;
                                    CONN_x.Setup.conditions.values{nsub}{idx}{session_count}{2}=spmfile.SPM.Sess(nses).U(ncondition).dur;
                                else,
                                    rt=conn_get_rt(nsub,nses,options.nset);
                                    CONN_x.Setup.conditions.values{nsub}{idx}{session_count}{1}=(spmfile.SPM.Sess(nses).U(ncondition).ons-0)*rt;
                                    CONN_x.Setup.conditions.values{nsub}{idx}{session_count}{2}=spmfile.SPM.Sess(nses).U(ncondition).dur*rt;
                                end
                            end
                        end
                    end
                    
                    if options.addcovariates
                        if isfield(spmfile.SPM,'Sess')&&isfield(spmfile.SPM.Sess(nses),'C')&&~isempty(spmfile.SPM.Sess(nses).C.C)
                            name='SPM covariates';
                            idx=strmatch(name,CONN_x.Setup.l1covariates.names,'exact');
                            if isempty(idx), idx=length(CONN_x.Setup.l1covariates.names); CONN_x.Setup.l1covariates.names{end+1}=' '; end
                            CONN_x.Setup.l1covariates.names{idx}=name;
                            CONN_x.Setup.l1covariates.files{nsub}{idx}{session_count}={'[raw values]',[],spmfile.SPM.Sess(nses).C.C};
                        else conn_disp(['warning subject ',num2str(nsub),': SPM.Sess.C not found (no SPM covariates imported)',' in file ',files{ifile}]);
                        end
                    end
                    if options.addrealignment
                        for remov=0:10,if conn_existfile(conn_prepend('rp_',conn_prepend(-remov,filename1),'.txt')); break; end; end
                        if remov<10
                            name='realignment';
                            idx=strmatch(name,CONN_x.Setup.l1covariates.names,'exact');
                            if isempty(idx), idx=length(CONN_x.Setup.l1covariates.names); CONN_x.Setup.l1covariates.names{end+1}=' '; end
                            CONN_x.Setup.l1covariates.names{idx}=name;
                            CONN_x.Setup.l1covariates.files{nsub}{idx}{session_count}=conn_file(conn_prepend('rp_',conn_prepend(-remov,filename1),'.txt'));
                        else conn_disp(['warning subject ',num2str(nsub),': file ',conn_prepend('rp_',filename1,'.txt'),' not found (no realignment covariates imported)',' in file ',files{ifile}]);
                        end
                    end
                    if options.addartfiles
                        for remov=0:10,if conn_existfile(conn_prepend('art_regression_outliers_',conn_prepend(-remov,filename1),'.mat')); break; end; end
                        if remov<10
                            name='ART covariates';
                            idx=strmatch(name,CONN_x.Setup.l1covariates.names,'exact');
                            if isempty(idx), idx=length(CONN_x.Setup.l1covariates.names); CONN_x.Setup.l1covariates.names{end+1}=' '; end
                            CONN_x.Setup.l1covariates.names{idx}=name;
                            CONN_x.Setup.l1covariates.files{nsub}{idx}{session_count}=conn_file(conn_prepend('art_regression_outliers_',conn_prepend(-remov,filename1),'.mat'));
                        else conn_disp(['warning subject ',num2str(nsub),': file ',conn_prepend('art_regression_outliers_',filename1,'.mat'),' not found (no ART covariates imported)',' in file ',files{ifile}]);
                        end
                    end
                end
            end
        end
        if options.addfunctional||options.addconditions||options.addcovariates||options.addrealignment||options.addartfiles, CONN_x.Setup.nsessions(nsub)=session_count; end
    end
end
if ~err&&options.addconditions, % fills possible empty conditions for each subject/session
    nconditions=length(CONN_x.Setup.conditions.names)-1;
    for nsub2=1:nsub,
        nsess=CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsub2));
        for ncondition=1:nconditions
            for nses=1:nsess
                if numel(CONN_x.Setup.conditions.values{nsub2})<ncondition||numel(CONN_x.Setup.conditions.values{nsub2}{ncondition})<nses||numel(CONN_x.Setup.conditions.values{nsub2}{ncondition}{nses})<1,
                    CONN_x.Setup.conditions.values{nsub2}{ncondition}{nses}{1}=[];
                end
                if numel(CONN_x.Setup.conditions.values{nsub2})<ncondition||numel(CONN_x.Setup.conditions.values{nsub2}{ncondition})<nses||numel(CONN_x.Setup.conditions.values{nsub2}{ncondition}{nses})<2,
                    CONN_x.Setup.conditions.values{nsub2}{ncondition}{nses}{2}=[];
                end
            end
        end
    end
end
	%catch,
	%	conn_disp(['warning: importing from ',CONN_x.Setup.spm{nsub}{1},': unexpected error (stopped importing for this subject)']); 
	%	err=err+1;
	%	conn_disp(lasterr);
	%end
if isfield(CONN_x,'gui')&&isnumeric(CONN_x.gui)&&CONN_x.gui, 
    if ~err, conn_msgbox([num2str(numel(SUBJECTS)),' subjects imported with no errors'],'Done',true);
    else conn_msgbox({['Import finished with ',num2str(err),' errors'],'see log or Matlab command window for details'},'WARNING!',true); 
    end
end
if nargout,out=CONN_x;end
end

