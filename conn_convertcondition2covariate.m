function out=conn_convertcondition2covariate(varargin)
global CONN_x;

if nargin<1||isempty(varargin{1}), return; end
nconditions=varargin;
donotremove=false;
donotapply=false;
usecovariateweighting=false;
subjects=1:CONN_x.Setup.nsubjects;
if iscell(nconditions)
    if ~isempty(nconditions)&&isequal(nconditions{1},'-DONOTAPPLYSUBJECTS')
        donotapply=true;
        donotremove=true;
        subjects=nconditions{2};
        nconditions=nconditions(3:end);
        usecovariateweighting=true;
    elseif ~isempty(nconditions)&&isequal(nconditions{1},'-DONOTAPPLY')
        donotapply=true;
        donotremove=true;
        nconditions=nconditions(2:end);
    elseif ~isempty(nconditions)&&isequal(nconditions{1},'-DONOTREMOVE')
        donotremove=true;
        nconditions=nconditions(2:end);
    end
    if ischar(nconditions{1})
        [ok,nconditions]=ismember(nconditions,CONN_x.Setup.conditions.names(1:end-1)); 
    else
        nconditions=cell2mat(nconditions);
    end
end
nconditions0=length(CONN_x.Setup.conditions.names);
nconditions=setdiff(nconditions,[0 nconditions0]);
out={};
RT={};
for ncondition=nconditions(:)'
    names=['Effect of ',CONN_x.Setup.conditions.names{ncondition}];
    nl1covariates=numel(CONN_x.Setup.l1covariates.names)-1;
    nl1covariate=nl1covariates+1;
    
    % adds first-level covariate
    for nsub=subjects(:)',
        nsess=CONN_x.Setup.nsessions(min(length(CONN_x.Setup.nsessions),nsub));
        for nses=1:nsess,
            if numel(RT)<nsub||numel(RT{nsub})<nses, RT{nsub}{nses}=conn_get_rt(nsub,nses); end
            rt=RT{nsub}{nses}/10;
            hrf=spm_hrf(rt);
            offs=ceil(100/rt);
            onset=CONN_x.Setup.conditions.values{nsub}{ncondition}{nses}{1};
            durat=CONN_x.Setup.conditions.values{nsub}{ncondition}{nses}{2};
            if 0,%numel(CONN_x.Setup.conditions.values{nsub}{ncondition}{nses})>2, val=CONN_x.Setup.conditions.values{nsub}{ncondition}{nses}{3}; % placeholder for block/event weighting
            else val=ones(size(onset)); 
            end
            if donotapply||(numel(CONN_x.Setup.nscans)>=nsub&&numel(CONN_x.Setup.nscans{nsub})>=nses)
                nscans=0;
                try, nscans=CONN_x.Setup.nscans{nsub}{nses}; end
                if donotapply&&nscans==0, 
                    nscans=100;
                    if ~isempty(onset), nscans=max(nscans,ceil(max(onset)/RT{nsub}{nses})+10); end
                    if ncondition==nconditions(1), conn_disp('fprintf','warning: undefined number of scans for subject %d session %d. Assumming %d scans. Please enter functional data first to avoid this warning\n',nsub,nses,nscans); end
                end
                x=zeros(offs+ceil(nscans*RT{nsub}{nses}/rt),1);
                if length(durat)>=1, 
                    for n1=1:length(onset), 
                        tdurat=max(rt,min(offs*rt+RT{nsub}{nses}*nscans-onset(n1),durat(min(length(durat),n1))));
                        in=offs+round(1+onset(n1)/rt+(0:tdurat/rt-1));
                        x(in(in>0))=val(n1); 
                    end
                end
                if CONN_x.Setup.acquisitiontype==1,
                    x=convn(x,hrf);
                end
                x=mean(reshape(x(offs+(1:10*nscans)),[10,nscans]),1)';%x=x(1+10*(0:nscans-1));
                if usecovariateweighting&&CONN_x.Setup.conditions.param(ncondition)>0
                    %try
                        covfilename=CONN_x.Setup.l1covariates.files{nsub}{CONN_x.Setup.conditions.param(ncondition)}{nses}{1};
                        switch(covfilename),
                            case '[raw values]',
                                covdata=CONN_x.Setup.l1covariates.files{nsub}{CONN_x.Setup.conditions.param(ncondition)}{nses}{3};
                            otherwise,
                                covdata=conn_loadtextfile(covfilename,false);
                        end
                        x=x.*max(0,sum(covdata,2));
                    %end
                end
                if donotapply, out{nsub}{ncondition}{nses}=x;
                else CONN_x.Setup.l1covariates.files{nsub}{nl1covariate}{nses}={'[raw values]',{sprintf('size [%s]',num2str(size(x))),'[raw values]'},x};
                end
            end
        end
    end
    if ~donotapply
        CONN_x.Setup.l1covariates.names{nl1covariate}=names;
        CONN_x.Setup.l1covariates.names{nl1covariate+1}=' ';
    end
end

% removes conditions
if ~donotremove
    nconditions=setdiff(1:nconditions0,nconditions);
    CONN_x.Setup.conditions.names=CONN_x.Setup.conditions.names(nconditions);
    nconditions=setdiff(nconditions,nconditions0);
    for n1=1:length(CONN_x.Setup.conditions.values), CONN_x.Setup.conditions.values{n1}={CONN_x.Setup.conditions.values{n1}{nconditions}}; end
    CONN_x.Setup.conditions.model=CONN_x.Setup.conditions.model(nconditions);
    CONN_x.Setup.conditions.param=CONN_x.Setup.conditions.param(nconditions);
    CONN_x.Setup.conditions.filter=CONN_x.Setup.conditions.filter(nconditions);
end

