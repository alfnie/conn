function conn_removesession(nsess,nsubs)
global CONN_x;

if nargin<2||isempty(nsubs), nsubs=1:CONN_x.Setup.nsubjects; end
for nsub=nsubs(:)'
    nses=setdiff(1:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsub)),nsess);
    CONN_x.Setup.functional{nsub}=CONN_x.Setup.functional{nsub}(nses);
    CONN_x.Setup.structural{nsub}=CONN_x.Setup.structural{nsub}(nses);
    try, for nalt=1:numel(CONN_x.Setup.secondarydataset), if isfield(CONN_x.Setup.secondarydataset(nalt),'functionals_explicit'), CONN_x.Setup.secondarydataset(nalt).functionals_explicit{nsub}=CONN_x.Setup.secondarydataset(nalt).functionals_explicit{nsub}(nses); end; end
    end
    try, CONN_x.Setup.unwarp_functional{nsub}=CONN_x.Setup.unwarp_functional{nsub}(nses);
    end
    try, CONN_x.Setup.coregsource_functional{nsub}=CONN_x.Setup.coregsource_functional{nsub}(nses);
    end
    %CONN_x.Setup.spm;
    %CONN_x.Setup.dicom;
    try, for nroi=1:numel(CONN_x.Setup.rois.files{nsub}), CONN_x.Setup.rois.files{nsub}{nroi}=CONN_x.Setup.rois.files{nsub}{nroi}(nses); end
    end
    try, for ncon=1:numel(CONN_x.Setup.conditions.values{nsub}), CONN_x.Setup.conditions.values{nsub}{ncon}=CONN_x.Setup.conditions.values{nsub}{ncon}(nses); end
    end
    try, for ncov=1:numel(CONN_x.Setup.l1covariates.files{nsub}), CONN_x.Setup.l1covariates.files{nsub}{ncov}=CONN_x.Setup.l1covariates.files{nsub}{ncov}(nses); end
    end
    %CONN_x.Setup.l2covariates.values;
    CONN_x.Setup.nscans{nsub}=CONN_x.Setup.nscans{nsub}(nses);
    CONN_x.Setup.nsessions(nsub)=numel(nses);
    %CONN_x.Setup.RT;
end