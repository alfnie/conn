% creates QC-FC plots separately within each site (for a multi-site study)

foldername = 'validsubjects';                        % output QA-report name
groupvariable = 'QC_validsubjects'; %'AllSubjects';                  % 2nd-level covariate name containing "valid subjects" group (e.g. 'AllSubjects')
sitevariable = {'MSU','NFB','PTSP','RHY','STUT'};%{'site_1','site_2','site_3'};    % 2nd-level covariate names containing SITES subject-groups
% conn load /myfolder/myproject.mat;            % uncomment this line to load a different conn project (by default it works on the currently-open project)  

QCvariables = {'QC_InvalidScans','QC_MeanMotion','QC_ProportionValidScans'};    % name of QC variables to use for QC-FC plots
DOCOMBI = true;                   % computes QC-FC plots combining all sites
DOSITES = true;                   % computes QC-FC plots separately for each site


[filepath,filename,fileext] = fileparts(conn_projectmanager('projectfile')); % project folder
filepath = conn_fullfile(filepath,filename);
[covs,covs_names]=conn_module('get','l2covariates');
idx=find(ismember(covs_names,groupvariable)); valid=covs(:,idx);
if DOSITES, sites=[]; for k=1:numel(sitevariable), idx=find(ismember(covs_names,sitevariable{k})); sites=[sites, covs(:,idx)]; end; end

if DOCOMBI
    conn_batch(...
        'QA.plots',         {'QA_DENOISE FC-QC'},...
        'QA.l2covariates',  QCvariables,...
        'subjects',         find(valid),...
        'QA.foldername',    fullfile(filepath,'results','qa',sprintf('QA_%s_all',foldername)) );
end

if DOSITES
    for k=1:numel(sitevariable),
        conn_batch(...
            'QA.plots',         {'QA_DENOISE FC-QC'},...
            'QA.l2covariates',  QCvariables,...
            'subjects',         find(sites(:,k)&valid),...
            'QA.foldername',    fullfile(filepath,'results','qa',sprintf('QA_%s_site%d',foldername,k)) );
    end
end

