function score=conn_qascores(option, folderout, subjects, QC_variables, Control_variables, SubjectExclusionCriteria)
% internal function
% conn_qascore('DataValidity', folderout, subjects)
% conn_qascore('DataQuality', folderout, subjects, QC_variables, Control_variables)
% conn_qascore('DataSensitivity', folderout, subjects, QC_variables)

global CONN_x;

score=[];
ONEPERDAY=true; % note: QC score plots are saved in a different directory each day (to keep records of previous computations while saving storage space). Set ONEPERDAY=false if you prefer QC score plots to be saved in a different directory each time they are computed
switch(lower(option)) % Data Validity score
    case 'datavalidity',
        if nargin<2||isempty(folderout)
            if ONEPERDAY, folderout=fullfile(CONN_x.folders.qa,['QA_GUIrequest_DataValidity_',datestr(now,'yyyy_mm_dd')]);
            else folderout=fullfile(CONN_x.folders.qa,['QA_GUIrequest_DataValidity_',datestr(now,'yyyy_mm_dd_HHMMSSFFF')]);
            end
        end
        if nargin<3||isempty(subjects), subjects=1:CONN_x.Setup.nsubjects; end
        try
            conn_process('qaplots',conn_server('util_localfile',folderout),{'QA_DENOISE histogram'},subjects);
            %conn_batch(...
            %    'QA.plots',         {'QA_DENOISE histogram'},...
            %    'QA.foldername',    folderout,...
            %    'subjects',subjects);
            QC_PeakFC=conn_module('get','l2covariates','QC_PeakFC');
            QC_IqrFC=conn_module('get','l2covariates','QC_IqrFC');
            QC_PeakFC=QC_PeakFC(subjects);
            QC_IqrFC=QC_IqrFC(subjects); 
            k=1.348980; %spm_invNcdf(.75)-spm_invNcdf(.25);
            DataValidityScore=exp(-k*abs(mean(QC_PeakFC)/mean(QC_IqrFC)));
            conn_savematfile(fullfile(folderout,'DataValidityScore.mat'),'subjects','QC_PeakFC','QC_IqrFC','DataValidityScore'),
            CONN_x.Preproc.qa.folders{1}=folderout;
            CONN_x.Preproc.qa.inputs{1}={subjects};
            CONN_x.Preproc.qa.options1=rmfield(CONN_x.Preproc,{'variables','qa'});
            CONN_x.Preproc.qa.DataValidityScore=DataValidityScore;
            score=DataValidityScore;
        end

    case 'dataquality'
        if nargin<2||isempty(folderout)
            if ONEPERDAY, folderout=fullfile(CONN_x.folders.qa,['QA_GUIrequest_DataQuality_',datestr(now,'yyyy_mm_dd')]);
            else folderout=fullfile(CONN_x.folders.qa,['QA_GUIrequest_DataQuality_',datestr(now,'yyyy_mm_dd_HHMMSSFFF')]);
            end
        end
        if nargin<3||isempty(subjects), subjects=1:CONN_x.Setup.nsubjects; end
        if nargin<4||isempty(QC_variables),
            [nill,QC_variables]=conn_module('get','l2covariates','^QC_');
            QC_variables=QC_variables(cellfun('length',regexp(QC_variables,'^(QC_InvalidScans|QC_ProportionValidScans|QC_MeanMotion)$'))>0);
        end
        if nargin<5||isempty(Control_variables), Control_variables={}; end
        if isempty(QC_variables)
            conn_disp('Unable to compute Data Quality score. Expected the following QC_ 2nd-level covariates: QC_InvalidScans|QC_ProportionValidScans|QC_MeanMotion');
        else
            try
                conn_process('qaplots',conn_server('util_localfile',folderout),{'QA_DENOISE FC-QC'},subjects,[],[],QC_variables,[],[],Control_variables);
                %conn_batch(...
                %    'QA.plots',         {'QA_DENOISE FC-QC'},...
                %    'QA.l2covariates',  QC_variables,...
                %    'QA.foldername',    folderout,...
                %    'subjects',subjects);
                DataQualityScore=cellfun(@(name)conn_jsonread(fullfile(folderout,sprintf('QA_DENOISE_QC-FC.measure%s.json',name)),'AfterDenoising_PercentOverlap'),QC_variables);
                PercentSignificant=cell2mat(cellfun(@(name)conn_jsonread(fullfile(folderout,sprintf('QA_DENOISE_QC-FC.measure%s.json',name)),'AfterDenoising_PercentSignificantEdges'),QC_variables,'uni',0));
                conn_savematfile(fullfile(folderout,'DataQualityScore.mat'),'subjects','QC_variables','Control_variables','DataQualityScore','PercentSignificant'),
                CONN_x.Preproc.qa.folders{2}=folderout;
                CONN_x.Preproc.qa.inputs{2}={subjects, QC_variables, Control_variables};
                CONN_x.Preproc.qa.options2=rmfield(CONN_x.Preproc,{'variables','qa'});
                CONN_x.Preproc.qa.DataQualityScore=DataQualityScore;
                score=min(DataQualityScore);
            end
        end

    case 'datasensitivity'
        if nargin<2||isempty(folderout)
            if ONEPERDAY, folderout=fullfile(CONN_x.folders.qa,['QA_GUIrequest_DataSensitivity_',datestr(now,'yyyy_mm_dd')]);
            else folderout=fullfile(CONN_x.folders.qa,['QA_GUIrequest_DataSensitivity_',datestr(now,'yyyy_mm_dd_HHMMSSFFF')]);
            end
        end
        if nargin<3||isempty(subjects), subjects=1:CONN_x.Setup.nsubjects; end
        if nargin<4||isempty(QC_variables)
            [nill,QC_variables]=conn_module('get','l2covariates','^QC_');
            QC_variables=QC_variables(cellfun('length',regexp(QC_variables,'^(QC_ProportionValidScans|QC_MeanMotion|QC_MeanGSchange|QC_NORM_struct|QC_DOF|QC_PeakFC|QC_StdFC)$'))>0);
        end
        if isempty(QC_variables)
            conn_disp('Unable to compute Data Sensitivity. Expected the following QC_ 2nd-level covariates: QC_ProportionValidScans|QC_MeanMotion|QC_MeanGSchange|QC_NORM_struct|QC_DOF|QC_PeakFC|QC_StdFC');
        else
            if nargin<5, SubjectExclusionCriteria=[]; end
            try
                conn_process('qaplots',conn_server('util_localfile',folderout),{'QA_COV'},subjects,[],[],QC_variables);
                %conn_batch(...
                %    'QA.plots',         {'QA_COV'},...
                %    'QA.l2covariates',  QC_variables,...
                %    'QA.foldername',    folderout,...
                %    'subjects',subjects);
                QC_OutlierScore=conn_module('get','l2covariates','QC_OutlierScore');
                if isequal(SubjectExclusionCriteria,'mild') % note: default behavior is to remove severe outliers (QC_OutlierScore>3)
                    QC_ValidSubjects=QC_OutlierScore<=1.5;
                    conn_importl2covariate({'QC_ValidSubjects','QC_OutlierSubjects','ExcludeOutlierSubjects'},{QC_ValidSubjects,~QC_ValidSubjects,0./~QC_ValidSubjects},0);
                elseif ~isempty(SubjectExclusionCriteria) % changes QC_ValidSubjects, QC_OutlierSubjects, and ExcludeOutlierSubjects to new participant exclusion criteria (fixed number of excluded subjects, instead of default OutlierScore>3 threshold)
                    [nill,idx]=sort(QC_OutlierScore);
                    QC_ValidSubjects=~isnan(QC_OutlierScore);
                    QC_ValidSubjects(idx(numel(idx)-SubjectExclusionCriteria+1:end))=false;
                    conn_importl2covariate({'QC_ValidSubjects','QC_OutlierSubjects','ExcludeOutlierSubjects'},{QC_ValidSubjects,~QC_ValidSubjects,0./~QC_ValidSubjects},0);
                end
                QC_ValidSubjects=conn_module('get','l2covariates','QC_ValidSubjects');
                QC_DOF=conn_module('get','l2covariates','QC_DOF_WelchSatterthwaite');
                if isempty(QC_DOF), QC_DOF=conn_module('get','l2covariates','QC_DOF'); end
                Nsubjects=sum(QC_ValidSubjects>0);
                DataSensitivityScore=1-spm_Ncdf(1.645-0.1003*sqrt(sum(QC_DOF(QC_ValidSubjects>0)-3))); % power in fixed-effects analysis to detect an average r=0.1 effect or higher
                DataSensitivityScore_singlesubject=1-spm_Ncdf(1.645-0.1003*sqrt(mean(QC_DOF(QC_ValidSubjects>0)-3))); % power in single-subject fixed-effects analysis
                conn_savematfile(fullfile(folderout,'DataSensitivityScore.mat'),'subjects','QC_variables','DataSensitivityScore','DataSensitivityScore_singlesubject','Nsubjects','QC_OutlierScore','QC_ValidSubjects','QC_DOF'),
                CONN_x.Preproc.qa.folders{3}=folderout;
                CONN_x.Preproc.qa.inputs{3}={subjects, QC_variables};
                CONN_x.Preproc.qa.options3=rmfield(CONN_x.Preproc,{'variables','qa'});
                CONN_x.Preproc.qa.DataSensitivityScore=DataSensitivityScore;
                score=DataSensitivityScore;
            end
        end

    case 'combined', %{'datavaliditydataquality','dataqualitydatavalidity'}
        if nargin<2||isempty(folderout)
            if ONEPERDAY, folderout=fullfile(CONN_x.folders.qa,['QA_GUIrequest_combined_',datestr(now,'yyyy_mm_dd')]);
            else folderout=fullfile(CONN_x.folders.qa,['QA_GUIrequest_combined_',datestr(now,'yyyy_mm_dd_HHMMSSFFF')]);
            end
        end
        if nargin<3||isempty(subjects), subjects=1:CONN_x.Setup.nsubjects; end
        if nargin<4||isempty(QC_variables),
            [nill,QC_variables]=conn_module('get','l2covariates','^QC_');
            QC_variables=QC_variables(cellfun('length',regexp(QC_variables,'^(QC_InvalidScans|QC_ProportionValidScans|QC_MeanMotion)$'))>0);
        end
        if nargin<5||isempty(Control_variables), Control_variables={}; end
        if isempty(QC_variables)
            conn_disp('Unable to compute Data Quality score. Expected the following QC_ 2nd-level covariates: QC_InvalidScans|QC_ProportionValidScans|QC_MeanMotion');
        else
            try
                conn_process('qaplots',conn_server('util_localfile',folderout),{'QA_DENOISE histogram','QA_DENOISE FC-QC'},subjects,[],[],QC_variables,[],[],Control_variables);

                QC_PeakFC=conn_module('get','l2covariates','QC_PeakFC');
                QC_IqrFC=conn_module('get','l2covariates','QC_IqrFC');
                QC_PeakFC=QC_PeakFC(subjects);
                QC_IqrFC=QC_IqrFC(subjects);
                k=1.348980; %spm_invNcdf(.75)-spm_invNcdf(.25);
                DataValidityScore=exp(-k*abs(mean(QC_PeakFC)/mean(QC_IqrFC)));
                conn_savematfile(fullfile(folderout,'DataValidityScore.mat'),'subjects','QC_PeakFC','QC_IqrFC','DataValidityScore'),
                CONN_x.Preproc.qa.folders{1}=folderout;
                CONN_x.Preproc.qa.inputs{1}={subjects};
                CONN_x.Preproc.qa.options1=rmfield(CONN_x.Preproc,{'variables','qa'});
                CONN_x.Preproc.qa.DataValidityScore=DataValidityScore;

                DataQualityScore=cellfun(@(name)conn_jsonread(fullfile(folderout,sprintf('QA_DENOISE_QC-FC.measure%s.json',name)),'AfterDenoising_PercentOverlap'),QC_variables);
                PercentSignificant=cell2mat(cellfun(@(name)conn_jsonread(fullfile(folderout,sprintf('QA_DENOISE_QC-FC.measure%s.json',name)),'AfterDenoising_PercentSignificantEdges'),QC_variables,'uni',0));
                conn_savematfile(fullfile(folderout,'DataQualityScore.mat'),'subjects','QC_variables','Control_variables','DataQualityScore','PercentSignificant'),
                CONN_x.Preproc.qa.folders{2}=folderout;
                CONN_x.Preproc.qa.inputs{2}={subjects, QC_variables, Control_variables};
                CONN_x.Preproc.qa.options2=rmfield(CONN_x.Preproc,{'variables','qa'});
                CONN_x.Preproc.qa.DataQualityScore=DataQualityScore;
                
                score=[DataValidityScore min(DataQualityScore)];
            end
        end
end
end