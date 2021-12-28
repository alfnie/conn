function [ok,matlabbatch,outputfiles,job_id]=conn_setup_preproc(STEPS,varargin)
% CONN_SETUP_PREPROC
% Runs individual preprocessing steps
%
% conn_setup_preproc(steps)
% runs preprocessing pipeline (default_*) or one/multiple individual preprocessing steps (structural_* and functional_*). Valid step names are (enter as cell array to run multiple sequential steps):
%
% conn_setup_preproc(steps,'param1_name',param1_value,'param2_name',param2_value,...)
% defines additional non-default values for parameters specific to individual steps
%
% conn_setup_preproc('steps')
% returns the full list of valid preprocessing-step names
%
% see "help conn_batch" for details about available preprocessing steps as well as a full list of additional preprocessing parameters
%


global CONN_x CONN_gui;
PREFERSPM8OVERSPM12=false; % set to true if you prefer to use SPM8 procedures over SPM12 ones (when SPM12 is installed)
ALLSETSPERMISSIONS=false;  % set to true if you want to allow dataset-1 or above preprocessing steps to import new ROIs, first-level covariates, or second-level covariates into your project
if isdeployed, spmver12=true;
else spmver12=str2double(regexp(spm('ver'),'SPM(\d+)','tokens','once'))>=12;
end
if isfield(CONN_gui,'font_offset'),font_offset=CONN_gui.font_offset; else font_offset=0; end
%if isfield(CONN_x,'pobj')&&isfield(CONN_x.pobj,'readonly')&&CONN_x.pobj.readonly, error('This procedure cannot be run while in view-only mode. Please re-load your project to enable edits'); end
if ~nargin, STEPS=''; varargin={'multiplesteps',1}; end
options=varargin;
steps={'default_mni','default_mnifield','default_mnidirectfield','default_ss','default_ssfield','default_ssnl',...
    'functional_surface_coreg&resample',...
    'structural_manualorient','structural_center','structural_segment',...
    'structural_normalize','structural_segment&normalize','structural_normalize_preservemasks',...
    'structural_manualspatialdef', ...
    'functional_removescans','functional_manualorient','functional_center','functional_centertostruct'...
    'functional_slicetime','functional_bandpass','functional_regression','functional_realign','functional_realign&unwarp',...
    'functional_realign&unwarp&fieldmap','functional_art','functional_coregister_affine_reslice',...
    'functional_segment',...
    'functional_manualspatialdef',...
    'functional_smooth','functional_motionmask',...
    'functional_segment&normalize_indirect','functional_normalize_indirect','functional_normalize_indirect_preservemasks', ...
    'functional_segment&normalize_direct','functional_normalize_direct', ...
    'functional_realign_noreslice', ...
    'functional_coregister_nonlinear', ...
    'functional_coregister_affine_noreslice', ...
    'functional_label', ...
    'functional_label_as_original', ...
    'functional_label_as_subjectspace', ...
    'functional_label_as_mnispace', ...
    'functional_label_as_surfacespace', ...
    'functional_label_as_smoothed', ...
    'functional_load', ...
    'functional_load_from_original', ...
    'functional_load_from_subjectspace', ...
    'functional_load_from_mnispace', ...
    'functional_load_from_surfacespace', ...
    'functional_load_from_smoothed', ...
    'functional_smooth_masked',...
    'functional_surface_resample', ...
    'functional_surface_smooth', ...
    'functional_vdm_apply', ...
    'functional_vdm_create' ...
    };
%'functional_normalize','functional_segment&normalize',...
steps_names={'<HTML><b>default preprocessing pipeline</b> for volume-based analyses (direct normalization to MNI-space)</HTML>','<HTML><b>preprocessing pipeline</b> for volume-based analyses (indirect normalization to MNI-space) when FieldMaps are available</HTML>','<HTML><b>preprocessing pipeline</b> for volume-based analyses (direct normalization to MNI-space) when FieldMaps are available</HTML>','<HTML><b>preprocessing pipeline</b> for surface-based analyses (in subject-space)</HTML>','<HTML><b>preprocessing pipeline</b> for surface-based analyses (in subject-space) when FieldMaps are available</HTML>','<HTML><b>preprocessing pipeline</b> for surface-based analyses (in subject-space) using nonlinear coregistration</HTML>',...
    'functional Direct Coregistration to structural without reslicing followed by Resampling of functional data at the location of FreeSurfer subject-specific structural cortical surface (converts volume- to surface- level data)', ...
    'structural Manual transformation (rotation/flip/translation/affine of structural volumes)','structural Center to (0,0,0) coordinates (translation)','structural Segmentation (Grey/White/CSF tissue estimation)',...
    'structural Normalization (MNI space normalization)','structural Segmentation & Normalization (simultaneous Grey/White/CSF segmentation and MNI normalization)','structural Normalization preserving Grey/White/CSF masks (MNI space normalization of structural, applying same transformation to existing Grey/White/CSF masks)',...
    'structural Manual deformation (non-linear transformation of structural volumes)', ...
    'functional Removal of initial scans (disregard initial functional scans)','functional Manual transformation (rotation/flip/translation/affine of functional volumes)','functional Center to (0,0,0) coordinates (translation)','functional Center to structural coordinates (translation)'...
    'functional Slice-Timing correction (STC; correction for inter-slice differences in acquisition time)','functional Band-pass filtering (temporal filtering of BOLD data)','functional Regression of temporal components (keep residuals of linear model to BOLD timeseries)','functional Realignment (subject motion estimation and correction)','functional Realignment & unwarp (subject motion estimation and correction)',...
    'functional Realignment & unwarp & susceptibility distortion correction (subject motion estimation and correction)','functional Outlier detection (ART-based identification of outlier scans for scrubbing)','functional Direct Coregistration to structural (rigid body transformation)',...
    'functional Segmentation (Grey/White/CSF segmentation)',...
    'functional Manual deformation (non-linear transformation of functional volumes)',...
    'functional Smoothing (spatial convolution with Gaussian kernel)','functional Motion-mask estimation (BOLD signal derivative wrt movement parameters)',...
    'functional Indirect Segmentation & Normalization (coregister functional/structural; structural segmentation & normalization; apply same deformation field to functional)', ...
    'functional Indirect Normalization (coregister functional/structural; structural normalization; apply same deformation field to functional)',...
    'functional Indirect Normalization preserving Grey/White/CSF masks (coregister functional/structural; structural normalization; apply same deformation field to functional and to existing Grey/White/CSF masks)',...
    'functional Direct Segmentation & Normalization (simultaneous Grey/White/CSF segmentation and MNI normalization)',...
    'functional Direct Normalization (MNI space normalization)', ...
    'functional Realignment without reslicing (subject motion estimation and correction)', ...
    'functional Indirect Coregistration to structural (non-linear transformation)', ...
    'functional Direct Coregistration to structural without reslicing (rigid body transformation)', ...
    'functional Label current functional files as new secondary dataset (custom label)', ...
    'functional Label current functional files as "original data"', ...
    'functional Label current functional files as "subject-space data"', ...
    'functional Label current functional files as "mni-space data"', ...
    'functional Label current functional files as "surface-space data"', ...
    'functional Label current functional files as "smoothed data"', ...
    'functional Load functional data from previously labeled dataset (custom label)', ...
    'functional Load functional data from "original data" dataset', ...
    'functional Load functional data from "subject-space data" dataset', ...
    'functional Load functional data from "mni-space data" dataset', ...
    'functional Load functional data from "surface-space data" dataset', ...
    'functional Load functional data from "smoothed data" dataset', ...
    'functional Masked Smoothing (spatial convolution with Gaussian kernel restricted to voxels within Grey Matter mask)', ...
    'functional Resampling of functional data at the location of FreeSurfer subject-specific structural cortical surface (converts volume- to surface- level data)', ...
    'functional Smoothing of surface-level functional data (spatial diffusion on surface tessellation)', ...
    'functional Susceptibility Distortion Correction using voxel-displacement maps (VDM)', ...
    'functional Creation of voxel-displacement map (VDM) for Susceptibility Distortion Correction' ...
    };
%'functional Normalization (MNI space normalization)','functional Segmentation & Normalization (simultaneous Grey/White/CSF segmentation and MNI normalization)',...
steps_descr={{'INPUT: structural&functional volumes','OUTPUT (all in MNI-space): skull-stripped normalized structural volume, Grey/White/CSF normalized masks, realigned slice-time corrected normalized smoothed functional volumes, subject movement ''realignment'' and ''scrubbing'' 1st-level covariate'},{'INPUT: structural&functional&FieldMap volumes',  'OUTPUT (all in MNI-space): skull-stripped normalized structural volume, Grey/White/CSF normalized masks, realigned&unwarp slice-time corrected normalized smoothed functional volumes, subject movement ''realignment'' and ''scrubbing'' 1st-level covariate'},{'INPUT: structural&functional&FieldMap volumes',  'OUTPUT (all in MNI-space): skull-stripped normalized structural volume, Grey/White/CSF normalized masks, realigned&unwarp slice-time corrected normalized smoothed functional volumes, subject movement ''realignment'' and ''scrubbing'' 1st-level covariate'},{'INPUT: structural&functional volumes','OUTPUT (all in subject-space): skull-stripped structural volume, Grey/White/CSF masks, realigned slice-time corrected coregistered functional volumes, subject movement ''realignment'' and ''scrubbing'' 1st-level covariate'},{'INPUT: structural&functional&VDM volumes','OUTPUT (all in subject-space): skull-stripped structural volume, Grey/White/CSF masks, realigned&unwarp slice-time corrected coregistered functional volumes, subject movement ''realignment'' and ''scrubbing'' 1st-level covariate'},{'INPUT: structural&functional volumes','OUTPUT (all in subject-space): skull-stripped structural volume, Grey/White/CSF masks, realigned slice-time corrected coregistered functional volumes, subject movement ''realignment'' and ''scrubbing'' 1st-level covariate'},...
    {'INPUT: functional data (volume files); structural volume; FreeSurfer-processed structural volume','OUTPUT: functional data (surface files)'}, ...
    {'INPUT: structural volume','OUTPUT: structural volume (same files re-oriented, not resliced)'}, {'INPUT: structural volume','OUTPUT: structural volume (same files translated, not resliced)'}, {'INPUT: structural volume','OUTPUT: skull-stripped structural volume, Grey/White/CSF masks (in same space as structural)'},...
    {'INPUT: structural volume','OUTPUT: skull-stripped normalized structural volume (in MNI space)'},{'INPUT: structural volume','OUTPUT: skull-stripped normalized structural volume, normalized Grey/White/CSF masks (all in MNI space)'},{'INPUT: structural volume; Grey/White/CSF masks (in same space as structural)','OUTPUT: skull-stripped normalized structural volume, normalized Grey/White/CSF masks (all in MNI space)'},...
    {'INPUT: structural volume; user-defined spatial deformation file (e.g. y_#.nii file)','OUTPUT: resampled structural volumes'}, ...
    {'INPUT: functional volumes','OUTPUT: temporal subset of functional volumes; temporal subset of first-level covariates (if already defined)'},{'INPUT: functional volumes','OUTPUT: functional volumes (same files re-oriented, not resliced)'},{'INPUT: functional volumes','OUTPUT: functional volumes (same files translated, not resliced)'},{'INPUT: structural and functional volumes','OUTPUT: functional volumes (same files translated, not resliced)'}, ...
    {'INPUT: functional volumes','OUTPUT: slice-timing corrected functional volumes'},{'INPUT: functional volumes','OUTPUT: band-pass filtered functional volumes'},{'INPUT: functional volumes; first-level covariates','OUTPUT: functional volumes with selected covariates regressed-out'},{'INPUT: functional volumes','OUTPUT: realigned functional volumes, mean functional image, subject movement ''realignment'' 1st-level covariate'},{'INPUT: functional volumes','OUTPUT: realigned&unwarp functional volumes, mean functional image, subject movement ''realignment'' 1st-level covariate'},...
    {'INPUT: functional volumes & VDM maps','OUTPUT: realigned&unwarp functional volumes, mean functional image, subject movement ''realignment'' 1st-level covariate'},{'INPUT: functional volumes, realignment parameters','OUTPUT: outlier scans 1st-level covariate, mean functional image, QA 2nd-level covariates'},{'INPUT: structural and mean functional volume (or first functional)','OUTPUT: coregistered functional volumes'},...
    {'INPUT: mean functional volume (or first functional)','OUTPUT: Grey/White/CSF masks (in same space as functional volume)'},...
    {'INPUT: functional volumes; user-defined spatial deformation file (e.g. y_#.nii file)','OUTPUT: resampled functional volumes'},...
    {'INPUT: functional volumes','OUTPUT: smoothed functional volumes'},{'INPUT: functional volumes','OUTPUT: motion masks'},...
    {'INPUT: structural volume; functional volumes','OUTPUT: skull-stripped normalized structural volume, normalized Grey/White/CSF masks; normalized functional volumes (all in MNI space)'},...
    {'INPUT: structural volume; functional volumes','OUTPUT: skull-stripped normalized structural volume, normalized functional volumes (all in MNI space)'},...
    {'INPUT: structural volume; functional volumes; Grey/White/CSF ROIs (in same space as structural volumes)','OUTPUT: skull-stripped normalized structural volume, normalized functional volumes, normalized Grey/White/CSF masks (all in MNI space)'},...
    {'INPUT: mean functional volume (or first functional)','OUTPUT: normalized functional volumes, normalized Grey/White/CSF masks'},...
    {'INPUT: mean functional volume (or first functional)','OUTPUT: normalized functional volumes'}, ...
    {'INPUT: functional volumes','OUTPUT: realigned functional volumes (same files re-oriented, not resliced), mean functional image, subject movement ''realignment'' 1st-level covariate'}, ...
    {'INPUT: structural and mean functional volume (or first functional)','OUTPUT: functional volumes coregistered to structural (direct normalization to MNI space + inverse deformation field transformation); Grey/White/CSF masks (in same space as functional volume)'}, ...
    {'INPUT: structural and mean functional volume (or first functional)','OUTPUT: functional volumes (all functional volumes are coregistered but not resliced)'}, ...
    {'INPUT: functional volumes','OUTPUT: none (one of the secondary datasets will point to current version of functional volumes)'}, ...
    {'INPUT: functional volumes','OUTPUT: none ("original data" secondary dataset will point to current version of functional volumes)'}, ...
    {'INPUT: functional volumes','OUTPUT: none ("subject-space data" secondary dataset will point to current version of functional volumes)'}, ...
    {'INPUT: functional volumes','OUTPUT: none ("mni-space data" secondary dataset will point to current version of functional volumes)'}, ...
    {'INPUT: functional volumes','OUTPUT: none ("surface-space data" secondary dataset will point to current version of functional volumes)'}, ...
    {'INPUT: functional volumes','OUTPUT: none ("smoothed data" datasets will point to current version of functional volumes)'}, ...
    {'INPUT: functional volumes (in secondary dataset)','OUTPUT: functional volumes (current functional volumes -typically the Primary dataset- will point to same files as one of the secondary datasets)'}, ...
    {'INPUT: functional volumes (in secondary dataset)','OUTPUT: functional volumes (current functional volumes -typically the Primary dataset- will point to same files as "original data" secondary dataset)'}, ...
    {'INPUT: functional volumes (in secondary dataset)','OUTPUT: functional volumes (current functional volumes -typically the Primary dataset- will point to same files as "subject-space data" secondary dataset)'}, ...
    {'INPUT: functional volumes (in secondary dataset)','OUTPUT: functional volumes (current functional volumes -typically the Primary dataset- will point to same files as "mni-space data" secondary dataset)'}, ...
    {'INPUT: functional volumes (in secondary dataset)','OUTPUT: functional volumes (current functional volumes -typically the Primary dataset- will point to same files as "surface-space data" secondary dataset)'}, ...
    {'INPUT: functional volumes (in secondary dataset)','OUTPUT: functional volumes (current functional volumes -typically the Primary dataset- will point to same files as "smoothed data" secondary dataset)'}, ...
    {'INPUT: functional volumes; Grey Matter ROI','OUTPUT: smoothed functional volumes'}, ...
    {'INPUT: functional data (volume files coregistered to structural); FreeSurfer-processed structural volume','OUTPUT: functional data (surface files)'}, ...
    {'INPUT: functional data (surface files)','OUTPUT: smoothed functional data (surface files)'}, ...
    {'INPUT: functional volumes & VDM maps','OUTPUT: Distortion corrected functional volumes'}, ...
    {'INPUT: functional volumes; double-echo FieldMap acquisition files (e.g. Magnitude+PhaseDiff volumes) in "fmap" dataset','OUTPUT: SPM VDM maps in "vdm" dataset'} ...
    };
%{'INPUT: mean functional volume (or first functional)','OUTPUT: normalized functional volumes'},{'INPUT: mean functional volume (or first functional)','OUTPUT: normalized functional volumes, normalized Grey/White/CSF masks '},...
steps_index=num2cell(1:numel(steps));
steps_combinedpipelines={...
    {'functional_label_as_original','functional_realign&unwarp','functional_center','functional_slicetime','functional_art','functional_segment&normalize_direct','functional_label_as_mnispace','structural_center','structural_segment&normalize','functional_smooth','functional_label_as_smoothed'},...
    {'functional_label_as_original','functional_vdm_create','functional_realign&unwarp&fieldmap','functional_center','functional_slicetime','functional_art','structural_center','functional_segment&normalize_indirect','functional_label_as_mnispace','functional_smooth','functional_label_as_smoothed'},...
    {'functional_label_as_original','functional_vdm_create','functional_realign&unwarp&fieldmap','functional_center','functional_slicetime','functional_art','functional_segment&normalize_direct','functional_label_as_mnispace','structural_center','structural_segment&normalize','functional_smooth','functional_label_as_smoothed'},...
    {'functional_label_as_original','functional_realign&unwarp','functional_slicetime','functional_art','functional_coregister_affine_noreslice','functional_label_as_subjectspace','functional_surface_resample','functional_label_as_surfacespace','functional_surface_smooth','functional_label_as_smoothed','structural_segment'},...
    {'functional_label_as_original','functional_vdm_create','functional_realign&unwarp&fieldmap','functional_slicetime','functional_art','functional_coregister_affine_noreslice','functional_label_as_subjectspace','functional_surface_resample','functional_label_as_surfacespace','functional_surface_smooth','functional_label_as_smoothed','structural_segment'},...
    {'functional_label_as_original','functional_realign&unwarp','functional_slicetime','functional_art','functional_coregister_nonlinear','functional_label_as_subjectspace','functional_surface_resample','functional_label_as_surfacespace','functional_surface_smooth','functional_label_as_smoothed'},...
    {'functional_coregister_affine_noreslice','functional_surface_resample'} ...
    };
for n=1:numel(steps_combinedpipelines),
    [ok,idx]=ismember(steps_combinedpipelines{n},steps);
    if ~all(ok), error('preprocessing step names have changed'); end
    steps_index{n}=idx(:)';
end
if nargin>0&&ischar(STEPS)&&strcmp(STEPS,'steps'), ok=sort(steps); return; end
steps_pipelines=cellfun('length',steps_index)>1;
steps_default=cellfun('length',regexp(steps,'^default_'))>0;
dogui=false;
sets=0;
subjects=1:CONN_x.Setup.nsubjects;
sessions=1:max(CONN_x.Setup.nsessions);
doimport=true;
typeselect='';
multiplesteps=iscell(STEPS)&numel(STEPS)>1;
showallsteps=false;
voxelsize_anat=1;
voxelsize_func=2;
boundingbox=[-90,-126,-72;90,90,108]; % default bounding-box
interp=[];
fwhm=[];
diffusionsteps=[];
sliceorder=[];
sliceorder_select=[];
label={};
load_label={};
ta=[];
unwarp=[];
removescans=[];
bp_filter=[];
bp_keep0=1;
reg_names={};
reg_dimensions=[];
reg_deriv=[];
reg_filter=[];
reg_detrend=1;
reg_lag=[];
reg_lagmax=8;
reg_skip=0;
reorient=[];
respatialdef=[];
coregtomean=1;
rtm=0;
coregsource={};
applytofunctional=false;
tpm_template=[];
affreg=[];
tpm_ngaus=[];
vdm_et1=[]; % eg. 2.84, 4.37;
vdm_et2=[]; % eg. 5.30, 6.83
vdm_ert=[]; % eg. 37.6
vdm_blip=[];% eg. -1
vdm_type=[];
vdm_fmap=[];
art_thresholds=[];
art_useconservative=1;
art_global_thresholds=[9 5 3];
art_motion_thresholds=[2 .9 .5];
art_global_threshold=art_global_thresholds(1+art_useconservative); % default art scan-to-scan global signal z-value thresholds
art_motion_threshold=art_motion_thresholds(1+art_useconservative); % default art scan-to-scan composite motion mm thresholds
art_use_diff_motion=1;
art_use_diff_global=1;
art_use_norms=1;
art_force_interactive=0;
art_drop_flag=0;
art_gui_display=true;
parallel_profile=[];
parallel_N=0;
functional_template=fullfile(fileparts(which('spm')),'templates','EPI.nii');
if ~conn_existfile(functional_template), functional_template=fullfile(fileparts(which('spm')),'toolbox','OldNorm','EPI.nii'); end
structural_template=fullfile(fileparts(which('spm')),'templates','T1.nii');
if ~conn_existfile(structural_template), structural_template=fullfile(fileparts(which('spm')),'toolbox','OldNorm','T1.nii'); end
selectedstep=1;
if ~isempty(STEPS)&&(ischar(STEPS)||(iscell(STEPS)&&numel(STEPS)==1))
    STEPS=regexprep(char(STEPS),{'^default_mniphase$','^default_mnidirectphase$','^default_ssphase$'},{'default_mnifield','default_mnidirectfield','default_ssfield'});
    [ok,idx]=ismember(lower(STEPS),steps(steps_pipelines));
    if ok, STEPS=steps(steps_index{idx}); selectedstep=idx;
    else
        lSTEPS=regexprep(lower(STEPS),'^run_|^update_|^interactive_','');
        if ~isempty(regexp(char(lSTEPS),'^functional_label_as_'))||~isempty(regexp(char(lSTEPS),'^functional_load_from_'))||ismember(lSTEPS,steps), STEPS=cellstr(STEPS);
        elseif conn_existfile(STEPS), load(STEPS,'STEPS','coregtomean'); if isempty(coregtomean), coregtomean=1; end
        else error('STEP name %s is not a valid preprocessing step or an existing preprocessing-pipeline file',STEPS);
        end
    end
elseif ~isempty(STEPS)
    if ischar(STEPS), STEPS=cellstr(STEPS); end
    STEPS=regexprep(STEPS,{'^default_mniphase$','^default_mnidirectphase$','^default_ssphase$'},{'default_mnifield','default_mnidirectfield','default_ssfield'});
    [ok,idx]=ismember(lower(STEPS),steps(steps_pipelines));
    nSTEPS={};
    for n1=1:numel(STEPS)
        if ok(n1), nSTEPS=[nSTEPS steps(steps_index{idx(n1)})]; selectedstep=idx(n1);
        else nSTEPS=[nSTEPS STEPS(n1)];
        end
    end
    STEPS=nSTEPS;
end
ok=0;
for n1=1:2:numel(options)-1,
    switch(lower(options{n1}))
        case 'select',
            typeselect=lower(char(options{n1+1}));
        case 'showallsteps',
            showallsteps=options{n1+1};
        case 'multiplesteps',
            multiplesteps=options{n1+1};
        case 'fwhm',
            fwhm=options{n1+1};
        case 'diffusionsteps',
            diffusionsteps=options{n1+1};
        case 'sliceorder',
            sliceorder=options{n1+1};
            if iscell(sliceorder), sliceorder=[sliceorder{:}]; end
        case 'label',
            label=options{n1+1};
            if ~iscell(label), label={label}; end
        case 'load_label',
            load_label=options{n1+1};
            if ~iscell(load_label), load_label={load_label}; end
        case 'ta',
            ta=options{n1+1};
        case 'unwarp',  % note: deprecated over CONN_x.Setup.unwarp_functional field
            unwarp=options{n1+1};
        case 'bp_filter',
            bp_filter=options{n1+1};
        case 'bp_keep0',
            bp_keep0=options{n1+1};
        case 'reg_names',
            reg_names=options{n1+1};
        case 'reg_dimensions',
            reg_dimensions=options{n1+1};
        case 'reg_deriv',
            reg_deriv=options{n1+1};
        case 'reg_filter',
            reg_filter=options{n1+1};
        case 'reg_detrend',
            reg_detrend=options{n1+1};
        case 'reg_lag',
            reg_lag=options{n1+1};
        case 'reg_lagmax',
            reg_lagmax=options{n1+1};
        case 'reg_skip',
            reg_skip=options{n1+1};
        case 'removescans',
            removescans=options{n1+1};
        case 'applytofunctional',
            applytofunctional=options{n1+1};
        case 'coregtomean',
            coregtomean=options{n1+1};
        case 'rtm',
            rtm=options{n1+1};
        case 'coregsource', % note: deprecated over CONN_x.Setup.coregsource_functional field
            coregsource=options{n1+1};
        case 'reorient',
            reorient=options{n1+1};
        case 'respatialdef',
            respatialdef=options{n1+1};
        case 'art_thresholds',
            art_thresholds=options{n1+1};
        case 'sets',
            sets=options{n1+1};
        case 'subjects',
            subjects=options{n1+1};
        case 'sessions',
            sessions=options{n1+1};
        case 'voxelsize',
            voxelsize_anat=options{n1+1};
            voxelsize_func=options{n1+1};
        case 'voxelsize_anat',
            voxelsize_anat=options{n1+1};
        case 'voxelsize_func',
            voxelsize_func=options{n1+1};
        case 'boundingbox',
            boundingbox=options{n1+1};
        case 'interp'
            interp=options{n1+1};
        case 'doimport',
            doimport=options{n1+1};
        case 'dogui',
            dogui=options{n1+1};
        case {'functional_template','template_functional','template_func'}
            functional_template=char(options{n1+1});
        case {'structural_template','template_structural','template_anat'}
            structural_template=char(options{n1+1});
        case 'usespm8methods',
            PREFERSPM8OVERSPM12=options{n1+1};
        case 'affreg'
            affreg=char(options{n1+1});
        case 'tpm_template',
            tpm_template=options{n1+1};
        case 'tpm_ngaus',
            tpm_ngaus=options{n1+1};
        case 'vdm_et1',
            vdm_et1=options{n1+1};
        case 'vdm_et2',
            vdm_et2=options{n1+1};
        case 'vdm_ert',
            vdm_ert=options{n1+1};
        case 'vdm_blip',
            vdm_blip=options{n1+1};
        case 'vdm_type',
            vdm_type=options{n1+1};
        case 'vdm_fmap',
            vdm_fmap=options{n1+1};
        case 'parallel_profile'
            parallel_profile=char(options{n1+1});
        case 'parallel_N'
            parallel_N=options{n1+1};
        otherwise
            error(['unrecognized option ',options{n1}]);
    end
end
if isfield(CONN_x,'pobj')&&isstruct(CONN_x.pobj)&&isfield(CONN_x.pobj,'subjects'), subjects=CONN_x.pobj.subjects; end % this field overwrites user-defined options
if ~isempty(sets)&&ischar(sets), sets=conn_datasetlabel(sets,'error'); end

if ~nargin||isempty(STEPS)||dogui,
    dogui=true;
    if showallsteps, idx=1:numel(steps);
    else idx=find(cellfun('length',regexp(steps,'^structural|^functional|^default')));
    end
    if ~isempty(typeselect)
        switch(typeselect)
            case 'structural', idx=find(cellfun('length',regexp(steps,'^structural'))&~steps_pipelines);
            case 'functional', idx=find(cellfun('length',regexp(steps,'^functional'))&~steps_pipelines);
        end
    end
    steps0=steps;
    steps=steps(idx);
    steps_names=steps_names(idx);
    steps_descr=steps_descr(idx);
    steps_pipelines=steps_pipelines(idx);
    steps_default=steps_default(idx);
    steps_index=steps_index(idx);
    [nill,steps_order]=sort(steps_names);
    steps_order=[sort(steps_order(steps_default(steps_order))) steps_order(~steps_default(steps_order))];
    scalefig=1+multiplesteps;
    dlg.steps=steps;
    dlg.steps_names=steps_names;
    dlg.steps_descr=steps_descr;
    dlg.steps_index=steps_index;
    dlg.steps_order=steps_order;
    dlg.fig=figure('units','norm','position',[.2,.3,.5+.2*(1|multiplesteps),.6],'menubar','none','numbertitle','off','name','SPM data preprocessing step','color',1*[1 1 1]);
    if multiplesteps,
        uicontrol('style','frame','units','norm','position',[0,.57,1,.43],'backgroundcolor',.9*[1 1 1],'foregroundcolor',.9*[1 1 1]);
        uicontrol('style','frame','units','norm','position',[.05,.2,.9,.33],'backgroundcolor',1*[1 1 1],'foregroundcolor',.8*[1 1 1]);
        %uicontrol('style','frame','units','norm','position',[.025,.025,.95,.55],'backgroundcolor',1*[1 1 1],'foregroundcolor',.75*[1 1 1],'fontsize',9+font_offset);
    end
    htm0=uicontrol('style','text','units','norm','position',[.05,.9,.85,.05],'backgroundcolor',1*[1 1 1],'foregroundcolor','k','horizontalalignment','left','string','Select individual preprocessing step:','fontweight','bold','fontsize',9+font_offset);
    dlg.m0=uicontrol('style','popupmenu','units','norm','position',[.05,.85,.85,.05],'string',steps_names(steps_order),'value',find(ismember(steps_order,selectedstep)),'backgroundcolor',1*[1 1 1],'foregroundcolor','k','tooltipstring','Select a data preprocessing step','callback',@(varargin)conn_setup_preproc_update,'fontsize',9+font_offset);
    if multiplesteps,
        set(htm0,'string','List of all available preprocessing steps:');
        set(dlg.m0,'tooltipstring','Select a data preprocessing step or pipeline and click ''Add'' to add it to your data preprocessing pipeline');
    end
    dlg.m6=uicontrol('style','text','units','norm','position',[.05,.725,.85,.1],'max',2,'string','','backgroundcolor',1*[1 1 1],'enable','inactive','horizontalalignment','left','fontsize',9+font_offset);
    dlg.m4=uicontrol('style','checkbox','units','norm','position',[.05,.65,.85,.05],'value',~coregtomean,'string','First functional volume as reference','backgroundcolor',1*[1 1 1],'tooltipstring','<HTML>Uses firts functional volume as reference in coregistration/normalization step <br/> - if unchecked coregistration/normalization uses mean-volume as reference instead<br/> - note: mean volume is created during realignment</HTML>','visible','off','fontsize',9+font_offset);
    %dlg.m3=uicontrol('style','checkbox','units','norm','position',[.1,.5,.8/scalefig,.05],'value',applytofunctional,'string','Apply structural deformation field to functional data as well','backgroundcolor',1*[1 1 1],'tooltipstring','Apply structural deformation field computed during structural normalization/segmentation step to coregistered functional data as well','visible','off','fontsize',9+font_offset);
    dlg.m2=uicontrol('style','popupmenu','units','norm','position',[.05,.55,.85,.05],'value',1,'string',{'Run process and import results to CONN project','Run process only (do not import results)','Interactive SPM batch editor only (do not run process)'}','backgroundcolor',1*[1 1 1],'fontsize',9+font_offset);
    dlg.m1=uicontrol('style','checkbox','units','norm','position',[.05,.08,.3,.05],'value',1,'string','Process all subjects','backgroundcolor',1*[1 1 1],'tooltipstring','Apply this preprocessing to all subjects in your curent CONN project','callback',@(varargin)conn_setup_preproc_update,'fontsize',9+font_offset);
    dlg.m5=uicontrol('style','listbox','units','norm','position',[.35,.11,.15,.08],'max',2,'string',arrayfun(@(n)sprintf('Subject%d',n),1:CONN_x.Setup.nsubjects,'uni',0),'backgroundcolor',1*[1 1 1],'tooltipstring','Select subjects','visible','off','fontsize',9+font_offset);
    dlg.m1b=uicontrol('style','checkbox','units','norm','position',[.05,.03,.3,.05],'value',1,'string','Process all sessions','backgroundcolor',1*[1 1 1],'tooltipstring','Apply this preprocessing to all sessions in your curent CONN project','callback',@(varargin)conn_setup_preproc_update,'fontsize',9+font_offset);
    dlg.m5b=uicontrol('style','listbox','units','norm','position',[.35,.03,.15,.08],'max',2,'string',arrayfun(@(n)sprintf('Session%d',n),1:max(CONN_x.Setup.nsessions),'uni',0),'backgroundcolor',1*[1 1 1],'tooltipstring','Select sessions','visible','off','fontsize',9+font_offset);
    if any(cellfun('length',regexp(dlg.steps_names,'^functional'))), dlg.m10=uicontrol('style','popupmenu','units','norm','position',[.05,.13,.3,.05],'value',1+sets,'string',[{'Process primary dataset'},arrayfun(@(n)sprintf('Process secondary dataset #%d %s',n,regexprep(CONN_x.Setup.secondarydataset(n).label,'(.+)','($1)')),1:numel(CONN_x.Setup.secondarydataset),'uni',0)],'tooltipstring','<HTML>Apply this preprocessing to selected functional dataset in your curent CONN project (the same dataset will hold the OUTPUT files of each preprocessing step)<br/> - note: only when preprocessing the primary dataset (dataset-0) non-essential OUTPUTS of each preprocessing step (ROIs, first- and <br/>second- level covariates) will be automatically imported into your CONN project</HTML>','visible','on','fontsize',9+font_offset);
    else dlg.m10=[];
    end
    [tstr,tidx]=conn_jobmanager('profiles');
    tnull=find(strcmp('Null profile',conn_jobmanager('profiles')));
    tlocal=find(strcmp('Background process (Unix,Mac)',tstr),1);
    tvalid=setdiff(1:numel(tstr),tnull);
    tstr=cellfun(@(x)sprintf('distributed processing (run on %s)',x),tstr,'uni',0);
    if 1, tvalid=tidx; if isunix&&~isempty(tlocal)&&~ismember(tlocal,tvalid), tvalid=[tvalid(:)' tlocal]; end
    elseif 1, tvalid=tidx; % show only default scheduler
    else tstr{tidx}=sprintf('<HTML><b>%s</b></HTML>',tstr{tidx});
    end
    toptions=[{'local processing (run on this computer)' 'queue/script it (save as scripts to be run later)'} tstr(tvalid)];
    if CONN_gui.isremote
        info=conn_remotely('info');
        if isfield(info,'host')&&~isempty(info.host), tnameserver=info.host;
        elseif isfield(info,'remote_ip')&&~isempty(info.remote_ip), tnameserver=info.remote_ip;
        else tnameserver='CONN server';
        end
        toptions=regexprep(toptions,'\<run on (this computer)?',['run on ',tnameserver,' ']);
    end    
    dlg.m9=uicontrol('style','popupmenu','units','norm','position',[.55,.12,.40,.05],'string',toptions,'value',1,'backgroundcolor',1*[1 1 1],'fontsize',8+CONN_gui.font_offset);
    if multiplesteps, dlg.m11=uicontrol('style','pushbutton','units','norm','position',[.55,.04,.2,.07],'string','Start','tooltipstring','Accept changes and run data preprocessing pipeline','callback','set(gcbf,''userdata'',0); uiresume(gcbf)','fontsize',9+font_offset,'fontweight','bold');
    else              dlg.m11=uicontrol('style','pushbutton','units','norm','position',[.55,.04,.2,.07],'string','Start','tooltipstring','Accept changes and run data preprocessing step','callback','set(gcbf,''userdata'',0); uiresume(gcbf)','fontsize',9+font_offset);
    end
    dlg.m12=uicontrol('style','pushbutton','units','norm','position',[.75,.04,.2,.07],'string','Cancel','callback','delete(gcbf)','fontsize',9+font_offset);
    if multiplesteps
        set([htm0 dlg.m0],'visible','off');
        dlg.m0b=uicontrol('style','text','units','norm','position',[.07,.45,.86,.06],'backgroundcolor',1*[1 1 1],'foregroundcolor',0*[1 1 1],'horizontalalignment','left','string','','fontweight','bold','fontsize',9+font_offset);
        set(dlg.m6,'position',[.07,.26,.86,.18],'backgroundcolor',1*[1 1 1]);
        set(dlg.m4,'position',[.07,.215,.86,.04],'backgroundcolor',1*[1 1 1]);
        set(dlg.m2,'visible','off');%'string',{'Run process and import results to CONN project'});
        %set(dlg.m3,'position',get(dlg.m3,'position')-[0 .075 0 0]);
        %set(dlg.m4,'position',get(dlg.m4,'position')-[0 .075 0 0]);
        set(dlg.fig,'name','SPM data preprocessing pipeline');
        uicontrol('style','text','units','norm','position',[.05,.915,.85,.05],'backgroundcolor',.9*[1 1 1],'foregroundcolor','k','horizontalalignment','left','string','Data preprocessing pipeline:','fontweight','bold','fontsize',11+font_offset);
        dlg.m7=uicontrol('style','listbox','units','norm','position',[.05,.59,.78,.325],'max',2,'string',{},'backgroundcolor',.9*[1 1 1],'tooltipstring','Define sequence of preprocessing steps','fontsize',9+font_offset,'callback','dlg=get(gcbo,''userdata''); str=get(gcbo,''string''); val=get(gcbo,''value''); if numel(val)==1, idx=find(strcmp(dlg.steps_names(dlg.steps_order),str{val})); if numel(idx)==1, set(dlg.m0,''value'',idx); feval(get(dlg.m0,''callback'')); end; end');
        dlg.m8a=uicontrol('style','pushbutton','units','norm','position',[.84,.87,.11,.045],'string','Add','fontweight','bold','tooltipstring','Adds new data preprocessing step to this list','callback',@conn_setup_preproc_update_add,'fontsize',9+font_offset);
        dlg.m8b=uicontrol('style','pushbutton','units','norm','position',[.84,.825,.11,.045],'string','Remove','tooltipstring','Removes selected preprocessing step from this list','callback','dlg=get(gcbo,''userdata''); str=get(dlg.m7,''string''); str=str(setdiff(1:numel(str),get(dlg.m7,''value''))); set(dlg.m7,''string'',str,''value'',[]); feval(get(dlg.m0,''callback'')); ','fontsize',9+font_offset);
        dlg.m8g=uicontrol('style','pushbutton','units','norm','position',[.84,.78,.11,.045],'string','Clear','tooltipstring','Removes all preprocessing steps from this list','callback','dlg=get(gcbo,''userdata''); set(dlg.m7,''string'',{},''value'',[]); feval(get(dlg.m0,''callback'')); ','fontsize',9+font_offset);
        dlg.m8c=uicontrol('style','pushbutton','units','norm','position',[.84,.735,.11,.045],'string','Move up','tooltipstring','Moves selected preprocessing step up in this list','callback','dlg=get(gcbo,''userdata''); str=get(dlg.m7,''string''); val=get(dlg.m7,''value''); idx=1:numel(str); idx(val)=min(idx(val))-1.5; [nill,idx]=sort(idx); str=str(idx); set(dlg.m7,''string'',str,''value'',find(rem(nill,1)~=0));','fontsize',9+font_offset);
        dlg.m8d=uicontrol('style','pushbutton','units','norm','position',[.84,.69,.11,.045],'string','Move down','tooltipstring','Moves selected preprocessing step down this list','callback','dlg=get(gcbo,''userdata''); str=get(dlg.m7,''string''); val=get(dlg.m7,''value''); idx=1:numel(str); idx(val)=max(idx(val))+1.5; [nill,idx]=sort(idx); str=str(idx); set(dlg.m7,''string'',str,''value'',find(rem(nill,1)~=0));','fontsize',9+font_offset);
        dlg.m8e=uicontrol('style','pushbutton','units','norm','position',[.84,.635,.11,.045],'string','Save','tooltipstring','Saves this data preprocessing pipeline list for future use','callback',@conn_setup_preproc_save,'fontsize',9+font_offset);
        dlg.m8f=uicontrol('style','pushbutton','units','norm','position',[.84,.59,.11,.045],'string','Load','tooltipstring','Loads data preprocessing pipeline list from file','callback',@conn_setup_preproc_load,'fontsize',9+font_offset);
        set([dlg.m7 dlg.m8a dlg.m8b dlg.m8c dlg.m8d dlg.m8e dlg.m8f dlg.m8g],'userdata',dlg);
    else dlg.m7=[]; dlg.m0b=[];
    end
    set([dlg.m0 dlg.m1 dlg.m1b],'userdata',dlg);
    %if isempty(STEPS)&&multiplesteps, STEPS=steps(steps_index{1}); end
    if ~isempty(STEPS)
        [tok,idx]=ismember(STEPS,steps);
        if multiplesteps, set(dlg.m7,'string',steps_names(idx(tok>0))');
        else set(dlg.m0, 'value',find(ismember(steps_order,idx(tok>0)),1));
        end
    end
    conn_setup_preproc_update(dlg.m0);
    if multiplesteps,
        conn_setup_preproc_load(dlg.m8f);
        if isempty(get(dlg.m7,'string')), conn_setup_preproc_update_add(dlg.m8a); end
    end
    uiwait(dlg.fig);
    
    if ~ishandle(dlg.fig), return; end
    pressedok=get(dlg.fig,'userdata');
    if isempty(pressedok), return; end
    if multiplesteps
        STEPS=get(dlg.m7,'string');
        [tok,idx]=ismember(STEPS,steps_names);
        STEPS=steps(idx(tok>0));
    else
        STEPS=steps(dlg.steps_order(get(dlg.m0,'value')));
        %idx0=find(steps_pipelines);
        %[ok,idx]=ismember(lower(STEPS),steps(idx0));
        %if ok, STEPS=steps0(steps_index{idx0(idx)}); end
    end
    %STEP_name=steps_names{get(dlg.m0,'value')};
    %if any(ismember(STEPS,{'structural_segment&normalize','structural_normalize'})), applytofunctional=get(dlg.m3,'value'); end
    if any(cellfun('length',regexp(STEPS,'^functional_coregister|^functional_normalize|functional_segment|functional_center'))), coregtomean=~get(dlg.m4,'value'); end
    if ~get(dlg.m1,'value'), subjects=get(dlg.m5,'value'); end
    if ~get(dlg.m1b,'value'), sessions=get(dlg.m5b,'value'); end
    if ~isempty(dlg.m10), sets=get(dlg.m10,'value')-1; end
    dorun=get(dlg.m2,'value');
    doparallel=get(dlg.m9,'value');
    if multiplesteps, conn_setup_preproc_save(dlg.m8f); end
    delete(dlg.fig);
    switch(dorun)
        case 1, STEPS=cellfun(@(x)['run_',x],STEPS,'uni',0); doimport=true;
        case 2, STEPS=cellfun(@(x)['run_',x],STEPS,'uni',0); doimport=false;
        case 3, STEPS=cellfun(@(x)['interactive_',x],STEPS,'uni',0); doimport=false;
        case 4, STEPS=cellfun(@(x)['update_',x],STEPS,'uni',0); doimport=true;
    end
    
    if doparallel>1
        if doparallel==2, parallel_profile=find(strcmp('Null profile',conn_jobmanager('profiles')));
        else parallel_profile=tvalid(doparallel-2);
            if conn_jobmanager('ispending')
                answ=conn_questdlg({'There are previous pending jobs associated with this project','This job cannot be submitted until all pending jobs finish',' ','Would you like to queue this job for later?','(pending jobs can be seen at Tools.Cluster/HPC.View pending jobs'},'Warning','Queue','Cancel','Queue');
                if isempty(answ)||strcmp(answ,'Cancel'), ok=false; end
                parallel_profile=find(strcmp('Null profile',conn_jobmanager('profiles')));
            end
        end
        if numel(subjects)>1,
            answer=inputdlg(sprintf('Number of parallel jobs? (1-%d)',numel(subjects)),'',1,{num2str(1)});
            if isempty(answer)||isempty(str2num(answer{1})), return; end
            parallel_N=str2num(answer{1});
        else parallel_N=1;
        end
    end
end

lSTEPS=regexprep(lower(STEPS),'^run_|^update_|^interactive_','');
sliceorder_select_options={'ascending','descending','interleaved (middle-top)','interleaved (bottom-up)','interleaved (top-down)','interleaved (Siemens)','interleaved (Philips)','BIDS'};
sliceorder_select_options_extended={'ascending (e.g. 1,2,3...9,10)','descending (e.g. 10,9,8...2,1)','interleaved middle-top (e.g. 10,5,9,4...6 1)','interleaved bottom-up (e.g. 1,3,5...8,10)','interleaved top-down (e.g. 10,8,6...3,1)','interleaved (Siemens) (e.g. 2,4,6...7,9)','interleaved (Philips) (e.g. 1,4,7...6,9)','BIDS (from functional .json metadata)'};
if any(ismember('functional_slicetime',lSTEPS))
    if ischar(sliceorder),
        [slok,sliceorder_select]=ismember(sliceorder,sliceorder_select_options);
        if ~slok, conn_disp(sprintf('Warning: incorrect sliceorder name %s',sliceorder)); sliceorder_select=[]; end
        sliceorder=[];
    end
    if isempty(sliceorder)&&(isempty(sliceorder_select)||dogui)
        [sliceorder_select,tok] = listdlg('PromptString','Select slice order:','ListSize',[400 200],'SelectionMode','single','InitialValue',sliceorder_select,'ListString',[sliceorder_select_options_extended,{'manually define','do not know (skip slice timing correction)'}]);
        if isempty(sliceorder_select), return; end
        if sliceorder_select==numel(sliceorder_select_options)+1
            sliceorder=inputdlg(['Slice order? (enter slice indexes from z=1 -first slice in image- to z=? -last slice- in the order they were acquired). Alternatively enter acquisition time of each slice in milliseconds (e.g. for multiband sequences). Press Cancel to enter this information at a later point (e.g. separately for each subject)'],'conn_setup_preproc',1,{' '});
            if ~isempty(sliceorder), sliceorder=str2num(sliceorder{1});
            else sliceorder=[];
            end
        elseif sliceorder_select==numel(sliceorder_select_options)+2
            sliceorder_select=[]; STEPS=STEPS(~ismember(lSTEPS,'functional_slicetime'));
        end
    end
end

if any(ismember('functional_removescans',lSTEPS))
    if isempty(removescans)||dogui
        if isempty(removescans), removescans=0; end
        removescans=inputdlg('Enter number of initial scans to remove','conn_setup_preproc',1,{num2str(removescans)});
        if isempty(removescans), return; end
        removescans=str2num(removescans{1});
    end
end

if any(ismember('functional_bandpass',lSTEPS))
    if isempty(bp_filter)||dogui
        if isempty(bp_filter), bp_filter=[0.01 0.10]; end
        bp_filter=inputdlg('Enter band-pass filter thresholds in Hz','conn_setup_preproc',1,{num2str(bp_filter)});
        if isempty(bp_filter), return; end
        bp_filter=str2num(bp_filter{1});
    end
end

if any(ismember('functional_regression',lSTEPS))
    if isempty(reg_names)||dogui
        temp_reg_names=CONN_x.Setup.l1covariates.names(1:end-1);
        if isempty(reg_names), answ=~cellfun('length',regexp(temp_reg_names,'^QC_'));
        else answ=reshape(ismember(temp_reg_names,reg_names),1,[]);
        end
        if isempty(reg_deriv), tansw=~cellfun('length',regexpi(temp_reg_names,'^effect of|realign|movement|motion'));
        else tansw=reshape(reg_deriv==0,1,[]);
        end
        answ=reshape([answ&tansw;answ&~tansw],1,[]);
        temp_reg_names=reshape([temp_reg_names;cellfun(@(x)sprintf('%s + 1st order temporal derivative',x),temp_reg_names,'uni',0)],1,[]);
        temp_reg_names=[{'time (detrending)'},temp_reg_names];
        if isempty(reg_detrend), answ=[true, answ];
        else answ=[reg_detrend, answ];
        end
        answ=listdlg('Promptstring','Select model regressors','selectionmode','multiple','liststring',temp_reg_names,'initialvalue',find(answ),'ListSize',[320 300]);
        if numel(answ)>=1,
            reg_detrend=any(answ==1);
            answ=answ(answ>1);
            uansw=unique(2*floor(answ/2));
            reg_names=temp_reg_names(uansw);
            reg_deriv=double(ismember(uansw+1,answ));
        else return;
        end
    end
end

if dogui&&any(ismember(lSTEPS,{'functional_vdm_create'}))
    thfig=figure('units','norm','position',[.4,.4,.35,.3],'color',1*[1 1 1],'name','VDM create settings','numbertitle','off','menubar','none');
    ht1=uicontrol('style','popupmenu','units','norm','position',[.1,.85,.8,.1],'string',arrayfun(@(n)sprintf('Fieldmap location: secondary dataset #%d %s',n,regexprep(CONN_x.Setup.secondarydataset(n).label,'(.+)','($1)')),1:numel(CONN_x.Setup.secondarydataset),'uni',0),'value',1,'backgroundcolor',1*[1 1 1],'tooltipstring','defines location of available fieldmap-sequence files');
    ht2=uicontrol('style','popupmenu','units','norm','position',[.1,.75,.8,.1],'string',{'Fieldmap type: automatically determine','Fieldmap type: Magnitude,Phasediff files','Fieldmap type: Real1,Imag1,Real2,Imag2 files','Fieldmap type: Pre-computed fieldmap file (Hz)'},'value',1,'backgroundcolor',1*[1 1 1],'tooltipstring','defines type of available fieldmap-sequence files');
    ht3=uicontrol('style','checkbox','units','norm','position',[.1,.64,.8,.1],'string','Read double-echo timing from BIDS / .json files','value',1,'backgroundcolor',1*[1 1 1],'tooltipstring','use information in .json sidecar files to estimate EchoTime and EPI Total Readout Time values');
    ht4=[];ht5=[];ht6=[];
    ht4a=uicontrol('style','text','units','norm','position',[.1,.5,.6,.1],'string','Fieldmap''s Short Echo Time (in ms)','horizontalalignment','left','backgroundcolor',1*[1 1 1],'enable','off');
    ht4=uicontrol('style','edit','units','norm','position',[.7,.5,.2,.1],'string',num2str(vdm_et1),'tooltipstring','defines Echo Time (in ms units) of first dual-echo acquisition (leave empty to import from .json / BIDS file)','enable','off');
    ht5a=uicontrol('style','text','units','norm','position',[.1,.4,.6,.1],'string','Fieldmap''s Long Echo Time (in ms)','horizontalalignment','left','backgroundcolor',1*[1 1 1],'enable','off');
    ht5=uicontrol('style','edit','units','norm','position',[.7,.4,.2,.1],'string',num2str(vdm_et2),'tooltipstring','defines Echo Time (in ms units) of second dual-echo acquisition (leave empty to import from .json / BIDS file)','enable','off');
    ht6a=uicontrol('style','text','units','norm','position',[.1,.3,.6,.1],'string','Functional''s Total Readout Time (in ms)','horizontalalignment','left','backgroundcolor',1*[1 1 1],'enable','off');
    ht6=uicontrol('style','edit','units','norm','position',[.7,.3,.2,.1],'string',num2str(vdm_ert),'tooltipstring','defines EPI Total Readout Time (in ms units) of functional data (note: equal to 1000/BandwidthPerPixelPhaseEncode;  leave empty to import from .json / BIDS file)','enable','off');
    ht7a=uicontrol('style','text','units','norm','position',[.1,.2,.6,.1],'string','Functional''s Blip direction (+1,-1,S,R)','horizontalalignment','left','backgroundcolor',1*[1 1 1],'enable','off');
    if isempty(vdm_blip), tvdm_blip=-1; else tvdm_blip=vdm_blip; end
    if ~ischar(tvdm_blip), tvdm_blip=num2str(tvdm_blip); end
    ht7=uicontrol('style','edit','units','norm','position',[.7,.2,.2,.1],'string',tvdm_blip,'tooltipstring','defines k-space traversal blip direction: +1 for positive direction, -1 for negative direction, leave empty or set to ''S'' to derive this information from the PhaseEncodingDirection field in .json/BIDS file, set to ''R'' to try the reverse direction of ''S''','enable','off');
    uicontrol('style','pushbutton','string','OK','units','norm','position',[.1,.01,.38,.15],'callback','uiresume');
    uicontrol('style','pushbutton','string','Cancel','units','norm','position',[.51,.01,.38,.15],'callback','delete(gcbf)');
    onoff={'on','off'};
    if isempty(vdm_fmap), vdm_fmap=conn_datasetlabel('fmap'); end
    if ischar(vdm_fmap), vdm_fmap=conn_datasetlabel(vdm_fmap); end
    if isempty(vdm_fmap), vdm_fmap=1; end
    set(ht1,'value',vdm_fmap);
    if isempty(vdm_type), set(ht2,'value',1);
    else set(ht2,'value',1+vdm_type); set([ht4 ht4a ht5 ht5a],'visible',onoff{1+(get(ht2,'value')==4)});
    end
    if isempty(vdm_et1)&&isempty(vdm_et2)&&isempty(vdm_ert)&&isempty(vdm_blip), set(ht3,'value',1);
    else set(ht3,'value',0); set([ht4 ht4a ht5 ht5a ht6 ht6a ht6 ht6a ht7 ht7a],'enable','on');
    end
    set(ht2,'userdata',[],'callback',@(varargin)set([ht4 ht4a ht5 ht5a],'visible',onoff{1+(get(ht2,'value')==4)}));
    set(ht3,'userdata',[],'callback',@(varargin)set([ht4 ht4a ht5 ht5a ht6 ht6a ht7 ht7a],'enable',onoff{1+get(ht3,'value')}));
    uiwait(thfig);
    if ~ishandle(thfig), return; end
    vdm_fmap=get(ht1,'value');
    vdm_type=get(ht2,'value')-1; if ~vdm_type, vdm_type=[]; end
    if get(ht3,'value'), vdm_et1=[]; vdm_et2=[]; vdm_ert=[]; vdm_blip=[];
    else
        temp=get(ht4,'string'); if isempty(temp), vdm_et1=temp; elseif ~isempty(str2num(temp)), vdm_et1=str2num(temp); else error('unable to interpret vdm_et1 input %s',temp); end
        temp=get(ht5,'string'); if isempty(temp), vdm_et2=temp; elseif ~isempty(str2num(temp)), vdm_et2=str2num(temp); else error('unable to interpret vdm_et2 input %s',temp); end
        temp=get(ht6,'string'); if isempty(temp), vdm_ert=temp; elseif ~isempty(str2num(temp)), vdm_ert=str2num(temp); else error('unable to interpret vdm_ert input %s',temp); end
        temp=get(ht7,'string'); if isempty(temp)||isequal(temp,'s')||isequal(temp,'S'), vdm_blip=[]; elseif ~isempty(str2num(temp)), vdm_blip=str2num(temp); elseif isequal(temp,'r')||isequal(temp,'R'), vdm_blip=0; else error('unable to interpret vdm_blip input %s',temp); end
    end
    delete(thfig);
    drawnow;
end

if any(ismember({'structural_manualorient','functional_manualorient'},lSTEPS))
    if isempty(reorient)||dogui
        ntimes=sum(ismember(lSTEPS,{'structural_manualorient','functional_manualorient'}));
        reorient={};
        opts={'translation to 0/0/0 coordinates',nan;
            '90-degree rotation around x-axis (x/y/z to x/-z/y)',[1 0 0;0 0 1;0 -1 0];
            '90-degree rotation around x-axis (x/y/z to x/z/-y)',[1 0 0;0 0 -1;0 1 0];
            '90-degree rotation around y-axis (x/y/z to -z/y/x)',[0 0 1;0 1 0;-1 0 0];
            '90-degree rotation around y-axis (x/y/z to z/y/-x)',[0 0 -1;0 1 0;1 0 0];
            '90-degree rotation around z-axis (x/y/z to y/-x/z)',[0 -1 0;1 0 0;0 0 1];
            '90-degree rotation around z-axis (x/y/z to -y/x/z)',[0 1 0;-1 0 0;0 0 1];
            '180-degree rotation around x-axis (x/y/z to x/-y/-z)',[1 0 0;0 -1 0;0 0 -1];
            '180-degree rotation around y-axis (x/y/z to -x/y/-z)',[-1 0 0;0 1 0;0 0 -1];
            '180-degree rotation around z-axis (x/y/z to -x/-y/z)',[-1 0 0;0 -1 0;0 0 1];
            'clockwise rotation around x-axis (arbitrary angle)',@(a)[1 0 0;0 cos(a) sin(a);0 -sin(a) cos(a)];
            'clockwise rotation around y-axis (arbitrary angle)',@(a)[cos(a) 0 sin(a);0 1 0;-sin(a) 0 cos(a)];
            'clockwise rotation around z-axis (arbitrary angle)',@(a)[cos(a) sin(a) 0;-sin(a) cos(a) 0;0 0 1];
            'non-rigid reflection along x-axis (x/y/z/ to -x/y/z)', [-1 0 0;0 1 0;0 0 1];
            'non-rigid reflection along y-axis (x/y/z/ to x/-y/z)', [1 0 0;0 -1 0;0 0 1];
            'non-rigid reflection along z-axis (x/y/z/ to x/y/-z)', [1 0 0;0 1 0;0 0 -1];
            'arbitrary affine transformation matrix (manually define 4x4 matrix)', 1;
            'arbitrary affine transformation matrix (load 4x4 matrix from file)', 2};
        for ntime=1:ntimes
            if ntimes>1 [treorient,tok] = listdlg('PromptString',sprintf('Select re-orientation transformation for STEP %d/%d:',ntime,ntimes),'ListSize',[300 200],'SelectionMode','single','ListString',opts(:,1));
            else [treorient,tok] = listdlg('PromptString','Select re-orientation transformation:','ListSize',[300 200],'SelectionMode','single','ListString',opts(:,1));
            end
            if isempty(treorient), return; end
            reorient{ntime}=opts{treorient,2};
            if isequal(reorient{ntime},1)
                answ=inputdlg('Enter affine transformation matrix (4x4 values)','conn_setup_preproc',1,{mat2str(eye(4))});
                if isempty(answ), return; end
                answ=str2num(answ{1});
                reorient{ntime}=answ;
            elseif isequal(reorient{ntime},2)
                [tfilename1,tfilename2]=uigetfile('*.mat','Select file',pwd);
                if ~ischar(tfilename1), return; end
                filename=fullfile(tfilename2,tfilename1);
                reorient{ntime}=filename;
            elseif isequal(reorient{ntime},3)
                [tfilename1,tfilename2]=uigetfile('*.nii','Select file',pwd);
                if ~ischar(tfilename1), return; end
                filename=fullfile(tfilename2,tfilename1);
                reorient{ntime}=filename;
            elseif isa(reorient{ntime},'function_handle'),
                answ=inputdlg('Angular rotation (in degrees)','conn_setup_preproc',1,{num2str(90)});
                if isempty(answ), return; end
                answ=str2num(answ{1});
                reorient{ntime}=reorient{ntime}(answ/180*pi);
            end
        end
    end
end

if any(ismember('functional_art',lSTEPS))
    if ~isempty(art_thresholds)
        art_global_threshold=art_thresholds(1);
        art_motion_threshold=art_thresholds(2);
        if numel(art_thresholds)>=3, art_use_diff_global=art_thresholds(3); end
        if numel(art_thresholds)>=4, art_use_diff_motion=art_thresholds(4); end
        if numel(art_thresholds)>=5, art_use_norms=art_thresholds(5); end
        if numel(art_thresholds)>=6, art_force_interactive=art_thresholds(6); end
        if numel(art_thresholds)>=7&&~isnan(art_thresholds(7)), art_motion_threshold(2)=art_thresholds(7); end
        if numel(art_thresholds)>=8, art_drop_flag=art_thresholds(8); end
    end
    if isempty(art_thresholds)||dogui
        thfig=figure('units','norm','position',[.4,.4,.3,.4],'color',1*[1 1 1],'name','Functional outlier detection settings','numbertitle','off','menubar','none');
        ht0=uicontrol('style','popupmenu','units','norm','position',[.05,.8,.9,.1],'string',{'Use liberal settings (99th percentiles in normative sample)','Use intermediate settings (97th percentiles in normative sample)','Use conservative settings (95th percentiles in normative sample)','Edit settings','Edit settings interactively (ART gui)'},'value',1+art_useconservative,'backgroundcolor',1*[1 1 1]);
        ht1a=uicontrol('style','text','units','norm','position',[.05,.7,.9,.05],'string','Global-signal z-value threshold','backgroundcolor',1*[1 1 1]);
        ht1=uicontrol('style','edit','units','norm','position',[.05,.6,.9,.1],'string',num2str(art_global_threshold));
        ht2a=uicontrol('style','text','units','norm','position',[.05,.5,.9,.05],'string','Subject-motion mm threshold','backgroundcolor',1*[1 1 1]);
        ht2=uicontrol('style','edit','units','norm','position',[.05,.4,.9,.1],'string',num2str(art_motion_threshold));
        ht3a=uicontrol('style','checkbox','units','norm','position',[.05,.3,.4,.05],'string','Use diff global','value',art_use_diff_global,'backgroundcolor',1*[1 1 1],'tooltipstring','Global-signal threshold based on scan-to-scan changes in global BOLD signal');
        ht3b=uicontrol('style','checkbox','units','norm','position',[.05,.25,.4,.05],'string','Use abs global','value',~art_use_diff_global,'backgroundcolor',1*[1 1 1],'tooltipstring','Global-signal threshold based on absolute global BOLD signal values');
        ht3c=uicontrol('style','checkbox','units','norm','position',[.05,.20,.4,.05],'string','Drop first scan(s)','value',art_drop_flag>0,'backgroundcolor',1*[1 1 1],'userdata',art_drop_flag,'tooltipstring','Flags first scan(s) in each session for removal');
        ht4a=uicontrol('style','checkbox','units','norm','position',[.55,.3,.4,.05],'string','Use diff motion','value',art_use_diff_motion,'backgroundcolor',1*[1 1 1],'tooltipstring','Subject-motion threshold based on scan-to-scan changes in motion parameters');
        ht4b=uicontrol('style','checkbox','units','norm','position',[.55,.25,.4,.05],'string','Use abs motion','value',~art_use_diff_motion,'backgroundcolor',1*[1 1 1],'tooltipstring','Subject-motion threshold based on absolute motion parameter values');
        ht5=uicontrol('style','checkbox','units','norm','position',[.55,.2,.9,.05],'string','Use comp motion','value',art_use_norms,'backgroundcolor',1*[1 1 1],'tooltipstring','Subject-motion threshold based on composite motion measure');
        uicontrol('style','pushbutton','string','OK','units','norm','position',[.1,.01,.38,.10],'callback','uiresume');
        uicontrol('style','pushbutton','string','Cancel','units','norm','position',[.51,.01,.38,.10],'callback','delete(gcbf)');
        set([ht1a ht1 ht2a ht2 ht3a ht4a ht3b ht3c ht4b ht5],'enable','off');
        set(ht0,'callback','h=get(gcbo,''userdata''); switch get(gcbo,''value''), case 1, set(h.handles,''enable'',''off''); set(h.handles(2),''string'',num2str(h.default{1}(1))); set(h.handles(4),''string'',num2str(h.default{2}(1))); set(h.handles([5:6 10]),''value'',1); set(h.handles([7 9]),''value'',0); case 2, set(h.handles,''enable'',''off''); set(h.handles(2),''string'',num2str(h.default{1}(2))); set(h.handles(4),''string'',num2str(h.default{2}(2))); set(h.handles([5:6 10]),''value'',1); set(h.handles([7 9]),''value'',0); case 3, set(h.handles,''enable'',''off''); set(h.handles(2),''string'',num2str(h.default{1}(3))); set(h.handles(4),''string'',num2str(h.default{2}(3))); set(h.handles([5:6 10]),''value'',1); set(h.handles([7 9]),''value'',0); case 4, set(h.handles,''enable'',''on''); case 5, set(h.handles,''enable'',''off''); end;','userdata',struct('handles',[ht1a ht1 ht2a ht2 ht3a ht4a ht3b ht3c ht4b ht5],'default',{{art_global_thresholds, art_motion_thresholds}}));
        %@(varargin)set([ht1a ht1 ht2a ht2 ht3a ht4a ht3b ht4b ht5],'enable',subsref({'on','off'},struct('type','{}','subs',{{1+(get(gcbo,'value')~=3)}}))));
        set(ht5,'callback','h=get(gcbo,''userdata''); temp=str2num(get(h.handles(4),''string'')); if get(gcbo,''value''), set(h.handles(3),''string'',''Subject-motion mm threshold''); temp=temp(1); else, set(h.handles(3),''string'',''Subject-motion translation/rotation thresholds [mm, rad]''); if numel(temp)<2, temp=[temp .02]; end; end; set(h.handles(4),''string'',mat2str(temp));','userdata',struct('handles',[ht1a ht1 ht2a ht2 ht3a ht4a ht3b ht3c ht4b ht5],'default',{{art_global_thresholds, art_motion_thresholds}}));
        set(ht3a,'callback',@(varargin)set(ht3b,'value',~get(gcbo,'value')));
        set(ht3b,'callback',@(varargin)set(ht3a,'value',~get(gcbo,'value')));
        set(ht3c,'callback','v=get(gcbo,''value''); if v, v=str2double(inputdlg({''Number of initial scans to remove''},'''',1,{num2str(get(gcbo,''userdata''))})); if isempty(v), v=0; end; end; set(gcbo,''value'',v>0); if v>0, set(gcbo,''userdata'',v); end');
        set(ht4a,'callback',@(varargin)set(ht4b,'value',~get(gcbo,'value')));
        set(ht4b,'callback',@(varargin)set(ht4a,'value',~get(gcbo,'value')));
        uiwait(thfig);
        if ~ishandle(thfig), return; end
        art_global_threshold=str2num(get(ht1,'string'));
        temp=str2num(get(ht2,'string'));
        art_motion_threshold=temp;
        art_use_diff_global=get(ht3a,'value');
        art_use_diff_motion=get(ht4a,'value');
        art_use_norms=get(ht5,'value');
        if get(ht3c,'value'), art_drop_flag=get(ht3c,'userdata'); else art_drop_flag=0; end
        art_force_interactive=get(ht0,'value')==5;
        delete(thfig);
        drawnow;
        if numel(art_motion_threshold)<2, art_thresholds=[art_global_threshold(1) art_motion_threshold(1) art_use_diff_global(1) art_use_diff_motion(1) art_use_norms(1) art_force_interactive(1) nan art_drop_flag(1)];
        else art_thresholds=[art_global_threshold(1) art_motion_threshold(1) art_use_diff_global(1) art_use_diff_motion(1) art_use_norms(1) art_force_interactive(1) art_motion_threshold(2) art_drop_flag(1)];
        end
        %answ=inputdlg({'Enter scan-to-scan global signal z-value threshold','Enter scan-to-scan composite motion mm threshold'},'conn_setup_preproc',1,{num2str(art_global_threshold),num2str(art_motion_threshold)});
        %if isempty(answ), return; end
        %art_global_threshold=str2num(answ{1});
        %art_motion_threshold=str2num(answ{2});
    end
end

if any(ismember({'structural_manualspatialdef','functional_manualspatialdef'},lSTEPS))
    if isempty(respatialdef)||dogui
        DOSPM12=~PREFERSPM8OVERSPM12&spmver12; %SPM12/SPM8
        ntimes=sum(ismember(lSTEPS,{'structural_manualspatialdef','functional_manualspatialdef'}));
        if isempty(respatialdef), respatialdef={};
        elseif ischar(respatialdef), respatialdef=cellstr(respatialdef);
        end
        for ntime=1:ntimes
            if numel(respatialdef)>=ntime, topt=respatialdef(ntime);
            else topt={};
            end
            if DOSPM12,
                str='Select deformation field volume (e.g. y_*.nii)'; conn_disp(str);
                [tfilename1,tfilename2]=uigetfile('*.nii',str,topt{:});
            else
                str='Select transformation file (e.g. *_seg_sn.mat)'; conn_disp(str);
                [tfilename1,tfilename2]=uigetfile('*.mat',str,topt{:});
            end
            if ~ischar(tfilename1), return; end
            filename=fullfile(tfilename2,tfilename1);
            respatialdef{ntime}=filename;
        end
    end
end

if dogui&&any(ismember(lSTEPS,{'structural_normalize','structural_normalize_preservemasks','structural_segment&normalize','structural_segment','functional_normalize','functional_segment&normalize','functional_segment','functional_segment&normalize_direct','functional_segment&normalize_indirect','functional_normalize_indirect','functional_normalize_indirect_preservemasks','functional_normalize_direct','structural_manualspatialdef','functional_manualspatialdef'}))
    DOSPM12=~PREFERSPM8OVERSPM12&spmver12; %SPM12/SPM8
    thfig=figure('units','norm','position',[.4,.4,.3,.2],'color',1*[1 1 1],'name','Segment/Normalize/Resample settings','numbertitle','off','menubar','none');
    if DOSPM12||any(ismember(lSTEPS,{'structural_segment','structural_segment&normalize','functional_segment','functional_segment&normalize'}))
        ht3=uicontrol('style','checkbox','units','norm','position',[.1,.75,.8,.1],'string','Use default Tissue Probability Maps','value',1,'backgroundcolor',1*[1 1 1],'tooltipstring','defines TPM file used by normalization/segmentation routine');
        ht4=[];ht5=[];
    else
        ht3=[];
        ht4=uicontrol('style','checkbox','units','norm','position',[.1,.85,.8,.1],'string','Use default structural template','value',1,'backgroundcolor',1*[1 1 1]);
        ht5=uicontrol('style','checkbox','units','norm','position',[.1,.75,.8,.1],'string','Use default functional template','value',1,'backgroundcolor',1*[1 1 1]);
    end
    ht1a=uicontrol('style','text','units','norm','position',[.1,.55,.6,.1],'string','Structurals target resolution (in mm)','horizontalalignment','left','backgroundcolor',1*[1 1 1]);
    ht1=uicontrol('style','edit','units','norm','position',[.7,.55,.2,.1],'string',num2str(voxelsize_anat),'tooltipstring','defines voxel-size of volumes created when resampling the structural volumes to the desired target space (e.g. MNI)');
    ht2a=uicontrol('style','text','units','norm','position',[.1,.45,.6,.1],'string','Functionals target resolution (in mm)','horizontalalignment','left','backgroundcolor',1*[1 1 1]);
    ht2=uicontrol('style','edit','units','norm','position',[.7,.45,.2,.1],'string',num2str(voxelsize_func),'tooltipstring','defines voxel-size of volumes created when resampling the functional volumes to the desired target space (e.g. MNI)');
    ht0a=uicontrol('style','text','units','norm','position',[.1,.35,.6,.1],'string','Bounding box (in mm)','horizontalalignment','left','backgroundcolor',1*[1 1 1]);
    ht0=uicontrol('style','edit','units','norm','position',[.7,.35,.2,.1],'string',mat2str(boundingbox),'tooltipstring','<HTML>defines bounding box of resampled volumes<br/> - enter a 2x3 matrix with minimum xyz values in the top row and maximum xyz values in the bottom row</HTML>');
    uicontrol('style','pushbutton','string','OK','units','norm','position',[.1,.01,.38,.15],'callback','uiresume');
    uicontrol('style','pushbutton','string','Cancel','units','norm','position',[.51,.01,.38,.15],'callback','delete(gcbf)');
    if ~any(ismember(lSTEPS,{'structural_normalize','structural_normalize_preservemasks','structural_segment&normalize','structural_segment','functional_segment&normalize_indirect','functional_normalize_indirect','functional_normalize_indirect_masks','structural_manualspatialdef'})), set([ht1a ht1 ht4],'enable','off'); end
    if ~any(ismember(lSTEPS,{'functional_normalize','functional_segment&normalize','functional_segment','functional_segment&normalize_direct','functional_normalize_direct','functional_segment&normalize_indirect','functional_normalize_indirect','functional_normalize_indirect_preservemasks','functional_manualspatialdef'})), set([ht2a ht2 ht5],'enable','off'); end
    if all(ismember(lSTEPS,{'structural_manualspatialdef','functional_manualspatialdef'})), set([ht3 ht4 ht5],'visible','off'); set(thfig,'name','Resample settings'); end
    set(ht3,'userdata',[],'callback','if ~get(gcbo,''value''), [t1,t0]=uigetfile(''*.nii;*.img'',''Select TPM file''); if ischar(t1), set(gco,''userdata'',fullfile(t0,t1)); else set(gcbo,''value'',1); end; end');
    set(ht4,'userdata',[],'callback','if ~get(gcbo,''value''), [t1,t0]=uigetfile(''*.nii;*.img'',''Select template file''); if ischar(t1), set(gco,''userdata'',fullfile(t0,t1)); else set(gcbo,''value'',1); end; end');
    set(ht5,'userdata',[],'callback','if ~get(gcbo,''value''), [t1,t0]=uigetfile(''*.nii;*.img'',''Select template file''); if ischar(t1), set(gco,''userdata'',fullfile(t0,t1)); else set(gcbo,''value'',1); end; end');
    uiwait(thfig);
    if ~ishandle(thfig), return; end
    temp=str2num(get(ht1,'string')); if ~isempty(temp), voxelsize_anat=temp; end
    temp=str2num(get(ht2,'string')); if ~isempty(temp), voxelsize_func=temp; end
    temp=str2num(get(ht0,'string')); if ~isempty(temp), boundingbox=temp; end
    if ~isempty(ht3), val=get(ht3,'value'); if ~val, tpm_template=get(ht3,'userdata'); end; end
    if ~isempty(ht4), val=get(ht4,'value'); if ~val, structural_template=get(ht4,'userdata'); end; end
    if ~isempty(ht5), val=get(ht5,'value'); if ~val, functional_template=get(ht5,'userdata'); end; end
    delete(thfig);
    drawnow;
end

if any(ismember('functional_smooth',lSTEPS))||any(ismember('functional_smooth_masked',lSTEPS))
    if isempty(fwhm)||dogui
        if isempty(fwhm), fwhm=8; end
        fwhm=inputdlg('Enter smoothing kernel FWHM (in mm)','conn_setup_preproc',1,{num2str(fwhm)});
        if isempty(fwhm), return; end
        fwhm=str2num(fwhm{1});
    end
end

if any(ismember('functional_surface_smooth',lSTEPS))
    if isempty(diffusionsteps)||dogui
        if isempty(diffusionsteps), diffusionsteps=40; end
        diffusionsteps=inputdlg('Enter number of diffusion steps for smoothing','conn_setup_preproc',1,{num2str(diffusionsteps)});
        if isempty(diffusionsteps), return; end
        diffusionsteps=str2num(diffusionsteps{1});
    end
end

if any(ismember('functional_label',lSTEPS))
    if isempty(label)||dogui
        nl=sum(ismember(lSTEPS,'functional_label'));
        if nl>1,
            if numel(label)~=nl, label=arrayfun(@(n)sprintf('Label%d',n),1:nl,'uni',0); end
            label=inputdlg(repmat({'Enter functional label'},1,nl),'conn_setup_preproc',1,label);
        else
            if isempty(label), label={datestr(now)}; end
            label=inputdlg('Enter functional label  (arbitrary description)','conn_setup_preproc',1,label);
        end
        if isempty(label), return; end
    end
end

if any(ismember('functional_load',lSTEPS))
    if isempty(load_label)||dogui
        nl=sum(ismember(lSTEPS,'functional_load'));
        if isempty(load_label), load_label=cell(1,nl); end
        for il=1:nl
            str=[{'Load from primary dataset'}, arrayfun(@(n)sprintf('Load from secondary dataset #%d %s',n,regexprep(CONN_x.Setup.secondarydataset(n).label,'(.+)','($1)')),1:numel(CONN_x.Setup.secondarydataset),'uni',0)];
            if isempty(load_label{il}), load_label{il}=0;
            elseif ischar(load_label{il}), load_label{il}=conn_datasetlabel(load_label{il},'error');
            end
            [jl,tok] = listdlg('PromptString',['Select secondary dataset ',num2str(il),'/',num2str(nl)],'SelectionMode','single','InitialValue',load_label{il}+1,'ListSize',[400 200],'ListString',str);
            if isempty(jl), return; end
            load_label{il}=jl-1;
        end
    end
end

if any(cellfun('length',regexp(lSTEPS,'^functional_label_as_'))&~ismember(lSTEPS,{'functional_label_as_original', 'functional_label_as_subjectspace', 'functional_label_as_mnispace', 'functional_label_as_surfacespace', 'functional_label_as_smoothed'}))
    labelsnewidx=find(cellfun('length',regexp(lSTEPS,'^functional_label_as_'))>0&~ismember(lSTEPS,{'functional_label_as_original', 'functional_label_as_subjectspace', 'functional_label_as_mnispace', 'functional_label_as_surfacespace', 'functional_label_as_smoothed'}));
    labelsoldidx=find(ismember(lSTEPS,'functional_label'));
    labelsnew=regexprep(lSTEPS(labelsnewidx),'^functional_label_as_','');
    lSTEPS(labelsnewidx)=regexprep(lSTEPS(labelsnewidx),'functional_label_as_.*','functional_label');
    STEPS(labelsnewidx)=regexprep(STEPS(labelsnewidx),'functional_label_as_.*','functional_label');
    if isempty(label)
        label=labelsnew;
    else
        label=[label(:)', labelsnew(:)'];
        [nill,idx]=sort([labelsoldidx(:)', labelsnewidx(:)']);
        label=label(idx);
    end
end

if any(cellfun('length',regexp(lSTEPS,'^functional_load_from_'))&~ismember(lSTEPS,{'functional_load_from_original', 'functional_load_from_subjectspace', 'functional_load_from_mnispace', 'functional_load_from_surfacespace', 'functional_load_from_smoothed'}))
    labelsnewidx=find(cellfun('length',regexp(lSTEPS,'^functional_load_from_'))>0&~ismember(lSTEPS,{'functional_load_from_original', 'functional_load_from_subjectspace', 'functional_load_from_mnispace', 'functional_load_from_surfacespace', 'functional_load_from_smoothed'}));
    labelsoldidx=find(ismember(lSTEPS,'functional_load'));
    labelsnew=regexprep(lSTEPS(labelsnewidx),'^functional_load_from_','');
    lSTEPS(labelsnewidx)=regexprep(lSTEPS(labelsnewidx),'functional_load_from_.*','functional_load');
    STEPS(labelsnewidx)=regexprep(STEPS(labelsnewidx),'functional_load_from_.*','functional_load');
    if isempty(load_label)
        load_label=labelsnew;
    else
        load_label=[load_label(:)', labelsnew(:)'];
        [nill,idx]=sort([labelsoldidx(:)', labelsnewidx(:)']);
        load_label=load_label(idx);
    end
end


% loginfo=struct('subjects',subjects,'steps',STEPS,...
%     'fwhm',fwhm,'sliceorder',sliceorder,'ta',ta,'unwarp',unwarp,'removescans',removescans,'applytofunctional',applytofunctional,...
%     'coregtomean',coregtomean,'reorient',reorient,'art_thresholds',art_thresholds,'voxelsize',voxelsize,'boundingbox',boundingbox,...
%     'doimport',doimport,'dogui',0,'functional_template',functional_template,'structural_template',structural_template,...
%     'tpm_template',tpm_template,'tpm_ngaus',tpm_ngaus);
if parallel_N>0,
    if ~isempty(parallel_profile), conn_jobmanager('setprofile',parallel_profile); end
    conn save;
    if isempty(sliceorder)&&~isempty(sliceorder_select), sliceorder=sliceorder_select_options{sliceorder_select}; end
    info=conn_jobmanager('submit','setup_preprocessing',subjects,parallel_N,[],...
        STEPS,...
        'sessions',sessions,'sets',sets,'fwhm',fwhm,'label',label,'load_label',load_label,'sliceorder',sliceorder,'ta',ta,'unwarp',unwarp,'removescans',removescans,'applytofunctional',applytofunctional,...
        'coregtomean',coregtomean,'rtm',rtm,'coregsource',coregsource,'reorient',reorient,'respatialdef',respatialdef,'art_thresholds',art_thresholds,'voxelsize_anat',voxelsize_anat,'voxelsize_func',voxelsize_func,'boundingbox',boundingbox,'interp',interp,'diffusionsteps',diffusionsteps,...
        'doimport',doimport,'dogui',0,'functional_template',functional_template,'structural_template',structural_template,...
        'affreg',affreg,'tpm_template',tpm_template,'tpm_ngaus',tpm_ngaus,'vdm_et1',vdm_et1,'vdm_et2',vdm_et2,'vdm_ert',vdm_ert,'vdm_blip',vdm_blip,'vdm_type',vdm_type,'bp_filter',bp_filter,'bp_keep0',bp_keep0,'reg_names',reg_names,'reg_dimensions',reg_dimensions,'reg_deriv',reg_deriv,'reg_filter',reg_filter,'reg_detrend',reg_detrend,'reg_lag',reg_lag,'reg_lagmax',reg_lagmax,'reg_skip',reg_skip);
    if isequal(parallel_profile,find(strcmp('Null profile',conn_jobmanager('profiles')))),
        ok=0;
    elseif dogui
        conn_jobmanager(info);
        ok=0;
    else
        [nill,finished]=conn_jobmanager('waitfor',info);
        if finished==2, ok=1+doimport;
        else ok=3;
        end
    end
    return;
elseif conn_projectmanager('inserver'), 
    if isempty(sliceorder)&&~isempty(sliceorder_select), sliceorder=sliceorder_select_options{sliceorder_select}; end
    conn_process('setup_preprocessing',...
        STEPS,...
        'sessions',sessions,'sets',sets,'fwhm',fwhm,'label',label,'load_label',load_label,'sliceorder',sliceorder,'ta',ta,'unwarp',unwarp,'removescans',removescans,'applytofunctional',applytofunctional,...
        'coregtomean',coregtomean,'rtm',rtm,'coregsource',coregsource,'reorient',reorient,'respatialdef',respatialdef,'art_thresholds',art_thresholds,'voxelsize_anat',voxelsize_anat,'voxelsize_func',voxelsize_func,'boundingbox',boundingbox,'interp',interp,'diffusionsteps',diffusionsteps,...
        'doimport',doimport,'dogui',0,'functional_template',functional_template,'structural_template',structural_template,...
        'affreg',affreg,'tpm_template',tpm_template,'tpm_ngaus',tpm_ngaus,'vdm_et1',vdm_et1,'vdm_et2',vdm_et2,'vdm_ert',vdm_ert,'vdm_blip',vdm_blip,'vdm_type',vdm_type,'bp_filter',bp_filter,'bp_keep0',bp_keep0,'reg_names',reg_names,'reg_dimensions',reg_dimensions,'reg_deriv',reg_deriv,'reg_filter',reg_filter,'reg_detrend',reg_detrend,'reg_lag',reg_lag,'reg_lagmax',reg_lagmax,'reg_skip',reg_skip);
    return;
else
    if ~isfield(CONN_x,'SetupPreproc')||~isfield(CONN_x.SetupPreproc,'log'), CONN_x.SetupPreproc.log={}; end
    try, spmver=spm('version'); catch, spmver=spm('ver'); end
    CONN_x.SetupPreproc.log{end+1}={'timestamp',datestr(now),'ver_CONN',conn('ver'),'ver_SPM',spmver,'steps',...
        STEPS,...
        'subjects',subjects,'sessions',sessions,'sets',sets,'fwhm',fwhm,'label',label,'load_label',load_label,'sliceorder',sliceorder,'sliceorder_select',sliceorder_select,'ta',ta,'unwarp',unwarp,'removescans',removescans,'applytofunctional',applytofunctional,...
        'coregtomean',coregtomean,'rtm',rtm,'coregsource',coregsource,'reorient',reorient,'respatialdef',respatialdef,'art_thresholds',art_thresholds,'voxelsize_anat',voxelsize_anat,'voxelsize_func',voxelsize_func,'boundingbox',boundingbox,'interp',interp,'diffusionsteps',diffusionsteps,...
        'doimport',doimport,'dogui',0,'functional_template',functional_template,'structural_template',structural_template,...
        'affreg',affreg,'tpm_template',tpm_template,'tpm_ngaus',tpm_ngaus,'vdm_et1',vdm_et1,'vdm_et2',vdm_et2,'vdm_ert',vdm_ert,'vdm_blip',vdm_blip,'vdm_type',vdm_type,'bp_filter',bp_filter,'bp_keep0',bp_keep0,'reg_names',reg_names,'reg_dimensions',reg_dimensions,'reg_deriv',reg_deriv,'reg_filter',reg_filter,'reg_detrend',reg_detrend,'reg_lag',reg_lag,'reg_lagmax',reg_lagmax,'reg_skip',reg_skip};
end
job_id={};

for iSTEP=1:numel(STEPS)
    matlabbatch={};
    outputfiles={};
    STEP=STEPS{iSTEP};
    idx=find(strcmpi(regexprep(lower(STEP),'^run_|^update_|^interactive_',''),steps));
    if ~isempty(idx), STEP_name=steps_names{idx(1)};
    else STEP_name='process';
    end
    ok=0;
    
    hmsg=[];
    if dogui, hmsg=conn_msgbox({['Preparing ',STEP_name],'Please wait...'},'');
    else conn_disp(['Preparing ',STEP_name,'. Please wait...']);
    end
    switch(regexprep(lower(STEP),'^run_|^update_|^interactive_',''))
        case 'functional_removescans'
            if removescans>0, conn_disp('fprintf','removing first %d scans\n',removescans);
            elseif removescans<0, conn_disp('fprintf','removing last %d scans\n',-removescans);
            end
            for isubject=1:numel(subjects),
                nsubject=subjects(isubject);
                %matlabbatch{end+1}.removescans.data={};
                nsess=CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject));
                for nses=1:nsess,
                    if ismember(nses,sessions)
                        filename=conn_get_functional(nsubject,nses,sets);
                        if isempty(filename), error('Functional data not yet defined for subject %d session %d',nsubject,nses); end
                        temp=cellstr(filename);
                        if numel(temp)==1,
                            temp=cellstr(conn_expandframe(temp{1}));
                        end
                        %matlabbatch{end}.removescans.data{end+1}=temp;
                        outputfiles{isubject}{nses}{1}=char(temp(max(0,removescans)+1:end+min(0,removescans)));
                        if ~sets||ALLSETSPERMISSIONS,
                            nl1covariates=length(CONN_x.Setup.l1covariates.names)-1;
                            for nl1covariate=1:nl1covariates,
                                try
                                    covfilename=CONN_x.Setup.l1covariates.files{nsubject}{nl1covariate}{nses}{1};
                                    switch(covfilename),
                                        case '[raw values]',
                                            data=CONN_x.Setup.l1covariates.files{nsubject}{nl1covariate}{nses}{3};
                                        otherwise,
                                            data=conn_loadtextfile(covfilename,false);
                                            %if isstruct(data), tempnames=fieldnames(data); data=data.(tempnames{1}); end
                                    end
                                    if size(data,1)==numel(temp)
                                        data=data(max(0,removescans)+1:end+min(0,removescans),:);
                                        outputfiles{isubject}{nses}{1+nl1covariate}=data;
                                    end
                                catch, conn_disp('fprintf','warning: problem cropping covariate %s subject %d session %d\n',CONN_x.Setup.l1covariates.names{nl1covariate},nsubject,nses);
                                end
                            end
                        end
                    end
                end
            end
            
        case 'functional_bandpass'
            conn_disp('fprintf','bandpass filtering (%.4f-%.4f Hz)\n',bp_filter(1),bp_filter(2));
            if bp_keep0, conn_disp('fprintf','maintaining average BOLD signal\n'); end
            for isubject=1:numel(subjects),
                nsubject=subjects(isubject);
                nsess=CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject));
                for nses=1:nsess,
                    if ismember(nses,sessions)
                        filename=conn_get_functional(nsubject,nses,sets);
                        if isempty(filename), error('Functional data not yet defined for subject %d session %d',nsubject,nses); end
                        filein=cellstr(filename);
                        fileout=conn_prepend('b',filein);
                        if numel(fileout)>1, fileout=fileout(1); end
                        Vin=spm_vol(char(filein));
                        Vout=Vin;
                        for nt=1:numel(Vout), Vout(nt).fname=char(fileout); Vout(nt).pinfo=[1;0;0]; Vout(nt).descrip='band-pass filtered'; end
                        %Vout=struct('fname',char(fileout),...
                        %    'mat',Vin(1).mat,...
                        %    'dim',Vin(1).dim,...
                        %    'n',[1,1],...
                        %    'pinfo',[1;0;0],...
                        %    'dt',[spm_type('float32'),spm_platform('bigend')],...
                        %    'descrip','band-pass filtered');
                        %Vout=repmat(Vout,[numel(Vin),1]);for nt=1:numel(Vin),Vout(nt).n=[nt,1];end
                        Vout=spm_create_vol(Vout);
                        tr=conn_get_rt(nsubject,nses,sets);
                        [gridx,gridy]=ndgrid(1:Vin(1).dim(1),1:Vin(1).dim(2));
                        xyz0=[gridx(:),gridy(:)]';
                        for slice=1:Vin(1).dim(3)
                            xyz=[xyz0; slice+zeros(1,size(xyz0,2)); ones(1,size(xyz0,2))];
                            y=spm_get_data(Vin(:)',xyz);
                            maskout=isnan(y);
                            y(maskout)=0;
                            my=sum(y,1)./max(1,sum(~maskout,1));
                            y(maskout)=my(ceil(find(maskout)./size(y,1)));
                            y=conn_filter(tr,bp_filter,y);
                            if bp_keep0, y=y+repmat(my,size(y,1),1); end
                            y=permute(reshape(y,[numel(Vin),Vin(1).dim(1:2)]),[2,3,1]);
                            for nt=1:numel(Vin),
                                t=y(:,:,nt);
                                Vout(nt)=spm_write_plane(Vout(nt),t,slice);
                            end
                        end
                        outputfiles{isubject}{nses}{1}=char(fileout);
                    end
                end
            end
            
        case 'functional_regression'
            if isempty(bp_filter)&&any(reg_filter), error('Band-pass filter not specified (missing #bp_filter field)'); end
            conn_disp('fprintf','regression of temporal component (%s)\n',sprintf('%s ',reg_names{:}));
            conn_disp('fprintf','dimensions = %s; derivatives = %s; filtered = %s; lags(0-%ss) = %s; skip = %d; detrend = %d\n',mat2str(reg_dimensions),mat2str(reg_deriv),mat2str(reg_filter),mat2str(reg_lagmax),mat2str(reg_lag),reg_skip,reg_detrend);
            for isubject=1:numel(subjects),
                nsubject=subjects(isubject);
                nsess=CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject));
                for nses=1:nsess,
                    if ismember(nses,sessions)
                        filename=conn_get_functional(nsubject,nses,sets);
                        if isempty(filename), error('Functional data not yet defined for subject %d session %d',nsubject,nses); end
                        filein=cellstr(filename);
                        fileout=conn_prepend('d',filein);
                        if numel(fileout)>1, fileout=fileout(1); end
                        Vin=spm_vol(char(filein));
                        Vout=Vin;
                        for nt=1:numel(Vout), Vout(nt).fname=char(fileout); Vout(nt).pinfo=[1;0;0]; Vout(nt).descrip='linear regressed'; end
                        if any(reg_lag), Vlag=struct('fname',conn_prepend('lag_',filein{1}),'mat',Vin(1).mat,'dim',Vin(1).dim,'n',[1,1],'pinfo',[1;0;0],'dt',[spm_type('float32'),spm_platform('bigend')],'descrip','lag (samples)'); end % note: lag in target wrt regressors
                        if ~reg_skip,
                            Vout=spm_create_vol(Vout);
                            if any(reg_lag), Vlag=spm_create_vol(Vlag); end
                        end
                        lagidx=[];
                        lagmax=[];
                        if any(reg_filter), rt=conn_get_rt(nsubject,nses,sets); end
                        X=[ones(numel(Vin),1)]; Xnames={'session (1)'};
                        if reg_detrend, X=[X,linspace(-1,1,numel(Vin))']; Xnames{end+1}='detrend (1)'; end
                        entercovariates=X;
                        reg_done=false(size(reg_names));
                        for nl1covariate=1:numel(reg_names)
                            icov=find(strcmp(CONN_x.Setup.l1covariates.names(1:end-1),reg_names{nl1covariate}));
                            if ~isempty(icov) % first-level covariates
                                reg_done(nl1covariate)=true;
                                cfilename=CONN_x.Setup.l1covariates.files{nsubject}{icov}{nses}{1};
                                assert(~isempty(cfilename),'covariate %s has not been defined for subject %d sessions %d',reg_names{nl1covariate},nsubject,nses);
                                switch(cfilename),
                                    case '[raw values]',
                                        data=CONN_x.Setup.l1covariates.files{nsubject}{icov}{nses}{3};
                                    otherwise,
                                        data=conn_loadtextfile(cfilename,false);
                                end
                                assert(size(data,1)==numel(Vin),'mismatched dimensions; functional data has %d timepoints, covariate %s has %d timepoints',numel(Vin),reg_names{nl1covariate},size(data,1));
                                if numel(reg_dimensions)>=nl1covariate, data=data(:,1:min(size(data,2),reg_dimensions(nl1covariate))); end
                                entercovariates=cat(2,entercovariates,data-repmat(mean(data,1),size(data,1),1)); % note: same as ROI entercovariates behavior                                
                                if numel(reg_deriv)>=nl1covariate&&reg_deriv(nl1covariate)>0, ddata=convn(cat(1,data(1,:),data,data(end,:)),[1;0;-1],'valid'); data=[data, ddata]; end
                                if numel(reg_deriv)>=nl1covariate&&reg_deriv(nl1covariate)>1, data=[data, convn(cat(1,ddata(1,:),ddata,ddata(end,:)),[1;0;-1],'valid')]; end
                                if numel(reg_filter)>=nl1covariate&&reg_filter(nl1covariate)>0, data=conn_filter(rt,bp_filter,data); end
                                if numel(reg_lag)>=nl1covariate&&reg_lag(nl1covariate)>0, lagidx=[lagidx, size(X,2)+(1:size(data,2))]; end
                                if nnz(data~=0&data~=1), X=cat(2,X,data-repmat(mean(data,1),size(data,1),1)); Xnames{end+1}=sprintf('%s (%d)',reg_names{nl1covariate},size(data,2)); % note: 0/1 covariates not centered
                                else X=cat(2,X,data); Xnames{end+1}=sprintf('%s (%d)',reg_names{nl1covariate},size(data,2));
                                end
                            end
                        end
                        Vsource=[];
                        for nl1covariate=1:numel(reg_names)
                            icov=find(strcmp(CONN_x.Setup.l1covariates.names(1:end-1),reg_names{nl1covariate}));
                            if isempty(icov) % ROIs
                                reg_done(nl1covariate)=true;
                                nroi=find(strcmp(CONN_x.Setup.rois.names(1:end-1),reg_names{nl1covariate}));
                                assert(~isempty(nroi),'unable to find first-level covariate or ROI named %s',reg_names{nl1covariate});
                                if isempty(Vsource)
                                    Vsource=CONN_x.Setup.functional{nsubject}{nses}{1};
                                    clear VsourceUnsmoothed;
                                    for nalt=1:numel(CONN_x.Setup.secondarydataset)
                                        VsourceUnsmoothed{nalt}=conn_get_functional(nsubject,nses,nalt,true);
                                    end
                                end
                                if CONN_x.Setup.analysisunits==1, scalinglevel='roi'; else scalinglevel='none'; end
                                if (nroi>3&&~CONN_x.Setup.rois.sessionspecific(nroi))||(nroi<=3&&~CONN_x.Setup.structural_sessionspecific), nsesstemp=1; else nsesstemp=nsess; end
                                temp=conn_maskserode(nsubject,nroi,false);
                                Vmask=temp{nsubject}{nroi}{min(nses,nsesstemp)};
                                %Vmask=CONN_x.Setup.rois.files{nsubject}{nroi}{min(nses,nsesstemp)}{1};
                                if CONN_x.Setup.rois.mask(nroi), mask=CONN_x.Setup.rois.files{nsubject}{1}{min(nses,nsesstemp)}{1}; else, mask=''; end
                                if CONN_x.Setup.rois.multiplelabels(nroi)&&numel(CONN_x.Setup.rois.files{nsubject}{nroi}{min(nses,nsesstemp)}{3})<=1, level='clusters'; else, level='rois'; end
                                if CONN_x.Setup.rois.unsmoothedvolumes(nroi), Vsourcethis=VsourceUnsmoothed{CONN_x.Setup.rois.unsmoothedvolumes(nroi)}; else, Vsourcethis=Vsource; end
                                if conn_surf_dimscheck(CONN_x.Setup.rois.files{nsubject}{nroi}{min(nses,nsesstemp)}{3})&&~conn_surf_dimscheck(Vsourcethis), fsanatomical=CONN_x.Setup.structural{nsubject}{min(nses,nsesstemp)}{1}; else fsanatomical=''; end
                                outputtype='none'; filenamerex='';
                                if CONN_x.Setup.rois.dimensions{nroi}>1,        % average&pca
                                    if CONN_x.Setup.rois.weighted(nroi), data=conn_rex(Vsourcethis,Vmask,'summary_measure','weighted eigenvariate','dims',CONN_x.Setup.rois.dimensions{nroi},'conjunction_mask',mask,'level',level,'scaling','none','select_clusters',0,'covariates',entercovariates,'fsanatomical',fsanatomical,'output_type',outputtype,'output_rex',filenamerex);
                                    else data=conn_rex(Vsourcethis,Vmask,'summary_measure','eigenvariate','dims',CONN_x.Setup.rois.dimensions{nroi},'conjunction_mask',mask,'level',level,'scaling',scalinglevel,'select_clusters',0,'covariates',entercovariates,'fsanatomical',fsanatomical,'output_type',outputtype,'output_rex',filenamerex);
                                    end
                                    if size(data,2)<CONN_x.Setup.rois.dimensions{nroi},
                                        conn_disp('fprintf','WARNING: not enough voxels or scans to extract %d dimensions from %s @ %s\n',CONN_x.Setup.rois.dimensions{nroi},Vsourcethis,Vmask);
                                        data=[data zeros(size(data,1),CONN_x.Setup.rois.dimensions{nroi}-size(data,2))];
                                    end
                                elseif CONN_x.Setup.rois.dimensions{nroi}==0||CONN_x.Setup.rois.weighted(nroi),   % weighted sum
                                    data=conn_rex(Vsourcethis,Vmask,'summary_measure','weighted sum','conjunction_mask',mask,'level',level,'scaling','none','select_clusters',0,'covariates',entercovariates,'fsanatomical',fsanatomical,'output_type',outputtype,'output_rex',filenamerex);
                                else                                            % average
                                    data=conn_rex(Vsourcethis,Vmask,'summary_measure','mean','conjunction_mask',mask,'level',level,'scaling',scalinglevel,'select_clusters',0,'covariates',entercovariates,'fsanatomical',fsanatomical,'output_type',outputtype,'output_rex',filenamerex);
                                end
                                [data,ok]=conn_nan(data);
                                assert(size(data,1)==numel(Vin),'mismatched dimensions; functional data has %d timepoints, covariate %s has %d timepoints',numel(Vin),reg_names{nl1covariate},size(data,1));
                                if numel(reg_dimensions)>=nl1covariate, data=data(:,1:min(size(data,2),reg_dimensions(nl1covariate))); end
                                if numel(reg_deriv)>=nl1covariate&&reg_deriv(nl1covariate)>0, ddata=convn(cat(1,data(1,:),data,data(end,:)),[1;0;-1],'valid'); data=[data, ddata]; end
                                if numel(reg_deriv)>=nl1covariate&&reg_deriv(nl1covariate)>1, data=[data, convn(cat(1,ddata(1,:),ddata,ddata(end,:)),[1;0;-1],'valid')]; end
                                if numel(reg_filter)>=nl1covariate&&reg_filter(nl1covariate)>0, data=conn_filter(rt,bp_filter,data); end
                                if numel(reg_lag)>=nl1covariate&&reg_lag(nl1covariate)>0, lagidx=[lagidx, size(X,2)+(1:size(data,2))]; end
                                X=cat(2,X,data-repmat(mean(data,1),size(data,1),1)); Xnames{end+1}=sprintf('%s (%d)',reg_names{nl1covariate},size(data,2));
                            end
                        end
                        if ~reg_skip
                            [gridx,gridy]=ndgrid(1:Vin(1).dim(1),1:Vin(1).dim(2));
                            xyz0=[gridx(:),gridy(:)]';
                            iX=pinv(X'*X);
                            if ~isempty(lagidx), rt=conn_get_rt(nsubject,nses,sets); end
                            for slice=1:Vin(1).dim(3)
                                xyz=[xyz0; slice+zeros(1,size(xyz0,2)); ones(1,size(xyz0,2))];
                                y=spm_get_data(Vin(:)',xyz);
                                maskout=isnan(y);
                                y(maskout)=0;
                                my=sum(y,1)./max(1,sum(~maskout,1));
                                y(maskout)=my(ceil(find(maskout)./size(y,1)));
                                if ~isempty(lagidx)
                                    if isempty(lagmax), lagmax=reg_lagmax/rt; end
                                    if reg_detrend, X(:,lagidx)=X(:,lagidx)-X(:,1:2)*(pinv(X(:,1:2)'*X(:,1:2))*X(:,1:2)'*X(:,lagidx)); end
                                    [yfit,nill,lag]=conn_lagregress(X,y,'select',lagidx,'maxdn',lagmax,'omit',1);
                                    y=y-yfit;
                                    Vlag=spm_write_plane(Vlag,reshape(lag,Vin(1).dim(1:2)),slice);
                                else
                                    b=iX*X'*y;
                                    b(1,:)=0; % note: keep constant term
                                    y=y-X*b;
                                end
                                y=permute(reshape(y,[numel(Vin),Vin(1).dim(1:2)]),[2,3,1]);
                                for nt=1:numel(Vin),
                                    t=y(:,:,nt);
                                    Vout(nt)=spm_write_plane(Vout(nt),t,slice);
                                end
                            end
                        end
                        conn_savetextfile(conn_prepend('dp_',filein{1},'.txt'),X);
                        conn_savematfile(conn_prepend('dp_',filein{1},'.mat'),'X','Xnames');
                        if ~reg_skip, outputfiles{isubject}{nses}{1}=char(fileout);
                        else outputfiles{isubject}{nses}{1}=char(filein);
                        end
                        outputfiles{isubject}{nses}{2}=conn_prepend('dp_',filein{1},'.txt');
                    end
                end
            end
            
        case 'functional_manualorient'
        case 'structural_manualorient'
        case 'functional_center'
        case 'functional_centertostruct'
        case 'structural_center'
        case {'functional_label','functional_label_as_original', 'functional_label_as_subjectspace', 'functional_label_as_mnispace', 'functional_label_as_surfacespace', 'functional_label_as_smoothed'}
        case {'functional_load','functional_load_from_original', 'functional_load_from_subjectspace', 'functional_load_from_mnispace', 'functional_load_from_surfacespace', 'functional_load_from_smoothed'}
            
        case 'functional_manualspatialdef'
            if iscell(respatialdef), trespatialdef=respatialdef{1}; respatialdef=respatialdef(2:end);
            else trespatialdef=respatialdef;
            end
            DOSPM12=~PREFERSPM8OVERSPM12&spmver12; %SPM12/SPM8
            if DOSPM12
                matlabbatch{end+1}.spm.spatial.normalise.write.woptions.bb=boundingbox;
                matlabbatch{end}.spm.spatial.normalise.write.woptions.vox=voxelsize_func.*[1 1 1];
                if ~isempty(interp), matlabbatch{end}.spm.spatial.normalise.write.woptions.interp=interp; end
            else
                matlabbatch{end+1}.spm.spatial.normalise.write.roptions.bb=boundingbox;
                matlabbatch{end}.spm.spatial.normalise.write.roptions.vox=voxelsize_func.*[1 1 1];
                if ~isempty(interp), matlabbatch{end}.spm.spatial.normalise.write.roptions.interp=interp; end
            end
            jsubject=0;
            for isubject=1:numel(subjects), % normalize write
                nsubject=subjects(isubject);
                nsess=intersect(sessions,1:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)));
                for nses=nsess(:)'
                    jsubject=jsubject+1;
                    if DOSPM12, matlabbatch{end}.spm.spatial.normalise.write.subj(jsubject).def={trespatialdef};
                    else        matlabbatch{end}.spm.spatial.normalise.write.subj(jsubject).matname={trespatialdef};
                    end
                    matlabbatch{end}.spm.spatial.normalise.write.subj(jsubject).resample={};
                    filename=conn_get_functional(nsubject,nses,sets);
                    if isempty(filename), error('Functional data not yet defined for subject %d session %d',nsubject,nses); end
                    temp=cellstr(filename);
                    if coregtomean, % keeps mean image in same space in case it is required later
                        [xtemp,failed]=conn_setup_preproc_meanimage(temp{1});
                        if ~isempty(xtemp),
                            xtemp={xtemp};
                            matlabbatch{end}.spm.spatial.normalise.write.subj(jsubject).resample=cat(1,matlabbatch{end}.spm.spatial.normalise.write.subj(jsubject).resample,xtemp);
                        end
                    end
                    matlabbatch{end}.spm.spatial.normalise.write.subj(jsubject).resample=cat(1,matlabbatch{end}.spm.spatial.normalise.write.subj(jsubject).resample,temp);
                    outputfiles{isubject}{nses}{1}=char(conn_prepend('w',temp));
                end
            end
            if ~jsubject, matlabbatch=matlabbatch(1:end-1); end
            
        case 'structural_manualspatialdef'
            if iscell(respatialdef), trespatialdef=respatialdef{1}; respatialdef=respatialdef(2:end);
            else trespatialdef=respatialdef;
            end
            DOSPM12=~PREFERSPM8OVERSPM12&spmver12; %SPM12/SPM8
            if DOSPM12
                matlabbatch{end+1}.spm.spatial.normalise.write.woptions.bb=boundingbox;
                matlabbatch{end}.spm.spatial.normalise.write.woptions.vox=voxelsize_anat.*[1 1 1];
                if ~isempty(interp), matlabbatch{end}.spm.spatial.normalise.write.woptions.interp=interp;
                else matlabbatch{end}.spm.spatial.normalise.write.woptions.interp=1;
                end
            else
                matlabbatch{end+1}.spm.spatial.normalise.write.roptions.bb=boundingbox;
                matlabbatch{end}.spm.spatial.normalise.write.roptions.vox=voxelsize_anat.*[1 1 1];
                if ~isempty(interp), matlabbatch{end}.spm.spatial.normalise.write.roptions.interp=interp;
                else matlabbatch{end}.spm.spatial.normalise.write.roptions.interp=1;
                end
            end
            jsubject=0;
            for isubject=1:numel(subjects), % normalize write
                nsubject=subjects(isubject);
                if CONN_x.Setup.structural_sessionspecific, nsess_struct=intersect(sessions,1:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)));
                else nsess_struct=1;
                end
                for nses=nsess_struct(:)'
                    if isempty(CONN_x.Setup.structural{nsubject}{nses}{1}), error(sprintf('No structural file defined for Subject %d Session %d. Please select a structural file',nsubject,nses));
                    elseif numel(CONN_x.Setup.structural{nsubject}{nses}{3})>1, error(sprintf('Multiple structural files found for Subject %d Session %d. Please select a single structural file',nsubject,nses));
                    end
                    outputfiles{isubject}{nses}{1}=CONN_x.Setup.structural{nsubject}{nses}{1};
                    if ismember(nses,sessions)
                        jsubject=jsubject+1;
                        if DOSPM12, matlabbatch{end}.spm.spatial.normalise.write.subj(jsubject).def={trespatialdef};
                        else        matlabbatch{end}.spm.spatial.normalise.write.subj(jsubject).matname={trespatialdef};
                        end
                        matlabbatch{end}.spm.spatial.normalise.write.subj(jsubject).resample=outputfiles{isubject}{nses}(1);
                    end
                    outputfiles{isubject}{nses}{1}=conn_prepend('w',outputfiles{isubject}{nses}{1});
                end
            end
            if ~jsubject, matlabbatch=matlabbatch(1:end-1); end
            
        case 'structural_segment'
            if ~PREFERSPM8OVERSPM12&&spmver12 %SPM12
                matlabbatch{end+1}.spm.spatial.preproc.channel.vols={};
                jsubject=0;
                for isubject=1:numel(subjects),
                    nsubject=subjects(isubject);
                    if CONN_x.Setup.structural_sessionspecific, nsess=CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)); else nsess=1; end
                    for nses=1:nsess
                        if ismember(nses,sessions)
                            jsubject=jsubject+1;
                            if isempty(CONN_x.Setup.structural{nsubject}{nses}{1}), error(sprintf('No structural file defined for Subject %d Session %d. Please select a structural file',nsubject,nses));
                            elseif numel(CONN_x.Setup.structural{nsubject}{nses}{3})>1, error(sprintf('Multiple structural files found for Subject %d Session %d. Please select a single structural file',nsubject,nses));
                            end
                            matlabbatch{end}.spm.spatial.preproc.channel.vols{jsubject}=CONN_x.Setup.structural{nsubject}{nses}{1};
                            outputfiles{isubject}{nses}{1}=CONN_x.Setup.structural{nsubject}{nses}{1};
                            outputfiles{isubject}{nses}{2}=conn_prepend('c1',CONN_x.Setup.structural{nsubject}{nses}{1},'.nii'); % note: fix SPM12 issue converting .img to .nii
                            outputfiles{isubject}{nses}{3}=conn_prepend('c2',CONN_x.Setup.structural{nsubject}{nses}{1},'.nii');
                            outputfiles{isubject}{nses}{4}=conn_prepend('c3',CONN_x.Setup.structural{nsubject}{nses}{1},'.nii');
                        end
                    end
                end
                if ~isempty(tpm_template), matlabbatch{end}.spm.spatial.preproc.tissue=conn_setup_preproc_tissue(tpm_template,tpm_ngaus,subjects,sessions); end
                if ~isempty(affreg), matlabbatch{end}.spm.spatial.preproc.warp.affreg=affreg; end
                matlabbatch{end}.spm.spatial.preproc.channel.vols=reshape(matlabbatch{end}.spm.spatial.preproc.channel.vols,[],1);
                matlabbatch{end}.spm.spatial.preproc.warp.write=[1 1];
                if ~jsubject, matlabbatch=matlabbatch(1:end-1); end
            else % SPM8
                matlabbatch{end+1}.spm.spatial.preproc.data={};
                jsubject=0;
                for isubject=1:numel(subjects),
                    nsubject=subjects(isubject);
                    if CONN_x.Setup.structural_sessionspecific, nsess=CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)); else nsess=1; end
                    for nses=1:nsess
                        if ismember(nses,sessions)
                            jsubject=jsubject+1;
                            if isempty(CONN_x.Setup.structural{nsubject}{nses}{1}), error(sprintf('No structural file defined for Subject %d Session %d. Please select a structural file',nsubject,nses));
                            elseif numel(CONN_x.Setup.structural{nsubject}{nses}{3})>1, error(sprintf('Multiple structural files found for Subject %d Session %d. Please select a single structural file',nsubject,nses));
                            end
                            matlabbatch{end}.spm.spatial.preproc.data{jsubject}=CONN_x.Setup.structural{nsubject}{nses}{1};
                            outputfiles{isubject}{nses}{1}=CONN_x.Setup.structural{nsubject}{nses}{1};
                            outputfiles{isubject}{nses}{2}=conn_prepend('c1',CONN_x.Setup.structural{nsubject}{nses}{1});
                            outputfiles{isubject}{nses}{3}=conn_prepend('c2',CONN_x.Setup.structural{nsubject}{nses}{1});
                            outputfiles{isubject}{nses}{4}=conn_prepend('c3',CONN_x.Setup.structural{nsubject}{nses}{1});
                        end
                    end
                end
                if ~isempty(tpm_template),
                    if ~isempty(subjects)&&~isempty(sessions)&&(isnumeric(tpm_template)||(ischar(tpm_template)&&size(tpm_template,1)==1&&~isempty(conn_datasetlabel(tpm_template))))
                        error('unsupported subject-specific TPM in SPM8; please upgrade to SPM12 instead')
                    else temp=cellstr(conn_expandframe(tpm_template));
                    end
                    if isempty(tpm_ngaus), tpm_ngaus=[2 2 2 4]; end % grey/white/CSF (+other implicit)
                    matlabbatch{end}.spm.spatial.preproc.opts.tpm=temp;
                    matlabbatch{end}.spm.spatial.preproc.opts.ngaus=ngaus(1:numel(temp)+1);
                end
                matlabbatch{end}.spm.spatial.preproc.roptions.bb=boundingbox; %
                matlabbatch{end}.spm.spatial.preproc.output.GM=[0,0,1];
                matlabbatch{end}.spm.spatial.preproc.output.WM=[0,0,1];
                matlabbatch{end}.spm.spatial.preproc.output.CSF=[0,0,1];
                if ~jsubject, matlabbatch=matlabbatch(1:end-1); end
            end
            jsubject=0;
            for isubject=1:numel(subjects),
                nsubject=subjects(isubject);
                if CONN_x.Setup.structural_sessionspecific, nsess=CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)); else nsess=1; end
                for nses=1:nsess
                    if ismember(nses,sessions)
                        jsubject=jsubject+1;
                        matlabbatch{end+1}.spm.util.imcalc.expression='(i2+i3+i4).*i1';
                        matlabbatch{end}.spm.util.imcalc.input=reshape(outputfiles{isubject}{nses}(1:4),[],1);
                        matlabbatch{end}.spm.util.imcalc.output=conn_prepend('c0',CONN_x.Setup.structural{nsubject}{nses}{1});
                        matlabbatch{end}.spm.util.imcalc.options.dtype=spm_type('float32');
                        outputfiles{isubject}{nses}{1}=conn_prepend('c0',CONN_x.Setup.structural{nsubject}{nses}{1});
                        matlabbatch{end+1}.spm.util.imcalc.expression='(i2+i3+i4)';
                        matlabbatch{end}.spm.util.imcalc.input=reshape(outputfiles{isubject}{nses}(1:4),[],1);
                        matlabbatch{end}.spm.util.imcalc.output=conn_prepend('c0mask',CONN_x.Setup.structural{nsubject}{nses}{1});
                        matlabbatch{end}.spm.util.imcalc.options.dtype=spm_type('float32');
                    end
                end
            end
            
        case {'structural_normalize','structural_normalize_preservemasks','functional_normalize_indirect','functional_normalize_indirect_preservemasks'}
            DOSPM12=~PREFERSPM8OVERSPM12&spmver12; %SPM12/SPM8
            if strcmp(regexprep(lower(STEP),'^run_|^update_|^interactive_',''),'functional_normalize_indirect')||strcmp(regexprep(lower(STEP),'^run_|^update_|^interactive_',''),'functional_normalize_indirect_preservemasks')
                jsubject=0;
                for isubject=1:numel(subjects), % coregister
                    nsubject=subjects(isubject);
                    if CONN_x.Setup.structural_sessionspecific, nsess_struct=intersect(sessions,1:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)));
                    else nsess_struct=1;
                    end
                    for nses_struct=nsess_struct(:)'
                        if CONN_x.Setup.structural_sessionspecific, nsess_func=nses_struct;
                        else nsess_func=intersect(sessions,1:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)));
                        end
                        if ~isempty(nsess_func)
                            jsubject=jsubject+1;
                            if isempty(CONN_x.Setup.structural{nsubject}{nses_struct}{1}), error(sprintf('No structural file defined for Subject %d Session %d. Please select a structural file',nsubject,nses_struct));
                            elseif numel(CONN_x.Setup.structural{nsubject}{nses_struct}{3})>1, error(sprintf('Multiple structural files found for Subject %d Session %d. Please select a single structural file',nsubject,nses_struct));
                            end
                            if ismember(nses_struct,sessions), reffile=CONN_x.Setup.structural{nsubject}{nses_struct}{1};
                            else
                                if DOSPM12, [matfile,nill,reffile]=conn_setup_preproc_meanimage(CONN_x.Setup.structural{nsubject}{nses_struct}{1},'norm_spm12');
                                else [matfile,nill,reffile]=conn_setup_preproc_meanimage(CONN_x.Setup.structural{nsubject}{nses_struct}{1},'norm_spm8');
                                end
                            end
                            filename=conn_get_functional(nsubject,nsess_func(1),sets);
                            if isempty(filename), error('Functional data not yet defined for subject %d session %d',nsubject,nsess_func(1)); end
                            temp=cellstr(filename);
                            if numel(temp)==1, ttemp=cellstr(conn_expandframe(temp{1})); else ttemp=temp; end
                            if coregtomean==2
                                if ~isempty(coregsource)&&iscell(coregsource)&&numel(coregsource)>=isubject
                                    xtemp={coregsource{isubject}};
                                elseif numel(CONN_x.Setup.coregsource_functional)>=nsubject
                                    xtemp=CONN_x.Setup.coregsource_functional{nsubject}(1);
                                else error('missing coregsource info');
                                end
                            elseif coregtomean,
                                [xtemp,failed]=conn_setup_preproc_meanimage(temp{1});
                                if isempty(xtemp),  errmsg=['Error preparing files for coregistration. Mean functional file not found (associated with functional data ',failed,' ; run either realignment or ART first, or select the option "first functional volume as reference" if you prefer to use the first functional volume instead of the mean functional volume as reference)']; conn_disp(errmsg); error(errmsg); end
                                xtemp={xtemp};
                            else xtemp=ttemp(1);
                            end
                            matlabbatch{end+1}.spm.spatial.coreg.estimate.source=xtemp;
                            matlabbatch{end}.spm.spatial.coreg.estimate.ref={reffile};
                            if 0,%coregtomean, matlabbatch{end}.spm.spatial.coreg.estimate.other=xtemp;
                            else matlabbatch{end}.spm.spatial.coreg.estimate.other={};
                            end
                            for nsestrue=nsess_func(:)'
                                filename=conn_get_functional(nsubject,nsestrue,sets);
                                if isempty(filename), error('Functional data not yet defined for subject %d session %d',nsubject,nsestrue); end
                                temp=cellstr(filename);
                                if numel(temp)==1, ttemp=cellstr(conn_expandframe(temp{1})); else ttemp=temp; end
                                matlabbatch{end}.spm.spatial.coreg.estimate.other=cat(1,matlabbatch{end}.spm.spatial.coreg.estimate.other,ttemp);
                            end
                            %                             if strcmp(regexprep(lower(STEP),'^run_|^update_|^interactive_',''),'functional_normalize_indirect_preservemasks')
                            %                                 fmask1=CONN_x.Setup.rois.files{nsubject}{1}{nses_struct}{1};
                            %                                 if isempty(fmask1), error('Grey Matter ROI data not yet defined for subject %d session %d',nsubject,nses_struct); end
                            %                                 fmask2=CONN_x.Setup.rois.files{nsubject}{2}{nses_struct}{1};
                            %                                 if isempty(fmask2), error('White Matter ROI data not yet defined for subject %d session %d',nsubject,nses_struct); end
                            %                                 fmask3=CONN_x.Setup.rois.files{nsubject}{3}{nses_struct}{1};
                            %                                 if isempty(fmask3), error('CSF ROI data not yet defined for subject %d session %d',nsubject,nses_struct); end
                            %                                 temp=[cellstr(fmask1);cellstr(fmask2);cellstr(fmask3)];
                            %                                 outputfiles{isubject}{nses_struct}{3}=temp;
                            %                                 matlabbatch{end}.spm.spatial.coreg.estimate.other=cat(1,matlabbatch{end}.spm.spatial.coreg.estimate.other,temp(:));
                            %                             end
                        end
                    end
                end
                if ~jsubject, matlabbatch=matlabbatch(1:end-1); end
            end
            if DOSPM12
                %note: structural_template disregarded (using tissue probability maps instead)
                matlabbatch{end+1}.spm.spatial.normalise.estwrite.woptions.bb=boundingbox;
                matlabbatch{end}.spm.spatial.normalise.estwrite.woptions.vox=voxelsize_anat.*[1 1 1];
                if ~isempty(interp), matlabbatch{end}.spm.spatial.normalise.estwrite.woptions.interp=interp;
                else matlabbatch{end}.spm.spatial.normalise.estwrite.woptions.interp=1;
                end
                if ~isempty(tpm_template), [nill,matlabbatch{end}.spm.spatial.normalise.estwrite.eoptions.tpm]=conn_setup_preproc_tissue(tpm_template,tpm_ngaus,subjects,sessions); end
                if ~isempty(affreg), matlabbatch{end}.spm.spatial.normalise.estwrite.eoptions.affreg=affreg; end
            else
                %note: tissue probability maps disregarded (using structural template instead)
                matlabbatch{end+1}.spm.spatial.normalise.estwrite.roptions.bb=boundingbox;
                matlabbatch{end}.spm.spatial.normalise.estwrite.roptions.vox=voxelsize_anat.*[1 1 1];
                matlabbatch{end}.spm.spatial.normalise.estwrite.eoptions.template={structural_template};
                if ~isempty(interp), matlabbatch{end}.spm.spatial.normalise.estwrite.eoptions.interp=interp;
                else matlabbatch{end}.spm.spatial.normalise.estwrite.eoptions.interp=1;
                end
            end
            jsubject=0;
            for isubject=1:numel(subjects), % structural normalize
                nsubject=subjects(isubject);
                if CONN_x.Setup.structural_sessionspecific, nsess=CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)); else nsess=1; end
                for nses=1:nsess
                    if ismember(nses,sessions)
                        jsubject=jsubject+1;
                        if isempty(CONN_x.Setup.structural{nsubject}{nses}{1}), error(sprintf('No structural file defined for Subject %d Session %d. Please select a structural file',nsubject,nses));
                        elseif numel(CONN_x.Setup.structural{nsubject}{nses}{3})>1, error(sprintf('Multiple structural files found for Subject %d Session %d. Please select a single structural file',nsubject,nses));
                        end
                        if DOSPM12, matlabbatch{end}.spm.spatial.normalise.estwrite.subj(jsubject).vol={CONN_x.Setup.structural{nsubject}{nses}{1}};
                        else        matlabbatch{end}.spm.spatial.normalise.estwrite.subj(jsubject).source={CONN_x.Setup.structural{nsubject}{nses}{1}};
                        end
                        matlabbatch{end}.spm.spatial.normalise.estwrite.subj(jsubject).resample={CONN_x.Setup.structural{nsubject}{nses}{1}};
                        if DOSPM12,
                            outputfiles{isubject}{nses}{1}=CONN_x.Setup.structural{nsubject}{nses}{1};
                            outputfiles{isubject}{nses}{5}=conn_prepend('y_',CONN_x.Setup.structural{nsubject}{nses}{1},'.nii');
                        else
                            outputfiles{isubject}{nses}{1}=CONN_x.Setup.structural{nsubject}{nses}{1};
                            outputfiles{isubject}{nses}{5}=conn_prepend('',CONN_x.Setup.structural{nsubject}{nses}{1},'_seg_sn.mat');
                        end
                    else
                        if DOSPM12, [outputfiles{isubject}{nses}{5},nill,outputfiles{isubject}{nses}{1}]=conn_setup_preproc_meanimage(CONN_x.Setup.structural{nsubject}{nses}{1},'norm_spm12');
                        else [outputfiles{isubject}{nses}{5},nill,outputfiles{isubject}{nses}{1}]=conn_setup_preproc_meanimage(CONN_x.Setup.structural{nsubject}{nses}{1},'norm_spm8');
                        end
                    end
                end
            end
            if ~jsubject, matlabbatch=matlabbatch(1:end-1); end
            
            doapplyfunctional=applytofunctional||strcmp(regexprep(lower(STEP),'^run_|^update_|^interactive_',''),'functional_normalize_indirect')||strcmp(regexprep(lower(STEP),'^run_|^update_|^interactive_',''),'functional_normalize_indirect_preservemasks');
            if doapplyfunctional % functional resample
                if DOSPM12
                    matlabbatch{end+1}.spm.spatial.normalise.write.woptions.bb=boundingbox;
                    matlabbatch{end}.spm.spatial.normalise.write.woptions.vox=voxelsize_func.*[1 1 1];
                    if ~isempty(interp), matlabbatch{end}.spm.spatial.normalise.write.woptions.interp=interp; end
                else
                    matlabbatch{end+1}.spm.spatial.normalise.write.roptions.bb=boundingbox;
                    matlabbatch{end}.spm.spatial.normalise.write.roptions.vox=voxelsize_func.*[1 1 1];
                    if ~isempty(interp), matlabbatch{end}.spm.spatial.normalise.write.roptions.interp=interp; end
                end
                jsubject=0;
                for isubject=1:numel(subjects), % normalize write
                    nsubject=subjects(isubject);
                    if CONN_x.Setup.structural_sessionspecific, nsess_struct=intersect(sessions,1:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)));
                    else nsess_struct=1;
                    end
                    for nses=nsess_struct(:)'
                        if CONN_x.Setup.structural_sessionspecific, nsess_func=nses;
                        else nsess_func=intersect(sessions,1:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)));
                        end
                        if any(ismember(nsess_func,sessions))
                            jsubject=jsubject+1;
                            if DOSPM12, matlabbatch{end}.spm.spatial.normalise.write.subj(jsubject).def=outputfiles{isubject}{nses}(5);
                            else        matlabbatch{end}.spm.spatial.normalise.write.subj(jsubject).matname=outputfiles{isubject}{nses}(5);
                            end
                            matlabbatch{end}.spm.spatial.normalise.write.subj(jsubject).resample={};
                            %outputfiles{isubject}{nses}=outputfiles{isubject}{nses}(1);
                            for nsestrue=nsess_func(:)'
                                if ismember(nsestrue,sessions)
                                    filename=conn_get_functional(nsubject,nsestrue,sets);
                                    if isempty(filename), error('Functional data not yet defined for subject %d session %d',nsubject,nsestrue); end
                                    temp=cellstr(filename);
                                    if coregtomean, % keeps mean image in same space in case it is required later
                                        [xtemp,failed]=conn_setup_preproc_meanimage(temp{1});
                                        if ~isempty(xtemp),
                                            xtemp={xtemp};
                                            matlabbatch{end}.spm.spatial.normalise.write.subj(jsubject).resample=cat(1,matlabbatch{end}.spm.spatial.normalise.write.subj(jsubject).resample,xtemp);
                                        end
                                    end
                                    matlabbatch{end}.spm.spatial.normalise.write.subj(jsubject).resample=cat(1,matlabbatch{end}.spm.spatial.normalise.write.subj(jsubject).resample,temp);
                                    outputfiles{isubject}{nsestrue}{2}=char(conn_prepend('w',temp));
                                end
                            end
                        end
                    end
                end
                if ~jsubject, matlabbatch=matlabbatch(1:end-1); end
            end
            
            if DOSPM12
                matlabbatch{end+1}.spm.spatial.normalise.write.woptions.bb=boundingbox;
                matlabbatch{end}.spm.spatial.normalise.write.woptions.vox=voxelsize_anat.*[1 1 1];
                if ~isempty(interp), matlabbatch{end}.spm.spatial.normalise.write.woptions.interp=interp;
                else matlabbatch{end}.spm.spatial.normalise.write.woptions.interp=1;
                end
            else
                matlabbatch{end+1}.spm.spatial.normalise.write.roptions.bb=boundingbox;
                matlabbatch{end}.spm.spatial.normalise.write.roptions.vox=voxelsize_anat.*[1 1 1];
                if ~isempty(interp), matlabbatch{end}.spm.spatial.normalise.write.roptions.interp=interp;
                else matlabbatch{end}.spm.spatial.normalise.write.roptions.interp=1;
                end
            end
            jsubject=0;
            for isubject=1:numel(subjects), % normalize write
                nsubject=subjects(isubject);
                if CONN_x.Setup.structural_sessionspecific, nsess_struct=intersect(sessions,1:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)));
                else nsess_struct=1;
                end
                for nses=nsess_struct(:)'
                    if CONN_x.Setup.structural_sessionspecific, nsess_func=nses;
                    else nsess_func=intersect(sessions,1:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)));
                    end
                    if ismember(nses,sessions),%||(doapplyfunctional&&any(ismember(nsess_func,sessions)))
                        jsubject=jsubject+1;
                        if DOSPM12, matlabbatch{end}.spm.spatial.normalise.write.subj(jsubject).def=outputfiles{isubject}{nses}(5);
                        else        matlabbatch{end}.spm.spatial.normalise.write.subj(jsubject).matname=outputfiles{isubject}{nses}(5);
                        end
                        matlabbatch{end}.spm.spatial.normalise.write.subj(jsubject).resample={};
                        matlabbatch{end}.spm.spatial.normalise.write.subj(jsubject).resample=outputfiles{isubject}{nses}(1);
                        %outputfiles{isubject}{nses}=outputfiles{isubject}{nses}(1);
                    end
                    outputfiles{isubject}{nses}{1}=conn_prepend('w',outputfiles{isubject}{nses}{1});
                    if strcmp(regexprep(lower(STEP),'^run_|^update_|^interactive_',''),'functional_normalize_indirect_preservemasks')||strcmp(regexprep(lower(STEP),'^run_|^update_|^interactive_',''),'structural_normalize_preservemasks')
                        if ismember(nses,sessions)
                            fmask1=CONN_x.Setup.rois.files{nsubject}{1}{nses}{1};
                            if isempty(fmask1), error('Grey Matter ROI data not yet defined for subject %d session %d',nsubject,nses); end
                            fmask2=CONN_x.Setup.rois.files{nsubject}{2}{nses}{1};
                            if isempty(fmask2), error('White Matter ROI data not yet defined for subject %d session %d',nsubject,nses); end
                            fmask3=CONN_x.Setup.rois.files{nsubject}{3}{nses}{1};
                            if isempty(fmask3), error('CSF ROI data not yet defined for subject %d session %d',nsubject,nses); end
                            temp=[cellstr(fmask1);cellstr(fmask2);cellstr(fmask3)];
                            outputfiles{isubject}{nses}{3}=temp;
                            matlabbatch{end}.spm.spatial.normalise.write.subj(jsubject).resample=cat(1,matlabbatch{end}.spm.spatial.normalise.write.subj(jsubject).resample,outputfiles{isubject}{nses}{3});
                            outputfiles{isubject}{nses}{3}=conn_prepend('w',outputfiles{isubject}{nses}{3});
                        end
                    end
                end
            end
            if ~jsubject, matlabbatch=matlabbatch(1:end-1); end
            
        case {'structural_segment&normalize','functional_segment&normalize_indirect'}
            DOSPM12=~PREFERSPM8OVERSPM12&spmver12; %SPM12/SPM8
            if strcmp(regexprep(lower(STEP),'^run_|^update_|^interactive_',''),'functional_segment&normalize_indirect')
                jsubject=0;
                for isubject=1:numel(subjects), % coregister
                    nsubject=subjects(isubject);
                    if CONN_x.Setup.structural_sessionspecific, nsess_struct=intersect(sessions,1:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)));
                    else nsess_struct=1;
                    end
                    for nses_struct=nsess_struct(:)' %note: any structural targets for in-sessions functionals
                        if CONN_x.Setup.structural_sessionspecific, nsess_func=nses_struct;
                        else nsess_func=intersect(sessions,1:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)));
                        end
                        if ~isempty(nsess_func)
                            jsubject=jsubject+1;
                            if isempty(CONN_x.Setup.structural{nsubject}{nses_struct}{1}), error(sprintf('No structural file defined for Subject %d Session %d. Please select a structural file',nsubject,nses_struct));
                            elseif numel(CONN_x.Setup.structural{nsubject}{nses_struct}{3})>1, error(sprintf('Multiple structural files found for Subject %d Session %d. Please select a single structural file',nsubject,nses_struct));
                            end
                            if ismember(nses_struct,sessions), reffile=CONN_x.Setup.structural{nsubject}{nses_struct}{1};
                            else
                                if DOSPM12, [matfile,nill,reffile]=conn_setup_preproc_meanimage(CONN_x.Setup.structural{nsubject}{nses_struct}{1},'norm_spm12');
                                else [matfile,nill,reffile]=conn_setup_preproc_meanimage(CONN_x.Setup.structural{nsubject}{nses_struct}{1},'norm_spm8');
                                end
                            end
                            filename=conn_get_functional(nsubject,nsess_func(1),sets);
                            if isempty(filename), error('Functional data not yet defined for subject %d session %d',nsubject,nsess_func(1)); end
                            temp=cellstr(filename);
                            if numel(temp)==1, ttemp=cellstr(conn_expandframe(temp{1})); else ttemp=temp; end
                            if coregtomean==2
                                if ~isempty(coregsource)&&iscell(coregsource)&&numel(coregsource)>=isubject
                                    xtemp={coregsource{isubject}};
                                elseif numel(CONN_x.Setup.coregsource_functional)>=nsubject
                                    xtemp=CONN_x.Setup.coregsource_functional{nsubject}(1);
                                else error('missing coregsource info');
                                end
                            elseif coregtomean,
                                [xtemp,failed]=conn_setup_preproc_meanimage(temp{1});
                                if isempty(xtemp),  errmsg=['Error preparing files for coregistration. Mean functional file not found (associated with functional data ',failed,' ; run either realignment or ART first, or select the option "first functional volume as reference" -#coregtomean field- if you prefer to use the first functional volume instead of the mean functional volume as reference)']; conn_disp(errmsg); error(errmsg); end
                                xtemp={xtemp};
                            else xtemp=ttemp(1);
                            end
                            matlabbatch{end+1}.spm.spatial.coreg.estimate.source=xtemp;
                            matlabbatch{end}.spm.spatial.coreg.estimate.ref={reffile};
                            if 0,%coregtomean, matlabbatch{end}.spm.spatial.coreg.estimate.other=xtemp;
                            else matlabbatch{end}.spm.spatial.coreg.estimate.other={};
                            end
                            for nsestrue=nsess_func(:)'
                                filename=conn_get_functional(nsubject,nsestrue,sets);
                                if isempty(filename), error('Functional data not yet defined for subject %d session %d',nsubject,nsestrue); end
                                temp=cellstr(filename);
                                if numel(temp)==1, ttemp=cellstr(conn_expandframe(temp{1})); else ttemp=temp; end
                                matlabbatch{end}.spm.spatial.coreg.estimate.other=cat(1,matlabbatch{end}.spm.spatial.coreg.estimate.other,ttemp);
                            end
                        end
                    end
                end
                if ~jsubject, matlabbatch=matlabbatch(1:end-1); end
            end
            if DOSPM12, matlabbatch{end+1}.spm.spatial.preproc.channel.vols={};
            else  matlabbatch{end+1}.spm.spatial.preproc.data={};
            end
            jsubject=0;
            for isubject=1:numel(subjects), % structural segment&normalize
                nsubject=subjects(isubject);
                if CONN_x.Setup.structural_sessionspecific, nsess=CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)); else nsess=1; end
                for nses=1:nsess % all structurals (in-sessions or not)
                    if ismember(nses,sessions)
                        jsubject=jsubject+1;
                        if isempty(CONN_x.Setup.structural{nsubject}{nses}{1}), error(sprintf('No structural file defined for Subject %d Session %d. Please select a structural file',nsubject,nses));
                        elseif numel(CONN_x.Setup.structural{nsubject}{nses}{3})>1, error(sprintf('Multiple structural files found for Subject %d Session %d. Please select a single structural file',nsubject,nses));
                        end
                        if DOSPM12,
                            matlabbatch{end}.spm.spatial.preproc.channel.vols{jsubject}=CONN_x.Setup.structural{nsubject}{nses}{1};
                            outputfiles{isubject}{nses}{1}=CONN_x.Setup.structural{nsubject}{nses}{1};
                            outputfiles{isubject}{nses}{5}=conn_prepend('y_',CONN_x.Setup.structural{nsubject}{nses}{1},'.nii');  % note: fix SPM12 issue converting .img to .nii
                        else
                            matlabbatch{end}.spm.spatial.preproc.data{jsubject}=CONN_x.Setup.structural{nsubject}{nses}{1};
                            outputfiles{isubject}{nses}{1}=CONN_x.Setup.structural{nsubject}{nses}{1};
                            outputfiles{isubject}{nses}{5}=conn_prepend('',CONN_x.Setup.structural{nsubject}{nses}{1},'_seg_sn.mat');
                        end
                    else
                        if DOSPM12, [outputfiles{isubject}{nses}{5},nill,outputfiles{isubject}{nses}{1}]=conn_setup_preproc_meanimage(CONN_x.Setup.structural{nsubject}{nses}{1},'norm_spm12');
                        else [outputfiles{isubject}{nses}{5},nill,outputfiles{isubject}{nses}{1}]=conn_setup_preproc_meanimage(CONN_x.Setup.structural{nsubject}{nses}{1},'norm_spm8');
                        end
                    end
                    if DOSPM12,
                        outputfiles{isubject}{nses}{2}=conn_prepend('c1',outputfiles{isubject}{nses}{1},'.nii'); % note: fix SPM12 issue converting .img to .nii
                        outputfiles{isubject}{nses}{3}=conn_prepend('c2',outputfiles{isubject}{nses}{1},'.nii');
                        outputfiles{isubject}{nses}{4}=conn_prepend('c3',outputfiles{isubject}{nses}{1},'.nii');
                    else
                        outputfiles{isubject}{nses}{2}=conn_prepend('c1',outputfiles{isubject}{nses}{1});
                        outputfiles{isubject}{nses}{3}=conn_prepend('c2',outputfiles{isubject}{nses}{1});
                        outputfiles{isubject}{nses}{4}=conn_prepend('c3',outputfiles{isubject}{nses}{1});
                    end
                end
            end
            if DOSPM12
                if ~isempty(tpm_template), matlabbatch{end}.spm.spatial.preproc.tissue=conn_setup_preproc_tissue(tpm_template,tpm_ngaus,subjects,sessions); end
                if ~isempty(affreg), matlabbatch{end}.spm.spatial.preproc.warp.affreg=affreg; end
                matlabbatch{end}.spm.spatial.preproc.warp.write=[1 1];
                matlabbatch{end}.spm.spatial.preproc.channel.vols=reshape(matlabbatch{end}.spm.spatial.preproc.channel.vols,[],1);
                if ~jsubject, matlabbatch=matlabbatch(1:end-1); end
            else
                if ~isempty(tpm_template),
                    if ~isempty(subjects)&&~isempty(sessions)&&(isnumeric(tpm_template)||(ischar(tpm_template)&&size(tpm_template,1)==1&&~isempty(conn_datasetlabel(tpm_template))))
                        error('unsupported subject-specific TPM in SPM8; please upgrade to SPM12 instead')
                    else temp=cellstr(conn_expandframe(tpm_template));
                    end
                    if isempty(tpm_ngaus), tpm_ngaus=[2 2 2 4]; end % grey/white/CSF (+other implicit)
                    matlabbatch{end}.spm.spatial.preproc.opts.tpm=temp;
                    matlabbatch{end}.spm.spatial.preproc.opts.ngaus=ngaus(1:numel(temp)+1);
                end
                matlabbatch{end}.spm.spatial.preproc.roptions.bb=boundingbox;
                matlabbatch{end}.spm.spatial.preproc.output.GM=[0,0,1];
                matlabbatch{end}.spm.spatial.preproc.output.WM=[0,0,1];
                matlabbatch{end}.spm.spatial.preproc.output.CSF=[0,0,1];
                if ~jsubject, matlabbatch=matlabbatch(1:end-1); end
            end
            
            doapplyfunctional=applytofunctional||strcmp(regexprep(lower(STEP),'^run_|^update_|^interactive_',''),'functional_segment&normalize_indirect');
            if doapplyfunctional % resample functional
                if DOSPM12
                    matlabbatch{end+1}.spm.spatial.normalise.write.woptions.bb=boundingbox;
                    matlabbatch{end}.spm.spatial.normalise.write.woptions.vox=voxelsize_func.*[1 1 1];
                    if ~isempty(interp), matlabbatch{end}.spm.spatial.normalise.write.woptions.interp=interp; end
                else
                    matlabbatch{end+1}.spm.spatial.normalise.write.roptions.bb=boundingbox;
                    matlabbatch{end}.spm.spatial.normalise.write.roptions.vox=voxelsize_func.*[1 1 1];
                    if ~isempty(interp), matlabbatch{end}.spm.spatial.normalise.write.roptions.interp=interp; end
                end
                jsubject=0;
                for isubject=1:numel(subjects), % normalize write functional
                    nsubject=subjects(isubject);
                    if CONN_x.Setup.structural_sessionspecific, nsess_struct=intersect(sessions,1:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)));
                    else nsess_struct=1;
                    end
                    for nses=nsess_struct(:)' %note: any structural targets for in-sessions functionals
                        if CONN_x.Setup.structural_sessionspecific, nsess_func=nses;
                        else nsess_func=intersect(sessions,1:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)));
                        end
                        if any(ismember(nsess_func,sessions))
                            jsubject=jsubject+1;
                            if DOSPM12, matlabbatch{end}.spm.spatial.normalise.write.subj(jsubject).def=outputfiles{isubject}{nses}(5);
                            else        matlabbatch{end}.spm.spatial.normalise.write.subj(jsubject).matname=outputfiles{isubject}{nses}(5);
                            end
                            matlabbatch{end}.spm.spatial.normalise.write.subj(jsubject).resample={};
                            for nsestrue=nsess_func(:)'
                                if ismember(nsestrue,sessions)
                                    filename=conn_get_functional(nsubject,nsestrue,sets);
                                    if isempty(filename), error('Functional data not yet defined for subject %d session %d',nsubject,nsestrue); end
                                    temp=cellstr(filename);
                                    if coregtomean, % keeps mean image in same space in case it is required later
                                        [xtemp,failed]=conn_setup_preproc_meanimage(temp{1});
                                        if ~isempty(xtemp),
                                            xtemp={xtemp};
                                            matlabbatch{end}.spm.spatial.normalise.write.subj(jsubject).resample=cat(1,matlabbatch{end}.spm.spatial.normalise.write.subj(jsubject).resample,xtemp);
                                        end
                                    end
                                    matlabbatch{end}.spm.spatial.normalise.write.subj(jsubject).resample=cat(1,matlabbatch{end}.spm.spatial.normalise.write.subj(jsubject).resample,temp);
                                    outputfiles{isubject}{nsestrue}{6}=char(conn_prepend('w',temp));
                                end
                            end
                        end
                    end
                end
                if ~jsubject, matlabbatch=matlabbatch(1:end-1); end
            end
            
            if DOSPM12
                matlabbatch{end+1}.spm.spatial.normalise.write.woptions.bb=boundingbox;
                matlabbatch{end}.spm.spatial.normalise.write.woptions.vox=voxelsize_anat.*[1 1 1];
                if ~isempty(interp), matlabbatch{end}.spm.spatial.normalise.write.woptions.interp=interp;
                else matlabbatch{end}.spm.spatial.normalise.write.woptions.interp=1;
                end
            else
                matlabbatch{end+1}.spm.spatial.normalise.write.roptions.bb=boundingbox;
                matlabbatch{end}.spm.spatial.normalise.write.roptions.vox=voxelsize_anat.*[1 1 1];
                if ~isempty(interp), matlabbatch{end}.spm.spatial.normalise.write.roptions.interp=interp;
                else matlabbatch{end}.spm.spatial.normalise.write.roptions.interp=1;
                end
            end
            jsubject=0;
            for isubject=1:numel(subjects), % normalize write structural
                nsubject=subjects(isubject);
                if CONN_x.Setup.structural_sessionspecific, nsess_struct=intersect(sessions,1:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)));
                else nsess_struct=1;
                end
                for nses=nsess_struct(:)' %note: any structural targets for in-sessions functionals
                    if CONN_x.Setup.structural_sessionspecific, nsess_func=nses;
                    else nsess_func=intersect(sessions,1:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)));
                    end
                    if ismember(nses,sessions),%||(doapplyfunctional&&any(ismember(nsess_func,sessions)))
                        jsubject=jsubject+1;
                        if DOSPM12, matlabbatch{end}.spm.spatial.normalise.write.subj(jsubject).def=outputfiles{isubject}{nses}(5);
                        else        matlabbatch{end}.spm.spatial.normalise.write.subj(jsubject).matname=outputfiles{isubject}{nses}(5);
                        end
                        matlabbatch{end}.spm.spatial.normalise.write.subj(jsubject).resample={};
                        if ismember(nses,sessions)
                            matlabbatch{end}.spm.spatial.normalise.write.subj(jsubject).resample=outputfiles{isubject}{nses}(1:4)';
                        end
                        %outputfiles{isubject}{nses}=outputfiles{isubject}{nses}(1:4);
                    end
                end
                if CONN_x.Setup.structural_sessionspecific, nsess=CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)); else nsess=1; end
                for nses=1:nsess
                    outputfiles{isubject}{nses}{1}=conn_prepend('w',outputfiles{isubject}{nses}{1});
                    outputfiles{isubject}{nses}{2}=conn_prepend('w',outputfiles{isubject}{nses}{2});
                    outputfiles{isubject}{nses}{3}=conn_prepend('w',outputfiles{isubject}{nses}{3});
                    outputfiles{isubject}{nses}{4}=conn_prepend('w',outputfiles{isubject}{nses}{4});
                end
            end
            if ~jsubject, matlabbatch=matlabbatch(1:end-1); end
            
            jsubject=0;
            for isubject=1:numel(subjects), % imcalc
                nsubject=subjects(isubject);
                if CONN_x.Setup.structural_sessionspecific, nsess=CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)); else nsess=1; end
                for nses=1:nsess
                    if ismember(nses,sessions)
                        jsubject=jsubject+1;
                        matlabbatch{end+1}.spm.util.imcalc.expression='(i2+i3+i4).*i1';
                        matlabbatch{end}.spm.util.imcalc.input=reshape(outputfiles{isubject}{nses}(1:4),[],1);
                        matlabbatch{end}.spm.util.imcalc.output=conn_prepend('wc0',CONN_x.Setup.structural{nsubject}{nses}{1});
                        matlabbatch{end}.spm.util.imcalc.options.dtype=spm_type('float32');
                        matlabbatch{end+1}.spm.util.imcalc.expression='(i2+i3+i4)';
                        matlabbatch{end}.spm.util.imcalc.input=reshape(outputfiles{isubject}{nses}(1:4),[],1);
                        matlabbatch{end}.spm.util.imcalc.output=conn_prepend('wc0mask',CONN_x.Setup.structural{nsubject}{nses}{1});
                        matlabbatch{end}.spm.util.imcalc.options.dtype=spm_type('float32');
                    end
                    outputfiles{isubject}{nses}{1}=conn_prepend('wc0',conn_prepend(-1,outputfiles{isubject}{nses}{1}));
                end
            end
            
        case {'functional_coregister_nonlinear'}
            % functional segment&normalize
            DOSPM12=~PREFERSPM8OVERSPM12&spmver12; %SPM12/SPM8
            if DOSPM12, matlabbatch{end+1}.spm.spatial.preproc.channel.vols={};
            else matlabbatch{end+1}.spm.spatial.preproc.data={};
            end
            for isubject=1:numel(subjects),
                nsubject=subjects(isubject);
                filename=conn_get_functional(nsubject,sessions(1),sets);
                if isempty(filename), error('Functional data not yet defined for subject %d session %d',nsubject,sessions(1)); end
                temp=cellstr(filename);
                if numel(temp)==1, ttemp=cellstr(conn_expandframe(temp{1})); else ttemp=temp; end
                if coregtomean==2
                    if ~isempty(coregsource)&&iscell(coregsource)&&numel(coregsource)>=isubject
                        xtemp={coregsource{isubject}};
                    elseif numel(CONN_x.Setup.coregsource_functional)>=nsubject
                        xtemp=CONN_x.Setup.coregsource_functional{nsubject}(1);
                    else error('missing coregsource info');
                    end
                elseif coregtomean,
                    [xtemp,failed]=conn_setup_preproc_meanimage(temp{1});
                    if isempty(xtemp),  errmsg=['Error preparing files for normalization. Mean functional file not found (associated with functional data ',failed,' ; run either realignment or ART first, or select the option "first functional volume as reference" -#coregtomean field- if you prefer to use the first functional volume instead of the mean functional volume as reference)']; conn_disp(errmsg); error(errmsg); end
                    xtemp={xtemp};
                else xtemp=ttemp(1);
                end
                if DOSPM12,
                    matlabbatch{end}.spm.spatial.preproc.channel.vols{isubject}=xtemp{1};
                    outputfiles{isubject}{1}=xtemp{1};
                    outputfiles{isubject}{2}=conn_prepend('y_',xtemp{1},'.nii');
                else
                    matlabbatch{end}.spm.spatial.preproc.data{isubject}=xtemp{1};
                    outputfiles{isubject}{1}=xtemp{1};
                    outputfiles{isubject}{2}=conn_prepend('',xtemp{1},'_seg_sn.mat');
                end
            end
            if DOSPM12,
                if ~isempty(tpm_template), matlabbatch{end}.spm.spatial.preproc.tissue=conn_setup_preproc_tissue(tpm_template,tpm_ngaus,subjects,sessions); end
                if ~isempty(affreg), matlabbatch{end}.spm.spatial.preproc.warp.affreg=affreg; end
                matlabbatch{end}.spm.spatial.preproc.channel.vols=reshape(matlabbatch{end}.spm.spatial.preproc.channel.vols,[],1);
                matlabbatch{end}.spm.spatial.preproc.warp.write=[1 1];
                matlabbatch{end+1}.spm.spatial.normalise.write.woptions.bb=boundingbox;
                matlabbatch{end}.spm.spatial.normalise.write.woptions.vox=voxelsize_func.*[1 1 1];
                if ~isempty(interp), matlabbatch{end}.spm.spatial.normalise.write.woptions.interp=interp; end
            else
                if ~isempty(tpm_template),
                    if ~isempty(subjects)&&~isempty(sessions)&&(isnumeric(tpm_template)||(ischar(tpm_template)&&size(tpm_template,1)==1&&~isempty(conn_datasetlabel(tpm_template))))
                        error('unsupported subject-specific TPM in SPM8; please upgrade to SPM12 instead')
                    else temp=cellstr(conn_expandframe(tpm_template));
                    end
                    if isempty(tpm_ngaus), tpm_ngaus=[2 2 2 4]; end % grey/white/CSF (+other implicit)
                    matlabbatch{end}.spm.spatial.preproc.opts.tpm=temp;
                    matlabbatch{end}.spm.spatial.preproc.opts.ngaus=ngaus(1:numel(temp)+1);
                end
                matlabbatch{end}.spm.spatial.preproc.roptions.bb=boundingbox;
                matlabbatch{end}.spm.spatial.preproc.output.GM=[0,0,1];
                matlabbatch{end}.spm.spatial.preproc.output.WM=[0,0,1];
                matlabbatch{end}.spm.spatial.preproc.output.CSF=[0,0,1];
                matlabbatch{end+1}.spm.spatial.normalise.write.roptions.bb=boundingbox;
                matlabbatch{end}.spm.spatial.normalise.write.roptions.vox=voxelsize_func.*[1 1 1];
                if ~isempty(interp), matlabbatch{end}.spm.spatial.normalise.write.roptions.interp=interp; end
            end
            for isubject=1:numel(subjects),
                nsubject=subjects(isubject);
                if DOSPM12, matlabbatch{end}.spm.spatial.normalise.write.subj(isubject).def=outputfiles{isubject}(2);
                else        matlabbatch{end}.spm.spatial.normalise.write.subj(isubject).matname=outputfiles{isubject}(2);
                end
                matlabbatch{end}.spm.spatial.normalise.write.subj(isubject).resample=outputfiles{isubject}(1);
                
                nsess=CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject));
                outputfiles{isubject}=repmat({},1,nsess);
                for nses=1:nsess
                    if ismember(nses,sessions)
                        filename=conn_get_functional(nsubject,nses,sets);
                        if isempty(filename), error('Functional data not yet defined for subject %d session %d',nsubject,nses); end
                        temp=cellstr(filename);
                        if numel(temp)==1, ttemp=cellstr(conn_expandframe(temp{1})); else ttemp=temp; end
                        matlabbatch{end}.spm.spatial.normalise.write.subj(isubject).resample=cat(1,matlabbatch{end}.spm.spatial.normalise.write.subj(isubject).resample,ttemp);
                        outputfiles{isubject}{nses}{7}=char(conn_prepend('w',temp));
                    end
                end
            end
            
            % structural segment/normalize
            if DOSPM12, matlabbatch{end+1}.spm.spatial.preproc.channel.vols={};
            else  matlabbatch{end+1}.spm.spatial.preproc.data={};
            end
            jsubject=0;
            for isubject=1:numel(subjects),
                nsubject=subjects(isubject);
                if CONN_x.Setup.structural_sessionspecific, nsess=CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)); else nsess=1; end
                for nses=1:nsess
                    if ismember(nses,sessions)
                        jsubject=jsubject+1;
                        if isempty(CONN_x.Setup.structural{nsubject}{nses}{1}), error(sprintf('No structural file defined for Subject %d Session %d. Please select a structural file',nsubject,nses));
                        elseif numel(CONN_x.Setup.structural{nsubject}{nses}{3})>1, error(sprintf('Multiple structural files found for Subject %d Session %d. Please select a single structural file',nsubject,nses));
                        end
                        if DOSPM12, matlabbatch{end}.spm.spatial.preproc.channel.vols{jsubject}=CONN_x.Setup.structural{nsubject}{nses}{1};
                        else matlabbatch{end}.spm.spatial.preproc.data{jsubject}=CONN_x.Setup.structural{nsubject}{nses}{1};
                        end
                    end
                    if DOSPM12
                        outputfiles{isubject}{nses}{1}=CONN_x.Setup.structural{nsubject}{nses}{1};
                        outputfiles{isubject}{nses}{2}=conn_prepend('c1',CONN_x.Setup.structural{nsubject}{nses}{1},'.nii'); % note: fix SPM12 issue converting .img to .nii
                        outputfiles{isubject}{nses}{3}=conn_prepend('c2',CONN_x.Setup.structural{nsubject}{nses}{1},'.nii');
                        outputfiles{isubject}{nses}{4}=conn_prepend('c3',CONN_x.Setup.structural{nsubject}{nses}{1},'.nii');
                        outputfiles{isubject}{nses}{5}=conn_prepend('y_',CONN_x.Setup.structural{nsubject}{nses}{1},'.nii');
                        outputfiles{isubject}{nses}{6}=conn_prepend('iy_',CONN_x.Setup.structural{nsubject}{nses}{1},'.nii');
                    else
                        outputfiles{isubject}{nses}{1}=CONN_x.Setup.structural{nsubject}{nses}{1};
                        outputfiles{isubject}{nses}{2}=conn_prepend('c1',CONN_x.Setup.structural{nsubject}{nses}{1});
                        outputfiles{isubject}{nses}{3}=conn_prepend('c2',CONN_x.Setup.structural{nsubject}{nses}{1});
                        outputfiles{isubject}{nses}{4}=conn_prepend('c3',CONN_x.Setup.structural{nsubject}{nses}{1});
                        outputfiles{isubject}{nses}{5}=conn_prepend('',CONN_x.Setup.structural{nsubject}{nses}{1},'_seg_sn.mat');
                        outputfiles{isubject}{nses}{6}=conn_prepend('',CONN_x.Setup.structural{nsubject}{nses}{1},'_seg_inv_sn.mat');
                    end
                end
            end
            if DOSPM12
                if ~isempty(tpm_template), matlabbatch{end}.spm.spatial.preproc.tissue=conn_setup_preproc_tissue(tpm_template,tpm_ngaus,subjects,sessions); end
                if ~isempty(affreg), matlabbatch{end}.spm.spatial.preproc.warp.affreg=affreg; end
                matlabbatch{end}.spm.spatial.preproc.warp.write=[1 1];
                matlabbatch{end}.spm.spatial.preproc.channel.vols=reshape(matlabbatch{end}.spm.spatial.preproc.channel.vols,[],1);
                if ~jsubject, matlabbatch=matlabbatch(1:end-1); end
                matlabbatch{end+1}.spm.spatial.normalise.write.woptions.bb=boundingbox;
                matlabbatch{end}.spm.spatial.normalise.write.woptions.vox=voxelsize_anat.*[1 1 1];
                if ~isempty(interp), matlabbatch{end}.spm.spatial.normalise.write.woptions.interp=interp;
                else matlabbatch{end}.spm.spatial.normalise.write.woptions.interp=1;
                end
            else
                if ~isempty(tpm_template),
                    if ~isempty(subjects)&&~isempty(sessions)&&(isnumeric(tpm_template)||(ischar(tpm_template)&&size(tpm_template,1)==1&&~isempty(conn_datasetlabel(tpm_template))))
                        error('unsupported subject-specific TPM in SPM8; please upgrade to SPM12 instead')
                    else temp=cellstr(conn_expandframe(tpm_template));
                    end
                    if isempty(tpm_ngaus), tpm_ngaus=[2 2 2 4]; end % grey/white/CSF (+other implicit)
                    matlabbatch{end}.spm.spatial.preproc.opts.tpm=temp;
                    matlabbatch{end}.spm.spatial.preproc.opts.ngaus=ngaus(1:numel(temp)+1);
                end
                matlabbatch{end}.spm.spatial.preproc.roptions.bb=boundingbox; %
                matlabbatch{end}.spm.spatial.preproc.output.GM=[0,0,1];
                matlabbatch{end}.spm.spatial.preproc.output.WM=[0,0,1];
                matlabbatch{end}.spm.spatial.preproc.output.CSF=[0,0,1];
                if ~jsubject, matlabbatch=matlabbatch(1:end-1); end
                matlabbatch{end+1}.spm.spatial.normalise.write.roptions.bb=boundingbox;
                matlabbatch{end}.spm.spatial.normalise.write.roptions.vox=voxelsize_anat.*[1 1 1];
                if ~isempty(interp), matlabbatch{end}.spm.spatial.normalise.write.roptions.interp=interp;
                else matlabbatch{end}.spm.spatial.normalise.write.roptions.interp=1;
                end
            end
            jsubject=0;
            for isubject=1:numel(subjects),
                nsubject=subjects(isubject);
                if CONN_x.Setup.structural_sessionspecific, nsess=CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)); else nsess=1; end
                for nses=1:nsess
                    if CONN_x.Setup.structural_sessionspecific, nsesstrue=nses;
                    else nsesstrue=1:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject));
                    end
                    if ismember(nses,sessions)||any(ismember(nsesstrue,sessions))
                        if ismember(nses,sessions),
                            jsubject=jsubject+1; % write struct normed (optional)
                            if DOSPM12, matlabbatch{end}.spm.spatial.normalise.write.subj(jsubject).def=outputfiles{isubject}{nses}(5);
                            else        matlabbatch{end}.spm.spatial.normalise.write.subj(jsubject).matname=outputfiles{isubject}{nses}(5);
                            end
                            matlabbatch{end}.spm.spatial.normalise.write.subj(jsubject).resample=outputfiles{isubject}{nses}(1:4)';
                        end
                        if any(ismember(nsesstrue,sessions))
                            jsubject=jsubject+1; % write func inv normed
                            if DOSPM12, matlabbatch{end}.spm.spatial.normalise.write.subj(jsubject).def=outputfiles{isubject}{nses}(6);
                            else        matlabbatch{end}.spm.spatial.normalise.write.subj(jsubject).matname=outputfiles{isubject}{nses}(6);
                            end
                            matlabbatch{end}.spm.spatial.normalise.write.subj(jsubject).resample={};%outputfiles{isubject}{nses}(1:4)';
                            %outputfiles{isubject}{nses}=outputfiles{isubject}{nses}(1:4);
                            for nsestrue=nsesstrue
                                if ismember(nsestrue,sessions)
                                    temp=cellstr(outputfiles{isubject}{nsestrue}{7});
                                    %if isempty(CONN_x.Setup.functional{nsubject}{nsestrue}{1}), error('Functional data not yet defined for subject %d session %d',nsubject,nsestrue); end
                                    %temp=cellstr(CONN_x.Setup.functional{nsubject}{nsestrue}{1});
                                    if coregtomean, % keeps mean image in same space in case it is required later
                                        [xtemp,failed]=conn_setup_preproc_meanimage(temp{1});
                                        if ~isempty(xtemp),
                                            xtemp={xtemp};
                                            matlabbatch{end}.spm.spatial.normalise.write.subj(jsubject).resample=cat(1,matlabbatch{end}.spm.spatial.normalise.write.subj(jsubject).resample,xtemp);
                                        end
                                    end
                                    matlabbatch{end}.spm.spatial.normalise.write.subj(jsubject).resample=cat(1,matlabbatch{end}.spm.spatial.normalise.write.subj(jsubject).resample,temp);
                                    outputfiles{isubject}{nsestrue}{7}=char(conn_prepend('w',temp));
                                end
                            end
                        end
                    end
                    %outputfiles{isubject}{nses}{1}=conn_prepend('w',outputfiles{isubject}{nses}{1});
                    %outputfiles{isubject}{nses}{2}=conn_prepend('w',outputfiles{isubject}{nses}{2});
                    %outputfiles{isubject}{nses}{3}=conn_prepend('w',outputfiles{isubject}{nses}{3});
                    %outputfiles{isubject}{nses}{4}=conn_prepend('w',outputfiles{isubject}{nses}{4});
                end
            end
            if ~jsubject, matlabbatch=matlabbatch(1:end-1); end
            jsubject=0;
            for isubject=1:numel(subjects),
                nsubject=subjects(isubject);
                if CONN_x.Setup.structural_sessionspecific, nsess=CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)); else nsess=1; end
                for nses=1:nsess
                    if ismember(nses,sessions)
                        jsubject=jsubject+1;
                        matlabbatch{end+1}.spm.util.imcalc.expression='(i2+i3+i4).*i1';
                        matlabbatch{end}.spm.util.imcalc.input=reshape(conn_prepend('w',outputfiles{isubject}{nses}(1:4)),[],1);
                        matlabbatch{end}.spm.util.imcalc.output=conn_prepend('wc0',CONN_x.Setup.structural{nsubject}{nses}{1});
                        matlabbatch{end}.spm.util.imcalc.options.dtype=spm_type('float32');
                        matlabbatch{end+1}.spm.util.imcalc.expression='(i2+i3+i4)';
                        matlabbatch{end}.spm.util.imcalc.input=reshape(conn_prepend('w',outputfiles{isubject}{nses}(1:4)),[],1);
                        matlabbatch{end}.spm.util.imcalc.output=conn_prepend('wc0mask',CONN_x.Setup.structural{nsubject}{nses}{1});
                        matlabbatch{end}.spm.util.imcalc.options.dtype=spm_type('float32');
                        jsubject=jsubject+1;
                        matlabbatch{end+1}.spm.util.imcalc.expression='(i2+i3+i4).*i1';
                        matlabbatch{end}.spm.util.imcalc.input=reshape(outputfiles{isubject}{nses}(1:4),[],1);
                        matlabbatch{end}.spm.util.imcalc.output=conn_prepend('c0',CONN_x.Setup.structural{nsubject}{nses}{1});
                        matlabbatch{end}.spm.util.imcalc.options.dtype=spm_type('float32');
                        matlabbatch{end+1}.spm.util.imcalc.expression='(i2+i3+i4)';
                        matlabbatch{end}.spm.util.imcalc.input=reshape(outputfiles{isubject}{nses}(1:4),[],1);
                        matlabbatch{end}.spm.util.imcalc.output=conn_prepend('c0mask',CONN_x.Setup.structural{nsubject}{nses}{1});
                        matlabbatch{end}.spm.util.imcalc.options.dtype=spm_type('float32');
                    end
                    outputfiles{isubject}{nses}{1}=conn_prepend('c0',CONN_x.Setup.structural{nsubject}{nses}{1});
                end
            end
            
        case 'functional_slicetime'
            sliceorder_all=sliceorder;
            if ~iscell(sliceorder_all),sliceorder_all={sliceorder_all}; end
            for isubject=1:numel(subjects),
                nsubject=subjects(isubject);
                if isubject<=numel(sliceorder_all), sliceorder=sliceorder_all{min(numel(sliceorder_all),isubject)}; end
                nsess=CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject));
                for nses=1:nsess
                    if ismember(nses,sessions)
                        matlabbatch{end+1}.spm.temporal.st.scans={};
                        [filename,cfile]=conn_get_functional(nsubject,nses,sets,[],'cfile');
                        if isempty(filename), error('Functional data not yet defined for subject %d session %d',nsubject,nses); end
                        nslice=cfile{3}(1).dim(3);
                        temp=cellstr(filename);
                        if numel(temp)==1,
                            temp=cellstr(conn_expandframe(temp{1}));
                        end
                        matlabbatch{end}.spm.temporal.st.scans{end+1}=temp;
                        matlabbatch{end}.spm.temporal.st.tr=conn_get_rt(nsubject,nses,sets);
                        matlabbatch{end}.spm.temporal.st.nslices=nslice;
                        if isempty(ta), matlabbatch{end}.spm.temporal.st.ta=matlabbatch{end}.spm.temporal.st.tr*(1-1/nslice);
                        else matlabbatch{end}.spm.temporal.st.ta=ta(min(numel(ta),nsubject));
                        end
                        while (numel(unique(sliceorder))~=nslice||max(sliceorder)~=nslice||min(sliceorder)~=1) && (numel(sliceorder)~=nslice||any(sliceorder<0|sliceorder>conn_get_rt(nsubject,nses,sets)*1000))
                            if isempty(sliceorder_select)
                                if ~isempty(sliceorder),
                                    conn_msgbox({['Subject ',num2str(nsubject),' Session ',num2str(nses), ' Incorrectly defined slice order vector'],[num2str(nslice),' slices']},'',2);
                                    sliceorder_select=[];
                                end
                                if isempty(sliceorder_select)
                                    [sliceorder_select,tok] = listdlg('PromptString',['Select slice order (subject ',num2str(nsubject),' session ',num2str(nses),'):'],'SelectionMode','single','ListString',[sliceorder_select_options,{'manually define'}]);
                                end
                                if isempty(sliceorder_select), return; end
                            end
                            switch(sliceorder_select)
                                case 1, sliceorder=1:nslice;        % ascending
                                case 2, sliceorder=nslice:-1:1;     % descending
                                case 3, sliceorder=round((nslice-(1:nslice))/2 + (rem((nslice-(1:nslice)),2) * (nslice - 1)/2)) + 1; % interleaved (middle-top)
                                case 4, sliceorder=[1:2:nslice 2:2:nslice]; % interleaved (bottom-up)
                                case 5, sliceorder=[nslice:-2:1, nslice-1:-2:1]; % interleaved (top-down)
                                case 6, sliceorder=[fliplr(nslice:-2:1) fliplr(nslice-1:-2:1)]; % interleaved (Siemens)
                                case 7, sliceorder=cell2mat(arrayfun(@(n)n:round(sqrt(nslice)):nslice,1:round(sqrt(nslice)),'uni',0)); % interleaved (Philips)
                                case 8, % BIDS json file
                                    str=conn_jsonread(matlabbatch{end}.spm.temporal.st.scans{1}{1},'SliceTiming');
                                    if isempty(str), conn_disp('fprintf','ERROR: SliceTiming information not found in json file associated with %s (new user-input required)\n',matlabbatch{end}.spm.temporal.st.scans{1}{1});
                                    else sliceorder=1000*reshape(str,1,[]); % (ms)
                                    end
                                case 9, % manually define
                                    sliceorder=1:nslice;
                                    sliceorder=inputdlg(['Slice order? (enter slice indexes from z=1 -first slice in image- to z=',num2str(nslice),' -last slice- in the order they were acquired). Alternatively enter acquisition time of each slice in milliseconds (e.g. for multiband sequences)'],'conn_setup_preproc',1,{sprintf('%d ',sliceorder)});
                                    if isempty(sliceorder), return;
                                    else sliceorder=str2num(regexprep(sliceorder{1},'[a-zA-Z]+',num2str(nslice)));
                                    end
                            end
                            if (numel(unique(sliceorder))~=nslice||max(sliceorder)~=nslice||min(sliceorder)~=1) && (numel(sliceorder)~=nslice||any(sliceorder<0|sliceorder>conn_get_rt(nsubject,nses,sets)*1000)), sliceorder_select=[]; end
                        end
                        if (numel(unique(sliceorder))~=nslice||max(sliceorder)~=nslice||min(sliceorder)~=1),
                            matlabbatch{end}.spm.temporal.st.so=sliceorder;
                            matlabbatch{end}.spm.temporal.st.refslice=mean(sliceorder); % slice timing (ms)
                        elseif 1, % note: convert to slice-timing (ms) syntax for SPM
                            [nill,isliceorder]=sort(sliceorder);
                            matlabbatch{end}.spm.temporal.st.so=1000*(isliceorder-1)*matlabbatch{end}.spm.temporal.st.ta/(nslice-1);
                            matlabbatch{end}.spm.temporal.st.refslice=mean(matlabbatch{end}.spm.temporal.st.so); % slice timing (ms)
                        else % keep slice-order syntax (back-compatibility)
                            matlabbatch{end}.spm.temporal.st.so=sliceorder;
                            matlabbatch{end}.spm.temporal.st.refslice=sliceorder(floor(nslice/2)); % slice order
                        end
                        if isempty(matlabbatch{end}.spm.temporal.st.scans), matlabbatch=matlabbatch(1:end-1); end
                        if ~isempty(sliceorder_select)&&all(sliceorder_select<9), sliceorder=[]; end
                        outputfiles{isubject}{nses}=char(conn_prepend('a',cellstr(filename)));
                    end
                end
            end
            
        case 'functional_realign_noreslice'
            for isubject=1:numel(subjects),
                nsubject=subjects(isubject);
                matlabbatch{end+1}.spm.spatial.realign.estwrite.data={};
                nsess=CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject));
                for nses=1:nsess
                    if ismember(nses,sessions)
                        filename=conn_get_functional(nsubject,nses,sets);
                        if isempty(filename), error('Functional data not yet defined for subject %d session %d',nsubject,nses); end
                        temp=cellstr(filename);
                        temp1=temp{1};
                        matlabbatch{end}.spm.spatial.realign.estwrite.data{end+1}=temp;
                        outputfiles{isubject}{nses}{1}=char(temp);
                        outputfiles{isubject}{nses}{2}=conn_prepend('rp_',temp1,'.txt');
                    end
                end
                matlabbatch{end}.spm.spatial.realign.estwrite.eoptions.rtm=rtm;
                matlabbatch{end}.spm.spatial.realign.estwrite.roptions.which=[0,1];
                if isempty(matlabbatch{end}.spm.spatial.realign.estwrite.data), matlabbatch=matlabbatch(1:end-1); end
            end
            
        case 'functional_realign'
            for isubject=1:numel(subjects),
                nsubject=subjects(isubject);
                matlabbatch{end+1}.spm.spatial.realign.estwrite.data={};
                nsess=CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject));
                for nses=1:nsess
                    if ismember(nses,sessions)
                        filename=conn_get_functional(nsubject,nses,sets);
                        if isempty(filename), error('Functional data not yet defined for subject %d session %d',nsubject,nses); end
                        temp=cellstr(filename);
                        temp1=temp{1};
                        matlabbatch{end}.spm.spatial.realign.estwrite.data{end+1}=temp;
                        outputfiles{isubject}{nses}{1}=char(conn_prepend('r',temp));
                        outputfiles{isubject}{nses}{2}=conn_prepend('rp_',temp1,'.txt');
                    end
                end
                matlabbatch{end}.spm.spatial.realign.estwrite.eoptions.rtm=rtm;
                matlabbatch{end}.spm.spatial.realign.estwrite.roptions.which=[2,1];
                if isempty(matlabbatch{end}.spm.spatial.realign.estwrite.data), matlabbatch=matlabbatch(1:end-1); end
            end
            
        case 'functional_realign&unwarp'
            for isubject=1:numel(subjects),
                nsubject=subjects(isubject);
                matlabbatch{end+1}.spm.spatial.realignunwarp.eoptions.rtm=rtm;
                nsess=CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject));
                jses=0;
                for nses=1:nsess
                    if ismember(nses,sessions)
                        jses=jses+1;
                        filename=conn_get_functional(nsubject,nses,sets);
                        if isempty(filename), error('Functional data not yet defined for subject %d session %d',nsubject,nses); end
                        temp=cellstr(filename);
                        if numel(temp)==1, ttemp=cellstr(conn_expandframe(temp{1})); else ttemp=temp; end
                        matlabbatch{end}.spm.spatial.realignunwarp.data(jses).scans=ttemp;
                        outputfiles{isubject}{nses}{1}=char(conn_prepend('u',temp));
                        outputfiles{isubject}{nses}{2}=conn_prepend('rp_',temp{1},'.txt');
                    end
                end
                if ~jses, matlabbatch=matlabbatch(1:end-1); end
            end
            
        case {'functional_realign&unwarp&fieldmap','functional_realign&unwarp&phasemap'}
            for isubject=1:numel(subjects),
                nsubject=subjects(isubject);
                matlabbatch{end+1}.spm.spatial.realignunwarp.eoptions.rtm=rtm;
                nsess=CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject));
                jses=0;
                for nses=1:nsess
                    if ismember(nses,sessions)
                        jses=jses+1;
                        filename=conn_get_functional(nsubject,nses,sets);
                        if isempty(filename), error('Functional data not yet defined for subject %d session %d',nsubject,nses); end
                        temp=cellstr(filename);
                        if numel(temp)==1, ttemp=cellstr(conn_expandframe(temp{1})); else ttemp=temp; end
                        if ~isempty(unwarp)&&iscell(unwarp)&&numel(unwarp)>=isubject&&numel(unwarp{isubject})>=nses
                            tmfile=unwarp{isubject}{nses};
                        elseif ~isempty(conn_datasetlabel('vdm'))
                            tmfile=conn_get_functional(nsubject,nses,'vdm');
                        elseif numel(CONN_x.Setup.unwarp_functional)>=nsubject&&numel(CONN_x.Setup.unwarp_functional{nsubject})>=nses
                            tmfile=CONN_x.Setup.unwarp_functional{nsubject}{nses}{1};
                        else
                            tmfile=conn_prepend('vdm',ttemp{1});
                            if ~conn_existfile(tmfile), tmfile=conn_prepend('vdm_',ttemp{1}); end
                            if ~conn_existfile(tmfile), tmfile=conn_prepend('vdm5_',ttemp{1}); end
                            if ~conn_existfile(tmfile)
                                tmfile=dir(fullfile(fileparts(ttemp{1}),'vdm*.nii'));
                                if numel(tmfile)==1, tmfile=fullfile(fileparts(ttemp{1}),tmfile(1).name);
                                else tmfile='';
                                end
                            end
                            if ~conn_existfile(tmfile)
                                tmfile=dir(fullfile(fileparts(ttemp{1}),'vdm*.img'));
                                if numel(tmfile)==1, tmfile=fullfile(fileparts(ttemp{1}),tmfile(1).name);
                                else tmfile='';
                                end
                            end
                        end
                        if isempty(tmfile),
                            conn_disp('fprintf','unable to find voxel-displacement file %s. Switching to manual selection\n',fullfile(fileparts(ttemp{1}),'vdm*'));
                            tmfile=spm_select(1,'^vdm.*',['SUBJECT ',num2str(nsubject),'SESSION ',num2str(nses),' Voxel-Displacement Map (vdm*)'],{tmfile},fileparts(ttemp{1}));
                        end
                        if isempty(tmfile),return;end
                        if isequal(tmfile,conn_prepend('vdm_',regexprep(ttemp{1},',\d+$',''))), outtmfile=tmfile;
                        elseif isequal(tmfile,conn_prepend('vdm5_',regexprep(ttemp{1},',\d+$',''))), outtmfile=tmfile;
                        else % note: forces single vdm file per session (duplicates file if necessary)
                            outtmfile=conn_prepend('vdm_',regexprep(ttemp{1},',\d+$',''));
                            if ispc, [ok,msg]=system(sprintf('copy "%s" "%s"',tmfile,outtmfile));
                            else [ok,msg]=system(sprintf('''cp'' -f ''%s'' ''%s''',tmfile,outtmfile));
                            end
                        end
                        matlabbatch{end}.spm.spatial.realignunwarp.data(jses).scans=ttemp;
                        matlabbatch{end}.spm.spatial.realignunwarp.data(jses).pmscan={outtmfile};
                        outputfiles{isubject}{nses}{1}=char(conn_prepend('u',temp));
                        outputfiles{isubject}{nses}{2}=conn_prepend('rp_',temp{1},'.txt');
                    end
                end
                if ~jses, matlabbatch=matlabbatch(1:end-1); end
            end
            
        case {'functional_vdm_apply'}
            PED=[];
            for isubject=1:numel(subjects),
                nsubject=subjects(isubject);
                matlabbatch{end+1}.spm.tools.fieldmap.applyvdm.roptions.which=[2 0]; % create only corrected functionals
                nsess=CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject));
                jses=0;
                for nses=1:nsess
                    if ismember(nses,sessions)
                        jses=jses+1;
                        filename=conn_get_functional(nsubject,nses,sets);
                        if isempty(filename), error('Functional data not yet defined for subject %d session %d',nsubject,nses); end
                        temp=cellstr(filename);
                        if numel(temp)==1, ttemp=cellstr(conn_expandframe(temp{1})); else ttemp=temp; end
                        if isempty(PED),
                            PED=conn_jsonread(ttemp{1},'PhaseEncodingDirection',false);
                            if isempty(PED)
                                conn_disp('fprintf','warning: unable to find PhaseEncodingDirection information in %s. Assuming Posterior-Anterior (j)\n',ttemp{1});
                                PED=2;
                            else
                                if iscell(PED), PED=char(PED); end
                                if isequal(PED,'i+')||isequal(PED,'i-')||isequal(PED,'i'), PED=1;
                                elseif isequal(PED,'j+')||isequal(PED,'j-')||isequal(PED,'j'), PED=2;
                                elseif isequal(PED,'k+')||isequal(PED,'k-')||isequal(PED,'k'), PED=3;
                                else
                                    if ~isempty(PED), conn_disp('fprintf','warning: unable to interpret PhaseEncodingDirection information in %s (%s). Assuming Posterior-Anterior (j)\n',ttemp{1},PED);
                                    else conn_disp('fprintf','warning: unable to find PhaseEncodingDirection information in %s. Assuming Posterior-Anterior (j)\n',ttemp{1});
                                    end
                                    PED=2;
                                end
                            end
                            matlabbatch{end}.spm.tools.fieldmap.applyvdm.roptions.pedir=PED;
                        end
                        if ~isempty(unwarp)&&iscell(unwarp)&&numel(unwarp)>=isubject&&numel(unwarp{isubject})>=nses
                            tmfile=unwarp{isubject}{nses};
                        elseif ~isempty(conn_datasetlabel('vdm'))
                            tmfile=conn_get_functional(nsubject,nses,'vdm');
                        elseif numel(CONN_x.Setup.unwarp_functional)>=nsubject&&numel(CONN_x.Setup.unwarp_functional{nsubject})>=nses
                            tmfile=CONN_x.Setup.unwarp_functional{nsubject}{nses}{1};
                        else
                            tmfile=conn_prepend('vdm',ttemp{1});
                            if ~conn_existfile(tmfile), tmfile=conn_prepend('vdm_',ttemp{1}); end
                            if ~conn_existfile(tmfile), tmfile=conn_prepend('vdm5_',ttemp{1}); end
                            if ~conn_existfile(tmfile)
                                tmfile=dir(fullfile(fileparts(ttemp{1}),'vdm*.nii'));
                                if numel(tmfile)==1, tmfile=fullfile(fileparts(ttemp{1}),tmfile(1).name);
                                else tmfile='';
                                end
                            end
                            if ~conn_existfile(tmfile)
                                tmfile=dir(fullfile(fileparts(ttemp{1}),'vdm*.img'));
                                if numel(tmfile)==1, tmfile=fullfile(fileparts(ttemp{1}),tmfile(1).name);
                                else tmfile='';
                                end
                            end
                        end
                        if isempty(tmfile),
                            conn_disp('fprintf','unable to find voxel-displacement file %s. Switching to manual selection\n',fullfile(fileparts(ttemp{1}),'vdm*'));
                            tmfile=spm_select(1,'^vdm.*',['SUBJECT ',num2str(nsubject),'SESSION ',num2str(nses),' Voxel-Displacement Map (vdm*)'],{tmfile},fileparts(ttemp{1}));
                        end
                        if isempty(tmfile),return;end
                        if isequal(tmfile,conn_prepend('vdm_',regexprep(ttemp{1},',\d+$',''))), outtmfile=tmfile;
                        elseif isequal(tmfile,conn_prepend('vdm5_',regexprep(ttemp{1},',\d+$',''))), outtmfile=tmfile;
                        else % note: forces single vdm file per session (duplicates file if necessary)
                            outtmfile=conn_prepend('vdm_',regexprep(ttemp{1},',\d+$',''));
                            if ispc, [ok,msg]=system(sprintf('copy "%s" "%s"',tmfile,outtmfile));
                            else [ok,msg]=system(sprintf('''cp'' -f ''%s'' ''%s''',tmfile,outtmfile));
                            end
                        end
                        matlabbatch{end}.spm.tools.fieldmap.applyvdm.data(jses).scans=ttemp;
                        matlabbatch{end}.spm.tools.fieldmap.applyvdm.data(jses).vdmfile={outtmfile};
                        outputfiles{isubject}{nses}{1}=char(conn_prepend('u',temp));
                    end
                end
                if ~jses, matlabbatch=matlabbatch(1:end-1); end
            end
            
        case 'functional_art'
            icov=find(strcmp(CONN_x.Setup.l1covariates.names(1:end-1),'realignment'));
            for isubject=1:numel(subjects),
                nsubject=subjects(isubject);
                matlabbatch{end+1}.art.P={};
                matlabbatch{end}.art.M={};
                matlabbatch{end}.art.global_threshold=art_global_threshold;
                matlabbatch{end}.art.motion_threshold=art_motion_threshold;
                matlabbatch{end}.art.use_diff_motion=art_use_diff_motion;
                matlabbatch{end}.art.use_diff_global=art_use_diff_global;
                matlabbatch{end}.art.use_norms=art_use_norms;
                matlabbatch{end}.art.drop_flag=art_drop_flag;
                matlabbatch{end}.art.gui_display=art_gui_display;
                nsess=CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject));
                for nses=1:nsess
                    if ismember(nses,sessions)
                        filename=conn_get_functional(nsubject,nses,sets);
                        if isempty(filename), error('Functional data not yet defined for subject %d session %d',nsubject,nses); end
                        temp=cellstr(filename);
                        temp1=temp{1};
                        matlabbatch{end}.art.P{end+1}=char(temp);
                        if isempty(icov),
                            for remov=0:10,if conn_existfile(conn_prepend('rp_',conn_prepend(-remov,temp1),'.txt')); break; end; end
                            if remov==10, errmsg=['Error preparing files for ART processing. No ''realignment'' covariate; alternative realignment parameters file ',conn_prepend('rp_',temp1,'.txt'),' not found']; conn_disp(errmsg); error(errmsg); end
                            matlabbatch{end}.art.M{end+1}=conn_prepend('rp_',conn_prepend(-remov,temp1),'.txt');
                            matlabbatch{end}.art.motion_file_type=0;
                        else
                            matlabbatch{end}.art.motion_file_type=0; % SPM-convention
                            cfilename=CONN_x.Setup.l1covariates.files{nsubject}{icov}{nses}{1};
                            assert(~isempty(cfilename),'covariate %s has not been defined for subject %d sessions %d',CONN_x.Setup.l1covariates.names{icov},nsubject,nses);
                            switch(cfilename),
                                case '[raw values]',
                                    matlabbatch{end}.art.M{end+1}=CONN_x.Setup.l1covariates.files{nsubject}{icov}{nses}{3};
                                otherwise,
                                    matlabbatch{end}.art.M{end+1}=cfilename;
                                    [nill,fname,fext]=fileparts(matlabbatch{end}.art.M{end});
                                    if isequal(lower(fext),'.mat'), temp=load(cfilename); if isstruct(temp), tempfieldname=fieldnames(temp); temp=temp.(tempfieldname{1}); end; matlabbatch{end}.art.M{end}=temp;
                                    elseif isequal(lower(fext),'.par'), matlabbatch{end}.art.motion_file_type=1;
                                    elseif isequal(lower(fext),'.txt')&&~isempty(regexp(lower(fname),'\.siemens$')), matlabbatch{end}.art.motion_file_type=2;
                                    elseif isequal(lower(fext),'.txt')&&~isempty(regexp(lower(fname),'\.deg$')), matlabbatch{end}.art.motion_file_type=3;
                                    end
                            end
                        end
                        outputfiles{isubject}{nses}{1}=conn_prepend('art_regression_outliers_',temp1,'.mat');
                        outputfiles{isubject}{nses}{2}=conn_prepend('art_regression_timeseries_',temp1,'.mat');
                        if nses==sessions(1), matlabbatch{end}.art.output_dir=fileparts(temp1); end
                    end
                end
                if isempty(matlabbatch{end}.art.P), matlabbatch=matlabbatch(1:end-1); end
            end
            
        case {'functional_coregister_affine_reslice'}
            jsubject=0;
            for isubject=1:numel(subjects),
                nsubject=subjects(isubject);
                if CONN_x.Setup.structural_sessionspecific, nsess_struct=intersect(sessions,1:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)));
                else nsess_struct=1;
                end
                for nses_struct=nsess_struct(:)'
                    if CONN_x.Setup.structural_sessionspecific, nsess_func=nses_struct;
                    else nsess_func=intersect(sessions,1:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)));
                    end
                    if ~isempty(nsess_func)
                        jsubject=jsubject+1;
                        filename=conn_get_functional(nsubject,nsess_func(1),sets);
                        if isempty(filename), error('Functional data not yet defined for subject %d session %d',nsubject,nsess_func(1)); end
                        temp=cellstr(filename);
                        if numel(temp)==1, ttemp=cellstr(conn_expandframe(temp{1})); else ttemp=temp; end
                        if coregtomean==2
                            if ~isempty(coregsource)&&iscell(coregsource)&&numel(coregsource)>=isubject
                                xtemp={coregsource{isubject}};
                            elseif numel(CONN_x.Setup.coregsource_functional)>=nsubject
                                xtemp=CONN_x.Setup.coregsource_functional{nsubject}(1);
                            else error('missing coregsource info');
                            end
                        elseif coregtomean,
                            [xtemp,failed]=conn_setup_preproc_meanimage(temp{1});
                            if isempty(xtemp),  errmsg=['Error preparing files for coregistration. Mean functional file not found (associated with functional data ',failed,' ; run either realignment or ART first, or select the option "first functional volume as reference" -#coregtomean field- if you prefer to use the first functional volume instead of the mean functional volume as reference)']; conn_disp(errmsg); error(errmsg); end
                            xtemp={xtemp};
                        else xtemp=ttemp(1);
                        end
                        matlabbatch{end+1}.spm.spatial.coreg.estwrite.source=xtemp;
                        if isempty(CONN_x.Setup.structural{nsubject}{nses_struct}{1}), error(sprintf('No structural file defined for Subject %d Session %d. Please select a structural file',nsubject,nses_struct));
                        elseif numel(CONN_x.Setup.structural{nsubject}{nses_struct}{3})>1, error(sprintf('Multiple structural files found for Subject %d Session %d. Please select a single structural file',nsubject,nses_struct));
                        end
                        matlabbatch{end}.spm.spatial.coreg.estwrite.ref=CONN_x.Setup.structural{nsubject}{nses_struct}(1);
                        if 0,%coregtomean, matlabbatch{end}.spm.spatial.coreg.estwrite.other=xtemp;
                        else matlabbatch{end}.spm.spatial.coreg.estwrite.other={};
                        end
                        for nsestrue=nsess_func(:)'
                            filename=conn_get_functional(nsubject,nsestrue,sets);
                            if isempty(filename), error('Functional data not yet defined for subject %d session %d',nsubject,nsestrue); end
                            temp=cellstr(filename);
                            if numel(temp)==1, ttemp=cellstr(conn_expandframe(temp{1})); else ttemp=temp; end
                            matlabbatch{end}.spm.spatial.coreg.estwrite.other=cat(1,matlabbatch{end}.spm.spatial.coreg.estwrite.other,ttemp);
                        end
                    end
                end
            end
            if ~jsubject, matlabbatch=matlabbatch(1:end-1); end
            
        case {'functional_coregister','functional_coregister_affine','functional_coregister_affine_noreslice'}
            jsubject=0;
            for isubject=1:numel(subjects),
                nsubject=subjects(isubject);
                if CONN_x.Setup.structural_sessionspecific, nsess_struct=intersect(sessions,1:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)));
                else nsess_struct=1;
                end
                for nses_struct=nsess_struct(:)'
                    if CONN_x.Setup.structural_sessionspecific, nsess_func=nses_struct;
                    else nsess_func=intersect(sessions,1:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)));
                    end
                    if ~isempty(nsess_func)
                        jsubject=jsubject+1;
                        filename=conn_get_functional(nsubject,nsess_func(1),sets);
                        if isempty(filename), error('Functional data not yet defined for subject %d session %d',nsubject,nsess_func(1)); end
                        temp=cellstr(filename);
                        if numel(temp)==1, ttemp=cellstr(conn_expandframe(temp{1})); else ttemp=temp; end
                        if coregtomean==2
                            if ~isempty(coregsource)&&iscell(coregsource)&&numel(coregsource)>=isubject
                                xtemp={coregsource{isubject}};
                            elseif numel(CONN_x.Setup.coregsource_functional)>=nsubject
                                xtemp=CONN_x.Setup.coregsource_functional{nsubject}(1);
                            else error('missing coregsource info');
                            end
                        elseif coregtomean,
                            [xtemp,failed]=conn_setup_preproc_meanimage(temp{1});
                            if isempty(xtemp),  errmsg=['Error preparing files for coregistration. Mean functional file not found (associated with functional data ',failed,' ; run either realignment or ART first, or select the option "first functional volume as reference" -#coregtomean field- if you prefer to use the first functional volume instead of the mean functional volume as reference)']; conn_disp(errmsg); error(errmsg); end
                            xtemp={xtemp};
                        else xtemp=ttemp(1);
                        end
                        matlabbatch{end+1}.spm.spatial.coreg.estimate.source=xtemp;
                        if isempty(CONN_x.Setup.structural{nsubject}{nses_struct}{1}), error(sprintf('No structural file defined for Subject %d Session %d. Please select a structural file',nsubject,nses_struct));
                        elseif numel(CONN_x.Setup.structural{nsubject}{nses_struct}{3})>1, error(sprintf('Multiple structural files found for Subject %d Session %d. Please select a single structural file',nsubject,nses_struct));
                        end
                        matlabbatch{end}.spm.spatial.coreg.estimate.ref=CONN_x.Setup.structural{nsubject}{nses_struct}(1);
                        if 0,%coregtomean, matlabbatch{end}.spm.spatial.coreg.estimate.other=xtemp;
                        else matlabbatch{end}.spm.spatial.coreg.estimate.other={};
                        end
                        for nsestrue=nsess_func(:)'
                            filename=conn_get_functional(nsubject,nsestrue,sets);
                            if isempty(filename), error('Functional data not yet defined for subject %d session %d',nsubject,nsestrue); end
                            temp=cellstr(filename);
                            if numel(temp)==1, ttemp=cellstr(conn_expandframe(temp{1})); else ttemp=temp; end
                            matlabbatch{end}.spm.spatial.coreg.estimate.other=cat(1,matlabbatch{end}.spm.spatial.coreg.estimate.other,ttemp);
                        end
                    end
                end
            end
            if ~jsubject, matlabbatch=matlabbatch(1:end-1); end
            
        case 'functional_segment'
            DOSPM12=~PREFERSPM8OVERSPM12&spmver12; %SPM12/SPM8
            if DOSPM12, matlabbatch{end+1}.spm.spatial.preproc.channel.vols={};
            else matlabbatch{end+1}.spm.spatial.preproc.data={};
            end
            for isubject=1:numel(subjects),
                nsubject=subjects(isubject);
                filename=conn_get_functional(nsubject,sessions(1),sets);
                if isempty(filename), error('Functional data not yet defined for subject %d session %d',nsubject,sessions(1)); end
                temp=cellstr(filename);
                if numel(temp)==1, ttemp=cellstr(conn_expandframe(temp{1})); else ttemp=temp; end
                if coregtomean==2
                    if ~isempty(coregsource)&&iscell(coregsource)&&numel(coregsource)>=isubject
                        xtemp={coregsource{isubject}};
                    elseif numel(CONN_x.Setup.coregsource_functional)>=nsubject
                        xtemp=CONN_x.Setup.coregsource_functional{nsubject}(1);
                    else error('missing coregsource info');
                    end
                elseif coregtomean,
                    [xtemp,failed]=conn_setup_preproc_meanimage(temp{1});
                    if isempty(xtemp),  errmsg=['Error preparing files for normalization. Mean functional file not found (associated with functional data ',failed,' ; run either realignment or ART first, or select the option "first functional volume as reference" -#coregtomean field- if you prefer to use the first functional volume instead of the mean functional volume as reference)']; conn_disp(errmsg); error(errmsg); end
                    xtemp={xtemp};
                else xtemp=ttemp(1);
                end
                if DOSPM12,
                    matlabbatch{end}.spm.spatial.preproc.channel.vols{isubject}=xtemp{1};
                    outputfiles{isubject}{1}=xtemp{1};
                    outputfiles{isubject}{2}=conn_prepend('c1',xtemp{1},'.nii'); % note: fix SPM12 issue converting .img to .nii
                    outputfiles{isubject}{3}=conn_prepend('c2',xtemp{1},'.nii');
                    outputfiles{isubject}{4}=conn_prepend('c3',xtemp{1},'.nii');
                else
                    matlabbatch{end}.spm.spatial.preproc.data{isubject}=xtemp{1};
                    outputfiles{isubject}{1}=xtemp{1};
                    outputfiles{isubject}{2}=conn_prepend('c1',xtemp{1});
                    outputfiles{isubject}{3}=conn_prepend('c2',xtemp{1});
                    outputfiles{isubject}{4}=conn_prepend('c3',xtemp{1});
                end
            end
            if DOSPM12,
                if ~isempty(tpm_template), matlabbatch{end}.spm.spatial.preproc.tissue=conn_setup_preproc_tissue(tpm_template,tpm_ngaus,subjects,sessions); end
                if ~isempty(affreg), matlabbatch{end}.spm.spatial.preproc.warp.affreg=affreg; end
                matlabbatch{end}.spm.spatial.preproc.channel.vols=reshape(matlabbatch{end}.spm.spatial.preproc.channel.vols,[],1);
                matlabbatch{end}.spm.spatial.preproc.warp.write=[1 1];
            else
                if ~isempty(tpm_template),
                    if ~isempty(subjects)&&~isempty(sessions)&&(isnumeric(tpm_template)||(ischar(tpm_template)&&size(tpm_template,1)==1&&~isempty(conn_datasetlabel(tpm_template))))
                        error('unsupported subject-specific TPM in SPM8; please upgrade to SPM12 instead')
                    else temp=cellstr(conn_expandframe(tpm_template));
                    end
                    if isempty(tpm_ngaus), tpm_ngaus=[2 2 2 4]; end % grey/white/CSF (+other implicit)
                    matlabbatch{end}.spm.spatial.preproc.opts.tpm=temp;
                    matlabbatch{end}.spm.spatial.preproc.opts.ngaus=ngaus(1:numel(temp)+1);
                end
                matlabbatch{end}.spm.spatial.preproc.roptions.bb=boundingbox; %
                matlabbatch{end}.spm.spatial.preproc.output.GM=[0,0,1];
                matlabbatch{end}.spm.spatial.preproc.output.WM=[0,0,1];
                matlabbatch{end}.spm.spatial.preproc.output.CSF=[0,0,1];
            end
            
        case {'functional_normalize','functional_normalize_direct'}
            DOSPM12=~PREFERSPM8OVERSPM12&spmver12; %SPM12/SPM8
            if DOSPM12
                %note: functional_template disregarded (using tissue probability maps instead)
                matlabbatch{end+1}.spm.spatial.normalise.estwrite.woptions.bb=boundingbox;
                matlabbatch{end}.spm.spatial.normalise.estwrite.woptions.vox=voxelsize_func.*[1 1 1];
                if ~isempty(interp), matlabbatch{end}.spm.spatial.normalise.estwrite.woptions.interp=interp; end
                if ~isempty(tpm_template), [nill,matlabbatch{end}.spm.spatial.normalise.estwrite.eoptions.tpm]=conn_setup_preproc_tissue(tpm_template,tpm_ngaus,subjects,sessions(1)); end
                if ~isempty(affreg), matlabbatch{end}.spm.spatial.normalise.estwrite.eoptions.affreg=affreg; end
            else
                %note: tissue probability maps disregarded (using functional_template instead)
                matlabbatch{end+1}.spm.spatial.normalise.estwrite.roptions.bb=boundingbox;
                matlabbatch{end}.spm.spatial.normalise.estwrite.roptions.vox=voxelsize_func.*[1 1 1];
                matlabbatch{end}.spm.spatial.normalise.estwrite.eoptions.template={functional_template};
                if ~isempty(interp), matlabbatch{end}.spm.spatial.normalise.estwrite.eoptions.interp=interp; end
            end
            for isubject=1:numel(subjects),
                nsubject=subjects(isubject);
                filename=conn_get_functional(nsubject,sessions(1),sets);
                if isempty(filename), error('Functional data not yet defined for subject %d session %d',nsubject,sessions(1)); end
                temp=cellstr(filename);
                if numel(temp)==1, ttemp=cellstr(conn_expandframe(temp{1})); else ttemp=temp; end
                if coregtomean==2
                    if ~isempty(coregsource)&&iscell(coregsource)&&numel(coregsource)>=isubject
                        xtemp={coregsource{isubject}};
                    elseif numel(CONN_x.Setup.coregsource_functional)>=nsubject
                        xtemp=CONN_x.Setup.coregsource_functional{nsubject}(1);
                    else error('missing coregsource info');
                    end
                elseif coregtomean,
                    [xtemp,failed]=conn_setup_preproc_meanimage(temp{1});
                    if isempty(xtemp),  errmsg=['Error preparing files for normalization. Mean functional file not found (associated with functional data ',failed,' ; run either realignment or ART first, or select the option "first functional volume as reference" -#coregtomean field- if you prefer to use the first functional volume instead of the mean functional volume as reference)']; conn_disp(errmsg); error(errmsg); end
                    xtemp={xtemp};
                else xtemp=ttemp(1);
                end
                if DOSPM12, matlabbatch{end}.spm.spatial.normalise.estwrite.subj(isubject).vol=xtemp;
                else        matlabbatch{end}.spm.spatial.normalise.estwrite.subj(isubject).source=xtemp;
                end
                nsess=CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject));
                if coregtomean, matlabbatch{end}.spm.spatial.normalise.estwrite.subj(isubject).resample=xtemp;
                else matlabbatch{end}.spm.spatial.normalise.estwrite.subj(isubject).resample={};
                end
                for nses=1:nsess
                    if ismember(nses,sessions)
                        filename=conn_get_functional(nsubject,nses,sets);
                        if isempty(filename), error('Functional data not yet defined for subject %d session %d',nsubject,nses); end
                        temp=cellstr(filename);
                        if numel(temp)==1, ttemp=cellstr(conn_expandframe(temp{1})); else ttemp=temp; end
                        matlabbatch{end}.spm.spatial.normalise.estwrite.subj(isubject).resample=cat(1,matlabbatch{end}.spm.spatial.normalise.estwrite.subj(isubject).resample,ttemp);
                        outputfiles{isubject}{nses}=char(conn_prepend('w',temp));
                    end
                end
            end
            
        case {'functional_segment&normalize','functional_segment&normalize_direct'}
            DOSPM12=~PREFERSPM8OVERSPM12&spmver12; %SPM12/SPM8
            if DOSPM12, matlabbatch{end+1}.spm.spatial.preproc.channel.vols={};
            else matlabbatch{end+1}.spm.spatial.preproc.data={};
            end
            for isubject=1:numel(subjects),
                nsubject=subjects(isubject);
                filename=conn_get_functional(nsubject,sessions(1),sets);
                if isempty(filename), error('Functional data not yet defined for subject %d session %d',nsubject,sessions(1)); end
                temp=cellstr(filename);
                if numel(temp)==1, ttemp=cellstr(conn_expandframe(temp{1})); else ttemp=temp; end
                if coregtomean==2
                    if ~isempty(coregsource)&&iscell(coregsource)&&numel(coregsource)>=isubject
                        xtemp={coregsource{isubject}};
                    elseif numel(CONN_x.Setup.coregsource_functional)>=nsubject
                        xtemp=CONN_x.Setup.coregsource_functional{nsubject}(1);
                    else error('missing coregsource info');
                    end
                elseif coregtomean,
                    [xtemp,failed]=conn_setup_preproc_meanimage(temp{1});
                    if isempty(xtemp),  errmsg=['Error preparing files for normalization. Mean functional file not found (associated with functional data ',failed,' ; run either realignment or ART first, or select the option "first functional volume as reference" -#coregtomean field- if you prefer to use the first functional volume instead of the mean functional volume as reference)']; conn_disp(errmsg); error(errmsg); end
                    xtemp={xtemp};
                else xtemp=ttemp(1);
                end
                if DOSPM12,
                    matlabbatch{end}.spm.spatial.preproc.channel.vols{isubject}=xtemp{1};
                    outputfiles{isubject}{1}=xtemp{1};
                    outputfiles{isubject}{2}=conn_prepend('c1',xtemp{1},'.nii'); % note: fix SPM12 issue converting .img to .nii
                    outputfiles{isubject}{3}=conn_prepend('c2',xtemp{1},'.nii');
                    outputfiles{isubject}{4}=conn_prepend('c3',xtemp{1},'.nii');
                    outputfiles{isubject}{5}=conn_prepend('y_',xtemp{1},'.nii');
                else
                    matlabbatch{end}.spm.spatial.preproc.data{isubject}=xtemp{1};
                    outputfiles{isubject}{1}=xtemp{1};
                    outputfiles{isubject}{2}=conn_prepend('c1',xtemp{1});
                    outputfiles{isubject}{3}=conn_prepend('c2',xtemp{1});
                    outputfiles{isubject}{4}=conn_prepend('c3',xtemp{1});
                    outputfiles{isubject}{5}=conn_prepend('',xtemp{1},'_seg_sn.mat');
                end
            end
            if DOSPM12,
                if ~isempty(tpm_template), matlabbatch{end}.spm.spatial.preproc.tissue=conn_setup_preproc_tissue(tpm_template,tpm_ngaus,subjects,sessions); end
                if ~isempty(affreg), matlabbatch{end}.spm.spatial.preproc.warp.affreg=affreg; end
                matlabbatch{end}.spm.spatial.preproc.channel.vols=reshape(matlabbatch{end}.spm.spatial.preproc.channel.vols,[],1);
                matlabbatch{end}.spm.spatial.preproc.warp.write=[1 1];
                matlabbatch{end+1}.spm.spatial.normalise.write.woptions.bb=boundingbox;
                matlabbatch{end}.spm.spatial.normalise.write.woptions.vox=voxelsize_func.*[1 1 1];
                if ~isempty(interp), matlabbatch{end}.spm.spatial.normalise.write.woptions.interp=interp; end
            else
                if ~isempty(tpm_template),
                    if ~isempty(subjects)&&~isempty(sessions)&&(isnumeric(tpm_template)||(ischar(tpm_template)&&size(tpm_template,1)==1&&~isempty(conn_datasetlabel(tpm_template))))
                        error('unsupported subject-specific TPM in SPM8; please upgrade to SPM12 instead')
                    else temp=cellstr(conn_expandframe(tpm_template));
                    end
                    if isempty(tpm_ngaus), tpm_ngaus=[2 2 2 4]; end % grey/white/CSF (+other implicit)
                    matlabbatch{end}.spm.spatial.preproc.opts.tpm=temp;
                    matlabbatch{end}.spm.spatial.preproc.opts.ngaus=ngaus(1:numel(temp)+1);
                end
                matlabbatch{end}.spm.spatial.preproc.roptions.bb=boundingbox; %
                matlabbatch{end}.spm.spatial.preproc.output.GM=[0,0,1];
                matlabbatch{end}.spm.spatial.preproc.output.WM=[0,0,1];
                matlabbatch{end}.spm.spatial.preproc.output.CSF=[0,0,1];
                matlabbatch{end+1}.spm.spatial.normalise.write.roptions.bb=boundingbox;
                matlabbatch{end}.spm.spatial.normalise.write.roptions.vox=voxelsize_func.*[1 1 1];
                if ~isempty(interp), matlabbatch{end}.spm.spatial.normalise.write.roptions.interp=interp; end
            end
            for isubject=1:numel(subjects),
                nsubject=subjects(isubject);
                if DOSPM12, matlabbatch{end}.spm.spatial.normalise.write.subj(isubject).def=outputfiles{isubject}(5);
                else        matlabbatch{end}.spm.spatial.normalise.write.subj(isubject).matname=outputfiles{isubject}(5);
                end
                matlabbatch{end}.spm.spatial.normalise.write.subj(isubject).resample=outputfiles{isubject}(1:4)';
                outputfiles{isubject}=outputfiles{isubject}(1:4);
                nsess=CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject));
                for nses=1:nsess
                    if ismember(nses,sessions)
                        filename=conn_get_functional(nsubject,nses,sets);
                        if isempty(filename), error('Functional data not yet defined for subject %d session %d',nsubject,nses); end
                        temp=cellstr(filename);
                        if numel(temp)==1, ttemp=cellstr(conn_expandframe(temp{1})); else ttemp=temp; end
                        matlabbatch{end}.spm.spatial.normalise.write.subj(isubject).resample=cat(1,matlabbatch{end}.spm.spatial.normalise.write.subj(isubject).resample,ttemp);
                        outputfiles{isubject}{4+nses}=char(conn_prepend('w',temp));
                    end
                end
                outputfiles{isubject}{1}=conn_prepend('w',outputfiles{isubject}{1});
                outputfiles{isubject}{2}=conn_prepend('w',outputfiles{isubject}{2});
                outputfiles{isubject}{3}=conn_prepend('w',outputfiles{isubject}{3});
                outputfiles{isubject}{4}=conn_prepend('w',outputfiles{isubject}{4});
            end
            
        case 'functional_smooth_masked'
            if size(fwhm,1)>1&&~iscell(fwhm), fwhm=num2cell(fwhm,2); end
            if iscell(fwhm), this_fwhm=fwhm{1}; fwhm=fwhm(2:end);
            else this_fwhm=fwhm;
            end
            if isempty(this_fwhm)
                this_fwhm=inputdlg('Enter smoothing FWHM (in mm)','conn_setup_preproc',1,{num2str(8)});
                if isempty(this_fwhm), return; end
                this_fwhm=str2num(this_fwhm{1});
            end
            for isubject=1:numel(subjects),
                nsubject=subjects(isubject);
                nsess=CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject));
                SVARIANT=1; % smoothing variant: 1: hard transitions (out-of-mask voxels are always unchanged); 2: smooth transitions (out-of-mask voxels are changed)
                LAMBDA=5;
                for nses=1:nsess
                    if ismember(nses,sessions)
                        filename=conn_get_functional(nsubject,nses,sets);
                        if isempty(filename), error('Functional data not yet defined for subject %d session %d',nsubject,nses); end
                        temp=cellstr(filename);
                        if CONN_x.Setup.structural_sessionspecific, nses_struct=nses;
                        else nses_struct=1;
                        end
                        fmask=CONN_x.Setup.rois.files{nsubject}{1}{nses_struct}{1};
                        if isempty(fmask), error('Grey Matter ROI data not yet defined for subject %d session %d',nsubject,nses); end
                        vmask=spm_vol(fmask);
                        vol=spm_vol(char(temp));
                        volout=vol;
                        %voloutmask=vol(1); voloutmask.fname=conn_prepend('p',vol(1).fname); voloutmask.mat=vol(1).mat; voloutmask.dim=vol(1).dim; voloutmask.pinfo=[1;0;0]; voloutmask.dt=[spm_type('float32') spm_platform('bigend')];
                        for n=1:numel(vol),volout(n).fname=conn_prepend('m',vol(n).fname); volout(n).mat=vol(1).mat; volout(n).dim=vol(1).dim; volout(n).pinfo=[1;0;0]; end
                        [gridx,gridy,gridz]=ndgrid(1:vol(1).dim(1),1:vol(1).dim(2),1:vol(1).dim(3));xyz=vol(1).mat*[gridx(:),gridy(:),gridz(:),ones(numel(gridx),1)]';
                        tempout=conn_prepend('m',temp); spm_unlink(tempout{:});
                        volout=spm_create_vol(volout);
                        vox=sqrt(sum(vol(1).mat(1:3,1:3).^2));
                        mask=max(0,min(1,double(reshape(spm_get_data(vmask,pinv(vmask.mat)*xyz),vol(1).dim))));
                        smask=zeros(size(mask));
                        spm_smooth(mask,smask,[1 1 1].*this_fwhm./vox);
                        mask0=(1-smask).^LAMBDA;
                        for n=1:numel(vol)
                            data=double(reshape(spm_get_data(vol(n),pinv(vol(n).mat)*xyz),vol(1).dim));
                            sdata=zeros(size(data));
                            spm_smooth(data.*mask,sdata,[1 1 1].*this_fwhm./vox);
                            switch(SVARIANT)
                                case 1, volout(n)=spm_write_vol(volout(n),(1-mask).*data + mask.*sdata./max(eps,smask));
                                case 2, volout(n)=spm_write_vol(volout(n),(sdata+data.*mask0)./max(eps,smask+mask0));
                                    %case 3, volout(n)=spm_write_vol(volout(n),(sdata.*smask+1/LAMBDA*data)./(smask+1/LAMBDA));
                            end
                        end
                        %voloutmask=spm_write_vol(voloutmask,mask);
                        outputfiles{isubject}{nses}=char(conn_prepend('m',temp));
                        %outputfiles{isubject}{nses}{2}=conn_prepend('p',vol(1).fname);
                    end
                end
            end
            
        case 'functional_smooth'
            if size(fwhm,1)>1&&~iscell(fwhm), fwhm=num2cell(fwhm,2); end
            if iscell(fwhm), this_fwhm=fwhm{1}; fwhm=fwhm(2:end);
            else this_fwhm=fwhm;
            end
            if isempty(this_fwhm)
                this_fwhm=inputdlg('Enter smoothing FWHM (in mm)','conn_setup_preproc',1,{num2str(8)});
                if isempty(this_fwhm), return; end
                this_fwhm=str2num(this_fwhm{1});
            end
            matlabbatch{end+1}.spm.spatial.smooth.fwhm=[1 1 1].*this_fwhm;
            matlabbatch{end}.spm.spatial.smooth.data={};
            for isubject=1:numel(subjects),
                nsubject=subjects(isubject);
                nsess=CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject));
                for nses=1:nsess
                    if ismember(nses,sessions)
                        filename=conn_get_functional(nsubject,nses,sets);
                        if isempty(filename), error('Functional data not yet defined for subject %d session %d',nsubject,nses); end
                        temp=cellstr(filename);
                        matlabbatch{end}.spm.spatial.smooth.data=cat(1,matlabbatch{end}.spm.spatial.smooth.data,temp);
                        outputfiles{isubject}{nses}=char(conn_prepend('s',temp));
                    end
                end
            end
            
        case 'functional_motionmask'
            for isubject=1:numel(subjects),
                nsubject=subjects(isubject);
                nsess=CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject));
                for nses=1:nsess
                    if ismember(nses,sessions)
                        filename=conn_get_functional(nsubject,nses,sets);
                        if isempty(filename), error('Functional data not yet defined for subject %d session %d',nsubject,nses); end
                        temp=cellstr(filename);
                        outputfiles{isubject}{nses}=temp;
                    end
                end
            end
            
        case 'functional_surface_resample'
            for isubject=1:numel(subjects),
                nsubject=subjects(isubject);
                nsess=CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject));
                for nses=1:nsess
                    if ismember(nses,sessions)
                        if CONN_x.Setup.structural_sessionspecific, nsess_struct=nses;
                        else nsess_struct=1;
                        end
                        filename=conn_get_functional(nsubject,nses,sets);
                        filestruct=deblank(CONN_x.Setup.structural{nsubject}{nsess_struct}{1});
                        assert(~isempty(filename), 'Functional data not yet defined for subject %d session %d',nsubject,nses);
                        assert(~isempty(filestruct), 'Structural data not yet defined for subject %d session %d',nsubject,nses);
                        assert(conn_checkFSfiles(filestruct),'No FreeSurfer data found for structural file %s',filestruct);
                        outputfiles{isubject}{nses}=conn_surf_resample(filename,filestruct);
                    end
                end
            end
            
        case 'functional_surface_smooth'
            if size(diffusionsteps,1)>1&&~iscell(diffusionsteps), diffusionsteps=num2cell(diffusionsteps,2); end
            if iscell(diffusionsteps), this_diffusionsteps=diffusionsteps{1}; diffusionsteps=diffusionsteps(2:end);
            else this_diffusionsteps=diffusionsteps;
            end
            if isempty(this_diffusionsteps)
                this_diffusionsteps=inputdlg('Enter number of diffusion steps for smoothing','conn_setup_preproc',1,{num2str(10)});
                if isempty(this_diffusionsteps), return; end
                this_diffusionsteps=str2num(this_diffusionsteps{1});
            end
            %matlabbatch{end+1}.functional_surface_smooth.diffusionsteps=this_diffusionsteps;
            %matlabbatch{end}.functional_surface_smooth.data={};
            for isubject=1:numel(subjects),
                nsubject=subjects(isubject);
                nsess=CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject));
                for nses=1:nsess
                    if ismember(nses,sessions)
                        filename=conn_get_functional(nsubject,nses,sets);
                        if isempty(filename), error('Functional data not yet defined for subject %d session %d',nsubject,nses); end
                        filenameout=conn_surf_smooth(filename,this_diffusionsteps);
                        outputfiles{isubject}{nses}=char(filenameout);
                        %matlabbatch{end}.functional_surface_smooth.data=cat(1,matlabbatch{end}.functional_surface_smooth.data,temp);
                        %outputfiles{isubject}{nses}=char(conn_prepend('s',temp));
                    end
                end
            end
            
        case 'functional_vdm_create'
            spm_jobman('initcfg');
            for isubject=1:numel(subjects),
                nsubject=subjects(isubject);
                nsess=CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject));
                VDM=[];
                for nses=1:nsess
                    if ismember(nses,sessions)
                        filename=conn_get_functional(nsubject,nses,sets);
                        if isempty(filename), error('Functional data not yet defined for subject %d session %d',nsubject,nses); end
                        if 1 % reference functional = 1st volume
                            filename=cellstr(filename);
                            if numel(filename)==1,filename=cellstr(conn_expandframe(filename{1}));end
                            filename=filename{1};
                            newfilename=regexprep(conn_prepend('ref_',filename),'\,\d+$','');
                            vol=spm_vol(filename);
                            val=spm_read_vols(vol(1));
                            mat0=vol(1).mat;
                            dim0=vol(1).dim;
                            newvol=struct('mat',vol(1).mat,'dim',vol(1).dim,'fname',newfilename,'pinfo',[1;0;0],'n',[1,1],'dt',[spm_type('float32') spm_platform('bigend')],'descrip','');
                            try, conn_setup_preproc_filedelete(newfilename); end
                            spm_write_vol(newvol,val);
                            filename=newfilename;
                        end
                        if isempty(vdm_fmap), vdm_fmap='fmap'; end
                        if ischar(vdm_fmap)&&isempty(conn_datasetlabel(vdm_fmap)), error('Dataset %s (fieldmap acquisition files) not found',vdm_fmap); end
                        fmap=conn_get_functional(nsubject,nses,vdm_fmap);
                        if ~isempty(fmap)
                            fmap=cellstr(fmap);
                            if numel(fmap)==1,fmap=cellstr(conn_expandframe(fmap{1})); end
                            ET1=vdm_et1; ET2=vdm_et2; ERT=vdm_ert; BLIP=vdm_blip;
                            if isequal(BLIP,0)||isequal(BLIP,'R')||isequal(BLIP,'r'), reverseBLIP=true; BLIP=[]; else reverseBLIP=false; end
                            if isequal(vdm_type,1)||(isempty(vdm_type)&&(numel(fmap)==2||numel(fmap)==3)), % Magnitude1+PhaseDiff or Magnitude1+Magnitude2+PhaseDiff
                                fmap=[fmap(end) fmap(1)]; % note: sorts as PhaseDiff+Magnitude for SPM FieldMap_create
                                scphase=FieldMap('Scale',fmap{1});
                                fmap{1}=scphase.fname; % scaled phase
                                vol=spm_vol(char(fmap));
                                phaseabs=spm_read_vols(vol);
                                pmean=angle(sum(sum(sum(phaseabs(:,:,:,2).*exp(1i*phaseabs(:,:,:,1)),1),2),3));
                                pdemeaned=angle(exp(1i*(phaseabs(:,:,:,1)-pmean))); % note: demean phase info to avoid global shifts
                                newvol=struct('mat',vol(1).mat,'dim',vol(1).dim,'fname',fmap{1},'pinfo',[1;0;0],'n',[1,1],'dt',[spm_type('float32') spm_platform('bigend')],'descrip','');
                                try, conn_setup_preproc_filedelete(fmap{1}); end
                                spm_write_vol(newvol,pdemeaned);
                                if isempty(ERT), ERT=1000*conn_jsonread(filename,'TotalReadoutTime'); end
                                if isempty(ERT), ERT=1000./conn_jsonread(filename,'BandwidthPerPixelPhaseEncode'); end
                                if nses==1&&isempty(ERT), conn_disp('fprintf','warning: unable to find TotalReadoutTime or BandwidthPerPixelPhaseEncode information in %s\n',filename); end
                                if isempty(ET1), ET1=1000*conn_jsonread(fmap{1},'EchoTime1'); end
                                if isempty(ET1), ET1=1000*conn_jsonread(fmap{2},'EchoTime1'); end
                                if isempty(ET1), ET1=1000*conn_jsonread(fmap{2},'EchoTime'); end
                                %if isempty(ET1), ET1=1000*conn_jsonread(fmap{2},'EchoTime'); end
                                if nses==1&&isempty(ET1), conn_disp('fprintf','warning: unable to find EchoTime1 or EchoTime information in %s\n',fmap{2}); end %ET1=4.37;
                                if isempty(ET2), ET2=1000*conn_jsonread(fmap{1},'EchoTime2'); end
                                if isempty(ET2), ET2=1000*conn_jsonread(fmap{2},'EchoTime2'); end
                                if isempty(ET2), ET2=ET1+1000*conn_jsonread(fmap{1},'EchoTimeDifference'); end %ETDIFF=2.46
                                if isempty(ET2), ET2=ET1+1000*conn_jsonread(fmap{2},'EchoTimeDifference'); end
                                if isempty(ET2), ET2=1000*conn_jsonread(fmap{1},'EchoTime'); end
                                if nses==1&&isempty(ET2), conn_disp('fprintf','warning: unable to find EchoTime2 or EchoTimeDifference information in %s\n',fmap{1}); end
                                if isempty(BLIP),
                                    BLIP=conn_jsonread(filename,'PhaseEncodingDirection',false);
                                    if iscell(BLIP), BLIP=char(BLIP); end
                                    if isequal(BLIP,'i+')||isequal(BLIP,'i'), BLIP=sign([0 1 0 0]*mat0*[1 0 0 0]');
                                    elseif isequal(BLIP,'i-'), BLIP=sign([0 1 0 0]*mat0*[-1 0 0 0]');
                                    elseif isequal(BLIP,'j+')||isequal(BLIP,'j'), BLIP=sign([0 1 0 0]*mat0*[0 1 0 0]');
                                    elseif isequal(BLIP,'j-'), BLIP=sign([0 1 0 0]*mat0*[0 -1 0 0]');
                                    elseif isequal(BLIP,'k+')||isequal(BLIP,'k'), BLIP=sign([0 1 0 0]*mat0*[0 0 1 0]');
                                    elseif isequal(BLIP,'k-'), BLIP=sign([0 1 0 0]*mat0*[0 0 -1 0]');
                                    elseif ~isempty(BLIP), error('unable to interpret PhaseEncodingDirection %s (expected ''j+'' or ''j-'' directions)',BLIP);
                                    end
                                    if reverseBLIP, BLIP=-BLIP; end
                                end
                                if nses==1&&isempty(BLIP), conn_disp('fprintf','warning: unable to find PhaseEncodingDirection information in %s\n',filename); end
                                if ~isempty(ET1)&&~isempty(ET2)&&~isempty(ERT)&&~isempty(BLIP)
                                    conn_disp('fprintf','Creating vdm file for subject %d session %d...\n',nsubject,nses);
                                    if ET1>ET2, [ET1,ET2]=deal(ET2,ET1); end
                                    pm_defaults;
                                    pm_def.sessname='session'; pm_def.SHORT_ECHO_TIME=ET1; pm_def.LONG_ECHO_TIME=ET2; pm_def.TOTAL_EPI_READOUT_TIME=ERT; pm_def.INPUT_DATA_FORMAT='PM'; pm_def.EPI_BASED_FIELDMAPS=0; pm_def.K_SPACE_TRAVERSAL_BLIP_DIR=BLIP; pm_def.MASKBRAIN=1; pm_def.match_vdm=1; %pm_def.write_unwarped=1;
                                    conn_disp('fprintf','   Phase : %s\n   Magnitude : %s\n   reference : %s\n',fmap{1},fmap{2},filename);
                                    conn_setup_preproc_disp(pm_def,'   options');
                                    vfmap=spm_vol(char(fmap)); vfmap(3:end)=[]; % note: disregards any additional magnitude images
                                    VDM = FieldMap_create(vfmap,{filename},pm_def); %[ET1,ET2,0,ERT,-1]
                                    try
                                        tvol=spm_vol(char(VDM{1}.fname));
                                        tdat=spm_read_vols(tvol);
                                        tvol.fname=conn_prepend('reverse_',tvol.fname);
                                        spm_write_vol(tvol,-tdat);
                                    end
                                else error('insufficient information for vdm creation. Skipping subject %d session %d...\n',nsubject,nses);
                                end
                            elseif isequal(vdm_type,3)||(isempty(vdm_type)&&numel(fmap)==1), % FieldMap [note: work in progress; needs further testing]
                                units=conn_jsonread(fmap{1},'Units',false); 
                                newfmap1=conn_prepend('',fmap{1},['_session',num2str(nses),'.nii']);
                                if isempty(units)||~ischar(units),
                                    conn_disp('fprintf','units of FieldMap %s not found, assuming Hz\n',fmap{1});
                                    ct=1;
                                else
                                    switch(lower(units))
                                        case 'hz', ct=1;
                                        case 'rad/s', ct=1/2/pi;
                                        case 'tesla', ct=2.6752219e8/2/pi;
                                        otherwise, error('unrecognized units %s found in FieldMap %s (expected Hz, rad/s, or Tesla)',units,fmap{1});
                                    end
                                end
                                if 1
                                    vol=spm_vol(fmap{1});
                                    val=ct*spm_read_vols(vol);
                                    fmap{1}=conn_prepend('sc',fmap{1});
                                    newvol=struct('mat',vol(1).mat,'dim',vol(1).dim,'fname',newfmap1,'pinfo',[1;0;0],'n',[1,1],'dt',[spm_type('float32') spm_platform('bigend')],'descrip','');
                                    try, conn_setup_preproc_filedelete(newfmap1); end
                                    spm_write_vol(newvol,val);
                                end
                                fmap{1}=newfmap1; % avoids potential issues when using the same file across different sessions
                                if isempty(ERT), ERT=1000*conn_jsonread(filename,'TotalReadoutTime'); end
                                if isempty(ERT), ERT=1000./conn_jsonread(filename,'BandwidthPerPixelPhaseEncode'); end
                                if isempty(ERT), 
                                    PED=conn_jsonread(filename,'PhaseEncodingDirection',false); 
                                    if iscell(PED), PED=char(PED); end
                                    if ~isempty(PED)&&ischar(PED)&&any(PED(1)=='ijk'), ERT=1000*conn_jsonread(filename,'EffectiveEchoSpacing')*dim0(PED(1)-'h'); end
                                end
                                if nses==1&&isempty(ERT), conn_disp('fprintf','warning: unable to find TotalReadoutTime information (or BandwidthPerPixelPhaseEncode or EffectiveEchoSpacing) in %s\n',filename); end
                                if isempty(BLIP),
                                    PED=conn_jsonread(filename,'PhaseEncodingDirection',false);
                                    if iscell(PED), PED=char(PED); end
                                    if isequal(PED,'i+')||isequal(PED,'i'), BLIP=sign([0 1 0 0]*mat0*[1 0 0 0]');
                                    elseif isequal(PED,'i-'), BLIP=sign([0 1 0 0]*mat0*[-1 0 0 0]');
                                    elseif isequal(PED,'j+')||isequal(PED,'j'), BLIP=sign([0 1 0 0]*mat0*[0 1 0 0]');
                                    elseif isequal(PED,'j-'), BLIP=sign([0 1 0 0]*mat0*[0 -1 0 0]');
                                    elseif isequal(PED,'k+')||isequal(PED,'k'), BLIP=sign([0 1 0 0]*mat0*[0 0 1 0]');
                                    elseif isequal(PED,'k-'), BLIP=sign([0 1 0 0]*mat0*[0 0 -1 0]');
                                    elseif ~isempty(PED), error('unable to interpret PhaseEncodingDirection %s (expected ''j+'' or ''j-'' directions)',PED);
                                    else BLIP=[];
                                    end
                                    if reverseBLIP, BLIP=-BLIP; end
                                end
                                if nses==1&&isempty(BLIP), conn_disp('fprintf','warning: unable to find PhaseEncodingDirection information in %s\n',filename); end
                                if ~isempty(ERT)&&~isempty(BLIP)
                                    conn_disp('fprintf','Creating vdm file for subject %d session %d...\n',nsubject,nses);
                                    pm_defaults;
                                    pm_def.sessname='session'; pm_def.TOTAL_EPI_READOUT_TIME=ERT; pm_def.EPI_BASED_FIELDMAPS=0; pm_def.K_SPACE_TRAVERSAL_BLIP_DIR=BLIP; pm_def.MASKBRAIN=1; pm_def.match_vdm=0; %pm_def.match_vdm=1; %pm_def.write_unwarped=1;
                                    conn_disp('fprintf','   FieldMap : %s\n   ref : %s\n',fmap{1},filename);
                                    conn_setup_preproc_disp(pm_def,'   options');
                                    VDM = FieldMap_create(char(fmap),{filename},pm_def); %[ET1,ET2,0,ERT,-1]                                    
                                    try
                                        tvol=spm_vol(char(VDM{1}.fname));
                                        tdat=spm_read_vols(tvol);
                                        tvol.fname=conn_prepend('reverse_',tvol.fname);
                                        spm_write_vol(tvol,-tdat);
                                    end
                                else error('insufficient information for vdm creation. Skipping subject %d session %d...\n',nsubject,nses);
                                end
                            elseif isequal(vdm_type,2)||(isempty(vdm_type)&&numel(fmap)==4), % Real1+Imag1+Real2+Imag2 [note: work in progress; needs further testing]
                                if isempty(ERT), ERT=1000*conn_jsonread(filename,'TotalReadoutTime'); end
                                if isempty(ERT), ERT=1000./conn_jsonread(filename,'BandwidthPerPixelPhaseEncode'); end
                                if nses==1&&isempty(ERT), conn_disp('fprintf','warning: unable to find TotalReadoutTime or BandwidthPerPixelPhaseEncode information in %s\n',filename); end
                                if isempty(ET1), ET1=1000*conn_jsonread(fmap{1},'EchoTime'); end
                                if isempty(ET1), ET1=1000*conn_jsonread(fmap{2},'EchoTime'); end
                                if nses==1&&isempty(ET1), conn_disp('fprintf','warning: unable to find EchoTime information in %s or %s\n',fmap{1},fmap{2}); end %ET1=4.37;
                                if isempty(ET2), ET2=1000*conn_jsonread(fmap{3},'EchoTime'); end
                                if isempty(ET2), ET2=1000*conn_jsonread(fmap{4},'EchoTime'); end
                                if nses==1&&isempty(ET2), conn_disp('fprintf','warning: unable to find EchoTime information in %s or %s\n',fmap{3},fmap{4}); end %ET1=4.37;
                                if isempty(BLIP),
                                    BLIP=conn_jsonread(filename,'PhaseEncodingDirection',false);
                                    if iscell(BLIP), BLIP=char(BLIP); end
                                    if isequal(BLIP,'i+')||isequal(BLIP,'i'), BLIP=sign([0 1 0 0]*mat0*[1 0 0 0]');
                                    elseif isequal(BLIP,'i-'), BLIP=sign([0 1 0 0]*mat0*[-1 0 0 0]');
                                    elseif isequal(BLIP,'j+')||isequal(BLIP,'j'), BLIP=sign([0 1 0 0]*mat0*[0 1 0 0]');
                                    elseif isequal(BLIP,'j-'), BLIP=sign([0 1 0 0]*mat0*[0 -1 0 0]');
                                    elseif isequal(BLIP,'k+')||isequal(BLIP,'k'), BLIP=sign([0 1 0 0]*mat0*[0 0 1 0]');
                                    elseif isequal(BLIP,'k-'), BLIP=sign([0 1 0 0]*mat0*[0 0 -1 0]');
                                    elseif ~isempty(BLIP), error('unable to interpret PhaseEncodingDirection %s (expected ''j+'' or ''j-'' directions)',BLIP);
                                    end
                                    if reverseBLIP, BLIP=-BLIP; end
                                end
                                if nses==1&&isempty(BLIP), conn_disp('fprintf','warning: unable to find PhaseEncodingDirection information in %s\n',filename); end
                                if ~isempty(ET1)&&~isempty(ET2)&&~isempty(ERT)&&~isempty(BLIP)
                                    conn_disp('fprintf','Creating vdm file for subject %d session %d...\n',nsubject,nses);
                                    pm_defaults;
                                    pm_def.sessname='session'; pm_def.SHORT_ECHO_TIME=ET1; pm_def.LONG_ECHO_TIME=ET2; pm_def.TOTAL_EPI_READOUT_TIME=ERT; pm_def.INPUT_DATA_FORMAT='RI'; pm_def.EPI_BASED_FIELDMAPS=0; pm_def.K_SPACE_TRAVERSAL_BLIP_DIR=BLIP; pm_def.MASKBRAIN=1; pm_def.match_vdm=1; %pm_def.write_unwarped=1;
                                    conn_disp('fprintf','   Real1 : %s\n   Imag1 : %s\n   Real2 : %s\n   Imag2 : %s\n   ref : %s\n',fmap{1},fmap{2},fmap{3},fmap{4},filename);
                                    conn_setup_preproc_disp(pm_def,'   options');
                                    VDM = FieldMap_create(char(fmap),{filename},pm_def); %[ET1,ET2,0,ERT,-1]
                                    try
                                        tvol=spm_vol(char(VDM{1}.fname));
                                        tdat=spm_read_vols(tvol);
                                        tvol.fname=conn_prepend('reverse_',tvol.fname);
                                        spm_write_vol(tvol,-tdat);
                                    end
                                else error('insufficient information for vdm creation. Skipping subject %d session %d...\n',nsubject,nses);
                                end
                            else error('type of fieldmap sequence files could not be determined from number of available files (%d) in dataset %s for subject %d session %d',numel(fmap),mat2str(vdm_fmap),nsubject,nses);
                            end
                        end
                        if isempty(VDM), outputfiles{isubject}{nses}='';
                        else outputfiles{isubject}{nses}=char(VDM{1}.fname);
                        end
                    end
                end
            end
            
        otherwise
            error(['unrecognized option ',STEP]);
    end
    
    
    if dogui&&ishandle(hmsg), delete(hmsg); end
    hmsg=[];
    newmatlabbatch={}; newin=[];
    for n=1:numel(matlabbatch) % fix to allow subject-specific TPMs
        if isfield(matlabbatch{n},'spm')&&isfield(matlabbatch{n}.spm,'spatial')&&isfield(matlabbatch{n}.spm.spatial,'preproc')&&isfield(matlabbatch{n}.spm.spatial.preproc,'tissue')&&iscell(matlabbatch{n}.spm.spatial.preproc.tissue)
            for n2=1:numel(matlabbatch{n}.spm.spatial.preproc.tissue)
                newmatlabbatch{end+1}=matlabbatch{n};
                newmatlabbatch{end}.spm.spatial.preproc.tissue=matlabbatch{n}.spm.spatial.preproc.tissue{n2};
                newmatlabbatch{end}.spm.spatial.preproc.channel.vols=matlabbatch{n}.spm.spatial.preproc.channel.vols(n2);
                newin(end+1)=n;
            end
        elseif isfield(matlabbatch{n},'spm')&&isfield(matlabbatch{n}.spm,'spatial')&&isfield(matlabbatch{n}.spm.spatial,'normalise')&&isfield(matlabbatch{n}.spm.spatial.normalise,'estwrite')&&isfield(matlabbatch{n}.spm.spatial.normalise.estwrite,'eoptions')&&isfield(matlabbatch{n}.spm.spatial.normalise.estwrite.eoptions,'tpm')&&size(matlabbatch{n}.spm.spatial.normalise.estwrite.eoptions.tpm,2)>1
            for n2=1:size(matlabbatch{n}.spm.spatial.normalise.estwrite.eoptions.tpm,2)
                newmatlabbatch{end+1}=matlabbatch{n};
                newmatlabbatch{end}.spm.spatial.normalise.estwrite.eoptions.tpm=matlabbatch{n}.spm.spatial.normalise.estwrite.eoptions.tpm(:,n2);
                newmatlabbatch{end}.spm.spatial.normalise.estwrite.subj=matlabbatch{n}.spm.spatial.normalise.estwrite.subj(n2);
                newin(end+1)=n;
            end
        end
    end
    if ~isempty(newmatlabbatch)
        keepin=setdiff(1:numel(matlabbatch),newin);
        newin=[keepin newin];
        matlabbatch=[matlabbatch(keepin), newmatlabbatch];
        [nill,idxin]=sort(newin);
        matlabbatch=matlabbatch(idxin);
    end
    
    
    if strncmp(lower(STEP),'interactive_',numel('interactive_'))
        doimport=false;
        if any(strcmpi(regexprep(lower(STEP),'^run_|^update_|^interactive_',''),{'functional_art'}))
            for n=1:numel(matlabbatch)
                conn_art('sess_file',matlabbatch{n}.art);
            end
        elseif any(strcmpi(regexprep(lower(STEP),'^run_|^update_|^interactive_',''),{'functional_removescans','functional_bandpass','functional_regression','functional_manualorient','structural_manualorient','functional_center','functional_centertostruct','structural_center','functional_motionmask','functional_label','functional_label_as_original', 'functional_label_as_subjectspace', 'functional_label_as_mnispace', 'functional_label_as_surfacespace', 'functional_label_as_smoothed','functional_load','functional_load_from_original', 'functional_load_from_subjectspace', 'functional_load_from_mnispace', 'functional_load_from_surfacespace', 'functional_load_from_smoothed','functional_surface_smooth','functional_surface_resample','functional_vdm_create'}))
        elseif ~isempty(matlabbatch)
            spm_jobman('initcfg');
            try, spm_get_defaults('mat.format','-v7.3'); end
            try
                job_id=spm_jobman('interactive',matlabbatch);
                % outputs=cfg_util('getAllOutputs', job_id)
            catch
                ok=-1;
            end
        end
    else %if strncmp(lower(STEP),'run_',numel('run_'))
        if dogui, hmsg=conn_msgbox({['Performing ',STEP_name],'Please wait...'},'');
        else conn_disp(['Performing ',STEP_name,'. Please wait...']);
        end
        if any(strcmpi(regexprep(lower(STEP),'^run_|^interactive_',''),{'functional_art'}))
            conn_disp('fprintf','\nART preprocessing job\n');
            conn_setup_preproc_disp(matlabbatch);
            for n=1:numel(matlabbatch)
                if art_force_interactive, h=conn_art('sess_file',matlabbatch{n}.art);
                else h=conn_art('sess_file',matlabbatch{n}.art,'visible','off');
                end
                if strcmp(get(h,'name'),'art'), %close(h);
                elseif strcmp(get(gcf,'name'),'art'), h=gcf;%close(gcf);
                else h=findobj(0,'name','art'); %close(h);
                end
                if art_force_interactive, uiwait(h);
                else
                    %try
                    %    if isfield(matlabbatch{n}.art,'output_dir')
                    %        saveas(h,fullfile(matlabbatch{n}.art.output_dir,'art_screenshot.fig'));
                    %    end
                    %end
                    try
                        if isfield(matlabbatch{n}.art,'output_dir')
                            conn_print(h,fullfile(matlabbatch{n}.art.output_dir,'art_screenshot.jpg'),'-nogui');
                        end
                        close(h);
                    end
                end
            end
        elseif any(strcmpi(regexprep(lower(STEP),'^run_|^update_|^interactive_',''),{'functional_motionmask'}))
            for isubject=1:numel(outputfiles),
                for nses=1:numel(outputfiles{isubject})
                    if ismember(nses,sessions)
                        outputfiles{isubject}{nses}=conn_computeMaskMovement(outputfiles{isubject}{nses});
                    end
                end
            end
        elseif any(strcmpi(regexprep(lower(STEP),'^run_|^update_|^interactive_',''),{'functional_removescans','functional_bandpass','functional_regression','functional_manualorient','structural_manualorient','functional_center','functional_centertostruct','structural_center','functional_motionmask','functional_label','functional_label_as_original', 'functional_label_as_subjectspace', 'functional_label_as_mnispace', 'functional_label_as_surfacespace', 'functional_label_as_smoothed','functional_load','functional_load_from_original', 'functional_load_from_subjectspace', 'functional_load_from_mnispace', 'functional_load_from_surfacespace', 'functional_load_from_smoothed'}))
        elseif strncmp(lower(STEP),'update_',numel('update_'))
        elseif ~isempty(matlabbatch)
            spm_jobman('initcfg');
            try, spm_get_defaults('mat.format','-v7.3'); end
            conn_disp('fprintf','\nSPM preprocessing job\n');
            conn_setup_preproc_disp(matlabbatch);
            debugskip=false;
            if ~debugskip
                warning('off','MATLAB:RandStream:ActivatingLegacyGenerators');
                job_id=spm_jobman('run',matlabbatch);
                warning('on','MATLAB:RandStream:ActivatingLegacyGenerators');
            end
        end
        if dogui&&ishandle(hmsg), delete(hmsg);
        else conn_disp(['Done ',STEP_name]);
        end
        ok=1;
    end
    if ishandle(hmsg), delete(hmsg); end
    
    if ok>=0&&doimport
        if dogui, hmsg=conn_msgbox({'Importing results to CONN project','Please wait...'},'');
        else conn_disp(['Importing results to CONN project. Please wait...']);
        end
        switch(regexprep(lower(STEP),'^run_|^update_|^interactive_',''))
            case 'functional_removescans'
                for isubject=1:numel(subjects),
                    nsubject=subjects(isubject);
                    for nses=1:numel(outputfiles{isubject})
                        if ismember(nses,sessions)
                            nV=conn_set_functional(nsubject,nses,sets,outputfiles{isubject}{nses}{1});
                            if ~sets
                                for nl1covariate=1:numel(outputfiles{isubject}{nses})-1
                                    if size(outputfiles{isubject}{nses}{1+nl1covariate},1)>0
                                        CONN_x.Setup.l1covariates.files{nsubject}{nl1covariate}{nses}={'[raw values]',[],outputfiles{isubject}{nses}{1+nl1covariate}};
                                    end
                                end
                            end
                        end
                    end
                end
                
            case {'functional_bandpass'}
                for isubject=1:numel(subjects),
                    nsubject=subjects(isubject);
                    for nses=1:numel(outputfiles{isubject})
                        if ismember(nses,sessions)
                            nV=conn_set_functional(nsubject,nses,sets,outputfiles{isubject}{nses}{1});
                        end
                    end
                end
                
            case {'functional_regression'}
                if ~sets||ALLSETSPERMISSIONS,
                    icov=find(strcmp(CONN_x.Setup.l1covariates.names(1:end-1),'QC_regressors'));
                    if isempty(icov),
                        icov=numel(CONN_x.Setup.l1covariates.names);
                        CONN_x.Setup.l1covariates.names{icov}='QC_regressors';
                        CONN_x.Setup.l1covariates.names{icov+1}=' ';
                    end
                else icov=[];
                end
                for isubject=1:numel(subjects),
                    nsubject=subjects(isubject);
                    for nses=1:numel(outputfiles{isubject})
                        if ismember(nses,sessions)
                            nV=conn_set_functional(nsubject,nses,sets,outputfiles{isubject}{nses}{1});
                            if ~sets||ALLSETSPERMISSIONS, CONN_x.Setup.l1covariates.files{nsubject}{icov}{nses}=conn_file(outputfiles{isubject}{nses}{2}); end
                        end
                    end
                end
                
            case 'functional_manualorient'
                if iscell(reorient), treorient=reorient{1}; reorient=reorient(2:end);
                else treorient=reorient;
                end
                if ischar(treorient)
                    R=load(treorient,'-mat');
                    fR=fieldnames(R);
                    treorient=R.(fR{1});
                end
                filename={};
                for isubject=1:numel(subjects),
                    nsubject=subjects(isubject);
                    nsess=CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject));
                    translation=[];
                    for nses=1:nsess
                        if ismember(nses,sessions)
                            [filename{nsubject}{nses},cfile]=conn_get_functional(nsubject,nses,sets,[],'cfile');
                            if isempty(filename{nsubject}{nses}), error('Functional data not yet defined for subject %d session %d',nsubject,nses); end
                            temp=cellstr(filename{nsubject}{nses});
                            if numel(temp)==1,
                                temp=cellstr(conn_expandframe(temp{1}));
                            end
                            if coregtomean, % keeps mean image in same space in case it is required later
                                [xtemp,failed]=conn_setup_preproc_meanimage(temp{1});
                                if ~isempty(xtemp), temp=[{xtemp};temp]; end
                            end
                            M=cell(1,numel(temp));
                            for n=1:numel(temp)
                                M{n}=spm_get_space(temp{n});
                                if isnan(treorient)
                                    if isempty(translation)
                                        translation=-M{n}(1:3,1:3)*cfile{3}(1).dim'/2 - M{n}(1:3,4);
                                        try, R=eye(4);R(1:3,4)=translation; save(conn_prepend('centering_',temp{n},'.mat'),'R','-mat'); R=inv(R); save(conn_prepend('icentering_',temp{n},'.mat'),'R','-mat'); end
                                    end
                                    M{n}(1:3,4)=M{n}(1:3,4)+translation;
                                else
                                    if size(treorient,1)==4, R=treorient;
                                    else R=[treorient zeros(3,1); zeros(1,3) 1];
                                    end
                                    M{n}=R*M{n};
                                    if n==1, try, save(conn_prepend('reorient_',temp{n},'.mat'),'R','-mat'); R=inv(R); save(conn_prepend('ireorient_',temp{n},'.mat'),'R','-mat'); end; end
                                end
                            end
                            for n=1:numel(temp)
                                spm_get_space(temp{n},M{n});
                            end
                        end
                    end
                end
                for isubject=1:numel(subjects),
                    nsubject=subjects(isubject);
                    nsess=CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject));
                    for nses=1:nsess
                        if ismember(nses,sessions)
                            nV=conn_set_functional(nsubject,nses,sets,filename{nsubject}{nses});
                        end
                    end
                end
                
            case 'structural_manualorient'
                SAVETODIFFERENTFILE=true;
                if iscell(reorient), treorient=reorient{1}; reorient=reorient(2:end);
                else treorient=reorient;
                end
                if ischar(treorient)
                    R=load(treorient,'-mat');
                    fR=fieldnames(R);
                    treorient=R.(fR{1});
                end
                jsubject=0;
                for isubject=1:numel(subjects),
                    nsubject=subjects(isubject);
                    if CONN_x.Setup.structural_sessionspecific, nsess=CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)); else nsess=1; end
                    for nses=1:nsess
                        if ismember(nses,sessions)
                            jsubject=jsubject+1;
                            if isempty(CONN_x.Setup.structural{nsubject}{nses}{1}), error(sprintf('No structural file defined for Subject %d Session %d. Please select a structural file',nsubject,nses));
                            elseif numel(CONN_x.Setup.structural{nsubject}{nses}{3})>1, error(sprintf('Multiple structural files found for Subject %d Session %d. Please select a single structural file',nsubject,nses));
                            end
                            temp=CONN_x.Setup.structural{nsubject}{nses}{1};
                            M=spm_get_space(temp);
                            if isnan(treorient)
                                translation=-M(1:3,1:3)*CONN_x.Setup.structural{nsubject}{nses}{3}(1).dim'/2 - M(1:3,4);
                                M(1:3,4)=M(1:3,4)+translation;
                                try, R=eye(4);R(1:3,4)=translation; save(conn_prepend('centering_',temp,'.mat'),'R','-mat'); R=inv(R); save(conn_prepend('icentering_',temp,'.mat'),'R','-mat'); end
                            else
                                if size(treorient,1)==4, R=treorient;
                                else R=[treorient zeros(3,1); zeros(1,3) 1];
                                end
                                M=R*M;
                                try, save(conn_prepend('reorient_',temp,'.mat'),'R','-mat'); R=inv(R); save(conn_prepend('ireorient_',temp,'.mat'),'R','-mat'); end
                            end
                            if SAVETODIFFERENTFILE
                                a=spm_vol(temp);
                                b=spm_read_vols(a);
                                temp=conn_prepend('c',temp);
                                a.fname=regexprep(temp,',\d+$','');
                                spm_write_vol(a,b);
                            end
                            spm_get_space(temp,M);
                            [CONN_x.Setup.structural{nsubject}{nses},nV]=conn_file(temp);
                        end
                    end
                    if ~CONN_x.Setup.structural_sessionspecific, CONN_x.Setup.structural{nsubject}(2:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)))=CONN_x.Setup.structural{nsubject}(1); end
                end
                
            case {'functional_manualspatialdef'}
                for isubject=1:numel(subjects),
                    nsubject=subjects(isubject);
                    nsess=intersect(1:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)),sessions);
                    for nses=nsess(:)'
                        conn_set_functional(nsubject,nses,sets,outputfiles{isubject}{nses}{1});
                    end
                end
                
            case {'structural_manualspatialdef'}
                for isubject=1:numel(subjects),
                    nsubject=subjects(isubject);
                    nsess=intersect(1:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)),sessions);
                    for nses=nsess(:)'
                        if CONN_x.Setup.structural_sessionspecific, nses_struct=nses;
                        else nses_struct=1;
                        end
                        CONN_x.Setup.structural{nsubject}{nses}=conn_file(outputfiles{isubject}{nses_struct}{1});
                    end
                end
                
            case 'functional_label_as_original'
                conn_datasetcopy(sets,'original data',subjects);
            case 'functional_label_as_subjectspace'
                conn_datasetcopy(sets,'subject-space data',subjects);
            case 'functional_label_as_mnispace'
                conn_datasetcopy(sets,'mni-space data',subjects);
            case 'functional_label_as_surfacespace'
                conn_datasetcopy(sets,'surface-space data',subjects);
            case 'functional_label_as_smoothed'
                conn_datasetcopy(sets,'smoothed data',subjects);
            case 'functional_label',
                if iscell(label), this_label=label{1}; label=label(2:end);
                else this_label=label;
                end
                conn_datasetcopy(sets,this_label,subjects);
            case 'functional_load_from_original'
                conn_datasetcopy('original data',sets,subjects);
            case 'functional_load_from_subjectspace'
                conn_datasetcopy('subject-space data',sets,subjects);
            case 'functional_load_from_mnispace'
                conn_datasetcopy('mni-space data',sets,subjects);
            case 'functional_load_from_surfacespace'
                conn_datasetcopy('surface-space data',sets,subjects);
            case 'functional_load_from_smoothed'
                conn_datasetcopy('smoothed data',sets,subjects);
            case 'functional_load',
                if iscell(load_label), this_label=load_label{1}; load_label=load_label(2:end);
                else this_label=load_label;
                end
                conn_datasetcopy(this_label,sets,subjects);
                
            case 'functional_center'
                treorient=nan;
                filename={};
                for isubject=1:numel(subjects),
                    nsubject=subjects(isubject);
                    nsess=CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject));
                    translation=[];
                    for nses=1:nsess
                        if ismember(nses,sessions)
                            [filename{nsubject}{nses},cfile]=conn_get_functional(nsubject,nses,sets,[],'cfile');
                            if isempty(filename{nsubject}{nses}), error('Functional data not yet defined for subject %d session %d',nsubject,nses); end
                            temp=cellstr(filename{nsubject}{nses});
                            if numel(temp)==1,
                                temp=cellstr(conn_expandframe(temp{1}));
                            end
                            if coregtomean, % keeps mean image in same space in case it is required later
                                [xtemp,failed]=conn_setup_preproc_meanimage(temp{1});
                                if ~isempty(xtemp), temp=[{xtemp};temp]; end
                            end
                            M=cell(1,numel(temp));
                            for n=1:numel(temp)
                                M{n}=spm_get_space(temp{n});
                                if isempty(translation)
                                    translation=-M{n}(1:3,1:3)*cfile{3}(1).dim'/2 - M{n}(1:3,4);
                                    conn_disp('fprintf','Functional centering translation x/y/z = %s (Subject %d)\n',mat2str(translation'),nsubject);
                                    try, R=eye(4);R(1:3,4)=translation; save(conn_prepend('centering_',temp{n},'.mat'),'R','-mat'); R=inv(R); save(conn_prepend('icentering_',temp{n},'.mat'),'R','-mat'); end
                                end
                                M{n}(1:3,4)=M{n}(1:3,4)+translation;
                            end
                            for n=1:numel(temp)
                                spm_get_space(temp{n},M{n});
                            end
                        end
                    end
                end
                for isubject=1:numel(subjects),
                    nsubject=subjects(isubject);
                    nsess=CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject));
                    for nses=1:nsess
                        if ismember(nses,sessions)
                            nV=conn_set_functional(nsubject,nses,sets,filename{nsubject}{nses});
                        end
                    end
                end
                
            case 'functional_centertostruct'
                treorient=nan;
                filename={};
                for isubject=1:numel(subjects),
                    nsubject=subjects(isubject);
                    nsess=CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject));
                    translation=[];
                    for nses=1:nsess
                        if ismember(nses,sessions)
                            if CONN_x.Setup.structural_sessionspecific, nsess_struct=nses;
                            else nsess_struct=1;
                            end
                            if isempty(CONN_x.Setup.structural{nsubject}{nsess_struct}{1}), error(sprintf('No structural file defined for Subject %d Session %d. Please select a structural file',nsubject,nsess_struct));
                            elseif numel(CONN_x.Setup.structural{nsubject}{nsess_struct}{3})>1, error(sprintf('Multiple structural files found for Subject %d Session %d. Please select a single structural file',nsubject,nsess_struct));
                            end
                            cfile_struct=CONN_x.Setup.structural{nsubject}{nsess_struct};
                            [filename{nsubject}{nses},cfile]=conn_get_functional(nsubject,nses,sets,[],'cfile');
                            if isempty(filename{nsubject}{nses}), error('Functional data not yet defined for subject %d session %d',nsubject,nses); end
                            temp=cellstr(filename{nsubject}{nses});
                            if numel(temp)==1,
                                temp=cellstr(conn_expandframe(temp{1}));
                            end
                            if coregtomean, % keeps mean image in same space in case it is required later
                                [xtemp,failed]=conn_setup_preproc_meanimage(temp{1});
                                if ~isempty(xtemp), temp=[{xtemp};temp]; end
                            end
                            M=cell(1,numel(temp));
                            for n=1:numel(temp)
                                M{n}=spm_get_space(temp{n});
                                if isempty(translation)
                                    M0=spm_get_space(cfile_struct{1});
                                    translation=-M{n}(1:3,1:3)*cfile{3}(1).dim'/2 - M{n}(1:3,4) +M0(1:3,1:3)*cfile_struct{3}(1).dim'/2 + M0(1:3,4);
                                    conn_disp('fprintf','Functional centering translation x/y/z = %s (Subject %d)\n',mat2str(translation'),nsubject);
                                    try, R=eye(4);R(1:3,4)=translation; save(conn_prepend('centeringtostruct_',temp{n},'.mat'),'R','-mat'); R=inv(R); save(conn_prepend('icentering_',temp{n},'.mat'),'R','-mat'); end
                                end
                                M{n}(1:3,4)=M{n}(1:3,4)+translation;
                            end
                            for n=1:numel(temp)
                                spm_get_space(temp{n},M{n});
                            end
                        end
                    end
                end
                for isubject=1:numel(subjects),
                    nsubject=subjects(isubject);
                    nsess=CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject));
                    for nses=1:nsess
                        if ismember(nses,sessions)
                            nV=conn_set_functional(nsubject,nses,sets,filename{nsubject}{nses});
                        end
                    end
                end
                
            case 'structural_center'
                SAVETODIFFERENTFILE=true;
                treorient=nan;
                jsubject=0;
                for isubject=1:numel(subjects),
                    nsubject=subjects(isubject);
                    if CONN_x.Setup.structural_sessionspecific, nsess=CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)); else nsess=1; end
                    for nses=1:nsess
                        if ismember(nses,sessions)
                            jsubject=jsubject+1;
                            if isempty(CONN_x.Setup.structural{nsubject}{nses}{1}), error(sprintf('No structural file defined for Subject %d Session %d. Please select a structural file',nsubject,nses));
                            elseif numel(CONN_x.Setup.structural{nsubject}{nses}{3})>1, error(sprintf('Multiple structural files found for Subject %d Session %d. Please select a single structural file',nsubject,nses));
                            end
                            temp=CONN_x.Setup.structural{nsubject}{nses}{1};
                            M=spm_get_space(temp);
                            translation=-M(1:3,1:3)*CONN_x.Setup.structural{nsubject}{nses}{3}(1).dim'/2 - M(1:3,4);
                            M(1:3,4)=M(1:3,4)+translation;
                            conn_disp('fprintf','Structural centering translation x/y/z = %s (Subject %d)\n',mat2str(translation'),nsubject);
                            try, R=eye(4);R(1:3,4)=translation; save(conn_prepend('centering_',temp,'.mat'),'R','-mat'); R=inv(R); save(conn_prepend('icentering_',temp,'.mat'),'R','-mat'); end
                            %M(1:3,4)=-M(1:3,1:3)*CONN_x.Setup.structural{nsubject}{nses}{3}(1).dim'/2;
                            if SAVETODIFFERENTFILE
                                a=spm_vol(temp);
                                b=spm_read_vols(a);
                                temp=conn_prepend('c',temp);
                                a.fname=regexprep(temp,',\d+$','');
                                spm_write_vol(a,b);
                            end
                            spm_get_space(temp,M);
                            [CONN_x.Setup.structural{nsubject}{nses},nV]=conn_file(temp);
                        end
                    end
                    if ~CONN_x.Setup.structural_sessionspecific, CONN_x.Setup.structural{nsubject}(2:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)))=CONN_x.Setup.structural{nsubject}(1); end
                end
                
            case 'structural_segment'
                for isubject=1:numel(subjects),
                    nsubject=subjects(isubject);
                    nsess=intersect(1:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)),sessions);
                    for nses=nsess(:)'
                        if CONN_x.Setup.structural_sessionspecific, nses_struct=nses;
                        else nses_struct=1;
                        end
                        CONN_x.Setup.structural{nsubject}{nses}=conn_file(outputfiles{isubject}{nses_struct}{1});
                        if ~sets||ALLSETSPERMISSIONS
                            CONN_x.Setup.rois.files{nsubject}{1}{nses}=conn_file(outputfiles{isubject}{nses_struct}{2});
                            CONN_x.Setup.rois.files{nsubject}{2}{nses}=conn_file(outputfiles{isubject}{nses_struct}{3});
                            CONN_x.Setup.rois.files{nsubject}{3}{nses}=conn_file(outputfiles{isubject}{nses_struct}{4});
                        end
                    end
                end
                
            case {'structural_segment&normalize','functional_segment&normalize_indirect'}
                for isubject=1:numel(subjects),
                    nsubject=subjects(isubject);
                    nsess=intersect(1:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)),sessions);
                    for nses=nsess(:)'
                        if CONN_x.Setup.structural_sessionspecific, nses_struct=nses;
                        else nses_struct=1;
                        end
                        CONN_x.Setup.structural{nsubject}{nses}=conn_file(outputfiles{isubject}{nses_struct}{1});
                        if ~sets||ALLSETSPERMISSIONS
                            CONN_x.Setup.rois.files{nsubject}{1}{nses}=conn_file(outputfiles{isubject}{nses_struct}{2});
                            CONN_x.Setup.rois.files{nsubject}{2}{nses}=conn_file(outputfiles{isubject}{nses_struct}{3});
                            CONN_x.Setup.rois.files{nsubject}{3}{nses}=conn_file(outputfiles{isubject}{nses_struct}{4});
                        end
                        if applytofunctional||strcmp(regexprep(lower(STEP),'^run_|^update_|^interactive_',''),'functional_segment&normalize_indirect'),
                            conn_set_functional(nsubject,nses,sets,outputfiles{isubject}{nses}{6});
                        end
                    end
                end
                
            case 'functional_coregister_nonlinear'
                for isubject=1:numel(subjects),
                    nsubject=subjects(isubject);
                    nsess=intersect(1:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)),sessions);
                    for nses=nsess(:)'
                        if CONN_x.Setup.structural_sessionspecific, nses_struct=nses;
                        else nses_struct=1;
                        end
                        %CONN_x.Setup.structural{nsubject}{nses}=conn_file(outputfiles{isubject}{nses_struct}{1});
                        if ~sets||ALLSETSPERMISSIONS
                            CONN_x.Setup.rois.files{nsubject}{1}{nses}=conn_file(outputfiles{isubject}{nses_struct}{2});
                            CONN_x.Setup.rois.files{nsubject}{2}{nses}=conn_file(outputfiles{isubject}{nses_struct}{3});
                            CONN_x.Setup.rois.files{nsubject}{3}{nses}=conn_file(outputfiles{isubject}{nses_struct}{4});
                        end
                        conn_set_functional(nsubject,nses,sets,outputfiles{isubject}{nses}{7});
                    end
                end
                
            case {'structural_normalize','structural_normalize_preservemasks','functional_normalize_indirect','functional_normalize_indirect_preservemasks'}
                for isubject=1:numel(subjects),
                    nsubject=subjects(isubject);
                    nsess=intersect(1:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)),sessions);
                    for nses=nsess(:)'
                        if CONN_x.Setup.structural_sessionspecific, nses_struct=nses;
                        else nses_struct=1;
                        end
                        CONN_x.Setup.structural{nsubject}{nses}=conn_file(outputfiles{isubject}{nses_struct}{1});
                        if applytofunctional||strcmp(regexprep(lower(STEP),'^run_|^update_|^interactive_',''),'functional_normalize_indirect'),
                            conn_set_functional(nsubject,nses,sets,outputfiles{isubject}{nses}{2});
                        end
                        if strcmp(regexprep(lower(STEP),'^run_|^update_|^interactive_',''),'functional_normalize_indirect_preservemasks')||strcmp(regexprep(lower(STEP),'^run_|^update_|^interactive_',''),'structural_normalize_preservemasks')
                            if ismember(nses_struct,sessions)
                                CONN_x.Setup.rois.files{nsubject}{1}{nses}=conn_file(outputfiles{isubject}{nses_struct}{3}{1});
                                CONN_x.Setup.rois.files{nsubject}{2}{nses}=conn_file(outputfiles{isubject}{nses_struct}{3}{2});
                                CONN_x.Setup.rois.files{nsubject}{3}{nses}=conn_file(outputfiles{isubject}{nses_struct}{3}{3});
                            end
                        end
                    end
                end
                
            case 'functional_segment'
                if ~sets||ALLSETSPERMISSIONS
                    for isubject=1:numel(subjects),
                        nsubject=subjects(isubject);
                        nsess=intersect(1:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)),sessions);
                        for nses=nsess(:)'
                            CONN_x.Setup.rois.files{nsubject}{1}{nses}=conn_file(outputfiles{isubject}{2});
                            CONN_x.Setup.rois.files{nsubject}{2}{nses}=conn_file(outputfiles{isubject}{3});
                            CONN_x.Setup.rois.files{nsubject}{3}{nses}=conn_file(outputfiles{isubject}{4});
                        end
                    end
                end
                
            case {'functional_segment&normalize','functional_segment&normalize_direct'}
                for isubject=1:numel(subjects),
                    nsubject=subjects(isubject);
                    nsess=intersect(1:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)),sessions);
                    for nses=nsess(:)'
                        if ~sets||ALLSETSPERMISSIONS
                            CONN_x.Setup.rois.files{nsubject}{1}{nses}=conn_file(outputfiles{isubject}{2});
                            CONN_x.Setup.rois.files{nsubject}{2}{nses}=conn_file(outputfiles{isubject}{3});
                            CONN_x.Setup.rois.files{nsubject}{3}{nses}=conn_file(outputfiles{isubject}{4});
                        end
                        conn_set_functional(nsubject,nses,sets,outputfiles{isubject}{4+nses});
                    end
                end
                
            case {'functional_slicetime','functional_normalize','functional_normalize_direct','functional_smooth','functional_smooth_masked','functional_surface_smooth','functional_surface_resample'}
                for isubject=1:numel(subjects),
                    nsubject=subjects(isubject);
                    nsess=intersect(1:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)),sessions);
                    for nses=nsess(:)'
                        conn_set_functional(nsubject,nses,sets,outputfiles{isubject}{nses});
                    end
                end
                
            case 'functional_vdm_create'
                if ~sets||ALLSETSPERMISSIONS,
                    for isubject=1:numel(subjects),
                        nsubject=subjects(isubject);
                        nsess=intersect(1:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)),sessions);
                        for nses=nsess(:)'
                            CONN_x.Setup.unwarp_functional{nsubject}{nses}=conn_file(outputfiles{isubject}{nses});
                            conn_set_functional(nsubject,nses,'vdm',outputfiles{isubject}{nses});
                        end
                    end
                end
                
            case 'functional_art'
                if ~sets||ALLSETSPERMISSIONS,
                    icov=find(strcmp(CONN_x.Setup.l1covariates.names(1:end-1),'QC_timeseries'));
                    if isempty(icov),
                        icov=numel(CONN_x.Setup.l1covariates.names);
                        CONN_x.Setup.l1covariates.names{icov}='QC_timeseries';
                        CONN_x.Setup.l1covariates.names{icov+1}=' ';
                    end
                    for isubject=1:numel(subjects),
                        nsubject=subjects(isubject);
                        nsess=intersect(1:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)),sessions);
                        for nses=nsess(:)'
                            CONN_x.Setup.l1covariates.files{nsubject}{icov}{nses}=conn_file(outputfiles{isubject}{nses}{2});
                        end
                    end
                    icov0=icov;
                    icov=find(strcmp(CONN_x.Setup.l1covariates.names(1:end-1),'scrubbing'));
                    if isempty(icov),
                        icov=numel(CONN_x.Setup.l1covariates.names);
                        CONN_x.Setup.l1covariates.names{icov}='scrubbing';
                        CONN_x.Setup.l1covariates.names{icov+1}=' ';
                    end
                    for isubject=1:numel(subjects),
                        nsubject=subjects(isubject);
                        nsess=intersect(1:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)),sessions);
                        for nses=nsess(:)'
                            CONN_x.Setup.l1covariates.files{nsubject}{icov}{nses}=conn_file(outputfiles{isubject}{nses}{1});
                        end
                    end
                    MEANSOVERVALIDONLY=true; % note: switch to computing MeanMotion and MeanGSchange only over valid scans
                    y1=zeros(CONN_x.Setup.nsubjects,1);y2=zeros(CONN_x.Setup.nsubjects,1);y3=nan(CONN_x.Setup.nsubjects,1);y4=zeros(CONN_x.Setup.nsubjects,1);y5=nan(CONN_x.Setup.nsubjects,1);y6=zeros(CONN_x.Setup.nsubjects,1);yok=true;
                    for isubject=1:numel(subjects),
                        nsubject=subjects(isubject);
                        nsess=CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject));
                        for nses=1:nsess
                            try
                                temp=load(CONN_x.Setup.l1covariates.files{nsubject}{icov}{nses}{1});
                                validscans=~any(temp.R~=0,2);
                                y1(nsubject)=y1(nsubject)+sum(validscans,1);        % ValidScans
                                y2(nsubject)=y2(nsubject)+sum(sum(temp.R~=0));      % InvalidScans
                                try,
                                    temp=load(CONN_x.Setup.l1covariates.files{nsubject}{icov0}{nses}{1});
                                catch
                                    temp=load(regexprep(CONN_x.Setup.l1covariates.files{nsubject}{icov}{nses}{1},'art_regression_outliers_(and_movement_)?','art_regression_timeseries_'));
                                end
                                y3(nsubject)=max(y3(nsubject), max(abs(temp.R(:,end)),[],1) );      % MaxMotion
                                y5(nsubject)=max(y5(nsubject), max(abs(temp.R(:,end-1)),[],1) );    % MaxGSchange
                                if MEANSOVERVALIDONLY
                                    y4(nsubject)=y4(nsubject)+sum(abs(temp.R(validscans,end)),1);   % MeanMotion
                                    y6(nsubject)=y6(nsubject)+sum(abs(temp.R(validscans,end-1)),1); % MeanGSchange
                                else
                                    y4(nsubject)=y4(nsubject)+mean(abs(temp.R(:,end)),1);           % MeanMotion
                                    y6(nsubject)=y6(nsubject)+mean(abs(temp.R(:,end-1)),1);         % MeanGSchange
                                end
                            catch
                                yok=false;
                            end
                        end
                        if MEANSOVERVALIDONLY
                            y4(nsubject)=y4(nsubject)/y1(nsubject);
                            y6(nsubject)=y6(nsubject)/y1(nsubject);
                        else
                            y4(nsubject)=y4(nsubject)/nsess;
                            y6(nsubject)=y6(nsubject)/nsess;
                        end
                    end
                    if yok,
                        str_global=sprintf(' (outliers threshold = %s)',mat2str(art_global_threshold));
                        str_motion=sprintf(' (outliers threshold = %s)',mat2str(art_motion_threshold));
                        conn_importl2covariate({'QC_ValidScans','QC_InvalidScans','QC_MaxMotion','QC_MeanMotion','QC_MaxGSchange','QC_MeanGSchange'},{y1,y2,y3,y4,y5,y6},0,subjects,{'CONN Quality Assurance: Number of valid (non-outlier) scans','CONN Quality Assurance: Number of outlier scans',['CONN Quality Assurance: Largest motion observed',str_motion],['CONN Quality Assurance: Average motion observed (disregarding outlier scans)',str_motion],['CONN Quality Assurance: Largest global BOLD signal changes observed',str_global],['CONN Quality Assurance: Average global BOLD signal changes observed (disregarding outlier scans)',str_global]});
                    end
                end
                
            case {'functional_realign','functional_realign&unwarp','functional_realign&unwarp&fieldmap','functional_realign&unwarp&phasemap','functional_realign_noreslice'}
                icov=find(strcmp(CONN_x.Setup.l1covariates.names(1:end-1),'realignment'));
                if isempty(icov)&&(~sets||ALLSETSPERMISSIONS),
                    icov=numel(CONN_x.Setup.l1covariates.names);
                    CONN_x.Setup.l1covariates.names{icov}='realignment';
                    CONN_x.Setup.l1covariates.names{icov+1}=' ';
                end
                for isubject=1:numel(subjects),
                    nsubject=subjects(isubject);
                    nsess=intersect(1:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)),sessions);
                    for nses=nsess(:)'
                        conn_set_functional(nsubject,nses,sets,outputfiles{isubject}{nses}{1});
                        if ~sets||ALLSETSPERMISSIONS
                            CONN_x.Setup.l1covariates.files{nsubject}{icov}{nses}=conn_file(outputfiles{isubject}{nses}{2});
                        end
                    end
                end
                
            case {'functional_vdm_apply'}
                for isubject=1:numel(subjects),
                    nsubject=subjects(isubject);
                    nsess=intersect(1:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)),sessions);
                    for nses=nsess(:)'
                        conn_set_functional(nsubject,nses,sets,outputfiles{isubject}{nses}{1});
                    end
                end
                
                
            case {'functional_coregister','functional_coregister_affine','functional_coregister_affine_noreslice','functional_coregister_affine_reslice'}
                filename={};
                for isubject=1:numel(subjects),
                    nsubject=subjects(isubject);
                    nsess=intersect(1:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)),sessions);
                    for nses=nsess(:)'
                        filename{nsubject}{nses}=conn_get_functional(nsubject,nses,sets);
                        if isempty(filename{nsubject}{nses}), error('Functional data not yet defined for subject %d session %d',nsubject,nses); end
                    end
                end
                for isubject=1:numel(subjects),
                    nsubject=subjects(isubject);
                    nsess=intersect(1:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)),sessions);
                    for nses=nsess(:)'
                        conn_set_functional(nsubject,nses,sets,filename{nsubject}{nses});
                    end
                end
                
            case 'functional_motionmask'
                if ~sets||ALLSETSPERMISSIONS
                    iroi=find(strcmp(CONN_x.Setup.rois.names(1:end-1),'MotionMask'),1);
                    if isempty(iroi),
                        iroi=numel(CONN_x.Setup.rois.names);
                        CONN_x.Setup.rois.names{iroi}='MotionMask';
                        CONN_x.Setup.rois.dimensions{iroi}=1;
                        CONN_x.Setup.rois.mask(iroi)=0;
                        CONN_x.Setup.rois.subjectspecific(iroi)=1;
                        CONN_x.Setup.rois.sessionspecific(iroi)=1;
                        CONN_x.Setup.rois.multiplelabels(iroi)=1;
                        CONN_x.Setup.rois.regresscovariates(iroi)=0;
                        CONN_x.Setup.rois.unsmoothedvolumes(iroi)=1;
                        CONN_x.Setup.rois.weighted(iroi)=1;
                        CONN_x.Setup.rois.names{iroi+1}=' ';
                    end
                    for isubject=1:numel(subjects),
                        nsubject=subjects(isubject);
                        nsess=intersect(1:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)),sessions);
                        for nses=nsess(:)'
                            CONN_x.Setup.rois.files{nsubject}{iroi}{nses}=conn_file(outputfiles{isubject}{nses});
                        end
                    end
                end
        end
        if dogui&&ishandle(hmsg), delete(hmsg); end
        ok=2;
        if isfield(CONN_x,'filename')&&~isempty(CONN_x.filename), conn save; end
    end
    
    if ok<0, return; end
end
if ~dogui, conn_disp('Done'); end
end


function conn_setup_preproc_update(hdl)
if ~nargin, hdl=gcbo; end
dlg=get(hdl,'userdata');
val=get(dlg.m0,'value');
val2=dlg.steps_order(val);
str=get(dlg.m0,'string');
%if any(ismember(cat(1,str(val),get(dlg.m7,'string')),{'structural Segmentation & Normalization','structural Normalization'})),
%    if any(ismember(str(val),{'structural Segmentation & Normalization','structural Normalization'})),
%        set(dlg.m3,'visible','on');
%    else set(dlg.m3,'visible','off');
%    end
%if any(ismember(cat(1,str(val),get(dlg.m7,'string')),{'functional Coregistration to structural','functional Normalization','functional Segmentation & Normalization','functional Segmentation'})),
%if any(ismember(str(val),{'functional Coregistration to structural','functional Normalization','functional Segmentation & Normalization','functional Segmentation'})),
if any(cellfun('length',regexp(str(val),'^functional Coregistration|^functional Normalization|^functional Segmentation|^functional Center|^functional Direct|^functional Indirect'))),
    set(dlg.m4,'visible','on');
else set(dlg.m4,'visible','off');
end
set(dlg.m6,'string',dlg.steps_descr{val2});
if ~isempty(dlg.m0b),
    nm7=numel(get(dlg.m7,'string'));
    vm7=get(dlg.m7,'value');
    if ~nm7,
        set([dlg.m6 dlg.m0b dlg.m4],'visible','off');
        set([dlg.m8b dlg.m8c dlg.m8d dlg.m8g],'enable','off');
    else
        if nm7==1||numel(vm7)~=1, set(dlg.m0b,'string',regexprep(dlg.steps_names{val2},'<.*?>',''));
        else set(dlg.m0b,'string',[sprintf('Step %d/%d: ',vm7,nm7) regexprep(dlg.steps_names{val2},'<.*?>','')]);
        end
        set([dlg.m6 dlg.m0b],'visible','on');
        set([dlg.m8b dlg.m8c dlg.m8d dlg.m8g],'enable','on');
    end
end
if get(dlg.m1,'value'), set(dlg.m5,'visible','off');
else set(dlg.m5,'visible','on');
end
if get(dlg.m1b,'value'), set(dlg.m5b,'visible','off');
else set(dlg.m5b,'visible','on');
end
if ~isempty(dlg.m7)&&isempty(get(dlg.m7,'string')), set(dlg.m11,'enable','off'); else set(dlg.m11,'enable','on'); end
end

function conn_setup_preproc_update_add(hdl,varargin)
if ~nargin, hdl=gcbo; end
dlg=get(hdl,'userdata');
val=listdlg('liststring',dlg.steps_names(dlg.steps_order),'selectionmode','single','promptstring',{'Select pipeline (in bold) or individual preprocessing step (in regular font)'},'ListSize',[800 400]);
if ~isempty(val)
    ival=dlg.steps_order(val);
    val=dlg.steps_index{ival};
    set(dlg.m7,'string',cat(1,get(dlg.m7,'string'),dlg.steps_names(val)'));
    tidx=find(strcmp(dlg.steps(dlg.steps_order),dlg.steps{val(1)}));
    if numel(tidx)==1, set(dlg.m0,'value',tidx); end
    conn_setup_preproc_update(hdl); %feval(get(dlg.m0,'callback'));
end
end

function conn_setup_preproc_save(hdl,varargin)
global CONN_x;
if ~nargin, hdl=gcbo; end
dlg=get(hdl,'userdata');
STEPS=get(dlg.m7,'string');
[tok,idx]=ismember(STEPS,dlg.steps_names);
STEPS=dlg.steps(idx(tok>0));
coregtomean=1;
if any(cellfun('length',regexp(STEPS,'^functional_coregister|^functional_normalize|functional_segment|functional_center'))), coregtomean=~get(dlg.m4,'value'); end
if nargin==1
    CONN_x.SetupPreproc.steps=STEPS;
    CONN_x.SetupPreproc.coregtomean=coregtomean;
elseif nargin==2&&ischar(varargin{1})&&conn_existfile(varargin{1})
    filename=varargin{1};
    save(filename,'STEPS','coregtomean');
    conn_disp('fprintf','Data preprocessing pipeline saved to %s\n',filename);
else
    [outputpathfilename, outputpathpathname]=uiputfile('*.mat','Save data preprocessing pipeline',fullfile(fileparts(which('conn')),'utils','preprocessingpipelines'));
    if ~ischar(outputpathfilename), return; end
    filename=fullfile(outputpathpathname,outputpathfilename);
    save(filename,'STEPS','coregtomean');
    conn_disp('fprintf','Data preprocessing pipeline saved to %s\n',filename);
end
end

function conn_setup_preproc_load(hdl,varargin)
global CONN_x
if ~nargin, hdl=gcbo; end
dlg=get(hdl,'userdata');
if nargin==1
    if ~isfield(CONN_x,'SetupPreproc')||~isfield(CONN_x.SetupPreproc,'steps')||~isfield(CONN_x.SetupPreproc,'coregtomean'), return; end
    STEPS=CONN_x.SetupPreproc.steps;
    coregtomean=CONN_x.SetupPreproc.coregtomean;
    filename='';
elseif nargin==2&&ischar(varargin{1})&&conn_existfile(varargin{1})
    filename=varargin{1};
    load(filename,'STEPS','coregtomean');
else
    [outputpathfilename,outputpathpathname]=uigetfile('*.mat','Select file',fullfile(fileparts(which('conn')),'utils','preprocessingpipelines'));
    if ~ischar(outputpathfilename), return; end
    filename=fullfile(outputpathpathname,outputpathfilename);
    load(filename,'STEPS','coregtomean');
end

if ~exist('STEPS','var'), conn_msgbox(sprintf('Problem loading file %s. Incorrect format',filename),'',2);
else
    [tok,idx]=ismember(STEPS,dlg.steps);
    if ~all(tok), conn_disp('Warning: some preprocessing steps do not have valid names'); end
    steps=dlg.steps(idx(tok>0));
    set(dlg.m7,'string',dlg.steps_names(idx(tok>0)));
    if exist('coregtomean','var'),
        if isempty(coregtomean), coregtomean=1; end
        set(dlg.m4,'value',~coregtomean);
    end
    if isempty(steps), set(dlg.m7,'value',[]);
    else
        set(dlg.m7,'value',1);
        tidx=find(strcmp(dlg.steps(dlg.steps_order),steps{1}));
        if numel(tidx)==1, set(dlg.m0,'value',tidx); conn_setup_preproc_update(dlg.m0); end
    end
end
end

function conn_setup_preproc_filedelete(filename)
if ispc, [ok,nill]=system(sprintf('del "%s"',filename));
else [ok,nill]=system(sprintf('rm ''%s''',filename));
end
end

function  [tissue,tpm]=conn_setup_preproc_tissue(tpm_template,tpm_ngaus,subjects,sessions)
global CONN_x;
if isnumeric(tpm_template)||(ischar(tpm_template)&&size(tpm_template,1)==1&&~isempty(conn_datasetlabel(tpm_template)))
    tpm_expanded={};
    tpm={};
    jsubject=0;
    for isubject=1:numel(subjects),
        nsubject=subjects(isubject);
        if CONN_x.Setup.structural_sessionspecific, nsess=CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsubject)); else nsess=1; end
        for nses=1:nsess
            if ismember(nses,sessions)
                jsubject=jsubject+1;
                tfile=conn_get_functional(nsubject,nses,tpm_template);
                tpm_expanded{end+1}=cellstr(conn_expandframe(tfile));
                tpm{end+1}=tfile;
            end
        end
    end
else
    tpm_expanded={cellstr(conn_expandframe(tpm_template))};
    tpm=reshape(cellstr(tpm_template),[],1);
end
tissue={};
for ifile=1:numel(tpm_expanded),
    file=tpm_expanded{ifile};
    if isempty(tpm_ngaus), tpm_ngaus=[1 1 2 3 4 2]; end % grey/white/CSF/bone/soft/air
    if numel(tpm_ngaus)<numel(file), tpm_ngaus=[tpm_ngaus(:)' 4+zeros(1,numel(file)-numel(tpm_ngaus))]; end
    for n=1:numel(file)
        tissue{ifile}(n)=struct('tpm',{file(n)},'ngaus',tpm_ngaus(min(n,numel(tpm_ngaus))),'native',[1 0],'warped',[0 0]);
    end
end
if numel(tissue)==1, tissue=tissue{1}; end % regular case, single tpm file 
end








