function varargout=conn_batch(varargin)
% CONN BATCH batch functionality for connectivity toolbox
% 
% Defines experiment information and/or run processing steps programmatically
% 
% CONN_BATCH syntax:
% 
% 1) conn_batch(BATCH);
%    where BATCH is a structure (with fields defined in the section below)
%    e.g. 
%       clear BATCH;
%       BATCH.Setup.RT=2;
%       conn_batch(BATCH);
% 
% 2) conn_batch('fieldname1',fieldvalue1,'fieldname2',fieldvalue2,...) 
%    where 'fieldname*' are individual BATCH structure fields
%    e.g. 
%       conn_batch('Setup.RT',2); 
% 
% 3) conn_batch(batchfilename)
%    where batchfilename is a .mat file containing a batch structure
%    or a .json file containing a batch structure
%    or a .m file containing a Matlab script
%    e.g.
%       conn_batch('mybatchfile.mat');
%
% note1: in standalone releases use syntax (from system-prompt):
%     conn batch batchfilename    : runs a batch file (.m or .mat)
%     conn batch "matlabcommands" : runs one or several matlab commands 
% note2: syntax conn_batch({BATCH1 BATCH2 ...}) processes sequentially multiple batch structures (equivalent to 
%       conn_batch(BATCH1); conn_batch(BATCH2); ... e.g. useful when defining/running multiple first- or 
%       second- level analyses)
%
%__________________________________________________________________________________________________________________
% 
% BATCH structure fields:
% 
%  filename          : conn_*.mat project file (defaults to currently open project)
%  subjects          : Subset of subjects to run processing steps or define parameters for (defaults to all subjects)
%  parallel          : Parallelization options (defaults to local procesing / no parallelization)
%  Setup             : Information/processes regarding experiment Setup and Preprocessing
%  Denoising         : Information/processes regarding Denoising step
%  Analysis          : Information/processes regarding first-level analyses
%  Results           : Information/processes regarding second-level analyses/results
%  QA                : Information/processes regarding Quality Assurance plots
% 
% 
% BATCH.parallel DEFINES PARALLELIZATION OPTIONS (applies to any Setup/Setup.preprocessing/Denoising/Analysis/QA steps) %!
% 
%    parallel.N                         : number of parallel jobs; 0 to run locally ([0])
%    parallel.profile                   : (optional) name of parallelization profile 
%                                           if undefined CONN uses the default parallelization profile defined in GUI.Tools.GridSettings
%                                           see "conn_jobmanager profiles" for a list of all available profiles
%                                           see GUI Tools.GridSettings for additional information and to add/edit profiles
%                                           use the profile name 'Null profile' to queue this job (queued/scripted jobs are prepared but
%                                            not submitted; see GUI.Tools.SeePendingJobs to submit next queued job)
%    parallel.cmd_submitoptions         : (optional) alternative value for parallelization profile 'in-line' additional submit-settings 
%                                           defaults to the chosen parallelization profile value for this field
%    parallel.cmd_submitoptions_infile   : (optional) alternative value for parallelization profile 'in-file' additional submit-settings
%                                           defaults to the chosen parallelization profile value for this field
%    parallel.cmd_rundeployed            : (optional) aternative value for profile 'nodes use pre-compiled CONN only' setting
%                                           defaults to the chosen parallelization profile value for this field
%    parallel.cmd_checkstatus_automatic  : (optional) aternative value for profile 'check jobs status automatically' setting
%                                           defaults to the chosen parallelization profile value for this field
%    parallel.immediatereturn            : (optional) 1/0 : 1 returns control to Matlab without waiting for parallel job to finish ([0])
% 
% BATCH.Setup DEFINES EXPERIMENT SETUP AND PERFORMS INITIAL DATA EXTRACTION AND/OR PREPROCESSING STEPS %!
% 
%    Setup.isnew                        : 1/0 set to 1 when defining a new conn project in order to initialize it, disregarding any
%                                         previous information stored in this project [0]
%    Setup.done                         : 1/0: 0 defines fields only; 1 defines fields and runs SETUP processing steps afterwards [0]
%    Setup.overwrite                    : (for done=1) 1/0 overwrites target files if they exist [1]
% 
%    Setup.nsubjects                    : Number of subjects
%    Setup.RT                           : Repetition time (seconds) [NaN: read from sidecar json file if available]
%    Setup.acquisitiontype              : 1/0: Continuous acquisition of functional volumes [1] 
% 
%    Setup.functionals                  : functionals{nsub}{nses} char array of functional volume files (dataset-0; for voxel-level analyses; 
%                                         see secondarydatasets below) 
%    Setup.structurals                  : structurals{nsub} char array of structural volume files 
%                                         OR structurals{nsub}{nses} char array of anatomical session-specific volume files 
%    Setup.secondarydatasets            : define one or several additional functional datasets (e.g. fmap or vdm files for susceptibility 
%                                         distorition, alternative functional files for ROI-level timeseries extraction, etc.) 
%                                         [default secondarydatasets equal to struct('functionals_type',2)]
%       Setup.secondarydatasets.functionals_label : label of secondary functional dataset (e.g. 'fmap', 'vdm', ...)
%       Setup.secondarydatasets.functionals_explicit : functionals_explicit{nsub}{nses} char array of volume files (only for functionals_type==4 only) 
%       Setup.secondarydatasets.functionals_type : alternative option for defining files in secondary functional dataset (1-4) [2]: 
%                                          1: same files as functional data
%                                          2: same files as functional data field after removing leading 's' from filename
%                                          3: other (same as functional data field but using alternative filename-change rule; 
%                                             see functionals_rule above and help conn_rulebasedfilename); 
%                                          4: other (explicitly specify the functional volume files using functionals_explicit field above)
%       Setup.secondarydatasets.functionals_rule : (only for functionals_type==3 only) regexprep(filename,functionals_rule{2},functionals_rule{3})  
%                                         converts filenames in 'Setup.functionals' field to filenames that will be used when extracting
%                                          BOLD signal ROI timeseries (if functionals_rule{1}==2 filename is interpreted as a full path;  
%                                         if functionals_rule{1}==1 filename is interpreted as only the file *name* -no file path, no 
%                                         file extension-)    
%    Setup.add                          : 1/0; use 0 (default) to define the full set of subjects in your experiment; use 1 to define an 
%                                         additional set of subjects (to be added to any already-existing subjects in your project) [0]
%                                         When using Setup.add=1, the following fields are expected to contain the information for the new
%                                          /added subjects *only*: Setup.functionals, Setup.structurals, Setup.functionals_explicit, 
%                                          Setup.vdm_functionals, Setup.fmap_functionals, Setup.coregsource_functionals, Setup.spmfiles, 
%                                          Setup.masks.Grey/White/CSF, Setup.rois.files, Setup.conditions.onsets/durations, Setup.covariates.files
%                                          Setup.subjects.effects, Setup.subjects.groups
%                                         When using Setup.add=1 in combination with Setup.done, Setup.preprocessing, Denoising.done, and/or 
%                                          Analysis.done only the new/added subjects will be processed
%                                         When using Setup.add=1 the BATCH.subjects field is disregarded/overwritten to point to the new/added 
%                                          subjects only
%                                         note: Setup.add cannot be used in combination with any of the Setup.rois.add, Setup.conditions.add, or 
%                                          Setup.covariates.add options within the same batch structure
% 
%    Setup.masks                        : defines gray/white/CSF masks for each subject
%      Setup.masks.Grey                 : masks.Grey{nsub} char array of grey matter mask volume file [defaults to Grey mask from structural] 
%      Setup.masks.White                : masks.White{nsub} char array of white matter mask volume file [defaults to White mask from structural] 
%      Setup.masks.CSF                  : masks.CSF{nsub} char array of CSF mask volume file [defaults to CSF mask from structural] 
%                                       : each of these fields can also be defined as a double cell array for session-specific files (e.g. 
%                                          mask.Grey{nsub}{nses} grey matter file for subject nsub and session nses)
%                                       : each of these fields can also be defined as a structure with fields files/dimensions/etc. 
%                                         (same as 'Setup.rois' below).
%    Setup.rois                         : defines arbitrary ROIs
%      Setup.rois.names                 : rois.names{nroi} char array of ROI name [defaults to ROI filename]
%      Setup.rois.files                 : rois.files{nroi}{nsub}{nses} char array of roi file (rois.files{nroi}{nsub} char array of roi file, 
%                                         to use the same roi for all sessions; or rois.files{nroi} char array of roi file, to use the same 
%                                         roi for all subjects)
%      Setup.rois.dimensions            : rois.dimensions{nroi} number of ROI dimensions - # temporal components to extract from ROI [1] (set 
%                                         to 1 to extract the average timeseries within ROI voxels; set to a number greater than 1 to extract 
%                                         additional PCA timeseries within ROI voxels 
%      Setup.rois.weighted              : rois.weighted(nroi) 1/0 to use weighted average/PCA computation when extracting temporal components 
%                                         from each ROI (BOLD signals are weighted by the ROI mask value at each voxel)
%      Setup.rois.multiplelabels        : rois.multiplelabels(nroi) 1/0 to indicate roi file contains multiple labels/ROIs (default: set to 
%                                         1 if there exist an associated .txt or .xls file with the same filename and in the same folder as 
%                                         the roi file)
%      Setup.rois.labelfiles            : rois.labelfiles{nroi} to import label file from alternative location
%      Setup.rois.mask                  : rois.mask(nroi) 1/0 to mask with grey matter voxels [0] 
%      Setup.rois.regresscovariates     : rois.regresscovariates(nroi) 1/0 to regress known first-level covariates before computing PCA 
%                                         decomposition of BOLD signal within ROI [1 if dimensions>1; 0 otherwise] 
%      Setup.rois.dataset               : rois.dataset(nroi) index n to Secondary Dataset #n identifying the version of functional data 
%                                         coregistered  to this ROI to extract BOLD timeseries from [1] (set to 0 to extract BOLD signal 
%                                         from functional data instead; secondary datasets may be identified by their index or by their 
%                                         label -see 'functional_label' preprocessing step)
%      Setup.rois.add                   : 1/0; use 0 (default) to define the full set of ROIs to be used in your analyses; use 1 to define 
%                                         an additional set of ROIs (to be added to any already-existing ROIs in your project) [0]
%  
%    Setup.conditions                   : defines analysis conditions (all functional connectivity measures will be computed separately for
%                                         each defined condition)
%      Setup.conditions.names           : conditions.names{ncondition} char array condition name
%      Setup.conditions.onsets          : conditions.onsets{ncondition}{nsub}{nses} vector of condition onsets (in seconds)
%      Setup.conditions.durations       : conditions.durations{ncondition}{nsub}{nses} vector of condition durations (in seconds)
%      Setup.conditions.param           : conditions.param(ncondition) temporal modulation (0 for no temporal modulation; positive index to 
%                                         first-level covariate for other temporal interactions) 
%      Setup.conditions.filter          : conditions.filter{ncondition} temporal/frequency decomposition ([] for no decomposition; [low high] 
%                                         for fixed band-pass frequency filter; [N] for filter bank decompositoin with N frequency filters; 
%                                         [Duration Onsets] in seconds for sliding-window decomposition where Duration is a scalar and Onsets 
%                                         is a vector of two or more sliding-window onset values) 
%      Setup.conditions.missingdata     : 1/0 Allow subjects with missing condition data (empty onset/duration fields in *all* of the
%                                         sessions) [0] 
%      Setup.conditions.model           : (optional) conditions.model{ncondition} cell array characterizing definition of secondary conditions, 
%                                         as a function of other regular conditions. Format {@fun, condition_name1, condition_name2, ...}
%                                         where @fun is a function handle characterizing a function fun(x) with x a matrix of functional
%                                         connectivity values for one or multiple conditions (rows) and one or multiple datapoints (columns)
%                                         Example {@(x)mean(x,1),'Time1','Time2','Time3'} will compute the mean of the three Time# conditions
%                                         Alternative forms of conditions.model{ncondition}{1} values (functional definition)
%                                          'avg','std','min','max'   : to compute average/standard-deviation/minimum/maximum across selected conditions
%                                          'lin',G       : to compute regressor coefficients in linear model X=G*B fitting all selected conditions
%                                             where G is a [conditions x regressors] design matrix (note: for subject-specific values define
%                                             G as a char array using Matlab notation and enter second-level covariate names as values). When
%                                             G contains multiple columns the first column defines the effect of interest (the rest are
%                                             estimated but simply disregarded; not taken into second-level analyses)
%                                          @fun    : (single argument function) to compute an arbitrary function fun(X) of condition 
%                                             values X (conditions x datapoints matrix) for each subject
%                                          @fun    : (two or three arguments function for subject-specific functions) to compute an 
%                                             arbitrary function fun(X,c,cnames) as a function of not only the condition values X but also
%                                             of second-level covariate values C (1 x covariates vector) for all second-level covariates 
%                                             in your project (covariate names will be listed as the third argument of fun for reference)
%      Setup.conditions.importfile      : (optional) Alternatively, importfile is a char or cell array pointing to a '*.txt','*.csv', or 
%                                         BIDS- '*.tsv' file containing conditions names/onsets/durations information (see help 
%                                         conn_importcondition)
%      Setup.conditions.importfile_options : (for conditions.importfile procedure only) cell array containing additional options to pass 
%                                            to conn_importcondition when importing condition info (see help conn_importcondition)
%      Setup.conditions.add             : 1/0; use 0 (default) to define the full set of conditions to be used in your analyses; use 1 to 
%                                         define an additional set of conditions (to be added to any already-existing conditions in your 
%                                         project) [0]
%
%    Setup.covariates                   : defines 1st-level covariates (timeseries for each subject/session; these may be later used as potential
%                                         factors during denoising, as seed timeseries for 1st-level analyses, as modulatory terms, etc.)
%      Setup.covariates.names           : covariates.names{ncovariate} char array of first-level covariate name
%      Setup.covariates.files           : covariates.files{ncovariate}{nsub}{nses} char array of covariate file 
%      Setup.covariates.add             : 1/0; use 0 (default) to define the full set of covariates to be used in your analyses; use 1 to 
%                                         define an additional set of covariates (to be added to any already-existing covariates in your 
%                                         project) [0]
%  
%    Setup.subjects                     : defines 2nd-level covariates (arbitrary continuous/categorical/ordinal data for each subject)
%      Setup.subjects.effects           : subjects.effects{neffect} vector of size [nsubjects,1] defining second-level effects
%      Setup.subjects.effect_names      : subjects.effect_names{neffect} char array of second-level covariate name
%      Setup.subjects.effect_descrip    : (optional) subjects.effect_descrip{neffect} char array of effect description (long name; for  
%                                         display purposes only)
%      Setup.subjects.add               : 1/0; use 0 (default) to define the full set of covariates to be used in your analyses; use 1 to 
%                                         define an additional set of covariates (to be added to any already-existing covariates in your 
%                                         project) [0]
%  
%    Setup.subjects                     : alternative way of defining 2nd-level covariates (subject groups only)
%      Setup.subjects.groups            : subjects.group vector of size [nsubjects,1] (with values from 1 to ngroup) defining subject groups
%      Setup.subjects.group_names       : subjects.group_names{ngroup} char array of second-level group name
%      Setup.subjects.group_descrip     : (optional) subjects.group_descrip{neffect} char array of group description (long name; for display 
%                                         purposes only)
%      Setup.subjects.add               : 1/0; use 0 (default) to define the full set of covariates to be used in your analyses; use 1 to 
%                                         define an additional set of covariates (to be added to any already-existing covariates in your 
%                                         project) [0]
%  
%    Setup.analyses                     : Vector of index to analysis types (1: ROI-to-ROI; 2: Seed-to-voxel; 3: Voxel-to-voxel; 4: Dynamic 
%                                         FC) Defaults to vector [1,2,3,4] (all analyses)
%    Setup.voxelmask                    : Analysis mask (voxel-level analyses): 1: Explicit mask (brainmask.nii); 2: Implicit mask 
%                                         (subject-specific) [1] 
%    Setup.voxelmaskfile                : Explicit mask file (only when voxelmask=1) [fullfile(fileparts(which('spm')),'apriori',
%                                        'brainmask.nii')] 
%    Setup.voxelresolution              : Analysis space (voxel-level analyses): 1: Volume-based template (SPM; default 2mm isotropic 
%                                         or same as explicit mask if specified); 2: Same as structurals; 3: Same as functionals; 
%                                         4: Surface-based template (Freesurfer) [1] 
%    Setup.analysisunits                : BOLD signal units: 1: PSC units (percent signal change); 2: raw units [1] 
%    Setup.outputfiles                  : Optional output files (outputfiles(1): 1/0 creates confound beta-maps; outputfiles(2): 1/0 creates 
%                                         confound-corrected timeseries; outputfiles(3): 1/0 creates seed-to-voxel r-maps) ;outputfiles(4): 
%                                         1/0 creates seed-to-voxel p-maps) ;outputfiles(5): 1/0 creates seed-to-voxel FDR-p-maps); 
%                                         outputfiles(6): 1/0 creates ROI-extraction REX files; [0,0,0,0,0,0] 
%    Setup.spmfiles                     : Optionally, spmfiles{nsub} is a char array pointing to the 'SPM.mat' source file to extract Setup 
%                                         information from for each subject (use alternatively spmfiles{nsub}{nses} for session-specific 
%                                         SPM.mat files) 
%    Setup.spmfiles_options             : (for Setup.spmfiles procedure) Cell array containing additional options to pass to conn_importspm 
%                                          when importing experiment info from spmfiles (see help conn_importspm)
%    Setup.vdm_functionals              : (for Setup.preprocessing.steps=='realign&unwarp&fieldmap') vdm_functionals{nsub}{nses} char 
%                                         array of voxel-displacement volumes (vdm* file; explicitly entering these volumes here superceeds CONN's 
%                                         default option to search for/use vdm* files in same directory as functional data) 
%    Setup.fmap_functionals             : (for Setup.preprocessing.steps=='vdm_create') fmap_functionals{nsub}{nses} char 
%                                         array of fieldmap sequence files (magnitude1+phasediff or real1+imag1+real2+imag2 or fieldmap (Hz) volumes)
%    Setup.coregsource_functionals      : (for Setup.preprocessing.steps=='functional_coregister/segment/normalize') 
%                                         coregsource_functionals{nsub} char array of source volume for coregistration/normalization/
%                                         segmentation (used only when preprocessing "coregtomean" field is set to 2, user-defined source 
%                                         volumes are used in this case instead of either the first functional volume (coregtomean=0) or the 
%                                         mean functional volume (coregtomean=1) for coregistration/normalization/segmentation) 
%    Setup.localcopy                    : (for Setup.structural, Setup.functional, Setup.secondarydatasets, and Setup.rois) 1/0 : copies structural/
%                                         functional files into conn_*/data/BIDS folder before importing into CONN [0]
%    Setup.binary_threshold             : (for BOLD extraction from Grey/White/CSF ROIs) Threshold value # for binarizing Grey/White/CSF 
%                                         masks [.5 .5 .5] 
%    Setup.binary_threshold_type        : (for BOLD extraction from Grey/White/CSF ROIs) 1: absolute threshold (keep voxels with values
%                                         above x); 2: percentile threshold (keep x% of voxels with the highest values) [1 1 1] 
%    Setup.exclude_grey_matter          : (for BOLD extration from White/CSF ROIs) threhsold for excluding Grey matter voxels (nan for no 
%                                         threshold) [nan nan nan]
%    Setup.erosion_steps                : (for BOLD extraction from Grey/White/CSF ROIs) integer numbers are interpreted as erosion kernel 
%                                         size for Grey/White/CSF mask erosion after binarization; non-integer numbers are interpreted as 
%                                         percentile voxels kept after erosion [0 1 1]
%    Setup.erosion_neighb               : (for BOLD extraction from Grey/White/CSF ROIs; only when using integer erosion_steps/ kernel sizes, 
%                                         this field is disregarded otherwise) Neighborhood size for Grey/White/CSF mask erosion after
%                                         binarization (a voxel is eroded if there are more than masks_erosion_neighb zeros within the 
%                                         (2*masks_erosionsteps+1)^3-neighborhood of each voxel) [1 1 1]
% 
%  
% BATCH.Setup.preprocessing PERFORMS DATA PREPROCESSING STEPS (realignment/slicetiming/coregistration/segmentation/normalization/smoothing) %!
%
%  Setup.preprocessing.steps         : List of data preprocessing steps (cell array containing a subset of the following step names, in 
%                                     the desired order; e.g. {'functional_realign','functional_art'}):
%                      PIPELINES:
%                        'default_mni'                           : default MNI-space preprocessing pipeline
%                        'default_mnifield'                      : same as default_mni but with vdm/fieldmap information and indirect normalization
%                        'default_mnidirectfield'                : same as default_mni but with vdm/fieldmap information and direct normalization
%                        'default_ss'                            : default subject-space preprocessing pipeline
%                        'default_ssfield'                       : same as default_ss but with vdm/fieldmap information
%                        'default_ssnl'                          : same as default_ss but with non-linear coregistration
%                      INDIVIDUAL STRUCTURAL STEPS:
%                        'structural_center'                     : centers structural data to origin (0,0,0) coordinates
%                        'structural_manualorient'               : applies user-defined affine transformation to structural data
%                        'structural_manualspatialdef'           : applies user-defined spatial deformation to structural data
%                        'structural_mask'                       : masks structural data using inclusive or exclusive mask
%                        'structural_segment&normalize'          : structural unified normalization and segmentation 
%                        'structural_segment&normalize&lesion'   : structural unified normalization and segmentation with lesion mask
%                                                                   (normalizes structural data and creates a modified TPM for functional
%                                                                   normalization that includes the lesion as an added tissue class)
%                        'structural_normalize'                  : structural normalization to MNI space (without segmentation)
%                        'structural_normalize_preservemasks'    : structural normalization to MNI space with user-defined Grey/White/CSF masks 
%                                                                   (normalizes structural data and applies same transformation to user-defined 
%                                                                   Grey/White/CSF mask ROIs)
%                        'structural_segment'                    : structural segmentation (Grey/White/CSF tissue classes)
%                      INDIVIDUAL FUNCTIONAL (or combined functional/structural) STEPS:
%                        'functional_art'                        : functional identification of outlier scans (from motion 
%                                                                   displacement and global signal changes)
%                        'functional_bandpass'                   : functional band-pass filtering
%                        'functional_center'                     : centers functional data to origin (0,0,0) coordinates
%                        'functional_coregister_affine'          : functional affine coregistration to structural volumes
%                        'functional_coregister_nonlinear'       : functional non-linear coregistration to structural volumes
%                        'functional_label'                      : labels current functional files (to list of Secondary Datasets)
%                        'functional_load'                       : assigns current functional files (from list of Secondary Datasets)
%                        'functional_manualorient'               : applies user-defined affine transformation to functional data
%                        'functional_manualspatialdef'           : applies user-defined spatial deformation to functional data
%                        'functional_mask'                       : masks functional data using inclusive or exclusive mask
%                        'functional_motionmask'                 : creates functional motion masks (mean BOLD signal spatial 
%                                                                   derivatives wrt motion parameters)
%                        'functional_normalize_direct'           : functional direct normalization
%                        'functional_normalize_indirect'         : functional indirect normalization (coregister to structural; 
%                                                                   normalize structural; apply same transform to functionals)
%                        'functional_normalize_indirect_preservemasks': functional indirect normalization with user-defined Grey/White/CSF masks 
%                                                                   (coregister to structural; normalize structural; apply same transformation to 
%                                                                   functionals as well as to Grey/White/CSF masks)
%                        'functional_realign'                    : functional realignment
%                        'functional_realign_noreslice'          : functional realignment without reslicing 
%                                                                   (applies transform to source header files)
%                        'functional_realign&unwarp'             : functional_realignment + unwarp 
%                                                                   (removes motion-by-inhomogeneity interactions)
%                        'functional_realign&unwarp&fieldmap'    : functional_realignemnt&unwarp + distortion correction 
%                                                                   (corrects static inhomogeneity distortions)
%                        'functional_regression'                 : removal of user-defined temporal components from BOLD timeseries (keeps 
%                                                                   residuals of linear regression model)
%                        'functional_removescans'                : removes user-defined number of initial scans from functional
%                        'functional_roiextract'                 : extraction of ROI timeseries (compute BOLD timeseres within ROI)
%                        'functional_segment'                    : functional segmentation (Grey/White/CSF tissue classes)
%                        'functional_segment&normalize_direct'   : functional direct unified normalization and segmentation
%                        'functional_segment&normalize_indirect' : functional indirect unified normalization and segmentation
%                                                                   (coregister to structural; normalize and segment structural; 
%                                                                   apply same transformation to functionals)
%                        'functional_slicetime'                  : functional slice-timing correction
%                        'functional_smooth'                     : functional spatial smoothing
%                        'functional_smooth_masked'              : functional spatial masked-smoothing (spatial convolution with Gaussian kernel 
%                                                                   restricted to voxels within custom functional mask
%                        'functional_vdm_create'                 : creation of vdm (voxel-displacement-map) from fieldmap dataset (reads 'fmap' 
%                                                                   secondary functional dataset containing magnitude and phasediff images and 
%                                                                   creates 'vdm' secondary functional dataset containing voxel-displacement map)
%                       SURFACE FUNCTIONAL STEPS:
%                        'functional_surface_resample'           : resample functional data at the location of FreeSurfer subject-specific 
%                                                                   structural cortical surface
%                        'functional_surface_smooth'             : functional spatial diffusion of surface data 
%                        'functional_surface_coreg&resample'     : coregister&resample functional data at the location of FreeSurfer subject-specific 
%                                                                   structural cortical surface
%                      If steps is left empty or unset a gui will prompt the user to specify the desired preprocessing pipeline 
%                      If steps points to an existing preprocessing-pipeline file (e.g. saved from GUI) the corresponding 
%                       preprocessing-pipeline will be run
%   
%      Setup.preprocessing.affreg          : (normalization) affine registration before normalization ['mni']
%      Setup.preprocessing.art_thresholds  : (functional_art) ART thresholds for identifying outlier scans 
%                                            art_thresholds(1): threshold value for global-signal (z-value; default 5) 
%                                            art_thresholds(2): threshold value for subject-motion (mm; default .9) 
%                                           additional options: 
%                                            art_thresholds(3): 1/0 global-signal threshold based on scan-to-scan changes
%                                                               in global-BOLD measure (default 1) 
%                                            art_thresholds(4): 1/0 subject-motion threshold based on scan-to-scan changes 
%                                                               in subject-motion measure (default 1) 
%                                            art_thresholds(5): 1/0 subject-motion threhsold based on composite-movement 
%                                                               measure (default 1) 
%                                            art_thresholds(6): 1/0 force interactive mode (ART gui) (default 0) 
%                                            art_thresholds(7): [only when art_threshold(5)=0] subject-motion threshold 
%                                                               based on rotation measure 
%                                            art_thresholds(8): N number of initial scans to be flagged for removal 
%                                                               (default 0)
%                                          note: when art_threshold(5)=0, art_threshold(2) defines the threshold based on the translation 
%                                           measure, and art_threhsold(7) defines the threshold based on the rotation measure; otherwise 
%                                           art_threshold(2) defines the (single) threshold based on the composite-motion measure 
%                                          note: the default art_thresholds(1:2) [5 .9] values correspond to the "intermediate" 
%                                           (97th percentile) settings; to use the "conservative" (95th percentile) settings use 
%                                           [3 .5]; to use the "liberal" (99th percentile) settings use [9 2] values instead
%                                          note: art needs subject-motion files to estimate possible outliers. If a 'realignment' 
%                                           first-level covariate exists it will load the subject-motion parameters from that first-
%                                           level covariate; otherwise it will look for a rp_*.txt file (SPM format) in the same 
%                                           folder as the functional data
%                                          note: subject-motion files can be in any of the following formats: a) *.txt file (SPM 
%                                           format; three translation parameters in mm followed by pitch/roll/yaw in radians); 
%                                           b) *.par (FSL format; three Euler angles in radians followed by translation parameters 
%                                           in mm); c) *.siemens.txt (Siemens MotionDetectionParameter.txt format); d) *.deg.txt (same 
%                                           as SPM format but rotations in degrees instead of radians)
%      Setup.preprocessing.boundingbox     : (normalization) target bounding box for resliced volumes (mm) [-90,-126,-72;90,90,108] 
%      Setup.preprocessing.bp_filter       : (functional_bandpass, functional_regression) Low- and High- frequency thresholds (in Hz)
%      Setup.preprocessing.bp_keep0        : (functional_bandpass) 0: removes average BOLD signal (freq=0Hz component); 1: keeps average BOLD signal 
%                                            in output independent of band-pass filter values; [1]
%      Setup.preprocessing.coregtomean     : (functional_coregister/segment/normalize) 0: use first volume; 1: use mean volume (computed during 
%                                            realignment); 2: use user-defined source volume (see Setup.coregsource_functionals field) [1]
%      Setup.preprocessing.diffusionsteps  : (surface_smooth) number of diffusion steps
%      Setup.preprocessing.fwhm            : (functional_smooth) Smoothing factor (mm) [8]
%      Setup.preprocessing.interp          : (normalization) target voxel interpolation method (0:nearest neighbor; 1:trilinear; 2 or higher:n-order 
%                                            spline) [4]
%      Setup.preprocessing.label           : (functional_label) label of secondary dataset (note: the following functional step names do not require an 
%                                            explicit label field: 'functional_label_as_original', 'functional_label_as_subjectspace', 
%                                            'functional_label_as_mnispace', 'functional_label_as_surfacespace', 'functional_label_as_smoothed')
%      Setup.preprocessing.load_label      : (functional_load) label of secondary dataset (note: the following functional step names do not require an 
%                                            explicit continue field: 'functional_load_from_original', 'functional_load_from_subjectspace', 
%                                            'functional_load_from_mnispace', 'functional_load_from_surfacespace', 'functional_load_from_smoothed')
%      Setup.preprocessing.mask_names_anat : (strucutral_mask) list of ROI names (if multiple ROIs, the intersection of all ROIs will be used as mask)
%      Setup.preprocessing.mask_names_func : (functional_mask, functional_smooth_masked) list of ROI names (if multiple ROIs, the intersection of all ROIs will be used as mask)
%      Setup.preprocessing.mask_inclusive_anat : (strucutral_mask) 1: inclusive ROI mask (keep voxels inside ROI); 0: exclusive ROI mask (keep voxels outside ROI) [1]
%      Setup.preprocessing.mask_inclusive_func : (functional_mask, functional_smooth_masked) 1: inclusive ROI mask (keep voxels inside ROI); 0: exclusive ROI mask (keep voxels outside ROI) [1]
%      Setup.preprocessing.reg_names       : (functional_regression) list of first-level covariates to use as model regressors / design matrix (valid 
%                                            entries are first-level covariate names or ROI names)
%      Setup.preprocessing.reg_dimensions  : (functional_regression) list of maximum number of dimensions (one value for each model regressor in 
%                                            reg_names)
%      Setup.preprocessing.reg_deriv       : (functional_regression) list of 0/1/2 values (one value for each model regressor in reg_names): add 
%                                            first- or second- order derivatives to each model regressor
%      Setup.preprocessing.reg_filter      : (functional_regression) list of 0/1 values (one value for each model regressor in reg_names):  
%                                            band-pass filter individual model regressors (filter specified in bp_filter field)
%      Setup.preprocessing.reg_lag         : (functional_regression) list of 0/1 values (one value for each model regressor in reg_names):
%                                            removes lagged regressor (estimating optimal lag between each regressor and the BOLD signal at each voxel)
%      Setup.preprocessing.reg_lagmax      : (functional_regression) maximum lag (in seconds); when estimating optimal lag (see reg_lag) lags
%                                            between -reg_lagmax and +reg_lagmax are considered
%      Setup.preprocessing.reg_detrend     : (functional_regression) 1: adds a linear/detrending term to model regressors [1]
%      Setup.preprocessing.reg_skip        : (functional_regression) 1: does not create output functional files, only creates session-specific 
%                                            dp_*.txt files with covariate timeseries to be included later in an arbitrary first-level model [0]
%      Setup.preprocessing.removescans     : (functional_removescans) number of initial scans to remove
%      Setup.preprocessing.reorient        : (functional/structural_manualorient) 3x3 or 4x4 transformation matrix or filename containing
%                                            corresponding matrix
%      Setup.preprocessing.respatialdef    : (functional/structural_manualspatialdef) nifti deformation file (e.g. y_*.nii or *seg_sn.mat files)
%      Setup.preprocessing.roi_names       : (functional_roiextract) list of ROI names
%                                              additional 1st-level covariate names may be included (to be regressed-out from ROI timeseries)
%      Setup.preprocessing.roi_dimensions  : (functional_roiextract) list of maximum number of dimensions (one value for each entry in roi_names)
%      Setup.preprocessing.roi_deriv       : (functional_roiextract) list of 0/1/2 values (one value for each entry in roi_names): adds 
%                                            first- or second- order derivatives to each extracted timeseries
%      Setup.preprocessing.roi_filter      : (functional_roiextract) list of 0/1 values (one value for each entry in roi_names):  
%                                            band-pass filter individual timeseries (filter specified in bp_filter field)
%      Setup.preprocessing.roi_detrend     : (functional_roiextract) 1: detrends extracted BOLD timeseries [0]
%      Setup.preprocessing.roi_scale       : (functional_roiextract) 1: scales extracted BOLD timeseries to PSC units (within each ROI) [1]
%      Setup.preprocessing.rtm             : (functional_realign) 0: use first volume; 1: use mean volume [0]
%      Setup.preprocessing.rmask           : (functional_realign) 1: applies implicit masking (voxels outside of field of view in >=1 image are set to NaN); [1]
%      Setup.preprocessing.sliceorder      : (functional_slicetime) acquisition order (vector of indexes; 1=first slice in image; note: use cell
%                                             array for subject-specific vectors)
%                                            alternatively sliceorder may also be defined as one of the following strings: 'ascending',
%                                             'descending','interleaved (middle-top)','interleaved (bottom-up)','interleaved (top-down)',
%                                             'interleaved (Siemens)','interleaved (Philips)','BIDS' (this option reads slice timing information 
%                                             from .json files)
%                                            alternatively sliceorder may also be defined as a vector containing the acquisition time in 
%                                             milliseconds for each slice (e.g. for multi-band sequences) 
%      Setup.preprocessing.ta              : (functional_slicetime) acquisition time (TA) in seconds (used to determine slice times when 
%                                             sliceorder is defined by a vector of slice indexes; note: use vector for subject-specific 
%                                             values). Defaults to (1-1/nslices)*TR where nslices is the number of slices
%      Setup.preprocessing.template_structural: (structural_normalize SPM8 only) anatomical template file for approximate coregistration 
%                                             [spm/template/T1.nii]
%      Setup.preprocessing.template_functional: (functional_normalize SPM8 only) functional template file for normalization 
%                                             [spm/template/EPI.nii]
%      Setup.preprocessing.tpm_template    : (any segment/normalize option in SPM12) tissue probability map [spm/tpm/TPM.nii]
%                                            alternatively, location of subject-specific TPM files (secondary functional dataset number or name ['tpm'])
%      Setup.preprocessing.tpm_ngaus       : (structural_segment, structural_segment&normalize in SPM8&SPM12) number of gaussians for each 
%                                             tissue probability map
%      Setup.preprocessing.tpm_structlesion: (structural_segment&normalize&lesion) name of ROI containing a structural-lesion mask
%                                             (the lesion mask is expected to be coregistered with the structural, as part of structural normalization 
%                                              a new TPM template will be created with the lesion as an added tissue class)
%      Setup.preprocessing.vdm_et1         : (functional_vdm_create) ET1 (Echo Time first echo in fieldmap sequence) 
%                                             (default [] : read from .json file / BIDS)
%      Setup.preprocessing.vdm_et2         : (functional_vdm_create) ET2 (Echo Time second echo in fieldmap sequence)
%                                             (default [] : read from .json file / BIDS)
%      Setup.preprocessing.vdm_ert         : (functional_vdm_create) ERT (Effective Readout Time in funcional data)
%                                             (default [] : read from .json file / BIDS)
%      Setup.preprocessing.vdm_blip        : (functional_vdm_create) k-space traversal blip direction along the y-axis following SPM convention 
%                                             (i.e. positive for P>A, negative for A>P)
%                                            use +1 or -1 to specify this value explicitly
%                                            leave empty to read from .json file /BIDS PhaseEncodingDirection field and use the formula 
%                                               BLIP = sign([0 1 0 0]*vol.mat*[i j k 0]'; 
%                                               e.g. PhaseEncodingDirection='j+', mat=[-1 0 0;0 1 0;0 0 1] => BLIP=sign([0 1 0 0]*mat*[0 1 0 0])=+1
%                                            use 0 to reverse the sign from the above formula
%      Setup.preprocessing.vdm_type        : (functional_vdm_create only) type of fieldmap sequence files 
%                                               []  : automatically detect
%                                               1   : magnitude+phasediff (or magnitude1+magnitude2+phasediff)
%                                               2   : real1+imag1+real2+imag2
%                                               3   : fieldmapHz (e.g. --fout output from TOPUP)
%      Setup.preprocessing.vdm_fmap        : (functional_vdm_create only) location of fieldmap sequence files (secondary functional dataset number
%                                             or label containing fieldmap sequence files) ['fmap']
%      Setup.preprocessing.voxelsize_anat  : (structural normalization) target voxel size for resliced volumes (mm) [1]
%      Setup.preprocessing.voxelsize_func  : (functional normalization) target voxel size for resliced volumes (mm) [2]
%      Setup.preprocessing.sessions        : defines functional sessions to preprocess [1:max # of sessions]
%      Setup.preprocessing.sets            : defines functional dataset to preprocess (0 for functional data; [1-N] or labels for Secondary 
%                                             Datasets) [0]
%  
%  
% BATCH.Denoising PERFORMS DENOISING STEPS (confound removal & filtering) %!
% 
%    Denoising.done                 : 1/0: 0 defines fields only; 1 runs DENOISING processing steps [0]
%    Denoising.overwrite            : (for done=1) 1/0: overwrites target files if they exist [1]
%    Denoising.filter               : vector with two elements specifying band pass filter: low-frequency & high-frequency cutoffs (Hz)
%    Denoising.detrending           : 0/1/2/3: BOLD times-series polynomial detrending order (0: no detrending; 1: linear detrending; 
%                                     ... 3: cubic detrending) 
%    Denoising.despiking            : 0/1/2: temporal despiking with a hyperbolic tangent squashing function (1:before regression; 
%                                     2:after regression) [0] 
%    Denoising.regbp                : 1/2: order of band-pass filtering step (1 = RegBP: regression followed by band-pass; 2 = Simult: 
%                                     simultaneous regression&band-pass) [1] 
%    Denoising.confounds            : Cell array of confound names (alternatively see 'confounds.names' below)
% 
%    Denoising.confounds            : alternatively confounds can be a structure with fields
%      Denoising.confounds.names    : confounds.names{nconfound} char array of confound name (confound names can be: 'Grey Matter',
%                                     'White Matter','CSF',any ROI name, any covariate name, or 'Effect of *' where * represents 
%                                      any condition name])
%                                     note: use confounds.names={} to specify CONN's default set of confounds, or confounds.names={''}
%                                      to indicate no confounds at all
%      Denoising.confounds.dimensions : confounds.dimensions{nconfound} number of confound dimensions [defaults to using all dimensions 
%                                     available for each confound variable]
%      Denoising.confounds.deriv    : confounds.deriv{nconfound} include temporal derivatives up to n-th order of each effect (0 for 
%                                     raw timeseries, 1 for raw+firstderivative timeseries, etc.) [0|1]
%      Denoising.confounds.power    : confounds.power{nconfound} include powers up to n-th order of each effect (1 for linear effects, 
%                                     2 for linear+quadratic effect, etc.) [1]
%      Denoising.confounds.filter   : (for regbp==1) confounds.filter{nconfound} band-pass filter confound regressors before entering 
%                                     in regression equation [0]
%  
%  
% BATCH.Analysis PERFORMS FIRST-LEVEL ANALYSES (ROI-to-ROI and seed-to-voxel) %!
% 
%    Analysis.done                  : 1/0: 0 defines fields only; 1 runs ANALYSIS processing steps [0]
%    Analysis.overwrite             : (for done=1) 1/0: overwrites target files if they exist [1]
%    Analysis.name                  : analysis name (identifying each set of independent analysis)
%                                     (alternatively sequential index identifying each set of independent analyses [1])
%                      
%    Analysis.measure               : connectivity measure used, 1 = 'correlation (bivariate)', 2 = 'correlation (semipartial)', 3 = 
%                                     'regression (bivariate)', 4 = 'regression (multivariate)'; [1] 
%    Analysis.weight                : within-condition weight, 1 = 'none', 2 = 'hrf', 3 = 'hanning'; [2] 
%    Analysis.modulation            : temporal modulation, 0 = standard weighted GLM analyses; 1 = gPPI analyses of condition-specific 
%                                     temporal modulation factor, or a string for PPI analyses of other temporal modulation factor 
%                                     (same for all conditions; valid strings are ROI names and 1st-level covariate names)'; [0] 
%    Analysis.conditions            : (for modulation==1 only) list of task condition names to be simultaneously entered in gPPI 
%                                     model (leave empty for default 'all existing conditions') [] 
%    Analysis.type                  : analysis type, 1 = 'ROI-to-ROI', 2 = 'Seed-to-Voxel', 3 = 'all'; [3] 
%    Analysis.sources               : Cell array of sources names (seeds) (source names can be: any ROI name) (if this variable does 
%                                     not exist the toolbox will perform the analyses for all of the existing ROIs which are not 
%                                     defined as confounds in the Denoising step) 
%                                     note: partial ROI name matches are accepted (e.g. using 'networks.Language' will match all
%                                     ROIs starting with that token)
% 
%    Analysis.sources               : alternatively sources can be a structure with fields
%      Analysis.sources.names       : sources.names{nsource} char array of source names (seeds)
%      Analysis.sources.dimensions  : sources.dimensions{nsource} number of source dimensions [1]
%      Analysis.sources.deriv       : sources.deriv{nsource} number of derivatives for each dimension [0]
%      Analysis.sources.fbands      : sources.fbands{nsource} number of frequency bands for each dimension [1]
% 
% 
% BATCH.vvAnalysis PERFORMS FIRST-LEVEL ANALYSES (voxel-to-voxel) %!
%
%    vvAnalysis.done                : 1/0: 0 defines fields only; 1 runs ANALYSIS processing steps [0]
%    vvAnalysis.overwrite           : (for done=1) 1/0: overwrites target files if they exist [1]
%    vvAnalysis.name                : analysis name (identifying each set of independent analysis)
%                                    (alternatively sequential index identifying each set of independent analyses [1])
%  
%    vvAnalysis.measures            : voxel-to-voxel measure name (type 'conn_v2v measurenames' for a list of default measures) (if 
%                                     this variable does not exist the toolbox will perform the analyses for all of the default 
%                                     voxel-to-voxel measures) 
%                           'group-PCA'             : Principal Component Analysis of BOLD timeseries
%                           'group-ICA'             : Independent Component Analysis of BOLD timeseries
%                           'group-MVPA'/'MCOR'     : MultiVoxel Pattern Analysis of connectivity patterns (MCOR)
%                           'IntrinsicConnectivity' : Intrinsic Connectivity Contrast (pICC0)
%                           'LocalCorrelation'      : Integrated Local Correlation (ILC,LCOR)     
%                           'InterHemisphericCorrelation' : Inter-hemispheric Correlation (IHC)
%                           'GlobalCorrelation'     : Integrated Global Correlation (IGC,GCOR)   
%                           'RadialCorrelation'     : Radial Correlation Contrast (RCC)
%                           'RadialSimilarity'      : Radial Similarity Contrast (RSC)
%                           'ALFF'                  : Amplitude of Low Frequency Fluctuations
%                           'fALFF'                 : fractional ALFF
%
%    vvAnalysis.measures            : alternatively voxel-to-voxel measures can be a structure with fields
%      vvAnalysis.measures.names: measures.names voxel-to-voxel measure name (see 'conn_v2v measurenames' for a list of valid 
%                                     measure names)
%      vvAnalysis.measures.factors  : (for group-PCA, group-ICA, group-MVPA) number of group-level components to estimate
%      vvAnalysis.measures.kernelsupport : (for ILC, RCC) local support (FWHM mm) of smoothing kernel [8]
%      vvAnalysis.measures.norm     : (for ILC,ICC,IHC,RCC,RSC,ALFF,fALFF) 0/1 normalize values to z-scores [1]
%      vvAnalysis.measures.mask     : (for group-PCA, group-ICA, group-MVPA) optional mask for group-level component estimation 
%                                     (e.g. masked ICA)
%      vvAnalysis.measures.options  : (for group-ICA) optional ICA method : string containing 'GICA1' or 'GICA3' for choice of ICA back-
%                                     projection method; string containing 'tanh','gauss', or 'pow3' for ICA estimation method (G1/G2/G3)
%      vvAnalysis.measures.dimensions : number of subject-level dimensions to retain (subject-level dimensionality reduction) [64]
%  
%  
% BATCH.dynAnalysis PERFORMS FIRST-LEVEL ANALYSES (dynamic connectivity) %!
%
%    dynAnalysis.done               : 1/0: 0 defines fields only; 1 runs ANALYSIS processing steps [0]
%    dynAnalysis.overwrite          : (for done=1) 1/0: overwrites target files if they exist [1]
%    dynAnalysis.name               : analysis name (identifying each set of independent analysis)
%                                     (alternatively sequential index identifying each set of independent analyses [1])
%  
%    dynAnalysis.sources            : Cell array of sources names (seeds) (source names can be: any ROI name) (if this variable does 
%                                     not exist the toolbox will perform the analyses for all of the existing ROIs which are not 
%                                     defined as confounds in the Denoising step) 
%    dynAnalysis.factors            : Number of group-level dynamic components to estimate [20]
%    dynAnalysis.window             : Length of temporal windows (FWHM in seconds) [30]
%  
%  
% BATCH.Results PERFORMS SECOND-LEVEL ANALYSES (ROI-to-ROI, Seed-to-Voxel, Voxel-to-Voxel analyses) %!
% 
%    Results.name                   : analysis name (identifying each set of first-level independent analysis)
%                                     (alternatively sequential index identifying each set of first-level independent analyses [1])
%    Results.display                : 1/0 display results [1]
%    Results.saveas                 : (optional) name to save between-subjects/between_conditions contrast
%    Results.foldername             : (optional) alternative folder name to store the results
% 
%    Results.between_subjects
%      Results.between_subjects.effect_names : cell array of second-level effect names
%      Results.between_subjects.contrast : contrast vector (same size as effect_names)
%  
%    Results.between_conditions [defaults to multiple analyses, one per condition]
%      Results.between_conditions.effect_names : cell array of condition names (as in Setup.conditions.names)
%      Results.between_conditions.contrast : contrast vector (same size as effect_names)
%  
%    Results.between_sources  [if analysis includes multiple sources or multiple measures defaults to one analysis per source/measure]
%      Results.between_sources.effect_names : cell array of source names 
%                                     (as in Analysis.regressors for seed-to-voxel analyses), for multivariate sources 
%                                     they are appended with _N_M -where N is an index ranging from 1 to 1+derivative order, and M 
%                                     is an index ranging from 1 to the number of dimensions specified for each ROI; for example 
%                                     ROINAME_2_3 corresponds to the first derivative of the third PCA component extracted from the 
%                                     roi ROINAME) 
%                                     (as in Analysis.measures for voxel-to-voxel analyses)
%      Results.between_sources.contrast : contrast vector (same size as effect_names)
%  
%  
%
% BATCH.QA PERFORMS QUALITY ASSURANCE PLOTS %!
% 
%    QA.foldername                  : output folder where QA plots will be stored [results/qa/QA_#date#]
%    QA.plots                       : list of QA plots to create (cell array of labels below)
%                        QA_NORM structural     (1) : structural data + outline of MNI TPM template
%                        QA_NORM functional     (2) : mean functional data + outline of MNI TPM template
%                        QA_NORM rois           (3) : ROI data + outline of MNI TPM template  
%                        QA_REG functional      (10): display mean functional data + structural data overlay
%                        QA_REG structural      (4) : structural data + outline of ROI
%                        QA_REG functional      (5) : mean functional data + outline of ROI
%                        QA_REG mni             (6) : reference MNI structural template + outline of ROI
%                        QA_COREG functional    (7) : display same single-slice (z=0) across multiple sessions/datasets
%                        QA_TIME functional     (8) : display same single-slice (z=0) across all timepoints within each session
%                        QA_TIMEART functional  (9) : display same single-slice (z=0) across all timepoints within each session together with 
%                                                     ART timeseries (global signal changes and framewise displacement)
%                        QA_DENOISE histogram  (11) : histogram of voxel-to-voxel correlation values (before and after denoising)
%                        QA_DENOISE timeseries (12) : BOLD signal traces before and after denoising
%                        QA_DENOISE FC-QC      (13) : histogram of FC-QC associations; between-subject correlation between QC (Quality Control) 
%                                                     and FC (Functional Connectivity) measures
%                        QA_DENOISE scatterplot(14) : scatterplot of FC (Functional Connectivity r coeff) vs. distance (mm)
%                        QA_SPM design         (21) : SPM review design matrix (from SPM.mat files only)
%                        QA_SPM contrasts      (22) : SPM review contrast specification (from SPM.mat files only)
%                        QA_SPM results        (23) : SPM review contrast effect-size (from SPM.mat files only)
%                        QA_COV                (31) : histogram display of second-level variables
%
%   QA.rois                         : (only for plots==3,4,5,6) ROI numbers to include (defaults to WM -roi#2-)
%   QA.sets                         : (only for plots==2,7,8,9,10) functional dataset number (defaults to dataset-0)
%   QA.l2covariates                 : (only for plots==13,31) l2 covariate names (defaults to all QC_*)
%   QA.l1contrasts                  : (only for plots==23) l1 contrast name (defaults to first contrast)
%   QA.conditions                   : (only for plots==11,13,14,15) FC & QC-FC plots aggregate across sesssions where 
%                                      the selected conditions are present (defaults to all sessions)
%__________________________________________________________________________________________________________________
% 
% See 
%   conn_batch_workshop_nyu.m 
%   conn_batch_workshop_nyu_parallel.m 
%   conn_batch_humanconnectomeproject.m 
% for additional information and examples of use.
% 

%$

global CONN_x;
varargout={};
if ~nargin, help(mfilename); return; end

if nargin==1&&~ischar(varargin{1}), batch=varargin{1}; %batch(BATCH) syntax
elseif nargin==1&&ischar(varargin{1}), 
    if ~isempty(regexp(varargin{1},'\.mat$'))
        data=conn_loadmatfile(varargin{1},'-mat'); tnames=fieldnames(data); batch=data.(tnames{1}); %batch(batchfilename.mat) syntax
    elseif ~isempty(regexp(varargin{1},'\.m$')) %batch(batchfilename.m) syntax
        conn_batch_eval(varargin{1});
        return;
    elseif ~isempty(regexp(varargin{1},'\.json$'))
        batch=conn_jsonread(varargin{1}); %batch(batchfilename.json) syntax
        tnames=fieldnames(batch); if numel(tnames)==1&&strcmpi(tnames{1},'batch'), batch=batch.(tnames{1}); end
    else %batch(Matlabcommand) syntax
        evalin('base',varargin{1});
        return;
    end
elseif nargin>=1&&isa(varargin{1},'function_handle') %batch(@fcn,params)
    feval(varargin{:});
    return;
else %batch(fieldname,fieldvalue,...) syntax
    batch=[];
    for n=1:2:nargin-1
        str=regexp(varargin{n},'\.','split');
        batch=setfield(batch,str{:},varargin{n+1});
    end
end
if iscell(batch),%batch({BATCH1, BATCH2...}) syntax
    for nbatch=1:numel(batch),conn_batch(batch{nbatch});end
    return;
elseif numel(batch)>1, %batch([BATCH1 BATCH2]) syntax
    for nbatch=1:numel(batch),conn_batch(batch(nbatch));end
    return;
end

%% NEW step
if isfield(batch,'New'), % obsolete functionality / for backwards compatibility only
    if isfield(batch,'parallel')&&isfield(batch.parallel,'N')&&batch.parallel.N>0, error('BATCH.New option is not compatible with parallel processing. Please use the newer BATCH.Setup.preprocessing fields instead'); return; end
    conn_disp('Use of BATCH.New for preprocessing is no longer supported and may become obsolete in future releases. Please modify your scripts using BATCH.New to use instead BATCH.Setup.preprocessing (see doc conn_batch for details)'); 
    OPTIONS=struct('RT',2,'FWHM',8,'VOX',2,'CONN_DISPLAY',0,'STRUCTURAL_TEMPLATE',fullfile(fileparts(which('spm')),'templates','T1.nii'),'FUNCTIONAL_TEMPLATE',fullfile(fileparts(which('spm')),'templates','EPI.nii'),'SO',[],'UNWARP',[]);
    if ~conn_existfile(OPTIONS.FUNCTIONAL_TEMPLATE), OPTIONS.FUNCTIONAL_TEMPLATE=fullfile(fileparts(which('spm')),'toolbox','OldNorm','EPI.nii'); end
    if ~conn_existfile(OPTIONS.STRUCTURAL_TEMPLATE), OPTIONS.STRUCTURAL_TEMPLATE=fullfile(fileparts(which('spm')),'toolbox','OldNorm','T1.nii'); end
    if isfield(batch,'filename')&&~isempty(batch.filename),OPTIONS.CONN_NAME=batch.filename; end
    if isfield(batch.New,'center')&&~isempty(batch.New.center),OPTIONS.CENTER=batch.New.center;end
    if isfield(batch.New,'reorient')&&~isempty(batch.New.reorient),OPTIONS.REORIENT=batch.New.reorient;end
    if isfield(batch.New,'RT')&&~isempty(batch.New.RT),OPTIONS.RT=batch.New.RT;end; 
    if isfield(batch.New,'fwhm')&&~isempty(batch.New.fwhm),OPTIONS.FWHM=batch.New.fwhm;end
    if isfield(batch.New,'FWHM')&&~isempty(batch.New.FWHM),OPTIONS.FWHM=batch.New.FWHM;end
    if isfield(batch.New,'VOX')&&~isempty(batch.New.VOX),OPTIONS.VOX=batch.New.VOX;end
    if isfield(batch.New,'sliceorder')&&~isempty(batch.New.sliceorder),OPTIONS.SO=batch.New.sliceorder;end
    if isfield(batch.New,'unwarp')&&~isempty(batch.New.unwarp),OPTIONS.UNWARP=batch.New.unwarp;end
    if isfield(batch.New,'removescans')&&~isempty(batch.New.removescans),OPTIONS.removescans=batch.New.removescans;end
    if isfield(batch.New,'coregtomean')&&~isempty(batch.New.coregtomean),OPTIONS.coregtomean=batch.New.coregtomean;end
    if isfield(batch.New,'applytofunctional')&&~isempty(batch.New.applytofunctional),OPTIONS.applytofunctional=batch.New.applytofunctional;end
    if isfield(batch.New,'art_thresholds')&&~isempty(batch.New.art_thresholds),OPTIONS.art_thresholds=batch.New.art_thresholds;end
    if isfield(batch.New,'steps')&&~isempty(batch.New.steps),OPTIONS.STEPS=batch.New.steps;end
    if isfield(batch.New,'template_structural')&&~isempty(batch.New.template_structural),OPTIONS.STRUCTURAL_TEMPLATE=batch.New.template_structural;end
    if isfield(batch.New,'template_functional')&&~isempty(batch.New.template_functional),OPTIONS.FUNCTIONAL_TEMPLATE=batch.New.template_functional;end
    if isfield(batch.New,'tpm_template'),OPTIONS.tpm_template=batch.New.tpm_template;end
    if isfield(batch.New,'tpm_structlesion'),OPTIONS.tpm_structlesion=batch.New.tpm_structlesion;end
    if isfield(batch.New,'tpm_ngaus'),OPTIONS.tpm_ngaus=batch.New.tpm_ngaus;end
    if isfield(batch.New,'functionals')&&~isempty(batch.New.functionals),
        OPTIONS.FUNCTIONAL_FILES=batch.New.functionals;
    end
    if isfield(batch.New,'structurals')&&~isempty(batch.New.structurals),
        OPTIONS.STRUCTURAL_FILES=batch.New.structurals;
    end
    conn_setup_wizard(OPTIONS);
end

PAR_CMD={};
PAR_ARG={};

%% SETUP step
if isfield(batch,'Setup'),
    if isfield(batch,'filename'),
        if ~isempty(batch.filename)&&ischar(batch.filename), [nill,nill,extt]=fileparts(batch.filename); if isempty(extt), batch.filename=[batch.filename,'.mat']; end; end
        if (isfield(batch.Setup,'isnew')&&batch.Setup.isnew)||~conn_existfile(batch.filename),
            conn init;                   % initializes CONN_x structure
            %CONN_x.Setup.RT=nan;
            CONN_x.filename=batch.filename;
            if conn_existfile(CONN_x.filename), conn_jobmanager('cleardmat'); end
        else
            CONN_x.filename=batch.filename;
            CONN_x.gui=0;
            conn load;                      % loads existing conn_* project
            CONN_x.gui=1;
        end
    end
    
    if ~isfield(batch.Setup,'overwrite'),batch.Setup.overwrite='Yes';end
    if isscalar(batch.Setup.overwrite)&&~isstruct(batch.Setup.overwrite)&&ismember(double(batch.Setup.overwrite),[1 89 121]), batch.Setup.overwrite='Yes'; end
    if isfield(batch.Setup,'add')&&batch.Setup.add, 
        SUBJECTS=CONN_x.Setup.nsubjects+(1:batch.Setup.nsubjects); 
        CONN_x.Setup.nsubjects=conn_merge(CONN_x.Setup.nsubjects,CONN_x.Setup.nsubjects+batch.Setup.nsubjects); 
    else
        if isfield(batch.Setup,'nsubjects')&&~isempty(batch.Setup.nsubjects),
            if batch.Setup.nsubjects~=CONN_x.Setup.nsubjects, CONN_x.Setup.nsubjects=conn_merge(CONN_x.Setup.nsubjects,batch.Setup.nsubjects); end
        end
        if isfield(batch,'subjects'), SUBJECTS=batch.subjects; 
        else SUBJECTS=1:CONN_x.Setup.nsubjects;
        end
    end
    if isfield(batch.Setup,'spmfiles')&&~isempty(batch.Setup.spmfiles),
        CONN_x.gui=struct('overwrite',batch.Setup.overwrite);
        if isfield(batch.Setup,'spmfiles_options')&&~isempty(batch.Setup.spmfiles_options),args=batch.Setup.spmfiles_options; 
        else args={};
        end
        conn_importspm(batch.Setup.spmfiles,args{:},'subjects',SUBJECTS);
        CONN_x.gui=1;
    end
    if isfield(batch.Setup,'RT')&&~isempty(batch.Setup.RT),CONN_x.Setup.RT(SUBJECTS)=batch.Setup.RT;
    elseif isfield(batch.Setup,'rt')&&~isempty(batch.Setup.rt),CONN_x.Setup.RT(SUBJECTS)=batch.Setup.rt;
    end
    if isfield(batch.Setup,'acquisitiontype')&&~isempty(batch.Setup.acquisitiontype),
        CONN_x.Setup.acquisitiontype=1+(batch.Setup.acquisitiontype~=1);
    end
    if isfield(batch.Setup,'analyses'),
        CONN_x.Setup.steps=accumarray(batch.Setup.analyses(:),1,[4,1])';
    end
    if isfield(batch.Setup,'voxelmask')&&~isempty(batch.Setup.voxelmask),
        CONN_x.Setup.analysismask=batch.Setup.voxelmask;
    end
    if isfield(batch.Setup,'voxelmaskfile')&&~isempty(batch.Setup.voxelmaskfile),
        CONN_x.Setup.explicitmask=conn_file(batch.Setup.voxelmaskfile);
    end
    if isfield(batch.Setup,'voxelresolution')&&~isempty(batch.Setup.voxelresolution),
        CONN_x.Setup.spatialresolution=batch.Setup.voxelresolution;
    end
    if isfield(batch.Setup,'analysisunits')&&~isempty(batch.Setup.analysisunits),
        CONN_x.Setup.analysisunits=batch.Setup.analysisunits;
    end
    if isfield(batch.Setup,'outputfiles'),
        CONN_x.Setup.outputfiles=batch.Setup.outputfiles;
    end
    if isfield(batch.Setup,'surfacesmoothing'),
        conn_disp('warning: surfacesmoothing field is no longer supported. Please use explicit ''surface_smoothing'' preprocessing step instead to smooth surface-level functional data');
        %CONN_x.Setup.surfacesmoothing=batch.Setup.surfacesmoothing;
    end
    if isfield(batch.Setup,'functionals')&&~isempty(batch.Setup.functionals),
        localcopy=false; if isfield(batch.Setup,'localcopy')&&batch.Setup.localcopy, localcopy=true; end
        localcopy_reduce=false; if isfield(batch.Setup,'localcopy_reduce')&&batch.Setup.localcopy_reduce, localcopy_reduce=true; end
        if localcopy, conn_updatefolders; end
        if ~iscell(batch.Setup.functionals), batch.Setup.functionals={batch.Setup.functionals}; end
        for isub=1:numel(SUBJECTS),
            nsub=SUBJECTS(isub);
            if ~iscell(batch.Setup.functionals{isub}), batch.Setup.functionals{isub}={batch.Setup.functionals{isub}}; end
            CONN_x.Setup.nsessions(nsub)=length(batch.Setup.functionals{isub});
            for nses=1:CONN_x.Setup.nsessions(nsub),
                if localcopy, [nill,nill,nV]=conn_importvol2bids(batch.Setup.functionals{isub}{nses},nsub,nses,'func',[],[],[],localcopy_reduce);
                else
                    [CONN_x.Setup.functional{nsub}{nses},nV]=conn_file(batch.Setup.functionals{isub}{nses});
                    CONN_x.Setup.nscans{nsub}{nses}=nV;
                end
                
            end
        end
    end
%     if isfield(batch.Setup,'roiextract'), deal(batch.Setup.secondarydatasets.roiextract=batch.Setup.roiextract; end
%     if isfield(batch.Setup,'roiextract_rule'), batch.Setup.secondarydatasets.roiextract_rule=batch.Setup.roiextract_rule; end
%     if isfield(batch.Setup,'roiextract_functionals'), batch.Setup.secondarydatasets.roiextract_functionals=batch.Setup.roiextract_functionals; end
    if isfield(batch.Setup,'roifunctionals')&&~isfield(batch.Setup,'secondarydatasets'), batch.Setup.secondarydatasets=batch.Setup.roifunctionals; end
    if isfield(batch.Setup,'secondarydatasets')
        localcopy=false; if isfield(batch.Setup,'localcopy')&&batch.Setup.localcopy, localcopy=true; end
        localcopy_reduce=false; if isfield(batch.Setup,'localcopy_reduce')&&batch.Setup.localcopy_reduce, localcopy_reduce=true; end
        for nalt=1:numel(batch.Setup.secondarydatasets)
            if iscell(batch.Setup.secondarydatasets), tsecondarydataset=batch.Setup.secondarydatasets{nalt};
            else tsecondarydataset=batch.Setup.secondarydatasets(nalt);
            end
            ialt=nalt;
            if isfield(tsecondarydataset,'functionals_rule')
                CONN_x.Setup.secondarydataset(ialt).functionals_rule=tsecondarydataset.functionals_rule;
            elseif isfield(tsecondarydataset,'roiextract_rule')
                CONN_x.Setup.secondarydataset(ialt).functionals_rule=tsecondarydataset.roiextract_rule;
            end
            if isfield(tsecondarydataset,'functionals_label')
                CONN_x.Setup.secondarydataset(ialt).label=tsecondarydataset.functionals_label;
            end
            if isfield(tsecondarydataset,'functionals_explicit')||isfield(tsecondarydataset,'roiextract_functionals')
                if isfield(tsecondarydataset,'functionals_explicit'), temp=tsecondarydataset.functionals_explicit;
                else temp=tsecondarydataset.roiextract_functionals;
                end
                if ~isempty(temp)
                    if ~iscell(temp), temp={temp}; end
                    for isub=1:numel(SUBJECTS),
                        nsub=SUBJECTS(isub);
                        %CONN_x.Setup.nsessions(nsub)=length(batch.Setup.functionals_explicit{nsub});
                        if ~isempty(temp{isub})
                            if ~iscell(temp{isub}), temp{isub}={temp}; end
                            for nses=1:CONN_x.Setup.nsessions(nsub),
                                if localcopy,
                                    switch(char(CONN_x.Setup.secondarydataset(ialt).label))
                                        case 'fmap', args={'fmap','fmap'}; % copies dataset 'fmap' as .../fmap/*fmap.nii file
                                        case 'vdm',  args={'fmap','vdm'};  % copies dataset 'vdm' as .../fmap/*vdm.nii file
                                        otherwise, args={CONN_x.Setup.secondarydataset(ialt).label,[]}; % copies other datasets as ../dataset<label>/*unknown.nii file
                                    end
                                    [nill,nill,nV]=conn_importvol2bids(temp{isub}{nses},nsub,nses,args{:},[],[],localcopy_reduce);
                                else conn_set_functional(nsub,nses,ialt,temp{isub}{nses});
                                end
                                %[CONN_x.Setup.secondarydataset(ialt).functionals_explicit{nsub}{nses},nV]=conn_file(temp{isub}{nses});
                            end
                        end
                    end
                end
            end
            if isfield(tsecondarydataset,'functionals_type')
                CONN_x.Setup.secondarydataset(ialt).functionals_type=tsecondarydataset.functionals_type;
            elseif isfield(tsecondarydataset,'roiextract')
                CONN_x.Setup.secondarydataset(ialt).functionals_type=tsecondarydataset.roiextract;
            end
        end
    end
    if isfield(batch.Setup,'vdm_functionals')&&~isfield(batch.Setup,'unwarp_functionals'), batch.Setup.unwarp_functionals=batch.Setup.vdm_functionals; end
    if isfield(batch.Setup,'unwarp_functionals')
        localcopy=false; if isfield(batch.Setup,'localcopy')&&batch.Setup.localcopy, localcopy=true; end
        localcopy_reduce=false; if isfield(batch.Setup,'localcopy_reduce')&&batch.Setup.localcopy_reduce, localcopy_reduce=true; end
        for isub=1:numel(SUBJECTS),
            nsub=SUBJECTS(isub);
            %CONN_x.Setup.nsessions(nsub)=length(batch.Setup.unwarp_functionals{nsub});
            if ~iscell(batch.Setup.unwarp_functionals{isub}), batch.Setup.unwarp_functionals{isub}={batch.Setup.unwarp_functionals{isub}}; end
            for nses=1:CONN_x.Setup.nsessions(nsub),
                fname=batch.Setup.unwarp_functionals{isub}{min(numel(batch.Setup.unwarp_functionals{isub}),nses)};
                if localcopy, [nill,nill,nV]=conn_importvol2bids(fname,nsub,nses,'fmap','vdm',[],[],localcopy_reduce);
                else
                    [CONN_x.Setup.unwarp_functional{nsub}{nses},nV]=conn_file(fname);
                    conn_set_functional(nsub,nses,'vdm',fname);
                end
                %CONN_x.Setup.nscans{nsub}{nses}=nV;
            end
        end
    end
    if isfield(batch.Setup,'fmap_functionals')
        localcopy=false; if isfield(batch.Setup,'localcopy')&&batch.Setup.localcopy, localcopy=true; end
        localcopy_reduce=false; if isfield(batch.Setup,'localcopy_reduce')&&batch.Setup.localcopy_reduce, localcopy_reduce=true; end
        for isub=1:numel(SUBJECTS),
            nsub=SUBJECTS(isub);
            %CONN_x.Setup.nsessions(nsub)=length(batch.Setup.fmap_functionals{nsub});
            if ~iscell(batch.Setup.fmap_functionals{isub}), batch.Setup.fmap_functionals{isub}={batch.Setup.fmap_functionals{isub}}; end
            for nses=1:CONN_x.Setup.nsessions(nsub),
                fname=batch.Setup.fmap_functionals{isub}{min(numel(batch.Setup.fmap_functionals{isub}),nses)};
                if localcopy, [nill,nill,nV]=conn_importvol2bids(fname,nsub,nses,'fmap','fmap',[],[],localcopy_reduce);
                else
                    %[CONN_x.Setup.fmap_functional{nsub}{nses},nV]=conn_file(fname);
                    conn_set_functional(nsub,nses,'fmap',fname);
                end
                %CONN_x.Setup.nscans{nsub}{nses}=nV;
            end
        end
    end
    if isfield(batch.Setup,'coregsource_functionals')
        for isub=1:numel(SUBJECTS),
            nsub=SUBJECTS(isub);
            [CONN_x.Setup.coregsource_functional{nsub},nV]=conn_file(batch.Setup.coregsource_functionals{isub});
        end
    end
    
    allfields=fieldnames(batch.Setup);
    idx=cellfun('length',regexp(allfields,'_functionals$'))>0;
    idx=idx(~ismember(allfields(idx),{'fmap_functionals','vdm_functionals','unwarp_functionals','coregsource_functionals','roiextract_functionals'})); % *_functionals interpreted as other secondary datasets
    for nidx=1:numel(idx) % '*_functionals'
        setupfieldname=allfields{idx(nidx)};
        datasetname=regexprep(setupfieldname,'_functionals$','');
        localcopy=false; if isfield(batch.Setup,'localcopy')&&batch.Setup.localcopy, localcopy=true; end
        localcopy_reduce=false; if isfield(batch.Setup,'localcopy_reduce')&&batch.Setup.localcopy_reduce, localcopy_reduce=true; end
        for isub=1:numel(SUBJECTS),
            nsub=SUBJECTS(isub);
            if ~iscell(batch.Setup.(setupfieldname){isub}), batch.Setup.(setupfieldname){isub}={batch.Setup.(setupfieldname){isub}}; end
            for nses=1:CONN_x.Setup.nsessions(nsub),
                fname=batch.Setup.(setupfieldname){isub}{min(numel(batch.Setup.(setupfieldname){isub}),nses)};
                if localcopy, [nill,nill,nV]=conn_importvol2bids(fname,nsub,nses,datasetname,[],[],[],localcopy_reduce);
                else
                    %[CONN_x.Setup.fmap_functional{nsub}{nses},nV]=conn_file(fname);
                    conn_set_functional(nsub,nses,datasetname,fname);
                end
                %CONN_x.Setup.nscans{nsub}{nses}=nV;
            end
        end
    end
    
    if isfield(batch.Setup,'cwthreshold')
        if ~isfield(batch.Setup,'binary_threshold_type'), batch.Setup.binary_threshold_type=[1 1 1]; end
        if ~isfield(batch.Setup,'exclude_grey_matter'), batch.Setup.exclude_grey_matter=[nan nan nan]; end
        if numel(batch.Setup.cwthreshold)==1
            if ~isfield(batch.Setup,'binary_threshold'), batch.Setup.binary_threshold=[.5 batch.Setup.cwthreshold batch.Setup.cwthreshold]; end
        elseif numel(batch.Setup.cwthreshold)==2
            if ~isfield(batch.Setup,'binary_threshold'), batch.Setup.binary_threshold=[.5 batch.Setup.cwthreshold(1) batch.Setup.cwthreshold(1)]; end
            if ~isfield(batch.Setup,'erosion_steps'), batch.Setup.erosion_steps=[0 batch.Setup.cwthreshold(2) batch.Setup.cwthreshold(2)]; end
        elseif numel(batch.Setup.cwthreshold)==4
            if ~isfield(batch.Setup,'binary_threshold'), batch.Setup.binary_threshold=[.5 batch.Setup.cwthreshold([1,3])]; end
            if ~isfield(batch.Setup,'erosion_steps'), batch.Setup.erosion_steps=[0 batch.Setup.cwthreshold([2,4])]; end
        else
            error('unexpected batch.Setup.cwthreshold field values. Please use batch.Setup.binary_threshold, batch.Setup.erosion_steps, and batch.Setup.erosion_neighb fields instead');
        end
    end
    if isfield(batch.Setup,'binary_threshold')
        CONN_x.Setup.erosion.binary_threshold(1:numel(batch.Setup.binary_threshold))=batch.Setup.binary_threshold;
    end
    if isfield(batch.Setup,'binary_threshold_type')
        CONN_x.Setup.erosion.binary_threshold_type(1:numel(batch.Setup.binary_threshold_type))=batch.Setup.binary_threshold_type;
    end
    if isfield(batch.Setup,'exclude_grey_matter')
        CONN_x.Setup.erosion.exclude_grey_matter(1:numel(batch.Setup.exclude_grey_matter))=batch.Setup.exclude_grey_matter;
    end
    if isfield(batch.Setup,'erosion_steps')
        CONN_x.Setup.erosion.erosion_steps(1:numel(batch.Setup.erosion_steps))=batch.Setup.erosion_steps;
    end
    if isfield(batch.Setup,'erosion_neighb')
        CONN_x.Setup.erosion.erosion_neighb(1:numel(batch.Setup.erosion_neighb))=batch.Setup.erosion_neighb;
    end

    if isfield(batch.Setup,'structurals')&&~isempty(batch.Setup.structurals),
        localcopy=false; if isfield(batch.Setup,'localcopy')&&batch.Setup.localcopy, localcopy=true; end
        localcopy_reduce=false; if isfield(batch.Setup,'localcopy_reduce')&&batch.Setup.localcopy_reduce, localcopy_reduce=true; end
        if localcopy, conn_updatefolders; end
        CONN_x.Setup.structural_sessionspecific=0; 
        for isub=1:numel(SUBJECTS),
            nsub=SUBJECTS(isub);
            nsess=CONN_x.Setup.nsessions(nsub);
            temp=batch.Setup.structurals{isub};
            if ischar(temp), temp={temp}; end
            if numel(temp)>1, CONN_x.Setup.structural_sessionspecific=1; end
            if ~CONN_x.Setup.structural_sessionspecific, nsesstemp=1; else nsesstemp=nsess; end
            for nses=1:nsess,
                if nses>nsesstemp, CONN_x.Setup.structural{nsub}{nses}=CONN_x.Setup.structural{nsub}{nsesstemp};
                elseif localcopy, conn_importvol2bids(temp{min(numel(temp),nses)},nsub,nses,'anat',[],[],[],localcopy_reduce);
                else CONN_x.Setup.structural{nsub}{nses}=conn_file(temp{min(numel(temp),nses)});
                end
            end
        end
    end
    if isfield(batch.Setup,'masks'),
        localcopy=false; if isfield(batch.Setup,'localcopy')&&batch.Setup.localcopy, localcopy=true; end 
        localcopy_reduce=false; if isfield(batch.Setup,'localcopy_reduce')&&batch.Setup.localcopy_reduce, localcopy_reduce=true; end
        masks={'Grey','White','CSF'};
        if ~isstruct(batch.Setup.masks), batch.Setup.masks=struct('Grey',{batch.Setup.masks}); end
        for nmask=1:length(masks),
            if isfield(batch.Setup.masks,masks{nmask})&&~isempty(batch.Setup.masks.(masks{nmask})),
                if ~isstruct(batch.Setup.masks.(masks{nmask})),
                    if isequal(SUBJECTS,1:CONN_x.Setup.nsubjects)
                        subjectspecific=0;
                        sessionspecific=0;
                    else
                        subjectspecific=CONN_x.Setup.rois.subjectspecific(nmask);
                        sessionspecific=CONN_x.Setup.rois.sessionspecific(nmask);
                    end
                    temp1=batch.Setup.masks.(masks{nmask});
                    if ischar(temp1), temp1={temp1}; end
                    if numel(temp1)>1||CONN_x.Setup.nsubjects==1, subjectspecific=1; end
                    for isub=1:numel(SUBJECTS),
                        nsub=SUBJECTS(isub);
                        temp2=temp1{min(numel(temp1),isub)};
                        if ischar(temp2), temp2={temp2}; end
                        if numel(temp2)>1, sessionspecific=1; end
                        for nses=1:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsub))
                            if localcopy&&subjectspecific, 
                                CONN_x.Setup.rois.subjectspecific(nmask)=subjectspecific;
                                CONN_x.Setup.rois.sessionspecific(nmask)=sessionspecific;
                                if sessionspecific, conn_importvol2bids(temp2{min(numel(temp2),nses)},nsub,nses,'roi',[],[],[],localcopy_reduce,nmask);
                                elseif nses==1, conn_importvol2bids(temp2{min(numel(temp2),nses)},nsub,[],'roi',[],[],[],localcopy_reduce,nmask);
                                else CONN_x.Setup.rois.files{nsub}{nmask}{nses}=CONN_x.Setup.rois.files{nsub}{nmask}{1};
                                end
                            else CONN_x.Setup.rois.files{nsub}{nmask}{nses}=conn_file(temp2{min(numel(temp2),nses)});
                            end
                        end
                    end
                    CONN_x.Setup.rois.subjectspecific(nmask)=subjectspecific;
                    CONN_x.Setup.rois.sessionspecific(nmask)=sessionspecific;
                else,
                    if isequal(SUBJECTS,1:CONN_x.Setup.nsubjects)
                        subjectspecific=0;
                        sessionspecific=0;
                    else
                        subjectspecific=CONN_x.Setup.rois.subjectspecific(nmask);
                        sessionspecific=CONN_x.Setup.rois.sessionspecific(nmask);
                    end
                    temp1=batch.Setup.masks.(masks{nmask}).files;
                    if ischar(temp1), temp1={temp1}; end
                    if numel(temp1)>1||CONN_x.Setup.nsubjects==1, subjectspecific=1; end
                    for isub=1:numel(SUBJECTS),
                        nsub=SUBJECTS(isub);
                        temp2=temp1{min(numel(temp1),isub)};
                        if ischar(temp2), temp2={temp2}; end
                        if numel(temp2)>1, sessionspecific=1; end
                        for nses=1:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsub))
                            if localcopy&&subjectspecific, 
                                CONN_x.Setup.rois.subjectspecific(nmask)=subjectspecific;
                                CONN_x.Setup.rois.sessionspecific(nmask)=sessionspecific;
                                if sessionspecific, conn_importvol2bids(temp2{min(numel(temp2),nses)},nsub,nses,'roi',[],[],[],localcopy_reduce,nmask);
                                elseif nses==1, conn_importvol2bids(temp2{min(numel(temp2),nses)},nsub,[],'roi',[],[],[],localcopy_reduce,nmask);
                                else CONN_x.Setup.rois.files{nsub}{nmask}{nses}=CONN_x.Setup.rois.files{nsub}{nmask}{1};
                                end
                            else CONN_x.Setup.rois.files{nsub}{nmask}{nses}=conn_file(temp2{min(numel(temp2),nses)});
                            end
                        end
                    end
                    if isfield(batch.Setup.masks.(masks{nmask}),'dimensions')&&~iscell(batch.Setup.masks.(masks{nmask}).dimensions), batch.Setup.masks.(masks{nmask}).dimensions=num2cell(batch.Setup.masks.(masks{nmask}).dimensions); end
                    if isfield(batch.Setup.masks.(masks{nmask}),'dimensions'), CONN_x.Setup.rois.dimensions{nmask}=batch.Setup.masks.(masks{nmask}).dimensions; if CONN_x.Setup.rois.dimensions{nmask}==0, CONN_x.Setup.rois.weighted(nmask)=1; end; end
                    if isfield(batch.Setup.masks.(masks{nmask}),'regresscovariates'), CONN_x.Setup.rois.regresscovariates(nmask)=batch.Setup.masks.(masks{nmask}).regresscovariates; end
                    if isfield(batch.Setup.masks.(masks{nmask}),'weighted'), CONN_x.Setup.rois.weighted(nmask)=batch.Setup.masks.(masks{nmask}).weighted; end
                    if isfield(batch.Setup.masks.(masks{nmask}),'dataset'), CONN_x.Setup.rois.unsmoothedvolumes(nmask)=batch.Setup.masks.(masks{nmask}).dataset; 
                    elseif isfield(batch.Setup.masks.(masks{nmask}),'roiextract'), CONN_x.Setup.rois.unsmoothedvolumes(nmask)=batch.Setup.masks.(masks{nmask}).roiextract; 
                    end
                    CONN_x.Setup.rois.subjectspecific(nmask)=subjectspecific;
                    CONN_x.Setup.rois.sessionspecific(nmask)=sessionspecific;
                end
            end
        end
    end
%     if isfield(batch.Setup,'masks')&&isfield(batch.Setup.masks,'Grey')&&~isempty(batch.Setup.masks.Grey),for nsub=1:CONN_x.Setup.nsubjects,CONN_x.Setup.rois.files{nsub}{1}=conn_file(batch.Setup.masks.Grey{nsub});end; end
%     if isfield(batch.Setup,'masks')&&isfield(batch.Setup.masks,'White')&&~isempty(batch.Setup.masks.White),for nsub=1:CONN_x.Setup.nsubjects,CONN_x.Setup.rois.files{nsub}{2}=conn_file(batch.Setup.masks.White{nsub});end; end
%     if isfield(batch.Setup,'masks')&&isfield(batch.Setup.masks,'CSF')&&~isempty(batch.Setup.masks.CSF),for nsub=1:CONN_x.Setup.nsubjects,CONN_x.Setup.rois.files{nsub}{3}=conn_file(batch.Setup.masks.CSF{nsub});end; end
    if isfield(batch.Setup,'rois'),%&&~isempty(batch.Setup.rois),
        localcopy=false; if isfield(batch.Setup,'localcopy')&&batch.Setup.localcopy, localcopy=true; end 
        localcopy_reduce=false; if isfield(batch.Setup,'localcopy_reduce')&&batch.Setup.localcopy_reduce, localcopy_reduce=true; end
        if ~isstruct(batch.Setup.rois), 
            temp=batch.Setup.rois;
            batch.Setup.rois=struct;
            batch.Setup.rois.files=temp;
            for n1=1:length(temp), ttemp=temp{n1}; while iscell(ttemp), ttemp=ttemp{1}; end; [nill,name,nameext]=fileparts(ttemp); batch.Setup.rois.names{n1}=name;end; 
        end
        if isfield(batch.Setup.rois,'add')&&batch.Setup.rois.add, 
            if isfield(batch.Setup,'isnew')&&batch.Setup.isnew, conn importrois; end
            n0=length(CONN_x.Setup.rois.names)-1;
        else
            n0=3;
            try, if isfield(batch.Setup.rois,'names')&&numel(batch.Setup.rois.names)>=3&&isequal(regexprep(batch.Setup.rois.names(1:3),'\s',''),regexprep(CONN_x.Setup.rois.names(1:3),'\s','')), n0=0; end; end
        end
        if isfield(batch.Setup.rois,'dataset')&&ischar(batch.Setup.rois.dataset), batch.Setup.rois.dataset={batch.Setup.rois.dataset}; end
        if isfield(batch.Setup.rois,'dataset')&&iscell(batch.Setup.rois.dataset), 
            labels=batch.Setup.rois.dataset;
            batch.Setup.rois.dataset=zeros(size(batch.Setup.rois.dataset));
            ok=cellfun(@isnumeric,labels);
            batch.Setup.rois.dataset(ok)=cell2mat(labels(ok));
            for nok=reshape(find(~ok),1,[])
                [tidx,tisnew]=conn_datasetlabel(labels{nok});
                if ~tisnew, batch.Setup.rois.dataset(nok)=tidx; ok(nok)=true; end
            end
            %[ok(~ok),batch.Setup.rois.dataset(~ok)]=ismember(labels(~ok),{CONN_x.Setup.secondarydataset.label}); 
            %if any(~ok), [ok(~ok),batch.Setup.rois.dataset(~ok)]=ismember(lower(labels(~ok)),lower({CONN_x.Setup.secondarydataset.label})); end
            assert(all(ok),'could not find a Secondary Dataset label match to %s',sprintf('%s ',labels{~ok}));
        end
        for n1=1:length(batch.Setup.rois.files),
            if isequal(SUBJECTS,1:CONN_x.Setup.nsubjects)
                subjectspecific=0; 
                sessionspecific=0;
            else
                subjectspecific=CONN_x.Setup.rois.subjectspecific(n0+n1);
                sessionspecific=CONN_x.Setup.rois.sessionspecific(n0+n1);
            end            
            temp1=batch.Setup.rois.files{n1};
            if ischar(temp1), temp1={temp1}; end
            if numel(temp1)>1||CONN_x.Setup.nsubjects==1, subjectspecific=1; end
            for isub=1:numel(SUBJECTS),
                nsub=SUBJECTS(isub);
                temp2=temp1{min(numel(temp1),isub)};
                if ischar(temp2), temp2={temp2}; end
                if numel(temp2)>1, sessionspecific=1; end
                for nses=1:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsub))
                    [nill,name,nameext]=fileparts(temp2{min(numel(temp2),nses)});
                    %[V,str,icon]=conn_getinfo(batch.Setup.rois.files{n1}{nsub});
                    %CONN_x.Setup.rois.files{nsub}{n0+n1}={batch.Setup.rois.files{n1}{nsub},str,icon};
                    if localcopy&&subjectspecific, 
                        if ~isfield(batch.Setup.rois,'names')||length(batch.Setup.rois.names)<n1||isempty(batch.Setup.rois.names{n1}), batch.Setup.rois.names{n1}=name; end 
                        CONN_x.Setup.rois.names{n0+n1}=batch.Setup.rois.names{n1};
                        CONN_x.Setup.rois.subjectspecific(n0+n1)=subjectspecific;
                        CONN_x.Setup.rois.sessionspecific(n0+n1)=sessionspecific;
                        if sessionspecific, conn_importvol2bids(temp2{min(numel(temp2),nses)},nsub,nses,'roi',[],[],[],localcopy_reduce,n0+n1);
                        elseif nses==1, conn_importvol2bids(temp2{min(numel(temp2),nses)},nsub,[],'roi',[],[],[],localcopy_reduce,n0+n1);
                        else CONN_x.Setup.rois.files{nsub}{n0+n1}{nses}=CONN_x.Setup.rois.files{nsub}{n0+n1}{1};
                        end
                    else CONN_x.Setup.rois.files{nsub}{n0+n1}{nses}=conn_file(temp2{min(numel(temp2),nses)});
                    end
                    if isfield(batch.Setup.rois,'labelfiles')&&numel(batch.Setup.rois.labelfiles)>=n1
                        [nill,nill,nameext2]=fileparts(batch.Setup.rois.labelfiles{n1});
                        if ~ismember(nameext2,{'.txt','.csv','.xls'}), nameext2='.txt'; end
                        tmfile=conn_prepend('',CONN_x.Setup.rois.files{nsub}{n0+n1}{nses}{1},nameext2);
                        try, conn_fileutils('copyfile',batch.Setup.rois.labelfiles{n1},tmfile); end
                        %if ispc, [ok,msg]=system(sprintf('copy "%s" "%s"',batch.Setup.rois.labelfiles{n1},tmfile));
                        %else [ok,msg]=system(sprintf('cp ''%s'' ''%s''',batch.Setup.rois.labelfiles{n1},tmfile));
                        %end
                    end
                end
            end
            if ~isfield(batch.Setup.rois,'names')||length(batch.Setup.rois.names)<n1||isempty(batch.Setup.rois.names{n1}), batch.Setup.rois.names{n1}=name; end % note: need to set names first in case localcopy==1
            if isfield(batch.Setup.rois,'dimensions')&&~iscell(batch.Setup.rois.dimensions), batch.Setup.rois.dimensions=num2cell(batch.Setup.rois.dimensions); end
            if ~isfield(batch.Setup.rois,'dimensions')||length(batch.Setup.rois.dimensions)<n1||isempty(batch.Setup.rois.dimensions{n1}), batch.Setup.rois.dimensions{n1}=1; end
            if ~isfield(batch.Setup.rois,'mask')||length(batch.Setup.rois.mask)<n1, batch.Setup.rois.mask(n1)=0; end
            if ~isfield(batch.Setup.rois,'multiplelabels')||length(batch.Setup.rois.multiplelabels)<n1, batch.Setup.rois.multiplelabels(n1)=(strcmp(nameext,'.img')|strcmp(nameext,'.nii')|strcmp(nameext,'.mgz'))&(conn_existfile(conn_prepend('',CONN_x.Setup.rois.files{1}{n0+n1}{1}{1},'.txt'))|conn_existfile(conn_prepend('',CONN_x.Setup.rois.files{1}{n0+n1}{1}{1},'.csv'))|conn_existfile(conn_prepend('',CONN_x.Setup.rois.files{1}{n0+n1}{1}{1},'.xls'))); end
            if ~isfield(batch.Setup.rois,'regresscovariates')||length(batch.Setup.rois.regresscovariates)<n1, batch.Setup.rois.regresscovariates(n1)=double(batch.Setup.rois.dimensions{n1}>1); end
            if ~isfield(batch.Setup.rois,'weighted')||length(batch.Setup.rois.weighted)<n1, batch.Setup.rois.weighted(n1)=double(batch.Setup.rois.dimensions{n1}==0); end
            if ~isfield(batch.Setup.rois,'dataset')||length(batch.Setup.rois.dataset)<n1, 
                if isfield(batch.Setup.rois,'roiextract')&&length(batch.Setup.rois.roiextract)>=n1, batch.Setup.rois.dataset(n1)=batch.Setup.rois.roiextract(n1); 
                else batch.Setup.rois.dataset(n1)=1; 
                end
            end
            CONN_x.Setup.rois.names{n0+n1}=batch.Setup.rois.names{n1}; CONN_x.Setup.rois.names{n0+n1+1}=' ';
            CONN_x.Setup.rois.dimensions{n0+n1}=batch.Setup.rois.dimensions{n1};
            CONN_x.Setup.rois.mask(n0+n1)=batch.Setup.rois.mask(n1);
            CONN_x.Setup.rois.subjectspecific(n0+n1)=subjectspecific;
            CONN_x.Setup.rois.sessionspecific(n0+n1)=sessionspecific;
            CONN_x.Setup.rois.multiplelabels(n0+n1)=batch.Setup.rois.multiplelabels(n1);
            CONN_x.Setup.rois.regresscovariates(n0+n1)=batch.Setup.rois.regresscovariates(n1);
            CONN_x.Setup.rois.weighted(n0+n1)=batch.Setup.rois.weighted(n1);
            CONN_x.Setup.rois.unsmoothedvolumes(n0+n1)=batch.Setup.rois.dataset(n1);
        end
        for nsub=1:CONN_x.Setup.nsubjects,% disregards other existing rois
            CONN_x.Setup.rois.files{nsub}=CONN_x.Setup.rois.files{nsub}(1:n0+length(batch.Setup.rois.files)); 
        end
        CONN_x.Setup.rois.names=CONN_x.Setup.rois.names(1:n0+length(batch.Setup.rois.files)+1);
        CONN_x.Setup.rois.names{end}=' ';
        CONN_x.Setup.rois.dimensions=CONN_x.Setup.rois.dimensions(1:n0+length(batch.Setup.rois.files));
        CONN_x.Setup.rois.mask=CONN_x.Setup.rois.mask(1:n0+length(batch.Setup.rois.files));
        CONN_x.Setup.rois.subjectspecific=CONN_x.Setup.rois.subjectspecific(1:n0+length(batch.Setup.rois.files));
        CONN_x.Setup.rois.sessionspecific=CONN_x.Setup.rois.sessionspecific(1:n0+length(batch.Setup.rois.files));
        CONN_x.Setup.rois.multiplelabels=CONN_x.Setup.rois.multiplelabels(1:n0+length(batch.Setup.rois.files));
        CONN_x.Setup.rois.regresscovariates=CONN_x.Setup.rois.regresscovariates(1:n0+length(batch.Setup.rois.files));
        CONN_x.Setup.rois.weighted=CONN_x.Setup.rois.weighted(1:n0+length(batch.Setup.rois.files));
        CONN_x.Setup.rois.unsmoothedvolumes=CONN_x.Setup.rois.unsmoothedvolumes(1:n0+length(batch.Setup.rois.files));
    elseif isfield(batch.Setup,'isnew')&&batch.Setup.isnew,
        conn importrois;
    end

    if isfield(batch.Setup,'conditions')&&~isempty(batch.Setup.conditions)&&isfield(batch.Setup.conditions,'importfile'),
        if ~isfield(batch.Setup.conditions,'add')||~batch.Setup.conditions.add,CONN_x.Setup.conditions.names={' '};end
        opts={}; if isfield(batch.Setup.conditions,'importfile_options'), opts=batch.Setup.conditions.importfile_options; end
        if ~iscell(opts), opts={opts}; end
        conn_importcondition(batch.Setup.conditions.importfile,'subjects',SUBJECTS,opts{:});
    end
    if isfield(batch.Setup,'conditions')&&~isempty(batch.Setup.conditions)&&isfield(batch.Setup.conditions,'names'),
        if isfield(batch.Setup.conditions,'add')&&batch.Setup.conditions.add, nl0=numel(CONN_x.Setup.conditions.names)-1; 
        else nl0=0;
        end 
        nconditions=numel(batch.Setup.conditions.names);
        CONN_x.Setup.conditions.names(nl0+(1:nconditions+1))={batch.Setup.conditions.names{:},' '};
        for isub=1:numel(SUBJECTS),
            nsub=SUBJECTS(isub);
            for nconditions=1:length(batch.Setup.conditions.names),
                for nses=1:CONN_x.Setup.nsessions(nsub),
                    if isfield(batch.Setup.conditions,'onsets')
                        % note: accepts [ncond x nsub x nses] matrix with NaN indicating [] or inf values
                        if ~iscell(batch.Setup.conditions.onsets), batch.Setup.conditions.onsets=arrayfun(@(n)shiftdim(batch.Setup.conditions.onsets(n,:,:),1),1:size(batch.Setup.conditions.onsets,1),'uni',0); end
                        if ~iscell(batch.Setup.conditions.onsets{nconditions}), batch.Setup.conditions.onsets{nconditions}=arrayfun(@(n)shiftdim(batch.Setup.conditions.onsets{nconditions}(n,:,:),1),1:size(batch.Setup.conditions.onsets{nconditions},1),'uni',0); end
                        if ~iscell(batch.Setup.conditions.onsets{nconditions}{isub}), batch.Setup.conditions.onsets{nconditions}{isub}=arrayfun(@(n)shiftdim(batch.Setup.conditions.onsets{nconditions}{isub}(n,:,:),1),1:size(batch.Setup.conditions.onsets{nconditions}{isub},1),'uni',0); end
                        if ~iscell(batch.Setup.conditions.durations), batch.Setup.conditions.durations=arrayfun(@(n)shiftdim(batch.Setup.conditions.durations(n,:,:),1),1:size(batch.Setup.conditions.durations,1),'uni',0); end
                        if ~iscell(batch.Setup.conditions.durations{nconditions}), batch.Setup.conditions.durations{nconditions}=arrayfun(@(n)shiftdim(batch.Setup.conditions.durations{nconditions}(n,:,:),1),1:size(batch.Setup.conditions.durations{nconditions},1),'uni',0); end
                        if ~iscell(batch.Setup.conditions.durations{nconditions}{isub}), batch.Setup.conditions.durations{nconditions}{isub}=arrayfun(@(n)shiftdim(batch.Setup.conditions.durations{nconditions}{isub}(n,:,:),1),1:size(batch.Setup.conditions.durations{nconditions}{isub},1),'uni',0); end
                        t1=batch.Setup.conditions.onsets{nconditions}{isub}{nses};
                        t2=batch.Setup.conditions.durations{nconditions}{isub}{nses};
                        if isempty(t1)||all(~isfinite(t1)), t1=[]; t2=[];
                        elseif isempty(t2), t2=repmat(inf,size(t1));
                        elseif any(~isfinite(t2)), t2(~isfinite(t2))=inf;
                        end
                        %batch.Setup.conditions.onsets{nconditions}{isub}{nses}=reshape(t1,1,[]);`
                        %batch.Setup.conditions.durations{nconditions}{isub}{nses}=reshape(t2,1,[]);
                        CONN_x.Setup.conditions.values{nsub}{nl0+nconditions}{nses}={t1,t2};
                    else
                        CONN_x.Setup.conditions.values{nsub}{nl0+nconditions}{nses}={0,inf};
                    end
                end
            end
        end
        if isfield(batch.Setup.conditions,'model')
            CONN_x.Setup.conditions.model(nl0+(1:nconditions))=batch.Setup.conditions.model;
        else
            CONN_x.Setup.conditions.model(nl0+(1:nconditions))=cell(1,nconditions);
        end
        if isfield(batch.Setup.conditions,'param')
            CONN_x.Setup.conditions.param(nl0+(1:nconditions))=batch.Setup.conditions.param;
        else
            CONN_x.Setup.conditions.param(nl0+(1:nconditions))=zeros(1,nconditions);
        end
        if isfield(batch.Setup.conditions,'filter')
            CONN_x.Setup.conditions.filter(nl0+(1:nconditions))=batch.Setup.conditions.filter;
        else
            CONN_x.Setup.conditions.filter(nl0+(1:nconditions))=cell(1,nconditions);
        end
    end
    if isfield(batch.Setup,'conditions')&&isfield(batch.Setup.conditions,'missingdata'), CONN_x.Setup.conditions.missingdata=batch.Setup.conditions.missingdata; end
    
    if isfield(batch.Setup,'covariates')&&~isempty(batch.Setup.covariates),
        if isfield(batch.Setup.covariates,'add')&&batch.Setup.covariates.add, nl0=numel(CONN_x.Setup.l1covariates.names)-1; 
        else nl0=0;
        end
        ncovariates=numel(batch.Setup.covariates.names);
        CONN_x.Setup.l1covariates.names(nl0+(1:ncovariates+1))={batch.Setup.covariates.names{:},' '};
        for isub=1:numel(SUBJECTS),
            nsub=SUBJECTS(isub);
            for nl1covariates=1:length(batch.Setup.covariates.files),
                for nses=1:CONN_x.Setup.nsessions(nsub),
                    if ~iscell(batch.Setup.covariates.files{nl1covariates}), batch.Setup.covariates.files{nl1covariates}={batch.Setup.covariates.files{nl1covariates}}; end
                    if ~iscell(batch.Setup.covariates.files{nl1covariates}{isub}), batch.Setup.covariates.files{nl1covariates}{isub}={batch.Setup.covariates.files{nl1covariates}{isub}}; end
                    CONN_x.Setup.l1covariates.files{nsub}{nl0+nl1covariates}{nses}=conn_file(batch.Setup.covariates.files{nl1covariates}{isub}{nses});
                end
            end
        end
    end
    if isfield(batch.Setup,'subjects')&&~isempty(batch.Setup.subjects),
        if ~isfield(batch.Setup.subjects,'add')||batch.Setup.subjects.add==0
            CONN_x.Setup.l2covariates.names={' '};
            CONN_x.Setup.l2covariates.descrip={};
            CONN_x.Setup.l2covariates.values=repmat({{}},[CONN_x.Setup.nsubjects,1]);
        end
        if isfield(batch.Setup.subjects,'group_names')&&~isempty(batch.Setup.subjects.group_names),
            for ngroup=1:length(batch.Setup.subjects.group_names),
                idx=strmatch(batch.Setup.subjects.group_names{ngroup},CONN_x.Setup.l2covariates.names,'exact');
                if isempty(idx),
                    nl2covariates=length(CONN_x.Setup.l2covariates.names);
                    CONN_x.Setup.l2covariates.names{nl2covariates}=batch.Setup.subjects.group_names{ngroup};
                    CONN_x.Setup.l2covariates.names{nl2covariates+1}=' ';
                else, nl2covariates=idx;end
                for isub=1:numel(SUBJECTS),
                    nsub=SUBJECTS(isub);
                    CONN_x.Setup.l2covariates.values{nsub}{nl2covariates}=(batch.Setup.subjects.groups(isub)==ngroup);
                end
                if isfield(batch.Setup.subjects,'group_descrip'), CONN_x.Setup.l2covariates.descrip{nl2covariates}=batch.Setup.subjects.group_descrip{ngroup};
                elseif isfield(batch.Setup.subjects,'descrip'), CONN_x.Setup.l2covariates.descrip{nl2covariates}=batch.Setup.subjects.descrip{ngroup};
                else CONN_x.Setup.l2covariates.descrip{nl2covariates}='';
                end
            end
        end
        if isfield(batch.Setup.subjects,'effect_names')&&~isempty(batch.Setup.subjects.effect_names),
            for neffect=1:length(batch.Setup.subjects.effect_names),
                idx=strmatch(batch.Setup.subjects.effect_names{neffect},CONN_x.Setup.l2covariates.names,'exact');
                if isempty(idx),
                    nl2covariates=length(CONN_x.Setup.l2covariates.names);
                    CONN_x.Setup.l2covariates.names{nl2covariates}=batch.Setup.subjects.effect_names{neffect};
                    CONN_x.Setup.l2covariates.names{nl2covariates+1}=' ';
                else, nl2covariates=idx;end
                for isub=1:numel(SUBJECTS),
                    nsub=SUBJECTS(isub);
                    CONN_x.Setup.l2covariates.values{nsub}{nl2covariates}=batch.Setup.subjects.effects{neffect}(isub);
                end
                if isfield(batch.Setup.subjects,'effect_descrip'), CONN_x.Setup.l2covariates.descrip{nl2covariates}=batch.Setup.subjects.effect_descrip{neffect};
                elseif isfield(batch.Setup.subjects,'descrip'), CONN_x.Setup.l2covariates.descrip{nl2covariates}=batch.Setup.subjects.descrip{neffect};
                else CONN_x.Setup.l2covariates.descrip{nl2covariates}='';
                end
            end
        end
    end
    
    if isfield(batch.Setup,'preprocessing')
        if isfield(batch.Setup.preprocessing,'steps'),%&&~isempty(batch.Setup.preprocessing.steps),
            OPTIONS={'dogui',0,'fwhm',8,'art_thresholds',[5 .9]}; %note: art_useconservative=0 uses [9 2]; art_useconservative=1 uses [5 .9]; art_useconservative=2 uses [3 .5]
            if isempty(batch.Setup.preprocessing.steps)&&~isfield(batch.Setup.preprocessing,'multiplesteps'), batch.Setup.preprocessing.multiplesteps=1; end
            steps=batch.Setup.preprocessing.steps;
            batch.Setup.preprocessing=rmfield(batch.Setup.preprocessing,'steps');
            options=fieldnames(batch.Setup.preprocessing);
            for n=1:numel(options)
                OPTIONS=[OPTIONS {options{n}, batch.Setup.preprocessing.(options{n})}];
            end
            if isfield(batch,'parallel')&&isfield(batch.parallel,'N')&&batch.parallel.N>0,
                PAR_CMD{end+1}='setup_preprocessing'; PAR_ARG{end+1}={[],steps,OPTIONS{:}};
            else
                conn_process('setup_preprocessing',steps,'subjects',SUBJECTS,OPTIONS{:});
                CONN_x.gui=1;
                if isfield(CONN_x,'filename')&&~isempty(CONN_x.filename), conn save; end
            end
        else
            pnames=fieldnames(batch.Setup.preprocessing);
            if ~isempty(pnames), conn_disp('fprintf','warning: no preprocessing steps entered. Skipping %s information\n',sprintf('%s ',pnames{:})); end
        end
        %conn_setup_preproc(steps,OPTIONS{:});
    end
    
    if isfield(batch.Setup,'done')&&batch.Setup.done,
        if ~conn_projectmanager('inserver')&&~(isfield(batch,'parallel')&&isfield(batch.parallel,'N')&&batch.parallel.N>0), conn save; end
        if ~isfield(batch.Setup,'overwrite'), batch.Setup.overwrite='Yes'; end
        if isscalar(batch.Setup.overwrite)&&~isstruct(batch.Setup.overwrite)&&ismember(double(batch.Setup.overwrite),[1 89 121]), batch.Setup.overwrite='Yes'; end
        CONN_x.gui=struct('overwrite',batch.Setup.overwrite,'subjects',SUBJECTS);
        if isfield(batch,'parallel')&&isfield(batch.parallel,'N')&&batch.parallel.N>0, 
            %CONN_x.isready(2)=1;
            PAR_CMD{end+1}='Setup'; PAR_ARG{end+1}={CONN_x.gui};
        else
            conn_process Setup;
            CONN_x.gui=1;
            conn save;
        end
    else,
        if isfield(batch,'filename'), conn save; end
    end
else
    if isfield(batch,'filename'),
        if ~isempty(batch.filename)&&ischar(batch.filename), [nill,nill,extt]=fileparts(batch.filename); if isempty(extt), batch.filename=[batch.filename,'.mat']; end; end
        CONN_x.filename=batch.filename;
        CONN_x.gui=0;
        conn load;                      % loads existing conn_* project
        CONN_x.gui=1;
    end
    if isfield(batch,'subjects'), SUBJECTS=batch.subjects;
    else SUBJECTS=1:CONN_x.Setup.nsubjects;
    end
end

%% DENOISING step
if isfield(batch,'Preprocessing')&&~isfield(batch,'Denoising'),batch.Denoising=batch.Preprocessing; end
if isfield(batch,'Denoising'),
%     if isfield(batch,'filename'),
%         CONN_x.filename=batch.filename;
%         CONN_x.gui=0;
%         conn load;                      % loads existing conn_* project
%         CONN_x.gui=1;
%     end
    if isfield(batch.Denoising,'filter')&&~isempty(batch.Denoising.filter),
        CONN_x.Preproc.filter=batch.Denoising.filter;          % frequency filter (band-pass values, in Hz)
    end
    if isfield(batch.Denoising,'detrending')&&~isempty(batch.Denoising.detrending),
        CONN_x.Preproc.detrending=batch.Denoising.detrending;          
    end
    if isfield(batch.Denoising,'despiking')&&~isempty(batch.Denoising.despiking),
        CONN_x.Preproc.despiking=batch.Denoising.despiking;          
    end
    if isfield(batch.Denoising,'regbp')&&~isempty(batch.Denoising.regbp),
        CONN_x.Preproc.regbp=batch.Denoising.regbp;          
    end
    if isfield(batch.Denoising,'confounds')
        if isempty(batch.Denoising.confounds),
            batch.Denoising.confounds=struct;
            batch.Denoising.confounds.names={};
        else
            if ~isstruct(batch.Denoising.confounds),
                temp=batch.Denoising.confounds;
                batch.Denoising.confounds=struct;
                batch.Denoising.confounds.names=cellstr(temp);
            end
        end
        CONN_x.Preproc.confounds.names=batch.Denoising.confounds.names;
        if isfield(batch.Denoising.confounds,'dimensions')&&~isempty(batch.Denoising.confounds.dimensions), CONN_x.Preproc.confounds.dimensions=batch.Denoising.confounds.dimensions; else CONN_x.Preproc.confounds.dimensions={}; end
        if isfield(batch.Denoising.confounds,'deriv')&&~isempty(batch.Denoising.confounds.deriv), CONN_x.Preproc.confounds.deriv=batch.Denoising.confounds.deriv; else CONN_x.Preproc.confounds.deriv={}; end
        if isfield(batch.Denoising.confounds,'power')&&~isempty(batch.Denoising.confounds.power), CONN_x.Preproc.confounds.power=batch.Denoising.confounds.power; else CONN_x.Preproc.confounds.power={}; end
        if isfield(batch.Denoising.confounds,'filter')&&~isempty(batch.Denoising.confounds.filter), CONN_x.Preproc.confounds.filter=batch.Denoising.confounds.filter; else CONN_x.Preproc.confounds.filter={}; end
    end
    
    if isfield(batch.Denoising,'done')&&batch.Denoising.done,
        if ~conn_projectmanager('inserver')&&~(isfield(batch,'parallel')&&isfield(batch.parallel,'N')&&batch.parallel.N>0), conn save; end
        if ~isfield(batch.Denoising,'overwrite'), batch.Denoising.overwrite='Yes'; end
        if isscalar(batch.Denoising.overwrite)&&~isstruct(batch.Denoising.overwrite)&&ismember(double(batch.Denoising.overwrite),[1 89 121]), batch.Denoising.overwrite='Yes'; end
        CONN_x.gui=struct('overwrite',batch.Denoising.overwrite,'subjects',SUBJECTS);
        if isfield(batch,'parallel')&&isfield(batch.parallel,'N')&&batch.parallel.N>0, 
            PAR_CMD{end+1}='Denoising'; PAR_ARG{end+1}={CONN_x.gui};
        else
            conn_process Denoising;
            CONN_x.gui=1;
            conn save;
        end
    else
        if isfield(batch.Denoising,'confounds')&&~isempty(batch.Denoising.confounds), conn_process setup_updatedenoising; end
        if isfield(batch,'filename'), conn save; end
    end
end

%% ANALYSIS step
if isfield(batch,'Analysis')&&isfield(batch.Analysis,'measures')&&~isfield(batch,'vvAnalysis'), batch.vvAnalysis=batch.Analysis; batch=rmfield(batch,'Analysis'); end

if isfield(batch,'Analysis'),
%     if isfield(batch,'filename'),
%         CONN_x.filename=batch.filename;
%         CONN_x.gui=0;
%         conn load;                      % loads existing conn_* project
%         CONN_x.gui=1;
%     end
    if isfield(batch.Analysis,'sources')||~isfield(batch.Analysis,'measures'),
        if isfield(batch.Analysis,'name'), batch.Analysis.analysis_number=batch.Analysis.name; end
        if ~isfield(batch.Analysis,'analysis_number')||isempty(batch.Analysis.analysis_number),batch.Analysis.analysis_number=1; end
        if ~isfield(batch.Analysis,'modulation')||isempty(batch.Analysis.modulation),batch.Analysis.modulation=0; end
        if ~isfield(batch.Analysis,'measure')||isempty(batch.Analysis.measure),batch.Analysis.measure=1; end
        if ~isfield(batch.Analysis,'weight')||isempty(batch.Analysis.weight),batch.Analysis.weight=2; end
        if ~isfield(batch.Analysis,'type')||isempty(batch.Analysis.type),batch.Analysis.type=3; end
        if ~isfield(batch.Analysis,'conditions'),batch.Analysis.conditions=[]; end
        if iscell(batch.Analysis.analysis_number), batch.Analysis.analysis_number=char(batch.Analysis.analysis_number); end
        if ischar(batch.Analysis.analysis_number)
            ianalysis=strmatch(batch.Analysis.analysis_number,{CONN_x.Analyses.name},'exact');
            if isempty(ianalysis), 
                ianalysis=numel(CONN_x.Analyses)+1;
                CONN_x.Analyses(ianalysis)=struct(...
                 'name','SBC_01',...
                 'sourcenames',{{}},...
                 'regressors',	struct('names',{{}},'types',{{}},'deriv',{{}},'fbands',{{}},'dimensions',{{}}),...
                 'type',3,...
                 'measure',1,...
                 'modulation',0,...
                 'conditions',[],...
                 'weight',2,...
                 'sources',{{}});
                %CONN_x.Analyses(ianalysis)=CONN_x.Analyses(end);
                CONN_x.Analyses(ianalysis).name=batch.Analysis.analysis_number; 
            end
            batch.Analysis.analysis_number=ianalysis;
        end
        if ischar(batch.Analysis.weight)
            temp=strmatch(lower(batch.Analysis.weight),{'none','hrf','hanning'});
            if isempty(temp), error('unrecognized Analysis.weight %s',batch.Analysis.weight); end
            batch.Analysis.weight=temp;
        end
        if ischar(batch.Analysis.measure)
            temp=strmatch(lower(batch.Analysis.measure),{'correlation (bivariate)','correlation (semipartial)','regression (bivariate)','regression (multivariate)'});
            if isempty(temp), error('unrecognized Analysis.measure %s',batch.Analysis.measure); end
            batch.Analysis.measure=temp;
        end
        if ischar(batch.Analysis.type)
            temp=strmatch(lower(batch.Analysis.type),{'roi-to-roi', 'seed-to-voxel', 'all'});
            if isempty(temp), error('unrecognized Analysis.type %s',batch.Analysis.type); end
            batch.Analysis.type=temp;
        end
        CONN_x.Analysis=batch.Analysis.analysis_number;
        CONN_x.Analyses(CONN_x.Analysis).modulation=batch.Analysis.modulation;
        CONN_x.Analyses(CONN_x.Analysis).measure=batch.Analysis.measure;
        CONN_x.Analyses(CONN_x.Analysis).weight=batch.Analysis.weight;
        CONN_x.Analyses(CONN_x.Analysis).type=batch.Analysis.type;
        CONN_x.Analyses(CONN_x.Analysis).conditions=batch.Analysis.conditions;
        if ~isfield(batch.Analysis,'sources')||isempty(batch.Analysis.sources),
            CONN_x.Analyses(CONN_x.Analysis).regressors.names={};
        else
            batch.Analysis.sources
            if ~isstruct(batch.Analysis.sources),
                CONN_x.Analyses(CONN_x.Analysis).regressors.names=batch.Analysis.sources;
                CONN_x.Analyses(CONN_x.Analysis).regressors.dimensions=repmat({1},size(batch.Analysis.sources));
                CONN_x.Analyses(CONN_x.Analysis).regressors.deriv=repmat({0},size(batch.Analysis.sources));
                CONN_x.Analyses(CONN_x.Analysis).regressors.types=repmat({'roi'},size(batch.Analysis.sources));
                CONN_x.Analyses(CONN_x.Analysis).regressors.fbands=repmat({1},size(batch.Analysis.sources));
            else
                CONN_x.Analyses(CONN_x.Analysis).regressors.names=batch.Analysis.sources.names;
                CONN_x.Analyses(CONN_x.Analysis).regressors.dimensions=batch.Analysis.sources.dimensions;
                CONN_x.Analyses(CONN_x.Analysis).regressors.deriv=batch.Analysis.sources.deriv;
                CONN_x.Analyses(CONN_x.Analysis).regressors.types=repmat({'roi'},size(batch.Analysis.sources.names));
                if isfield(batch.Analysis.sources,'fbands')
                    CONN_x.Analyses(CONN_x.Analysis).regressors.fbands=batch.Analysis.sources.fbands;
                else
                    CONN_x.Analyses(CONN_x.Analysis).regressors.fbands=repmat({1},1,numel(batch.Analysis.sources.names));
                end
            end
        end
    end
    
    if isfield(batch.Analysis,'done')&&batch.Analysis.done,
        if ~conn_projectmanager('inserver')&&~(isfield(batch,'parallel')&&isfield(batch.parallel,'N')&&batch.parallel.N>0), conn save; end
        if ~isfield(batch.Analysis,'overwrite'), batch.Analysis.overwrite='Yes'; end
        if isscalar(batch.Analysis.overwrite)&&~isstruct(batch.Analysis.overwrite)&&ismember(double(batch.Analysis.overwrite),[1 89 121]), batch.Analysis.overwrite='Yes'; end
        CONN_x.gui=struct('overwrite',batch.Analysis.overwrite,'subjects',SUBJECTS);
        if isfield(batch,'parallel')&&isfield(batch.parallel,'N')&&batch.parallel.N>0, 
            PAR_CMD{end+1}='Analyses_seedandroi'; 
            if ~isfield(batch.Analysis,'analysis_number'), PAR_ARG{end+1}={CONN_x.gui}; 
            else PAR_ARG{end+1}={CONN_x.gui, batch.Analysis.analysis_number}; 
            end
        else
            if ~isfield(batch.Analysis,'analysis_number'), conn_process('Analyses_seedandroi');
            else conn_process('Analyses_seedandroi',batch.Analysis.analysis_number);
            end
            CONN_x.gui=1;
            conn save;
        end
%     elseif isfield(batch.Analysis,'sources')||~isfield(batch.Analysis,'measures'),
%         if ~isfield(batch.Analysis,'analysis_number'), conn_process('analyses_seedsetup');
%         else conn_process('analyses_seedsetup',batch.Analysis.analysis_number);
%         end
    else
        if isfield(batch,'filename'), conn save; end
    end
end

%% ANALYSIS step
if isfield(batch,'vvAnalysis'),
%     if isfield(batch,'filename'),
%         CONN_x.filename=batch.filename;
%         CONN_x.gui=0;
%         conn load;                      % loads existing conn_* project
%         CONN_x.gui=1;
%     end
    if isfield(batch.vvAnalysis,'measures')
        if isfield(batch.vvAnalysis,'name'), batch.vvAnalysis.analysis_number=batch.vvAnalysis.name; end
        if ~isfield(batch.vvAnalysis,'analysis_number')||isempty(batch.vvAnalysis.analysis_number),batch.vvAnalysis.analysis_number=1; end
        if iscell(batch.vvAnalysis.analysis_number), batch.vvAnalysis.analysis_number=char(batch.vvAnalysis.analysis_number); end
        if ischar(batch.vvAnalysis.analysis_number)
            ianalysis=strmatch(batch.vvAnalysis.analysis_number,{CONN_x.vvAnalyses.name},'exact');
            if isempty(ianalysis), 
                ianalysis=numel(CONN_x.vvAnalyses)+1;
                CONN_x.vvAnalyses(ianalysis)=struct(...
                 'name','V2V_01',...
                 'measurenames',{{}},...
                 'variables',  conn_v2v('measures'),...
                 'regressors', conn_v2v('empty'),...
                 'measures',{{}},...
                 'options','',...
                 'mask',[]);
%                 CONN_x.vvAnalyses(ianalysis)=CONN_x.vvAnalyses(end);
                CONN_x.vvAnalyses(ianalysis).name=batch.vvAnalysis.analysis_number; 
            end
            batch.vvAnalysis.analysis_number=ianalysis;
        end
        CONN_x.vvAnalysis=batch.vvAnalysis.analysis_number;
        if isempty(batch.vvAnalysis.measures),
            CONN_x.vvAnalyses(CONN_x.vvAnalysis).regressors.names={};
        elseif ~isstruct(batch.vvAnalysis.measures),
            batch.vvAnalysis.measures=cellstr(batch.vvAnalysis.measures);
            CONN_x.vvAnalyses(CONN_x.vvAnalysis).regressors.names=batch.vvAnalysis.measures;
            CONN_x.vvAnalyses(CONN_x.vvAnalysis).regressors.measuretype=repmat({[]},size(batch.vvAnalysis.measures));
            CONN_x.vvAnalyses(CONN_x.vvAnalysis).regressors.global=repmat({[]},size(batch.vvAnalysis.measures));
            CONN_x.vvAnalyses(CONN_x.vvAnalysis).regressors.localsupport=repmat({[]},size(batch.vvAnalysis.measures));
            CONN_x.vvAnalyses(CONN_x.vvAnalysis).regressors.deriv=repmat({[]},size(batch.vvAnalysis.measures));
            CONN_x.vvAnalyses(CONN_x.vvAnalysis).regressors.norm=repmat({[]},size(batch.vvAnalysis.measures));
            CONN_x.vvAnalyses(CONN_x.vvAnalysis).regressors.filename=repmat({''},size(batch.vvAnalysis.measures));
            CONN_x.vvAnalyses(CONN_x.vvAnalysis).regressors.dimensions_in=repmat({[]},size(batch.vvAnalysis.measures));
            CONN_x.vvAnalyses(CONN_x.vvAnalysis).regressors.dimensions_out=repmat({[]},size(batch.vvAnalysis.measures));
            CONN_x.vvAnalyses(CONN_x.vvAnalysis).options='';
            CONN_x.vvAnalyses(CONN_x.vvAnalysis).mask=[];
        else
            if ~isfield(batch.vvAnalysis.measures,'measureclass'),batch.vvAnalysis.measures.measureclass=repmat({[]},size(batch.vvAnalysis.measures.names)); end
            if ~isfield(batch.vvAnalysis.measures,'type'),batch.vvAnalysis.measures.type=repmat({[]},size(batch.vvAnalysis.measures.names)); end
            if ~isfield(batch.vvAnalysis.measures,'kernelsupport'),batch.vvAnalysis.measures.kernelsupport=repmat({[]},size(batch.vvAnalysis.measures.names)); end
            if ~isfield(batch.vvAnalysis.measures,'kernelshape'),batch.vvAnalysis.measures.kernelshape=repmat({[]},size(batch.vvAnalysis.measures.names)); end
            if ~isfield(batch.vvAnalysis.measures,'dimensions'),batch.vvAnalysis.measures.dimensions=repmat({[]},size(batch.vvAnalysis.measures.names)); end
            if ~isfield(batch.vvAnalysis.measures,'factors'),batch.vvAnalysis.measures.factors=repmat({[]},size(batch.vvAnalysis.measures.names)); end
            if ~isfield(batch.vvAnalysis.measures,'norm'),batch.vvAnalysis.measures.norm=repmat({[]},size(batch.vvAnalysis.measures.names)); end
            if ~isfield(batch.vvAnalysis.measures,'dimensions_in'),batch.vvAnalysis.measures.dimensions_in=batch.vvAnalysis.measures.dimensions; end
            if ~isfield(batch.vvAnalysis.measures,'dimensions_out'),batch.vvAnalysis.measures.dimensions_out=batch.vvAnalysis.measures.factors; end
            if ~isfield(batch.vvAnalysis.measures,'options'),batch.vvAnalysis.measures.options=''; end
            if ~isfield(batch.vvAnalysis.measures,'mask'),batch.vvAnalysis.measures.mask=[]; end
            if ~iscell(batch.vvAnalysis.measures.names), batch.vvAnalysis.measures.names={batch.vvAnalysis.measures.names}; end
            if ~iscell(batch.vvAnalysis.measures.measureclass), batch.vvAnalysis.measures.measureclass={batch.vvAnalysis.measures.measureclass}; end
            if ~iscell(batch.vvAnalysis.measures.kernelsupport), batch.vvAnalysis.measures.kernelsupport={batch.vvAnalysis.measures.kernelsupport}; end
            if ~iscell(batch.vvAnalysis.measures.type), batch.vvAnalysis.measures.type={batch.vvAnalysis.measures.type}; end
            if ~iscell(batch.vvAnalysis.measures.kernelshape), batch.vvAnalysis.measures.kernelshape={batch.vvAnalysis.measures.kernelshape}; end
            if ~iscell(batch.vvAnalysis.measures.norm), batch.vvAnalysis.measures.norm={batch.vvAnalysis.measures.norm}; end
            if ~iscell(batch.vvAnalysis.measures.dimensions_in), batch.vvAnalysis.measures.dimensions_in={batch.vvAnalysis.measures.dimensions_in}; end
            if ~iscell(batch.vvAnalysis.measures.dimensions_out), batch.vvAnalysis.measures.dimensions_out={batch.vvAnalysis.measures.dimensions_out}; end
            CONN_x.vvAnalyses(CONN_x.vvAnalysis).regressors.names=batch.vvAnalysis.measures.names;
            CONN_x.vvAnalyses(CONN_x.vvAnalysis).regressors.measuretype=batch.vvAnalysis.measures.measureclass;
            CONN_x.vvAnalyses(CONN_x.vvAnalysis).regressors.localsupport=batch.vvAnalysis.measures.kernelsupport;
            CONN_x.vvAnalyses(CONN_x.vvAnalysis).regressors.global=batch.vvAnalysis.measures.type;
            CONN_x.vvAnalyses(CONN_x.vvAnalysis).regressors.deriv=batch.vvAnalysis.measures.kernelshape;
            CONN_x.vvAnalyses(CONN_x.vvAnalysis).regressors.norm=batch.vvAnalysis.measures.norm;
            CONN_x.vvAnalyses(CONN_x.vvAnalysis).regressors.filename=repmat({''},size(batch.vvAnalysis.measures.names));
            CONN_x.vvAnalyses(CONN_x.vvAnalysis).regressors.dimensions_in=batch.vvAnalysis.measures.dimensions_in;
            CONN_x.vvAnalyses(CONN_x.vvAnalysis).regressors.dimensions_out=batch.vvAnalysis.measures.dimensions_out;
            CONN_x.vvAnalyses(CONN_x.vvAnalysis).options=batch.vvAnalysis.measures.options;
            if isempty(batch.vvAnalysis.measures.mask), CONN_x.vvAnalyses(CONN_x.vvAnalysis).mask=batch.vvAnalysis.measures.mask;
            else CONN_x.vvAnalyses(CONN_x.vvAnalysis).mask=conn_file(batch.vvAnalysis.measures.mask);
            end
        end
    end
    
    if isfield(batch.vvAnalysis,'done')&&batch.vvAnalysis.done,
        if ~conn_projectmanager('inserver')&&~(isfield(batch,'parallel')&&isfield(batch.parallel,'N')&&batch.parallel.N>0), conn save; end
        if ~isfield(batch.vvAnalysis,'overwrite'), batch.vvAnalysis.overwrite='Yes'; end
        if isscalar(batch.vvAnalysis.overwrite)&&~isstruct(batch.vvAnalysis.overwrite)&&ismember(double(batch.vvAnalysis.overwrite),[1 89 121]), batch.vvAnalysis.overwrite='Yes'; end
        CONN_x.gui=struct('overwrite',batch.vvAnalysis.overwrite,'subjects',SUBJECTS);
        if isfield(batch,'parallel')&&isfield(batch.parallel,'N')&&batch.parallel.N>0, 
            PAR_CMD{end+1}='Analyses_vv'; 
            if ~isfield(batch.vvAnalysis,'analysis_number'), PAR_ARG{end+1}={CONN_x.gui}; 
            else PAR_ARG{end+1}={CONN_x.gui, batch.vvAnalysis.analysis_number}; 
            end
        else
            if ~isfield(batch.vvAnalysis,'analysis_number'), conn_process('Analyses_vv');
            else conn_process('Analyses_vv',batch.vvAnalysis.analysis_number);
            end
            CONN_x.gui=1;
            conn save;
        end
%     elseif isfield(batch.vvAnalysis,'measures')
%         if ~isfield(batch.vvAnalysis,'analysis_number'), conn_process('Analyses_vvsetup');
%         else conn_process('Analyses_vvsetup',batch.vvAnalysis.analysis_number);
%         end
    else
        if isfield(batch,'filename'), conn save; end
    end
end

%% ANALYSIS step
if isfield(batch,'dynAnalysis'),
%     if isfield(batch,'filename'),
%         CONN_x.filename=batch.filename;
%         CONN_x.gui=0;
%         conn load;                      % loads existing conn_* project
%         CONN_x.gui=1;
%     end
    if isfield(batch.dynAnalysis,'name'), batch.dynAnalysis.analysis_number=batch.dynAnalysis.name; end
    if ~isfield(batch.dynAnalysis,'analysis_number')||isempty(batch.dynAnalysis.analysis_number),batch.dynAnalysis.analysis_number=1; end
    if iscell(batch.dynAnalysis.analysis_number), batch.dynAnalysis.analysis_number=char(batch.dynAnalysis.analysis_number); end
    if ischar(batch.dynAnalysis.analysis_number)
        ianalysis=strmatch(batch.dynAnalysis.analysis_number,{CONN_x.dynAnalyses.name},'exact');
        if isempty(ianalysis),
            ianalysis=numel(CONN_x.dynAnalyses)+1;
            CONN_x.dynAnalyses(ianalysis)=struct(...
                 'name','DYN_01',...
                 'regressors', struct('names',{{}}),...
                 'variables', struct('names',{{}}),...
                 'Ncomponents',20,...
                 'condition',[],...
                 'analyses',3,...
                 'window',30,...
                 'output',[1 1 0]);
            %CONN_x.dynAnalyses(ianalysis)=CONN_x.dynAnalyses(end);
            CONN_x.dynAnalyses(ianalysis).name=batch.dynAnalysis.analysis_number;
        end
        batch.dynAnalysis.analysis_number=ianalysis;
    end
    CONN_x.dynAnalysis=batch.dynAnalysis.analysis_number;

    if isfield(batch.dynAnalysis,'sources'), CONN_x.dynAnalyses(CONN_x.dynAnalysis).regressors.names=batch.dynAnalysis.sources; end
    if isfield(batch.dynAnalysis,'factors'), CONN_x.dynAnalyses(CONN_x.dynAnalysis).Ncomponents=batch.dynAnalysis.factors; end
    if isfield(batch.dynAnalysis,'window'), CONN_x.dynAnalyses(CONN_x.dynAnalysis).window=batch.dynAnalysis.window; end
    
    if isfield(batch.dynAnalysis,'condition'), CONN_x.dynAnalyses(CONN_x.dynAnalysis).condition=batch.dynAnalysis.condition; end
    if isfield(batch.dynAnalysis,'done_step1only')&&batch.DynAnalysis.done_step1, stepname='Analyses_dyn_step1'; batch.dynAnalysis.done=true;
    elseif isfield(batch.dynAnalysis,'done_step2only')&&batch.DynAnalysis.done_step1, stepname='Analyses_dyn_step2'; batch.dynAnalysis.done=true;
    else stepname='Analyses_dyn';
    end
    if isfield(batch.dynAnalysis,'done')&&batch.dynAnalysis.done,
        if ~conn_projectmanager('inserver')&&~(isfield(batch,'parallel')&&isfield(batch.parallel,'N')&&batch.parallel.N>0), conn save; end
        if ~isfield(batch.dynAnalysis,'overwrite'), batch.dynAnalysis.overwrite='Yes'; end
        if isscalar(batch.dynAnalysis.overwrite)&&~isstruct(batch.dynAnalysis.overwrite)&&ismember(double(batch.dynAnalysis.overwrite),[1 89 121]), batch.dynAnalysis.overwrite='Yes'; end
        CONN_x.gui=struct('overwrite',batch.dynAnalysis.overwrite,'subjects',SUBJECTS);
        if isfield(batch,'parallel')&&isfield(batch.parallel,'N')&&batch.parallel.N>0, 
            PAR_CMD{end+1}=stepname; 
            if ~isfield(batch.dynAnalysis,'analysis_number'), PAR_ARG{end+1}={CONN_x.gui}; 
            else PAR_ARG{end+1}={CONN_x.gui, batch.dynAnalysis.analysis_number}; 
            end
        else
            if ~isfield(batch.dynAnalysis,'analysis_number'), conn_process(stepname);
            else conn_process(stepname,batch.dynAnalysis.analysis_number);
            end
            CONN_x.gui=1;
            conn save;
        end
%     elseif isfield(batch.dynAnalysis,'sources')
%         if ~isfield(batch.dynAnalysis,'analysis_number'), conn_process('Analyses_dynsetup');
%         else conn_process('Analyses_dynsetup',batch.dynAnalysis.analysis_number);
%         end
    else
        if isfield(batch,'filename'), conn save; end
    end
end

if isfield(batch,'QA'),
    if isfield(batch.QA,'foldername'), qafolder=batch.QA.foldername;
    else qafolder=fullfile(CONN_x.folders.qa,['QA_',datestr(now,'yyyy_mm_dd_HHMMSSFFF')]); ;
    end
    if isfield(batch.QA,'plots'), procedures=batch.QA.plots; 
    else procedures=[];
    end
    if isfield(batch.QA,'rois'), validrois=batch.QA.rois; 
    else validrois=[];
    end
    if isfield(batch.QA,'sets'), validsets=batch.QA.sets; 
    else validsets=[];
    end
    if isfield(batch.QA,'l2covariates'), nl2covariates=batch.QA.l2covariates; 
    else nl2covariates=[];
    end
    if isfield(batch.QA,'l1contrasts'), nl1contrasts=batch.QA.l1contrasts; 
    else nl1contrasts=[];
    end
    if isfield(batch.QA,'conditions'), validconditions=batch.QA.conditions; 
    else validconditions=[];
    end
    if isfield(batch,'parallel')&&isfield(batch.parallel,'N')&&batch.parallel.N>0,
        PAR_CMD{end+1}='qaplots';
        PAR_ARG{end+1}={[],qafolder,procedures,SUBJECTS,validrois,validsets,nl2covariates,nl1contrasts,validconditions};
    else
        conn_process('qaplots',qafolder,procedures,SUBJECTS,validrois,validsets,nl2covariates,nl1contrasts,validconditions);
        CONN_x.gui=1;
    end
end

% parallel submission
if isfield(batch,'parallel')&&isfield(batch.parallel,'N')&&batch.parallel.N>0
    conn save;
    if isfield(batch.parallel,'profile'), conn_jobmanager('setprofile',batch.parallel.profile);
    else conn_jobmanager('setprofile'); % if not specified, use default
    end
    if ~isempty(PAR_CMD),
        pnames=fieldnames(batch.parallel);
        pnames=pnames(~ismember(pnames,{'N','profile','immediatereturn'}));
        for n=1:numel(pnames), conn_jobmanager('options',pnames{n},batch.parallel.(pnames{n})); end
        info=conn_jobmanager('submit',PAR_CMD,SUBJECTS,batch.parallel.N,PAR_ARG);
        if ~isfield(batch.parallel,'immediatereturn')||~batch.parallel.immediatereturn, conn_jobmanager('waitfor',info); end
    end
    if isfield(batch,'Results')||isfield(batch,'vvResults'), 
        clear batchtemp;
        if isfield(batch,'Results'), batchtemp.Results=batch.Results;
        else batchtemp.vvResults=batch.vvResults;
        end
        if isfield(batch,'filename'), batchtemp.filename=batch.filename; else, batchtemp.filename=CONN_x.filename; end
        info=conn_jobmanager('submit',{'batch'},[],1,{{[],batchtemp}});
        if ~isfield(batch.parallel,'immediatereturn')||~batch.parallel.immediatereturn, conn_jobmanager('waitfor',info); end
        return
    end
end

%% RESULTS step
isvv=false; 
if isfield(batch,'Results')&&isfield(batch.Results,'between_measures'), isvv=true; end
if isfield(batch,'vvResults'), isvv=true; batch.Results=batch.vvResults; batch=rmfield(batch,'vvResults'); end
if isfield(batch,'Results'),
    if isfield(batch,'filename'),
        CONN_x.filename=batch.filename;
        CONN_x.gui=0;
        conn load;                      % loads existing conn_* project
        CONN_x.gui=1;
    end
    if isfield(batch.Results,'name'), batch.Results.analysis_number=batch.Results.name; end
    if ~isfield(batch.Results,'analysis_number')||isempty(batch.Results.analysis_number), batch.Results.analysis_number=1; end
    if ischar(batch.Results.analysis_number)
        if isvv, ianalysis=strmatch(batch.Results.analysis_number,{CONN_x.vvAnalyses.name},'exact');
        else
            ianalysis=strmatch(batch.Results.analysis_number,{CONN_x.Analyses.name},'exact');
            if isempty(ianalysis), isvv=true; ianalysis=strmatch(batch.Results.analysis_number,{CONN_x.vvAnalyses.name},'exact'); end
        end
        if isempty(ianalysis), error('unrecognized analysis %s',batch.Results.analysis_number); end
        batch.Results.analysis_number=ianalysis;
    end
    if isvv&&isfield(batch.Results,'between_sources')&&~isfield(batch.Results,'between_measures'), batch.Results.between_measures=batch.Results.between_sources; batch.Results=rmfield(batch.Results,'between_sources'); end
    if isvv, CONN_x.vvAnalysis=batch.Results.analysis_number;
    else CONN_x.Analysis=batch.Results.analysis_number;
    end
    if isfield(batch.Results,'foldername'),CONN_x.Results.foldername=batch.Results.foldername;else CONN_x.Results.foldername=''; end

    if isfield(batch.Results,'between_subjects')&&~isempty(batch.Results.between_subjects),
        for neffect=1:length(batch.Results.between_subjects.effect_names),
            idx=strmatch(batch.Results.between_subjects.effect_names{neffect},CONN_x.Setup.l2covariates.names,'exact');
            if isempty(idx), 
                if isfield(batch.Results.between_subjects,'effects')
                    nl2covariates=length(CONN_x.Setup.l2covariates.names);
                    CONN_x.Setup.l2covariates.names{nl2covariates}=batch.Results.between_subjects.effect_names{neffect};
                    CONN_x.Setup.l2covariates.names{nl2covariates+1}=' ';
                    for nsub=1:CONN_x.Setup.nsubjects,
                        CONN_x.Setup.l2covariates.values{nsub}{nl2covariates}=batch.Results.between_subjects.effects{neffect}(nsub);
                    end
                    CONN_x.Setup.l2covariates.descrip{nl2covariates}='';
                elseif any(batch.Results.between_subjects.effect_names{neffect}=='*')&&all(ismember(regexp(batch.Results.between_subjects.effect_names{neffect},'\s*\*\s*','split'),CONN_x.Setup.l2covariates.names(1:end-1)))
                    [nill,idx]=ismember(regexp(batch.Results.between_subjects.effect_names{neffect},'\s*\*\s*','split'),CONN_x.Setup.l2covariates.names(1:end-1));
                    nl2covariates=length(CONN_x.Setup.l2covariates.names);
                    CONN_x.Setup.l2covariates.names{nl2covariates}=batch.Results.between_subjects.effect_names{neffect};
                    CONN_x.Setup.l2covariates.names{nl2covariates+1}=' ';
                    for nsub=1:CONN_x.Setup.nsubjects,
                        CONN_x.Setup.l2covariates.values{nsub}{nl2covariates}=prod([CONN_x.Setup.l2covariates.values{nsub}{idx}]);
                    end
                    CONN_x.Setup.l2covariates.descrip{nl2covariates}='';
                else
                    error(['unknown subject effect ',batch.Results.between_subjects.effect_names{neffect}]); return;
                end
            end
        end
        CONN_x.Results.xX.nsubjecteffects=zeros(1,length(batch.Results.between_subjects.effect_names));
        for neffect=1:length(batch.Results.between_subjects.effect_names),
            idx=strmatch(batch.Results.between_subjects.effect_names{neffect},CONN_x.Setup.l2covariates.names,'exact');
            if isempty(idx), error(['unknown subject effect ',batch.Results.between_subjects.effect_names{neffect}]); return;
            else, CONN_x.Results.xX.nsubjecteffects(neffect)=idx(1); end
        end
        assert(isfield(batch.Results.between_subjects,'contrast'),'missing field batch.Results.between_subjects.contrast');
        assert(size(batch.Results.between_subjects.contrast,2)==length(batch.Results.between_subjects.effect_names),'number of columns in between_subjects.contrast (%d) must match number of elements in between_subjects.effect_names (%d)',size(batch.Results.between_subjects.contrast,2),length(batch.Results.between_subjects.effect_names));
        [CONN_x.Results.xX.nsubjecteffects,idx]=sort(CONN_x.Results.xX.nsubjecteffects); % note: resorts list of subject effects (bug fix Gregor Lichtner)
        CONN_x.Results.xX.csubjecteffects=batch.Results.between_subjects.contrast(:,idx);
        CONN_x.Results.xX.nsubjecteffectsbyname=CONN_x.Setup.l2covariates.names(CONN_x.Results.xX.nsubjecteffects);
       
        if ~isfield(batch.Results,'between_conditions')||isempty(batch.Results.between_conditions),
            if 1, %isfield(batch.Results,'done')&&batch.Results.done
                clear batchtemp;
                if isfield(batch,'filename'), batchtemp.filename=batch.filename; else, batchtemp.filename=CONN_x.filename; end
                if isvv
                    batchtemp.vvResults=batch.Results;
                    for ncondition=1:length(CONN_x.Setup.conditions.names)-1,
                        batchtemp.vvResults.between_conditions.effect_names={CONN_x.Setup.conditions.names{ncondition}};
                        batchtemp.vvResults.between_conditions.contrast=[1];
                        conn_batch(batchtemp);
                    end
                else
                    batchtemp.Results=batch.Results;
                    for ncondition=1:length(CONN_x.Setup.conditions.names)-1,
                        batchtemp.Results.between_conditions.effect_names={CONN_x.Setup.conditions.names{ncondition}};
                        batchtemp.Results.between_conditions.contrast=[1];
                        conn_batch(batchtemp);
                    end
                end
            end
        else
            CONN_x.Results.xX.nconditions=zeros(1,length(batch.Results.between_conditions.effect_names));
            for neffect=1:length(batch.Results.between_conditions.effect_names),
                idx=strmatch(batch.Results.between_conditions.effect_names{neffect},CONN_x.Setup.conditions.names,'exact');
                if isempty(idx), error(['unknown condition ',batch.Results.between_conditions.effect_names{neffect}]); return;
                else, CONN_x.Results.xX.nconditions(neffect)=idx(1); end
            end
            CONN_x.Results.xX.cconditions=batch.Results.between_conditions.contrast;
            CONN_x.Results.xX.nconditionsbyname=CONN_x.Setup.conditions.names(CONN_x.Results.xX.nconditions);

            if isfield(batch.Results,'saveas')&&~isempty(batch.Results.saveas)
                conn_contrastmanager('add',0,batch.Results.saveas);
            end
            if ~isvv&&any(CONN_x.Analyses(CONN_x.Analysis).type==[1,3]), % &&isfield(batch.Results,'done')&&batch.Results.done
                CONN_x.gui=struct('overwrite','Yes');
                if isfield(batch.Results,'display'), CONN_x.gui.display_results=batch.Results.display; end
                conn_process('results_roi');
                CONN_x.gui=1;
                CONN_x.Results.foldername=[];
                %conn save;
            end
            
            if isvv
                if ~isfield(batch.Results,'between_measures')||isempty(batch.Results.between_measures),
                    if 1,%isfield(batch.Results,'done')&&batch.Results.done
                        clear batchtemp;
                        if isfield(batch,'filename'), batchtemp.filename=batch.filename; else, batchtemp.filename=CONN_x.filename; end
                        batchtemp.Results=batch.Results;
                        for nmeasure=1:length(CONN_x.vvAnalyses(CONN_x.vvAnalysis).measures),
                            batchtemp.Results.between_measures.effect_names=conn_v2v('cleartext',CONN_x.vvAnalyses(CONN_x.vvAnalysis).measures(nmeasure));
                            batchtemp.Results.between_measures.contrast=[1];
                            conn_batch(batchtemp);
                        end
                    end
                elseif isfield(batch.Results,'between_measures'),
                    CONN_x.Results.xX.nmeasures=zeros(1,length(batch.Results.between_measures.effect_names));
                    CONN_x.Results.xX.nmeasuresbyname=cell(1,length(batch.Results.between_measures.effect_names));
                    for neffect=1:length(batch.Results.between_measures.effect_names),
                        idx=strmatch(batch.Results.between_measures.effect_names{neffect},conn_v2v('cleartext',CONN_x.vvAnalyses(CONN_x.vvAnalysis).measures),'exact');
                        if isempty(idx), error(['unknown measure ',batch.Results.between_measures.effect_names{neffect}]); return;
                        else
                            CONN_x.Results.xX.nmeasuresbyname(neffect)=conn_v2v('cleartext',CONN_x.vvAnalyses(CONN_x.vvAnalysis).measures(idx(1)));
                            CONN_x.Results.xX.nmeasures(neffect)=idx(1);
                        end
                    end
                    CONN_x.Results.xX.cmeasures=batch.Results.between_measures.contrast;
                    %conn save;
                    
                    if 1, %isfield(batch.Results,'done')&&batch.Results.done,
                        CONN_x.gui=struct('overwrite','Yes');
                        if isfield(batch.Results,'display'), CONN_x.gui.display_results=batch.Results.display; end
                        conn_process('results_voxel','dosingle','voxel-to-voxel');
                        CONN_x.gui=1;
                        CONN_x.Results.foldername=[];
                        %conn save;
                    end
                end
            else
                if any(CONN_x.Analyses(CONN_x.Analysis).type==[2,3]) && (~isfield(batch.Results,'between_sources')||isempty(batch.Results.between_sources)),
                    if 1, %isfield(batch.Results,'done')&&batch.Results.done
                        clear batchtemp;
                        if isfield(batch,'filename'), batchtemp.filename=batch.filename; else, batchtemp.filename=CONN_x.filename; end
                        batchtemp.Results=batch.Results;
                        for nsource=1:length(CONN_x.Analyses(CONN_x.Analysis).sources),
                            batchtemp.Results.between_sources.effect_names={CONN_x.Analyses(CONN_x.Analysis).sources{nsource}};
                            batchtemp.Results.between_sources.contrast=[1];
                            conn_batch(batchtemp);
                        end
                    end
                elseif isfield(batch.Results,'between_sources')&&any(CONN_x.Analyses(CONN_x.Analysis).type==[2,3]),
                    CONN_x.Results.xX.nsources=zeros(1,length(batch.Results.between_sources.effect_names));
                    CONN_x.Results.xX.nsourcesbyname=cell(1,length(batch.Results.between_sources.effect_names));
                    for neffect=1:length(batch.Results.between_sources.effect_names),
                        idx=strmatch(batch.Results.between_sources.effect_names{neffect},CONN_x.Analyses(CONN_x.Analysis).sources,'exact');
                        if isempty(idx), idx=strmatch(batch.Results.between_sources.effect_names{neffect},CONN_x.Analyses(CONN_x.Analysis).sources); end
                        if isempty(idx), error(['unknown source ',batch.Results.between_sources.effect_names{neffect}]); return;
                        elseif numel(idx)>1, error(['multiple possible matches to source ',batch.Results.between_sources.effect_names{neffect}]); return;
                        else
                            CONN_x.Results.xX.nsources(neffect)=idx(1);
                            CONN_x.Results.xX.nsourcesbyname(neffect)=CONN_x.Analyses(CONN_x.Analysis).sources(idx(1));
                        end
                    end
                    CONN_x.Results.xX.csources=batch.Results.between_sources.contrast;
                    %conn save;
                    
                    if 1, %isfield(batch.Results,'done')&&batch.Results.done,
                        CONN_x.gui=struct('overwrite','Yes');
                        if isfield(batch.Results,'display'), CONN_x.gui.display_results=batch.Results.display; end
                        conn_process('results_voxel','dosingle','seed-to-voxel');
                        CONN_x.gui=1;
                        CONN_x.Results.foldername=[];
                        %conn save;
                    end
                end
            end
        end
    end
end

if isfield(batch,'filename'), try, conn_projectmanager('tag',''); end; end

end

function conn_batch_eval(filename)
str=conn_fileutils('fileread',filename);
str=regexprep(str,'\s*function .*?\n(.*?)(end[\s\n]*)?$','$1');
evalin('base',str);
end


