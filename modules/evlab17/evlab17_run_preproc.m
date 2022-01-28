function fileout=evlab17_run_preproc(varargin)
% EVLAB17_RUN_PREPROC runs preprocessing pipeline
%   evlab17_run_preproc('/myfolder/run_preproc.cfg'); 
%      loads preprocessing info from /myfolder/run_preproc.cfg file
%      and runs associated pipeline
%
%   evlab17_run_preproc('/myfolder/datafiles.cfg','pipeline_preproc_DefaultMNI.cfg'); 
%      loads data info from /myfolder/datafiles.cfg file
%      and preprocessing info from preproc_pipeline_DefaultMNI.cfg
%      and runs associated pipeline
%
%   evlab17_run_preproc('/myfolder/run_preproc.cfg', ..., [], fieldname1, fieldvalue1, ...)
%      loads preprocessing info from run_preproc.cfg file
%      and additional preprocessing info from fieldname/fieldvalue pairs
%      and runs associated pipeline
%
%  OPTIONS (entered in .cfg file with fieldnames preceeded by #, or entered as argument pairs to evlab17_run_preproc)
%
%  DATA DEFINITIONS:
%      functionals     : list of functional files, if appropriate
%                       enter only one file per run/session
%                       if starting from NIFTI files, enter here the functional files (e.g. /mydata/nii/984000-7.nii)
%                       if starting from DICOM files (see 'dicoms' field below), enter here the runs that contain functional data
%                           files may be specified by DICOM series numbers (e.g. entering 7 identifies files named *-7.nii output by the dicom converter step)
%                           files may be specified by filenames without extension (e.g. entering 984000-7 identifies the 984000-7.nii files generated from 984000-7-*.dcm files)
%                           enter simply * to indicate all DICOM series
%                           enter simply ? (or do not include the 'functionals' field) to let users select functional files interactively (GUI)
%      structurals     : list of structural files, if appropriate
%                       enter either a single file (for session-invariant structurals) or one file per run (for session-specific structurals)
%                       if starting from NIFTI files, enter here the functional files (e.g. /mydata/nii/984000-3.nii)
%                       if starting from DICOM files (see 'dicoms' field below), enter here the filenames/runs that contain structural data
%                           files may be specified by DICOM series numbers (e.g. entering 3 identifies files named *-3.nii output by the dicom converter step)
%                           files may be specified by filenames without extension (e.g. 984000-3)
%                           enter simply * to indicate all DICOM series
%                           enter simply ? to let user select functional files interactively (GUI)
%      dicoms          : list of dicom files, if appropriate
%                           for <root>-<slice#>.dcm file naming convention: enter one file per run (only first file -1.dcm from each dicom series); output files will be named <root>.nii
%                           for other file naming conventions: enter all dicom files; output files will be broken down by DICOM series and renamed as run-<series#>.nii
%                       files may be specified explicitly (e.g. /mydata/dicoms/984000-1-1.dcm)
%                       files may also include wildcards (e.g. /mydata/dicoms/*-1.dcm)
%      vdm_functionals: list if voxel-displacement maps, if appropriate (for preprocessing.steps=='functional_realign&unwarp&phasemap')
%                       enter a single file or one file per run/session (vdm* file; explicitly entering these volumes here superceeds CONN's 
%                       default option to search for/use vdm* files in same directory as functional data) 
%      fmap_functionals: list if fieldmap volumes, if appropriate (for preprocessing.steps=='functional_vdm_create')
%                       enter either a) magnitude1+phasediff images; b) real1+imag1+real2+imag2; or c) fieldmap (in Hz) fieldmap acquisition sequence volumes
%                       (note: use fmap_functionals_ss fieldname if entering session/run-specific sets of fieldmap files)
%      rois            : list of roi files, if appropriate
%
%    (additional data-definition options:)
%      RT              : repetition time in seconds
%      functionals_path: (when entering nifti files as functionals) folder where nifti files are stored [default: full paths are explicitly defined in functionals field]
%      structurals_path: (when entering nifti files as structurals) folder where dicom files are stored [default: full paths are explicitly defined in structurals field]
%      rois_path: (when entering nifti files as rois) folder where roi files are stored [default: full paths are explicitly defined in rois field]
%      dicoms_path: (when entering dicom files as structurals and/or functionals) folder where dicom files are stored [full paths are explicitly defined in dicoms field]
%      dicoms_outputfolder: (when entering dicom files as structurals and/or functionals) folder where nii files are stored (see 'folderout' options in conn_dcm2nii)  ['../nii']
%      dicoms_outputoverwrite: (when entering dicom files as structurals and/or functionals) 1/0 overwrite existing .nii files if they exist [0]
%
%  PREPROCESSING-STEPS DEFINITIONS:
%      steps           : list of preprocessing steps (in order) to be performed. Available steps are:
%        functional_art                        : functional identification of outlier scans (from motion displacement and global signal changes)
%        functional_bandpass                   : functional band-pass filtering
%        functional_center                     : centers functional data to origin (0,0,0) coordinates
%        functional_centertostruct             : centers functional data to approximate structural-volume coordinates
%        functional_coregister_affine_reslice  : functional affine coregistration to structural volumes
%        functional_coregister_affine_noreslice: functional affine coregistration to structural volumes without reslicing (applies transformation to source header files)
%        functional_coregister_nonlinear       : functional non-linear coregistration to structural volumes
%        functional_label                      : labels current functional files as part of list of Secondary Datasets
%        functional_manualorient               : applies user-defined affine transformation to functional data
%        functional_manualspatialdef           : applies user-defined spatial deformation to functional data
%        functional_motionmask                 : creates functional motion masks (mean BOLD signal spatial derivatives wrt motion parameters)
%        functional_normalize_direct           : functional direct normalization
%        functional_normalize_indirect         : functional indirect normalization (coregister to structural; normalize structural; apply same transformation to functionals)
%        functional_normalize_indirect_preservemasks: functional indirect normalization with user-defined Grey/White/CSF masks (coregister to structural; normalize structural; apply same transformation to functionals as well as to Grey/White/CSF masks)
%        functional_realign                    : functional realignment
%        functional_realign_noreslice          : functional realignment without reslicing (applies transformation to source header files)
%        functional_realign&unwarp             : functional realignment & unwarp (motion-by-inhomogeneity interactions)
%        functional_realign&unwarp&fieldmap    : functional realignemnt & unwarp & inhomogeneity correction (from vdm/phasemap files)
%        functional_regression                 : removal of user-defined temporal components from BOLD timeseries (keeps residuals of linear regression model)
%        functional_removescans                : removes user-defined number of initial scans from functional data
%        functional_segment                    : functional segmentation (Gray/White/CSF tissue classes)
%        functional_segment&normalize_direct   : functional direct unified normalization and segmentation
%        functional_segment&normalize_indirect : functional indirect unified normalization and segmentation (coregister to structural; normalize and segment structural; apply same transformation to functionals)
%        functional_slicetime                  : functional slice-timing correction
%        functional_smooth                     : functional spatial smoothing (spatial convolution with Gaussian kernel)
%        functional_smooth_masked              : functional spatial masked-smoothing (spatial convolution with Gaussian kernel restricted to voxels within Grey Matter mask)
%        functional_surface_coreg&resample     : coregister&resample functional data at the location of FreeSurfer subject-specific structural cortical surface
%        functional_surface_resample           : resample functional data at the location of FreeSurfer subject-specific structural cortical surface
%        functional_surface_smooth             : functional spatial diffusion of surface data
%        functional_vdm_create                 : creation of vdm (voxel-displacement-map) from fieldmap dataset (reads 'fmap' secondary functional dataset containing magnitudfe and phasediff images and creates 'vdm' secondary functional dataset containing voxel-displacement map)
%        structural_center                     : centers structural data to origin (0,0,0) coordinates
%        structural_manualorient               : applies user-defined affine transformation to structural data
%        structural_manualspatialdef           : applies user-defined spatial deformation to structural data
%        structural_segment&normalize          : structural unified normalization and segmentation 
%        structural_normalize                  : structural normalization to MNI space (without segmentation)
%        structural_normalize_preservemasks    : structural normalization to MNI space with user-defined Grey/White/CSF masks (normalizes structural data and applies same transformation to user-defined Grey/White/CSF mask ROIs)
%        structural_segment                    : structural segmentation (Grey/White/CSF tissue classes)
%
%   (additional preprocessing-step options are required only when including specific preprocessing steps, default values are specified below between brackets)
%      fwhm            : (functional_smooth) Smoothing factor (mm) [8]
%      coregtomean     : (functional_coregister/segment/normalize) 0: use first volume; 1: use mean volume (computed during realignment); 2: use user-defined source volume (see Setup.coregsource_functionals field) [1]
%      sliceorder      : (functional_slicetime) acquisition order (vector of indexes; 1=first slice in image; note: use cell array for subject-specific vectors)
%                       alternatively sliceorder may also be defined as one of the following strings: 'ascending','descending','interleaved (middle-top)','interleaved (bottom-up)','interleaved (top-down)','interleaved (Siemens)','BIDS'
%                       alternatively sliceorder may also be defined as a vector containing the acquisition time in milliseconds for each slice (e.g. for multi-band sequences)
%      ta              : (functional_slicetime) acquisition time (TA) in seconds (used to determine slice times when sliceorder is defined by a vector of slice indexes; note: use vector for subject-specific values). Defaults to (1-1/nslices)*TR where nslices is the number of slices
%      art_thresholds  : (functional_art) ART thresholds for identifying outlier scans
%                                            art_thresholds(1): threshold value for global-signal (z-value; default 5)
%                                            art_thresholds(2): threshold value for subject-motion (mm; default .9)
%                        additional options: art_thresholds(3): 1/0 global-signal threshold based on scan-to-scan changes in global-BOLD measure (default 1)
%                                            art_thresholds(4): 1/0 subject-motion threshold based on scan-to-scan changes in subject-motion measure (default 1)
%                                            art_thresholds(5): 1/0 subject-motion threhsold based on composite-movement measure (default 1)
%                                            art_thresholds(6): 1/0 force interactive mode (ART gui) (default 0)
%                                            art_thresholds(7): [only when art_threshold(5)=0] subject-motion threshold based on rotation measure
%                                            art_thresholds(8): N number of initial scans to be flagged for removal (default 0)
%                            note: when art_threshold(5)=0, art_threshold(2) defines the threshold based on the translation measure, and art_threhsold(7) defines the threshold based on the rotation measure; otherwise art_threshold(2) defines the (single) threshold based on the composite-motion measure
%                            note: the default art_thresholds(1:2) [5 .9] values correspond to the "intermediate" (97th percentile) settings, to use the "conservative" (95th percentile) settings use [3 .5], to use the "liberal" (99th percentile) settings use [9 2] values instead
%                            note: art needs subject-motion files to estimate possible outliers. If a 'realignment' first-level covariate exists it will load the subject-motion parameters from that first-level covariate; otherwise it will look for a rp_*.txt file (SPM format) in the same folder as the functional data
%                                  subject-motion files can be in any of the following formats: a) *.txt file (SPM format; three translation parameters in mm followed by pitch/roll/yaw in radians); b) *.par (FSL format; three Euler angles in radians followed by translation parameters in mm); c) *.siemens.txt (Siemens MotionDetectionParameter.txt format); d) *.deg.txt (same as SPM format but rotations in degrees instead of radians)
%      removescans     : (functional_removescans) number of initial scans to remove
%      reorient        : (functional/structural_manualorient) 3x3 or 4x4 transformation matrix or filename containing corresponding matrix
%      respatialdef    : (functional/structural_manualspatialdef) nifti deformation file (e.g. y_*.nii or *seg_sn.mat files)
%      voxelsize_anat  : (normalization) target voxel size for resliced anatomical volumes (mm) [1]
%      voxelsize_func  : (normalization) target voxel size for resliced functional volumes (mm) [2]
%      boundingbox     : (normalization) target bounding box for resliced volumes (mm) [-90,-126,-72;90,90,108]
%      interp          : (normalization) target voxel interpolation method (0:nearest neighbor; 1:trilinear; 2 or higher:n-order spline) [4]
%      template_anat   : (structural_normalize SPM8 only) anatomical template file for approximate coregistration [spm/template/T1.nii]
%      template_func   : (functional_normalize SPM8 only) functional template file for normalization [spm/template/EPI.nii]
%      affreg          : (normalization) affine registration before normalization ['mni']
%      tpm_template    : (structural_segment, structural_segment&normalize in SPM8, and any segment/normalize option in SPM12) tissue probability map [spm/tpm/TPM.nii]
%      tpm_ngaus       : (structural_segment, structural_segment&normalize in SPM8&SPM12) number of gaussians for each tissue probability map
%      diffusionsteps  : (surface_smooth) number of diffusion steps
%      vdm_et1         : (functional_vdm_create) ET1 (Echo Time first echo in fieldmap sequence) (default [] : read from .json file / BIDS)
%      vdm_et2         : (functional_vdm_create) ET2 (Echo Time second echo in fieldmap sequence) (default [] : read from .json file / BIDS)
%      vdm_ert         : (functional_vdm_create) ERT (Effective Readout Time in funcional data) (default [] : read from .json file / BIDS)
%      vdm_blip        : (functional_vdm_create) k-space traversal blip direction (+1 or -1; default -1)
%      vdm_type        : (functional_vdm_create only) type of fieldmap sequence files ([]: automatically detect; 1: magnitude+phasediff (or magnitude1+magnitude2+phasediff); 2: real1+imag1+real2+imag2; 3: fieldmapHz)
%      vdm_fmap        : (functional_vdm_create only) location of fieldmap sequence files (secondary functional dataset number or label containing fieldmap sequence files) ['fmap']
%      bp_filter       : (functional_bandpass) Low- and High- frequency thresholds (in Hz)
%      bp_keep0        : (functional_bandpass) 0: removes average BOLD signal (freq=0Hz component); 1: keeps average BOLD signal in output independent of band-pass filter values; [1]
%      reg_names       : (functional_regression) list of first-level covariates to use as model regressors / design matrix (valid entries are first-level covariate names)
%      reg_dimensions  : (functional_regression) maximum number of dimensions for each model regressor (inf to use all dimensions) [inf]
%      reg_deriv       : (functional_regression) 0: none; 1: adds first-order derivatves; 2: adds second-order derivatives to each model regressor [0]
%      reg_skip        : (functional_regression) 1: does not create output functional files, only creates session-specific dp_*.txt files with covariate timeseries to be included later in an arbitrary first-level model [0]
%      label           : (functional_label) label of secondary dataset [datestr(now)]
%
%   (general options, default values are specified below between brackets)
%      dogui           : 1/0 launch GUI interface to edit/modify preprocessing the pipeline before starting
%      localcopy       : 1/0 copy original functional/structural data to local BIDS folder [0]
%      parallel.N      : (parallelization options) number of parallel jobs [0]
%      parallel.profile: (parallelization options) cluster computing profile (see "conn_jobmanager profiles" to see a list of available parallelization profiles) ['Grid Engine']
%      qa_plots        : (segmentation and/or normalization steps) 1/0 to create quality assurance plots [0]
%      qa_plist        : (segmentation and/or normalization steps) list of quality assurance plots to create (default value: [] for all QA displays)
%      qa_folder       : (segmentation and/or normalization steps) output directory for QA plots

%      qa_parallel     : (model review displays) 1/0 to create plots in parallel (using background or cluster computing resources) [0]
%      qa_profile      : (model review displays) cluster computing profile for QA plots (see "conn_jobmanager profiles" to see a list of available parallelization profiles) [background process]
%

evlab17_module init;
fileout=[];

% loads .cfg files
options=struct;
for n=1:nargin
    if isempty(varargin{n}), break; end
    filename=varargin{n};
    if ischar(filename)
        if isempty(dir(filename)),
            if ~isempty(dir(fullfile(fileparts(which(mfilename)),filename))), filename=fullfile(fileparts(which(mfilename)),filename);
            else
                fprintf('warning: file %s not found\n',filename);
                filename=which(filename);
            end
        end
        fprintf('loading file %s\n',filename);
    end
    options=conn_loadcfgfile(filename,options);
end
for n=n+1:2:nargin-1
    fieldname=regexp(varargin{n},'\.','split');
    fieldvalue=varargin{n+1};
    options=setfield(options,fieldname{:},fieldvalue);
end
fprintf('%s options:\n',mfilename);
disp(options);
options0=options;
options0.arguments=varargin;

% interprets info
tag=datestr(now,'yyyy_mm_dd_HHMMSSFFF');
fields={};
functional_runs=[];
structural_runs=[];
functional_selected=[];
structural_selected=[];
functional_files={};
structural_files={};
fpaths={};
niifolder=[];
dicom_options={};
isdataset=isfield(options,'dataset');
qa_plots=evlab17_module('default','qa_plots');
if isdataset
    fileout=conn_prepend('',options.dataset,'.mat');
    %qafolder=fullfile(fileparts(fileout),['QA_',tag]);
    %%qafolder=fullfile(conn_prepend('',fileout,''),'results','qa',['QA_',tag]);
    options=rmfield(options,'dataset');
else
    fileout=fullfile(pwd,['evlab17_',tag,'.mat']); 
    %qafolder=fullfile(niifolder,['QA_',tag]);
end

if isfield(options,'structurals')
    if ischar(options.structurals), options.structurals=cellstr(options.structurals); end
    if isequal(options.structurals,{'*'})
        structural_runs='*';
    elseif isequal(options.structurals,{'?'})
        structural_runs=[];
    elseif isnumeric(options.structurals)   % number indicate session numbers
        structural_runs=options.structurals;
    else
        %[fpaths2,fnames2]=cellfun(@fileparts,options.structurals,'uni',0);
        %partialfilenames=cellfun('length',fpaths2)==0;
        partialfilenames=cellfun(@(x)isempty(regexp(x,'^[\\\/]')),options.structurals);
        allfiles={};                    % strings indicate filenames (with potential wildcards)
        for n=1:numel(options.structurals)
            if partialfilenames(n)&&isfield(options,'structurals_path')
                if isempty(regexp(char(options.structurals_path),'^[\\\/]')), options.structurals_path=conn_fullfile(conn_prepend('',fileout,''),char(options.structurals_path)); end
                partialfilenames(n)=false;
                files=conn_dir(fullfile(char(options.structurals_path),options.structurals{n}),'-R');
                if isempty(files), error('file %s not found',fullfile(char(options.structurals_path),options.structurals{n})); end
                allfiles=cat(1,allfiles,conn_sortfilenames(cellstr(files)));
            elseif partialfilenames(n)&&~isfield(options,'structurals_path'), 
                allfiles=cat(1,allfiles,options.structurals(n));
            else
                files=conn_dir(options.structurals{n},'-R');
                if isempty(files), error('file %s not found',options.structurals{n}); end
                allfiles=cat(1,allfiles,conn_sortfilenames(cellstr(files)));
            end
        end
        fields{end+1}='structurals';
        fields{end+1}={allfiles};
        structural_selected=allfiles;
        fprintf('structural files:\n');
        disp(char(allfiles));
        structural_files=allfiles;
        if ~isempty(structural_files), niifolder=fileparts(structural_files{1}); end
    end
    options=rmfield(options,'structurals');
end
if isfield(options,'structurals_path'), options=rmfield(options,'structurals_path'); end
if isfield(options,'functionals')
    if ischar(options.functionals), options.functionals=cellstr(options.functionals); end
    if isequal(options.functionals,{'*'})
        functional_runs=inf;
    elseif isequal(options.functionals,{'?'})
        functional_runs=[];
    elseif isnumeric(options.functionals)   % number indicate session numbers
        functional_runs=options.functionals;
    else
        %[fpaths2,fnames2]=cellfun(@fileparts,options.functionals,'uni',0);
        %partialfilenames=cellfun('length',fpaths2)==0;
        partialfilenames=cellfun(@(x)isempty(regexp(x,'^[\\\/]')),options.functionals);
        allfiles={};                    % strings indicate filenames (with potential wildcards)
        for n=1:numel(options.functionals)
            if partialfilenames(n)&&isfield(options,'functionals_path')
                if isempty(regexp(char(options.functionals_path),'^[\\\/]')), options.functionals_path=conn_fullfile(conn_prepend('',fileout,''),char(options.functionals_path)); end
                partialfilenames(n)=false;
                files=conn_dir(fullfile(char(options.functionals_path),options.functionals{n}),'-R');
                if isempty(files), error('file %s not found',fullfile(char(options.functionals_path),options.functionals{n})); end
                allfiles=cat(1,allfiles,conn_sortfilenames(cellstr(files)));
            elseif partialfilenames(n)&&~isfield(options,'functionals_path'), 
                allfiles=cat(1,allfiles,options.functionals(n));
            else
                files=conn_dir(options.functionals{n},'-R');
                if isempty(files), error('file %s not found',options.functionals{n}); end
                allfiles=cat(1,allfiles,conn_sortfilenames(cellstr(files)));
            end
        end
        fields{end+1}='functionals';
        fields{end+1}={allfiles};
        functional_selected=allfiles;
        fprintf('functional files:\n');
        disp(char(allfiles));
        functional_files=allfiles;
        if ~isempty(functional_files), niifolder=fileparts(functional_files{1}); end
    end
    options=rmfield(options,'functionals');
end
if isfield(options,'vdm_functionals')&&~isfield(options,'unwarp_functionals'), options.unwarp_functionals=options.vdm_functionals; options=rmfield(options,'vdm_functionals'); end
if isfield(options,'vdm_functionals_ss')&&~isfield(options,'unwarp_functionals_ss'), options.unwarp_functionals_ss=options.vdm_functionals_ss; options=rmfield(options,'vdm_functionals_ss'); end
if isfield(options,'unwarp_functionals_ss'), options.unwarp_functionals=options.unwarp_functionals_ss; end
if isfield(options,'unwarp_functionals')
    if ischar(options.unwarp_functionals), options.unwarp_functionals={options.unwarp_functionals}; end
    if isfield(options,'functionals_path')
        %[fpaths2,fnames2]=cellfun(@fileparts,options.unwarp_functionals,'uni',0);
        %partialfilenames=cellfun('length',fpaths2)==0;
        partialfilenames=cellfun(@(x)isempty(regexp(x,'^[\\\/]')),options.unwarp_functionals);
        allfiles={};                    % strings indicate filenames (with potential wildcards)
        for n=1:numel(options.unwarp_functionals)
            if partialfilenames(n)&&isfield(options,'functionals_path')
                if isempty(regexp(char(options.functionals_path),'^[\\\/]')), options.functionals_path=conn_fullfile(conn_prepend('',fileout,''),char(options.functionals_path)); end
                partialfilenames(n)=false;
                files=conn_dir(fullfile(char(options.functionals_path),options.unwarp_functionals{n}),'-R');
                if isempty(files), error('file %s not found',fullfile(char(options.functionals_path),options.unwarp_functionals{n})); end
                allfiles=cat(1,allfiles,conn_sortfilenames(cellstr(files)));
            elseif partialfilenames(n)&&~isfield(options,'functionals_path'),
                allfiles=cat(1,allfiles,options.unwarp_functionals(n));
            else
                files=conn_dir(options.unwarp_functionals{n},'-R');
                if isempty(files), error('file %s not found',options.unwarp_functionals{n}); end
                allfiles=cat(1,allfiles,conn_sortfilenames(cellstr(files)));
            end
        end
    else allfiles=options.unwarp_functionals;
    end
    unwarp_files=allfiles;
    options.unwarp_functionals={allfiles};
end
if isfield(options,'fmap_functionals_ss'), options.fmap_functionals=options.fmap_functionals_ss; end
if isfield(options,'fmap_functionals')
    if ischar(options.fmap_functionals), options.fmap_functionals={options.fmap_functionals}; end
    if isfield(options,'functionals_path')
        %[fpaths2,fnames2]=cellfun(@fileparts,options.fmap_functionals,'uni',0);
        %partialfilenames=cellfun('length',fpaths2)==0;
        partialfilenames=cellfun(@(x)isempty(regexp(x,'^[\\\/]')),options.fmap_functionals);
        allfiles={};                    % strings indicate filenames (with potential wildcards)
        for n=1:numel(options.fmap_functionals)
            if partialfilenames(n)&&isfield(options,'functionals_path')
                if isempty(regexp(char(options.functionals_path),'^[\\\/]')), options.functionals_path=conn_fullfile(conn_prepend('',fileout,''),char(options.functionals_path)); end
                partialfilenames(n)=false;
                files=conn_dir(fullfile(char(options.functionals_path),options.fmap_functionals{n}),'-R');
                if isempty(files), error('file %s not found',fullfile(char(options.functionals_path),options.fmap_functionals{n})); end
                allfiles=cat(1,allfiles,conn_sortfilenames(cellstr(files)));
            elseif partialfilenames(n)&&~isfield(options,'functionals_path'),
                allfiles=cat(1,allfiles,options.fmap_functionals(n));
            else
                files=conn_dir(options.fmap_functionals{n},'-R');
                if isempty(files), error('file %s not found',options.fmap_functionals{n}); end
                allfiles=cat(1,allfiles,conn_sortfilenames(cellstr(files)));
            end
        end
    else allfiles=options.fmap_functionals;
    end
    fmap_files=allfiles;
    options.fmap_functionals={{char(allfiles)}}; % same set of fmap files for all runs 
end
if isfield(options,'functionals_path'), options=rmfield(options,'functionals_path'); end
if isfield(options,'rois_ss'), options.rois=options.rois_ss; end
if isfield(options,'rois')
    if ischar(options.rois), options.rois={options.rois}; end
    if iscell(options.rois), options.rois=struct('files',{options.rois}); % interpreted as one file per ROI
    elseif isstruct(options.rois)
        fnames=fieldnames(options.rois);
        newrois=struct('names',{{}},'files',{{}});
        for n=1:numel(fnames)
            fstruct=options.rois.(fnames{n});
            if ischar(fstruct), fstruct={fstruct}; end
            if iscell(fstruct), fstruct=struct('files',{fstruct}); end % interpreted as one file per run/session
            newrois.names{n}=fnames{n};
            newrois.files{n}{1}=cellstr(fstruct.files);
            if isfield(fstruct,'names'), newrois.names{n}=fstruct.names; end
            if isfield(fstruct,'dimensions'), newrois.dimensions{n}=fstruct.dimensions; end
            if isfield(fstruct,'weighted'), newrois.weighted(n)=fstruct.weighted; end
            if isfield(fstruct,'multiplelabels'), newrois.multiplelabels(n)=fstruct.multiplelabels; end
            if isfield(fstruct,'mask'), newrois.mask(n)=fstruct.mask; end
            if isfield(fstruct,'regresscovariates'), newrois.regresscovariates(n)=fstruct.regresscovariates; end
            if isfield(fstruct,'dataset'), newrois.dataset{n}=fstruct.dataset; end
            if isfield(fstruct,'add'), newrois.add=fstruct.add; end
        end
        options.rois=newrois;
    end
    if isfield(options,'rois_path')
        for nrois=1:numel(options.rois.files)
            roifiles=options.rois.files{nrois};
            if iscell(roifiles)&&numel(roifiles)==1, roifiles=roifiles{1}; end
            if ~iscell(roifiles), roifiles={roifiles}; end
            partialfilenames=cellfun(@(x)isempty(regexp(x,'^[\\\/]')),roifiles);
            roi_files={};                    % strings indicate filenames (with potential wildcards)
            for n=1:numel(roifiles)
                if partialfilenames(n)&&isfield(options,'rois_path')
                    if isempty(regexp(char(options.rois_path),'^[\\\/]')), options.rois_path=conn_fullfile(conn_prepend('',fileout,''),char(options.rois_path)); end
                    partialfilenames(n)=false;
                    files=conn_dir(fullfile(char(options.rois_path),roifiles{n}),'-R');
                    if isempty(files), error('file %s not found',fullfile(char(options.rois_path),roifiles{n})); end
                    roi_files=cat(1,roi_files,conn_sortfilenames(cellstr(files)));
                elseif partialfilenames(n)&&~isfield(options,'rois_path'),
                    roi_files=cat(1,roi_files,roifiles(n));
                else
                    files=conn_dir(roifiles{n},'-R');
                    if isempty(files), error('file %s not found',roifiles{n}); end
                    roi_files=cat(1,roi_files,conn_sortfilenames(cellstr(files)));
                end
            end
            options.rois.files{nrois}{1}=roi_files; 
        end
    end
    if ~isfield(options.rois,'add')&&isfield(options,'isnew')&&options.isnew, options.rois.add=1; end % note: fix to keep CONN's default ROIs
end
if isfield(options,'rois_path'), options=rmfield(options,'rois_path'); end
if isfield(options,'dicoms_outputfolder')
    dicom_options{end+1}='folderout';
    dicom_options{end+1}=char(options.dicoms_outputfolder);
    options=rmfield(options,'dicoms_outputfolder');
    niifolder=char(options.dicoms_outputfolder);
end
if isfield(options,'dicoms_outputoverwrite')
    dicom_options{end+1}='overwrite';
    dicom_options{end+1}=options.dicoms_outputoverwrite;
    options=rmfield(options,'dicoms_outputoverwrite');
end
if isfield(options,'dicoms')
    allfiles={};
    for n=1:numel(options.dicoms)
        if isfield(options,'dicoms_path')
            if isempty(regexp(char(options.dicoms_path),'^[\\\/]')), options.dicoms_path=conn_fullfile(conn_prepend('',fileout,''),char(options.dicoms_path)); end
            files=conn_dir(fullfile(char(options.dicoms_path),options.dicoms{n}),'-R');
            if isempty(files), error('file %s not found',fullfile(char(options.dicoms_path),options.dicoms{n})); end
        else
            files=conn_dir(options.dicoms{n},'-R');
            if isempty(files), error('file %s not found',options.dicoms{n}); end
        end
        allfiles=cat(1,allfiles,conn_sortfilenames(cellstr(files)));
    end
    %[outstruct,files,fdescrip,ftypes,frt]=conn_dcmconvert(allfiles,'folderout','../nii','overwrite',false,dicom_options{:});
    [files,ftypes,frt,fdescrip]=conn_dcm2nii(allfiles,dicom_options{:});
    if isfield(options,'RT')
        if any(frt(~isnan(frt))~=options.RT), fprintf('WARNING: Inconsistent RT values found in dicom (%s) compared to RT field entered in .cfg file (%s). The value in .cfg file takes precedence',mat2str(frt),mat2str(options.RT)); end
    else
        frt=unique(frt(~isnan(frt))); 
        if numel(frt)==1, options.RT=frt; end
    end
    [fpaths,fnames]=cellfun(@fileparts,files,'uni',0);
    niifolder=fpaths{1};
    %%%%
    if ~isempty(functional_selected)
        [fpaths2,fnames2]=cellfun(@fileparts,functional_selected,'uni',0);
        partialfilenames=find(cellfun('length',fpaths2)==0);
        for n=1:numel(partialfilenames)
            idx=find(strcmp(fnames2{partialfilenames(n)},fnames));
            if numel(idx)==1, functional_selected(partialfilenames(n))=files(idx); end
        end
        if ~isempty(partialfilenames)
            idx=2*find(strcmp(fields(1:2:end-1),'functionals'),1)-1;
            fields{idx+1}={functional_selected};
            niifolder=fileparts(functional_selected{1});
            disp('functional files:');
            disp(char(functional_selected));
            functional_files=functional_selected;
        end
    else
        if isempty(functional_runs)&&isfield(options,'steps')&&any(cellfun('length',regexp(options.steps,'^functional_'))>0)
            [functional_runs,tok] = listdlg('PromptString','Select FUNCTIONAL run(s)','SelectionMode','multiple','ListString',cellfun(@(a,b)sprintf('%s : %s',a,b),fnames,fdescrip,'uni',0),'InitialValue',find(ftypes>10));
        elseif ~isempty(functional_runs)
            if isinf(functional_runs), functional_runs=1:numel(fnames);
            else
                for n=1:numel(functional_runs),
                    idx=find(cellfun('length',regexp(fnames,sprintf('-%d$',functional_runs(n)))));
                    if numel(idx)~=1, disp(char(fnames)); error('unable to identify functional run %d',functional_runs(n)); end
                    functional_runs(n)=idx;
                end
            end
        end
        if ~isempty(functional_runs)
            fields{end+1}='functionals';
            fields{end+1}={files(functional_runs)};
            niifolder=fpaths{functional_runs(1)};
            disp('functional runs:');
            disp(char(files(functional_runs)));
            functional_files=files(functional_runs);
        end
    end
    if ~isempty(structural_selected)
        [fpaths2,fnames2]=cellfun(@fileparts,structural_selected,'uni',0);
        partialfilenames=find(cellfun('length',fpaths2)==0);
        for n=1:numel(partialfilenames)
            idx=find(strcmp(fnames2{partialfilenames(n)},fnames));
            if numel(idx)==1, structural_selected(partialfilenames(n))=files(idx); end
        end
        if ~isempty(partialfilenames)
            idx=2*find(strcmp(fields(1:2:end-1),'structurals'),1)-1;
            fields{idx+1}={structural_selected};
            niifolder=fileparts(structural_selected{1});
            disp('structural files:');
            disp(char(structural_selected));
            structural_files=structural_selected;
        end
    else
        if isempty(structural_runs)&&isfield(options,'steps')&&any(cellfun('length',regexp(options.steps,'^structural_'))>0)
            ok=false;
            showonly=1:numel(files); showonly(functional_runs)=[];
            while ~ok
                [structural_runs,tok] = listdlg('PromptString','Select STRUCTURAL run(s)','SelectionMode','multiple','ListString',cellfun(@(a,b)sprintf('%s : %s',a,b),fnames(showonly),fdescrip(showonly),'uni',0),'InitialValue',find(ftypes(showonly)==1,1));
                if ~isempty(structural_runs), structural_runs=showonly(structural_runs); end
                ok=ismember(numel(structural_runs),[0 1 numel(functional_runs)]);
            end
        elseif ~isempty(structural_runs)
            if isequal(structural_runs,'*'), structural_runs=1:numel(fnames);
            else
                for n=1:numel(structural_runs),
                    idx=find(cellfun('length',regexp(fnames,sprintf('-%d$',structural_runs(n)))));
                    if numel(idx)~=1, disp(char(fnames)); error('unable to identify structural run %d',structural_runs(n)); end
                    structural_runs(n)=idx;
                end
            end
        end
        if ~isempty(structural_runs)
            fields{end+1}='structurals';
            fields{end+1}={files(structural_runs)};
            niifolder=fpaths{structural_runs(1)};
            disp('structural files:');
            disp(char(files(structural_runs)));
            structural_files=files(structural_runs);
        end
    end
    options=rmfield(options,'dicoms');
end
if isfield(options,'dicoms_path'), options=rmfield(options,'dicoms_path'); end
assert(numel(functional_files)<=1||numel(structural_files)<=1||numel(functional_files)==numel(structural_files),'mismatched number of functional and structural runs (%d functional; %d structural). The number of structural volumes should be either one (session-independent structural) or equal to the number of functional sessions (session-specific structurals)',numel(functional_files),numel(structural_files));
if isfield(options,'unwarp_functionals_ss') % one set of fmap files per session
    assert(numel(unwarp_files)==numel(functional_files),'mismatched number of unwarp_functionals_ss files');
    unwarp_files=reshape(unwarp_files,[],numel(functional_files));
    options.unwarp_functionals={arrayfun(@(n)char(unwarp_files(:,n)),1:numel(functional_files),'uni',0)};
    options=rmfield(options,'unwarp_functionals_ss'); 
end
if isfield(options,'fmap_functionals_ss') % one set of fmap files per session
    assert(~rem(numel(fmap_files),numel(functional_files)),'mismatched number of fmap_functionals_ss files');
    fmap_files=reshape(fmap_files,[],numel(functional_files));
    options.fmap_functionals={arrayfun(@(n)char(fmap_files(:,n)),1:numel(functional_files),'uni',0)};
    options=rmfield(options,'fmap_functionals_ss'); 
end
if isempty(niifolder)
    if ~isempty(functional_files), niifolder=fileparts(functional_files{1}); 
    elseif ~isempty(structural_files), niifolder=fileparts(structural_files{1}); 
    end
    if isempty(niifolder), niifolder=pwd; end
end
if ~isdataset, fileout=fullfile(niifolder,['evlab17_',tag,'.mat']); end
if isfield(options,'parallel')||isfield(options,'localcopy'), options.filename=fileout; end

names=fieldnames(options); 
namesQA=names(cellfun('length',regexp(names,'^qa_'))>0);
optionsQA=struct('qa_plist','preprocessing');
for n=1:numel(namesQA)
    optionsQA.(namesQA{n})=options.(namesQA{n});
    options=rmfield(options,namesQA{n});
end
if isfield(options,'parallel'), 
    qaoptions.parallel=options.parallel; 
    if qa_plots, options.parallel.immediatereturn=false; end
end
if isfield(options,'subjects'), qaoptions.subjects=options.subjects; end
namesDesign=names(ismember(names,{'design','runs','files','path'})); % keep this info for future reference
optionsDesign=struct;
for n=1:numel(namesDesign)
    optionsDesign.(namesDesign{n})=options.(namesDesign{n});
    options=rmfield(options,namesDesign{n});
end
names=fieldnames(options); 
for n=1:numel(names), % all other options are passed directly to conn_module
    fields{end+1}=names{n};
    fields{end+1}=options.(names{n});
    options=rmfield(options,names{n});
end

if ~isempty(fields), 
    pwd0=pwd;
    cd(niifolder);
    evlab17_module('preprocessing',fields{:}); 
    evlab17_module('setinfo','preprocessing',fields);
    if isfield(optionsDesign,'design')||isfield(optionsDesign,'files'), evlab17_module('setinfo','design',optionsDesign); end
    evlab17_module('setinfo','run_preproc_input',options0);
    evlab17_module('save',fileout);
    fprintf('EVLAB17 dataset info saved in %s\n',fileout);
    if evlab17_module('inconnfolders'), conn_fixpermissions(fileout,[],true); % dataset folder
    else conn_fixpermissions(fileparts(fileout)); % nii folder
    end
    if qa_plots, evlab17_run_qa(optionsQA,[],'dataset',fileout); end
    cd(pwd0);
end

end




