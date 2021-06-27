function fileout=evlab17_run_results(varargin)
% EVLAB17_RUN_RESULTS runs second-level analyses
%
%  OPTIONS (entered in .cfg file with fieldnames preceeded by #, or entered as argument pairs to evlab17_run_results)
%
%      data            : list of nifti files entered into second-level analysis (Nsubjects x Nmeasures)
%                           note: when entering multiple files per subject (e.g. repeated measures) enter first all files (one per subject) for measure#1, followed by all files for measure#2, etc.
%                        alternatively list of SPM.mat files containing first-level analyses (Nsubjects x 1, or Nsubjects x Nmeasures)
%                        alternatively list of folder names containing SPM.mat first-level analyses (Nsubjects x 1, or Nsubjects x Nmeasures)
%      design          : design matrix (Neffects x Nsubjects) 
%                           enter one row for each modeled effect (across subjects)
%                           each row should contain one value/number per subject
%      contrast_between: between-subjects contrast vector/matrix (Nc1 x Neffects) 
%      contrast_within : within-subjects contrast vector/matrix (Nc2 x Nmeasures)
%      mask            : (optional) analysis mask; default no masking
%      contrast_names  : (optional, only when entering SPM.mat files in #data field) list of contrast names to select from first-level analysis files (Nmeasures x 1)
%      data_labels     : (optional) labels of columns of data matrix
%      design_labels   : (optional) labels of columns of design matrix
%      analysistype    : (optional) analysis type 1: include both parametric and non-parametric stats; 2: include only parametric stats (Random Field Theory assumptions); 3: include only non-parametric stats (permutation/randomization analyses)
%      folder          : (optional) folder where analysis are stored; default current folder
%

% note: this functionality moved to conn_module('glm',...) (2020 alfnie)

evlab17_module init silent;
fileout=[];
spmfolder=conn_module('glm',varargin{:});
conn_fixpermissions(spmfolder);

end




