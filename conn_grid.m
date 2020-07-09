% 
% Using CONN in a computer cluster environment
% 
% CONN grid computing options allow you to process your subjects in parallel, with each subject, 
% or each of several smaller groups of subjects, being processed independently by a different 
% computer node. This can significantly increase the speed of your analyses, allowing you to 
% process hundreds of subjects in the time it would take to process just one or a few subjects.
% 
% Configuration/Installation: %!
%
% To use this functionality you need to have access to a computer that is part of a distributed
% cluster environment. Using or installing CONN on a cluster environment does not require administrator 
% privileges. Each node in your cluster simply needs to have access to the folders containing:
% 
% * option 1) SPM, CONN, and Matlab (and associated Matlab licenses)
% * option 2) standalone CONN, and Matlab Runtime (MCR, freely available; no licenses required)
% * as well as the data folders containing your conn project and data
% 
% CONN supports natively the following cluster computing job schedulers:
% 
% * Grid Engine (Sun/Oracle Grid Engine, Open Grid Scheduler, or compatible system)
% * PBS/Torque (Portable Batch System, or compatible)
% * LFS (Platform Load Sharing Facility, or compatible)
% * Slurm (Simple Linux Utility for Resource Management, or compatible)
% 
% CONN will work out-of-the-box in many institutional clusters simply by selecting the appropriate
% job scheduler that handles job submission in your cluster environment. In addition, profiles for 
% other schedulers, or system-specific settings (e.g to increase your system default walltime settings, 
% enter optional project or account ids, etc.), can be easily created/edited in 
% CONN Tools.Cluster/HPCsettings if necessary.
% 
% To get started simply install CONN on a computer that is part of your distributed cluster environment.
% Then launch CONN GUI interactively and select Tools.Cluster/HPCsettings. There you can select and test 
% your system job scheduler. If in doubt simply try one of the four default configurations (GridEngine, 
% PBS, LSF, or SLURM) and click on 'Test profile' (each test may take up to a few minutes). If you see a 
% 'failed' message, click 'Cancel' and try a different profile. If you see a 'Test finished correctly'
% message, select this as your default profile and save your profile configuration options for future 
% sessions and/or for other users. This process needs to be performed just once. 
%
% Usage: %!
%
% When running any of CONN's analysis steps (preprocessing, Setup, Denoising, and First-level) you can 
% select between running them locally on your current computer or running them in parallel. 
%
% If using CONN's GUI, before pressing the 'Start' button to start your desired analysis step, simply 
% change the field that reads 'local processing' to 'distributed processing'. You will then be prompted
% to select the number of jobs that you want to submit (typically you may divide your analyses into as 
% many as one job per subject; keep in mind that your cluster environment may have some recommendations 
% or limitations regarding the maximum number of jobs that may be run simultaneously, or the maximum-time 
% that individual jobs will be allowed to run). CONN will then automatically:
% 
%    a) break down your analysis into the selected number of jobs
%    b) submit these jobs to your cluster (using your default cluster configuration options)
%    c) track these jobs progress and let you re-submit any failed jobs if necessary
%    and d) merge the results back into your conn project when they are finished
% 
% If using CONN's batch processing options, simply enter the desired number of jobs/nodes in the 
% batch.parallel.N field as part of your conn_batch fields. See help conn_batch for additional details
%
