function varargout = art(varargin)
% art - module for automatic and manual detection and removal of outliers.
%
% art
% Lauches the art GUI. This will prompt the user to enter the functional
% volumes and motion pararameters (SPM, FSL, or Siemens supported) for one
% subject (one or multiple sessions). It then displays four graphs:
%
% The top graph is the global brain activation mean as a function of
% time, with all of the identified outlier scans marked. Optionally it also
% overlais the task-related regressors (SPM designs only) and break down of
% number of identified outliers per condition. 
%
% The second graph shows the timeseries derived from the global BOLD signal
% and used to identify potential outlier scans (either absolute or
% scan-to-scan differences in the global BOLD signal normalized to
% z-scores) 
%
% The third graph shows the timeseries derived from the subject motion
% parameters and used to identify potential outlier scans (either the
% absolute or scan-to-scan differences in the individual motion parameters
% -x,y,z translation parameters and roll,pitch,yaw rotation parameters-, or
% the absolute or scan-to-scan differences in a single 'composite motion'
% measure -maximal movement of any voxel within the brain bounding box-)
%
% Using default threshold values for each of the bottom three graphs
% we define outliers as points that exceed the threshold in at least
% one of the global signal or motion parameter graphs. The thresholds are
% shown as horizontal black lines in each of the graphs.
%
% Points which are identified as outliers, are indicated by a vertical
% black line in the graph that corresponds to the outlying
% parameter(s). For example, the if the absolute value of the Y motion
% parameter for time t=17 is above the motion threshold, it is
% identified as an outlier and indicated by a black vertical line at
% t=17 in the third graph. The union of all outliers is indicated by
% red vertical lines on the top graph. The list of outliers is also
% displayed in the editable text box below the graphs.  The current
% values of the thresholds are displayed by the side of the
% corresponding graphs. These values may can be changed by the user
% either by pressing the up/down buttons, which increment/decrement
% the current value by 10%, or by specifying a new value in the text
% box.
%
% In Addition, the user can manually add or remove points from the
% list of outliers by editting the list. Note that the list is only
% updated once the curser points outside the text box (i.e. click the
% mouse somewhere outside the text box). Since any changes made by the
% user are overridden once the thresholds are updated, it is
% recommended to do any manual changes as the last step before saving.
%
% In the 'options' section of the GUI, the user can select to display the
% task-related design, the spectra of the global BOLD signal motion or
% task-related parameters, the level of task-correlated motion and
% task-correlated global BOLD signal, as well as the corrected analysis
% mask (without the influence of the identified outliers)
%
% By default art generates the following output files:
%  Regressor files (one per session, stored in the same folders as the
%   functional volumes, and named art_regression_outliers_*.mat and 
%   art_regression_outliers_and_movement_*.mat). This files can be entered
%   as covariates in the first-level analyses in order to effectively
%   remove the identified outlier scans from further analyses
%  Analysis mask (one file named art_mask.img defining the analysis mask
%   after disregarding outlier scans). This file can be entered as an
%   explicit analysis mask in the first-level analyses in order to
%   avoid any influences of the outlier scans on the implicit analysis mask
%   computation (on SPM you will also need to modify the defaults in order
%   to skip the implicit masking operation, e.g. set defaults.mask.thresh =
%   -inf) 
%
% Pressing the save button lets the user choose wheter to save the
% list of identified outliers, motion statistics, graphs, outlier
% regressors, or Analysis mask. 
%
%
% art('sess_file','filename.cfg');
% Uses the configuration file filename.cfg to define the art analysis
% Options: see example.cfg file for more information
%
% art('sess_file',batch)
% Uses the structure batch to define the art analysis
% Options:
%   batch.P                   : batch.P{nses} [char] functional filename(s) for session nses
%   batch.M                   : batch.M{nses} [char] realignment filename for session nses
%   batch.global_threshold    : global BOLD signal threshold (z-score)
%   batch.motion_threshold    : motion threshold(s)
%   batch.use_diff_motion     : 1/0 use scan-to-scan differences in motion parameters
%   batch.use_diff_global     : 1/0 use scan-to-scan differences in global BOLD signal
%   batch.use_norms           : 1/0 use motion composite measure
%   batch.drop_flag           : number of initial scans to flag as outliers (removal of initial scans)
%   batch.motion_file_type    : indicates type of realignment file (0: SPM rp_*.txt file; 1: FSL .par file; 2: Siemens .txt file; 3: .txt SPM-format but rotation parameters in degrees)
%   batch.close               : 1/0 close gui
%   batch.print               : 1/0 print gui
%   batch.output_dir          : directory for output files (default same folder as first-session functional files)
%
% See also art_batch
%

% ----------------------------------------------------------------------
% - Added voxel-wise SNR and variability displays; save composite motion
%   measure; clean-up code -Alfonso 06/11  
% - Added analysis mask option, modified pow spec display and other GUI display
%   options, last update 9/25/09 - Sue, Alfonso and Darren
% - if multiple sessions are specified, standard deviations are calculated
%   within sessions
%   oliver hinds 2008-04-23
% - added support for reading siemens motion paramter file format
%   oliver hinds 2008-04-23
% - added ability to read file names and sesions from config file
%   oliver hinds 2008-04-23
% - tiny fix to make Matlab 6.5 compatible - Shay 5/14/07
% - added "signal-task correlation" - 5/11/07
% - added "motion-task correlation" and "show spectrum"
% - minor GUI changes to support Windows and open a large graph
%   in a separate window. Also fixed starnge motion params filename
%   bug. Shay Mozes, 5/2/2007
% - fixed bug in display of motion outlier on the graph, Shay Mozes 4/30/2007
% - superimpose task conditions on the z-graph, Shay Mozes, 4/24/2007
% - added support for SPM5, Shay Mozes, 4/2007
% - now supporting FSL .par format, 4/9/2007
% - new GUI and features Shay Mozes, 2006
% + Mar. 2007 from art_global.m, by Paul Mazaika, April 2004.
%   from artdetect4.m, by Jeff Cooper, Nov. 2002
%   from artdetect3.m, by Sue Whitfield artdetect.m
% Sue Whitfield 2000



%% ---------------   GUI initialization  -----------------------------
% GUIDE GUI default creation and callback support
% --------------------------------------------------------------------
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @art_OpeningFcn, ...
    'gui_OutputFcn',  @art_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
end



%% ---------------   Main process -----------------------------
% Main functionality of art (this process runs once for each dataset)
% --------------------------------------------------------------------

% -----------------------------------------------------------------------
% ART_OPENINGFCN
% This function is called before opening the gui.
% It executes the main functionality of art, extracting the global signal
% and movement parameters, defining the measures used to identify outliers,
% printing movement statistics, and calling the Update* functions to
% identify outliers and generate the corresponding plots
% -----------------------------------------------------------------------
function art_OpeningFcn(hObject, eventdata, handles, varargin)
% -----------------------
%Initialize
% -----------------------
warning('OFF','all')
%find spm version
if isdeployed, spm_ver='SPM12';
else spm_ver = spm('ver');
end
switch(spm_ver),
    case 'SPM99', spm_ver=1;
    case 'SPM2', spm_ver=2;
    case 'SPM5', spm_ver=5;
    case {'SPM8b','SPM8'}, spm_ver=8;
    case {'SPM12b','SPM12'}, spm_ver=12;
    otherwise, art_disp(['Warning! unrecognized SPM version ',spm_ver]); spm_ver=12;
end

%clear data from previous sessions
try
    setappdata(handles.showDesign,'SPM',[]);
    setappdata(handles.showDesign,'sessions',[]);
    setappdata(handles.showDesign,'SPMfile',[]);
    setappdata(handles.zthresh,'g',[]);
    setappdata(handles.mvthresh,'mv_data',[]);
    setappdata(handles.mvthresh,'drop_flag',[]);
    setappdata(handles.mvthresh,'altval',[]);
    setappdata(handles.rtthresh,'altval',[]);
    setappdata(handles.savefile,'path',[]);
    setappdata(handles.savefile,'datafiles',[]);
    setappdata(handles.savefile,'stats_file',[]);
    setappdata(handles.savefile,'analyses',[]);
    setappdata(handles.mvthresh,'mv_stats',[]);
    setappdata(handles.zthresh,'zoutliers',[]);
    setappdata(handles.mvthresh,'mv_norm_outliers',[]);
    setappdata(handles.mvthresh,'mv_x_outliers',[]);
    setappdata(handles.mvthresh,'mv_y_outliers',[]);
    setappdata(handles.mvthresh,'mv_z_outliers',[]);
    setappdata(handles.rtthresh,'rt_norm_outliers',[]);
    setappdata(handles.rtthresh,'rt_p_outliers',[]);
    setappdata(handles.rtthresh,'rt_r_outliers',[]);
    setappdata(handles.rtthresh,'rt_y_outliers',[]);
catch %#ok<*CTCH>
end

% ------------------------
% Default values for outliers
% ------------------------
z_thresh = 3.0;             %global signal threshold
mvmt_thresh = 0.5;          %absolute subject motion threshold
rotat_thresh = .05;         %absolute subject rotation threshold
mvmt_diff_thresh = 1.0;     %scan-to-scan subject motion threshold
rotat_diff_thresh = .02;    %scan-to-scan subject rotation threshold

output_path='';
sess_file='';
stats_file='';


% look for args in varargin, 
for i=1:2:numel(varargin)
    if strcmp(varargin{i}, 'sess_file')
        sess_file = varargin{i+1};
    elseif strcmp(varargin{i}, 'stats_file')
        stats_file = varargin{i+1};
    elseif strcmp(varargin{i}, 'output_path')
        output_path = varargin{i+1};
    elseif i==numel(varargin),
        if exist(varargin{i},'file'),
            [filepath,filename,fileext]=fileparts(varargin{i});
            if strcmp(fileext,'.cfg'), sess_file=varargin{i}; end
        end
    end
end
setappdata(handles.savefile,'stats_file',stats_file);


% ------------------------
% Collect files
% ------------------------
if ~isempty(sess_file) % read config
    [num_sess,global_type_flag,drop_flag,gui_display,motionFileType,motion_threshold,global_threshold,use_diff_motion,use_diff_global,use_norms,SPMfile,mask_file,output_dir,do_close,do_print,P,M] = ...
        read_art_sess_file(sess_file); 
else
    motion_threshold=[];global_threshold=[];SPMfile=[];mask_file=[];output_dir='';do_close=false;do_print=false;
    use_diff_motion=1;use_diff_global=1;use_norms=1;
    num_sess = spm_input('How many sessions?',1,'n',1,1);
    
    global_type_flag = spm_input('Which global mean to use?', 1, 'm', ...
        'Regular | User Mask',...
        [1 2], 1);
    motionFileType = spm_input('Select type of motion params file.',1,'m',...
        ' txt(SPM) | par(FSL) | txt(Siemens)', ...
        [0 1 2], 0);
    
    P=cell(1,num_sess);
    M=cell(1,num_sess);
    for i = 1:num_sess
        switch spm_ver
            case {1,2}
                P{i} = spm_get(Inf,'.img',['Select functional volumes for session'  num2str(i) ':']);
            case {5,8,12}
                P{i} = spm_select(Inf,'image',['Select functional volumes for session'  num2str(i) ':']);
                %P{i} = spm_select(Inf,'.*\.nii|.*\.img',['Select functional volumes for session'  num2str(i) ':']);
        end
        if motionFileType == 0 %SPM format
            switch spm_ver
                case {1,2}
                    mvmt_file = spm_get(1,'.txt',['Select movement params file for session' num2str(i) ':']);
                case {5,8,12}
                    mvmt_file = spm_select(1,'^.*\.txt$',['Select movement params file for session' num2str(i) ':']);
            end
            M{i} =load(mvmt_file);
        elseif motionFileType == 1 %FSL format
            switch spm_ver
                case {1,2}
                    mvmt_file = spm_get(1,'.par',['Select movement params file for session' num2str(i) ':']);
                case {5,8,12}
                    mvmt_file = spm_select(1,'^.*\.par$',['Select movement params file for session' num2str(i) ':']);
            end
            M{i} =load(mvmt_file);
        elseif motionFileType == 2 % Siemens MotionDetectionParameter.txt
            switch spm_ver
                case {1,2}
                    mvmt_file = spm_get(1,'.txt',['Select movement params file for session' num2str(i) ':']);
                case {5,8,12}
                    mvmt_file = spm_select(1,'^.*\.txt$',['Select movement params file for session' num2str(i) ':']);
            end
            M{i} = read_siemens_motion_parm_file(mvmt_file);
        end
        output_path = fileparts(mvmt_file);
    end
    drop_flag = 0;
    gui_display = true; % placeholder future use skipping Java functionality (alfnie)
end

if global_type_flag==2,
    if isempty(mask_file),
        mask_file = spm_select(1, '.*\.nii|.*\.img', 'Select mask image in functional space');
    end
    mask=spm_vol(mask_file);
end
setappdata(handles.showDesign,'sessions',-(1:num_sess));
setappdata(handles.showDesign,'SPMfile',SPMfile);
if ~isempty(SPMfile), temp=load(SPMfile); setappdata(handles.showDesign,'SPM',temp.SPM); clear temp; end
datafiles=cell(1,length(P));
for i=1:length(P),
    if ~isempty(P{i}),
        datafiles{i}=strtrim(P{i}(1,:));
    else
        datafiles{i}=[];
    end;
end % <alfnie>: keep filenames of functional data (first scan per session only)
if ~isempty(motion_threshold), 
    mvmt_thresh=motion_threshold(1); mvmt_diff_thresh=motion_threshold(1); 
    if numel(motion_threshold)>1, rotat_thresh=motion_threshold(2); rotat_diff_thresh=motion_threshold(2); end
end
if ~isempty(global_threshold), z_thresh=global_threshold; end
%if drop_flag, disp('warning: explicitly dropping 1st scan no longer supported. Edit the ''all outliers'' box to specify any number of additional outlier scans'); end 
% dropflag extended to add first N scans as outliers in savefile_Callback_SaveRegressor (alfnie 09/15)

mv_data = [];
for i = 1:length(M)
    mv_data = vertcat(mv_data,M{i});
end

%translate to SPM format: x y z (in mm) pitch roll yaw (in radians)
if motionFileType == 1, %FSL order of fields is three Euler angles (x,y,z in radians) then three translation params (x,y,z in mm).
    tmp = mv_data(:,1:3);
    mv_data(:,1:3) = mv_data(:,4:6);
    mv_data(:,4:6) = tmp;
elseif motionFileType == 3, %x y z (in mm) pitch roll jaw (in degrees)
    mv_data(:,4:6)=mv_data(:,4:6)/180*pi;
end


% ---------------------------------------
% Compute Global signal and analysis mask
% ---------------------------------------
g = cell(1,num_sess); % g is a cell array of the global mean for each scan in each session
gsigma=g;gmean=g;dgsigma=g;dgmean=g;
maskscan={};
maskscan_files=0;
art_mask_temporalfile=['art_mask_temporalfile',datestr(now,'yyyymmddHHMMSSFFF'),'.mat'];
Data_Sum=0;
Data_SumSquared=0;
VY=cell(1,num_sess);
cumdisp;
notcoregistered=false;
for sess=1:num_sess
    art_disp('fprintf','%-4s: ',['Mapping files for session ' num2str(sess) '...']);
    VY{sess}     = spm_vol(P{sess});
    art_disp('fprintf','%3s\n','...done')
    
    switch spm_ver
        case {1,2}
            if any(any(diff(cat(1,VY{sess}.dim),1,1),1)&[1,1,1,0])
                error('images do not all have the same dimensions')
            end
        case {5,8,12}
            if any(any(diff(cat(1,VY{sess}.dim),1,1),1))
                error('images do not all have the same dimensions')
            end
    end
    if ~notcoregistered
        for n=1:numel(VY{sess}),
            if any(any(VY{sess}(n).mat~=VY{1}(1).mat)), notcoregistered=true; break; end
        end
    end
    
    % --------------------------------------------------
    % Extracts global-signal mask and (optionally) compute analysis mask
    % note: scan-specific global-signal mask = voxels above mean(all voxels)/8
    %       scan-specific analysis mask = voxels above 0.8*mean(global-signal mask)
    %       global-signal mask = intersection of {scan-specific global-signal mask}
    %       analysis mask = intersection of {scan-specific analysis mask}
    % --------------------------------------------------
    nscans = numel(VY{sess});
    g{sess} = zeros(nscans,4);
    art_disp('fprintf','%-4s: %3s','Calculating globals...',' ')
    if sess==1
        VY1inv=pinv(VY{sess}(1).mat);
        [tempx,tempy,tempz]=ind2sub(VY{sess}(1).dim(1:3),1:prod(VY{sess}(1).dim(1:3)));
        xyz_voxel=[tempx(:),tempy(:),tempz(:),ones(numel(tempx),1)]';
        xyz=VY{sess}(1).mat*xyz_voxel;
    end
    if global_type_flag==1  % regular mean : Global-conjunction (uses conjunction of individual scan masks; individual scan mask are defined as voxels above mean/8 for each scan; see art_maskglobal_scan)
        Mask=ones(VY{sess}(1).dim(1:3));
        for i = 1:nscans,
            if notcoregistered, temp=spm_get_data(VY{sess}(i),pinv(VY{sess}(i).mat)*xyz);
            else temp=reshape(spm_get_data(VY{sess}(i),xyz_voxel),VY{sess}(i).dim); %temp=spm_read_vols(VY{sess}(i));
            end
            [maskscan{end+1},masktemp]=art_maskglobal_scan(temp,VY{sess}(i),VY{sess}(1),VY1inv); %#ok<AGROW>
            if numel(maskscan)>=100, % stop maskcan from growing too much
                maskscan_files=maskscan_files+1;
                try
                    save(fullfile(output_dir,[art_mask_temporalfile(1:end-4),num2str(maskscan_files),'.mat']),'maskscan','-v7.3');
                catch
                    art_disp('warning: unable to write to ',output_dir,' folder. Writing output files to ',pwd,' instead.');
                    output_dir=pwd;
                    save(fullfile(output_dir,[art_mask_temporalfile(1:end-4),num2str(maskscan_files),'.mat']),'maskscan','-v7.3');
                end
                maskscan={};
            end
            Mask(masktemp)=0;
            cumdisp([num2str(i),'/',num2str(nscans)]);
        end
    elseif global_type_flag==2  % user-defined mask
        Mask=spm_get_data(mask,pinv(mask.mat)*xyz);
    end
    idxMask=find(Mask);

    % --------------------------------------------------
    % computes global signal
    % --------------------------------------------------
    badmask=(numel(idxMask)<numel(Mask)/10) && global_type_flag~=2;
    for i = 1:nscans,
        if notcoregistered, temp=reshape(spm_get_data(VY{sess}(i),pinv(VY{sess}(i).mat)*xyz),VY{1}(1).dim);
        else temp=reshape(spm_get_data(VY{sess}(i),xyz_voxel),VY{sess}(i).dim); %temp=spm_read_vols(VY{sess}(i));
        end
        if badmask, g{sess}(i) = spm_global(VY{sess}(i));
        else g{sess}(i)=mean(temp(idxMask));
        end
        Data_Sum=Data_Sum+temp;
        Data_SumSquared=Data_SumSquared+temp.^2;
        %imagesc(Data_Sum(:,:,30));colorbar;drawnow;
    end
    art_disp('fprintf','...done');cumdisp;
    
    % --------------------------------------------------
    % Compute derived signals for outlier identification
    % --------------------------------------------------
    % g{sess} columns:
    % 1: global signal
    % 2: standardized global signal
    % 3: scan-to-scan differences in global signal
    % 4: standardized scan-to-scan differences in global signal
    
    gsigma{sess} = .7413*diff(prctile(g{sess}(:,1),[25,75]));gsigma{sess}(gsigma{sess}==0)=1; % robus standard-deviation
    gmean{sess} = median(g{sess}(:,1)); % robust mean
    g{sess}(:,2)=(g{sess}(:,1)-gmean{sess})/max(eps,gsigma{sess}); % z-score
    g{sess}(2:end,3)=diff(g{sess}(:,1),1,1);
    dgsigma{sess} = .7413*diff(prctile(g{sess}(:,3),[25,75]));dgsigma{sess}(dgsigma{sess}==0)=1;
    dgmean{sess} = median(g{sess}(:,3));
    g{sess}(2:end,4)=(g{sess}(2:end,3)-dgmean{sess})/max(eps,dgsigma{sess});
    z_thresh = 0.1*round(z_thresh*10);
end

VY=cat(1,VY{:}); VY1=VY(1); %#ok<NASGU>
try
    save(fullfile(output_dir,art_mask_temporalfile),'maskscan','maskscan_files','VY1','VY','notcoregistered','xyz','xyz_voxel','Data_Sum','Data_SumSquared','-v7.3'); 
catch
    art_disp(['warning: unable to write to ',output_dir,' folder. Writing output files to ',pwd,' instead.']);
    output_dir=pwd;
    save(fullfile(output_dir,art_mask_temporalfile),'maskscan','maskscan_files','VY1','VY','notcoregistered','xyz','xyz_voxel','Data_Sum','Data_SumSquared','-v7.3'); 
end
set(handles.figure1,'closerequestfcn',['try,if ispc,[nill,ok]=system(''del "',fullfile(output_dir,[art_mask_temporalfile(1:end-4),'"*.mat']),''');else [nill,ok]=system(''rm ''''',fullfile(output_dir,[art_mask_temporalfile(1:end-4),'''''*.mat']),'''); end; end; delete(gcbf);']);

%update text fields
set(handles.data_stdv,'String',num2str(cat(2,gsigma{:}),'%0.1f '));
set(handles.zthresh,'String',num2str(z_thresh));
set(handles.mvthresh,'String',num2str(mvmt_thresh));
set(handles.rtthresh,'String',num2str(rotat_thresh));

% ------------------------------------------------------------------------
% Compute Movement parameters and derived signals for outlier identification
% ------------------------------------------------------------------------
% mv_data columns:
% 1-6   : raw movement parameters
% 7     : euclidean norm of raw movement parameters
% 8-13  : scan-to-scan differences in raw movement parameters
% 14-31 : 6 control points trajectories (placed on center of faces of bounding box; x,y,z coordinates for each control point)
% 32    : composite measure: euclidean norm of control point trajectories
% 33-50 : scan-to-scan differences in 6 control points trajectories
% 51    : composite measure: max scan-to-scan movement across the 6 control points

mv_data=[mv_data,zeros([size(mv_data,1),51-size(mv_data,2)])];
respos=diag([70,70,75]);resneg=diag([-70,-110,-45]);
res=[respos,zeros(3,1),zeros(3,4),zeros(3,4),eye(3),zeros(3,1); % 6 control points: [+x,+y,+z,-x,-y,-z];
    zeros(3,4),respos,zeros(3,1),zeros(3,4),eye(3),zeros(3,1);
    zeros(3,4),zeros(3,4),respos,zeros(3,1),eye(3),zeros(3,1);
    resneg,zeros(3,1),zeros(3,4),zeros(3,4),eye(3),zeros(3,1);
    zeros(3,4),resneg,zeros(3,1),zeros(3,4),eye(3),zeros(3,1);
    zeros(3,4),zeros(3,4),resneg,zeros(3,1),eye(3),zeros(3,1);];
for i=1:size(mv_data,1)
    temp=spm_matrix([1*mv_data(i,1:3),mv_data(i,4:6)]); temp=temp(:)';
    mv_data(i,14:31)=temp*res';
end
cur_sess_start=0;
for sess=1:num_sess
    n=length(g{sess}(:,1));
    mv_data(cur_sess_start+(1:n),7) = sqrt(sum(abs(mv_data(cur_sess_start+(1:n),1:3)).^2,2));
    mv_data(cur_sess_start+(2:n),8:13)  = diff(mv_data(cur_sess_start+(1:n),1:6),1,1);
    mv_data(cur_sess_start+(1:n),32)=sqrt(mean(abs(detrend(mv_data(cur_sess_start+(1:n),14:31),'constant')).^2,2));
    mv_data(cur_sess_start+(2:n),33:50)=diff(mv_data(cur_sess_start+(1:n),14:31),1,1);
    mv_data(cur_sess_start+(2:n),51)=max(sqrt(sum(reshape(abs(mv_data(cur_sess_start+(2:n),33:50)).^2,[n-1,3,6]),2)),[],3);
    cur_sess_start = cur_sess_start + n;
end

%save application data for use in callbacks
setappdata(handles.zthresh,'g',g);
setappdata(handles.mvthresh,'mv_data',mv_data);
setappdata(handles.mvthresh,'altval',num2str(mvmt_diff_thresh));
setappdata(handles.rtthresh,'altval',num2str(rotat_diff_thresh));
setappdata(handles.mvthresh,'drop_flag',drop_flag);
setappdata(handles.savefile,'path',output_path);
setappdata(handles.savefile,'dir',output_dir);
setappdata(handles.savefile,'datafiles',datafiles);
setappdata(handles.savefile,'spm_ver',spm_ver);
setappdata(handles.savefile,'mv_data_raw',M);
setappdata(handles.savefile,'art_mask_temporalfile',art_mask_temporalfile);

if ~isempty(SPMfile), set(handles.showDesign,'Value',get(handles.showDesign,'Max')); end
set(handles.norms,'Value',use_norms);
set(handles.diff1,'value',use_diff_global);
set(handles.diff2,'value',use_diff_motion);
if use_diff_global||use_diff_motion
    diffsglobalandmotion_Callback(hObject, [], handles, 2*use_diff_global-1, 2*use_diff_motion-1);
else
    UpdateGlobal(hObject, eventdata, handles,1.0);
    UpdateMovement(hObject, eventdata, handles,1.0);
    UpdateRotation(hObject, eventdata, handles,1.0);
    UnionOutliers(hObject, eventdata, handles);
    UpdateSummaryplot(hObject, eventdata, handles);
end

idx=str2num(get(handles.all_outliers, 'String')); %#ok<*ST2NM>
for sess=1:num_sess
    art_disp('fprintf','\nSession %d global statistics -  mean: %7.4f stdv: %7.4f',sess,gmean{sess},gsigma{sess});
end
art_disp('fprintf','\n');
art_disp('fprintf','Outlier detection: %d identified outliers\n',length(idx));

% ------------------------
% Compute and print statistics of movement
%--------------------------
mv_data = getappdata(handles.mvthresh,'mv_data');
mv_stats = [mean(abs(mv_data)); std(abs(mv_data)); max(abs(mv_data)) ];
setappdata(handles.mvthresh,'mv_stats',mv_stats);
art_disp('fprintf','\n\nStatistics of movement data:\n\n');
art_disp('fprintf','%5s%10s%10s%10s%11s%10s%9s%10s\n',' ','x','y','z',' pitch','roll','yaw','norm');
art_disp('fprintf','%7s%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n','mean ',mv_stats(1,1:3),mv_stats(1,4:6),mv_stats(1,32));
art_disp('fprintf','%7s%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n','stdv ',mv_stats(2,1:3),mv_stats(2,4:6),mv_stats(2,32));
art_disp('fprintf','%7s%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n\n\n','max ',mv_stats(3,1:3),mv_stats(3,4:6),mv_stats(3,32));
% BEGIN ohinds 2008-04-23: save stats to file
if ~isempty(stats_file)
    if isempty(fileparts(stats_file)),stats_file=fullfile(output_dir,stats_file);end
    fp = fopen(stats_file,'wt');
    if fp ~= -1
        art_disp('fprintf','saving global motion stats to %s\n',stats_file);
        fprintf(fp,'%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n',mv_stats(1,1:3),mv_stats(1,4:6),mv_stats(1,32));
        fprintf(fp,'%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n',mv_stats(2,1:3),mv_stats(2,4:6),mv_stats(2,32));
        fprintf(fp,'%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n\n\n',mv_stats(3,1:3),mv_stats(3,4:6),mv_stats(3,32));
        fclose(fp);
    end
end
% END ohinds 2008-04-23: save stats to file

% Choose default command line output for art
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);
set(handles.figure1,'resizefcn','art(''showOptions_Callback'',gcbo,[],guidata(gcbo))');

% --------------------------------------------------------
% Saves regressor matrix with outliers & new analysis mask
%---------------------------------------------------------
savefile_Callback_SaveRegressor(handles);
try
savefile_Callback_SaveMask(handles);
catch
    art_disp('Warning. Error encountered during ART implicit-mask file creation. Skipping this step');
end
try % fix Matlab14 new graphic objects issue (uipanel shown on top of uicontrol objects)
    h1=findobj(handles.figure1,'type','uipanel');
    h2=findobj(handles.figure1,'type','uicontrol');
    x1=cell2mat(get(h1,'position'));
    x2=cell2mat(get(h2,'position'));
    for n=1:size(x2,1)
        xpos=repmat(x2(n,1:2),[size(x1,1),1]);
        [md,n0]=max(sign(min(min(xpos-x1(:,1:2),x1(:,1:2)+x1(:,3:4)-xpos),[],2)),[],1);
        if md>0, set(h2(n),'parent',h1(n0)); end
    end
    set(h2,'parent',handles.figure1);
end
try, if do_print, saveas(handles.figure1,fullfile(output_dir,'art_screenshot.fig')); end; end
try, if do_close, close(handles.figure1); end; end

end




%% ---------------   GUI callback functions --------------------------
% Processing steps for each of the art GUI options 
% --------------------------------------------------------------------

% -----------------------------------------------------------------------
% UPDATEGLOBAL
% This function identifies outliers based on the BOLD globalsignal and
% generates the corresponding plot (gui 2nd plot from the top)
% -----------------------------------------------------------------------
function UpdateGlobal(hObject, eventdata, handles, incr) %#ok<INUSL>

%get data
z_thresh = str2num(get(handles.zthresh,'String'));
g = getappdata(handles.zthresh,'g');
num_sess = length(g);

%calc new outliers
% BEGIN ohinds 2008-04-23: plot zscores
%axes(handles.zvalue); % (avoids turning figure visible unnecessarily; alfnie 2017)
cla(handles.zvalue);
hold(handles.zvalue,'on');
cur_sess_start=1;
z_thresh = z_thresh*incr; 
idxind=2;
out_idx = cell(1,num_sess);
for sess=1:num_sess
    if get(handles.diff1,'value'),
        out_idx{sess} = cur_sess_start+(find(abs(g{sess}(:,idxind)) > z_thresh|abs([g{sess}(2:end,idxind);0]) > z_thresh))'-1; 
    else
        out_idx{sess} = cur_sess_start+(find(abs(g{sess}(:,idxind)) > z_thresh))'-1; 
    end
    %update plot
    plot(cur_sess_start:cur_sess_start+size(g{sess},1)-1, g{sess}(:,idxind),'parent',handles.zvalue);
    cur_sess_start = cur_sess_start + size(g{sess},1);
end
out_idx=cat(2,out_idx{:});
set(handles.zvalue,'xlim',[0,cur_sess_start]);
l=ylabel(handles.zvalue,'global mean\newline     [std]');
set(l,'VerticalAlignment','bottom','horizontalalignment','center');
set(handles.zvalue,'XTickLabel',[]);

thresh_x = 1:cur_sess_start-1;
thresh_y = z_thresh*ones(1,length(thresh_x));
line(thresh_x, thresh_y, 'Color', 'black','parent',handles.zvalue);
line(thresh_x, -1*thresh_y, 'Color', 'black','parent',handles.zvalue);

%update text
set(handles.zthresh,'String',num2str(z_thresh));
setappdata(handles.zthresh,'zoutliers',out_idx);
hold(handles.zvalue,'off');
% END ohinds 2008-04-23: plot zscores

%plot outliers
axes_lim = get(handles.zvalue, 'YLim');
axes_height = axes_lim;
for i = 1:length(out_idx)
    line((out_idx(i)*ones(1, length(axes_height))), axes_height, 'Color', 'black','parent',handles.zvalue);
end
end


% -----------------------------------------------------------------------
% UPDATEMOVEMENT
% This function identifies outliers based on the subject movement
% parameters (either translation parameters or composite motion) and
% generates the corresponding plot (gui 3rd plot from the top) 
% -----------------------------------------------------------------------
function UpdateMovement(hObject, eventdata, handles, incr) %#ok<INUSL>

%get data
mvmt_thresh = str2num(get(handles.mvthresh,'String'));
mv_data = getappdata(handles.mvthresh,'mv_data');

%calc new outliers
mvmt_thresh = mvmt_thresh*incr;

if get(handles.diff2,'value'),
    out_mvmt_idx = (find(abs(mv_data(:,1:3)) > mvmt_thresh | [abs(mv_data(2:end,1:3));[0,0,0]] > mvmt_thresh ))';
else
    out_mvmt_idx = (find(abs(mv_data(:,1:3)) > mvmt_thresh))';
end
out_mvmt_idx_X=out_mvmt_idx(out_mvmt_idx<=size(mv_data,1));
out_mvmt_idx_Y=out_mvmt_idx(out_mvmt_idx>size(mv_data,1)&out_mvmt_idx<=2*size(mv_data,1))-size(mv_data,1);
out_mvmt_idx_Z=out_mvmt_idx(out_mvmt_idx>2*size(mv_data,1))-2*size(mv_data,1);

%find norm outliers
if (get(handles.norms,'Value') == get(handles.norms,'Max'))
    normv=mv_data(:,32);
    if get(handles.diff2,'value')
        out_mvmt_idx_norm = find(normv>mvmt_thresh|[normv(2:end);0]>mvmt_thresh);
    else
        out_mvmt_idx_norm = find(normv>mvmt_thresh);
    end
    out_mvmt_idx_norm = out_mvmt_idx_norm';
    setappdata(handles.mvthresh,'mv_norm_outliers',out_mvmt_idx_norm);
end

%update text
set(handles.mvthresh,'String',num2str(mvmt_thresh));
setappdata(handles.mvthresh,'mv_x_outliers',out_mvmt_idx_X);
setappdata(handles.mvthresh,'mv_y_outliers',out_mvmt_idx_Y);
setappdata(handles.mvthresh,'mv_z_outliers',out_mvmt_idx_Z);

%axes(handles.mvmtGraph); % (avoids turning figure visible unnecessarily; alfnie 2017)
cla(handles.mvmtGraph);
if (get(handles.norms,'Value') == get(handles.norms,'Max'))
    plot(normv,'parent',handles.mvmtGraph);
    set(handles.mvmtGraph,'xlim',[0,length(normv)+1]);
else
    plot(mv_data(:,1:3),'parent',handles.mvmtGraph);
    set(handles.mvmtGraph,'xlim',[0,size(mv_data,1)+1]);
end
l=ylabel(handles.mvmtGraph,'movement \newline   [mm]');
set(l,'VerticalAlignment','bottom','horizontalalignment','center');
set(handles.mvmtGraph,'XTickLabel',[]);
if (get(handles.norms,'Value') ~= get(handles.norms,'Max'))
    legend(handles.mvmtGraph,'x', 'y', 'z','Location','East');
end
h = handles.mvmtGraph;
set(h,'Ygrid','on');

thresh_mv_x = 1:size(mv_data,1);
thresh_mv_y = mvmt_thresh*ones(1,size(mv_data,1));
line(thresh_mv_x, thresh_mv_y, 'Color', 'black','parent',handles.mvmtGraph);
if ~(get(handles.norms,'Value') == get(handles.norms,'Max'))
    line(thresh_mv_x, -1*thresh_mv_y, 'Color', 'black','parent',handles.mvmtGraph);
end

axes_lim = get(handles.mvmtGraph, 'YLim');
axes_height = axes_lim;
if ~(get(handles.norms,'Value') == get(handles.norms,'Max'))
    for i = 1:length(out_mvmt_idx_X)
        line((out_mvmt_idx_X(i)*ones(1, length(axes_height))), axes_height, 'Color', 'black','parent',handles.mvmtGraph);
    end
    for i = 1:length(out_mvmt_idx_Y)
        line((out_mvmt_idx_Y(i)*ones(1, length(axes_height))), axes_height, 'Color', 'black','parent',handles.mvmtGraph);
    end
    for i = 1:length(out_mvmt_idx_Z)
        line((out_mvmt_idx_Z(i)*ones(1, length(axes_height))), axes_height, 'Color', 'black','parent',handles.mvmtGraph);
    end
else
    for i = 1:length(out_mvmt_idx_norm)
        line((out_mvmt_idx_norm(i)*ones(1, length(axes_height))), axes_height, 'Color', 'black','parent',handles.mvmtGraph);
    end
end
end


% -----------------------------------------------------------------------
% UPDATEROTATION
% This function identifies outliers based on the subject rotation
% parameters and generates the corresponding plot (gui 4th plot from the top)  
% (note: only available when not using composite movement measures)
% -----------------------------------------------------------------------
function UpdateRotation(hObject, eventdata, handles, incr) %#ok<INUSL>

%get data
rotat_thresh = str2num(get(handles.rtthresh,'String'));
mv_data = getappdata(handles.mvthresh,'mv_data');

%calc new outliers
rotat_thresh = rotat_thresh*incr;

out_rotat_idx = (find(abs(mv_data(:,4:6)) > rotat_thresh))';
out_rotat_idx_p=out_rotat_idx(out_rotat_idx<=size(mv_data,1));
out_rotat_idx_r=out_rotat_idx(out_rotat_idx>size(mv_data,1)&out_rotat_idx<=2*size(mv_data,1))-size(mv_data,1);
out_rotat_idx_y=out_rotat_idx(out_rotat_idx>2*size(mv_data,1))-2*size(mv_data,1);
%update text
set(handles.rtthresh,'String',num2str(rotat_thresh));
setappdata(handles.rtthresh,'rt_p_outliers',out_rotat_idx_p);
setappdata(handles.rtthresh,'rt_r_outliers',out_rotat_idx_r);
setappdata(handles.rtthresh,'rt_y_outliers',out_rotat_idx_y);

%if composite measure
if (get(handles.norms,'Value') == get(handles.norms,'Max'))
    setappdata(handles.rtthresh,'rt_norm_outliers',[]);
    legend(handles.rotatGraph,'off');
    cla(handles.rotatGraph);
    set([handles.rotatGraph,handles.rt_up,handles.rt_down,handles.rtthresh,handles.text13],'visible','off')
    set(handles.axes_mask,'position',get(handles.axes_mask,'position').*[1,1,1,0]+[0,0,0,.40]);
    set(handles.all_outliers,'position',get(handles.all_outliers,'position').*[1,1,1,0]+[0,0,0,.29]);
    return
end
set(handles.axes_mask,'position',get(handles.axes_mask,'position').*[1,1,1,0]+[0,0,0,.25]);
set(handles.all_outliers,'position',get(handles.all_outliers,'position').*[1,1,1,0]+[0,0,0,.14]);

%axes(handles.rotatGraph); % (avoids turning figure visible unnecessarily; alfnie 2017)
cla(handles.rotatGraph);
plot(mv_data(:,4:6),'parent',handles.rotatGraph);
set(handles.rotatGraph,'xlim',[0,size(mv_data,1)+1]);
l=ylabel(handles.rotatGraph,'rotation \newline  [rad]');
set(l,'VerticalAlignment','Bottom');
set(handles.rotatGraph,'XTickLabel',[]);
legend(handles.rotatGraph,'pitch', 'roll', 'yaw', 'Location', 'East');
h = handles.rotatGraph;
set(h,'Ygrid','on');

thresh_rt_x = 1:length(mv_data);
thresh_rt_y = rotat_thresh*ones(1,length(mv_data));
y_lim = get(handles.rotatGraph, 'YLim');
line(thresh_rt_x, thresh_rt_y, 'Color', 'black','parent',handles.rotatGraph);
line(thresh_rt_x, -1*thresh_rt_y, 'Color', 'black','parent',handles.rotatGraph);
for i = 1:length(out_rotat_idx_p)
    line((out_rotat_idx_p(i)*ones(1, 2)), y_lim, 'Color', 'black','parent',handles.rotatGraph);
end
for i = 1:length(out_rotat_idx_r)
    line((out_rotat_idx_r(i)*ones(1, 2)), y_lim, 'Color', 'black','parent',handles.rotatGraph);
end
for i = 1:length(out_rotat_idx_y)
    line((out_rotat_idx_y(i)*ones(1, 2)), y_lim, 'Color', 'black','parent',handles.rotatGraph);
end
set([h,handles.rt_up,handles.rt_down,handles.rtthresh,handles.text13],'visible','on')
end


% -----------------------------------------------------------------------
% UPDATESUMMARYPLOT
% This function generates the summary plot (gui 1st plot from the top)
% displaying the global signal and all of the identified outliers, and
% optionally the design matrix information and corresponding breakdown of
% outliers by condition
% -----------------------------------------------------------------------
function UpdateSummaryplot(hObject, eventdata, handles)
g = getappdata(handles.zthresh,'g');
num_sess = length(g);
tmps = get(handles.all_outliers,'String');
if ~isempty(tmps)
    nstrings = size(tmps,1);
    idx=cell(1,nstrings);
    for i=1:nstrings
        idx{i}=round(str2num(tmps(i,:)));
    end
    idx=cat(2,idx{:});
    set(handles.all_outliers, 'String', int2str(idx));
else
    idx = [];
end

%plot global mean
%axes(handles.globalMean); % (avoids turning figure visible unnecessarily; alfnie 2017)
cla(handles.globalMean);
% BEGIN ohinds 2008-04-23: plot and print global mean
hold(handles.globalMean,'on');

cur_sess_start=1;
rng_mean=0;rng_minmax=[-inf,-inf];
for sess=1:num_sess
    %rng{sess} = range(g{sess}(:,1));
    rng_mean=rng_mean+mean(g{sess}(:,1));
    rng_minmax=max(rng_minmax,[-min(g{sess}(:,1)),max(g{sess}(:,1))]);
    plot(cur_sess_start:cur_sess_start+length(g{sess}(:,1))-1, g{sess}(:,1),'parent',handles.globalMean);
    %ylabstr = sprintf('%s %f (%d)', ylabstr, rng{sess}, sess); % ohinds: can't put the range for all sessions on the ylabel, not enough room
    cur_sess_start = cur_sess_start + length(g{sess}(:,1));
end
rng_mean=rng_mean/num_sess;
set(handles.globalMean,'xlim',[0,cur_sess_start],'ylim',sort((rng_minmax.*[-1 1])*[1.1,-.1;-.1,1.1])+[0 eps]);
ylabel(handles.globalMean,'mean image\newlineintensity');
xlabel(handles.globalMean,'scans');
% END ohinds 2008-04-23: plot global mean

y_lim = get(handles.globalMean, 'YLim');
cur_sess_start=1;
for sess=1:num_sess
    patch(cur_sess_start+[0,0,(length(g{sess}(:,1))-1)*[1,1]],[ylim,fliplr(ylim)],-ones(1,4),.9+.05*rem(sess,2)*[1,1,1],'edgecolor','none','parent',handles.globalMean);
    if num_sess>1
        text(cur_sess_start+(length(g{sess}(:,1))-1)/2,ylim*[-.1;1.1],['Session ',num2str(sess)],'horizontalalignment','center','parent',handles.globalMean);
    end
    cur_sess_start = cur_sess_start + length(g{sess}(:,1));
end
for i = 1:length(idx)
    line((idx(i)*ones(1, 2)), y_lim, 'Color', 'red','parent',handles.globalMean);
end
hold(handles.globalMean,'off');
analyses=getappdata(handles.savefile,'analyses');
analyses.outliers.scans=idx;
setappdata(handles.savefile,'analyses',analyses);

%show design (moved to global plot <alfnie> 2009-01)
if (get(handles.showDesign,'Value') == get(handles.showDesign,'Max'))
    [SPM,design,names] = get_design(handles);
    stats_file=getappdata(handles.showDesign,'SPMfile');
    hold(handles.globalMean,'on');
    colors = {'k:','b:','r:','g:','c:','m:','y:'};
    h=plot(1,nan,'.','markersize',1,'parent',handles.globalMean);
    for i=1:size(design,2)
        h(i+1)=plot(1:size(design,1) , rng_mean+sum(rng_minmax)/2*design(:,i),colors{mod(i,5)+1},'MarkerSize',4,'parent',handles.globalMean);
    end
    % computes number of outliers per condition <alfnie> 2009-01
    out_idx=round(idx(idx>0));
    if cur_sess_start-1~=size(design,1),
        art_disp(['warning: incorrect number of scans (design matrix: ',num2str(size(design,1)),' ; functional data: ',num2str(cur_sess_start-1),')']);
        outliers_per_condition=length(out_idx);
    else
        outliers_per_condition=[length(out_idx),sum(abs(design(out_idx,:)>0),1); size(design,1),sum(abs(design(:,:)>0),1)];
    end
    if size(outliers_per_condition,2)==length(names)+1,
        legendnames={['Total :',num2str(outliers_per_condition(1,1)),' outlier scans (',num2str(100*outliers_per_condition(1,1)/max(eps,outliers_per_condition(2,1)),'%0.0f'),'%)']};
        for i=1:length(names), legendnames{i+1}=[names{i},' :',num2str(outliers_per_condition(1,i+1),'%0.0f'),' outlier scans (',num2str(100*outliers_per_condition(1,i+1)/max(eps,outliers_per_condition(2,i+1)),'%0.0f'),'%)']; end
        legend(h,legendnames{:});
    end
    analyses=getappdata(handles.savefile,'analyses');
    analyses.outliers.condition_effects=outliers_per_condition(1,2:end);
    analyses.outliers.condition_names=names;
    setappdata(handles.savefile,'analyses',analyses);
    [statsfile_path,statsfile_name] = fileparts(stats_file); if isempty(statsfile_path),statsfile_path=pwd;end;
    stats_file_outliers=fullfile(statsfile_path,[statsfile_name,'_outliers.txt']);
    if ~isempty(stats_file_outliers)
        art_disp('fprintf','Number of outliers\n');
        art_disp('fprintf','%10s ','Total');
        art_disp('fprintf','%10s ',names{:});
        art_disp('fprintf','\n');
        art_disp('fprintf','%10.0f ',outliers_per_condition(1,:));
        art_disp('fprintf','\n');
        art_disp('fprintf',' %9.1f%%',100*outliers_per_condition(1,:)./max(eps,outliers_per_condition(2,:)));
        art_disp('fprintf','\n');
        fp=fopen(stats_file_outliers,'w');
        art_disp('fprintf','saving outlier statistics to %s\n',stats_file_outliers);
        fprintf(fp,'%10s ','Total');
        fprintf(fp,'%10s ',names{:});
        fprintf(fp,'\n');
        fprintf(fp,'%10.0f ',outliers_per_condition(1,:));
        fprintf(fp,'\n');
        fprintf(fp,' %9.1f%%',100*outliers_per_condition(1,:)./max(eps,outliers_per_condition(2,:)));
        fprintf(fp,'\n');
        fclose(fp); 
    end
    hold(handles.globalMean,'off')
    
else
    legend(handles.globalMean,'off');
end
if (get(handles.showOptions,'Value') > 1)
    showOptions_Callback(hObject, eventdata, handles);
else
    showOptions_Callback(hObject, eventdata, []);
end
end


% -----------------------------------------------------------------------
% UNIONOUTLIERS
% This function retrieves all of the outliers identified from the global
% signal and movement parameters and updates the list of all outliers
% -----------------------------------------------------------------------
function UnionOutliers(hObject, eventdata, handles) %#ok<INUSL>

%get data
idx = getappdata(handles.zthresh,'zoutliers');
if ~(get(handles.norms,'Value') == get(handles.norms,'Max'))
    idx = [idx , getappdata(handles.mvthresh,'mv_x_outliers')];
    idx = [idx , getappdata(handles.mvthresh,'mv_y_outliers')];
    idx = [idx , getappdata(handles.mvthresh,'mv_z_outliers')];
    idx = [idx , getappdata(handles.rtthresh,'rt_p_outliers')];
    idx = [idx , getappdata(handles.rtthresh,'rt_r_outliers')];
    idx = [idx , getappdata(handles.rtthresh,'rt_y_outliers')];
else
    idx = [idx , getappdata(handles.rtthresh,'rt_norm_outliers')];
    idx = [idx , getappdata(handles.mvthresh,'mv_norm_outliers')];
end
idx = unique(idx);
%update data
set(handles.all_outliers, 'String', int2str(idx));
end


% -----------------------------------------------------------------------
% DIFFSGLOBALANDMOTION_CALLBACK
% This function executes when any of the 'Use diff' checkboxes are toggled
% It toggles between 'absolute' and 'scan-to-scan' measures for any of the
% global signal or subject movement parameters, and update the
% corresponding plots
% -----------------------------------------------------------------------
function diffsglobalandmotion_Callback(hObject, eventdata, handles,option_global,option_motion)
if nargin<4||isempty(option_global), option_global=0; end
if nargin<5||isempty(option_motion), option_motion=0; end

if option_motion>0 % diff_motion
    %switch thresholds
    tmp = get(handles.mvthresh,'String');
    set(handles.mvthresh,'String',getappdata(handles.mvthresh,'altval'));
    setappdata(handles.mvthresh,'altval',tmp);
    tmp = get(handles.rtthresh,'String');
    set(handles.rtthresh,'String',getappdata(handles.rtthresh,'altval'));
    setappdata(handles.rtthresh,'altval',tmp);
    
    %switch data used
    mv_data = getappdata(handles.mvthresh,'mv_data');
    tmp = mv_data(:,1:6);
    mv_data(:,1:6) = mv_data(:,8:13);
    mv_data(:,8:13) = tmp;
    tmp = mv_data(:,14:32);
    mv_data(:,14:32) = mv_data(:,33:51);
    mv_data(:,33:51) = tmp;
    setappdata(handles.mvthresh,'mv_data',mv_data);
end

if option_global>0 %diff_global
    g = getappdata(handles.zthresh,'g');
    for n1=1:length(g), g{n1}(:,[2,4])=g{n1}(:,[4,2]); end
    setappdata(handles.zthresh,'g',g);
end

if option_global
    UpdateGlobal(hObject, eventdata, handles,1.0);
end
if option_motion
    UpdateMovement(hObject, eventdata, handles,1.0);
    UpdateRotation(hObject, eventdata, handles,1.0);
end
UnionOutliers(hObject, eventdata, handles);
UpdateSummaryplot(hObject, eventdata, handles);
end


% -----------------------------------------------------------------------
% SHOWCORR_CALLBACK
% Computes and displays the movement-task correlations
% -----------------------------------------------------------------------
function showCorr_Callback(hObject, eventdata, handles) %#ok<INUSL>
if (get(handles.showCorr,'Value') == get(handles.showCorr,'Max'))
    %display correlations
    [SPM,design,names] = get_design(handles);
    mv_data = getappdata(handles.mvthresh,'mv_data');
    sessions = getappdata(handles.showDesign,'sessions');
    f = figure;
    setappdata(handles.showCorr,'figure',f);
    nrows=zeros(1,length(sessions)+1);
    cm=cell(1,length(sessions));
    for sess=1:length(sessions),
        s = sessions(sess);
        rows = SPM.Sess(s).row;
        cols = SPM.Sess(s).col(1:length(SPM.Sess(s).U)); % extracts only effects of interest (no covariates) from design matrix
        nrows(sess+1)=length(rows);
        
        %create partial matrix to correlate (we only want to correlate with the motion parameters within each session). NOTE: This may cause weird behaviour in weird designs... (note: author please clarify)
        part = [SPM.xX.X(rows,cols) mv_data(sum(nrows(1:sess))+(1:nrows(sess+1)),1:6)];
        cm{sess} = corrcoef(part);
        a = subplot(length(sessions),1,sess);
        
        imagesc(cm{sess}(1:end-6,end-5:end),[-1,1]);
        colorbar;
        set(a,'XTickLabel',{'x','y','z','pitch','roll','yaw'});
        set(a,'YTick',1:length(cols));
        set(a,'YTickLabel',names);
        title(sprintf('Session %d',s));
        analyses=getappdata(handles.savefile,'analyses');
        analyses.motion_task_correlation(sess)=struct('r',cm{sess}(1:end-6,end-5:end),'rows',{get(a,'yticklabel')},'cols',{get(a,'xticklabel')});
        setappdata(handles.savefile,'analyses',analyses);
        
    end
else
    f = getappdata(handles.showCorr,'figure');
    if ishandle(f)
        close(f);
    end
end
end

% -----------------------------------------------------------------------
% SHOW_SIGNAL_CORR_CALLBACK
% Computes and displays the BOLD signal-task correlations
% -----------------------------------------------------------------------
function show_signal_corr_Callback(hObject, eventdata, handles) %#ok<INUSL>
if (get(handles.sigCorr,'Value') == get(handles.sigCorr,'Max'))
    %display correlations
    SPM = get_design(handles);
    g = getappdata(handles.zthresh,'g');
    sessions = getappdata(handles.showDesign,'sessions');
    f = figure;
    setappdata(handles.sigCorr,'figure',f);
    cm=cell(1,length(sessions));
    for sess=1:length(sessions)
        s = sessions(sess);
        rows = SPM.Sess(s).row;
        cols = SPM.Sess(s).col(1:length(SPM.Sess(s).U)); % extracts only effects of interest (no covariates) from design matrix
        
        %create partial matrix to correlate (we only want to correlate with the motion parameters within each session). NOTE: This may cause weird behaviour in weird designs... (note: author please clarify)
        part = [SPM.xX.X(rows,cols) g{sess}(:,1)];
        cm{sess} = corrcoef(part);
        a = subplot(length(sessions),1,sess);
        
        imagesc(cm{sess}(end:end,1:end-1),[-1,1]);
        colorbar;
        names=cat(2,SPM.Sess(s).U(1:length(cols)).name);
        set(a,'XTick',1:length(cols));
        set(a,'XTickLabel',names);
        set(a,'YTick',1);
        set(a,'YTickLabel','mean activation');
        title(sprintf('Session %d',s));
        analyses=getappdata(handles.savefile,'analyses');
        analyses.signal_task_correlation(sess)=struct('r',cm{sess}(end,1:end-1),'rows',{get(a,'yticklabel')},'cols',{get(a,'xticklabel')});
        setappdata(handles.savefile,'analyses',analyses);
    end
else
    f = getappdata(handles.sigCorr,'figure');
    if ishandle(f)
        close(f);
    end
end
end

% -----------------------------------------------------------------------
% SHOWSPEC_CALLBACK
% Computes and displays the spectrum of BOLD signal, movement, and task
% regressors
% -----------------------------------------------------------------------
function showSpec_Callback(hObject, eventdata, handles) %#ok<INUSL>
if (get(handles.showSpec,'Value') == get(handles.showSpec,'Max'))
    %get data and compute power spectrum
    SPM = get_design(handles);
    mv_data = getappdata(handles.mvthresh,'mv_data');
    g = getappdata(handles.zthresh,'g');
    sessions = getappdata(handles.showDesign,'sessions');
    f = figure;
    setappdata(handles.showSpec,'figure',f);
    
    %sampling freq.
    sf = 1/SPM.xY.RT;
    
    nrows=zeros(1,length(sessions)+1);
    for sess=1:length(sessions),
        s = sessions(sess);
        rows = SPM.Sess(s).row;
        cols = SPM.Sess(s).col(1:length(SPM.Sess(s).U)); % extracts only effects of interest (no covariates) from design matrix
        nrows(sess+1)=length(rows);
        
        %create partial design matrix which only contains relevant data
        %for curent session
        data = [SPM.xX.X(rows,cols),mv_data(sum(nrows(1:sess))+(1:nrows(sess+1)),1:6),g{sess}(:,1)]; %alfnie 08/2009: added global signal
        
        cf = sf/2; %Nyquist freq.
        n = size(data,1);
        
        %this is done in a loop (and not in matrix ops) since dct encounters memory problems for large matrices.
        hold on
        n=n*5;freqs = (0:cf/n:cf-cf/n)';f = zeros(5*size(data,1),size(data,2));F=f;%alfnie 08/2009: resample
        for i=1:size(data,2)
            %%calculate dct
            temp=(abs(fft(detrend(data(:,i)).*hanning(size(data,1)),2*5*size(data,1))).^2);%alfnie 08/2009: plot spectral densities
            f(:,i) = temp(1:end/2);
            %normalize
            F(:,i) = f(:,i)/(sum(abs(f(:,i))./freqs([2,2:end])));F(1,i)=nan; % normalization (plots show same area in log-freq space)
        end
        
        a = subplot(length(sessions), 1, sess);
        hs = plot(freqs,F,'-'); axis tight; set(gca,'xlim',[sf/size(data,1),cf],'xscale','log','yscale','lin');set(gcf,'color','w');%.9412*[1,1,1])
        names=cat(2,SPM.Sess(s).U(1:length(cols)).name);
        names(end+1:end+7) = {'x','y','z','pitch','roll','yaw','BOLD'};
        l = legend(names,'Location','EastOutside');
        pos = get(l,'Position');
        set(l,'Visible','off');
        
        
        for i = 1:length(names)
            color = get(hs(i),'Color');
            box = uicontrol('Style','checkbox','String',names(i),'ForegroundColor', color,'backgroundcolor','w','Callback',{@showSpec_Callback_setvisibility},'Value',1,'UserData',hs(i),'Units','normalized');
            tmppos = get(box,'Position');
            tmppos(1:2) = pos(1:2);
            tmppos(3)=1-tmppos(1);
            set(box,'Position',tmppos);
            pos(2) = pos(2)+ 0.035;
            if (i > length(names)-7)&&i<length(names)
                set(box,'Value',0);
                showSpec_Callback_setvisibility(box,0);
            end
        end
        
        title(sprintf('Session %d',s));
        ylabel('Power density (normalized)');
        %try finding highpass freq. in SPM.xX.K(s).Hparam
        try
            cutoff = 1/SPM.xX.K(s).HParam;
        catch
            art_disp('fprintf','no highpass cutoff frequency found in SPM.mat, using default (128).\n');
            cutoff = 1/128;
        end
        %draw cutoff freq. DRG (2009-08-25) added to show cutoff frequency more clearly
        x_lim = xlim(a);
        l = patch([x_lim(1) x_lim(1) cutoff cutoff],[ylim(a) fliplr(ylim(a))],-ones(1,4),[.8 .8 .8],'edgecolor','none');
        xlabel(sprintf('Frequency [Hz], cutoff=1/%i',1/cutoff));%DRG (2009-08-25) This line was moved from ~15 lines up.
        set(a,'UserData',l);
        
        analyses=getappdata(handles.savefile,'analyses');
        analyses.motion_task_spectra(sess)=struct('Power',f,'rows',freqs,'cols',{names});
        setappdata(handles.savefile,'analyses',analyses);
        
    end
    
else
    f = getappdata(handles.showSpec,'figure');
    if ishandle(f)
        close(f);
    end
end
end

%toggle visibility for spectrum graph
function showSpec_Callback_setvisibility(handle,tmp) %#ok<INUSD>
h = get(handle,'UserData');
if (get(handle,'Value') == get(handle,'Max'))
    set(h(1),'Visible','on');
else
    set(h(1),'Visible','off');
end
% DRG (2009-08-25) DON'T NEED TO DO THIS. LIMITS ARE ALWAYS THIS SAME ANYWAYS BECAUSE DATA IS THERE, BUT JUST HIDDEN. redraw the cutoff: a = ancestor(h(1),'axes');l = get(a,'UserData');set(l,'Visible','off');ylim = get(a,'YLim');set(l,'YData',ylim,'Visible','on');
end

% -----------------------------------------------------------------------
% SHOWOPTIONS_CALLBACK
% Computes and displays the analysis mask, voxel-wise variance and SNR,
% based on current list of outliers 
% -----------------------------------------------------------------------
function showOptions_Callback(hObject, eventdata, handles, rescale_clim) %#ok<INUSL>
persistent plotdata
if isempty(handles), return; end
if isempty(plotdata)||~isfield(plotdata,'figurehandle')||~isequal(plotdata.figurehandle,handles.figure1), % note: keeps only one file at a time in memory (to speed up processing while working on the same art project, while avoiding running out of memory when having too many art windows open)
    output_dir=getappdata(handles.savefile,'dir');
    art_mask_temporalfile=getappdata(handles.savefile,'art_mask_temporalfile');
    plotdata=load(fullfile(output_dir,art_mask_temporalfile));%,'maskscan','VY1','VY','notcoregistered','xyz','Data_Sum','Data_SumSquared');
    plotdata.Ma=spm_read_vols(plotdata.VY1);
    plotdata.Ma=plotdata.Ma/max(plotdata.Ma(:));
    plotdata.figurehandle=handles.figure1;
    temp=plotdata.VY1.mat(1:3,1:3);
    dir_order=zeros(1,3);
    [nill,dir_order(3)]=max(abs(temp(3,:)));dir_order(3)=dir_order(3)*sign(temp(3,dir_order(3)));
    [nill,dir_order(1)]=max(abs(temp(1,:)).*((1:3)~=dir_order(3)));dir_order(1)=dir_order(1)*sign(temp(1,dir_order(1)));
    dir_order(2)=setdiff(1:3,abs(dir_order([1,3]))); dir_order(2)=dir_order(2)*sign(temp(2,dir_order(2)));
    plotdata.dirorder=dir_order; % directions of storage dimensions most similar to spatial directions (x,y,z)
    first=1;
else
    first=0;
end
if nargin<4, rescale_clim=0; end

if isempty(handles)||get(handles.showOptions,'Value')==1,if ~isempty(handles),set([handles.axes_mask;get(handles.axes_mask,'children')],'visible','off'); set([handles.text_all_outliers,handles.all_outliers],'visible','on'); end; return;
elseif strcmp(get(handles.axes_mask,'visible'),'off'), set(get(handles.axes_mask,'children'),'visible','on'); set([handles.axes_mask,handles.text_all_outliers,handles.all_outliers],'visible','off'); end
option=get(handles.showOptions,'value');
if first, temp=get(handles.axes_mask,'children'); delete(temp(ishandle(temp))); end

% computes analysis mask
out_idx=round(str2num(get(handles.all_outliers, 'String')));
Mask=ones(plotdata.VY1.dim);
art_mask_temporalfile=getappdata(handles.savefile,'art_mask_temporalfile');
nbase=0;
for n0=[1:plotdata.maskscan_files 0]
    if n0, load(fullfile(output_dir,[art_mask_temporalfile(1:end-4),num2str(n0),'.mat']),'maskscan');
    else maskscan=plotdata.maskscan;
    end
    out_idx_infile=find(ismember(nbase+(1:numel(maskscan)),out_idx));
    for n1=setdiff(1:length(maskscan),out_idx_infile),Mask(maskscan{n1})=0;end
    nbase=nbase+numel(maskscan);
end
%for n1=setdiff(1:length(plotdata.maskscan),out_idx),Mask(plotdata.maskscan{n1})=0;end

if option>2
    % computes Var/SNR
    hw=[]; if numel(out_idx)>100, try, hw=waitbar(0,'updating Variance/SNR plots. Please wait'); end; end
    N=numel(plotdata.VY);
    Data_Sum=plotdata.Data_Sum;
    Data_SumSquared=plotdata.Data_SumSquared;
    for n1=1:numel(out_idx),
        i=out_idx(n1);
        if plotdata.notcoregistered, temp=reshape(spm_get_data(plotdata.VY(i),pinv(plotdata.VY(i).mat)*plotdata.xyz),plotdata.VY(1).dim);
        else temp=reshape(spm_get_data(plotdata.VY(i),plotdata.xyz_voxel),plotdata.VY(i).dim); end; %temp=spm_read_vols(plotdata.VY(i)); end
        Data_Sum=Data_Sum-temp;
        Data_SumSquared=Data_SumSquared-temp.^2;
        if numel(out_idx)>100&&~isempty(hw), waitbar(n1/numel(out_idx),hw); end
    end
    Data_Sum=Data_Sum/max(eps,N-numel(out_idx));
    Data_SumSquared=Data_SumSquared/max(eps,N-numel(out_idx));
    Data_Std=reshape(sqrt(max(0,Data_SumSquared-Data_Sum.^2)),plotdata.VY1.dim)*(N-numel(out_idx))/max(eps,N-numel(out_idx)-1);
    Data_SNR=Data_Sum./max(eps,Data_Std);
    if numel(out_idx)>100&&~isempty(hw), close(hw); end
end

switch(option)
    case 2, b=Mask.*(1+plotdata.Ma);cscale=nan;
    case 3, b=Mask.*Data_Std; cscale=max(b(:)); 
    case 4, b=Mask.*Data_SNR; cscale=max(b(:)); 
end

% generates display
b=permute(b,abs(plotdata.dirorder));
for n1=1:3,if plotdata.dirorder(n1)<0, b=b(end:-1:1,:,:); end; b=permute(b,[2,3,1]); end; % turn to ~spatial (xyz)
b=permute(b(:,end:-1:1,:),[2,1,3]);
sb=any(any(b,1),2);
b=b(:,:,sb);
slices=1:size(b,3);
set(handles.axes_mask,'units','points');
size1=get(handles.axes_mask,'position');
set(handles.axes_mask,'units','normalized');
size1=(size1(3)/size(b,2))/(size1(4)/size(b,1));
  nhoriz=1:length(slices);
  nverti=ceil(length(slices)./nhoriz);
  area=(min(size1./nhoriz,1./nverti).^2);  
[nill,nhoriz]=max(area);
temp=reshape(b(:,:,slices),[size(b,1),size(b,2)*length(slices)]);
temp2=[];
for n1=1:ceil(length(slices)/nhoriz),
    temp2=cat(1,temp2,[temp(:,size(b,2)*(n1-1)*nhoriz+1:size(b,2)*min(n1*nhoriz,length(slices))),zeros(size(b,1),max(0,size(b,2)*(n1*nhoriz-length(slices))))]);
end
if isnan(cscale),
    temp2=[temp2,zeros(size(temp2,1),20)];
else
    temp2=[temp2,zeros(size(temp2,1),10),linspace(cscale,0,size(temp2,1))'*ones(1,10)];
end
if first||~isfield(plotdata,'h')||~ishandle(plotdata.h),
    axes(handles.axes_mask);
    plotdata.h=imagesc(-temp2,'parent',handles.axes_mask);
    colormap(handles.axes_mask,gray);
    if isnan(cscale), set(handles.axes_mask,'clim',[-2,0],'ytick',[],'visible','off'); 
    else set(handles.axes_mask,'ytick',1,'yticklabel',{num2str(cscale,'%0.1f')},'visible','on'); end
    set(handles.axes_mask,'xtick',[],'xcolor','w','yaxislocation','right','box','off');
    axis(handles.axes_mask,'equal','tight');
    if option==3, plotdata.h2=text(size(temp2,2)+8,size(temp2,1)/2,'standard deviation','rotation',90,'fontsize',8,'horizontalalignment','center');
    elseif option==4, plotdata.h2=text(size(temp2,2)+8,size(temp2,1)/2,'SNR','rotation',90,'fontsize',8,'horizontalalignment','center');
    else plotdata.h2=text(size(temp2,2)+8,size(temp2,1)/2,'','rotation',90,'fontsize',8,'horizontalalignment','center');
    end
else
    set(plotdata.h,'cdata',-temp2);
    if isnan(cscale), set(handles.axes_mask,'clim',[-2,0],'ytick',[],'visible','off'); end
    if ~isnan(cscale)&&rescale_clim, set(handles.axes_mask,'clim',[-cscale,0],'ytick',1,'yticklabel',{num2str(cscale,'%0.1f')},'visible','on'); end
    axis(handles.axes_mask,'equal','tight');
    if option==3, set(plotdata.h2,'string','standard deviation');
    elseif option==4, set(plotdata.h2,'string','SNR');
    else set(plotdata.h2,'string','');
    end
    set(plotdata.h2,'position',[size(temp2,2)+8,size(temp2,1)/2,0]);
end
end


% -----------------------------------------------------------------------
% SAVEFILE_CALLBACK
% Saves art output files (motion statistics, graphs, SPM regressors, or
% Analysis mask)
% -----------------------------------------------------------------------
function savefile_Callback(hObject, eventdata, handles) %#ok<INUSL>

[S,v] = listdlg('PromptString','What would you like to save?',...
    'SelectionMode','multiple',...
    'ListSize',[160,120], ...
    'ListString',{'Outliers','Motion Statistics','Graphs','SPM regressors','Analysis mask','Voxel-wise Variability','Voxel-wise SNR','Mean functional'});

if v==0
    return;
end
%get path
tmpdir = pwd;
pathname = '.'; %getappdata(hObject, 'path');
cd(pathname);
for s=S(:)',
    switch s
        case 1 %save outliers
            out_idx = round(str2num(get(handles.all_outliers, 'String'))); %#ok<NASGU>
            %ask user to choose filename
            filter = {'*.mat';'*.txt'};
            %ext = {'.mat';'.txt'};
            [filename, pathname, filteridx] = uiputfile( filter,'Save outliers as:');
            
            %save according to file format
            switch filteridx
                %binary MAT file
                case 1
                    filename = strcat(char(pathname), char(filename));
                    save(filename ,'out_idx', '-mat');
                    %txt file
                case 2
                    filename = strcat(char(pathname), char(filename));
                    save(filename,'out_idx','-ascii');
            end
            
        case 2 %save motion statistics
            mv_stats = getappdata(handles.mvthresh,'mv_stats'); %#ok<NASGU>
            analyses=getappdata(handles.savefile,'analyses'); %#ok<NASGU>
            
            %save statistics to .mat file
            %mv_stats has 7 columns corresponding to x y z pitch roll yaw norm
            %and 3 rows corresponding to mean, stdv and max of the absolute values of
            %the movement parameters
            
            %ask user to choose filename
            filter = {'*.mat'; '*.txt'};
            [filename, pathname, filteridx] = uiputfile( filter,'Save motion statistics as');
            
            %save according to file format
            switch filteridx
                case 1 %binary MAT file
                    filename = strcat(char(pathname), char(filename));
                    save(filename ,'mv_stats','analyses','-mat');
                case 2 %txt file
                    filename = strcat(char(pathname), char(filename));
                    save(filename,'mv_stats','analyses','-ascii');
            end
            
        case 3 %save graphs
            %ask user to choose filename
            filter = {'*.jpg';'*.eps';'*.fig'};
            [filename, pathname] = uiputfile( filter,'Save figure as:');
            
            filename = strcat(char(pathname), char(filename));
            saveas(gcf,filename);
            
        case 4, % save SPM regressors
            savefile_Callback_SaveRegressor(handles);
            
        case {5,6,7,8}, % save mask
            savefile_Callback_SaveMask(handles,s-4);
    end
end

cd(tmpdir);
end

% -----------------------------------------------------------------------
% SAVEFILE_CALLBACK_SAVEREGRESSORS
% saves SPM regressor files 
% One regressor file per session, named art_regression_outliers_*.mat, and
% stored in the original location of each functional series (or if not
% possible to current directory)
% Regressor matrices contains 1's at the location of each outlier. This
% implements outlier removal in SPM when these regressor files are used as
% covariates 
% -----------------------------------------------------------------------
function savefile_Callback_SaveRegressor(handles)
g = getappdata(handles.zthresh,'g');
drop_flag=getappdata(handles.mvthresh,'drop_flag');
M = getappdata(handles.savefile,'mv_data_raw');
mv_data = getappdata(handles.mvthresh,'mv_data');
num_sess = length(g);
out_idx = round(str2num(get(handles.all_outliers, 'String')));
datafiles=getappdata(handles.savefile,'datafiles');

cur_idx=0;
for j=1:num_sess,
    idx=find(out_idx>cur_idx&out_idx<=cur_idx+size(g{j},1));
    lidx=out_idx(idx)-cur_idx;
    if drop_flag, lidx=union(lidx, 1:min(size(g{j},1),drop_flag)); end
    if size(g{j},1)~=size(M{j},1), error('mismatch number of dimensions (%d scans from BOLD signal; %d timepoints from movement parameters)',size(g{j},1),size(M{j},1)); end
    R1=zeros(size(g{j},1),length(lidx));
    R1(lidx,:)=eye(length(lidx));
    RT=cat(2,abs(g{j}(:,2)),mv_data(cur_idx+(1:size(g{j},1)),32)); % threshold timeseries (global&motion)
    RB=cat(2,M{j},mv_data(cur_idx+(1:size(g{j},1)),32)); % extended motion
    R2=cat(2,R1,RB); % outliers and extended motion
    [datafiles_path,datafiles_name] = fileparts(datafiles{j});
    art_disp(['Saving SPM regressor file ',fullfile(datafiles_path,['art_regression_outliers_',datafiles_name,'.mat']),' and ',fullfile(datafiles_path,['art_regression_outliers_and_movement_',datafiles_name,'.mat'])]);
    try
        R=R1;save(fullfile(datafiles_path,['art_regression_outliers_',datafiles_name,'.mat']),'R','-mat'); %#ok<NASGU>
        R=R2;save(fullfile(datafiles_path,['art_regression_outliers_and_movement_',datafiles_name,'.mat']),'R','-mat'); %#ok<NASGU>
        R=RT;save(fullfile(datafiles_path,['art_regression_timeseries_',datafiles_name,'.mat']),'R','-mat'); %#ok<NASGU>
        R=[min(RB,[],1);mean(RB,1);max(RB,[],1)];save(fullfile(datafiles_path,['art_movement_stats_',datafiles_name,'.mat']),'R','-mat'); %#ok<NASGU>
    catch
        try
            tdatafiles_path=datafiles_path; tdatafiles_path(tdatafiles_path==filesep)='_';
            R=R1;save(fullfile('./',[tdatafiles_path,'_art_regression_outliers_',datafiles_name,'.mat']),'R','-mat'); %#ok<NASGU>
            R=R2;save(fullfile('./',[tdatafiles_path,'_art_regression_outliers_and_movement_',datafiles_name,'.mat']),'R','-mat'); %#ok<NASGU>
            R=RT;save(fullfile('./',[tdatafiles_path,'_art_regression_timeseries_',datafiles_name,'.mat']),'R','-mat'); %#ok<NASGU>
            R=[min(RB,[],1);mean(RB,1);max(RB,[],1)];save(fullfile('./',[tdatafiles_path,'_art_movement_stats_',datafiles_name,'.mat']),'R','-mat'); %#ok<NASGU>
        catch
            [filename, pathname] =uiputfile({'*.mat'},['Save session ',num2str(j),' outliers regressor:'],fullfile(datafiles_path,['art_regression_outliers_',datafiles_name,'.mat']));
            filename = fullfile(char(pathname), char(filename));
            R=R1;save(filename ,'R','-mat'); %#ok<NASGU>
            [filename, pathname] =uiputfile({'*.mat'},['Save session ',num2str(j),' outliers&motion regressor:'],fullfile(datafiles_path,['art_regression_outliers_and_movement_',datafiles_name,'.mat']));
            filename = fullfile(char(pathname), char(filename));
            R=R2;save(filename ,'R','-mat'); %#ok<NASGU>
            [filename, pathname] =uiputfile({'*.mat'},['Save session ',num2str(j),' motion regressor:'],fullfile(datafiles_path,['art_regression_timeseries_',datafiles_name,'.mat']));
            filename = fullfile(char(pathname), char(filename));
            R=RT;save(filename ,'R','-mat'); %#ok<NASGU>
            [filename, pathname] =uiputfile({'*.mat'},['Save session ',num2str(j),' stats:'],fullfile(datafiles_path,['art_movement_stats_',datafiles_name,'.mat']));
            filename = fullfile(char(pathname), char(filename));
            R=[min(RB,[],1);mean(RB,1);max(RB,[],1)];save(filename ,'R','-mat'); %#ok<NASGU>
        end
    end
    
    cur_idx=cur_idx+size(g{j},1);
end
end

% -----------------------------------------------------------------------
% SAVEFILE_CALLBACK_SAVEMASK
% saves Analysis mask after disregarding outlier scans.
% Filename 'art_mask.img' saved to current directory or location specified
% by 'output_dir' field in .cfg file
% -----------------------------------------------------------------------
function savefile_Callback_SaveMask(handles,option)
if nargin<2, option=[1 4]; end
output_dir=getappdata(handles.savefile,'dir');

art_mask_temporalfile=getappdata(handles.savefile,'art_mask_temporalfile');
plotdata=load(fullfile(output_dir,art_mask_temporalfile));%,'maskscan','VY1','VY','notcoregistered','xyz','Data_Sum','Data_SumSquared');
out_idx=round(str2num(get(handles.all_outliers, 'String')));
datafiles=getappdata(handles.savefile,'datafiles');
[datafiles_path,datafiles_name,datafiles_ext] = spm_fileparts(datafiles{1});
% computes Analysis Mask
Mask=ones(plotdata.VY1.dim);
nbase=0;
for n0=[1:plotdata.maskscan_files 0]
    if n0, load(fullfile(output_dir,[art_mask_temporalfile(1:end-4),num2str(n0),'.mat']),'maskscan');
    else maskscan=plotdata.maskscan;
    end
    out_idx_infile=find(ismember(nbase+(1:numel(maskscan)),out_idx));
    for n1=setdiff(1:length(maskscan),out_idx_infile),Mask(maskscan{n1})=0;end
    nbase=nbase+numel(maskscan);
end
%for n1=setdiff(1:length(plotdata.maskscan),out_idx),Mask(plotdata.maskscan{n1})=0;end
if any(option>1)
    % computes Var/SNR
    hw=[]; if numel(out_idx)>100, try, hw=waitbar(0,'updating Variance/SNR plots. Please wait'); end; end
    N=numel(plotdata.VY);
    Data_Sum=plotdata.Data_Sum;
    Data_SumSquared=plotdata.Data_SumSquared;
    for n1=1:numel(out_idx),
        i=out_idx(n1);
        if plotdata.notcoregistered, temp=reshape(spm_get_data(plotdata.VY(i),pinv(plotdata.VY(i).mat)*plotdata.xyz),plotdata.VY(1).dim);
        else temp=reshape(spm_get_data(plotdata.VY(i),plotdata.xyz_voxel),plotdata.VY(i).dim); end %temp=spm_read_vols(plotdata.VY(i)); end
        Data_Sum=Data_Sum-temp;
        Data_SumSquared=Data_SumSquared-temp.^2;
        if numel(out_idx)>100&&~isempty(hw), waitbar(n1/numel(out_idx),hw); end
    end
    Data_Sum=Data_Sum/max(eps,N-numel(out_idx));
    Data_SumSquared=Data_SumSquared/max(eps,N-numel(out_idx));
    Data_Std=reshape(sqrt(max(0,Data_SumSquared-Data_Sum.^2)),plotdata.VY1.dim)*(N-numel(out_idx))/max(eps,N-numel(out_idx)-1);
    Data_SNR=Data_Sum./max(eps,Data_Std);
    if numel(out_idx)>100&&~isempty(hw), close(hw); end
end
for noption=1:numel(option)
    ffilename='';
    switch option(noption)
        case 1, b=Mask; filename=['art_mask_',datafiles_name,datafiles_ext]; filetype='uint8'; filedescription='analysis mask';
        case 2, b=Data_Std; filename=['art_ResStd_',datafiles_name,datafiles_ext]; filetype='float32'; filedescription='Voxel-wise variability';
        case 3, b=Mask.*Data_SNR; filename=['art_SNR_',datafiles_name,datafiles_ext]; filetype='float32'; filedescription='Voxel-wise SNR';
        case 4, b=Data_Sum; filename=['art_mean_',datafiles_name,datafiles_ext]; filetype='float32'; filedescription='Voxel-wise mean';
    end
    V=struct('mat',plotdata.VY1.mat,'dim',plotdata.VY1.dim,'pinfo',[1;0;0],'dt',[spm_type(filetype) spm_platform('bigend')],'descrip',filedescription);
    if ~isempty(ffilename), V.fname=ffilename;
    else V.fname=fullfile(output_dir,filename);
    end
    spm_write_vol(V,b);
    art_disp(['New ',filedescription,' file saved to ',V.fname]);
end
end

% -----------------------------------------------------------------------
% READ_ART_SESS_FILE
% reads a .cfg file (art session file) defining the art processing options
% ohinds 2008-04-23
% -----------------------------------------------------------------------
function [num_sess,global_type_flag,drop_flag,gui_display,motionFileType,motion_threshold,global_threshold,use_diff_motion,use_diff_global,use_norms,SPMfile,mask_file,output_dir,do_close,do_print,P,M] = read_art_sess_file(sess_file)

loadfromfile=~isstruct(sess_file);
if loadfromfile&&~exist(sess_file,'file')
    error(['session file ' sess_file ' cant be opened']);
end

num_sess = 1;
global_type_flag = 1;
drop_flag = 0;
gui_display = 1;
motionFileType = 0;
image_dir = ''; 
motion_dir = '';
auto_motion_fname = 0;
motion_threshold=[];
global_threshold=[];
use_diff_motion=1;
use_diff_global=1;
use_norms=1;
SPMfile=[];
mask_file=[];
output_dir='';
do_close=false;
do_print=false;

if loadfromfile
    fp = fopen(sess_file);
    % read each param
    s = fscanf(fp,'%s',1);
    while ~strcmp(s,'end')
        if ~isempty(s) && s(1) == '#'
            % skips until end of line
            fgetl(fp);
        elseif strcmp(s,'sessions:')
            % num_sessions
            num_sess = fscanf(fp,'%d',1);
        elseif strcmp(s,'global_mean:')
            % global_type_flag
            global_type_flag = fscanf(fp,'%d',1);
        elseif strcmp(s,'drop_flag:')
            % drop_flag
            drop_flag = fscanf(fp,'%d',1);
        elseif strcmp(s,'gui_display:')
            gui_display = fscanf(fp,'%d',1);
        elseif strcmp(s,'motion_file_type:')
            motionFileType = fscanf(fp,'%d',1);
        elseif strcmp(s,'image_dir:')
            image_dir = fscanf(fp,'%s',1);
        elseif strcmp(s,'motion_dir:')
            motion_dir = fscanf(fp,'%s',1);
        elseif strcmp(s,'motion_fname_from_image_fname:')
            auto_motion_fname = str2num(fscanf(fp,'%s',1));
        elseif strcmp(s,'motion_threshold:')
            motion_threshold = fscanf(fp,'%f',2);
        elseif strcmp(s,'global_threshold:')
            global_threshold = fscanf(fp,'%f',1);
        elseif strcmp(s,'spm_file:')
            SPMfile = fscanf(fp,'%s',1);
        elseif strcmp(s,'use_diff_motion:')
            use_diff_motion = fscanf(fp,'%d',1);
        elseif strcmp(s,'use_diff_global:')
            use_diff_global = fscanf(fp,'%d',1);
        elseif strcmp(s,'use_norms:')
            use_norms = fscanf(fp,'%d',1);
        elseif strcmp(s,'mask_file:')
            mask_file = fscanf(fp,'%s',1);
        elseif strcmp(s,'output_dir:')
            output_dir = art_fullfile(fscanf(fp,'%s',1));
        elseif strcmp(s,'close:')
            do_close = fscanf(fp,'%d',1);
        elseif strcmp(s,'print:')
            do_print = fscanf(fp,'%d',1);
        end
        s = fscanf(fp,'%s',1);
    end
    
    M = {};
    P = {};
    % read the filenames
    s = fscanf(fp,'%s',1);
    while(~strcmp(s,'end'))
        if strcmp(s,'session')
            sess = fscanf(fp,'%d',1);
            type = fscanf(fp,'%s',1);
            
            
            % set up P
            if size(P,2) < sess
                P{sess} = {}; %#ok<AGROW>
            end
        elseif numel(s)>0&&s(1) == '#'
            % skips until end of line
            fgetl(fp);
        elseif strcmp(type,'image')
            if any(s=='?'),
                idx=find(s=='?');
                ns=length(idx);
                for sn=0:10^ns-1,
                    st=s;st(idx)=num2str(sn,['%0',num2str(ns),'d']);
                    if ~isempty(dir(fullfile(image_dir,st))),
                        P{sess}{end+1} = fullfile(image_dir,st); %#ok<AGROW>
                    end
                end
                s(idx)=num2str(1,['%0',num2str(ns),'d']);
            else
                P{sess}{end+1} = fullfile(image_dir,s); %#ok<AGROW>
            end
            
            if auto_motion_fname && length(P{sess})<=1
                tmotion_dir=image_dir;
                if motionFileType == 2
                    M{sess} = read_siemens_motion_parm_file(strprepend('',fullfile(tmotion_dir,s),'.txt')); %#ok<AGROW>
                elseif motionFileType == 0
                    M{sess}=[]; %#ok<AGROW>
                    for n1=0:5,
                        if ~isempty(dir(strprepend('rp_',strprepend(-n1,fullfile(tmotion_dir,s)),'.txt'))),
                            M{sess} = load(strprepend('rp_',strprepend(-n1,fullfile(tmotion_dir,s)),'.txt')); %#ok<AGROW>
                            break;
                        end
                    end
                    if isempty(M{sess}),error(['No motion file found: ',strprepend('rp_',fullfile(tmotion_dir,s),'.txt'),' or similar']); end
                else
                    M{sess} = load(strprepend('',fullfile(tmotion_dir,s),'.par')); %#ok<AGROW>
                end
            end
            
        elseif strcmp(type,'movement') || strcmp(type,'motion')
            if motionFileType == 2
                M{sess} = read_siemens_motion_parm_file(fullfile(motion_dir,s)); %#ok<AGROW>
            else
                M{sess} = load(fullfile(motion_dir,s)); %#ok<AGROW>
            end
        end
        s = fscanf(fp,'%s',1);
    end
    
    for i=1:numel(P)
        P{i} = char(P{i}); %#ok<AGROW> %make_spm_file_matrix(P{i});
    end
    
    fclose(fp);
else
    % num_sess,global_type_flag,drop_flag,motionFileType,motion_threshold,
    % global_threshold,use_diff_motion,use_diff_global,use_norms,SPMfile,
    % mask_file,output_dir,P,M
    % P{nses}{nfile} (files); M{nses} (motion files)
    
    if isfield(sess_file,'num_sess'), num_sess=sess_file.num_sess; end
    if isfield(sess_file,'global_type_flag'), global_type_flag=sess_file.global_type_flag; end
    if isfield(sess_file,'drop_flag'), drop_flag=sess_file.drop_flag; end
    if isfield(sess_file,'gui_display'), gui_display=sess_file.gui_display; end
    if isfield(sess_file,'motionFileType'), motionFileType=sess_file.motionFileType; end
    if isfield(sess_file,'motion_threshold'), motion_threshold=sess_file.motion_threshold; end
    if isfield(sess_file,'global_threshold'), global_threshold=sess_file.global_threshold; end
    if isfield(sess_file,'use_diff_motion'), use_diff_motion=sess_file.use_diff_motion; end
    if isfield(sess_file,'use_diff_global'), use_diff_global=sess_file.use_diff_global; end
    if isfield(sess_file,'use_norms'), use_norms=sess_file.use_norms; end
    if isfield(sess_file,'SPMfile'), SPMfile=sess_file.SPMfile; end
    if isfield(sess_file,'mask_file'), mask_file=sess_file.mask_file; end
    if isfield(sess_file,'output_dir'), output_dir=sess_file.output_dir; end
    if isfield(sess_file,'close'), do_close=sess_file.close; end
    if isfield(sess_file,'print'), do_print=sess_file.print; end
    if isfield(sess_file,'P'), P=sess_file.P; end
    if isfield(sess_file,'M'), M=sess_file.M; end

    if isfield(sess_file,'sessions'), num_sess=sess_file.sessions; end
    if isfield(sess_file,'global_mean'), global_type_flag=sess_file.global_mean; end
    if isfield(sess_file,'motion_file_type'), motionFileType=sess_file.motion_file_type; end
    if isfield(sess_file,'motion_fname_from_image_fname'), motion_threshold=sess_file.motion_fname_from_image_fname; end
    if isfield(sess_file,'spm_file'), SPMfile=sess_file.spm_file; end

    % note: fix to avoid eval (8/2014)
%     convertnames={'sessions','num_sess'; 'global_mean','global_type_flag'; 'motion_file_type','motionFileType'; 'motion_fname_from_image_fname','motion_threshold'; 'spm_file','SPMfile'}; 
%     fields=fieldnames(sess_file);
%     for n=1:numel(fields)
%         [ok,idx]=ismember(fields{n},convertnames(:,1));
%         idx=find(idx,1);
%         if ~isempty(idx), eval([convertnames(idx,2),'=sess_file.',fields{n}]);
%         else eval([fields{n},'=sess_file.',fields{n}]);
%         end
%     end
    switch(motionFileType)
        case 2, for n=1:numel(M),M{n}=read_siemens_motion_parm_file(M{n}); end
        otherwise, for n=1:numel(M),M{n}=load(M{n}); end
    end
    num_sess=numel(M);
end
end



%% ---------------   Other GUI functions -------------------------------
% The following functions handle simple gui callback events
% ----------------------------------------------------------------------

% --- Executes on button press in z_up.
function z_up_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
UpdateGlobal(hObject, eventdata, handles, 1.05);
UnionOutliers(hObject, eventdata, handles)
UpdateSummaryplot(hObject, eventdata, handles);
end

% --- Executes on button press in z_down.
function z_down_Callback(hObject, eventdata, handles)
UpdateGlobal(hObject, eventdata, handles,1/1.05);
UnionOutliers(hObject, eventdata, handles);
UpdateSummaryplot(hObject, eventdata, handles);
end

% --- Executes on zthreshold.
function zthresh_Callback(hObject, eventdata, handles)
UpdateGlobal(hObject, eventdata, handles,1.0);
UnionOutliers(hObject, eventdata, handles);
UpdateSummaryplot(hObject, eventdata, handles);
end

% --- Executes on button press in mv_up.
function mv_up_Callback(hObject, eventdata, handles)
UpdateMovement(hObject, eventdata, handles,1.05);
UnionOutliers(hObject, eventdata, handles);
UpdateSummaryplot(hObject, eventdata, handles);
end

% --- Executes on button press in mv_down.
function mv_down_Callback(hObject, eventdata, handles)
UpdateMovement(hObject, eventdata, handles,1/1.05);
UnionOutliers(hObject, eventdata, handles);
UpdateSummaryplot(hObject, eventdata, handles);
end

% --- Executes on update of mvthresh.
function mvthresh_Callback(hObject, eventdata, handles)
UpdateMovement(hObject, eventdata, handles,1.0);
UnionOutliers(hObject, eventdata, handles);
UpdateSummaryplot(hObject, eventdata, handles);
end

% --- Executes on button press in rt_up.
function rt_up_Callback(hObject, eventdata, handles)
UpdateRotation(hObject, eventdata, handles,1.05)
UnionOutliers(hObject, eventdata, handles);
UpdateSummaryplot(hObject, eventdata, handles);
end

% --- Executes on button press in rt_down.
function rt_down_Callback(hObject, eventdata, handles)
UpdateRotation(hObject, eventdata, handles,1/1.05)
UnionOutliers(hObject, eventdata, handles);
UpdateSummaryplot(hObject, eventdata, handles);
end

% --- Executes on update of mvthresh.
function rtthresh_Callback(hObject, eventdata, handles)
UpdateRotation(hObject, eventdata, handles,1.0)
UnionOutliers(hObject, eventdata, handles);
UpdateSummaryplot(hObject, eventdata, handles);
end

% --- Executes when comp motion checkbox is changed
function norms_Callback(hObject, eventdata, handles)
UpdateMovement(hObject, eventdata, handles,1.0);
UpdateRotation(hObject, eventdata, handles,1.0);
UnionOutliers(hObject, eventdata, handles);
UpdateSummaryplot(hObject, eventdata, handles);
end

% --- Executes when editing the all_outliers list
function all_outliers_Callback(hObject, eventdata, handles)
UpdateSummaryplot(hObject, eventdata, handles);
end

% --- Executes during object creation, after setting all properties.
function all_outliers_CreateFcn(hObject, eventdata, handles) %#ok<INUSD>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Outputs from this function are returned to the command line.
function varargout = art_OutputFcn(hObject, eventdata, handles) %#ok<INUSL>
varargout{1} = handles.output;
end

% --- Executes during object creation, after setting all properties.
function zthresh_CreateFcn(hObject, eventdata, handles) %#ok<INUSD>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes during object creation, after setting all properties.
function mvthresh_CreateFcn(hObject, eventdata, handles) %#ok<INUSD>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes during object creation, after setting all properties.
function rtthresh_CreateFcn(hObject, eventdata, handles) %#ok<INUSD>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes when show design checkbox is changed
function showDesign_Callback(hObject, eventdata, handles)
UpdateSummaryplot(hObject, eventdata, handles)
end


%% ---------------   Utility functions -------------------------------
% General utility functions (not specific to ART)
% --------------------------------------------------------------------

% -----------------------------------------------------------------------
% ART_MASKGLOBAL_SCAN
% Computes the analysis mask for a single volume
% -----------------------------------------------------------------------
function [idxremove_analysis,idxremove_globalsignal]=art_maskglobal_scan(data,VYi,VY1,VY1inv)
data1=data>mean(data(~isnan(data)))/8; % global-signal mask
idxremove_globalsignal=find(data1~=true);
data=(data>0.80*mean(data(data1>0))); % analysis mask
idxremove_analysis=find(data~=true);
if any(any(VYi.mat~=VY1.mat)), % resample global-signal and analysis mask voxels to VY1 space
    [tempx,tempy,tempz]=ind2sub(VYi.dim(1:3),idxremove_analysis);
    xyz=round(VY1inv*VYi.mat*[tempx(:),tempy(:),tempz(:),ones(numel(tempx),1)]');
    idxremove_analysis=sub2ind(VY1.dim(1:3),max(1,min(VY1.dim(1),xyz(1,:))),max(1,min(VY1.dim(2),xyz(2,:))),max(1,min(VY1.dim(3),xyz(3,:))));
    [tempx,tempy,tempz]=ind2sub(VYi.dim(1:3),idxremove_globalsignal);
    xyz=round(VY1inv*VYi.mat*[tempx(:),tempy(:),tempz(:),ones(numel(tempx),1)]');
    idxremove_globalsignal=sub2ind(VY1.dim(1:3),max(1,min(VY1.dim(1),xyz(1,:))),max(1,min(VY1.dim(2),xyz(2,:))),max(1,min(VY1.dim(3),xyz(3,:))));
end
idxremove_analysis=uint32(idxremove_analysis(:));
idxremove_globalsignal=uint32(idxremove_globalsignal(:));
end

% -----------------------------------------------------------------------
% GET_DESIGN
% Imports SPM design matrix information from SPM.mat file
% -----------------------------------------------------------------------
function [SPM,design,names] = get_design(handles)
SPM = getappdata(handles.showDesign,'SPM');
SPMbak=SPM;
if (isempty(SPM))
    if isdeployed, spm_ver='SPM12';
    else spm_ver = spm('ver');
    end
    switch(spm_ver),
        case 'SPM99', spm_ver=1;
        case 'SPM2', spm_ver=2;
        case 'SPM5', spm_ver=5;
        case {'SPM8b','SPM8'}, spm_ver=8;
        case {'SPM12b','SPM12'}, spm_ver=12;
        otherwise, art_disp(['Warning! unrecognized SPM version ',spm_ver]); spm_ver=8;
    end
    switch spm_ver
        case {1,2}
            tmpfile = spm_get(1,'.mat','Select design matrix:');
        case {5,8,12}
            tmpfile = spm_select(1,'^.*\.mat$','Select design matrix:');
    end
    if isempty(tmpfile), error('No design matrix selected'); end
    load(tmpfile);
    setappdata(handles.showDesign,'SPMfile',tmpfile);
    setappdata(handles.showDesign,'SPM',SPM);
end
sessions = getappdata(handles.showDesign,'sessions');
if isempty(sessions)||any(sessions<0),
    if length(sessions)==length(SPM.Sess) && all(sessions==-(1:length(sessions))),
        sessions = abs(sessions);
    elseif numel(SPM.Sess)<numel(sessions)
        setappdata(handles.showDesign,'SPM',SPMbak);
        error('Incorrect number of sessions in SPM.mat file');
    else
        tmpsess = inputdlg('What session(s) to use? (e.g. 1 or [1,2])','',1,{['[',num2str(abs(sessions)),']']});
        sessions = str2num(char(tmpsess));
    end
    setappdata(handles.showDesign,'sessions',sessions);
end
rows = [];
cols = [];
names={};
for s = sessions
    rows = [rows SPM.Sess(s).row]; %#ok<AGROW>
    cols = [cols SPM.Sess(s).col(1:length(SPM.Sess(s).U))]; %#ok<AGROW> % extracts only effects of interest (no covariates) from design matrix
    names=cat(1,names,SPM.Sess(s).U(:).name);
end

design = SPM.xX.X(rows,cols);
end

% -----------------------------------------------------------------------
% READ_SIEMENS_MOTION_PARM_FILE
% reads a siemens motion detection parameter file
% Oliver Hinds <ohinds@mit.edu>
% 2007-07-23
% -----------------------------------------------------------------------
function mp = read_siemens_motion_parm_file(fname)

mp = [];
% open the file
fp = fopen(fname);
if fp == -1
    error('coulnd''t open motion parm file.');
end

% read each parameter
i = 1;
while(~feof(fp))
    % read the motion header
    fscanf(fp,'%s',6);    
    if feof(fp)
        break;
    end    
    if i == 1
        fscanf(fp,'%s',5);
    end    
    fscanf(fp,'%s',7);
    for j=1:6
        fscanf(fp,'%s',4);
        mp(i,j) = fscanf(fp,'%f',1); %#ok<AGROW>
    end    
    fscanf(fp,'%s',1);
    i=i+1;
end

% siemens keeps their params in a different order than spm does, and their rotations in degrees. fix it
m = mp;
mp(:,4:6) = m(:,4:6)*pi/180;
mp(:,1)   = m(:,2);
mp(:,2)   = m(:,1);
mp(:,3)   =-m(:,3);
end

% -----------------------------------------------------------------------
% MAKE_SPM_FILE_MATRIX
% Take a cell array and pad appropritately to make a matrix
% -----------------------------------------------------------------------
function m = make_spm_file_matrix(p)

mx = -1;
for i=1:numel(p)
    if size(p{i},2) > mx
        mx = size(p{i},2);
    end
end

m = char(32*ones(numel(p),mx));
for i=1:numel(p)
    m(i,mx-length(p{i})+1:end) = p{i};
end
end

% -----------------------------------------------------------------------
% ZSCORE
% Computes standardized z-scores
% -----------------------------------------------------------------------
function z = zscore(x)
%ZSCORE Standardized z score.
% z=zscore(x);
%
stdx=std(x);stdx(stdx==0)=1;
z=(x-mean(x))./stdx;
end

% -----------------------------------------------------------------------
% RANGE
% Computes sample range
% -----------------------------------------------------------------------
function d=range(x)
%RANGE  Sample range.
%d=range(x);
%
d=max(x)-min(x);
end

% -----------------------------------------------------------------------
% STRPREPEND
% pre-pend filename with string
% -----------------------------------------------------------------------
function fileout=strprepend(str1,file,str2)

[fpath,ffile,fext]=fileparts(file);
if nargin<3, str2=fext; end
if ~ischar(str1),ffile=ffile(1+abs(str1):end);str1=''; end
if ~ischar(str2),ffile=ffile(1:end-abs(str2));str2=''; end
fileout=fullfile(fpath,[str1,ffile,str2]);
end

% -----------------------------------------------------------------------
% CUMDISP
% Persistent display
% -----------------------------------------------------------------------
function cumdisp(txt)
% CUMDISP persistent disp
% cumdisp; initializes persistent display
% cumdisp(text); displays persistent text
%
persistent oldtxt;
if nargin<1,
    oldtxt='';
    fprintf(1,'\n');
else
    fprintf(1,[repmat('\b',[1,length(oldtxt)]),txt]);
    oldtxt=sprintf(txt);
end
end

function art_disp(varargin)
persistent isconn
if isempty(isconn), isconn=~isempty(which('conn_disp')); end
if isconn, conn_disp(varargin{:});
elseif nargin>0&&ischar(varargin{1})&&strcmp(varargin{1},'fprintf'), fprintf(varargin{2:end});
else disp(varargin{:});
end
end
% -----------------------------------------------------------------------
% PRCTILE
% computes smaple percentile 
% -----------------------------------------------------------------------
function z=prctile(x,p)
nx=length(x);
z=zeros(size(p));
sx=sort(x);
q = [0,100*(0.5:(nx-0.5))./nx,100]';
xx = [sx(1);sx(:);sx(end)];
z(:) = interp1q(q,xx,p(:));
end

% -----------------------------------------------------------------------
% HANNING
% Outputs Hann window
% -----------------------------------------------------------------------
function w=hanning(n)

if ~rem(n,2),%even
    w = .5*(1 - cos(2*pi*(1:n/2)'/(n+1)));
    w=[w;flipud(w)];
else %odd
    w = .5*(1 - cos(2*pi*(1:(n+1)/2)'/(n+1)));
    w = [w; flipud(w(1:end-1))];
end
end

% -----------------------------------------------------------------------
% ART_FULLFILE
% Builds absolute full filename from parts
% -----------------------------------------------------------------------
function filename=art_fullfile(varargin)

if ~nargin, filename=pwd; return; 
elseif nargin==1, filename=varargin{1}; 
else filename=fullfile(varargin{:});
end
if isempty(filename), filename=pwd; return; end
[filename_path,filename_name,filename_ext]=fileparts(filename);
if isempty(filename_path),
    filename_path=pwd;
else
    cwd=pwd;
    cd(filename_path);
    filename_path=pwd;
    cd(cwd);
end
filename=fullfile(filename_path,[filename_name,filename_ext]);
end

