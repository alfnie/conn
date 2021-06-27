function evlab17_roiresultsplots(filename,varargin)
% EVLAB17_ROIRESULTSPLOTS displays second-level ROI-based analysis results
%   evlab17_roiresultsplots('/myfolder/REX.mat');
%      displays second-level ROI-based results stored in /myfolder
%   evlab17_roiresultsplots('atlas.cfg');
%      displays second-level ROI-based results defined in atlas.cfg file (cfg file needs to contain #folder field specifying target folder and #roi_id specifying target filename)
%

evlab17_module init silent;

% interprets info
results_path='';
results_name='';
if isdir(filename), results_path=filename; filename=[]; end
if ~isempty(filename)
    [t1,t2,ext]=spm_fileparts(filename);
    if strcmp(ext,'.cfg')
        options=conn_loadcfgfile(filename);
        assert(isfield(options,'folder'),'missing #folder information in %s',filename);
        results_path=options.folder;
        if isfield(options,'roi_id'), results_name=sprintf('REX_%s.mat',options.roi_id); end
    else 
        results_path=t1;
        results_name=[t2,ext];
    end
end
if isempty(results_name)
    [tfilename,tpathname]=uigetfile({'REX*.mat','REX file (REX*.mat)'; '*',  'All Files (*)'},'Select ROI results file');
    if ~ischar(tfilename)||isempty(tfilename), return; end
    results_path=tpathname;
    results_name=tfilename;
end
fprintf('Loading file %s\n',fullfile(results_path,results_name));
conn_rex('results',fullfile(results_path,results_name),'output_type','saverex','gui',0); 


end



