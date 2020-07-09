function conn_assignfslutlabels(iroi)
% assigns conn/utils/surf/FreeSurferColorLUT.txt labels file to ROI (e.g. for aparc+aseg.mgz files)
% conn_assignfslutlabels(iroi)
%  iroi : ROI name or ROI index in CONN project
%

OVERWRITE=false;

global CONN_x;
if ischar(iroi), iroi=find(strcmp(iroi,CONN_x.Setup.rois.names(1:end-1)),1); end
if isempty(iroi)||iroi<1||iroi>=numel(CONN_x.Setup.rois.names), error('Out of range ROI index %d',iroi); end
files=conn_prepend('',CONN_x.Setup.rois.files{iroi},'.txt');
file0=fullfile(fileparts(which(mfilename)),'utils','surf','FreeSurferColorLUT.txt');
for n1=1:numel(files),
    if ~iscell(files{n1}), files{n1}={files{n1}}; end
    for n2=1:numel(files{n1})
        if ~iscell(files{n1}{n2}), files{n1}{n2}={files{n1}{n2}}; end
        if OVERWRITE||~conn_existfile(files{n1}{n2}{1})
            if ispc, [nill,nill]=system(sprintf('copy "%s" "%s"',file0,files{n1}{n2}{1}));
            else     [nill,nill]=system(sprintf('cp ''%s'' ''%s''',file0,files{n1}{n2}{1}));
            end
        end
    end
end

