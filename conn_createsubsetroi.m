function [fileout,Nrois]=conn_createsubsetroi(filein, fileout, rois)
% CONN_CREATESUBSETROI creates nifti ROI file/atlas selecting from an atlas file a subset of regions
% 
% conn_createsubsetroi(filein, fileout, rois)
%    filein     : input mask filename
%    fileout    : output ROI filename (by default [filein].ROI.nii)
%    rois       : cell array with list of ROIs to keep in output file (see conn_roilabels for a list of ROIs in an atlas file)
%                 (note: each element may in turn be a nested cell array with multiple ROIs to merge multiple ROIs into one)
%

if ischar(filein), filein=cellstr(filein); end
if size(filein,1)>1, filein=reshape(filein,1,[]); end
if nargin<3||isempty(rois)
    rois={};
    for nin=1:numel(filein)
        trois=conn_roilabels(filein{nin});
        if ~isempty(trois), rois=[rois, trois]; end
    end
    if isempty(rois), conn_msgbox({sprintf('No ROI labels found in file %s',sprintf('%s ',filein{:})),'Please select an atlas input file (e.g. with a sidecar .txt file containing ROI labels)'},'',2); return; end
    value=listdlg('liststring',rois,'selectionmode','multiple','initialvalue',[],'promptstring',{'Select one ROI to create a new individual-ROI mask file','(selecting multiple ROIs will create a mask combining them)'},'ListSize',[500 200]);
    if isempty(value), return; end
    rois={rois(value)};
end
if nargin<2||isempty(fileout), fileout=conn_prepend('',filein,'.subset.nii'); end
if any(conn_server('util_isremotefile',[{fileout} filein])), 
    [fileout,Nrois]=conn_server('run',mfilename,conn_server('util_localfile',fileout),conn_server('util_localfile',filein), rois); 
    fileout=conn_server('util_remotefile',fileout);
    return
else
    filein=conn_server('util_localfile',filein);
    fileout=conn_server('util_localfile',fileout);
end

for nin=1:numel(filein)
    fprintf('reading file %s\n',filein{nin});
    if nin==1, [mask{nin},vol]=conn_vol_read(filein{nin});
    else mask{nin}=conn_vol_read(filein{nin},filein{1}); % resamples to first filein file
    end
    original_rois{nin}=conn_roilabels(filein{nin});
    if isempty(original_rois{nin}), [nill,tname,nill]=fileparts(filein{nin}); original_rois{nin}={tname}; end
end

if ischar(rois), rois=cellstr(rois); end

Nrois=0;
output_vol=zeros(size(mask{1}));
output_rois={};
for nroi=1:numel(rois) % each new ROI
    thisroi=rois{nroi};
    if ~iscell(thisroi), thisroi={thisroi}; end % enter cell array with multiple ROIs to merge multiple ROIs into one
    idx=false(size(mask{1}));
    str=[];
    for n=1:numel(thisroi) % multiple elements will be merged into one
        for nin=1:numel(filein) % check on all input files until a match is found
            iroi=find(strcmp(original_rois{nin},thisroi{n})); ok=numel(iroi)>0; 
            if ~ok, iroi=find(strncmp(original_rois{nin},thisroi{n},numel(thisroi{n}))); ok=numel(iroi)>0; end
            if ok, ok=nnz(ismember(mask{nin},iroi))>0; end
            if ok, break; end
        end
        if ~ok, error('unable to find ROI %s in %s',thisroi{n},sprintf('%s ',filein{:})); end
        if numel(iroi)>1, 
            fprintf('Warning: %d ROI names matched to %s (combining all)\n',numel(iroi),thisroi{n}); 
        end
        idx=idx|ismember(mask{nin},iroi);
        if n==numel(thisroi), str=[str,thisroi{n}];
        else str=[str,thisroi{n},'+'];
        end
    end
    Nrois=Nrois+1;
    output_vol(idx)=Nrois;
    output_rois{end+1}=sprintf('%s\n',str);
end

vol=struct('fname',fileout, 'mat',vol.mat, 'dim',vol.dim, 'dt', [spm_type('uint16') spm_platform('bigend')], 'pinfo',[1;0;0],'descrip',''); 
try, spm_unlink(fileout); end
vol=spm_write_vol(vol,output_vol);
conn_fileutils('filewrite_raw',conn_prepend('',fileout,'.txt'), output_rois);
end
