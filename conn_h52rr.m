function varargout=conn_h52rr(filename,saveout)
% CONN_H52RR
% converts .h5 timeseries data (timepoints x ROIs)
% into .mtx.nii connectivity matrix (ROIs x ROIs)
% and stores it in an associated .mtx.nii file
%
% e.g. (creates an output task-rest_bold_Atlas_hp2000_clean_GSR_parcellated.mtx.nii file)
% conn_h52rr('task-rest_bold_Atlas_hp2000_clean_GSR_parcellated.h5');
%
% e.g. (without creating an output file)
% [R,data]=conn_h52rr('task-rest_bold_Atlas_hp2000_clean_GSR_parcellated.h5');
%

if nargin<2||isempty(saveout), saveout=nargout==0; end

varargout={};
info=h5info(filename);
DatasetName=info.Datasets(1).Name;
data=h5read(filename,['/',DatasetName]); % timepoints x ROIs
R=atanh(corr(data)); % ROIs x ROIs
R(isnan(R))=0;
R(1:size(R,1)+1:end)=0;
if saveout, conn_mtx_write(conn_prepend('',filename,'.mtx.nii'),R); end
if nargout>0, varargout={R,data}; end

