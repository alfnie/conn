function [data,names,coords,samples] = conn_mtx_read(filename)
% CONN_MTX_READ reads matrix data from nifti file
%
% data = conn_mtx_read(filename) reads matrix numeric data from filename
%  filename : input filename *.mtx.nii
%  data     : output data [nrois, nrois, nsamples] 3d matrix
%
% [data ,names, coords, labels] = conn_mtx_read(filename) outputs additional information from sidecar .json file
%  names    : list of ROI names {1,nrois} cell array (default ROI#1,ROI#2,...)
%  coords   : list of ROI xyz coordinates {1,nrois} matrix (default [0 0 0],[0 0 0],...)
%  samples  : list of sample labels {1,nsamples} cell array (default Sample#1,Sample#2,...)
%

if ~nargin, help(mfilename); return; end

vol = spm_vol(filename);
data = permute(spm_read_vols(vol),[1,2,4,3]);
info = conn_jsonread(filename);
if isfield(info,'names'), names=info.names; else names={}; end
if isfield(info,'coords'), coords=info.coords; else coords={}; end
if isfield(info,'samples'), samples=info.samples; else samples={}; end

