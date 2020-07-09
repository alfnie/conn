function conn_setup_preproc_reorientcombine(varargin)
% conn_setup_preproc_reorientcombine combines multiple sequential reorient_*.mat files into a single file
%
% example: 
%  conn_setup_preproc_reorientcombine reorient_out.mat reorient_step1.mat reorient_step2.mat  reorient_step3.mat
%

R0=eye(4);
for n=2:nargin
    load(varargin{n},'R');
    if all(size(R)==3), R=[R zeros(3,1); zeros(1,3) 1]; end
    R0=R*R0;
end
R=R0;
save(varargin{1},'R','-mat');

