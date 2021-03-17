function fileout=conn_surf_vol2surf(filein,FSfolder)
% CONN_SURF_VOL2SURF resamples MNI-space volume at [lh|rh].[pial|mid|white].surf surfaces
%
% conn_surf_vol2surf('file.nii')
% creates 'file.surf.nii' file with sapled values at each vertex averaged across the three surfaces and smoothed (10 steps of iterative diffusion smoothing)
%

if nargin<2, FSfolder=''; end
if any(conn_server('util_isremotefile',filein)), fileout=conn_server('util_remotefile',conn_server('run',mfilename,conn_server('util_localfile',filein),FSfolder)); return; end
if isempty(FSfolder), FSfolder=fullfile(fileparts(which('conn')),'utils','surf'); end

tV=conn_surf_extract(filein,{fullfile(FSfolder,'lh.pial.surf'),fullfile(FSfolder,'rh.pial.surf'),fullfile(FSfolder,'lh.mid.surf'),fullfile(FSfolder,'rh.mid.surf'),fullfile(FSfolder,'lh.white.surf'),fullfile(FSfolder,'rh.white.surf')},[],10); % use this for smoother display
%tV=conn_surf_extract(filein,{fullfile(FSfolder,'lh.mid.surf'),fullfile(FSfolder,'rh.mid.surf')});                                                                                                                                          % use this for accurate values
tV=cell2mat(tV);
maskV=isnan(tV)|(tV==0);
tV(maskV)=0;
tV=sum(reshape(tV,[],2,3),3)./max(1,sum(reshape(~maskV,[],2,3),3));

fileout=conn_prepend('',filein,'.surf.nii');
conn_surf_write(fileout,tV(:));
fprintf('created file %s\n',fileout);

