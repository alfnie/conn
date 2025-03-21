function write_curv(fname, curv, fnum)
% write_curv(fname, curv, fnum)
%
% writes a curvature vector into a binary file
%				fname - name of file to write to
%				curv  - vector of curvatures
%				fnum  - # of faces in surface.
%


%
% write_curv.m
%
% Original Author: Bruce Fischl
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:55:10 $
%    $Revision: 1.2 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA). 
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
%

if nargin<3||isempty(fnum), fnum=327680; end % alfnie added default (res=8)
if any(conn_server('util_isremotefile',fname)), conn_server('run',mfilename,conn_server('util_localfile',fname), curv, fnum); return; end
fname=conn_server('util_localfile',fname);

% open it as a big-endian file
fid = fopen(fname, 'wb', 'b') ;
vnum = length(curv) ;
NEW_VERSION_MAGIC_NUMBER = 16777215;
conn_freesurfer_fwrite3(fid, NEW_VERSION_MAGIC_NUMBER ) ;
fwrite(fid, vnum,'int32') ;
fwrite(fid, fnum,'int32') ;
fwrite(fid, 1, 'int32');
fwrite(fid, curv, 'float') ;
fclose(fid) ;

