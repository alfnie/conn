function [vertex_coords, faces] = conn_freesurfer_read_surf(fname)
%
% [vertex_coords, faces] = read_surf(fname)
% reads a the vertex coordinates and face lists from a surface file
% note that reading the faces from a quad file can take a very long
% time due to the goofy format that they are stored in. If the faces
% output variable is not specified, they will not be read so it
% should execute pretty quickly.
%


%
% read_surf.m
%
% Original Author: Bruce Fischl
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:55:10 $
%    $Revision: 1.3 $
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


%fid = fopen(fname, 'r') ;
%nvertices = fscanf(fid, '%d', 1);
%all = fscanf(fid, '%d %f %f %f %f\n', [5, nvertices]) ;
%curv = all(5, :)' ;

% open it as a big-endian file


%QUAD_FILE_MAGIC_NUMBER =  (-1 & 0x00ffffff) ;
%NEW_QUAD_FILE_MAGIC_NUMBER =  (-3 & 0x00ffffff) ;

TRIANGLE_FILE_MAGIC_NUMBER =  16777214 ;
QUAD_FILE_MAGIC_NUMBER =  16777215 ;

fid = fopen(fname, 'rb', 'b') ;
if (fid < 0)
    str = sprintf('could not open curvature file %s.', fname) ;
    error(str) ;
end
magic = conn_freesurfer_fread3(fid) ;

if(magic == QUAD_FILE_MAGIC_NUMBER)
    vnum = conn_freesurfer_fread3(fid) ;
    fnum = conn_freesurfer_fread3(fid) ;
    vertex_coords = fread(fid, vnum*3, 'int16') ./ 100 ;
    if (nargout > 1)
        for i=1:fnum
            for n=1:4
                faces(i,n) = conn_freesurfer_fread3(fid) ;
            end
        end
    end
elseif (magic == TRIANGLE_FILE_MAGIC_NUMBER)
    fgets(fid) ;
    fgets(fid) ;
    vnum = fread(fid, 1, 'int32') ;
    fnum = fread(fid, 1, 'int32') ;
    vertex_coords = fread(fid, vnum*3, 'float32') ;
    faces = fread(fid, fnum*3, 'int32') ;
    faces = reshape(faces, 3, fnum)' ;
end

vertex_coords = reshape(vertex_coords, 3, vnum)' ;

try
    tcoords=reshape(vertex_coords(1+faces,:),[],3,3); 
    if median(sum((tcoords(:,1,:)-repmat(mean(tcoords(:,1,:),1),size(tcoords,1),1)).*(cross(tcoords(:,2,:)-tcoords(:,1,:),tcoords(:,3,:)-tcoords(:,1,:))),3),1)>0, 
        faces=fliplr(faces); 
    end
end
fclose(fid) ;
