function [vertex_coords, faces] = write_surf(fname,vertex_coords,faces)

TRIANGLE_FILE_MAGIC_NUMBER =  16777214 ;
fid = fopen(fname, 'wb', 'b') ;
if (fid < 0)
    str = sprintf('could not create curvature file %s.', fname) ;
    error(str) ;
end
conn_freesurfer_fwrite3(fid,TRIANGLE_FILE_MAGIC_NUMBER);
fclose(fid);
fid=fopen(fname,'at','b');
fprintf(fid,'\n');
fprintf(fid,'\n');
fclose(fid);
fid=fopen(fname,'ab','b');
vnum = size(vertex_coords,1) ;
fnum = size(faces,1);
fwrite(fid,vnum,'int32');
fwrite(fid,fnum,'int32');
fwrite(fid,vertex_coords','float32');
fwrite(fid,faces','int32');
fclose(fid);
