function fh=conn_meshroi_display(filename)

fh=conn_mesh_display('',{'',filename});
fh('brain_transparency',0);
fh('sub_transparency',0);
