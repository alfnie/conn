function fh=conn_createiconroi(filename)

fh=conn_mesh_display('',{'',filename});
fh('brain_transparency',0);
fh('sub_transparency',0);
fh('zoomout');
fh('background',[0 0 0]);
fh('bookmark',conn_prepend('',filename,'.icon.jpg'));