function filename=conn_surf_resample(filename,FS_folder)
% conn_surf_resample 
% resample functional data at the location of FreeSurfer subject-specific structural cortical surface 
%

tfilename=cellstr(filename);
refinfo=[];
tfileout={};
for n=1:numel(tfilename)
    [nill,tfileout{n},nill,refinfo]=conn_surf_extract(tfilename{n},[],FS_folder,0,false,true,refinfo);
end
if ischar(filename), filename=char(tfileout);
else filename=tfileout;
end



