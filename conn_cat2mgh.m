function conn_cat2mgh(filename)
% CONN_CAT2MGH converts CAT12 output files to MGH format (only files necessary by CONN)
% conn_cat2mgh(filename)
%   filename : T1 image processed by CAT12
%         
% example:
% conn_cat2mgh('/data/T1.nii')
%   inputs (expected files)
%       /data/T1.nii
%       /data/surf/?h.central.T1.gii
%       /data/surf/?h.sphere.reg.T1.gii
%       /data/surf/?h.thickness.T1 files
%   outputs:
%      /data/surf/?h.pial
%      /data/surf/?h.white
%      /data/surf/?h.sphere.reg
%      /data/mri/T1.nii
%

if any(conn_server('util_isremotefile',filename)), conn_server('run',mfilename,conn_server('util_localfile',filename)); return; 
else filename=conn_server('util_localfile',filename);
end 
if iscell(filename)
    cellfun(@conn_cat2mgh,filename,'uni',0);
else
    [filepath,filename,fileext]=fileparts(conn_fullfile(filename));
    vol=spm_vol(fullfile(filepath,[filename,fileext]));
    data=spm_read_vols(vol);
    a.vox2ras1=vol.mat;
    a.volsize=vol.dim([2 1 3]);
    a.volres = sqrt(sum(vol.mat(:,1:3).^2,1));
    a.vox2ras0=conn_freesurfer_vox2ras_1to0(vol.mat);
    a.tkrvox2ras=conn_freesurfer_vox2ras_tkreg(a.volsize,a.volres);
    
    vout=struct('fname',fullfile(filepath,'mri','T1.nii'),'mat',vol.mat,'dim',vol.dim,'n',[1,1],'pinfo',[1;0;0],'dt',vol.dt,'descrip','FreeSurfer conn_cat2mgh');
    spm_write_vol(vout,data);
    
    for lhrh={'lh','rh'}
        data=gifti(fullfile(filepath,'surf',[lhrh{1},'.central.',filename,'.gii']));
        vertices=a.tkrvox2ras(1:3,:)*pinv(a.vox2ras0)*[double(data.vertices), ones(size(data.vertices,1),1)]'; %to mgh coordinates (xyzvol' = a.vox2ras0*pinv(a.tkrvox2ras)*xyzsurf')
        vertices=vertices(1:3,:)';
        
        th=conn_freesurfer_read_curv(fullfile(filepath,'surf',[lhrh{1},'.thickness.',filename]));
        vn=vertexNormal(triangulation(double(data.faces),vertices));
        conn_freesurfer_write_surf(fullfile(filepath,'surf',[lhrh{1},'.white']),vertices-repmat(th/2,[1 3]).*vn,data.faces-1);
        conn_freesurfer_write_surf(fullfile(filepath,'surf',[lhrh{1},'.pial']),vertices+repmat(th/2,[1 3]).*vn,data.faces-1);

        data=gifti(fullfile(filepath,'surf',[lhrh{1},'.sphere.reg.',filename,'.gii']));
        conn_freesurfer_write_surf(fullfile(filepath,'surf',[lhrh{1},'.sphere.reg']),vertices,data.faces-1);
    end
end
end










