function filename=conn_mgh2nii(filename)
% CONN_MGH2NII converts MGH .mgh/.mgz files to NIFTI .nii format
%
% filename = conn_mgh2nii(filename)
%

if any(conn_server('util_isremotefile',filename)), filename=conn_server('util_remotefile',conn_server('run',mfilename,conn_server('util_localfile',filename))); return; end
ischarfilename=ischar(filename);
filename=cellfun(@strtrim,cellstr(filename),'uni',0);
filenameout=regexprep(filename,'\.mgh\s*$|\.mgz\s*$','.nii');
redo=~cellfun(@conn_existfile,filenameout);
if ~any(redo), filename=filenameout; if ischarfilename, filename=char(filename); end; return; end

fprintf('converting mgh files to nifti format...');
[pathname,name,ext]=spm_fileparts(filename{1});
if strcmp(ext,'.mgz')
    filename(redo)=gunzip(filename(redo));
    if ispc, [ok,msg]=cellfun(@(x)system(sprintf('move "%s" "%s.mgh"',x,x)),filename(redo),'uni',0);
    else     [ok,msg]=cellfun(@(x)system(sprintf('mv -f ''%s'' ''%s.mgh''',x,x)),filename(redo),'uni',0);
    end
    if ~all(cellfun(@(x)isequal(x,0),ok)), error(['error converting mgh to nifti format. Please check file permissions in folder ',pathname]); end
    filename(redo)=cellfun(@(x)sprintf('%s.mgh',x),filename(redo),'uni',0);
end
for n=find(redo(:)')
    a=conn_freesurfer_MRIread(filename{n});
    b=permute(a.vol,[2,1,3,4]);
    dt=[spm_type('float32') spm_platform('bigend')];
    if ~nnz(rem(b,1)~=0)
        if ~nnz(b<0),
            if ~nnz(b>255),         dt(1)=spm_type('uint8');
            elseif ~nnz(b>65535),   dt(1)=spm_type('uint16');
            end
        end
    end
    V=struct('mat',a.vox2ras1,'dim',a.volsize([2 1 3]),'dt',dt,'fname',filenameout{n},'pinfo',[1;0;0]);
    try
        sb=[size(b) 1 1 1 1];
        if prod(sb(1:3))==prod(conn_surf_dims(8)) && nnz(sb(1:3)~=1)==1
            b=reshape(b,[conn_surf_dims(8) sb(4:end)]);
            V.dim=conn_surf_dims(8);
            V.mat=eye(4);
        end
        if size(b,4)>1
            V=repmat(V,[size(b,4),1]);for nh=1:size(b,4),V(nh).n=[nh,1];end
            V=spm_create_vol(V);
            for nh=1:size(b,4), V(nh)=spm_write_vol(V(nh),b(:,:,:,nh)); end
        else
            spm_write_vol(V,b);
        end
    catch
        error(['error writing nifti file. Please check file permissions in folder ',pathname]); 
    end
    [tpath,tname,text]=spm_fileparts(filenameout{n});
    if ~isempty(regexp(tname,'^aparc.*\+aseg$'))
        tfile=conn_prepend('',filenameout{n},'.txt');
        file0=fullfile(fileparts(which(mfilename)),'utils','surf','FreeSurferColorLUT.txt');
        if ~conn_existfile(tfile), try, conn_fileutils('copyfile',file0,tfile); end; end
    end
end
filename=filenameout; 
if ischarfilename, filename=char(filename); end; 
fprintf('done\n');

