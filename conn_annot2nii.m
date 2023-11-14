function niifilename = conn_annot2nii(varargin)
% CONN_ANNOT2NII lh.filename.annot rh.filename.annot
% converts freesurfer .annot file into .nii/.txt surface nifti ROI file
%
% filename = conn_annot2nii(filename)
%


FORCEREDO=false; % forces creation of target file even if it already exists

if nargin==1&&ischar(varargin{1}), roifiles=varargin;
elseif nargin==1&&iscell(varargin{1}), roifiles=varargin{1};
else roifiles=varargin;
end
if any(conn_server('util_isremotefile',roifiles)), niifilename=conn_server('util_remotefile',conn_server('run',mfilename,conn_server('util_localfile',roifiles))); return; 
else roifiles=conn_server('util_localfile',roifiles);
end

if numel(roifiles)==1
    [file_path,file_name,file_ext]=fileparts(roifiles{1});
    if strncmp(file_name,'lh.',3), roifiles=[roifiles, {fullfile(file_path,[regexprep(file_name,'^lh\.','rh.') file_ext])}]; 
    elseif strncmp(file_name,'rh.',3), roifiles=[roifiles, {fullfile(file_path,[regexprep(file_name,'^rh\.','lh.') file_ext])}]; 
    end
end

log=struct('name',{{}},'hem',[],'data',{{}},'labels',{{}},'folder',{{}});
for nfile=1:numel(roifiles)
    % gets info from .annot file
    roifile=deblank(roifiles{nfile});
    [temp_vert,temp_label,temp_table]=conn_freesurfer_read_annotation(roifile,0);
    [nill,temp_rois]=ismember(temp_label,temp_table.table(:,5));
    temp_colors=temp_table.table(:,1:3)/255;
    names_rois=temp_table.struct_names;
    
    [file_path,file_name,file_ext]=fileparts(roifile);
    if isempty(file_path), file_path=pwd; end
    if strcmp(file_name(1:2),'lh'), lhrh=1; 
    elseif strcmp(file_name(1:2),'rh'), lhrh=2; 
    else error('unrecognized file naming convention (lh.* rh.* filenames expected)');
    end
    
    log.name{end+1}=file_name(4:end);
    log.hem(end+1)=lhrh;
    log.data{end+1}=temp_rois;
    log.labels{end+1}=names_rois;
    log.folder{end+1}=file_path;
end

% creates associated .nii / .txt files
niifilename={};
mismatchedsize=[];
for nfile1=1:numel(log.name)
    for nfile2=nfile1+1:numel(log.name)
        if strcmp(log.name{nfile1},log.name{nfile2})&&isequal(sort(log.hem([nfile1 nfile2])),[1 2])
            ifile=[nfile1,nfile2];
            fname=fullfile(file_path,[log.name{ifile(1)},'.surf.nii']);
            [nill,idx]=sort(log.hem(ifile));
            ifile=ifile(idx);
            dim0=conn_surf_dims(8);
            dim=dim0.*[1 1 2];
            converted=true;
            FS_folder=[];
            if (numel(log.data{ifile(1)})~=prod(dim0)||numel(log.data{ifile(2)})~=prod(dim0))&&isequal(log.folder{ifile(1)},log.folder{ifile(2)}) % converts subject-space to fsaverage
                [fpath1,fname1,fext1]=fileparts(log.folder{ifile(1)});
                if strcmp(fname1,'label')&&conn_existfile(fullfile(fpath1,'surf','lh.sphere.reg'))&&conn_existfile(fullfile(fpath1,'surf','rh.sphere.reg')), FS_folder=fullfile(fpath1,'surf');
                elseif conn_existfile(fullfile(log.folder{ifile(1)},'lh.sphere.reg'))&&conn_existfile(fullfile(log.folder{ifile(1)},'rh.sphere.reg')), FS_folder=log.folder{ifile(1)};
                end
                if ~isempty(FS_folder)
                    xyz_ref1=conn_freesurfer_read_surf(fullfile(FS_folder,'lh.sphere.reg'));
                    xyz_ref2=conn_freesurfer_read_surf(fullfile(FS_folder,'rh.sphere.reg'));
                    if size(xyz_ref1,1)==numel(log.data{ifile(1)})&&size(xyz_ref2,1)==numel(log.data{ifile(2)}), converted=false; end
                end
            end
            if ~converted||numel(log.data{ifile(1)})+numel(log.data{ifile(2)})==prod(dim)
                niifilename{end+1}=fname;
                if FORCEREDO||~conn_existfile(fname),
                    if ~converted
                        [xyz_sphere,sphere2ref,ref2sphere]=conn_surf_sphere([],xyz_ref1);
                        log.data{ifile(1)}=log.data{ifile(1)}(ref2sphere);
                        [xyz_sphere,sphere2ref,ref2sphere]=conn_surf_sphere([],xyz_ref2);
                        log.data{ifile(2)}=log.data{ifile(2)}(ref2sphere);
                    end
                    assert(numel(log.data{ifile(1)})==numel(log.data{ifile(2)}),'Annotation file defined in subject-space (%d&%d vertices). Only .annot files in freesurfer fsaverage tessellation allowed (unless ?h.sphere.reg data available)',numel(log.data{ifile(1)}),numel(log.data{ifile(2)}));
                    data=[log.data{ifile(1)}(:) log.data{ifile(2)}(:)];
                    names_rois=log.labels{ifile(1)}; none=find(strncmp('None',names_rois,4)); data(ismember(data(:,1),none),1)=0;
                    names_rois=log.labels{ifile(2)}; none=find(strncmp('None',names_rois,4)); data(ismember(data(:,2),none),2)=0;
                    data(data(:,2)>0,2)=numel(log.labels{ifile(1)})+data(data(:,2)>0,2);
                    if isequal(log.labels{ifile(1)},log.labels{ifile(2)})
                        names_rois=[cellfun(@(x)[x ' (L)'],names_rois,'uni',0); cellfun(@(x)[x ' (R)'],names_rois,'uni',0)];
                    else 
                        %names_rois=[log.labels{ifile(1)}; log.labels{ifile(2)}];
                        names_rois=[cellfun(@(x)[x ' (L)'],regexprep(log.labels{ifile(1)},'^l\.','','ignorecase'),'uni',0); cellfun(@(x)[x ' (R)'],regexprep(log.labels{ifile(2)},'^r\.','','ignorecase'),'uni',0)];
                    end
                    V=struct('mat',eye(4),'dim',dim,'pinfo',[1;0;0],'fname',fname,'dt',[spm_type('uint16') spm_platform('bigend')]);
                    spm_write_vol(V,reshape(data,dim));
                    fprintf('Created file %s\n',fname);
                    fname=fullfile(file_path,[log.name{ifile(1)},'.surf.txt']);
                    fh=fopen(fname,'wt');
                    for n=1:max(data(:))
                        fprintf(fh,'%s\n',names_rois{n});
                    end
                    fclose(fh);
                    try, conn_surf_surf2vol(fname,[],FS_folder,.5); end
                end
                %fprintf('Created file %s\n',fname);
            else mismatchedsize=numel(log.data{ifile(1)})+numel(log.data{ifile(2)});
            end
        end
    end
end
if isempty(niifilename)&&isempty(mismatchedsize), error('nii file not created'); 
elseif isempty(niifilename), error('nii file not created. Number of vertices (%d) does not match freesurfer fsaverage tessellation (%d)',mismatchedsize,prod(dim)); 
end
niifilename=char(niifilename);
