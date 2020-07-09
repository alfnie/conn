function [data,fileout,fileout2,refinfo]=conn_surf_extract(filename,surfnames,FS_folder,smooth,DODISP,DOSAVE,refinfo,refidxcl)

if nargin<3, FS_folder=''; end
if nargin<4||isempty(smooth), smooth=0; end
if nargin<5||isempty(DODISP), DODISP=false; end
if nargin<6||isempty(DOSAVE), DOSAVE=false; end
if nargin<7||isempty(refinfo), refinfo={}; refinfoinit=true; 
else refinfoinit=false; 
end
if nargin<8, refidxcl='all'; end
fileout='';
fileout2='';
if isempty(FS_folder) % extracts from single surface file (assumes all files in MNI space)
    if DODISP
        [file_path,file_name,file_ext]=fileparts(surfnames);
        if isempty(regexp([file_name,file_ext],'^lh|^rh'))||~isempty(regexp([file_name,file_ext],'subcortical|cerebellum'))||any(strcmp([file_name,file_ext],{'lh.pial.surf','rh.pial.surf'}))
            surfnames={surfnames};
        else
            surfnames={fullfile(fileparts(which(mfilename)),'surf',[file_name(1:3),'pial.surf']),surfnames};
        end
        for n1=1:numel(surfnames),[tfile_path,tfile_name,tfile_ext]=fileparts(surfnames{n1}); surfnames_redux{n1}=[tfile_name,tfile_ext]; end
        hfig=figure('units','norm','position',[.4,.5,.3,.2],'color','w','name','Extract values at surface coordinates','numbertitle','off','menubar','none');
        uicontrol('style','text','units','norm','position',[.1,.8,.8,.1],'string','Extract values at coordinates:','backgroundcolor','w','horizontalalignment','left');
        ht1=uicontrol('style','popupmenu','units','norm','position',[.1,.7,.8,.1],'string',surfnames_redux,'value',1);
        uicontrol('style','text','units','norm','position',[.1,.45,.8,.1],'string','Smoothing level:','backgroundcolor','w','horizontalalignment','left');
        ht2=uicontrol('style','edit','units','norm','position',[.1,.35,.8,.1],'string',num2str(smooth),'tooltipstring','Enter radius (number of vertices) for surface-smoothing of extracted values');
        uicontrol('style','pushbutton','string','OK','units','norm','position',[.1,.01,.38,.15],'callback','uiresume');
        uicontrol('style','pushbutton','string','Cancel','units','norm','position',[.51,.01,.38,.15],'callback','delete(gcbf)');
        uiwait(hfig);
        if ishandle(hfig)
            surfnames=surfnames{get(ht1,'value')};
            smooth=str2num(get(ht2,'string'));
            delete(hfig);
        else
            data=[];
            return;
        end
    end
    voldata=[];
    if isstruct(filename)&&isfield(filename,'fname'), vol=filename;
    elseif isstruct(filename)&&isfield(filename,'data'), vol=filename.vol; voldata=filename.data;
    else  vol=spm_vol(filename);
    end
    if ~iscell(surfnames),surfnames={surfnames};end
    data={};
    for i=1:numel(surfnames)
        if isstruct(surfnames{i})
            xyz=surfnames{i}.vertices;
            faces=surfnames{i}.faces;
        else
            [xyz,faces]=conn_freesurfer_read_surf(surfnames{i});
            faces=faces+1;
        end
        %ijk=pinv(vol.mat)*a.vox2ras0*pinv(a.tkrvox2ras)*[xyz,ones(size(xyz,1),1)]';
        % note: a.vox2ras0*pinv(a.tkrvox2ras) = eye(4)
        ijk=pinv(vol(1).mat)*[xyz,ones(size(xyz,1),1)]';
        if isempty(voldata), data_ref=max(spm_get_data(vol,ijk),[],1); % reads from file
        else % reads from 3d matrix directly
            valid=all(ijk(1:3,:)>=1,1)&ijk(1,:)<=size(voldata,1)&ijk(2,:)<=size(voldata,2)&ijk(3,:)<=size(voldata,3);
            data_ref=zeros(size(voldata,4),size(ijk,2));
            for ndata=1:size(voldata,4), data_ref(ndata,valid)=voldata([1,size(voldata,1),size(voldata,1)*size(voldata,2)]*(round(ijk(1:3,valid))-1)+1); end
            data_ref=max(data_ref,[],1);
        end
        data_ref=data_ref(:);
        if smooth>0
            data_ref(isnan(data_ref))=0;
            A=double(sparse(repmat(faces,3,1),repmat(faces,1,3), 1)>0);
            A=double(A|speye(size(A,1)));
            A=A*sparse(1:size(A,1),1:size(A,1),1./sum(A,2));
            mdata_ref=double(data_ref~=0);
            for n=1:smooth,
                data_ref=A*data_ref; 
                mdata_ref=A*mdata_ref; 
            end
            data_ref(mdata_ref<.5)=0;
        end
        data{i}=data_ref;
    end
    if numel(surfnames)==1, data=data{1}; end
    
else % extracts from arbitrary freesurfer surfaces (subject-specific coordinates)
    if refinfoinit
        if isstruct(FS_folder)
            vol=FS_folder;
            FS_folder=fileparts(vol.fname);
            [temp1,temp2]=spm_fileparts(FS_folder);
            if strcmp(temp2,'mri')||strcmp(temp2,'anat'), FS_folder=temp1; end
        else
            names={'T1.nii','brain.nii','T1.mgh','brain.mgh','T1.mgz','brain.mgz'};
            if ~(ischar(FS_folder)&&isdir(FS_folder))
                [FS_folder,names_name,names_ext]=spm_fileparts(FS_folder);
                if isempty(FS_folder), FS_folder=pwd; end
                %names={[names_name,names_ext]};
            end
            [temp1,temp2]=spm_fileparts(FS_folder);
            if strcmp(temp2,'mri')||strcmp(temp2,'anat'), FS_folder=temp1; end
            valid=cellfun(@(name)~isempty(dir(fullfile(FS_folder,'mri',name))),names);
            if ~any(valid),
                if ~nargout, error('missing anatomical data');
                else
                    conn_disp('missing anatomical data'); data=[]; fileout=''; return;
                end
            end
            ivalid=find(valid,1);
            if ~isempty(regexp(names{ivalid},'\.mgh\s*$|\.mgz\s*$')), tname=conn_mgh2nii(fullfile(FS_folder,'mri',names{ivalid})); 
            else tname=fullfile(FS_folder,'mri',names{ivalid});
            end
            vol=spm_vol(tname);
            %a=MRIread(fullfile(FS_folder,'mri',names{find(valid,1)}),true); % note: when using spm_vol, a.tkrvox2ras=vox2ras_tkreg(a.volsize,a.volres)
        end
        a.vox2ras1=vol.mat;
        a.volsize=vol.dim([2 1 3]);
        a.volres = sqrt(sum(vol.mat(:,1:3).^2,1));
        a.vox2ras0=conn_freesurfer_vox2ras_1to0(vol.mat);
        a.tkrvox2ras=conn_freesurfer_vox2ras_tkreg(a.volsize,a.volres);
    end
    alpha=.05:.1:.95;
    DOMEAN=1; % 1:mean; 2:median; 3:mode
    if nargin<2||isempty(surfnames), surfnames={'.white','.pial','.sphere.reg'}; end
    resolution=8;
    hems={'lh','rh'};
    voldata=[];
    if isstruct(filename)&&isfield(filename,'fname')
        vol=filename;
        [file_path,file_name,file_ext,file_num]=spm_fileparts(vol.fname);
    elseif isstruct(filename)&&isfield(filename,'data'), 
        vol=filename.vol; 
        voldata=filename.data;
        [file_path,file_name,file_ext,file_num]=spm_fileparts(vol.fname);
    else
        vol=spm_vol(filename);
        [file_path,file_name,file_ext,file_num]=spm_fileparts(filename);
    end
    data=[];
    dataraw=[];
    alpha=permute(alpha(:),[2,3,1]);
    for hem=1:2,
        shem=hems{hem};
        if refinfoinit
            [xyz1,faces]=conn_freesurfer_read_surf(fullfile(FS_folder,'surf',[shem,surfnames{1}]));
            faces=faces+1;
            if numel(surfnames)>=2&&~isempty(surfnames{2})
                xyz2=conn_freesurfer_read_surf(fullfile(FS_folder,'surf',[shem,surfnames{2}]));
                xyz_data=conn_bsxfun(@times,alpha,xyz1')+conn_bsxfun(@times,1-alpha,xyz2');
            else
                xyz_data=xyz1';
            end
            refinfo{hem}.xyz=a.vox2ras0*pinv(a.tkrvox2ras)*[xyz_data(:,:);ones(1,size(xyz_data,2)*size(xyz_data,3))];
            refinfo{hem}.sxyz=[size(xyz_data,1),size(xyz_data,2),size(xyz_data,3)];
            if numel(surfnames)>=3&&~isempty(surfnames{3}) % resample at sphere reference grid
                xyz_ref=conn_freesurfer_read_surf(fullfile(FS_folder,'surf',[shem,surfnames{3}]));
                [xyz_sphere,sphere2ref,ref2sphere]=conn_surf_sphere(resolution,xyz_ref);
                refinfo{hem}.xyz_sphere=xyz_sphere;
                refinfo{hem}.ref2sphere=ref2sphere;
                refinfo{hem}.sphere2ref=sphere2ref;
            end
        end
        if ~isequal(refidxcl,'all')
            if hem==1, trefidxcl=refidxcl(refidxcl<=numel(refinfo{hem}.ref2sphere));
            else       trefidxcl=refidxcl(refidxcl>numel(refinfo{hem}.ref2sphere))-numel(refinfo{hem}.ref2sphere);
            end
            if numel(surfnames)>=3&&~isempty(surfnames{3}) % resample at sphere reference grid
                trefidxcl=refinfo{hem}.ref2sphere(trefidxcl);
            end
            ijkmm=refinfo{hem}.xyz(:,conn_bsxfun(@plus,trefidxcl(:),refinfo{hem}.sxyz(2)*(0:refinfo{hem}.sxyz(3)-1)));
            ijk=pinv(vol(1).mat)*ijkmm;
        else
            ijkmm=refinfo{hem}.xyz;
            ijk=pinv(vol(1).mat)*ijkmm;
        end
        voldone=false; nvol=1; data_ref_log={};
        while ~voldone
            if isempty(voldata)&&numel(vol)==1, 
                data_ref=spm_get_data(vol,ijk); 
                voldone=true; % reads from file
            elseif isempty(voldata),
                data_ref=spm_get_data(vol(nvol),pinv(vol(nvol).mat)*ijkmm); % reads from file
                voldone=(nvol>=numel(vol));
            else % reads from 3d matrix directly
                valid=all(ijk(1:3,:)>=1,1)&ijk(1,:)<=size(voldata,1)&ijk(2,:)<=size(voldata,2)&ijk(3,:)<=size(voldata,3);
                data_ref=zeros(size(voldata,4),size(ijk,2));
                for ndata=1:size(voldata,4), data_ref(ndata,valid)=voldata([1,size(voldata,1),size(voldata,1)*size(voldata,2)]*(round(ijk(1:3,valid))-1)+1); end
                voldone=true;
            end
            if numel(alpha)>1
                if DOMEAN==1 % mean ignoring nan's
                    data_ref=reshape(data_ref,[],refinfo{hem}.sxyz(3))';
                    idata_ref=isnan(data_ref)|(data_ref==0);
                    data_ref(idata_ref)=0;
                    data_ref=sum(data_ref,1)./max(eps,sum(~idata_ref,1));
                elseif DOMEAN==2 % median ignoring nan's
                    data_ref=sort(reshape(data_ref,[refinfo{hem}.sxyz(2),refinfo{hem}.sxyz(3)]),2);
                    idata_ref=max(1,ceil(sum(~isnan(data_ref),2)/2));
                    data_ref=data_ref((1:size(data_ref,1))'+size(data_ref,1)*(idata_ref-1))';
                else
                    data_ref=reshape(data_ref,[],refinfo{hem}.sxyz(3))';
                    data_ref=mode(data_ref,1);
                end
            end
            if numel(surfnames)>=3&&~isempty(surfnames{3}) % resample at sphere reference grid
                if isequal(refidxcl,'all')
                    data_ref=data_ref(:,refinfo{hem}.ref2sphere);
                end
                faces=refinfo{hem}.xyz_sphere.faces;
            end
            if isempty(voldata)&&numel(vol)>1
                data_ref_log{nvol}=data_ref;
                nvol=nvol+1;
            end
        end
        if ~isempty(data_ref_log), data_ref=cat(1,data_ref_log{:}); end
        dataraw_ref=data_ref;
        if smooth>0 % smooths on surface
            if ~isequal(refidxcl,'all'), error('Surface smoothing option not supported when extracting only data from a subset of vertices'); end
            data_ref(isnan(data_ref))=0;
            mask=double(data_ref~=0);
            if refinfoinit||~isfield(refinfo{hem},'A')||isempty(refinfo{hem}.A)
                A=double(sparse(repmat(faces,3,1),repmat(faces,1,3), 1)>0);
                A=double(A|speye(size(A,1)));
                %A=A*sparse(1:size(A,1),1:size(A,1),1./sum(A,2));
                refinfo{hem}.A=A;
            else
                A=refinfo{hem}.A;
            end
            for n=1:smooth,
                %mask=double(data_ref~=0);
                data_ref=(data_ref*A')./max(eps,mask*A'); 
            end
        end
%         % removes outliers
%         data_localmsk=ones(size(data_ref));
%         data_localavg=data_ref;
%         data_localvar=data_ref.^2;
%         for n=1:32,
%             mask=double(data_localavg~=0)*A';
%             data_localavg=(data_localavg*A')./max(eps,mask);
%             data_localvar=(data_localvar*A')./max(eps,mask);
%             data_localmsk=(data_localmsk*A')./max(eps,mask);
%         end
%         data_localavg=data_localavg./max(eps,data_localmsk);
%         data_localvar=data_localvar./max(eps,data_localmsk);
%         data_localvar=data_localvar-data_localavg.^2;
%         if refinfoinit||~isfield(refinfo{hem},'A2')||isempty(refinfo{hem}.A2)
%             A2=A;
%             for n=1:8,A2=A2*A2;end
%             A2=A2*sparse(1:size(A2,1),1:size(A2,1),1./sum(A2,2));
%             refinfo{hem}.A2=A2;
%         else
%             A2=refinfo{hem}.A2;
%         end
%         for n=1:smooth,
%             mask=double(data_ref~=0);
%             data_ref=(data_ref*A')./max(eps,mask*A');
%         end
        % saves data
        if 0,%numel(surfnames)>=3&&~isempty(surfnames{3})&&DOSAVE % saves .lh / .rh individual files
            if ~isequal(refidxcl,'all'), error('Volume-save option not supported when extracting only data from a subset of vertices'); end
            dim=conn_surf_dims(resolution);
            if resolution~=8, sext=['.surf',num2str(resolution),file_ext];
            else sext=['.surf',file_ext];
            end
            if ~isempty(file_num), sext=['.',regexprep(file_num,'\D',''),sext]; end
            newvol=struct('fname',fullfile(file_path,[file_name,'.',shem,sext]),...
                'mat',eye(4),...
                'dim',dim,...
                'pinfo',[1;0;0],...
                'dt',[spm_type('float32'),spm_platform('bigend')],...
                'descrip','surface data');
            spm_write_vol(newvol,reshape(dataraw_ref,dim));
            %conn_disp(['Created file ',newvol.fname]);
            if smooth>0
                newvol=struct('fname',fullfile(file_path,[file_name,'.',shem,'.smooth',num2str(smooth),sext]),...
                    'mat',eye(4),...
                    'dim',dim,...
                    'pinfo',[1;0;0],...
                    'dt',[spm_type('float32'),spm_platform('bigend')],...
                    'descrip','surface data');
                spm_write_vol(newvol,reshape(data_ref,dim));
                %conn_disp(['Created file ',newvol.fname]);
            end
        end
        data=[data,data_ref];
        dataraw=[dataraw,dataraw_ref];
    end
    if numel(surfnames)>=3&&~isempty(surfnames{3})&&DOSAVE
        if ~isequal(refidxcl,'all'), error('Volume-save option not supported when extracting only data from a subset of vertices'); end
        dim=conn_surf_dims(resolution);
        dim=dim.*[1 1 2];
        if resolution~=8, sext=['.surf',num2str(resolution),file_ext];
        else sext=['.surf',file_ext];
        end
        if ~isempty(file_num), sext=['.',regexprep(file_num,'\D',''),sext]; end
        newvol=struct('fname',fullfile(file_path,[file_name,sext]),...
            'mat',eye(4),...
            'dim',dim,...
            'pinfo',[1;0;0],...
            'dt',[spm_type('float32'),spm_platform('bigend')],...
            'descrip','surface data');
        if ~isempty(dir(newvol.fname)), try, spm_unlink(newvol.fname); end; end
        if size(dataraw,1)>1
            newvol=repmat(newvol,[size(dataraw,1),1]);for nh=1:size(dataraw,1),newvol(nh).n=[nh,1];end
            newvol=spm_create_vol(newvol);
            for nh=1:size(dataraw,1), newvol(nh)=spm_write_vol(newvol(nh),reshape(dataraw(nh,:),dim)); end
        else spm_write_vol(newvol,reshape(dataraw,dim));
        end
        %conn_disp(['Created file ',newvol.fname]);
        fileout=newvol.fname;
        fileout2=newvol.fname;
        if smooth>0
            newvol=struct('fname',fullfile(file_path,[file_name,'.smooth',num2str(smooth),sext]),...
                'mat',eye(4),...
                'dim',dim,...
                'pinfo',[1;0;0],...
                'dt',[spm_type('float32'),spm_platform('bigend')],...
                'descrip','surface data');
            if ~isempty(dir(newvol.fname)), try, spm_unlink(newvol.fname); end; end
            if size(dataraw,1)>1
                newvol=repmat(newvol,[size(data,1),1]);for nh=1:size(data,1),newvol(nh).n=[nh,1];end
                newvol=spm_create_vol(newvol);
                for nh=1:size(data,1), newvol(nh)=spm_write_vol(newvol(nh),reshape(data(nh,:),dim)); end
            else spm_write_vol(newvol,reshape(data,dim));
            end
            %conn_disp(['Created file ',newvol.fname]);
            fileout=newvol.fname;
        end
    end
end
end

