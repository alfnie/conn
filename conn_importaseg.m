function varargout=conn_importaseg(filename,filelabels,justchecking)

global CONN_x;
if ~nargin
    for nsub=1:CONN_x.Setup.nsubjects
        for nses=1:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsub))
            filename=CONN_x.Setup.structural{nsub}{nses}{1};
            if conn_importaseg(fileparts(filename),[],true)
                filenames=conn_importaseg(fileparts(filename));
                for nseg=1:3
                    CONN_x.Setup.rois.files{nsub}{nseg}{nses}=conn_file(filenames{nseg});
                end
            end
        end
    end
    return
end

if nargin<3||isempty(justchecking), justchecking=false; end
if nargin<2||isempty(filelabels),  filelabels=fullfile(fileparts(which(mfilename)),'utils','surf','FreeSurferColorLUT.txt'); end
if ~justchecking&&any(conn_server('util_isremotefile',filename)), varargout={conn_server('run',mfilename,conn_server('util_localfile',filename),filelabels,justchecking)}; return; end
filename=conn_server('util_localfile',filename);

if conn_fileutils('isdir',filename), filename=fullfile(filename,'aseg.mgz'); end
[file_path,file_name,file_ext,file_num]=spm_fileparts(filename);
filenames=arrayfun(@(n)fullfile(file_path,sprintf('c%d_%s.img',n,file_name)),1:3,'uni',0);
ok1=conn_existfile(filename);
ok2=ok1&all(cellfun(@(x)conn_existfile(x),filenames));

if justchecking
    varargout={ok1,ok2};
else
    varargout={filenames};
    if ~ok1, error(['File ',filename,' not found']); 
    elseif ok2, 
        conn_disp('fprintf','Note: Files %s already exist. Skipping conversion\n',sprintf('%s ',filenames{:}));
    else
        [nill,nill,filename]=conn_file(filename)
        a=spm_vol(filename);
        b=spm_read_vols(a);
        idxvoxels=b>0;
        [ub,nill,iub]=unique(b(idxvoxels));
        
        [id,PU]=textread(filelabels,'%s%s%*[^\n]','delimiter',' \t'); % FreeSurfer *LUT.txt or equivalent format (ROI_NUMBER ROI_LABEL)
        id0=str2double(id);
        idxnull=find(isnan(id0));
        if ~numel(idxnull)||numel(strmatch('#',id(idxnull)))==numel(idxnull)
            idxvalid=find(~isnan(id0));
            id=id0(idxvalid);
            PU=PU(idxvalid);
            b=round(b);ub=round(ub);
            XYZnames=cell(1,numel(ub));
            for nub=1:numel(ub),
                nid1=find(id==ub(nub),1);
                if ~isempty(nid1), XYZnames{nub}=deblank(PU{nid1});
                else XYZnames{nub}=['undefined-',num2str(ub(nub))];end
            end
            b(idxvoxels)=iub;
        else error(['File ',filelabels,' format not recognized']);
        end
        
        Labels={{'Left-Cerebral-Cortex','Right-Cerebral-Cortex'},...
            {'Left-Cerebral-White-Matter','Right-Cerebral-White-Matter'},...
            {'Left-Lateral-Ventricle','Left-Inf-Lat-Vent','3rd-Ventricle','4th-Ventricle','CSF','Right-Lateral-Ventricle','Right-Inf-Lat-Vent','5th-Ventricle'}};
        for nreg=1:3,
            [ok,idx]=ismember(Labels{nreg},XYZnames);
            idx=idx(ok);
            b_reg{nreg}=ismember(b,idx);
            
            V=struct('mat',a.mat,'dim',a.dim,'dt',[spm_type('uint8') spm_platform('bigend')],'fname',filenames{nreg});
            spm_write_vol(V,b_reg{nreg});
            conn_disp(['Saved file ',V.fname,' (',num2str(nnz(b_reg{nreg})),' voxels)']);
        end
    end
end
