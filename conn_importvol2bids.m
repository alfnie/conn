function [out,changed,nV]=conn_importvol2bids(filename,nsub,nses,fmod,ftype,docopy,dodisp,dosinglesubjectreduce,roinumber,keepnames)
% CONN_IMPORTVOL2BIDS imports functional/anatomical file into CONN/BIDS directory
% conn_importvol2bids(filename,nsub,nses,fmod,[froot])
%   fmod   : file modality (anat, func, dwi, fmap, roi)
%   froot  : filename root (defaults: T1w for 'anat'; task-rest_bold for 'func')

% note: for secondary datasets fmod uses the format "dataset(\d+)(.*)" where $1=nset and $2=fmod (e.g. dataset1func)

global CONN_x;
SOFTLINK=false; % change to "true" to use soft links instead of copying files (note: this feature is not yet fully tested)

if nargin<=1
    if nargin==1, dosinglesubjectreduce=filename;
    else dosinglesubjectreduce=false;
    end
    changed=false;
    out={};
    V=[];
    for nsub=1:CONN_x.Setup.nsubjects
        nsess=CONN_x.Setup.nsessions(min(length(CONN_x.Setup.nsessions),nsub));
        if ~CONN_x.Setup.structural_sessionspecific, nsesstemp=1; else nsesstemp=nsess; end
        for nses=1:nsesstemp,
            if ~isempty(CONN_x.Setup.structural{nsub}{nses}{1})
                [out{1}{nsub}{nses},tchanged]=conn_importvol2bids(CONN_x.Setup.structural{nsub}{nses}{1},nsub,nses,'anat',[],[],true,dosinglesubjectreduce);
                changed=changed|tchanged;
            end
        end
        for nses=nsesstemp+1:nsess,
            CONN_x.Setup.structural{nsub}{nses}=CONN_x.Setup.structural{nsub}{nsesstemp};
            out{1}{nsub}{nses}=CONN_x.Setup.structural{nsub}{nses}{1};
        end
        for nses=1:nsess,
            if ~isempty(CONN_x.Setup.functional{nsub}{nses}{1})
                [out{2}{nsub}{nses},tchanged]=conn_importvol2bids(CONN_x.Setup.functional{nsub}{nses}{1},nsub,nses,'func',[],[],true,dosinglesubjectreduce);
                changed=changed|tchanged;
            end
        end
        for nses=1:nsess,
            if numel(CONN_x.Setup.unwarp_functional)>=nsub&&numel(CONN_x.Setup.unwarp_functional{nsub})>=nses
                if ~isempty(CONN_x.Setup.unwarp_functional{nsub}{nses}{1})
                    [out{3}{nsub}{nses},tchanged]=conn_importvol2bids(CONN_x.Setup.unwarp_functional{nsub}{nses}{1},nsub,nses,'vdm','vdm',[],true,dosinglesubjectreduce);
                    changed=changed|tchanged;
                end
            end
        end
        if ~isempty(conn_datasetlabel('vdm'))
            for nses=1:nsess,
                fname=conn_get_functional(nsub,nses,'vdm');
                if ~isempty(fname)
                    [out{4}{nsub}{nses},tchanged]=conn_importvol2bids(fname,nsub,nses,'vdm','vdm',[],true,dosinglesubjectreduce);
                    changed=changed|tchanged;
                end
            end
        end
        if ~isempty(conn_datasetlabel('fmap'))
            for nses=1:nsess,
                fname=conn_get_functional(nsub,nses,'fmap');
                if ~isempty(fname)
                    [out{5}{nsub}{nses},tchanged]=conn_importvol2bids(fname,nsub,nses,'fmap','fmap',[],true,dosinglesubjectreduce);
                    changed=changed|tchanged;
                end
            end
        end
        for nroi=1:numel(CONN_x.Setup.rois.names)-1
            if CONN_x.Setup.rois.subjectspecific(nroi) %||nsub==1 
                if ~CONN_x.Setup.rois.sessionspecific(nroi), nsesstemp=1; else nsesstemp=nsess; end
                for nses=1:nsesstemp,
                    if ~isempty(CONN_x.Setup.rois.files{nsub}{nroi}{nses}{1})
                        if CONN_x.Setup.rois.sessionspecific(nroi), saveasnses=nses;
                        else saveasnses=[];
                        end                        
                        [out{6}{nsub}{nses}{nroi},tchanged]=conn_importvol2bids(CONN_x.Setup.rois.files{nsub}{nroi}{nses}{1},nsub,saveasnses,'roi',[],[],true,dosinglesubjectreduce,nroi);
                        changed=changed|tchanged;
                    end
                end
                for nses=nsesstemp+1:nsess,
                    CONN_x.Setup.rois.files{nsub}{nroi}{nses}=CONN_x.Setup.rois.files{nsub}{nroi}{nsesstemp};
                    out{6}{nsub}{nses}{nroi}=CONN_x.Setup.rois.files{nsub}{nroi}{nses}{1};
                end
            %else % note: imports only one sample (subject #1) for subject-independent rois?
            %    for nses=1:nsess,
            %        CONN_x.Setup.rois.files{nsub}{nroi}{nses}=CONN_x.Setup.rois.files{1}{nroi}{nses};
            %        out{6}{nsub}{nses}{nroi}=CONN_x.Setup.rois.files{nsub}{nroi}{nses}{1};
            %    end
            end
        end
    end
    return
end

SESSASRUNS=true; % change to false for replicating older-versions structure
changed=false;
nV=[];
if nargin<10||isempty(keepnames), keepnames=false; end
if nargin<9||isempty(roinumber), roinumber=[]; end
if nargin<8||isempty(dosinglesubjectreduce), dosinglesubjectreduce=false; end
if nargin<7||isempty(dodisp), dodisp=false; end
if nargin<6||isempty(docopy), docopy=true; end
if nargin<5||isempty(ftype), ftype=[]; end
nset=0; 
if keepnames
    if ~iscell(filename), filename=cellstr(filename); end
    idx=regexp(filename{1},'[\\\/]derivatives[\\\/]([^\\\/]+[\\\/])?sub-');
    if ~isempty(idx), % derivatives/*/sub-   to derivatives/conn/sub-
        newfilepath=fullfile(filename{1}(1:idx(end)-1),'derivatives','conn');
        idx=regexp(filename{1},'[\\\/]derivatives[\\\/]([^\\\/]+[\\\/])?sub-','end');
        newfilename=filename{1}(idx(end)-3:end);
    else % sub- to derivatives/conn/sub-
        idx=regexp(filename{1},'[\\\/]sub-[^\\\/]+[\\\/]');
        if ~isempty(idx),
            newfilepath=fullfile(filename{1}(1:idx(end)-1),'derivatives','conn');
            newfilename=filename{1}(idx(end)+1:end);
        else
            [fpath,fname,fext]=fileparts(filename{1});
            if isempty(fpath), fpath=pwd; end
            newfilepath=fpath;
            newfilename=[fname,fext];
        end
    end
    newfilename=regexprep(newfilename,'\.gz$','');
    out=fullfile(newfilepath,newfilename);
    [fpath,nill,nill]=fileparts(out);
    conn_fileutils('mkdir',fpath);
else
    if dosinglesubjectreduce, BIDSfolder=fileparts(fileparts(CONN_x.folders.bids));
    elseif ismember(fmod,{'func','anat','dwi','fmap','roi','rois'}), BIDSfolder=fullfile(CONN_x.folders.bids,'dataset');  % root directory
    elseif ~isempty(regexp(fmod,'^dataset\d*')), nset=regexp(fmod,'^dataset(\d*)','tokens','once'); nset=str2double([nset{:}]); if isnan(nset), nset=0; end; 
        BIDSfolder=fullfile(CONN_x.folders.bids,regexp(fmod,'^dataset\d*','match','once')); 
        if isequal(CONN_x.Setup.secondarydataset(nset).label,'fmap'), fmod='fmap'; 
        else fmod=regexprep(fmod,'^dataset\d*','');
        end
    else BIDSfolder=fullfile(CONN_x.folders.bids,'derivatives'); % derivatives
    end
    if SESSASRUNS&&~isempty(nses), frun=sprintf('run-%02d_',nses(1));
    else frun='';
    end
    if isequal(fmod,'rois'), fmod='roi'; end % backcompat fix
    if isempty(ftype)
        switch(fmod)
            case 'func', ftype='bold'; if isempty(frun), newfilename='task-rest_bold.nii'; else newfilename=sprintf('%sbold.nii',frun); end
            case 'anat', ftype='T1w'; newfilename=sprintf('%sT1w.nii',frun);
            case 'dwi',  ftype='dwi'; newfilename=sprintf('%sdwi.nii',frun);
            case 'vdm', ftype='vdm'; newfilename=sprintf('%svdm.nii',frun);
            case 'fmap', ftype='fmap'; newfilename=sprintf('%sfmap.nii',frun);
            case 'roi', ftype='roi';
                if ~isempty(roinumber)&& numel(CONN_x.Setup.rois.names)>=roinumber, newfilename=sprintf('%sroi-%s.nii',frun,regexprep(CONN_x.Setup.rois.names{roinumber},'[^\w\d_]',''));
                else newfilename=sprintf('%sroi.nii',frun);
                end
            otherwise,   ftype='unknown'; newfilename=sprintf('%sunknown.nii',frun);
        end
    else newfilename=sprintf('%s%s.nii',frun,ftype);
    end
    if dosinglesubjectreduce>1, subjectid=fmod; if isempty(subjectid), subjectid='data'; end; fsub=subjectid;
    elseif dosinglesubjectreduce, [nill,subjectid,nill]=fileparts(CONN_x.filename); fsub=subjectid;
    elseif ~isempty(nsub),    fsub=sprintf('sub-%04d',nsub(1));
    else fsub='';
    end
    if ~SESSASRUNS&&~isempty(nses), fses=sprintf('ses-%02d',nses(1)); else fses=''; end
    newfilepath=BIDSfolder; conn_fileutils('mkdir',newfilepath);
    if ~dosinglesubjectreduce
        if ~isempty(fsub),      conn_fileutils('mkdir',newfilepath,fsub); newfilepath=fullfile(newfilepath,fsub); end
        if ~isempty(fses),      conn_fileutils('mkdir',newfilepath,fses); newfilepath=fullfile(newfilepath,fses);  end
    end
    if ~isempty(fmod),      conn_fileutils('mkdir',newfilepath,fmod); newfilepath=fullfile(newfilepath,fmod); end
    if ~isempty(fses),      newfilename=sprintf('%s_%s',fses,newfilename); end
    if ~isempty(fsub),      newfilename=sprintf('%s_%s',fsub,newfilename); end
    out=fullfile(newfilepath,newfilename);
end

if iscell(filename), filename=char(filename); end
[nill,nill,fext]=fileparts(out);
[tfileroot,tfileext1,tfileext2,tfilenum]=conn_fileparts(filename(1,:));
exts={};
if strcmp(fext,'.nii')&&((size(filename,1)>1) || (strcmp(tfileext1,'.img')&&isempty(tfileext2))), % convert .img to .nii, or multiple .nii to single 4d .nii
    if docopy
        f=conn_dir(conn_prepend('',out,'.*'),'-R');
        if ~isempty(f)
            f=cellstr(f);
            conn_fileutils('spm_unlink',f{:});
        end
        changed=true;
        if size(filename,1)>1 % if merging multiple compressed files, decompress first to temporal files in target folder
            filename=cellstr(filename);
            for n1=1:numel(filename)
                [nill,tfileext1,tfileext2,tfilenum]=conn_fileparts(filename{n1});
                if isequal(tfileext2,'.gz') 
                    tfilename=conn_prepend('',out,['.part',num2str(n1),tfileext1,tfileext2]);
                    conn_fileutils('copyfile',filename{n1},tfilename);
                    filename{n1}=conn_gz2nii(tfilename);
                    conn_fileutils('spm_unlink',tfilename);
                end
            end
            filename=char(filename);
        end
        a=conn_fileutils('spm_vol',filename);
        conn_fileutils('spm_file_merge',a,out);
        filename=out;
    end
else % copy/link file
    out=conn_prepend('',fullfile(newfilepath,newfilename),[tfileext1,tfileext2,tfilenum]); % keep extension of input file (note: possible .gz included)
    if docopy && ~isequal(filename, out)
        f=conn_dir(conn_prepend('',out,'.*'),'-R');
        if ~isempty(f)
            f=cellstr(f);
            conn_fileutils('spm_unlink',f{:});
        end
        changed=true;
        exts={[tfileext1,tfileext2]};
        if strcmp(tfileext1,'.img'), exts=[exts, {'.hdr'}]; end
    end
end
if docopy&&changed
    switch(fmod) % possible associated files
        case 'dwi',  exts=[exts {'.bval','.bvec','.json'}];
        case 'roi', exts=[exts {'.txt','.csv','.xls','.info','.icon.jpg','.json'}];
        otherwise,   exts=[exts {'.json'}];
    end
    for nexts=1:numel(exts) % copy/link original and additional files if needed
        tfilename=[tfileroot,exts{nexts}];
        if conn_existfile(tfilename)
            outfile=conn_prepend('',fullfile(newfilepath,newfilename),exts{nexts});
            if SOFTLINK&~ispc, [ok,nill]=system(['ln -fs ''',tfilename,''' ''',outfile,'''']);
            else try, conn_fileutils('copyfile',tfilename,outfile); end
            end
            %if ispc, [ok,nill]=system(['copy "',tfilename,'" "',outfile,'"']);
            %elseif SOFTLINK, [ok,nill]=system(['ln -fs ''',tfilename,''' ''',outfile,'''']);
            %else, [ok,nill]=system(['cp ''',tfilename,''' ''',outfile,'''']);
            %end
        end
    end
    if strcmp(fmod,'func')
        tfilename=conn_prepend('vdm',filename); % copy/link original and vdm files if needed
        if ~conn_existfile(tfilename), tfilename=conn_prepend('vdm_',filename); end
        if ~conn_existfile(tfilename), tfilename=conn_prepend('vdm5_',filename); end
        if ~conn_existfile(tfilename)
            tfilename=conn_dirn(fullfile(fileparts(filename),'vdm*.nii'));
            if numel(tfilename)==1, tfilename=fullfile(fileparts(filename),tfilename(1).name);
            else tfilename='';
            end
        end
        if ~conn_existfile(tfilename)
            tfilename=conn_dirn(fullfile(fileparts(filename),'vdm*.img'));
            if numel(tfilename)==1, tfilename=fullfile(fileparts(filename),tfilename(1).name);
            else tfilename='';
            end
        end
        if conn_existfile(tfilename)
            if ~isempty(regexp(tfilename,'\.img$'))
                outfile=conn_prepend('vdm_',fullfile(newfilepath,newfilename),'.img');
                if SOFTLINK&~ispc, [ok,nill]=system(['ln -fs ''',tfilename,''' ''',outfile,'''']);
                else try, conn_fileutils('copyfile',tfilename,outfile); end
                end
%                 if ispc, [ok,nill]=system(['copy "',tfilename,'" "',outfile,'"']);
%                 elseif SOFTLINK, [ok,nill]=system(['ln -fs ''',tfilename,''' ''',outfile,'''']);
%                 else, [ok,nill]=system(['cp ''',tfilename,''' ''',outfile,'''']);
%                 end
                outfile=conn_prepend('vdm_',fullfile(newfilepath,newfilename),'.hdr');
                tfilename=conn_prepend('',tfilename,'.hdr');
                if SOFTLINK&~ispc, [ok,nill]=system(['ln -fs ''',tfilename,''' ''',outfile,'''']);
                else try, conn_fileutils('copyfile',tfilename,outfile); end
                end
%                 if ispc, [ok,nill]=system(['copy "',tfilename,'" "',outfile,'"']);
%                 elseif SOFTLINK, [ok,nill]=system(['ln -fs ''',tfilename,''' ''',outfile,'''']);
%                 else, [ok,nill]=system(['cp ''',tfilename,''' ''',outfile,'''']);
%                 end
            else
                outfile=conn_prepend('vdm_',fullfile(newfilepath,newfilename));
                if SOFTLINK&~ispc, [ok,nill]=system(['ln -fs ''',tfilename,''' ''',outfile,'''']);
                else try, conn_fileutils('copyfile',tfilename,outfile); end
                end
%                 if ispc, [ok,nill]=system(['copy "',tfilename,'" "',outfile,'"']);
%                 elseif SOFTLINK, [ok,nill]=system(['ln -fs ''',tfilename,''' ''',outfile,'''']);
%                 else, [ok,nill]=system(['cp ''',tfilename,''' ''',outfile,'''']);
%                 end
            end
        end
    end
    if strcmp(fmod,'anat')
        [ok,nill,nill,fsfiles]=conn_checkFSfiles(filename); % copy/link original and additional FS files if needed
        if ok
            for nexts=1:numel(fsfiles) 
                tfilename=fsfiles{nexts};
                [tfilename_path,tfilename_name,tfilename_ext]=spm_fileparts(tfilename);
                [nill,tfilename_path1,tfilename_path2]=fileparts(tfilename_path);
                tfilename_path=[tfilename_path1 tfilename_path2];
                conn_fileutils('mkdir',fileparts(newfilepath),tfilename_path);
                outfile=fullfile(fileparts(newfilepath),tfilename_path,[tfilename_name tfilename_ext]);
                if SOFTLINK&~ispc, [ok,nill]=system(['ln -fs ''',tfilename,''' ''',outfile,'''']);
                else try, conn_fileutils('copyfile',tfilename,outfile); end
                end
%                 if ispc, [ok,nill]=system(['copy "',tfilename,'" "',outfile,'"']);
%                 elseif SOFTLINK, [ok,nill]=system(['ln -fs ''',tfilename,''' ''',outfile,'''']);
%                 else, [ok,nill]=system(['cp ''',tfilename,''' ''',outfile,'''']);
%                 end
            end
            if conn_importaseg(fileparts(filename),[],true), % import gray/white/csf masks from aseg.mgz
                try,
                    tfilename=conn_importaseg(fileparts(filename));
                    for nseg=1:3
                        CONN_x.Setup.rois.files{nsub(end)}{nseg}{nses(end)}=conn_file(tfilename{nseg});
                    end
                end
            end
        end
    end
end
if docopy
    if isempty(nses), nses=1; end
    if strcmp(fmod,'func')&&~conn_existfile(conn_prepend('',out,'.json')),
        try % if none exists, initializes .json with TR info
            conn_fileutils('spm_jsonwrite',conn_prepend('',out,'.json'),struct('RepetitionTime',conn_get_rt(nsub(end),nses(end))));
        end
    end
    if changed&&dodisp, fprintf('Created file %s from %s\n',out,filename(1,:)); end
    switch(fmod) % initializes CONN_x functional/structural fields
        case 'func'
            nV=conn_set_functional(nsub(end),nses(end),nset,out); % note: when nsubs = [nsubs1 nsubs2] or nses = [nses1 nses2] file_bids(nsub1,nses1) is assigned to nsub2,nses2 (in order to allow subject- and session- independent data)
%             if ~nset, 
%                 [CONN_x.Setup.functional{nsub(end)}{nses(end)},nV]=conn_file(out); 
%             else 
%                 CONN_x.Setup.secondarydataset(nset).functionals_type=4; 
%                 [CONN_x.Setup.secondarydataset(nset).functionals_explicit{nsub(end)}{nses(end)},nV]=conn_file(out);
%             end
        case 'anat'
            [CONN_x.Setup.structural{nsub(end)}{nses(end)},nV]=conn_file(out);
        case 'roi',
            if ~isempty(roinumber), CONN_x.Setup.rois.files{nsub(end)}{roinumber}{nses(end)}=conn_file(out); end
        case 'vdm'
            if isequal(ftype,'vdm'), [CONN_x.Setup.unwarp_functional{nsub(end)}{nses(end)},nV]=conn_file(out); end
            conn_set_functional(nsub(end),nses(end),ftype,out);
        otherwise
            conn_set_functional(nsub(end),nses(end),ftype,out);
    end
end
end

function varargout=conn_fileparts(filename)
varargout=regexp(deblank(filename(1,:)),'^(.*?)(\.[^\.]*?)(\.gz)?(\,\d+)?$','tokens','once');
if numel(varargout)<nargout, varargout=[varargout, repmat({''},1,nargout-numel(varargout))]; end
end
