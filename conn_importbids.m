function [ok,ERR,WRN]=conn_importbids(filenames,varargin);
global CONN_x;

options=struct('type','functional',...
               'localcopy',false,...
               'copytoderiv',false,...
               'nset',0,...
               'bidsname','',...
               'sessions',0,...
               'subjects_id',[],...
               'subjects',[]);
for n1=1:2:nargin-2, if ~isfield(options,lower(varargin{n1})), error('unknown option %s',lower(varargin{n1})); else options.(lower(varargin{n1}))=varargin{n1+1}; end; end
if ~isempty(options.subjects), nsubs=options.subjects; else nsubs=1:CONN_x.Setup.nsubjects; end
if ~isempty(options.bidsname), bidsname=options.bidsname;
elseif ~options.nset,          bidsname='func';
else                           bidsname=sprintf('dataset%dfunc',options.nset);
end
if isequal(filenames,'all')&&isequal(options.type,'conditions'), filenames=conn_module('get','functionals'); end % note: syntax "conn_importbids all type conditions"

ok=false;
ERR={};
WRN={};
file_name=[];
for isub=1:numel(nsubs),
    if isempty(options.subjects_id), filename=filenames{isub};
    else
        if isempty(file_name), [nill,file_name,nill]=cellfun(@fileparts,filenames,'uni',0); end
        id=options.subjects_id{isub};
        select=cellfun('length',regexp(file_name,['^sub-',id,'$|^sub-',id,'_']))>0;
        filename=filenames(select);
        if isempty(filename), conn_disp('fprintf','warning: no match found for subject sub-%s\n',id); end
    end
    nsub=nsubs(isub);
    if isempty(options.sessions), nses=[];
    elseif iscell(options.sessions), nses=options.sessions{isub};
    else nses=options.sessions;
    end
    if isequal(nses,0), % match number of sessions to input data
        if isequal(options.type,'functional')
            if numel(CONN_x.Setup.nsessions)==1&&CONN_x.Setup.nsubjects>1, CONN_x.Setup.nsessions=CONN_x.Setup.nsessions+zeros(1,CONN_x.Setup.nsubjects); end
            CONN_x.Setup.nsessions(nsub)=numel(filename);
        end
        nses=[];
    end
    nsess=CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsub));
    if isempty(nses), nsesstemp=1:nsess;
    else nsesstemp=nses;
    end
    switch(options.type)
        case 'functional'
            if numel(filename)~=numel(nsesstemp), ERR{end+1}=sprintf('Subject %d functional import skipped : mismatched number of selected files (%d) and number of target sessions (%d)',nsub,numel(filename),numel(nsesstemp));
            else
                if ~options.nset, 
                    if numel(CONN_x.Setup.RT)<nsub, CONN_x.Setup.RT=CONN_x.Setup.RT(min(numel(CONN_x.Setup.RT), 1:nsub)); end
                    CONN_x.Setup.RT(nsub)=nan;
                end
                for n2=1:numel(nsesstemp)
                    nses=nsesstemp(n2);
                    if options.copytoderiv
                        [out,nill,nV]=conn_importvol2bids(filename{n2},nsub,nses,'func',[],[],[],[],[],true);
                    elseif options.localcopy,
                        [out,nill,nV]=conn_importvol2bids(filename{n2},nsub,nses,bidsname);
                        %if ~options.nset, CONN_x.Setup.nscans{nsub}{nses}=nV; end
                    else
                        out=filename{n2};
                        nV=conn_set_functional(nsub,nses,options.nset,out);
                    end
                    conn_disp('fprintf','functional %s imported to subject %d session %d\n',filename{n2},nsub,nses);
                    if ~options.nset, % create task-* condition names
                    end
                    if ~options.nset, % search for fmriprep confound regressors file 
                        fname=conn_prepend('',regexprep(filename{n2},'(_space-[^\._\\\/]*)?(_res-[^\._\\\/]*)?_desc-[^\.\\\/]*\.nii(\.gz)?$',''),'_desc-confounds_timeseries.tsv');
                        if ~conn_existfile(fname), fname=conn_prepend('',regexprep(filename{n2},'(_mod-[^\._\\\/]*)?(_echo-[^\._\\\/]*)?(_recording-[^\._\\\/]*)?(_proc-[^\._\\\/]*)?(_space-[^\._\\\/]*)?(_res-[^\._\\\/]*)?_desc-[^\.\\\/]*\.nii(\.gz)?$',''),'_desc-confounds_timeseries.tsv'); end %
                        if ~conn_existfile(fname), fname=conn_prepend('',regexprep(filename{n2},'(_space-[^\._\\\/]*)?(_res-[^\._\\\/]*)?_desc-[^\.\\\/]*\.nii(\.gz)?$',''),'_desc-confounds_regressors.tsv'); end % back-compatibility
                        if conn_existfile(fname)
                            data=conn_loadtextfile(fname);
                            if isstruct(data)
                                conn_disp('fprintf','regressors %s imported to subject %d session %d\n',fname,nsub,nses);
                                names=reshape(fieldnames(data),1,[]);
                                enames=regexprep(names,{'.*','(_?\d+|_[xyz])$'},{'QC_$0',''});
                                enames=[enames,regexprep(enames,{'^QC_(trans|rot)$','^QC_motion_outlier$','^QC_(std_dvars|framewise_displacement)$'},{'realignment','scrubbing','QC_timeseries'})];
                                [unames,nill,uidx]=unique(enames);
                                [nill,tidx]=sort(cellfun('length',regexp(unames,'^QC_')));
                                for tn1=reshape(tidx,1,[])
                                    tname=unames{tn1};
                                    R=[]; for tn2=reshape(unique(1+rem(find(uidx==tn1)-1,numel(names))),1,[]), R=cat(2,R,data.(names{tn2})); end; 
                                    try, nanfirst=isnan(R(1,:))&all(~isnan(R(2:end,:)),1); if any(nanfirst), R(1,nanfirst)=0; end; end % note: fix issue with first timepoint in dvars and framewise_displacement variables having a NaN value (converts to 0's)
                                    idx=strmatch(tname,CONN_x.Setup.l1covariates.names,'exact');
                                    if isempty(idx), idx=length(CONN_x.Setup.l1covariates.names); CONN_x.Setup.l1covariates.names{end+1}=' '; end                                    
                                    try, if size(R,2)>1, out2=conn_prepend('',out,sprintf('.%s.mat',tname)); conn_savematfile(out2,'R'); R=out2; end; end % note: saves multivariate timeseries to files (to keep conn_*.mat filesize small)
                                    CONN_x.Setup.l1covariates.names{idx}=tname;
                                    CONN_x.Setup.l1covariates.files{nsub}{idx}{nses}=conn_file(R);
                                end
                            else
                                conn_disp('fprintf','warning: unexpected format of file %s (subject %d session %d). Skipping import\n',fname,nsub,nses);
%                             else
%                                 name='QC_fmriprep';
%                                 idx=strmatch(name,CONN_x.Setup.l1covariates.names,'exact');
%                                 if isempty(idx), idx=length(CONN_x.Setup.l1covariates.names); CONN_x.Setup.l1covariates.names{end+1}=' '; end
%                                 CONN_x.Setup.l1covariates.names{idx}=name;
%                                 CONN_x.Setup.l1covariates.files{nsub}{idx}{nses}=conn_file(fname);
                            end
                        end
                    end
                    if ~options.nset, % search for fmap/*_magnitude.nii and fmap/*_phasediff.nii or fmap/*_real1.nii fmap/*_real2.nii fmap/*_imag1.nii fmap/*_imag2.nii fmap
                        enames={'phasediff','magnitude','magnitude1','magnitude2','real1','real2','imag1','imag2','fieldmap','fmap'};
                        fnames={};
                        for ne=1:numel(enames)
                            fname=regexprep(filename{n2},'func(\/|\\)([^\\\/]*)_bold\.nii(\.gz)?$',['fmap$1$2_',enames{ne},'.nii']);
                            if ~conn_existfile(fname), fname=regexprep(filename{n2},'func(\/|\\)([^\\\/]*)_bold\.nii(\.gz)?$',['fmap$1$2_',enames{ne},'.nii.gz']); end
                            if ~conn_existfile(fname), tname=regexprep(filename{n2},'func(\/|\\)([^\\\/]*)_bold\.nii(\.gz)?$',['fmap$1*_',enames{ne},'.nii']); tdir=conn_dirn(tname); if numel(tdir)==1, fname=fullfile(fileparts(tname),tdir.name); end; end
                            if ~conn_existfile(fname), tname=regexprep(filename{n2},'func(\/|\\)([^\\\/]*)_bold\.nii(\.gz)?$',['fmap$1*_',enames{ne},'.nii.gz']); tdir=conn_dirn(tname); if numel(tdir)==1, fname=fullfile(fileparts(tname),tdir.name); end; end
                            if conn_existfile(fname), fnames{ne}=fname; else fnames{ne}=''; end
                        end
                        v=cellfun('length',fnames);
                        valid={[2,1],[3,4,1],[3,1],[5,6,7,8],9,10};
                        for nv=1:numel(valid), if all(v(valid{nv})), fnames=fnames(valid{nv}); break; end; end
                        if all(cellfun('length',fnames))&&~options.nset
                            if options.copytoderiv
                                [nill,nill,nV]=conn_importvol2bids(char(fnames),nsub,nses,'fmap','fmap',[],[],[],[],true);
                            elseif options.localcopy,
                                [nill,nill,nV]=conn_importvol2bids(char(fnames),nsub,nses,'fmap','fmap');
                            else
                                nV=conn_set_functional(nsub,nses,'fmap',char(fnames));
                            end
                            conn_disp('fprintf','fmap %s imported to subject %d session %d\n',nsub,nses);
                        end
                    end
                end
            end
        case 'structural'
            if isempty(filename), ok=true; %return; end
            else
                if isequal(nsesstemp,1:nsess), CONN_x.Setup.structural_sessionspecific=numel(filename)>1&numel(filename)==nsess; end
                if CONN_x.Setup.structural_sessionspecific
                    if numel(filename)~=numel(nsesstemp), ERR{end+1}=sprintf('Subject %d structural import skipped : mismatched number of selected files (%d) and number of target sessions (%d)',nsub,numel(filename),numel(nsesstemp));
                    else
                        for n2=1:numel(nsesstemp)
                            nses=nsesstemp(n2);
                            if options.copytoderiv, conn_importvol2bids(filename{n2},nsub,nses,'anat',[],[],[],[],[],true);
                            elseif options.localcopy, conn_importvol2bids(filename{n2},nsub,nses,'anat');
                            else CONN_x.Setup.structural{nsub}{nses}=conn_file(filename{n2});
                            end
                        end
                    end
                else
                    if numel(filename)~=1,
                        conn_disp('fprintf','warning: subject %d listed multiple strutural files (%d), importing first file only\n',nsub,numel(filename));
                        WRN{end+1}=sprintf('warning: subject %d listed multiple strutural files (%d), importing first file only\n',nsub,numel(filename));
                    end
                    for n2=1:numel(nsesstemp)
                        nses=nsesstemp(n2);
                        if options.copytoderiv, conn_importvol2bids(filename{1},nsub,[1,nses],'anat',[],[],[],[],[],true);
                        elseif options.localcopy, conn_importvol2bids(filename{1},nsub,[1,nses],'anat');
                        else CONN_x.Setup.structural{nsub}{nses}=conn_file(filename{1});
                        end
                        conn_disp('fprintf','structural %s imported to subject %d session %d\n',filename{1},nsub,nses);
                    end
                end
                for n2=1:numel(nsesstemp) % search for fmriprep gray/white/csf files
                    nses=nsesstemp(n2);
                    if CONN_x.Setup.structural_sessionspecific, n3=n2; else n3=1; end
                    f=reshape({regexprep(filename{n3},'desc-preproc_T1w\.nii(\.gz)?$','label-GM_probseg.nii'),regexprep(filename{n3},'desc-preproc_T1w\.nii(\.gz)?$','label-GM_probseg.nii.gz'),...
                        regexprep(filename{n3},'desc-preproc_T1w\.nii(\.gz)?$','label-WM_probseg.nii'),regexprep(filename{n3},'desc-preproc_T1w\.nii(\.gz)?$','label-WM_probseg.nii.gz'),...
                        regexprep(filename{n3},'desc-preproc_T1w\.nii(\.gz)?$','label-CSF_probseg.nii'),regexprep(filename{n3},'desc-preproc_T1w\.nii(\.gz)?$','label-CSF_probseg.nii.gz')},2,[]);
                    ef=reshape(conn_existfile(f),size(f))&cellfun(@(x)~isequal(filename{n3},x),f);
                    for nmask=1:3,
                        if any(ef(:,nmask)),
                            localcopy_reduce=false;
                            temp2=f{find(ef(:,nmask),1),nmask};
                            if options.copytoderiv,
                                if CONN_x.Setup.structural_sessionspecific, conn_importvol2bids(temp2,nsub,nses,'roi',[],[],[],localcopy_reduce,nmask,true);
                                elseif nses==1, conn_importvol2bids(temp2,nsub,[],'roi',[],[],[],localcopy_reduce,nmask,true);
                                else CONN_x.Setup.rois.files{nsub}{nmask}{nses}=CONN_x.Setup.rois.files{nsub}{nmask}{1};
                                end
                            elseif options.localcopy,
                                if CONN_x.Setup.structural_sessionspecific, conn_importvol2bids(temp2,nsub,nses,'roi',[],[],[],localcopy_reduce,nmask);
                                elseif nses==1, conn_importvol2bids(temp2,nsub,[],'roi',[],[],[],localcopy_reduce,nmask);
                                else CONN_x.Setup.rois.files{nsub}{nmask}{nses}=CONN_x.Setup.rois.files{nsub}{nmask}{1};
                                end
                            else
                                CONN_x.Setup.rois.files{nsub}{nmask}{nses}=conn_file(temp2);
                            end
                            conn_disp('fprintf','ROI%d %s imported to subject %d session %d\n',nmask,temp2,nsub,nses);
                        end
                    end
                end
            end
        case 'conditions'
            if isempty(filename), ok=true; %return; end
            else
                if numel(filename)~=numel(nsesstemp), ERR{end+1}=sprintf('Subject %d conditions import skipped : mismatched number of selected files (%d) and number of target sessions (%d)',nsub,numel(filename),numel(nsesstemp));
                else
                    for n2=1:numel(nsesstemp)
                        nses=nsesstemp(n2);
                        tname=filename{n2};
                        fname=conn_prepend('',regexprep(tname,'_bold\.nii(\.gz)?$','.nii'),'_events.tsv');
                        if ~conn_existfile(fname), fname=conn_prepend('',regexprep(tname,'(_space-[^\._\\\/]*)?(_res-[^\._\\\/]*)?(_desc-[^\._\\\/]*)_bold\.nii(\.gz)?$','$3.nii'),'_events.tsv'); end
                        if ~conn_existfile(fname), fname=conn_prepend('',regexprep(tname,'(_space-[^\._\\\/]*)?(_res-[^\._\\\/]*)?(_desc-[^\._\\\/]*)_bold\.nii(\.gz)?$','.nii'),'_events.tsv'); end
                        if ~conn_existfile(fname)&&~isempty(regexp(filename{n2},'[\\\/]derivatives[\\\/]fmriprep([\\\/]sub-)')) % search condition info in raw-data directory (for fmriprep import)
                            tname=regexprep(filename{n2},'[\\\/]derivatives[\\\/]fmriprep([\\\/]sub-)','$1');
                            fname=conn_prepend('',regexprep(tname,'_bold\.nii(\.gz)?$','.nii'),'_events.tsv');
                            if ~conn_existfile(fname), fname=conn_prepend('',regexprep(tname,'(_space-[^\._\\\/]*)?(_res-[^\._\\\/]*)?(_desc-[^\._\\\/]*)_bold\.nii(\.gz)?$','$3.nii'),'_events.tsv'); end
                            if ~conn_existfile(fname), fname=conn_prepend('',regexprep(tname,'(_space-[^\._\\\/]*)?(_res-[^\._\\\/]*)?(_desc-[^\._\\\/]*)_bold\.nii(\.gz)?$','.nii'),'_events.tsv'); end
                        end
                        if conn_existfile(fname),
                            conn_importcondition({fname},'subjects',nsub,'sessions',nses,'breakconditionsbysession',false,'deleteall',false);
                            conn_disp('fprintf','conditions in %s imported\n',fname);
                            %                     elseif ~isempty(regexp(filename{n2},'_task-rest_[^\\\/]*$'))
                            %                         cname=char(regexp(filename{n2},'(_ses-[^\._\\\/]*)?_task-rest_[^\\\/]*$','tokens','once'));
                            %                         if isempty(cname), cname='rest';
                            %                         else cname=[regexprep(cname,'^_ses-',''),'_','rest'];
                            %                         end
                            %                         conn_importcondition(struct('conditions',{{cname}},'onsets',0,'durations',inf),'subjects',nsub,'sessions',nses,'breakconditionsbysession',false,'deleteall',false);
                            %                         conn_disp('fprintf','condition %s imported as %s\n',filename{n2},cname);
                        elseif ~isempty(regexp(filename{n2},'_task-([^_\\\/]+)_[^\\\/]*$'))
                            cname=regexp(filename{n2},'(_ses-[^_\\\/]*)?_task-([^_\\\/]+)_[^\\\/]*$','tokens','once');
                            if numel(cname)==2&&~isempty(cname{1}), cname=[regexprep(cname{1},'^_ses-',''),'_',cname{2}];
                            elseif numel(cname)==2&&isempty(cname{1}), cname=cname{2};
                            else cname=char(cname);
                            end
                            conn_importcondition(struct('conditions',{{cname}},'onsets',0,'durations',inf),'subjects',nsub,'sessions',nses,'breakconditionsbysession',false,'deleteall',false);
                            conn_disp('fprintf','warning: condition file %s not found. Imported as %s\n',fname,cname);
                            WRN{end+1}=sprintf('warning: condition file %s not found. Imported as %s\n',fname,cname);
                        elseif ~isempty(regexp(filename{n2},'_ses-([^_\\\/]+)_[^\\\/]*$'))
                            cname=char(regexp(filename{n2},'_ses-([^_\\\/]+)_[^\\\/]*$','tokens','once'));
                            conn_importcondition(struct('conditions',{{cname}},'onsets',0,'durations',inf),'subjects',nsub,'sessions',nses,'breakconditionsbysession',false,'deleteall',false);
                            conn_disp('fprintf','warning: condition file %s not found. Imported as %s\n',fname,cname);
                            WRN{end+1}=sprintf('warning: condition file %s not found. Imported as %s\n',fname,cname);
                        else
                            conn_disp('fprintf','warning: condition file %s not found. Skipped\n',fname);
                            WRN{end+1}=sprintf('warning: condition file %s not found. Skipped\n',fname);
                        end
                    end
                end
            end
        otherwise, error('unrecognized type %s (structural|functional|conditions)',options.type);
    end
end
ok=isempty(ERR);

end


% if nargin<1||isempty(filepath), filepath=CONN_x.Setup.bids{1}; end
% 
% subjects=conn_dir(fullfile(filepath,'sub-*','-cell','-dir','-R'));
% files_task=conn_dir(fullfile(filepath,'*_events.tsv'),'-cell','-sort');
% 
% files_func=cat(2,reshape(conn_dir(fullfile(filepath,'sub-*_bold.nii'),'-cell'),1,[]),reshape(conn_dir(fullfile(filepath,'sub-*_bold.nii.gz'),'-cell'),1,[]));
% [files_func_path,files_func_name,files_func_ext]=cellfun(@fileparts,files_func,'uni',0);
% [files_func_path1,files_func_path2]=cellfun(@fileparts,files_func_path,'uni',0);
% files_func_uname=regexprep(files_func_name,'sub-\d+_?|run-\d+_?|ses-\d+_?|.nii$','');
% files_func_valid=strcmp(files_func_path2,'func');
% 
% files_anat=cat(2,reshape(conn_dir(fullfile(filepath,'sub-*.nii'),'-cell'),1,[]),reshape(conn_dir(fullfile(filepath,'sub-*.nii.gz'),'-cell'),1,[]));
% [files_anat_path,files_anat_name,files_anat_ext]=cellfun(@fileparts,files_anat,'uni',0);
% [files_anat_path1,files_anat_path2]=cellfun(@fileparts,files_func_path,'uni',0);
% files_anat_uname=regexprep(files_anat_name,'sub-\d+_?|run-\d+_?|ses-\d+_?|.nii$','');
% files_anat_valid=strcmp(files_func_path2,'anat');
