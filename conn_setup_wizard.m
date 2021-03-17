function conn_setup_wizard(OPTIONS)
% obsolete fcn. tbr

% CONN_SETUP_WIZARD
% This script will take raw functional and anatomical volumes, perform
% initial preprocessing analyses (realignment, slicetime correction, 
% coregistration, segmentation, normalisation, and smoothing), and initialize 
% a conn project with this information. 
%

cwd=pwd;
DEFAULTS=struct('STEPS',[],'RT',2,'FWHM',8,'VOX',2,'NSUBJECTS',1,'NSESSIONS',1,'FUNCTIONAL_FILES',[],'STRUCTURAL_FILES',[],'CENTER',1,'REORIENT',eye(4),'CONN_NAME','','CONN_SAVE',1,'CONN_DISPLAY',1,'STRUCTURAL_TEMPLATE',fullfile(fileparts(which('spm')),'templates','T1.nii'),'FUNCTIONAL_TEMPLATE',fullfile(fileparts(which('spm')),'templates','EPI.nii'),'SO',[],'UNWARP',[],...
    'removescans',[],'reorient',[],'coregtomean',true,'applytofunctional',false,'voxelsize',2,'boundingbox',[-90,-126,-72;90,90,108],'fwhm',[],'sliceorder',[],'art_thresholds',[],'tpm_template',[],'tpm_ngaus',[],'unwarp',[]);
if isempty(dir(DEFAULTS.FUNCTIONAL_TEMPLATE)), DEFAULTS.FUNCTIONAL_TEMPLATE=fullfile(fileparts(which('spm')),'toolbox','OldNorm','EPI.nii'); end
if isempty(dir(DEFAULTS.STRUCTURAL_TEMPLATE)), DEFAULTS.STRUCTURAL_TEMPLATE=fullfile(fileparts(which('spm')),'toolbox','OldNorm','T1.nii'); end

if nargin<1, OPTIONS=[]; end

if isempty(which('spm_get_defaults'))||(~isdeployed&&strcmp(spm('ver'),'SPM5')), spm('defaults','fmri');
else spm_get_defaults; % note: in spm8 (and above?) use this instead to avoid overwriting the defaults if those have been dynamically changed
end

if ~isfield(OPTIONS,'STEPS'),OPTIONS.STEPS=DEFAULTS.STEPS;end

if isempty(OPTIONS.STEPS), OPTIONS.STEPS={};
elseif ~iscell(OPTIONS.STEPS), OPTIONS.STEPS={OPTIONS.STEPS}; 
end
if isempty(OPTIONS.STEPS)||any(cellfun('length',regexp(OPTIONS.STEPS,'^run_|^update_|^interactive_|^functional_|^structural_|^default_'))) % new preprocessing pipeline
    % initializes conn project
    if ~isfield(OPTIONS,'CONN_DISPLAY'),OPTIONS.CONN_DISPLAY=DEFAULTS.CONN_DISPLAY;end
    if ~isfield(OPTIONS,'CONN_NAME'),
        OPTIONS.CONN_NAME=DEFAULTS.CONN_NAME;
        [filename,pathname]=uiputfile({'*.mat','conn-project files (conn_*.mat)';'*','All Files (*)'},'Select New CONN project filename:',OPTIONS.CONN_NAME);
%         [filename,pathname]=uiputfile('conn_*.mat','Save new experiment data',OPTIONS.CONN_NAME);
        if isequal(filename,0), return; end
        OPTIONS.CONN_NAME=fullfile(pathname,filename);
    elseif ~isempty(dir(OPTIONS.CONN_NAME)),
        Ransw=conn_questdlg([OPTIONS.CONN_NAME,' project already exists, Overwrite?'],'warning','Yes','No','No');
        if isempty(Ransw)||strcmp(Ransw,'No'), return; end
    end
    if ~isfield(OPTIONS,'RT'),
        OPTIONS.RT=DEFAULTS.RT;
        answ=inputdlg('Enter Repetition-Time (in seconds)','conn_setup_wizard',1,{num2str(OPTIONS.RT)});
        if isempty(answ), return; end
        OPTIONS.RT=str2num(answ{1});
    end
    if ~isfield(OPTIONS,'voxelsize'), OPTIONS.voxelsize=DEFAULTS.voxelsize; end
    if ~isfield(OPTIONS,'boundingbox'), OPTIONS.boundingbox=DEFAULTS.boundingbox; end
    if isfield(OPTIONS,'FWHM'), OPTIONS.fwhm=OPTIONS.FWHM; end
    if ~isfield(OPTIONS,'fwhm'), OPTIONS.fwhm=DEFAULTS.fwhm; end
    if isfield(OPTIONS,'SO'), OPTIONS.sliceorder=OPTIONS.SO; end
    if isfield(OPTIONS,'unwarp'), OPTIONS.unwarp=OPTIONS.UNWARP; end
    if ~isfield(OPTIONS,'sliceorder'), OPTIONS.sliceorder=DEFAULTS.sliceorder; end
    if ~isfield(OPTIONS,'unwarp'), OPTIONS.unwarp=DEFAULTS.unwarp; end
    if ~isfield(OPTIONS,'removescans'), OPTIONS.removescans=DEFAULTS.removescans; end
    if ~isfield(OPTIONS,'reorient'), OPTIONS.reorient=DEFAULTS.reorient; end
    if ~isfield(OPTIONS,'coregtomean'), OPTIONS.coregtomean=DEFAULTS.coregtomean; end
    if ~isfield(OPTIONS,'applytofunctional'), OPTIONS.applytofunctional=DEFAULTS.applytofunctional; end
    if ~isfield(OPTIONS,'art_thresholds'), OPTIONS.art_thresholds=DEFAULTS.art_thresholds; end
    if isfield(OPTIONS,'FUNCTIONAL_FILES'), 
        OPTIONS.NSUBJECTS=length(OPTIONS.FUNCTIONAL_FILES);
    elseif ~isfield(OPTIONS,'NSUBJECTS'),   
        OPTIONS.NSUBJECTS=DEFAULTS.NSUBJECTS;
        answ=inputdlg('Enter number of subjects','conn_setup_wizard',1,{num2str(OPTIONS.NSUBJECTS)});
        if isempty(answ), return; end
        OPTIONS.NSUBJECTS=str2num(answ{1});
    end
    
    clear batch;
    batch.filename=OPTIONS.CONN_NAME;
    batch.Setup.isnew=1;
    batch.Setup.nsubjects=OPTIONS.NSUBJECTS;
    batch.Setup.RT=OPTIONS.RT;
    batch.Setup.done=0;
    if any(ismember(OPTIONS.STEPS,{'default_ss','default_ssphase'})),  
        batch.Setup.voxelmask=2;        % implicit mask
        batch.Setup.voxelresolution=4;  % (4): surface-based (3):same as functional
        %batch.Setup.rois.functionals_type=1;  % extract ROI from same files
        batch.Setup.rois={};            % no default ROIs (they are in MNI space)
    end
    conn_batch(batch);

    % loads raw functional/anatomical files
    if ~isfield(OPTIONS,'FUNCTIONAL_FILES'),
        OPTIONS.FUNCTIONAL_FILES=cell(OPTIONS.NSUBJECTS,1);
        if ~isfield(OPTIONS,'NSESSIONS'),OPTIONS.NSESSIONS=DEFAULTS.NSESSIONS;end
        OPTIONS.NSESSIONS=OPTIONS.NSESSIONS(min(length(OPTIONS.NSESSIONS),1:OPTIONS.NSUBJECTS));
        for nsubject=1:OPTIONS.NSUBJECTS,
            temp=inputdlg(['Subject ',num2str(nsubject),': Enter number of sessions'],'conn_setup_wizard',1,{num2str(OPTIONS.NSESSIONS(min(length(OPTIONS.NSESSIONS),nsubject)))});
            if isempty(temp), return; end
            temp=str2num(temp{1});
            OPTIONS.NSESSIONS(nsubject)=temp;
            OPTIONS.FUNCTIONAL_FILES{nsubject}=cell(OPTIONS.NSESSIONS(min(length(OPTIONS.NSESSIONS),nsubject)),1);
            for nsession=1:OPTIONS.NSESSIONS(min(length(OPTIONS.NSESSIONS),nsubject)),
                OPTIONS.FUNCTIONAL_FILES{nsubject}{nsession}=cellstr(spm_select(Inf,'\.img$|\.nii$',['SUBJECT ',num2str(nsubject),' SESSION ',num2str(nsession),' functional volumes'],OPTIONS.FUNCTIONAL_FILES{nsubject}{nsession}));
                if isempty(OPTIONS.FUNCTIONAL_FILES{nsubject}{nsession}{1}),return;end
            end
        end
    end
    if ~isfield(OPTIONS,'STRUCTURAL_FILES'),
        OPTIONS.STRUCTURAL_FILES=cell(OPTIONS.NSUBJECTS,1);
        for nsubject=1:OPTIONS.NSUBJECTS,
            OPTIONS.STRUCTURAL_FILES{nsubject}=cellstr(spm_select(1,'\.img$|\.nii$',['SUBJECT ',num2str(nsubject),' structural volume'],OPTIONS.STRUCTURAL_FILES{nsubject}));
            if isempty(OPTIONS.STRUCTURAL_FILES{nsubject}{1}),return;end
        end
    else,
        for nsub=1:length(OPTIONS.STRUCTURAL_FILES),if ~iscell(OPTIONS.STRUCTURAL_FILES{nsub}),OPTIONS.STRUCTURAL_FILES{nsub}=cellstr(OPTIONS.STRUCTURAL_FILES{nsub});end;end
    end
    if ~isfield(OPTIONS,'STRUCTURAL_TEMPLATE'),OPTIONS.STRUCTURAL_TEMPLATE=DEFAULTS.STRUCTURAL_TEMPLATE;end
    if ~isfield(OPTIONS,'FUNCTIONAL_TEMPLATE'),OPTIONS.FUNCTIONAL_TEMPLATE=DEFAULTS.FUNCTIONAL_TEMPLATE;end
    if ~isfield(OPTIONS,'tpm_template'),OPTIONS.tpm_template=DEFAULTS.tpm_template;end
    if ~isfield(OPTIONS,'tpm_ngaus'),OPTIONS.tpm_ngaus=DEFAULTS.tpm_ngaus;end
    if iscell(OPTIONS.STRUCTURAL_TEMPLATE),OPTIONS.STRUCTURAL_TEMPLATE=char(OPTIONS.STRUCTURAL_TEMPLATE);end
    if iscell(OPTIONS.FUNCTIONAL_TEMPLATE),OPTIONS.FUNCTIONAL_TEMPLATE=char(OPTIONS.FUNCTIONAL_TEMPLATE);end
    for nsub=1:length(OPTIONS.FUNCTIONAL_FILES),for nses=1:length(OPTIONS.FUNCTIONAL_FILES{nsub}),if ~iscell(OPTIONS.FUNCTIONAL_FILES{nsub}{nses}),OPTIONS.FUNCTIONAL_FILES{nsub}{nses}=cellstr(OPTIONS.FUNCTIONAL_FILES{nsub}{nses});end;end;end;
%     if ~isfield(OPTIONS,'FWHM'),OPTIONS.FWHM=[]; end
%     if ~isfield(OPTIONS,'SO'),OPTIONS.SO=DEFAULTS.SO;end
    
    clear batch;
    batch.filename=OPTIONS.CONN_NAME;
    batch.Setup.isnew=0;
    batch.Setup.functionals=OPTIONS.FUNCTIONAL_FILES; 
    batch.Setup.structurals=OPTIONS.STRUCTURAL_FILES;
    batch.Setup.conditions.names{1}='rest';
    for nsubject=1:length(OPTIONS.FUNCTIONAL_FILES),for nsess=1:length(OPTIONS.FUNCTIONAL_FILES{nsubject}),
            batch.Setup.conditions.onsets{1}{nsubject}{nsess}=0;
            batch.Setup.conditions.durations{1}{nsubject}{nsess}=inf;
        end;end
    batch.Setup.done=0;
    conn_batch(batch);
    conn save;
    
    % runs preprocessing steps
    if isempty(OPTIONS.STEPS)
        opts={'MNI-space template (default)','Subject-space template','Blank template (manually define)'};
        answ=conn_questdlg('Select preprocessing pipeline template:','',opts{:},opts{1});
        if isempty(answ), return; end
        [nill,answ]=ismember(answ,opts);
        switch(answ)
            case 1, OPTIONS.STEPS='default_mni';
            case 2, OPTIONS.STEPS='default_ss'; 
            case 3, OPTIONS.STEPS=''; OPTIONS.CONN_DISPLAY=1; 
        end
    end
    conn_setup_preproc(OPTIONS.STEPS,'multiplesteps',1,'dogui',OPTIONS.CONN_DISPLAY,'structural_template',OPTIONS.STRUCTURAL_TEMPLATE,'functional_template',OPTIONS.FUNCTIONAL_TEMPLATE,'tpm_template',OPTIONS.tpm_template,'tpm_ngaus',OPTIONS.tpm_ngaus,...
       'voxelsize',OPTIONS.voxelsize, 'boundingbox',OPTIONS.boundingbox, 'fwhm',OPTIONS.fwhm, 'sliceorder',OPTIONS.sliceorder, 'unwarp',OPTIONS.unwarp,'removescans',OPTIONS.removescans,'reorient',OPTIONS.reorient,'coregtomean',OPTIONS.coregtomean,'applytofunctional',OPTIONS.applytofunctional,'art_thresholds',OPTIONS.art_thresholds);
%     if isempty(OPTIONS.STEPS)
%         steps=[];
%     elseif any(ismember(OPTIONS.STEPS,{'defaultSS'})),  
%         steps='defaultSS';
%     else
%         steps='defaultMNI';
%     end
%     conn_setup_preproc(steps,'multiplesteps',1,'dogui',OPTIONS.CONN_DISPLAY,'fwhm',OPTIONS.FWHM,'sliceorder',OPTIONS.SO,'structural_template',OPTIONS.STRUCTURAL_TEMPLATE,'functional_template',OPTIONS.FUNCTIONAL_TEMPLATE);
    conn_disp(['****************************************']);
    conn_disp(['Finished all preprocessing steps']);
    conn_disp(['****************************************']);
    conn save;
    if OPTIONS.CONN_DISPLAY,
        conn
        conn('load',OPTIONS.CONN_NAME);
        conn gui_setup
    end
    
    
else % old preprocessing pipeline
    conn_disp('Warning: this preprocessing pipeline will be discontinued in future releases. See help conn_batch for new syntax');
    %if ~isfield(OPTIONS,'CONN_SAVE'),OPTIONS.CONN_SAVE=DEFAULTS.CONN_SAVE;end
    if ~isfield(OPTIONS,'CONN_DISPLAY'),OPTIONS.CONN_DISPLAY=DEFAULTS.CONN_DISPLAY;end
    if ~isfield(OPTIONS,'CENTER'),OPTIONS.CENTER=DEFAULTS.CENTER;end
    if ~isfield(OPTIONS,'REORIENT'),OPTIONS.REORIENT=DEFAULTS.REORIENT;end
    STEPS={'segmentation','slicetiming','realignment','coregistration','normalization','smoothing','initialization'};
    if isempty(OPTIONS.STEPS),
        if OPTIONS.CONN_DISPLAY,
            answ=listdlg('Promptstring','Select preprocessing steps','selectionmode','multiple','liststring',STEPS,'initialvalue',1:length(STEPS));
            OPTIONS.STEPS={STEPS{answ}};
        else, OPTIONS.STEPS=STEPS; end
    else,
        OPTIONS.STEPS=lower(OPTIONS.STEPS);
    end
    
    if ~isempty(strmatch('initialization',OPTIONS.STEPS,'exact')),OPTIONS.CONN_SAVE=1;else,OPTIONS.CONN_SAVE=0;end
    if OPTIONS.CONN_SAVE&&~isfield(OPTIONS,'CONN_NAME'),OPTIONS.CONN_NAME=DEFAULTS.CONN_NAME;[filename,pathname]=uiputfile('conn_*.mat','Save new experiment data',OPTIONS.CONN_NAME);OPTIONS.CONN_NAME=fullfile(pathname,filename);end;
    if ~isfield(OPTIONS,'RT'),OPTIONS.RT=DEFAULTS.RT;OPTIONS.RT=inputdlg('Enter Repetition-Time (in seconds)','conn_setup_wizard',1,{num2str(OPTIONS.RT)});OPTIONS.RT=str2num(OPTIONS.RT{1});end
    if ~isfield(OPTIONS,'FWHM'),OPTIONS.FWHM=DEFAULTS.FWHM;if ~isempty(strmatch('smoothing',OPTIONS.STEPS,'exact')),OPTIONS.FWHM=inputdlg('Enter smoothing FWHM (in mm)','conn_setup_wizard',1,{num2str(OPTIONS.FWHM)});OPTIONS.FWHM=str2num(OPTIONS.FWHM{1});end;end
    if ~isfield(OPTIONS,'VOX'),OPTIONS.VOX=DEFAULTS.VOX;if ~isempty(strmatch('coregistration',OPTIONS.STEPS,'exact'))||~isempty(strmatch('normalization',OPTIONS.STEPS,'exact')),OPTIONS.VOX=inputdlg('Enter resampled voxel size (in mm)','conn_setup_wizard',1,{num2str(OPTIONS.VOX)});OPTIONS.VOX=str2num(OPTIONS.VOX{1});end;end
    if ~isfield(OPTIONS,'SO'),OPTIONS.SO=DEFAULTS.SO;end
    save connsetup_wizard_job_setup.mat OPTIONS
    if ~isfield(OPTIONS,'FUNCTIONAL_FILES'),
        if ~isfield(OPTIONS,'NSUBJECTS'),OPTIONS.NSUBJECTS=DEFAULTS.NSUBJECTS;OPTIONS.NSUBJECTS=inputdlg('Enter number of subjects','conn_setup_wizard',1,{num2str(OPTIONS.NSUBJECTS)});OPTIONS.NSUBJECTS=str2num(OPTIONS.NSUBJECTS{1});end
        OPTIONS.FUNCTIONAL_FILES=cell(OPTIONS.NSUBJECTS,1);
        if ~isfield(OPTIONS,'NSESSIONS'),OPTIONS.NSESSIONS=DEFAULTS.NSESSIONS;end
        OPTIONS.NSESSIONS=OPTIONS.NSESSIONS(min(length(OPTIONS.NSESSIONS),1:OPTIONS.NSUBJECTS));
        for nsubject=1:OPTIONS.NSUBJECTS,
            temp=inputdlg(['Subject ',num2str(nsubject),' Enter number of sessions'],'conn_setup_wizard',1,{num2str(OPTIONS.NSESSIONS(min(length(OPTIONS.NSESSIONS),nsubject)))});temp=str2num(temp{1});OPTIONS.NSESSIONS(nsubject)=temp;
            OPTIONS.FUNCTIONAL_FILES{nsubject}=cell(OPTIONS.NSESSIONS(min(length(OPTIONS.NSESSIONS),nsubject)),1);
            for nsession=1:OPTIONS.NSESSIONS(min(length(OPTIONS.NSESSIONS),nsubject)),
                OPTIONS.FUNCTIONAL_FILES{nsubject}{nsession}=cellstr(spm_select(Inf,'\.img$|\.nii$',['SUBJECT ',num2str(nsubject),' SESSION ',num2str(nsession),' functional volumes'],OPTIONS.FUNCTIONAL_FILES{nsubject}{nsession}));
                if isempty(OPTIONS.FUNCTIONAL_FILES{nsubject}{nsession}{1}),return;end
            end
        end
    end
    
    for nsub=1:length(OPTIONS.FUNCTIONAL_FILES),for nses=1:length(OPTIONS.FUNCTIONAL_FILES{nsub}),if ~iscell(OPTIONS.FUNCTIONAL_FILES{nsub}{nses}),OPTIONS.FUNCTIONAL_FILES{nsub}{nses}=cellstr(OPTIONS.FUNCTIONAL_FILES{nsub}{nses});end;end;end;
    %OPTIONS.XFUNCTIONAL_FILES=OPTIONS.FUNCTIONAL_FILES;for nsub=1:length(OPTIONS.FUNCTIONAL_FILES),for nses=1:length(OPTIONS.FUNCTIONAL_FILES{nsub}),[tempa,tempb,tempc]=fileparts(OPTIONS.FUNCTIONAL_FILES{nsub}{nses}{1}); if length(OPTIONS.FUNCTIONAL_FILES{nsub}{nses})==1&&strcmp(tempc,'.nii'),OPTIONS.XFUNCTIONAL_FILES{nsub}{nses}=cfg_getfile('ExtFPList',tempa,['^',tempb,tempc],1:1e4); end;end;end
    OPTIONS.XFUNCTIONAL_FILES=OPTIONS.FUNCTIONAL_FILES;for nsub=1:length(OPTIONS.FUNCTIONAL_FILES),for nses=1:length(OPTIONS.FUNCTIONAL_FILES{nsub}),[tempa,tempb,tempc]=fileparts(OPTIONS.FUNCTIONAL_FILES{nsub}{nses}{1}); if length(OPTIONS.FUNCTIONAL_FILES{nsub}{nses})==1&&strcmp(tempc,'.nii'),OPTIONS.XFUNCTIONAL_FILES{nsub}{nses}=cellstr(spm_select('ExtFPList',tempa,['^',tempb,tempc],1:1e4)); end;end;end
    OPTIONS.NSUBJECTS=length(OPTIONS.FUNCTIONAL_FILES);
    if ~isfield(OPTIONS,'STRUCTURAL_FILES'),
        OPTIONS.STRUCTURAL_FILES=cell(OPTIONS.NSUBJECTS,1);
        for nsubject=1:OPTIONS.NSUBJECTS,
            OPTIONS.STRUCTURAL_FILES{nsubject}=cellstr(spm_select(1,'\.img$|\.nii$',['SUBJECT ',num2str(nsubject),' structural volume'],OPTIONS.STRUCTURAL_FILES{nsubject}));
            if isempty(OPTIONS.STRUCTURAL_FILES{nsubject}{1}),return;end
        end
    else,
        for nsub=1:length(OPTIONS.STRUCTURAL_FILES),if ~iscell(OPTIONS.STRUCTURAL_FILES{nsub}),OPTIONS.STRUCTURAL_FILES{nsub}=cellstr(OPTIONS.STRUCTURAL_FILES{nsub});end;end
    end
    if ~isfield(OPTIONS,'STRUCTURAL_TEMPLATE'),OPTIONS.STRUCTURAL_TEMPLATE=DEFAULTS.STRUCTURAL_TEMPLATE;if ~isempty(strmatch('coregistration',OPTIONS.STEPS,'exact'))||~isempty(strmatch('normalization',OPTIONS.STEPS,'exact')),OPTIONS.STRUCTURAL_TEMPLATE=spm_select(1,'\.img$|\.nii$','anatomical template',{OPTIONS.STRUCTURAL_TEMPLATE});end;end
    if ~isfield(OPTIONS,'FUNCTIONAL_TEMPLATE'),OPTIONS.FUNCTIONAL_TEMPLATE=DEFAULTS.FUNCTIONAL_TEMPLATE;if ~isempty(strmatch('coregistration',OPTIONS.STEPS,'exact'))||~isempty(strmatch('normalization',OPTIONS.STEPS,'exact')),OPTIONS.FUNCTIONAL_TEMPLATE=spm_select(1,'\.img$|\.nii$','functional template',{OPTIONS.FUNCTIONAL_TEMPLATE});end;end
    if iscell(OPTIONS.STRUCTURAL_TEMPLATE),OPTIONS.STRUCTURAL_TEMPLATE=char(OPTIONS.STRUCTURAL_TEMPLATE);end
    if iscell(OPTIONS.FUNCTIONAL_TEMPLATE),OPTIONS.FUNCTIONAL_TEMPLATE=char(OPTIONS.FUNCTIONAL_TEMPLATE);end
    save connsetup_wizard_job_setup.mat OPTIONS
    
    %% SPM preprocessing batch process
    for nsubject=1:OPTIONS.NSUBJECTS,
        conn_disp(['********************************']);
        conn_disp(['Preparing subject ',num2str(nsubject)]);
        conn_disp(['********************************']);
        matlabbatch={};
        prefix={[],[],[]};
        % Coregister Anatomical
        if ~isempty(strmatch('coregistration',OPTIONS.STEPS,'exact'))||~isempty(strmatch('normalization',OPTIONS.STEPS,'exact')),
            if OPTIONS.CENTER||any(any(OPTIONS.REORIENT~=eye(4))),
                a=spm_vol(OPTIONS.STRUCTURAL_FILES{nsubject}{1});b=spm_read_vols(a);
                a.mat=OPTIONS.REORIENT*a.mat;
                if OPTIONS.CENTER, a.mat(1:3,4)=-a.mat(1:3,1:3)*a.dim'/2; end
                spm_write_vol(a,b);
            end
            matlabbatch{end+1}.spm.spatial.coreg.estimate.ref={OPTIONS.STRUCTURAL_TEMPLATE};
            matlabbatch{end}.spm.spatial.coreg.estimate.source=OPTIONS.STRUCTURAL_FILES{nsubject};
        end
        % Segment Anatomical
        if ~isempty(strmatch('segmentation',OPTIONS.STEPS,'exact')),
            matlabbatch{end+1}.spm.spatial.preproc.data=OPTIONS.STRUCTURAL_FILES{nsubject};
            matlabbatch{end}.spm.spatial.preproc.output.GM=[1,1,1];
            matlabbatch{end}.spm.spatial.preproc.output.WM=[1,1,1];
            matlabbatch{end}.spm.spatial.preproc.output.CSF=[1,1,1];
            prefix{2}=['m',prefix{2}];
        end
        % Skull-stripped anatomical
        if ~isempty(strmatch('coregistration',OPTIONS.STEPS,'exact')),
            matlabbatch{end+1}.spm.util.imcalc.expression='(i1+i2+i3).*i4';
            matlabbatch{end}.spm.util.imcalc.input=cat(1,conn_prepend('c1',OPTIONS.STRUCTURAL_FILES{nsubject}),conn_prepend('c2',OPTIONS.STRUCTURAL_FILES{nsubject}),conn_prepend('c3',OPTIONS.STRUCTURAL_FILES{nsubject}),OPTIONS.STRUCTURAL_FILES{nsubject});
            matlabbatch{end}.spm.util.imcalc.output=conn_prepend('c0',OPTIONS.STRUCTURAL_FILES{nsubject}{1});
            matlabbatch{end}.spm.util.imcalc.options.dtype=spm_type('float32');
        end
        % Corrects slice-timing Functional
        if ~isempty(strmatch('slicetiming',OPTIONS.STEPS,'exact')),
            nslices=zeros(1,length(OPTIONS.XFUNCTIONAL_FILES{nsubject}));
            for nses=1:length(OPTIONS.XFUNCTIONAL_FILES{nsubject}),
                tempvol=spm_vol(OPTIONS.XFUNCTIONAL_FILES{nsubject}{nses}{1});
                nslices(nses)=tempvol.dim(3);
            end
            uniquenslices=unique(nslices);
            so=OPTIONS.SO;
            if iscell(so), so=so{min(numel(so),nsubject)}; end
            for nslice=uniquenslices(:)',
                nses=find(nslices==nslice);
                while numel(unique(so))~=nslice||max(so)~=nslice||min(so)~=1,
                    if isempty(so), so=1:nslice; end
                    so=inputdlg(['Subject ',num2str(nsubject),', ','Sessions ',num2str(nses(:)'),': Acquisition order? (1=first slice in image; ',num2str(nslice),'=last slice)'],'conn_setup_wizard',1,{num2str(so)});so=str2num(so{1});
                end
                matlabbatch{end+1}.spm.temporal.st.scans=OPTIONS.XFUNCTIONAL_FILES{nsubject}(nses);
                matlabbatch{end}.spm.temporal.st.tr=OPTIONS.RT(min(numel(OPTIONS.RT),nsubject));
                matlabbatch{end}.spm.temporal.st.nslices=nslice;
                matlabbatch{end}.spm.temporal.st.ta=OPTIONS.RT(min(numel(OPTIONS.RT),nsubject))*(1-1/nslice);
                matlabbatch{end}.spm.temporal.st.refslice=floor(nslice/2);
                matlabbatch{end}.spm.temporal.st.so=so;
            end
            prefix{1}=['a',prefix{1}];
        end
        prefix1=prefix{1};
        % Realign Functional
        if ~isempty(strmatch('realignment',OPTIONS.STEPS,'exact')),
            matlabbatch{end+1}.spm.spatial.realign.estwrite.data=conn_prepend(prefix{1},OPTIONS.XFUNCTIONAL_FILES{nsubject});
            matlabbatch{end}.spm.spatial.realign.estwrite.eoptions.rtm=0;
            matlabbatch{end}.spm.spatial.realign.estwrite.roptions.which=[2,1];
            prefix{1}=['r',prefix{1}];
        end
        % Coregister Functional
        if ~isempty(strmatch('coregistration',OPTIONS.STEPS,'exact')),
            matlabbatch{end+1}.spm.spatial.coreg.estimate.ref=conn_prepend('c0',OPTIONS.STRUCTURAL_FILES{nsubject});
            matlabbatch{end}.spm.spatial.coreg.estimate.source=conn_prepend(['mean',prefix1],{OPTIONS.XFUNCTIONAL_FILES{nsubject}{1}{1}});
            temp=conn_prepend(prefix{1},OPTIONS.XFUNCTIONAL_FILES{nsubject});
            matlabbatch{end}.spm.spatial.coreg.estimate.other=cat(1,temp{:});
        end
        % Normalize Structural/Functional
        if ~isempty(strmatch('normalization',OPTIONS.STEPS,'exact')),
            matlabbatch{end+1}.spm.spatial.normalise.write.subj.matname=conn_prepend('',OPTIONS.STRUCTURAL_FILES{nsubject},'_seg_sn.mat');
            matlabbatch{end}.spm.spatial.normalise.write.subj.resample=cat(1,conn_prepend(prefix{2},OPTIONS.STRUCTURAL_FILES{nsubject}));
            matlabbatch{end}.spm.spatial.normalise.write.roptions.vox=OPTIONS.VOX*[1 1 1];
            %matlabbatch{end}.spm.spatial.normalise.write.roptions.bb=[-90,-126,-72;90,90,108];
            
            matlabbatch{end+1}.spm.spatial.normalise.estwrite.eoptions.template={OPTIONS.FUNCTIONAL_TEMPLATE};
            matlabbatch{end}.spm.spatial.normalise.estwrite.subj.source=conn_prepend(['mean',prefix1],{OPTIONS.XFUNCTIONAL_FILES{nsubject}{1}{1}});
            temp=conn_prepend(prefix{1},OPTIONS.XFUNCTIONAL_FILES{nsubject});
            matlabbatch{end}.spm.spatial.normalise.estwrite.subj.resample=cat(1,temp{:},conn_prepend(['mean',prefix1],{OPTIONS.XFUNCTIONAL_FILES{nsubject}{1}{1}}));
            matlabbatch{end}.spm.spatial.normalise.estwrite.roptions.vox=OPTIONS.VOX*[1 1 1];
            %matlabbatch{end}.spm.spatial.normalise.estwrite.roptions.bb=[-90,-126,-72;90,90,108];
            prefix{1}=['w',prefix{1}];
            prefix{2}=['w',prefix{2}];
            prefix{3}=['w',prefix{3}];
        elseif ~isempty(strmatch('normalization_old',OPTIONS.STEPS,'exact')),
            matlabbatch{end+1}.spm.spatial.normalise.write.subj.matname=conn_prepend('',OPTIONS.STRUCTURAL_FILES{nsubject},'_seg_sn.mat');
            temp=conn_prepend(prefix{1},OPTIONS.XFUNCTIONAL_FILES{nsubject});
            matlabbatch{end}.spm.spatial.normalise.write.subj.resample=cat(1,temp{:},conn_prepend('mean',{OPTIONS.XFUNCTIONAL_FILES{nsubject}{1}{1}}),conn_prepend(prefix{2},OPTIONS.STRUCTURAL_FILES{nsubject}));
            matlabbatch{end}.spm.spatial.normalise.write.roptions.vox=OPTIONS.VOX*[1 1 1];
            prefix{1}=['w',prefix{1}];
            prefix{2}=['w',prefix{2}];
            prefix{3}=['w',prefix{3}];
        elseif ~isempty(strmatch('coregistration',OPTIONS.STEPS,'exact')),
            t=load(char(conn_prepend('',OPTIONS.STRUCTURAL_FILES{nsubject},'_seg_sn.mat')));t.Tr=[];save(char(conn_prepend('',OPTIONS.STRUCTURAL_FILES{nsubject},'_seg_sn_aff.mat')),'-struct','t');
            matlabbatch{end+1}.spm.spatial.normalise.write.subj.matname=conn_prepend('',OPTIONS.STRUCTURAL_FILES{nsubject},'_seg_sn_aff.mat');
            temp=conn_prepend(prefix{1},OPTIONS.XFUNCTIONAL_FILES{nsubject});
            matlabbatch{end}.spm.spatial.normalise.write.subj.resample=cat(1,temp{:},conn_prepend(['mean',prefix1],{OPTIONS.XFUNCTIONAL_FILES{nsubject}{1}{1}}),conn_prepend(prefix{2},OPTIONS.STRUCTURAL_FILES{nsubject}));
            matlabbatch{end}.spm.spatial.normalise.write.roptions.vox=OPTIONS.VOX*[1 1 1];
            matlabbatch{end}.spm.spatial.normalise.write.roptions.prefix='r';
            prefix{1}=['r',prefix{1}];
            prefix{2}=['r',prefix{2}];
        end
        % Smooth Functional
        if ~isempty(strmatch('smoothing',OPTIONS.STEPS,'exact')),
            temp=conn_prepend(prefix{1},OPTIONS.XFUNCTIONAL_FILES{nsubject});
            matlabbatch{end+1}.spm.spatial.smooth.data=cat(1,temp{:});
            matlabbatch{end}.spm.spatial.smooth.fwhm=OPTIONS.FWHM*[1,1,1];
            prefix{1}=['s',prefix{1}];
        end
        
        save(['setup_conn_wizard_job',num2str(nsubject),'.mat'],'matlabbatch','prefix');
    end
    spm_jobman('initcfg');
    for nsubject=1:OPTIONS.NSUBJECTS,
        conn_disp(['********************************']);
        conn_disp(['Preprocessing subject ',num2str(nsubject)]);
        conn_disp(['********************************']);
        load(['setup_conn_wizard_job',num2str(nsubject),'.mat'],'matlabbatch','prefix');
        if ~isempty(matlabbatch),
            spm_jobman('run',matlabbatch);
        end
    end
    
    %% CONN Setup definition batch process
    if OPTIONS.CONN_SAVE,
        conn_disp(['********************************']);
        conn_disp(['Initializing conn project ',OPTIONS.CONN_NAME]);
        conn_disp(['********************************']);
        clear batch;
        batch.filename=OPTIONS.CONN_NAME;
        if ~isempty(dir(batch.filename)),
            Ransw=conn_questdlg([OPTIONS.CONN_NAME,' project already exists, Overwrite?'],'warning','Yes','No','No');
            if strcmp(Ransw,'Yes'), batch.Setup.isnew=1; else, batch.Setup.isnew=0; end
        else,batch.Setup.isnew=1;end
        batch.Setup.nsubjects=OPTIONS.NSUBJECTS;
        if ~isempty(strmatch('coregistration',OPTIONS.STEPS,'exact'))&&isempty(strmatch('normalization',OPTIONS.STEPS,'exact')),batch.Setup.normalized=0;else,batch.Setup.normalized=1;end
        batch.Setup.RT=OPTIONS.RT;
        batch.Setup.functionals=conn_prepend(prefix{1},OPTIONS.FUNCTIONAL_FILES);
        batch.Setup.structurals=conn_prepend(prefix{2},OPTIONS.STRUCTURAL_FILES);
        if ~isempty(strmatch('segmentation',OPTIONS.STEPS,'exact')),
            batch.Setup.masks.Grey=conn_prepend([prefix{3},'c1'],OPTIONS.STRUCTURAL_FILES);
            batch.Setup.masks.White=conn_prepend([prefix{3},'c2'],OPTIONS.STRUCTURAL_FILES);
            batch.Setup.masks.CSF=conn_prepend([prefix{3},'c3'],OPTIONS.STRUCTURAL_FILES);
        end
        if ~isempty(strmatch('realignment',OPTIONS.STEPS,'exact')),
            batch.Setup.covariates.names{1}='realignment';
            for nsubject=1:length(OPTIONS.FUNCTIONAL_FILES),
                for nsess=1:length(OPTIONS.FUNCTIONAL_FILES{nsubject}),
                    batch.Setup.covariates.files{1}{nsubject}{nsess}=conn_prepend(['rp_',prefix1],OPTIONS.FUNCTIONAL_FILES{nsubject}{nsess}{1},'.txt');
                end
            end
        end
        batch.Setup.conditions.names{1}='rest';
        for nsubject=1:length(OPTIONS.FUNCTIONAL_FILES),for nsess=1:length(OPTIONS.FUNCTIONAL_FILES{nsubject}),
                batch.Setup.conditions.onsets{1}{nsubject}{nsess}=0;
                batch.Setup.conditions.durations{1}{nsubject}{nsess}=inf;
            end;end
        batch.Setup.done=0;
        conn_batch(batch);
        
        if OPTIONS.CONN_DISPLAY,
            conn
            conn('load',OPTIONS.CONN_NAME);
            conn gui_setup
        end
    end
    
    conn_disp(['********************************']);
    conn_disp(['Done']);
    conn_disp(['********************************']);
end


