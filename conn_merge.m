
function value=conn_merge(value0,value, copyfiles, disregardcurrent, softlink, mergeinfo, skipchecks)
% CONN_MERGE 
%
% internal function:  merge conn*.mat projects
%
% conn_merge(filenames) merges one or more projects with current project
%   additional parameters: (see code for details)
%   conn_merge(filenames,[], copyfiles [=true], disregardcurrent [=false], softlink [=false], mergeinfo [=true], skipchecks [=false])
% conn_merge([],n) changes number of subjects in current project (adds or removes subjects)
% conn_merge([],[0 idx]) keeps subjects indexes by idx in current project (removes others)
%

if nargin<3||isempty(copyfiles), copyfiles=true; end % set to false if merged project should only contain project definitions but no connectivity datafiles
                                                     % set to 2 if files should not be instantly copied but instead a .sh script containing the copy commands should be generated 
if nargin<4||isempty(disregardcurrent), disregardcurrent=false; end % set to true to disregard current project info when merging multiple projects
if nargin<5||isempty(softlink), softlink=false; end % set to true if merged project should contain symbolic links to original connectivity datafiles instead of a copy of the original datafiles
if nargin<6||isempty(mergeinfo), mergeinfo=true; end % set to false if merged project definitions should not be changed (only copy connectivity datafiles)
if nargin<7||isempty(skipchecks), skipchecks=false; end % set to true if consistency checks between projects are skipped

DEBUG=false;
MAXMEM=50; % maximum number of simultaneous project files that may be loaded at a time
global CONN_x;

if nargin>0&&ischar(value0), % merge multiple project files (value0 = filenames)
    filenames=value0;
    if ~iscell(filenames), filenames=cellstr(filenames); end
    if numel(filenames)>MAXMEM&&~disregardcurrent
        for nbatch=1:MAXMEM:numel(filenames)
            value = conn_merge(char(filenames(nbatch:min(numel(filenames),nbatch+MAXMEM-1))));
        end
        return
    end
    if disregardcurrent, value0=0; 
    else value0=CONN_x.Setup.nsubjects; 
    end
    idx0=1:value0; % index to original subjects
    value=value0;
    other{1}.CONN_x=CONN_x;
    otheridx1=ones(1,value0);
    otheridx2=1:value0;
    otheridx3=ones(1,value0);
    n1=0;
    for n1a=1:size(filenames,1),
        conn_disp('fprintf','loading info from project %s\n',filenames{n1a});
        temp=load(deblank(filenames{n1a}),'CONN_x','-mat');
        if ~isfield(temp,'CONN_x'), conn_disp('fprintf','warning: invalid format, disregarding file %s\n',filenames{n1a});
        else
            n1=n1+1;
            other{1+n1}=temp;
            if isfield(other{1+n1}.CONN_x,'pobj')&&~other{1+n1}.CONN_x.pobj.holdsdata, % .dmat extended projects terminal nodes
                other{1+n1}.CONN_x.filename=conn_fullfile(conn_projectmanager('parentfile',deblank(filenames{n1a}),other{1+n1}.CONN_x.pobj));
                other{1+n1}.CONN_x=conn_updatefolders(other{1+n1}.CONN_x);
                newsubs=other{1+n1}.CONN_x.pobj.subjects;
                thissubs=newsubs;
                idx0(ismember(idx0,newsubs))=[]; % allow overwrite individual subjects info
                copyfiles=false;
            else
                other{1+n1}.CONN_x.filename=deblank(filenames{n1a});
                other{1+n1}.CONN_x=conn_updatefolders(other{1+n1}.CONN_x);
                newsubs=value+(1:other{1+n1}.CONN_x.Setup.nsubjects);
                thissubs=1:numel(newsubs);
            end
            %nsubs=other{1+n1}.CONN_x.Setup.nsubjects;
            %value=value+nsubs;
            %otheridx1=[otheridx1,1+n1+zeros(1,nsubs)];   % index to conn project
            %otheridx2=[otheridx2,1:nsubs];               % index to subject within conn project
            nsubs=numel(newsubs);
            if max(newsubs)>value, value=max(newsubs); end
            otheridx1(newsubs)=1+n1;                     % index to conn project
            otheridx2(newsubs)=thissubs;                 % index to subject within conn project
            otheridx3=[otheridx3,ones(1,nsubs)];         % 1/0 indicating whether subject data exists
        end
    end
    idx1=1:value; % index to final subjects
else % change number of subjects in current project (value0 = original number of subjects; value = new number of subjects)
    if nargin<1||isempty(value0), value0=CONN_x.Setup.nsubjects; end
    idx0=1:value0; % index to original subjects
    if numel(value)>1, idx1=value(value>0); value=numel(idx1);
    elseif numel(value)==1&&value>=value0, idx1=1:value;
    else
        names=[repmat('Subject ',[value0,1]),num2str((1:value0)')];
        Ransw=1:value0; %zeros(1,value0-value+1);
        while length(Ransw)>0&&length(Ransw)~=value0-value,
            Ransw=listdlg('name',['Removing subjects!'],'PromptString',['Select ',num2str(value0-value),' subjects to remove'],'ListString',names,'SelectionMode','multiple');
        end
        if ~isempty(Ransw), idx1=setdiff(1:value0,Ransw); % index to final subjects
        else value=value0; return; 
        end
        Ransw=conn_questdlg({'This will permanently remove any connectivity data for the selected subjects.','Are you sure you want to continue?'},'','Yes','No','No');
        if ~strcmp(Ransw,'Yes'), value=value0; return; end
    end
    other{1}.CONN_x=CONN_x;
    otheridx1=ones(1,value);
    otheridx2=min(value0,idx1);
    otheridx3=idx1<=value0;
end
if mergeinfo
    TXT={};
    if length(other)>1&&~isempty(otheridx1), % check same definitions across projects
        REF=otheridx1(1); % first project with >0 subjects
        for n1=REF:length(other),
            if numel(other{n1}.CONN_x.Setup.nsessions)==1, other{n1}.CONN_x.Setup.nsessions=other{n1}.CONN_x.Setup.nsessions+zeros(1,other{n1}.CONN_x.Setup.nsubjects); end
            if numel(other{n1}.CONN_x.Setup.RT)==1, other{n1}.CONN_x.Setup.RT=other{n1}.CONN_x.Setup.RT+zeros(1,other{n1}.CONN_x.Setup.nsubjects); end
        end
        if ~skipchecks
            for n1=REF+1:length(other),
                %if other{REF}.CONN_x.Setup.RT~=other{n1}.CONN_x.Setup.RT, TXT{end+1}=['Mismatched TR with project ',other{n1}.CONN_x.filename]; end
                if length(other{REF}.CONN_x.Setup.rois.names)~=length(other{n1}.CONN_x.Setup.rois.names) || any(any(size(strvcat(other{REF}.CONN_x.Setup.rois.names{:}))~=size(strvcat(other{n1}.CONN_x.Setup.rois.names{:})))) || any(any((strvcat(other{REF}.CONN_x.Setup.rois.names{:}))~=(strvcat(other{n1}.CONN_x.Setup.rois.names{:})))), TXT{end+1}=['Mismatched ROI names with project ',other{n1}.CONN_x.filename]; end
                if length(other{REF}.CONN_x.Setup.conditions.names)~=length(other{n1}.CONN_x.Setup.conditions.names) || any(any(size(strvcat(other{REF}.CONN_x.Setup.conditions.names{:}))~=size(strvcat(other{n1}.CONN_x.Setup.conditions.names{:})))) || any(any((strvcat(other{REF}.CONN_x.Setup.conditions.names{:}))~=(strvcat(other{n1}.CONN_x.Setup.conditions.names{:})))), TXT{end+1}=['Mismatched condition names with project ',other{n1}.CONN_x.filename]; end
                if length(other{REF}.CONN_x.Setup.l1covariates.names)~=length(other{n1}.CONN_x.Setup.l1covariates.names) || any(any(size(strvcat(other{REF}.CONN_x.Setup.l1covariates.names{:}))~=size(strvcat(other{n1}.CONN_x.Setup.l1covariates.names{:})))) || any(any((strvcat(other{REF}.CONN_x.Setup.l1covariates.names{:}))~=(strvcat(other{n1}.CONN_x.Setup.l1covariates.names{:})))), TXT{end+1}=['Mismatched first-level covariate names with project ',other{n1}.CONN_x.filename]; end
                if length(other{REF}.CONN_x.Setup.l2covariates.names)~=length(other{n1}.CONN_x.Setup.l2covariates.names) || any(any(size(strvcat(other{REF}.CONN_x.Setup.l2covariates.names{:}))~=size(strvcat(other{n1}.CONN_x.Setup.l2covariates.names{:})))) || any(any((strvcat(other{REF}.CONN_x.Setup.l2covariates.names{:}))~=(strvcat(other{n1}.CONN_x.Setup.l2covariates.names{:})))), TXT{end+1}=['Mismatched second-level covariate names with project ',other{n1}.CONN_x.filename]; end
            end
            if ~isempty(TXT),
                conn_disp(strvcat(TXT{:}));
                value=value0;
                return;
            end
        end
        if REF>1
            CONN_x.Setup.rois=other{REF}.CONN_x.Setup.rois;
            CONN_x.Setup.conditions=other{REF}.CONN_x.Setup.conditions;
            CONN_x.Setup.l1covariates=other{REF}.CONN_x.Setup.l1covariates;
            CONN_x.Setup.l2covariates=other{REF}.CONN_x.Setup.l2covariates;
        end
    end
    if length(CONN_x.Setup.nsessions)==1&&CONN_x.Setup.nsubjects>1, CONN_x.Setup.nsessions=CONN_x.Setup.nsessions+zeros(1,CONN_x.Setup.nsubjects); end
    if length(CONN_x.Setup.RT)==1&&CONN_x.Setup.nsubjects>1, CONN_x.Setup.RT=CONN_x.Setup.RT+zeros(1,CONN_x.Setup.nsubjects); end
    for n1=1:length(other), % consistency check for source/measure/condition _list info
        if ~isfield(CONN_x,'Analyses')&&isfield(other{n1}.CONN_x,'Analyses'), CONN_x.Analyses=other{n1}.CONN_x.Analyses; end
        if isfield(CONN_x,'Analyses')&&isfield(other{n1}.CONN_x,'Analyses')
            for ianalysis=1:numel(CONN_x.Analyses)
                sourcenames1=CONN_x.Analyses(ianalysis).sourcenames;
                if numel(other{n1}.CONN_x.Analyses)>=ianalysis
                    sourcenames2=other{n1}.CONN_x.Analyses(ianalysis).sourcenames;
                    if ~isequal(sourcenames1,sourcenames2)&&~isempty(sourcenames1)&&~isempty(sourcenames2)&&~isequal(sourcenames1(1:min(numel(sourcenames1),numel(sourcenames2))),sourcenames2(1:min(numel(sourcenames1),numel(sourcenames2)))), TXT{end+1}=['Incompatible source names with project ',other{n1}.CONN_x.filename]; end
                    if numel(sourcenames2)>numel(sourcenames1), CONN_x.Analyses(ianalysis).sourcenames=sourcenames2; end
                end
            end
        end
        if ~isfield(CONN_x,'vvAnalyses')&&isfield(other{n1}.CONN_x,'vvAnalyses'), CONN_x.vvAnalyses=other{n1}.CONN_x.vvAnalyses; end
        if isfield(CONN_x,'vvAnalyses')&&isfield(other{n1}.CONN_x,'vvAnalyses')
            if isfield(CONN_x.vvAnalyses,'measurenames')&&isfield(other{n1}.CONN_x.vvAnalyses,'measurenames')
                for ianalysis=1:numel(CONN_x.vvAnalyses)
                    measurenames1=CONN_x.vvAnalyses(ianalysis).measurenames;
                    if numel(other{n1}.CONN_x.vvAnalyses)>=ianalysis
                        measurenames2=other{n1}.CONN_x.vvAnalyses(ianalysis).measurenames;
                        if ~isequal(measurenames1,measurenames2)&&~isempty(measurenames1)&&~isempty(measurenames2)&&~isequal(measurenames1(1:min(numel(measurenames1),numel(measurenames2))),measurenames2(1:min(numel(measurenames1),numel(measurenames2)))), TXT{end+1}=['Incompatible voxel-to-voxel measure names with project ',other{n1}.CONN_x.filename]; end
                        if numel(measurenames2)>numel(measurenames1), CONN_x.vvAnalyses(ianalysis).measurenames=measurenames2; end
                    end
                end
            end
        end
        if isfield(CONN_x.Setup.conditions,'allnames')&&isfield(other{n1}.CONN_x.Setup.conditions,'allnames')
            allnames1=CONN_x.Setup.conditions.allnames;
            allnames2=other{n1}.CONN_x.Setup.conditions.allnames;
            if ~isequal(allnames1,allnames2)&&~isempty(allnames1)&&~isempty(allnames2)&&~isequal(allnames1(1:min(numel(allnames1),numel(allnames2))),allnames2(1:min(numel(allnames1),numel(allnames2)))), TXT{end+1}=['Incompatible condition names with project ',other{n1}.CONN_x.filename]; end
            if numel(allnames2)>numel(allnames1), CONN_x.Setup.conditions.allnames=allnames2; end
        end
    end
    if ~isempty(TXT),
        conn_disp(strvcat(TXT{:}));
        value=value0;
        return;
    end
    for nsub=setdiff(idx1,idx0), % new subjects only
        copyfromproject=otheridx1(nsub);
        copyfromsubject=otheridx2(nsub);
        if ~copyfromsubject, copyfromsubject=1; end
        CONN_x.Setup.functional{nsub}=other{copyfromproject}.CONN_x.Setup.functional{copyfromsubject};
        CONN_x.Setup.structural{nsub}=other{copyfromproject}.CONN_x.Setup.structural{copyfromsubject};
        try, for nalt=1:numel(other{copyfromproject}.CONN_x.Setup.secondarydataset), if isfield(other{copyfromproject}.CONN_x.Setup.secondarydataset(nalt),'functionals_explicit')&&~isempty(other{copyfromproject}.CONN_x.Setup.secondarydataset(nalt).functionals_explicit), CONN_x.Setup.secondarydataset(nalt).functionals_explicit{nsub}=other{copyfromproject}.CONN_x.Setup.secondarydataset(nalt).functionals_explicit{copyfromsubject}; end; end
        end
        try, if isfield(other{copyfromproject}.CONN_x.Setup,'unwarp_functional')&&~isempty(other{copyfromproject}.CONN_x.Setup.unwarp_functional), CONN_x.Setup.unwarp_functional{nsub}=other{copyfromproject}.CONN_x.Setup.unwarp_functional{copyfromsubject}; end
        end
        try, if isfield(other{copyfromproject}.CONN_x.Setup,'coregsource_functional')&&~isempty(other{copyfromproject}.CONN_x.Setup.coregsource_functional), CONN_x.Setup.coregsource_functional{nsub}=other{copyfromproject}.CONN_x.Setup.coregsource_functional{copyfromsubject}; end
        end
        CONN_x.Setup.spm{nsub}=other{copyfromproject}.CONN_x.Setup.spm{copyfromsubject};
        CONN_x.Setup.dicom{nsub}=other{copyfromproject}.CONN_x.Setup.dicom{copyfromsubject};
        CONN_x.Setup.rois.files{nsub}=other{copyfromproject}.CONN_x.Setup.rois.files{copyfromsubject};
        try, CONN_x.Setup.conditions.values{nsub}=other{copyfromproject}.CONN_x.Setup.conditions.values{copyfromsubject}; end
        CONN_x.Setup.l1covariates.files{nsub}=other{copyfromproject}.CONN_x.Setup.l1covariates.files{copyfromsubject};
        CONN_x.Setup.l2covariates.values{nsub}=other{copyfromproject}.CONN_x.Setup.l2covariates.values{copyfromsubject};
        CONN_x.Setup.nscans{nsub}=other{copyfromproject}.CONN_x.Setup.nscans{copyfromsubject};
        CONN_x.Setup.nsessions(nsub)=other{copyfromproject}.CONN_x.Setup.nsessions(min(numel(other{copyfromproject}.CONN_x.Setup.nsessions),copyfromsubject));
        CONN_x.Setup.RT(nsub)=other{copyfromproject}.CONN_x.Setup.RT(min(numel(other{copyfromproject}.CONN_x.Setup.RT),copyfromsubject));
    end
    if isempty(idx1), tempidx1=1;
    else tempidx1=idx1;
    end
    CONN_x.Setup.functional={CONN_x.Setup.functional{tempidx1}};
    CONN_x.Setup.structural={CONN_x.Setup.structural{tempidx1}};
    try, for nalt=1:numel(CONN_x.Setup.secondarydataset), if isfield(CONN_x.Setup.secondarydataset(nalt),'functionals_explicit'), CONN_x.Setup.secondarydataset(nalt).functionals_explicit={CONN_x.Setup.secondarydataset(nalt).functionals_explicit{tempidx1}}; end; end
    end
    try, CONN_x.Setup.unwarp_functional={CONN_x.Setup.unwarp_functional{tempidx1}};
    end
    try, CONN_x.Setup.coregsource_functional={CONN_x.Setup.coregsource_functional{tempidx1}};
    end
    CONN_x.Setup.spm={CONN_x.Setup.spm{tempidx1}};
    CONN_x.Setup.dicom={CONN_x.Setup.dicom{tempidx1}};
    CONN_x.Setup.rois.files={CONN_x.Setup.rois.files{tempidx1}};
    try, CONN_x.Setup.conditions.values={CONN_x.Setup.conditions.values{tempidx1}}; end
    CONN_x.Setup.l1covariates.files={CONN_x.Setup.l1covariates.files{tempidx1}};
    CONN_x.Setup.l2covariates.values={CONN_x.Setup.l2covariates.values{tempidx1}};
    CONN_x.Setup.nscans={CONN_x.Setup.nscans{tempidx1}};
    CONN_x.Setup.nsessions=CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),tempidx1));
    CONN_x.Setup.RT=CONN_x.Setup.RT(min(numel(CONN_x.Setup.RT),tempidx1));
end

% rename/copy analysis files
if copyfiles && (value<value0||any(otheridx3(setdiff(idx1,idx0)))) && ~isempty(CONN_x.filename),
    if copyfiles>1, DEBUG_FH=fopen(conn_prepend('conn_merge_',CONN_x.filename,'_copyfiles.sh'),'wt'); fclose(DEBUG_FH); end
    for nother=1:length(other),
        filepath=other{nother}.CONN_x.folders.data;
        filepathresults1=other{nother}.CONN_x.folders.preprocessing;
        filepathresults2={};for ianalysis=1:length(other{nother}.CONN_x.vvAnalyses), filepathresults2{ianalysis}=fullfile(other{nother}.CONN_x.folders.firstlevel_vv,other{nother}.CONN_x.vvAnalyses(ianalysis).name); end
        filepathresults3={};for ianalysis=1:length(other{nother}.CONN_x.Analyses), filepathresults3{ianalysis}=fullfile(other{nother}.CONN_x.folders.firstlevel,other{nother}.CONN_x.Analyses(ianalysis).name); end
        filepathresults4={};for ianalysis=1:length(other{nother}.CONN_x.dynAnalyses), filepathresults4{ianalysis}=fullfile(other{nother}.CONN_x.folders.firstlevel_dyn,other{nother}.CONN_x.dynAnalyses(ianalysis).name); end
        filepathresults5=other{nother}.CONN_x.folders.qa;
        filenameskey='Subject';
        filenamesall{nother}={};
        filesubjsall{nother}=zeros(0,4); % Subject number, Index to number in filename, Number of characters in number, Path number
        filepathall{nother}=unique({filepath,filepathresults1,filepathresults2{:},filepathresults3{:},filepathresults4{:},filepathresults5});
        for n0=1:length(filepathall{nother}),
            files=dir(fullfile(filepathall{nother}{n0},['*',filenameskey,'*']));
            for n1=1:length(files),
                k=strfind(lower(files(n1).name),lower(filenameskey));
                if ~isempty(k),
                    k=k(1);
                    m=files(n1).name(k+length(filenameskey):end);
                    [nill,idx]=min((m>='0'&m<='9'));if ~nill,m=m(1:idx-1);end
                    if ~isempty(m)&&length(str2num(m))==1,
                        filesubjsall{nother}=cat(1,filesubjsall{nother},[str2num(m),k+length(filenameskey),length(m),n0]);
                        filenamesall{nother}{end+1}=fullfile(filepathall{nother}{n0},files(n1).name);
                    end
                end
            end
        end
    end
    for nsub=1:max(value0,value),
        if nsub>length(idx1),% remove from other{1} subject# idx0(nsub)
            if ~isempty(filesubjsall{1})
                idx=find(filesubjsall{1}(:,1)==idx0(nsub));
                for n1=1:length(idx),
                    if ispc, [ok,nill]=mysystem(['del "',deblank(filenamesall{1}{idx(n1)}),'"']);
                    else, [ok,nill]=mysystem(['rm ''',deblank(filenamesall{1}{idx(n1)}),'''']);
                    end
                end
            end
        elseif nsub>length(idx0), % new subject 
            if otheridx3(nsub),%CONN_x{nsub} <- other{otheridx1(nsub)}.CONN_x{otheridx2(nsub)};
                idx=find(filesubjsall{otheridx1(nsub)}(:,1)==otheridx2(nsub));
                for n1=1:length(idx),
                    [filepath,filename,fileext]=fileparts(deblank(filenamesall{otheridx1(nsub)}{idx(n1)}));
                    newfilename=[filename(1:filesubjsall{otheridx1(nsub)}(idx(n1),2)-1),num2str(idx1(nsub),['%0',num2str(filesubjsall{otheridx1(nsub)}(idx(n1),3)),'d']),filename(filesubjsall{otheridx1(nsub)}(idx(n1),2)+filesubjsall{otheridx1(nsub)}(idx(n1),3):end)];
                    %filepath=filepathall{otheridx1(nsub)}{filesubjsall{otheridx1(nsub)}(idx(n1),4)};
                    newfilepath=filepathall{1}{filesubjsall{otheridx1(nsub)}(idx(n1),4)};
                    if ispc, [ok,nill]=mysystem(['copy "',fullfile(filepath,[filename,fileext]),'" "',fullfile(newfilepath,[newfilename,fileext]),'"']);
                    elseif softlink, [ok,nill]=mysystem(['ln -fs ''',fullfile(filepath,[filename,fileext]),''' ''',fullfile(newfilepath,[newfilename,fileext]),'''']);
                    else, [ok,nill]=mysystem(['cp ''',fullfile(filepath,[filename,fileext]),''' ''',fullfile(newfilepath,[newfilename,fileext]),'''']);
                    end
                end
            end
        elseif idx0(nsub)~=idx1(nsub), % subject renamed
            if ~isempty(filesubjsall{1})
                idx=find(filesubjsall{1}(:,1)==idx0(nsub));
                for n1=1:length(idx),
                    [filepath,filename,fileext]=fileparts(deblank(filenamesall{1}{idx(n1)}));
                    newfilename=[filename(1:filesubjsall{1}(idx(n1),2)-1),num2str(idx1(nsub),['%0',num2str(filesubjsall{1}(idx(n1),3)),'d']),filename(filesubjsall{1}(idx(n1),2)+filesubjsall{1}(idx(n1),3):end)];
                    tmp=strmatch(fullfile(filepath,[newfilename,fileext]),filenamesall{1},'exact');if isempty(tmp),newfilename=[filename(1:filesubjsall{1}(idx(n1),2)-1),num2str(idx1(nsub),['%',num2str(filesubjsall{1}(idx(n1),3)),'d']),filename(filesubjsall{1}(idx(n1),2)+filesubjsall{1}(idx(n1),3):end)];end;if isempty(strmatch(fullfile(filepath,[newfilename,fileext]),filenamesall{1},'exact')),conn_disp('warning, non existing target file'); end
                    if ispc, [ok,nill]=mysystem(['copy "',fullfile(filepath,[newfilename,fileext]),'" "',fullfile(filepath,[filename,fileext]),'"']);
                    else, [ok,nill]=mysystem(['cp ''',fullfile(filepath,[newfilename,fileext]),''' ''',fullfile(filepath,[filename,fileext]),'''']);
                    end
                end
            end
        end
    end
end
if ~nargout, CONN_x.Setup.nsubjects=value; end

    function [sys1,sys2]=mysystem(str)
        sys1=[];
        sys2=[];
        switch(DEBUG)
            case 0, 
                if copyfiles>1, DEBUG_FH=fopen(conn_prepend('conn_merge_',CONN_x.filename,'_copyfiles.sh'),'at'); fprintf(DEBUG_FH,'%s\n',str); fclose(DEBUG_FH);
                else [sys1,sys2]=system(str); disp(str);
                end
            case 1, disp(str);
        end
    end
end