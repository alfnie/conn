function out=conn_importcondition(filename,varargin)
% CONN_IMPORTCONDITION import condition onsets/durations/name information from file
%
% FILE FORMATS ACCEPTED:
% FILE FORMAT 1: (CONN single-file condition format; *.csv or *.txt extensions)
%     .csv file (comma-separated fields)
%     one header line, and five columns: condition_name, subject_number, session_number, onsets, durations
%
% Example syntax:
%   conn_importcondition('myfile.csv');
%
% Example file contents (myfile.csv):
%   condition_name, subject_number, session_number, onsets, durations
%   task, 1, 1, 0 50 100 150, 25
%   task, 2, 1, 0 50 100 150, 25
%   rest, 1, 1, 25 75 125 175, 25
%   rest, 2, 1, 25 75 125 175, 25
%
% note: leaving the subject_number or session_number empty in a row
% indicates all subjects or all sessions: e.g.
%   condition_name, subject_number, session_number, onsets, durations
%   task, , , 0 50 100 150, 25
%   rest, , , 25 75 125 175, 25
%
% note: enter 0 in 'onsets' column, and inf in 'durations' column to indicate 
% that a condition is present during the entire session: e.g.
%   condition_name, subject_number, session_number, onsets, durations
%   rest, , , 0, inf
%
%
% FILE FORMAT 2: (BIDS .tsv one file per subject/session condition format; *_events.tsv extension)
%     .tsv file (tab-separated fields)
%     one header line, and three columns: onset, duration, trial_type
%     (notes: any additional columns are disregarded; columns may appear in any order; "trial_type" column alternative names: "event_type", "condition")
%
% Example syntax:
%   conn_importcondition('*_task_events.tsv')
%     (use this syntax if the *_events.tsv files are stored in the same folder as the Dataset-0 functional data)
%   conn_importcondition({{'sub-01/ses-01/func/sub-01_task_events.tsv','sub-01/ses-02/func/sub-01_task_events.tsv'},{'sub-02/ses-01/func/sub-02_task_events.tsv','sub-02/ses-02/func/sub-02_task_events.tsv'}});
%     (use this syntax when explicitly defining the location of each *_events.tsv file)
%
% Example file contents (sub-01/ses-01/func/sub-01_task_events.tsv):
%   onset	duration	trial_type
%   0       25          task
%   25      25          rest
%   50      25          task
%   75      25          rest
%   100     25          task
%   125     25          rest
%   150     25          task
%   175     25          rest
%
% additional options:
% conn_importcondition(...,'breakconditionsbysession',true);
%   automatically breaks each condition into session-specific conditions
% conn_importcondition(...,'deleteall',true);
%   deletes any existing condition(s) before importing the new condition-information file
%

 
global CONN_x;

options=struct('breakconditionsbysession',false,'deleteall',false,'dogui',false,'subjects',1:CONN_x.Setup.nsubjects,'sessions',[]);
for n1=1:2:nargin-1, if ~isfield(options,lower(varargin{n1})), error('unknown option %s',lower(varargin{n1})); else options.(lower(varargin{n1}))=varargin{n1+1}; end; end

% read text data
if ~nargin||isempty(filename),
    fileoptions={'CONN-legacy (single .txt or .csv file)','BIDS-compatible (single *_events.tsv file)','BIDS-compatible (one *_events.tsv file in each subject/session folder)'};
    tooltipstrings={'<HTML>Single .txt or .csv text file with subject/session/condition/onset/duration information<br/>see <i>help conn_importcondition</i> for additional file-format information</HTML>',...
        '<HTML>Single .tsv file with condition/onset/duration information (common across all subjects/sessions)<br/>see <i>help conn_importcondition</i> for additional file-format information</HTML>',...
        '<HTML>Multiple .tsv files with condition/onset/duration information (one file for each subject&session; *_events.tsv files located in same folders and matchd filenames as functional data)<br/>see <i>help conn_importcondition</i> for additional file-format information</HTML>'};
    answ=conn_questdlg('','Import task/condition information from:',fileoptions{[1:numel(fileoptions) 1]},'tooltipstring',tooltipstrings{:});
    if isempty(answ), return; end
    if isequal(answ,fileoptions{1})
        [tfilename,tpathname]=uigetfile({'*.txt','text files (*.txt)'; '*.csv','CSV-files (*.csv)'; '*',  'All Files (*)'},'Select data file');
        if ~ischar(tfilename)||isempty(tfilename), return; end
        filename=fullfile(tpathname,tfilename);
        filetype=1;
    elseif isequal(answ,fileoptions{2})
        [tfilename,tpathname]=uigetfile({'*.tsv','TSV-files (*.tsv)'; '*',  'All Files (*)'},'Select data file');
        if ~ischar(tfilename)||isempty(tfilename), return; end
        filename=fullfile(tpathname,tfilename);
        filetype=2;
    else
        filename='*_events.tsv';
        filetype=6;
    end
    options.dogui=true;
else
    if isstruct(filename), filetype=5;
    elseif iscell(filename), filetype=4;
    elseif any(filename=='*'), filetype=3;
    else
        [nill,nill,fileext]=fileparts(filename);
        if strcmp(fileext,'.tsv'), filetype=2;
        else filetype=1;
        end
    end
end
switch(filetype)
    case 1, % CONN-legacy
        [conditions,nsubs,nsess,onsets,durations]=textread(filename, '%s%d%d%s%s','delimiter',',','headerlines',1);
    case 2, % BIDS single-file
        [onsets,durations,conditions]=conn_importcondition_readbids(filename);
        nsubs=zeros(size(conditions));
        nsess=zeros(size(conditions));
    case {3,4}, % BIDS one file per subject/session
        if filetype==3 % look for tsv files in functional directory
            cwd=fileparts(CONN_x.Setup.functional{1}{1}{1});
            fmask=filename;
            f=dir(fullfile(cwd,fmask));
            if isempty(f), error('no %s file found in %s',fmask,cwd); end
            f={f.name};
            nf=1;
            if numel(f)>1,
                nf=listdlg('PromptString','Select events file:','ListString',f,'SelectionMode','single','ListSize',[200 200]);
                if isempty(nf), return; end
            end
            f0=f{nf};
            filename={};
            for isub=1:numel(options.subjects),
                nsub=options.subjects(isub);
                nsess=CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsub));
                for nses=1:nsess
                    cwd=fileparts(CONN_x.Setup.functional{nsub}{nses}{1});
                    f=dir(fullfile(cwd,fmask));
                    if isempty(f), error('No %s file found in %s',fmask,cwd); end
                    f={f.name};
                    nf=1;
                    if numel(f)>1,
                        d=conn_wordld(f0,f);
                        [nill,nf]=min(d);
                    end
                    filename{isub}{nses}=fullfile(cwd,f{nf});
                    fprintf('Subject %d Session %d: %s\n',nsub,nses,filename{isub}{nses});
                end
            end
        end
        conditions={}; onsets=[]; durations=[]; nsess=[]; nsubs=[];
        for isub=1:numel(filename)
            nsub=options.subjects(isub);
            if ~iscell(filename{isub}), filename{isub}={filename{isub}}; allsess=isempty(options.sessions);
            else allsess=false; 
            end
            for nses=1:numel(filename{isub})
                [tonsets,tdurations,tconditions]=conn_importcondition_readbids(filename{isub}{nses});
                tnsubs=nsub+zeros(size(tconditions));
                if allsess, tnsess=zeros(size(tconditions));
                elseif ~isempty(options.sessions), tnsess=options.sessions(nses)+zeros(size(tconditions));
                else tnsess=nses+zeros(size(tconditions));
                end
                conditions=cat(1,conditions, tconditions);
                onsets=cat(1,onsets, tonsets);
                durations=cat(1,durations, tdurations);
                nsubs=cat(1,nsubs, tnsubs);
                nsess=cat(1,nsess,tnsess);
            end
        end
    case 5, % manually-defined condition info
        conditions=filename.conditions;
        onsets=filename.onsets;
        durations=filename.durations;
        if isfield(filename,'nsubs'), nsubs=filename.nsubs;
        else nsubs=options.subjects;
        end
        if isfield(filename,'nsess'), nsess=filename.nsess;
        else nsess=options.sessions;
        end
        if ~iscell(conditions), conditions={conditions}; end
        if numel(conditions)==1&&numel(onsets)>1, conditions=repmat(conditions,size(onsets)); end
        if numel(onsets)==1&&numel(conditions)>1, onsets=onsets+zeros(size(conditions)); end
        if numel(durations)==1&&numel(conditions)>1, durations=durations+zeros(size(conditions)); end
        if numel(nsubs)==1&&numel(conditions)>1, nsubs=nsubs+zeros(size(conditions)); end
        if numel(nsess)==1&&numel(conditions)>1, nsess=nsess+zeros(size(conditions)); end
    case 6, % BIDS multiple files per subject/session
        [ok,err,wrn]=conn_importbids('all','type','conditions','subjects',options.subjects,'sessions',options.sessions);
        if options.dogui&&isfield(CONN_x,'gui')&&isnumeric(CONN_x.gui)&&CONN_x.gui,
            if isempty(wrn), conn_msgbox('conditions imported with no errors','Done',1); 
            else conn_msgbox({['conditions imported with ',num2str(numel(wrn)),' warnings'],'see log for details'},'Done',1); 
            end
        end
        return
end
[names,nill,nconds]=unique(conditions);

% fills-in condition info
if any(nsubs>CONN_x.Setup.nsubjects), error('Subject number in file exceeds number of subjects in current CONN project'); end
for nsub=intersect(options.subjects,1:max(nsubs))
    if any(nsess(nsubs==nsub)>CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsub))), error('Session number in file exceeds number of sessions (subject %d)',nsub); end
end
if options.deleteall, CONN_x.Setup.conditions.names={' '}; end
for isub=1:numel(options.subjects),
    nsub=options.subjects(isub);
    for nses=1:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsub))
        for ncond=1:max(nconds),
            data=find((nsubs==0|nsubs==nsub)&(nsess==0|nsess==nses)&nconds==ncond);
            if ~isempty(data)
                name=names{ncond};
                if numel(data)>1&&iscell(onsets), error('multiple rows for condition %s subject %d session %d',name,nsub,nses); end
                if options.breakconditionsbysession, name=sprintf('%s_Session%d',name,nses); end
                idx=strmatch(name,CONN_x.Setup.conditions.names,'exact');
                if isempty(idx),
                    idx=length(CONN_x.Setup.conditions.names);
                    CONN_x.Setup.conditions.names{end+1}=' ';
                end
                CONN_x.Setup.conditions.model{idx}=[];
                CONN_x.Setup.conditions.param(idx)=0;
                CONN_x.Setup.conditions.filter{idx}=[];
                CONN_x.Setup.conditions.names{idx}=name;
                if iscell(onsets)
                    t1=str2num(onsets{data}); if isempty(t1)&&~isempty(deblank(t1)), error('incorrect syntax in onset field %s (condition %s, subject %d, session %d)',onsets{data}, name, nsub, nses); end
                    CONN_x.Setup.conditions.values{nsub}{idx}{nses}{1}=t1;
                    t1=str2num(durations{data}); if isempty(t1)&&~isempty(deblank(t1)), error('incorrect syntax in duration field %s (condition %s, subject %d, session %d)',onsets{data}, name, nsub, nses); end
                    CONN_x.Setup.conditions.values{nsub}{idx}{nses}{2}=t1;
                else
                    t1=onsets(data);CONN_x.Setup.conditions.values{nsub}{idx}{nses}{1}=t1(:)';
                    t1=durations(data);CONN_x.Setup.conditions.values{nsub}{idx}{nses}{2}=t1(:)';
                end
            end
        end
    end
end

% fills possible empty conditions for each subject/session
for isub=1:numel(options.subjects),
    nsub=options.subjects(isub);
    for ncondition=1:length(CONN_x.Setup.conditions.names)-1
        for nses=1:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsub))
            if numel(CONN_x.Setup.conditions.values{nsub})<ncondition||numel(CONN_x.Setup.conditions.values{nsub}{ncondition})<nses||numel(CONN_x.Setup.conditions.values{nsub}{ncondition}{nses})<1,
                CONN_x.Setup.conditions.values{nsub}{ncondition}{nses}{1}=[];
            end
            if numel(CONN_x.Setup.conditions.values{nsub})<ncondition||numel(CONN_x.Setup.conditions.values{nsub}{ncondition})<nses||numel(CONN_x.Setup.conditions.values{nsub}{ncondition}{nses})<2,
                CONN_x.Setup.conditions.values{nsub}{ncondition}{nses}{2}=[];
            end
        end
    end
end
if options.dogui&&isfield(CONN_x,'gui')&&isnumeric(CONN_x.gui)&&CONN_x.gui, 
    conn_msgbox([num2str(numel(names)),' conditions imported with no errors'],'Done',1);
end
if nargout,out=CONN_x;end
end



