function SPM=conn_updatespmfilepaths(filename,rule_SWD2PWD,rule_CONN2CONN,match1,match2)
% SPM_UPDATESPMFILEPATHS updates filepaths in SPM.mat file if data has moved to a different folder
%
% Basic syntax:
% conn_updatespmfilepaths(spmfile)  Finds any occurrences of the folder specified in SPM.swd and changes those to the current location of the SPM.mat file (or pwd if entering SPM structure manually)
%                                   Finds any occurrence of CONN specific folders (conn/rois, conn/utils/otherrois, conn/utils/surf, and conn/modules/el/root.subjects) and changes those to the location of the same folders in the current machine (fileparts(which('conn')))
%   spmfile                     : char array with full path to SPM.mat file 
%                                 (note: CONN_UPDATESPMFILEPATHS will create a copy of SPM.mat as backup_SPM.mat before making any changes to this file)
% e.g.
%   conn_updatespmfilepaths('/myanalyses/SPM.mat');
%
% Alternative syntax:
% SPM = conn_updatespmfilepaths(SPM)
%    SPM                        : SPM structure within SPM.mat file
% e.g. 
%   load(spmfile,'SPM');
%   SPM = conn_updatespmfilepaths(SPM);
%   save(spmfile,'SPM');
%
% Additional options:
% conn_updatespmfilepaths(spmfile, rule_SWD2PWD, rule_CONN2CONN, root_searchstring, root_replacestring)
%   rule_SWD2PWD                : 0/1 specifies whether CONN_UPDATESPMFILEPATHS will find any occurrences of the folder specified in SPM.swd and change those to the current location of the SPM.mat file (or pwd if entering SPM structure manually)
%   rule_CONN2CONN              : 0/1 specifies whether CONN_UPDATESPMFILEPATHS will find any occurrence of CONN specific folders and change those to the location of these folders in the current machine (fileparts(which('conn')))
%   root_*                      : explicitly defines other search/replace patterns as: [root_searchstring][body_filename] -> [root_replacestring][body_filename]
%      root_searchstring        :  search string (char array or cell array)
%      root_replacestring       :  replacement string (char array or cell array)
%      e.g. [/disk1/mydata][/subject01/session01/myfile.img] -> [/disk2/users/me/data][/subject01/session01/myfile.img]
% e.g.
%   conn_updatespmfilepaths('/myanalyses/SPM.mat', false, false, '/disk1/mydata','/disk2/users/me/data');
%


connfolder=fileparts(which('conn'));
if nargin<2||isempty(rule_SWD2PWD), rule_SWD2PWD=true; end
if nargin<3||isempty(rule_CONN2CONN), rule_CONN2CONN=true; end
if nargin<4||isempty(match1), match1={}; end
if nargin<5||isempty(match2), match2={}; end
if ischar(rule_SWD2PWD), rule_SWD2PWD=str2num(rule_SWD2PWD); end
if ischar(rule_CONN2CONN), rule_CONN2CONN=str2num(rule_CONN2CONN); end
if ischar(match1), match1=cellstr(match1); end
if ischar(match2), match2=cellstr(match2); end
assert(numel(match1)==numel(match2), 'mismatch number of search/replace pairs');
if ischar(filename) % input filename
    filename=conn_fullfile(filename);
    data=load(filename,'-mat');
    assert(isfield(data,'SPM'),'missing SPM structure in file %s',filename);
    [filepath,nill,nill]=fileparts(filename);
else % input SPM structure
    data.SPM=filename;
    filename='';
    filepath=pwd;
end
if rule_SWD2PWD&&isfield(data.SPM,'swd')&&~isempty(data.SPM.swd)&&~isempty(filepath)&&~isequal(data.SPM.swd,filepath) % adds SWD2PWD rule to match1/match2 pairs
    match1{end+1}=data.SPM.swd;
    match2{end+1}=filepath;
    data.SPM.swd=filepath;
end
changes=[];
for n1=1:numel(match1)
    changes=cat(2,changes, struct('key',match1{n1},'fullname2',match2{n1},'m1',numel(match1{n1}),'m2',numel(match2{n1})));
end
% updates the following SPM fields
update={...
    'xY.VY',...
    'xY.P',...
    'xVol.VRpv',...
    'Vbeta',...
    'VResMS',...
    'VM'};
for nupdate=1:length(update),
    str=regexp(update{nupdate},'\.','split');
    temp=[];
    try, temp=getfield(data.SPM,str{:}); end
    if ~isempty(temp),
        temp=conn_updatespmfilepaths_internal(temp);
        if ~isempty(temp), data.SPM=setfield(data.SPM,str{:},temp);
        else break;
        end
    end
end
% updates SPM.SPM.xCon(#).Vcon and SPM.SPM.xCon(#).Vspm fields
if isfield(data.SPM,'xCon')
    for ncon=1:numel(data.SPM.xCon)
        if isfield(data.SPM.xCon,'Vcon')
            temp=data.SPM.xCon(ncon).Vcon;
            temp=conn_updatespmfilepaths_internal(temp);
            if ~isempty(temp), data.SPM.xCon(ncon).Vcon=temp;
            else break;
            end
        end
        if isfield(data.SPM.xCon,'Vspm')
            temp=data.SPM.xCon(ncon).Vspm;
            temp=conn_updatespmfilepaths_internal(temp);
            if ~isempty(temp), data.SPM.xCon(ncon).Vspm=temp;
            else break;
            end
        end
    end
end

if isempty(filename)
    SPM=data.SPM;
else
    if ~conn_existfile(conn_prepend('backup_',filename)), conn_fileutils('copyfile',filename,conn_prepend('backup_',filename)); end
    conn_savematfile(filename,'-struct',data,'-v7.3');
end


    function var=conn_updatespmfilepaths_internal(var)
        if ~isfield(var,'fname')&&~ischar(var), return; end % only structure or char array
        if isfield(var,'fname')&&numel(var)>1, % fname structure array
            for nvar=1:numel(var),
                newvar=conn_updatespmfilepaths_internal(var(nvar));
                if ~isempty(newvar), var(nvar)=newvar; elseif ~isempty(var(nvar)), var=[]; break; end
            end
        else
            if isfield(var,'fname'), varfilename=var.fname;
            else varfilename=var;
            end
            if ~isempty(varfilename),
                varfilename=regexprep(cellstr(varfilename),'(^\s+)|(\s+$)','');
                if ~isempty(changes), 
                    varfilename=conn_updatefilepaths_change(varfilename,changes); 
                end
                if rule_CONN2CONN&&~isempty(connfolder)
                    varfilename=regexprep(varfilename, '^.*\/conn(\/modules\/el\/root\.subjects\/.+)$', [connfolder,'$1']);
                    varfilename=regexprep(varfilename, '^.*\/conn(\/rois\/[^\/]+)$', [connfolder,'$1']);
                    varfilename=regexprep(varfilename, '^.*\/conn(\/utils\/otherrois\/[^\/]+)$', [connfolder,'$1']);
                    varfilename=regexprep(varfilename, '^.*\/conn(\/utils\/surf\/[^\/]+)$', [connfolder,'$1']);
                end
                if isfield(var,'fname'), var.fname=char(varfilename);
                else var=char(varfilename);
                end
            end
        end
    end
end

function fname=conn_updatefilepaths_change(fname,changes)
waschar=false;
if ischar(fname), waschar=true; fname=cellstr(fname); end
if ~isempty(changes)
    for nch=numel(changes):-1:1
        ok=strncmp(changes(nch).key,fname,changes(nch).m1);
        if any(ok),
            fname(ok)=cellfun(@(x)[changes(nch).fullname2(1:changes(nch).m2),x(changes(nch).m1+1:end)],fname(ok),'uni',0);
        end
    end
end
if waschar, fname=char(fname); end
end




