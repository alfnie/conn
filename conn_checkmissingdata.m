function MD = conn_checkmissingdata(state,nconditions, nsources, force)
global CONN_x;

if nargin<1||isempty(state), state=2; end
if nargin<2||isempty(nconditions), nconditions=[]; end
if nargin<3||isempty(nsources), nsources=[]; end
if nargin<4||isempty(force), force=false; end
if isequal(state,'missinglabel'), MD=cellfun(@conn_checkmissingdata_ismissing,nconditions); return; end % internal use only

if state==2
    filepathresults1=fullfile(CONN_x.folders.firstlevel,CONN_x.Analyses(CONN_x.Analysis).name);
else
    filepathresults1=fullfile(CONN_x.folders.firstlevel_vv,CONN_x.vvAnalyses(CONN_x.vvAnalysis).name);
end
if nargin<2||isempty(nconditions)
    nconditions=CONN_x.Results.xX.nconditions;
end
icondition=[];isnewcondition=[];for ncondition=nconditions(:)',[icondition(ncondition),isnewcondition(ncondition)]=conn_conditionnames(CONN_x.Setup.conditions.names{ncondition}); end
switch state
    case 2,
        if nargin<3||isempty(nsources), nsources=CONN_x.Results.xX.nsources; end
        sources=CONN_x.Analyses(CONN_x.Analysis).sources;
    case 3,
        if nargin<3||isempty(nsources), nsources=CONN_x.Results.xX.nmeasures; end
        sources=CONN_x.vvAnalyses(CONN_x.vvAnalysis).measures;
end

MD=true(CONN_x.Setup.nsubjects,1);
list_tfilename={};
list_MD=[];
for nsub=1:CONN_x.Setup.nsubjects,
    for n0=1:numel(nconditions)
        for n1=1:min(1,numel(nsources))
            nroi=nsources(n1);
            roiname=sources{nroi};
            ncondition=nconditions(n0);
            tfilename='';
            switch state
                case 2, % seed-to-voxel
                    [iroi,isnew]=conn_sourcenames(roiname,'-');
                    if ~isnew, tfilename=fullfile(filepathresults1,['BETA_Subject',num2str(nsub,CONN_x.opt.fmt1),'_Condition',num2str(icondition(ncondition),'%03d'),'_Source',num2str(iroi,'%03d'),'.nii']); end
                    if ~isempty(tfilename)
                        if conn_projectmanager('inserver'), list_tfilename{end+1}=tfilename; list_MD(end+1)=nsub;
                        elseif conn_checkmissingdata_ismissing(tfilename), MD(nsub)=false; 
                        end
                    end
                case 3, % voxel-to-voxel
                    [iroi,isnew,ncomp]=conn_v2v('match_extended',roiname);
                    if ~isnew, tfilename=fullfile(filepathresults1,['BETA_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'_Measure',num2str(iroi,'%03d'),'_Component',num2str(ncomp,'%03d'),'.nii']); end
                    if ~isempty(tfilename)
                        if conn_projectmanager('inserver'), list_tfilename{end+1}=tfilename; list_MD(end+1)=nsub;
                        elseif conn_checkmissingdata_ismissing(tfilename), MD(nsub)=false; 
                        end
                    end
            end
        end
    end
end
if conn_projectmanager('inserver'), 
    ismissing=conn_server('run',mfilename,'missinglabel',conn_server('util_localfile',list_tfilename));
    for n=1:numel(ismissing), if ismissing(n), MD(list_MD(n))=false; end; end
end
end

function ok=conn_checkmissingdata_ismissing(tfilename)
ok=false;
try,
    try, a=nifti(tfilename);
    catch, a=spm_vol(tfilename);
    end
    try, if isequal(a.descrip,'CONNlabel:MissingData'), ok=true; end; end
end
end
