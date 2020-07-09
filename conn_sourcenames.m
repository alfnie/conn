function [iroi,isnew]=conn_sourcenames(name,flag)
global CONN_x;
if ~isfield(CONN_x,'Analysis'), CONN_x.Analysis=1; end
ianalysis=CONN_x.Analysis;
if nargin<2, flag='-'; end
if ~isfield(CONN_x.Analyses(ianalysis),'sourcenames'),
    CONN_x.Analyses(ianalysis).sourcenames={};
end
idx=strmatch(name,CONN_x.Analyses(ianalysis).sourcenames,'exact');
if ~isempty(idx),
    iroi=idx(1);
    isnew=0;
else,
    iroi=length(CONN_x.Analyses(ianalysis).sourcenames)+1;
    if flag=='+',
        CONN_x.Analyses(ianalysis).sourcenames{iroi}=name;
        filesourcenames=fullfile(CONN_x.folders.firstlevel,CONN_x.Analyses(ianalysis).name,'_list_sources.mat');
        filesourcenames=conn_projectmanager('projectfile',filesourcenames,CONN_x.pobj,'.mat');
        sourcenames=CONN_x.Analyses(ianalysis).sourcenames;
        save(filesourcenames,'sourcenames');
    end
    isnew=1;
end
end
