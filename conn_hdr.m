function varargout=conn_hdr(filename)
% CONN_HDR Displays setup information in conn*.mat file
%

if nargin<1 || isempty(filename),
    [filename,filepath]=uigetfile('conn*.mat');
    if ~ischar(filename), return; end
    filename=fullfile(filepath,filename);
    load(filename,'CONN_x');
elseif isstruct(filename),
    CONN_x=filename;
    filename='';
else
    load(filename,'CONN_x');
end
txt={};
nl1covariates=length(CONN_x.Setup.l1covariates.names)-1;
nl2covariates=length(CONN_x.Setup.l2covariates.names)-1;
nconditions=length(CONN_x.Setup.conditions.names)-1;
nrois=length(CONN_x.Setup.rois.names)-1;
txt{end+1}=sprintf('Project %s',filename);
txt{end+1}=sprintf('%d subjects',CONN_x.Setup.nsubjects);
txt{end+1}=sprintf('%d ROIs',nrois);
txt{end+1}=sprintf('%d conditions',nconditions);
txt{end+1}=sprintf('%d first-level covariates',nl1covariates);
txt{end+1}=sprintf('%d second-level covariates',nl2covariates);
for nsub=1:CONN_x.Setup.nsubjects,
    nsess=CONN_x.Setup.nsessions(min(length(CONN_x.Setup.nsessions),nsub));
    txt{end+1}=sprintf('------------------------------------');
    txt{end+1}=sprintf('Subject %d',nsub);
    for nl2covariate=1:nl2covariates,
        value=CONN_x.Setup.l2covariates.values{nsub}{nl2covariate};
        names=CONN_x.Setup.l2covariates.names{nl2covariate};
        txt{end+1}=sprintf('Second-level covariate "%s": %f',names,value);
    end
    txt{end+1}=sprintf('%d sessions',nsess);
    for nses=1:nsess,
        nscans=CONN_x.Setup.nscans{nsub}{nses};
        txt{end+1}=sprintf('Session %d (%d scans)',nses,nscans);
        txt{end+1}=sprintf('Anatomical volume: %s %s',CONN_x.Setup.structural{nsub}{nses}{1},CONN_x.Setup.structural{nsub}{nses}{2}{end});
        if size(CONN_x.Setup.functional{nsub}{nses}{1},1)>1,txt{end+1}=sprintf('Functional volumes : %s ... %s %s',CONN_x.Setup.functional{nsub}{nses}{1}(1,:),CONN_x.Setup.functional{nsub}{nses}{1}(end,:),CONN_x.Setup.functional{nsub}{nses}{2}{end});
        else, txt{end+1}=sprintf('Functional volumes: %s %s',CONN_x.Setup.functional{nsub}{nses}{1},CONN_x.Setup.functional{nsub}{nses}{2}{end}); end
        for nroi=1:nrois,
            names=CONN_x.Setup.rois.names{nroi};
            txt{end+1}=sprintf('ROI "%s": %s',names,CONN_x.Setup.rois.files{nsub}{nroi}{nses}{1});
        end
        for nl1covariate=1:nl1covariates,
            filename=CONN_x.Setup.l1covariates.files{nsub}{nl1covariate}{nses}{1};
            names=CONN_x.Setup.l1covariates.names{nl1covariate};
            if strcmp(filename,'[raw values]'),
                data=CONN_x.Setup.l1covariates.files{nsub}{nl1covariate}{nses}{3};
                txt{end+1}=sprintf('First-level covariate "%s" (%d scans, %d dimensions)',names,size(data,1),size(data,2));
            else
                txt{end+1}=sprintf('First-level covariate "%s": %s',names,filename);
            end
        end
    end
end
if ~nargout,conn_disp(char(txt));else,varargout{1}=txt;end

