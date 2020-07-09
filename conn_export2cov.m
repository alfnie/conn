function [ok,tfilename]=conn_export2cov(nl2covariates,tfilename)
% conn_export2cov exports 2nd-level covariate to file
% [ok,tfilename]=conn_export2cov(nl2covariates,tfilename)
%

global CONN_x;

ok=false;
if nargin<1||isempty(nl2covariates)
    nl2covariates=listdlg('liststring',CONN_x.Setup.l2covariates.names(1:end-1),'selectionmode','multiple','initialvalue',[],'promptstring',{'Select covariates to save'},'ListSize',[500 200]);
    if isempty(nl2covariates), return; end
end
if nargin<2||isempty(tfilename)
    [tfilename,tpathname]=uiputfile({'*.mat','MAT-files (*.mat)'; '*.txt','text files (*.txt)'; '*.csv','CSV-files (*.csv)'; '*',  'All Files (*)'},'Output data to file:');
    if ~ischar(tfilename)||isempty(tfilename), return; end
    tfilename=fullfile(tpathname,tfilename);
end
[nill,nill,tfileext]=fileparts(tfilename);
nl2covariates(nl2covariates==numel(CONN_x.Setup.l2covariates.names))=[];
tt=[];
for il2covariate=1:numel(nl2covariates),
    if length(CONN_x.Setup.l2covariates.descrip)<nl2covariates(il2covariate), CONN_x.Setup.l2covariates.descrip{nl2covariates(il2covariate)}=''; end
    t=[];
    for nsub=1:CONN_x.Setup.nsubjects,
        if length(CONN_x.Setup.l2covariates.values)<nsub, CONN_x.Setup.l2covariates.values{nsub}={}; end
        if length(CONN_x.Setup.l2covariates.values{nsub})<nl2covariates(il2covariate), CONN_x.Setup.l2covariates.values{nsub}{nl2covariates(il2covariate)}=nan; end
        t=cat(2,t,CONN_x.Setup.l2covariates.values{nsub}{nl2covariates(il2covariate)});
    end
    tt=cat(1,tt,t);
end
names=CONN_x.Setup.l2covariates.names(nl2covariates);
descrip=CONN_x.Setup.l2covariates.descrip(nl2covariates);
switch(tfileext)
    case '.mat'
        data=tt.';
        save(tfilename,'data','names','descrip');
    otherwise,
        if strcmp(tfileext,'.txt'), names=regexprep(names,'\s','');
        else                        names=regexprep(names,'\,','');
        end
        fh=fopen(tfilename,'wt');
        for n1=1:numel(names),
            if isempty(names{n1}), names{n1}='-'; end
            fprintf(fh,'%s',names{n1});
            if n1<numel(names)&&strcmp(tfileext,'.csv'), fprintf(fh,','); elseif n1<numel(names), fprintf(fh,' '); else fprintf(fh,'\n'); end
        end
        for n2=1:size(tt,2),
            for n1=1:size(tt,1),
                fprintf(fh,'%f',tt(n1,n2));
                if n1<size(tt,1)&&strcmp(tfileext,'.csv'), fprintf(fh,','); elseif n1<size(tt,1), fprintf(fh,' '); else fprintf(fh,'\n'); end
            end
        end
        fclose(fh);
end
ok=true;
end
