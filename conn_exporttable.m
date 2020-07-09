function conn_exporttable(filename,varargin)
% conn_exporttable
% exports ROI-level results table

global CONN_h;
if ~nargin||~ischar(filename),
    [filename,filepath]=uiputfile({'*.txt','*.txt (text file)';'*.csv','*.csv (comma-separated table)';'*.mat','*.mat (matlab file)'},'Save stats as');    
    if ~ischar(filename), return; end
    filename=fullfile(filepath,filename);
end
[filepath,filename,fileext]=fileparts(filename);
%P=conn_fdr(CONN_h.menus.m_results.roiresults.p,2);
switch(fileext),
    case '.csv',
        fh=fopen(fullfile(filepath,[filename,fileext]),'wt');
        if size(CONN_h.menus.m_results.roiresults.dof,2)>1
            fprintf(fh,'%s,%s,%s,%s,%s\n','Targets','beta',[CONN_h.menus.m_results.roiresults.statsname,'(',num2str(CONN_h.menus.m_results.roiresults.dof(1)),',',num2str(CONN_h.menus.m_results.roiresults.dof(2)),')'],'p-unc','p-FDR');
        else
            fprintf(fh,'%s,%s,%s,%s,%s\n','Targets','beta',[CONN_h.menus.m_results.roiresults.statsname,'(',num2str(CONN_h.menus.m_results.roiresults.dof(1)),')'],'p-unc','p-FDR');
        end
        for n1=1:numel(CONN_h.menus.m_results.roiresults.idx),
            n2=CONN_h.menus.m_results.roiresults.idx(n1);
            tmp=CONN_h.menus.m_results.roiresults.names2{n2};if length(tmp)>36,tmp=[tmp(1:31),'*',tmp(end-3:end)]; end;
            fprintf(fh,'%s,%10.6f,%10.6f,%10.6f,%10.6f\n',tmp,CONN_h.menus.m_results.roiresults.h(n2),CONN_h.menus.m_results.roiresults.F(n2),CONN_h.menus.m_results.roiresults.p(n2),CONN_h.menus.m_results.roiresults.P(n2));%P(n2));
        end
        fclose(fh);
    case '.txt',
        fh=fopen(fullfile(filepath,[filename,fileext]),'wt');
        if size(CONN_h.menus.m_results.roiresults.dof,2)>1
            fprintf(fh,'%-36s%10s%10s%12s%12s\n','Targets','beta',[CONN_h.menus.m_results.roiresults.statsname,'(',num2str(CONN_h.menus.m_results.roiresults.dof(1)),',',num2str(CONN_h.menus.m_results.roiresults.dof(2)),')'],'p-unc','p-FDR');
        else
            fprintf(fh,'%-36s%10s%10s%12s%12s\n','Targets','beta',[CONN_h.menus.m_results.roiresults.statsname,'(',num2str(CONN_h.menus.m_results.roiresults.dof(1)),')'],'p-unc','p-FDR');
        end
        for n1=1:numel(CONN_h.menus.m_results.roiresults.idx),
            n2=CONN_h.menus.m_results.roiresults.idx(n1);
            tmp=CONN_h.menus.m_results.roiresults.names2{n2};if length(tmp)>36,tmp=[tmp(1:31),'*',tmp(end-3:end)]; end;
            fprintf(fh,'%-36s%10.2f%10.2f%12f%12f\n',tmp,CONN_h.menus.m_results.roiresults.h(n2),CONN_h.menus.m_results.roiresults.F(n2),CONN_h.menus.m_results.roiresults.p(n2),CONN_h.menus.m_results.roiresults.P(n2));
        end
        fclose(fh);
    case '.mat',
        roiresults.names={CONN_h.menus.m_results.roiresults.names2{CONN_h.menus.m_results.roiresults.idx}};
        roiresults.beta=CONN_h.menus.m_results.roiresults.h(CONN_h.menus.m_results.roiresults.idx);
        roiresults.T=CONN_h.menus.m_results.roiresults.F(CONN_h.menus.m_results.roiresults.idx);
        roiresults.p_unc=CONN_h.menus.m_results.roiresults.p(CONN_h.menus.m_results.roiresults.idx);
        roiresults.p_fdr=CONN_h.menus.m_results.roiresults.P(CONN_h.menus.m_results.roiresults.idx);
        roiresults.Y=CONN_h.menus.m_results.roiresults.y(:,CONN_h.menus.m_results.roiresults.idx);
        roiresults.X=CONN_h.menus.m_results.roiresults.xX.X;
        save(fullfile(filepath,[filename,fileext]),'roiresults');
    case '.txt',
end
