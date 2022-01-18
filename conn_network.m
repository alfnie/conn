function results=conn_network(filename,rois,measures,normalization,thr)
% CONN_NETWORK network analyses gui
%
% conn_network(filename,rois,measures,normalization,thr);
%  filename:        resultsROI_Condition*.mat file (found in /results/firstlevel/)
%  rois:            indexes of ROIs to compute ROI-to-ROI connectivity matrix
%  measures:        1: Efficiency/cost measures; 2: PathDistance/Clustering measures
%  normalization:   0: uses raw (Fisher-transformed) correlation coefficients; 1: normalizes across subjects (transform to z-scores); 2: normalizes using cost
%  thr:             threshold value to compute adjacency matrix
%
% conn_network (without arguments) launches a gui to select the above parameters 
%

posstr=1;

if nargin<1||isempty(filename),
    global CONN_x;
    if ~isempty(CONN_x)&&isfield(CONN_x,'folders')&&isfield(CONN_x.folders,'firstlevel'),pathname=fullfile(CONN_x.folders.firstlevel,CONN_x.Analyses(1).name);
    else pathname=pwd; end
    [filename, pathname] = conn_fileutils('uigetfile','.mat', 'Select a condition file',fullfile(pathname,'resultsROI_Condition*.mat'));
    filename=fullfile(pathname,filename);
end

% loads connectivity values
if iscell(filename), Z=filename{1};names=filename{2};if numel(filename)>2,filename=filename{3};else filename=[]; end
else conn_loadmatfile(filename,'Z','names'); end
[nROI,nROI2,nSubjects]=size(Z);
if nargin<2||isempty(rois),
    [rois,ok]=listdlg('PromptString','Select ROIs:',...
        'SelectionMode','multiple',...
        'ListString',char(names),...
        'initialvalue',1:nROI);
    if ~ok,return;end
end
names={names{rois}};Z=Z(rois,rois,:);nROI=length(rois);

% chooses set of measures
if nargin<3||isempty(measures),
    measures=1;
    if isnumeric(posstr)&&isempty(findobj(0,'tag','Interactive')), spm('CreateIntWin'); end;
    measures=spm_input('Select measures',posstr,'m','Global efficienty/Local efficiency/Cost|Path distance/Clustering coefficient/Degree',[],measures);posstr='+1';
end
% normalizes across subjects
if nargin<4||isempty(normalization),
    normalization=3;
    if isnumeric(posstr)&&isempty(findobj(0,'tag','Interactive')), spm('CreateIntWin'); end;
    normalization=spm_input('Inter-subjects normalization',posstr,'m','none|z-scores|cost',[],normalization);posstr='+1';
    normalization=normalization-1;
end
ZZ=zeros(size(Z));
for n=1:size(Z,3),
    z=Z(:,:,n);
    z(1:size(z,1)+1:end)=nan; 
    idx=find(~isnan(z));
    switch(normalization)
        case 0,
        case 1,z=(z-mean(z(idx)))/std(z(idx));
        case 2,[st,sidx]=sort(z(idx)); z(idx(sidx))=((1:numel(idx))-.5)'/numel(idx);
    end
    ZZ(:,:,n)=z;
end
N=size(ZZ,1);

%explores variable thresholding
if nargin<5||isempty(thr),    
    if measures==1,
        ht=conn_waitbar(0,'Computing variable threshold range. Please wait...');
        switch(normalization)
            case 0,THR=0:.025:.75;
            case 1,THR=0:.1:4;
            case 2,THR=0:.01:.50;
        end
        %THR=linspace(0,.99*min(max(max(ZZ,[],1),[],2),[],3),20);
        net_d=zeros(size(ZZ,3),numel(THR));net_c=net_d;net_k=net_d;
        net_d_rand=zeros(1,numel(THR));net_c_rand=net_d_rand;net_k_rand=net_d_rand;
        net_d_latt=zeros(1,numel(THR));net_c_latt=net_d_latt;net_k_latt=net_d_latt;
%         idx=find(triu(ones(N,N)-eye(N)));
        for nthr=1:numel(THR),
            for n=1:size(ZZ,3),
                thr=THR(nthr);
                if normalization==2, C=(ZZ(:,:,n)>=1-thr); else C=(ZZ(:,:,n)>=thr); end
                %                 t=ZZ(:,:,n);st=sort(t(idx),'descend');
                %                 C=(t>=st(max(1,round(thr*numel(st)))));
                [nill,nill,nill,d,c,k]=conn_network_efficiency(C);
                net_c(n,nthr)=c;
                net_d(n,nthr)=d;
                net_k(n,nthr)=k;
            end
            thr=THR(nthr);
            [nill,nill,nill,d,c,k]=conn_network_efficiency([N,net_k(:,nthr)'],2e1);
            net_c_rand(nthr)=c;
            net_d_rand(nthr)=d;
            net_k_rand(nthr)=k;
            [nill,nill,nill,d,c,k]=conn_network_efficiency([N,net_k(:,nthr)'],-2e1);
            net_c_latt(nthr)=c;
            net_d_latt(nthr)=d;
            net_k_latt(nthr)=k;
            conn_waitbar(nthr/numel(THR),ht);
        end
        conn_waitbar('close',ht);
        %plots network measures
        x={net_d,net_c};%y={net_k};
        x0={{net_d_rand,net_d_latt},{net_c_rand,net_c_latt}};%y0={net_k_rand};
        switch(normalization)
            case 0,y={THR}; y0={THR};ylabels={{'\fontsize{16}correlation (Fisher-transformed)'}};
            case 1,y={THR}; y0={THR};ylabels={{'\fontsize{16}z-score'}};
            case 2,y={net_k}; y0={net_k_rand};ylabels={{'\fontsize{16}Cost (K)'}};
        end
        
        xlabels={{'\fontsize{16}\color{black}Global efficiency'},{'\fontsize{16}\color{black}Local efficiency'}};
        figure('name','Network theory: explore variable threshold range','numbertitle','off','units','norm','position',[.3 .3 .6 .6]);clf;
        for n1=1:numel(x),
            thr=mean(y{1},1);
            mnet_d=mean(x{n1},1);
            snet_d=std(x{n1},0,1);
            n=size(x{n1},1);
            m=size(x{n1},2);
            alpha=.05;
            subplot(numel(x),2,2*n1-1);
            hold on;
            h=patch([thr,fliplr(thr)],[mnet_d+2*snet_d,fliplr(mnet_d-2*snet_d)],'k');
            %h=patch([thr,fliplr(thr)],[mnet_d+spm_invTcdf(1-alpha/m/2,n-1)*snet_d/sqrt(n),fliplr(mnet_d-spm_invTcdf(1-alpha/m/2,n-1)*snet_d/sqrt(n))],'k');
            set(h,'facecolor',.75*[1,1,1],'edgecolor','none');
            h1=plot(thr,mnet_d,'k-','linewidth',2);
            h2=[];linetypes={'--',':','-.','.-'}; linecolors={.5*[1,1,1],.5*[1,1,1],.5*[1,1,1],.5*[1,1,1]}; for n2=1:numel(x0{n1}),
                thr0=mean(y0{1},1);
                mnet_d0=mean(x0{n1}{n2},1);
                h2(n2)=plot(thr0,mnet_d0,[linetypes{n2}],'linewidth',2,'color',linecolors{n2});
            end
            hold off;
            set(gcf,'color','w');set(gca,'fontsize',16,'xcolor',.75*[1,1,1],'ycolor',.75*[1,1,1]);
            ylabel(xlabels{n1});xlabel(ylabels{1});
            if n1==1, legend([h1,h2],{'Data','Random graph','Lattice'}); end
        end
        h=subplot(1,2,2); set(h,'position',get(h,'position')*[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 .25 0 .5]);
        mnet_d_1=mnet_d-mean(x0{1}{2},1); %title('Difference with lattice')
        mnet_d_2=mnet_d-mean(x0{n1}{1},1); %title('Difference with random');
        mnet_d=mnet_d_1+mnet_d_2; title('GE_{data}-GE_{lattice} + LE_{data}-LE_{random}');
        snet_d_1=std(x{1},0,1);
        snet_d_2=std(x{2},0,1);
        snet_d=sqrt(snet_d_1.^2+snet_d_2.^2);
        hold on;
        %             alllinesY=inf;
        %             for n2=1:numel(x0{n1}),
        %                 mnet_d0=mean(x0{n1}{n2},1);
        %                 alllinesY=min(alllinesY,mnet_d0);
        %             end
        %             mnet_d=mnet_d-alllinesY;
        h=patch([thr,fliplr(thr)],[mnet_d+2*snet_d,fliplr(mnet_d-2*snet_d)],'k');
        set(h,'facecolor',.75*[1,1,1],'edgecolor','none');
        h1=plot(thr,mnet_d,'k-','linewidth',2);
        hold off;
        set(gcf,'color','w');set(gca,'fontsize',16,'xcolor',.75*[1,1,1],'ycolor',.75*[1,1,1]);
        %ylabel(xlabels{n1});
        xlabel(ylabels{1});
        
        %clf; plot(mean(x{1},1)-x0{1}{2},mean(x{2},1)-x0{2}{1},'.-'); text(mean(x{1},1)-x0{1}{2},mean(x{2},1)-x0{2}{1},arrayfun(@(x)num2str(x,'%.2f'),thr,'uni',0)); hold off;
        
        if isempty(filename), return; end
        if isnumeric(posstr)&&isempty(findobj(0,'tag','Interactive')), spm('CreateIntWin'); end;
        switch(normalization)
            case 0,thr=spm_input('threshold (Z-score)',posstr,'r',mat2str(.5),[1,1]);posstr='+1';
            case 1,thr=spm_input('threshold (z-score)',posstr,'r',mat2str(1.5),[1,1]);posstr='+1';
            case 2,thr=spm_input('threshold (cost)',posstr,'r',mat2str(.15),[1,1]);posstr='+1';
        end
        if isempty(thr),return;end
    else
    end
end

h=conn_msgbox('computing network measures, please wait...','conn_network');
if measures==1
    measure_names=conn_network_measures;
    nmeasures=numel(measure_names);
    NET_values=zeros(size(ZZ,3),nROI,nmeasures);
    net_values=zeros(size(ZZ,3),nmeasures);
else
    net_d=zeros(size(ZZ,3),1);net_c=net_d;net_k=net_d;
    net_rc=net_c;net_rd=net_d;
    NET_d=zeros(size(ZZ,3),nROI);NET_c=NET_d;NET_k=NET_d;
end
for nthr=1:numel(thr),
    %computes adjacency matrix
    if measures==1,
        if normalization==2, ZZ_thr=(ZZ>=1-thr(nthr)); else ZZ_thr=(ZZ>=thr(nthr)); end
%         for n=1:size(ZZ,3),
%             if normalization==2, ZZ_thr(:,:,n)=(ZZ(:,:,n)>=1-thr(nthr)); else ZZ_thr(:,:,n)=ZZ(:,:,n)>=thr(nthr); end
%         end
    else
        ZZ_thr=ZZ>thr(nthr);
    end
    
    %compute network-level measures
    if measures==1,
        for n=1:size(ZZ,3),
            C=ZZ_thr(:,:,n);
            %[Eg,El,K,eg,el,k]=conn_network_efficiency(C);
            [nill,measure_values]=conn_network_measures(C);
            for nmeasure=1:nmeasures
                NET_values(n,:,nmeasure)=NET_values(n,:,nmeasure)+measure_values{nmeasure}';
                net_values(n,nmeasure)=net_values(n,nmeasure)+mean(measure_values{nmeasure}(~isnan(measure_values{nmeasure})));
            end
            %net_c(n)=net_c(n)+el;
            %net_d(n)=net_d(n)+eg;
            %net_k(n)=net_k(n)+k;
            %NET_c(n,:)=NET_c(n,:)+El';
            %NET_d(n,:)=NET_d(n,:)+Eg';
            %NET_k(n,:)=NET_k(n,:)+K';
        end
    else
        for n=1:size(ZZ,3),
            C=ZZ_thr(:,:,n);
            [c,ac,nac]=conn_network_clustering(C,5e2);
            [d,ad,nad]=conn_network_mindist(C,5e2);
            [k,ak]=conn_network_degree(C);
            net_rc(n)=net_rc(n)+ac;
            net_rd(n)=net_rd(n)+ad;
            net_c(n)=net_c(n)+nac;
            net_d(n)=net_d(n)+nad;
            net_k(n)=net_k(n)+ak;
            NET_c(n,:)=NET_c(n,:)+c';
            %     NET_d(n,:)=1./(sum((d>0).*(1./max(eps,d)),1)./sum(d>0,1));
            NET_d(n,:)=NET_d(n,:)+sum((d>0&~isinf(d)).*min(1e10,d),1)./sum(d>0&~isinf(d),1);
            NET_k(n,:)=NET_k(n,:)+k(:)';
        end
    end
end
if measures==1
    NET_values=NET_values/numel(thr);
    net_values=net_values/numel(thr);
else
    net_d=net_d/numel(thr);net_c=net_c/numel(thr);net_k=net_k/numel(thr);
    net_rc=net_rc/numel(thr);net_rd=net_rd/numel(thr);
    NET_d=NET_d/numel(thr);NET_c=NET_c/numel(thr);NET_k=NET_k/numel(thr);
end
if numel(thr)>1,
    %computes average adjacency matrix
    ZZ_thr=0;
    for nthr=1:numel(thr),
        if measures==1,
            if normalization==2, ZZ_thr=ZZ_thr+(ZZ>=1-thr(nthr)); else ZZ_thr=ZZ_thr+(ZZ>=thr(nthr)); end
        else
            ZZ_thr=ZZ_thr+(ZZ>thr(nthr));
        end
    end
    ZZ_thr=ZZ_thr/numel(thr);
end
if ishandle(h), close(h); end

if ~isempty(filename),
    [pathname,filename,fileext]=fileparts(filename);
    if strcmp(fileext,'.mat') % create adjacency matrices .mat file
        A=ZZ_thr; ROInames=names;
        conn_savematfile(fullfile(pathname,[filename,fileext]),'A','ROInames');
        fprintf(1,'Adjacency matrices saved as: %s\n',fullfile(pathname,[filename,fileext]));
    elseif strcmp(fileext,'.dl') % creates .dl UCINET format files
        for n=1:size(ZZ,3),
            filename_out=fullfile(pathname,[filename,'_Subject',num2str(n,'%03d'),'.dl']);
            fh=fopen(filename_out,'wt');
            if isequal(fh,-1),filename_out=fullfile(pwd,[filename,'_Sample',num2str(n,'%03d'),'.dl']);fh=fopen(filename_out,'wt');end
            fprintf(fh,'dl format=fullmatrix,n=%d\n',nROI);
            fprintf(fh,'labels:\n');
            for n1=1:nROI,name=names{n1};name(name==' ')=[];name=name(1:min(numel(name),18));fprintf(fh,'%s',name);if n1==nROI,fprintf(fh,'\n');else,fprintf(fh,',');end;end
            fprintf(fh,'data:\n');
            for n1=1:nROI,fprintf(fh,'%d ',ZZ_thr(n1,:,n)>.5);fprintf(fh,'\n');end
            fclose(fh);
            fprintf(1,'Output .dl file created : %s\n',filename_out);
        end
    else
        % creates network-level measures file
        filename_out=fullfile(pathname,[filename,'.csv']);
        %fh=fopen(filename_out,'wt');
        %if isequal(fh,-1),filename_out=fullfile(pwd,[filename,'.csv']);fh=fopen(filename_out,'wt');end
        fh={};
        if measures==1,
            tnames=names;for n1=1:numel(names),tnames{n1}(tnames{n1}==',')=';';end
            fh{end+1}=sprintf([',Network',repmat(',',[1,nmeasures])]);for n1=1:nROI,fh{end+1}=sprintf([',%s',repmat(',',[1,nmeasures-1])],tnames{n1});end;fh{end+1}=sprintf('\n');
            fh{end+1}=sprintf('Subject/Sample #,');
            for nmeasure=1:nmeasures,fh{end+1}=sprintf('%s,',measure_names{nmeasure});end
            for n1=1:nROI,
                for nmeasure=1:nmeasures,fh{end+1}=sprintf(',%s',measure_names{nmeasure});end
            end;
            fh{end+1}=sprintf('\n');
            for n=1:size(ZZ,3),
                fh{end+1}=sprintf('%d,',n);
                for nmeasure=1:nmeasures,fh{end+1}=sprintf('%f,',net_values(n,nmeasure));end
                for n1=1:nROI,
                    for nmeasure=1:nmeasures,fh{end+1}=sprintf(',%f',NET_values(n,n1,nmeasure));end
                end
                fh{end+1}=sprintf('\n');
            end
        else
            fh{end+1}=sprintf(',Global,,,,,');for n1=1:nROI,fh{end+1}=sprintf(',%s,,',names{n1});end;fh{end+1}=sprintf('\n');
            fh{end+1}=sprintf('Subject/Sample #,Average path distance,Normalized average path distance,Average clustering coefficient,Normalized average clustering coefficient,Average degree,');for n1=1:nROI,fh{end+1}=sprintf(',Average path distance,Average clustering coefficient,degree');end;fh{end+1}=sprintf('\n');
            for n=1:size(ZZ,3),
                fh{end+1}=sprintf('%d,%f,%f,%f,%f,%f,',n,net_rd(n),net_d(n),net_rc(n),net_c(n),net_k(n));
                for n1=1:nROI,fh{end+1}=sprintf(',%f,%f,%f',NET_d(n,n1),NET_c(n,n1),NET_k(n,n1));end;fh{end+1}=sprintf('\n');
            end
        end
        %fclose(fh);
        try, conn_fileutils('filewrite_raw',filename_out,fh);
        catch, filename_out=fullfile(pwd,[filename,'.csv']); conn_fileutils('filewrite_raw',filename_out,fh);
        end
        fprintf(1,'Network measures saved as: %s\n',filename_out);
    end
end
if measures==1,
    %results=struct('type',measures,'Z',ZZ,'Z_thr',ZZ_thr,'thr',thr,'rois',{names},'GlobalEfficiency',net_d,'LocalEfficiency',net_c,'Cost',net_k,'GlobalEfficiency_roi',NET_d,'LocalEfficiency_roi',NET_c,'Cost_roi',NET_k);
    for nmeasure=1:nmeasures
        structmeasures.([measure_names{nmeasure},'_roi'])=NET_values(:,:,nmeasure);
        structmeasures.(measure_names{nmeasure})=net_values(:,nmeasure);
    end
    results=struct('type',measures,'Z',ZZ,'Z_thr',ZZ_thr,'thr',thr,'rois',{names},'measures',structmeasures,'measure_names',{measure_names});
else
    results=struct('type',measures,'Z',ZZ,'Z_thr',ZZ_thr,'thr',thr,'rois',{names},'pathlength',net_rd,'clustering',net_rc,'norm_pathlength',net_d,'norm_clustering',net_c,'degree',net_k,'pathlength_roi',NET_d,'clustering_roi',NET_c,'degree_roi',NET_k);
    results.norm_pathlength_roi=diag(results.norm_pathlength./results.pathlength)*results.pathlength_roi;
    results.norm_clustering_roi=diag(results.norm_clustering./results.clustering)*results.clustering_roi;
end
if ~isempty(filename)&&strcmp(fileext,'.csv')
    filename_out=fullfile(pathname,[filename,'.networkmeasures.mat']);
    conn_savematfile(filename_out,'results');
    fprintf(1,'Network measures saved as: %s\n',filename_out);
end
