function conn_qaplots_con(fileSPM,qafolder,fname,dpires)
% qa plot for SPM contrasts

DOORTHOG=true;

[pwd1,nill]=fileparts(fileSPM);
[nill,pwd1name]=fileparts(pwd1);
SPM=struct; conn_loadmatfile(fileSPM,'SPM');
if nargin<2||isempty(qafolder), qafolder=pwd1; end
if nargin<3||isempty(fname), fname='QA_SPM_contrasts.jpg'; end
if nargin<4||isempty(dpires), dpires='-r150'; end

if (isfield(SPM,'xCon')&&~isempty(SPM.xCon)) || (isfield(SPM,'xConOriginal')&&~isempty(SPM.xConOriginal))
    sX=spm_sp('Set',SPM.xX.X);
    if sX.rk>0, opp=sX.v(:,[1:sX.rk])*sX.v(:,[1:sX.rk])';
    else opp=zeros( size(sX.X,2) );
    end
    estimablecols=max(abs(opp-speye(size(opp,1))),[],1)<= sX.tol;
    iXX=sX.v(:,[1:sX.rk])*diag(1./sX.ds(1:sX.rk).^2)*sX.v(:,[1:sX.rk])';dof=size(SPM.xX.X,1)-sX.rk; %iXX=pinv(SPM.xX.X'*SPM.xX.X);dof=size(SPM.xX.X,1)-rank(SPM.xX.X);
    if DOORTHOG
        orthcols=zeros(1,size(SPM.xX.X,2));
        for n=1:size(SPM.xX.X,2)
            x=SPM.xX.X(:,[1:n-1,n+1:size(SPM.xX.X,2)]);
            y=SPM.xX.X(:,n);
            yorth=y-x*pinv(x'*x)*x'*y;
            orthcols(n)=100*(norm(yorth)/max(eps,norm(y))).^2;
        end
    end
    sizecols=estimablecols'./max(eps,sqrt(diag(iXX)));  % C*pinv(X'*X)*(X'*(X*C')) / sqrt(C*pinv(X'*X)*C')
    pwd0=pwd;
    try, cd(qafolder); end
    for nplot=1,%:2
        if nplot==1 && (isfield(SPM,'xCon')&&~isempty(SPM.xCon)), ncons=numel(SPM.xCon); connames={SPM.xCon.name};
        elseif nplot==2 && (isfield(SPM,'xConOriginal')&&~isempty(SPM.xConOriginal)), ncons=numel(SPM.xConOriginal); connames={SPM.xConOriginal.name};
        else ncons=0;
        end
        if ncons
            estimablecons=false(1,ncons);sizecons=zeros(1,ncons);
            A=[];iA=[];jA=[];
            for nc=1:ncons
                if nplot==1, c=SPM.xCon(nc).c;
                elseif nplot==2, c=SPM.xConOriginal(nc).c;
                end
                if size(c,2)==1
                    jA=[jA;nc];
                    idx=find(any(c,2));
                    iA=[iA;setdiff(idx(:),iA)];
                    [i1,i2]=ismember(iA,idx);
                    A(nc,i1)=sum(c(idx(i2(i1)),:),2)';
                    estimablecons(nc)=all(all(abs(opp*c - c) <= sX.tol));
                    k=max([sum(c(c>0)) -sum(c(c<0)) abs(sum(c))]);
                    sizecons(nc)=estimablecons(nc)*k/max(eps,sqrt(c'*iXX*c)); %sizecons(nc)=sqrt((c'*c)/(c'*iXX*c));
                    if DOORTHOG
                        x=SPM.xX.X*null(c');
                        y=SPM.xX.X*c;
                        yorth=y-x*pinv(x'*x)*x'*y;
                        orthcons(nc)=100*(norm(yorth)/max(eps,norm(y))).^2;
                    end
                end
            end
            [nill,idx]=sort(0*sum(A,1)'+1e-6*iA(:)); iA=iA(idx);A=A(:,idx);
            A=A(jA,:); Aall=sum(A,2); Apos=sum(A.*(A>0),2); Aneg=sum(A.*(A<0),2);
            ytickval=unique(round(linspace(1,size(A,1),20))); xtickval=unique(round(linspace(1,size(A,2),20)));
            
            hfig=figure('units','norm','position',[.1 .3 .8 .6],'color','w','menubar','none');
            cmap=.25*0+.75*jet(255);cmap(128,:)=[1 1 1];colormap(cmap);
            axes('units','norm','position',[.15 .4 .25 .5]);
            if nnz(A)<1e2, conn_menu_plotmatrix(.8*A/max(abs(A(:))),'colormapcdata',round(128+127*(A/max(abs(A(:))))),'colormap',cmap,'shape','hdo');
            else imagesc(A/max(abs(A(:)))); axis tight; colormap(cmap); box off
            end
            set(gca,'ydir','reverse');
            if nnz(A)<1e2, axis equal; end
            %conn_menu_plotmatrix(.8*A/max(abs(A(:))),'scalesigned',1); set(gca,'ydir','reverse');
            %h=imagesc(A);
            %hold on;[i,j]=find(A); plot(repmat(j(:)',5,1)+repmat([-.5 -.5 .5 .5 -.5]',1,numel(j)), repmat(i(:)',5,1)+repmat([-.5 .5 .5 -.5 -.5]',1,numel(i)), 'k-'); hold off;
            set(gca,'clim',max(.01,max(abs(A(:))))*[-1 1],'ytick',ytickval,'yticklabel',[],'xtick',xtickval,'xticklabel',[],'xlim',[.5 size(A,2)+.5], 'ylim',[.5 size(A,1)+.5]);
            hold on; text(xtickval,1.05*(size(A,1)+.5)*ones(1,numel(xtickval)),regexprep(SPM.xX.name(iA(xtickval)),'\*bf\(1\)$',''),'rotat',-90,'horizontalalignment','left','interpreter','none'); hold off;
            txtnill=regexprep(SPM.xX.name(iA(xtickval)),'\*bf\(1\)$',''); [txtnill{estimablecols(iA(xtickval))}]=deal(''); hold on; text(xtickval,1.05*(size(A,1)+.5)*ones(1,numel(xtickval)),txtnill,'rotat',-90,'horizontalalignment','left','interpreter','none','backgroundcolor','y'); hold off;
            hold on; text(.5-.05*ones(1,numel(ytickval))*size(A,2),ytickval,connames(jA(ytickval)),'horizontalalignment','right','interpreter','none'); hold off;
            txtnill=connames(jA(ytickval)); [txtnill{estimablecons(jA(ytickval))}]=deal(''); hold on; text(.5-.05*ones(1,numel(ytickval))*size(A,2),ytickval,txtnill,'horizontalalignment','right','interpreter','none','backgroundcolor','k','color','w'); hold off;
            txt=arrayfun(@(a,b,c)sprintf('warning: scale = %s',mat2str(.1*round(10*setdiff(unique([a,b,c]),[0])))),Aall,Apos,Aneg,'uni',0);
            txtnill=ismember(round(Aall*1e3)/1e3,[0,1])&ismember(round(Apos*1e3)/1e3,[0,1])&ismember(round(-Aneg*1e3)/1e3,[0,1]);
            if any(txtnill), txt(txtnill)=repmat({''},1,nnz(txtnill)); end
            hold on; h=text(mean(xtickval)*ones(1,size(A,1)),1:size(A,1),txt,'color','k','rotat',0,'horizontalalignment','center','interpreter','none','fontweight','bold'); hold off; set(h,'fontsize',ceil(get(h(1),'fontsize')*.8));
            grid on;
            hax=gca;
            h=colorbar; try, set(h,'color',.5*[1 1 1],'box','off'); end
            hold on; text(mean(xtickval),.5-.1*(size(A,1)+.5),'CONTRAST DEFINITIONS','horizontalalignment','center','fontsize',13,'fontweight','bold'); hold off;
            
            axes('units','norm','position',[.5 .4 .20 .2]);
            %axes('units','norm','position',[.65 .65 .3 .25]);
            h=bar(sizecols(iA)); set(h,'facecolor',.85*[1 1 1]);
            ylim=1.1*max([1,max(sizecols(iA(xtickval))),max(sizecons(jA(ytickval)))]); %ylim=1.1*max(1,max(sizecols(iA(xtickval))));
            hold on; text(xtickval,-.05*ylim*ones(1,numel(xtickval)),regexprep(SPM.xX.name(iA(xtickval)),'\*bf\(1\)$',''),'rotat',-90,'horizontalalignment','left','interpreter','none'); hold off;
            txtnill=regexprep(SPM.xX.name(iA(xtickval)),'\*bf\(1\)$',''); [txtnill{estimablecols(iA(xtickval))}]=deal(''); hold on; text(xtickval,-.05*ylim*ones(1,numel(xtickval)),txtnill,'rotat',-90,'horizontalalignment','left','interpreter','none','backgroundcolor','k','color','w'); hold off;
            set(gca,'xlim',[0,max(xtickval)+1],'xtick',[],'ylim',[0 ylim],'box','off','xcolor',.5*[1 1 1],'ycolor',.5*[1 1 1]);
            grid on;
            hold on; text(mean(xtickval),ylim*1.1,'EFFECT ESTIMABILITY','horizontalalignment','center','fontsize',13,'fontweight','bold'); hold off;
            h=ylabel('T-statistic @ SNR=1'); set(h,'fontsize',10,'rotation',90,'color',.5*[1 1 1]);
            
            axes('units','norm','position',[.75 .4 .20 .2]);
            %axes('units','norm','position',[.2 .2 .3 .25]);
            h=bar(sizecons(jA)); set(h,'facecolor',.85*[1 1 1]);
            ylim=1.1*max([1,max(sizecols(iA(xtickval))),max(sizecons(jA(ytickval)))]); %ylim=1.1*max(1,max(sizecons(jA(ytickval))));
            hold on; text(ytickval,-.05*ylim*ones(1,numel(ytickval)),connames(jA(ytickval)),'rotat',-90,'horizontalalignment','left','interpreter','none'); hold off;
            txtnill=connames(jA(ytickval)); [txtnill{estimablecons(jA(ytickval))}]=deal(''); hold on; text(ytickval,-.05*ylim*ones(1,numel(ytickval)),txtnill,'rotat',-90,'horizontalalignment','left','interpreter','none','backgroundcolor','k','color','w'); hold off;
            set(gca,'xlim',[0,max(ytickval)+1],'xtick',[],'ylim',[0 ylim],'box','off','xcolor',.5*[1 1 1],'ycolor',.5*[1 1 1]);
            grid on;
            hold on; text(mean(ytickval),ylim*1.1,'CONTRAST ESTIMABILITY','horizontalalignment','center','fontsize',13,'fontweight','bold'); hold off;
            
            if DOORTHOG
                axes('units','norm','position',[.5 .7 .20 .2]);
                %axes('units','norm','position',[.65 .65 .3 .25]);
                h=bar(orthcols(iA)); set(h,'facecolor',.85*[1 1 1]);
                ylim=1.1*100;%max(1,max(orthcols(iA(xtickval))));
                %hold on; text(xtickval,-.05*ylim*ones(1,numel(xtickval)),regexprep(SPM.xX.name(iA(xtickval)),'\*bf\(1\)$',''),'rotat',-90,'horizontalalignment','left','interpreter','none'); hold off;
                %txtnill=regexprep(SPM.xX.name(iA(xtickval)),'\*bf\(1\)$',''); [txtnill{estimablecols(iA(xtickval))}]=deal(''); hold on; text(xtickval,-.05*ylim*ones(1,numel(xtickval)),txtnill,'rotat',-90,'horizontalalignment','left','interpreter','none','backgroundcolor','k','color','w'); hold off;
                set(gca,'xlim',[0,max(xtickval)+1],'xtick',[],'ylim',[0 ylim],'box','off','xcolor',.5*[1 1 1],'ycolor',.5*[1 1 1]);
                grid on;
                hold on; text(mean(xtickval),ylim*1.1,'EFFECT ORTHOGONALITY','horizontalalignment','center','fontsize',13,'fontweight','bold'); hold off;
                h=ylabel('percent variance'); set(h,'fontsize',10,'rotation',90,'color',.5*[1 1 1]);
                try, ytick=get(gca,'ytick'); yticklabel=cellstr(get(gca,'yticklabel')); set(gca,'ytick',ytick,'yticklabel',cellfun(@(x)[x,'%'],yticklabel,'uni',0)); end
                
                axes('units','norm','position',[.75 .7 .20 .2]);
                %axes('units','norm','position',[.2 .2 .3 .25]);
                h=bar(orthcons(jA)); set(h,'facecolor',.85*[1 1 1]);
                ylim=1.1*100;%max(1,max(orthcons(ytickval)));
                %hold on; text(ytickval,-.05*ylim*ones(1,numel(ytickval)),connames(ytickval),'rotat',-90,'horizontalalignment','left','interpreter','none'); hold off;
                %txtnill=connames(ytickval); [txtnill{estimablecons(ytickval)}]=deal(''); hold on; text(ytickval,-.05*ylim*ones(1,numel(ytickval)),txtnill,'rotat',-90,'horizontalalignment','left','interpreter','none','backgroundcolor','k','color','w'); hold off;
                set(gca,'xlim',[0,max(ytickval)+1],'xtick',[],'ylim',[0 ylim],'box','off','xcolor',.5*[1 1 1],'ycolor',.5*[1 1 1]);
                grid on;
                hold on; text(mean(ytickval),ylim*1.1,'CONTRAST ORTHOGONALITY','horizontalalignment','center','fontsize',13,'fontweight','bold'); hold off;
                try, ytick=get(gca,'ytick'); yticklabel=cellstr(get(gca,'yticklabel')); set(gca,'ytick',ytick,'yticklabel',cellfun(@(x)[x,'%'],yticklabel,'uni',0)); end
            end
            
            if any(~estimablecols(iA))||any(~estimablecons(jA))
                axes('units','norm','position',[.05 .05 .1 .1]);
                hold on; text(1,.5,'estimable contrasts/effects','horizontalalignment','left'); hold off;
                hold on; text(1,-.5,'non-estimable contrasts/effects','horizontalalignment','left','backgroundcolor','k','color','w'); hold off;
                set(gca,'xlim',[0 2],'ylim',[-1 1],'xtick',[],'ytick',[]); axis off; %,'box','on');
            end
            
            %axes('units','norm','position',[.35 .80 .5 .05]);
            %imagesc(repmat(reshape(estimablecons,1,[]),[1,1,3]));
            %set(gca,'xtick',[],'ytick',[],'box','on','position',[.35 .80 0 .05]+get(hax,'position').*[0 0 1 0]);
            %hax1=gca;
            %h=title('Contrasts estimability'); set(h,'fontweight','normal');
            %h=xlabel('(white = ok; black = non-estimable)'); set(h,'fontsize',ceil(get(h,'fontsize')/2));
            %axes('units','norm','position',[.05 .35 .02 .3]);
            %imagesc(repmat(reshape(estimablecols(iA),[],1),[1,1,3]));
            %set(gca,'ytick',[],'xtick',[],'box','on');
            %hax2=gca;
            %h=ylabel('Parameters estimability'); set(h,'fontweight','normal');
            
            %if nplot==1, fname=sprintf('QA_SPM_contrasts_%s.subject%03d.jpg',pwd1name,nsub); %fname='QA_SPM_contrasts.jpg';
            %elseif nplot==2, fname=sprintf('QA_SPM_contrasts_alloriginalcontrasts_%s.subject%03d.jpg',pwd1name,nsub); %fname='QA_SPM_contrasts_alloriginalcontrasts.jpg';
            %end
            
            conn_print(hfig,fullfile(qafolder,fname),'-nogui','-noerror',dpires,'-nopersistent');
            conn_args={'image_display',fname};
            opts={'forcecd'};
            conn_savematfile(conn_prepend('',fullfile(qafolder,fname),'.mat'),'conn_args','opts','-v7.3');
            
            %fh=fopen(conn_prepend('',fullfile(qafolder,fname),'.txt'),'wt');
            fh={};
            if all(estimablecons(jA)), fh{end+1}=sprintf('All Contrasts are estimable\n');
            else fh{end+1}=sprintf('List of non-estimable contrasts:\n'); for nc=find(~estimablecons(jA)), fh{end+1}=sprintf('%s\n',connames{jA(nc)}); end
            end
            if all(estimablecols(iA)), fh{end+1}=sprintf('All Design Matrix columns are estimable\n');
            else fh{end+1}=sprintf('List of non-estimable effects:\n'); for nc=find(~estimablecols(iA)), fh{end+1}=sprintf('%s\n',SPM.xX.name{iA(nc)}); end
            end
            %fclose(fh);
            conn_fileutils('filewrite_raw',conn_prepend('',fullfile(qafolder,fname),'.txt'),fh);
            info=struct('contrasts',{connames(jA)},'contrast_orthogonality',orthcons(jA),'contrast_estimability',sizecons(jA),...
                'effects',{SPM.xX.name(iA)},'effect_orthogonality',orthcols(iA),'effect_estimability',sizecols(iA));
            conn_fileutils('spm_jsonwrite',conn_prepend('',fullfile(qafolder,fname),'.json'),info,struct('indent',' '));
            if ishandle(hfig), delete(hfig); end
        end
    end
    try, (pwd0); end
else conn_disp('fprintf','No contrasts defined. Skipping file %s\n',fileSPM);
end
