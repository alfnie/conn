function fh=conn_icaexplore
global CONN_x CONN_h CONN_gui;

if isempty(CONN_gui)||~isfield(CONN_gui,'font_offset'), conn_font_init; end
filepathresults=fullfile(CONN_x.folders.firstlevel_vv,CONN_x.vvAnalyses(CONN_x.vvAnalysis).name);
ICAPCA='ICA';
if ~conn_existfile(fullfile(filepathresults,[ICAPCA,'.Timeseries.mat'])), ICAPCA='PCA'; end
if ~conn_existfile(fullfile(filepathresults,[ICAPCA,'.Timeseries.mat'])), conn_msgbox('Sorry, this option is not available until the first-level ICA analyses have been re-run','',2); return; end
%set(CONN_h.screen.hfig,'pointer','watch');drawnow
[data,conditions,weights]=deal({}); conn_loadmatfile(fullfile(filepathresults,[ICAPCA,'.Timeseries.mat']),'data','conditions','weights');
if conn_existfile(fullfile(filepathresults,[ICAPCA,'.Maps.nii'])), filevol=fullfile(filepathresults,[ICAPCA,'.Maps.nii']); vol=conn_fileutils('spm_localvol',filevol); 
else filevol=fullfile(filepathresults,[ICAPCA,'.ROIs.nii']); vol=conn_fileutils('spm_localvol',filevol);
end
B=spm_read_vols(vol);
[maxB,imaxB]=max(B,[],4);
[x,y,z]=ndgrid(1:vol(1).dim(1),1:vol(1).dim(2),1:vol(1).dim(3));
xyz=vol(1).mat*[x(:) y(:) z(:) ones(numel(x),1)]';
imat=pinv(vol(1).mat);
volref=CONN_gui.refs.canonical.V;
Bref=reshape(spm_get_data(volref,pinv(volref.mat)*xyz),vol(1).dim(1:3)).^2;
names=arrayfun(@(n)sprintf('%s_%d',ICAPCA,n),1:size(B,4),'uni',0);
labels=repmat({''},1,numel(names));
if conn_existfile(fullfile(filepathresults,[ICAPCA,'.ROIs.txt'])),
    try
        trefnames=regexp(fileread(fullfile(filepathresults,[ICAPCA,'.ROIs.txt'])),'[\n\r]*','split');
        trefnames=trefnames(cellfun('length',trefnames)>0);
        if numel(trefnames)==numel(labels), labels=trefnames; end
        for n=1:numel(labels), if strcmp(labels{n},num2str(n)), labels{n}=''; end; end
    end
end
TR=conn_get_rt;
Ncomponents=size(B,4);
KURx=zeros(1,Ncomponents);
SKWx=zeros(1,Ncomponents);
for n=1:Ncomponents,
    b=B(:,:,:,n);
    mask=~isnan(b);
    mb=mean(b(mask));
    sb=std(b(mask),0);
    KURx(n)=mean(((b(mask)-mb)/sb).^4);
    SKWx(n)=mean(((b(mask)-mb)/sb).^3);
end
CCx=nan(1,Ncomponents);
VARt=zeros(numel(data{1}),Ncomponents);
ZCt=zeros(numel(data{1}),Ncomponents);
validc=zeros(numel(data{1}),1);
for n1=1:numel(data) % subjects
    for n2=1:numel(data{n1}) % conditions
        if ~isempty(data{n1}{n2})
            mask=weights{n1,n2}{1}>0;
            b=data{n1}{n2}(mask,:)*sqrt(nnz(mask));
            remove=weights{n1,n2}{5}(mask)==1;
            VARt(n2,:)=VARt(n2,:)+std(b,0,1);
            mb=mean(b,1);
            b=conn_bsxfun(@minus,b,mb);
            b(remove,:)=0;
            ZCt(n2,:)=ZCt(n2,:)+((mean((b(1:end-1,:)<=0&b(2:end,:)>0)|(b(1:end-1,:)>=0&b(2:end,:)<0),1))/TR(min(numel(TR),n1)))/2;
            validc(n2)=validc(n2)+1;
        end
    end
end
VARt=conn_bsxfun(@rdivide,VARt,validc);
ZCt=conn_bsxfun(@rdivide,ZCt,validc);

boffset=[0 0 0 0];
conn_menu('frame2',boffset+[.045,.35,.91,.53],'');%'Component loadings');
poslist=boffset+[.07 .40 .40 .38];
ht3=conn_menu('listbox2',poslist,'',names,'<HTML>Select component(s) for display</HTML>',@(varargin)conn_icaexplore_update([0 1 1]));
for n=6:20, set(ht3,'string',repmat(' ',1,6*n),'fontname','monospaced','fontsize',8+CONN_gui.font_offset); if get(ht3,'extent')*[0 0 1 0]'>poslist(3), break; end; end
ht3fieldsize=sprintf('%d',n-1);ht3fieldsize2=sprintf('%d',2*(n-1));
%htemp=uicontrol('style','pushbutton','unit','norm','position',[.17,.85,.3,.05],'string','test','fontsize',8+CONN_gui.font_offset,'fontname','monospaced');set(htemp,'units','characters'); temp=get(htemp,'position'); set(htemp,'position',[temp(1:2) 13 max(1,temp(4))],'units','norm'); temp=get(htemp,'position'); delete(htemp);
%conn_menu('text2',boffset+[.17+.28/5 .85 2*.28/5 .04],'','spatial');
%conn_menu('text2',boffset+[.17+3*.28/5 .85 2*.28/5 .04],'','temporal');
temp=boffset+[.07 .78 .38/6 .04];
ht3b(1)=conn_menu('pushbutton2',temp+[0*temp(3) 0 0 0],'','network','Select network(s) for display',@(varargin)conn_icaexplore_update([1 0 0]));
ht3b(2)=conn_menu('pushbutton2',temp+[2*temp(3) 0 0 0],'','kurtosis','Spatial kurtosis (click to sort)',@(varargin)conn_icaexplore_update([2 0 0]));
ht3b(3)=conn_menu('pushbutton2',temp+[3*temp(3) 0 0 0],'','skewness','Spatial skewness (click to sort)',@(varargin)conn_icaexplore_update([3 0 0]));
ht3b(4)=conn_menu('pushbutton2',temp+[4*temp(3) 0 0 0],'','variability','Temporal component timeseries standard deviation for the selected condition and averaged across all subjects (click to sort)',@(varargin)conn_icaexplore_update([4 0 0]));
ht3b(5)=conn_menu('pushbutton2',temp+[5*temp(3) 0 0 0],'','frequency','Temporal component timeseries frequency (Hz) for the selected condition and averaged across all subjects (click to sort)',@(varargin)conn_icaexplore_update([5 0 0]));

% ht3b=uicontrol('style','text','units','norm','position',boffset+[.27,.85,.10,.04],'string','','fontname','default','fontsize',8+CONN_gui.font_offset,'backgroundcolor',CONN_gui.backgroundcolorA,'foregroundcolor',CONN_gui.fontcolorA); 
ht7=conn_menu('popup2',boffset+[.65 .84 .15 .04],'',{'Spatial components','ICA parcellation'},'<HTML> - <i>Spatial components</i> display spatial component scores for each selected ICA component (networks)<br/> - <i>ICA parcellation</i> displays the ICA network (among selected components) with the highest loading for each voxel</HTML>',@(varargin)conn_icaexplore_update([0 1 0]));
posimage=[.52,.40,.41,.44];
ht4=conn_menu('image2',boffset+posimage,'','','',@conn_icaexplore_mtncallback);
uicontrol('style','text','units','norm','position',boffset+[posimage(1)+posimage(3)/2-.070,posimage(2)-1*.059,.070,.045],'string','threshold','fontname','default','fontsize',8+CONN_gui.font_offset,'backgroundcolor',CONN_gui.backgroundcolor,'foregroundcolor',CONN_gui.fontcolorA); 
ht6=conn_menu('popup2',boffset+[.07,.36,.30,.04],'',{'<HTML><i> - ICA tools:</i></HTML>','Compute spatial match to template','Flip sign of individual networks/components', 'Label individual networks/components', 'Create ICA parcellation ROI file'},'<HTML> - <i>Spatial correlation</i> computes the spatial correlation and dice coefficients between spatial component scores of each ICA network and a user-defined reference file/mask<br/> - <i>Flip sign</i> flips the spatial component positive/negative loadings for selected components<br/> - <i>Label</i> adds user-defined labels to identify each network/component<br/> - <i>ICA parcellation</i> creates an ROI file identifying the ICA network with the highest spatial component score for each voxel</HTML>',@(varargin)conn_icaexplore_tools);
nfacshown=1:numel(names);
nfacselected=nfacshown;

conn_menu('frame2',boffset+[.045,.07,.91,.24],'');%'Component timeseries');
ht1=conn_menu('listbox2',boffset+[.79,.08,.075,.15],'Subjects',arrayfun(@(n)sprintf('Subject %d',n),1:CONN_x.Setup.nsubjects,'uni',0),'Select subject(s) for display',@(varargin)conn_icaexplore_update([0 0 1]));
ht2=conn_menu('listbox2',boffset+[.87,.08,.075,.15],'Conditions',conditions,'<HTML>Select condition(s) for display</HTML>',@(varargin)conn_icaexplore_update([0 0 1]));
ht5=conn_menu('image2',boffset+[.06,.09,.71,.17],'Temporal components');
conn_menu('update',ht4,{permute(Bref,[2,1,3]),permute(B(:,:,:,1),[2,1,3,4]),permute(abs(B(:,:,:,1)),[2,1,3,4])},{struct('mat',vol(1).mat,'dim',vol(1).dim),[]});
set(ht4.h10,'string',num2str(2));
conn_menu('updatethr',[],[],ht4.h10);
set(ht1,'max',2,'value',1);
set(ht2,'max',2,'value',1);
set(ht3,'max',2,'value',nfacselected);

sortby=0;
showdisp=0;
conn_icaexplore_update([1 1 1]);
fh=@conn_icaexplore_update;
%set(CONN_h.screen.hfig,'pointer','arrow');

    function conn_icaexplore_update(doreset)
        if nargin<1, doreset=[0 1 1]; end
        try, if isfield(CONN_h,'menus')&&isfield(CONN_h.menus,'waiticonObj'), CONN_h.menus.waiticonObj.start; end; end
        nsub=get(ht1,'value');
        ncond=get(ht2,'value');
        ifac=get(ht3,'value');
        nfacselected=nfacshown(ifac);
        
        if doreset(1)
            if doreset(1)==abs(sortby), sortby=-sortby;
            else sortby=doreset(1);
            end
            switch(abs(sortby))
                case 1, m=1:numel(names);
                case 2, m=KURx;
                case 3, m=SKWx;
                case 4, m=mean(VARt(ncond,:),1); doreset(2)=doreset(2)|doreset(3);
                case 5, m=mean(ZCt(ncond,:),1);  doreset(2)=doreset(2)|doreset(3);
            end
            if sortby>0, [nill,idx]=sort(m);
            else         [nill,idx]=sort(m,'descend');
            end
            nfacshown=idx;
            invidx=idx; invidx(idx)=1:numel(idx);
            ifac=sort(invidx(nfacselected));
            nfacselected=nfacshown(ifac);
            doreset(2)=1;
        end
        
        str=cellfun(@(a,b,c,d,e,f)sprintf(regexprep('%-Xs%-Ys%-Ys%-Ys%-Ys',{'Y','X'},{ht3fieldsize,ht3fieldsize2}),[a ' ' f],num2str(b),regexprep(num2str(c),'NaN','-'),num2str(d),num2str(e)),names,num2cell(KURx),num2cell(SKWx),num2cell(mean(VARt(ncond,:),1)),num2cell(mean(ZCt(ncond,:),1)),labels,'uni',0);
        set(ht3,'string',str(nfacshown),'value',ifac,'fontname','monospaced','fontsize',8+CONN_gui.font_offset);
            
        if doreset(2)
            val=get(ht7,'value');
            switch(val)
                case 1,
                    if conn_surf_dimscheck(vol), %if isequal(vol.dim,conn_surf_dims(8).*[1 1 2]),
                        if numel(nfacselected)>1&&~CONN_h.menus.m_results_surfhires
                            t1=reshape(B(:,:,[1,size(B,3)/2+1],nfacselected),size(B,1)*size(B,2),2,1,[]);
                            t2=abs(t1);
                            t1=t1(CONN_gui.refs.surf.default2reduced,:,:,:);
                            t2=t2(CONN_gui.refs.surf.default2reduced,:,:,:);
                            conn_menu('update',ht4,{CONN_gui.refs.surf.defaultreduced,t1,t2},{vol(1),1});
                        else
                            t1=reshape(B(:,:,:,nfacselected),size(CONN_gui.refs.surf.default(1).vertices,1),2,1,[]);
                            t2=abs(t1);
                            conn_menu('update',ht4,{CONN_gui.refs.surf.default,t1,t2},{vol(1),1});
                        end
                    else
                        conn_menu('update',ht4,{permute(Bref,[2,1,3]),permute(B(:,:,:,nfacselected),[2,1,3,4]),permute(abs(B(:,:,:,nfacselected)),[2,1,3,4])},{struct('mat',vol(1).mat,'dim',vol(1).dim),[]});
                    end
                    if val~=showdisp
                        maxb=max(abs(B(:)));
                        set(ht4.h9,'string',num2str(maxb));
                        conn_menu('updatecscale',[],[],ht4.h9);
                    end
                case 2,
                    %[maxB,tempimaxB]=max(abs(B(:,:,:,nfacselected)),[],4);
                    [maxB,tempimaxB]=max(max(0,B(:,:,:,nfacselected)),[],4);
                    if conn_surf_dimscheck(vol), %if isequal(vol.dim,conn_surf_dims(8).*[1 1 2]),
                        t1=reshape(tempimaxB,size(CONN_gui.refs.surf.default(1).vertices,1),2,1,[]);
                        t2=reshape(maxB,size(CONN_gui.refs.surf.default(1).vertices,1),2,1,[]);
                        conn_menu('update',ht4,{CONN_gui.refs.surf.default,t1,t2},{vol(1),1});
                    else
                        conn_menu('update',ht4,{permute(Bref,[2,1,3]),permute(tempimaxB,[2,1,3,4]),permute(maxB,[2,1,3,4])},{struct('mat',vol(1).mat,'dim',vol(1).dim),[]});
                    end
                    if val~=showdisp||~isequal(tempimaxB,imaxB)
                        maxb=max(abs(tempimaxB(:)));
                        set(ht4.h9,'string',num2str(maxb));
                        conn_menu('updatecscale',[],[],ht4.h9);
                    end
                    imaxB=tempimaxB;
            end
            showdisp=val;
        end
        
        if doreset(3)
            %x=[];
            x=0;
            for isub=1:numel(nsub)
                xt=[];
                for icond=1:numel(ncond)
                    if ~isempty(data{nsub(isub)}{ncond(icond)}), 
                        xt=[xt data{nsub(isub)}{ncond(icond)}(:,nfacselected)];
                    end
                end
                %x=[x xt];
                x=x+xt/numel(nsub); % average across subjects
            end
            conn_menu('update',ht5,x);
        end
        try, if isfield(CONN_h,'menus')&&isfield(CONN_h.menus,'waiticonObj'), CONN_h.menus.waiticonObj.stop; end; end
        
        return;

    end

    function conn_icaexplore_tools(varargin)
        if nargin>0, option=varargin{1};
        else option=get(ht6,'value')-1;
            set(ht6,'value',1);
        end
        switch(option)
            case 1, conn_icaexplore_spatialcorr(varargin{2:end});
            case 2, conn_icaexplore_flipsigns(varargin{2:end});
            case 3, conn_icaexplore_label(varargin{2:end});
            case 4, conn_icaexplore_parcellation(varargin{2:end});
        end
    end
        
    function conn_icaexplore_spatialcorr(varargin)
       Y=reshape(B,[],size(B,4));
       maskY=~isnan(Y); 
       Y(isnan(Y))=0; 
       v=[];
       maskv=[];
       refnames=[];
       hfig=figure('unit','norm','position',[.2 .2 .5 .5],'name','ICA match to template','numbertitle','off','menubar','none','color','w');
       hselect=uicontrol('units','norm','position',[.2 .925 .6 .05],'style','pushbutton','string','Select template file','callback',@(varargin)conn_icaexplore_spatialcorr_selecttemplate);
       hmenu=uicontrol('units','norm','position',[.2 .875 .6 .05],'style','popupmenu','string',{'Spatial correlation (correlation coefficient)','Spatial overlap of suprathreshold areas (dice coefficient)'},'value',1,'callback',@conn_icaexplore_spatialcorr_update);
       hax=axes('units','norm','position',[.2 .2 .6 .6]);
       try
           conn_icaexplore_spatialcorr_selecttemplate(fullfile(fileparts(which(mfilename)),'utils','otherrois','networksonly.nii'));
       end
       
        function conn_icaexplore_spatialcorr_selecttemplate(varargin)
            if nargin>0, filename=varargin{1};
            else
                [file_name,file_path]=conn_fileutils('uigetfile','*.nii;*.img','Select reference file/mask',pwd);
                if isequal(file_name,0), return; end
                filename=fullfile(file_path,file_name);
            end
            set(hselect,'tooltipstring',filename);
            maskV=conn_fileutils('spm_vol',filename);
            v=conn_fileutils('spm_get_data',maskV,pinv(maskV(1).mat)*xyz);
            maskv=~isnan(v);
            v(isnan(v))=0;
            if size(v,1)==1&&max(v)>1&&~any(rem(v,1)), v=full(sparse(v(v>0),find(v>0),1,max(v),numel(v))); end
            refnames=arrayfun(@(n)sprintf('reference volume #%d',n),1:size(v,1),'uni',0);
            if conn_existfile(conn_prepend('',filename,'.txt')),
                try
                    trefnames=regexp(fileread(conn_prepend('',filename,'.txt')),'[\n\r]*','split');
                    trefnames=trefnames(cellfun('length',trefnames)>0);
                    if numel(trefnames)==size(v,1), refnames=trefnames; end
                end
            end
            conn_icaexplore_spatialcorr_update;
        end
        function conn_icaexplore_spatialcorr_update(varargin)
            switch(get(hmenu,'value'))
                case 1, 
                    tY=maskY.*conn_bsxfun(@minus,Y,sum(Y.*maskY,1)./sum(maskY,1));
                    tv=maskv.*conn_bsxfun(@minus,v,sum(v.*maskv,2)./sum(maskv,2));
                    r=(tv*tY)./conn_bsxfun(@times,sqrt(sum(tv.^2,2)),sqrt(sum(tY.^2,1)));
                    fprintf('Correlation matrix (rows are reference volumes; columns are %s components)\n',ICAPCA);
                case 2,
                    thr=str2num(get(ht4.h10,'string'));
                    tY=maskY.*(Y>thr);
                    if sum(~isnan(unique(v)))>2, tv=maskv.*(v>thr);
                    else tv=maskv.*(v==max(v(:)));
                    end
                    r=2*(tv*tY)./conn_bsxfun(@plus,sum(tv,2),sum(tY,1));
                    fprintf('Dice coefficient matrix (rows are reference volumes; columns are %s components)\n',ICAPCA);
                
            end
            for ni=1:size(r,1), fprintf('%f ',r(ni,:)); fprintf('\n'); end
            if size(r,2)>3
                [maxr,imax]=sort((r),2,'descend');
                for ni=1:size(r,1),
                    conn_disp('fprintf','best three matches to %s are %s_%d (r=%f), %s_%d (r=%f), and %s_%d (r=%f)\n',refnames{ni},ICAPCA,imax(ni,1),r(ni,imax(ni,1)),ICAPCA,imax(ni,2),r(ni,imax(ni,2)),ICAPCA,imax(ni,3),r(ni,imax(ni,3)));
                end
            elseif size(r,2)>1,
                [maxr,imax]=max((r),[],2);
                for ni=1:size(r,1),
                    conn_disp('fprintf','best match to %s is %s_%d (r=%f)\n',refnames{ni},ICAPCA,imax(ni),r(ni,imax(ni)));
                end
            end
            if size(r,1)>3
                [maxr,imax]=sort((r),1,'descend');
                for ni=1:size(r,2),
                    conn_disp('fprintf','best three matches to %s_%d are %s (r=%f), %s (r=%f), and %s (r=%f)\n',ICAPCA,ni,refnames{imax(1,ni)},r(imax(1,ni),ni),refnames{imax(2,ni)},r(imax(2,ni),ni),refnames{imax(3,ni)},r(imax(3,ni),ni));
                end
            elseif size(r,1)>1,
                [maxr,imax]=max((r),[],1);
                for ni=1:size(r,2),
                    conn_disp('fprintf','best match to %s_%d is %s (r=%f)\n',ICAPCA,ni,refnames{imax(ni)},r(imax(ni),ni));
                end
            end
            cla(hax);
            rplot=[(max(0,r)) nan(size(r,1),1) max(max(0,r),[],2); nan(1,size(r,2)+2); max(max(0,r),[],1) nan(1,2)];
            [nill,hpatch]=conn_menu_plotmatrix(rplot,'parent',hax);
            set(hax,'ydir','reverse','xtick',[1:size(r,2) size(r,2)+2],'ytick',[1:size(r,1) size(r,1)+2],'xticklabel',[names {'BestMatch'}],'yticklabel',[refnames {'BestMatch'}],'xlim',[-1 size(rplot,2)+.5],'ylim',[.5 size(rplot,1)+2])
            try, set(hax,'xticklabelrotation',90,'yticklabelrotation',0,'TickLabelInterpreter','none'); 
            catch, set(hax,'xticklabel',[arrayfun(@num2str,1:numel(refnames),'uni',0) {'','BestMatch'}]); 
            end
            try, set(hax,'gridColor',.5*[1 1 1]); end
            axis(hax,'equal');
            grid(hax,'on');
        end
    end

    function conn_icaexplore_flipsigns(varargin)
        if nargin>0, iflip=varargin{1};
        else
            nset=listdlg('name','Flip signs','PromptString','Select component(s) to flip','ListString',names,'SelectionMode','multiple','ListSize',[200 200]);
            if ~isempty(nset)
                B(:,:,:,nset)=-B(:,:,:,nset);
                [maxB,imaxB]=max(B,[],4);
                SKWx(nset)=-SKWx(nset);
                for n1=1:numel(data) % subjects
                    for n2=1:numel(data{n1}) % conditions
                        if ~isempty(data{n1}{n2})
                            data{n1}{n2}(:,nset)=-data{n1}{n2}(:,nset);
                        end
                    end
                end
                conn_icaexplore_update;
                answ=conn_questdlg({'Saving these changes stores the (currently displayed) group-level spatial maps and component timeseries for future reference','(note: subject-level backprojected spatial maps or any imported group-ICA ROIs are NOT changed by this operation;','you may flip the sign manually in any second-level analysis by changing the sign of the between-source contrast elements)',' ','Not saving these changes does not modify any of the stored maps or timeseries.','The sign-flip operation effect is still viewable in the current ''summary'' display, but ','switching to a different tab or clicking again on the ''sumary'' button will revert to the original (saved) maps',' ','Save these changes now?'},'','Yes','Not now','Not now');
                if ~(isempty(answ)||strcmp(answ,'Not now')), 
                    conn_savematfile(fullfile(filepathresults,[ICAPCA,'.Timeseries.mat']),'data','-append');
                    try, conn_fileutils('deletefile',fullfile(filepathresults,[ICAPCA,'.Maps.nii'])); end
                    vol=spm_create_vol(vol);
                    for n1=1:size(B,4), vol(n1).fname=fullfile(filepathresults,[ICAPCA,'.Maps.nii']); vol(n1)=spm_write_vol(vol(n1),B(:,:,:,n1)); end
                    if conn_server('util_isremotefile',filevol), conn_cache('push',filevol); end
                end
            end
        end
    end


    function conn_icaexplore_label(varargin)
       hfig=figure('unit','norm','position',[.2 .2 .3 .3],'name','ICA network labels','numbertitle','off','menubar','none','color','w');
       hmenu=uicontrol('units','norm','position',[.2 .7 .6 .2],'style','popupmenu','string',names,'value',1,'callback',@(varargin)conn_icaexplore_label_update('select'));
       hedit=uicontrol('units','norm','position',[.2 .4 .6 .2],'style','edit','string','','callback',@(varargin)conn_icaexplore_label_update('edit'));
       uicontrol('units','norm','position',[.2 .6 .6 .1],'style','text','string','Label:');
       hok=uicontrol('style','pushbutton','units','norm','position',[.55,.025,.2,.1],'string','OK','tooltipstring','Accept changes','callback','uiresume(gcbf)');
       hcancel=uicontrol('style','pushbutton','units','norm','position',[.75,.025,.2,.1],'string','Cancel','callback','delete(gcbf)');
       newlabels=labels;
       conn_icaexplore_label_update('select');
       set(hfig,'handlevisibility','on','hittest','on');
       uiwait(hfig);
       if ~ishandle(hfig), return; end
       delete(hfig);
       labels=newlabels;
       %fh=fopen(fullfile(filepathresults,[ICAPCA,'.ROIs.txt']),'wt');
       fh={};
       for nl=1:numel(labels), 
           if isempty(labels{nl}), fh{end+1}=sprintf('%d\n',nl);
           else fh{end+1}=sprintf('%s\n',labels{nl}(labels{nl}>=32));
           end
       end
       %fclose(fh);
       conn_fileutils('filewrite_raw',fullfile(filepathresults,[ICAPCA,'.ROIs.txt']),fh);
       conn_icaexplore_update;
       
        function conn_icaexplore_label_update(str)
            switch(str)
                case 'select'
                    ilabel=get(hmenu,'value');
                    set(hedit,'string',newlabels{ilabel});
                case 'edit'
                    ilabel=get(hmenu,'value');
                    newlabels{ilabel}=get(hedit,'string');
            end
        end
    end

    function conn_icaexplore_parcellation(varargin)
       if nargin>0, filename=varargin{1}; 
       else
           minSize=0; % note: remove clusters with fewer than 10 voxels?
           selectedcomponents=nfacselected;
           allcomponentnames=labels;
           isemptynames=cellfun('length',allcomponentnames)==0;
           allcomponentnames(isemptynames)=names(isemptynames);
           selectedcomponents=listdlg('name','ICA Parcellation','PromptString','Select component(s) to include','ListString',allcomponentnames,'SelectionMode','multiple','InitialValue',selectedcomponents,'ListSize',[200 200]);
           if isempty(selectedcomponents), return; end
           [file_name,file_path]=uiputfile('*.nii;*.img','Select parcellation filename',pwd);
           if isequal(file_name,0), return; end
           filename=fullfile(file_path,file_name);
           [file_path,file_name,file_ext]=fileparts(filename);
           if isequal(file_path,0), return; end
           if ~any(strcmp(file_ext,{'.img','.nii'})), filename=fullfile(file_path,[file_name,'.nii']); end
           %[maxB,imaxB]=max(abs(B(:,:,:,selectedcomponents)),[],4);
           [maxB,imaxB]=max(max(0,B(:,:,:,selectedcomponents)),[],4);
           e0=struct('fname',filename,'descrip','conn ICA parcellation file','mat',vol(1).mat,'dim',vol(1).dim,'n',[1,1],'pinfo',[1;0;0],'dt',[spm_type('uint16'),spm_platform('bigend')]);
           try, conn_fileutils('spm_unlink',e0.fname); end
           thr=str2num(get(ht4.h10,'string'));
           txyz=round(imat(1:3,:)*xyz);
           posneg={'+','-'};
           nc=0;
           l3=zeros(size(maxB));
           %fh=fopen(conn_prepend('',filename,'.txt'),'wt');
           fh={};
           for np=1:numel(selectedcomponents)
               tB=B(:,:,:,selectedcomponents(np));
               l2=maxB>thr & imaxB==np;
               [nl,l]=conn_clusters(l2,txyz);
               [nill,idxnl]=sort(nl,'descend');
               for ni2=1:numel(nl)
                   if nl(idxnl(ni2))>=minSize
                       nc=nc+1;
                       tmask=l==idxnl(ni2);
                       l3(tmask)=nc;
                       if isempty(labels{selectedcomponents(np)}), tname=names{selectedcomponents(np)};
                       else tname=labels{selectedcomponents(np)};
                       end
                       fh{end+1}=sprintf('%s%s (%d,%d,%d) n=%d\n',tname,posneg{1+(mean(tB(tmask))<0)},round(mean(xyz(1,tmask))),round(mean(xyz(2,tmask))),round(mean(xyz(3,tmask))),nl(idxnl(ni2)));
                   end
               end
           end
           %fclose(fh);
           conn_fileutils('filewrite_raw',conn_prepend('',filename,'.txt'),fh);
           spm_write_vol(e0,l3);
           conn_msgbox(sprintf('ICA parcellation file %s saved\n',filename),'',true);
       end        
    end

    function [str0,str]=conn_icaexplore_mtncallback(varargin)
        txyz=varargin{1};
        if nargin>1, ifac=varargin{2}; 
        else ifac=[]; 
        end
        txyz=round(imat*[reshape(txyz(1:3),[],1);1]);
        txyz=max(1,min(vol(1).dim(1:3),txyz(1:3)'));
        tb=reshape(B(txyz(1),txyz(2),txyz(3),:),1,[]);
        tb(isnan(tb))=0;
        [maxb,idx]=sort(abs(tb),'descend');
        str0={};
        str={};
        try
            if any(maxb~=0)
                if showdisp==1
                    tf=sprintf('%s_%d (%.3f)',ICAPCA,nfacselected(ifac),tb(nfacselected(ifac)));
                    str0=[str0 {tf}];
                elseif showdisp==2
                    tf=sprintf('%s_%d (%.3f)',ICAPCA,nfacselected(imaxB(txyz(1),txyz(2),txyz(3))),tb(nfacselected(imaxB(txyz(1),txyz(2),txyz(3)))));
                    str0=[str0 {tf}];
                end
                hf=arrayfun(@(a,b)sprintf(' %s_%d (%.3f)',ICAPCA,a,b),idx(1:min(numel(idx),3)),tb(idx(1:min(numel(idx),3))),'uni',0);
                str=[str {' ','Highest component scores for this voxel:'} hf];
            end
        end
    end
        

end

