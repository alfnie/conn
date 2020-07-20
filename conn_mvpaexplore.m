function fh=conn_mvpaexplore
global CONN_x CONN_h CONN_gui;

if isempty(CONN_gui)||~isfield(CONN_gui,'font_offset'), CONN_gui.font_offset=0; end
filepathresults=fullfile(CONN_x.folders.firstlevel_vv,CONN_x.vvAnalyses(CONN_x.vvAnalysis).name);
[iroi,isnew,ncomp]=cellfun(@(x)conn_v2v('match_extended',x),CONN_x.vvAnalyses(CONN_x.vvAnalysis).measures);
if any(isnew), conn_msgbox('Sorry, this option is not available until the associated first-level MVPA analyses have been run','',2); return; end
volname=arrayfun(@(iroi,ncomp)fullfile(filepathresults,sprintf('PCAcov_Measure%03d_Component%03d.nii',iroi,ncomp)),iroi,ncomp,'uni',0);
ICAPCA='MVPA';
%XR=arrayfun(@(iroi,ncomp)fullfile(fullfile(CONN_x.folders.firstlevel_vv,CONN_x.vvAnalyses(CONN_x.vvAnalysis).name),['BETA_Subject',num2str(nsubs,'%03d'),'_Condition',num2str(CONN_h.menus.m_analyses.icondition(nconditions),'%03d'),'_Measure',num2str(iroi,'%03d'),'_Component',num2str(ncomp,'%03d'),'.nii']),iroi,ncomp,'uni',0);

vol=spm_vol(char(volname));
B=100*spm_read_vols(vol); % cumulative variance
mask=B>0;
prctB=[];
for n=1:size(B,4),b=B(:,:,:,n); temp=histc(b(b>0),0:1:100); prctB(:,n)=temp; end

[x,y,z]=ndgrid(1:vol(1).dim(1),1:vol(1).dim(2),1:vol(1).dim(3));
xyz=vol(1).mat*[x(:) y(:) z(:) ones(numel(x),1)]';
imat=pinv(vol(1).mat);
volref=CONN_gui.refs.canonical.V;
Bref=reshape(spm_get_data(volref,pinv(volref.mat)*xyz),vol(1).dim(1:3)).^2;
names=arrayfun(@(n)sprintf('%s_%d',ICAPCA,n),1:size(B,4),'uni',0);
labels=repmat({''},1,numel(names));
TR=conn_get_rt;
Ncomponents=size(B,4);
VARx=reshape(sum(sum(sum(B,1),2),3)./sum(sum(sum(mask,1),2),3),1,[]);

boffset=[0 0 0 0];
conn_menu('frame2',boffset+[.045,.08,.91,.80],'');
poslist=boffset+[.17 .40 .30 .38];
ht3=conn_menu('listbox2',poslist,'',names,'<HTML>Select component(s) for display</HTML>',@(varargin)conn_mvpaexplore_update([0 1 1]));
for n=6:20, set(ht3,'string',repmat(' ',1,4*n),'fontname','monospaced','fontsize',8+CONN_gui.font_offset); if get(ht3,'extent')*[0 0 1 0]'>poslist(3), break; end; end
ht3fieldsize=sprintf('%d',n-1);ht3fieldsize2=sprintf('%d',2*(n-1));
temp=boffset+[.17 .78 .28/4 .04];
ht3b(1)=conn_menu('pushbutton2',temp+[0*temp(3) 0 0 0],'','component','Select component(s) for display',@(varargin)conn_mvpaexplore_update([1 0 0]));
ht3b(2)=conn_menu('pushbutton2',temp+[2*temp(3) 0 0 0],'','variance (%)','Percent variance in connectivity profiles at each voxel explained by each component (click to sort)',@(varargin)conn_mvpaexplore_update([2 0 0]));
ht3b(3)=conn_menu('pushbutton2',temp+[3*temp(3) 0 0 0],'','cumulative (%)','Cumulative percent variance in connectivity profiles at each voxel explained by the first N components (click to sort)',@(varargin)conn_mvpaexplore_update([2 0 0]));

posimage=[.52,.40,.41,.44];
ht4=conn_menu('image2',boffset+posimage,'Percent variance explained by first N components (cumulative)','');%,'',@conn_mvpaexplore_mtncallback);
uicontrol('style','text','units','norm','position',boffset+[posimage(1)+posimage(3)/2-.070,posimage(2)-1*.059,.070,.045],'string','threshold','fontname','default','fontsize',8+CONN_gui.font_offset,'backgroundcolor',CONN_gui.backgroundcolor,'foregroundcolor',CONN_gui.fontcolorA); 
%ht6=conn_menu('popup2',boffset+[.07,.36,.30,.04],'',{'<HTML><i> - MVPA tools:</i></HTML>','Compute spatial match to template','Flip sign of individual networks/components', 'Label individual networks/components', 'Create ICA parcellation ROI file'},'<HTML> - <i>Spatial correlation</i> computes the spatial correlation and dice coefficients between spatial component scores of each ICA network and a user-defined reference file/mask<br/> - <i>Flip sign</i> flips the spatial component positive/negative loadings for selected components<br/> - <i>Label</i> adds user-defined labels to identify each network/component<br/> - <i>ICA parcellation</i> creates an ROI file identifying the ICA network with the highest spatial component score for each voxel</HTML>',@(varargin)conn_mvpaexplore_tools);
nfacshown=1:numel(names);
nfacselected=nfacshown;

% conn_menu('frame2',boffset+[.15,.07,.80,.24],'');%'Component timeseries');
% ht1=conn_menu('listbox2',boffset+[.79,.08,.075,.15],'Subjects',arrayfun(@(n)sprintf('Subject %d',n),1:CONN_x.Setup.nsubjects,'uni',0),'Select subject(s) for display',@(varargin)conn_mvpaexplore_update([0 0 1]));
% ht2=conn_menu('listbox2',boffset+[.87,.08,.075,.15],'Conditions',conditions,'<HTML>Select condition(s) for display</HTML>',@(varargin)conn_mvpaexplore_update([0 0 1]));

ht5=conn_menu('hist',boffset+[.62,.12,.21,.19],''); %'Cumulative percent variance (%)');
%ht5=conn_menu('image2',boffset+[.62,.09,.21,.22],''); %'Cumulative percent variance (%)');
conn_menu('update',ht4,{permute(Bref,[2,1,3]),permute(B(:,:,:,1),[2,1,3,4]),permute(abs(B(:,:,:,1)),[2,1,3,4])},{struct('mat',vol(1).mat,'dim',vol(1).dim),[]});
set(ht4.h10,'string',num2str(90));
conn_menu('updatethr',[],[],ht4.h10);
% set(ht1,'max',2,'value',1);
% set(ht2,'max',2,'value',1);
set(ht3,'string',names,'max',2,'value',nfacselected);
%conn_menu('updateplotstack',ht5,prctB);

sortby=0;
showdisp=0;
conn_mvpaexplore_update([1 1 1]);
fh=@conn_mvpaexplore_update;
%set(CONN_h.screen.hfig,'pointer','arrow');

    function conn_mvpaexplore_update(doreset)
        if nargin<1, doreset=[0 1 1]; end
        try, if isfield(CONN_h,'menus')&&isfield(CONN_h.menus,'waiticonObj'), CONN_h.menus.waiticonObj.start; end; end
%         nsub=get(ht1,'value');
%         ncond=get(ht2,'value');
        ifac=get(ht3,'value');
        nfacselected=nfacshown(ifac);
        
        if doreset(1)
            if doreset(1)==abs(sortby), sortby=-sortby;
            else sortby=doreset(1);
            end
            switch(abs(sortby))
                case 1, m=1:numel(names);
                case 2, m=VARx;
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
        
        str=cellfun(@(a,b,c,d)sprintf(regexprep('%-Xs%-Ys%-Ys',{'Y','X'},{ht3fieldsize,ht3fieldsize2}),[a ' ' d],num2str(b),num2str(c)),names,num2cell(diff([0 VARx])),num2cell(VARx),labels,'uni',0);
        set(ht3,'string',str(nfacshown),'value',ifac,'fontname','monospaced','fontsize',8+CONN_gui.font_offset);
            
        if doreset(2)
            val=1; %get(ht7,'value');
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
                    conn_menu('updatehist',ht5,{[0,0:100,100],[0,prctB(:,max(nfacselected))',0],sprintf('mean = %.2f',VARx(max(nfacselected)))});
                    %conn_menu('updateplotstack',ht5,prctB(:,nfacselected));

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
        
%         if doreset(3)
%             %x=[];
%             x=0;
%             for isub=1:numel(nsub)
%                 xt=[];
%                 for icond=1:numel(ncond)
%                     if ~isempty(data{nsub(isub)}{ncond(icond)}), 
%                         xt=[xt data{nsub(isub)}{ncond(icond)}(:,nfacselected)];
%                     end
%                 end
%                 %x=[x xt];
%                 x=x+xt/numel(nsub); % average across subjects
%             end
%             conn_menu('update',ht5,x);
%         end
        try, if isfield(CONN_h,'menus')&&isfield(CONN_h.menus,'waiticonObj'), CONN_h.menus.waiticonObj.stop; end; end
        
        return;

    end


    function [str0,str]=conn_mvpaexplore_mtncallback(varargin)
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

