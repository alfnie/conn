function hfig=conn_displayroi(option,varargin)
% internal function ROI-to-ROI results display
%

global CONN_x CONN_gui;
if ~nargin, option='init'; end
hfig=[];
if ~ischar(option), % gui-callback
    if ishandle(option)&&isequal(get(option,'type'),'figure'), hfig=option;
    else hfig=gcbf;
    end
    if isempty(hfig), hfig=gcf; end
    option=varargin{2}; 
    varargin=varargin(3:end); 
    margin=nargin-2; 
else
    margin=nargin;
    if isempty(regexpi(option,'^init')), hfig=gcf; end
end

init=false;
initxy=false;
DOSCALEF=true;
switch(lower(option)),
    case {'init','initfile','initspm'}
        init=true; 
        ncon=1;
        if isempty(CONN_gui)||~isfield(CONN_gui,'font_offset'), conn_font_init; end
        if margin>1,  
            data.source=0; %if numel(varargin)>=2, data.source=varargin{2}; else data.source=[]; end
            data.thres=1;
            data.side=3;
            data.initfile='';
            if strcmpi(option,'initfile')
                if conn_fileutils('isdir',varargin{1})
                    if conn_existfile(fullfile(varargin{1},'ROI.mat')), varargin{1}=fullfile(varargin{1},'ROI.mat');
                    elseif conn_existfile(fullfile(varargin{1},'SPM.mat')), varargin{1}=fullfile(varargin{1},'SPM.mat');
                    else error('unable to find ROI.mat or SPM.mat file in results directory %s',varargin{1});
                    end
                end
                [fpath,fname,fext]=fileparts(varargin{1});
                if strcmp([fname,fext],'SPM.mat'), option='initspm'; end
            end
            %h=conn_msgbox('computing ROI-level results, please wait...','conn_displayroi');
            if strcmpi(option,'initfile'), 
                data.initfile=varargin{1};
                results=conn_loadmatfile(varargin{1});
                results=results.ROI;
                filepath=fileparts(varargin{1});
                if isempty(data.source), data.source=0; end
                if numel(varargin)>=2, ncon=varargin{2}; end
                if numel(varargin)>=3, data.thres=varargin{3}; end
                if numel(varargin)>=4, data.side=varargin{4}; end
            elseif strcmpi(option,'initspm'), 
                data.initfile=varargin{1};
                filepath=fileparts(varargin{1});
                if isempty(filepath), filepath=pwd; end
                if isempty(data.source), data.source=0; end
                SPM=struct; conn_loadmatfile(varargin{1},'SPM');
                SPM_h=permute(SPM.xX_multivariate.h,[3,4,1,2]);
                if size(SPM_h,3)>1||size(SPM_h,4)>1, SPM_h=sqrt(sum(abs(SPM_h(:,:,:)).^2,3)); end
                SPM_F=permute(SPM.xX_multivariate.F,[3,4,1,2]);
                if isfield(SPM.xX_multivariate,'p'), SPM_p=SPM.xX_multivariate.p;
                else
                    SPM_p=nan(size(SPM_F)); 
                    idxvalid=~(isnan(SPM_F)|SPM_F==0); 
                    switch(SPM.xX_multivariate.statsname)
                        case 'T',  SPM_p(idxvalid)=spm_Tcdf(-SPM_F(idxvalid),SPM.xX_multivariate.dof);
                        case 'F',  SPM_p(idxvalid)=1-spm_Fcdf(SPM_F(idxvalid),SPM.xX_multivariate.dof(1),SPM.xX_multivariate.dof(2));
                        case 'X',  SPM_p(idxvalid)=1-spm_Xcdf(SPM_F(idxvalid),SPM.xX_multivariate.dof);
                    end
                end
                SPM_F(SPM_F==0)=nan;
                nrois=size(SPM_F,1);
                [fpath,fname,fext]=fileparts(SPM.xY.VY(1).fname);
                if isempty(fpath), fpath=filepath; end
                info=[];
                try, info=conn_jsonread(fullfile(fpath,[fname,fext])); end
                if isempty(info)
                    info.names=arrayfun(@(n)sprintf('ROI#%04d',n),1:nrois,'uni',0);
                    info.coords=repmat({[0 0 0]},nrois,1);
                end
                if ~iscell(info.coords), info.coords=num2cell(info.coords,2); end % [nroisx3]
                %z=permute(randn([nrois,nrois,size(SPM.xX_multivariate.Zfiles)]),[3,2,4,1]); 
                vol=conn_fileutils('spm_vol',char(SPM.xX_multivariate.Zfiles));
                z=conn_fileutils('spm_read_vols',vol);
                z=permute(reshape(z,[nrois,nrois,size(SPM.xX_multivariate.Zfiles)]),[3,2,4,1]); % subjects x rois (targets) x conditions x rois (seeds)
                % subjects x rois x conditions
                validrois=any(all(all(~isnan(z),1),3),4); % note: eliminates ROIs that have no valid connectivity data with any other ROI
                if ~all(validrois)
                    fprintf('warning: disregarding the following ROIs from this group-level analysis due to missing-data %s\n',sprintf('%s ',info.names{~validrois}));
                    z=z(:,validrois,:,validrois);
                    SPM_h=SPM_h(validrois,validrois,:,:);
                    SPM_F=SPM_F(validrois,validrois,:,:);
                    SPM_p=SPM_p(validrois,validrois,:,:);
                    info.names=info.names(validrois);
                    info.coords=info.coords(validrois);
                end
                results=struct(...
                    'xX', struct('isSurface',SPM.xX.isSurface,'isMtx',SPM.xX.isMtx,'SelectedSubjects',SPM.xX.SelectedSubjects,'name',{SPM.xX_multivariate.Xnames},'X',SPM.xX_multivariate.X), ...
                    'data',z,...
                    'h', SPM_h,...
                    'F', SPM_F,...
                    'p', SPM_p,...
                    'dof', SPM.xX_multivariate.dof,...
                    'statsname',SPM.xX_multivariate.statsname,...
                    'c', SPM.xX_multivariate.C,...
                    'c2', SPM.xX_multivariate.M,...
                    'ynames',{SPM.xX_multivariate.Ynames}, ...
                    'names',{info.names},...
                    'names2',{info.names},...
                    'xyz', {info.coords},...
                    'xyz2',{info.coords} ...
                    ); 
                if numel(varargin)>=2, ncon=varargin{2}; end
                if numel(varargin)>=3, data.thres=varargin{3}; end
            else
                [results,filepath]=conn_process(varargin{:});
                
            end
            if isempty(filepath), data.defaultfilepath=pwd;
            else data.defaultfilepath=filepath;
            end
            if ~isempty(data.initfile)&&isempty(fileparts(data.initfile)), data.initfile=fullfile(data.defaultfilepath,data.initfile); end
            %close(h);
        else
            data.thres=1;
            data.side=3;
            [filename,filepath]=conn_fileutils('uigetfile','*ROI*.mat');
            if ~ischar(filename), return; end
            results=conn_loadmatfile(fullfile(filepath,filename));results=results.ROI;
            data.initfile=fullfile(filepath,filename);
            data.source=0;
            if isempty(filepath), data.defaultfilepath=pwd;
            else data.defaultfilepath=filepath;
            end
        end
        hmsginit=conn_msgbox('Initializing. Please wait...','conn_displayroi',-1);
        data.roifile=fullfile(data.defaultfilepath,'ROI.mat');
        if isempty(data.initfile), data.initfile=data.roifile; end
        if ~isfield(results(1),'data')||numel(results)>1
            if 1 % note: disregard any pre-computed mvpa stats for consistency
                data.MVPAh=[];
                data.MVPAF=[];
                data.MVPAp=[];
                data.MVPAdof=[];
                data.MVPAstatsname=[];
            else
                data.MVPAh=cat(1,results.MVPAh);
                data.MVPAF=cat(1,results.MVPAF);
                data.MVPAp=cat(1,results.MVPAp);
                temp={results.MVPAdof};
                if any(cellfun('length',temp)>1), temp=cellfun(@(x)[ones(1,max(0,2-length(x))),x(:)'],temp,'uni',0); end
                data.MVPAdof=cell2mat(temp(:));
                data.MVPAstatsname=results(1).MVPAstatsname;
            end
            data.h=cat(1,results.h);
            data.F=cat(1,results.F);
            data.p=cat(1,results.p);
            data.dof=cat(1,results.dof);
            data.statsname=results(1).statsname;            
            temp=cat(4,results.y); % subjects x rois (targets) x nconditions x rois (seeds)
            data.results=results(1);
            data.results.data=temp;
            clear temp;
        else
            data.MVPAh=[];
            data.MVPAF=[];
            data.MVPAp=[];
            data.MVPAdof=[];
            data.MVPAstatsname=[];
            data.h=results.h;
            data.F=results.F;
            data.p=results.p;
            data.dof=repmat(results.dof,numel(data.F),1);
            data.statsname=results.statsname;            
            data.results=results;
        end
        if isempty(data.MVPAdof), data.MVPAdofstr={};
        elseif ~any(any(diff(data.MVPAdof,1,1),1),2),
            if size(data.MVPAdof,2)==1, data.MVPAdofstr=repmat({['(',num2str(data.MVPAdof(1)),')']},numel(data.MVPAdof),1);
            else data.MVPAdofstr=repmat({['(',num2str(data.MVPAdof(1,1)),',',num2str(data.MVPAdof(1,2)),')']},size(data.MVPAdof,1),1);
            end
        else
            if size(data.MVPAdof,2)==1, data.MVPAdofstr=arrayfun(@(dof)['(',num2str(dof(1)),')'],data.MVPAdof,'uni',0);
            else data.MVPAdofstr=cellfun(@(dof)['(',num2str(dof(1)),',',num2str(dof(2)),')'],num2cell(data.MVPAdof,2),'uni',0);
            end
        end
        if size(data.F,2)>=size(data.F,1)&&max(max(abs(data.F(:,1:size(data.F,1))-data.F(:,1:size(data.F,1))')))<1e-10, data.issymmetric=true; else data.issymmetric=false; end
        if strcmp(data.statsname,'T')&&all(data.dof(:)==data.dof(1)), data.p2=nan(size(data.p));data.p2(~isnan(data.F))=spm_Tcdf(data.F(~isnan(data.F)),data.dof(1)); 
        else data.p2=1-data.p; 
        end
        %if size(data.dof,2)>1&&~any(diff(data.dof(:,1))), data.statsname=[data.statsname,'(',num2str(data.dof(1)),')']; end
        if ~any(any(diff(data.dof,1,1),1),2),
            if size(data.dof,2)==1, data.dofstr=repmat({['(',num2str(data.dof(1)),')']},numel(data.dof),1);
            else data.dofstr=repmat({['(',num2str(data.dof(1,1)),',',num2str(data.dof(1,2)),')']},size(data.dof,1),1);
            end
        else
            if size(data.dof,2)==1, data.dofstr=arrayfun(@(dof)['(',num2str(dof(1)),')'],data.dof,'uni',0);
            else data.dofstr=cellfun(@(dof)['(',num2str(dof(1)),',',num2str(dof(2)),')'],num2cell(data.dof,2),'uni',0);
            end
        end
        data.names=results(1).names;
        data.names=regexprep(data.names,{'_1_1$','^rs\.','^rsREL\.','^aal\.'},'');
        data.namesreduced=regexprep(data.names,{'^BA\.(\d+) \(([LR])\)\. .*','^\((-?\d+),(-?\d+),(-?\d+)\)$','^SLrois\.|^aal\.|^atlas\.|^networks\.','\s\(([LlRr])\)','([^\(\)]*[^\.])\s*\(.+\)\s*$'},{'$1$2','($1 $2 $3)','',' ${lower($1)}','$1'});
        data.names2=results(1).names2;
        data.names2=regexprep(data.names2,{'_1_1$','^rs\.','^rsREL\.','^aal\.'},'');
        data.names2reduced=regexprep(data.names2,{'^BA\.(\d+) \(([LR])\)\. .*','^\((-?\d+),(-?\d+),(-?\d+)\)$','^SLrois\.|^aal\.|^atlas\.|^networks\.','\s\(([LlRr])\)','([^\(\)]*[^\.])\s*\(.+\)\s*$'},{'$1$2','($1 $2 $3)','',' ${lower($1)}','$1'});
        data.xyz=cat(1,results(1).xyz{:});
        data.xyz2=cat(1,results(1).xyz2{:});
        data.displaytheserois=1:length(data.names);
        data.tfceZ=[];
        data.tfceZpeaks=[];
        data.tfceZd=[];
        data.cMVPAh=[];
        data.cMVPAF=[];
        data.cMVPAp=[];
        
        data.thr=.01;
        data.thrtype=1;
        data.mvpathr=.05;
        data.mvpathrtype=6;
        data.thres_defaults={ {.05,1,.05,32}, {.01,1,.05,6}, {.05,4,0,1}, {.05,2,0,1}, {.01,1,.05,30}, {.001,1,.05,15} };
        if ~isempty(data.thres), 
            if iscell(data.thres)
                [data.thr,data.thrtype,data.mvpathr,data.mvpathrtype]=deal(data.thres{:});
                data.thres=numel(data.thres_defaults)+1;
            else
                [data.thr,data.thrtype,data.mvpathr,data.mvpathrtype]=deal(data.thres_defaults{data.thres}{:});
                data.side=3;
            end
        end
        data.mvpasortresultsby=2;        
        data.mvpathrtype_all={'none',...
            'cluster-level p-uncorrected (SPC TFCE score)',       'cluster-level p-FDR corrected (SPC TFCE score)',         'cluster-level p-FWE corrected (SPC TFCE score)',...
            'cluster-level p-uncorrected (SPC mass/intensity)',   'cluster-level p-FDR corrected (SPC mass/intensity)',     'cluster-level p-FWE corrected (SPC mass/intensity)',...
            'cluster-level p-uncorrected (SPC size)',             'cluster-level p-FDR corrected (SPC size)',               'cluster-level p-FWE corrected (SPC size)',...
            'network-level p-uncorrected (NBS TFCE score)',       'network-level p-FDR corrected (NBS TFCE score)',         'network-level p-FWE corrected (NBS TFCE score)',...
            'network-level p-uncorrected (NBS mass/intensity)',   'network-level p-FDR corrected (NBS mass/intensity)',     'network-level p-FWE corrected (NBS mass/intensity)',...
            'network-level p-uncorrected (NBS size)',             'network-level p-FDR corrected (NBS size)',               'network-level p-FWE corrected (NBS size)',...
            'ROI-level p-uncorrected (ROI TFCE score)',           'ROI-level p-FDR corrected (ROI TFCE score)',             'ROI-level p-FWE corrected (ROI TFCE score)',...
            'ROI-level p-uncorrected (ROI mass/intensity)',       'ROI-level p-FDR corrected (ROI mass/intensity)',         'ROI-level p-FWE corrected (ROI mass/intensity)',...
            'ROI-level p-uncorrected (ROI size)',                 'ROI-level p-FDR corrected (ROI size)',                   'ROI-level p-FWE corrected (ROI size)',...
            'ROI-level p-uncorrected (multivariate omnibus test)',             'ROI-level p-FDR corrected (multivariate omnibus test)',...
            'cluster-level p-uncorrected (multivariate omnibus test)',         'cluster-level p-FDR corrected (multivariate omnibus test)'};
        data.mvpathrtype_isnonparam=ismember(1:numel(data.mvpathrtype_all),[2:28]);
        data.mvpathrtype_iscluster=ismember(1:numel(data.mvpathrtype_all),[2:10,31:32]);
        data.mvpathrtype_isnetwork=ismember(1:numel(data.mvpathrtype_all),11:19);
        data.mvpathrtype_isroi=ismember(1:numel(data.mvpathrtype_all),[20:28,29:30]);
        data.mvpathrtype_isce=ismember(1:numel(data.mvpathrtype_all),[2:4,11:13,20:22]);
        data.mvpathrtype_ismass=ismember(1:numel(data.mvpathrtype_all),[5:7,14:16,23:25]);
        data.mvpathrtype_issize=ismember(1:numel(data.mvpathrtype_all),[8:10,17:19,26:28]);
        data.mvpathrtype_isftest=ismember(1:numel(data.mvpathrtype_all),[29:30,31:32]);
        data.mvpathrtype_shown=[14:16,31:32,5:7,29:30,23:25,1]; %[15,24,6,1]; 
        %data.mvpaside=3;
        data.PERM=[];
        data.iPERM=[];
        data.clusters=[];
        data.view=0;
        data.proj=[];
        data.x=[];
        data.y=[];
        data.z=[];
        data.bgz=0;
        data.maxz=[];
        data.display='connectivity';
        %data.displayreduced=0;
        %data.displaytheserois=1:length(data.names2);
        data.displayreduced=1;
        data.displaylabels=1;
        data.displaybrains=1;
        data.display3d=0;
        data.displaygui=1;
        data.pausegui=1;
        data.mvpaenablethr=1;
        data.enablethr=1;
        data.displayconnectionstats=0;
        data.displayroilabelsinstats=0;
        data.displayallmeasures=0;
        data.visible='on';
        data.plotposition={[.01,.27,.53,.60],[.01,.05,.88,.9]}; 
        data.plotconnoptions.menubar=false;
        data.plotconnoptions.LINEWIDTH=2;
        data.plotconnoptions.LINESTYLEMTX=0;
        data.plotconnoptions.DOFFSET=.35;
        data.plotconnoptions.BSCALE=.25;
        data.plotconnoptions.BTRANS=.10;
        data.plotconnoptions.LTRANS=1;
        data.plotconnoptions.RSCALE=.75;
        data.plotconnoptions.LCOLOR=3;
        data.plotconnoptions.LCOLORSCALE=1;
        data.plotconnoptions.LCURVE=2;
        data.plotconnoptions.LBUNDL=.5;
        data.plotconnoptions.FONTSIZE=max(4,[0,1]+4+CONN_gui.font_offset);
        data.plotconnoptions.FONTANGLE=0;
        if 1, data.plotconnoptions.BCOLOR=.975*[1,1,1];
        elseif isfield(CONN_gui,'backgroundcolor'), data.plotconnoptions.BCOLOR=CONN_gui.backgroundcolor;
        else data.plotconnoptions.BCOLOR=[0.12 0.126 0.132];
        end
        data.plotconnoptions.NPLOTS=12;
        data.plotconnoptions.Projections={[0,-1,0;0,0,1;-1,0,0],[1,0,0;0,0,1;0,1,0],[1,0,0;0,1,0;0,0,1]};
        data.plotconnoptions.nprojection=1;
        data.plotconnoptions.Projections_axes={{'y','z'},{'x','z'},{'x','y'}};
        FSfolder=fullfile(fileparts(which('conn')),'utils','surf');
        rend(1)=reducepatch(conn_surf_readsurf(fullfile(FSfolder,'lh.pial.surf')),.02,'fast');
        rend(2)=reducepatch(conn_surf_readsurf(fullfile(FSfolder,'rh.pial.surf')),.02,'fast');
        %[xyz,faces]=read_surf(fullfile(FSfolder,'lh.cortex.surf'));
        %rend(1)=reducepatch(struct('vertices',xyz,'faces',faces+1),.02,'fast');
        %[xyz,faces]=read_surf(fullfile(FSfolder,'rh.cortex.surf'));
        %rend(2)=reducepatch(struct('vertices',xyz,'faces',faces+1),.02,'fast');
        data.plotconnoptions.rende=struct('vertices',cat(1,rend.vertices),'faces',[rend(1).faces; size(rend(1).vertices,1)+rend(2).faces]);
        data.xy2=200*[cos(2*pi*(0:numel(data.displaytheserois)-1)'/numel(data.displaytheserois)),sin(2*pi*(0:numel(data.displaytheserois)-1)'/numel(data.displaytheserois))]; 
        data.xy2_clusters=[];
        if isfield(CONN_gui,'refs')&&isfield(CONN_gui.refs,'canonical')&&isfield(CONN_gui.refs.canonical,'filename')&&~isempty(CONN_gui.refs.canonical.filename)
            filename=CONN_gui.refs.canonical.filename;
        else
            filename=fullfile(fileparts(which('spm')),'canonical','avg152T1.nii');
        end
        data.ref=spm_vol(filename);
        color1=data.plotconnoptions.BCOLOR;
        color2=color1; %.975*[1 1 1];
        color3=color2+(.5-color2)*.15; %.9*[1 1 1];
        foregroundcolor=.5*.8+.2*(1-round(color1));
        
        %if all(data.plotconnoptions.BCOLOR>.8), color2=data.plotconnoptions.BCOLOR; 
        %else color2=.94*[1,1,1];
        %end
        %color1=[.5/6,1/6,2/6];
        %color2=[.5/6,1/6,2/6];
        hmsg=[];%figure('units','norm','position',[.01,.1,.98,.8],'numbertitle','off','name','ROI second-level results. Initializing...','color',color1,'colormap',gray,'menubar','none','toolbar','none','interruptible','off');
        h0=get(0,'screensize');
        hfig=figure('visible','off','renderer','opengl','units','pixels','position',[h0(3)-.75*h0(3)+2,h0(4)-.9*h0(4)-48,.75*h0(3)-2*2,.9*h0(4)]);
        %h0=get(0,'screensize'); h0=h0(1,3:4)-h0(1,1:2)+1; h0=h0/max(1,max(abs(h0))/2000);
        %minheight=500;
        %hfig=figure('visible','off','renderer','opengl','units','pixels','position',[0*72+1,h0(2)-max(minheight,.5*h0(1))-48,h0(1)-0*72-1,max(minheight,.5*h0(1))]);
        data.hfig=hfig;
        set(hfig,'units','norm','numbertitle','off','name',['ROI second-level results explorer ',data.defaultfilepath],'color',color1,'colormap',gray,'menubar','none','toolbar','none','interruptible','off','tag','conn_displayroi','keypressfcn',@conn_displayroi_keypress,'windowbuttondownfcn',@(varargin)conn_display_windowbuttonmotionfcn('down'),'windowbuttonupfcn',@(varargin)conn_display_windowbuttonmotionfcn('up'),'visible','on'); 
        %uicontrol('style','frame','units','norm','position',[.0,.95,.5,.05],'backgroundcolor',color2,'foregroundcolor',color2);
        hframe1=uicontrol('style','frame','units','norm','position',[0,0,1,.27],'backgroundcolor',color2,'foregroundcolor',color2,'parent',data.hfig);
        hframe2=uicontrol('style','frame','units','norm','position',[0,.87,1,.13],'backgroundcolor',color3,'foregroundcolor',color3,'parent',data.hfig);
        hframe3=0;%uicontrol('style','frame','units','norm','position',[.62,.35,.33,.45],'backgroundcolor',1*[1 1 1],'foregroundcolor',.85*[1,1,1],'parent',data.hfig);
        huicontrol_cthr=uicontrol('style','popupmenu','units','norm','position',[.20,.965,.605,.03],'string',{...
            'standard settings for cluster-based inferences #1: Functional Network Connectivity',...
            'standard settings for cluster-based inferences #2: Spatial Pairwise Clustering',...
            'standard settings for cluster-based inferences #3: Threshold Free Cluster Enhancement',...
            'alternative settings for connection-based inferences: parametric univariate statistics ',...
            'alternative settings for ROI-based inferences: parametric multivariate statistics',...
            'alternative settings for network-based inferences: Network Based Statistics',...
            '<HTML><i>customize (advanced Family-Wise Error control settings)</i></HTML>'},'tag','highlight','fontname','arial','fontsize',8+CONN_gui.font_offset,'callback',{@conn_displayroi,'fwec.option'},'value',data.thres,'tooltipstring','Select false-positive control method','backgroundcolor',.9*[1,1,1]);
        huicontrol_cthr0=uicontrol('style','text','units','norm','position',[.03,.925,.17,.03],'fontsize',8+CONN_gui.font_offset,'string','connection threshold: p < ','horizontalalignment','right','fontweight','bold','foregroundcolor',1-color3,'backgroundcolor',color3,'interruptible','off','parent',data.hfig);
        huicontrol_cthr1=uicontrol('style','edit','units','norm','position',[.20,.925,.10,.03],'fontsize',8+CONN_gui.font_offset,'string',num2str(data.thr),'foregroundcolor',1-color3,'backgroundcolor',color3,'interruptible','off','callback',{@conn_displayroi,'fwec.connectionlevel.value'},'tooltipstring','Connection-level threshold value (false-positive threshold value for individual connections)','parent',data.hfig);
        huicontrol_cthr2=uicontrol('style','popupmenu','units','norm','position',[.325,.915,.25,.04],'fontsize',8+CONN_gui.font_offset,'string',{'p-uncorrected','p-FDR corrected','p-FDR corrected (TFCE)','p-FWE corrected (TFCE)','F/T/X stat'},'foregroundcolor',1-color3,'backgroundcolor',color3,'tooltipstring','<HTML>False-positive control type for individual connections</HTML>','interruptible','off','callback',{@conn_displayroi,'fwec.connectionlevel.type'},'value',max(1,min(5, data.thrtype)),'parent',data.hfig);
        huicontrol_cthr3=uicontrol('style','popupmenu','units','norm','position',[.605,.915,.20,.04],'fontsize',8+CONN_gui.font_offset,'string',{'positive contrast (one-sided)','negative contrast (one-sided)','two-sided'},'foregroundcolor',1-color3,'backgroundcolor',color3,'tooltipstring','Analysis results directionality','interruptible','off','callback',{@conn_displayroi,'fwec.connectionlevel.side'},'value',data.side,'parent',data.hfig);
        huicontrol_ccthr0=uicontrol('style','text','units','norm','position',[.03,.885,.17,.03],'fontsize',8+CONN_gui.font_offset,'string','cluster threshold: p < ','horizontalalignment','right','fontweight','bold','foregroundcolor',1-color3,'backgroundcolor',color3,'fontweight','bold','interruptible','off','parent',data.hfig);
        huicontrol_ccthr1=uicontrol('style','edit','units','norm','position',[.20,.885,.10,.03],'fontsize',8+CONN_gui.font_offset,'string',num2str(data.mvpathr),'foregroundcolor',1-color3,'backgroundcolor',color3,'interruptible','off','callback',{@conn_displayroi,'fwec.clusterlevel.value'},'tooltipstring','<HTML>Cluster-level threshold value (false-positive threshold value for individual clusters/groups of connections)','parent',data.hfig);
        huicontrol_ccthr2=uicontrol('style','popupmenu','units','norm','position',[.325,.875,.48,.04],'fontsize',8+CONN_gui.font_offset,'string',data.mvpathrtype_all(data.mvpathrtype_shown),'value',find(data.mvpathrtype_shown==data.mvpathrtype),'foregroundcolor',1-color3,'backgroundcolor',color3,'tooltipstring',...
            ['<HTML>Type of cluster- or ROI- level false-positive control',...
            '<br/> <br/> - choose <i>network</i> measures for <b>non-parametric network-level inferences</b> (NBS: Network Based Statistics, Zalesky et al. 2010)<br/> Networks represent maximal subgraphs of suprathreshold-connected ROIs (groups of ROIs and suprathreshold effects/connections among them)<br/> Network size and Network mass measures both represent measures of degree/cost of these subgraphs (i.e. number and strength of suprathreshold effects/connections within each graph) <br/> Network TFCE scores represent a combined measure of network size and mass, defined as the Threshold Free Cluster Enhancement score for the chosen support section (Smith and Nichols 2009) <br/> Multiple comparison correction is implemented at the network-level (FWE/FDR across multiple networks). Network-level inferences remain valid when used in combination with arbitrary (e.g. p-uncorrected) connection thresholds<br/> e.g. <b>connection-level threshold p &#60 0.01 (p-uncorrected) & cluster-threshold p &#60 0.05 (network p-FDR corrected)</b>',...
            '<br/> <br/> - choose <i>cluster</i> measures for <b>parametric cluster-level inferences</b> (FNC: Functional Network Connetivity, Jafri et al. 2008)<br/> Clusters represent groups of effects/connections within- or between- networks (networks here refers to groups of ROIs)<br/> multivariate omnibus-test represents a multivariate measure characterizing the strength of all effects/connections within a network or between two networks<br/> Multiple comparison correction is implemented at the cluster-level (FWE/FDR across multiple clusters). Cluster-level inferences remain valid when used in combination with arbitrary (e.g. p-uncorrected) connection thresholds<br/> e.g. <b>connection-level threshold p &#60 0.05 (p-uncorrected) & cluster-threshold p &#60 0.05 (cluster p-FDR corrected)</b>',...
            '<br/> <br/> - choose <i>cluster</i> measures for <b>non-parametric cluster-level inferences</b> (SPC: Spatial Pairwise Clustering, Zalesky et al. 2012)<br/> Clusters represent groups of suprathreshold effects/connections<br/> Cluster size and cluster mass measures both represent measures of degree/cost of these clusters (i.e. number and strength of suprathreshold effects/connections within each cluster)<br/> Cluster TFCE scores represent a combined measure of cluster size and mass, defined as the Threshold Free Cluster Enhancement score for the chosen support section (Smith and Nichols 2009)  <br/> Multiple comparison correction is implemented at the cluster-level (FWE/FDR across multiple clusters). Cluster-level inferences remain valid when used in combination with arbitrary (e.g. p-uncorrected) connection thresholds<br/> e.g. <b>connection-level threshold p &#60 0.01 (p-uncorrected) & cluster-threshold p &#60 0.05 (cluster p-FDR corrected)</b>',...
            '<br/> <br/> - choose <i>ROI</i> for <b>parametric or non-parametric ROI-level inferences</b><br/> multivariate omnibus-test represents a multivariate measure characterizing the strength of all effects/connections from each ROI<br/> ROI size and ROI mass measures both represent measures of degree/cost of each ROI (i.e. number and strength of suprathreshold effects/connections from each ROI)<br/> ROI TFCE scores represent a combined measure of ROI size and mass, defined as the Threshold Free Cluster Enhancement score for the chosen support section (Smith and Nichols 2009) <br/> Multiple comparison correction is implemented at the ROI-level (FWE/FDR across multiple ROIs). ROI-level inferences are invariant to the choice of connection-level threshold<br/> e.g. <b>connection-level threshold p &#60 0.01 (p-uncorrected) & cluster-threshold p &#60 0.05 (ROI p-FDR corrected)</b>',...
            '<br/> <br/> - choose <i>none</i> for <b>parametric connection-level inferences</b> <br/> This option is only appropriately corrected for multiple comparisons when used in combination with p-FDR corrected connection thresholds (but not with p-uncorrected connection thresholds)<br/> e.g. <b>connection-level threshold p &#60 0.05 (p-FDR corrected) & cluster-threshold none</b><br/></HTML>'],...
            'interruptible','off','callback',{@conn_displayroi,'fwec.clusterlevel.type'},'value',max([1,find(data.mvpathrtype_shown==data.mvpathrtype,1)]),'parent',data.hfig);
        huicontrol_cthr4=uicontrol('style','checkbox','units','norm','position',[.88,.87,.12,.03],'string','show details','tag','highlight','value',0,'fontsize',8+CONN_gui.font_offset,'horizontalalignment','right','foregroundcolor',1-color3,'backgroundcolor',color3,'callback',{@conn_displayroi,'advancedthr'},'tooltipstring','Displays advanced thresholding options'); 
        hhelp=uicontrol('style','pushbutton','units','norm','position',[.81,.965,.02,.03],'fontsize',7+CONN_gui.font_offset,'foregroundcolor',foregroundcolor,'backgroundcolor',color2,'string','?','tag','highlight','tooltipstring','<HTML>Documentation about available methods of statistical inference</HTML>','interruptible','off','callback',@(varargin)conn('gui_help','url','http://www.conn-toolbox.org/fmri-methods/cluster-level-inferences'),'parent',data.hfig);
        
        %huicontrol_ccthr3=uicontrol('style','popupmenu','units','norm','position',[.66,.35,.33,.04],'fontsize',8+CONN_gui.font_offset,'string',{'threshold seed ROIs (F-test)','threshold seed ROIs (NBS; by intensity)','threshold seed ROIs (NBS; by size)','threshold networks (NBS; by intensity)','threshold networks (NBS; by size)'},'foregroundcolor',1-color2,'backgroundcolor',color2,'fontweight','bold','horizontalalignment','right','value',data.mvpathrmeasure,'tooltipstring','Threshold individual seed ROIs or individual networks (subsets of connected ROIs)','interruptible','off','callback',{@conn_displayroi,'mvpathrmeasure'},'parent',data.hfig);
        data.handles=[...
            huicontrol_cthr,...
            huicontrol_cthr1,...
            huicontrol_cthr2,...
            huicontrol_cthr3,...
            uicontrol('style','text','units','norm','position',[.15,.875,.70,.07],'fontsize',7+CONN_gui.font_offset,'string','','foregroundcolor',.5*[1 1 1],'backgroundcolor',color3,'horizontalalignment','center','parent',data.hfig),...
            huicontrol_cthr4, ... %0, ...%uicontrol('style','listbox','units','norm','position',[.57,.31,.41,.39],'fontsize',8+CONN_gui.font_offset,'string',' ','max',2,'foregroundcolor',.9-.8*color2,'backgroundcolor',color2,'tooltipstring','Selecting one or several seed ROIs limits the analyses only to the connectivity between the selected seeds and all of the ROIs in the network','interruptible','off','callback',{@conn_displayroi,'list1'},'keypressfcn',@conn_menu_search,'visible','off','parent',data.hfig),...
            uicontrol('style','text','units','norm','position',[.05,.22,.90,.03],'fontsize',7+CONN_gui.font_offset,'string',sprintf('%-24s  %-20s  %+12s  %+12s  %+12s','Analysis Unit','Statistic','p-unc','p-FDR','p-FWE'),'foregroundcolor',.5*[1 1 1],'backgroundcolor',color2,'fontname','monospaced','horizontalalignment','left','parent',data.hfig),...
            uicontrol('style','listbox','units','norm','position',[.05,.07,.90,.15],'fontsize',7+CONN_gui.font_offset,'string',' ','tag','highlight','fontname','monospaced','foregroundcolor',round(1-color2),'backgroundcolor',color2,'tooltipstring','Statistics for each connection, ROI, or cluster. Right-click to export table to .txt file','max',2,'interruptible','off','callback',{@conn_displayroi,'list2'},'keypressfcn',@conn_menu_search,'parent',data.hfig),...
            uicontrol('style','checkbox','units','norm','position',[.05,.04,.20,.03],'fontsize',8+CONN_gui.font_offset,'string','display extended stats','tag','highlight','foregroundcolor',1-color2,'backgroundcolor',color2,'interruptible','off','callback',{@conn_displayroi,'displayconnectionstats'},'value',data.displayconnectionstats,'tooltipstring','Check to display additional and post-hoc statistics for all suprathreshold units','parent',data.hfig),... %uicontrol('style','text','units','norm','position',[.61,.50,.38,.04],'fontsize',8+CONN_gui.font_offset,'string','Define thresholds:','foregroundcolor',color2,'backgroundcolor',1-.5*color2,'fontweight','bold','horizontalalignment','left','parent',data.hfig),...
            uicontrol('style','text','units','norm','position',[.64,.76,.14,.03],'string','Display&Print','horizontalalignment','center','fontweight','bold','fontname','arial','fontsize',9+CONN_gui.font_offset,'foregroundcolor',foregroundcolor,'backgroundcolor',color2,'parent',data.hfig),... %uicontrol('style','pushbutton','units','norm','position',[.84,.05,.14,.04],'fontsize',8+CONN_gui.font_offset,'string','non-parametric stats','callback',{@conn_displayroi,'enableperm'},'tooltipstring','Enables permutation-test based statistics (Cluster and ROI size/mass statistics)','parent',data.hfig),...
            uicontrol('style','text','units','norm','position',[.80,.76,.14,.03],'string','Tools','horizontalalignment','center','fontweight','bold','fontname','arial','fontsize',9+CONN_gui.font_offset,'foregroundcolor',foregroundcolor,'backgroundcolor',color2,'parent',data.hfig),... %uicontrol('style','text','units','norm','position',[.61,.95,.38,.04],'fontsize',8+CONN_gui.font_offset,'string','Define connectivity matrix:','foregroundcolor',color2,'backgroundcolor',1-.5*color2,'fontweight','bold','horizontalalignment','left','parent',data.hfig),...
            uicontrol('style','pushbutton','units','norm','position',[.64,.40,.29,.04],'fontsize',7+CONN_gui.font_offset,'foregroundcolor',foregroundcolor,'backgroundcolor',color2,'string',sprintf('Analysis of %d connections among %d ROIs',numel(data.names)*(numel(data.names)-1)/2*(1+~data.issymmetric),numel(data.names)),'tag','highlight','tooltipstring','<HTML>Defines subset of ROIs to include in these analyses (among all the sources selected in the first-level analysis definition)<br/>note: this choice affects all inferences</HTML>','interruptible','off','callback',{@conn_displayroi,'roi.select'},'parent',data.hfig),...
            uicontrol('style','checkbox','units','norm','position',[.25,.04,.20,.03],'fontsize',8+CONN_gui.font_offset,'string','display extended roi labels','tag','highlight','foregroundcolor',1-color2,'backgroundcolor',color2,'interruptible','off','callback',{@conn_displayroi,'displayroilabelstats'},'value',data.displayroilabelsinstats,'visible','off','tooltipstring','Check to include complete labels when describing ROIs','parent',data.hfig),... %uicontrol('style','pushbutton','units','norm','position',[.61,.10,.38,.04],'fontsize',8+CONN_gui.font_offset,'string','Select all','tooltipstring','Looks at the connectivity between all ROIs in the network','interruptible','off','callback',{@conn_displayroi,'selectall'},'parent',data.hfig),...
            uicontrol('style','pushbutton','units','norm','position',[.64,.36,.29,.04],'fontsize',7+CONN_gui.font_offset,'foregroundcolor',foregroundcolor,'backgroundcolor',color2,'string','ROIs sorted using hierarchical clustering','tag','highlight','tooltipstring','<HTML>Defines ROIs order and clusters<br/>note: this choice affects all cluster-based inferences (but not connection- ROI- or network-based inferences)</HTML>','interruptible','off','callback',{@conn_displayroi,'roi.order'},'parent',data.hfig),...
            huicontrol_ccthr1,...
            huicontrol_ccthr2,...
            uicontrol('style','pushbutton','units','norm','position',[.80,.68,.13,.04],'fontsize',8+CONN_gui.font_offset,'string','Export mask','tag','highlight','foregroundcolor',foregroundcolor,'backgroundcolor',color2,'callback',{@conn_displayroi,'export_mask'},'tooltipstring','Exports list of suprathreshold connections in ROI-to-ROI connectivity matrix','parent',data.hfig),... %0,...%uicontrol('style','popupmenu','units','norm','position',[.68,.34,.15,.04],'fontsize',8+CONN_gui.font_offset,'string',{'connection-level results','seed-level results','network-level results'},'foregroundcolor',1-color2,'backgroundcolor',color2,'tooltipstring','Criteria for sorting results in statistics table','interruptible','off','callback',{@conn_displayroi,'mvpasort'},'value',data.mvpasortresultsby),...
            uicontrol('style','pushbutton','units','norm','position',[.80,.72,.13,.04],'fontsize',8+CONN_gui.font_offset,'string','Import values','tag','highlight','foregroundcolor',foregroundcolor,'backgroundcolor',color2,'callback',{@conn_displayroi,'import_values'},'tooltipstring','Imports individual connectivity values for each suprathreshold connection and for each subject into CONN toolbox as second-level covariates','parent',data.hfig),... %uicontrol('style','checkbox','units','norm','position',[.61,.36,.03,.03],'fontsize',8+CONN_gui.font_offset,'foregroundcolor',1-color2,'backgroundcolor',color2,'interruptible','off','callback',{@conn_displayroi,'mvpaenablethr'},'value',data.mvpaenablethr,'tooltipstring','Enable Seed/Network threshold','parent',data.hfig)...
            uicontrol('style','pushbutton','units','norm','position',[.80,.64,.13,.04],'fontsize',8+CONN_gui.font_offset,'string','Export data','tag','highlight','foregroundcolor',foregroundcolor,'backgroundcolor',color2,'callback',{@conn_displayroi,'export_data'},'tooltipstring','Exports all connectvity values (entire ROI-to-ROI matrix) for each individual subject and condition to matrix NIFTI file','parent',data.hfig),... 
            hframe1,...%uicontrol('style','pushbutton','units','norm','position',[0,0,.10,.025],'fontsize',8+CONN_gui.font_offset,'backgroundcolor',get(hfig,'color'),'string','display options','tooltipstring','Controls the way functional results are displayed (right-click on figure to get this same menu and remove this button)','callback','set(findobj(gcbf,''type'',''uicontextmenu'',''tag'',''conn_displayroi_plot''),''visible'',''on'')'),...
            uicontrol('style','text','units','norm','position',[.05,0,.90,.03],'string','','foregroundcolor',1-color2,'backgroundcolor',color2,'parent',data.hfig),...
            huicontrol_cthr0,...
            huicontrol_ccthr0,...
            hframe2,...
            hframe3, ...
            uicontrol('style','pushbutton','units','norm','position',[.71,.695,.06,.05],'string','Ring Print','tag','highlight','fontname','arial','fontweight','bold','fontsize',8+CONN_gui.font_offset,'callback',{@conn_displayroi,'ring_print'},'tooltipstring','<HTML><b>Ring print</b><br/>Prints current results on ROI-ring display</HTML>','parent',data.hfig), ...
            uicontrol('style','pushbutton','units','norm','position',[.64,.530,.06,.05],'string','Glass display','tag','highlight','fontname','arial','fontweight','bold','fontsize',8+CONN_gui.font_offset,'callback',{@conn_displayroi,'glass_view'},'tooltipstring','<HTML><b>Glass display</b><br/>Displays current results on 3d glass-brain</HTML>','parent',data.hfig), ...
            uicontrol('style','pushbutton','units','norm','position',[.71,.530,.06,.05],'string','Glass print','tag','highlight','fontname','arial','fontweight','bold','fontsize',8+CONN_gui.font_offset,'callback',{@conn_displayroi,'glass_print'},'tooltipstring','<HTML><b>Glass print</b><br/>Prints current results on 3d glass-brain</HTML>','parent',data.hfig), ...
            uicontrol('style','pushbutton','units','norm','position',[.64,.475,.06,.05],'string','Plot effects','tag','highlight','fontname','arial','fontsize',8+CONN_gui.font_offset,'callback',{@conn_displayroi,'cluster_view'},'tooltipstring','<HTML><b>Plot effects</b><br/>Explores/displays average effect sizes within each suprathreshold cluster or connection</HTML>'),...
            uicontrol('style','pushbutton','units','norm','position',[.71,.475,.06,.05],'string','Plot design','tag','highlight','fontname','arial','fontsize',8+CONN_gui.font_offset,'callback',{@conn_displayroi,'plot_design'},'tooltipstring','<HTML><b>Plot design</b><br/>Displays General Linear Model design matrix and additional details</HTML>'), ...
            uicontrol('style','pushbutton','units','norm','position',[.80,.58,.13,.04],'string','Bookmark','tag','highlight','foregroundcolor',foregroundcolor,'backgroundcolor',color2,'fontname','arial','fontsize',8+CONN_gui.font_offset,'callback',{@conn_displayroi,'bookmark'},'tooltipstring','Bookmark this second-level results explorer view'),...
            uicontrol('style','pushbutton','units','norm','position',[.80,.54,.13,.04],'string','Open folder','tag','highlight','foregroundcolor',foregroundcolor,'backgroundcolor',color2,'fontname','arial','fontsize',8+CONN_gui.font_offset,'callback',{@conn_displayroi,'openfolder'},'tooltipstring','Open folder containing current second-level results files'), ...
            uicontrol('style','pushbutton','units','norm','position',[.80,.50,.13,.04],'string','Graphic options','tag','highlight','foregroundcolor',foregroundcolor,'backgroundcolor',color2,'fontname','arial','fontsize',8+CONN_gui.font_offset,'callback',{@conn_displayroi,'menubar'},'tooltipstring','Show/hide menubar for advanced display/graphic options'), ...
            uicontrol('style','pushbutton','units','norm','position',[.64,.640,.06,.05],'string','Matrix display','tag','highlight','fontname','arial','fontweight','bold','fontsize',8+CONN_gui.font_offset,'callback',{@conn_displayroi,'matrix_view'},'tooltipstring','<HTML><b>Matrix display</b><br/>Displays current results on ROI-to-ROI matrix display</HTML>','parent',data.hfig), ...
            uicontrol('style','pushbutton','units','norm','position',[.71,.640,.06,.05],'string','Matrix print','tag','highlight','fontname','arial','fontweight','bold','fontsize',8+CONN_gui.font_offset,'callback',{@conn_displayroi,'matrix_print'},'tooltipstring','<HTML><b>Matrix print</b><br/>Prints current results on ROI-to-ROI matrix display</HTML>','parent',data.hfig), ...
            uicontrol('style','pushbutton','units','norm','position',[.64,.695,.06,.05],'string','Ring Display','tag','highlight','fontname','arial','fontweight','bold','fontsize',8+CONN_gui.font_offset,'callback',{@conn_displayroi,'ring_view'},'tooltipstring','<HTML><b>Ring display</b><br/>Displays current results on ROI-ring display</HTML>','parent',data.hfig), ...
            uicontrol('style','pushbutton','units','norm','position',[.64,.585,.06,.05],'string','Lines display','tag','highlight','fontname','arial','fontweight','bold','fontsize',8+CONN_gui.font_offset,'callback',{@conn_displayroi,'lines_view'},'tooltipstring','<HTML><b>Connections display</b><br/>Displays current results on ROI-to-ROI connections display</HTML>','parent',data.hfig), ...
            uicontrol('style','pushbutton','units','norm','position',[.71,.585,.06,.05],'string','Lines print','tag','highlight','fontname','arial','fontweight','bold','fontsize',8+CONN_gui.font_offset,'callback',{@conn_displayroi,'lines_print'},'tooltipstring','<HTML><b>Connections print</b><br/>Prints current results on ROI-to-ROI connections display</HTML>','parent',data.hfig) ...
            ];
        uiwrap(hhelp);
        uiwrap(data.handles([8,  12,14,  17,18,19,31,32,33]));
        bp=[36 26 27 28 34 35 29 30 37 38 ];
        bp_isprint=[0 1 0 1 0 1 0 0 0 1];
        temp=imread(fullfile(fileparts(which(mfilename)),sprintf('conn_vproject_icon%02d.jpg',0))); temp=double(temp); printmask=round(temp/255);
        for n1=1:numel(bp),
            set(data.handles(bp(n1)),'units','pixel'); pt=get(data.handles(bp(n1)),'position'); set(data.handles(bp(n1)),'units','norm');
            temp=imread(fullfile(fileparts(which(mfilename)),sprintf('conn_displayroi_icon%02d.jpg',n1))); temp=double(temp); temp=temp/255; temp=max(0,min(1,(temp).^.5)); ft=min(size(temp,1)/ceil(pt(4)),size(temp,2)/ceil(pt(3))); if any(n1==[1,2]), ft=0.95*ft; elseif any(n1==[9,10]), ft=.90*ft; else ft=.70*ft; end;
            maxtemp=1;%mode(round(temp(:)*100))/100;
            if maxtemp<.5, temp=1-temp; maxtemp=1-maxtemp; end
            temp=max(0,min(.95, .75*temp+.25*temp/maxtemp.*repmat(shiftdim(color1,-1),[size(temp,1),size(temp,2),1,size(temp,4)]) ));
            if bp_isprint(n1)
                if ismember(n1,[2,4,6,10]), tempprintmask=printmask(ceil(size(printmask,1)/4)+(1:ceil(size(printmask,1)/2)),ceil(size(printmask,2)/4)+(1:ceil(size(printmask,2)/2))); else tempprintmask=printmask; end
                temp=.75*mean(color1)+.25*temp;
                if size(temp,1)>size(temp,2), temp=temp(1:size(temp,2),:,:); end
                temp(:,ceil(size(temp,2)/2+(1:size(temp,1))-size(temp,1)/2),:)=max(0,temp(:,ceil(size(temp,2)/2+(1:size(temp,1))-size(temp,1)/2),:)+(1-2*mean(color1))*repmat(.5*tempprintmask(round(linspace(1,size(tempprintmask,1),size(temp,1))),round(linspace(1,size(tempprintmask,2),size(temp,1)))),[1,1,3]));
            end
            tempr1=round(1:ft/10:size(temp,1)); tempr1=tempr1(1:floor(numel(tempr1)/10)*10);
            tempr2=round(1:ft/10:size(temp,2)); tempr2=tempr2(1:floor(numel(tempr2)/10)*10);
            temp=permute(mean(mean(reshape(temp(tempr1,tempr2,:),[10,numel(tempr1)/10,10,numel(tempr2)/10,size(temp,3)]),1),3),[2,4,5,1,3]);
            str=get(data.handles(bp(n1)),'string'); set(data.handles(bp(n1)),'cdata',temp,'string','');
        end
        if data.thrtype==5, set(data.handles(22),'string',['connection threshold: ',data.statsname,' > ']);
        elseif data.mvpathr==1, set(data.handles(22),'string','connection threshold: ');
        else set(data.handles(22),'string','connection threshold: p < ');
        end
        set(data.handles([2,3,4,15,16,22,23]),'visible','off');
        hc1=uicontextmenu('parent',hfig);
        uimenu(hc1,'Label','Export table','callback',@(varargin)conn_exportlist(data.handles(8),'',get(data.handles(7),'string')),'tag','donotdelete');
        %hc2=uimenu(hc1,'Label','Sort rows by','tag','donotdelete');
        %uimenu(hc2,'Label','Connections','callback',@(varargin)conn_displayroi('mvpasort',1),'tag','donotdelete');
        %uimenu(hc2,'Label','Seeds','callback',@(varargin)conn_displayroi('mvpasort',2),'tag','donotdelete');
        %uimenu(hc2,'Label','Clusters','callback',@(varargin)conn_displayroi('mvpasort',3),'tag','donotdelete');
        set(data.handles(8),'uicontextmenu',hc1);
        %set(data.handles(7),'string',sprintf('%-6s %-6s  %6s  %6s  %4s  %8s  %8s','Seed','Target','beta',[data.MVPAstatsname,'/',data.statsname],'dof','p-unc','p-FDR'));
        %uicontrol('style','text','units','norm','position',[.53,.50,.45,.04],'string','second-level analysis results','foregroundcolor','b','backgroundcolor',color2,'fontname','monospaced','fontweight','bold','horizontalalignment','center');
        %uicontrol('style','text','units','norm','position',[.03,.01,.05,.04],'string','view:  ','foregroundcolor',1-color1,'backgroundcolor',color1,'fontweight','bold','horizontalalignment','right');
        %uicontrol('style','text','units','norm','position',[.35,.01,.02,.04],'string','q','fontname','symbol','horizontalalignment','right','foregroundcolor',1-color1,'backgroundcolor',color1,'fontweight','bold','horizontalalignment','right');
        %uicontrol('style','text','units','norm','position',[.20,.01,.02,.04],'string','z','horizontalalignment','right','foregroundcolor',1-color1,'backgroundcolor',color1,'fontweight','bold','horizontalalignment','right');
        rcircle=sign([sin(linspace(0,2*pi,64)'),cos(linspace(0,2*pi,64))'])*diag([5,5]);
        data.plotaxes=[];%axes('units','norm','position',[.01,.08,.58,.84],'visible','off','parent',data.hfig);
        data.legendaxes=axes('parent',data.hfig);set(data.legendaxes,'units','norm','position',[.90,.82,.09,.05],'xtick',[],'ytick',[],'box','on','xcolor',.5*[1,1,1],'ycolor',.5*[1,1,1]);axis(data.legendaxes,'equal');set(data.legendaxes,'xlim',[95-20,200]);axis(data.legendaxes,'off');
        cmap=jet(256);%cmap=cmap(32:224,:)*.8;
        %cmap=cmap.^repmat(.1+.9*abs(linspace(1,-1,size(cmap,1)))',1,size(cmap,2));
        data.legend=[patch(100+rcircle(:,1),0+rcircle(:,2),'w','edgecolor','none','facecolor',[1,.5,.5],'parent',data.legendaxes),...
                     patch(150+rcircle(:,1),0+rcircle(:,2),'w','edgecolor','none','facecolor',[.5,.5,1],'parent',data.legendaxes),...
                     text(0,0,'','horizontalalignment','left','fontsize',6+CONN_gui.font_offset,'color',.5*[1,1,1],'parent',data.legendaxes),... % ROI-to-ROI effects:
                     text(100+10,0,'Positive','horizontalalignment','left','fontsize',6+CONN_gui.font_offset,'color',.5*[1,1,1],'parent',data.legendaxes),...
                     text(150+10,0,'Negative','horizontalalignment','left','fontsize',6+CONN_gui.font_offset,'color',.5*[1,1,1],'parent',data.legendaxes),...
                     text(100-5,0,'Negative','horizontalalignment','left','fontsize',6+CONN_gui.font_offset,'color',.5*[1,1,1],'horizontalalignment','right','parent',data.legendaxes),...
                     text(150+5,0,'Positive','horizontalalignment','left','fontsize',6+CONN_gui.font_offset,'color',.5*[1,1,1],'horizontalalignment','left','parent',data.legendaxes)]; 
        for n=1:size(cmap,1), data.legend=[data.legend patch(100+(n+[0 0 1 1])*50/size(cmap,1),[-5,5,5,-5],'w','edgecolor','none','facecolor',cmap(n,:),'tag','conn_displayroi_plotlegendcont','parent',data.legendaxes)]; end        
        set(data.legend,'visible','off');
        set(data.legend,'buttondownfcn',@(varargin)conn_displayroi(hfig,[],'display.colorbar.limit'));
        set(hfig,'userdata',data);
        conn_displayroi(hfig,[],'fwec.option',[],'immediatereturn');
        data=get(hfig,'userdata');
        
        if data.displayreduced
            conn_displayroi(hfig,[],'displayreduced');
            set(hfig,'visible','on');
            if ishandle(hmsg), delete(hmsg); end
            if ishandle(hmsginit), delete(hmsginit); end
            return;
        else
            set(hfig,'visible','on');
            if ishandle(hmsg), delete(hmsg); end
            if ishandle(hmsginit), delete(hmsginit); end
        end
        
        
    case 'advancedthr'
        data=get(hfig,'userdata');
        if margin>1, advanced=varargin{1}; set(data.handles(6),'value',advanced);
        else advanced=get(data.handles(6),'value');
        end
        if advanced
            conn_displayroi(hfig,[],'fwec.option',numel(data.thres_defaults)+1,'immediatereturn');
        else
            ok=false;for method=1:numel(data.thres_defaults), if isequal({data.thr,data.thrtype,data.mvpathr,data.mvpathrtype},data.thres_defaults{method})&&data.side==3, ok=true; break; end; end;
            if ~ok,
                conn_displayroi(hfig,[],'fwec.option',numel(data.thres_defaults)+1,'immediatereturn');
            else
                conn_displayroi(hfig,[],'fwec.option',method,'immediatereturn');
                data.mvpathrtype=method;
            end
        end
        return
        
    case 'fwec.option'
        data=get(hfig,'userdata');
        if margin>1&&~isempty(varargin{1}), value=varargin{1}; set(data.handles(1),'value',value);
        else  value=get(data.handles(1),'value');
        end
        data.thres=value;
        advanced=get(data.handles(6),'value');
        if value<=numel(data.thres_defaults)
            [data.thr,data.thrtype,data.mvpathr,data.mvpathrtype]=deal(data.thres_defaults{data.thres}{:}); 
            data.side=3;
            set(data.handles(2),'string',num2str(data.thr));
            set(data.handles(3),'value',data.thrtype);
            set(data.handles(15),'string',num2str(data.mvpathr));
            set(data.handles(16),'value',find(data.mvpathrtype_shown==data.mvpathrtype));
            if data.thrtype==5, set(data.handles(22),'string',['connection threshold: ',data.statsname,' > ']);
            elseif data.mvpathr==1, set(data.handles(22),'string','connection threshold: ');
            else set(data.handles(22),'string','connection threshold: p < ');
            end
            switch(value)
                case 1, tstr3='parametric multivariate statistics (cluster-level inferences, Functional Network Connectivity, Jafri et al., 2008)';
                case 2, tstr3='non-parametric statistics (cluster-level inferences, Spatial Pairwise Clustering, Zalesky et al., 2012)';
                case 3, tstr3='non-parametric statistics (cluster-level inferences, Threshold Free Cluster Enhancement, Smith and Nichols 2007)';
                case 4, tstr3='parametric univariate statistics (connection-level inferences, FDR corrected, Benjamini & Hochberg, 1995)';
                case 5, tstr3='parametric multivariate statistics (ROI-level inferences, FDR corrected, Benjamini & Hochberg, 1995)';
                case 6, tstr3='non-parametric statistics (network-level inferences, Network Based Statistics, Zalesky et al., 2010)';
            end
            tstr1=cellstr(get(data.handles(3),'string'));
            tstr2=cellstr(get(data.handles(16),'string'));
            tstr4=sprintf('%s %s %s',get(data.handles(22),'string'),get(data.handles(2),'string'),tstr1{get(data.handles(3),'value')});
            tstr5=sprintf('%s %s %s',get(data.handles(23),'string'),get(data.handles(15),'string'),tstr2{get(data.handles(16),'value')});
            if value==3||value==4, set(data.handles(5),'string',{tstr3,sprintf('%s',tstr4)});
            else set(data.handles(5),'string',{tstr3,sprintf('%s; %s',tstr5,tstr4)});
            end
            if advanced,
                set(data.handles(5),'string','','visible','off');
                set(data.handles([2,3,4,15,16,22,23]),'visible','on');
                if data.thrtype==3||data.thrtype==4, set(data.handles([15,16,23]),'visible','off'); end
            else
                set(data.handles(5),'visible','on');
                set(data.handles([2,3,4,15,16,22,23]),'visible','off');
            end
        else
            set(data.handles(5),'string','','visible','off');
            set(data.handles([2,3,4,15,16,22,23]),'visible','on');
            if data.thrtype==3||data.thrtype==4, set(data.handles([15,16,23]),'visible','off'); end
        	set(hfig,'userdata',data); return;
        end
        if margin>2&&isequal(varargin{2},'immediatereturn'), set(hfig,'userdata',data); return; end

    case 'fwec.connectionlevel.value',
        data=get(hfig,'userdata');
        if margin>1, value=varargin{1}; set(data.handles(2),'string',num2str(value))
        else         value=str2num(get(data.handles(2),'string'));
        end
        if ~isempty(value), data.thr=max(0,value); end
        data.visible='on';
    case 'fwec.connectionlevel.type',
        data=get(hfig,'userdata');
        if margin>1, value=varargin{1}; set(data.handles(3),'value',value); 
        else value=max(1,min(5,get(data.handles(3),'value')));
        end
        docont=true;
        if data.thrtype==1&&value==5&&data.side==3, data.thr=spm_invTcdf(1-data.thr/2,data.dof(1)); set(data.handles(2),'string',num2str(data.thr)); docont=false;
        elseif data.thrtype==1&&value==5, data.thr=spm_invTcdf(1-data.thr,data.dof(1)); set(data.handles(2),'string',num2str(data.thr)); docont=false;
        elseif data.thrtype==5&&value==1&&data.side==3, data.thr=1-spm_Tcdf(data.thr,data.dof(1)); data.thr=2*min(data.thr,1-data.thr); set(data.handles(2),'string',num2str(data.thr)); docont=false;
        elseif data.thrtype==5&&value==1, data.thr=1-spm_Tcdf(data.thr,data.dof(1)); set(data.handles(2),'string',num2str(data.thr)); docont=false;
        elseif data.thrtype<3&&value==5&&data.thr<1, data.thr=3; set(data.handles(2),'string',num2str(data.thr)); 
        elseif data.thrtype==5&&value==2, data.thr=.05; set(data.handles(2),'string',num2str(data.thr)); 
        elseif data.thrtype==5&&value==1, data.thr=.001; set(data.handles(2),'string',num2str(data.thr)); 
        end
        data.thrtype=value;
        if data.thrtype==5, set(data.handles(22),'string',['connection threshold: ',data.statsname,' > ']);
        elseif data.mvpathr==1, set(data.handles(22),'string','connection threshold: ');
        else set(data.handles(22),'string','connection threshold: p < ');
        end        
        if data.thrtype==3||data.thrtype==4, set(data.handles([15,16,23]),'visible','off');
        else set(data.handles([15,16,23]),'visible','on');
        end
        data.visible='on';
        if ~docont, set(hfig,'userdata',data); return; end
    case 'fwec.connectionlevel.side',
        data=get(hfig,'userdata');
        if margin>1, value=varargin{1}; set(data.handles(4),'value',value);
        else value=get(data.handles(4),'value');
        end
        data.side=value;
        %data.plotconnoptions.LCOLOR=1+(data.side<3);
        data.visible='on';
    case 'fwec.clusterlevel.value',
        data=get(hfig,'userdata');
        if margin>1, value=varargin{1}; set(data.handles(15),'string',num2str(value));
        else
            str=get(data.handles(15),'string');
            value=str2num(str);
        end
        if ~isempty(value),
            data.mvpathr=max(0,value);
        end
        data.visible='on';
    case 'fwec.clusterlevel.type',
        data=get(hfig,'userdata');
        if margin>1, value=varargin{1}; set(data.handles(16),'value',value);
        else         value=get(data.handles(16),'value');
        end
        data.mvpathrtype=data.mvpathrtype_shown(value);
        data.visible='on';
        
    case 'plot_design'
        data=get(hfig,'userdata');
        Y=repmat({''},size(data.results(1).xX.X,1),size(data.results(1).c2,2));
        idx=find(data.results(1).xX.SelectedSubjects);
        if isfield(data.results(1),'ynames'), for n1=1:numel(idx), for n2=1:size(data.results(1).c2,2), Y{n1,n2}=sprintf('subject %d %s',idx(n1),data.results(1).ynames{n2}); end; end
        else for n1=1:numel(idx), for n2=1:size(data.results(1).c2,2), Y{n1,n2}=sprintf('subject %d condition %d',idx(n1),n2); end; end
        end
        conn_displaydesign(data.results(1).xX.X,Y,data.results(1).c,data.results(1).c2,data.results(1).xX.name,true);
        return
        
    case {'import_values','cluster_view'}
        data=get(hfig,'userdata');        
        filename={fullfile(data.defaultfilepath,'results.clusters.nii')};
        tfilename=conn_displayroi_selectfiles(filename,data.initfile);
        if isempty(tfilename), return; 
        elseif iscell(tfilename), CLUSTERSFROMFILE='';
        else CLUSTERSFROMFILE=tfilename; 
        end
        hmsginit=conn_msgbox('Loading data. Please wait...','conn_displayroi',-1);
        selectedsubjects=data.results(1).xX.SelectedSubjects;
        if isfield(data.results(1),'ynames'), names_conditions=data.results(1).ynames;
        else names_conditions={};
        end
       
        y={};name={};
        y2=[];name2={};
        y3=[];name3={};
        txt1='';
        txt2='';
        if ~isempty(CLUSTERSFROMFILE)
            [Cdata,Cnames,Ccoords,Csamples] = conn_mtx_read(CLUSTERSFROMFILE);
            [ok,Cidx]=ismember(Cnames,data.results(1).names2);
            if ~all(ok), fprintf('warning: unable to find ROI %s. Skipping\n',sprintf('%s ',Cnames{~ok})); Cdata(~ok,~ok)=0; end
            for value=1:max(Cdata(:))
                [isource,itarget]=find(Cdata==value);
                if isempty(y3), y2(end+1)=0; else y2(end+1)=size(y3,3); end
                if numel(Csamples)>=value, name2{end+1}=Csamples{value}; else name2{end+1}=sprintf('cluster %d',value); end
                for ni=1:numel(isource)
                    ty=permute(data.results(1).data(:,Cidx(itarget(ni)),:,Cidx(isource(ni))),[1,3,2]);   % (seed) [subjects x targets x conditions]
                    if ~isempty(selectedsubjects)&&~rem(size(ty,1),nnz(selectedsubjects)) % fill-in with NaN for missing data
                        tty=nan(size(ty,1)/nnz(selectedsubjects)*numel(selectedsubjects),size(ty,2));
                        tty(repmat(logical(selectedsubjects),size(ty,1)/nnz(selectedsubjects),1),:)=ty;
                        ty=tty;
                        if isempty(names_conditions), names_conditions=arrayfun(@(n)sprintf('condition %d',n),1:size(ty,2),'uni',0); end
                        y3=cat(3,y3,ty);
                        name3{end+1}=sprintf('%s connection between %s and %s',name2{end},Cnames{isource(ni)},Cnames{itarget(ni)});
                        for n=1:size(ty,2)
                            y{end+1}=ty(:,n);
                            name{end+1}=sprintf('%s at %s',name3{end},names_conditions{n});
                        end
                    end
                end
            end
        else
            for value=1:size(data.list2,1)
                isource=[];
                itarget=[];
                if data.list2(value,2)>0, % connection
                    txt2=sprintf('connectivity between %s and %s',data.names2{data.list2(value,1)},data.names2{data.list2(value,2)});
                    isource=data.list2(value,1);
                    itarget=data.list2(value,2);
                elseif data.list2(value,1)>0, % seed
                    txt1=data.list2txt{value};
                    if isempty(y3), y2(end+1)=0; else y2(end+1)=size(y3,3); end
                    name2{end+1}=regexp(txt1,'^(Cluster|Network|ROI) \d+\/\d+','match','once');
                else % cluster
                    txt1=data.list2txt{value};
                    if isempty(y3), y2(end+1)=0; else y2(end+1)=size(y3,3); end
                    name2{end+1}=regexp(txt1,'^(Cluster|Network|ROI) \d+\/\d+','match','once');
                end
                if ~isempty(isource)
                    ty=permute(data.results(1).data(:,itarget,:,isource),[1,3,2]);   % (seed) [subjects x targets x conditions]
                    if ~isempty(selectedsubjects)&&~rem(size(ty,1),nnz(selectedsubjects)) % fill-in with NaN for missing data
                        tty=nan(size(ty,1)/nnz(selectedsubjects)*numel(selectedsubjects),size(ty,2));
                        tty(repmat(logical(selectedsubjects),size(ty,1)/nnz(selectedsubjects),1),:)=ty;
                        ty=tty;
                        if isempty(names_conditions), names_conditions=arrayfun(@(n)sprintf('condition %d',n),1:size(ty,2),'uni',0); end
                        y3=cat(3,y3,ty);
                        name3{end+1}=sprintf('%s %s',regexp(txt1,'^(Cluster|Network|ROI) \d+\/\d+','match','once'),regexprep(txt2,'\s+',' '));
                        for n=1:size(ty,2)
                            y{end+1}=ty(:,n);
                            name{end+1}=sprintf('%s %s at %s',regexp(txt1,'^(Cluster|Network|ROI) \d+\/\d+','match','once'),regexprep(txt2,'\s+',' '),names_conditions{n});
                            %                         name{end+1}=regexprep(sprintf('conn between %s and %s at %s',...
                            %                             data.names2reduced{isource},...
                            %                             data.names2reduced{itarget},...
                            %                             names_conditions{n}),'\s+',' ');
                        end
                    end
                end
            end
        end
        if isempty(y3)
        elseif strcmpi(option,'import_values')
            if ~isfield(CONN_x,'filename')||isempty(CONN_x.filename)||~isfield(CONN_x,'Setup')||~isfield(CONN_x.Setup,'nsubjects')||isempty(CONN_x.Setup.nsubjects)||any(rem(cellfun(@numel,y),CONN_x.Setup.nsubjects)), % no conn project loaded
                try, [tfilename,tfilepath]=uiputfile('*.mat','Save data as',fileparts(tfilename));
                catch, [tfilename,tfilepath]=uiputfile('*.mat','Save data as');
                end
                if ischar(tfilename), conn_savematfile(fullfile(tfilepath,tfilename),'name','y'); fprintf('data saved to file %s\n',fullfile(tfilepath,tfilename)); end
            elseif get(data.handles(9),'value')||isempty(y2), % one variable per connection
                conn_importl2covariate(name,y);
            else % one variable per cluster
                breaks=[y2,size(y3,3)];
                y2=zeros([size(y3,1),size(y3,2),numel(breaks)-1]);
                for n=1:numel(breaks)-1,
                    y2(:,:,n)=mean(y3(:,:,breaks(n)+1:breaks(n+1)),3);
                end
                y3={};name3={};
                for n=1:size(y2,2),
                    y3=[y3,num2cell(permute(y2(:,n,:),[1,3,2]),1)];
                    name3=[name3, cellfun(@(s)sprintf('%s at %s',s,names_conditions{n}),name2,'uni',0)];
                end
                conn_importl2covariate(name3,y3);
            end
        else
            PLOTBETA=true; % true: displays size of regressor effects; false: displays size of contrast effects
            if get(data.handles(9),'value')||isempty(y2), % one plot per connection
                conn_rex('test',data.results(1).xX,reshape(y3(selectedsubjects,:,:),[nnz(selectedsubjects)*size(y3,2),size(y3,3)]),data.results(1).c,names_conditions,name3,[],[],true,data.results(1).c2,[],true,PLOTBETA);
            else % one plot per cluster
                breaks=[y2,size(y3,3)];
                y2=zeros([size(y3,1),size(y3,2),numel(breaks)-1]);
                for n=1:numel(breaks)-1,
                    y2(:,:,n)=mean(y3(:,:,breaks(n)+1:breaks(n+1)),3);
                end
                conn_rex('test',data.results(1).xX,reshape(y2(selectedsubjects,:,:),[nnz(selectedsubjects)*size(y2,2),size(y2,3)]),data.results(1).c,names_conditions,name2,[],[],true,data.results(1).c2,[],true,PLOTBETA);
            end
        end
        if ishandle(hmsginit), delete(hmsginit); end
        return
                             
    case 'export_data'
        DOSORT=true; % set to false to keep original ROI order; set to true to use ROI order from these analyses
        data=get(hfig,'userdata');
        if margin>1, tfilename=varargin{1};
        else
            [tfilename,tpathname]=uiputfile({'*.nii','NIFTI files (*.nii)'; '*',  'All Files (*)'},'Output data to file:',data.defaultfilepath);
            if ischar(tfilename)&&~isempty(tfilename),
                tfilename=fullfile(tpathname,tfilename);
            else return; 
            end
        end
        if isfield(data.results,'data'), Y=permute(data.results(1).data,[1,3,4,2]);
        else Y=permute(cat(4,data.results.y),[1,3,4,2]);
        end
        if DOSORT, idx=data.displaytheserois(data.displaytheserois<=size(Y,3));
        else idx=sort(data.displaytheserois(data.displaytheserois<=size(Y,3)));
        end
        Y=Y(:,:,idx,idx);
        Y=permute(Y,[3,4,1,2]); 
        ColumnNames=data.names2(idx);
        ColumnGroups=data.clusters(idx);
        
        SampleNames=repmat({''},size(Y,3),size(Y,4));
        idx=find(data.results(1).xX.SelectedSubjects);
        if isfield(data.results(1),'ynames'), for n1=1:numel(idx), for n2=1:size(Y,4), SampleNames{n1,n2}=sprintf('subject %d %s',idx(n1),data.results(1).ynames{n2}); end; end        
        else for n1=1:numel(idx), for n2=1:size(Y,4), SampleNames{n1,n2}=sprintf('subject %d measure %d',idx(n1),n2); end; end        
        end
        conn_mtx_write(tfilename,Y(:,:,:),ColumnNames, data.xyz2(idx), SampleNames);
        conn_disp('fprintf','Connectivity matrix data saved in %s\n',tfilename);
        return
        
    case {'export','export_mask'}
        data=get(hfig,'userdata');
        if margin>1, tfilename=varargin{1};
        else
            [tfilename,tpathname]=uiputfile({'*.nii','NIFI files (*.nii)'; '*.txt','text files (*.txt)'; '*.csv','CSV-files (*.csv)'; '*.mat','MAT-files (*.mat)'; '*',  'All Files (*)'},'Output mask to file:',fullfile(data.defaultfilepath,'results.nii'));
            if ischar(tfilename)&&~isempty(tfilename),
                tfilename=fullfile(tpathname,tfilename);
            else return
            end
        end
        if ~isempty(tfilename)
            [nill,nill,tfileext]=fileparts(tfilename);
            z=data.CM_z;
            R=z(data.displaytheserois(data.displaytheserois<=size(z,1)),data.displaytheserois);
            R(isnan(R))=0;
            z=data.CM_z0;
            R_unthresholded=z(data.displaytheserois(data.displaytheserois<=size(z,1)),data.displaytheserois);
            R_unthresholded(isnan(R))=0;
            ColumnNames=data.names2(data.displaytheserois);
            ColumnGroups=data.clusters(data.displaytheserois);
            R_clusters_unthresholded=data.CLUSTER_labels(data.displaytheserois(data.displaytheserois<=size(z,1)),data.displaytheserois);
            R_clusters=R_clusters_unthresholded;R_clusters(R==0)=0;
            %R_clusters=R_clusters_unthresholded;R_clusters(~ismember(R_clusters,unique(R_clusters(R~=0))))=0;
            
            if ~isempty(data.CLUSTER_selected)
                [ok,idx1]=ismember(R_clusters,data.CLUSTER_selected);
                R_clusters(~ok)=0;
                R_clusters(ok)=idx1(ok);
                R_clusters_names=data.CLUSTER_selected_names;
            else R_clusters_names={};
            end
            switch(tfileext)
                case '.nii'
                    conn_mtx_write(tfilename,R,ColumnNames, data.xyz2(data.displaytheserois));
                    conn_mtx_write(conn_prepend('',tfilename,'.orig.nii'),R_unthresholded,ColumnNames,data.xyz2(data.displaytheserois));
                    conn_mtx_write(conn_prepend('',tfilename,'.mask.nii'),double(R~=0),ColumnNames,data.xyz2(data.displaytheserois));
                    conn_disp('fprintf','Thresholded connectivity matrix saved in %s\n',tfilename);
                    if ~isempty(R_clusters_names), 
                        for n=1:numel(R_clusters_names)
                            idx=find(strncmp(data.list2txt(data.list2visible),R_clusters_names{n},numel(R_clusters_names{n})),1);
                            if ~isempty(idx),R_clusters_names{n}=data.list2txt{data.list2visible(idx)}; end
                        end
                        conn_mtx_write(conn_prepend('',tfilename,'.clusters.orig.nii'),R_clusters_unthresholded,ColumnNames,data.xyz2(data.displaytheserois),R_clusters_names); 
                        conn_mtx_write(conn_prepend('',tfilename,'.clusters.nii'),R_clusters,ColumnNames,data.xyz2(data.displaytheserois),R_clusters_names); 
                        conn_disp('fprintf','Suprathreshold connectivity clusters saved in %s\n',conn_prepend('',tfilename,'.clusters.nii'));
                    end
                case '.txt'
                    %fh=fopen(tfilename,'wt');
                    fh={};
                    for nt=1:numel(ColumnNames), fh{end+1}=sprintf('%s\t',ColumnNames{nt}); end; fh{end+1}=sprintf('\n');
                    for nt=1:numel(ColumnNames), for nt2=1:numel(ColumnNames), fh{end+1}=sprintf('%f\t',R(nt,nt2)); end; fh{end+1}=sprintf('\n'); end
                    %fclose(fh); 
                    conn_fileutils('filewrite_raw',tfilename, fh);
                    conn_disp('fprintf','Thresholded connectivity matrix saved in %s\n',tfilename);
                case '.csv'
                    %fh=fopen(tfilename,'wt');
                    fh={};
                    for nt=1:numel(ColumnNames), fh{end+1}=sprintf(',%s',ColumnNames{nt}); end; fh{end+1}=sprintf('\n');
                    for nt=1:numel(ColumnNames),
                        fh{end+1}=sprintf('%s',ColumnNames{nt});
                        for nt2=1:numel(ColumnNames), fh{end+1}=sprintf(',%f',R(nt,nt2)); end; fh{end+1}=sprintf('\n');
                    end
                    conn_fileutils('filewrite_raw',tfilename, fh);
                    %fclose(fh); 
                    conn_disp('fprintf','Thresholded connectivity matrix saved in %s\n',tfilename);
                case '.mat', 
                    conn_savematfile(tfilename,'R','R_unthresholded','ColumnNames','ColumnGroups'); 
                    conn_disp('fprintf','Thresholded connectivity matrix saved in %s\n',tfilename);
            end
        end
        return
        
    case 'openfolder'
        data=get(hfig,'userdata');
        conn_fileutils('cd',data.defaultfilepath);
        try
            if ispc, [nill,nill]=system(sprintf('start "%s"',data.defaultfilepath));
            else [nill,nill]=system(sprintf('open ''%s''',data.defaultfilepath));
            end
        end
        return
        
    case 'bookmark',
        data=get(hfig,'userdata');
        tfilename=[];
        descr='';
        if isfield(CONN_gui,'slice_display_skipbookmarkicons'), SKIPBI=CONN_gui.slice_display_skipbookmarkicons;
        else SKIPBI=false;
        end
        conn_args={'displayroi','initfile',data.initfile};
        opts={};%{'forcecd'};
        [fullfilename,tfilename,descr]=conn_bookmark('save',...
            tfilename,...
            descr,...
            conn_args,...
            opts);
        if isempty(fullfilename), return; end
        if ~SKIPBI,
            tht=conn_msgbox('Printing bookmark icon. Please wait...','',-1);
            conn_print(gcbf,conn_prepend('',fullfilename,'.jpg'),'-nogui','-r50','-nopersistent');
            if ishandle(tht), delete(tht); end
        end
        return;
    
    case 'mvpaextend'
        data=get(hfig,'userdata');
        if margin>1, value=varargin{1}; set(data.handles(16),'value',find(mvpathrtype_shown==value,1));
        else         value=data.mvpathrtype_shown(get(data.handles(16),'value'));
        end
        if isequal(data.mvpathrtype_shown,1:numel(data.mvpathrtype_all)), data.mvpathrtype_shown=[14:16,5:7,23:25,29:30,1]; %[15,24,6,1]; %[1,3,9,21]; 
        else data.mvpathrtype_shown=1:numel(data.mvpathrtype_all);
        end
        docont=false;
        if ~any(value==data.mvpathrtype_shown), value=data.mvpathrtype_shown(1); docont=true; end
        set(data.handles(16),'string',data.mvpathrtype_all(data.mvpathrtype_shown),'value',find(data.mvpathrtype_shown==value,1));
        data.mvpathrtype=value;
        if ~docont, set(hfig,'userdata',data); return; end
        data.visible='on';
%     case 'mvpaside',
%         hfig=gcbf;
%         data=get(hfig,'userdata');
%         value=get(data.handles(17),'value');
%         data.mvpaside=value;
%         data.visible='on';
    case 'mvpasort',
        data=get(hfig,'userdata');
        data.mvpasortresultsby=varargin{1};
        data.visible='on';
    case {'displayreduced','roi.select'}
        data=get(hfig,'userdata');
        if strcmpi(option,'roi.select'), data.displayreduced=2;
        elseif margin>1, data.displayreduced=varargin{1}; 
        end
        olddisplaytheserois=data.displaytheserois;
        oldclusters=data.clusters;
        switch(data.displayreduced),
            case 0,
                data.displaytheserois=1:length(data.names2);
            case 1,
                data.displaytheserois=1:length(data.names);
            case 2,
                if margin>2&&~isempty(varargin{2}),
                    initial=varargin{2};
                    listrois=[];
                    for n1=1:length(initial),
                        idx=strmatch(initial{n1},data.names,'exact');
                        if isempty(idx),
                            idx=strmatch(initial{n1},data.names); % allows partial-name matches
                        end
                        if isempty(idx), fprintf('warning: unable to find ROI %s. Skipping\n',initial{n1}); end
                        listrois=[listrois,idx(:)'];
                    end
                    if numel(listrois)>1, data.displaytheserois=listrois;
                    else
                        if numel(listrois)==1, conn_disp('Please select more than one ROI'); end
                        return;
                    end
                else
                    idxresortv=1:numel(data.names);
                    temp=regexp(data.names,'BA\.(\d*) \(L\)','tokens'); itemp=~cellfun(@isempty,temp); idxresortv(itemp)=-2e6+cellfun(@(x)str2double(x{1}),temp(itemp));
                    temp=regexp(data.names,'BA\.(\d*) \(R\)','tokens'); itemp=~cellfun(@isempty,temp); idxresortv(itemp)=-1e6+cellfun(@(x)str2double(x{1}),temp(itemp));
                    [nill,idxresort]=sort(idxresortv);
                    [nill,tidx]=ismember(data.displaytheserois,idxresort);
                    answ=listdlg('Promptstring','Select ROIs to include in this group-analysis:','selectionmode','multiple','liststring',data.names2(idxresort),'initialvalue',sort(tidx),'ListSize',[520 300]);
                    if numel(answ)>1, data.displaytheserois=sort(idxresort(answ));
                    else
                        if numel(answ)==1, conn_disp('Please select more than one ROI'); end
                        return;
                    end
                end
        end
        
        new1displaytheserois=olddisplaytheserois(ismember(olddisplaytheserois,data.displaytheserois)); % existing ones
        new2displaytheserois=data.displaytheserois(~ismember(data.displaytheserois,olddisplaytheserois)); % new ones
        data.displaytheserois=[reshape(new1displaytheserois,1,[]),reshape(new2displaytheserois,1,[])];
        data.source=data.source(data.source==0 | data.source<=length(data.displaytheserois));if isempty(data.source),data.source=1;end
        %results=conn_process('results_roi',data.displaytheserois);
        if margin<=1, h=conn_msgbox('Updating ROI-level results. Please wait...','conn_displayroi',-1); 
        else h=[]; 
        end
        for nresults=1:size(data.results(1).data,4), %numel(data.results) % note: obsolete?
            domvpa=data.displaytheserois;
            ndims=min(4,ceil(sqrt(size(data.results(1).data,1))/4));
            ndims=max(1,min(min(numel(domvpa),size(data.results(1).data,2)), ndims ));
            if ndims<numel(domvpa)
                y=data.results(1).data(:,domvpa,:,nresults);
                y(:,any(any(isnan(y),1),3),:)=[]; % subjects x rois x conditions
                sy=[size(y),1,1];
                y=reshape(permute(y,[1,3,2]),sy(1)*sy(3),sy(2)); % (subjects x conditions) x rois
                [Q,D,R]=svd(y,0);
                ndims=max(1,min(size(R,2),ndims));
                d=D(1:size(D,1)+1:size(D,1)*min(size(D))).^2;
                ndims=min([ndims,find(cumsum(d)/sum(d)>.95,1)]); % 95 percent variance
                y=y*R(:,1:ndims);
                MVPAy=permute(reshape(y,[sy(1),sy(3),ndims]),[1,3,2]);
                %data.results(nresults).MVPApcacov=d(1:ndims)/sum(d);
            else
                y=data.results(1).data(:,domvpa,:,nresults);
                y=y(:,~any(any(isnan(y),1),3),:);
                MVPAy=y;
                %data.results(nresults).MVPApcacov=[];
            end
            [dataresults(nresults).MVPAh,dataresults(nresults).MVPAF,dataresults(nresults).MVPAp,dataresults(nresults).MVPAdof,dataresults(nresults).MVPAstatsname]=conn_glm(data.results(1).xX.X,MVPAy(:,:),data.results(1).c,kron(data.results(1).c2,eye(size(MVPAy,2))));
            if isequal(dataresults(nresults).MVPAstatsname,'T'), 
                %data.results(nresults).MVPAstatsname='F'; data.results(nresults).MVPAdof=[1,data.results(nresults).MVPAdof]; data.results(nresults).MVPAF=data.results(nresults).MVPAF.^2;
                dataresults(nresults).MVPAp=2*min(dataresults(nresults).MVPAp,1-dataresults(nresults).MVPAp); 
            end
        end
        if ishandle(h), close(h); end
        data.MVPAF=cat(1,dataresults.MVPAF);
        data.MVPAp=cat(1,dataresults.MVPAp);
        temp={dataresults.MVPAdof};
        if any(cellfun('length',temp)>1), temp=cellfun(@(x)[ones(1,max(0,2-length(x))),x(:)'],temp,'uni',0); end
        data.MVPAdof=cell2mat(temp(:));
        %data.MVPAdof=cat(1,data.results.MVPAdof);
        data.MVPAstatsname=dataresults(1).MVPAstatsname;
        %data.MVPApcacov=cat(1,data.results.MVPApcacov);
        %if size(data.MVPAdof,2)>1&&~any(diff(data.MVPAdof(:,1))), data.MVPAstatsname=[data.MVPAstatsname,'(',num2str(data.MVPAdof(1)),')']; end
        if ~any(any(diff(data.MVPAdof,1,1),1),2),
            if size(data.MVPAdof,2)==1, data.MVPAdofstr=repmat({['(',num2str(data.MVPAdof(1)),')']},numel(data.MVPAdof),1);
            else data.MVPAdofstr=repmat({['(',num2str(data.MVPAdof(1,1)),',',num2str(data.MVPAdof(1,2)),')']},size(data.MVPAdof,1),1);
            end
        else
            if size(data.MVPAdof,2)==1, data.MVPAdofstr=arrayfun(@(dof)['(',num2str(dof(1)),')'],data.MVPAdof,'uni',0);
            else data.MVPAdofstr=cellfun(@(dof)['(',num2str(dof(1)),',',num2str(dof(2)),')'],num2cell(data.MVPAdof,2),'uni',0);
            end
        end
        %set(data.handles(7),'string',sprintf('%-6s %-6s  %6s  %6s  %4s  %8s  %8s','Seed','Target','beta',[data.MVPAstatsname,'/',data.statsname],'dof','p-unc','p-FDR'));
        if ~isequal(olddisplaytheserois,data.displaytheserois)||~isequal(oldclusters,data.clusters), %if ~isempty(new2displaytheserois) % skip recomputing clusters when only removing ROIs?
            data.xy2=zeros(length(data.names2),2); data.xy2(data.displaytheserois,:)=200*[cos(2*pi*(0:numel(data.displaytheserois)-1)'/numel(data.displaytheserois)),sin(2*pi*(0:numel(data.displaytheserois)-1)'/numel(data.displaytheserois))];
            data.xy2_clusters=[];
            data.clusters=[];
            data.names_clusters={};
            if isfield(data,'clusters_options')&&~isempty(data.clusters_options), data=conn_displayroi_clusters(data); end
        end
        if ~isequal(olddisplaytheserois,data.displaytheserois)||~isequal(oldclusters,data.clusters), 
            data.PERM=[]; 
            data.tfceZ=[];
            data.cMVPAF=[];
            %f=conn_dir(conn_displayroi_simfilename(data.roifile,'all'),'-R','-cell');
            %if ~isempty(f), conn_fileutils('spm_unlink',f{:}); end
        end
        if ishandle(h), close(h); end
        data.proj=[];data.x=[];data.y=[];data.z=[];
        data.bgz=0;
        data.visible='on';

    case {'roi.order.export','roi.order.save'}
        if margin>1&&~isempty(varargin{1}), answ=varargin{1};
        else answ=conn_questdlg('','Save ROI order:','Save ROI order/groups to file','Save ROI order/groups to clipboard','Save ROI order/groups to file');
        end
        try, answ=regexprep(answ,'order/cluster','order/group'); end
        data=get(hfig,'userdata');
        ROIconfiguration=struct('xy2',data.xy2,'displaytheserois',data.displaytheserois,'xy2_clusters',data.xy2_clusters,'clusters',data.clusters,'names2',{data.names2},'names_clusters',{data.names_clusters});
        if isequal(answ,'Save ROI order/groups to file')
            if margin>2&&~isempty(varargin{2}), tfilename=varargin{2};
            else
                [tfilename,tfilepath]=uiputfile('connROIorder.mat','Save ROI order as');
                if isequal(tfilename,0), return; end
                tfilename=fullfile(tfilepath,tfilename);
            end
            conn_savematfile(tfilename,'ROIconfiguration','-mat');
        elseif isequal(answ,'Save ROI order/groups to clipboard')
            assignin('base','ROIconfiguration',ROIconfiguration);
        end
        return
        
    case {'roi.order','roi.order.import','roi.order.importfromfile','roi.order.load','roi.order.loadfromfile'}
        data=get(hfig,'userdata');
        exstr='';
        Answ={'Use hierarchical clustering method (default)',...
            'Use CONN atlas apriori order/groups (atlas ROIs only)',...
            'Use CONN networks apriori order/groups (network ROIs only)',...
            'Load ROI order/groups from file',...
            'Load ROI order/groups from clipboard',...
            'Save ROI order/groups to file' ,...
            'Save ROI order/groups to clipboard',...
            'Manually define ROI order/groups'};
%             'Edit cluster labels'};
        option=regexprep(option,{'load/import|import/load|import','save/export|export/save|export'},{'load','save'},'ignorecase');
        if strcmpi(option,'roi.order.loadfromfile')
            varargin=[Answ(4),varargin{:}];
            margin=margin+1;
        end
        if margin>1, 
            answ=varargin{1};
            if isnumeric(answ), answ=Answ{answ}; end
        else
            if strcmpi(option,'roi.order')
                answ=conn_questdlg('','Define ROI order:',Answ{:},Answ{1});
            else
                answ=conn_questdlg('','Load ROI order:',Answ{2:5},Answ{2});
            end
        end
        if isequal(answ,Answ{1})
            conn_displayroi(hfig,[],'clusters','hc');
            return
        elseif isequal(answ,Answ{2})
            if ~conn_existfile(fullfile(fileparts(which('conn')),'rois','atlas.groups.mat')), conn_msgbox('Unable to find atlas.groups.mat file. Please update to latest release of CONN and try again','',2); return; 
            else
                load(fullfile(fileparts(which('conn')),'rois','atlas.groups.mat'),'ROIconfiguration','-mat');
                exstr='ROIs sorted using CONN atlas apriori order/groups';
                data.plotconnoptions.DOFFSET=.35;
            end
        elseif isequal(answ,Answ{3})
            if ~conn_existfile(fullfile(fileparts(which('conn')),'rois','networks.groups.mat')), conn_msgbox('Unable to find networks.groups.mat file. Please update to latest release of CONN and try again','',2); return; 
            else
                load(fullfile(fileparts(which('conn')),'rois','networks.groups.mat'),'ROIconfiguration','-mat');
                exstr='ROIs sorted using CONN networks apriori order/groups';
                data.plotconnoptions.DOFFSET=.70;
            end
        elseif isequal(answ,Answ{4})
            if margin>2&&~isempty(varargin{2}), tfilename=varargin{2};
            else 
                [tfilename,tfilepath]=conn_fileutils('uigetfile','connROIorder.mat','Load ROI order from');
                if isequal(tfilename,0), return; end
                tfilename=fullfile(tfilepath,tfilename);
            end
            if ~isempty(tfilename)
                ROIconfiguration=struct; conn_loadmatfile(tfilename,'ROIconfiguration','-mat');
                exstr='ROIs sorted manually (from file)';
            end
        elseif isequal(answ,Answ{5})
            try, ROIconfiguration=evalin('base','ROIconfiguration');
                exstr='ROIs sorted manually (from clipboard)';
            catch, conn_msgbox('Unable to import ROI configuration information. Please use ''Save ROI order/groups to clipboard'' first from this or a different ROI second-level results window','',2); return; 
            end
        elseif isequal(answ,Answ{6})||isequal(answ,Answ{7})
            conn_displayroi(hfig,[],'roi.order.save',answ,varargin{2:end});
            return
        elseif isequal(answ,Answ{8})
            ROIconfiguration=struct('xy2',data.xy2,'displaytheserois',data.displaytheserois,'xy2_clusters',data.xy2_clusters,'clusters',data.clusters,'names2',{data.names2},'names_clusters',{data.names_clusters});
            tfilename=conn_roiclusters(ROIconfiguration,[],[],[],fullfile(data.defaultfilepath,'ROIorder.mat'));
            if isempty(tfilename), return; end
            conn_loadmatfile(tfilename,'ROIconfiguration','-mat');
            exstr='ROIs sorted manually';
            %conn_displayroi(hfig,[],'display.groups.labels')
        else return;
        end
        olddisplaytheserois=data.displaytheserois;
        oldclusters=data.clusters;
        extnames2=ROIconfiguration.names2(ROIconfiguration.displaytheserois);
        ok=ismember(extnames2,data.names);
        for n1=reshape(find(~ok),1,[])
            idx=strmatch(extnames2{n1},data.names); % allows partial-name matches
            if numel(idx)==1, extnames2{n1}=data.names{idx}; end
        end
        [ok,idx]=ismember(data.names,extnames2);
        if ~nnz(ok), 
            conn_msgbox('Unable to import ROI configuration information. No matching ROIs','',2); return; end
        data.displaytheserois=find(ok);
        data.xy2(:)=nan;
        if isfield(ROIconfiguration,'xy2'), 
            data.xy2(data.displaytheserois,:)=ROIconfiguration.xy2(ROIconfiguration.displaytheserois(idx(ok)),:);
        else
            [idxok,tidx]=sort(idx(ok));
            clusters=1+cumsum(diff(ROIconfiguration.clusters(ROIconfiguration.displaytheserois(idxok)))~=0);
            a=-pi+2*pi*((1:numel(data.displaytheserois))'+[1;clusters(:)]-1)/(numel(data.displaytheserois)+max(clusters)); 
            data.xy2(data.displaytheserois(tidx),:)=200*[cos(a) sin(a)];
%             [clusters,tidx]=sort(ROIconfiguration.clusters(ROIconfiguration.displaytheserois(idx(ok))));
%             a=-pi+2*pi*((1:numel(data.displaytheserois))'+clusters-1)/(numel(data.displaytheserois)+max(clusters)); 
%             data.xy2(data.displaytheserois(tidx),:)=200*[cos(a) sin(a)];
        end
        data.clusters(:)=0;
        data.clusters(data.displaytheserois)=ROIconfiguration.clusters(ROIconfiguration.displaytheserois(idx(ok)));
        if isfield(ROIconfiguration,'xy2_clusters'), 
            data.xy2_clusters=data.xy2;
            data.xy2_clusters(data.displaytheserois,:)=ROIconfiguration.xy2_clusters(ROIconfiguration.displaytheserois(idx(ok)),:);
        elseif ~isempty(data.clusters) % note: automatic fill-in if missing
            data.xy2_clusters=data.xy2;
            mxy2_clusters=[]; for n1=1:max(data.clusters), if any(data.clusters==n1), mxy2_clusters(n1,:)=mean(data.xy2(data.clusters==n1,:),1); end; end
            data.xy2_clusters(data.clusters>0,:)=mxy2_clusters(data.clusters(data.clusters>0),:);
        else data.xy2_clusters=[];
        end
        if isfield(ROIconfiguration,'names_clusters')&&~isempty(ROIconfiguration.names_clusters), data.names_clusters=ROIconfiguration.names_clusters;
        else data.names_clusters={};
        end        
        if isfield(ROIconfiguration,'xyz2')&&~isempty(ROIconfiguration.xyz2), % note: automatic fill-in if missing
            data.xyz2(:)=nan;
            data.xyz2(data.displaytheserois,:)=ROIconfiguration.xyz2(ROIconfiguration.displaytheserois(idx(ok)),:);
        end
        %data.displaytheserois=ROIconfiguration.displaytheserois;
        %data.xy2=ROIconfiguration.xy2; 
        %data.clusters=ROIconfiguration.clusters;
        [nill,tidx]=sort(mod(pi+angle(data.xy2(data.displaytheserois,:)*[1;1i]),2*pi));data.displaytheserois=data.displaytheserois(tidx);
        data.proj=[];data.x=[];data.y=[];data.z=[];
        data.bgz=0;
        data.visible='on';
        try, if ~isempty(exstr), set(data.handles(14),'string',exstr); end; end

        
        if ~isequal(olddisplaytheserois,data.displaytheserois)||~isequal(oldclusters,data.clusters), 
            data.displayreduced=2;
            %set(data.handles(12),'value',data.displayreduced+1);
            data.PERM=[]; 
            data.tfceZ=[];
            data.cMVPAF=[];
            data.source=data.source(data.source==0 | data.source<=length(data.displaytheserois));if isempty(data.source),data.source=1;end
            %results=conn_process('results_roi_seed',data.displaytheserois);
            h=conn_msgbox('updating ROI-level results, please wait...','conn_displayroi');
            %f=conn_dir(conn_displayroi_simfilename(data.roifile,'all'),'-R','-cell');
            %if ~isempty(f), conn_fileutils('spm_unlink',f{:}); end
            for nresults=1:size(data.results(1).data,4), %numel(data.results) % note: obsolete?
                domvpa=data.displaytheserois;
                ndims=min(4,ceil(sqrt(size(data.results(1).data,1))/4)); % rule of thumb for number of dimensions
                ndims=max(1,min(min(numel(domvpa),size(data.results(1).data,2)), ndims ));
                if ndims<numel(domvpa)
                    y=data.results(1).data(:,domvpa,:,nresults);
                    y(:,any(any(isnan(y),1),3),:)=[]; % subjects x rois x conditions
                    sy=[size(y),1,1];
                    y=reshape(permute(y,[1,3,2]),sy(1)*sy(3),sy(2)); % (subjects x conditions) x rois
                    [Q,D,R]=svd(y,0);
                    ndims=max(1,min(size(R,2),ndims));
                    d=D(1:size(D,1)+1:size(D,1)*min(size(D))).^2;
                    ndims=min([ndims,find(cumsum(d)/sum(d)>.95,1)]); % 95 percent variance
                    y=y*R(:,1:ndims);
                    MVPAy=permute(reshape(y,[sy(1),sy(3),ndims]),[1,3,2]);
                    %data.results(nresults).MVPApcacov=d(1:ndims)/sum(d);
                else
                    y=data.results(1).data(:,domvpa,:,nresults);
                    y=y(:,~any(any(isnan(y),1),3),:);
                    MVPAy=y;
                    %data.results(nresults).MVPApcacov=[];
                end
                [dataresults(nresults).MVPAh,dataresults(nresults).MVPAF,dataresults(nresults).MVPAp,dataresults(nresults).MVPAdof,dataresults(nresults).MVPAstatsname]=conn_glm(data.results(1).xX.X,MVPAy(:,:),data.results(1).c,kron(data.results(1).c2,eye(size(MVPAy,2))));
                if isequal(dataresults(nresults).MVPAstatsname,'T'),
                    %data.results(nresults).MVPAstatsname='F'; data.results(nresults).MVPAdof=[1,data.results(nresults).MVPAdof]; data.results(nresults).MVPAF=data.results(nresults).MVPAF.^2;
                    dataresults(nresults).MVPAp=2*min(dataresults(nresults).MVPAp,1-dataresults(nresults).MVPAp);
                end
            end
            if ishandle(h), close(h); end
            data.MVPAF=cat(1,dataresults.MVPAF);
            data.MVPAp=cat(1,dataresults.MVPAp);
            temp={dataresults.MVPAdof};
            if any(cellfun('length',temp)>1), temp=cellfun(@(x)[ones(1,max(0,2-length(x))),x(:)'],temp,'uni',0); end
            data.MVPAdof=cell2mat(temp(:));
            %data.MVPAdof=cat(1,data.results.MVPAdof);
            data.MVPAstatsname=dataresults(1).MVPAstatsname;
            %data.MVPApcacov=cat(1,data.results.MVPApcacov);
            %if size(data.MVPAdof,2)>1&&~any(diff(data.MVPAdof(:,1))), data.MVPAstatsname=[data.MVPAstatsname,'(',num2str(data.MVPAdof(1)),')']; end
            if ~any(any(diff(data.MVPAdof,1,1),1),2),
                if size(data.MVPAdof,2)==1, data.MVPAdofstr=repmat({['(',num2str(data.MVPAdof(1)),')']},numel(data.MVPAdof),1);
                else data.MVPAdofstr=repmat({['(',num2str(data.MVPAdof(1,1)),',',num2str(data.MVPAdof(1,2)),')']},size(data.MVPAdof,1),1);
                end
            else
                if size(data.MVPAdof,2)==1, data.MVPAdofstr=arrayfun(@(dof)['(',num2str(dof(1)),')'],data.MVPAdof,'uni',0);
                else data.MVPAdofstr=cellfun(@(dof)['(',num2str(dof(1)),',',num2str(dof(2)),')'],num2cell(data.MVPAdof,2),'uni',0);
                end
            end
        end
        
    case {'ring_view','ring_print','default_print'}
        data=get(hfig,'userdata');
        if margin>1, options={varargin{1},'-nogui'};
        else options={fullfile(data.defaultfilepath,'print01.jpg')};
        end
        hc=data.plotconnoptions.BCOLOR; %[.95 .95 .9];
        hfig=figure('units','norm','position',[.4 .25 .6 .7],'color',hc,'menubar','none','name','matrix display','numbertitle','off');
        hax=copyobj(data.plotaxes,hfig);
        hax2=copyobj(data.legendaxes,hfig);
        set(hax,'units','norm','position',[.01,.05,.88,.9],'color',hc);
        %set(findobj(hax,'tag','plot_brainbackground'),'facecolor',hc);
        if strcmp(lower(option),'ring_print')||strcmp(lower(option),'default_print')
            conn_print(options{:});
            close(hfig);
        else
            set(hfig,'userdata',struct('hfig',hfig,'hax',hax));
            hc1=uimenu(hfig,'Label','Effects');
            hc2=uimenu(hc1,'Label','fontsize');
            uimenu(hc2,'Label','increase labels fontsize','callback','h=findobj(gcbf,''type'',''text''); s=max(1,cell2mat(get(h,''fontsize''))+1); for ss=reshape(unique(s),1,[]), set(h(s==ss),''fontsize'',ss); end');
            uimenu(hc2,'Label','decrease labels fontsize','callback','h=findobj(gcbf,''type'',''text''); s=max(1,cell2mat(get(h,''fontsize''))-1); for ss=reshape(unique(s),1,[]), set(h(s==ss),''fontsize'',ss); end');
            uimenu(hc2,'Label','set labels fontsize','callback','h=findobj(gcbf,''type'',''text''); s=conn_menu_inputdlg(''Enter fontsize'',''conn_displayroi'',1,{num2str(round(mean(cell2mat(get(h,''fontsize'')))))}); if ~isempty(s), s=str2num(s{1}); if ~isempty(s), set(h,''fontsize'',s); end; end');
            hc2=uimenu(hc1,'Label','background');
            uimenu(hc2,'Label','white background','callback','data=get(gcbf,''userdata''); h=findobj(gcbf,''type'',''text''); nc0=get(data.hfig,''color''); set([data.hfig data.hax],''color'',[1 1 1]); hc=cell2mat(get(h,''color'')); for nc=unique(hc,''rows'')'', idx=all(bsxfun(@eq,nc'',hc),2); set(h(idx),''color'',max(0,min(1,[1 1 1]-abs(nc''-nc0)))); end');
            uimenu(hc2,'Label','light background','callback','data=get(gcbf,''userdata''); h=findobj(gcbf,''type'',''text''); nc0=get(data.hfig,''color''); set([data.hfig data.hax],''color'',[.95 .95 .9]); hc=cell2mat(get(h,''color'')); for nc=unique(hc,''rows'')'', idx=all(bsxfun(@eq,nc'',hc),2); set(h(idx),''color'',max(0,min(1,[.95 .95 .9]-abs(nc''-nc0)))); end');
            uimenu(hc2,'Label','dark background','callback','data=get(gcbf,''userdata''); h=findobj(gcbf,''type'',''text''); nc0=get(data.hfig,''color''); set([data.hfig data.hax],''color'',[.11 .11 .11]); hc=cell2mat(get(h,''color'')); for nc=unique(hc,''rows'')'', idx=all(bsxfun(@eq,nc'',hc),2); set(h(idx),''color'',max(0,min(1,[.11 .11 .11]+abs(nc''-nc0)))); end');
            uimenu(hc2,'Label','black background','callback','data=get(gcbf,''userdata''); h=findobj(gcbf,''type'',''text''); nc0=get(data.hfig,''color''); set([data.hfig data.hax],''color'',[0 0 0]); hc=cell2mat(get(h,''color'')); for nc=unique(hc,''rows'')'', idx=all(bsxfun(@eq,nc'',hc),2); set(h(idx),''color'',max(0,min(1,[0 0 0]+abs(nc''-nc0)))); end');
            hc1=uimenu(hfig,'Label','Print');
            uimenu(hc1,'Label','current view','callback',@(varargin)conn_print(options{:}));
        end
        return

        
    case {'glass_view','glass_print','matrix_view','matrix_print','lines_view','lines_print','glass0_view','glass0_print'}
        if margin>1, options={varargin{1},'-nogui'};
        else         options={};
        end
        data=get(hfig,'userdata');
        datax=data.xyz2(:,1);%*data.proj(:,1);
        datay=data.xyz2(:,2);%*data.proj(:,2);
        dataz=data.xyz2(:,3);%*data.proj(:,3);
        weight1=cumpatch({hfig,1},'maskout');
        weight2=cumpatch({hfig,2},'maskout');
        weight3=cumpatch({hfig,3},'maskout');  
        if isempty(weight3), return; end
        idxkeep_source=weight3(1,weight3(2,:)>0);
        idxkeep=weight2(1,weight2(2,:)>0);
        idxkeep_conn1=weight1(1,weight1(end,:)>0);
        idxkeep_conn2=weight1(2,weight1(end,:)>0);
        if isempty(idxkeep), idxkeep=union(idxkeep_conn1,idxkeep_conn2); end
        ok=ismember(data.displaytheserois,idxkeep); idxkeep=data.displaytheserois(ok); % note: resorts
        c=repmat({[1,1,1]},[numel(idxkeep),1]);
        for n1=1:numel(idxkeep_source),idxc=find(idxkeep==idxkeep_source(n1));if ~isempty(idxc), c(idxc,:)=repmat({[.75,.75,.75]},[numel(idxc),1]); end; end
        %for n1=1:numel(data.source),idxc=find(idxkeep==data.source(n1));c(idxc,:)=repmat({[.25,.25,.25]},[numel(idxc),1]); end
        %c=mat2cell(data.plot_cmap(round(1+(size(data.plot_cmap,1)-1)*(1+data.plot_K(idxkeep))/2),:),ones(numel(idxkeep),1),3);
        %z=data.plot_z(data.displaytheserois(idxkeep(data.displaytheserois(idxkeep)<=size(data.plot_z,1))),data.displaytheserois(idxkeep)).*sign(data.h(data.displaytheserois(idxkeep(data.displaytheserois(idxkeep)<=size(data.plot_z,1))),data.displaytheserois(idxkeep)));
        %z(~ismember(idxkeep,idxkeep_source),:)=nan;
        % ring placeholder xyz/2
        if strcmp(lower(option),'matrix_view')||strcmp(lower(option),'matrix_print')
            z=zeros(numel(idxkeep));
            for n1=1:numel(idxkeep_conn1), 
                i=find(idxkeep==idxkeep_conn1(n1)); j=find(idxkeep==idxkeep_conn2(n1)); 
                z(i,j)=data.CM_z(idxkeep(i),idxkeep(j)).*sign(data.h(idxkeep(i),idxkeep(j))); 
            end
            z0=data.CM_z0(idxkeep,idxkeep);
            if data.issymmetric, 
                [i,j]=find(z~=0&z'==0);
                z(j+size(z,1)*(i-1))=z(i+size(z,1)*(j-1));
                [i,j]=find(z0~=0&z0'==0);
                z0(j+size(z0,1)*(i-1))=z0(i+size(z0,1)*(j-1));
            end
            znames=data.names2reduced(idxkeep);
            if max(data.clusters(idxkeep))>1, 
                tidx=find(diff(data.clusters(idxkeep))~=0);
                if isempty(data.names_clusters), toptions={tidx+.5, []}; 
                else toptions={tidx+.5, data.names_clusters(data.clusters(idxkeep([tidx(:)' numel(idxkeep)])))};
                end
            else toptions={[],[]};
            end
            if ~isempty(data.CLUSTER_selected)
                tcl=data.CLUSTER_labels(idxkeep,idxkeep);
                [ok,idx1]=ismember(tcl,data.CLUSTER_selected);
                tcl(~ok)=0;
                tcl(ok)=idx1(ok);
                toptions{end+1}=tcl;
                toptions{end+1}=data.CLUSTER_selected_names;
            else toptions=[toptions, {[],[]}];
            end
            toptions{end+1}=[];
            toptions{end+1}=znames;
            z(1:size(z,1)+1:numel(z))=nan;
            z0(1:size(z,1)+1:numel(z0))=nan;
            fh=conn_montage_display(...
                cat(4,z,z0),...
                {'suprathreshold group-level results','raw (unthresholded) group-level results'},...
                'matrix',...
                [],[],toptions{:});
            if ~isempty(data.maxz), fh('colorscale','rescale',data.maxz); end
            if ~isempty(regexp(lower(option),'print$')), 
                fh('background',[1 1 1]);
                fh('print',options{:}); fh('close'); 
            end
        elseif strcmp(lower(option),'lines_view')||strcmp(lower(option),'lines_print')
            z=zeros(numel(idxkeep));
            tcl=data.CLUSTER_labels(idxkeep,idxkeep);
            for n1=1:numel(idxkeep_conn1), 
                i=find(idxkeep==idxkeep_conn1(n1)); j=find(idxkeep==idxkeep_conn2(n1)); 
                z(i,j)=data.CM_z(idxkeep(i),idxkeep(j)).*sign(data.h(idxkeep(i),idxkeep(j))); 
            end
            tidx=idxkeep;
            if data.issymmetric, 
                [i,j]=find(z~=0&z'==0);
                z(j+size(z,1)*(i-1))=z(i+size(z,1)*(j-1));
                %z=triu(z);
                %p=amd(z~=0+eye(size(z,1))); %[p,q,r,s]=dmperm((z~=0)+eye(size(z,1)));
                %z=z(p,p); %z=triu(z(p,p));
                %tidx=idxkeep(p);
            end
            hstruct=conn_table_display(z, ...
                'colorsign',[],...
                'showtril',[],...
                'rlabel',data.names2reduced(tidx),...
                'clabel',data.names2reduced(tidx),...
                'rfontsize',data.plotconnoptions.FONTSIZE(1),...
                'cfontsize',data.plotconnoptions.FONTSIZE(1),...
                'vfontsize',max(1,data.plotconnoptions.FONTSIZE(1)-2),...
                'dfontsize',.5*data.plotconnoptions.FONTSIZE(1));
            if strcmp(lower(option),'lines_print')
                set([hstruct.hfig hstruct.hax],'color',[1 1 1]);
                conn_print(options{:});
                close(hstruct.hfig);
            end
        else % glass
            z=nan(numel(idxkeep));
            tcl=data.CLUSTER_labels(idxkeep,idxkeep);
            for n1=1:numel(idxkeep_conn1), 
                i=find(idxkeep==idxkeep_conn1(n1)); j=find(idxkeep==idxkeep_conn2(n1)); 
                z(i,j)=data.plot_z(idxkeep(i),idxkeep(j)).*sign(data.h(idxkeep(i),idxkeep(j))); 
            end
            fh=conn_mesh_display('','',[],...
                struct('sph_names',{data.names2reduced(idxkeep)},'sph_xyz',[datax(idxkeep),datay(idxkeep),dataz(idxkeep)],...
                'sph_r',3*ones(numel(idxkeep),1),...
                'sph_shapes',{data.names2(idxkeep)},...
                'sph_c',{c}),...%{repmat({[.9,.9,.9]},[1,numel(idxkeep)])}), ...
                z,... %+1*sign(z).*tcl, ...
                [], .2, [0,-1e-8,1],[],data.defaultfilepath);
            if strcmp(lower(option),'glass0_view')||strcmp(lower(option),'glass0_print')
                fh('brain',4);
                fh('brain_transparency',0);
                fh('sub_transparency',0);
                fh('mask_transparency',.05);
                %fh('material',[.1 1 1 .25 0]);
                fh('axis','on');
                try
                    nprojection=data.plotconnoptions.nprojection;
                    if nprojection<=0, nprojection=3; end
                    switch(nprojection)
                        case 1, fh('view',[-1,0,0]);
                        case 2, fh('view',[0,-1,0]);
                        case 3, fh('view',[0,-1e-8,1]);
                    end
                end
            else
                fh('brain',4);
                fh('brain_transparency',0);
                fh('sub_transparency',0);
                fh('mask_transparency',.15);
                fh('material',[]);
                fh('axis','on');
                fh('roi_color',.75*rand(numel(idxkeep),3));
                fh('view',[0,-1e-8,1]);
                fh('roi_shape','real');
                fh('roi_transparency',.15);
            end
            if ~isempty(regexp(lower(option),'print$')), 
                fh('background',[1 1 1]);
                fh('print',3,options{:}); fh('close'); 
            end
        end
        return
        
    case 'display.projection',
        data=get(hfig,'userdata');
        data.proj=varargin{1};
        data.x=[];data.y=[];data.z=[];
        data.bgz=0;
        
    case 'menubar',
        data=get(hfig,'userdata');
        data.plotconnoptions.menubar=data.plotconnoptions.menubar==0;
        
    case 'selectall',
        data=get(hfig,'userdata');
        data.source=0;
        data.bgz=0;
        %set(data.handles(11),'value',max(0,min(1,data.bgz/200+.5)));
        if data.view>0, data.bgimage=spm_get_data(data.ref,pinv(data.ref.mat)*[data.proj(:,1:3)*[data.bgx(:),data.bgy(:),data.bgz+zeros(prod(size(data.bgx)),1)]';ones(1,prod(size(data.bgx)))]); end
        
    case 'list1', % -obsolete-
        data=get(hfig,'userdata');
        if margin>1, value=varargin{1}; set(data.handles(6),'value',value);
        else         value=get(data.handles(6),'value');
        end
        if all(value>0&value<=length(data.displaytheserois)), 
            data.source=value;
            data.bgz=mean(data.z(data.displaytheserois(data.source)));
            %set(data.handles(11),'value',max(0,min(1,data.bgz/200+.5)));
            if data.view>0, data.bgimage=spm_get_data(data.ref,pinv(data.ref.mat)*[data.proj(:,1:3)*[data.bgx(:),data.bgy(:),data.bgz+zeros(prod(size(data.bgx)),1)]';ones(1,prod(size(data.bgx)))]); end
            %set(data.refaxes,'cdata',convn(convn(reshape(data.bgimage,size(data.bgx)),conn_hanning(5),'same'),conn_hanning(5)','same'));
        else
            data.source=0; 
            data.bgz=0;
            %set(data.handles(11),'value',max(0,min(1,data.bgz/200+.5)));
            if data.view>0, data.bgimage=spm_get_data(data.ref,pinv(data.ref.mat)*[data.proj(:,1:3)*[data.bgx(:),data.bgy(:),data.bgz+zeros(prod(size(data.bgx)),1)]';ones(1,prod(size(data.bgx)))]); end
            %set(data.refaxes,'cdata',convn(convn(reshape(data.bgimage,size(data.bgx)),conn_hanning(5),'same'),conn_hanning(5)','same'));
        end
        data.visible='on';
        
    case {'list2','list2clear'}
        data=get(hfig,'userdata');
        if strcmpi(option,'list2clear'),
            values=[];
            set(data.handles(8),'value',values);
        elseif margin>1, values=varargin{1}; set(data.handles(8),'value',values);
        else values=get(data.handles(8),'value');
        end
        if ~data.displayconnectionstats, 
            if all(values<=numel(data.list2visible)), values=data.list2visible(values); 
            else values=size(data.list2,1)+1;
            end
        end
        %if data.view==0, props={'facealpha',.10,'facealpha',min(.99,max(.001,abs(data.plotconnoptions.LTRANS)))};%{'facealpha',.10,'facealpha',min(.99,max(.001,abs(data.plotconnoptions.LTRANS)))};
        %else, props={'visible','off','visible','on'};
        %end
        if isempty(values)||any(values==0)||any(values>size(data.list2,1)), 
            if data.displayconnectionstats, values=1:size(data.list2,1); 
            else values=data.list2visible;
            end
        end
        if isempty(values)||any(values>size(data.list2,1))
            %set(cat(2,data.plotsadd2{cellfun('length',data.plotsadd2)>0}),props{3:4});%'linewidth',data.plotconnoptions.LINEWIDTH,'edgealpha',data.plotconnoptions.LTRANS);
            %set(cat(2,data.plotsadd3{cellfun('length',data.plotsadd3)>0}),props{3:4});
        else
            %set(cat(2,data.plotsadd2{:}),props{1:2});%'linewidth',data.plotconnoptions.LINEWIDTH/2,'edgealpha',data.plotconnoptions.LTRANS/2);
            %set(cat(2,data.plotsadd3{:}),props{1:2});
            mask=zeros(1,max(max(data.list2(:,1:2)))); % mask for ROIs
            maskc=[]; % mask for connections
            maskr=[]; % mask for highlighted seeds
            for value=values(:)'
                if value>0,%&&value<=length(data.plotsadd2),%&&~isempty(data.plotsadd2{value}),
                    if data.list2(value,2)>0,     mask(data.list2(value,1:2))=1; maskc=[maskc data.list2(value,1:2)];                   % individual connections
                    elseif data.list2(value,1)>0, mask(data.list2(value,1))=1; maskr=[maskr data.list2(value,1)];                       % seed
                    else
                        mask(data.list2(find(data.list2(:,3)==data.list2(value,3)&data.list2(:,1)>0),1))=1;
                        mask(data.list2(find(data.list2(:,3)==data.list2(value,3)&data.list2(:,2)>0),2))=1; 
                        maskc=[maskc reshape(data.list2(find(data.list2(:,3)==data.list2(value,3)&data.list2(:,1)>0&data.list2(:,2)>0),1:2)',1,[])];   % cluster
                    %else                          mask(data.list2(find(data.list2(:,3)==data.list2(value,3)&data.list2(:,1)>0),1))=1; maskr=[maskr reshape(data.list2(find(data.list2(:,3)==data.list2(value,3)&data.list2(:,1)>0),1),1,[])];   % network
                    end
                    %if data.list2(value,2)>0,     n2=find(data.list2(:,1)==data.list2(value,1)&data.list2(:,2)==data.list2(value,2)|data.list2(:,1)==data.list2(value,2)&data.list2(:,2)==data.list2(value,1));
                    %elseif data.list2(value,1)>0, n2=find(data.list2(:,1)==data.list2(value,1)|data.list2(:,2)==data.list2(value,1)); 
                    %else                          n2=find(data.list2(:,3)==data.list2(value,3));
                    %end
                    %set(cat(2,data.plotsadd2{n2}),props{3:4});%'linewidth',data.plotconnoptions.LINEWIDTH,'edgealpha',data.plotconnoptions.LTRANS);
                    %nall=unique(data.list2(n2,1:2));
                    %mask(nall(nall>0))=1;
                    %set(cat(2,data.plotsadd3{nall(nall>0)}),props{3:4});
                    if 0 % stat-color
                        cmap=jet(256);%cmap=cmap(32:224,:)*.8;
                        %cmap=cmap.^repmat(.1+.9*abs(linspace(1,-1,size(cmap,1)))',1,size(cmap,2));
                        n1n2=data.list2(n2,1:2);
                        n1n2valid=find(all(n1n2>0,2));
                        j=data.F(n1n2(n1n2valid,1)+size(data.F,1)*(n1n2(n1n2valid,2)-1));
                        J=max([eps;abs(j)]);
                        %for nt=1:numel(n1n2valid),
                        %    set(data.plotsadd2{nt},'edgecolor',cmap(ceil(size(cmap,1)/2)+round(floor(size(cmap,1)/2)*max(-1,min(1,j(nt)/J))),:));
                        %end
                    end
                end
            end
            if numel(values)==1&&all(values>0)
                if data.list2(value,2)>0, txt=sprintf('connectivity between %s and %s',data.names2{data.list2(value,1)},data.names2{data.list2(value,2)});
                elseif data.list2(value,1)>0, txt=sprintf('connectivity with %s',data.names2{data.list2(value,1)});
                else txt=sprintf('cluster comprising %d ROIs and %d connections among them',nnz(mask),ceil(nnz(data.list2(:,3)==data.list2(value,3)&data.list2(:,1)>0&data.list2(:,2)>0))); 
                %else txt=sprintf('cluster comprising %d ROIs and %d connections among them',nnz(mask),ceil(nnz(data.list2(:,3)==data.list2(value,3)&data.list2(:,1)>0&data.list2(:,2)>0)/2)); % note: display only (assuming bidirectional)
                end
                set(data.handles(21),'string',txt);
            else
                set(data.handles(21),'string','');
            end
            if isempty(maskc), 
                v=cumpatch({hfig,1},'mask',[1:numel(mask);  mask]); 
                w=zeros(size(v));
                if ~isempty(maskr)&&all(maskr<=numel(v)), v(maskr)=max(1,max(v(maskr))); w(maskr)=1; end 
                if isequal(size(v),size(w)), z=1i*w+v; else z=w; end
                cumpatch({hfig,2},'mask',[1:numel(mask); z]);
                cumpatch({hfig,3},'mask',[1:numel(mask); z]);
            else
                cumpatch({hfig,1},'mask',maskc);
                w=zeros(size(maskc));
                if ~isempty(maskr)&&all(maskr<=numel(mask)), mask(maskr)=max(1,max(mask(maskr))); w(maskr)=1; end
                %if isequal(size(v),size(w)), z=w+v; else z=w; end
                cumpatch({hfig,2},'mask',[1:numel(mask); mask]);
                cumpatch({hfig,3},'mask',[1:numel(mask); mask]);
            end
            
%             if 0
%                 n1=data.list2(value,1);
%                 n2=data.list2(value,2);
%                 y=cat(3,data.results(n1).y);
%                 if ~n2, y=y(:,data.displaytheserois);
%                     yc=data.names2(data.displaytheserois);
%                     MVPAy=cat(3,data.results(n1).MVPAy);
%                 else    y=y(:,n2);
%                     yc=data.names2(n2);
%                     MVPAy=[];
%                 end
%                 data.exploreplot.X=data.results(1).xX.X;
%                 data.exploreplot.c=data.results(1).c;
%                 data.exploreplot.c2=data.results(1).c2;
%                 data.exploreplot.effects=data.results(1).xX.X*data.results(1).c';
%                 data.exploreplot.y=y;
%                 data.exploreplot.y_fit=data.results(1).xX.X*(pinv(data.results(1).xX.X)*y);
%                 data.exploreplot.yc=yc;
%                 data.exploreplot.MVPAy=MVPAy;
%                 data.exploreplot.MVPAy_fit=data.results(1).xX.X*(pinv(data.results(1).xX.X)*MVPAy);
%                 %kron(data.results(nresults).c2,eye(size(data.results(nresults).MVPAy,2)))
%             end
        end
        %if get(data.handles(6),'listboxtop')>size(get(data.handles(6),'string'),1), set(data.handles(6),'listboxtop',1); end
        if get(data.handles(8),'listboxtop')>size(get(data.handles(8),'string'),1), set(data.handles(8),'listboxtop',1); end
        uicontrol(data.handles(8));
        return;
        data.visible='on';
    
    case {'view','view-axial','view-coronal','view-sagittal','view-ring'} % - obsolete - 
        data=get(hfig,'userdata');
        switch(lower(option))
            case 'view-ring', data.view=0; data.displaybrains=0;
            case 'view-axial', data.view=1; data.displaybrains=1;
            case 'view-coronal', data.view=2; data.displaybrains=1;
            case 'view-sagittal', data.view=3; data.displaybrains=1;
        end
        %data.view=get(data.handles(9),'value');
        data.proj=[];data.x=[];data.y=[];data.z=[];
        data.bgz=0;
%     case {'displayefffectsize-on','displayefffectsize-off','displayefffectsize-none'}
%         hfig=gcbf;
%         data=get(hfig,'userdata');
%         switch(lower(option))
%             case 'displayefffectsize-on', data.displayeffectsize=1;
%             case 'displayefffectsize-off', data.displayeffectsize=0;
%             case 'displayefffectsize-none', data.displayeffectsize=-1;
%         end
%         data.proj=[];data.x=[];data.y=[];data.z=[];
%         data.bgz=0;
    case 'display3d'
        data=get(hfig,'userdata');
        if 1, % ring placeholder 0
            data.view=1;
            data.proj=[];data.x=[];data.y=[];data.z=[];
            data.bgz=0;
        end
        data.display3d=1;
    case 'displaygui',
        data=get(hfig,'userdata');
        if (isempty(varargin)&&~data.displaygui)||(~isempty(varargin)&&strcmpi(varargin{1},'on'));
            if isfield(data,'displaygui_hcontrols'),
                data.displaygui=true;
                hcontrols=data.displaygui_hcontrols;
                set(hcontrols(ishandle(hcontrols)),'visible','on');
                set(data.plotaxes,'units','norm','position',data.plotposition{1});
            end
        else
            data.displaygui=false;
            hcontrols=data.handles(data.handles~=0);
            hcontrols=findobj(hcontrols,'flat','visible','on');
            hcontrols=hcontrols(ishandle(hcontrols));
            set(hcontrols,'visible','off');
            data.displaygui_hcontrols=hcontrols;
            set(data.plotaxes,'units','norm','position',data.plotposition{2});
        end            
%         data.proj=[];data.x=[];data.y=[];data.z=[];
%         data.bgz=0;
        set(hfig,'userdata',data);
        return
    case 'pausegui',
        data=get(hfig,'userdata');
        data.pausegui=1-(data.pausegui>0);
        set(hfig,'userdata',data);
        conn_display_windowbuttonmotionfcn('pause',data.pausegui);
        return
    case 'print'
        data=get(hfig,'userdata');
        if margin>1, options={varargin{1},'-nogui'};
        else options={fullfile(data.defaultfilepath,'print01.jpg')};
        end
        hcontrols=data.handles(data.handles~=0);
        hcontrols=findobj(hcontrols,'flat','visible','on');
        hcontrols=hcontrols(ishandle(hcontrols));
        set(hcontrols,'visible','off');
        set(data.plotaxes,'units','norm','position',data.plotposition{2});
        try, conn_print(options{:}); end
        set(hcontrols,'visible','on');
        set(data.plotaxes,'units','norm','position',data.plotposition{1});
        return
    case 'changebackground',
        data=get(hfig,'userdata');
        if margin>1, filename=varargin{1};
        else         filename=spm_select(1,'\.img$|\.nii$',['Select background anatomical image'],{},fileparts(data.ref.fname));
        end
        data.ref=spm_vol(filename);
        data.proj=[];data.x=[];data.y=[];data.z=[];
        data.bgz=0;
    case {'display.colorbar.limits','display.colorbar.limit'}
        data=get(hfig,'userdata');
        if numel(varargin)>=1, data.maxz=max(abs(varargin{1}));
        else 
            val=data.maxz;
            val=conn_menu_inputdlg({'Enter new colorbar limit:'},'Rescale colorbar',1,{num2str(val)});
            if isempty(val), return; end
            val=str2num(val{1}); 
            data.maxz=max(abs(val));
        end
        data.x=[];data.y=[];data.z=[];
        data.bgz=0;        
    case 'refresh',
        data=get(hfig,'userdata');
        data.x=[];data.y=[];data.z=[];
        data.bgz=0;
    case 'display.viewcycle'
        data=get(hfig,'userdata');
        data.plotconnoptions.nprojection=1+mod(data.plotconnoptions.nprojection,3);
    case {'displaytype','display.type'}
        value=varargin{1};
        %if numel(varargin)>=2, hfig=varargin{1}; value=varargin{2};
        %else hfig=gcbf; value=varargin{1};
        %end
        %if isempty(hfig), hfig=gcf; end
        if ischar(value),value=str2num(value); end
        if isempty(value), return; end
        data=get(hfig,'userdata');
        data.plotconnoptions.LINESTYLEMTX=value;
        data.x=[];data.y=[];data.z=[];
        data.bgz=0;
    case 'display.labelsize'
        data=get(hfig,'userdata');
        data.plotconnoptions.DOFFSET=varargin{1};
        initxy=true;        
    case 'display.fontsize'
        data=get(hfig,'userdata');
        data.plotconnoptions.FONTSIZE=varargin{1};
    case 'display.backgroundcolor'
        data=get(hfig,'userdata');
        data.plotconnoptions.BCOLOR=varargin{1};
    case {'displayoptions','display.options'}
        data=get(hfig,'userdata');
        answer=conn_menu_inputdlg({'Brain display size (0,inf)','Brain display orientation (0=automatic; 1=sagittal; 2=coronal; 3=axial)','Brain display contrast (0,1)','Connectivity display (0: lines; 1: matrix)','Connectivity lines width (-inf,inf; negative for proportional to stats)','Connectivity lines transparency (-1,1; negative for proportional to stats)','Connectivity lines curvature (-inf,inf)','Connectivity lines bundling (0,inf)','ROI sphere size (0,inf)','Space reserved for labels (0,inf)','Fontsize for labels (pts)','Rotation of labels (degrees)','Image background color (rgb)'},'display options',1,...
            {num2str(data.plotconnoptions.BSCALE),num2str(data.plotconnoptions.nprojection),num2str(data.plotconnoptions.BTRANS),num2str(data.plotconnoptions.LINESTYLEMTX),num2str(data.plotconnoptions.LINEWIDTH),num2str(data.plotconnoptions.LTRANS),num2str(data.plotconnoptions.LCURVE),num2str(data.plotconnoptions.LBUNDL),num2str(data.plotconnoptions.RSCALE),num2str(data.plotconnoptions.DOFFSET),num2str(data.plotconnoptions.FONTSIZE),num2str(data.plotconnoptions.FONTANGLE),num2str(data.plotconnoptions.BCOLOR)});
        try
            data.plotconnoptions.BSCALE=str2num(answer{1});
            data.plotconnoptions.nprojection=str2num(answer{2});
            data.plotconnoptions.BTRANS=str2num(answer{3});
            data.plotconnoptions.LINESTYLEMTX=str2num(answer{4});
            data.plotconnoptions.LINEWIDTH=str2num(answer{5});
            data.plotconnoptions.LTRANS=str2num(answer{6});
            data.plotconnoptions.LCURVE=str2num(answer{7});
            data.plotconnoptions.LBUNDL=str2num(answer{8});
            data.plotconnoptions.RSCALE=str2num(answer{9});
            data.plotconnoptions.DOFFSET=str2num(answer{10});
            data.plotconnoptions.FONTSIZE=str2num(answer{11});
            data.plotconnoptions.FONTANGLE=str2num(answer{12});
            data.plotconnoptions.BCOLOR=str2num(answer{13});
            initxy=true;
        catch
        end
    case 'enablethr' % - obsolete - 
        data=get(hfig,'userdata');
        value=varargin{1};
        data.enablethr=value;
        if ischar(data.enablethr), data.enablethr=strcmpi(data.enablethr,'on'); end
    case 'mvpaenablethr' % - obsolete - 
        value=varargin{1};
        data=get(hfig,'userdata');
        data.mvpaenablethr=value;
        if ischar(data.mvpaenablethr), data.mvpaenablethr=strcmpi(data.mvpaenablethr,'on'); end
%     case {'mvpaenablethr','enablethr'}
%         hfig=gcbf;
%         data=get(hfig,'userdata');
%         data.enablethr=get(data.handles(19),'value');
%         data.mvpaenablethr=get(data.handles(18),'value');
    case 'displayroilabelstats'
        data=get(hfig,'userdata');
        if margin>1, value=varargin{1}; set(data.handles(13),'value',value);
        else value=get(data.handles(13),'value');
        end
        data.displayroilabelsinstats=value;        
        data.displayconnectionstats=1;
        if data.displayroilabelsinstats, txt=[regexprep(data.list2txt,'\\\\(\d*)\\\\(.*?)\\\\','$2') {' '}];
        else txt=[regexprep(data.list2txt,'\\\\(\d*)\\\\(.*?)\\\\','$1') {' '}];
        end
        set(data.handles(8),'string',txt ,'value',max(1,min(numel(data.list2txt)+1, get(data.handles(8),'value'))));
        if get(data.handles(8),'listboxtop')>numel(txt), set(data.handles(8),'listboxtop',1); end
        set(hfig,'userdata',data);
        return
    case 'displayconnectionstats'
        data=get(hfig,'userdata');
        if margin>1, value=varargin{1}; set(data.handles(9),'value',value);
        else         value=get(data.handles(9),'value');
        end
        str0=cellstr(get(data.handles(8),'string'));
        val0=get(data.handles(8),'value');
        data.displayconnectionstats=value;        
        if data.displayconnectionstats,
            if data.displayroilabelsinstats, txt=[regexprep(data.list2txt,'\\\\(\d*)\\\\(.*?)\\\\','$2') {' '}];
            else txt=[regexprep(data.list2txt,'\\\\(\d*)\\\\(.*?)\\\\','$1') {' '}];
            end
            set(data.handles(8),'string',txt ,'value',max(1,min(numel(data.list2txt)+1, get(data.handles(8),'value'))));
            try, [ok,idx]=ismember(str0(val0),txt(1:end-1)); if any(ok), set(data.handles(8),'value',idx(ok)); end; end
            set(data.handles(13),'visible','on');
        else
            if data.displayroilabelsinstats, txt=[regexprep(data.list2txt(data.list2visible),'\\\\(\d*)\\\\(.*?)\\\\','$2') {' '}];
            else txt=[regexprep(data.list2txt(data.list2visible),'\\\\(\d*)\\\\(.*?)\\\\','$1') {' '}];
            end
            set(data.handles(8),'string',txt ,'value',max(1,min(numel(data.list2visible)+1, get(data.handles(8),'value'))));
            try, [ok,idx]=ismember(str0(val0),txt(1:end-1)); if any(ok), set(data.handles(8),'value',idx(ok)); end; end
            set(data.handles(13),'visible','off');
        end
        if get(data.handles(8),'listboxtop')>numel(txt), set(data.handles(8),'listboxtop',1); end
        set(hfig,'userdata',data);
        return
    case {'edgecolors1','edgecolors2','edgecolors3','edgecolors4','edgecolors5'}
        data=get(hfig,'userdata');
        if strcmp(option,'edgecolors1'), data.plotconnoptions.LCOLOR=1;
        elseif strcmp(option,'edgecolors2'),  data.plotconnoptions.LCOLOR=2;
        elseif strcmp(option,'edgecolors3'),  data.plotconnoptions.LCOLOR=3;
        elseif strcmp(option,'edgecolors4'),  data.plotconnoptions.LCOLORSCALE=data.plotconnoptions.LCOLORSCALE/1.5;
        elseif strcmp(option,'edgecolors5'),  data.plotconnoptions.LCOLORSCALE=data.plotconnoptions.LCOLORSCALE*1.5;
        end
    case {'edgewidths1','edgewidths2','edgewidths3','edgewidths4'}
        data=get(hfig,'userdata');
        if strcmp(option,'edgewidths1'), data.plotconnoptions.LINEWIDTH=abs(data.plotconnoptions.LINEWIDTH);
        elseif strcmp(option,'edgewidths2'), data.plotconnoptions.LINEWIDTH=-abs(data.plotconnoptions.LINEWIDTH);
        elseif strcmp(option,'edgewidths3'), data.plotconnoptions.LINEWIDTH=2*data.plotconnoptions.LINEWIDTH;
        elseif strcmp(option,'edgewidths4'), data.plotconnoptions.LINEWIDTH=1/2*data.plotconnoptions.LINEWIDTH;
        end
    case {'edgeopacity1','edgeopacity2','edgeopacity3','edgeopacity4'}
        data=get(hfig,'userdata');
        if strcmp(option,'edgeopacity1'), data.plotconnoptions.LTRANS=abs(data.plotconnoptions.LTRANS);
        elseif strcmp(option,'edgeopacity2'), data.plotconnoptions.LTRANS=-abs(data.plotconnoptions.LTRANS);
        elseif strcmp(option,'edgeopacity3'), data.plotconnoptions.LTRANS=2*data.plotconnoptions.LTRANS;
        elseif strcmp(option,'edgeopacity4'), data.plotconnoptions.LTRANS=1/2*data.plotconnoptions.LTRANS;
        end
    case {'labelsoff','labelson','labelspartial','labelsfull'}
        data=get(hfig,'userdata'); 
        data.displaylabels=2*strcmpi(option,'labelsfull')+strcmpi(option,'labelson')+.5*strcmpi(option,'labelspartial');
        if data.displaylabels==0, set(findobj(hfig,'tag','textstring','-or','tag','textstringpartial'),'visible','off'); 
        elseif data.displaylabels==1, set(findobj(hfig,'tag','textstring','-or','tag','textstringpartial'),'visible','on');
        elseif data.displaylabels>1, set(findobj(hfig,'tag','textstring','-or','tag','textstringpartial'),'visible','on','color',1-data.plotconnoptions.BCOLOR);
        else set(findobj(hfig,'tag','textstringpartial'),'visible','off');set(findobj(hfig,'tag','textstring'),'visible','on');
        end
        set(hfig,'userdata',data); 
        return
    case {'labels1','labels2','labels3'}
        data=get(hfig,'userdata');
        if strcmp(option,'labels1'), data.plotconnoptions.FONTSIZE=data.plotconnoptions.FONTSIZE+1;
        elseif strcmp(option,'labels2'), data.plotconnoptions.FONTSIZE=max(1,data.plotconnoptions.FONTSIZE-1);
        elseif strcmp(option,'labels3'), s=conn_menu_inputdlg('Enter fontsize','conn_displayroi',1,{num2str(data.plotconnoptions.FONTSIZE(1))}); if ~isempty(s), s=str2num(s{1}); if ~isempty(s), data.plotconnoptions.FONTSIZE(1)=s(1); end; end;
        end
        h=findobj(hfig,'tag','textstring','-or','tag','textstringpartial');
        if ~isempty(h), set(h,'fontsize',data.plotconnoptions.FONTSIZE(1)); end
        set(hfig,'userdata',data);
        return
    case {'labelsedit','display.rois.labels'}
        data=get(hfig,'userdata');
        name=data.names2reduced;
        %name=get(state.handles.sphplots_txt,'string');
        ok=true;
        thfig=dialog('units','norm','position',[.3,.4,.6,.4],'windowstyle','normal','name','ROI labels','color','w','resize','on');
        uicontrol(thfig,'style','text','units','norm','position',[.1,.85,.8,.10],'string',sprintf('New ROI label names (%d)',numel(name)),'backgroundcolor','w','fontsize',9+CONN_gui.font_offset,'fontweight','bold');
        ht1=uicontrol(thfig,'style','edit','units','norm','position',[.1,.30,.8,.55],'max',2,'string',name,'fontsize',8+CONN_gui.font_offset,'horizontalalignment','left','tooltipstring','manually edit the ROI labels');
        ht2=uicontrol(thfig,'style','edit','units','norm','position',[.1,.20,.8,.1],'max',1,'string','','fontsize',8+CONN_gui.font_offset,'horizontalalignment','left','tooltipstring','enter Matlab command for fast-editing all ROIs simultaneously (str is input variable cell array; ouput is cell array; e.g. "lower(str)")','callback','ht1=get(gcbo,''userdata''); set(ht1,''string'',feval(inline(get(gcbo,''string''),''str''),get(ht1,''string'')))','userdata',ht1);
        uicontrol(thfig,'style','pushbutton','string','Apply','units','norm','position',[.1,.01,.38,.10],'callback','uiresume','fontsize',8+CONN_gui.font_offset);
        uicontrol(thfig,'style','pushbutton','string','Cancel','units','norm','position',[.51,.01,.38,.10],'callback','delete(gcbf)','fontsize',8+CONN_gui.font_offset);
        while ok
            uiwait(thfig);
            ok=ishandle(thfig);
            if ok,
                newname=get(ht1,'string');
                if numel(newname)~=numel(name), conn_msgbox(sprintf('Number of labels entered (%d) does not match expected value (%d)',numel(newname),numel(name)),'',2);
                else
                    delete(thfig);
                    data.names2reduced=newname;
                    ok=false;
                end
            else return;
            end
        end
    case {'groupsedit','display.groups.labels'}
        data=get(hfig,'userdata');
        name=data.names_clusters;
        if isempty(name), 
            name=arrayfun(@(n)sprintf('Group #%d (%s)',n,sprintf('%s ',data.names2reduced{data.clusters==n})),1:max(data.clusters),'uni',0); 
            %name=regexprep(name,'^(.{128})(.+)$','$1 ...)');
        end
        %name=get(state.handles.sphplots_txt,'string');
        ok=true;
        thfig=dialog('units','norm','position',[.3,.4,.6,.4],'windowstyle','normal','name','Group labels','color','w','resize','on');
        uicontrol(thfig,'style','text','units','norm','position',[.1,.85,.8,.10],'string',sprintf('New Group label names (%d)',numel(name)),'backgroundcolor','w','fontsize',9+CONN_gui.font_offset,'fontweight','bold');
        ht1=uicontrol(thfig,'style','edit','units','norm','position',[.1,.30,.8,.55],'max',2,'string',name,'fontsize',8+CONN_gui.font_offset,'horizontalalignment','left','tooltipstring','manually edit the Group labels');
        ht2=uicontrol(thfig,'style','edit','units','norm','position',[.1,.20,.8,.1],'max',1,'string','','fontsize',8+CONN_gui.font_offset,'horizontalalignment','left','tooltipstring','enter Matlab command for fast-editing all Groups simultaneously (str is input variable cell array; ouput is cell array; e.g. "lower(str)")','callback','ht1=get(gcbo,''userdata''); set(ht1,''string'',feval(inline(get(gcbo,''string''),''str''),get(ht1,''string'')))','userdata',ht1);
        uicontrol(thfig,'style','pushbutton','string','Apply','units','norm','position',[.1,.01,.24,.10],'callback','uiresume','fontsize',8+CONN_gui.font_offset);
        uicontrol(thfig,'style','pushbutton','string','Delete all','units','norm','position',[.38,.01,.24,.10],'fontsize',8+CONN_gui.font_offset,'callback','ht1=get(gcbo,''userdata''); set(ht1,''string'',regexprep(get(ht1,''string''),''.*'',''''))','userdata',ht1);
        uicontrol(thfig,'style','pushbutton','string','Cancel','units','norm','position',[.66,.01,.24,.10],'callback','delete(gcbf)','fontsize',8+CONN_gui.font_offset);
        while ok
            uiwait(thfig);
            ok=ishandle(thfig);
            if ok,
                newname=get(ht1,'string');
                if numel(newname)~=numel(name), conn_msgbox(sprintf('Number of labels entered (%d) does not match expected value (%d)',numel(newname),numel(name)),'',2);
                else
                    delete(thfig);
                    data.names_clusters=newname;
                    ok=false;
                end
            else return;
            end
        end
        initxy=true;
       
    case {'brainsoff','brainson','brainssingle','brainsmany'}
        data=get(hfig,'userdata');
        switch(lower(option))
            case 'brainson', data.displaybrains=2;
            case 'brainsmany', data.displaybrains=1.5;
            case 'brainssingle', data.displaybrains=1;
            otherwise, data.displaybrains=0;
        end
        initxy=true;
%         data.displaybrains=2*strcmpi(option,'brainson')+strcmpi(option,'brainspartial');
%     case 'slider1',
%         hfig=gcbf;
%         data=get(hfig,'userdata');
%         value=get(data.handles(10),'value');
%         ang=(value-.5)*pi;
%         switch(data.view),
%             case 1, data.proj=[1,0,0;0,cos(ang),-sin(ang);0,sin(ang),cos(ang)];
%             case 2, data.proj=[1,0,0;0,sin(ang),-cos(ang);0,cos(ang),sin(ang)];
%             case 3, data.proj=[0,sin(ang),cos(ang);1,0,0;0,cos(ang),-sin(ang)];
%         end
%         data.x=[];data.y=[];data.z=[];
    case 'slider2',
        data=get(hfig,'userdata');
        if margin>1, value=varargin{1}; set(data.handles(11),'value',value);
        else         value=get(data.handles(11),'value');
        end
        data.bgz=(value-.5)*200;
        data.bgimage=spm_get_data(data.ref,pinv(data.ref.mat)*[data.proj(:,1:3)*[data.bgx(:),data.bgy(:),data.bgz+zeros(prod(size(data.bgx)),1)]';ones(1,prod(size(data.bgx)))]);
        set(data.refaxes,'cdata',reshape(data.bgimage,size(data.bgx)));%convn(convn(reshape(data.bgimage,size(data.bgx)),conn_hanning(5),'same'),conn_hanning(5)','same'));
        set(hfig,'userdata',data);
        return;

    case 'close'
        delete(hfig);
        return
        
    case 'clusters',
        options=varargin;
        data=get(hfig,'userdata');
        olddisplaytheserois=data.displaytheserois;
        oldclusters=data.clusters;
        data=conn_displayroi_clusters(data,options{:});
        if ~isequal(olddisplaytheserois,data.displaytheserois)||~isequal(oldclusters,data.clusters), 
            data.PERM=[]; 
            data.tfceZ=[];
            data.cMVPAF=[];
        end
        data.proj=[];data.x=[];data.y=[];data.z=[];
        data.bgz=0;
    otherwise
        conn_disp('fprintf','unrecognized conn_displayroi option %s\n',lower(option));
        return
%     case 'cluster',
%         hfig=gcf;
%         data=get(hfig,'userdata');
%         ncluster=varargin{1};
%         data.display='connectivity';
%         data.visible='on';
%         data.source=find(data.clusters==ncluster);
%         data.bgz=mean(data.z(data.source));
%         %set(data.handles(11),'value',max(0,min(1,data.bgz/200+.5)));
%         data.bgimage=spm_get_data(data.ref,pinv(data.ref.mat)*[data.proj(:,1:3)*[data.bgx(:),data.bgy(:),data.bgz+zeros(prod(size(data.bgx)),1)]';ones(1,prod(size(data.bgx)))]);
%         set(data.handles(6),'value',data.source);
%         %conn_disp(['Cluster #',num2str(ncluster)]);conn_disp(strvcat(data.names{data.source}));
end






% selects optimal view
if isempty(data.view),
    data.view=0;
end
% projector associated with view
if isempty(data.proj),
    switch(data.view),
        case 1, data.proj=[1,0,0;0,1,0;0,0,1];
        case 2, data.proj=[1,0,0;0,0,1;0,1,0];
        case 3, data.proj=[0,0,1;1,0,0;0,1,0];
        case 0, data.proj=[1,0,0;0,1,0;0,0,1];
    end
end
% projects coordinates and background image
scale=1;
if isempty(data.x)||isempty(data.y),
    if ~data.view
        data.x=data.xy2*data.proj(1:2,1);data.y=data.xy2*data.proj(1:2,2);data.z=zeros(size(data.xy2,1),1);
    else
        data.x=data.xyz2*data.proj(:,1);data.y=data.xyz2*data.proj(:,2);data.z=data.xyz2*data.proj(:,3);
        lim=[1,1,1;data.ref.dim];refminmax=sort([lim((dec2bin(0:7)-'0'+1)+repmat([0,2,4],[8,1])),ones(8,1)]*data.ref.mat(1:3,:)'*data.proj(:,1:2));
        [data.bgx,data.bgy]=meshgrid(refminmax(1,1):scale:refminmax(end,1),refminmax(1,2):scale:refminmax(end,2));
        data.bgimage=spm_get_data(data.ref,pinv(data.ref.mat)*[data.proj(:,1:3)*[data.bgx(:),data.bgy(:),data.bgz+zeros(prod(size(data.bgx)),1)]';ones(1,prod(size(data.bgx)))]);
    end
    initxy=true;
end
%if isequal(data.source,1:numel(data.displaytheserois)), data.source=0; end
N=length(data.names);
N2=length(data.names2);
switch(data.display),
    case 'connectivity',
        %cmap=jet(64);cmap=cmap(8:56,:);
        cmap=jet(256);%cmap=cmap(32:224,:)*.8;
        %cmap=cmap.^repmat(.1+.9*abs(linspace(1,-1,size(cmap,1)))',1,size(cmap,2));
        %maxdataclusters=max(data.clusters);
        %if ~data.displaygui||(numel(data.source)==1&&data.source>0), 
        %    set(data.handles([14 15 16 18]),'visible','off');
        %else
        %    set(data.handles([14 15 16 18]),'visible','on');
        %end
%         if ~data.PERMenabled,
%             %set(data.handles(14),'string','threshold seed ROIs (F-test)','value',1);
%             set(data.handles(10),'string','Enable permutation tests');
%             %data.mvpathrmeasure=1;
%         else
%             %set(data.handles(14),'string',{'threshold seed ROIs (F-test)','threshold seed ROIs (NBS; by intensity)','threshold seed ROIs (NBS; by size)','threshold networks (NBS; by intensity)','threshold networks (NBS; by size)'});
%             set(data.handles(10),'string','Disable permutation tests');
%         end
        %if data.mvpathrmeasure==1, set(data.handles(16),'string',{'p-uncorrected','p-FDR'},'value',min(2,data.mvpathrtype)); data.mvpathrtype=min(2,data.mvpathrtype);
        %elseif data.mvpathrmeasure>3, set(data.handles(16),'string',{'p-uncorrected','p-FWE'},'value',min(2,data.mvpathrtype)); data.mvpathrtype=min(2,data.mvpathrtype);
        %else set(data.handles(16),'string',{'p-uncorrected','p-FDR','p-FWE'});
        %end
        %if data.mvpathrmeasure==1, data.mvpathrtype=1+(data.mvpathrtype>1); set(data.handles(16),'value',data.mvpathrtype); 
        %elseif data.mvpathrmeasure>3, data.mvpathrtype=1+2*(data.mvpathrtype>1); set(data.handles(16),'value',data.mvpathrtype); 
        %end
        %if data.mvpaenablethr, set(data.handles([14:16]),'enable','on');
        %else set(data.handles([14:16]),'enable','off');
        %end
        %if data.enablethr, set(data.handles(1:4),'enable','on');
        %else set(data.handles(1:4),'enable','off');
        %end
        if isequal(data.statsname,'T')
            set(data.handles(4),'value',data.side); 
            if data.enablethr, set(data.handles(4),'enable','on'); end
        else
            set(data.handles(4),'enable','off','value',3);
        end
        if data.mvpathrtype==1,
            set(data.handles(15),'enable','off'); 
        else
            set(data.handles(15),'enable','on'); 
        end
        
        %figure(hfig);
        set(hfig,'pointer','watch');%drawnow;
        %hcontrols=findobj(hfig,'enable','on');
        hcontrols=findobj(hfig,'enable','on','-not','style','edit');
        hcontrols=hcontrols(ishandle(hcontrols));
        set(hcontrols,'enable','off');
        th1=axes('units','norm','position',data.plotposition{1},'parent',data.hfig);th2=patch([0 0 1 1],[0 1 1 0],'k','edgecolor','none','facecolor',get(hfig,'color'),'facealpha',.5,'parent',th1);set(th1,'xlim',[0 1],'ylim',[0 1],'visible','off'); 
        %figure(data.hfig); set(0,'CurrentFigure',data.hfig); axes(th1); 
        drawnow; delete([th1 th2]);
        
        %set(data.handles(17),'value',data.mvpasortresultsby);
        %mvpaenabled=data.mvpaenablethr&&~(numel(data.source)==1&&data.source>0);

        if ~isfield(data,'clusters_options')||isempty(data.clusters_options), 
            olddisplaytheserois=data.displaytheserois;
            oldclusters=data.clusters;
            data.clusters_options=struct('type','hc','groups',nan,'param',.05);
            data=conn_displayroi_clusters(data); 
            if ~isequal(olddisplaytheserois,data.displaytheserois)||~isequal(oldclusters,data.clusters),
                data.PERM=[];
                data.tfceZ=[];
                data.cMVPAF=[];
            end
            if ~data.view
                data.x=data.xy2*data.proj(1:2,1);data.y=data.xy2*data.proj(1:2,2);data.z=zeros(size(data.xy2,1),1);
            else
                data.x=data.xyz2*data.proj(:,1);data.y=data.xyz2*data.proj(:,2);data.z=data.xyz2*data.proj(:,3);
                lim=[1,1,1;data.ref.dim];refminmax=sort([lim((dec2bin(0:7)-'0'+1)+repmat([0,2,4],[8,1])),ones(8,1)]*data.ref.mat(1:3,:)'*data.proj(:,1:2));
                [data.bgx,data.bgy]=meshgrid(refminmax(1,1):scale:refminmax(end,1),refminmax(1,2):scale:refminmax(end,2));
                data.bgimage=spm_get_data(data.ref,pinv(data.ref.mat)*[data.proj(:,1:3)*[data.bgx(:),data.bgy(:),data.bgz+zeros(prod(size(data.bgx)),1)]';ones(1,prod(size(data.bgx)))]);
            end
            initxy=true;
        end
        
        % computes stat threshold
        if isequal(data.statsname,'T')
            switch(data.side),
                case 1,p=data.p; Fthr=data.F;
                case 2,p=data.p2; Fthr=-data.F;
                case 3,p=2*min(data.p,data.p2); Fthr=abs(data.F);
            end
        else
            p=data.p;
            Fthr=data.F;
        end
        p(setdiff(1:N,data.displaytheserois),:)=nan;
        p(:,setdiff(1:N2,data.displaytheserois))=nan;
        if ~data.mvpathrtype_isftest(data.mvpathrtype), 
            if any(~data.source), P=reshape(conn_fdr(p(:)),size(p));
            else P=nan(size(p)); tempidx=data.source(data.source>0&data.source<=length(data.displaytheserois)); temp=p(intersect(1:N,data.displaytheserois(tempidx)),:); temp(:)=conn_fdr(temp(:)); P(intersect(1:N,data.displaytheserois(tempidx)),:)=temp; 
            end
        else P=conn_fdr(p,2);
        end
        data.P=P;
        data.Ppos=conn_fdr(p,2);
        mvpap=data.MVPAp;
        mvpap(setdiff(1:N,data.displaytheserois))=nan;
        if any(~data.source), mvpaP=conn_fdr(mvpap); 
        else mvpaP=nan(size(mvpap)); tempidx=data.source(data.source>0&data.source<=length(data.displaytheserois)); temp=mvpap(intersect(1:N,data.displaytheserois(tempidx))); temp(:)=conn_fdr(temp(:)); mvpaP(intersect(1:N,data.displaytheserois(tempidx)))=temp; end;
        data.MVPAP=mvpaP;
        
        THR=data.thr;
        switch(data.thrtype)
            case 1, THR_TYPE=1; %p-unc
            case 2, THR_TYPE=3; %p-fdr
            case 3, THR_TYPE=5; THR=0; %p-fdr TCFE
            case 4, THR_TYPE=5; THR=0; %p-fwe TFCE
            case 5, THR_TYPE=4; %T/F/X
        end
        SIDE=data.side;
        if isempty(data.PERM)||~any(data.PERM.Pthr==THR&data.PERM.Pthr_type==THR_TYPE&data.PERM.Pthr_side==SIDE)
            if conn_existfile(conn_displayroi_simfilename(data.roifile,THR_TYPE,THR,data.displaytheserois))||((data.mvpathrtype_isnonparam(data.mvpathrtype)||data.thrtype==3||data.thrtype==4)&&conn_displayroi_randomise(data,THR_TYPE,THR,init~=1)),
                try, 
                    data.PERM=conn_loadmatfile(conn_displayroi_simfilename(data.roifile,THR_TYPE,THR,data.displaytheserois)); 
                    if ~isfield(data.PERM,'VERSION'), data.PERM.VERSION=0; end
                    if data.PERM.VERSION<2, data.PERM=[]; conn_fileutils('spm_unlink',conn_displayroi_simfilename(data.roifile,THR_TYPE,THR,data.displaytheserois)); end % note: disregard older versions                    
                end
            end
        end
        if isempty(data.PERM)||~any(data.PERM.Pthr==THR&data.PERM.Pthr_type==THR_TYPE&data.PERM.Pthr_side==SIDE)
            if data.mvpathrtype_isnonparam(data.mvpathrtype)||data.thrtype==3||data.thrtype==4
                conn_msgbox({'Unable to compute non-parametric statistics. Please try again later.','Switching to parametric statistics'},'',2);
                data.thres=1; % note: assume default#1 is parametric
                set(data.handles(1),'value',data.thres);
                set(hfig,'userdata',data);
                conn_displayroi(hfig,[],'fwec.option',[],'immediatereturn');
                data=get(hfig,'userdata');
                THR=data.thr;
                switch(data.thrtype)
                    case 1, THR_TYPE=1; %p-unc
                    case 2, THR_TYPE=3; %p-fdr
                    case 3, THR_TYPE=5; THR=0; %p-fdr TCFE
                    case 4, THR_TYPE=5; THR=0; %p-fwe TFCE
                    case 5, THR_TYPE=4; %T/F/X
                end
                SIDE=data.side;
            end
        end
%         if data.mvpathrtype_isnonparam(data.mvpathrtype)&&(isempty(data.PERM)||~any(data.PERM.Pthr==data.thr&data.PERM.Pthr_type==data.thrtype&data.PERM.Pthr_side==data.side))
%             % update permutation tests
%             niterations=10000;
%             Y=permute(cat(4,data.results.y),[1,3,4,2]);
%             Y=Y(:,:,data.displaytheserois(data.displaytheserois<=N),data.displaytheserois);
%             tthr=repmat(data.thr,1,3);
%             tthrtype=repmat(data.thrtype,1,3);
%             thrside=1:3;
%             %if data.thrtype==1, tthr=tthr(tthrtype==1);tthrside=tthrside(tthrtype==1);tthrtype=tthrtype(tthrtype==1); end
%             %if ~any(tthr==data.thr&tthrtype==data.thrtype&tthrside==data.side), tthr=[data.thr, tthr]; tthrtype=[data.thrtype, tthrtype]; tthrside=[data.side, tthrside]; end
%             try
%                 data.PERM=conn_randomise(data.results(1).xX.X,Y,data.results(1).c,data.results(1).c2,tthr,tthrtype,tthrside,niterations,data.PERM,[],'matrix');
%             end
%         end
        
        if data.enablethr
            switch(data.thrtype)
                case 1, show=p<=data.thr;
                case 2, show=P<=data.thr;
                case 5, show=Fthr>=data.thr;
            end
        else
            show=~isnan(p);
            THR_TYPE=nan;
        end
        if (data.thrtype==3||data.thrtype==4) % TFCE FDR/FWE
            data.mvpathrtype=1; % note: force no cluster thresholding
            if isempty(data.tfceZ)
                Z=data.F(data.displaytheserois,data.displaytheserois);
                Z=(Z+Z')/2; Z(~triu(ones(size(Z)),1))=0; % note: tfce uses symmetric mtx
                if isequal(data.statsname,'T')
                    [tfceZpos,tfceZpospeaks,tfceZposd]=conn_tfce(abs(Z).*(Z>0),'Hmin',1); %T stat positive-sided
                    [tfceZneg,tfceZnegpeaks,tfceZnegd]=conn_tfce(abs(Z).*(Z<0),'Hmin',1); %T stat negative-sided
                    data.tfceZ=tfceZpos.*(Z>0)-tfceZneg.*(Z<0);
                    data.tfceZpeaks=(tfceZpospeaks&(Z>0))|(tfceZnegpeaks&(Z<0));
                    data.tfceZd=tfceZposd.*(Z>0)-tfceZnegd.*(Z<0);
                else
                    [data.tfceZ,data.tfceZpeaks,data.tfceZd]=conn_tfce(sqrt(abs(Z)),'Hmin',1); %F stat (keep same scale as two-tailed T)
                end
                temp=zeros(size(p)); temp(data.displaytheserois,data.displaytheserois)=data.tfceZ; data.tfceZ=temp;
                temp=zeros(size(p)); temp(data.displaytheserois,data.displaytheserois)=data.tfceZd; data.tfceZd=temp;
                temp=false(size(p)); temp(data.displaytheserois,data.displaytheserois)=data.tfceZpeaks; data.tfceZpeaks=temp;
            end
            if data.enablethr&&~isempty(data.PERM)&&any(data.PERM.Pthr==THR&data.PERM.Pthr_type==THR_TYPE&data.PERM.Pthr_side==data.side)
                if isequal(data.statsname,'T')|data.side==3, ttZ=abs(data.tfceZ); tfce_peak_mask=data.tfceZpeaks; ttZb=abs(data.tfceZd); %tfce_cluster_mask=tfceZclusters;
                elseif data.side==1,         ttZ=data.tfceZ; tfce_peak_mask=data.tfceZpeaks&(data.tfceZ>0); ttZb=data.tfceZd; %tfce_cluster_mask=tfceZclusters.*(tZ>0);
                elseif data.side==2,         ttZ=-data.tfceZ; tfce_peak_mask=data.tfceZpeaks&(data.tfceZ<0); ttZb=-data.tfceZd; %tfce_cluster_mask=tfceZclusters.*(tZ<0);
                else error('incorrect option ''side''');
                end
                iPERM=find(data.PERM.Pthr==THR&data.PERM.Pthr_type==THR_TYPE&data.PERM.Pthr_side==data.side);
                tfce_max_dist=[];
                tfce_peak_p_unc=0;
                NiPERM=0;
                for niPERM=1:numel(iPERM),
                    iperm=iPERM(niPERM);
                    kiPERM=size(data.PERM.Dist_Voxel_statmax{iperm},1);
                    tfce_max_dist=[tfce_max_dist;data.PERM.Dist_Voxel_statmax{iperm}];                    
                    if nnz(data.PERM.Hist_Voxel_stat{iperm})<2, tfce_peak_p_unc=tfce_peak_p_unc+kiPERM*double(1+round(data.PERM.maxT*100*ttZb(tfce_peak_mask)/sqrt(prod(data.PERM.model.dims)/prod(data.PERM.model.dims(1:2))))<=find(data.PERM.Hist_Voxel_stat{iperm}));
                    else tfce_peak_p_unc=tfce_peak_p_unc+kiPERM*max(0,min(1,interp1(find(data.PERM.Hist_Voxel_stat{iperm}),flipud(cumsum(flipud(nonzeros(data.PERM.Hist_Voxel_stat{iperm})))),1+round(data.PERM.maxT*100*ttZb(tfce_peak_mask)/sqrt(prod(data.PERM.model.dims)/prod(data.PERM.model.dims(1:2)))),'linear','extrap')));
                    end
                    NiPERM=NiPERM+kiPERM;
                % note: Hist_Voxel_stat*sqrt(prod(data.PERM.model.dims)/prod(data.PERM.model.dims(1:2)))
                end
                tfce_peak_p_unc=tfce_peak_p_unc/NiPERM; % peak p-unc
                tfce_peak_p_FDR=conn_fdr(tfce_peak_p_unc); % peak p-FDR
                if data.thrtype==3, % TFCE FDR  
                    temp0=ttZb(tfce_peak_mask);
                    temp1=max(temp0(tfce_peak_p_FDR>data.thr)); if isempty(temp1), temp1=0; end
                    temp2=min(temp0(tfce_peak_p_FDR<=data.thr)); if isempty(temp2), temp2=inf; end
                    thr=(temp1+temp2)/2;
                elseif data.thr>1, thr=0; 
                else thr=interp1([0 ((1:numel(tfce_max_dist))-0.5)/numel(tfce_max_dist) 1],[min(tfce_max_dist),reshape(sort(tfce_max_dist),1,[]),max(tfce_max_dist)],1-max(0,min(1,data.thr))); % TFCE FWE  
                end
                data.ttZ=ttZ;
                data.ttZb=ttZb;
                data.ttZbthr=thr;
                show=ttZb>thr;
                V=double(show);
                V(show)=abs(ttZb(show));
            else
                show=false(size(p));
                V=double(show);
            end
            tdatatype='matrix';
        else
            show(setdiff(1:N,data.displaytheserois),:)=0;
            show(:,setdiff(1:N2,data.displaytheserois))=0;
            V=double(show);
            if ~isempty(data.PERM)&&data.PERM.VERSION<1, DOSCALEF=false; end
            if ~DOSCALEF||~isequal(data.statsname,'T'), V(show)=abs(data.F(show));
            else V(show)=abs(data.F(show)).^2;
            end
            % cluster-size & cluster-mass for each cluster
            %[nill,tidx]=sort(mod(pi+angle(data.xy2(data.displaytheserois,:)*[1;1i]),2*pi));data.displaytheserois=data.displaytheserois(tidx);
            if data.mvpathrtype_isnetwork(data.mvpathrtype), tdatatype='network';
            else tdatatype='matrix';
            end
        end
        
        if isempty(data.cMVPAF)&&~isempty(data.clusters)
            numc=max(data.clusters(data.displaytheserois));
            ndims=min(4,ceil(sqrt(size(data.results(1).data,1))/4)); % rule of thumb for number of dimensions
            data.cMVPAF=nan(numc,numc);
            data.cMVPAp=nan(numc,numc);
            data.cMVPAstatsname=repmat({''},[numc,numc]);
            data.cMVPAdofstr=repmat({''},[numc,numc]);
            data.cMVPAlabel=zeros(numc,numc);
            for num1=1:numc,
                domvpa1=data.displaytheserois(data.clusters(data.displaytheserois)==num1);
                for num2=1:num1
                    domvpa2=data.displaytheserois(data.clusters(data.displaytheserois)==num2);
                    y=[];
                    for nresults1=reshape(domvpa1,1,[])
                        ty=data.results(1).data(:,domvpa2,:,nresults1);
                        ty(:,any(any(isnan(ty),1),3),:)=[]; % subjects x rois x conditions
                        sy=[size(ty),1,1];
                        y=cat(2,y,reshape(permute(ty,[1,3,2]),[sy(1)*sy(3),sy(2)])); %(subjectsxconditions) x rois
                    end
                    if ~isempty(y)
                        if size(y,2)<=size(y,1), [Q,D,R]=svd(y,0);
                        else [R,D,Q]=svd(y',0);
                        end
                        ndims2=max(1,min(size(R,2),ndims));
                        d=D(1:size(D,1)+1:size(D,1)*min(size(D))).^2;
                        ndims2=min([ndims2,find(cumsum(d)/sum(d)>.95,1)]); % 95 percent variance
                        y=y*R(:,1:ndims2);
                        y=permute(reshape(y,[sy(1),sy(3),ndims2]),[1,3,2]); % subjects x components x conditions
                        data.cMVPApcacov=d(1:ndims2)/sum(d);
                        [data.cMVPAh,data.cMVPAF(num1,num2),data.cMVPAp(num1,num2),data.cMVPAdof,data.cMVPAstatsname{num1,num2}]=conn_glm(data.results(1).xX.X,y(:,:),data.results(1).c,kron(data.results(1).c2,eye(size(y,2))));
                        if ~nnz(diff(y(:,:),1)), data.cMVPAp(num1,num2)=NaN; data.cMVPAF(num1,num2)=0; end % note: skips constant data
                        if isequal(data.cMVPAstatsname{num1,num2},'T'),
                            data.cMVPAstatsname{num1,num2}='F';
                            data.cMVPAF(num1,num2)=data.cMVPAF(num1,num2).^2;
                            data.cMVPAp(num1,num2)=2*min(data.cMVPAp(num1,num2),1-data.cMVPAp(num1,num2));
                            data.cMVPAdofstr{num1,num2}=['(1,',num2str(data.cMVPAdof(1)),')'];
                        else
                            data.cMVPAdofstr{num1,num2}=['(',num2str(data.cMVPAdof(1)),',',num2str(data.cMVPAdof(2)),')'];
                        end
                        data.cMVPAlabel(num1,num2)=1;
                    end
                end                
            end
            data.cMVPAlabel(data.cMVPAlabel>0)=1:nnz(data.cMVPAlabel>0);
            data.cMVPAP=data.cMVPAp; data.cMVPAP(:)=conn_fdr(data.cMVPAp(:));
        end
        if data.mvpathrtype_isroi(data.mvpathrtype)
            CLUSTER_labels=zeros(size(show));
            CLUSTER_labels(data.displaytheserois,data.displaytheserois)=repmat((1:numel(data.displaytheserois))',1,numel(data.displaytheserois));
            nclL=accumarray(CLUSTER_labels(CLUSTER_labels>0),1);
        elseif data.mvpathrtype_isftest(data.mvpathrtype)&&data.mvpathrtype_iscluster(data.mvpathrtype)
            CLUSTER_labels=zeros(size(show));
            temp=max(data.cMVPAlabel,data.cMVPAlabel');
            CLUSTER_labels(data.displaytheserois,data.displaytheserois)=temp(data.clusters(data.displaytheserois),data.clusters(data.displaytheserois));
            %numc=size(data.cMVPAF,1);
            %CLUSTER_labels(data.displaytheserois,data.displaytheserois)=repmat(reshape(data.clusters(data.displaytheserois),[],1),1,numel(data.displaytheserois))+numc*(repmat(reshape(data.clusters(data.displaytheserois),[],1),1,numel(data.displaytheserois))'-1);
            nclL=accumarray(CLUSTER_labels(CLUSTER_labels>0),1);
        elseif isequal(data.statsname,'T')&&strcmp(tdatatype,'matrix') % note: remove '&&strcmp(tdatatype,'matrix')' to have network-level stats separate for positive/negative effects
            CLUSTER_labels1=zeros(size(show));
            [nclL1,CLUSTER_labels1(data.displaytheserois,data.displaytheserois)]=conn_clusters(show(data.displaytheserois,data.displaytheserois)&(data.F(data.displaytheserois,data.displaytheserois)>0),tdatatype);
            CLUSTER_labels2=zeros(size(show));
            [nclL2,CLUSTER_labels2(data.displaytheserois,data.displaytheserois)]=conn_clusters(show(data.displaytheserois,data.displaytheserois)&(data.F(data.displaytheserois,data.displaytheserois)<0),tdatatype);
            CLUSTER_labels=CLUSTER_labels1;
            CLUSTER_labels(CLUSTER_labels1==0&CLUSTER_labels2>0)=CLUSTER_labels2(CLUSTER_labels1==0&CLUSTER_labels2>0)+numel(nclL1);
            nclL=[reshape(nclL1,[],1);reshape(nclL2,[],1)];
        else
            CLUSTER_labels=zeros(size(show));
            [nclL,CLUSTER_labels(data.displaytheserois,data.displaytheserois)]=conn_clusters(show(data.displaytheserois,data.displaytheserois),tdatatype);
        end
        %[nclL,CLUSTER_labels]=conn_clusters(show);
        data.CLUSTER_labels=CLUSTER_labels;
        CLroi_labels=max(CLUSTER_labels,[],2);
        mask=CLUSTER_labels>0;
        mclL=accumarray(CLUSTER_labels(mask),V(mask),[max([0,max(CLUSTER_labels(mask))]),1]);
        if ~DOSCALEF, sclL=mclL.^3./max(eps,nclL).^2.5/3; % cluster TFCE-score
        else  sclL=mclL.^(3/2)./max(eps,nclL)/3; % cluster TFCE-score
        end
        % seed-size & seed-mass for each seed
        nsdL=sum(show,2);
        msdL=sum(V,2);
        if ~DOSCALEF, ssdL=msdL.^3./max(eps,nsdL).^2.5/3; % seed TFCE-score
        else ssdL=msdL.^(3/2)./max(eps,nsdL)/3; % seed TFCE-score
        end
        if (data.thrtype==3||data.thrtype==4), maxZ=accumarray(CLUSTER_labels(mask),V(mask),[max([0,max(CLUSTER_labels(mask))]),1],@max); end
        if ~isempty(data.PERM)&&any(data.PERM.Pthr==THR&data.PERM.Pthr_type==THR_TYPE&data.PERM.Pthr_side==data.side)
            data.iPERM=find(data.PERM.Pthr==THR&data.PERM.Pthr_type==THR_TYPE&data.PERM.Pthr_side==data.side);
            PERMp_cluster_size_unc=0;PERMp_seed_size_unc=0;
            PERMp_cluster_size_FWE=0;PERMp_seed_size_FWE=0;
            PERMp_cluster_mass_unc=0;PERMp_seed_mass_unc=0;
            PERMp_cluster_mass_FWE=0;PERMp_seed_mass_FWE=0;
            PERMp_cluster_score_unc=0;PERMp_seed_score_unc=0;
            PERMp_cluster_score_FWE=0;PERMp_seed_score_FWE=0;
            if (data.thrtype==3||data.thrtype==4) % TFCE FDR/FWE
                %peaks(n1)=nnz(tfce_peak_mask(clusters{n1}));
                PERMp_tfce_FWE=zeros(size(maxZ));PERMp_tfce_FDR=PERMp_tfce_FWE;PERMp_tfce_unc=PERMp_tfce_FWE;
                for ncl=1:numel(maxZ)
                    PERMp_tfce_FWE(ncl)=mean(tfce_max_dist>maxZ(ncl));
                    tidx=1:numel(tfce_peak_p_FDR);
                    tidx(CLUSTER_labels(tfce_peak_mask)~=ncl)=[];
                    if ~isempty(tidx),
                        PERMp_tfce_FDR(ncl)=min(tfce_peak_p_FDR(tidx));
                        PERMp_tfce_unc(ncl)=min(tfce_peak_p_unc(tidx));
                    end
                end
            else
                NiPERM=0;
                for niPERM=1:numel(data.iPERM),
                    iperm=data.iPERM(niPERM);
                    kiPERM=size(data.PERM.Dist_Cluster_sizemax{iperm},1);
                    if nnz(data.PERM.Hist_Seed_size{iperm})<2, PERMp_seed_size_unc=PERMp_seed_size_unc+kiPERM*double(1+nsdL<=find(data.PERM.Hist_Seed_size{iperm}));
                    else PERMp_seed_size_unc=PERMp_seed_size_unc+kiPERM*max(0,min(1,interp1(find(data.PERM.Hist_Seed_size{iperm}),flipud(cumsum(flipud(nonzeros(data.PERM.Hist_Seed_size{iperm})))),1+nsdL,'linear','extrap')));
                    end
                    PERMp_seed_size_FWE=PERMp_seed_size_FWE+kiPERM*mean(conn_bsxfun(@ge,data.PERM.Dist_Seed_sizemax{iperm}',nsdL),2);
                    if nnz(data.PERM.Hist_Seed_mass{iperm})<2, PERMp_seed_mass_unc=PERMp_seed_mass_unc+kiPERM*double(1+round(data.PERM.maxT*msdL)<=find(data.PERM.Hist_Seed_mass{iperm}));
                    else PERMp_seed_mass_unc=PERMp_seed_mass_unc+kiPERM*max(0,min(1,interp1(find(data.PERM.Hist_Seed_mass{iperm}),flipud(cumsum(flipud(nonzeros(data.PERM.Hist_Seed_mass{iperm})))),1+round(data.PERM.maxT*msdL),'linear','extrap')));
                    end
                    PERMp_seed_mass_FWE=PERMp_seed_mass_FWE+kiPERM*mean(conn_bsxfun(@ge,data.PERM.Dist_Seed_massmax{iperm}',msdL),2);
                    if nnz(data.PERM.Hist_Seed_score{iperm})<2, PERMp_seed_score_unc=PERMp_seed_score_unc+kiPERM*double(1+round(data.PERM.maxT*ssdL)<=find(data.PERM.Hist_Seed_score{iperm}));
                    else PERMp_seed_score_unc=PERMp_seed_score_unc+kiPERM*max(0,min(1,interp1(find(data.PERM.Hist_Seed_score{iperm}),flipud(cumsum(flipud(nonzeros(data.PERM.Hist_Seed_score{iperm})))),1+round(data.PERM.maxT*ssdL),'linear','extrap')));
                    end
                    PERMp_seed_score_FWE=PERMp_seed_score_FWE+kiPERM*mean(conn_bsxfun(@ge,data.PERM.Dist_Seed_scoremax{iperm}',ssdL),2);
                    if data.mvpathrtype_isnetwork(data.mvpathrtype),
                        if nnz(data.PERM.Hist_Network_size{iperm})<2, PERMp_cluster_size_unc=PERMp_cluster_size_unc+kiPERM*double(1+nclL<=find(data.PERM.Hist_Network_size{iperm}));
                        else PERMp_cluster_size_unc=PERMp_cluster_size_unc+kiPERM*max(0,min(1,interp1(find(data.PERM.Hist_Network_size{iperm}),flipud(cumsum(flipud(nonzeros(data.PERM.Hist_Network_size{iperm})))),1+nclL,'linear','extrap')));
                        end
                        PERMp_cluster_size_FWE=PERMp_cluster_size_FWE+kiPERM*mean(conn_bsxfun(@ge,data.PERM.Dist_Network_sizemax{iperm}',nclL),2);
                        if nnz(data.PERM.Hist_Network_mass{iperm})<2, PERMp_cluster_mass_unc=PERMp_cluster_mass_unc+kiPERM*double(1+round(data.PERM.maxT*mclL)<=find(data.PERM.Hist_Network_mass{iperm}));
                        else PERMp_cluster_mass_unc=PERMp_cluster_mass_unc+kiPERM*max(0,min(1,interp1(find(data.PERM.Hist_Network_mass{iperm}),flipud(cumsum(flipud(nonzeros(data.PERM.Hist_Network_mass{iperm})))),1+round(data.PERM.maxT*mclL),'linear','extrap')));
                        end
                        PERMp_cluster_mass_FWE=PERMp_cluster_mass_FWE+kiPERM*mean(conn_bsxfun(@ge,data.PERM.Dist_Network_massmax{iperm}',mclL),2);
                        if nnz(data.PERM.Hist_Network_score{iperm})<2, PERMp_cluster_score_unc=PERMp_cluster_score_unc+kiPERM*double(1+round(data.PERM.maxT*sclL)<=find(data.PERM.Hist_Network_score{iperm}));
                        else PERMp_cluster_score_unc=PERMp_cluster_score_unc+kiPERM*max(0,min(1,interp1(find(data.PERM.Hist_Network_score{iperm}),flipud(cumsum(flipud(nonzeros(data.PERM.Hist_Network_score{iperm})))),1+round(data.PERM.maxT*sclL),'linear','extrap')));
                        end
                        PERMp_cluster_score_FWE=PERMp_cluster_score_FWE+kiPERM*mean(conn_bsxfun(@ge,data.PERM.Dist_Network_scoremax{iperm}',sclL),2);
                    else
                        if nnz(data.PERM.Hist_Cluster_size{iperm})<2, PERMp_cluster_size_unc=PERMp_cluster_size_unc+kiPERM*double(1+nclL<=find(data.PERM.Hist_Cluster_size{iperm}));
                        else PERMp_cluster_size_unc=PERMp_cluster_size_unc+kiPERM*max(0,min(1,interp1(find(data.PERM.Hist_Cluster_size{iperm}),flipud(cumsum(flipud(nonzeros(data.PERM.Hist_Cluster_size{iperm})))),1+nclL,'linear','extrap')));
                        end
                        PERMp_cluster_size_FWE=PERMp_cluster_size_FWE+kiPERM*mean(conn_bsxfun(@ge,data.PERM.Dist_Cluster_sizemax{iperm}',nclL),2);
                        if nnz(data.PERM.Hist_Cluster_mass{iperm})<2, PERMp_cluster_mass_unc=PERMp_cluster_mass_unc+kiPERM*double(1+round(data.PERM.maxT*mclL)<=find(data.PERM.Hist_Cluster_mass{iperm}));
                        else PERMp_cluster_mass_unc=PERMp_cluster_mass_unc+kiPERM*max(0,min(1,interp1(find(data.PERM.Hist_Cluster_mass{iperm}),flipud(cumsum(flipud(nonzeros(data.PERM.Hist_Cluster_mass{iperm})))),1+round(data.PERM.maxT*mclL),'linear','extrap')));
                        end
                        PERMp_cluster_mass_FWE=PERMp_cluster_mass_FWE+kiPERM*mean(conn_bsxfun(@ge,data.PERM.Dist_Cluster_massmax{iperm}',mclL),2);
                        if nnz(data.PERM.Hist_Cluster_score{iperm})<2, PERMp_cluster_score_unc=PERMp_cluster_score_unc+kiPERM*double(1+round(data.PERM.maxT*sclL)<=find(data.PERM.Hist_Cluster_score{iperm}));
                        else PERMp_cluster_score_unc=PERMp_cluster_score_unc+kiPERM*max(0,min(1,interp1(find(data.PERM.Hist_Cluster_score{iperm}),flipud(cumsum(flipud(nonzeros(data.PERM.Hist_Cluster_score{iperm})))),1+round(data.PERM.maxT*sclL),'linear','extrap')));
                        end
                        PERMp_cluster_score_FWE=PERMp_cluster_score_FWE+kiPERM*mean(conn_bsxfun(@ge,data.PERM.Dist_Cluster_scoremax{iperm}',sclL),2);
                    end
                    NiPERM=NiPERM+kiPERM;
                end
                PERMp_cluster_size_unc=PERMp_cluster_size_unc/NiPERM;PERMp_seed_size_unc=PERMp_seed_size_unc/NiPERM;
                PERMp_cluster_size_FWE=PERMp_cluster_size_FWE/NiPERM;PERMp_seed_size_FWE=PERMp_seed_size_FWE/NiPERM;
                PERMp_cluster_mass_unc=PERMp_cluster_mass_unc/NiPERM;PERMp_seed_mass_unc=PERMp_seed_mass_unc/NiPERM;
                PERMp_cluster_mass_FWE=PERMp_cluster_mass_FWE/NiPERM;PERMp_seed_mass_FWE=PERMp_seed_mass_FWE/NiPERM;
                PERMp_cluster_score_unc=PERMp_cluster_score_unc/NiPERM;PERMp_seed_score_unc=PERMp_seed_score_unc/NiPERM;
                PERMp_cluster_score_FWE=PERMp_cluster_score_FWE/NiPERM;PERMp_seed_score_FWE=PERMp_seed_score_FWE/NiPERM;
                PERMp_seed_size_unc(setdiff(1:N,data.displaytheserois))=nan;
                PERMp_seed_size_FWE(setdiff(1:N,data.displaytheserois))=nan;
                temp=PERMp_seed_size_unc;
                if any(~data.source), tempP=conn_fdr(temp);
                else tempP=nan(size(temp)); tempidx=data.source(data.source>0&data.source<=length(data.displaytheserois)); temp=temp(intersect(1:N,data.displaytheserois(tempidx))); temp(:)=conn_fdr(temp(:)); tempP(intersect(1:N,data.displaytheserois(tempidx)))=temp; end;
                PERMp_seed_size_FDR=tempP;
                PERMp_seed_mass_unc(setdiff(1:N,data.displaytheserois))=nan;
                PERMp_seed_mass_FWE(setdiff(1:N,data.displaytheserois))=nan;
                temp=PERMp_seed_mass_unc;
                if any(~data.source), tempP=conn_fdr(temp);
                else tempP=nan(size(temp)); tempidx=data.source(data.source>0&data.source<=length(data.displaytheserois)); temp=temp(intersect(1:N,data.displaytheserois(tempidx))); temp(:)=conn_fdr(temp(:)); tempP(intersect(1:N,data.displaytheserois(tempidx)))=temp; end;
                PERMp_seed_mass_FDR=tempP;
                PERMp_seed_score_unc(setdiff(1:N,data.displaytheserois))=nan;
                PERMp_seed_score_FWE(setdiff(1:N,data.displaytheserois))=nan;
                temp=PERMp_seed_score_unc;
                if any(~data.source), tempP=conn_fdr(temp);
                else tempP=nan(size(temp)); tempidx=data.source(data.source>0&data.source<=length(data.displaytheserois)); temp=temp(intersect(1:N,data.displaytheserois(tempidx))); temp(:)=conn_fdr(temp(:)); tempP(intersect(1:N,data.displaytheserois(tempidx)))=temp; end;
                PERMp_seed_score_FDR=tempP;
                
                PERMp_cluster_size_FDR=conn_fdr(PERMp_cluster_size_unc);
                PERMp_cluster_mass_FDR=conn_fdr(PERMp_cluster_mass_unc);
                PERMp_cluster_score_FDR=conn_fdr(PERMp_cluster_score_unc);
                %conn_disp('mass');conn_disp(char(data.names(find(PERMp_seed_mass_FWE<=.05))));
                %conn_disp('size');conn_disp(char(data.names(find(PERMp_seed_size_FWE<=.05))));
                %conn_disp('score');conn_disp(char(data.names(find(PERMp_seed_score_FWE<=.05))));
            end
        else
%             if data.mvpathrtype_isnonparam(data.mvpathrtype)||data.thrtype==3||data.thrtype==4
%                 conn_msgbox({'Unable to compute non-parametric statistics. Please try again later.','Switching to parametric statistics'},'',2);
%                 data.thres=1; % note: assume default#1 is parametric
%                 set(data.handles(1),'value',data.thres);
%                 set(hfig,'userdata',data);
%                 conn_displayroi('fwec.option',hfig,'immediatereturn');
%                 data=get(hfig,'userdata');
% %                 [data.thr,data.thrtype,data.mvpathr,data.mvpathrtype]=deal(data.thres_defaults{data.thres}{:});
% %                 data.side=3;
% %                 set(data.handles(2),'string',num2str(data.thr));
% %                 set(data.handles(3),'value',data.thrtype);
% %                 set(data.handles(15),'string',num2str(data.mvpathr));
% %                 set(data.handles(16),'value',find(data.mvpathrtype_shown==data.mvpathrtype));
% %                 if data.thrtype==3, set(data.handles(22),'string',['connection threshold: ',data.statsname,' > ']);
% %                 elseif data.mvpathr==1, set(data.handles(22),'string','connection threshold: ');
% %                 else set(data.handles(22),'string','connection threshold: p < ');
% %                 end
% %                 set(data.handles([2,3,4,15,16,22,23]),'visible','off');
%             end
            data.iPERM=[];
        end
        if data.thrtype==3 % tfce
            seedmask=zeros(size(data.MVPAp));
            netmask=PERMp_tfce_FDR;
        elseif data.thrtype==4 % tfce
            seedmask=zeros(size(data.MVPAp));
            netmask=PERMp_tfce_FWE;
        elseif 1,%mvpaenabled,
            switch(data.mvpathrtype)
                case 1, % none
                    seedmask=zeros(size(data.MVPAp));
                    netmask=zeros(size(nclL));
                case {2,11}, % cluster score uncorrected
                    seedmask=zeros(size(data.MVPAp));
                    netmask=PERMp_cluster_score_unc;
                case {3,12}, % cluster score FDR-corrected
                    seedmask=zeros(size(data.MVPAp));
                    netmask=PERMp_cluster_score_FDR;
                case {4,13}, % cluster score FWE-corrected
                    seedmask=zeros(size(data.MVPAp));
                    netmask=PERMp_cluster_score_FWE;
                case {5,14}, % cluster intensity uncorrected
                    seedmask=zeros(size(data.MVPAp));
                    netmask=PERMp_cluster_mass_unc;
                case {6,15}, % cluster intensity FDR-corrected
                    seedmask=zeros(size(data.MVPAp));
                    netmask=PERMp_cluster_mass_FDR;
                case {7,16}, % cluster intensity FWE-corrected
                    seedmask=zeros(size(data.MVPAp));
                    netmask=PERMp_cluster_mass_FWE;
                case {8,17}, % cluster size uncorrected
                    seedmask=zeros(size(data.MVPAp));
                    netmask=PERMp_cluster_size_unc;
                case {9,18}, % cluster size FDR-corrected
                    seedmask=zeros(size(data.MVPAp));
                    netmask=PERMp_cluster_size_FDR;
                case {10,19}, % cluster size FWE-corrected
                    seedmask=zeros(size(data.MVPAp));
                    netmask=PERMp_cluster_size_FWE;
                case 20, % seed score uncorrected
                    seedmask=PERMp_seed_score_unc;
                    netmask=zeros(size(nclL));
                case 21, % seed score FDR-corrected
                    seedmask=PERMp_seed_score_FDR;
                    netmask=zeros(size(nclL));
                case 22, % seed score FWE-corrected
                    seedmask=PERMp_seed_score_FWE;
                    netmask=zeros(size(nclL));
                case 23, % seed intensity uncorrected
                    seedmask=PERMp_seed_mass_unc;
                    netmask=zeros(size(nclL));
                case 24, % seed intensity FDR-corrected
                    seedmask=PERMp_seed_mass_FDR;
                    netmask=zeros(size(nclL));
                case 25, % seed intensity FWE-corrected
                    seedmask=PERMp_seed_mass_FWE;
                    netmask=zeros(size(nclL));
                case 26, % seed size uncorrected
                    seedmask=PERMp_seed_size_unc;
                    netmask=zeros(size(nclL));
                case 27, % seed size FDR-corrected
                    seedmask=PERMp_seed_size_FDR;
                    netmask=zeros(size(nclL));
                case 28, % seed size FWE-corrected
                    seedmask=PERMp_seed_size_FWE;
                    netmask=zeros(size(nclL));
                case 29, % seed F-test uncorrected
                    seedmask=data.MVPAp;
                    netmask=zeros(size(nclL));
                case 30, % seed F-test FDR-corrected
                    seedmask=data.MVPAP;
                    netmask=zeros(size(nclL));
                case 31, % cluster F-test uncorrected
                    seedmask=zeros(size(data.MVPAp));
                    netmask=data.cMVPAp(data.cMVPAlabel>0);
                case 32, % cluster F-test FDR-corrected
                    seedmask=zeros(size(data.MVPAp));
                    netmask=data.cMVPAP(data.cMVPAlabel>0);
            end
        else
            seedmask=zeros(size(data.MVPAp));
            netmask=zeros(size(nclL));
        end
        %if data.displayreduced, p=p(:,1:N); end
        %p=p(intersect(1:N,data.displaytheserois),data.displaytheserois);
        seedmask(setdiff(1:N,data.displaytheserois))=nan;
%         if isfield(data,'displayeffectsize')&&data.displayeffectsize>0
%             z=abs(data.h);
%         elseif isfield(data,'displayeffectsize')&&data.displayeffectsize<0
%             z=ones(size(data.h));
%         else
        if data.thrtype==3||data.thrtype==4, z=abs(data.tfceZ); data.CM_z0=data.tfceZ; % tfce
        else z=abs(data.F); data.CM_z0=data.F; 
        end
%         end
        z(isnan(z))=0;
        seedz=data.MVPAF;
%         maxz=max(abs(z(:)));
%         if maxz>0, z(z>0)=z(z>0)/maxz; end
        if data.enablethr
            switch(data.thrtype),
                case 1, z(p>data.thr)=0;
                case 2, z(P>data.thr)=0;
                case {3,4}, %z(data.ttZb<data.ttZbthr)=0;
                case 5, z(Fthr<data.thr)=0;
            end
        end
        %if ~isfield(data,'mvpathrtype'), data.mvpathrtype=2; end
        %if ~isfield(data,'mvpathr'), data.mvpathr=.05; end
        netz=ones(size(netmask));
        %data.CLUSTER_selected=[];
        if data.thrtype==3||data.thrtype==4 % tfce
            netz(netmask>data.thr|isnan(netmask))=0;
            z(~ismember(CLUSTER_labels,find(netz)))=0;
            seedz(~isnan(seedmask))=1;
            %data.CLUSTER_selected=find(netz);
        elseif data.mvpathrtype>1, 
            if data.mvpathrtype_iscluster(data.mvpathrtype)|data.mvpathrtype_isnetwork(data.mvpathrtype), % cluster
                netz(netmask>data.mvpathr|isnan(netmask))=0;
                z(~ismember(CLUSTER_labels,find(netz)))=0;
                seedz(~isnan(seedmask))=1;
                %seedz(~ismember(CLroi_labels,find(netz)))=0;
                %data.CLUSTER_selected=find(netz);
            elseif data.mvpathrtype_isroi(data.mvpathrtype) % seed
                seedz(seedmask>data.mvpathr|isnan(seedmask))=0;
                z(seedz==0|isnan(seedz),:)=0;
                netz=seedz(data.displaytheserois);
                %data.CLUSTER_selected=find(seedz(data.displaytheserois));
            end
        else
            seedz(~isnan(seedmask))=1; % connection 
        end
        %if max(z(:))>0, z=ceil(z/max(z(:))*3); end
        z(isnan(p))=nan;
        %seedz(isnan(mvpap))=nan;
        if ~any(data.source==0), 
            z(setdiff(1:size(z,1),data.displaytheserois(data.source)),:)=nan; 
            %seedz(setdiff(1:numel(seedz),data.displaytheserois(data.source)))=nan; 
        end
        if ~data.mvpathrtype_isroi(data.mvpathrtype),seedz=max(max(z,[],2),max(z(:,1:size(z,1)),[],1)');end
        data.CM_z=z;
        data.CM_seedz=seedz;
        maxz=max(abs(z(:)));
        if maxz>0, z(z>0)=z(z>0)/maxz; end
        if isfield(data,'maxz')&&~isempty(data.maxz), emaxz=data.maxz; 
        else emaxz=maxz;
        end
        
        if data.plotconnoptions.LCOLOR==1, 
            if isequal(data.statsname,'T'), set(data.legend(1:5),'visible','on'); 
            else set(data.legend([1 3 4]),'visible','on');set(data.legend([2 5]),'visible','off');
            end
            set(data.legend(6:end),'visible','off');
        elseif data.plotconnoptions.LCOLOR==2, set(data.legend,'visible','off');
        elseif data.plotconnoptions.LCOLOR==3, 
            set(data.legend(1:5),'visible','off'); set(data.legend([3 6:end]),'visible','on');
            set(data.legend(7),'string',num2str(emaxz,'%.2f'));
            if isequal(data.statsname,'T'), 
                set(data.legend(6),'string',num2str(-emaxz,'%.2f'));
                wtemp=linspace(-1,1,size(cmap,1));
            else 
                set(data.legend(6),'string','0');
                wtemp=linspace(0,1,size(cmap,1));
            end
            for ntemp=1:size(cmap,1), set(data.legend(7+ntemp),'facecolor',cmap(max(1,min(size(cmap,1), round(size(cmap,1)/2+size(cmap,1)/2*sign(wtemp(ntemp))*abs(wtemp(ntemp))^data.plotconnoptions.LCOLORSCALE) )) ,:)); end
        end
        %set(data.handles(21),'backgroundcolor',data.plotconnoptions.BCOLOR,'foregroundcolor',1-data.plotconnoptions.BCOLOR);
        ntemp=nnz(data.displaytheserois<=numel(data.names));
        set(data.handles(12),'string',sprintf('Analysis of %d connections among %d ROIs',ntemp*(ntemp-1)/2*(1+~data.issymmetric),ntemp));
        % text lists
        if ~isfield(CONN_gui,'parse_html'), CONN_gui.parse_html={'<HTML><FONT color=rgb(100,100,100)>','</FONT></HTML>'}; end
        txt1={};for n1=1:length(data.displaytheserois),
            txt1{end+1}=(sprintf('%-s (%d)',data.names2{data.displaytheserois(n1)},n1));
            if data.displaytheserois(n1)>numel(data.names), txt1{end}=[CONN_gui.parse_html{1},txt1{end},CONN_gui.parse_html{2}]; end
        end;
        %txt1{end+1}=' ';
        txt1=char(txt1);%strvcat(txt1{:});
        %tp=[];sort2=[];txt2={};for n1=1:N,for n2=[1:n1-1,n1+1:N2+data.displayreduced*(N-N2)],
        tp3=[];sort3=[];txt3={};index3=[];
        parse_html1={'',''};%%regexprep(CONN_gui.parse_html,{'<HTML>','</HTML>','<FONT color=rgb\(\d+,\d+,\d+\)>','</FONT>'},{'<HTML><pre>','</pre></HTML>','<b>','</b>'});
        parse_html2={'',''};%regexprep(CONN_gui.parse_html,{'<HTML>','</HTML>','<FONT color=rgb\(\d+,\d+,\d+\)>'},{'<HTML><pre>','</pre></HTML>','<FONT color=rgb(128,0,0)>'});
        parse_html3={'',''};%regexprep(CONN_gui.parse_html,{'<HTML>','</HTML>','<FONT color=rgb\(\d+,\d+,\d+\)>'},{'<HTML><pre>','</pre></HTML>','<FONT color=rgb(0,0,128)>'});
        tp4=[];txt4={};index4=[];
        if data.thrtype==3||data.thrtype==4, sortmeasure=-maxZ;
        elseif data.mvpathrtype_isftest(data.mvpathrtype)&&data.mvpathrtype_iscluster(data.mvpathrtype), sortmeasure=data.cMVPAp(data.cMVPAlabel>0)-1e-10*abs(data.cMVPAF(data.cMVPAlabel>0));
        elseif ~data.mvpathrtype_isnonparam(data.mvpathrtype), sortmeasure=-mclL;
        elseif data.mvpathrtype_isce(data.mvpathrtype), sortmeasure=PERMp_cluster_score_unc-1e-10*sclL;
        elseif data.mvpathrtype_ismass(data.mvpathrtype), sortmeasure=PERMp_cluster_mass_unc-1e-10*sclL;
        else sortmeasure=PERMp_cluster_size_unc-1e-10*sclL;
        end
        data.displayallmeasures=data.displayconnectionstats;    % note: set to 1 to display all size/mass/score measures simultaneously     
        [nill,tidx]=sort(sortmeasure);
        data.CLUSTER_selected=[];
        data.CLUSTER_selected_names={};
        for n1=1:numel(tidx),
            nb1=tidx(n1);
            if netz(nb1)>0 % cluster info
                if data.thrtype==3||data.thrtype==4, tstr='Cluster';
                elseif data.mvpathrtype_iscluster(data.mvpathrtype),tstr='Cluster';
                elseif data.mvpathrtype_isnetwork(data.mvpathrtype),tstr='Network';
                else tstr='ROI';
                end
                tstr=sprintf('%s %d/%d',tstr,n1,numel(netz));
                if data.mvpathrtype_iscluster(data.mvpathrtype)||data.mvpathrtype_isnetwork(data.mvpathrtype)||(data.thrtype==3||data.thrtype==4)
                    data.CLUSTER_selected(n1)=nb1;
                    data.CLUSTER_selected_names{n1}=tstr;
                end
                if data.displayallmeasures, tstr2='';
                else tstr2=tstr;
                end
                if data.mvpathrtype_isftest(data.mvpathrtype)&&data.mvpathrtype_iscluster(data.mvpathrtype)
                    nb2=find(data.cMVPAlabel==nb1);
                    txt4{end+1}=sprintf('%-24s  %-20s  %12.6f  %12.6f',tstr,[data.cMVPAstatsname{nb2},data.cMVPAdofstr{nb2},' = ',num2str(data.cMVPAF(nb2),'%0.2f')],data.cMVPAp(nb2),data.cMVPAP(nb2));
                    tp4=cat(1,tp4,data.cMVPAF(nb2));
                    index4(end+1)=nb1;
                elseif isempty(data.iPERM),
                    if data.thrtype==3||data.thrtype==4,
                        txt4{end+1}=sprintf('%-24s  %-20s',tstr,['TFCE = ',num2str(maxZ(nb1),'%0.2f')]);
                        tp4=cat(1,tp4,maxZ(nb1));
                        index4(end+1)=nb1;
                    else
                        if data.mvpathrtype_isce(data.mvpathrtype)||data.displayallmeasures
                            txt4{end+1}=sprintf('%-24s  %-20s',tstr,['Score = ',num2str(sclL(nb1),'%0.2f')]);
                            tp4=cat(1,tp4,sclL(nb1));
                            index4(end+1)=nb1;
                        end
                        if data.mvpathrtype_ismass(data.mvpathrtype)||data.displayallmeasures
                            txt4{end+1}=sprintf('%-24s  %-20s',tstr2,['Mass = ',num2str(mclL(nb1),'%0.2f')]);
                            tp4=cat(1,tp4,mclL(nb1));
                            index4(end+1)=nb1;
                        end
                        if data.mvpathrtype_issize(data.mvpathrtype)||data.displayallmeasures
                            txt4{end+1}=sprintf('%-24s  %-20s',tstr2,['Size = ',num2str(nclL(nb1),'%d')]);
                            tp4=cat(1,tp4,mclL(nb1));
                            index4(end+1)=nb1;
                        end
                    end
                else
                    if data.thrtype==3||data.thrtype==4,
                        txt4{end+1}=sprintf('%-24s  %-20s  %12.6f  %12.6f  %12.6f',tstr,['TFCE = ',num2str(maxZ(nb1),'%0.2f')],PERMp_tfce_unc(nb1),PERMp_tfce_FDR(nb1),PERMp_tfce_FWE(nb1));
                        tp4=cat(1,tp4,maxZ(nb1));
                        index4(end+1)=nb1;
                    else
                        if data.mvpathrtype_isce(data.mvpathrtype)||data.displayallmeasures
                            txt4{end+1}=sprintf('%-24s  %-20s  %12.6f  %12.6f  %12.6f',tstr,['Score = ',num2str(sclL(nb1),'%0.2f')],PERMp_cluster_score_unc(nb1),PERMp_cluster_score_FDR(nb1),PERMp_cluster_score_FWE(nb1));
                            tp4=cat(1,tp4,sclL(nb1));
                            index4(end+1)=nb1;
                        end
                        if data.mvpathrtype_ismass(data.mvpathrtype)||data.displayallmeasures
                            txt4{end+1}=sprintf('%-24s  %-20s  %12.6f  %12.6f  %12.6f',tstr2,['Mass = ',num2str(mclL(nb1),'%0.2f')],PERMp_cluster_mass_unc(nb1),PERMp_cluster_mass_FDR(nb1),PERMp_cluster_mass_FWE(nb1));
                            tp4=cat(1,tp4,mclL(nb1));
                            index4(end+1)=nb1;
                        end
                        if data.mvpathrtype_issize(data.mvpathrtype)||data.displayallmeasures
                            txt4{end+1}=sprintf('%-24s  %-20s  %12.6f  %12.6f  %12.6f',tstr2,['Size = ',num2str(nclL(nb1),'%d')],PERMp_cluster_size_unc(nb1),PERMp_cluster_size_FDR(nb1),PERMp_cluster_size_FWE(nb1));
                            tp4=cat(1,tp4,mclL(nb1));
                            index4(end+1)=nb1;
                        end
                    end
                end
            end
        end
        
        if ~data.mvpathrtype_isnonparam(data.mvpathrtype), sortmeasure=mvpaP-1e-10*abs(data.MVPAF);
        elseif data.mvpathrtype_isce(data.mvpathrtype), sortmeasure=PERMp_seed_score_unc-1e-10*ssdL;
        elseif data.mvpathrtype_ismass(data.mvpathrtype), sortmeasure=PERMp_seed_mass_unc-1e-10*msdL;
        else sortmeasure=PERMp_seed_size_unc-1e-10*nsdL;
        end
        %switch(data.mvpathrtype) % sorting seeds for display
        %    case {1,20,21}, sortmeasure=mvpaP-1e-10*abs(data.MVPAF);
        %    case {2,3,4,8,9,10,14,15,16}, sortmeasure=PERMp_seed_mass_unc-1e-10*msdL;
        %    case {5,6,7,11,12,13,17,18,19}, sortmeasure=PERMp_seed_size_unc-1e-10*nsdL;
        %end
        sortedroinumbers=zeros(size(sortmeasure));
        [nill,tidx]=sort(sortmeasure(data.displaytheserois));
        sortedroinumbers(data.displaytheserois(tidx))=1:numel(tidx);
        for nt1=1:numel(tidx),
            na1=tidx(nt1);
            n1=data.displaytheserois(na1);
            if n1<=N&&seedz(n1)>0, % seed info
                %if data.displayroilabelsinstats, tstr=sprintf('%s',data.names2reduced{n1}); % (%d) na1 
                %else tstr=sprintf('ROI %d/%d',nt1,numel(tidx));
                %end
                tstr0=sprintf('ROI %d/%d',sortedroinumbers(n1),numel(tidx));
                tstr=sprintf('%s\\\\\\\\ %-32s\\\\',tstr0,data.names2reduced{n1}); 
                tstr=[tstr repmat(' ',1,max(0,24-numel(tstr0)))];
                if data.mvpathrtype_isroi(data.mvpathrtype)
                    data.CLUSTER_selected(nt1)=na1;
                    data.CLUSTER_selected_names{nt1}=tstr0;
                end
                %if numel(tstr)>24, tstr=[tstr(1:24),'..']; end
                if data.displayallmeasures, tstr2='';
                else tstr2=tstr;
                end
                if ~data.mvpathrtype_isnonparam(data.mvpathrtype)
                    %txt3{end+1}=sprintf('%-6s %6s  %6s  %6.2f  %4d  %12.6f  %12.6f',['(',num2str(na1),')'],'*   ','',data.MVPAF(n1),data.MVPAdof(n1,end),data.MVPAp(n1),data.MVPAP(n1));
                    txt3{end+1}=sprintf('%-24s  %-20s  %12.6f  %12.6f',tstr,[data.MVPAstatsname,data.MVPAdofstr{n1},' = ',num2str(data.MVPAF(n1),'%0.2f')],data.MVPAp(n1),data.MVPAP(n1));
                    txt3{end}=[parse_html1{1},txt3{end},parse_html1{2}];
                    tp3=cat(1,tp3,sortmeasure(n1)); sort3=cat(1,sort3,n1);
                    index3(end+1)=na1;
                else
                    if isempty(data.iPERM)%||~permenabled
                        if data.mvpathrtype_isce(data.mvpathrtype)||data.displayallmeasures
                            txt3{end+1}=sprintf('%-24s  %-20s',tstr,['Score = ',num2str(ssdL(n1),'%0.2f')]);
                            txt3{end}=[parse_html1{1},txt3{end},parse_html1{2}];
                            tp3=cat(1,tp3,sortmeasure(n1)); sort3=cat(1,sort3,n1);
                            index3(end+1)=na1;
                        end
                        if data.mvpathrtype_ismass(data.mvpathrtype)||data.displayallmeasures
                            txt3{end+1}=sprintf('%-24s  %-20s',tstr2,['Mass = ',num2str(msdL(n1),'%0.2f')]);
                            txt3{end}=[parse_html1{1},txt3{end},parse_html1{2}];
                            tp3=cat(1,tp3,sortmeasure(n1)); sort3=cat(1,sort3,n1);
                            index3(end+1)=na1;
                        end
                        if data.mvpathrtype_issize(data.mvpathrtype)||data.displayallmeasures
                            txt3{end+1}=sprintf('%-24s  %-20s',tstr2,['Size = ',num2str(nsdL(n1),'%d')]);
                            txt3{end}=[parse_html1{1},txt3{end},parse_html1{2}];
                            tp3=cat(1,tp3,sortmeasure(n1)); sort3=cat(1,sort3,n1);
                            index3(end+1)=na1;
                        end
                    else
                        if data.mvpathrtype_isce(data.mvpathrtype)||data.displayallmeasures
                            txt3{end+1}=sprintf('%-24s  %-20s  %12.6f  %12.6f  %12.6f',tstr,['Score = ',num2str(ssdL(n1),'%0.2f')],PERMp_seed_score_unc(n1),PERMp_seed_score_FDR(n1),PERMp_seed_score_FWE(n1));
                            txt3{end}=[parse_html1{1},txt3{end},parse_html1{2}];
                            tp3=cat(1,tp3,sortmeasure(n1)); sort3=cat(1,sort3,n1);
                            index3(end+1)=na1;
                        end
                        if data.mvpathrtype_ismass(data.mvpathrtype)||data.displayallmeasures
                            txt3{end+1}=sprintf('%-24s  %-20s  %12.6f  %12.6f  %12.6f',tstr2,['Mass = ',num2str(msdL(n1),'%0.2f')],PERMp_seed_mass_unc(n1),PERMp_seed_mass_FDR(n1),PERMp_seed_mass_FWE(n1));
                            txt3{end}=[parse_html1{1},txt3{end},parse_html1{2}];
                            tp3=cat(1,tp3,sortmeasure(n1)); sort3=cat(1,sort3,n1);
                            index3(end+1)=na1;
                        end
                        if data.mvpathrtype_issize(data.mvpathrtype)||data.displayallmeasures
                            txt3{end+1}=sprintf('%-24s  %-20s  %12.6f  %12.6f  %12.6f',tstr2,['Size = ',num2str(nsdL(n1),'%d')],PERMp_seed_size_unc(n1),PERMp_seed_size_FDR(n1),PERMp_seed_size_FWE(n1));
                            txt3{end}=[parse_html1{1},txt3{end},parse_html1{2}];
                            tp3=cat(1,tp3,sortmeasure(n1)); sort3=cat(1,sort3,n1);
                            index3(end+1)=na1;
                        end
                        %nsdL(n1),PERMp_seed_size_unc(n1),PERMp_seed_size_FDR(n1),PERMp_seed_size_FWE(n1));
                    end
                end
            end
        end
        [nill,idxsort3]=sort(tp3); sort3=sort3(idxsort3,:); txt3=txt3(idxsort3); index3=index3(idxsort3);
        %[nill,idxsort3]=sort(tp3-10*(1+tp4(1+CLroi_labels(index3)))); sort3=sort3(idxsort3,:); txt3=txt3(idxsort3); index3=index3(idxsort3);
        
        tp2=[];sort2=[];txt2={};index2=[];
        sortmeasure=P-1e-10*abs(data.F); % sorting connections for display
        for na1=1:length(data.displaytheserois),
            n1=data.displaytheserois(na1);
            for na2=1:length(data.displaytheserois),
                n2=data.displaytheserois(na2);
                if n1<=N&&n1~=n2&&z(n1,n2)>0, % connection info
                    if data.mvpathrtype_isroi(data.mvpathrtype)||~isfield(data,'issymmetric')||~data.issymmetric||na1<na2
                        index2(end+1)=N*(n2-1)+n1; % indexes to z matrix
                        if nnz(z>0)<=1e5 % skip connection-level stats if above 1e5
                            %txt2{end+1}=(sprintf('%-6s %-6s  %6.2f  %6.2f  %4d  %12.6f  %12.6f',['(',num2str(na1),')'],['(',num2str(na2),')'],data.h(n1,n2),data.F(n1,n2),data.dof(n1,end),p(n1,n2),P(n1,n2)));
                            %if data.displayroilabelsinstats, tname1=sprintf('(%s)',data.names2reduced{n1}); else tname1=sprintf('(%d)',sortedroinumbers(n1)); end %na1); end
                            %if data.displayroilabelsinstats, tname2=sprintf('(%s)',data.names2reduced{n2}); else tname2=sprintf('(%d)',sortedroinumbers(n2)); end %na2); end
                            tname1=sprintf('\\\\%03d\\\\%-32s\\\\',sortedroinumbers(n1),data.names2{n1});
                            tname2=sprintf('\\\\%03d\\\\%-32s\\\\',sortedroinumbers(n2),data.names2{n2});
                            tname3=repmat(' ',1,max(0,24-numel(sprintf(' Connection %03d-%03d',sortedroinumbers(n1),sortedroinumbers(n2)))));
                            if 0,%data.mvpathrtype_isftest(data.mvpathrtype), txt2{end+1}=sprintf('%-24s  %-20s  %12.6f  %12.6f',[' ',tname1,'-',tname2,tname3],sprintf('%s%s = %.2f',data.statsname,data.dofstr{n1},data.F(n1,n2)),p(n1,n2),data.Ppos(n1,n2)); % P(n1,n2)); % note: a posteriori p-fdr values (for each ROI seed)
                            else txt2{end+1}=sprintf('%-24s  %-20s  %12.6f  %12.6f',[' Connection ',tname1,'-',tname2,tname3],sprintf('%s%s = %.2f',data.statsname,data.dofstr{n1},data.F(n1,n2)),p(n1,n2),P(n1,n2));
                            end
                            
                            %if numel(data.names2reduced{n1})<=5, tname1=data.names2reduced{n1}; else tname1=['(',num2str(na1),')']; end
                            %if numel(data.names2reduced{n2})<=5, tname2=data.names2reduced{n2}; else tname2=['(',num2str(na2),')']; end
                            %txt2{end+1}=sprintf('%-20s  %-20s  %12.6f  %12.6f',['  conn ',tname1,'-',tname2],[data.statsname,'(',num2str(data.dof(n1,end)),') = ',num2str(data.F(n1,n2),'%0.2f')],p(n1,n2),P(n1,n2));
                            if data.h(n1,n2)>=0, txt2{end}=[parse_html2{1},txt2{end},parse_html2{2}];
                            else                 txt2{end}=[parse_html3{1},txt2{end},parse_html3{2}];
                            end
                        else txt2{end+1}='--';
                        end
                        tp2=cat(1,tp2,sortmeasure(n1,n2)); sort2=cat(1,sort2,[n1,n2]);
                    end
                end
            end
        end
        [nill,idxsort2]=sort(tp2); sort2=sort2(idxsort2,:); txt2=txt2(idxsort2); index2=index2(idxsort2);
        
        data.list2=[];
        data.list2txt={};
        data.list2visible=[];
        done2=ones(size(index2));
        done3=ones(size(index3));
        done4=ones(size(index4));
        if data.thrtype==3|data.thrtype==4|data.mvpathrtype_iscluster(data.mvpathrtype)|data.mvpathrtype_isnetwork(data.mvpathrtype), % cluster-level     %data.mvpasortresultsby==3,%data.PERMenabled&&mvpaenabled&&data.mvpathrmeasure>3
            for curcl=1:numel(txt4)
                if done4(curcl)>0
                    i=find(index4==index4(curcl));
                    if ~isempty(i)
                        data.list2txt=cat(2,data.list2txt,txt4(i));
                        data.list2=cat(1,data.list2,[zeros(numel(i),2),index4(curcl)+zeros(numel(i),1)]);
                        data.list2visible=cat(1,data.list2visible,ones(numel(i),1));
                        done4(i)=0;
                    end
                    idx2=find(ismember(index2,find(CLUSTER_labels==index4(curcl))));
                    if ~isempty(idx2),
                        done2(idx2)=0;
                        data.list2txt=cat(2,data.list2txt,txt2(idx2));
                        data.list2=cat(1,data.list2,[sort2(idx2,:),index4(curcl)+zeros(numel(idx2),1)]);
                        data.list2visible=cat(1,data.list2visible,zeros(numel(idx2),1));
                    end
%                     for na1=find(CLroi_labels(data.displaytheserois(index3))==index4(curcl))'
%                         if done3(na1)>0
%                             i=find(index3==index3(na1));
%                             if ~isempty(i)
%                                 done3(i)=0;
%                                 data.list2txt=cat(2,data.list2txt,{sprintf(' %s (%d)',data.names2reduced{sort3(i(1))},index3(na1))});
%                                 data.list2=cat(1,data.list2,[sort3(i(1)),0,index4(curcl)]);
%                                 data.list2visible=cat(1,data.list2visible,0);
%                                 data.list2txt=cat(2,data.list2txt,txt3(i));
%                                 data.list2=cat(1,data.list2,[sort3(i),zeros(numel(i),1),index4(curcl)+zeros(numel(i),1)]);
%                                 data.list2visible=cat(1,data.list2visible,zeros(numel(i),1));
%                             end
%                             idx2=find(index2==index3(na1)&done2);
%                             if ~isempty(idx2)
%                                 done2(idx2)=0;%index2(idx2)=0;
%                                 data.list2txt=cat(2,data.list2txt,txt2(idx2));
%                                 data.list2=cat(1,data.list2,[sort2(idx2,:),index4(curcl)+zeros(numel(idx2),1)]);
%                                 data.list2visible=cat(1,data.list2visible,zeros(numel(idx2),1));
%                             end
%                         end
%                     end
                end
            end
        elseif data.mvpathrtype_isroi(data.mvpathrtype), % ROI-level     %data.mvpasortresultsby==2
            for na1=1:numel(txt3)
                if done3(na1)>0
                    i=find(index3==index3(na1));
                    if ~isempty(i)
                        %data.list2txt=cat(2,data.list2txt,{sprintf(' %s (%d)',data.names2reduced{sort3(i(1))},index3(na1))});
                        %data.list2=cat(1,data.list2,[sort3(i(1)),0,0]);
                        %data.list2visible=cat(1,data.list2visible,1);
                        done3(i)=0;
                        data.list2txt=cat(2,data.list2txt,txt3(i));
                        data.list2=cat(1,data.list2,[sort3(i),zeros(numel(i),2)]);
                        data.list2visible=cat(1,data.list2visible,ones(numel(i),1));
                    end
                    idx2=find(ceil(index2/N)==data.displaytheserois(index3(na1))&done2);
                    if ~isempty(idx2),
                        done2(idx2)=0;
                        data.list2txt=cat(2,data.list2txt,txt2(idx2));
                        data.list2=cat(1,data.list2,[sort2(idx2,:),zeros(numel(idx2),1)]);
                        data.list2visible=cat(1,data.list2visible,zeros(numel(idx2),1));
                    end
                end
            end
        else
            data.list2txt=txt2;
            data.list2=[sort2,zeros(size(sort2,1),1)];
            data.list2visible=ones(size(sort2,1),1);
            done2(:)=0;
        end
        idx2=find(done2>0);
        if ~isempty(idx2)
            data.list2txt=cat(2,data.list2txt,txt2(idx2));
            index2(idx2)=0;
            data.list2=cat(1,data.list2,[sort2(idx2,:),zeros(numel(idx2),1)]);
            data.list2visible=cat(1,data.list2visible,ones(numel(idx2),1));
        end
        data.list2visible=find(data.list2visible);
        
%         [nill,idxsort]=sort(tp2); data.list2=sort2(idxsort,:); txt2=strvcat(txt2{idxsort},' ');
%         
        %if isequal(data.source,0), set(data.handles(6),'string',txt1,'value',max(1,min(size(txt1,1), 1:numel(data.displaytheserois))));
        %else set(data.handles(6),'string',txt1,'value',max(1,min(size(txt1,1), max(1,unique(data.source))))); %get(data.handles(6),'value'))));
        %end
        %if get(data.handles(6),'listboxtop')>size(get(data.handles(6),'string'),1), set(data.handles(6),'listboxtop',1); end
        if isempty(data.list2txt), set(data.handles(8),'string',{' no significant results '} ,'value',[]);
        else
            if data.displayroilabelsinstats, tstr=regexprep(data.list2txt,'\\\\(\d*)\\\\(.*?)\\\\','$2');
            else tstr=regexprep(data.list2txt,'\\\\(\d*)\\\\(.*?)\\\\','$1');
            end
            if data.displayconnectionstats, set(data.handles(8),'string',[tstr {' '}] ,'value',max(1,min(numel(data.list2txt)+1, get(data.handles(8),'value'))));
            else set(data.handles(8),'string',[tstr(data.list2visible) {' '}] ,'value',max(1,min(numel(data.list2visible)+1, get(data.handles(8),'value'))));
            end
        end
        if get(data.handles(8),'listboxtop')>size(get(data.handles(8),'string'),1), set(data.handles(8),'listboxtop',1); end
        if ~isfield(data,'displaylabels'),data.displaylabels=0;end
        if ~isfield(data,'displaybrains'),data.displaybrains=0;end
%         if ~isfield(data,'clusters_options')||isempty(data.clusters_options), 
%             data.clusters_options=struct('type','hc','groups',nan,'param',.05);
%             data=conn_displayroi_clusters(data); 
%             if ~data.view
%                 data.x=data.xy2*data.proj(1:2,1);data.y=data.xy2*data.proj(1:2,2);data.z=zeros(size(data.xy2,1),1);
%             else
%                 data.x=data.xyz2*data.proj(:,1);data.y=data.xyz2*data.proj(:,2);data.z=data.xyz2*data.proj(:,3);
%                 lim=[1,1,1;data.ref.dim];refminmax=sort([lim((dec2bin(0:7)-'0'+1)+repmat([0,2,4],[8,1])),ones(8,1)]*data.ref.mat(1:3,:)'*data.proj(:,1:2));
%                 [data.bgx,data.bgy]=meshgrid(refminmax(1,1):scale:refminmax(end,1),refminmax(1,2):scale:refminmax(end,2));
%                 data.bgimage=spm_get_data(data.ref,pinv(data.ref.mat)*[data.proj(:,1:3)*[data.bgx(:),data.bgy(:),data.bgz+zeros(prod(size(data.bgx)),1)]';ones(1,prod(size(data.bgx)))]);
%             end
%         end
        
        % plots
        EPS=1e-0;
        SINGLEBRAIN=data.displaybrains<=1;        
        LABELONSIGNONLY=data.displaylabels<1;
        LABELONALL=data.displaylabels>1;
        SQUAREALLROIS=false;
        OFFSET=.15;
        EMPH=2+2*(any(~data.source)|numel(data.source)>1);
        
        
%         figure(hfig);set(hfig,'pointer','watch');%drawnow;
%         th1=axes('units','norm','position',[.03,.06,.45,.88]);th2=patch([0 0 1 1],[0 1 1 0],'k');set(th2,'facealpha',.5);set(th1,'xlim',[0 1],'ylim',[0 1],'visible','off'); drawnow; delete([th1 th2]);
        if ~data.view, set(hfig,'color',data.plotconnoptions.BCOLOR); else set(hfig,'color','w'); end
        rcircle=(1+.25*(data.displaylabels>0)*(data.view>0))*[sin(linspace(0,2*pi,64)'),cos(linspace(0,2*pi,64))']*diag([5,5]);
        rtriang=[-1-1i*.5;0;-1+1i*.5];
        rsquare=[0,-.5;0,.5;.05,.5;.05,-.5]*diag([20,min(50,1.0*2*pi*200/length(data.displaytheserois)/1.125)]);
        rsquarelong=[0,-.5;0,.5;1,.5;1,-.5]*diag([20,min(50,1.0*2*pi*200/length(data.displaytheserois)/1.125)]);
        %rsquareshort=([.5+.3*sin(linspace(0,2*pi,16)') .5*cos(linspace(0,2*pi,16)')])*diag([20,min(50,1.0*2*pi*200/length(data.displaytheserois)/1.125)]);
        rsquareshort=([0,-.5;0,.5;.8,.5;.8,-.5])*diag([20,min(50,1.0*2*pi*200/length(data.displaytheserois)/1.125)]);        
        temp=200+rsquare(:,1)+1i*rsquare(:,2);temp([3,4])=temp([2,1])*real(temp(4))./real(temp(2));rsquare=[real(temp)-200, imag(temp)];
        temp=200+rsquareshort(:,1)+1i*rsquareshort(:,2);temp([3,4])=temp([2,1])*real(temp(4))./real(temp(2));rsquareshort=[real(temp)-200, imag(temp)];
        temp=200+rsquarelong(:,1)+1i*rsquarelong(:,2);temp([3,4])=temp([2,1])*real(temp(4))./real(temp(2));rsquarelong=[real(temp)-200, imag(temp)];
        if ishandle(data.plotaxes), delete(data.plotaxes); end
        h=findobj(hfig,'tag','conn_displayroi_plot'); if ~isempty(h),delete(h); end
        if data.displaygui, data.plotaxes=axes('units','norm','position',data.plotposition{1},'visible','off','parent',data.hfig);
        else data.plotaxes=axes('units','norm','position',data.plotposition{2},'visible','off','parent',data.hfig);
        end
        if isfield(data,'axeslim')&&~isempty(data.axeslim), set(data.plotaxes,'xlim',data.axeslim{1},'ylim',data.axeslim{2}); end
        if data.view
            lim=[1,1,1;data.ref.dim];refminmax=sort([lim((dec2bin(0:7)-'0'+1)+repmat([0,2,4],[8,1])),ones(8,1)]*data.ref.mat(1:3,:)'*data.proj(:,1:2));
            temp=reshape(data.bgimage,size(data.bgx)); %convn(convn(reshape(data.bgimage,size(data.bgx)),conn_hanning(5),'same'),conn_hanning(5)','same'); 
            temp(isnan(temp))=0;
            temp=round(1+(1-.2*temp/max(temp(:)))*(size(get(hfig,'colormap'),1)-1));
            data.refaxes=image(refminmax(1,1):scale:refminmax(end,1),refminmax(1,2):scale:refminmax(end,2),temp,'parent',data.plotaxes);hold(data.plotaxes,'on');
            set(data.refaxes,'cdatamapping','direct');
        else
            if 0,%~SQUAREALLROIS&&data.displaybrains<2, % ring circle reference
                xy0=200*(cos(2*pi*linspace(0,1,1e3))'+1i*sin(2*pi*linspace(0,1,1e3))');
                patch(1.08*real(xy0),1.08*imag(xy0),-20+zeros(size(xy0)),'w','facecolor',.1+.8*data.plotconnoptions.BCOLOR,'edgecolor',.1+.8*data.plotconnoptions.BCOLOR,'parent',data.plotaxes);
                patch(1*real(xy0),1*imag(xy0),-10+zeros(size(xy0)),'w','facecolor',data.plotconnoptions.BCOLOR,'edgecolor',data.plotconnoptions.BCOLOR,'parent',data.plotaxes);
            end
            data.refaxes=data.plotaxes;
            set(data.plotaxes,'color',data.plotconnoptions.BCOLOR,'xcolor',data.plotconnoptions.BCOLOR,'ycolor',data.plotconnoptions.BCOLOR,'buttondownfcn',@conn_displayroi_menubuttondownfcn,'userdata',data.hfig);
        end
        data.buttondown=struct('h1',data.plotaxes);set(data.plotaxes,'tag','conn_displayroi_plot');
        if ~data.plotconnoptions.menubar||isempty(findobj(hfig,'type','uimenu','-and','-not','tag','donotdelete'))
            if data.plotconnoptions.menubar, hc1=hfig;delete(findobj(hfig,'type','uimenu','-and','-not','tag','donotdelete'));%uimenu(hfig,'Label','Display options');
            else hc1=uicontextmenu;delete(findobj(hfig,'type','uimenu','-and','-not','tag','donotdelete'));
            end
            %set(hc1,'tag','conn_displayroi_plot','callback',@(varargin)set(data.handles(20),'visible','off'));
            if data.plotconnoptions.menubar
                if data.view>0
                    ht=uimenu(hc1,'Label','View');
                    uimenu(ht,'Label','connectome ring','callback',{@conn_displayroi,'view-ring'});
                    uimenu(ht,'Label','axial view (x-y)','callback',{@conn_displayroi,'view-axial'});
                    uimenu(ht,'Label','coronal view (x-z)','callback',{@conn_displayroi,'view-coronal'});
                    uimenu(ht,'Label','sagittal view (y-z)','callback',{@conn_displayroi,'view-sagittal'});
                    uimenu(ht,'Label','3d display','callback',{@conn_displayroi,'display3d'});
                    ht=uimenu(hc1,'Label','ROIs');
                    uimenu(ht,'Label','show ROI labels','callback',{@conn_displayroi,'labelson'});
                    uimenu(ht,'Label','hide ROI labels','callback',{@conn_displayroi,'labelsoff'});
                    uimenu(ht,'Label','increase labels fontsize','callback',{@conn_displayroi,'labels1'});
                    uimenu(ht,'Label','decrease labels fontsize','callback',{@conn_displayroi,'labels2'});
                    uimenu(ht,'Label','labels: set labels fontsize','callback',{@conn_displayroi,'labels3'});
                    uimenu(ht,'Label','edit ROI labels','callback',{@conn_displayroi,'labelsedit'});
                    ht=uimenu(hc1,'Label','Connections');
                    uimenu(ht,'Label','color: positive/negative = red/blue','callback',{@conn_displayroi,'edgecolors1'});
                    uimenu(ht,'Label','color: proportional to stats (rgb colormap)','callback',{@conn_displayroi,'edgecolors3'});
                    uimenu(ht,'Label','color: increase rgb colormap contrast','callback',{@conn_displayroi,'edgecolors4'});
                    uimenu(ht,'Label','color: decrease rgb colormap contrast','callback',{@conn_displayroi,'edgecolors5'});
                    uimenu(ht,'Label','thickness: fixed width','callback',{@conn_displayroi,'edgewidths1'},'separator','on');
                    uimenu(ht,'Label','thickness: proportional to stats','callback',{@conn_displayroi,'edgewidths2'});
                    uimenu(ht,'Label','thickness: increase','callback',{@conn_displayroi,'edgewidths3'});
                    uimenu(ht,'Label','thickness: decrease','callback',{@conn_displayroi,'edgewidths4'});
                    %uimenu(ht,'Label','arrow-widths scaled by T-values','callback',{@conn_displayroi,'displayefffectsize-off'});
                    %uimenu(ht,'Label','arrow-widths scaled by beta-values','callback',{@conn_displayroi,'displayefffectsize-on'});
                    %uimenu(ht,'Label','fixed arrow-widths','callback',{@conn_displayroi,'displayefffectsize-none'});
                    ht=uimenu(hc1,'Label','Options');
                    uimenu(ht,'Label','Change background anatomical image','callback',{@conn_displayroi,'changebackground'});
                    %ht=uimenu(hc1,'Label','Menubar');
                    %uimenu(ht,'Label','on','callback',['set(gcbf,''menubar'',''figure'');']);
                    %uimenu(ht,'Label','off','callback',['set(gcbf,''menubar'',''none'');']);
                    uimenu(ht,'Label','Refresh GUI','callback',{@conn_displayroi,'refresh'});
                    uimenu(ht,'Label','Pause/Resume GUI','callback',{@conn_displayroi,'pausegui'});
                    uimenu(ht,'Label','Show/Hide extended cluster threshold options','callback',{@conn_displayroi,'mvpaextend'});
                    uimenu(ht,'Label','Show/Hide GUI','callback',{@conn_displayroi,'displaygui'});
                    uimenu(ht,'Label','Print (high-res)','callback',{@conn_displayroi,'print'});
                else
                    %ht=uimenu(hc1,'Label','View');
                    %uimenu(ht,'Label','connectome ring','callback',{@conn_displayroi,'view-ring'});
                    %uimenu(ht,'Label','axial view (x-y)','callback',{@conn_displayroi,'view-axial'});
                    %uimenu(ht,'Label','coronal view (x-z)','callback',{@conn_displayroi,'view-coronal'});
                    %uimenu(ht,'Label','sagittal view (y-z)','callback',{@conn_displayroi,'view-sagittal'});
                    %uimenu(ht,'Label','3d display','callback',{@conn_displayroi,'display3d'});
                    ht=uimenu(hc1,'Label','ROIs');
                    %ht=uimenu(hc1,'Label','ROIs display order');
                    %%uimenu(ht,'Label','order: export current order to file or to other figure','callback',{@conn_displayroi,'roi.order.export'});
                    %%uimenu(ht,'Label','order: import ROI order from file or from other figure','callback',{@conn_displayroi,'roi.order.import'});
                    %%uimenu(ht,'Label','order: change to list order','callback',{@conn_displayroi,'clusters','none'});
                    %ht1=uimenu(ht,'Label','Define new ROI ordering criterium');
                    %%uimenu(ht,'Label','order: change to hierarchical clustering alg. (default; contiguous ROIs are functionally similar)','callback',{@conn_displayroi,'clusters','hc'});
                    %%uimenu(ht,'Label','order ROIs by networks (minimum degree algorithm)','callback',{@conn_displayroi,'clusters','amd'});
                    %uimenu(ht,'Label','order: change to minimum degree alg. (contiguous ROIs are in the same subgraph)','callback',{@conn_displayroi,'clusters','hcnet'});
                    %uimenu(ht,'Label','order: change to reverse Cuthill-McKee alg. (minimize supra-threshold connection lengths)','callback',{@conn_displayroi,'clusters','symrcm'});
                    %ht=uimenu(hc1,'Label','ROI labels');
                    %%uimenu(ht,'Label','view: show ROIs in ring only','callback',{@conn_displayroi,'brainsoff'},'separator','on');
                    %%uimenu(ht,'Label','view: show ROIs in reference brain displays only','callback',{@conn_displayroi,'brainson'});
                    %%uimenu(ht,'Label','view: show ROIs in ring and in reference brain displays','callback',{@conn_displayroi,'brainssingle'});
                    uimenu(ht,'Label','view: show reference brain displays','callback',{@conn_displayroi,'brainsmany'});
                    uimenu(ht,'Label','view: show reference brain display','callback',{@conn_displayroi,'brainssingle'});
                    %uimenu(ht,'Label','view: show group-reference brain displays','callback',{@conn_displayroi,'brainson'});
                    uimenu(ht,'Label','view: hide reference brain display(s)','callback',{@conn_displayroi,'brainsoff'});
                    uimenu(ht,'Label','labels: show ROI labels','callback',{@conn_displayroi,'labelson'});%,'separator','on');
                    %uimenu(ht,'Label','labels: show ROI labels for relevant ROIs only','callback',{@conn_displayroi,'labelspartial'});
                    uimenu(ht,'Label','labels: hide ROI labels','callback',{@conn_displayroi,'labelsoff'});
                    uimenu(ht,'Label','labels: increase labels fontsize','callback',{@conn_displayroi,'labels1'});
                    uimenu(ht,'Label','labels: decrease labels fontsize','callback',{@conn_displayroi,'labels2'});
                    uimenu(ht,'Label','labels: set labels fontsize','callback',{@conn_displayroi,'labels3'});
                    uimenu(ht,'Label','labels: edit ROI labels','callback',{@conn_displayroi,'labelsedit'});
                    uimenu(ht,'Label','labels: edit Group labels','callback',{@conn_displayroi,'groupsedit'});
                    ht=uimenu(hc1,'Label','Connections');
                    uimenu(ht,'Label','type: display connections as lines','callback',{@conn_displayroi,'displaytype',0});
                    uimenu(ht,'Label','type: display connectivity matrix','callback',{@conn_displayroi,'displaytype',1});
                    uimenu(ht,'Label','type: display connectivity polar matrix','callback',{@conn_displayroi,'displaytype',2});
                    uimenu(ht,'Label','color: positive/negative = red/blue','callback',{@conn_displayroi,'edgecolors1'},'separator','on');
                    uimenu(ht,'Label','color: ring colorwheel','callback',{@conn_displayroi,'edgecolors2'});
                    uimenu(ht,'Label','color: proportional to stats (rgb colormap)','callback',{@conn_displayroi,'edgecolors3'});
                    uimenu(ht,'Label','color: increase rgb colormap contrast','callback',{@conn_displayroi,'edgecolors4'});
                    uimenu(ht,'Label','color: decrease rgb colormap contrast','callback',{@conn_displayroi,'edgecolors5'});
                    uimenu(ht,'Label','thickness: fixed width','callback',{@conn_displayroi,'edgewidths1'},'separator','on');
                    uimenu(ht,'Label','thickness: proportional to stats','callback',{@conn_displayroi,'edgewidths2'});
                    uimenu(ht,'Label','thickness: increase','callback',{@conn_displayroi,'edgewidths3'});
                    uimenu(ht,'Label','thickness: decrease','callback',{@conn_displayroi,'edgewidths4'});
                    uimenu(ht,'Label','opacity: fixed opacity','callback',{@conn_displayroi,'edgeopacity1'},'separator','on');
                    uimenu(ht,'Label','opacity: proportional to stats','callback',{@conn_displayroi,'edgeopacity2'});
                    uimenu(ht,'Label','opacity: increase','callback',{@conn_displayroi,'edgeopacity3'});
                    uimenu(ht,'Label','opacity: decrease','callback',{@conn_displayroi,'edgeopacity4'});
                    ht=uimenu(hc1,'Label','Options');
                    uimenu(ht,'Label','Advanced display options','callback',{@conn_displayroi,'displayoptions'});
                    %ht=uimenu(hc1,'Label','Menubar');
                    %uimenu(ht,'Label','on','callback',['set(gcbf,''menubar'',''figure'');']);
                    %uimenu(ht,'Label','off','callback',['set(gcbf,''menubar'',''none'');']);
                    uimenu(ht,'Label','Refresh GUI','callback',{@conn_displayroi,'refresh'});
                    %uimenu(ht,'Label','Pause/Resume interactive ring display','callback',{@conn_displayroi,'pausegui'});
                    uimenu(ht,'Label','Switch axial/sagittal/coronal display','callback',{@conn_displayroi,'display.viewcycle'});
                    uimenu(ht,'Label','Show/Hide advanced thresholding options','callback',{@conn_displayroi,'mvpaextend'});
                    uimenu(ht,'Label','Show/Hide all GUI menus','callback',{@conn_displayroi,'displaygui'});
                    %uimenu(ht,'Label','Print (high-res)','callback',{@conn_displayroi,'print'});
                end
                if ~data.plotconnoptions.menubar
                    set(data.refaxes,'uicontextmenu',hc1);
                    set(hfig,'uicontextmenu',hc1);
                end
            end
        end
        
        if ~data.view&&isfield(data,'displaybrains')&&data.displaybrains>=2
            if isempty(data.clusters), NPLOTS=data.plotconnoptions.NPLOTS;
            else NPLOTS=max(data.clusters); 
            end
            if SINGLEBRAIN, offset=OFFSET+.5;
            else offset=OFFSET+data.plotconnoptions.DOFFSET*(data.displaylabels>0);
            end
            xy=data.x+1i*data.y;
            mr=mean(data.plotconnoptions.rende.vertices,1);
            data.xb=zeros(size(data.x));
            data.yb=zeros(size(data.y));
            data.zb=zeros(size(data.y));
            data.bclusters=zeros(size(data.x));
            for n1=1:NPLOTS,
                if isempty(data.clusters), idx=find(abs(angle(xy.*exp(-1i*2*pi/NPLOTS*n1)))<2*pi/NPLOTS/2);
                    %idx=find(max(1,min(NPLOTS,ceil((angle(xy)+pi)/2/pi*NPLOTS)))==n1);
                else idx=find(data.clusters==n1); 
                end
                if ~isempty(idx)
                    mx=mean(xy(idx));
                    px=cos(angle(mx));
                    py=sin(angle(mx));
                    if ~data.plotconnoptions.nprojection
                        cprojection=cellfun(@(x)sum(std(data.xyz2(idx,:)*x(1:2,:)',1,1).^2),data.plotconnoptions.Projections);
                        cprojection(2)=nan;
                        [nill,nprojection]=max(cprojection);
                    else
                        nprojection=data.plotconnoptions.nprojection;
                    end
                    p=data.plotconnoptions.Projections{nprojection}'*(1+offset)*data.plotconnoptions.BSCALE;
                    if nprojection==1&&mean(data.xyz2(idx,1)>0)>.5, p(:,1)=-p(:,1); end
                    xy2=[data.xyz2(idx,:)-mr(ones(numel(idx),1),:),ones(numel(idx),1)]*[p;(1+offset)*(1+.75*data.plotconnoptions.BSCALE)*200*px,(1+offset)*(1+.75*data.plotconnoptions.BSCALE)*200*py,0];
                    data.xb(idx)=xy2(:,1);
                    data.yb(idx)=xy2(:,2);
                    data.zb(idx)=xy2(:,3);
                    data.bclusters(idx)=n1;
                end
            end
        end
        
        %data.plotsadd2=cell(1,size(data.list2,1));
        %data.plotsadd3=cell(1,numel(data.names2));
        idxtext=[];idxtexthl=[];
        hold(data.plotaxes,'on');
        K=[];KK=[];%J=0;
        for na1=1:length(data.displaytheserois),%size(z,2),%N2,
            n1=data.displaytheserois(na1);
            temp1=data.h(z(:,n1)>0,n1);
            temp2=data.F(z(:,n1)>0,n1);
            if n1<=size(z,1), temp1=[temp1;data.h(n1,z(n1,:)>0)']; temp2=[temp2;data.F(n1,z(n1,:)>0)']; end
            if isempty(temp1), K(na1)=0;
            elseif data.plotconnoptions.LCOLOR==1,
                K(na1)=sum(temp1>0)-sum(temp1<0);
            else
                %K(na1)=mean(sign(temp1).*abs(temp1).^data.plotconnoptions.LCOLORSCALE);
                K(na1)=mean(sign(temp1).*abs(temp2).^data.plotconnoptions.LCOLORSCALE);
            end
            if ~(data.mvpaenablethr&&n1<=numel(seedz)&&seedz(n1)>0||data.enablethr&&(any(z(:,n1)>0)||data.enablethr&&n1<=size(z,1)&&any(z(n1,:)))), K(na1)=nan; end
            KK(n1)=K(na1);
            %J=max([J;abs(data.h(z(:,n1)>0&(1:size(z,1))'~=n1,n1))]);
            %if na1<=size(z,1), K(na1)=K(na1)+(sum(data.h(n1,z(n1,:)>0)>0,2)'-sum(data.h(n1,z(n1,:)>0)<0,2)');
            %end
        end
        KKscale=max(eps,max(abs(K)));
        Kscale=max(eps,max(abs(K)));
        KK=KK/KKscale;
        K=K/Kscale;
        %K(K>0)=K(K>0)/max(eps,max(K));
        %K(K<0)=K(K<0)/max(eps,max(-K));
        semic0=pi/3;
        semic1=linspace(-1,1,64)+1i*(cos(linspace(-semic0,semic0,64))-cos(semic0));
        semic2=linspace(-1,1,64)-1i*(cos(linspace(-semic0,semic0,64))-cos(semic0));
        ssemic1=[semic1,fliplr(semic1)];
        semic1factor1=imag(semic1).^.1;
        semic1factor2=(imag(ssemic1)*2).^.05;
        dwidthfactor1=([.5,min(5,1*abs(rsquare(1,2)))]*[abs(real(semic1)).^2;1-abs(real(semic1)).^2])/2;
        rsemic1=struct('vertices',zeros(2*numel(semic1),3), 'faces',[reshape([1:numel(semic1)-1; 2:numel(semic1)],1,[])' reshape([2:numel(semic1); 2*numel(semic1)-1:-1:numel(semic1)+1],1,[])' reshape([2*numel(semic1):-1:numel(semic1)+2;2*numel(semic1):-1:numel(semic1)+2],1,[])']);
        sf=ones(1,length(data.displaytheserois));
        markthese=zeros(length(data.names2),1);
        if data.displaybrains==2
            datax=data.xb;
            datay=data.yb;
            dataz=data.zb;
        else
            datax=data.x;
            datay=data.y;
            dataz=data.z;
        end
        cumpatch({data.hfig,1},'init'); % lines 
        cumpatch({data.hfig,2},'init'); % brain black selector
        cumpatch({data.hfig,3},'init'); % ring black selector
%         set(data.hfig,'windowbuttonmotionfcn',[]);
        ringsquares=struct('x',[],'y',[],'index',[],'cluster',[],'h',[],'names',{{}},'gca',data.plotaxes,'gcf',data.hfig,'gct',data.handles(21));
        [nill,datarank]=sort(angle(data.xy2(data.displaytheserois,1)+1i*data.xy2(data.displaytheserois,2)));datarank(datarank)=1:numel(datarank);
        if nnz(z>0)>1e5, LINESTYLEMTX=max(1,data.plotconnoptions.LINESTYLEMTX);
        else LINESTYLEMTX=data.plotconnoptions.LINESTYLEMTX;
        end
        if LINESTYLEMTX>1
            datax(data.displaytheserois)=200*cos(-pi+2*pi*max(0,datarank-.5)/numel(datarank));
            datay(data.displaytheserois)=200*sin(-pi+2*pi*max(0,datarank-.5)/numel(datarank));
        end
        
        for na1=1:length(data.displaytheserois),%size(z,2),%N2,
            n1=data.displaytheserois(na1);
            if (~data.view&&SQUAREALLROIS)||(0&any(na1==data.source)||(n1<=size(z,1)&&(any(z(n1,:)>0)&&data.enablethr||n1<=numel(seedz)&&seedz(n1)>0&&data.mvpaenablethr))||any(z(:,n1)>0)&&data.enablethr),%((n1<=size(z,1)||~data.displayreduced)&&any(z(:,n1)>0)),
                if data.mvpathrtype>1, markthese(n1)=1; end
                %if n1<=size(z,1), sf(na1)=max(sf(na1),max(z(n1,:))); end
                %if n1<=size(z,2), sf(na1)=max(sf(na1),max(z(:,n1))); end
                %if any(na1==data.source), sf(na1)=1; else sf(na1)=sf(na1)/3; end
                if data.view>0
                    if (n1<=size(z,1)&&any(z(n1,:)>0))||any(z(:,n1)>0),%||data.z(n1)>0
                        h=patch(datax(n1)+rcircle(:,1)*max(.5,1+1e-3*data.z(n1))*sf(na1),datay(n1)+rcircle(:,2)*max(.5,1+1e-3*data.z(n1))*sf(na1),max(1,5*data.z(n1))*EPS+200+zeros(size(rcircle,1),1),'w','parent',data.plotaxes);
                        k=K(na1);%(sum(data.h(z(:,n1)>0,n1)>0,1)-sum(data.h(z(:,n1)>0,n1)<0,1))/max(eps,sum(data.h(z(:,n1)>0,n1)>0,1)+sum(data.h(z(:,n1)>0,n1)<0,1));
                        if 1,%numel(data.source)==1&&data.source~=0
                            if isnan(k), set(h,'facecolor','none');
                            else set(h,'facecolor',cmap(round(1+(size(cmap,1)-1)*(1+k)/2),:));
                            end
                            if n1<=size(z,1)&&seedz(n1)>0, set(h,'linewidth',2);
                            elseif any(na1==data.source)
                                if any(z(:,n1)>0)||(n1<=size(z,1)&&any(z(n1,:)>0)), set(h,'linewidth',2);
                                else set(h,'facecolor','w','zdata',get(h,'zdata')-max(1,5*data.z(n1))*EPS-200);%,'edgecolor','k','linewidth',1);
                                end
                            else set(h,'edgecolor','w','linewidth',1);
                            end
                        else
                            set(h,'facecolor','w');
                            if any(na1==data.source), set(h,'edgecolor','k','linewidth',2); else, set(h,'edgecolor','k','linewidth',1); end
                        end
                        set(h,'buttondownfcn',@conn_displayroi_menubuttondownfcn,'userdata',data.hfig,'interruptible','off');
                        %data.plotsadd3{n1}=h;
                    end
                else
                    if 0,%data.displaybrains<2, % ring
%                         k=K(na1);
%                         zinc=200; %200;
%                         tx=datax(n1)+rsquare*[(1+0*data.MVPAF(min(size(z,1),n1))*(n1<=size(z,1)&&seedz(n1)>=0))*datax(n1)/200;-datay(n1)/200]*max(.5,1+1e-3*data.z(n1))*sf(na1);
%                         ty=datay(n1)+rsquare*[(1+0*data.MVPAF(min(size(z,1),n1))*(n1<=size(z,1)&&seedz(n1)>=0))*datay(n1)/200;datax(n1)/200]*max(.5,1+1e-3*data.z(n1))*sf(na1);
%                         tz=max(1,min(100,5*data.z(n1)))*EPS+zinc+zeros(size(rsquare,1),1);
%                         %h=patch(tx,ty,tz,1-get(data.hfig,'color'),'edgecolor','none');
%                         %if isnan(k), set(h,'facecolor','none'); else set(h,'facecolor',cmap(round(1+(size(cmap,1)-1)*(1+k)/2),:)); end
%                         %hold on;
%                         if isnan(k), facecolor=[1 1 1];facealpha=0;
%                         else facecolor=cmap(round(1+(size(cmap,1)-1)*(1+k)/2),:);facealpha=1;
%                         end
%                         tpatch=struct('vertices',[tx,ty,tz],'faces',1:numel(tx),'facevertexcdata',repmat(facecolor,numel(tx),1),'facevertexalphadata',repmat(facealpha,numel(tx),1),'coords',{data.names2(n1)},'index',n1);
%                         h=cumpatch({data.hfig,2},tpatch,'edgecolor','none','parent',data.plotaxes);
%                         if (n1<=size(z,1)&&seedz(n1)>0&&data.mvpathrtype>1) || (((n1<=size(z,1)&&any(z(n1,:)>0))||any(z(:,n1)>0))&&data.mvpathrtype==1), markthese(n1)=2; %set(h,'facecolor',1-get(data.hfig,'color')); 
%                         elseif any(na1==data.source), %set(h,'facecolor',.4+.2*get(data.hfig,'color'));
%                         else %set(h,'facecolor',.3+.4*get(data.hfig,'color'));
%                         end
%                         %data.plotsadd3{n1}=h;
                    else
                        if n1<=size(z,1)&&seedz(n1)>0, markthese(n1)=2; end
                    end
                end
                %if all(data.h(z(:,n1)>0,n1)>0), set(h,'facecolor',[1,0,0]);%*(.5+.5*k/max(eps,max(sum(z>0,1)))));
                %elseif all(data.h(z(:,n1)>0,n1)<0), k=sum(data.h(z(:,n1)>0,n1)<0,1);set(h,'facecolor',[0,0,1]);%*(.5+.5*k/max(eps,max(sum(z>0,1)))));
                %else, set(h,'facecolor',[0,1,0]); end
            end
            if data.view==0&&data.displaybrains<2, % ring
                zinc=200; %200;
                tx=datax(n1)+rsquarelong*[datax(n1)/200;-datay(n1)/200]*sf(na1);
                ty=datay(n1)+rsquarelong*[datay(n1)/200;datax(n1)/200]*sf(na1);
                ringsquares.x(:,end+1)=tx;
                ringsquares.y(:,end+1)=ty;
                ringsquares.index(:,end+1)=n1;
                ringsquares.names{end+1}=data.names2{n1};
                if isfield(data,'clusters')&&numel(data.clusters)>=n1, ringsquares.cluster(end+1)=data.clusters(n1); else ringsquares.cluster(end+1)=0; end
                tx=datax(n1)+rsquareshort*[datax(n1)/200;-datay(n1)/200]*sf(na1);
                ty=datay(n1)+rsquareshort*[datay(n1)/200;datax(n1)/200]*sf(na1);
                tz=2*zinc+zeros(size(rsquareshort,1),1);
                facealpha=1; facecolor=0*(.1+.8*data.plotconnoptions.BCOLOR)+1*round(1-data.plotconnoptions.BCOLOR);
                %if (~data.view&&SQUAREALLROIS)||(0&any(na1==data.source)||(n1<=size(z,1)&&(any(z(n1,:)>0)&&data.enablethr||n1<=numel(seedz)&&seedz(n1)>0&&data.mvpaenablethr))||any(z(:,n1)>0)&&data.enablethr),facealpha=1;
                %else facealpha=0;
                %end
                tpatch=struct('vertices',[tx,ty,tz],'faces',1:numel(tx),'facevertexcdata',repmat(facecolor,numel(tx),1),'facevertexalphadata',repmat(facealpha,numel(tx),1),'coords',{data.names2(n1)},'index',n1);
                h=cumpatch({data.hfig,3},tpatch,'edgecolor','none','parent',data.plotaxes);
                if ~data.view&&LINESTYLEMTX % diagonal reference
                    x=[datax([n1;n1]),datay([n1;n1])];
                    if LINESTYLEMTX==2||LINESTYLEMTX==4
                        x0=x/200;
                        x=max(0,datarank([na1;na1])-.5);
                        %x=(.5+angle(x*[1;1i])/pi/2);
                        fr1=.99*sqrt(max(0,x-.5)/numel(datarank));
                        fr2=.99*sqrt(max(0,x+.5)/numel(datarank));
                        %if data.issymmetric&datarank(na2)<datarank(na1), x=flipud(x); x0=flipud(x0); end
                        da=min(50/200,2*pi/length(data.displaytheserois))/2; %/1.125
                        tx1=(x0(1,1)-da*x0(1,2));ty1=(x0(1,2)+da*x0(1,1)); ntxty1=sqrt(max(eps,tx1.^2+ty1.^2)); tx1=tx1/ntxty1; ty1=ty1/ntxty1;
                        tx2=(x0(1,1)+da*x0(1,2)); ty2=(x0(1,2)-da*x0(1,1)); ntxty2=sqrt(max(eps,tx2.^2+ty2.^2)); tx2=tx2/ntxty2; ty2=ty2/ntxty2;
                        r1=fr1(2);
                        r2=fr2(2);
                        tsemic1=struct('vertices',200*[tx1*r1 ty1*r1; tx1*r2 ty1*r2; tx2*r2 ty2*r2; tx2*r1 ty2*r1],'faces',[1 2 3 4]);
                    else
                        x=max(0,datarank([na1;na1])-.5)/numel(datarank);
                        tx=.9*1.4142*(-100+200*x(1));
                        ty=.9*1.4142*(-100+200*x(2));
                        dx=.9*1.4142*200/numel(datarank)/2; %/1.125
                        dy=dx;
                        tsemic1=struct('vertices',[tx-dx ty-dy 110; tx-dx ty+dy 110; tx+dx ty+dy 110; tx+dx ty-dy 110],'faces',[1 2 3 4]);
                        %tsemic1.vertices(:,1)=-tsemic1.vertices(:,1); % flip ud
                        %tsemic1.vertices=tsemic1.vertices*[0.707106781186548 -0.707106781186547 0;0.707106781186547 0.707106781186548 0;0 0 1];
                        tsemic1.vertices(:,2)=-tsemic1.vertices(:,2); % flip ud
                    end
                    linetrans=1;
                    %tempc=data.plotconnoptions.BCOLOR;
                    %tsemic1.facevertexcdata=repmat(tempc,size(tsemic1.vertices,1),1);
                    %tsemic1.facevertexalphadata=repmat(linetrans,size(tsemic1.vertices,1),1);
                    %tsemic1.coords=[x];
                    %tsemic1.index=[n1;n1];
                    patch(tsemic1,'edgecolor','none','facecolor',1-data.plotconnoptions.BCOLOR,'parent',data.plotaxes); 
                end
            end
            idxtext=[idxtext;na1];
            if (0&any(na1==data.source)||(n1<=size(z,1)&&(any(z(n1,:)>0)&&data.enablethr||n1<=numel(seedz)&&seedz(n1)>0&&data.mvpaenablethr))||any(z(:,n1)>0)&&data.enablethr),%((n1<=size(z,1)||~data.displayreduced)&&any(z(:,n1)>0)),
                if data.mvpathrtype>1||(n1<=size(z,1)&&any(z(n1,:)>0))||any(z(:,n1)>0),%||data.z(n1)>0
                    idxtexthl=[idxtexthl;1];
                else idxtexthl=[idxtexthl;0];
                end
            else idxtexthl=[idxtexthl;0];
            end
            if n1<=N,
                for na2=1:length(data.displaytheserois), % connections
                    n2=data.displaytheserois(na2); %n2=[1:n1-1,n1+1:N2], 
                    %if n1~=n2&&((n1<n2&&z(n1,n2)>0)||(n1>n2&&z(n1,n2)>0&&n2<=N&&~(z(n2,n1)>0))), %(n2<=N||~data.displayreduced)&&(z(n1,n2)>0),
                    if z(n1,n2)>0,%&&n1~=n2, %(n2<=N||~data.displayreduced)&&(z(n1,n2)>0),
                        if data.issymmetric&&n2<n1&&z(n2,n1)>0 % avoids duplicated lines
                            %idxplotsadd1=find(data.list2(:,1)==n1&data.list2(:,2)==n2);
                            %idxplotsadd2=find(data.list2(:,1)==n2&data.list2(:,2)==n1);
                            %if ~isempty(idxplotsadd1)&&~isempty(idxplotsadd2), data.plotsadd2(idxplotsadd1)=data.plotsadd2(idxplotsadd2); end
                        else
                            x=[datax([n1;n2]),datay([n1;n2])];
                            dx=([data.x([n2;n1]),data.y([n2;n1])]-[data.x([n1;n2]),data.y([n1;n2])])/200;
                            if data.view>0 % obsolete
                                dx=dx./repmat(max(eps,sqrt(sum(abs(dx).^2,2))),[1,2]);
                                %h=patch(x(:,1)+dx(:,1).*5,x(:,2)+dx(:,2).*5,0*max(1,data.z([n1;n2]))*EPS+.5,'r-','linewidth',round(1+3*z(n1,n2)),'facecolor','none','edgecolor',[1,.5,.5]+0*[0,1,1]*z(n1,n2),'edgealpha',.5,'visible',data.visible);
                                linewidth=max(.1,abs(data.plotconnoptions.LINEWIDTH)*(4*(data.plotconnoptions.LINEWIDTH<0)*abs(z(n1,n2))+1*(data.plotconnoptions.LINEWIDTH>0))); % round(1+3*z(n1,n2))
                                h=plot3(x(:,1)+dx(:,1).*5,x(:,2)+dx(:,2).*5,0*max(1,data.z([n1;n2]))*EPS+.5,'r-','linewidth',linewidth,'color',[1,.5,.5]+0*[0,1,1]*z(n1,n2),'visible',data.visible);
                                %idxplotsadd2=find((data.list2(:,1)==n1&data.list2(:,2)==n2)|(data.list2(:,1)==n2&data.list2(:,2)==n1));
                                %if ~isempty(idxplotsadd2), data.plotsadd2(idxplotsadd2)={h}; end
                                %idxplotsadd2=find((data.list2(:,1)==n1&data.list2(:,2)==0)|(data.list2(:,1)==n2&data.list2(:,2)==0));
                                %for nt=1:numel(idxplotsadd2),data.plotsadd2{idxplotsadd2(nt)}=cat(2,data.plotsadd2{idxplotsadd2(nt)},h); end
                                %if data.h(n1,n2)<0, set(h,'color',[.5,.5,1]+0*[1,1,0]*z(n1,n2)); end
                                if data.plotconnoptions.LCOLOR==1
                                    if data.h(n1,n2)<0, set(h,'color',[.5,.5,1]); end
                                    %set(h,'edgecolor',cmap(ceil(size(cmap,1)/2)+round(floor(size(cmap,1)/2)*data.h(n1,n2)/J),:));
                                elseif data.plotconnoptions.LCOLOR==2
                                    tempc=(hsv2rgb([(1+angle(x(1))/pi)/2,1,1])+hsv2rgb([(1+angle(x(2))/pi)/2,1,1]))/2;
                                    set(h,'color',tempc);
                                else
                                    %tempc=cmap(max(1,min(size(cmap,1), round(size(cmap,1)/2+sign(data.h(n1,n2))*abs(data.h(n1,n2))^data.plotconnoptions.LCOLORSCALE*size(cmap,1)/2/Kscale))),:);
                                    tempc=cmap(max(1,min(size(cmap,1), round(size(cmap,1)/2+sign(data.h(n1,n2))*abs(data.F(n1,n2))^data.plotconnoptions.LCOLORSCALE*size(cmap,1)/2/Kscale))),:);
                                    set(h,'color',tempc);
                                end
                                %%set(h,'buttondownfcn',@conn_displayroi_menubuttondownfcn,'userdata',data);
                                %rangle=rtriang*exp(j*angle(diff(datax([n1;n2])+j*datay([n1;n2]),1,1)))*(1/2+1/2*z(n1,n2));
                                %h=patch(x(2,1)+dx(2,1)*7+5*real(rangle),x(2,2)+dx(2,2)*7+5*imag(rangle),0*max(1,data.z(n2))*EPS+.5+zeros(size(rangle)),'w'); set(h,'edgecolor','none','facecolor',[1,.5,.5]+0*[0,1,1]*z(n1,n2),'visible',data.visible);
                                %if data.h(n1,n2)<0, set(h,'facecolor',[.5,.5,1]+0*[1,1,0]*z(n1,n2)); end
                            elseif LINESTYLEMTX % matrix
                                if LINESTYLEMTX==2||LINESTYLEMTX==4 % angle
                                    x0=x/200;
                                    x=max(0,datarank([na1;na2])-.5);
                                    if data.issymmetric&x(2)>x(1), x=flipud(x); x0=flipud(x0); end
                                    %x0=x0./repmat(sqrt(sum(x0.^2,2)),1,2);
                                    %x=(.5+angle(x*[1;1i])/pi/2);
                                    fr1=.99*sqrt(max(0,x-.5)/numel(datarank));
                                    fr2=.99*sqrt(max(0,x+.5)/numel(datarank));
                                    %fr1=.99*log(1+max(0,x-.5)/numel(datarank))/log(2); 
                                    %fr2=.99*log(1+max(0,x+.5)/numel(datarank))/log(2); 
                                    %if data.issymmetric&datarank(na2)<datarank(na1), x=flipud(x); x0=flipud(x0); end
                                    da=min(50/200,2*pi/length(data.displaytheserois))/2; %/1.125
                                    tx1=(x0(1,1)-da*x0(1,2));ty1=(x0(1,2)+da*x0(1,1)); ntxty1=sqrt(max(eps,tx1.^2+ty1.^2)); tx1=tx1/ntxty1; ty1=ty1/ntxty1; 
                                    tx2=(x0(1,1)+da*x0(1,2)); ty2=(x0(1,2)-da*x0(1,1)); ntxty2=sqrt(max(eps,tx2.^2+ty2.^2)); tx2=tx2/ntxty2; ty2=ty2/ntxty2; 
                                    r1=fr1(2);
                                    r2=fr2(2);
                                    if data.issymmetric&&LINESTYLEMTX==2, 
                                        v1=[tx1*r1 ty1*r1; tx1*r2 ty1*r2; tx2*r2 ty2*r2; tx2*r1 ty2*r1];
                                        tx1=(x0(2,1)-da*x0(2,2));ty1=(x0(2,2)+da*x0(2,1)); ntxty1=sqrt(max(eps,tx1.^2+ty1.^2)); tx1=tx1/ntxty1; ty1=ty1/ntxty1;
                                        tx2=(x0(2,1)+da*x0(2,2)); ty2=(x0(2,2)-da*x0(2,1)); ntxty2=sqrt(max(eps,tx2.^2+ty2.^2)); tx2=tx2/ntxty2; ty2=ty2/ntxty2;
                                        r1=fr1(1);
                                        r2=fr2(1);
                                        tsemic1=struct('vertices',200*[v1; tx1*r1 ty1*r1; tx1*r2 ty1*r2; tx2*r2 ty2*r2; tx2*r1 ty2*r1],'faces',[1 2 3 4; 5 6 7 8]);
                                    else
                                        tsemic1=struct('vertices',200*[tx1*r1 ty1*r1; tx1*r2 ty1*r2; tx2*r2 ty2*r2; tx2*r1 ty2*r1],'faces',[1 2 3 4]);
                                    end
                                else % square
                                    x=max(0,datarank([na1;na2])-.5)/numel(datarank);
                                    if data.issymmetric&x(2)>x(1), x=flipud(x); end
                                    tx=.9*1.4142*(-100+200*x(1));
                                    ty=.9*1.4142*(-100+200*x(2));
                                    dx=.9*1.4142*200/numel(datarank)/2; %/1.125;
                                    dy=dx;
                                    if data.issymmetric&&LINESTYLEMTX==1, tsemic1=struct('vertices',[tx-dx ty-dy 100; tx-dx ty+dy 100; tx+dx ty+dy 100; tx+dx ty-dy 100;ty-dy tx-dx 100; ty-dy tx+dx 100; ty+dy tx+dx 100; ty+dy tx-dx 100],'faces',[1 2 3 4;5 6 7 8]);
                                    else tsemic1=struct('vertices',[tx-dx ty-dy 100; tx-dx ty+dy 100; tx+dx ty+dy 100; tx+dx ty-dy 100],'faces',[1 2 3 4]);
                                    end
                                    %tsemic1.vertices(:,1)=-tsemic1.vertices(:,1); % flip ud
                                    %tsemic1.vertices=tsemic1.vertices*[0.707106781186548 -0.707106781186547 0;0.707106781186547 0.707106781186548 0;0 0 1];
                                    tsemic1.vertices(:,2)=-tsemic1.vertices(:,2); % flip ud
                                end
                                linetrans=min(1,max(.001,abs(data.plotconnoptions.LTRANS)*((data.plotconnoptions.LTRANS<0)*abs(z(n1,n2))+(data.plotconnoptions.LTRANS>0))));
                                tempc=[1 .25 .25];
                                if data.plotconnoptions.LCOLOR==1
                                    if data.h(n1,n2)<0, tempc=[.25,.25,1]; end
                                elseif data.plotconnoptions.LCOLOR==2
                                    tempc=(hsv2rgb([x(1),1,1])+hsv2rgb([x(2),1,1]))/2;
                                else
                                    %tempc=cmap(max(1,min(size(cmap,1), round(size(cmap,1)/2+sign(data.h(n1,n2))*abs(data.h(n1,n2))^data.plotconnoptions.LCOLORSCALE*size(cmap,1)/2/Kscale))),:);
                                    tempc=cmap(max(1,min(size(cmap,1), round(size(cmap,1)/2+sign(data.h(n1,n2))*abs(data.CM_z0(n1,n2)/emaxz)^data.plotconnoptions.LCOLORSCALE*size(cmap,1)/2))),:);
                                end
                                tsemic1.facevertexcdata=repmat(tempc,size(tsemic1.vertices,1),1);
                                tsemic1.facevertexalphadata=repmat(linetrans,size(tsemic1.vertices,1),1);
                                if data.issymmetric&&(LINESTYLEMTX==1||LINESTYLEMTX==2), 
                                    tsemic1.coords=[x,flipud(x)];
                                    tsemic1.index=[n1 n2;n2 n1];
                                else
                                    tsemic1.coords=[x];
                                    tsemic1.index=[n1;n2];
                                end
                                h=cumpatch({data.hfig,1},tsemic1,'edgecolor','none','parent',data.plotaxes); %max(1,data.plotconnoptions.LINEWIDTH*(1+0*(n1<=numel(seedz)&seedz(min(numel(seedz),n1))>0&n2<=numel(seedz)&seedz(min(numel(seedz),n2))>0)))));
                            else % lines
                                x=x*[1;1i];
                                dx=dx(1,:)*[1;1i];
                                if data.displaybrains==2&&isfield(data,'bclusters')&&length(data.bclusters)>=max(n1,n2)&&data.bclusters(n1)==data.bclusters(n2), lcurve=0;
                                else lcurve=data.plotconnoptions.LCURVE;
                                end
                                %lcurve=lcurve*sign(data.h(n1,n2));
                                if 0,
                                    xtracurve=(1+2*exp(-abs(dx)*20));
                                    temp=(real(semic1)+xtracurve*lcurve*1i*imag(semic1)*((1-(.98+.01*rand)*abs(dx(1)/2).^2)));
                                    xt1=(x(1)+x(2))/2 + (x(2)-x(1))/2 * temp;
                                    xt2=(x(1)+x(2))/2 - (x(2)-x(1))/2 * temp(end:-1:1);
                                elseif 1
                                    xtracurve=exp(-abs(x(1)-x(2))/200)*semic1factor1;
                                    p1=.5-real(semic1)/2;
                                    xt1=x(1).*p1.^2 + x(2).*(1-p1).^2+(x(1)+x(2))/2*(2*p1.*(1-p1)).*xtracurve;
                                    xt2=xt1;
                                elseif 0 %lcurve<0,
                                    xt1=(x(1)+(1+real(semic1))/2*(x(2)-x(1)))./(1-lcurve/4*(1-abs(real(semic1)).^2));
                                    xt2=xt1;
                                else
                                    tlcurve=1+lcurve*abs(dx(1)/2)*.5;
                                    w0=(.5+real(semic1)/2).^tlcurve;
                                    w1=x(1)*w0;
                                    w2=x(2)*fliplr(w0);
                                    xt1=w1+w2;
                                    xt2=xt1;
                                end
                                if lcurve*min(abs(xt1))<lcurve*min(abs(xt2)), xt=xt1; else xt=xt2; end
                                %h=patch(1*[real(xt),fliplr(real(xt))],1*[imag(xt),fliplr(imag(xt))],[linspace(data.xyz2(n1,3),data.xyz2(n2,3),numel(xt)),linspace(data.xyz2(n2,3),data.xyz2(n1,3),numel(xt))],'r','edgecolor',[1,.25,.25],'facecolor','none','edgealpha',.05+.95*z(n1,n2).^EMPH,'linewidth',data.plotconnoptions.LINEWIDTH*(1+0*(n1<=numel(seedz)&seedz(min(numel(seedz),n1))>0&n2<=numel(seedz)&seedz(min(numel(seedz),n2))>0)));
                                linewidth=max(.1,abs(data.plotconnoptions.LINEWIDTH)*(4*(data.plotconnoptions.LINEWIDTH<0)*abs(z(n1,n2))+1*(data.plotconnoptions.LINEWIDTH>0)));
                                linetrans=min(1,max(.001,abs(data.plotconnoptions.LTRANS)*((data.plotconnoptions.LTRANS<0)*abs(z(n1,n2))+(data.plotconnoptions.LTRANS>0))));
                                %h=patch(1*[real(xt),fliplr(real(xt))],1*[imag(xt),fliplr(imag(xt))],(xy2(n1,3)+xy2(n2,3))/2+(xy2(n2,3)-xy2(n1,3))/2*real([semic1,fliplr(semic1)])+z(n1,n2)*EPS+.5+1000*abs(z(n1,n2))*[imag(semic1),fliplr(imag(semic1))],'r','edgecolor',[1,.25,.25],'facecolor','none','edgealpha',linetrans,'linewidth',linewidth); %max(1,data.plotconnoptions.LINEWIDTH*(1+0*(n1<=numel(seedz)&seedz(min(numel(seedz),n1))>0&n2<=numel(seedz)&seedz(min(numel(seedz),n2))>0)))));
                                tempc=[1 .25 .25];
                                if data.plotconnoptions.LCOLOR==1
                                    if data.h(n1,n2)<0, tempc=[.25,.25,1]; end
                                elseif data.plotconnoptions.LCOLOR==2
                                    tempc=(hsv2rgb([(1+angle(x(1))/pi)/2,1,1])+hsv2rgb([(1+angle(x(2))/pi)/2,1,1]))/2;
                                else
                                    %tempc=cmap(max(1,min(size(cmap,1), round(size(cmap,1)/2+sign(data.h(n1,n2))*abs(data.F(n1,n2))^data.plotconnoptions.LCOLORSCALE*size(cmap,1)/2/Kscale))),:);
                                    tempc=cmap(max(1,min(size(cmap,1), round(size(cmap,1)/2+sign(data.h(n1,n2))*abs(data.CM_z0(n1,n2)/emaxz)^data.plotconnoptions.LCOLORSCALE*size(cmap,1)/2))),:);
                                end
                                %dwidth=1i*(xt(end)-xt(1));dwidth=1*linewidth*dwidth/abs(dwidth);
                                dwidth=1i*[xt(2)-xt(1) xt(3:end)-xt(1:end-2) xt(end)-xt(end-1)]; 
                                dwidth=linewidth*exp(1i*angle(dwidth)).*dwidthfactor1;
                                %dwidth=linewidth*exp(1i*angle(dwidth)).*[abs(rsquare(1,2))*[1 .5] ones(1,numel(xt)-4) [.5 1]*abs(rsquare(1,2))]/2;
                                %dwidth=linewidth*exp(1i*angle(dwidth)).*ones(1,numel(xt))/2;
                                %dwidth=linewidth*exp(1i*angle(dwidth)).*[0 .5 ones(1,numel(xt)-4) .5 0]/2;
                                tsemic1=rsemic1;
                                tsemic1.vertices=[[real(xt+dwidth),real(xt(end:-1:1)-dwidth(end:-1:1))]' ...
                                                  [imag(xt+dwidth),imag(xt(end:-1:1)-dwidth(end:-1:1))]' ...
                                                  (+100*abs(z(n1,n2))*semic1factor2)'];
                                                  %[0*((dataz(n1)+dataz(n2))/2+(dataz(n2)-dataz(n1))/2*real(ssemic1))+0*z(n1,n2)*EPS+0*.5+100*abs(z(n1,n2))*(imag(ssemic1)*2).^.05]'];
                                tsemic1.facevertexcdata=repmat(tempc,size(tsemic1.vertices,1),1);
                                tsemic1.facevertexalphadata=repmat(linetrans,size(tsemic1.vertices,1),1);
                                if ~isempty(data.xy2_clusters), tsemic1.coords=[x;5*data.xy2_clusters([n1,n2],:)*[1;1i];1e6*sign(data.h(n1,n2))]; % note: for sign-specific smoothing
                                %if ~isempty(data.xy2_clusters), tsemic1.coords=[x;1e6*data.clusters([n1,n2]);1e6*sign(data.h(n1,n2))]; % note: for sign-specific smoothing
                                else tsemic1.coords=[x;1000*sign(data.h(n1,n2))]; % note: for sign-specific smoothing
                                end
                                tsemic1.index=[n1;n2];
                                h=cumpatch({data.hfig,1},tsemic1,'edgecolor','none','parent',data.plotaxes); %max(1,data.plotconnoptions.LINEWIDTH*(1+0*(n1<=numel(seedz)&seedz(min(numel(seedz),n1))>0&n2<=numel(seedz)&seedz(min(numel(seedz),n2))>0)))));
                                %h=patch(tsemic1,'edgecolor','none','facecolor',tempc,'facealpha',linetrans); %max(1,data.plotconnoptions.LINEWIDTH*(1+0*(n1<=numel(seedz)&seedz(min(numel(seedz),n1))>0&n2<=numel(seedz)&seedz(min(numel(seedz),n2))>0)))));
%                                 h=patch([real(xt),fliplr(real(xt))],...
%                                     [imag(xt),fliplr(imag(xt))],...
%                                     (dataz(n1)+dataz(n2))/2+(dataz(n2)-dataz(n1))/2*real([semic1,fliplr(semic1)])+0*z(n1,n2)*EPS+0*.5+1000*abs(z(n1,n2))*[imag(semic1),fliplr(imag(semic1))],...
%                                     'r','edgecolor',tempc,'facecolor','none','edgealpha',linetrans,'linewidth',linewidth); %max(1,data.plotconnoptions.LINEWIDTH*(1+0*(n1<=numel(seedz)&seedz(min(numel(seedz),n1))>0&n2<=numel(seedz)&seedz(min(numel(seedz),n2))>0)))));
                                %idxplotsadd2=find((data.list2(:,1)==n1&data.list2(:,2)==n2)|(data.list2(:,1)==n2&data.list2(:,2)==n1));
                                %if ~isempty(idxplotsadd2), data.plotsadd2(idxplotsadd2)={h}; end
                                %idxplotsadd2=find((data.list2(:,1)==n1&data.list2(:,2)==0)|(data.list2(:,1)==n2&data.list2(:,2)==0));
                                %for nt=1:numel(idxplotsadd2),data.plotsadd2{idxplotsadd2(nt)}=cat(2,data.plotsadd2{idxplotsadd2(nt)},h); end
                                %if ~isempty(data.clusters), wcolor=data.clusters([n1,n2])/maxdataclusters; wcolor(~wcolor)=[]; if any(wcolor), set(h,'edgecolor',cmap(ceil(wcolor(ceil(numel(wcolor)*rand))*size(cmap,1)),:)); end; end
                                %if ~isempty(data.clusters), wcolor=data.clusters([n1,n2])>0; if any(wcolor), set(h,'edgecolor',wcolor'/sum(wcolor)*cmap(ceil(data.clusters([n1,n2])*size(cmap,1)/maxdataclusters),:)); end; end
                                %set(h,'buttondownfcn',@conn_displayroi_menubuttondownfcn,'userdata',data);
                            end
                        end
                    end
                end
            end
        end

        bbox=[];
        if ~data.view&&(~isempty(data.clusters)||(isfield(data,'displaybrains')&&data.displaybrains)) % brain displays
            if isempty(data.clusters), NPLOTS=data.plotconnoptions.NPLOTS;
            else NPLOTS=max(data.clusters); 
            end
            if SINGLEBRAIN, offset=OFFSET+.5;
            else offset=OFFSET+data.plotconnoptions.DOFFSET*(data.displaylabels>0);
            end
            
            xy=data.x+1i*data.y;
            mr=mean(data.plotconnoptions.rende.vertices,1);
            colors=[.3+.4*get(data.hfig,'color');1-get(data.hfig,'color')];
            rende=[];
            if SINGLEBRAIN&&data.displaybrains % single-brain
                nprojection=data.plotconnoptions.nprojection;
                if nprojection<=0, nprojection=3; end
                p=data.plotconnoptions.Projections{nprojection}'*(1+offset)*data.plotconnoptions.BSCALE;
                rende=data.plotconnoptions.rende;
                rende.vertices=1*[1.05*detrend(rende.vertices,'constant'),ones(size(rende.vertices,1),1)]*[p;-400,+200*(1+max(0,data.plotconnoptions.DOFFSET-.35)),0];
                w=(-1+2*(mean(data.plotconnoptions.BCOLOR)>.5))*rende.vertices(:,3); w=sqrt(max(0,0+1*(w-min(w(:)))/max(eps,max(w(:))-min(w(:)))));
                tmin1=min(rende.vertices,[],1);tmax1=max(rende.vertices,[],1);
                h0=patch(struct('vertices',repmat(tmin1-.1*(tmax1-tmin1),[4,1])+1.2*[0;1;1;0]*[tmax1(1)-tmin1(1) 0 0]+1.2*[0;0;1;1]*[0 tmax1(2)-tmin1(2) 0],'faces',[1,2,3,4]),'edgecolor',1-data.plotconnoptions.BCOLOR,'facecolor','none','tag','plot_brainbackground','parent',data.plotaxes);
                h1=text(tmin1(1)-.25/2*(tmax1(1)-tmin1(1)+tmax1(2)-tmin1(2)),(tmin1(2)+tmax1(2))/2,tmin1(3),data.plotconnoptions.Projections_axes{nprojection}{2},'color',1-data.plotconnoptions.BCOLOR,'fontsize',data.plotconnoptions.FONTSIZE(1),'horizontalalignment','center','parent',data.plotaxes);
                h2=text((tmin1(1)+tmax1(1))/2,tmin1(2)-.25/2*(tmax1(1)-tmin1(1)+tmax1(2)-tmin1(2)),tmin1(3),data.plotconnoptions.Projections_axes{nprojection}{1},'color',1-data.plotconnoptions.BCOLOR,'fontsize',data.plotconnoptions.FONTSIZE(1),'horizontalalignment','center','parent',data.plotaxes);
                h=patch(rende,'edgecolor','none','facevertexcdata',(1-w)*(1-data.plotconnoptions.BCOLOR)+w*data.plotconnoptions.BCOLOR,'facealpha',data.plotconnoptions.BTRANS,'facecolor','inter','parent',data.plotaxes);
                %set([h1 h2],'tag','conn_displayroi_brain','buttondownfcn','conn_displayroi(''display.viewcycle'');','userdata',data.hfig);
                if isempty(bbox), bbox=[tmin1-(tmax1-tmin1).*[1 .1 0]; tmax1+(tmax1-tmin1).*[1 .1 0]];
                else bbox=[min(bbox(1,:),tmin1-(tmax1-tmin1)*[1 .1 0]); max(bbox(2,:),tmax1+(tmax1-tmin1).*[1 .1 0])];
                end
            end
            for n1=1:NPLOTS,
                if isempty(data.clusters), idx=find(abs(angle(xy.*exp(-1i*2*pi/NPLOTS*n1)))<2*pi/NPLOTS/2);
                    %idx=find(max(1,min(NPLOTS,ceil((angle(xy)+pi)/2/pi*NPLOTS)))==n1);
                else idx=find(data.clusters==n1); 
                end
                idx1=find(markthese(idx)~=0);
                if isempty(idx1)
                    if isfield(data,'displaybrains')&&isfield(data,'names_clusters')&&~isempty(data.names_clusters)&&~isempty(idx) % text labels
                        mx=mean(xy(idx));
                        if isempty(data.clusters), a0=2*pi/NPLOTS*n1;
                        else a0=angle(mx);
                        end
                        if data.displaybrains<2, efact=(1.11+data.plotconnoptions.DOFFSET);
                        else efact=(1+offset)*(1.11+data.plotconnoptions.DOFFSET+.75*data.plotconnoptions.BSCALE);
                        end
                        tx=efact*200*cos(a0);
                        ty=efact*200*sin(a0);
                        h=conn_menu_plottext(tx,ty,data.names_clusters{n1},'horizontalalignment','center','fontsize',data.plotconnoptions.FONTSIZE(end),'color',.05+.9*data.plotconnoptions.BCOLOR,'parent',data.plotaxes);
                        if isempty(bbox), bbox=1.15*[tx ty 0; tx ty 0];
                        else bbox=[min(bbox(1,:),1.15*[tx ty 0]); max(bbox(2,:),1.15*[tx ty 0])];
                        end
                        %set(h,'tag','textstring'); 
                    end
                else
                    mx=mean(xy(idx));
                    px=cos(angle(mx));
                    py=sin(angle(mx));
                    if isempty(data.clusters), a0=2*pi/NPLOTS*n1;
                    else a0=angle(mx);
                    end
                    if isfield(data,'displaybrains')&&data.displaybrains<2 % cluster reference markers
                        joffset=min(50,1.0*2*pi*200/length(data.displaytheserois)/1.125)/200/2;
%                         if ~isempty(data.names_clusters)&&isfield(data,'displaylabels')&&data.displaylabels, 
% %                             koffset=.11+data.plotconnoptions.DOFFSET;
% %                             ltemp=[(1+koffset-.02)*195,(1+koffset)*195*ones(1,62),(1+koffset-.02)*195,(1+.11+.02)*200,(1+.11)*200*ones(1,62),(1+.11+.02)*200].*exp(1i*(a0+[linspace(-joffset+min(angle(xy(idx).*exp(-1i*a0))),joffset+max(angle(xy(idx).*exp(-1i*a0))),64),linspace(joffset+max(angle(xy(idx).*exp(-1i*a0))),-joffset+min(angle(xy(idx).*exp(-1i*a0))),64)]));
% %                             patch(real(ltemp),imag(ltemp),0*1000+zeros(size(ltemp)),'k','facecolor',.05+.9*data.plotconnoptions.BCOLOR,'edgecolor','none','parent',data.plotaxes);
%                             koffset=.11+.95*data.plotconnoptions.DOFFSET;
%                             plot((1+koffset+.0*rem(n1,2))*200*[.98,ones(1,62),.98].*exp(1i*(a0+linspace(-joffset+min(angle(xy(idx).*exp(-1i*a0))),joffset+max(angle(xy(idx).*exp(-1i*a0))),64))),'k-','color',.7-.4*get(data.hfig,'color'),'linewidth',2,'parent',data.plotaxes);
%                         end
                        koffset=.11;
                        %if isfield(data,'displaybrains')&&data.displaybrains, koffset=offset; else koffset=.11; end
                        plot((1+koffset+.0*rem(n1,2))*200*[.98,ones(1,62),.98].*exp(1i*(a0+linspace(-joffset+min(angle(xy(idx).*exp(-1i*a0))),joffset+max(angle(xy(idx).*exp(-1i*a0))),64))),'k-','color',.7-.4*get(data.hfig,'color'),'linewidth',2,'parent',data.plotaxes);
                        ltemp=[(1+koffset)*195,(1+koffset)*195*ones(1,62),(1+koffset)*195,(1+.11)*180,(1+.11)*180*ones(1,62),(1+.11)*180].*exp(1i*(a0+[linspace(-joffset+min(angle(xy(idx).*exp(-1i*a0))),joffset+max(angle(xy(idx).*exp(-1i*a0))),64),linspace(joffset+max(angle(xy(idx).*exp(-1i*a0))),-joffset+min(angle(xy(idx).*exp(-1i*a0))),64)]));
                        patch(real(ltemp),imag(ltemp),0*1000+zeros(size(ltemp)),'k','facecolor',.1+.8*data.plotconnoptions.BCOLOR,'edgecolor','none','linewidth',2,'parent',data.plotaxes);
                        patch(real(ltemp),imag(ltemp),1000+zeros(size(ltemp)),'k','facecolor','none','edgecolor',round(data.plotconnoptions.BCOLOR),'linewidth',2,'parent',data.plotaxes);
                    end
                    if isfield(data,'displaybrains')&&isfield(data,'names_clusters')&&~isempty(data.names_clusters)  % text labels
                        if data.displaybrains<2, efact=(1.11+data.plotconnoptions.DOFFSET);
                        else efact=(1+offset)*(1.11+data.plotconnoptions.DOFFSET+.75*data.plotconnoptions.BSCALE);
                        end
                        tx=efact*200*cos(a0);
                        ty=efact*200*sin(a0);
                        h=conn_menu_plottext(tx,ty,data.names_clusters{n1},'horizontalalignment','center','fontsize',data.plotconnoptions.FONTSIZE(end),'color',1-data.plotconnoptions.BCOLOR,'parent',data.plotaxes);
                        if isempty(bbox), bbox=1.15*[tx ty 0; tx ty 0];
                        else bbox=[min(bbox(1,:),1.15*[tx ty 0]); max(bbox(2,:),1.15*[tx ty 0])];
                        end
                        %set(h,'tag','textstring');
                    end
                    if isfield(data,'displaybrains')&&data.displaybrains
                        if ~SINGLEBRAIN&&~data.plotconnoptions.nprojection
                            cprojection=cellfun(@(x)sum(std(data.xyz2(idx,:)*x(1:2,:)',1,1).^2),data.plotconnoptions.Projections);
                            cprojection(2)=nan;
                            [nill,nprojection]=max(cprojection);
                        else
                            nprojection=data.plotconnoptions.nprojection;
                        end
                        if nprojection<=0, nprojection=3; end
                        p=data.plotconnoptions.Projections{nprojection}'*(1+offset)*data.plotconnoptions.BSCALE;
                        if SINGLEBRAIN
                            xy2=1*[data.xyz2(idx,:)-mr(ones(numel(idx),1),:),ones(numel(idx),1)]*[p;-400,+200*(1+max(0,data.plotconnoptions.DOFFSET-.35)),0]; % roi positions in brain
                        else % multiple-brains
                            if nprojection==1&&mean(data.xyz2(idx,1)>0)>.5, p(:,1)=-p(:,1); end
                            rende=data.plotconnoptions.rende;
                            rende.vertices=[1.05*detrend(rende.vertices,'constant'),ones(size(rende.vertices,1),1)]*[p;(1+offset)*(1+.75*data.plotconnoptions.BSCALE)*200*px,(1+offset)*(1+.75*data.plotconnoptions.BSCALE)*200*py,0];
                            %h=patch(rende,'edgecolor','none','facecolor',.8-.6*get(data.hfig,'color'),'facealpha',data.plotconnoptions.BTRANS);
                            w=(-1+2*(mean(data.plotconnoptions.BCOLOR)>.5))*rende.vertices(:,3); w=sqrt(max(0,0+1*(w-min(w(:)))/max(eps,max(w(:))-min(w(:)))));
                            tmin1=min(rende.vertices,[],1);tmax1=max(rende.vertices,[],1);
                            %h0=patch(struct('vertices',repmat(tmin1-.1*(tmax1-tmin1),[4,1])+1.2*[0;1;1;0]*[tmax1(1)-tmin1(1) 0 0]+1.2*[0;0;1;1]*[0 tmax1(2)-tmin1(2) 0],'faces',[1,2,3,4]),'edgecolor',1-data.plotconnoptions.BCOLOR,'facecolor',data.plotconnoptions.BCOLOR,'parent',data.plotaxes);
                            h=patch(rende,'edgecolor','none','facevertexcdata',(1-w)*(1-data.plotconnoptions.BCOLOR)+w*data.plotconnoptions.BCOLOR,'facealpha',data.plotconnoptions.BTRANS,'facecolor','inter','parent',data.plotaxes);
                            %set(h,'tag','conn_displayroi_brain','buttondownfcn','conn_displayroi(''display.viewcycle'');','userdata',data.hfig);
                            if isempty(bbox), bbox=[tmin1-.11*(tmax1-tmin1); tmax1+.11*(tmax1-tmin1)];
                            else bbox=[min(bbox(1,:),tmin1-.11*(tmax1-tmin1)); max(bbox(2,:),tmax1+.11*(tmax1-tmin1))];
                            end
                            %xy2=[datax(idx) datay(idx) dataz(idx)]; % roi positions in connections (ring/brain)
                            xy2=[data.xyz2(idx,:)-mr(ones(numel(idx),1),:),ones(numel(idx),1)]*[p;(1+offset)*(1+.75*data.plotconnoptions.BSCALE)*200*px,(1+offset)*(1+.75*data.plotconnoptions.BSCALE)*200*py,0]; % roi positions in brain
                            
                        end
                        for n2=idx1(:)' % rois in brains
                            zinc=100; %100;
                            facecolor=round(1-data.plotconnoptions.BCOLOR); %cmap(round(1+(size(cmap,1)-1)*(1+k)/2),:)
                            tx=xy2(n2,1)+data.plotconnoptions.RSCALE*rcircle(:,1)*(1+offset)*data.plotconnoptions.BSCALE*(2+0*markthese(idx(n2)));
                            ty=xy2(n2,2)+data.plotconnoptions.RSCALE*rcircle(:,2)*(1+offset)*data.plotconnoptions.BSCALE*(2+0*markthese(idx(n2)));
                            tz=xy2(n2,3)+zeros(size(rcircle,1),1)+zinc*markthese(idx(n2));
                            ringsquares.x(:,end+1)=tx(round((1:size(ringsquares.x,1))*(numel(tx)+1)/(size(ringsquares.x,1)+1)));
                            ringsquares.y(:,end+1)=ty(round((1:size(ringsquares.y,1))*(numel(ty)+1)/(size(ringsquares.y,1)+1)));
                            ringsquares.index(:,end+1)=idx(n2);
                            ringsquares.names{end+1}=data.names2{idx(n2)};
                            ringsquares.cluster(end+1)=n1;
                            %h=patch(tx,ty,tz,...
                            %    'k','facecolor',facecolor,'edgecolor','none','parent',data.plotaxes); %1-get(data.hfig,'color'));
                            %set(h,'tag','conn_displayroi_roi');
                            facealpha=1; 
                            tpatch=struct('vertices',[tx,ty,tz],'faces',1:numel(tx),'facevertexcdata',repmat(facecolor,numel(tx),1),'facevertexalphadata',repmat(facealpha,numel(tx),1),'coords',{data.names2(idx(n2))},'index',idx(n2));
                            ttxt={['x,y,z = (',num2str(data.xyz2(idx(n2),1),'%1.0f'),',',num2str(data.xyz2(idx(n2),2),'%1.0f'),',',num2str(data.xyz2(idx(n2),3),'%1.0f'),') mm'],...
                               data.names2{idx(n2)}};
                            h=cumpatch({data.hfig,2},tpatch,'edgecolor','none','parent',data.plotaxes,'tag','conn_displayroi_roi','buttondownfcn',@conn_displayroi_menubuttondownfcn,'userdata',[{data.hfig} ttxt],'interruptible','off');
%                             k=KK(idx(n2));
%                             if isnan(k), set(h,'facecolor','none');
%                             %else set(h,'facecolor',cmap(round(1+(size(cmap,1)-1)*(1+k)/2),:));
%                             end
%                             set(h,'buttondownfcn',@conn_displayroi_menubuttondownfcn,'userdata',{data.hfig, idx(n2)},'interruptible','off');
                            %ttxt={['x,y,z = (',num2str(data.xyz2(idx(n2),1),'%1.0f'),',',num2str(data.xyz2(idx(n2),2),'%1.0f'),',',num2str(data.xyz2(idx(n2),3),'%1.0f'),') mm'],...
                            %    data.names2{idx(n2)}};
                            %set(h,'buttondownfcn',@conn_displayroi_menubuttondownfcn,'userdata',[{data.hfig} ttxt],'interruptible','off');
                            %data.plotsadd3{idx(n2)}=[data.plotsadd3{idx(n2)},h];
                        end
                        %if data.displaybrains<2
                        %    p1=(1+koffset+.0*rem(n1,2))*200*[0,10+ones(1,62),0].*exp(1i*(a0+[-joffset+min(angle(xy(idx).*exp(-1i*a0))),linspace(-joffset+min(angle(xy(idx).*exp(-1i*a0))),joffset+max(angle(xy(idx).*exp(-1i*a0))),62),joffset+max(angle(xy(idx).*exp(-1i*a0)))]));
                        %    h=patch(real(p1),imag(p1),-ones(size(p1)),'w','facecolor',.1+.8*data.plotconnoptions.BCOLOR,'edgecolor','none','parent',data.plotaxes);
                        %end
                    end
                end
            end
        end
                
        if ~LINESTYLEMTX&&data.displaybrains<2, cumpatch({data.hfig,1},'smooth&update',data.plotconnoptions.LBUNDL);
        else cumpatch({data.hfig,1},'update');
        end
        h1=cumpatch({data.hfig,2},'update');
        h2=cumpatch({data.hfig,3},'update');
        set([h1(:);h2(:)],'buttondownfcn',@conn_displayroi_menubuttondownfcn,'userdata',data.hfig,'interruptible','off');
        if ~isempty(markthese), 
            cumpatch({data.hfig,2},'mask',[1:numel(markthese);double(markthese'>1)]); 
            cumpatch({data.hfig,3},'mask',[1:numel(markthese);double(markthese'>1)]); 
        end
        data.plot_idxselected=idxtext(idxtexthl>0);
        data.plot_K=K;
        data.plot_cmap=cmap;
        data.plot_z=z;
        %%cumpatch({data.hfig,3},'mask',[1:N; zeros(1,N)]);

        if isfield(data,'displaylabels')&&data.displaylabels % text labels
            temp=data.names2reduced; for n1=1:numel(temp), if numel(temp{n1})>50, temp{n1}=temp{n1}(1:50); end; end
            tmvpaz=[seedz;0];
            if isfield(data,'names_clusters')&&~isempty(data.names_clusters)&&data.plotconnoptions.DOFFSET<=.01, efact=1.05; else efact=1; end
            h=text(efact*1.025*datax(data.displaytheserois(idxtext)),efact*1.025*datay(data.displaytheserois(idxtext)),...
                max(0,5*data.z(data.displaytheserois(idxtext)))*EPS+202,{temp{data.displaytheserois(idxtext)}},'parent',data.plotaxes);
                %(tmvpaz(min(numel(tmvpaz),data.displaytheserois(idxtext)))>0)+1002,{temp{data.displaytheserois(idxtext)}});
                %(seedz(data.displaytheserois(idxtext(data.displaytheserois(idxtext)<=size(z,1))))>0)+202,{temp{data.displaytheserois(idxtext)}});
            set(h,'tag','textstringpartial','clipping','off');
            if data.view>0
                set(h,'fontsize',data.plotconnoptions.FONTSIZE(1),'color','k','horizontalalignment','center','interpreter','none','fontweight','normal','backgroundcolor','none');
            else
                fontsize=data.plotconnoptions.FONTSIZE(1);
                %if isfield(data,'displaybrains')&&data.displaybrains, fontsize=max(4,fontsize-3); end
                set(h,'fontsize',fontsize,'color',.1*.5+.9*data.plotconnoptions.BCOLOR,'horizontalalignment','left','interpreter','none','fontweight','normal','backgroundcolor','none');
                if LABELONSIGNONLY, set(h,'visible','off'); end 
                set(h(LABELONALL|idxtexthl>0),'color',1-data.plotconnoptions.BCOLOR,'tag','textstring','visible','on');%,'fontweight','bold');
                %set(h(seedz(data.displaytheserois(idxtext(data.displaytheserois(idxtext)<=size(z,1))))>0),'color',1-get(data.hfig,'color'),'visible','on');%,'fontweight','bold');
                roundang=5; %45;
                for n1=1:numel(h),tpos=get(h(n1),'position');ang=angle(tpos*[1;1i;0])/pi*180+data.plotconnoptions.FONTANGLE;if abs(mod(ang,360)-180)<90, set(h(n1),'rotation',roundang*round(ang/roundang)+180,'position',tpos*1.10,'horizontalalignment','right'); else set(h(n1),'rotation',roundang*round(ang/roundang),'position',tpos*1.10); end; end
                cumpatch({data.hfig,3},'addfield',data.displaytheserois(idxtext),'text',h);
            end
            set(h,'buttondownfcn',@conn_displayroi_menubuttondownfcn,'userdata',data.hfig,'interruptible','off');            
        else
            data.displaylabels=0;
        end
        
        %%cumpatch({data.hfig,1},'update');
        conn_display_windowbuttonmotionfcn('init',ringsquares);
        set(data.hfig,'windowbuttonmotionfcn',@conn_display_windowbuttonmotionfcn);
        if data.pausegui, conn_display_windowbuttonmotionfcn('pause',data.pausegui); end
        
        %data.plotsadd2=data.plotsadd2(idxsort);
        hold(data.plotaxes,'off');
        set(data.plotaxes,'ydir','normal');
        axis(data.plotaxes,'equal'); set(data.plotaxes,'xtick',[],'ytick',[]); %axis(data.plotaxes,'off');
        if ~data.view&&initxy&&data.displaylabels, 
            oldtlim=[nan nan nan nan];
            for nlim=1:2
                tlim=findobj(data.plotaxes,'tag','textstring','-or','tag','textstringpartial','visible','on');
                if ~isempty(tlim)
                    %tic; tlim=get(tlim,'extent'); toc
                    %temp=get(tlim,'position');
                    %if iscell(temp) temp=cell2mat(temp); end
                    %[nill,temp]=sort(temp,1);
                    %tlim=get(tlim(temp([1 end],:)),'extent');
                    tlim=get(tlim,'extent');
                    if iscell(tlim) tlim=cell2mat(tlim(:)); end
                    tlim=[tlim(:,1:2);tlim(:,1:2)+tlim(:,3:4)]; %[get(data.plotaxes,'xlim')',get(data.plotaxes,'ylim')']];
                    if ~isempty(bbox), tlim=[tlim;bbox(:,1:2)]; end
                    newtlim=[min(tlim(:,1))-1,max(tlim(:,1))+1,min(tlim(:,2))-1,max(tlim(:,2))+1];
                    if max(oldtlim-newtlim)<.01, break; end
                    set(data.plotaxes,'xlim',newtlim(1:2),'ylim',newtlim(3:4));
                    oldtlim=newtlim;
                end
            end
            data.axeslim={get(data.plotaxes,'xlim'),get(data.plotaxes,'ylim')};
        else set(data.plotaxes,'xlim',data.axeslim{1},'ylim',data.axeslim{2});
        end
        %if ~data.view&&data.displaylabels&&data.displaybrains~=1, set(data.plotaxes,'xlim',(1.5+.5*(data.displaylabels>0)+.5*data.displaybrains)*[-200,200]); end
        if data.display3d, data.display3d=0; datadisplay3d=1;else datadisplay3d=0; end
        set(hfig,'userdata',data);
        set(hcontrols(ishandle(hcontrols)),'enable','on');
        set(hfig,'pointer','arrow');
        %try, uicontrol(data.handles(6)); end

        if datadisplay3d
            idxkeep=idxtext(idxtexthl>0);
            c=mat2cell(cmap(round(1+(size(cmap,1)-1)*(1+K(idxkeep))/2),:),ones(numel(idxkeep),1),3);
            for n1=1:numel(data.source),idxc=find(idxkeep==data.source(n1));c(idxc,:)=repmat({[.25,.25,.25]},[numel(idxc),1]); end
            % ring placeholder xyz/2
            conn_mesh_display('','',[],...
                struct('sph_names',{data.names2reduced(data.displaytheserois(idxkeep))},...
                 'sph_xyz',[data.x(data.displaytheserois(idxkeep)),data.y(data.displaytheserois(idxkeep)),data.z(data.displaytheserois(idxkeep))],...
                 'sph_r',3*ones(numel(idxkeep),1),...
                 'sph_shapes',{data.names2(data.displaytheserois(idxkeep))},...
                 'sph_c',{c}),...%{repmat({[.9,.9,.9]},[1,numel(idxkeep)])}), ...
                z(data.displaytheserois(idxkeep(data.displaytheserois(idxkeep)<=size(z,1))),data.displaytheserois(idxkeep)).*sign(data.h(data.displaytheserois(idxkeep(data.displaytheserois(idxkeep)<=size(z,1))),data.displaytheserois(idxkeep))), ...
                [], .2, [0,-1e-8,1],[],data.defaultfilepath);
        end
        
    case 'clusters',
        figure(hfig);set(hfig,'pointer','watch');drawnow;
        rcircle=1.3*[sin(linspace(0,2*pi,64)'),cos(linspace(0,2*pi,64))']*diag([5,5]);
        h=findobj('tag','conn_displayroi_plot');
        if ~isempty(h),delete(h); end
        h=axes('units','norm','position',[.03,.06,.55,.88]);
        lim=[1,1,1;data.ref.dim];refminmax=sort([lim((dec2bin(0:7)-'0'+1)+repmat([0,2,4],[8,1])),ones(8,1)]*data.ref.mat(1:3,:)'*data.proj(:,1:2));
        temp=reshape(data.bgimage,size(data.bgx)); %convn(convn(reshape(data.bgimage,size(data.bgx)),conn_hanning(5),'same'),conn_hanning(5)','same'); 
        temp(isnan(temp))=0;
        temp=round(1+(.7+.3*temp/max(temp(:)))*(size(get(hfig,'colormap'),1)-1));
        data.refaxes=image(refminmax(1,1):scale:refminmax(end,1),refminmax(1,2):scale:refminmax(end,2),temp);hold on;
        %data.refaxes=imagesc(refminmax(1,1):scale:refminmax(end,1),refminmax(1,2):scale:refminmax(end,2),convn(convn(reshape(data.bgimage,size(data.bgx)),conn_hanning(5),'same'),conn_hanning(5)','same'));hold on;
        h=gca;data.buttondown=struct('h1',h);set(h,'tag','conn_displayroi_plot');

        cmap=jet(max(1,max(data.clusters)));
        EPS=1;
        %data.plotsadd2={};
        idxtext=[];
        for n1=1:N,
            hold on;
            if length(data.clusters)>=n1&&data.clusters(n1)>0,
                h=patch(data.x(n1)+rcircle(:,1)*max(.5,1+1e-3*data.z(n1)),data.y(n1)+rcircle(:,2)*max(.5,1+1e-3*data.z(n1)),max(1,data.z(n1))*EPS+.10+zeros(size(rcircle,1),1),'w');
                set(h,'edgecolor','none','facecolor',cmap(data.clusters(n1),:));
                hold on;
                set(h,'buttondownfcn',@conn_displayroi_menubuttondownfcn,'userdata',data.hfig,'interruptible','off');
                idxtext=[idxtext;n1];
                %h=text(data.x(n1),data.y(n1),max(1,data.z(n1))*EPS+.10+2,num2str(n1));
                %set(h,'fontsize',9,'color','w','horizontalalignment','center','interpreter','none','fontweight','bold','backgroundcolor','none');
                %set(h,'buttondownfcn',@conn_displayroi_menubuttondownfcn,'userdata',data);
            end
        end
        hold off;
        h=text(data.x(data.displaytheserois(idxtext)),data.y(data.displaytheserois(idxtext)),max(1,data.z(data.displaytheserois(idxtext)))*EPS+.10+2,{data.names2reduced{data.displaytheserois(idxtext)}});        
        set(h,'fontsize',8+CONN_gui.font_offset,'color','w','horizontalalignment','center','interpreter','none','fontweight','normal','backgroundcolor','none');
        %set(h,'buttondownfcn',@conn_displayroi_menubuttondownfcn,'userdata',data.hfig);
        set(gca,'ydir','normal');
        axis equal; axis off;
        set(hfig,'userdata',data);
        set(hfig,'pointer','arrow','visible','on');
end
end


function data=conn_displayroi_clusters(data,varargin)
NCLUSTERS=[];
LAMBDAPOS=[];
NPLOTS=[];
MAXROIS=[];
DOASK=false;
if isfield(data,'clusters_options')&&~isempty(data.clusters_options)
    switch(data.clusters_options.type)
        case 'hc'
            TYPE='hc';
            NCLUSTERS=data.clusters_options.groups;
            LAMBDAPOS=data.clusters_options.param;
        case 'none'
            TYPE='none';
            try, NPLOTS=data.plotconnoptions.NPLOTS; end
        otherwise, %case {'hcnet','symrcm','amd'}
            TYPE=data.clusters_options.type;
            MAXROIS=data.clusters_options.groups;
    end
else
    TYPE='hc';
    NCLUSTERS=nan;
end
if ~isempty(varargin)
    DOASK=true;
    TYPE=varargin{1};
    switch(TYPE)
        case 'hc'
            if numel(varargin)>1, NCLUSTERS=varargin{2}; end
            if numel(varargin)>2, LAMBDAPOS=varargin{3}; end
        case 'none'
            if numel(varargin)>1, NPLOTS=varargin{2}; end
        otherwise, %case {'hcnet','symrcm','amd'}
            if numel(varargin)>1, MAXROIS=varargin{2}; end
    end
end
if strcmp(TYPE,'hc')
    if DOASK||isempty(LAMBDAPOS)||isempty(NCLUSTERS),
        if isempty(LAMBDAPOS), LAMBDAPOS=.05; end
        if isnan(NCLUSTERS), NCLUSTERS=[]; end
        answer=conn_menu_inputdlg({'Number of groups/clusters (leave empty for automatic cutoff)','Hierarchical Clustering criteria: -1=Labels; 0=Functional; 1=Positional'},'clustering parameters',1,{num2str(NCLUSTERS),num2str(LAMBDAPOS)});
        try
            NCLUSTERS=str2num(answer{1});
            if isempty(NCLUSTERS), NCLUSTERS=nan; end
            LAMBDAPOS=str2num(answer{2});
        catch
            return;
        end
    end
elseif strcmp(TYPE,'none')
    if DOASK||isempty(NPLOTS)
        if isempty(NPLOTS), NPLOTS=10; try, NPLOTS=data.plotconnoptions.NPLOTS; end; end
        answer=conn_menu_inputdlg({'Number of groups/clusters'},'clustering parameters',1,{num2str(NPLOTS)});
        try
            data.plotconnoptions.NPLOTS=str2num(answer{1});
        catch
            return;
        end
    end
elseif strcmp(TYPE,'hcnet'),
    if DOASK||isempty(MAXROIS),
        if isempty(MAXROIS), MAXROIS=inf; end
        answer=conn_menu_inputdlg({'Maximum number of ROIs per brain display / cluster'},'display parameters',1,{num2str(MAXROIS)});
        try, MAXROIS=str2num(answer{1}); end
    end
elseif strcmp(TYPE,'symrcm'),
    if DOASK||isempty(MAXROIS),
        if isempty(MAXROIS), MAXROIS=inf; end
        answer=conn_menu_inputdlg({'Maximum number of ROIs per brain display / cluster'},'display parameters',1,{num2str(MAXROIS)});
        try, MAXROIS=str2num(answer{1}); end
    end
elseif strcmp(TYPE,'amd')
    if DOASK||isempty(MAXROIS),
        MAXROIS=inf;
        answer=conn_menu_inputdlg({'Maximum number of ROIs per brain display /cluster (set to inf for one brain display per network)'},'display parameters',1,{num2str(MAXROIS)});
        try, MAXROIS=str2num(answer{1}); end
    end
end
if strcmp(TYPE,'none')
    data.displaytheserois=sort(data.displaytheserois);
    data.xy2=zeros(length(data.names2),2); data.xy2(data.displaytheserois,:)=200*[cos(2*pi*(0:numel(data.displaytheserois)-1)'/numel(data.displaytheserois)),sin(2*pi*(0:numel(data.displaytheserois)-1)'/numel(data.displaytheserois))];
    data.clusters=zeros(length(data.names2),1); data.clusters(data.displaytheserois)=1+floor(data.plotconnoptions.NPLOTS*(0:numel(data.displaytheserois)-1)/numel(data.displaytheserois));
    data.xy2_clusters=zeros(length(data.names2),2); data.xy2_clusters(data.displaytheserois,:)=200*[cos(2*pi*(.5+floor(data.plotconnoptions.NPLOTS*(0:numel(data.displaytheserois)-1)/numel(data.displaytheserois))')/data.plotconnoptions.NPLOTS),cos(2*pi*(.5+floor(data.plotconnoptions.NPLOTS*(0:numel(data.displaytheserois)-1)/numel(data.displaytheserois))')/data.plotconnoptions.NPLOTS)];
    data.names_clusters={};
    data.clusters_options=struct('type','none','groups',[]);
    exstr='ROIs sorted manually';
    
elseif strcmp(TYPE,'hc')
    %if isempty(which('dendrogram')), error('Sorry. This option requires matlab statistics toolbox'); end
    N=length(data.names);
    N2=length(data.names2);
    data.displaytheserois=sort(data.displaytheserois);
    i1=intersect(1:N,data.displaytheserois);
    i2=intersect(1:N2,data.displaytheserois);
    %X=data.h;
    if 1 % note: hierarchical clustering based on average connectivity data, not on specific second-level stats in each analysis
        for nresults=1:size(data.results(1).data,4)
            ty=data.results(1).data(:,:,:,nresults);
            ty(isnan(ty))=0;
            y=mean(mean(ty,1),3);
            if nresults==1, X=zeros(size(data.results(1).data,4),numel(y)); end
            X(nresults,:)=y(:)';
        end        
    elseif 1
        X=sign(data.h).*abs(data.F);
    else
        X=sign(data.h);
        if isequal(data.statsname,'T')
            switch(data.side),
                case 1,p=data.p; Fthr=data.F;
                case 2,p=data.p2; Fthr=-data.F;
                case 3,p=2*min(data.p,data.p2);  Fthr=abs(data.F);
            end
            %set(data.handles(4),'enable','on','value',data.side);
        else
            p=data.p;
            Fthr=data.F;
            %set(data.handles(4),'enable','off','value',3);
        end
        switch(data.thrtype),
            case 1, X(p>data.thr)=0;
            case 2, X(data.P>data.thr)=0;
            case {3,4}, X(data.ttZb<data.ttZbthr)=0;
            case 5, X(Fthr<data.thr)=0;                
        end
    end
    X(isnan(X))=0;
    X(1:size(X,1)+1:size(X,1)*size(X,1))=0;
    %             if size(data.results(1).xX.X,2)==1,
    %                 X(1:size(X,1)+1:size(X,1)*size(X,1))=1*max(max(abs(X(i1,i2))));
    %             else X(1:size(X,1)+1:size(X,1)*size(X,1))=0;
    %             end
    X=X(i1,i2,:);
    Xt=X;%[X;[1e2+zeros(1,N),zeros(1,N2-N)]];
    xyz2t=data.xyz2(i2,:);%[data.xyz2(i2,:),[1e2+zeros(N,1);zeros(N2-N,1)]];
    names2t=data.names2(i2);
    xyz2t=xyz2t*diag([2,1,1]); % bias x-dir
    
    [j1,j2,j3]=find(X);
    Y=sqrt(max(0,permute(sum(abs(conn_bsxfun(@minus,Xt,permute(Xt,[1 3 2]))).^2,1),[3 2 1])-sparse(j1,j2,j3.^2,size(X,2),size(X,2))-sparse(j2,j1,j3.^2,size(X,2),size(X,2)))); % removes diagonal elements of correlation matrix from distance computation
    if LAMBDAPOS>=0
        Y2=sqrt(max(0,permute(sum(abs(conn_bsxfun(@minus,xyz2t,permute(xyz2t,[3 2 1]))).^2,2),[1 3 2])));
    else
        Y2=conn_wordld(names2t,names2t);
    end
    %figure
    %Y = sqrt(pdist(Xt', 'euclidean').^2-squareform(sparse(j1,j2,j3.^2,size(X,2),size(X,2))+sparse(j2,j1,j3.^2,size(X,2),size(X,2)))); % removes diagonal elements of correlation matrix from distance computation
    %Y2 = pdist(xyz2t, 'euclidean');
    I=tril(ones(size(X,2)),-1);
    Y0=sqrt((1-abs(LAMBDAPOS))*Y.^2/max(eps,mean(Y(I>0).^2))+abs(LAMBDAPOS)*Y2.^2/max(eps,mean(Y2(I>0).^2)));
    Y0(eye(size(Y0))>0)=0;
    Y0=(Y0+Y0')/2;
    Y0b=Y0(I>0)';%sqrt((1-abs(LAMBDAPOS))*Yb.^2/mean(Y.^2)+abs(LAMBDAPOS)*Y2b.^2/mean(Y2b.^2));
    Z = conn_statslinkage(Y0b, 'co');
    data.clusters_X=X;
    semic0=1-1/size(X,2);
    idx=conn_statsoptimalleaforder(Z,Y0); %reference: Bar-Joseph, Z., Gifford, D.K., and Jaakkola, T.S. (2001). Fast optimal leaf ordering for hierarchical clustering. Bioinformatics 17, Suppl 1:S22?9. PMID: 11472989
    %[H,t,idx]=conn_statsdendrogram(Z,0,'labels',data.names2(i2),'orientation','left');
    if isnan(NCLUSTERS) % automatic cutoff number of clusters
        [nill,idxmax]=max((1:size(Z,1))'/size(Z,1)-Z(:,3)/max(Z(:,3))); %disp(size(Z,1)-idxmax+1);
        %[nill,idxmax]=max(log((1:size(Z,1))'/size(Z,1))-(Z(:,3)/max(Z(:,3))).^2); %disp(size(Z,1)-idxmax+1);
        nclusters=size(Z,1)-idxmax+1;
    else nclusters=NCLUSTERS;
    end
    T = conn_statscluster(Z, 'maxclust', nclusters);
    data.clusters=zeros(N2,1); data.clusters(i2)=T; 
    if isempty(setdiff(i2,data.displaytheserois(idx))), [nill,nill,data.clusters(data.displaytheserois(idx))]=unique(data.clusters(data.displaytheserois(idx)),'stable'); end
    %for n1=1:NCLUSTERS,conn_disp(' ');conn_disp(['Cluster #',num2str(n1)]);conn_disp(strvcat(data.names2{data.clusters==n1}));end
    TT=[]; for n1=1:min(size(X,2),nclusters), TT(:,n1) = conn_statscluster(Z, 'maxclust', n1); end;
    if size(TT,2)>1,
        %TT=TT(:,2:end); TT=mean(detrend(cumsum([ones(1,size(TT,2));diff(TT(idx,:),1,1)~=0],1)./repmat(2:min(size(X,2),nclusters),[size(TT,1),1]),'constant'),2)';
        wblank=.125;
        TT=TT(:,2:end); TT=mean((-.5+cumsum([ones(1,size(TT,2));diff(TT(idx,:),1,1)~=0],1))./repmat(2:min(size(X,2),nclusters),[size(TT,1),1]),2)'; TT=TT-(max(TT)+min(TT))/2;
        xy=exp(1i*(linspace(-pi*(semic0-1e-4-semic0*wblank),pi*(semic0-1e-4-semic0*wblank),numel(i2))+semic0*wblank*2*pi*TT));
        xy_clusters=exp(1i*(linspace(-pi*(semic0-1e-4-semic0*.95),pi*(semic0-1e-4-semic0*.95),numel(i2))+semic0*.95*2*pi*TT));
    else
        xy=exp(1i*(linspace(-pi*(semic0-1e-4),pi*(semic0-1e-4),numel(i2))));
        xy_clusters=xy;
    end
    data.xy2=zeros(length(data.names2),2);
    data.xy2(data.displaytheserois(idx),:)=200*[real(xy)',imag(xy)'];
    data.xy2_clusters=zeros(length(data.names2),2);
    data.xy2_clusters(data.displaytheserois(idx),:)=200*[real(xy_clusters)',imag(xy_clusters)'];
    data.displaytheserois=data.displaytheserois(idx);
    data.names_clusters={};
    data.clusters_options=struct('type','hc','groups',NCLUSTERS,'param',LAMBDAPOS);
    exstr='ROIs sorted using hierarchical clustering';
    
else
    N=length(data.names);
    N2=length(data.names2);
    data.displaytheserois=sort(data.displaytheserois);
    i1=intersect(1:N,data.displaytheserois);
    i2=intersect(1:N2,data.displaytheserois);
    if data.thrtype>1,
        p=data.P;
    else
        if isequal(data.statsname,'T')
            switch(data.side),
                case 1,p=data.p;
                case 2,p=data.p2;
                case 3,p=2*min(data.p,data.p2);
            end
        else
            p=data.p;
        end
    end
    mask=p(i1,i2)<=data.thr;
    [Nr1,Nr2]=size(mask);
    [i,j,v]=find(mask);
    M=sparse(i,j,v,max(Nr1,Nr2),max(Nr1,Nr2)); % handles non-square matrices
    M=M|M';
    M=M|speye(size(M,1));
    exstr='';
    if strcmp(TYPE,'hcnet'),
        [idx0,Labels]=conn_hcnet(p(i1,i2));
        docluster=true;
        TT=[];%fliplr(Labels(:,1:end-1));
        exstr='ROIs sorted using minimum degree';
    elseif strcmp(TYPE,'symrcm'),
        idx0=symrcm(M);
        docluster=true;
        TT=[];
        exstr='ROIs sorted using reverse Cuthill-McKee';
    elseif strcmp(TYPE,'amd')
        idx0=amd(M,struct('dense',sqrt(size(M,1))));
        docluster=true;
        TT=[];
        exstr='ROIs sorted using minimum degree';
    end
    if docluster
        [p,q,r,s]=dmperm(M(idx0,idx0));
        if ~isequal(p(:)',1:numel(p)),conn_disp('warning, resorted ROIs'); end
        p=idx0(p);
    else
        p=idx0;
        aM=full(sum(M(:,p),1)>1);
        if all(aM)
            r=[1,numel(p)+1];
        else
            p=[p(~aM),p(aM)];
            r=[1,sum(~aM)+1,numel(p)+1];
        end
    end
    n=diff(r);
    i=find(n>1);
    Nlabels=zeros(size(M,1),1); % node labels
    i0=0;
    for i1=1:numel(i)
        if n(i(i1))>MAXROIS
            k=ceil(n(i(i1))/MAXROIS);
            Nlabels(p(r(i(i1)):r(i(i1)+1)-1))=i0+ceil((1:n(i(i1)))/n(i(i1))*k);
            i0=i0+k;
        else
            Nlabels(p(r(i(i1)):r(i(i1)+1)-1))=i0+1;
            i0=i0+1;
        end
    end
    if any(Nlabels==0)
        j=sum(Nlabels==0);
        if j>MAXROIS
            k=ceil(j/MAXROIS);
            Nlabels(Nlabels==0)=i0+ceil((1:j)/j*k);
            i0=i0+k;
        else
            Nlabels(Nlabels==0)=i0+1;
            i0=i0+1;
        end
    end
    rankp=p;rankp(rankp)=1:numel(rankp);
    [nill,idx]=sortrows([Nlabels,rankp(:)]);
    data.clusters=zeros(N2,1);
    data.clusters(i2)=Nlabels+(min(Nlabels)==0);
    data.names_clusters={};
    data.clusters_options=struct('type',TYPE,'groups',MAXROIS);
    
    semic0=1-1/size(mask,2);
    if ~isempty(TT)
        %TT=mean(detrend(cumsum([ones(1,size(TT,2));diff(TT(idx,:),1,1)~=0],1)./repmat(2:size(TT,2)+1,[size(TT,1),1]),'constant'),2)';
        %xy=exp(1i*(linspace(-pi*(semic0-1e-4-semic0*.25/2),pi*(semic0-1e-4-semic0*.25/2),numel(i2))+semic0*.25*pi*TT));
        wblank=.125;
        TT=detrend(mean((-.5+cumsum([ones(1,size(TT,2));diff(TT(idx,:),1,1)~=0],1))./repmat(2:min(size(X,2),nclusters),[size(TT,1),1]),2)','constant');
        xy=exp(1i*(linspace(-pi*(semic0-1e-4-semic0*wblank),pi*(semic0-1e-4-semic0*wblank),numel(i2))+semic0*wblank*2*pi*TT));
        xy_clusters=exp(1i*(linspace(-pi*(semic0-1e-4-semic0*.95),pi*(semic0-1e-4-semic0*.95),numel(i2))+semic0*.95*2*pi*TT));
    else
        xy=exp(1i*(linspace(-pi*(semic0-1e-4),pi*(semic0-1e-4),numel(i2))));
        xy_clusters=xy;
    end
    data.xy2=zeros(length(data.names2),2);
    data.xy2(data.displaytheserois(idx),:)=200*[real(xy)',imag(xy)'];
    data.xy2_clusters=zeros(length(data.names2),2);
    data.xy2_clusters(data.displaytheserois(idx),:)=200*[real(xy_clusters)',imag(xy_clusters)'];
    data.displaytheserois=data.displaytheserois(idx);
end
try, if ~isempty(exstr), set(data.handles(14),'string',exstr); end; end
end

function h=cumpatch(npatch,option,varargin)
persistent id patchobj patchopt patchvertices patchmask hplot ischanged;
h=[];
if isempty(id),
    id={npatch};
    npatch=1;
else
    if ischar(option)&&strcmp(option,'init'), % deletes data from old figures
        n=find(arrayfun(@(x)ishandle(id{x}{1}),1:numel(id))>0); 
        if numel(n)<numel(id), id=id(n); patchobj=patchobj(n); patchopt=patchopt(n); patchvertices=patchvertices(n); patchmask=patchmask(n); hplot=hplot(n); ischanged=ischanged(n); end
    end 
    ipatch=[]; for n1=1:numel(id), if isequal(id{n1},npatch), ipatch=n1; break; end; end
    if isempty(ipatch), ipatch=numel(id)+1; end
    id{ipatch}=npatch; 
    npatch=ipatch; 
end
if ischar(option)
    switch(option)
        case 'init',   
            patchobj{npatch}=[];
        case 'maskout'
            if isempty(patchmask), h=[]; 
            elseif numel(patchobj)>=npatch&&~isempty(patchobj{npatch}), h=[patchobj{npatch}.index;patchmask{npatch}];
            else h=zeros(2,0);
            end
        case 'addfield'
            if numel(patchobj)>=npatch&&~isempty(patchobj{npatch})
                index=varargin{1};
                fieldname=varargin{2};
                fieldvalue=varargin{3};
                if isempty(index)
                    patchobj{npatch}.(fieldname)=fieldvalue;
                else
                    [ok,idx]=ismember(index,patchobj{npatch}.index);
                    if isfield(patchobj{npatch},fieldname), tempvalues=patchobj{npatch}.(fieldname); end
                    tempvalues(idx(ok))=fieldvalue(ok);
                    patchobj{npatch}.(fieldname)=tempvalues;
                end
            end
        case 'mask',
            if numel(patchobj)>=npatch&&~isempty(patchobj{npatch})
                if isempty(varargin)
                    if isempty(ischanged{npatch}), return; end
                    set(hplot{npatch},'facevertexalphadata',patchobj{npatch}.facevertexalphadata(:));
                    if isfield(patchobj{npatch},'text')
                        tcolor=2*(mean(mean(patchobj{npatch}.facevertexcdata,1))>=.5)-1;
                        [ka,nill,kb]=unique(mean(patchobj{npatch}.facevertexalphadata,1));
                        if isfield(patchobj{npatch},'textemphasis'), tbase=patchobj{npatch}.textemphasis; else tbase=1; end
                        for kn=1:numel(ka),set(patchobj{npatch}.text(kb==kn),'color',max(0,min(1, .5+tbase*tcolor*(-.4+.8*[1 1 1]*ka(kn))))); end
                    end
                    ischanged{npatch}=[];
                    if isfield(patchobj{npatch},'coords'), patchmask{npatch}=ones(1,size(patchobj{npatch}.index,2));
                    else patchmask{npatch}=[];
                    end
                else
                    index=varargin{1};
                    pairs=patchobj{npatch}.index;
                    N=size(pairs,2);
                    mask=zeros(1,size(pairs,2));
                    umask=zeros(1,size(pairs,2));
                    if isempty(index),
                    elseif size(index,1)<2 % specific connections
                        ok=ismember(sort(pairs,1)',sort(reshape(real(index),2,[]),1)','rows')';
                        mask(ok)=1;
                    else % all connections to/from individual ROIs
                        h=real(index(end,:));
                        for n=1:size(pairs,1)
                            [ok,idx]=ismember(pairs(n,:),index(1,:));
                            mask(ok)=max(mask(ok),real(index(2,idx(ok))));
                            umask(ok)=max(umask(ok),imag(index(2,idx(ok))));
                        end
                    end
                    facevertexcdata=patchobj{npatch}.facevertexcdata;
                    facevertexalphadata=patchobj{npatch}.facevertexalphadata;
                    vertices=patchobj{npatch}.vertices;
                    mask=max(0,min(1,mask));
                    facevertexalphadata=facevertexalphadata.*repmat(.5+.5*mask,[size(facevertexalphadata,1),1]);
                    if any(umask)&&numel(umask)*size(facevertexalphadata,1)==size(facevertexcdata,1), 
                        umask=reshape(repmat(max(0,min(1,umask)),[size(facevertexalphadata,1),1]),size(facevertexcdata,1),1);
                        facevertexcdata=repmat(1-umask,1,3).*facevertexcdata + umask*[1 0 0];
                        if numel(facevertexalphadata)==numel(umask), facevertexalphadata=facevertexalphadata.*reshape(0+1*umask,size(facevertexalphadata)); end
                    elseif size(vertices,2)==3&&numel(mask)*size(facevertexalphadata,1)==size(facevertexcdata,1), % pop-up selected connections
                        tmask=reshape(repmat(max(0,min(1,mask)),[size(facevertexalphadata,1),1]),size(facevertexcdata,1),1);
                        vertices(:,3)=vertices(:,3)+100*tmask;
                        tmask=repmat(tmask,1,size(facevertexcdata,2));
                        try, bg=mean(get(id{npatch}{1},'color'));
                        catch, bg=.975;
                        end
                        facevertexcdata=facevertexcdata.*tmask + (.025+.95*bg)*(1-tmask); % almost-background-color
                        %facevertexcdata=facevertexcdata.*tmask + repmat(mean(facevertexcdata,2),1,size(facevertexcdata,2)).*(1-tmask);
                    end
                    %set(hplot{npatch},'facevertexcdata',facevertexcdata,'facevertexalphadata',max(.001,facevertexalphadata(:)));
                    set(hplot{npatch},'facevertexcdata',facevertexcdata,'facevertexalphadata',.001+.998*facevertexalphadata(:),'vertices',vertices);
                    if isfield(patchobj{npatch},'text')
                        tcolor=2*(mean(mean(patchobj{npatch}.facevertexcdata,1))>=.5)-1;
                        [ka,nill,kb]=unique(mask);
                        if isfield(patchobj{npatch},'textemphasis'), tbase=patchobj{npatch}.textemphasis; else tbase=1; end
                        for kn=1:numel(ka),set(patchobj{npatch}.text(kb==kn),'color',max(0,min(1, .5+tbase*tcolor*(-.4+.8*[1 1 1]*ka(kn))))); end
                    end
                    patchmask{npatch}=mask;
                    if nargout>0,
                        try
                            if isfield(patchobj{npatch},'C')
                                if isempty(index), 
                                    h=ones(1,0); 
                                else
                                    w=patchobj{npatch}.C;
                                    k=full(sparse(1,index(1,:),real(index(2,:)),1,max(max(index(1,:)),size(w,1))));
                                    %k=zeros(1,size(w,1));k(index(1,:))=real(index(2,:));
                                    if numel(k)>size(w,2), w=[w zeros(size(w,1),numel(k)-size(w,2))]; end
                                    if numel(k)>size(w,1), w=[w; zeros(numel(k)-size(w,1),size(w,2))]; end
                                    k=k*w;
                                    h=min(1,max(h,k(index(1,:))));
                                end
                            else
                                for n=1:size(pairs,1)
                                    [ok,idx]=ismember(pairs(n,:),index(1,:));
                                    h=max(h,accumarray(idx(ok)',mask(ok),[numel(h),1],@max)');
                                end
                            end
                        end
                    end
                    ischanged{npatch}=index;
                end
            end
        case 'smooth&update',
            if numel(patchobj)>=npatch&&~isempty(patchobj{npatch})
                if ~isempty(varargin), lbundl=varargin{1}; end
                vertices=cat(1,patchobj{npatch}.vertices);
                if lbundl>0
                    coords=cat(2,patchobj{npatch}.coords);
                    N=size(coords,2);
                    vertices=permute(reshape(vertices,size(vertices,1)/N,N,size(vertices,2)),[2,1,3]);
                    if 1, % avoids duplicated lines
                        vertices=cat(1,vertices,vertices(:,[end/2+1:end,1:1:end/2],:));
                        if size(coords,1)<5, coords=cat(2,coords,coords([2,1,3:end],:));
                        else coords=cat(2,coords,coords([2,1,4,3,5:end],:));
                        end
                        N=2*N;
                    end
                    d1=sum(abs(coords).^2,1);
                    if N<1e4, 
                        H=exp(-(sqrt(max(0,repmat(d1,[N,1])+repmat(d1',[1,N])-2*real(coords'*coords)))/max(eps,200*lbundl)).^2);
                    elseif 0
                        H=speye(N);
                    else
                        tdmax=sqrt(-log(.001))*max(eps,200*lbundl);
                        i=[];j=[];k=[];
                        for n1=1:N, 
                            td=d1+d1(n1)-2*real(coords(:,n1)'*coords); 
                            ti=find(td<tdmax);
                            if ~isempty(ti)
                                j=[j ti];
                                i=[i n1+zeros(1,numel(ti))];
                                k=[k exp(-(td(ti)/max(eps,200*lbundl)).^2)];
                            end
                        end
                        H=sparse(i,j,k,N,N);
                    end
                    H=sparse(1:N,1:N,1./max(eps,sum(H,2)))*H;
                    newvertices=vertices;
                    newvertices(:,:)=full(H*vertices(:,:));
                    newvertices(:,:,3)=vertices(:,:,3);
                    if 1, % avoids duplicated lines
                        vertices=vertices(1:end/2,:,:);
                        newvertices=newvertices(1:end/2,:,:);
                        coords=coords(:,1:end/2);
                        N=N/2;
                    end
                    p1=(max(0,sqrt(max(0,sum(vertices(:,:,1:2).^2,3)))/200-.75)/.25).^2;
                    p2=0*repmat(exp(-(abs(coords(1,:)-coords(2,:))'/50).^2),[1,size(vertices,2)]);
                    p=repmat(min(1,p1+p2),[1,1,size(vertices,3)]);
                    vertices=vertices.*p + newvertices.*(1-p);
                    vertices=reshape(permute(vertices,[2,1,3]),[],size(vertices,3));
                end
                temp=struct('vertices',vertices, 'faces',cat(1,patchobj{npatch}.faces), 'facevertexcdata',cat(1,patchobj{npatch}.facevertexcdata), 'facevertexalphadata',cat(1,patchobj{npatch}.facevertexalphadata));
                hplot{npatch}=patch(temp,patchopt{npatch}{:},'facecolor','inter','facealpha','inter'); 
                if isfield(patchobj{npatch},'index'), temp.index=cat(2,patchobj{npatch}.index); temp.facevertexalphadata=reshape(temp.facevertexalphadata,[],size(temp.index,2)); end
                if isfield(patchobj{npatch},'index')&&size(temp.index,1)==2, temp.C=sparse([temp.index(1,:) temp.index(2,:)],[temp.index(2,:) temp.index(1,:)],1); end
                if isfield(patchobj{npatch},'coords'), temp.coords=cat(2,patchobj{npatch}.coords); end
                patchobj{npatch}=temp;
                if isfield(patchobj{npatch},'index'), patchmask{npatch}=ones(1,size(patchobj{npatch}.index,2));
                else patchmask{npatch}=[];
                end
                ischanged{npatch}=[];
                try, set(hplot{npatch},'alphadatamapping','none'); end
                h=hplot{npatch};
            end
        case 'update', 
            if numel(patchobj)>=npatch&&~isempty(patchobj{npatch})
                if 0
                    aLines=cat(3,patchobj{npatch}.vertices);
                    mLines=(aLines(1:end/2,:,:)+aLines(end:-1:end/2+1,:,:))/2;
                    [Lines,LineWidths]=conn_menu_bundle(mLines,[],[],5,true);
                    for n1=1:numel(patchobj{npatch}), patchobj{npatch}(n1).vertices=[Lines(:,:,n1,end); flipud(Lines(:,:,n1,end))]+aLines(:,:,n1)-[mLines(:,:,n1);flipud(mLines(:,:,n1))]; end
                end
                temp=struct('vertices',cat(1,patchobj{npatch}.vertices), 'faces',cat(1,patchobj{npatch}.faces), 'facevertexcdata',cat(1,patchobj{npatch}.facevertexcdata), 'facevertexalphadata',cat(1,patchobj{npatch}.facevertexalphadata));
                hplot{npatch}=patch(temp,patchopt{npatch}{:},'facecolor','inter','facealpha','inter'); 
                if isfield(patchobj{npatch},'index'), temp.index=cat(2,patchobj{npatch}.index); temp.facevertexalphadata=reshape(temp.facevertexalphadata,[],size(temp.index,2)); end
                if isfield(patchobj{npatch},'index')&&size(temp.index,1)==2, temp.C=sparse([temp.index(1,:) temp.index(2,:)],[temp.index(2,:) temp.index(1,:)],1); end
                if isfield(patchobj{npatch},'coords'), temp.coords=cat(2,patchobj{npatch}.coords); end
                patchobj{npatch}=temp;
                if isfield(patchobj{npatch},'index'), patchmask{npatch}=ones(1,size(patchobj{npatch}.index,2));
                else patchmask{npatch}=[];
                end
                ischanged{npatch}=[];
                try, set(hplot{npatch},'alphadatamapping','none'); end
                h=hplot{npatch};
            end
    end
else
    if numel(patchobj)<npatch||isempty(patchobj{npatch})
        patchobj{npatch}=option; 
        patchopt{npatch}=varargin;
        patchvertices{npatch}=size(option.vertices,1);
    else
        m=numel(patchobj{npatch});
        if isfield(option,'coords'), patchobj{npatch}(m+1)=struct('vertices',option.vertices,'faces',patchvertices{npatch}+option.faces,'facevertexcdata',option.facevertexcdata,'facevertexalphadata',option.facevertexalphadata,'coords',option.coords,'index',option.index);
        else patchobj{npatch}(m+1)=struct('vertices',option.vertices,'faces',patchvertices{npatch}+option.faces,'facevertexcdata',option.facevertexcdata,'facevertexalphadata',option.facevertexalphadata);
        end
        patchvertices{npatch}=patchvertices{npatch}+size(option.vertices,1);
    end
end
end

function conn_displayroi_menubuttondownfcn(varargin)
temp=get(gcbo,'userdata');
if iscell(temp)
    data=get(temp{1},'userdata');
    idx2=temp{2};
    okin=true;
%     set(data.handles(21),'string',[temp{2},'      ',temp{3}]);
else
    data=get(temp,'userdata');
    idx2=[];
end
if isempty(idx2)
    xyz=get(data.buttondown.h1,'currentpoint');
    x=xyz(1,1);y=xyz(1,2);z=0;
    [mind,idx]=min(sqrt(abs(data.x(data.displaytheserois)-x).^2+abs(data.y(data.displaytheserois)-y).^2)+1e-10*abs(data.z(data.displaytheserois)-z).^2);
    idx2=data.displaytheserois(idx);
    okin=sqrt(max(0,x.^2+y.^2))<216;
end
txt=data.names2{idx2};
if 0
    h=findobj('tag','conn_displayroi_menubuttondownfcn');if isempty(h), h=figure('units','pixels','position',[get(0,'pointerlocation')-[125,-200],250,40]);else, figure(h); end;
    set(h,'units','pixels','position',[get(0,'pointerlocation')-[125,-100],0,0]+[0,0,1,1].*get(h,'position'),'menubar','none','numbertitle','off','color','k','tag','conn_displayroi_menubuttondownfcn');
    clf(h);text(0,1,['x,y,z = (',num2str(data.xyz2(idx2,1),'%1.0f'),',',num2str(data.xyz2(idx2,2),'%1.0f'),',',num2str(data.xyz2(idx2,3),'%1.0f'),') mm'],'color','y','fontweight','bold','horizontalalignment','center','fontsize',9);
    text(0,0,['(',num2str(idx),') : ',txt],'color','y','fontweight','bold','horizontalalignment','center','fontsize',9,'interpreter','none');set(gca,'units','norm','position',[0,0,1,1],'xlim',[-1,1],'ylim',[-.5,1.5],'visible','off');
else
    if okin
        txt={['x,y,z = (',num2str(data.xyz2(idx2,1),'%1.0f'),',',num2str(data.xyz2(idx2,2),'%1.0f'),',',num2str(data.xyz2(idx2,3),'%1.0f'),') mm'],...
            [' : ',txt]};
        set(data.handles(21),'string',[txt{1},'      ',txt{2}]);
        if data.displayconnectionstats,
            idx2a=find(data.list2(:,1)==idx2);
            idx2b=find(data.list2(:,2)==idx2);
        elseif ~isempty(data.list2visible)
            idx2a=find(data.list2(data.list2visible,1)==idx2);
            idx2b=find(data.list2(data.list2visible,2)==idx2);
        else idx2a=[]; idx2b=[];
        end
        if ~isempty(idx2a)||~isempty(idx2b), set(data.handles(8),'value',[idx2a;idx2b]);end
        if ~isempty(idx2a), set(data.handles(8),'listboxtop',min(idx2a));
        elseif ~isempty(idx2b), set(data.handles(8),'listboxtop',min(idx2b));
        end
        %if ~isempty(idx2a)||~isempty(idx2b), conn_displayroi('list2'); end
    else
        set(data.handles(21),'string','');
        set(data.handles(8),'value',[]);
        %conn_displayroi('list2'); 
    end
end
%hc=get(0,'children');if length(hc)>0&&hc(1)~=h,hc=[h;hc(hc~=h)];set(0,'children',h); end
end


function conn_display_windowbuttonmotionfcn(option,varargin)
persistent data busy pause down;
if isempty(busy), busy=0; end
if isempty(pause), pause=0; end
if isempty(down), down=0; end
if isempty(data), ndata=1; 
else
    if nargin>0&&ischar(option)&&strcmp(option,'init'), 
        data=data(arrayfun(@(x)ishandle(data(x).gtf),1:numel(data))>0); % deletes data from old figures
        id=varargin{1}.gcf;
    else id=gcbf;
    end
    if isempty(data), ndata=1;
    else
        ndata=find(arrayfun(@(x)isequal(data(x).gtf,id),1:numel(data)),1);
        if isempty(ndata), ndata=numel(data)+1; end
    end
end
if nargin>0&&ischar(option)&&strcmp(option,'pause'), 
    pause=varargin{1};
elseif nargin>0&&ischar(option)&&strcmp(option,'init')
    busy=true;
    data(ndata).x=varargin{1}.x;
    data(ndata).y=varargin{1}.y;
    data(ndata).mx=mean(data(ndata).x(1:ceil(end/2),:),1);
    data(ndata).my=mean(data(ndata).y(1:ceil(end/2),:),1);
    if isempty(data(ndata).x), data(ndata).dx=[]; data(ndata).dy=[];
    else
        data(ndata).dx=data(ndata).x([2:end 1],:)-data(ndata).x;
        data(ndata).dy=data(ndata).y([2:end 1],:)-data(ndata).y;
    end
    data(ndata).index=varargin{1}.index;
    data(ndata).cluster=varargin{1}.cluster;
    data(ndata).h=varargin{1}.h;
    data(ndata).names=varargin{1}.names;
    data(ndata).gta=varargin{1}.gca;
    data(ndata).gtf=varargin{1}.gcf;
    data(ndata).gtt=varargin{1}.gct;
    data(ndata).highight.h1=findobj(varargin{1}.gcf,'tag','highlight');
    %data(ndata).highight.h2=findobj(gcbf,'tag','highlight_pointer');
    busy=false;
elseif numel(data)>=ndata&&~isempty(data(ndata).x)
    if ~isempty(gcbf)&&isequal(gcbf,gcf)&&~busy
        busy=true;
        if nargin>0&&ischar(option)&&strcmp(option,'down'), down=1;
        elseif nargin>0&&ischar(option)&&strcmp(option,'up'), down=0; busy=false; return;
        end
        hfig=data(ndata).gtf;
        set(hfig,'units','pixels');
        p1=get(0,'pointerlocation');
        p2=get(hfig,'position');
        p3=get(0,'screensize');
        p2(1:2)=p2(1:2)+p3(1:2)-1; % note: fix issue when connecting to external monitor/projector
        pos=(p1-p2(1:2));
        posfig=(p1-p2(1:2))./p2(3:4);
        set(hfig,'currentpoint',pos);
        pos=[0 0];
        try
            pos=get(data(ndata).gta,'currentpoint');
            pos=pos(1,1:2);
        end
        npos=sqrt(max(0,sum(pos.^2)));
        [nill,idx]=min(mean((pos(1)-data(ndata).x).^2+(pos(2)-data(ndata).y).^2,1));
        dlim=264;
        if npos<200||(npos>200*(1.11+.35)&&nill>16^2),idx=[]; end
        xlim=get(data(ndata).gta,'xlim');
        ylim=get(data(ndata).gta,'ylim');
        if numel(xlim)==2&&pos(1)>=xlim(1)&&pos(1)<=xlim(2)&&numel(ylim)==2&&pos(2)>=ylim(1)&&pos(2)<=ylim(2), set(hfig,'Pointer','crosshair'); 
        else set(hfig,'Pointer','arrow'); 
        end
        %if npos>=200&&npos<=216, [nill,idx]=min(mean((pos(1)-data(ndata).x).^2+(pos(2)-data(ndata).y).^2,1));
        %else idx=find(all(data(ndata).dx.*(pos(1)-data(ndata).x)+data(ndata).dy.*(pos(2)-data(ndata).y)>=0,1),1);
        %end
        if ~down&&pause
            if ~isempty(idx), % mouse inside an ROI (behavior when paused and mouse not clicked)
                set(data(ndata).gtt,'string',sprintf('connectivity with %s',data(ndata).names{idx}));%,'visible','on');
            elseif ~isempty(data(ndata).highight.h1)
                c0=mean(get(gcbf,'color')<.5);
                active=get(data(ndata).highight.h1,'position');
                if iscell(active), active=cell2mat(active(:)); end
                xpos=repmat(posfig,[size(active,1),1]);
                md=find(min(min(xpos-active(:,1:2),active(:,1:2)+active(:,3:4)-xpos),[],2)>-.001)';
                set(data(ndata).highight.h1,'foregroundcolor',[.4 .4 .4]+.2*c0);%,{'position'},xpos1);
                if ~isempty(md),
                    set(data(ndata).highight.h1(md),'foregroundcolor',[.0 .0 .0]+1*c0);%,{'position'},xpos1);
                end
            end
        elseif (down||~pause)&&isempty(idx)
            npos=sqrt(max(0,sum(pos.^2)));
            if npos<=216  % mouse inside ring (note: see rsquarelong)
                if npos>200, pos=pos/npos*200; end
                d=max(eps,sqrt(max(0,(pos(1)-data(ndata).mx).^2+(pos(2)-data(ndata).my).^2)));
                [mind,idx]=min(d);
                mind=200-npos; 
                %w=exp(-(10*(d/mind-1)).^2);
                w=exp(-((d-mind)*(10+mind)/max(eps,mind)/25).^2);
                set(data(ndata).gtt,'string','');%sprintf('FWHM = %d ROIs around %s',nnz(w>.5),data(ndata).names{idx}));%,'visible','on');
                v=cumpatch({data(ndata).gtf,1},'mask',[data(ndata).index; w]);
                if isequal(size(v),size(w)), z=1i*w+v; else z=w; end
                cumpatch({data(ndata).gtf,2},'mask',[data(ndata).index; z]);
                cumpatch({data(ndata).gtf,3},'mask',[data(ndata).index; z]);
            else %if numel(xlim)==2,%&&pos(1)>=xlim(1)&&pos(1)<=xlim(2) %&& numel(ylim)==2&&pos(2)>=ylim(1)&&pos(2)<=ylim(2) % mouse outside of ring
                conn_displayroi(hfig,[],'list2clear');
%                 w=ones(size(data(ndata).index));
%                 v=cumpatch({data(ndata).gtf,1},'mask',[data(ndata).index; w]);
%                 if isequal(size(v),size(w)), z=w+v; else z=w; end
%                 cumpatch({data(ndata).gtf,2},'mask',[data(ndata).index; z]);
%                 cumpatch({data(ndata).gtf,3},'mask',[data(ndata).index; z]);
%                 set(data(ndata).gtt,'string','');%,'visible','on');
%                 cumpatch({data(ndata).gtf,1},'mask');
%                 %cumpatch({data(ndata).gtf,2},'mask');
%                 cumpatch({data(ndata).gtf,3},'mask');%,[data(ndata).index; zeros(1,numel(data(ndata).index))]);
%                 %set(h,'facealpha',1);
            %else
            %    set(data(ndata).gtt,'string','','visible','on');
            end
        elseif (down||~pause)  
            if npos>216&&npos<=200*(1.11+.35), % mouse inside cluster (behavior when unpaused or mouse clicked) 
                w=data(ndata).cluster==data(ndata).cluster(idx);
                set(data(ndata).gtt,'string','');%sprintf('FWHM = %d ROIs around %s',nnz(w>.5),data(ndata).names{idx}));%,'visible','on');
                v=cumpatch({data(ndata).gtf,1},'mask',[data(ndata).index; w]);
                if isequal(size(v),size(w)), z=1i*w+v; else z=w; end
                cumpatch({data(ndata).gtf,2},'mask',[data(ndata).index; z]);
                cumpatch({data(ndata).gtf,3},'mask',[data(ndata).index; z]);
            else % mouse inside an ROI (behavior when unpaused or mouse clicked)
                set(data(ndata).gtt,'string',sprintf('connectivity with %s',data(ndata).names{idx}));%,'visible','on');
                w=data(ndata).index==data(ndata).index(idx);
                %w=(1:numel(data(ndata).index))==idx;
                v=cumpatch({data(ndata).gtf,1},'mask',[data(ndata).index; w]);
                if isequal(size(v),size(w)), z=1i*w+v; else z=w; end
                cumpatch({data(ndata).gtf,2},'mask',[data(ndata).index; z]);
                cumpatch({data(ndata).gtf,3},'mask',[data(ndata).index; z]);
            end
        end
        busy=false;
    end
end
end

function conn_displayroi_keypress(hdl,event)
%     option='';
%     %disp(event.Key)
%     if isequal(event.Key,'space')
%         conn_displayroi('pausegui');
%     end
end

function [filename_rois,filename_sources,viewrex]=conn_displayroi_selectfiles(filename_rois,filename_sources,viewrex)
global CONN_gui;
if nargin<2||isempty(filename_sources), filename_sources=''; end
if nargin<3||isempty(viewrex), viewrex=0; end
if isempty(CONN_gui)||~isfield(CONN_gui,'font_offset'), conn_font_init; end

filename_rois0=filename_rois;
filename_sources0=filename_sources;
thfig=dialog('units','norm','position',[.3,.4,.4,.25],'windowstyle','normal','name','REX interface','color','w','resize','on');
uicontrol(thfig,'style','text','units','norm','position',[.1,.75,.8,.20],'string',{'Explore model effects and connectivity values','within individual connections/clusters'},'backgroundcolor','w','fontsize',9+CONN_gui.font_offset,'fontweight','bold');
ht1=uicontrol(thfig,'style','popupmenu','units','norm','position',[.1,.50,.8,.20],'string',{'clusters of interest in current analysis','others clusters of interest (select exported mask/clusters file)'},'fontsize',8+CONN_gui.font_offset,'horizontalalignment','left','callback',@conn_displayroi_selectfiles_callback1,'tooltipstring','Define the clusters of interest');
%ht2=uicontrol(thfig,'style','popupmenu','units','norm','position',[.1,.30,.8,.20],'string',{'effect/activation/connectivity values in current analysis','other effect/activation/connectivity values (select second-level SPM.mat file)'},'fontsize',8+CONN_gui.font_offset,'horizontalalignment','left','callback',@conn_displayroi_selectfiles_callback2,'tooltipstring','Define the activation/connectivity values');
%if ~isempty(viewrex), ht3=uicontrol(thfig,'style','checkbox','units','norm','position',[.1,.3,.8,.1],'string','enable REX gui','value',0,'fontsize',8+CONN_gui.font_offset,'horizontalalignment','left','backgroundcolor','w','tooltipstring','Displays REX gui interface for additional options'); end
uicontrol(thfig,'style','pushbutton','string','OK','units','norm','position',[.1,.01,.38,.2],'callback','uiresume','fontsize',8+CONN_gui.font_offset);
uicontrol(thfig,'style','pushbutton','string','Cancel','units','norm','position',[.51,.01,.38,.2],'callback','delete(gcbf)','fontsize',8+CONN_gui.font_offset);
uiwait(thfig);
ok=ishandle(thfig);
if ok, 
    %if ~isempty(viewrex), viewrex=get(ht3,'value'); end
    delete(thfig);
else   filename_rois=[]; filename_sources=[]; 
end

    function conn_displayroi_selectfiles_callback1(varargin)
        if get(ht1,'value')==1
            filename_rois=filename_rois0;
        else
            tfilename=filename_rois;
            if iscell(tfilename)&&~isempty(tfilename), tfilename=tfilename{1}; end
            [tfilename,tpathname]=conn_fileutils('uigetfile','*.nii; *.img','Select mask/clusters file',tfilename);
            if ischar(tfilename), filename_rois=fullfile(tpathname,tfilename);
            else
                filename_rois=filename_rois0;
                set(ht1,'value',1);
            end
        end
    end
    function conn_displayroi_selectfiles_callback2(varargin)
        if get(ht2,'value')==1
            filename_sources=filename_sources0;
        else
            [tfilename,tpathname]=conn_fileutils('uigetfile','*.mat','Select SPM.mat file',filename_sources);
            if ischar(tfilename), filename_sources=fullfile(tpathname,tfilename);
            else
                filename_sources=filename_sources0;
                set(ht2,'value',1);
            end
        end
    end
end

function simfilename=conn_displayroi_simfilename(spmfile,THR_TYPE,THR,listrois)
if nargin==2&&isequal(THR_TYPE,'all')
    simfilename=fullfile(fileparts(spmfile),'nonparametricroi_p*.mat');
else
    indexfile=fullfile(fileparts(spmfile),'nonparametricroi_pindex.mat');
    listrois=listrois(:)'; %listrois=reshape(unique(listrois),1,[]);
    dosave=true;
    if ~conn_existfile(indexfile), 
        idx=1;
        listsrois={listrois};
    else
        listsrois={}; conn_loadmatfile(indexfile,'listsrois');
        idx=find(cellfun(@(x)isequal(x,listrois),listsrois),1);
        if isempty(idx), idx=numel(listsrois)+1; 
        else dosave=false;
        end
        listsrois{idx}=listrois;
    end    
    if dosave, conn_savematfile(indexfile,'listsrois'); end
    simfilename=char(arrayfun(@(a,b)fullfile(fileparts(spmfile),sprintf('nonparametricroi_p%d_%.8f_i%d.mat',a,b,idx)),THR_TYPE,THR,'uni',0));
    if ~isempty(conn_fileutils('dir',conn_prepend('parallel_*_',simfilename))), conn_process('results_nonparametric_collapse',conn_server('util_localfile_filesep',[],simfilename)); end
end
end

function uiwrap(h)
global CONN_gui;
if ~isfield(CONN_gui,'uicontrol_border'), CONN_gui.uicontrol_border=2; end
if isfield(CONN_gui,'isjava')&&~CONN_gui.isjava, return; end
for nh=1:numel(h)
    hpar=get(h(nh),'parent');
    bgcolor=get(hpar,'color');
    set(h(nh),'units','pixels');
    tpos=get(h(nh),'position');
    switch(get(h(nh),'style'))
        case 'pushbutton'
            htb=[uicontrol('style','frame','units','pixels','position',tpos+[tpos(3)-2*CONN_gui.uicontrol_border,0,2*CONN_gui.uicontrol_border-tpos(3),0],'foregroundcolor',bgcolor,'backgroundcolor',bgcolor,'units','norm','parent',hpar),...
                uicontrol('style','frame','units','pixels','position',tpos+[0,tpos(4)-CONN_gui.uicontrol_border,0,CONN_gui.uicontrol_border-tpos(4)],'foregroundcolor',bgcolor,'backgroundcolor',bgcolor,'units','norm','parent',hpar),...
                uicontrol('style','frame','units','pixels','position',tpos+[0,0,0,2*CONN_gui.uicontrol_border-tpos(4)],'foregroundcolor',bgcolor,'backgroundcolor',bgcolor,'units','norm','parent',hpar),...
                uicontrol('style','frame','units','pixels','position',tpos+[0,0,CONN_gui.uicontrol_border-tpos(3),0],'foregroundcolor',bgcolor,'backgroundcolor',bgcolor,'units','norm','parent',hpar)];
        case 'listbox',
            htb=[uicontrol('style','frame','units','pixels','position',tpos+[tpos(3)-2*CONN_gui.uicontrol_border,0,2*CONN_gui.uicontrol_border-tpos(3),0],'foregroundcolor',bgcolor,'backgroundcolor',bgcolor,'units','norm','parent',hpar),...
                uicontrol('style','frame','units','pixels','position',tpos+[0,tpos(4)-CONN_gui.uicontrol_border,0,CONN_gui.uicontrol_border-tpos(4)],'foregroundcolor',bgcolor,'backgroundcolor',bgcolor,'units','norm','parent',hpar),...
                uicontrol('style','frame','units','pixels','position',tpos+[0,-1,0,CONN_gui.uicontrol_border+2-tpos(4)],'foregroundcolor',bgcolor,'backgroundcolor',bgcolor,'units','norm','parent',hpar),...
                uicontrol('style','frame','units','pixels','position',tpos+[0,0,CONN_gui.uicontrol_border-tpos(3),0],'foregroundcolor',bgcolor,'backgroundcolor',bgcolor,'units','norm','parent',hpar)];
        case 'edit'
            htb=[uicontrol('style','frame','units','pixels','position',tpos+[tpos(3)-2*CONN_gui.uicontrol_border,0,2*CONN_gui.uicontrol_border-tpos(3),0],'foregroundcolor',bgcolor,'backgroundcolor',bgcolor,'units','norm','parent',hpar),...
                uicontrol('style','frame','units','pixels','position',tpos+[0,tpos(4)-CONN_gui.uicontrol_border,0,CONN_gui.uicontrol_border-tpos(4)],'foregroundcolor',bgcolor,'backgroundcolor',bgcolor,'units','norm','parent',hpar),...
                uicontrol('style','frame','units','pixels','position',tpos+[0,0,0,CONN_gui.uicontrol_border+1-tpos(4)],'foregroundcolor',bgcolor,'backgroundcolor',bgcolor,'units','norm','parent',hpar),...
                uicontrol('style','frame','units','pixels','position',tpos+[0,0,CONN_gui.uicontrol_border-tpos(3),0],'foregroundcolor',bgcolor,'backgroundcolor',bgcolor,'units','norm','parent',hpar)];        
    end
    set(h(nh),'units','norm');
end
end

function ok=conn_displayroi_randomise(data,THR_TYPE,THR,dogui)
ok=true;
simfilename=conn_displayroi_simfilename(data.roifile,THR_TYPE,THR,data.displaytheserois);
if THR==0, niters=1000;
else niters=max(1000,round(1/THR));
end
[fpath,fname,fext]=fileparts(data.initfile);
if strcmp([fname,fext],'SPM.mat')
    maskfile=fullfile(fileparts(data.roifile),'mask.nii');
    if ~conn_existfile(maskfile), maskfile=fullfile(fileparts(data.roifile),'mask.img'); end
    if conn_existfile(maskfile), mask=conn_fileutils('spm_read_vols',maskfile)>0;
    else mask=[];
    end
else mask=[];
end
SIDE=1:3;
THR=THR+[0 0 0];
THR_TYPE=THR_TYPE+[0 0 0];
N=length(data.names);
% update permutation tests
if isfield(data.results,'data'), Y=permute(data.results(1).data,[1,3,4,2]); 
else Y=permute(cat(4,data.results.y),[1,3,4,2]);
end
Y=Y(:,:,data.displaytheserois(data.displaytheserois<=N),data.displaytheserois);
if ~isempty(mask), mask=mask(data.displaytheserois(data.displaytheserois<=N),data.displaytheserois); end
try
    conn_randomise(data.results(1).xX.X,Y,data.results(1).c,data.results(1).c2,[],THR,THR_TYPE,SIDE,niters,simfilename,[],'matrix',mask);
    ok=true;
catch
   ok=false;
end
    % 'X','Y','c','m','THR','THR_TYPE','SIDE','niters','simfilename'
    % c=data.results(1).c;
    % m=data.results(1).c2;
    % X=data.results(1).xX.X;
    % conn_savematfile(tfilename,'X','Y','c','m','THR','THR_TYPE','SIDE','niters','simfilename');
    % conn_jobmanager('options','profile',parallel);
    % info=conn_jobmanager('submit','orphan_results_nonparametric',N,N,[],tfilename);
    % info=conn_jobmanager(info,'','donotupdate');
    % ok=true;
end




