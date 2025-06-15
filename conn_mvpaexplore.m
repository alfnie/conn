function fh=conn_mvpaexplore
global CONN_x CONN_h CONN_gui;

if isempty(CONN_gui)||~isfield(CONN_gui,'font_offset'), conn_font_init; end
filepathresults=fullfile(CONN_x.folders.firstlevel_vv,CONN_x.vvAnalyses(CONN_x.vvAnalysis).name);
[iroi,isnew,ncomp]=cellfun(@(x)conn_v2v('match_extended',x),CONN_x.vvAnalyses(CONN_x.vvAnalysis).measures);
if any(isnew), conn_msgbox('Sorry, this option is not available until the associated first-level MVPA analyses have been run','',2); return; end

if CONN_gui.usehighres, filename=fullfile(fileparts(which('conn')),'utils','surf','referenceT1_icbm.nii');
else filename=fullfile(fileparts(which('conn')),'utils','surf','referenceT1_trans.nii');
end
%filename=fullfile(fileparts(which('conn')),'utils','surf','referenceT1_icbm.nii');
[voldata,volref]=conn_vol_read(filename);
imat=pinv(volref.mat);

filepath=CONN_x.folders.preprocessing;
outcomenames=CONN_x.vvAnalyses(CONN_x.vvAnalysis).measures;
outcomeisource=[];for n1=1:length(outcomenames),
    [outcomeisource(n1),isnew,outcomencompsource(n1)]=conn_v2v('match_extended',outcomenames{n1});
    if isnew, error('Measure %s not found in global measures list. Please re-run first-level analyses',outcomenames{n1}); end
end
validsources=1:numel(outcomenames);
nconditions=length(CONN_x.Setup.conditions.names)-1;
icondition=[];isnewcondition=[];for ncondition=1:nconditions,[icondition(ncondition),isnewcondition(ncondition)]=conn_conditionnames(CONN_x.Setup.conditions.names{ncondition}); end
isvalidcondition=~isnewcondition;
for n1=numel(outcomeisource)
    isvalidcondition=isvalidcondition&conn_existfile(arrayfun(@(n)fullfile(CONN_x.folders.firstlevel_vv,CONN_x.vvAnalyses(CONN_x.vvAnalysis).name,['BETA_Subject',num2str(1,'%03d'),'_Condition',num2str(n,'%03d'),'_Measure',num2str(outcomeisource(n1),'%03d'),'_Component',num2str(outcomencompsource(n1),'%03d'),'.nii']),icondition,'uni',0)); 
end
validconditions=find(isvalidcondition);
isvalidsubject=conn_checkmissingdata(3,validconditions,validsources);
validsubjects=find(isvalidsubject);
iAllSubjects=find(ismember(CONN_x.Setup.l2covariates.names(1:end-1),{'AllSubjects'}),1);
%iAllSubjects=[]; % uncomment this line to use all subjects (instead of the AllSubjects covariate) for method = 1:5 (note: not 6) 
if ~isempty(iAllSubjects), AllSubjects=cell2mat(arrayfun(@(n)CONN_x.Setup.l2covariates.values{n}{iAllSubjects},reshape(validsubjects,[],1),'uni',0)); 
else AllSubjects=ones(numel(validsubjects),1); 
end
ncondition=validconditions(1);

V0=conn_vol(fullfile(filepath,['vvPC_Subject',num2str(validsubjects(1),'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'.mat']));
iV0mat=inv(V0.matdim.mat);
filename_B1=[]; % vvPC filenames
filename_S1={}; % MVPA score filenames
filename_S1vol=[];
filename_S1imat=[];
Nt=[];
xybak=[];
V1=[];          % vvPC volumes
X={};           % vvPC at the selected seed voxel/condition for each subject
Y={};           % vvPC at the selected target slice/condition for each subject
Z=[];           % MVPA scores at selected seed voxel/condition for each subject
thrZ=[];        % MVPA scores threshold (for high- vs. low- scoring subjects)
IDX=[];         % vvPC target voxel indices
XYZ=[0 50 28];  % seed voxel coordinates
S=[];           % structural at target slice
method=1;       % method (1:all subjects average, 2:low-score, 3:high-score, 4:high-low, 5:custom)
neig=1;         % eigenpattern number
W=[];           % between-subjects contrast
Wcustom=[];
covselected=[];
conselected=[];
txtmethod='';
nslice=round(iV0mat(3,:)*[XYZ(:);1]);
Sslice=round(imat(3,:)*[XYZ(:);1]);
dataview=3;
Si1=[];
Si2=[];
[filename_B1,V1,Nt]=conn_mvpaexplore_getinfo(0,filepath,icondition(ncondition),validsubjects);

%filename{nsub}=fullfile(filepathresults,['BETA_Subject',num2str(nsub,'%03d'),'_Condition',num2str(CONN_h.menus.m_results.icondition(nconditions(ncondition)),'%03d'),'_Measure',num2str(CONN_h.menus.m_results.outcomeisource(nsources(nsource)),'%03d'),'_Component',num2str(CONN_h.menus.m_results.outcomencompsource(nsources(nsource)),'%03d'),'.nii']);

% frame
boffset=[0 0 0 0];
conn_menu('frame',boffset+[.02,.15,.65,.74],'');
%conn_menu('frame2',boffset+[.045,.08,.91,.80],'');
% image left
posimage=[.08,.30,.16,.34];
ht2=conn_menu('image',boffset+posimage,'','','',@conn_mvpaexplore_mtncallback,@conn_mvpaexplore_click,'');
ht2title=conn_menu('pushbutton',boffset+[posimage(1),posimage(2)+posimage(4),posimage(3),.045],'','Seed','',@(varargin)conn_mvpaexplore_update('coordinates'));
set(ht2.h5b,'callback',@(varargin)conn_mvpaexplore_update('changeseedview')); % + symbol change view
set(ht2.h5,'callback',@(varargin)conn_mvpaexplore_update('changeseedslice'),'value',1); % slider change slice
hold(ht2.h1,'on'); ht2dot=plot(0,0,'color','k','marker','o','markerfacecolor','r','markersize',8+CONN_gui.font_offset,'parent',ht2.h1); hold(ht2.h1,'off'); 
try, addlistener(ht2.h5, 'ContinuousValueChange',@(varargin)conn_mvpaexplore_update('changeseedslice')); end
% image right
posimage=[.36,.20,.28,.54];
[ht4,ht4title]=conn_menu('image',boffset+posimage,' ','','','','','');
%uicontrol('style','text','units','norm','position',boffset+[posimage(1)+posimage(3)/2-.070,posimage(2)-1*.059,.070,.045],'string','threshold','fontname','default','fontsize',8+CONN_gui.font_offset,'backgroundcolor',CONN_gui.backgroundcolor,'foregroundcolor',CONN_gui.fontcolorA); 
ht4slider=conn_menu('slider',boffset+[posimage(1)+posimage(3),posimage(2),.015,posimage(4)],'','','z-slice',@(varargin)conn_mvpaexplore_update('slice'));
%try, addlistener(ht4slider, 'ContinuousValueChange',@(varargin)conn_mvpaexplore_update('slice')); end
set(ht4slider,'visible','off');
conn_menumanager('onregion',ht4slider,1,boffset+posimage+[0 0 .015 0]);
set(ht4slider,'min',1,'max',V0.matdim.dim(3),'sliderstep',min(.5,[1,10]/(V0.matdim.dim(3)-1)),'value',nslice);
%ht4slice=uicontrol('style','text','units','norm','position',boffset+[posimage(1)+posimage(3)/2-.059/2,posimage(2)-1*.045,.059,.045],'string',sprintf('z = %d mm',round(V0.matdim.mat(3,:)*[0;0;nslice;1])),'fontname','default','fontsize',8+CONN_gui.font_offset,'backgroundcolor',CONN_gui.backgroundcolor,'foregroundcolor',CONN_gui.fontcolorA); 
ht4slice=conn_menu('text',boffset+[posimage(1)+posimage(3)/2-.059/2,posimage(2)-1*.045,.059,.045],'',sprintf('z = %d mm',round(V0.matdim.mat(3,:)*[0;0;nslice;1])));
%ht4title=conn_menu('text2',boffset+[posimage(1),posimage(2)+posimage(4),posimage(3),.04],'',txtmethod);
% lists
[ht7, ht7title]=conn_menu('listbox2',boffset+[.86,.45,.09,.12],'fc-MVPA scores',arrayfun(@(n)sprintf('Eigenpattern #%d scores',n),1:numel(CONN_x.vvAnalyses(CONN_x.vvAnalysis).measures),'uni',0),'Select eigenpattern scores for display',@(varargin)conn_mvpaexplore_update('scores'));
ht5=conn_menu('listbox2',boffset+[.70,.45,.15,.12],'Subjects',{'All Subjects','High-scoring subjects','Low-scoring subjects','Low- & High- scoring subjects','High vs. Low between-subjects contrast','Custom groups or between-subjects contrast'},'<HTML>Select display option<br/> - Select <i><b>AllSubjects</b></i> to compute average connectivity across all subjects (note: the display will exclude subjects not in the <i>AllSubjects</i> second-level covariate group)<br/> - Select <i><b>low/high scoring subjects</b></i> to compute average connectivity separately within those subjects with low- vs. high- eigenvariate scores, characterizing <br/> principal axes of diversity in FC across subjects (for the selected eigenvariate) (note: the display will exclude subjects not in the <i>AllSubjects</i> group)<br/> - Select <i><b>custom group or between-subjects contrast</b></i> to define your own contrast across subjects (e.g. a different group, a between-group comparison, etc.)</HTML>',@(varargin)conn_mvpaexplore_update('subjects'));
ht6=conn_menu('listbox2',boffset+[.70,.26,.15,.12],'Conditions',CONN_x.Setup.conditions.names(validconditions),'<HTML>Select condition(s) for display</HTML>',@(varargin)conn_mvpaexplore_update('conditions'));
set([ht7,ht7title],'visible','off');
set([ht5],'max',2);
% buttons
ht11=conn_menu('pushbuttonblue2',boffset+[.70,.20,.07,.045],'','display 3D','creates whole-brain 3d view of current functional connectivity display',@(varargin)conn_mvpaexplore_update('display3d'));
ht13=conn_menu('pushbuttonblue2',boffset+[.77,.20,.07,.045],'','export data','<HTML>creates whole-brain NIFTI volumes with seed-to-voxel maps for each subject for the selected seed voxel</HTML>',@(varargin)conn_mvpaexplore_update('exportdata'));
ht12=conn_menu('pushbuttonblue2',boffset+[.84,.20,.07,.045],'','import values','<HTML>import fc-MVPA eigenpattern score values at selected seed for each subject as 2nd-level covariates</HTML>',@(varargin)conn_mvpaexplore_update('importvalues'));
set(ht12,'visible','off');
% hist
posimage=[.725,.73,.20,.10];
[ht21,ht22]=conn_menu('hist',boffset+posimage,'');
ht21title=conn_menu('text2',boffset+[posimage(1),posimage(2)-.07,posimage(3),.04],'','eigenpattern scores');
ht24=conn_menu('edit2',boffset+[posimage(1)+posimage(3)/2,posimage(2)+posimage(4)+.05,.06,.04],'','','<HTML>Select eigenpattern scores threshold dividing low- and high- scoring subjects<br/> - leave empty to specify the sample median (default)</HTML>',@(varargin)conn_mvpaexplore_update('threshold'));
ht25=uicontrol('style','frame','units','norm','position',boffset+[posimage(1)+posimage(3)/2-.01,posimage(2)+posimage(4)+.05-.01,.06+.02,.04+.02],'foregroundcolor',CONN_gui.backgroundcolor,'backgroundcolor',CONN_gui.backgroundcolor,'parent',CONN_h.screen.hfig);
set(ht25,'visible','on'); conn_menumanager('onregion',ht25,-1,boffset+[posimage(1)-.01,posimage(2)-.07,posimage(3)+.02,posimage(4)+.17]);
ht23=uicontrol('style','frame','units','norm','position',boffset+[posimage(1)-.01,posimage(2)-.07,posimage(3)+.02,posimage(4)+.12],'foregroundcolor',CONN_gui.backgroundcolor,'backgroundcolor',CONN_gui.backgroundcolor,'parent',CONN_h.screen.hfig);

conn_menu('updateimage',ht2,volref);
conn_menu('updateslider1',ht2,Sslice);
set(ht2.h5,'value',Sslice);
set(ht2.h5,'sliderstep',min(1,2*[1,10]/(get(ht2.h5,'max')-1)));

conn_mvpaexplore_click(XYZ);
conn_mvpaexplore_update refresh;

set(ht4.h9,'string','.20');
conn_menu('updatecscale',[],[],ht4.h9);
set(ht4.h10,'string','0');
conn_menu('updatethr',[],[],ht4.h10);

fh=@conn_mvpaexplore_update;

    function conn_mvpaexplore_click(pos,varargin)
        pos=reshape(pos,[],1);
        XYZ=pos(1:3);
        txyz=imat*[XYZ;1];
        set(ht2title,'string',sprintf('Seed (%d,%d,%d) mm',round(pos(1)),round(pos(2)),round(pos(3))));
        %data=get(ht2.h2,'userdata');
        switch(dataview)
            case 1, set(ht2dot,'xdata',volref.dim(1)+1-txyz(1),'ydata',volref.dim(3)+1-txyz(3));
            case 2, set(ht2dot,'xdata',volref.dim(2)+1-txyz(2),'ydata',volref.dim(3)+1-txyz(3));
            case 3, set(ht2dot,'xdata',volref.dim(1)+1-txyz(1),'ydata',volref.dim(2)+1-txyz(2));
        end
        conn_mvpaexplore_update seed;
    end

    function conn_mvpaexplore_update(option,varargin)
        switch(option)
            case 'importvalues'
                if ~isempty(Z)
                    Zallsubjects=nan(CONN_x.Setup.nsubjects,size(Z,2)); Zallsubjects(validsubjects,:)=Z;
                    conn_importl2covariate({sprintf('fc-MVPA: seed (%d,%d,%d) eigenpattern #%d scores during %s',XYZ(1),XYZ(2),XYZ(3),neig,CONN_x.Setup.conditions.names{ncondition})},{Zallsubjects},true,validsubjects,{sprintf('Imported eigenpattern scores from analysis %s',CONN_x.vvAnalyses(CONN_x.vvAnalysis).name)});
                end
                return
            case 'display3d'
                if size(W,1)>1, conn_msgbox('Sorry, only single-image display available. Please select in the ''subjects'' list an option displaying a single image','',2); 
                else
                    hm=conn_msgbox('Creating whole-brain data for 3D display, please wait','',-1);
                    for n1=1:size(W,1)
                        filenameout=fullfile(filepathresults,sprintf('connGUIimage%d.nii',n1));
                        w=zeros(V0.matdim.dim);
                        nv=iV0mat*[XYZ;1];
                        w(max(1,min(V0.matdim.dim(1),round(nv(1)))),max(1,min(V0.matdim.dim(2),round(nv(2)))),max(1,min(V0.matdim.dim(3),round(nv(3)))))=1;
                        conn_process('vv2rr',w(:)','style','vv2rv','saveas',filenameout,'validsubjects',validsubjects,'contrastsubjects',W(n1,:),'validconditions',ncondition,'contrastconditions',1);
                        tfh=conn_mesh_display(filenameout);
                        tfh('visible','off');
                        try, tfh('colorbar','rescale',str2num(get(ht4.h9,'string'))*[-1 1]); end
                        tfh('colormap','bluewhitered');
                        tfh('brain',2);
                        tfh('mask','off');
                        tfh('colorbar','on', txtmethod);
                        tfh('visible','on');
                    end
                    if ishandle(hm), delete(hm); end
                end
                return
            case 'exportdata'
                if size(W,1)>1, conn_msgbox('Sorry, only single-image display available. Please select in the ''subjects'' list an option displaying a single image','',2); 
                else
                    hm=conn_msgbox('Creating whole-brain data, please wait','',-1);
                    for n1=1:size(W,1)
                        filenameout=fullfile(filepathresults,sprintf('connGUIimage%d.nii',n1));
                        w=zeros(V0.matdim.dim);
                        nv=iV0mat*[XYZ;1];
                        w(max(1,min(V0.matdim.dim(1),round(nv(1)))),max(1,min(V0.matdim.dim(2),round(nv(2)))),max(1,min(V0.matdim.dim(3),round(nv(3)))))=1;
                        conn_process('vv2rr',w(:)','style','vv2rv','saveas',filenameout,'validsubjects',validsubjects,'contrastsubjects',[],'validconditions',ncondition,'contrastconditions',1);
                    end
                    if ishandle(hm), delete(hm); end
                end
                return
            case 'changeseedslice'
                %conn_menu('updateslider1',ht2);
                datan=get(ht2.h5,'value'); %data=get(ht2.h2,'userdata'); datan=data.n;
                pos=imat*[XYZ;1];
                switch(dataview)
                    case 1, pos(2)=datan; XYZ(2)=round(volref.mat(2,:)*pos); % coronal
                    case 2, pos(1)=volref.dim(1)+1-datan; XYZ(1)=round(volref.mat(1,:)*pos); % sagittal
                    case 3, pos(3)=datan; XYZ(3)=round(volref.mat(3,:)*pos); % axial
                end
                %conn_mvpaexplore_click(XYZ);
                %return
                %switch(dataview)
                %    case 1, set(ht2dot,'xdata',volref.dim(1)+1-pos(1),'ydata',volref.dim(3)+1-pos(3));
                %    case 2, set(ht2dot,'xdata',volref.dim(2)+1-pos(2),'ydata',volref.dim(3)+1-pos(3));
                %    case 3, set(ht2dot,'xdata',volref.dim(1)+1-pos(1),'ydata',volref.dim(2)+1-pos(2));
                %end
                set(ht2title,'string',sprintf('Seed (%d,%d,%d) mm',round(XYZ(1)),round(XYZ(2)),round(XYZ(3))));
                X={};
                Z=[];
                
            case {'changeseedview','coordinates'}
                pos=[];
                if strcmp(option,'changeseedview'), 
                    conn_menu('updateview',ht2); 
                    pos=XYZ;
                    data=get(ht2.h2,'userdata');
                    set(ht2.h5,'sliderstep',min(1,2*[1,10]/(data.x1.rdim(data.view)-1)));
                    dataview=data.view;
                else
                    answ=conn_menu_inputdlg({'seed XYZ coordinates (mm)'},'',1,{mat2str(XYZ(:)')},struct('Resize','on'));
                    if numel(answ)==1&&~isempty(answ{1}),
                        answ=str2num(answ{1});
                        if numel(answ)==3,
                            pos=reshape(answ(1:3),3,1);
                        end
                    end
                end
                if ~isempty(pos)
                    %data=get(ht2.h2,'userdata');
                    switch(dataview)
                        case 1, kslice=max(1,min(volref.dim(2), round([0 1 0 0]*imat*[pos;1]) )); % coronal
                        case 2, kslice=max(1,min(volref.dim(1), volref.dim(3)+1-round([1 0 0 0]*imat*[pos;1]) )); % sagittal
                        case 3, kslice=max(1,min(volref.dim(3), round([0 0 1 0]*imat*[pos;1]) )); % axial
                    end
                    conn_menu('updateslider1',ht2,kslice);
                    set(ht2.h5,'value',kslice);
                    conn_mvpaexplore_click(pos(:));
                end
                return
            case 'subjects'
                method=get(ht5,'value');
                if all(method==1|method==6), 
                    set([ht7,ht7title,ht12],'visible','off'); 
                    set(ht23,'visible','on');
                else 
                    set([ht7,ht7title,ht12],'visible','on'); 
                    set(ht23,'visible','off');
                end
                if any(method==6)
                    covselected=listdlg('liststring',CONN_x.Setup.l2covariates.names(1:end-1),'selectionmode','multiple','initialvalue',covselected,'promptstring','Select subject-effects:','ListSize',[300 200]);
                    if isempty(covselected), Wcustom=zeros(1,numel(validsubjects));
                    else
                        if numel(covselected)==1, conselected=1;
                        else
                            while 1
                                if isempty(conselected), conselected=eye(numel(covselected)); end
                                answ=conn_menu_inputdlg(sprintf('Between-subjects contrast (vector with %d values, or matrix with %d columns)',numel(covselected),numel(covselected)),'',1,mat2str(conselected));
                                if isempty(answ), return; end
                                conselected=str2num(answ{1});
                                if size(conselected,2)==numel(covselected), break; end
                                conselected=[];
                            end
                        end
                    end
                    Xcustom=cell2mat(arrayfun(@(n)cell2mat(CONN_x.Setup.l2covariates.values{n}(covselected)),reshape(validsubjects,[],1),'uni',0));
                    Icustom=~(all(Xcustom==0,2)|any(isnan(Xcustom),2));
                    Xcustom(~Icustom,:)=[];
                    Wcustom=zeros(size(conselected,1),numel(validsubjects));
                    Wcustom(:,Icustom')=conselected*inv(Xcustom'*Xcustom)*Xcustom';
                end                
                txtmethod='';
            case 'scores'
                neig=get(ht7,'value');
                txtmethod='';
                filename_S1={};
            case 'seed'
                X={};
                Z=[];
            case 'slice'
                nslice=max(1,min(V0.matdim.dim(3), round(get(ht4slider,'value')) ));
                set(ht4slice,'string',sprintf('z = %d mm',round(V0.matdim.mat(3,:)*[0;0;nslice;1])));
                Y={};S={};
            case 'conditions'
                ncondition=validconditions(get(ht6,'value'));
                [filename_B1,V1,Nt]=conn_mvpaexplore_getinfo(0,filepath,icondition(ncondition),validsubjects);
%                 tstr=num2str(icondition(ncondition),'%03d');
%                 filename_B1=arrayfun(@(nsub)fullfile(filepath,['vvPC_Subject',num2str(nsub,'%03d'),'_Condition',tstr,'.mat']), validsubjects,'uni',0);
                nv=iV0mat*[XYZ;1];
                set(CONN_h.screen.hfig,'pointer',CONN_gui.waiticon);drawnow
                [X,Y,IDX]=conn_mvpaexplore_getinfo(1,filename_B1,V1,Nt,nslice,nv);
                filename_S1={};
            case 'threshold'
                if ~isempty(Z)
                    thrZ=regexprep(get(ht24,'string'),'[^\+\-\.\d\(\)]','');
                    plotn=false; 
                    if ~isempty(thrZ), thrZ=str2num(thrZ); plotn=true; end
                    if isempty(thrZ)||isnan(thrZ), thrZ=median(Z(:)); set(ht24,'string',sprintf('thr = %.2f',thrZ)); end
                    [Z_pdf,Z_range]=conn_menu_plothist(Z(:)-thrZ,.25);
                    Z_rangeLow=Z_range(Z_range<=0);
                    Z_pdfLow=Z_pdf(Z_range<=0);
                    Z_rangeHigh=Z_range(Z_range>0);
                    Z_pdfHigh=Z_pdf(Z_range>0);
                    conn_menu('updatehist',ht21,{thrZ+[Z_rangeLow(1),Z_rangeLow',Z_rangeLow(end),0,Z_rangeHigh(1),Z_rangeHigh',Z_rangeHigh(end)],[0,0*Z_pdfLow',0,0,0,Z_pdfHigh',0],[0,Z_pdfLow',0,0,0,0*Z_pdfHigh',0]});
                    if plotn,
                        set(ht21.h6,'string',{'Low-scoring subjects ',sprintf('N=%d ',nnz(Z<=thrZ))},'horizontalalignment','right');
                        set(ht21.h7,'string',{'High-scoring subjects',sprintf('N=%d',nnz(Z>thrZ))},'horizontalalignment','left');
                    else
                        set(ht21.h6,'string','Low-scoring subjects ','horizontalalignment','right');
                        set(ht21.h7,'string','High-scoring subjects','horizontalalignment','left');
                    end
                    set(ht21.h2,'visible','off');
                end
        end
        if isempty(Y)||(isempty(filename_S1)&&any(ismember(method,[2,3,4,5])))
            set(CONN_h.screen.hfig,'pointer',CONN_gui.waiticon);drawnow
        end
        if isempty(X)
            nv=iV0mat*[XYZ;1];
            [X]=conn_mvpaexplore_getinfo(2,filename_B1,V1,Nt,nv);
        end
        if isempty(Y)
            [Y,IDX]=conn_mvpaexplore_getinfo(3,filename_B1,V1,Nt,nslice);
        end
        involrefspace=true;
        if isempty(S)
            if ~involrefspace % in vvPC resolution
                [tx,ty]=ndgrid(1:V0.matdim.dim(1),1:V0.matdim.dim(2));
                txyz=V0.matdim.mat*[tx(:) ty(:) nslice+zeros(numel(tx),1) ones(numel(tx),1)]';
                S=reshape(spm_get_data(volref,imat*txyz),V0.matdim.dim(1:2));
            else % in volref resolution
                [tx,ty]=ndgrid(1:volref.dim(1),1:volref.dim(2));
                Sslice=round([0 0 1 0]*imat*V0.matdim.mat*[0;0;nslice;1]);
                txyz=[tx(:) ty(:) Sslice+zeros(numel(tx),1) ones(numel(tx),1)]';
                S=reshape(spm_get_data(volref,txyz),volref.dim(1:2));
                txyz=iV0mat*volref.mat*[tx(:) ty(:) Sslice+zeros(numel(tx),1) ones(numel(tx),1)]';
                Si1=reshape(txyz(1,:),volref.dim(1:2));
                Si2=reshape(txyz(2,:),volref.dim(1:2));
            end
        end
        if isempty(txtmethod)
            if numel(method)>1, txtmethod='Functional connectivity with seed (r)';
            elseif method==1,   txtmethod='Functional connectivity with seed (r)';
            elseif method==2,   txtmethod=sprintf('Connectivity with seed (r) in subjects with high eigenpattern #%d scores',neig);
            elseif method==3,   txtmethod=sprintf('Connectivity with seed (r) in subjects with low eigenpattern #%d scores',neig);
            elseif method==4,   txtmethod=sprintf('Connectivity (r) in subjects with low (left) and high (right) eigenpattern #%d scores',neig);
            elseif method==5,   txtmethod='Difference in connectivity between high- and low-scoring subjects';
            elseif method==6,   txtmethod='Connectivity with seed - custom contrast';
            end
            set(ht4title,'string',txtmethod);
            set(ht21title,'string',sprintf('eigenpattern #%d scores',neig));
        end
        xy=0;
        W=[];
        for nmethod=1:numel(method)
            switch(method(nmethod))
                case 1, % average across AllSubjects
                    w=zeros(1,numel(validsubjects));
                    w(AllSubjects>0)=1/nnz(AllSubjects>0);
                    %w=ones(1,numel(validsubjects))/numel(validsubjects);
                case {2,3,4,5}, % low/high-score subjects
                    if isempty(filename_S1)
                        tstr=num2str(icondition(ncondition),'%03d');
                        for isub=1:numel(validsubjects)
                            nsub=validsubjects(isub);
                            filename_S1{isub}=fullfile(CONN_x.folders.firstlevel_vv,CONN_x.vvAnalyses(CONN_x.vvAnalysis).name,['BETA_Subject',num2str(nsub,'%03d'),'_Condition',tstr,'_Measure',num2str(outcomeisource(neig),'%03d'),'_Component',num2str(outcomencompsource(neig),'%03d'),'.nii']);
                        end
                        filename_S1vol=conn_fileutils('spm_vol',char(filename_S1));
                        filename_S1imat=inv(filename_S1vol(1).mat);
                        Z=[];
                        xybak=[];
                    end
                    if isempty(Z)
                        Z=conn_fileutils('spm_get_data',filename_S1vol,filename_S1imat*[XYZ;1]);
                        Z(~(AllSubjects>0))=nan;
                        Z=Z/max(eps,std(Z(~isnan(Z))));
                        thrZ=median(Z(~isnan(Z)));
                        set(ht24,'string',sprintf('thr = %.2f',thrZ));
                        [Z_pdf,Z_range]=conn_menu_plothist(Z(~isnan(Z))-thrZ,.25);
                        Z_rangeLow=Z_range(Z_range<=0);
                        Z_pdfLow=Z_pdf(Z_range<=0);
                        Z_rangeHigh=Z_range(Z_range>0);
                        Z_pdfHigh=Z_pdf(Z_range>0);
                        conn_menu('updatehist',ht21,{thrZ+[Z_rangeLow(1),Z_rangeLow',Z_rangeLow(end),0,Z_rangeHigh(1),Z_rangeHigh',Z_rangeHigh(end)],[0,0*Z_pdfLow',0,0,0,Z_pdfHigh',0],[0,Z_pdfLow',0,0,0,0*Z_pdfHigh',0]});
                        sz=sort(Z(:));set(ht21.h8,'xdata',sz,'ydata',mod(sz,.10*interp1(Z_range,Z_pdf,sz-thrZ)),'visible','on');
                        set(ht21.h6,'string','Low-scoring subjects ','horizontalalignment','right');
                        set(ht21.h7,'string','High-scoring subjects','horizontalalignment','left');
                        set(ht21.h2,'visible','off');
                    end
                    if method(nmethod)==2,      w=reshape((Z>thrZ)/max(eps,nnz(Z>thrZ)),1,[]); % high
                    elseif method(nmethod)==3,  w=reshape((Z<=thrZ)/max(eps,nnz(Z<=thrZ)),1,[]); % low
                    elseif method(nmethod)==4,  w=[reshape((Z<=thrZ)/max(eps,nnz(Z<=thrZ)),1,[]);reshape((Z>thrZ)/max(eps,nnz(Z>thrZ)),1,[])]; % low and high
                    else w=reshape((Z>thrZ)/max(eps,nnz(Z>thrZ))-(Z<=thrZ)/max(eps,nnz(Z<=thrZ)),1,[]); % high - low
                    end
                    %xy=xy/numel(validsubjects);
                    %if ~isempty(xybak)&&xybak(:)'*xy(:)<0, xy=-xy; end
                case 6
                    w=Wcustom;
            end
            W=cat(1,W,w);
        end
        xybak=xy;
        if nnz(W)==0
            conn_menu('update',ht4,{permute(S,[2,1,4,3]),[],[]},{V0.matdim,nslice});
        else
            xy=conn_mvpaexplore_getinfo(4,X,Y,W);
            if size(xy,1)==1
                xyMap=zeros(V0.matdim.dim(1:2));
                xyMask=zeros(V0.matdim.dim(1:2));
                xyMap(IDX)=xy;
                xyMask(IDX)=1;
            else
                xyMap=zeros([V0.matdim.dim(1:2),size(xy,1)]);
                xyMask=zeros([V0.matdim.dim(1:2),size(xy,1)]);
                for n1=1:size(xy,1), xyMap(prod(V0.matdim.dim(1:2))*(n1-1)+IDX)=xy(n1,:); end
                for n1=1:size(xy,1), xyMask(prod(V0.matdim.dim(1:2))*(n1-1)+IDX)=1; end
            end
            hfilt=conn_hanning(5); hfilt=hfilt/sum(hfilt); xyMap=xyMap.*xyMask+convn(convn(xyMap.*xyMask,hfilt,'same'),hfilt','same')./max(.5,convn(convn(xyMask,hfilt,'same'),hfilt','same')).*(1-xyMask);
            if ~involrefspace
                conn_menu('update',ht4,{permute(S,[2,1,4,3]),permute(xyMap,[2,1,4,3]),1+0*permute(xyMask,[2,1,4,3])},{V0.matdim,nslice});
            else
                for n1=1:size(xyMap,3), txy=interp2(xyMap(:,:,n1),Si2,Si1); if n1==1, NEWxyMap=zeros([size(txy,1),size(txy,2),size(xyMap,3)]); end; NEWxyMap(:,:,n1)=txy; end
                for n1=1:size(xyMask,3), txy=interp2(xyMask(:,:,n1),Si2,Si1); if n1==1, NEWxyMask=zeros([size(txy,1),size(txy,2),size(xyMask,3)]); end; NEWxyMask(:,:,n1)=txy; end
                NEWS=repmat(S,[1,1,size(xyMap,3)]);
                conn_menu('update',ht4,{permute(NEWS,[2,1,4,3]),permute(NEWxyMap,[2,1,4,3]),1+0*permute(NEWxyMask,[2,1,4,3])},{volref,Sslice});
            end
        end
        set(ht4.h10,'visible','off');

        try, if isfield(CONN_h,'menus')&&isfield(CONN_h.menus,'waiticonObj'), CONN_h.menus.waiticonObj.stop; end; end
        set(CONN_h.screen.hfig,'pointer','arrow');
        %drawnow;
        try, drawnow nocallbacks;
        catch, drawnow;
        end
    end

    function [str0,str]=conn_mvpaexplore_mtncallback(varargin)
        if CONN_gui.iscursordown, conn_mvpaexplore_click(varargin{:}); end
        str0={};
        str={};
    end
        

end

