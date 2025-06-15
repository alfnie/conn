function [dataplot,infoplot,data1plot]=conn_vproject(param,nonparam,views,projection,thres,side,parametric,res,box,select,threshold,data1plot,spmfile,voxeltovoxel,issurface)
% internal function voxel-level results display
%
global CONN_gui CONN_x;
if isempty(CONN_gui)||~isfield(CONN_gui,'font_offset'), conn_font_init; end
if isempty(CONN_gui)||~isfield(CONN_gui,'isremote'), CONN_gui.isremote=false; end
%CONN_VPROJECT Volume display
%

if numel(param)==1 && ishandle(param), % callbacks from UI objects
    init=0;
    doredraw=false;
    if ishandle(param)&&isequal(get(param,'type'),'figure'), GCF=param;
    else GCF=gcbf; 
    end
    if isempty(GCF), GCF=gcf; end
    if isstruct(nonparam), ARG=nonparam; end;
	OPTION=views; 
    if nargin>3, OPTION2=projection; else OPTION2=''; end
    if nargin>4, OPTION3=thres; else OPTION3=''; end
    if isempty(OPTION), return; end
    DATA=get(GCF,'userdata');
    a=DATA.a;
    b=DATA.b;
    c=DATA.c;
    d=DATA.d;
    e=DATA.e;
    f=DATA.f;
    g=DATA.g;
    SIM=DATA.SIM; 
    parametric=DATA.parametric;
    paramoptions=DATA.paramoptions;
    projection=DATA.projection;
    thres=DATA.thres;
    threshold=DATA.threshold;
    res=DATA.res;
    box=DATA.box;
    views=DATA.views;
    side=DATA.side;
    select=DATA.select;
    pointer=DATA.pointer;
    mat=DATA.mat;
    spmfile=DATA.spmfile;
    if isfield(DATA,'voxeltovoxel'), voxeltovoxel=DATA.voxeltovoxel; else voxeltovoxel=0; end
    if isfield(DATA,'issurface'), issurface=DATA.issurface; else issurface=0; end
    data1plot=[];
    if isstruct(OPTION),
    else,
        switch(lower(OPTION)),
            case 'buttondown',
                xlim=get(DATA.axes,'xlim');ylim=get(DATA.axes,'ylim');xpos=get(DATA.axes,'currentpoint');
                if xpos(1,1)>=xlim(1)&&xpos(1,1)<=xlim(2)&&xpos(1,2)>=ylim(1)&&xpos(1,2)<=ylim(2),
                    DATA.buttondown=get(0,'pointerlocation');
                    DATA.selectiontype=get(GCF,'selectiontype');
                    set(GCF,'windowbuttonmotionfcn',{@conn_vproject,'buttonmotion'},'busyaction','queue','userdata',DATA);
                    res=DATA.res*3/2;
                    box=1;
                else, DATA.buttondown=[]; set(GCF,'userdata',DATA); return; end
            case 'buttonup',
                set(GCF,'windowbuttonmotionfcn',[],'busyaction','cancel');
                p1=DATA.buttondown;
                if ~isempty(p1),
                    p2=get(0,'pointerlocation')-p1;
                    if strcmp(DATA.selectiontype,'extend'),
                        ang=-.005*(p2(2)+p2(1));projection([1,2],:)=[cos(ang),sin(ang);-sin(ang),cos(ang)]*projection([1,2],:);
                    else,
                        %ang=.01*p2(2);projection([2,3],:)=[cos(ang),sin(ang);-sin(ang),cos(ang)]*projection([2,3],:);
                        ang=-.01*p2(1);projection([1,3],:)=[cos(ang),sin(ang);-sin(ang),cos(ang)]*projection([1,3],:);
                    end
                    DATA.projection=projection;
                    set(GCF,'windowbuttonmotionfcn',[],'userdata',DATA);
                else, return; end
            case 'buttonmotion',
                p1=DATA.buttondown;
                if isempty(p1), return; end
                p2=get(0,'pointerlocation')-p1;
                if strcmp(DATA.selectiontype,'extend'),
                    ang=-.005*(p2(2)+p2(1));projection([1,2],:)=[cos(ang),sin(ang);-sin(ang),cos(ang)]*projection([1,2],:);
                else,
                    %ang=.01*p2(2);projection([2,3],:)=[cos(ang),sin(ang);-sin(ang),cos(ang)]*projection([2,3],:);
                    ang=-.01*p2(1);projection([1,3],:)=[cos(ang),sin(ang);-sin(ang),cos(ang)]*projection([1,3],:);
                end
                res=DATA.res*3/2;
                box=1;
            case 'resolution',
                res=str2num(get(DATA.handles(2),'string'));
                if isempty(res), res=DATA.res; end
                set(DATA.handles(2),'string',num2str(res));
                DATA.res=res;
                set(GCF,'userdata',DATA);
                init=-1;
            case 'threshold',
                threshold=str2num(get(DATA.handles(4),'string'));
                if isempty(threshold), threshold=DATA.threshold; end
                set(DATA.handles(4),'string',num2str(threshold));
                DATA.threshold=threshold;
                set(GCF,'userdata',DATA);
                init=-1;
            case 'keypress',
                if exist('ARG','var')&&isfield(ARG,'Character')
                    switch(lower(ARG.Character)),
                        case 'y',
                            DATA.rotation=DATA.rotation+[10,0,0]*(2*(length(ARG.Modifier)>0)-1);
                        case 'x',
                            DATA.rotation=DATA.rotation+[0,10,0]*(2*(length(ARG.Modifier)>0)-1);
                        case 'z',
                            DATA.rotation=DATA.rotation+[0,0,10]*(2*(length(ARG.Modifier)>0)-1);
                        otherwise,
                            return;
                    end
                    ang=-.01*DATA.rotation(1);projection([1,3],:)=[cos(ang),sin(ang);-sin(ang),cos(ang)]*projection([1,3],:);
                    ang=.01*DATA.rotation(2);projection([2,3],:)=[cos(ang),sin(ang);-sin(ang),cos(ang)]*projection([2,3],:);
                    ang=-.01*DATA.rotation(3);projection([1,2],:)=[cos(ang),sin(ang);-sin(ang),cos(ang)]*projection([1,2],:);
                    DATA.rotation=[0,0,0];
                    DATA.projection=projection;
                    set(GCF,'windowbuttonmotionfcn',[],'userdata',DATA);
                else
                    return;
                end
            case 'selectroi'
                if ~isempty(OPTION2), selectcluster=OPTION2;
                else selectcluster=get(DATA.handles(8),'value');
                end
                %if numel(selectcluster)>1, selectcluster=selectcluster(1); set(DATA.handles(8),'value',selectcluster); end
                if selectcluster>length(DATA.clusters), selectcluster=[]; end
                if isempty(selectcluster)||isequal(selectcluster,0)||any(selectcluster>length(DATA.clusters)), selectcluster=[]; select=[]; set(DATA.handles(8),'value',selectcluster); else, select=cat(1,DATA.clusters{selectcluster}); end
                DATA.select=select; DATA.selectcluster=selectcluster; set(GCF,'userdata',DATA);
                init=-2;
            case {'fwec.voxellevel','fwec.voxellevel.value','fwec.voxellevel.type'},
                if strcmp(OPTION,'fwec.voxellevel.value'), 
                    thres{1}=OPTION2;
                    set(DATA.handles(2),'string',num2str(thres{1})); 
                elseif strcmp(OPTION,'fwec.voxellevel.type'), 
                    thres{2}=OPTION2;
                    set(DATA.handles(3),'value',thres{2}); 
                elseif isempty(OPTION2)
                    temp=str2num(get(DATA.handles(2),'string')); 
                    if ~isempty(temp), thres{1}=temp; else, set(DATA.handles(2),'string',num2str(thres{1})); end
                    thres{2}=get(DATA.handles(3),'value');
                else
                    if iscell(OPTION2), thres(1:numel(OPTION2))=OPTION2;
                    else thres{1}=OPTION2;
                    end
                    set(DATA.handles(2),'string',num2str(thres{1}));
                    set(DATA.handles(3),'value',thres{2});
                end
                docont=true;

                if DATA.thres{2}==1&&thres{2}==5&&DATA.mat{6}=='F', thres{1}=spm_invFcdf(1-thres{1},DATA.mat{3}(1),DATA.mat{3}(2)); set(DATA.handles(2),'string',mat2str(thres{1})); docont=false;
                elseif DATA.thres{2}==1&&thres{2}==5&&DATA.side==3, thres{1}=spm_invTcdf(1-thres{1}/2,DATA.mat{3}(end)); set(DATA.handles(2),'string',mat2str(thres{1})); docont=false;
                elseif DATA.thres{2}==1&&thres{2}==5, thres{1}=spm_invTcdf(1-thres{1},DATA.mat{3}(end)); set(DATA.handles(2),'string',mat2str(thres{1})); docont=false;
                elseif DATA.thres{2}==5&&thres{2}==1&&DATA.mat{6}=='F', thres{1}=1-spm_Fcdf(thres{1},DATA.mat{3}(1),DATA.mat{3}(2)); set(DATA.handles(2),'string',mat2str(thres{1},4)); docont=false;
                elseif DATA.thres{2}==5&&thres{2}==1&&DATA.side==3,     thres{1}=1-spm_Tcdf(thres{1},DATA.mat{3}(end)); thres{1}=2*min(thres{1},1-thres{1}); set(DATA.handles(2),'string',mat2str(thres{1},4)); docont=false;
                elseif DATA.thres{2}==5&&thres{2}==1, thres{1}=1-spm_Tcdf(thres{1},DATA.mat{3}(end)); set(DATA.handles(2),'string',mat2str(thres{1},4)); docont=false;
                elseif (thres{2}==3||thres{2}==4)&&(thres{1}<0|thres{1}>1), thres{1}=max(0,min(1,thres{1})); set(DATA.handles(2),'string',num2str(thres{1}));
                elseif thres{2}<3&&(thres{1}<0|thres{1}>.5), thres{1}=max(0,min(.5,thres{1})); set(DATA.handles(2),'string',num2str(thres{1}));
                elseif thres{2}==5&&thres{1}<0, thres{1}=max(0,thres{1}); set(DATA.handles(2),'string',num2str(thres{1}));
                end
                if thres{2}<3, set(DATA.handles(1),'string','voxel threshold: p <'); 
                elseif (thres{2}==3||thres{2}==4), set(DATA.handles(1),'string','voxel threshold: p <'); 
                else set(DATA.handles(1),'string',['voxel threshold: ',DATA.mat{6},' >']); 
                end
                icp=[4 5 6  15]; if (thres{2}==3||thres{2}==4), set(DATA.handles(icp),'visible','off'); 
                else set(DATA.handles(icp),'visible','on'); 
                end
                DATA.thres=thres;
                %if ~docont, set(GCF,'userdata',DATA); return; end
                DATA.selectcluster=[];
                set(GCF,'userdata',DATA);
                init=-1;
                doredraw=true;
            case {'fwec.clusterlevel','fwec.clusterlevel.value','fwec.clusterlevel.type'},
                if strcmp(OPTION,'fwec.clusterlevel.value'), 
                    thres{3}=OPTION2;
                    set(DATA.handles(5),'string',num2str(thres{3})); 
                elseif strcmp(OPTION,'fwec.clusterlevel.type'), 
                    thres{4}=OPTION2;
                    thres4=thres{4}+(thres{4}>=6);
                    if DATA.parametric==1&&thres4>8, conn_msgbox('This threshold option is only available when using non-parametric statistics','',2); thres{4}=3;
                    elseif DATA.parametric==2&&ismember( thres4,[5 6 7]), conn_msgbox('This threshold option is only available when using parametric statistics','',2); thres{4}=3;
                    end
                    set(DATA.handles(6),'value',thres{4});
                elseif isempty(OPTION2), 
                    temp=str2num(get(DATA.handles(5),'string'));
                    if ~isempty(temp), thres{3}=temp; else, set(DATA.handles(5),'string',num2str(thres{3})); end
                    thres{4}=get(DATA.handles(6),'value');
                    thres4=thres{4}+(thres{4}>=6);
                    if DATA.parametric==1&&thres4>8, conn_msgbox('This threshold option is only available when using non-parametric statistics','',2); thres{4}=3;
                    elseif DATA.parametric==2&&ismember( thres4,[5 6 7]), conn_msgbox('This threshold option is only available when using parametric statistics','',2); thres{4}=3;
                    end
                    set(DATA.handles(6),'value',thres{4});
                else
                    if iscell(OPTION2), thres(2+(1:numel(OPTION2)))=OPTION2;
                    else thres{3}=OPTION2;
                    end
                    set(DATA.handles(5),'string',num2str(thres{3}));
                    set(DATA.handles(6),'value',thres{4});
                end
                DATA.selectcluster=[];
                if thres{4}==1||thres{4}==8, set(DATA.handles(4),'string','cluster threshold: k >'); 
                elseif DATA.parametric==1&&thres{4}>=5&&thres{4}<=6, set(DATA.handles(4),'string','peak threshold: p <'); 
                else set(DATA.handles(4),'string','cluster threshold: p <'); 
                end
                DATA.thres=thres;set(GCF,'userdata',DATA);
                init=-1;
            case {'export_mat'},
                if isempty(OPTION2)
                    [filename,filepath]=uiputfile('*.mat','Save structure as',fileparts(spmfile));
                else
                    [filepath,filename_name,filename_ext]=fileparts(OPTION2);
                    if isempty(filepath), filepath=fileparts(spmfile); end
                    if isempty(filepath), filepath=pwd; end
                    filename=[filename_name,filename_ext];
                end
                if ischar(filename),
                    data=rmfield(DATA,'handles');
                    data=rmfield(DATA,'axes');
                    conn_savematfile(fullfile(filepath,filename),'data');
                    conn_disp('fprintf','Results structure file saved as %s\n',fullfile(filepath,filename));
                end
            case {'export','export_mask','export_mask_selected'},
                if isempty(OPTION2)
                    %[filename,filepath]=uiputfile('*.nii;*.img','Save mask as',fileparts(spmfile));
                    [filename,filepath]=conn_fileutils('uiputfile','*.nii;*.img','Save mask as',fileparts(spmfile));
                else
                    [filepath,filename_name,filename_ext]=fileparts(OPTION2);
                    if isempty(filepath), filepath=fileparts(spmfile); end
                    if isempty(filepath), filepath=pwd; end
                    filename=[filename_name,filename_ext];
                end
                if ischar(filename),
                    [nill,filename_name,filename_ext]=fileparts(filename);
                    if isempty(filename_ext), filename=[filename,'.nii']; end
                    d0=cellstr(get(DATA.handles(11),'string'));d1=cellstr(get(DATA.handles(3),'string'));d2=cellstr(get(DATA.handles(6),'string'));
                    descrip1=sprintf('%s %s %s %s',get(DATA.handles(1),'string'),get(DATA.handles(2),'string'),d1{get(DATA.handles(3),'value')},d0{get(DATA.handles(11),'value')});
                    if (DATA.thres{2}==3||DATA.thres{2}==4), 
                        tempT=e; 
                        descrip2='none';
                        descrip3=sprintf('%s %s',get(DATA.handles(13),'string'),get(DATA.handles(14),'string'));
                    else 
                        tempT=d;
                        descrip2=sprintf('%s %s %s',get(DATA.handles(4),'string'),get(DATA.handles(5),'string'),d2{get(DATA.handles(6),'value')});
                        descrip3=sprintf('%s %s',get(DATA.handles(13),'string'),get(DATA.handles(14),'string'));
                    end
                    if strcmp(OPTION,'export_mask_selected')&&isfield(DATA,'select')&&~isempty(DATA.select), maskselected=zeros(size(tempT)); maskselected(DATA.select)=1;
                    else maskselected=ones(size(tempT));
                    end
                    descrip=sprintf('%s ; %s',descrip1,descrip2);                    
                    V=struct('mat',mat{1},'dim',size(b(:,:,:,1)),'dt',[spm_type('float32') spm_platform('bigend')],'fname',fullfile(filepath,filename),'descrip',descrip);
                    V=conn_fileutils('spm_write_vol',V,maskselected.*tempT.*double(b(:,:,:,end)>0));
                    try, conn_fileutils('spm_jsonwrite',conn_prepend('',fullfile(filepath,filename),'.json'),struct('description',sprintf('%s ; %s ; %s',descrip1,descrip2,descrip3)),struct('indent',' ')); end
                    conn_disp('fprintf','Results suprathreshold stats file (with voxel-level stats at each suprathreshold voxel) saved as %s\n',V.fname);
                    
                    try, 
                        conn_exportlist(DATA.handles(8),fullfile(filepath,[filename_name,'.Table.tsv']),get(DATA.handles(7),'string'),[],false,true); 
                        conn_disp('__nolog','fprintf','Results table saved as %s\n',fullfile(filepath,[filename_name,'.Table.tsv']));
                    end
%                         %fh=fopen(fullfile(filepath,[filename_name,'.Table.txt']),'at');
%                         str={};
%                         str{end+1}=sprintf('\n\n');
%                         for nbtemp=1:numel(DATA.clusters),
%                             str{end+1}=sprintf('\nCluster ');
%                             str{end+1}=sprintf('%+d ',[DATA.xyzpeak{nbtemp},1]*DATA.mat{1}(1:3,:)');
%                             str{end+1}=sprintf(' :\n');
%                             tstr=DATA.clusternames{nbtemp}.txt;
%                             for n=1:numel(tstr),str{end+1}=sprintf('%s\n',tstr{n});end
%                         end
%                         str{end+1}=sprintf('\nAll clusters combined :\n');
%                         tstr=DATA.clusternames{end}.txt;
%                         for n=1:numel(tstr),str{end+1}=sprintf('%s\n',tstr{n});end
%                         %fclose(fh);
%                         conn_fileutils('fileappend_raw',fullfile(filepath,[filename_name,'.Table.txt']),str);
%                     end

                    V=struct('mat',mat{1},'dim',size(b(:,:,:,1)),'pinfo',[1;0;0],'dt',[spm_type('uint32') spm_platform('bigend')],'fname',fullfile(filepath,[filename_name,'.ROIs.nii']),'descrip',descrip);
                    btemp=zeros(size(b(:,:,:,1)));
                    for nbtemp=1:numel(DATA.clusters), btemp(DATA.clusters{nbtemp})=nbtemp; end
                    V=conn_fileutils('spm_write_vol',V,btemp);
                    str={};
                    %fh=fopen(fullfile(filepath,[filename_name,'.ROIs.txt']),'wt');
                    for nbtemp=1:numel(DATA.clusters),
                        str{end+1}=[sprintf('Cluster_%d ',nbtemp), sprintf('%+03d ',[DATA.xyzpeak{nbtemp},1]*DATA.mat{1}(1:3,:)')];
                        str{end+1}=sprintf('\n');
                    end
                    %fclose(fh);
                    conn_fileutils('filewrite_raw',fullfile(filepath,[filename_name,'.ROIs.txt']),str);
                    conn_disp('__nolog','fprintf','Results ROI file (defining one ROI for each suprathreshold cluster) saved as %s\n',V.fname);

                    V=struct('mat',mat{1},'dim',size(b(:,:,:,1)),'dt',[spm_type('float32') spm_platform('bigend')],'fname',fullfile(filepath,[filename_name,'.Mask.nii']),'descrip',descrip);
                    V=conn_fileutils('spm_write_vol',V,maskselected.*double(b(:,:,:,end)>0));
                    %conn_disp('fprintf','Mask file saved as %s\n',V.fname);
                    try, conn_fileutils('spm_jsonwrite',fullfile(filepath,[filename_name,'.json']),struct('HeightThreshold',descrip1,'SizeThreshold',descrip2,'Stats',descrip3)); end
                    conn_disp('__nolog','fprintf','Results suprathreshold mask file (with 1/0 values identifying suprathreshold voxels) saved as %s\n',V.fname);
                    
                    V=struct('mat',mat{1},'dim',size(b(:,:,:,1)),'dt',[spm_type('float32') spm_platform('bigend')],'fname',fullfile(filepath,[filename_name,'.nonthr.nii']),'descrip',descrip);
                    V=conn_fileutils('spm_write_vol',V,tempT);
                    conn_disp('__nolog','fprintf','Results unthresholded stats file (with voxel-level stats at all voxels) saved as %s\n',V.fname);
                    
                    V=struct('mat',mat{1},'dim',size(b(:,:,:,1)),'pinfo',[1;0;0],'dt',[spm_type('uint32') spm_platform('bigend')],'fname',fullfile(filepath,[filename_name,'.anat.nii']),'descrip',descrip);
                    btemp=zeros(size(b(:,:,:,1)));
                    mbtemp=0; sizebtemp=[]; 
                    for nbtemp=1:numel(DATA.clusters), 
                        %tempsizebtemp=accumarray(DATA.clusternames{nbtemp}.anat(:),1)'; % measure to sort: size
                        tempsizebtemp=accumarray(DATA.clusternames{nbtemp}.anat(:),abs(tempT(DATA.clusters{nbtemp}(:))),[],@max)'; % measure to sort: maximum voxel-level stats
                        sizebtemp=[sizebtemp, [nbtemp+0*tempsizebtemp; 1:numel(tempsizebtemp); tempsizebtemp]]; mbtemp=mbtemp+numel(tempsizebtemp); 
                    end
                    if isempty(sizebtemp), idxbtemp=[];
                    else [nill,idxbtemp]=sort(-sizebtemp(end,:)); rankbtemp=idxbtemp; rankbtemp(idxbtemp)=1:numel(idxbtemp);  % sort output regions by measure above
                    end
                    mbtemp=0; for nbtemp=1:numel(DATA.clusters), btemp(DATA.clusters{nbtemp})=rankbtemp(mbtemp+DATA.clusternames{nbtemp}.anat); mbtemp=mbtemp+sum(sizebtemp(1,:)==nbtemp); end
                    V=conn_fileutils('spm_write_vol',V,btemp);
                    str={};
                    %fh=fopen(fullfile(filepath,[filename_name,'.ROIs.txt']),'wt');
                    for ibtemp=reshape(idxbtemp,1,[])
                        str{end+1}=[sprintf('cluster %d/%d peak %.1f: ',sizebtemp(1,ibtemp),numel(DATA.clusters),sizebtemp(3,ibtemp)),DATA.clusternames{sizebtemp(1,ibtemp)}.txt{sizebtemp(2,ibtemp)}];
                        str{end+1}=sprintf('\n');
                    end
                    %fclose(fh);
                    conn_fileutils('filewrite_raw',fullfile(filepath,[filename_name,'.anat.txt']),str);
                    conn_disp('__nolog','fprintf','Results anatomical reference atlas file (conjunction of anatomical reference atlas and suprathreshold clusters) saved as %s\n',V.fname);

                    str={sprintf('ROI\tVoxels\tPeak statistics\tMNI coordinates (mm)\n')};
                    for ibtemp=reshape(idxbtemp,1,[])
                        tn=str2double(regexp(DATA.clusternames{sizebtemp(1,ibtemp)}.txt{sizebtemp(2,ibtemp)},'^\d+','match','once'));
                        tstr=regexprep(DATA.clusternames{sizebtemp(1,ibtemp)}.txt{sizebtemp(2,ibtemp)},'^\d+ voxels \(\d+%\) covering ','');
                        tstr=regexp(tstr,' with center at ','split');
                        str{end+1}=sprintf('%s\t%d\t%.1f\t%s\n',tstr{1},tn(1),sizebtemp(3,ibtemp),tstr{2});
                    end
                    conn_fileutils('filewrite_raw',fullfile(filepath,[filename_name,'.anat.tsv']),str);
                    conn_disp('__nolog','fprintf','Results anatomical reference table saved as %s\n',fullfile(filepath,[filename_name,'.anat.tsv']));

                    if (DATA.thres{2}==3||DATA.thres{2}==4), 
                        tidx=find(DATA.f&b(:,:,:,end)>0);
                        [tx,ty,tz]=ind2sub(size(b(:,:,:,1)), tidx);
                        txyz=mat{1}*[tx(:) ty(:) tz(:) ones(numel(tx),1)]';
                        [nill,tidx]=sort(-abs(tempT(tidx))); % sorts TFCE scores in in descending order
                        txyz=txyz(:,tidx);
                        conn_createmniroi(fullfile(filepath,[filename_name,'.PEAKs.nii']),txyz(1:3,:)',5);
                    end
                    

%                     [peaks,peaks_idx]=conn_peaks({tempT,mat{1}},12);
%                     peaks_suprathreshold=b(peaks_idx+size(b,1)*size(b,2)*size(b,3))>0;
%                     %peaks_F=d(peaks_idx);
%                     %peaks_p=exp(-b(peaks_idx(peaks_suprathreshold)));
%                     %dof=mat{3};
%                     save(fullfile(filepath,[filename_name,'.PEAKS.mat']),'peaks','peaks_suprathreshold');%,'peaks_F','peaks_p','dof');
                    %fprintf('Peaks file saved as %s\n',fullfile(filepath,[filename_name,'.PEAKS.mat']));
                    %if ~isempty(C)
                    %    V=struct('mat',mat{1},'dim',size(C),'pinfo',[1;0;0],'dt',[spm_type('uint32') spm_platform('bigend')],'fname',fullfile(filepath,[filename_name,'.SEGs.img']),'descrip',descrip);
                    %    V=conn_fileutils('spm_write_vol',V,C);
                    %    conn_disp('fprintf','Segmentation file saved as %s\n',V.fname);
                    %    fh=fopen(fullfile(filepath,[filename_name,'.SEGs.txt']),'wt');
                    %    for nbtemp=1:numel(peaks_idx),
                    %        fprintf(fh,'%+d ',peaks(nbtemp,:));
                    %        fprintf(fh,'\n');
                    %    end
                    %    fclose(fh);
                    %end
                end
                return;
                
            case 'plot_design',
                hmsg=conn_msgbox('Initializing. Please wait...','',-1);
                SPM=struct; conn_loadmatfile(spmfile,'SPM');
                if isfield(SPM,'xX_multivariate')
                    if 0,%isfield(SPM.xX_multivariate,'Zcontr')
                        conn_displaydesign(SPM.xX_multivariate.X,SPM.xX_multivariate.Zfiles,SPM.xX_multivariate.C,SPM.xX_multivariate.Zcontr,SPM.xX_multivariate.Xnames,false);
                    else
                        conn_displaydesign(SPM.xX_multivariate.X,reshape({SPM.xY.VY.fname},size(SPM.xX_multivariate.X,1),[]),SPM.xX_multivariate.C,SPM.xX_multivariate.M,SPM.xX_multivariate.Xnames,false);
                    end
                else
                    conn_displaydesign(SPM.xX.X,reshape({SPM.xY.VY.fname},size(SPM.xX.X,1),[]),SPM.xCon(1).c',1,SPM.xX.name,false);
                end
                if ishandle(hmsg), delete(hmsg); end
                return;
                
            case 'polar_view'
                if ~isempty(OPTION2), options={OPTION2,'-nogui'};
                else options={};
                end
                filename=fullfile(fileparts(spmfile),'results.nii');
                conn_vproject(GCF,[],'export_mask',filename);
                filename=fullfile(fileparts(spmfile),'results.ROIs.nii');
                conn_polar_display(filename,[],[],[],'Cluster (x,y,z)');
                return;
                
            case {'surface_view','surface_print','surface_gui_print','default_print'}
                if ~isempty(OPTION2), options={OPTION2,'-nogui'};
                else options={};
                end
                filename=fullfile(fileparts(spmfile),'results.nii');
                conn_vproject(GCF,[],'export_mask_selected',filename);
                fh=conn_mesh_display(filename,'');
                fh('visible','off');
                if (DATA.thres{2}==3||DATA.thres{2}==4), tstr='TFCE';
                elseif numel(DATA.mat{3})>1&&isequal(DATA.mat{6},'T'), tstr=[DATA.mat{6},'(',num2str(DATA.mat{3}(end)),')'];
                elseif numel(DATA.mat{3})>1, tstr=[DATA.mat{6},'(',num2str(DATA.mat{3}(1)),',',num2str(DATA.mat{3}(2)),')'];
                else tstr=[DATA.mat{6},'(',num2str(DATA.mat{3}),')'];
                end
                fh('colorbar','on', tstr);
                fh('visible','on');
                if ~isempty(regexp(OPTION,'print$')), fh('background',[1 1 1]); fh('print',4,options{:}); fh('close'); end
                %fh('zoomin');
                return;

            case {'volume_view','volume_print'}
                if ~isempty(OPTION2), options={OPTION2,'-nogui'};
                else options={};
                end
                filename=fullfile(fileparts(spmfile),'results.nii');
                conn_vproject(GCF,[],'export_mask_selected',filename);
                if DATA.issurface, vfilename=conn_surf_surf2vol(filename); 
                else vfilename=filename;
                end
                fh=conn_mesh_display('',vfilename); % without surface projection
                fh('visible','off');
                %fh=conn_mesh_display(filename,vfilename); % with surface projection as well
                if DATA.issurface, fh('brain',1); fh('sub_transparency',0); 
                else fh('sub_transparency',.5);
                end
                fh('brain_color',[1 1 1]);
                fh('sub_color',[1 1 1]);
                fh('brain_transparency',.75);
                fh('act_transparency',.75);
                fh('act_color',[1 0 0]); fh('colormap',[1 0 0]); 
                fh('visible','on');
                if ~isempty(regexp(OPTION,'print$')), fh('background',[1 1 1]); fh('print',7,options{:}); fh('close'); end
                return;

            case {'glass_view','glass_print'}
                if ~isempty(OPTION2), options={OPTION2,'-nogui'};
                else options={};
                end
                filename=fullfile(fileparts(spmfile),'results.nii');
                conn_vproject(GCF,[],'export_mask_selected',filename);
                if DATA.issurface, filename=conn_surf_surf2vol(filename); end
                fh=conn_mesh_display('',filename,[],[],[],[],.5);
                fh('visible','off');
                %fh('background',[1 1 1]);
                if DATA.issurface, 
                    fh('brain',2); 
                    fh('sub_transparency',0);
                    fh('brain_transparency',.05);
                    fh('axis','on');
                elseif 0
                    fh('brain',4);
                    fh('brain_transparency',0);
                    fh('sub_transparency',0);
                    fh('mask_transparency',.05);
                    fh('mask_color',[1 1 1]);
                    %fh('material',[.1 1 1 .25 0]);
                    fh('axis','on');
                else
                    fh('brain',2);
                    fh('ref_file',filename,2);
                    fh('ref_pos',[80,-110,90]);
                    fh('brain_transparency',.02);
                    fh('sub_transparency',.02);
                    if 1
                        fh('material',[.1 1 1 .25 0]);
                        fh('ref_all','on');
                        fh('ref_material','dull');
                        fh('act_transparency',.0);
                    else 
                        fh('material',[]);
                    end
                    fh('axis','on');
                    if (DATA.thres{2}==3||DATA.thres{2}==4), tstr='TFCE';
                    elseif numel(DATA.mat{3})>1&&isequal(DATA.mat{6},'T'), tstr=[DATA.mat{6},'(',num2str(DATA.mat{3}(end)),')'];
                    elseif numel(DATA.mat{3})>1, tstr=[DATA.mat{6},'(',num2str(DATA.mat{3}(1)),',',num2str(DATA.mat{3}(2)),')'];
                    else tstr=[DATA.mat{6},'(',num2str(DATA.mat{3}),')'];
                    end
                    fh('colorbar','on', tstr);
                end
                fh('visible','on');
                if ~isempty(regexp(OPTION,'print$')), fh('background',[1 1 1]); fh('print',3,options{:}); fh('close'); end
                return;
                
            case 'spm_view'
                [spmfile_path,spmfile_name]=fileparts(spmfile);
                cwd=pwd;conn_fileutils('cd',spmfile_path);
                SPM=struct; conn_loadmatfile(spmfile,'SPM');
                spm defaults fmri;
                [hReg,xSPM,SPM] = spm_results_ui('setup',SPM);
                assignin('base','hReg',hReg);
                assignin('base','xSPM',xSPM);
                assignin('base','SPM',SPM);
                figure(spm_figure('GetWin','Interactive'));
                %cd(cwd);
                return;
                
            case {'slice_view','slice_print'}
                if ~isempty(OPTION2), options={OPTION2,'-nogui'};
                else options={};
                end
                [spmfile_path,spmfile_name]=fileparts(spmfile);
                if (DATA.thres{2}==3||DATA.thres{2}==4), tempT=DATA.e; tempstats='TFCE'; tempP=[];
                else tempT=DATA.d; tempstats=DATA.mat{6}; tempP=DATA.c;
                end
                if isfield(DATA,'select')&&~isempty(DATA.select), maskselected=zeros(size(tempT)); maskselected(DATA.select)=1;
                else maskselected=ones(size(tempT));
                end
                fh=conn_slice_display(struct('T',tempT,'p',tempP,'stats',tempstats,'dof',DATA.mat{3},'mat',DATA.mat{1},'supra',maskselected.*abs(tempT).*(DATA.b(:,:,:,2)>0),'clusters',struct('stats',{DATA.cluster_stats},'idx',{DATA.clusters}),'filename',spmfile),...
                    [],...
                    spmfile_path);
                if (DATA.thres{2}==3||DATA.thres{2}==4), tstr='TFCE';
                elseif numel(DATA.mat{3})>1&&isequal(DATA.mat{6},'T'), tstr=[DATA.mat{6},'(',num2str(DATA.mat{3}(end)),')'];
                elseif numel(DATA.mat{3})>1, tstr=[DATA.mat{6},'(',num2str(DATA.mat{3}(1)),',',num2str(DATA.mat{3}(2)),')'];
                else tstr=[DATA.mat{6},'(',num2str(DATA.mat{3}),')'];
                end
                fh('colorbar','on', tstr);
                fh('contour_transparency',1);
                fh('slice_transparency',.5);
                fh('multisliceset',1,16,8);
                fh('pointer_mm',[0 0 10]);
                fh('togglegui',1);
                if ~isempty(regexp(OPTION,'print$')), 
                    fh('background',[1 1 1]); 
                    fh('print',options{:});
                    fh('close'); 
                end
                return
                
            case 'cluster_view_old'
                filename=fullfile(fileparts(spmfile),'results.nii');
                conn_vproject(GCF,[],'export_mask',filename);
                filename=fullfile(fileparts(spmfile),'results.ROIs.nii');
                [tfilename,tspmfile,tviewrex]=conn_vproject_selectfiles(filename,spmfile,0);
                if isempty(tfilename), return; end
                [tspmfile_path,tspmfile_name]=fileparts(tspmfile);
                cwd=pwd;
                conn_fileutils('cd',tspmfile_path);
                
                if tviewrex
                    conn_rex(tspmfile,tfilename,'output_type','saverex','level','clusters','select_clusters',0,'s',[],'gui',1);%,'mstats',false);
                else
                    htfig=conn_msgbox('Loading connectivity values. Please wait','',-1);
                    conn_rex(tspmfile,tfilename,'output_type','saverex','level','clusters','select_clusters',0,'steps',{'extract','results'},'s',[],'gui',0);%,'mstats',false);
                    if ishandle(htfig), delete(htfig); end
                end
                cd(cwd);
                return;
                
            case {'network_view','network_print'}
                if ~isempty(OPTION2), options={OPTION2,'-nogui'};
                else options={};
                end
                if ~isfield(CONN_x,'filename')||isempty(CONN_x.filename), % not conn project loaded
                    conn_msgbox({'This option is only available when a CONN project has been loaded (e.g. from the main CONN gui)','Please load a CONN project and try again'},'',2); return;
                end
                validconditions=[];
                for ncondition=1:numel(CONN_x.Setup.conditions.names)-1,
                    [icondition,isnewcondition]=conn_conditionnames(CONN_x.Setup.conditions.names{ncondition});
                    if ~isnewcondition&&isempty(CONN_x.Setup.conditions.model{ncondition})&&conn_existfile(fullfile(CONN_x.folders.preprocessing,['vvPC_Subject',num2str(1,'%03d'),'_Condition',num2str(icondition,'%03d'),'.mat'])), validconditions=[validconditions, ncondition]; end
                end
                if isempty(validconditions), fprintf('warning: no conditions have been denoised in the voxel-to-voxel pipeline. Please repeat Denoising to enable these analyses'); return; end
                conditionoptions=cellfun(@(x)['Display connectivity during ',x, ' condition'],CONN_x.Setup.conditions.names(validconditions),'uni',0);
                if numel(validconditions)>1, conditionoptions=[conditionoptions, {'User-defined between-conditions contrast'}]; end
                thfig=dialog('units','norm','position',[.3,.3,.4,.3],'windowstyle','normal','name','Network display','color','w','resize','on');
                uicontrol(thfig,'style','text','units','norm','position',[.1,.75,.2,.10],'string','Seed:','backgroundcolor','w','fontsize',9+CONN_gui.font_offset,'horizontalalignment','left','fontweight','bold');
                ht0=uicontrol(thfig,'style','popupmenu','units','norm','position',[.3,.75,.6,.10],'string',{'Display connectivity with selected cluster in current analysis','Display connectivity with other seed (select mask/ROI file)'},'value',1,'fontsize',8+CONN_gui.font_offset,'userdata',[],'callback',@conn_vproject_callbackfcn0);
                uicontrol(thfig,'style','text','units','norm','position',[.1,.60,.2,.10],'string','Subjects:','backgroundcolor','w','fontsize',9+CONN_gui.font_offset,'horizontalalignment','left','fontweight','bold');
                ht1=uicontrol(thfig,'style','popupmenu','units','norm','position',[.3,.60,.6,.10],'string',{'Display same between-subjects contrast as in current analysis','Display average connectivity across all subjects','Display other contrast across subjects'},'value',2,'fontsize',8+CONN_gui.font_offset,'userdata',[],'callback',@conn_vproject_callbackfcn1);
                uicontrol(thfig,'style','text','units','norm','position',[.1,.45,.2,.10],'string','Conditions:','backgroundcolor','w','fontsize',9+CONN_gui.font_offset,'horizontalalignment','left','fontweight','bold');
                idxcondition=find(strcmp(CONN_x.Setup.conditions.names(validconditions),'rest'),1);
                if isempty(idxcondition), idxcondition=1; end
                ht2=uicontrol(thfig,'style','popupmenu','units','norm','position',[.3,.45,.6,.10],'string',conditionoptions,'value',idxcondition,'fontsize',8+CONN_gui.font_offset,'userdata',[],'callback',@(varargin)conn_vproject_callbackfcn2(validconditions));
                uicontrol(thfig,'style','text','units','norm','position',[.1,.30,.2,.10],'string','Display:','backgroundcolor','w','fontsize',9+CONN_gui.font_offset,'horizontalalignment','left','fontweight','bold');
                ht3=uicontrol(thfig,'style','popupmenu','units','norm','position',[.3,.30,.6,.10],'string',{'Surface display','Slice display'},'value',1,'fontsize',8+CONN_gui.font_offset);
                uicontrol(thfig,'style','pushbutton','string','Ok','units','norm','position',[.1,.01,.38,.15],'callback','uiresume','fontsize',8+CONN_gui.font_offset);
                uicontrol(thfig,'style','pushbutton','string','Cancel','units','norm','position',[.51,.01,.38,.15],'callback','delete(gcbf)','fontsize',8+CONN_gui.font_offset);
                uiwait(thfig);
                ok=ishandle(thfig);
                if ~ok, return; end
                val1=get(ht1,'value');
                switch(val1)
                    case 1,
                        SPM=struct; conn_loadmatfile(spmfile,'SPM');
                        X=SPM.xX_multivariate.X; % design
                        C=SPM.xX_multivariate.C; % between-subjects contrast
                        validsubjects=find(SPM.xX.SelectedSubjects);
                    case 2,
                        X=ones(CONN_x.Setup.nsubjects,1); % design
                        C=1;
                        validsubjects=1:CONN_x.Setup.nsubjects;
                    case 3,
                        tdata=get(ht1,'userdata');
                        X=cell2mat(arrayfun(@(n)cell2mat(CONN_x.Setup.l2covariates.values{n}(tdata.selected)),(1:CONN_x.Setup.nsubjects)','uni',0));
                        C=tdata.contrast;
                        validsubjects=1:CONN_x.Setup.nsubjects;
                        remove=all(X==0,2)|any(isnan(X),2);
                        X(remove,:)=[];
                        validsubjects(remove)=[];
                end
                val2=get(ht2,'value');
                if val2<=numel(validconditions)
                    ncond=validconditions(val2);
                    ccond=1;
                else
                    tdata=get(ht2,'userdata');
                    ncond=validconditions(tdata.selected);
                    ccond=tdata.contrast;
                end
                if nnz(ccond<0)&nnz(X*C'<0), dtxt='average connectivity with cluster contrasted across subjects and conditions (diff r)';
                elseif nnz(ccond<0), dtxt='average connectivity with cluster contrasted across conditions (diff r)';
                elseif nnz(X*C'<0), dtxt='average connectivity with cluster contrasted across subjects (diff r)';
                else dtxt='average connectivity with cluster (r)';
                end
                if get(ht0,'value')>1, option_roi=get(ht0,'userdata'); else option_roi=''; end
                option_display=get(ht3,'value');
                delete(thfig);
                assert(size(X,2)==size(C,2),'expected %d values in between-subjects contrast, found %d',size(X,2),size(C,2));
                assert(size(ccond,2)==numel(ncond),'expected %d values in between-conditions contrast, found %d',size(ccond,2),numel(ncond));
                alpha=X*pinv(X'*X)*C'; % alpha for N=inf

                hmsg=conn_msgbox('Computing seed-to-voxel maps for selected seed region. Please wait...','',-1);
                if isempty(option_roi),
                    filename=fullfile(fileparts(spmfile),'results.nii');
                    conn_vproject(GCF,[],'export_mask_selected',filename);
                    filenamein=fullfile(fileparts(spmfile),'results.Mask.nii');
                    filenameout=fullfile(fileparts(spmfile),'results.SeedtoVoxelMap.nii');
                else
                    filenamein=option_roi;
                    filenameout=conn_prepend('',option_roi,'.SeedtoVoxelMap.nii');
                end
                conn_process('vv2rr',filenamein,'style','vv2rv','validconditions',ncond,'contrastconditions',ccond,'validsubjects',validsubjects,'contrastsubjects',alpha','saveas',filenameout);
                if option_display==1
                    fh=conn_mesh_display(filenameout);
                    fh('visible','off');
                    fh('colorbar','rescale','symmetric');
                    fh('colormap','bluewhitered');
                    fh('colormap','darker');
                    fh('brain',2);
                    fh('mask','off');
                    fh('colorbar','on', dtxt);
                    fh('visible','on');
                    %fh('colorbar','rescale','symmetric');
                    %fh('material',[]);
                    if ishandle(hmsg), delete(hmsg); end
                    if ~isempty(regexp(OPTION,'print$')), fh('background',.95*[1 1 1]); fh('print',7,options{:}); fh('close'); end
                else
                    fh=conn_slice_display(filenameout);
                    fh('colorbar','rescale','symmetric');
                    fh('colormap','bluewhitered');
                    fh('act_transparency',.8);
                    fh('colorbar','on', dtxt);
                    %fh('colorbar','rescale','symmetric');
                    fh('multisliceset',1,16,8);
                    fh('pointer_mm',[0 0 10]);
                    fh('togglegui',1);
                    if ishandle(hmsg), delete(hmsg); end
                    if ~isempty(regexp(OPTION,'print$')),
                        fh('background',[1 1 1]);
                        fh('print',options{:});
                        fh('close');
                    end
                end
                return
                
            case {'cluster_import','cluster_view'}
                if ~isempty(OPTION2), 
                    tfilename=OPTION2;
                    tspmfile=spmfile;
                    tviewrex=false; 
                else 
                    filename=fullfile(fileparts(spmfile),'results.nii');
                    conn_vproject(GCF,[],'export_mask',filename);
                    filename=fullfile(fileparts(spmfile),'results.ROIs.nii');
                    if strcmpi(OPTION,'cluster_view')&&~conn_server('util_isremotefile',filename), tviewrex=false; else tviewrex=[]; end
                    [tfilename,tspmfile,tviewrex]=conn_vproject_selectfiles(filename,spmfile,tviewrex);
                end
                if isempty(tfilename), return; end
                [tspmfile_path,tspmfile_name]=fileparts(tspmfile);
                cwd=pwd;
                conn_fileutils('cd',tspmfile_path);
                
                if tviewrex
                    conn_rex(tspmfile,tfilename,'output_type','saverex','level','clusters','select_clusters',0,'s',[],'gui',1);%,'mstats',false);
                else
                    htfig=conn_msgbox('Loading connectivity values. Please wait','',-1);
                    [y,name,info]=conn_rex(tspmfile,tfilename,'level','clusters');
                    if ishandle(htfig), delete(htfig); end
                    name=regexprep(name,'^results\.ROIs\.?','cluster ');
                    if strcmpi(OPTION,'cluster_import')
                        %temp=conn_loadmatfile(tspmfile,'SPM');
                        if isfield(info.SPM.SPM.xX,'SelectedSubjects')&&~rem(size(y,1),nnz(info.SPM.SPM.xX.SelectedSubjects)) % fill-in with NaN for missing data
                            ty=nan(size(y,1)/nnz(info.SPM.SPM.xX.SelectedSubjects)*numel(info.SPM.SPM.xX.SelectedSubjects),size(y,2));
                            ty(repmat(logical(info.SPM.SPM.xX.SelectedSubjects),size(y,1)/nnz(info.SPM.SPM.xX.SelectedSubjects),1),:)=y;
                            y=ty;
                        end
                        if (~isfield(info,'extractcontrasts')||~info.extractcontrasts)&&isfield(info.SPM.SPM,'xX_multivariate')&&isfield(info.SPM.SPM.xX_multivariate,'Zfiles')
                            if ~rem(size(y,1),CONN_x.Setup.nsubjects)&&size(y,1)/CONN_x.Setup.nsubjects==numel(info.SPM.SPM.xX_multivariate.Znames)
                                name=reshape(cellfun(@(a,b)sprintf('%s %s',a,b),repmat(name(:)',numel(info.SPM.SPM.xX_multivariate.Znames),1),repmat(info.SPM.SPM.xX_multivariate.Znames(:),1,numel(name)),'uni',0),1,[]);
                                y=reshape(y,CONN_x.Setup.nsubjects,[],size(y,2));
                            end
                        else
                            if ~rem(size(y,1),CONN_x.Setup.nsubjects)&&size(y,1)/CONN_x.Setup.nsubjects==numel(info.SPM.SPM.xX_multivariate.Ynames)
                                name=reshape(cellfun(@(a,b)sprintf('%s %s',a,b),repmat(name(:)',numel(info.SPM.SPM.xX_multivariate.Ynames),1),repmat(info.SPM.SPM.xX_multivariate.Ynames(:),1,numel(name)),'uni',0),1,[]);
                                y=reshape(y,CONN_x.Setup.nsubjects,[],size(y,2));
                            end
                        end
                        if ~isfield(CONN_x,'filename')||isempty(CONN_x.filename)||~isfield(CONN_x,'Setup')||~isfield(CONN_x.Setup,'nsubjects')||isempty(CONN_x.Setup.nsubjects)||any(rem(size(y,1),CONN_x.Setup.nsubjects)), % not conn project loaded
                            [tfilename,tfilepath]=uiputfile('*.mat','Save data as',fileparts(tspmfile));
                            if ischar(tfilename), conn_savematfile(fullfile(tfilepath,tfilename),'name','y'); fprintf('data saved to file %s\n',fullfile(tfilepath,tfilename)); end
                        else
                            conn_importl2covariate(name,num2cell(y,1));
                        end
                    else
                        s=1:length(info.ROInames);
                        cname={};mcon=[]; if isfield(info,'mstats'), mstats=info.mstats; else mstats=true; end
                        try, cname=info.SPM.SPM.xX_multivariate.Znames; mcon=info.SPM.SPM.xX_multivariate.Zcontr; end
                        if isempty(cname), try, cname=info.SPM.SPM.xX_multivariate.Ynames; mcon=info.SPM.SPM.xX_multivariate.M; end; end
                        if isfield(info.SPM.SPM,'xX_multivariate')&&isfield(info.SPM.SPM.xX_multivariate,'dof')
                            c=info.SPM.SPM.xX_multivariate.C;
                            xX=info.SPM.SPM.xX_multivariate;
                            %if numel(mcon)<=1, mstats=false; end
                        else
                            if isfield(info,'ic')&&~isempty(info.ic), Ic=info.ic;
                            elseif isfield(info,'gui')&&~info.gui&&isfield(info.SPM.SPM,'xCon')&&numel(info.SPM.SPM.xCon)==1, Ic=1;
                            else [Ic,info.SPM.SPM.xCon] = spm_conman(info.SPM.SPM,'T|F',inf,'Select contrast','',1);
                            end
                            c=[info.SPM.SPM.xCon(Ic).c]';
                            if isempty(cname), cname={info.SPM.SPM.xCon(Ic).name}; end
                            xX=info.SPM.SPM.xX;
                            mcon=1;
                            %mstats=false;
                        end
                        if ~isfield(info,'ROIinfo'), info.ROIinfo=[]; end
                        if ~isfield(info,'gui'), info.gui=true; end
                        conn_rex('test',xX,info.ROIdata(:,s),c,cname,{info.ROInames{s}},s,info.ROIinfo,mstats,mcon,info.SPM.SPM,true);
                    end
                end
                cd(cwd);
                return;
                
            case {'reference_atlas','ref_atlas'}
                if nargin>3, DATA.refatlas=OPTION2;
                else
                    [tfilename,tpathname]=conn_fileutils('uigetfile','*.nii; *.img','Select background reference atlas',DATA.refatlas);
                    if ischar(tfilename), DATA.refatlas=fullfile(tpathname,tfilename); end
                    %DATA.refatlas=[];
                end
                doredraw=true;
                init=-1;
                DATA.selectcluster=[];
            case 'subjects_view'
                dataplot=conn_displaysubject(spmfile);
%                 [spmfile_path,spmfile_name]=fileparts(spmfile);
%                 cwd=pwd;cd(spmfile_path);
%                 conn_loadmatfile(spmfile,'SPM');
%                 idxsubjects=1:size(SPM.xX.X,1);
%                 if isfield(SPM.xX,'SelectedSubjects'), consubjects=find(SPM.xX.SelectedSubjects); consubjects=consubjects(1+rem(idxsubjects-1,numel(consubjects))); 
%                 else consubjects=idxsubjects;
%                 end
%                 datafiles=reshape({SPM.xY.VY.fname},size(SPM.xY.VY));
%                 
%                 ok=true;
%                 if isfield(SPM.xX,'isSurface')&&SPM.xX.isSurface, dispopts={'Surface display'};
%                 elseif isfield(SPM.xX,'SelectedSubjects'), dispopts={'Slice display (reference anatomical)','Slice display (own subject anatomical)','Surface display','Volume display','Glass display'};
%                 else dispopts={'Slice display','Surface display','Volume display','Glass display'};
%                 end
%                 thfig=dialog('units','norm','position',[.3,.4,.3,.3],'windowstyle','normal','name','Plot individual subject','color','w','resize','on');
%                 uicontrol(thfig,'style','text','units','norm','position',[.1,.85,.8,.10],'string','Display type:','horizontalalignment','left','backgroundcolor','w','fontsize',9+CONN_gui.font_offset,'fontweight','bold');
%                 ht1=uicontrol(thfig,'style','popup','units','norm','position',[.1,.75,.8,.10],'max',2,'string',dispopts,'fontsize',8+CONN_gui.font_offset,'horizontalalignment','left','tooltipstring','select type of display');
%                 uicontrol(thfig,'style','text','units','norm','position',[.1,.6,.8,.10],'string','Subject(s):','horizontalalignment','left','backgroundcolor','w','fontsize',9+CONN_gui.font_offset,'fontweight','bold');
%                 ht2=uicontrol(thfig,'style','listbox','units','norm','position',[.1,.4,.8,.20],'max',2,'string',arrayfun(@(n)sprintf('subject %d',n),consubjects,'uni',0),'fontsize',8+CONN_gui.font_offset,'horizontalalignment','left','tooltipstring','select individual subject(s) to display');
%                 uicontrol(thfig,'style','text','units','norm','position',[.1,.25,.40,.10],'string','Threshold:','horizontalalignment','left','backgroundcolor','w','fontsize',9+CONN_gui.font_offset,'fontweight','bold');
%                 ht3=uicontrol(thfig,'style','edit','units','norm','position',[.5,.25,.4,.10],'max',1,'string','.25','fontsize',8+CONN_gui.font_offset,'horizontalalignment','left','tooltipstring','select voxel-level threshold (voxels with absolute connectivity values below this threshold are not displayed / color coded)');
%                 uicontrol(thfig,'style','text','units','norm','position',[.1,.15,.40,.10],'string','Colormap range:','horizontalalignment','left','backgroundcolor','w','fontsize',9+CONN_gui.font_offset,'fontweight','bold');
%                 ht4=uicontrol(thfig,'style','edit','units','norm','position',[.5,.15,.4,.10],'max',1,'string','-1 1','fontsize',8+CONN_gui.font_offset,'horizontalalignment','left','tooltipstring','select colorbar range (range of connectivity values color coded)');
%                 uicontrol(thfig,'style','pushbutton','string','Ok','units','norm','position',[.1,.01,.38,.10],'callback','uiresume','fontsize',8+CONN_gui.font_offset);
%                 uicontrol(thfig,'style','pushbutton','string','Cancel','units','norm','position',[.51,.01,.38,.10],'callback','delete(gcbf)','fontsize',8+CONN_gui.font_offset);
%                 while ok
%                     uiwait(thfig);
%                     ok=ishandle(thfig);
%                     if ok,
%                         dispopt=dispopts{get(ht1,'value')};
%                         thr=str2num(get(ht3,'string'));
%                         vrange=str2num(get(ht4,'string'));
%                         nidx=get(ht2,'value');
%                         if ~isempty(thr)&&~isempty(vrange)
%                             ok=false;
%                             delete(thfig);
%                         end
%                     else dispopt=[];
%                     end
%                 end                
%                 if isempty(dispopt),dataplot=[];cd(cwd);return; end
%                 idxsubjects=idxsubjects(nidx);
%                 consubjects=consubjects(nidx);
%                 switch(dispopt)
%                     case 'Volume display'
%                         fhall={};
%                         for n=1:numel(idxsubjects)
%                             fh=conn_mesh_display([],datafiles(idxsubjects(n),:),[],[],[],thr);
%                             set(fh('figurehandle'),'name',['Subject ',mat2str(consubjects(n))]);
%                             fhall{end+1}=fh;
%                         end
%                         fh=fhall;
%                     case 'Surface display'
%                         fhall={};
%                         for n=1:numel(idxsubjects)
%                             fh=conn_mesh_display(datafiles(idxsubjects(n),:),[],[],[],[],thr);
%                             %fh('colormap','hot');
%                             fh('colorbar','rescale',vrange);
%                             fh('colorbar','on');
%                             set(fh('figurehandle'),'name',['Subject ',mat2str(consubjects(n))]);
%                             fhall{end+1}=fh;
%                         end
%                         fh=fhall;
%                     case {'Slice display','Slice display (reference anatomical)'},
%                         fh=conn_slice_display(datafiles(idxsubjects,:),'',...
%                             spmfile_path,thr);
%                         fh('contour_transparency',1);
%                         fh('colormap','hot');
%                         fh('colorbar','rescale',vrange);
%                         fh('colorbar','on');
%                         set(fh('figurehandle'),'name',['Subject(s) ',mat2str(consubjects(:)')]);
%                     case 'Slice display (own subject anatomical)'
%                         fhall={};
%                         for n=1:numel(idxsubjects)
%                             fh=conn_slice_display(datafiles(idxsubjects(n),:),CONN_x.Setup.structural{consubjects(n)}{1}{1},...
%                                 spmfile_path,thr);
%                             fh('contour_transparency',1);
%                             fh('colormap','hot');
%                             fh('colorbar','rescale',vrange);
%                             fh('colorbar','on');
%                             set(fh('figurehandle'),'name',['Subject ',mat2str(consubjects(n))]);
%                             fhall{end+1}=fh;
%                         end
%                         fh=fhall;
%                     case 'Glass display'
%                         fhall={};
%                         for n=1:numel(idxsubjects)
%                             fh=conn_mesh_display([],datafiles(idxsubjects(n),:),[],[],[],thr);
%                             set(fh('figurehandle'),'name',['Subject ',mat2str(consubjects(n))]);
%                             fh('brain',4);
%                             fh('background',[1 1 1]);
%                             fh('brain_transparency',0);
%                             fh('sub_transparency',0);
%                             fh('mask_transparency',.2);
%                             fh('material',[.1 1 1 .25 0]);
%                             fhall{end+1}=fh;
%                         end
%                         fh=fhall;
%                 end
%                 dataplot=fh;
%                 cd(cwd);
                return;
                
            case 'openfolder'
                conn_fileutils('cd',fileparts(spmfile));
                try
                    if ispc, [nill,nill]=system(sprintf('start "%s"',fileparts(spmfile)));
                    else [nill,nill]=system(sprintf('open ''%s''',fileparts(spmfile)));
                    end
                end
                return
                
            case 'bookmark',
                if ~isempty(OPTION2), 
                    if iscell(OPTION2), tfilename=OPTION2{1}; descr=OPTION2{2};
                    else tfilename=OPTION2; descr='';
                    end
                else
                    tfilename=[]; %fullfile(fileparts(spmfile),'results.bookmark.jpg');
                    descr='';
                end
                if isfield(CONN_gui,'slice_display_skipbookmarkicons'), SKIPBI=CONN_gui.slice_display_skipbookmarkicons;
                else SKIPBI=false;
                end
                conn_args={'display',spmfile,[],DATA.thres,DATA.side,DATA.parametric};
                opts={'forcecd'};
                [fullfilename,tfilename,descr]=conn_bookmark('save',...
                    tfilename,...
                    descr,...
                    conn_args,...
                    opts);
                if isempty(fullfilename), return; end
                if ~SKIPBI, 
                    tht=conn_msgbox('Printing bookmark icon. Please wait...','',-1);
                    conn_print(GCF,conn_prepend('',fullfilename,'.jpg'),'-nogui','-r50','-nopersistent'); 
                    if ishandle(tht), delete(tht); end
                end
                return;
            case 'connectome',
                optionDistancePeaks=12; % minimum distance between extracted peaks (mm)
                switch(DATA.side)
                    case 1, [peaks,peaks_idx]=conn_peaks({d,mat{1}},optionDistancePeaks);
                    case 2, [peaks,peaks_idx]=conn_peaks({-d,mat{1}},optionDistancePeaks);
                    case 3, [peaks,peaks_idx]=conn_peaks({abs(d),mat{1}},optionDistancePeaks);
                end
                if 0 % one peak per cluster
                    peaks_C=zeros(size(d));
                    for n1=1:numel(DATA.clusters), peaks_C(DATA.clusters{n1})=n1; end
                    withinClusterPeak=zeros(numel(DATA.clusters)+1,1);
                    withinClusterPeak(1+peaks_C(flipud(peaks_idx)))=flipud(peaks_idx);
                    otherwithinClusterPeak=find(peaks_C(peaks_idx)>0&~ismember(peaks_idx,withinClusterPeak(2:end)));
                    peaks_idx(otherwithinClusterPeak)=[];
                    peaks(otherwithinClusterPeak,:)=[];
                end
                peaks_suprathreshold=b(peaks_idx+size(b,1)*size(b,2)*size(b,3))>0;
                [spmfile_path,spmfile_name]=fileparts(spmfile);
                conn_savematfile(fullfile(spmfile_path,'PEAKS.mat'),'peaks','peaks_suprathreshold');
                conn_process('extract_connectome',peaks(peaks_suprathreshold,:),[peaks(peaks_suprathreshold,:);peaks(~peaks_suprathreshold,:)],-1);
                return
%                 ROI=conn_process('results_connectome',spmfile_path,-1);
%                 conn_savematfile(fullfile(spmfile_path,'ROI.mat'),'ROI');
%                 conn_displayroi('initfile','results_roi',0,fullfile(spmfile_path,'ROI.mat'));
                %conn_displayroi('initfile','results_connectome',0,spmfile_path,peaks_suprathreshold,fullfile(spmfile_path,'ROI.mat'));
                %conn_displayroi('init','results_connectome',0,spmfile_path,peaks_suprathreshold);
                
            case 'fwec.voxellevel.side',
                if isempty(OPTION2), side=get(DATA.handles(11),'value');
                else side=OPTION2;
                end
                if DATA.mat{6}=='T'
                    set(DATA.handles(11),'value',side);
                    if side~=DATA.side,
                        switch(side),
                            case 1, b=c;
                            case 2, b=1-c;
                            case 3, b=2*min(c,1-c);
                        end
                        %b(b==0)=nan;
                        b=-log(max(eps,b));
                        b(isnan(c))=nan;
                        DATA.side=side;
                        doredraw=true;
                        init=-1;
                        DATA.selectcluster=[];
                    end
                end
            case {'computenonparametric','computeparametric'}
                if isempty(OPTION2), skipgui=false;
                else skipgui=OPTION2;
                end
                if DATA.parametric==1, conn_msgbox('This option is only available when using non-parametric statistics','',2); return; end
                switch(DATA.thres{2}),
                    case 1, THR_TYPE=1; %'vox-unc',
                    case 2, THR_TYPE=3; %'fdr-all'
                    case 3, THR_TYPE=5;%'TFCE',
                    case 4, THR_TYPE=5;%'TFCE',
                    case 5, THR_TYPE=4;%'T/F/X stat',
                end
                THR=DATA.thres{1};
                SIDE=DATA.side;
                if isstruct(skipgui), guiparams=skipgui;
                else guiparams=~skipgui; 
                end
                if conn_vproject_randomise(spmfile,THR_TYPE,THR,guiparams,CONN_gui.isremote);
                    try, 
                        SIM=conn_loadmatfile(conn_vproject_simfilename(spmfile,THR_TYPE,THR)); 
                        if ~isfield(SIM,'VERSION'), SIM.VERSION=0; end
                        if SIM.VERSION<2, SIM=[]; conn_fileutils('spm_unlink',conn_vproject_simfilename(spmfile,THR_TYPE,THR)); end % note: disregard older versions
                    end
                end
                init=-1;
                DATA.SIM=SIM;
                DATA.selectcluster=[];
            case 'fwec.clusterlevel.parametric',
                if nnz(paramoptions)>1
                    if isempty(OPTION2), parametric=get(DATA.handles(15),'value');
                    else parametric=OPTION2;
                    end
                    thres4=DATA.thres{4}+(DATA.thres{4}>=6);
                    if parametric==1&&thres4>8, conn_msgbox('This threshold option is only available when using non-parametric statistics','',2); parametric=2;
                    elseif parametric==2&&ismember( thres4,[5 6 7]), conn_msgbox('This threshold option is only available when using parametric statistics','',2); parametric=1;
                    end
                    set(DATA.handles(15),'value',parametric);
                    if parametric~=DATA.parametric,
                        switch(parametric),
                            case 1,
                                a=DATA.param.backg;
                                b=DATA.param.logp;
                                c=DATA.param.p;
                                d=DATA.param.F;
                                mat=DATA.param.stats;
                                %set(DATA.handles(7),'string',sprintf('%15s%13s%13s%13s%13s%13s%13s','Clusters (x,y,z)','size','size p-FWE','size p-FDR','size p-unc','peak p-FWE','peak p-unc'));
                            case 2,
                                a=DATA.nonparam.backg;
                                b=DATA.nonparam.logp;
                                c=DATA.nonparam.p;
                                d=DATA.nonparam.F;
                                mat=DATA.nonparam.stats;
                                %set(DATA.handles(7),'string',sprintf('%15s%13s%13s%13s%13s%13s%13s%13s%13s','Clusters (x,y,z)','size','size p-FWE','size p-FDR','size p-unc','mass','mass p-FWE','mass p-FDR','mass p-unc'));
                        end
                        %if (DATA.thres{2}==3||DATA.thres{2}==4), set(DATA.handles(7),'string',sprintf('%15s%13s%13s%13s','Clusters (x,y,z)','size','peak TFCE','peak p-FWE')); end
                        if DATA.mat{6}=='T'
                            switch(DATA.side),
                                case 1, b=c;
                                case 2, b=1-c;
                                case 3, b=2*min(c,1-c);
                            end
                            b=-log(max(eps,b));
                            b(isnan(c))=nan;
                        end
                        DATA.a=a;DATA.b=b;DATA.c=c;DATA.d=d;DATA.mat=mat;DATA.parametric=parametric;
                        %                     %b(b==0)=nan;
                        %                     b=-log(max(eps,b));
                        init=-1;
                        doredraw=true;
                        DATA.selectcluster=[];
                    end
                end
            case 'advancedthr'
                if isempty(OPTION2), advanced=get(DATA.handles(40),'value');
                else advanced=OPTION2; set(DATA.handles(40),'value',advanced);
                end
                if advanced
                    conn_vproject(GCF,[],'fwec.option',4,'immediatereturn');
%                     icp=[1 2 3 4 5 6  11 15];
%                     set(DATA.handles(icp),'visible','on');
%                     if thres{2}==3||thres{2}==4, set(DATA.handles([4 5 6  15]),'visible','off'); end
%                     set(DATA.handles(33),'visible','off','string','');
                    return
                else
                    ok=false;for method=1:numel(DATA.thres_defaults), if isequal(DATA.thres,DATA.thres_defaults{method}(1:4))&&DATA.side==3, ok=true; break; end; end; 
                    if ~ok, 
                        conn_vproject(GCF,[],'fwec.option',4,'immediatereturn');
%                         set(DATA.handles(38),'value',true);
                    else
                        conn_vproject(GCF,[],'fwec.option',method);%,'immediatereturn');
%                         icp=[1 2 3 4 5 6  11 15]; 
%                         set(DATA.handles(icp),'visible','off');
%                         set(DATA.handles(33),'visible','on');
                    end
                    return
                end
                
            case 'fwec.option'
                if isempty(OPTION2), method=get(DATA.handles(32),'value');
                else method=OPTION2; set(DATA.handles(32),'value',method);
                end
                advanced=get(DATA.handles(40),'value');
                if method==1&&~DATA.paramoptions(1)&&DATA.issurface, conn_msgbox('Sorry, Random Field Theory parametric statistics not available in surface-based analyses','Option not available',2); set(DATA.handles(32),'value',4); return; 
                elseif method==1&&~DATA.paramoptions(1), conn_msgbox('Parametric analysis option not available. Please re-run second-level analyses','Incomplete analysis files',2); set(DATA.handles(32),'value',4); return; 
                elseif (method==2||method==3)&&~DATA.paramoptions(2), conn_msgbox('Non-parametric analysis option not available. Please re-run second-level analyses','Incomplete analysis files',2); set(DATA.handles(32),'value',4); return; 
                end
                icp=[1 2 3 4 5 6  11 15]; 
                switch(method)
                    case {1,2,3},
                        if method==1, thres=DATA.thres_defaults{1}(1:4); parametric=DATA.thres_defaults{1}{5};
                        elseif method==2, thres=DATA.thres_defaults{2}(1:4); parametric=DATA.thres_defaults{2}{5};
                        else thres=DATA.thres_defaults{3}(1:4); parametric=DATA.thres_defaults{3}{5};
                        end
                        side=3;
                        if advanced, 
                            set(DATA.handles(icp),'visible','on');
                            if thres{2}==3||thres{2}==4, set(DATA.handles([4 5 6 15]),'visible','off'); end
                        else set(DATA.handles(icp),'visible','off');
                        end
                        set(DATA.handles(1),'string','voxel threshold: p <'); 
                        set(DATA.handles(2),'string',num2str(thres{1}));
                        set(DATA.handles(3),'value',thres{2});
                        if method==3, set(DATA.handles(4),'string','cluster threshold: k >'); 
                        else set(DATA.handles(4),'string','cluster threshold: p <'); 
                        end
                        set(DATA.handles(5),'string',num2str(thres{3}));
                        set(DATA.handles(6),'value',thres{4});
                        set(DATA.handles(11),'value',side);
                        set(DATA.handles(15),'value',parametric);
                        if advanced, set(DATA.handles(33),'visible','off','string','');
                        else set(DATA.handles(33),'visible','on');
                        end
                        if parametric==1, 
                            a=DATA.param.backg;
                            b=DATA.param.logp;
                            c=DATA.param.p;
                            d=DATA.param.F;
                            mat=DATA.param.stats;
                            %set(DATA.handles(7),'string',sprintf('%15s%13s%13s%13s%13s%13s%13s','Clusters (x,y,z)','size','size p-FWE','size p-FDR','size p-unc','peak p-FWE','peak p-unc'));
                        else 
                            a=DATA.nonparam.backg;
                            b=DATA.nonparam.logp;
                            c=DATA.nonparam.p;
                            d=DATA.nonparam.F;
                            mat=DATA.nonparam.stats;
                            %set(DATA.handles(7),'string',sprintf('%15s%13s%13s%13s%13s%13s%13s%13s%13s','Clusters (x,y,z)','size','size p-FWE','size p-FDR','size p-unc','mass','mass p-FWE','mass p-FDR','mass p-unc'));
                        end
                        %if thres{2}==3, set(DATA.handles(7),'string',sprintf('%15s%13s%13s%13s','Clusters (x,y,z)','size','peak TFCE','peak p-FWE')); end
                        if mat{6}=='T'
                            switch(side),
                                case 1, b=c;
                                case 2, b=1-c;
                                case 3, b=2*min(c,1-c);
                            end
                            b=-log(max(eps,b));
                            b(isnan(c))=nan;
                        end
                        DATA.a=a;DATA.b=b;DATA.c=c;DATA.d=d;DATA.mat=mat;DATA.parametric=parametric;
                        DATA.side=side;
                        DATA.thres=thres;
                        set(GCF,'userdata',DATA);
                        if isequal(OPTION3, 'immediatereturn'), return; end
                        init=-1;
                        doredraw=true;
                        DATA.selectcluster=[];                        
                    case 4
                        set(DATA.handles(icp),'visible','on');
                        if thres{2}==3||thres{2}==4, set(DATA.handles([4 5 6  15]),'visible','off'); end
                        set(DATA.handles(33),'visible','off','string','');
                        return
                end
            case 'close'
                delete(GCF);
                return
                
            otherwise,
                conn_disp('fprintf','unrecognized conn_vproject option %s\n',lower(OPTION));
                return
        end
    end
else, %initialization
    init=1;
    doredraw=false;
    if nargin<15 || isempty(issurface), issurface=false; end
    if nargin<14 || isempty(voxeltovoxel), voxeltovoxel=0; end
    if nargin<13 || isempty(spmfile), spmfile=[]; end
    if nargin<12 || isempty(data1plot), data1plot=[]; end
    if nargin<11 || isempty(threshold), threshold=.5; end
    if nargin<10 || isempty(select), select=[]; end
    if nargin<9 || isempty(box), box=0; end
    if nargin<8 || isempty(res), res=1; end
    if nargin<7 || isempty(parametric), parametric=[]; end
    if nargin<6 || isempty(side), side=1; end
    if nargin<5 || isempty(thres), thres={.05,1,0,1}; end
    if nargin<4 || isempty(projection), projection=[-1,0,0;0,0,-1;0,-1,0]; end; 
    %if nargin<4 || isempty(projection), projection=[-0.3529,0.9347,-0.0429;0.0247,-0.0365,-0.9990;-0.9353,-0.3537,-0.0103]; end
    if nargin<3 || isempty(views), views='full'; end
%     if nargin<9 || isempty(threshold), if isempty(b), threshold=[nanmean(a(:))]; else, threshold=[nanmean(a(:)),2]; end; end
    pointer=[0 0 0];
    vproject_display;

    GCF=gcf;
    if ~strcmp(views,'none'), 
        clf; 
        paramoptions=[~isempty(param.F), ~isempty(nonparam.F)];
        if isempty(parametric)
            if paramoptions(1), parametric=1;
            elseif paramoptions(2), parametric=2;
            else error('no analyses found');
            end
        end
        if parametric==1
            a=param.backg;
            b=param.logp;
            c=param.p;
            d=param.F;
            mat=param.stats;
        else
            a=nonparam.backg;
            b=nonparam.logp;
            c=nonparam.p;
            d=nonparam.F;
            mat=nonparam.stats;
        end
        if mat{6}=='T'&&side~=1
            switch(side),
                case 1, b=c;
                case 2, b=1-c;
                case 3, b=2*min(c,1-c);
            end
            b=-log(max(eps,b));
            b(isnan(c))=nan;
        end
        e=[];f=[];g=[];
        DATA.a=a;
        DATA.b=b;
        DATA.c=c;
        DATA.d=d;
        DATA.e=e;
        DATA.f=f;
        DATA.g=g;
        DATA.projection=projection;
        DATA.thres_defaults={ {.001,1,.05,3, 1}, {.010,1,.05,9, 2}, {.05,4,0,1, 2} };
        if isnumeric(thres), parametric=DATA.thres_defaults{thres}{5}; thres=DATA.thres_defaults{thres}(1:4); side=3; end
        DATA.thres=thres; 
        DATA.threshold=threshold;
        DATA.res=res;
        DATA.box=box;
        DATA.views=views;
        DATA.rotation=[0,0,0];
        DATA.selectcluster=[];
        DATA.select=select;
        DATA.pointer=pointer;
        DATA.mat=mat;
        DATA.clusters={[]};
        DATA.side=side;
        DATA.spmfile=spmfile;
        DATA.voxeltovoxel=voxeltovoxel;
        DATA.issurface=issurface;
        if issurface, 
            surfparams=load(fullfile(fileparts(which(mfilename)),'utils','surf','surf_top.mat'),'A');
            DATA.datatype=kron(speye(2),surfparams.A); % explicit adjacency matrix
            surfparams={conn_surf_readsurf(fullfile(fileparts(which(mfilename)),'utils','surf','lh.white.surf')),conn_surf_readsurf(fullfile(fileparts(which(mfilename)),'utils','surf','lh.pial.surf')),conn_surf_readsurf(fullfile(fileparts(which(mfilename)),'utils','surf','rh.white.surf')),conn_surf_readsurf(fullfile(fileparts(which(mfilename)),'utils','surf','rh.pial.surf'))};
            DATA.surfcoords=[(surfparams{1}.vertices+surfparams{2}.vertices)/2;(surfparams{3}.vertices+surfparams{4}.vertices)/2];
        else
            DATA.datatype='volume';
            DATA.surfcoords=[];
        end
        DATA.refatlas=[];
        DATA.SIM=[];
        DATA.param=param;
        DATA.nonparam=nonparam;
        DATA.parametric=parametric;
        DATA.paramoptions=paramoptions;
        DATA.peakFDR=0;
    end
end

if 1, backgroundcolor=[.975 .975 .975];
elseif isfield(CONN_gui,'backgroundcolor'), backgroundcolor=CONN_gui.backgroundcolor;
else backgroundcolor=[0.12 0.126 0.132];
end
color3=backgroundcolor+(.5-backgroundcolor)*.15; %.9*[1 1 1];
foregroundcolor=.5*.8+.2*(1-round(backgroundcolor));
if init>0, % initialize
    %map=[gray(64).^2;.5*repmat([1,1,0],[64,1])+.5*((gray(64).^2)*diag([1,1,0]));.5*repmat([1,0,0],[64,1])+.5*((gray(64).^2)*diag([1,0,0]))];
    %map=[.5+.5*(gray(64).^2);.5*repmat([1,1,0],[64,1])+.5*((gray(64).^.5)*diag([1,1,0]));.5*repmat([0,0,1],[64,1])+.5*((gray(64).^.5)*diag([0,0,1]))];map(1,:)=1;
    map=[.25+.75*(gray(64).^2)*diag([1,1,1]);0*repmat([1,1,1],[64,1])+1*((gray(64).^1)*diag([1,0,0]));.25*repmat([1,1,0],[64,1])+.75*((gray(64).^.5)*diag([1,1,0]))];map(1,:)=.25;
    switch(views),
        case 'full',
            clf(GCF);
            if isempty(b),
                DATA.handles=[uicontrol('style','text','units','norm','position',[.6,.95,.1,.05],'string','resolution','backgroundcolor',figcolor,'foregroundcolor','b'),...
                    uicontrol('style','edit','units','norm','position',[.7,.95,.1,.05],'string',num2str(DATA.res),'callback',{@conn_vproject,'resolution'},'backgroundcolor',figcolor,'foregroundcolor',1-round(figcolor)),...
                    uicontrol('style','text','units','norm','position',[.8,.95,.1,.05],'string','threshold','backgroundcolor',figcolor,'foregroundcolor','b'),...
                    uicontrol('style','edit','units','norm','position',[.9,.95,.1,.05],'string',num2str(DATA.threshold),'callback',{@conn_vproject,'threshold'},'backgroundcolor',figcolor,'foregroundcolor','b')];
                DATA.axes=subplot(122); set(DATA.axes,'units','norm','position',[.55,.1,.4,.8],'visible','off');
                set(GCF,'name','results explorer','numbertitle','off','color',figcolor,'units','norm','position',[.1,.4,.8,.5],'interruptible','off','busyaction','cancel','colormap',map,'userdata',DATA); %,'windowbuttondownfcn',{@conn_vproject,'buttondown'},'windowbuttonupfcn',{@conn_vproject,'buttonup'},'windowbuttonmotionfcn',[],'keypressfcn',{@conn_vproject,'keypress'},
            else
                backgroundcolor2=backgroundcolor;
                uicontrol('style','frame','units','norm','position',[.0,.87,1,.13],'backgroundcolor',color3,'foregroundcolor',color3);
                %uicontrol('style','frame','units','norm','position',[.54,.42,.41,.42],'backgroundcolor',backgroundcolor2,'foregroundcolor',.85*[1,1,1]);
                %uicontrol('style','frame','units','norm','position',[.49,.39,.11,.41],'backgroundcolor',.9*[1,1,1],'foregroundcolor',.5*[1 1 1]);
                %dp=@(a,b)[.615+(a-1)*.065,.465+(b-1)*.055,.06,.05];
                dp1=@(a,b)[.68+(a-1)*.06,.51+(b-1)*.04,.05,.035];
                dp2=@(a,b)[.80+(a-1)*.15,.55+(b-1)*.04,.15,.04];
                set(GCF,'color',backgroundcolor);
                DATA.handles=[...
                    uicontrol('style','text','units','norm','position',[.05,.925,.15,.03],'string','voxel threshold: p <','fontname','arial','fontsize',8+CONN_gui.font_offset,'fontweight','bold','foregroundcolor',1-color3,'backgroundcolor',color3,'parent',GCF),...
                    uicontrol('style','edit','units','norm','position',[.20,.925,.10,.03],'string',num2str(thres{1}),'fontname','arial','fontsize',8+CONN_gui.font_offset,'callback',{@conn_vproject,'fwec.voxellevel'},'tooltipstring','Voxel-level threshold value: False-positive threshold value for individual voxels (based on T/F/X statistics)','foregroundcolor',1-color3,'backgroundcolor',color3,'parent',GCF),...
                    uicontrol('style','popupmenu','units','norm','position',[.325,.925,.25,.03],'string',{'p-uncorrected','p-FDR corrected','p-FDR corrected (TFCE)','p-FWE corrected (TFCE)','F/T/X stat'},'fontname','arial','fontsize',8+CONN_gui.font_offset,'callback',{@conn_vproject,'fwec.voxellevel'},'value',thres{2},'tooltipstring','False-positive control type for individual voxels','foregroundcolor',1-color3,'backgroundcolor',color3,'parent',GCF),...
                    uicontrol('style','text','units','norm','position',[.05,.885,.15,.03],'string','cluster threshold: p <','fontname','arial','fontsize',8+CONN_gui.font_offset,'fontweight','bold','foregroundcolor',1-color3,'backgroundcolor',color3,'parent',GCF),...
                    uicontrol('style','edit','units','norm','position',[.20,.885,.10,.03],'string',num2str(thres{3}),'fontname','arial','fontsize',8+CONN_gui.font_offset,'callback',{@conn_vproject,'fwec.clusterlevel'},'tooltipstring','Cluster-level threshold value: False-positive threshold value for individual clusters (based on cluster size/mass/peak)','foregroundcolor',1-color3,'backgroundcolor',color3,'parent',GCF),...
                    uicontrol('style','popupmenu','units','norm','position',[.325,.885,.25,.03],'string',{'cluster-size','cluster-size p-FWE corrected','cluster-size p-FDR corrected','cluster-size p-uncorrected','peak-voxel p-FWE corrected','peak-voxel p-FDR corrected','peak-voxel p-uncorrected','cluster-mass','cluster-mass p-FWE corrected','cluster-mass p-FDR corrected','cluster-mass p-uncorrected'},'fontname','arial','fontsize',7+CONN_gui.font_offset,'callback',{@conn_vproject,'fwec.clusterlevel'},'value',thres{4},'tooltipstring','False-positive control type for individual clusters','foregroundcolor',1-color3,'backgroundcolor',color3,'parent',GCF),...
                    uicontrol('style','text','units','norm','position',[.05,.30,.9,.03],'string',sprintf('%15s%13s%13s%13s%13s%13s%13s%15s','Cluster (x,y,z)','size','size p-FWE','size p-FDR','size p-unc','peak p-FWE','peak p-FDR','peak p-unc'),'fontname','arial','backgroundcolor',backgroundcolor,'foregroundcolor',foregroundcolor,'horizontalalignment','left','fontname','monospaced','fontsize',7+CONN_gui.font_offset,'parent',GCF),...
                    uicontrol('style','listbox','units','norm','position',[.05,.15,.9,.15],'string','','tag','highlight','max',2,'value',[],'backgroundcolor',backgroundcolor,'foregroundcolor',foregroundcolor,'horizontalalignment','left','fontname','monospaced','fontsize',7+CONN_gui.font_offset,'callback',{@conn_vproject,'selectroi'},'parent',GCF),...
                    uicontrol('style','pushbutton','units','norm','position',dp1(1,7),'string','Surface display','tag','highlight_image','fontname','arial','fontweight','bold','fontsize',8+CONN_gui.font_offset,'callback',{@conn_vproject,'surface_view'},'tooltipstring','<HTML><b>Surface display</b><br/>Displays selected clusters projected to cortical surface (show voxel-level stats at ICBM reference T1 surface)</HTML>','parent',GCF),...
                    uicontrol('style','text','units','norm','position',[.03,.80,.27,.04],'string','Voxels in selected-cluster','fontname','arial','fontsize',9+CONN_gui.font_offset,'backgroundcolor',backgroundcolor,'foregroundcolor',.75*[1 1 1],'horizontalalignment','left','max',2,'visible','off','value',[],'parent',GCF),...
                    uicontrol('style','popupmenu','units','norm','position',[.605,.925,.2,.03],'string',{'positive contrast (one-sided)','negative contrast (one-sided)','two-sided'},'value',side,'fontname','arial','fontsize',8+CONN_gui.font_offset,'callback',{@conn_vproject,'fwec.voxellevel.side'},'tooltipstring','Analysis results directionality','foregroundcolor',1-color3,'backgroundcolor',color3,'parent',GCF),...
                    uicontrol('style','edit','units','norm','position',[.05,.07,.9,.08],'string','','fontname','arial','fontsize',7+CONN_gui.font_offset,'visible','off','backgroundcolor',backgroundcolor,'foregroundcolor',foregroundcolor,'horizontalalignment','left','max',2,'value',[],'parent',GCF),...
                    uicontrol('style','text','units','norm','position',[.855,.930,.1,.03],'string','','fontname','arial','fontsize',8+CONN_gui.font_offset,'horizontalalignment','right','foregroundcolor',.5*[1 1 1],'backgroundcolor',color3,'parent',GCF),...
                    uicontrol('style','text','units','norm','position',[.855,.900,.1,.03],'string','','fontname','arial','fontsize',8+CONN_gui.font_offset,'horizontalalignment','right','foregroundcolor',.5*[1 1 1],'backgroundcolor',color3,'parent',GCF),...
                    uicontrol('style','popupmenu','units','norm','position',[.605,.885,.2,.03],'string',{'parametric stats','non-parametric stats'},'value',parametric,'fontname','arial','fontsize',8+CONN_gui.font_offset,'callback',{@conn_vproject,'fwec.clusterlevel.parametric'},'tooltipstring','Cluster-level statistics based on parametric analyses (Random Field Theory) or non-parametric analyses (permutation/randomization tests)','foregroundcolor',1-color3,'backgroundcolor',color3,'parent',GCF),...
                    uicontrol('style','pushbutton','units','norm','position',dp1(2,2),'string','Plot effects','tag','highlight_image','fontname','arial','fontsize',8+CONN_gui.font_offset,'callback',{@conn_vproject,'cluster_view'},'tooltipstring','<HTML><b>Plot effects</b><br/>Explores/displays average effect sizes within each significant cluster</HTML>','parent',GCF),...
                    uicontrol('style','pushbutton','units','norm','position',dp2(1,5),'string','Export mask','tag','highlight','foregroundcolor',foregroundcolor,'backgroundcolor',backgroundcolor2,'fontname','arial','fontsize',8+CONN_gui.font_offset,'callback',{@conn_vproject,'export_mask'},'tooltipstring','Exports mask of supra-threshold voxels','parent',GCF),...
                    uicontrol('style','pushbutton','units','norm','position',dp2(1,6),'string','Import values','tag','highlight','foregroundcolor',foregroundcolor,'backgroundcolor',backgroundcolor2,'fontname','arial','fontsize',8+CONN_gui.font_offset,'callback',{@conn_vproject,'cluster_import'},'tooltipstring','Imports average connectivity measure values within each significant cluster and for each subject into CONN toolbox as second-level covariates','parent',GCF),...
                    uicontrol('style','pushbutton','units','norm','position',dp2(1,2),'string','Bookmark','tag','highlight','foregroundcolor',foregroundcolor,'backgroundcolor',backgroundcolor2,'fontname','arial','fontsize',8+CONN_gui.font_offset,'callback',{@conn_vproject,'bookmark'},'tooltipstring','Bookmark this second-level results explorer view','parent',GCF),...
                    uicontrol('style','pushbutton','units','norm','position',dp1(1,6),'string','Volume display','tag','highlight_image','fontname','arial','fontweight','bold','fontsize',8+CONN_gui.font_offset,'callback',{@conn_vproject,'volume_view'},'tooltipstring','<HTML><b>Volume display</b><br/>Displays selected clusters on 3d brain (show cluster surfaces)</HTML>','parent',GCF),...
                    uicontrol('style','pushbutton','units','norm','position',dp1(1,0),'string','SPM display','tag','highlight_image','fontname','arial','fontweight','bold','fontsize',8+CONN_gui.font_offset,'callback',{@conn_vproject,'spm_view'},'tooltipstring','<HTML><b>SPM display</b><br/>Displays results in SPM</HTML>','parent',GCF),...
                    uicontrol('style','pushbutton','units','norm','position',dp1(1,4),'string','Glass display','tag','highlight_image','fontname','arial','fontweight','bold','fontsize',8+CONN_gui.font_offset,'callback',{@conn_vproject,'glass_view'},'tooltipstring','<HTML><b>Glass display</b><br/>Displays selected clusters on 3d glass-brain (show maximum-intensity-projection stats)</HTML>','parent',GCF),...
                    uicontrol('style','pushbutton','units','norm','position',dp1(1,5),'string','Slice display','tag','highlight_image','fontname','arial','fontweight','bold','fontsize',8+CONN_gui.font_offset,'callback',{@conn_vproject,'slice_view'},'tooltipstring','<HTML><b>Slice display</b><br/>Displays selected clusters on individual slices (show voxel-level stats at each slice)</HTML>','parent',GCF),...
                    uicontrol('style','pushbutton','units','norm','position',dp2(1,4),'string','non-parametric stats','tag','highlight','foregroundcolor',foregroundcolor,'backgroundcolor',backgroundcolor2,'fontname','arial','fontsize',8+CONN_gui.font_offset,'callback',{@conn_vproject,'computeparametric'},'tooltipstring','<HTML>Compute additional permutation/randomization simulations<br/>(only applicable when ''non-parametric stats'' is selected in the cluster-level threshold options above)</HTML>','parent',GCF),...
                    uicontrol('style','pushbutton','units','norm','position',dp1(1,1),'string','Plot subjects','tag','highlight_image','fontname','arial','fontsize',8+CONN_gui.font_offset,'callback',{@conn_vproject,'subjects_view'},'tooltipstring','<HTML><b>Plot subjects</b><br/>Explores/displays individual-subject maps</HTML>','parent',GCF),...
                    uicontrol('style','text','units','norm','position',[.68,.80,.12,.025],'string','Display&Print','backgroundcolor',backgroundcolor,'foregroundcolor',foregroundcolor,'horizontalalignment','center','fontweight','bold','fontname','arial','fontsize',9+CONN_gui.font_offset,'parent',GCF),...
                    uicontrol('style','text','units','norm','position',[.80,.80,.15,.025],'string','Tools','backgroundcolor',backgroundcolor,'foregroundcolor',foregroundcolor,'horizontalalignment','center','fontweight','bold','fontname','arial','fontsize',9+CONN_gui.font_offset,'parent',GCF),...
                    uicontrol('style','text','units','norm','position',[.075,.80,.4,.025],'string','Results','backgroundcolor',backgroundcolor,'foregroundcolor',foregroundcolor,'horizontalalignment','center','fontweight','bold','fontname','arial','fontsize',9+CONN_gui.font_offset,'parent',GCF),...
                    uicontrol('style','text','units','norm','position',dp2(1,3),'string','','backgroundcolor',backgroundcolor2,'foregroundcolor',.5*[1 1 1],'horizontalalignment','center','fontname','arial','fontsize',7+CONN_gui.font_offset), ... 
                    uicontrol('style','pushbutton','units','norm','position',dp1(2,1),'string','Plot design','tag','highlight_image','foregroundcolor',foregroundcolor,'fontname','arial','fontsize',8+CONN_gui.font_offset,'callback',{@conn_vproject,'plot_design'},'tooltipstring','<HTML><b>Plot design</b><br/>Displays General Linear Model design matrix and additional details</HTML>','parent',GCF),...
                    uicontrol('style','pushbutton','units','norm','position',dp2(1,1),'string','Open folder','tag','highlight','foregroundcolor',foregroundcolor,'backgroundcolor',backgroundcolor2,'fontname','arial','fontsize',8+CONN_gui.font_offset,'callback',{@conn_vproject,'openfolder'},'tooltipstring','Open folder containing current second-level results files','parent',GCF), ...
                    uicontrol('style','popupmenu','units','norm','position',[.20,.965,.605,.03],'string',{'standard settings for cluster-based inferences #1: Random Field Theory parametric statistics','standard settings for cluster-based inferences #2: Permutation/randomization analysis','standard settings for cluster-based inferences #3: Threshold Free Cluster Enhancement','<HTML><i>customize (advanced Family-Wise Error control settings)</i></HTML>'},'tag','highlight','foregroundcolor',foregroundcolor,'fontname','arial','fontsize',8+CONN_gui.font_offset,'callback',{@conn_vproject,'fwec.option'},'value',1,'tooltipstring','Select false-positive control method','backgroundcolor',.9*[1,1,1],'parent',GCF),...
                    uicontrol('style','text','units','norm','position',[.10,.88,.75,.07],'string','','horizontalalignment','center','backgroundcolor',color3,'foregroundcolor',.5*[1 1 1],'horizontalalignment','center','fontname','arial','fontsize',7+CONN_gui.font_offset,'parent',GCF), ... 
                    uicontrol('style','pushbutton','units','norm','position',dp1(2,7),'string','Surface print','tag','highlight_image','fontname','arial','fontweight','bold','fontsize',8+CONN_gui.font_offset,'callback',{@conn_vproject,'surface_print'},'tooltipstring','<HTML><b>Surface print</b><br/>Prints selected clusters projected to cortical surface (show voxel-level stats at ICBM reference T1 surface)</HTML>','parent',GCF),...
                    uicontrol('style','pushbutton','units','norm','position',dp1(2,6),'string','Volume print','tag','highlight_image','fontname','arial','fontweight','bold','fontsize',8+CONN_gui.font_offset,'callback',{@conn_vproject,'volume_print'},'tooltipstring','<HTML><b>Volume print</b><br/>Prints selected clusters on 3d brain (show cluster surfaces)</HTML>','parent',GCF),...
                    uicontrol('style','pushbutton','units','norm','position',dp1(2,5),'string','Slice print','tag','highlight_image','fontname','arial','fontweight','bold','fontsize',8+CONN_gui.font_offset,'callback',{@conn_vproject,'slice_print'},'tooltipstring','<HTML><b>Slice print</b><br/>Prints selected clusters on individual slices (show voxel-level stats at each slice)</HTML>','parent',GCF),...
                    uicontrol('style','pushbutton','units','norm','position',dp1(2,4),'string','Glass print','tag','highlight_image','fontname','arial','fontweight','bold','fontsize',8+CONN_gui.font_offset,'callback',{@conn_vproject,'glass_print'},'tooltipstring','<HTML><b>Glass print</b><br/>Prints selected clusters on 3d glass-brain (show maximum-intensity-projection stats)</HTML>','parent',GCF), ...
                    uicontrol('style','pushbutton','units','norm','position',dp1(1,3),'string','Connectivity display','tag','highlight_image','fontname','arial','fontweight','bold','fontsize',8+CONN_gui.font_offset,'callback',{@conn_vproject,'network_view'},'tooltipstring','<HTML><b>Network view</b><br/>Displays the connectivity pattern between the selected cluster and the rest of the brain<br/>(average seed-based connectivity map, with all voxels within the selected cluster as seeds)</HTML>','parent',GCF),...
                    uicontrol('style','pushbutton','units','norm','position',dp1(2,3),'string','Connectivity print','tag','highlight_image','fontname','arial','fontweight','bold','fontsize',8+CONN_gui.font_offset,'callback',{@conn_vproject,'network_print'},'tooltipstring','<HTML><b>Network print</b><br/>Prints the connectivity pattern between the selected cluster and the rest of the brain<br/>(average seed-based connectivity map, with all voxels within the selected cluster as seeds)</HTML>','parent',GCF), ...
                    uicontrol('style','checkbox','units','norm','position',[.88,.870,.12,.03],'string','show details','tag','highlight','value',0,'fontsize',8+CONN_gui.font_offset,'horizontalalignment','right','foregroundcolor',foregroundcolor,'backgroundcolor',.9*[1,1,1],'callback',{@conn_vproject,'advancedthr'},'tooltipstring','Displays advanced thresholding options','parent',GCF), ...
                    uicontrol('style','pushbutton','units','norm','position',[.81,.965,.02,.03],'fontsize',7+CONN_gui.font_offset,'foregroundcolor',foregroundcolor,'backgroundcolor',backgroundcolor2,'string','?','tag','highlight','tooltipstring','<HTML>Documentation about available methods of statistical inference</HTML>','interruptible','off','callback',@(varargin)conn('gui_help','url','http://www.conn-toolbox.org/fmri-methods/cluster-level-inferences'),'parent',GCF) , ...
                    uicontrol('style','pushbutton','units','norm','position',dp1(1,2),'string','Polar display','tag','highlight_image','fontname','arial','fontweight','bold','fontsize',8+CONN_gui.font_offset,'callback',{@conn_vproject,'polar_view'},'tooltipstring','<HTML><b>Polar display</b><br/>Displays overlap between each significant cluster and a set of canonical brain networks (Yeo&Buckner 2011 7-network parcellation)</HTML>','parent',GCF) ]; 
                    %uicontrol('style','text','units','norm','position',[.4,.30,.1,.025],'string','voxel p-cor','backgroundcolor','k','foregroundcolor','y','horizontalalignment','left'),...
                    %uicontrol('style','listbox','units','norm','position',[.2,.05,.1,.25],'string','','backgroundcolor','k','foregroundcolor','w','horizontalalignment','left'),...
                    %uicontrol('style','listbox','units','norm','position',[.3,.05,.1,.25],'string','','backgroundcolor','k','foregroundcolor','w','horizontalalignment','left'),...
                    %uicontrol('style','listbox','units','norm','position',[.4,.05,.1,.25],'string','','backgroundcolor','k','foregroundcolor','w','horizontalalignment','left')];
                uiwrap(DATA.handles([8,12,17,18,19,24,31,41]));
                if DATA.mat{6}~='T', set(DATA.handles(11),'value',3,'enable','off'); end; 
                bp=[9 20 22 23 21 16 25 30 38 34 35 37 36 39 42]; 
                bp_isprint=[0 0 0 0 0 0 0 0 0 1 1 1 1 1 0]; 
                temp=imread(fullfile(fileparts(which(mfilename)),sprintf('conn_vproject_icon%02d.jpg',0))); temp=double(temp); printmask=round(temp/255); 
                for n1=1:numel(bp),
                    set(DATA.handles(bp(n1)),'units','pixel'); pt=get(DATA.handles(bp(n1)),'position'); set(DATA.handles(bp(n1)),'units','norm'); 
                    temp=imread(fullfile(fileparts(which(mfilename)),sprintf('conn_vproject_icon%02d.jpg',n1))); temp=double(temp); temp=temp/255; temp=max(0,min(1,(temp).^.5)); ft=min(size(temp,1)/ceil(pt(4)),size(temp,2)/ceil(pt(3))); 
                    if ismember(n1,[1,2,10,11]), ft=.45*ft; elseif ismember(n1,[4,13]), ft=.50*ft; elseif ismember(n1,[9,14]), ft=.65*ft; elseif ismember(n1,[3,12]), ft=.9*ft; end; 
                    maxtemp=1;%mode(round(temp(:)*100))/100;
                    if maxtemp<.5, temp=1-temp; maxtemp=1-maxtemp; end
                    temp=max(0,min(1, .75*temp+.25*temp/maxtemp.*repmat(shiftdim(backgroundcolor,-1),[size(temp,1),size(temp,2),1,size(temp,4)]) ));
                    temp=min(.95,temp(round(1:ft:size(temp,1)),round(1:ft:size(temp,2)),:));
                    if bp_isprint(n1)
                        if ismember(n1,[3,12,14]), tempprintmask=printmask(ceil(size(printmask,1)/4)+(1:ceil(size(printmask,1)/2)),ceil(size(printmask,2)/4)+(1:ceil(size(printmask,2)/2))); else tempprintmask=printmask; end
                        temp=.75*mean(backgroundcolor2)+.25*temp;
                        temp(:,ceil(size(temp,2)/2+(1:size(temp,1))-size(temp,1)/2),:)=max(0,min(1, temp(:,ceil(size(temp,2)/2+(1:size(temp,1))-size(temp,1)/2),:)+(1-2*mean(backgroundcolor2))*repmat(.5*tempprintmask(round(linspace(1,size(tempprintmask,1),size(temp,1))),round(linspace(1,size(tempprintmask,2),size(temp,1)))),[1,1,3]) ));
                    end
                    str=get(DATA.handles(bp(n1)),'string'); set(DATA.handles(bp(n1)),'cdata',temp,'string',''); 
                end %uicontrol('units','pixel','position',pt.*[1 1 1 0]+[0 pt(4)-18-CONN_gui.font_offset 0 18+CONN_gui.font_offset],'style','text','string',str,'fontsize',9+CONN_gui.font_offset,'backgroundcolor','w','foregroundcolor','k'); end
                if ~isfield(DATA,'peakFDR')||~DATA.peakFDR, set(DATA.handles(6),'string',{'cluster-size','cluster-size p-FWE corrected','cluster-size p-FDR corrected','cluster-size p-uncorrected','peak-voxel p-FWE corrected','peak-voxel p-uncorrected','cluster-mass','cluster-mass p-FWE corrected','cluster-mass p-FDR corrected','cluster-mass p-uncorrected'}); end
                if DATA.thres{2}==3, set(DATA.handles(7),'string',sprintf('%-15s%13s%13s%13s%13s%13s%13s','Cluster (x,y,z)','size','peaks','TFCE','peak p-FWE','peak p-FDR','peak p-unc')); 
                elseif DATA.thres{2}==4, set(DATA.handles(7),'string',sprintf('%-15s%13s%13s%13s%13s%13s%13s','Cluster (x,y,z)','size','peaks','TFCE','peak p-FWE','peak p-FDR','peak p-unc')); 
                elseif DATA.parametric==1, set(DATA.handles(7),'string',sprintf('%-15s%13s%13s%13s%13s%13s%13s','Cluster (x,y,z)','size','size p-FWE','size p-FDR','size p-unc','peak p-FWE','peak p-unc'));
                elseif DATA.parametric==2, set(DATA.handles(7),'string',sprintf('%-15s%13s%13s%13s%13s%13s%13s%13s%13s','Cluster (x,y,z)','size','size p-FWE','size p-FDR','size p-unc','mass','mass p-FWE','mass p-FDR','mass p-unc'));
                end
                if nnz(DATA.paramoptions)<2, 
                    if DATA.paramoptions(1), set(DATA.handles(15),'string',{'parametric stats'},'value',1,'enable','off');
                    elseif DATA.paramoptions(2),  set(DATA.handles(15),'string',{'non-parametric stats'},'value',1,'enable','off');
                    else set(DATA.handles(15),'visible','off');
                    end
                end
                if DATA.issurface, set(DATA.handles([21 23 36]),'visible','off'); end
                ok=false;for method=1:numel(DATA.thres_defaults), if isequal(DATA.thres,DATA.thres_defaults{method}(1:4))&&DATA.side==3, ok=true; break; end; end; if ~ok, method=numel(DATA.thres_defaults)+1; end
                set(DATA.handles(32),'value',method);                
                icp=[1 2 3 4 5 6  11 15]; 
                if method<4
                    if DATA.thres{2}<3, set(DATA.handles(1),'string','voxel threshold: p <');
                    elseif (DATA.thres{2}==3||DATA.thres{2}==4), set(DATA.handles(1),'string','voxel threshold: p <');
                    else set(DATA.handles(1),'string',['voxel threshold: ',DATA.mat{6},' >']);
                    end
                    if method==3, set(DATA.handles(4),'string','cluster threshold: k >');
                    else set(DATA.handles(4),'string','cluster threshold: p <');
                    end
                end
                if method<4,%&&~get(DATA.handles(38),'value')
                    set(DATA.handles(icp),'visible','off');
                else
                    set(DATA.handles(icp),'visible','on');
                    if DATA.thres{2}==3||DATA.thres{2}==4, set(DATA.handles([4 5 6  15]),'visible','off'); end
                    set(DATA.handles(33),'visible','off','string','');
                end                
                hc1=uicontextmenu;
                uimenu(hc1,'Label','Export table','callback',@(varargin)conn_exportlist(DATA.handles(8),'',get(DATA.handles(7),'string')));
                set(DATA.handles(8),'uicontextmenu',hc1);
                hc1=uicontextmenu;
                uimenu(hc1,'Label','Export table','callback',@(varargin)conn_exportlist(DATA.handles(12)));
                uimenu(hc1,'Label','Change reference atlas','callback',{@conn_vproject,'ref_atlas'});
                set(DATA.handles(12),'uicontextmenu',hc1);
                %uicontrol('style','pushbutton','units','norm','position',[.775,.02,.1,.04],'string','Export table','fontname','arial','fontsize',8+CONN_gui.font_offset,'callback',@(varargin)conn_exportlist(DATA.handles(12)),'tooltipstring','Exports table'),...
                %hc1=uicontextmenu;
                %uimenu(hc1,'Label','Export stats','callback',{@conn_vproject,'export_stats'});
                %set(DATA.handles(8),'uicontextmenu',hc1);
                if voxeltovoxel
                    DATA.handles(10)=uicontrol('style','pushbutton','units','norm','position',[.875,.845,.1,.04],'string','post-hoc peak analyses','fontname','arial','fontsize',8+CONN_gui.font_offset,'callback',{@conn_vproject,'connectome'});
                end
                %if DATA.mat{6}~='T'||isempty(DATA.mat{2}), %enable Fcluster% 
                %    delete(DATA.handles(6)); DATA.handles(6)=uicontrol('style','popupmenu','units','norm','position',[.625,.895,.1,.03],'string',{'peak p-FDR corrected','extent k-value'},'fontname','arial','fontsize',8+CONN_gui.font_offset,'callback',{@conn_vproject,'fwec.clusterlevel'},'value',min(2,thres{4}));
                %end
                %if isempty(DATA.mat{2}), set(DATA.handles(6),'visible','off','value',4); DATA.thres{4}=4; DATA.thres{3}=10; thres=DATA.thres; set(DATA.handles(5),'string',num2str(thres{3})); set(DATA.handles(4),'string','extent threshold: cluster k ='); end
                DATA.axes=axes('units','norm','position',[.10,.46,.35,.33],'visible','off','tag','highlight_pointer');%[.55,.35,.4,.6],'visible','off');
                h0=get(0,'screensize');
                set(GCF,'name',['results explorer ',DATA.spmfile],'numbertitle','off','color',backgroundcolor,'units','pixels','position',[h0(3)-.75*h0(3)+2,h0(4)-.9*h0(4)-48,.75*h0(3)-2,.9*h0(4)],'interruptible','off','busyaction','cancel','colormap',map,'tag','conn_vproject','userdata',DATA,'windowbuttonmotionfcn',@uimotionfcn); %,'windowbuttondownfcn',{@conn_vproject,'buttondown'},'windowbuttonupfcn',{@conn_vproject,'buttonup'},'windowbuttonmotionfcn',[],'keypressfcn',{@conn_vproject,'keypress'},'colormap',map,'userdata',DATA); 
                figure(GCF);
            end
        case '3d',
            clf;
            DATA.handles=[uicontrol('style','text','units','norm','position',[0,0,.1,.05],'string','resolution','backgroundcolor','k','foregroundcolor','w'),...
                uicontrol('style','edit','units','norm','position',[.1,0,.1,.05],'string',num2str(DATA.res),'callback',{@conn_vproject,'resolution'},'backgroundcolor','k','foregroundcolor','w'),...
                uicontrol('style','text','units','norm','position',[.2,0,.1,.05],'string','threshold','backgroundcolor','k','foregroundcolor','w'),...
                uicontrol('style','edit','units','norm','position',[.3,0,.1,.05],'string',num2str(DATA.threshold),'callback',{@conn_vproject,'threshold'},'backgroundcolor','k','foregroundcolor','w')];
            DATA.axes=gca;set(DATA.axes,'visible','off');
            set(GCF,'name','Volume display','numbertitle','off','color','w','interruptible','off','busyaction','cancel','colormap',map,'userdata',DATA);%,'windowbuttondownfcn',{@conn_vproject,'buttondown'},'windowbuttonupfcn',{@conn_vproject,'buttonup'},'windowbuttonmotionfcn',[],'keypressfcn',{@conn_vproject,'keypress'},'colormap',map,'userdata',DATA);
        case 'orth',
            clf;
            DATA.handles=[uicontrol('style','text','units','norm','position',[0,0,.1,.05],'string','resolution','backgroundcolor','k','foregroundcolor','w'),...
                uicontrol('style','edit','units','norm','position',[.1,0,.1,.05],'string',num2str(DATA.res),'callback',{@conn_vproject,'resolution'},'backgroundcolor','k','foregroundcolor','w'),...
                uicontrol('style','text','units','norm','position',[.2,0,.1,.05],'string','threshold','backgroundcolor','k','foregroundcolor','w'),...
                uicontrol('style','edit','units','norm','position',[.3,0,.1,.05],'string',num2str(DATA.threshold),'callback',{@conn_vproject,'threshold'},'backgroundcolor','k','foregroundcolor','w')];
            DATA.axes=gca;set(DATA.axes,'visible','off');
            %DATA.axes=gca;
            set(GCF,'color','w','windowbuttondownfcn',[],'windowbuttonupfcn',[],'windowbuttonmotionfcn',[],'keypressfcn',[],'colormap',map,'userdata',DATA);
        case 'none',
            if ~nargout, DATA.axes=gca;set(DATA.axes,'visible','off'); end
            %DATA.axes=gca;
            %set(GCF,'color','w','userdata',DATA);
    end
end

% display
if strcmp(views,'full'),
    if init~=0, % redraws orth views
        hcontrols=findobj(GCF,'enable','on','-not','style','edit');
        hcontrols=hcontrols(ishandle(hcontrols));
        set(hcontrols,'enable','off');
        set(GCF,'pointer','watch'); 
        %if doredraw, drawnow; end
        if ~isempty(b),%length(threshold)>1,
            if isempty(DATA.clusters)||init~=-2, % compute/display stats
                bnew=b(:,:,:,1); %bnew(~a)=nan;%bnew(a<=threshold(1))=nan;
                [bnew,txt,xyzpeak,clusters,clustervol,clusternames,ithres,cluster_stats,doneparam,DATA.SIM,NiPERM,enew,fnew,gnew,statsOut]=vproject_display(bnew,threshold,thres,mat,DATA.side,DATA.d,DATA.e,DATA.f,DATA.g,DATA.peakFDR,DATA.parametric,DATA.SIM,DATA.spmfile,init,DATA.datatype,DATA.surfcoords,DATA.refatlas);
                b(:,:,:,2)=bnew; threshold(3)=0;%.5;
                %for n1=1:4, set(DATA.handles(8+n1),'string',strvcat(txt{n1})); end
                set(DATA.handles(8),'string',strvcat(txt,' '),'value',max(1,min(size(txt,1)+1,get(DATA.handles(8),'value'))));
                %if init>0, set(DATA.handles(8),'value',[]); end
                if isempty(NiPERM), set(DATA.handles(29),'string','','visible','off');
                else set(DATA.handles(29),'string',sprintf('(%d simulations)',NiPERM),'visible','on');
                end
                %if DATA.parametric==2&&~doneparam, set(DATA.handles(24),'visible','on');
                %else set(DATA.handles(24),'visible','off');
                %end
                %if DATA.parametric==2&&~doneparam, set(DATA.handles(24),'visible','on');
                %else set(DATA.handles(24),'visible','off');
                %end
                DATA.stdprojections=cell(1,4); 
                DATA.txt=txt;
                DATA.clusternames=clusternames;
                DATA.clustervol=clustervol;
                DATA.ithres=ithres;
                DATA.cluster_stats=cluster_stats;
                DATA.e=enew;
                DATA.f=fnew;
                DATA.g=gnew;
                DATA.cluster_stats_all=statsOut;
            else
                clusters=DATA.clusters;
                xyzpeak=DATA.xyzpeak;
                clusternames=DATA.clusternames;
                clustervol=DATA.clustervol;
            end
            if isfield(DATA,'selectcluster'), selectcluster=DATA.selectcluster;else, selectcluster=get(DATA.handles(8),'value'); end
            %if isempty(selectcluster), selectcluster=length(clusters)+1; end
            if all(selectcluster<=length(clusters)), select=cat(1,clusters{selectcluster}); clusternames=arrayfun(@(i)clusternames{i}.txt,selectcluster,'uni',0); clusternames=cat(2,clusternames{:}); titleclusternames=sprintf('Selected cluster %d/%d ',selectcluster,length(clusters));
            elseif (isempty(selectcluster)||any(selectcluster==length(clusters)+1))&&~isempty(clusternames), select=[]; clusternames=clusternames{length(clusters)+1}.txt; titleclusternames=''; %'All suprathreshold voxels:';
            else, select=[]; clusternames={}; titleclusternames=''; end
        end
        %M={[-1,0,0;0,0,-1;0,-1,0],[0,-1,0;0,0,-1;-1,0,0],[-1,0,0;0,1,0;0,0,1]};
        M={[1,0,0;0,0,-1;0,1,0],[0,-1,0;0,0,-1;1,0,0],[1,0,0;0,-1,0;0,0,-1]};
        if ~isfield(DATA,'stdprojections') || init==-1, DATA.stdprojections=cell(1,4); end
        tplot={};infoplot={};cplot={};
        for n1=1:3, 
            [tplot{n1},nill,nill,nill,b2,tres,DATA.stdprojections{n1},cplot{n1}]=vproject_view(a,b,c,clustervol,M{n1},threshold,res,select,DATA.stdprojections{n1},DATA.surfcoords);
            infoplot{n1}=struct('boundingbox',b2,'res',tres,'size',[size(tplot{n1},1),size(tplot{n1},2)],'projection',projection);
            %[tplot{n1},infoplot{n1},DATA.stdprojections{n1}]=conn_vproject(a,b,c,d,'none',M{n1},thres,res,0,select,mat,threshold,DATA.stdprojections{n1}); 
        end
        %tplot=[tplot{1},tplot{2};tplot{3},ones(size(tplot{3},1),size(tplot{2},2))];
        tplot=cat(1,cat(2,tplot{1},tplot{2}),cat(2,tplot{3},tplot{1}(1)*ones([size(tplot{3},1),size(tplot{2},2),size(tplot{2},3)])));
        cplot=cat(1,cat(2,cplot{1},cplot{2}),cat(2,cplot{3},cplot{1}(1)*ones([size(cplot{3},1),size(cplot{2},2),size(cplot{2},3)])));
        
        if ~isempty(b), h=DATA.axes; cla(h); %h=subplot(3,4,1,'parent',GCF); set(h,'units','norm','position',[.10,.46,.35,.33] ); % DATA.axes
        else, h=subplot(3,4,1,'parent',GCF); set(h,'units','norm','position',[.05,.1,.4,.8]);
        end
        %pres=.5;image(1:pres:size(tplot,2),1:pres:size(tplot,1),round(max((ceil(tplot(1:pres:end,1:pres:end)/64)-1)*64+1,convn(convn(tplot(1:pres:end,1:pres:end),conn_hanning(7)/4,'same'),conn_hanning(7)'/4,'same'))));axis equal; axis off;
        pres=.25;
        conf=conn_hanning(2/pres+1);conf=conf/sum(conf); 
        tplotidx1=1:pres:size(tplot,1); tplotidx1=[repmat(tplotidx1(1),1,(numel(conf)-1)/2), tplotidx1, repmat(tplotidx1(end),1,(numel(conf)-1)/2)];
        tplotidx2=1:pres:size(tplot,2); tplotidx2=[repmat(tplotidx2(1),1,(numel(conf)-1)/2), tplotidx2, repmat(tplotidx2(end),1,(numel(conf)-1)/2)];
        temp=tplot(round(tplotidx1),round(tplotidx2),:);
        if 0,%mean(backgroundcolor)<.5, 
            %tmask=find(all(diff(temp,1,3)==0,3));
            %for nt=1:3, temp(tmask+(nt-1)*size(temp,1)*size(temp,2))=(1-backgroundcolor(nt))+(2*backgroundcolor(nt)-1)*temp(tmask+(nt-1)*size(temp,1)*size(temp,2)); end
            %tmask=find(all(temp==1,3));
            %for nt=1:3, temp(tmask+(nt-1)*size(temp,1)*size(temp,2))=backgroundcolor(nt); end
            for nt=1:3, temp(:,:,nt)=(1-backgroundcolor(nt))+(2*backgroundcolor(nt)-1)*temp(:,:,nt); end
            temp=convn(convn(temp,conf,'valid'),conf','valid');
        else
            temp=.975*convn(convn(temp,conf,'valid'),conf','valid');
        end
        %backmask=find(all(temp==1,3));for n=1:3,temp(backmask+(n-1)*size(temp,1)*size(temp,2))=backgroundcolor(n); end
        hi=image(1:pres:size(tplot,2),1:pres:size(tplot,1),temp,'parent',h);axis(h,'equal','off');
        hold(h,'on'); hf=text(.95*size(tplot,2),.95*size(tplot,1),titleclusternames,'color','k','horizontalalignment','right','fontsize',5+CONN_gui.font_offset,'parent',h); hold(h,'off');
        set([hi,hf],'buttondownfcn',@conn_vproject_imcallback);
        %image(round(convn(convn(tplot(round(1:.5:end),round(1:.5:end),:),conn_hanning(3)/2,'same'),conn_hanning(3)'/2,'same')));axis equal; axis off;
        %image(tplot);axis equal; axis off;
        data1plot=DATA.stdprojections{4};
        if ~isempty(clusternames),set(DATA.handles(12),'string',clusternames,'value',[],'visible','on'); %set(DATA.handles(10),'string',titleclusternames); 
        else, set(DATA.handles(12),'string','','value',1,'visible','off'); %set(DATA.handles(10),'string',''); 
        end
        if (DATA.thres{2}==3||DATA.thres{2}==4), set(DATA.handles(13),'string',['TFCE > ',num2str(DATA.ithres(1),'%0.2f')]);
        elseif numel(DATA.mat{3})>1&&isequal(DATA.mat{6},'T')&&isequal(DATA.side,3), set(DATA.handles(13),'string',['|',DATA.mat{6},'(',num2str(DATA.mat{3}(end)),')|',' > ',num2str(DATA.ithres(1),'%0.2f')]);
        elseif numel(DATA.mat{3})>1&&isequal(DATA.mat{6},'T')&&isequal(DATA.side,2), set(DATA.handles(13),'string',['-',DATA.mat{6},'(',num2str(DATA.mat{3}(end)),')',' > ',num2str(DATA.ithres(1),'%0.2f')]);
        elseif numel(DATA.mat{3})>1&&isequal(DATA.mat{6},'T'), set(DATA.handles(13),'string',[DATA.mat{6},'(',num2str(DATA.mat{3}(end)),')',' > ',num2str(DATA.ithres(1),'%0.2f')]);
        elseif numel(DATA.mat{3})>1, set(DATA.handles(13),'string',[DATA.mat{6},'(',num2str(DATA.mat{3}(1)),',',num2str(DATA.mat{3}(2)),')',' > ',num2str(DATA.ithres(1),'%0.2f')]);
        else set(DATA.handles(13),'string',[DATA.mat{6},'(',num2str(DATA.mat{3}),')',' > ',num2str(DATA.ithres(1),'%0.2f')]);
        end
        if (DATA.thres{2}==3||DATA.thres{2}==4), set(DATA.handles(14),'string','');
        else set(DATA.handles(14),'string',['k ',8805,' ',num2str(DATA.ithres(2))]);
        end
        if DATA.paramoptions(1), set(DATA.handles(19),'visible','on');
        else set(DATA.handles(19),'visible','off');
        end
        tstr1=cellstr(get(DATA.handles(3),'string'));
        tstr2=cellstr(get(DATA.handles(6),'string'));
        if (DATA.thres{2}==3||DATA.thres{2}==4), 
            if DATA.thres{2}==3, tstr3='nonparametric statistics (Threshold Free Cluster Enhancement, Smith and Nichols 2007; Topological FDR, Chumbley et al. 2010)';
            else tstr3='nonparametric statistics (Threshold Free Cluster Enhancement, Smith and Nichols 2007)';
            end
            tstr4=sprintf('%s %s %s',get(DATA.handles(1),'string'),get(DATA.handles(2),'string'),tstr1{get(DATA.handles(3),'value')});
            set(DATA.handles(33),'string',{tstr3,sprintf('%s',tstr4)});
        else
            if parametric==1, tstr3='parametric statistics (Gaussian Random Field theory, Worsley et al. 1996)';
            else tstr3='nonparametric statistics (permutation/randomization analyses, Bullmore et al. 1999)';
            end
            tstr4=sprintf('%s %s %s',get(DATA.handles(1),'string'),get(DATA.handles(2),'string'),tstr1{get(DATA.handles(3),'value')});
            tstr5=sprintf('%s %s %s',get(DATA.handles(4),'string'),get(DATA.handles(5),'string'),tstr2{get(DATA.handles(6),'value')});
            set(DATA.handles(33),'string',{tstr3,sprintf('%s; %s',tstr5,tstr4)});
        end
        if DATA.thres{2}==3, set(DATA.handles(7),'string',sprintf('%-15s%13s%13s%13s%13s%13s%13s','Cluster (x,y,z)','size','peaks','TFCE','peak p-FWE','peak p-FDR','peak p-unc'));
        elseif DATA.thres{2}==4, set(DATA.handles(7),'string',sprintf('%-15s%13s%13s%13s%13s%13s%13s','Cluster (x,y,z)','size','peaks','TFCE','peak p-FWE','peak p-FDR','peak p-unc'));
        elseif DATA.parametric==1, set(DATA.handles(7),'string',sprintf('%-15s%13s%13s%13s%13s%13s%13s','Cluster (x,y,z)','size','size p-FWE','size p-FDR','size p-unc','peak p-FWE','peak p-unc'));
        elseif DATA.parametric==2, set(DATA.handles(7),'string',sprintf('%-15s%13s%13s%13s%13s%13s%13s%13s%13s','Cluster (x,y,z)','size','size p-FWE','size p-FDR','size p-unc','mass','mass p-FWE','mass p-FDR','mass p-unc'));
        end
        if ~isempty(b),
            DATA.b=b;DATA.threshold=threshold;DATA.xyzpeak=xyzpeak;DATA.infoplot=infoplot;DATA.select=select;DATA.clusters=clusters;DATA.cplot=cplot;
            %set(GCF,'userdata',DATA);
        end
        %views='none';
        set(hcontrols(ishandle(hcontrols)),'enable','on');
        set(GCF,'pointer','arrow');%,'userdata',DATA);
    end
    set(DATA.axes,'visible','off','tag','highlight_pointer');
    %subplot(122);conn_vproject(a,'none',projection,threshold,res,box); 
    %if ~isempty(b), DATA.axes=subplot(222); else, DATA.axes=subplot(122); end
end
if strcmp(views,'orth'),
%     %M={[-1,0,0;0,0,-1;0,-1,0],[0,-1,0;0,0,-1;-1,0,0],[-1,0,0;0,1,0;0,0,1]};
%     M={[1,0,0;0,0,-1;0,1,0],[0,1,0;0,0,-1;1,0,0],[1,0,0;0,1,0;0,0,-1]};
%     set(GCF,'pointer','watch'); drawnow
%     for n1=1:3, 
%         subplot(2,2,n1);
%         conn_vproject(a,b,c,d,'none',M{n1},threshold,res,box); 
%     end;
%     set(GCF,'pointer','arrow');
else, % none/full
    if ~nargout, set(GCF,'pointer','watch'); drawnow; end
    % volume plot
	[dataplot,idx,mask,B,b2,res,data1plot]=vproject_view(a,b,c,[],projection,threshold,res,select,data1plot,DATA.surfcoords);
    %if size(dataplot,3)>1, dataplot=dataplot(:,:,1)+dataplot(:,:,2); end
    infoplot=struct('boundingbox',b2,'res',res,'size',[size(dataplot,1),size(dataplot,2)],'projection',projection);
    DATA.infoplot{4}=infoplot;DATA.stdprojections{4}=data1plot; set(GCF,'userdata',DATA);
    if 0,%~nargout,  % note: do not display interactive single-display
        axes(DATA.axes);
        %image(dataplot);axis equal;axis off;
        pres=.25;conf=conn_hanning(2/pres+1);conf=conf/sum(conf); image(1:pres:size(dataplot,2),1:pres:size(dataplot,1),convn(convn(dataplot(round(1:pres:end),round(1:pres:end),:),conf,'same'),conf','same'));axis equal; axis off;
        %pres=.5;image(1:pres:size(dataplot,2),1:pres:size(dataplot,1),round(max((ceil(dataplot(1:pres:end,1:pres:end)/64)-1)*64+1,convn(convn(dataplot(1:pres:end,1:pres:end),conn_hanning(5)/3,'same'),conn_hanning(5)'/3,'same'))));axis equal; axis off;
        DATA.axes=gca;
        % bounding box
        if box,
            B3=B'*pinv(projection);
            B3=[1+(B3(:,1)-b2(1,1))/res, 1+(B3(:,2)-b2(1,2))/res, 1+(B3(:,3)-b2(1,3))/res];
            border={[1,3,7,5,1],[2,6,4,8,2],[1,3,2,8,1],[7,5,4,6,7]};
            [minBt3]=max(B3(:,3));
            for n0=1:length(border),
                for n1=1:length(border{n0})-1,
                    color=(B3(border{n0}(n1),3)==minBt3|B3(border{n0}(n1+1),3)==minBt3);
                    Bt1=linspace(B3(border{n0}(n1),1),B3(border{n0}(n1+1),1),100);
                    Bt2=linspace(B3(border{n0}(n1),2),B3(border{n0}(n1+1),2),100);
                    Bt3=linspace(B3(border{n0}(n1),3),B3(border{n0}(n1+1),3),100);
                    if color,
                        Bt1(Bt3>idx(max(1,min(prod(size(mask)), round(Bt2)+size(idx,1)*round(Bt1-1)))) & ...
                            mask(max(1,min(prod(size(mask)), round(Bt2)+size(idx,1)*round(Bt1-1)))) )=nan;
                    end
                    hold on; h=plot(Bt1,Bt2,'-'); set(h,'color',.25*[0,0,1],'linewidth',4-3*color);hold off;
                end
            end
        end
    end
    if ~nargout, set(GCF,'pointer','arrow'); end
end

    function conn_vproject_imcallback(varargin)
        tp1=get(0,'pointerlocation');
        tp2=get(GCF,'position');
        tp3=get(0,'screensize');
        tp2(1:2)=tp2(1:2)+tp3(1:2)-1; % note: fix issue when connecting to external monitor/projector
        tpos=(tp1-tp2(1:2));
        set(GCF,'currentpoint',tpos);
        if strcmp(get(gcbo,'type'),'axes'), th=gcbo;
        else th=get(gcbo,'parent');
        end
        if strcmp(get(th,'type'),'axes')
            xy=get(th,'currentpoint');
            xy=round(xy(1,1:2));
            if isfield(DATA,'cplot'),
                z=DATA.cplot(max(1,min(size(DATA.cplot,1), xy(2))),max(1,min(size(DATA.cplot,2),xy(1))));
                if z==0, z=[]; end
                zbak=get(DATA.handles(8),'value');
                if isequal(zbak,length(DATA.clusters)+1); zbak=[]; end
                if ~isequal(zbak,z)
                    set(DATA.handles(8),'value',z);
                    conn_vproject(GCF,[],'selectroi');
                end
            end
        end
    end
end

function [dataplot,idx,mask,B,b2,res,data1plot,clusteridplot]=vproject_view(a,b,c,clusterid,projection,threshold,res,select,data1plot,surfcoords)
%interpolate&project
if ~isempty(data1plot),
    mask=data1plot.mask;
    idx=data1plot.idx;
    trans=data1plot.trans;
    B=data1plot.B;
    b2=data1plot.b2;
else,
    if ~isempty(b), a=a+j*b(:,:,:,end); end
    if ~isempty(surfcoords)
        data1plot.recoords=round([surfcoords,ones(size(surfcoords,1),1)]*pinv([-2 0 0 92-10;0 2 0 -128+10;0 0 2 -74+20;0 0 0 1])');
        data1plot.recoords=data1plot.recoords(:,1:3);
        data1plot.recoordsmax=4+max(data1plot.recoords,[],1);
        data1plot.reindex=data1plot.recoords(:,1)+data1plot.recoordsmax(1)*(data1plot.recoords(:,2)-1)+data1plot.recoordsmax(1)*data1plot.recoordsmax(2)*(data1plot.recoords(:,3)-1);
        a=accumarray(data1plot.recoords,a(:),data1plot.recoordsmax,@mean);
        %b=accumarray([data1plot.recoords 1;data1plot.recoords 2],b(:));
    end
    bb={[1,size(a,2)]-(size(a,2)+1)/2,[1,size(a,2)]-(size(a,2)+1)/2,[1,size(a,3)]-(size(a,3)+1)/2};
    B=[];for n1=1:8,B=[B,[bb{1}(1+rem(n1,2));bb{2}(1+rem(floor(n1/2),2));bb{3}(1+rem(floor(n1/4),2))]]; end
    B2=B'*pinv(projection);
    b2=[min(B2,[],1);max(B2,[],1)];
    %res=res*((prod(b2(2,:)-b2(1,:))./prod(size(a))).^(1/3));
    if 1,%faster&avoids memory problems
        nz2=b2(1,3):res:b2(2,3);
        [x2,y2]=meshgrid(b2(1,1):res:b2(2,1),b2(1,2):res:b2(2,2)); % image plane index
        N=numel(x2);
        d2=[x2(:),y2(:)]*projection(1:2,:);
        d2=d2-repmat([bb{2}(1),bb{1}(1),bb{3}(1)],[N,1]); % 3d matrix index
        d2=reshape(d2,[size(x2),3]);
        trans=zeros(size(x2));
        idx=zeros([size(x2),1+(length(threshold)>1)]);
        idx1=1:N;idx2=[];
        for n1=1:length(nz2), % depth
            x3=round(d2(idx1)+nz2(n1)*projection(3,1));
            y3=round(d2(idx1+N)+nz2(n1)*projection(3,2));
            z3=round(d2(idx1+N*2)+nz2(n1)*projection(3,3));
            %a3=interp3(x1,y1,z1,a,x3,y3,z3,'nearest');
            %idxs=max(1,min(size(a,1), 1+round(y3-bb{2}(1)) )) + size(a,1)*max(0,min(size(a,2)-1, round(x3-bb{1}(1)) )) + size(a,1)*size(a,2)*max(0,min(size(a,3)-1, round(z3-bb{3}(1)) ));
            idxs=max(1,min(size(a,1), 1+y3 )) + size(a,1)*max(0,min(size(a,2)-1, x3 )) + size(a,1)*size(a,2)*max(0,min(size(a,3)-1, z3 ));
            a3=a(idxs);
            idx0=real(a3)>threshold(1)|imag(a3)>threshold(end);
            idx(idx1(idx0))=(length(nz2)+1-n1)/length(nz2);
            idx2=[idx2,idx1(idx0)];
            idx1=idx1(~idx0);
            if length(threshold)>1,
                x3=d2(idx2)+nz2(n1)*projection(3,1);
                y3=d2(idx2+N)+nz2(n1)*projection(3,2);
                z3=d2(idx2+N*2)+nz2(n1)*projection(3,3);
                %a3=interp3(x1,y1,z1,a,x3,y3,z3,'nearest');
                %idxt=max(1,min(size(a,1), 1+round(y3-bb{2}(1)) )) + size(a,1)*max(0,min(size(a,2)-1, round(x3-bb{1}(1)) )) + size(a,1)*size(a,2)*max(0,min(size(a,3)-1, round(z3-bb{3}(1)) ));
                idxt=max(1,min(size(a,1), 1+round(y3) )) + size(a,1)*max(0,min(size(a,2)-1, round(x3) )) + size(a,1)*size(a,2)*max(0,min(size(a,3)-1, round(z3) ));
                a3=a(idxt);
                if 0,
                    idx00=find(imag(a3)>threshold(end));
                    idx(N+idx2(idx00))=max(idx(N+idx2(idx00)),imag(a3(idx00))); %n1;
                    idxtrans=find(trans(idx2(idx00))==0);
                    trans(idx2(idx00(idxtrans)))=idxt(idx00(idxtrans));
                else,
                    idx00=find(imag(a3)>threshold(end));
                    idxtrans=find(~idx(N+idx2(idx00)));%idxtrans=find(imag(a3(idx00))>=idx(N+idx2(idx00)));
                    idx(N+idx2(idx00(idxtrans)))=(length(nz2)+1-n1)/length(nz2);%imag(a3(idx00(idxtrans)));
                    trans(idx2(idx00(idxtrans)))=idxt(idx00(idxtrans));
                end
                
                %idx2=idx2(~idx00);
            end
        end
        mask=(idx>0);
        %if any(any(mask(:,:,1)>0)), idx(~mask(:,:,1))=5*max(idx(mask(:,:,1)>0)); else, idx(~mask(:,:,1))=length(nz2); end
    else,
        [x2,y2,z2]=meshgrid(b2(1,1):res:b2(2,1),b2(1,2):res:b2(2,2),b2(1,3):res:b2(2,3));
        d=[x2(:),y2(:),z2(:)]*projection;
        x2=reshape(d(:,1),size(x2));
        y2=reshape(d(:,2),size(y2));
        z2=reshape(d(:,3),size(z2));
        %a2=interp3(x1,y1,z1,a,x2,y2,z2,'nearest');
        a2=a(max(1,min(size(a,1), 1+round(y2-bb{2}(1)) )) + size(a,1)*max(0,min(size(a,2)-1, round(x2-bb{1}(1)) )) + size(a,1)*size(a,2)*max(0,min(size(a,3)-1, round(z2-bb{3}(1)) )));
        [mask,idx]=max(a2>threshold,[],3);idx(~mask)=max(idx(mask>0));
    end
    data1plot.mask=mask;
    data1plot.idx=idx;
    data1plot.trans=trans;
    data1plot.B=B;
    data1plot.b2=b2;
end
conf1=conn_hanning(3)/2;conf1=conf1/sum(conf1);
conf2=conn_hanning(3)/2;conf2=conf2/sum(conf2);
if size(idx,3)==1, dataplot=-idx; dataplot=1+63*(dataplot-min(dataplot(:)))/(max(dataplot(:))-min(dataplot(:)));
else, 
    %dataplot1=-idx(:,:,1); dataplot1=1+63*(dataplot1-min(dataplot1(:)))/(max(dataplot1(:))-min(dataplot1(:)));
    %conjmask=mask(:,:,1)&mask(:,:,2);
    %dataplot2=-idx(:,:,2); dataplot2=1+63*(dataplot2-min(dataplot2(:)))/(max(dataplot2(:))-min(dataplot2(:)));
    %dataplot=dataplot1; dataplot(conjmask)=64+dataplot2(conjmask);
    if 1,
        dataplot1=idx(:,:,1);
        dataplot2=idx(:,:,2);
        %dataplot1=convn(convn(dataplot1,conf1,'same'),conf1','same');
        %dataplot2=convn(convn(dataplot2,conf2,'same'),conf2','same');
        k=.75;%.2;
        dataplotA=(dataplot2==0).*(k+(1-k)*dataplot1.^4) + (dataplot2>0).*(k+(1-k)*dataplot1.^4);
        dataplotA(dataplot1==0)=1; %.126; % 1 % background color
        dataplotB=(dataplot2==0).*dataplotA + (dataplot2>0).*(dataplotA.*max(0,min(.8,.4+.4*tanh(5*(dataplot1-dataplot2)))));
        %dataplotB=(dataplot2==0).*dataplotA + (dataplot2>0).*(dataplotA.*max(0,min(.75,tanh(1.5*(1*dataplot1-dataplot2)))));
        %dataplotB=(dataplot2==0).*dataplotA + (dataplot2>0).*(.25*dataplotA.*(dataplot2<=.9*dataplot1));
        %dataplotB=(dataplot2==0).*dataplotA + (dataplot2>0).*(0*(dataplot2>.95*dataplot1) + .5*dataplotA.*(dataplot2<=.95*dataplot1));
        dataplot=cat(3,dataplotA,dataplotB,dataplotB);
    else,
        dataplot1=idx(:,:,1);
        dataplot1=0+1*(dataplot1-min(dataplot1(~isinf(dataplot1))))/(max(dataplot1(~isinf(dataplot1)))-min(dataplot1(~isinf(dataplot1))));
        dataplot2=idx(:,:,2); if any(dataplot2(:)>0),mindataplot=.95*min(dataplot2(dataplot2>0));else,mindataplot=0;end;
        %dataplot2=dataplot2>mindataplot;
        dataplot2=max(0,dataplot2-mindataplot);
        dataplot2=convn(convn(dataplot2,conf1,'same'),conf1','same');
        dataplot2=0+1*dataplot2/max(eps,max(dataplot2(:)));
        %dataplot1=max(prctile(dataplot1(dataplot1<max(dataplot1(:))),5),dataplot1);dataplot2=max(0,dataplot1-dataplot2);
        %dataplot1=.75+.25*(dataplot1).^2;dataplot2=max(0,dataplot1-dataplot2);
        dataplot1=convn(convn(max(0,min(1,dataplot1)),conf2,'same'),conf2','same');
        dataplot1=0+1*(dataplot1).^4;
        %dataplotA=max(0,min(1, dataplot1+.25*(dataplot2>0)));dataplotB=max(0,min(1, dataplot1-dataplot2+.0*(dataplot2>0)));
        dataplotA=max(0,min(1, dataplot1.*(dataplot2==0)+(.25+dataplot1).*(dataplot2>0)));
        dataplotB=max(0,min(1, dataplot1.*(dataplot2==0)+(.25+dataplot1).*(.5-.5*dataplot2).*(dataplot2>0)));
        dataplot=cat(3,dataplotA,dataplotB,dataplotB);
        %idxback=find(all(dataplot<.1,3)); dataplot(idxback)=1;dataplot(idxback+size(dataplot,1)*size(dataplot,2))=1;dataplot(idxback+2*size(dataplot,1)*size(dataplot,2))=1;
        %dataplot=cat(3,dataplot1,dataplot2,dataplot2);
    end
end
if ~isempty(select),
    if ~isempty(surfcoords), select=data1plot.reindex(select); end
    [transa,transb,transc]=unique(trans(:));%trans=a(c);    
    [nill,idxp]=intersect(transa(:),select(:));
    d=logical(zeros(size(transa)));d(idxp)=1;
    e=find(d(transc));
    %%dataplot(e)=dataplot(e)+64;
    dataplot3=zeros([size(dataplot,1),size(dataplot,2)]);dataplot3(e)=1;
    mask0=.80*dataplot(:,:,1)+.20*dataplot(:,:,2);
    dataplot(:,:,1)=dataplot(:,:,1).*dataplot3 + mask0.*(1-dataplot3);
    dataplot(:,:,2)=0*dataplot(:,:,2).*dataplot3 + mask0.*(1-dataplot3);
    dataplot(:,:,3)=0*dataplot(:,:,3).*dataplot3 + mask0.*(1-dataplot3);
    %dataplot3=convn(convn(dataplot3,conf1,'same'),conf1,'same');e=find(dataplot3>0);
    %n=size(dataplot,1)*size(dataplot,2);
    %dataplot(e+1*n)=dataplot(e+1*n).*(1-dataplot3(e))+.5*dataplot(e+1*n).*(dataplot3(e));
    %dataplot(e+2*n)=dataplot(e+2*n).*(1-dataplot3(e))+.5*dataplot(e+2*n).*(dataplot3(e));
    %dataplot(e+0*n)=dataplot(e+0*n).*(1-dataplot3(e))+.5*dataplot(e+0*n).*(dataplot3(e));
end
if size(idx,3)~=1,
    if ~isempty(surfcoords), c=accumarray(data1plot.recoords,c(:),data1plot.recoordsmax,@mean); end
    n=size(dataplot,1)*size(dataplot,2);
    idxtemp1=find(dataplot2>0);
    idxtemp2=find(c(trans(idxtemp1))>.5);
    idxtemp=idxtemp1(idxtemp2);
    temp3=dataplot(idxtemp+2*n); 
    dataplot(idxtemp+2*n)=dataplot(idxtemp+0*n); 
    dataplot(idxtemp+0*n)=temp3;
    clusteridplot=zeros([size(dataplot,1),size(dataplot,2)]);
    if ~isempty(clusterid),
        if ~isempty(surfcoords), clusterid=accumarray(data1plot.recoords,clusterid(:),data1plot.recoordsmax,@min); end
        idxtemp1=find(dataplot2>0);
        idxtemp2=find(clusterid(trans(idxtemp1)));
        idxtemp=idxtemp1(idxtemp2);
        clusteridplot(idxtemp)=clusterid(trans(idxtemp1(idxtemp2)));
    end
end
end

function [anew,txt,xyzpeak,clusters,clustervol,clusternames,ithres,cluster_stats,doneall,SIM,NiPERM,tfceZ,tfceZpeaks,tfceZd,statsOut]=vproject_display(a,threshold,thres,mat,side,Z,tfceZ,tfceZpeaks,tfceZd,peakfdr,parametric,SIM,spmfile,init,datatype,surfcoords,refatlas)
global CONN_gui;
persistent refsrois;
if ~nargin, refsrois=[]; return; end
DOSCALEF=true;
if isempty(refsrois)||(isfield(refsrois,'filename')&&~isempty(refatlas)&&~isequal(refsrois.filename,refatlas))
    if isempty(refatlas)&&isfield(CONN_gui,'refs')&&isfield(CONN_gui.refs,'rois')&&isfield(CONN_gui.refs.rois,'filename')&&~isempty(CONN_gui.refs.rois.filename),
        refsrois=CONN_gui.refs.rois;
    else
        %filename=fullfile(fileparts(which('conn')),'utils','otherrois','Td.img');
        if isempty(refatlas), filename=fullfile(fileparts(which('conn')),'rois','atlas.nii');
        else filename=refatlas;
        end
        %filename=fullfile(fileparts(which('conn')),'utils','otherrois','BA.img');
        [filename_path,filename_name,filename_ext]=fileparts(filename);
        V=conn_fileutils('spm_vol',filename);
        [idxlabels,strlabels]=rex(filename,filename,'level','clusters','disregard_zeros',false); strlabels=regexprep(strlabels,['^',filename_name,'\.'],''); 
        %strlabels=textread(fullfile(filename_path,[filename_name,'.txt']),'%s','delimiter','\n');
        tempdata=conn_fileutils('spm_read_vols',V); if numel(V)>1, [nill,tempdata]=max(tempdata,[],4); tempdata(~nill)=0; idxlabels=1:numel(strlabels); end
        refsrois=struct('filename',filename,'filenameshort',filename_name,'V',V,'data',tempdata,'labels',{strlabels},'labelsidx',full(sparse(1,round(idxlabels(:)'),1:numel(idxlabels))));
    end
    if isempty(surfcoords)||conn_surf_dimscheck(refsrois.V)
        [x,y,z]=ndgrid(1:size(a,1),1:size(a,2),1:size(a,3));xyz=[x(:),y(:),z(:)];
    else
        xyz=round(surfcoords);
    end
    if ~isfield(refsrois,'labelsidx'), refsrois.labelsidx=1:numel(refsrois.labels); end
    refsrois.data=spm_get_data(refsrois.V,pinv(refsrois.V(1).mat)*mat{1}*[xyz,ones(size(xyz,1),1)]');
    if numel(refsrois.V)>1, [nill,refsrois.data]=max(refsrois.data,[],1); refsrois.data(~nill)=0; end
    maxrefsrois=max(refsrois.data(:)); refsrois.count=hist(refsrois.data(:),0:maxrefsrois);
    refsrois.labels=[{'not-labeled'},refsrois.labels(:)'];
    refsrois.labelsidx=[1,1+refsrois.labelsidx(:)'];
end

if nargin<5, side=1; end
if nargin<4, mat={eye(4),numel(a)}; end
idx0=find(~isnan(a));
p=exp(-a(idx0));
p0=a;p0(idx0)=p;
P=a;P(idx0)=conn_fdr(p(:));
ithres=[nan,nan];
doneall=true;

if thres{2}==3||thres{2}==4
    if isempty(tfceZ)
        if mat{6}=='T'
            [tfceZpos,tfceZpospeaks,tfceZposd]=conn_tfce(abs(Z).*(Z>0),'Hmin',1,'adjacency',datatype); %T stat positive-sided
            [tfceZneg,tfceZnegpeaks,tfceZnegd]=conn_tfce(abs(Z).*(Z<0),'Hmin',1,'adjacency',datatype); %T stat negative-sided
            tfceZ=tfceZpos.*(Z>0)-tfceZneg.*(Z<0);
            tfceZpeaks=(tfceZpospeaks&(Z>0))|(tfceZnegpeaks&(Z<0));
            tfceZd=tfceZposd.*(Z>0)-tfceZnegd.*(Z<0);
        else
            [tfceZ,tfceZpeaks,tfceZd]=conn_tfce(sqrt(abs(Z)),'Hmin',1,'adjacency',datatype); %F stat (keep same scale as two-tailed T)
        end
        %tfceZclusters=conn_watershed(abs(Z),'minH',1);
        %tfceZclusters=tfceZd;
    end
end

NiPERM=[]; 
tZ=Z;
switch(thres{2}),
    case 1,%'vox-unc',
        % voxels above height threshold
        maskvoxels=a>-log(thres{1});
        idxvoxels=find(maskvoxels);
        anew=zeros(size(a));anew(idxvoxels)=a(idxvoxels)+log(thres{1});%1;
        [xt,yt,zt]=ind2sub(size(a),idxvoxels);
        if ~isempty(idxvoxels),
            if mat{6}=='T'&&max(Z(:))>0&&min(Z(:))<0,
                [c1,C1]=conn_clusters(maskvoxels&Z>0,datatype); C1=C1(idxvoxels)';
                [c2,C2]=conn_clusters(maskvoxels&Z<0,datatype); C2=C2(idxvoxels)';
                C=C1;
                C(C1==0&C2>0)=C2(C1==0&C2>0)+numel(c1);
                c=[reshape(c1,[],1);reshape(c2,[],1)]';
            else
                [c,C]=conn_clusters(maskvoxels,datatype); C=C(idxvoxels)'; c=c';
                %[xt,yt,zt]=ind2sub(size(a),idxvoxels); C=spm_clusters([xt,yt,zt]'); c=hist(C,1:max(C));
            end
        end
        if mat{6}=='T'&&side<3, %T stat one-sided
            u=spm_invTcdf(1-thres{1},mat{3}(end));
            ithres(1)=u;
        elseif mat{6}=='T'&&side==3, %two-sided T stat
            u=spm_invFcdf(1-thres{1},1,mat{3}(end));
            ithres(1)=sqrt(u);
        elseif mat{6}=='F' %F stat
            u=spm_invFcdf(1-thres{1},mat{3}(1),mat{3}(2));
            ithres(1)=u;
        elseif mat{6}=='X' %X2 stat
            u=spm_invXcdf(1-thres{1},mat{3});
            ithres(1)=u;
        end
    case 2,%'vox-FDR',
        maskvoxels=P<thres{1};
        idxvoxels=find(maskvoxels);
        anew=zeros(size(a));anew(idxvoxels)=-log(P(idxvoxels))+log(thres{1});%1;
        [xt,yt,zt]=ind2sub(size(a),idxvoxels);
        if ~isempty(idxvoxels),
            if mat{6}=='T'&&max(Z(:))>0&&min(Z(:))<0,
                [c1,C1]=conn_clusters(maskvoxels&Z>0,datatype); C1=C1(idxvoxels)';
                [c2,C2]=conn_clusters(maskvoxels&Z<0,datatype); C2=C2(idxvoxels)';
                C=C1;
                C(C1==0&C2>0)=C2(C1==0&C2>0)+numel(c1);
                c=[reshape(c1,[],1);reshape(c2,[],1)]';
            else
                [c,C]=conn_clusters(maskvoxels,datatype); C=C(idxvoxels)'; c=c';
                %[xt,yt,zt]=ind2sub(size(a),idxvoxels);C=spm_clusters([xt,yt,zt]');c=hist(C,1:max(C)); % C: cluster per voxel; c: #voxels per clusters
            end
        end
        if isempty(idxvoxels), u=1; ithres(1)=inf;
        else 
            if mat{6}=='T'&&side<3, %T stat
                u=spm_invTcdf(1-exp(-min(a(idxvoxels))),mat{3}(end)); ithres(1)=u;
            elseif mat{6}=='T'&&side==3, %two-sided T stat
                u=spm_invFcdf(1-exp(-min(a(idxvoxels))),1,mat{3}(end)); ithres(1)=sqrt(u);
            elseif mat{6}=='F' %F stat
                u=spm_invFcdf(1-exp(-min(a(idxvoxels))),mat{3}(1),mat{3}(2)); ithres(1)=u;
            elseif mat{6}=='X' %X2 stat
                u=spm_invXcdf(1-exp(-min(a(idxvoxels))),mat{3}); ithres(1)=u;
            end
        end
    case {3,4,5},%'T/F/X stat',
        if thres{2}==3||thres{2}==4,
            if thres{2}==3, tZ=tfceZd; % peak TFCE
            else tZ=tfceZd;
            end
            THR_TYPE=5;%'TFCE',
            THR=0;%thres{1};
            SIDE=side;
            if isempty(SIM)||~any(SIM.Pthr==THR&SIM.Pthr_type==THR_TYPE&SIM.Pthr_side==SIDE)
                if conn_existfile(conn_vproject_simfilename(spmfile,THR_TYPE,THR))||conn_vproject_randomise(spmfile,THR_TYPE,THR,init~=1,CONN_gui.isremote),
                    try,
                        SIM=conn_loadmatfile(conn_vproject_simfilename(spmfile,THR_TYPE,THR));
                        if ~isfield(SIM,'VERSION'), SIM.VERSION=0; end
                        if SIM.VERSION<2, SIM=[]; conn_fileutils('spm_unlink',conn_vproject_simfilename(spmfile,THR_TYPE,THR)); end % note: disregard older versions
                    end
                end
            end
            if ~isempty(SIM)&&any(SIM.Pthr==THR&SIM.Pthr_type==THR_TYPE&SIM.Pthr_side==SIDE), iPERM=find(SIM.Pthr==THR&SIM.Pthr_type==THR_TYPE&SIM.Pthr_side==SIDE);
            else iPERM=[];
            end
            if mat{6}~='T'|side==3, ttZ=abs(tfceZ); tfce_peak_mask=tfceZpeaks; ttZb=abs(tfceZd); %tfce_cluster_mask=tfceZclusters;
            elseif side==1,         ttZ=tfceZ; tfce_peak_mask=tfceZpeaks&(tfceZ>0); ttZb=tfceZd; %tfce_cluster_mask=tfceZclusters.*(tZ>0);
            elseif side==2,         ttZ=-tfceZ; tfce_peak_mask=tfceZpeaks&(tfceZ<0); ttZb=-tfceZd; %tfce_cluster_mask=tfceZclusters.*(tZ<0);
            else error('incorrect option ''side''');
            end
            if ~isempty(iPERM)
                tfce_max_dist=[];
                tfce_peak_p_unc=0;
                NiPERM=0;
                for niPERM=1:numel(iPERM),
                    iperm=iPERM(niPERM);
                    kiPERM=size(SIM.Dist_Voxel_statmax{iperm},1);
                    tfce_max_dist=[tfce_max_dist;SIM.Dist_Voxel_statmax{iperm}];                    
                    if nnz(SIM.Hist_Voxel_stat{iperm})<2, tfce_peak_p_unc=tfce_peak_p_unc+kiPERM*double(1+round(SIM.maxT*100*ttZb(tfce_peak_mask)/sqrt(prod(SIM.model.dims)/prod(SIM.model.dims(1:2))))<=find(SIM.Hist_Voxel_stat{iperm}));
                    else tfce_peak_p_unc=tfce_peak_p_unc+kiPERM*max(0,min(1,interp1(find(SIM.Hist_Voxel_stat{iperm}),flipud(cumsum(flipud(nonzeros(SIM.Hist_Voxel_stat{iperm})))),1+round(SIM.maxT*100*ttZb(tfce_peak_mask)/sqrt(prod(SIM.model.dims)/prod(SIM.model.dims(1:2)))),'linear','extrap')));
                    end
                    NiPERM=NiPERM+kiPERM;
                % note: Hist_Voxel_stat*sqrt(prod(SIM.model.dims)/prod(SIM.model.dims(1:2)))
                end
                tfce_peak_p_unc=tfce_peak_p_unc/NiPERM; % peak p-unc
                tfce_peak_p_FDR=conn_fdr(tfce_peak_p_unc); % peak p-FDR
                %idx=find(tfce_peak_mask);
                %tfce_cluster_borders=tfce_cluster_mask==0;
                %tfce_cluster_mask=tfce_cluster_mask>0&ismember(tfce_cluster_mask,tfce_cluster_mask(idx(tfce_peak_p_FDR<=thres{1})));
                if thres{2}==3, % TFCE FDR  
                    temp0=ttZb(tfce_peak_mask);
                    temp1=max(temp0(tfce_peak_p_FDR>thres{1})); if isempty(temp1), temp1=0; end
                    temp2=min(temp0(tfce_peak_p_FDR<=thres{1})); if isempty(temp2), temp2=inf; end
                    thres{1}=(temp1+temp2)/2;
                elseif thres{1}>1, thres{1}=0; 
                else thres{1}=interp1([0 ((1:numel(tfce_max_dist))-0.5)/numel(tfce_max_dist) 1],[min(tfce_max_dist),reshape(sort(tfce_max_dist),1,[]),max(tfce_max_dist)],1-max(0,min(1,thres{1}))); % TFCE FWE  
                end
            else thres{1}=0; 
            end
        end
        % voxels above height threshold
        if mat{6}=='T'&&side==2, maskvoxels=tZ<-thres{1}; %T stat negative-sided
        elseif mat{6}=='T'&&side==3, maskvoxels=abs(tZ)>thres{1}; %T stat two-sided
        else maskvoxels=tZ>thres{1}; 
        end
        %if thres{2}==3||thres{2}==4, maskvoxels=maskvoxels&tfce_cluster_borders==0; end % note: force peak-level clusters
        idxvoxels=find(maskvoxels); 
        anew=zeros(size(a));anew(idxvoxels)=a(idxvoxels)-min(a(idxvoxels))+.01;
        [xt,yt,zt]=ind2sub(size(a),idxvoxels);
        if ~isempty(idxvoxels),
            if mat{6}=='T'&&max(Z(:))>0&&min(Z(:))<0,
                [c1,C1]=conn_clusters(maskvoxels&Z>0,datatype); C1=C1(idxvoxels)';
                [c2,C2]=conn_clusters(maskvoxels&Z<0,datatype); C2=C2(idxvoxels)';
                C=C1;
                C(C1==0&C2>0)=C2(C1==0&C2>0)+numel(c1);
                c=[reshape(c1,[],1);reshape(c2,[],1)]';
            else
                [c,C]=conn_clusters(maskvoxels,datatype); C=C(idxvoxels)'; c=c';
                %[xt,yt,zt]=ind2sub(size(a),idxvoxels); C=spm_clusters([xt,yt,zt]'); c=hist(C,1:max(C));
            end
        end
        u=thres{1};
        if mat{6}=='T'&&side<3, %T stat
            ithres(1)=thres{1};
        elseif mat{6}=='T'&&side==3, %two-sided T stat
            ithres(1)=abs(u);
            u=abs(u).^2;
        else %F/X stat
            ithres(1)=u;
        end
end
% if mat{6}=='T'&&side>2, %T stat two-sided
%     ithres(1)=spm_invTcdf(1-(1-spm_Tcdf(ithres(1),mat{3}(2)))/2,mat{3}(2));
% end

txt=[];xyzpeak={};clusters={};clusternames={};cluster_stats={};clustervol=zeros(size(a));statsOut=[];
if ~isempty(idxvoxels),
    if isempty(surfcoords), xyz=[xt,yt,zt];
    else xyz=round(surfcoords(idxvoxels,:)); 
    end
    xyzrois=xyz;
    try, if isempty(surfcoords)||conn_surf_dimscheck(refsrois.V), xyzrois=[xt,yt,zt]; end; end
    idx1=find(c>0);
    idxn=[];idx2={};idxz=[];for n1=1:length(idx1),idx2{n1}=find(C==(idx1(n1))); idxn(n1)=length(idx2{n1}); idxz(n1)=max(abs(tZ(idxvoxels(idx2{n1}))));end
    if thres{2}==3||thres{2}==4, [nill,sidxn]=sort(idxz(:));sidxn=flipud(sidxn);
    else [nill,sidxn]=sort(idxn(:));sidxn=flipud(sidxn);
    end
    xyzpeak=cell(1,length(idx1));
    clusters=cell(1,length(idx1));
    [k,cp,cP,kmass,cpmass,cPmass,kscore,cpscore,cPscore,minP,minp,maxZ,pPFWE,pPFDR,pPunc,peaks]=deal(nan(1,length(idx1)));
    if ~peakfdr, thres4=thres{4}+(thres{4}>=6);
    else thres4=thres{4};
    end
    %k=zeros(1,length(idx1));cp=k;cP=k;Ez=k;cPFDR=k; minP=k;minp=k;maxZ=k;pPFWE=k;
    if thres{2}~=3&&thres{2}~=4&&parametric==2 % non-parametric stats
        switch(thres{2}),
            case 1, THR_TYPE=1; %'vox-unc',
            case 2, THR_TYPE=3; %'fdr-all'
            case 3, THR_TYPE=5;%'TFCE',
            case 4, THR_TYPE=5;%'TFCE',
            case 5, THR_TYPE=4;%'T/F/X stat',
        end
        THR=thres{1};
        SIDE=side;
        if isempty(SIM)||~any(SIM.Pthr==THR&SIM.Pthr_type==THR_TYPE&SIM.Pthr_side==SIDE)
            if conn_existfile(conn_vproject_simfilename(spmfile,THR_TYPE,THR))||((thres4~=1&&thres4~=8)&&conn_vproject_randomise(spmfile,THR_TYPE,THR,init~=1,CONN_gui.isremote)),
                try, 
                    SIM=conn_loadmatfile(conn_vproject_simfilename(spmfile,THR_TYPE,THR)); 
                    if ~isfield(SIM,'VERSION'), SIM.VERSION=0; end
                    if SIM.VERSION<2, SIM=[]; conn_fileutils('spm_unlink',conn_vproject_simfilename(spmfile,THR_TYPE,THR)); end % note: disregard older versions
                end
            end
        end
        if ~isempty(SIM)&&any(SIM.Pthr==THR&SIM.Pthr_type==THR_TYPE&SIM.Pthr_side==SIDE), iPERM=find(SIM.Pthr==THR&SIM.Pthr_type==THR_TYPE&SIM.Pthr_side==SIDE);
        else iPERM=[];
        end
    end
    for n1=1:length(idx1),
        k(n1)=idxn(sidxn(n1));
        clusters{n1}=idxvoxels(idx2{sidxn(n1)});
        clustervol(clusters{n1})=n1;
        if mat{6}~='T'|side==3, [nill,idxminp]=max(abs(tZ(clusters{n1}))); 
        elseif side==1, [nill,idxminp]=max(tZ(clusters{n1}));
        elseif side==2, [nill,idxminp]=max(-tZ(clusters{n1}));
        end
        minp(n1)=p0(clusters{n1}(idxminp));
        %[minp(n1),idxminp]=min(p0(clusters{n1}));
        minP(n1)=P(clusters{n1}(idxminp));
        maxZ(n1)=tZ(clusters{n1}(idxminp));
        if side==2, maxZ(n1)=-maxZ(n1);
        elseif side==3, maxZ(n1)=abs(maxZ(n1)); 
        end
        if ~isempty(SIM)&&SIM.VERSION<1, DOSCALEF=false; end
        if ~DOSCALEF||mat{6}~='T', kmass(n1)=sum(abs(Z(clusters{n1})));
        else kmass(n1)=sum(abs(Z(clusters{n1})).^2);
        end
        if ~DOSCALEF, kscore(n1)=kmass(n1)^3/max(eps,k(n1))^2.5/3;
        else kscore(n1)=kmass(n1)^(3/2)/max(eps,k(n1))/3;
        end
        %kmass(n1)=sum(abs(Z(clusters{n1}))-ithres(1));
        %[minP(n1),idxminp]=min(P(clusters{n1}));
        %minp(n1)=exp(-max(a(clusters{n1})));
        xyzpeak{n1}=xyz(idx2{sidxn(n1)}(idxminp),:);
        if thres{2}~=3&&thres{2}~=4&&parametric==1 % parametric stats
            if ~isempty(mat{2}),
                if mat{6}=='T'&&side<3, %T stat
                    try
                        if isempty(mat{5}), [cP(n1),cp(n1)]=spm_P(1,k(n1)*mat{2}(end)/mat{4},u,mat{3},'T',mat{2},1,mat{4});
                        else, [cP(n1),cp(n1)]=spm_P(1,k(n1)*mat{5},u,mat{3},'T',mat{2},1,mat{4});
                            if isnan(cP(n1))&&~cp(n1),cP(n1)=0; end
                        end
                    catch
                        cP(n1)=nan; cp(n1)=nan;
                    end
                    try
                        pPFWE(n1)=spm_P(1,0,maxZ(n1),mat{3},'T',mat{2},1,mat{4});
                    catch
                        pPFWE(n1)=nan;
                    end
                elseif mat{6}=='T'|mat{6}=='F', % two-sided T-stat or F-stat
                    ok=false;
                    if isdeployed||min([12, str2num(regexprep(spm('ver'),'[^\d]',''))])>8,
                        try
                            if isempty(mat{5}), [cP(n1),cp(n1)]=spm_P(1,k(n1)*mat{2}(end)/mat{4},u,mat{3},'F',mat{2},1,mat{4});
                            else, [cP(n1),cp(n1)]=spm_P(1,k(n1)*mat{5},u,mat{3},'F',mat{2},1,mat{4});
                                if isnan(cP(n1))&&~cp(n1),cP(n1)=0; end
                            end
                            ok=true;
                        end
                    end
                    if ~ok
                        try
                            if isempty(mat{5}), [nill,cP(n1),nill,cp(n1)] =stat_thres(mat{2},mat{4},mat{7},[mat{3};inf,inf], .05, u, k(n1)*mat{2}(end)/mat{4});
                            else, [nill,cP(n1),nill,cp(n1)] =stat_thres(mat{2},mat{4},mat{7},[mat{3};inf,inf], .05, u, k(n1)*mat{5});
                                if isnan(cP(n1))&&~cp(n1),cP(n1)=0; end
                            end
                        catch
                            cP(n1)=nan; cp(n1)=nan;
                        end
                    end
                    try
                        if mat{6}=='T'&&side==3, %two-sided T stat
                            pPFWE(n1)=spm_P(1,0,maxZ(n1).^2,mat{3},'F',mat{2},1,mat{4});
                        else
                            pPFWE(n1)=spm_P(1,0,maxZ(n1),mat{3},'F',mat{2},1,mat{4});
                        end
                    catch
                        pPFWE(n1)=nan;
                    end
                else
                    cP(n1)=nan; cp(n1)=nan; pPFWE(n1)=nan;
                end
            else 
                cP(n1)=nan;cp(n1)=nan; pPFWE(n1)=nan;
            end
            cpmass(n1)=nan; cPmass(n1)=nan; 
            cpscore(n1)=nan; cPscore(n1)=nan; 
        elseif thres{2}~=3&&thres{2}~=4&&parametric==2 % non-parametric stats
            nclL=k(n1);
            mclL=kmass(n1);
            sclL=kscore(n1);
            if ~isempty(iPERM)
                PERMp_cluster_size_unc=0;
                PERMp_cluster_size_FWE=0;
                PERMp_cluster_mass_unc=0;
                PERMp_cluster_mass_FWE=0;
                PERMp_cluster_score_unc=0;
                PERMp_cluster_score_FWE=0;
                NiPERM=0;
                for niPERM=1:numel(iPERM),
                    iperm=iPERM(niPERM);
                    kiPERM=size(SIM.Dist_Cluster_sizemax{iperm},1);
                    if nnz(SIM.Hist_Cluster_size{iperm})<2, PERMp_cluster_size_unc=PERMp_cluster_size_unc+kiPERM*double(1+nclL<=find(SIM.Hist_Cluster_size{iperm}));
                    else PERMp_cluster_size_unc=PERMp_cluster_size_unc+kiPERM*max(0,min(1,interp1(find(SIM.Hist_Cluster_size{iperm}),flipud(cumsum(flipud(nonzeros(SIM.Hist_Cluster_size{iperm})))),1+nclL,'linear','extrap')));
                    end
                    %PERMp_cluster_size_FDR=conn_fdr(PERMp_cluster_size_unc);
                    PERMp_cluster_size_FWE=PERMp_cluster_size_FWE+kiPERM*mean(conn_bsxfun(@ge,SIM.Dist_Cluster_sizemax{iperm}',nclL),2);
                    if nnz(SIM.Hist_Cluster_mass{iperm})<2, PERMp_cluster_mass_unc=PERMp_cluster_mass_unc+kiPERM*double(1+round(SIM.maxT*mclL)<=find(SIM.Hist_Cluster_mass{iperm}));
                    else PERMp_cluster_mass_unc=PERMp_cluster_mass_unc+kiPERM*max(0,min(1,interp1(find(SIM.Hist_Cluster_mass{iperm}),flipud(cumsum(flipud(nonzeros(SIM.Hist_Cluster_mass{iperm})))),1+round(SIM.maxT*mclL),'linear','extrap')));
                    end
                    %PERMp_cluster_mass_FDR=conn_fdr(PERMp_cluster_mass_unc);
                    PERMp_cluster_mass_FWE=PERMp_cluster_mass_FWE+kiPERM*mean(conn_bsxfun(@ge,SIM.Dist_Cluster_massmax{iperm}',mclL),2);
                    if nnz(SIM.Hist_Cluster_score{iperm})<2, PERMp_cluster_score_unc=PERMp_cluster_score_unc+kiPERM*double(1+round(SIM.maxT*sclL)<=find(SIM.Hist_Cluster_score{iperm}));
                    else PERMp_cluster_score_unc=PERMp_cluster_score_unc+kiPERM*max(0,min(1,interp1(find(SIM.Hist_Cluster_score{iperm}),flipud(cumsum(flipud(nonzeros(SIM.Hist_Cluster_score{iperm})))),1+round(SIM.maxT*sclL),'linear','extrap')));
                    end
                    %PERMp_cluster_score_FDR=conn_fdr(PERMp_cluster_score_unc);
                    PERMp_cluster_score_FWE=PERMp_cluster_score_FWE+kiPERM*mean(conn_bsxfun(@ge,SIM.Dist_Cluster_scoremax{iperm}',sclL),2);
                    NiPERM=NiPERM+kiPERM;
                end
                cP(n1)=PERMp_cluster_size_FWE/NiPERM;
                cp(n1)=PERMp_cluster_size_unc/NiPERM;
                cPmass(n1)=PERMp_cluster_mass_FWE/NiPERM;
                cpmass(n1)=PERMp_cluster_mass_unc/NiPERM;
                cPscore(n1)=PERMp_cluster_score_FWE/NiPERM;
                cpscore(n1)=PERMp_cluster_score_unc/NiPERM;
            else
                doneall=false;
                cP(n1)=nan;cp(n1)=nan; cPmass(n1)=nan;cpmass(n1)=nan; cPscore(n1)=nan;cpscore(n1)=nan;
            end
            pPFWE(n1)=nan;
        elseif (thres{2}==3||thres{2}==4)&&parametric==2 % TFCE
            peaks(n1)=nnz(tfce_peak_mask(clusters{n1}));
            if ~isempty(iPERM)
                pPFWE(n1)=mean(tfce_max_dist>maxZ(n1));
                idx=1:numel(tfce_peak_p_FDR);
                idx(clustervol(tfce_peak_mask)~=n1)=[];
                if ~isempty(idx), 
                    pPFDR(n1)=min(tfce_peak_p_FDR(idx)); 
                    pPunc(n1)=min(tfce_peak_p_unc(idx)); 
                end
            end
        end
    end
    if peakfdr>0,
        try
            if peakfdr==1, % consider all peak voxels for peak-FDR correction
                if side==2, tZ=-Z(idxvoxels); else tZ=Z(idxvoxels); end
                mtZ=min(tZ);
                [maxN,maxZ,maxXYZ,maxA]=spm_max(1-mtZ+tZ,xyz'); 
                maxZ=maxZ+mtZ-1;
                [maxOK,maxI]=min(sum(abs(conn_bsxfun(@minus,permute(maxXYZ',[1,3,2]),permute(cell2mat(xyzpeak'),[3,1,2]))).^2,3),[],1);
            else % consider only onepeak voxel within each cluster for peak-FDR correction
                maxI=1:numel(maxZ);
            end
            Ez=nan(size(maxZ));
            if mat{6}=='T'&&side<3,
                [nill,nill,Eu] = spm_P_RF(1,0,u,mat{3},'T',mat{2},1);
                for n1=1:numel(maxZ)
                    [nill,nill,Ez(n1)] = spm_P_RF(1,0,maxZ(n1),mat{3},'T',mat{2},1);
                end
            elseif mat{6}=='T'&&side==3,
                [nill,nill,Eu] = spm_P_RF(1,0,u,mat{3},'F',mat{2},1);
                for n1=1:numel(maxZ)
                    [nill,nill,Ez(n1)] = spm_P_RF(1,0,maxZ(n1).^2,mat{3},'F',mat{2},1);
                end
            elseif mat{6}=='F',
                [nill,nill,Eu] = spm_P_RF(1,0,u,mat{3},'F',mat{2},1);
                for n1=1:numel(maxZ)
                    [nill,nill,Ez(n1)] = spm_P_RF(1,0,maxZ(n1),mat{3},'F',mat{2},1);
                end
            else Eu=Ez;
            end
            pPFDR=conn_fdr(Ez/Eu);
            pPFDR=pPFDR(maxI);
            maxZ=maxZ(maxI);
        catch
            pPFDR=nan(size(k));
        end
    end
    cPFDR=conn_fdr(cp); 
    cPFDRmass=conn_fdr(cpmass);
    cPFDRscore=conn_fdr(cpscore);
    for n1=1:length(idx1),
        %temp=['( ',sprintf('%+03.0f ',(mat{1}(1:3,:)*[xyzpeak{n1}';1])'),') '];
        temp=[sprintf('%+03.0f ',(mat{1}(1:3,:)*[xyzpeak{n1}';1])')];
        cluster_stats{end+1}=sprintf('%s  k = %d  p = %.6f',temp,k(n1),cp(n1));
        temp=[temp,repmat(' ',[1,max(0,15-length(temp))])];
        if peakfdr>0,
            txt=strvcat(txt,[...
                temp,...
                [sprintf('%13d',k(n1))],...
                [sprintf('%13f',cP(n1))],...
                [sprintf('%13f',cPFDR(n1))],...
                [sprintf('%13f',cp(n1))],...
                [sprintf('%13f',pPFWE(n1))],...
                [sprintf('%13f',pPFDR(n1))],...
                [sprintf('%13f',minp(n1))]]);
            statsOut=struct('k',k,'cP',cP,'cPFDR',cPFDR,'cp',cp,'pPFWE',pPFWE,'pPFDR',pPFDR,'minp',minp);
            %[sprintf('%15f',minP(n1))]]);
        elseif thres{2}==3 % TFCE fdr
            txt=strvcat(txt,[...
                temp,...
                [sprintf('%13d',k(n1))],...
                [sprintf('%13d',peaks(n1))],...
                [sprintf('%13.2f',maxZ(n1))],...
                [sprintf('%13f',pPFWE(n1))],...
                [sprintf('%13f',pPFDR(n1))],...
                [sprintf('%13f',pPunc(n1))]]);
            statsOut=struct('k',k,'peaks',peaks,'maxZ',maxZ,'pPFWE',pPFWE,'pPFDR',pPFDR,'pPunc',pPunc);
        elseif thres{2}==4 % TFCE fwe
            txt=strvcat(txt,[...
                temp,...
                [sprintf('%13d',k(n1))],...
                [sprintf('%13d',peaks(n1))],...
                [sprintf('%13.2f',maxZ(n1))],...
                [sprintf('%13f',pPFWE(n1))],...
                [sprintf('%13f',pPFDR(n1))],...
                [sprintf('%13f',pPunc(n1))]]);
            statsOut=struct('k',k,'peaks',peaks,'maxZ',maxZ,'pPFWE',pPFWE,'pPFDR',pPFDR,'pPunc',pPunc);
        elseif parametric==1 % param
            txt=strvcat(txt,[...
                temp,...
                [sprintf('%13d',k(n1))],...
                [sprintf('%13f',cP(n1))],...
                [sprintf('%13f',cPFDR(n1))],...
                [sprintf('%13f',cp(n1))],...
                [sprintf('%13f',pPFWE(n1))],...
                [sprintf('%13f',minp(n1))]]);
            statsOut=struct('k',k,'cPFWE',pPFWE,'cPFDR',cPFDR,'pPFWE',pPFWE,'minp',minp);
        else % nonparam
            txt=strvcat(txt,[...
                temp,...
                [sprintf('%13d',k(n1))],...
                [sprintf('%13f',cP(n1))],...
                [sprintf('%13f',cPFDR(n1))],...
                [sprintf('%13f',cp(n1))],...
                [sprintf('%13.2f',kmass(n1))],...
                [sprintf('%13f',cPmass(n1))],...
                [sprintf('%13f',cPFDRmass(n1))],...
                [sprintf('%13f',cpmass(n1))]]);
            statsOut=struct('k',k,'cP',cP,'cPFDR',cPFDR,'cp',cp,'kmass',kmass,'cPmass',cPmass,'cPFDRmass',cPFDRmass,'cpmass',cpmass);
        end
    end
    if thres{2}==3
        idxclu=find(peaks>0);
        idxrem=find(~(peaks>0));
    else
        switch(thres4),
            case 0,%none
                idxclu=1:numel(k);
                idxrem=[];
            case 1,%'k-value',
                idxclu=find(k>=thres{3});
                idxrem=find(~(k>=thres{3}));
            case 2,%'clu-FWE',
                idxclu=find(cP<thres{3});
                idxrem=find(~(cP<thres{3}));
            case 3,%'clu-FDR',
                idxclu=find(cPFDR<thres{3});
                idxrem=find(~(cPFDR<thres{3}));
            case 4,%'clu-unc',
                idxclu=find(cp<thres{3});
                idxrem=find(~(cp<thres{3}));
            case 5,%'peak-FWE',
                idxclu=find(pPFWE<thres{3});
                idxrem=find(~(pPFWE<thres{3}));
            case 6,%'peak-FDR',
                idxclu=find(pPFDR<thres{3});
                idxrem=find(~(pPFDR<thres{3}));
            case 7,%'peak-unc',
                idxclu=find(minp<thres{3});
                idxrem=find(~(minp<thres{3}));
            case 8,%'mass-value',
                idxclu=find(kmass>=thres{3});
                idxrem=find(~(kmass>=thres{3}));
            case 9,%'mass clu-FWE',
                idxclu=find(cPmass<thres{3});
                idxrem=find(~(cPmass<thres{3}));
            case 10,%'mass clu-FDR',
                idxclu=find(cPFDRmass<thres{3});
                idxrem=find(~(cPFDRmass<thres{3}));
            case 11,%'mass clu-unc',
                idxclu=find(cpmass<thres{3});
                idxrem=find(~(cpmass<thres{3}));
        end
    end
    if ~isempty(idxclu), ithres(2)=min(k(idxclu)); end
    for n1=1:length(idxrem), anew(clusters{idxrem(n1)})=0; end
    txt=txt(idxclu,:);
    if ~isempty(txt), txt=char(regexprep(cellstr(txt),'NaN','---')); end
    xyzpeak={xyzpeak{idxclu}};
    clusters={clusters{idxclu}};
    [nill,clustervol]=ismember(clustervol,idxclu);
    cluster_stats={cluster_stats{idxclu}};
    if ~isempty(idxclu),
        for n1=1:length(idxclu)+1,
            if n1<=length(idxclu),xyztemp=xyzrois(idx2{sidxn(idxclu(n1))},:);
            else, xyztemp=xyzrois(cat(2,idx2{sidxn(idxclu)}),:); end;
            xyztemp=mat{1}*[xyztemp,ones(size(xyztemp,1),1)]';
            v=spm_get_data(refsrois.V,pinv(refsrois.V(1).mat)*xyztemp);
            if numel(refsrois.V)>1, [nill,v]=max(v,[],1); v(~nill)=0; end
            [uv,nill,iv]=unique(v);
            clusternames{n1}.uvb=zeros(1,length(uv));for n2=1:length(uv),clusternames{n1}.uvb(n2)=sum(v==uv(n2));end
            clusternames{n1}.anat=iv;
            clusternames{n1}.xyz=zeros(3,length(uv)); 
            for n2=1:length(uv), 
                txyz=xyztemp(1:3,v==uv(n2));
                mxyz=mean(txyz,2);
                [nill,i]=min(sum(abs(txyz-repmat(mxyz,1,size(txyz,2))).^2,1));
                clusternames{n1}.xyz(:,n2)=round(txyz(1:3,i)); % return point within ROI closest to ROI-centroid
            end
            clusternames{n1}.uvc=refsrois.count(1+uv);
            clusternames{n1}.uvd={refsrois.labels{refsrois.labelsidx(1+uv)}};
            [nill,uvidx]=sort(-clusternames{n1}.uvb+1e10*(uv==0)); rankvidx=uvidx; rankvidx(uvidx)=1:numel(uvidx);
            clusternames{n1}.uvb=clusternames{n1}.uvb(uvidx);clusternames{n1}.uvc=clusternames{n1}.uvc(uvidx);clusternames{n1}.uvd={clusternames{n1}.uvd{uvidx}};clusternames{n1}.xyz=clusternames{n1}.xyz(:,uvidx);clusternames{n1}.anat=rankvidx(clusternames{n1}.anat);
            clusternames{n1}.txt=cell(1,length(uv)); for n2=1:length(uv),clusternames{n1}.txt{n2}=sprintf('%d voxels (%d%%) covering %0.0f%% of %s with center at (%+d,%+d,%+d)',clusternames{n1}.uvb(n2),round(100*clusternames{n1}.uvb(n2)/size(xyztemp,2)),clusternames{n1}.uvb(n2)/clusternames{n1}.uvc(n2)*100,clusternames{n1}.uvd{n2}, clusternames{n1}.xyz(1,n2),clusternames{n1}.xyz(2,n2),clusternames{n1}.xyz(3,n2)); end
            %clusternames{n1}.txt=cell(1,length(uv)); for n2=1:length(uv),clusternames{n1}.txt{n2}=[num2str(clusternames{n1}.uvb(n2)),' voxels (',num2str(round(100*clusternames{n1}.uvb(n2)/size(xyztemp,1))),'%) covering ',num2str(clusternames{n1}.uvb(n2)/clusternames{n1}.uvc(n2)*100,'%0.0f'),'% of ',refsrois.filenameshort,'.',clusternames{n1}.uvd{n2}]; end
            %clusternames{n1}.txt=strvcat(clusternames{n1}.txt{:});
        end
    end
%     clusternames={clusternames{idxclu,1},clusternames{idxclu(idxclulast),2}};
    %for n1=1:length(clusternames),disp(' ');disp(clusternames{n1}.txt);end
    if ~doneall&&(thres4~=1&&thres4~=8), txt=strvcat('Cluster threshold information unavailable yet. Please re-compute non-parametric statatistics',txt); end
end

%if params.plotlist, figure('color','w'); end
%h=[];xlim=[];ylim=[];
%disp(txt);

%     if params.permutation && ~isempty(params.L),
%         txt{end}=[txt{end},'  p=',num2str(mean(params.L(:,1)>=k),'%0.4f'),' (',num2str(length(unique(params.L(params.L(:,1)>=k,2)))/params.Ln,'%0.4f'),')'];
%     end
%     if params.plotlist,
%         %subplot(length(idx1),1,n1);
%         h(n1)=axes('units','norm','position',...
%             [.9,.5*(1-n1/max(length(idx1)+1,5)+.05/max(length(idx1)+1,5)),.09,.5*.9/max(length(idx1)+1,5)]);
%         mx1=mean(X1(:,idx(idx2{sidxn(n1)})),2);
%         mx2=mean(X2(:,idx(idx2{sidxn(n1)})),2);
%         plot(mx1,mx2,'.'); h0=ylabel(txt{end});set(h0,'fontsize',8+CONN_gui.font_offset,'rotation',0,'horizontalalignment','right')
%         params.output.clusters{n1}=[mx1,mx2];
%         xlim=cat(1,xlim,[min(mx1),max(mx1)]);ylim=cat(1,ylim,[min(mx2),max(mx2)]);
%     end
%   end
%if params.plotlist,
%    for n1=1:length(idx1),
%        set(h(n1),'xlim',[min(xlim(:,1))-.01,max(xlim(:,2))+.01],'ylim',[min(ylim(:,1))-.01,max(ylim(:,2))+.01]);
%        if n1~=length(idx1), set(h(n1),'xticklabel',[],'yticklabel',[],'fontsize',6); end
%    end
%end
end

function [peak_threshold, extent_threshold, peak_threshold_1, extent_threshold_1] = ...
   stat_thres(search_volume, num_voxels, fwhm, df, p_val_peak, ...
   cluster_threshold, p_val_extent, nconj, nvar, EC_file, expr)
%
% stat_thresh.m
% Modified version of STAT_THRESHOLD function by Keith Worsley. The original
% STAT_THRESHOLD function is part of the FMRISTAT packgage. The main use of the
% original version is to calculate peak and cluster thresholds, both corrected
% and uncorrected. Details on input and output arguments are found in the
% original STAT_THRESHOLD function available from Keith Worsley's web site
% at McGill University's Math & Stat department. 
%
% This stat_thresh.m function is a customized version to be called by a function
% spm_list_nS for non-stationarity correction of RFT-based cluster size test.
% The input and output of this function is therefore modified for producing cluster
% p-values (corrected) under non-stationarity. The modification includes:
%   -supressing the output from being displayed.
%   -the number of cluster p-values it can calculate has been increased to 500 clusters
%    (the default in the original version was 5).
%   -the p_val_extent is treated as extent, no matter how small it is.
%
% stat_thresh is called by spm_list_nS in the following format:
% [PEAK_P CLUSTER_P PEAK_P_1 CLUSTER_P_1] = 
%    stat_thresh(V_RESEL,NUM_VOX,1,[DF_ER;DF_RPV],ALPHA,CL_DEF_TH,CL_RESEL);
% PARAMETERS:
%    V_RESEL:      The seach volume in terms of resels. It is a 1x4 vector
%                  describing the topological characteristic of the search
%                  volume.
%    NUM_VOX:      The number of voxels in the search volume.
%    DF_ER:        Degrees of freedom of error
%    DF_RPV:       Degrees of freedom of RPV image estimation. Usually the same
%                  as the error df.
%    ALPHA:        The significance level of the peak (arbitrarily set to 0.05) 
%    CL_DEF_TH:    The cluster defining threshold. Can be entered in terms of
%                  a p-value (uncorrected) or a t-value.
%    CL_RESEL:     The cluster size in terms of resel
%
%    PEAK_P:       Peak p-value (FWE corrected). Not used for our purpose.
%    PEAK_P_1:     Peak p-value (uncorrected). Not used for our purpose.
%    CLUSTER_P:    Cluster p-value (FWE corrected)
%    CLUSTER_P_1:  Cluster p-value (uncorrected)
%    
%                           ----------------
%
% More etails on non-stationary cluster size test can be found in
%
% Worsley K J, Andermann M, Koulis T, MacDonald D and Evans A C
%   Detecting Changes in Nonisotropic Images
%   Human Brain Mapping 8: 98-101 (1999)
%
% Hayasaka S, Phan K L, Liberzon I, Worsley K J, and Nichols T E
%   Nonstationary cluster size inference with random-field and permutation methods
%   NeuroImage 22: 676-687 (2004)
%
%
%-----------------------------------------------------------------------------------
% Version 0.76b   Feb 19, 2007  by Satoru Hayasaka
%

% ############################################################################
% COPYRIGHT:   Copyright 2003 K.J. Worsley 
%              Department of Mathematics and Statistics,
%              McConnell Brain Imaging Center, 
%              Montreal Neurological Institute,
%              McGill University, Montreal, Quebec, Canada. 
%              keith.worsley@mcgill.ca , www.math.mcgill.ca/keith
%
%              Permission to use, copy, modify, and distribute this
%              software and its documentation for any purpose and without
%              fee is hereby granted, provided that this copyright
%              notice appears in all copies. The author and McGill University
%              make no representations about the suitability of this
%              software for any purpose.  It is provided "as is" without
%              express or implied warranty.
% ############################################################################

% ############################################################################
% UPDATES:
%
%          Variable nvar is rounded so that it is recognized as an integer.
%          Feb 19, 2007 by Satoru Hayasaka
% ############################################################################


% Defaults:
if nargin<1;  search_volume=[];  end
if nargin<2;  num_voxels=[];  end
if nargin<3;  fwhm=[];  end
if nargin<4;  df=[];  end
if nargin<5;  p_val_peak=[];  end
if nargin<6;  cluster_threshold=[];  end
if nargin<7;  p_val_extent=[];  end
if nargin<8;  nconj=[];  end
if nargin<9;  nvar=[];  end
if nargin<10;  EC_file=[];  end
if nargin<11;  expr=[];  end

if isempty(search_volume);  search_volume=1000000;  end
if isempty(num_voxels);  num_voxels=1000000;  end
if isempty(fwhm);  fwhm=0.0;  end
if isempty(df);  df=Inf;  end
if isempty(p_val_peak);  p_val_peak=0.05;  end
if isempty(cluster_threshold);  cluster_threshold=0.001;  end
if isempty(p_val_extent);  p_val_extent=0.05;  end
if isempty(nconj);  nconj=1;  end
if isempty(nvar);  nvar=1;  end

if size(fwhm,1)==1; fwhm(2,:)=fwhm; end
if size(fwhm,2)==1; scale=1; else scale=fwhm(1,2)/fwhm(1,1); fwhm=fwhm(:,1); end;
isscale=(scale>1); 

if length(num_voxels)==1; num_voxels(2,1)=1; end

if size(search_volume,2)==1
   radius=(search_volume/(4/3*pi)).^(1/3);
   search_volume=[ones(length(radius),1) 4*radius 2*pi*radius.^2 search_volume];
end
if size(search_volume,1)==1
   search_volume=[search_volume; [1 zeros(1,size(search_volume,2)-1)]];
end
lsv=size(search_volume,2);
fwhm_inv=all(fwhm>0)./(fwhm+any(fwhm<=0));
resels=search_volume.*repmat(fwhm_inv,1,lsv).^repmat(0:lsv-1,2,1);
invol=resels.*(4*log(2)).^(repmat(0:lsv-1,2,1)/2);
for k=1:2
   D(k,1)=max(find(invol(k,:)))-1;
end

% determines which method was used to estimate fwhm (see fmrilm or multistat): 
df_limit=4;

% max number of pvalues or thresholds to print:
% it can print out a ton of stuff! (the original default was 5)
nprint=500;

if length(df)==1; df=[df 0]; end
if size(df,1)==1; df=[df; Inf Inf]; end
if size(df,2)==1; df=[df [0; df(2,1)]]; end

% is_tstat=1 if it is a t statistic
is_tstat=(df(1,2)==0);
if is_tstat
   df1=1;
   df2=df(1,1);
else
   df1=df(1,1);
   df2=df(1,2);
end
if df2 >= 1000; df2=Inf; end
df0=df1+df2;

dfw1=df(2,1);
dfw2=df(2,2);
if dfw1 >= 1000; dfw1=Inf; end
if dfw2 >= 1000; dfw2=Inf; end

if length(nvar)==1; nvar(2,1)=df1; end
nvar = round(nvar); %-to make sure that nvar is integer!

if isscale & (D(2)>1 | nvar(1,1)>1 | df2<Inf)
   D
   nvar
   df2
   fprintf('Cannot do scale space.');
   return
end

Dlim=D+[scale>1; 0];
DD=Dlim+nvar-1;

% Values of the F statistic:
t=((1000:-1:1)'/100).^4;
% Find the upper tail probs cumulating the F density using Simpson's rule:
if df2==Inf
   u=df1*t;
   b=exp(-u/2-log(2*pi)/2+log(u)/4)*df1^(1/4)*4/100;
else  
   u=df1*t/df2;
   b=exp(-df0/2*log(1+u)+log(u)/4-betaln(1/2,(df0-1)/2))*(df1/df2)^(1/4)*4/100;
end
t=[t; 0];
b=[b; 0];
n=length(t);
sb=cumsum(b);
sb1=cumsum(b.*(-1).^(1:n)');
pt1=sb+sb1/3-b/3;
pt2=sb-sb1/3-b/3;
tau=zeros(n,DD(1)+1,DD(2)+1);
tau(1:2:n,1,1)=pt1(1:2:n);
tau(2:2:n,1,1)=pt2(2:2:n);
tau(n,1,1)=1;
tau(:,1,1)=min(tau(:,1,1),1);

% Find the EC densities:
u=df1*t;
 kk=(max(DD)-1+min(DD))/2;
 uu=conn_bsxfun(@power,u,0:.5:kk);
 [ii,jj]=ndgrid(0:kk,0:kk);
 gammalnii=gammalni(0:max(DD)+kk);
for d=1:max(DD)
   for e=0:min(min(DD),d)
      s1=0;
      cons=-((d+e)/2+1)*log(pi)+gammaln(d)+gammaln(e+1);
      for k=0:(d-1+e)/2
         %[i,j]=ndgrid(0:k,0:k);
         i=ii(1:k+1,1:k+1); j=jj(1:k+1,1:k+1);
         if df2==Inf
            q1=log(pi)/2-((d+e-1)/2+i+j)*log(2);
         else
            q1=(df0-1-d-e)*log(2)+gammaln((df0-d)/2+i)+gammaln((df0-e)/2+j) ...
               -gammalni(df0-d-e+i+j+k)-((d+e-1)/2-k)*log(df2);
         end
         %q2=cons-gammalni(i+1)-gammalni(j+1)-gammalni(k-i-j+1) ...
         %   -gammalni(d-k-i+j)-gammalni(e-k-j+i+1);
         q2=cons-gammalnii(1+max(0,i+1))-gammalnii(1+max(0,j+1))-gammalnii(1+max(0,k-i-j+1)) ...
            -gammalnii(1+max(0,d-k-i+j))-gammalnii(1+max(0,e-k-j+i+1));
         s2=sum(sum(exp(q1+q2)));
         if s2>0
            %s1=s1+(-1)^k*u.^((d+e-1)/2-k)*s2;
            s1=s1+(-1)^k*uu(:,1+2*((d+e-1)/2-k))*s2;
         end
      end
      if df2==Inf
         s1=s1.*exp(-u/2);
      else
         s1=s1.*exp(-(df0-2)/2*log(1+u/df2));
      end
      if DD(1)>=DD(2)
         tau(:,d+1,e+1)=s1;
         if d<=min(DD)
            tau(:,e+1,d+1)=s1;
         end
      else
         tau(:,e+1,d+1)=s1;      
         if d<=min(DD)
            tau(:,d+1,e+1)=s1;
         end
      end
   end
end

% For multivariate statistics, add a sphere to the search region:
a=zeros(2,max(nvar));
for k=1:2
   j=(nvar(k)-1):-2:0;
   a(k,j+1)=exp(j*log(2)+j/2*log(pi) ...
      +gammaln((nvar(k)+1)/2)-gammaln((nvar(k)+1-j)/2)-gammaln(j+1));
end
rho=zeros(n,Dlim(1)+1,Dlim(2)+1);
for k=1:nvar(1)
   for l=1:nvar(2)
      rho=rho+a(1,k)*a(2,l)*tau(:,(0:Dlim(1))+k,(0:Dlim(2))+l);
   end
end

if is_tstat
   t=[sqrt(t(1:(n-1))); -flipdim(sqrt(t),1)];
   rho=[rho(1:(n-1),:,:); flipdim(rho,1)]/2;
   for i=0:D(1)
      for j=0:D(2)
         rho(n-1+(1:n),i+1,j+1)=-(-1)^(i+j)*rho(n-1+(1:n),i+1,j+1);
      end
   end
   rho(n-1+(1:n),1,1)=rho(n-1+(1:n),1,1)+1;
   n=2*n-1;
end

% For scale space:
if scale>1
   kappa=D(1)/2;
   tau=zeros(n,D(1)+1);
   for d=0:D(1)
      s1=0;
      for k=0:d/2
         s1=s1+(-1)^k/(1-2*k)*exp(gammaln(d+1)-gammaln(k+1)-gammaln(d-2*k+1) ...
            +(1/2-k)*log(kappa)-k*log(4*pi))*rho(:,d+2-2*k,1);
      end
      if d==0
         cons=log(scale);
      else
         cons=(1-1/scale^d)/d;
      end
      tau(:,d+1)=rho(:,d+1,1)*(1+1/scale^d)/2+s1*cons;
   end
   rho(:,1:(D(1)+1),1)=tau;
end

if D(2)==0
   d=D(1);
   if nconj>1
      % Conjunctions:
      b=gamma(((0:d)+1)/2)/gamma(1/2);
      for i=1:d+1
         rho(:,i,1)=rho(:,i,1)/b(i);
      end
      m1=zeros(n,d+1,d+1);
      for i=1:d+1
         j=i:d+1;
         m1(:,i,j)=rho(:,j-i+1,1);
      end
      for k=2:nconj
         for i=1:d+1
            for j=1:d+1
               m2(:,i,j)=sum(rho(:,1:d+2-i,1).*m1(:,i:d+1,j),2);
            end
         end
         m1=m2;
      end
      for i=1:d+1
         rho(:,i,1)=m1(:,1,i)*b(i);
      end
   end
   
   if ~isempty(EC_file)
      if d<3
         rho(:,(d+2):4,1)=zeros(n,4-d-2+1);
      end
      fid=fopen(EC_file,'w');
      % first 3 are dimension sizes as 4-byte integers:
      fwrite(fid,[n max(d+2,5) 1],'int');
      % next 6 are bounding box as 4-byte floats: 
      fwrite(fid,[0 0 0; 1 1 1],'float');
      % rest are the data as 4-byte floats:
      if ~isempty(expr)
         eval(expr);
      end
      fwrite(fid,t,'float');
      fwrite(fid,rho,'float');
      fclose(fid);
   end
end

if all(fwhm>0)
   pval_rf=zeros(n,1);
   for i=1:D(1)+1
      for j=1:D(2)+1
         pval_rf=pval_rf+invol(1,i)*invol(2,j)*rho(:,i,j);
      end
   end
else
   pval_rf=Inf;
end

% Bonferroni 
pt=rho(:,1,1);
pval_bon=abs(prod(num_voxels))*pt;

% Minimum of the two:
pval=min(pval_rf,pval_bon);

tlim=1;
if p_val_peak(1) <= tlim
   peak_threshold=minterp1(pval,t,p_val_peak);
   if length(p_val_peak)<=nprint
      peak_threshold;
   end
else
   % p_val_peak is treated as a peak value:
   P_val_peak=interp1(t,pval,p_val_peak);
   peak_threshold=P_val_peak;
   if length(p_val_peak)<=nprint
      P_val_peak;
   end
end

if fwhm<=0 | any(num_voxels<0)
   peak_threshold_1=p_val_peak+NaN;
   extent_threshold=p_val_extent+NaN;
   extent_threshold_1=extent_threshold;
   return
end

% Cluster_threshold:
% ###-Changed so that cluster_threshold is considered as cluster extent no matter what.
if cluster_threshold > eps
   tt=cluster_threshold;
else
   % cluster_threshold is treated as a probability:
   tt=minterp1(pt,t,cluster_threshold);
   Cluster_threshold=tt;
end

d=D(1);
rhoD=interp1(t,rho(:,d+1,1),tt);
p=interp1(t,pt,tt);

% Pre-selected peak:

pval=rho(:,d+1,1)./rhoD;
if p_val_peak(1) <= tlim 
   peak_threshold_1=minterp1(pval,t, p_val_peak);
   if length(p_val_peak)<=nprint
      peak_threshold_1;
   end
else
   % p_val_peak is treated as a peak value:
   P_val_peak_1=interp1(t,pval,p_val_peak);
   peak_threshold_1=P_val_peak_1;
   if length(p_val_peak)<=nprint
      P_val_peak_1;
   end
end

if  D(1)==0 | nconj>1 | nvar(1)>1 | D(2)>0 | scale>1
    extent_threshold=p_val_extent+NaN;
    extent_threshold_1=extent_threshold;
    if length(p_val_extent)<=nprint
       extent_threshold;
       extent_threshold_1;
    end
    return
end

% Expected number of clusters:

% ###-Change tlim to a small number so that p_val_extent is considered as cluster extent
tlim = eps;

EL=invol(1,d+1)*rhoD;
cons=gamma(d/2+1)*(4*log(2))^(d/2)/fwhm(1)^d*rhoD/p;

if df2==Inf & dfw1==Inf
   if p_val_extent(1) <= tlim 
      pS=-log(1-p_val_extent)/EL;
      extent_threshold=(-log(pS)).^(d/2)/cons;
      pS=-log(1-p_val_extent);
      extent_threshold_1=(-log(pS)).^(d/2)/cons;
      if length(p_val_extent)<=nprint
         extent_threshold;
         extent_threshold_1;
      end
   else
      % p_val_extent is now treated as a spatial extent:
      pS=exp(-(p_val_extent*cons).^(2/d));
      P_val_extent=1-exp(-pS*EL);
      extent_threshold=P_val_extent;
      P_val_extent_1=1-exp(-pS);
      extent_threshold_1=P_val_extent_1;
      if length(p_val_extent)<=nprint
         P_val_extent;
         P_val_extent_1;
      end
   end
else
   % Find dbn of S by taking logs then using fft for convolution:
   ny=2^12;
   a=d/2;
   b2=a*10*max(sqrt(2/(min(df1+df2,dfw1))),1);
   if df2<Inf
      b1=a*log((1-(1-0.000001)^(2/(df2-d)))*df2/2);
   else
      b1=a*log(-log(1-0.000001));
   end
   dy=(b2-b1)/ny;
   b1=round(b1/dy)*dy;
   y=((1:ny)'-1)*dy+b1;
   numrv=1+(d+1)*(df2<Inf)+d*(dfw1<Inf)+(dfw2<Inf);
   f=zeros(ny,numrv);
   mu=zeros(1,numrv);
   if df2<Inf
      % Density of log(Beta(1,(df2-d)/2)^(d/2)):
      yy=exp(y./a)/df2*2;  
      yy=yy.*(yy<1);
      f(:,1)=(1-yy).^((df2-d)/2-1).*((df2-d)/2).*yy/a;
      mu(1)=exp(gammaln(a+1)+gammaln((df2-d+2)/2)-gammaln((df2+2)/2)+a*log(df2/2));
   else
      % Density of log(exp(1)^(d/2)):
      yy=exp(y./a);   
      f(:,1)=exp(-yy).*yy/a;
      mu(1)=exp(gammaln(a+1));
   end
   
   nuv=[];
   aav=[];
   if df2<Inf
      nuv=[df1+df2-d  df2+2-(1:d)];
      aav=[a repmat(-1/2,1,d)]; 
   end
   if dfw1<Inf
      if dfw1>df_limit
         nuv=[nuv dfw1-dfw1/dfw2-(0:(d-1))];
      else
         nuv=[nuv repmat(dfw1-dfw1/dfw2,1,d)];
      end
      aav=[aav repmat(1/2,1,d)];
   end   
   if dfw2<Inf
      nuv=[nuv dfw2];
      aav=[aav -a];
   end   
   
   for i=1:(numrv-1)
      nu=nuv(i);
      aa=aav(i);
      yy=y/aa+log(nu);
      % Density of log((chi^2_nu/nu)^aa):
      f(:,i+1)=exp(nu/2*yy-exp(yy)/2-(nu/2)*log(2)-gammaln(nu/2))/abs(aa);
      mu(i+1)=exp(gammaln(nu/2+aa)-gammaln(nu/2)-aa*log(nu/2));
   end
   % Check: plot(y,f); sum(f*dy,1) should be 1
      
   omega=2*pi*((1:ny)'-1)/ny/dy;
   shift=complex(cos(-b1*omega),sin(-b1*omega))*dy;
   prodfft=prod(fft(f),2).*shift.^(numrv-1);
   % Density of Y=log(B^(d/2)*U^(d/2)/sqrt(det(Q))):
   ff=real(ifft(prodfft));
   % Check: plot(y,ff); sum(ff*dy) should be 1
   mu0=prod(mu);
   % Check: plot(y,ff.*exp(y)); sum(ff.*exp(y)*dy.*(y<10)) should equal mu0   
   
   alpha=p/rhoD/mu0*fwhm(1)^d/(4*log(2))^(d/2);
   
   % Integrate the density to get the p-value for one cluster: 
   pS=cumsum(ff(ny:-1:1))*dy;
   pS=pS(ny:-1:1);
   % The number of clusters is Poisson with mean EL:
   pSmax=1-exp(-pS*EL);
   
   if p_val_extent(1) <= tlim 
      yval=minterp1(-pSmax,y,-p_val_extent);
      % Spatial extent is alpha*exp(Y) -dy/2 correction for mid-point rule:
      extent_threshold=alpha*exp(yval-dy/2);
      % For a single cluster:
      yval=minterp1(-pS,y,-p_val_extent);
      extent_threshold_1=alpha*exp(yval-dy/2);
      if length(p_val_extent)<=nprint
         extent_threshold;
         extent_threshold_1;
      end
   else
      % p_val_extent is now treated as a spatial extent:
      P_val_extent=interp1(y,pSmax,log(p_val_extent/alpha)+dy/2);
      extent_threshold=P_val_extent;
      % For a single cluster:
      P_val_extent_1=interp1(y,pS,log(p_val_extent/alpha)+dy/2);
      extent_threshold_1=P_val_extent_1;
      if length(p_val_extent)<=nprint
         P_val_extent;
         P_val_extent_1;
      end
   end
   
end
   
return
end

function x=gammalni(n)
i=find(n>=0);
x=Inf+n;
if ~isempty(i)
   x(i)=gammaln(n(i));
end
return
end

function iy=minterp1(x,y,ix)
% interpolates only the monotonically increasing values of x at ix
n=length(x);
mx=x(1);
my=y(1);
xx=x(1);
for i=2:n
   if x(i)>xx
      xx=x(i);
      mx=[mx xx];
      my=[my y(i)];
   end
end
iy=interp1(mx,my,ix);
return
end

function uiwrap(h)
global CONN_gui;
if ~isfield(CONN_gui,'uicontrol_border'), CONN_gui.uicontrol_border=2; end
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

function uimotionfcn(option,varargin)
h1=findobj(gcbf,'tag','highlight');
h2=findobj(gcbf,'tag','highlight_pointer');
h3=findobj(gcbf,'tag','highlight_image');
if ~isempty(h1)||~isempty(h2)
    p1=get(0,'pointerlocation');
    p2=get(gcbf,'position');
    p3=get(0,'screensize');
    p2(1:2)=p2(1:2)+p3(1:2)-1; % note: fix issue when connecting to external monitor/projector
    pos=(p1-p2(1:2))./p2(3:4);
    c1=get(gcbf,'color');
    c0=mean(c1<.5);
    if ~isempty(h1)
        active=get(h1,'position');
        if iscell(active), active=cell2mat(active(:)); end
        xpos=repmat(pos,[size(active,1),1]); 
        %xpos1=active; xpos1=round(xpos1*500)/500;
        %if size(xpos,1)>1, xpos1=num2cell(xpos1,2); end
        md=find(min(min(xpos-active(:,1:2),active(:,1:2)+active(:,3:4)-xpos),[],2)>-.001)';
        set(h1,'foregroundcolor',[.4 .4 .4]+.2*c0);%,{'position'},xpos1);
        if ~isempty(md), 
            %xpos1=active(md,:); xpos1=round(xpos1*500)/500+repmat([-.0002,-.0002,.0004,.0004],size(xpos1,1),1);
            %if size(xpos,1)>1, xpos1=num2cell(xpos1,2); end
            set(h1(md),'foregroundcolor',[.0 .0 .0]+1*c0);%,{'position'},xpos1); 
        end
    end
    if ~isempty(h2)
        active=get(h2,'position');
        if iscell(active), active=cell2mat(active(:)); end
        xpos=repmat(pos,[size(active,1),1]);
        md=find(min(min(xpos-active(:,1:2),active(:,1:2)+active(:,3:4)-xpos),[],2)>-.001)';
        if ~isempty(md), set(gcf,'Pointer','crosshair'); else set(gcf,'Pointer','arrow'); end
    end
    if ~isempty(h3)
        active=get(h3,'position');
        if iscell(active), active=cell2mat(active(:)); end
        xpos=repmat(pos,[size(active,1),1]); 
        md=find(min(min(xpos-active(:,1:2),active(:,1:2)+active(:,3:4)-xpos),[],2)>-.001)';
        if ~isempty(md), 
            for n1=1:numel(h3)
                cd=get(h3(n1),'cdata');
                mcd=max(max(cd,[],1),[],2);
                mcd=mcd/mcd(1);
                if any(n1==md), mcd(2:3)=mcd(2:3)*1.10; end
                if any(mcd~=1), set(h3(n1),'cdata',cd.*repmat(1./mcd,[size(cd,1),size(cd,2),1])); end
            end
        end
    end
end
end

function conn_vproject_callbackfcn0(varargin)
global CONN_x;
if get(gcbo,'value')>1,
    tdata=char(get(gcbo,'userdata'));
    [tfilename,tpathname]=conn_fileutils('uigetfile','*.nii; *.img','Select mask/ROI file',tdata);
    if ischar(tfilename), tdata=fullfile(tpathname,tfilename); set(gcbo,'userdata',tdata); 
    else set(gcbo,'value',1); 
    end
end
end
function conn_vproject_callbackfcn1(varargin)
global CONN_x;
if get(gcbo,'value')==3,
    data=get(gcbo,'userdata');
    if isempty(data), data.selected=1; data.contrast=1; set(gcbo,'userdata',data); end
    data.selected=listdlg('liststring',CONN_x.Setup.l2covariates.names(1:end-1),'selectionmode','multiple','initialvalue',data.selected,'promptstring','Select subject-effects:','ListSize',[300 200]);
    if isempty(data.selected), set(gcbo,'value',1); return; end
    if numel(data.selected)==1, data.contrast=1;
    else
        if size(data.contrast,2)~=numel(data.selected), data.contrast=ones(1,numel(data.selected))/numel(data.selected); end
        answ=conn_menu_inputdlg(sprintf('Between-subjects contrast (vector with %d values)',numel(data.selected)),'',1,mat2str(data.contrast));
        if isempty(answ), return; end
        data.contrast=str2num(answ{1});
    end
    set(gcbo,'userdata',data);
end
end
function conn_vproject_callbackfcn2(validconditions,varargin)
global CONN_x;
if get(gcbo,'value')>numel(validconditions),
    data=get(gcbo,'userdata');
    if isempty(data), data.selected=1; data.contrast=1; set(gcbo,'userdata',data); end
    data.selected=listdlg('liststring',CONN_x.Setup.conditions.names(validconditions),'selectionmode','multiple','initialvalue',data.selected,'promptstring','Select conditions:','ListSize',[300 200]);
    if isempty(data.selected), set(gcbo,'value',1); return; end
    if numel(data.selected)==1, data.contrast=1;
    else
        if size(data.contrast,2)~=numel(data.selected), data.contrast=ones(1,numel(data.selected))/numel(data.selected); end
        answ=conn_menu_inputdlg(sprintf('Between-conditions contrast (vector with %d values)',numel(data.selected)),'',1,mat2str(data.contrast));
        if isempty(answ), return; end
        data.contrast=str2num(answ{1});
    end
    set(gcbo,'userdata',data);
end
end


function [filename_rois,filename_sources,viewrex]=conn_vproject_selectfiles(filename_rois,filename_sources,viewrex)
global CONN_gui;
if isempty(CONN_gui)||~isfield(CONN_gui,'font_offset'), conn_font_init; end

filename_rois0=filename_rois;
filename_sources0=filename_sources;
thfig=dialog('units','norm','position',[.3,.4,.4,.25],'windowstyle','normal','name','REX interface','color','w','resize','on');
uicontrol(thfig,'style','text','units','norm','position',[.1,.75,.8,.20],'string',{'Explore model effects and activation/connectivity values','within individual ROIs/clusters'},'backgroundcolor','w','fontsize',9+CONN_gui.font_offset,'fontweight','bold');
ht1=uicontrol(thfig,'style','popupmenu','units','norm','position',[.1,.55,.8,.15],'string',{'clusters of interest in current analysis','others clusters of interest or ROIs (select exported mask/ROI file)'},'fontsize',8+CONN_gui.font_offset,'horizontalalignment','left','callback',@conn_vproject_selectfiles_callback1,'tooltipstring','Define the clusters of interest');
ht2=uicontrol(thfig,'style','popupmenu','units','norm','position',[.1,.40,.8,.15],'string',{'effect/activation/connectivity values in current analysis','other effect/activation/connectivity values (select second-level SPM.mat file)'},'fontsize',8+CONN_gui.font_offset,'horizontalalignment','left','callback',@conn_vproject_selectfiles_callback2,'tooltipstring','Define the activation/connectivity values');
if ~isempty(viewrex), ht3=uicontrol(thfig,'style','checkbox','units','norm','position',[.1,.3,.8,.1],'string','enable REX gui','value',0,'fontsize',8+CONN_gui.font_offset,'horizontalalignment','left','backgroundcolor','w','tooltipstring','Displays REX gui interface for additional options'); end
uicontrol(thfig,'style','pushbutton','string','OK','units','norm','position',[.1,.01,.38,.2],'callback','uiresume','fontsize',8+CONN_gui.font_offset);
uicontrol(thfig,'style','pushbutton','string','Cancel','units','norm','position',[.51,.01,.38,.2],'callback','delete(gcbf)','fontsize',8+CONN_gui.font_offset);
uiwait(thfig);
ok=ishandle(thfig);
if ok, 
    if ~isempty(viewrex), viewrex=get(ht3,'value'); end
    delete(thfig);
else   filename_rois=[]; filename_sources=[]; 
end

    function conn_vproject_selectfiles_callback1(varargin)
        if get(ht1,'value')==1
            filename_rois=filename_rois0;
        else
            [tfilename,tpathname]=conn_fileutils('uigetfile','*.nii; *.img','Select mask/ROI file',filename_rois);
            if ischar(tfilename), filename_rois=fullfile(tpathname,tfilename);
            else
                filename_rois=filename_rois0;
                set(ht1,'value',1);
            end
        end
    end
    function conn_vproject_selectfiles_callback2(varargin)
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

function simfilename=conn_vproject_simfilename(spmfile,THR_TYPE,THR)
if nargin==2&&isequal(THR_TYPE,'all')
    simfilename=fullfile(fileparts(spmfile),'nonparametric_p*.mat');
else
    simfilename=char(arrayfun(@(a,b)fullfile(fileparts(spmfile),sprintf('nonparametric_p%d_%.8f.mat',a,b)),THR_TYPE,THR,'uni',0));
    if ~isempty(conn_fileutils('dir',conn_prepend('parallel_*_',simfilename))), conn_process('results_nonparametric_collapse',conn_server('util_localfile_filesep',[],simfilename)); end
end
end

function ok=conn_vproject_randomise(spmfile,THR_TYPE,THR,dogui,isremote)
ok=true;
[tstr,tidx]=conn_jobmanager('profiles');
tnull=find(strcmp('Null profile',tstr));
tlocal=find(strcmp('Background process (Unix,Mac)',tstr),1);
tvalid=setdiff(1:numel(tstr),tnull);
tstr=cellfun(@(x)sprintf('distributed processing (run on %s)',x),tstr,'uni',0);
if 1, tvalid=tidx; if isunix&&~isempty(tlocal)&&~ismember(tlocal,tvalid), tvalid=[tvalid(:)' tlocal]; end % show default+local profiles
elseif 1, tvalid=tidx; % show only default profile
else tstr{tidx}=sprintf('<HTML><b>%s</b></HTML>',tstr{tidx}); % show all profiles
end
if THR==0, v2=1000;
else v2=max(1000,round(1/THR));
end

fh=figure('units','norm','position',[.4,.4,.3,.2],'menubar','none','numbertitle','off','name','compute non-parametric statistics','color','w');
h1=uicontrol('units','norm','position',[.1,.7,.4,.15],'style','text','string','# of new simulations: ','fontweight','bold','backgroundcolor',get(fh,'color'));
h2=uicontrol('units','norm','position',[.5,.7,.4,.15],'style','edit','string',num2str(v2),'tooltipstring','<HTML>Number of new data permutations/randomizations that will be evaluated in order to compute cluster-level statistics<br/> - note: if already previously computed, these new simulations will be added to any pre-existing ones (e.g. to increase the total number of simulations)</HTML>');
toptions=[{'local processing (run on this computer)'} tstr(tvalid)];
if isremote
    info=conn_remotely('info');
    if isfield(info,'host')&&~isempty(info.host), tnameserver=info.host;
    elseif isfield(info,'remote_ip')&&~isempty(info.remote_ip), tnameserver=info.remote_ip;
    else tnameserver='CONN server';
    end
    toptions=regexprep(toptions(2:end),'\<run on (this computer)?',['run on ',tnameserver,' ']);
end
h3=uicontrol('style','popupmenu','units','norm','position',[.1,.2,.8,.15],'string',toptions,'value',1,'callback',@conn_vproject_randomise_nowlater); 
%h3=uicontrol('style','popupmenu','units','norm','position',[.1,.2,.8,.15],'string',[{'local processing (run on this computer)' 'queue/script it (save as scripts to be run later)'} tstr(tvalid)],'value',1,'callback',@conn_vproject_randomise_nowlater); 
%h3=uicontrol('units','norm','position',[.1,.5,.8,.25],'style','popupmenu','string',{'Run simulations now','Run simulations later'},'callback',@conn_vproject_randomise_nowlater);
h4=uicontrol('units','norm','position',[.1,.5,.4,.15],'style','text','string','# of parallel jobs: ','fontweight','bold','backgroundcolor',get(fh,'color'),'visible','off');
h5=uicontrol('units','norm','position',[.5,.5,.4,.15],'style','edit','string','1','visible','off','tooltipString','Number of parallel jobs to divide the permutation/randomization analyses');
h6=uicontrol('units','norm','position',[.1,.025,.4,.175],'style','pushbutton','string','Start','callback','uiresume(gcbf)');
h7=uicontrol('units','norm','position',[.5,.025,.4,.175],'style','pushbutton','string','Cancel','callback','close(gcbf)');
if nargin>=4&&isstruct(dogui), % run in parallel without prompting gui
    set(h2,'string',num2str(dogui.number_of_simulations));
    set(h5,'string',num2str(dogui.number_of_parallel_jobs));
    assert(dogui.number_of_parallel_jobs>0 | ~isremote, 'Not possible to run this locally when working with remote projects, please specify a value greater than zero in ''number_of_parallel_jobs''');
    if dogui.number_of_parallel_jobs>0, set(h3,'value',1+~isremote);
    else set(h3,'value',1);
    end
    dogui=false;
end
if nargin<4||dogui, uiwait(fh); end
if ishandle(fh),
    v2=str2num(get(h2,'string'));
    v3=get(h3,'value');
    v5=get(h5,'string');
    close(fh);
    if isremote, v3=v3+1; end
    niters=v2;
    simfilename=conn_vproject_simfilename(spmfile,THR_TYPE,THR);
    maskfile=fullfile(fileparts(spmfile),'mask.nii');
    if niters<=0, return; end
    if ~conn_existfile(maskfile), maskfile=fullfile(fileparts(spmfile),'mask.img'); end
    if conn_existfile(maskfile), mask=conn_fileutils('spm_read_vols',conn_fileutils('spm_vol',maskfile)); 
    else mask=[]; 
    end
    SIDE=1:3;
    THR=THR+[0 0 0];
    THR_TYPE=THR_TYPE+[0 0 0];
    SPM=struct; conn_loadmatfile(spmfile,'SPM');
    if isfield(SPM,'xX_multivariate')
        X=SPM.xX_multivariate.X;
        c=SPM.xX_multivariate.C;
        m=SPM.xX_multivariate.M;
        if isfield(SPM.xX_multivariate,'type'), opt=SPM.xX_multivariate.type; 
        else opt=[];
        end
    else
        conn_disp('Warning: xX_multivariate design info not found. Assuming no covariance modeling design. First contrast only');
        ncon=1;
        X=SPM.xX.X;
        c=SPM.xCon(ncon).c';
        m=1;
        opt=[];
    end
    if isfield(SPM,'repeatedsubjects'), groupingsamples=SPM.repeatedsubjects; else groupingsamples=[]; end
    if v3==1, % now
        ht=conn_msgbox('Preparing data. Please wait...','results explorer',-1);
        try
            a=conn_fileutils('spm_vol',char(SPM.xY.Y));
        catch
            a=SPM.xY.VY;
        end
        if ~conn_existfile(a(1).fname)&&isfield(SPM,'swd')&&isempty(fileparts(a(1).fname)),for n=1:numel(a), if isempty(fileparts(a(n).fname)), a(n).fname=fullfile(SPM.swd,a(n).fname); end; end; end
        if ~conn_existfile(a(1).fname), conn_msgbox({sprintf('Unable to find file %s',a(1).fname),'Please re-compute second-level model and try again'},'Outdated file references',2); return; end
        y=[];
        if isempty(mask)
            y=conn_fileutils('spm_read_vols',a);
            mask=~any(isnan(y),4)&any(diff(y,1,4)~=0,4);
            y=reshape(y,size(y,1)*size(y,2)*size(y,3),size(y,4));
            y=y(mask,:);
        end
        fmask=find(mask);
        [i,j,k]=ind2sub(size(mask),fmask);
        xyz=[i(:) j(:) k(:)]';
        if isempty(y), y=conn_fileutils('spm_get_data',a,xyz)'; end
        y=permute(reshape(y,size(y,1),size(X,1),[]),[2,3,1]);
        if conn_surf_dimscheck(a),
            surfparams=load(fullfile(fileparts(which(mfilename)),'utils','surf','surf_top.mat'),'A');
            adj=kron(speye(2),surfparams.A);
            adj=adj(fmask,fmask);
        else adj=xyz;
        end
        if ishandle(ht), delete(ht); end 
        try
            conn_randomise(X,y,c,m,opt,THR,THR_TYPE,SIDE,niters,simfilename,[],adj,[],groupingsamples);
            ok=true;
        catch
            ok=false;
        end
    elseif 1 % distributed
        parallel=tvalid(v3-1); 
        %if v3==2, parallel=find(strcmp('Null profile',conn_jobmanager('profiles'))); 
        %elseif v3>2, parallel=tvalid(v3-2); 
        %end
        if isfield(SPM.xY,'Y'), Y=char(SPM.xY.Y);
        else Y=char({SPM.xY.VY.fname});
        end
        if ~isremote&&~conn_existfile(deblank(Y(1,:)))&&isfield(SPM,'swd')&&isempty(fileparts(deblank(Y(1,:)))),Y=cellstr(Y); for n=1:numel(Y), if isempty(fileparts(deblank(Y{n}))), Y{n}=fullfile(SPM.swd,Y{n}); end; end; Y=char(Y); end
        if ~isremote&&~conn_existfile(deblank(Y(1,:))), conn_msgbox({sprintf('Unable to find file %s',deblank(Y(1,:))),'Please re-compute second-level model and try again'},'Outdated file references',2); return; end
        tfilesep=conn_projectmanager('filesep');
        Y=conn_server('util_localfile_filesep',tfilesep,Y);
        simfilename=conn_server('util_localfile_filesep',tfilesep,simfilename);
        N=str2num(v5);
        tfilename=fullfile(fileparts(spmfile),'args_nonparam.mat');
        niters=max(1,ceil(niters/N));
        conn_savematfile(tfilename,'X','Y','c','m','THR','THR_TYPE','SIDE','niters','simfilename','mask','groupingsamples','opt');
        conn_jobmanager('options','profile',parallel);
        info=conn_jobmanager('submit','orphan_results_nonparametric',N,N,[],conn_server('util_localfile_filesep',tfilesep,tfilename));
        info=conn_jobmanager(info,'','donotupdate'); 
        ok=true; 
        %if v3>2, info=conn_jobmanager(info,'','donotupdate'); ok=true; 
        %else ok=false; 
        %end
        %conn_process('results_nonparametric_collapse',simfilename);
    else % later (OBSOLETE)
        [file2_path,file2_name,file2_ext]=fileparts(v5);
        if isfield(SPM.xY,'Y'), Y=char(SPM.xY.Y);
        else Y=char({SPM.xY.VY.fname});
        end
        if conn_existfile(fullfile(file2_path,[file2_name,'.m']))
            answ=conn_questdlg(sprintf('Overwrite %s file?',fullfile(file2_path,[file2_name,'.m'])),'warning','Yes','No','Yes');
            if strcmp(answ,'No'), ok=false; return; end
        end
        conn_savematfile(fullfile(file2_path,[file2_name,'.mat']),'X','Y','c','m','THR','THR_TYPE','SIDE','niters','simfilename','mask','groupingsamples','opt');
        conn_disp('fprintf','Created file %s\n',fullfile(file2_path,[file2_name,'.mat']));
        fh=fopen(fullfile(file2_path,[file2_name,'.m']),'wt');
        fprintf(fh,'opt=[];\n');
        fprintf(fh,'load %s;\n',fullfile(file2_path,[file2_name,'.mat']));
        fprintf(fh,'a=spm_vol(Y);\n');
        fprintf(fh,'y=[];\n');
        fprintf(fh,'if isempty(mask)\n');
        fprintf(fh,'    y=spm_read_vols(a);\n');
        fprintf(fh,'    mask=~any(isnan(y),4)&~all(y==0,4);\n');
        fprintf(fh,'    y=reshape(y,size(y,1)*size(y,2)*size(y,3),size(y,4));\n');
        fprintf(fh,'    y=y(mask,:);\n');
        fprintf(fh,'end\n');
        fprintf(fh,'[i,j,k]=ind2sub(size(mask),find(mask));\n');
        fprintf(fh,'xyz=[i(:) j(:) k(:)]'';\n');
        fprintf(fh,'if isempty(y), y=spm_get_data(a,xyz)''; end\n');
        fprintf(fh,'y=permute(reshape(y,size(y,1),size(X,1),[]),[2,3,1]);\n');
        fprintf(fh,'conn_randomise(X,y,c,m,opt,THR,THR_TYPE,SIDE,niters,simfilename,[],xyz,[],groupingsamples);\n');
        fclose(fh);
        conn_disp('fprintf','Created file %s\n',fullfile(file2_path,[file2_name,'.m']));
        ok=false;
    end
else
    ok=true;
end

    function conn_vproject_randomise_nowlater(varargin)
        v3=get(h3,'value');
        if v3==1, set([h4,h5],'visible','off'); else set([h4,h5],'visible','on'); end
    end
end



