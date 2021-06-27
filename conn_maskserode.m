function Vmask=conn_maskserode(validsubjects,validrois,REDO)
% conn_maskerode Grey/White/CSF mask thresholding/erosion

global CONN_x;

if nargin<1||isempty(validsubjects), validsubjects=1:CONN_x.Setup.nsubjects; end
if nargin<2||isempty(validrois), validrois=1:3; end
if nargin<3||isempty(REDO), REDO='yes'; end
if conn_projectmanager('inserver'), conn_process('maskserode',validsubjects,validrois,REDO); return; end

clear Vmask;
for nsub=validsubjects,
    nsess=CONN_x.Setup.nsessions(min(length(CONN_x.Setup.nsessions),nsub));
    for nroi=validrois,
        if (nroi>3&&~CONN_x.Setup.rois.sessionspecific(nroi))||(nroi<=3&&~CONN_x.Setup.structural_sessionspecific), nsesstemp=1; else nsesstemp=nsess; end
        if nroi>3&&~CONN_x.Setup.rois.subjectspecific(nroi), nsubstemp=1; else nsubstemp=nsub; end
        for nses=1:nsesstemp,
            Vmask{nsub}{nroi}{nses}=CONN_x.Setup.rois.files{nsubstemp}{nroi}{nses}{1};
        end
    end
    FORCERECOMPUTE=false;   % force recomputing eroded file
    SKIPRECOMPUTE=false;    % skip recomputing eroded file if it already exists
    QASTEP2=[0,1];          % QA volume variables: 1 after erosion; 0 before erosion
    for nroi=intersect(validrois,1:3),
        THR=CONN_x.Setup.erosion.binary_threshold(nroi);
        THRTYPE=CONN_x.Setup.erosion.binary_threshold_type(nroi);
        THRGM=CONN_x.Setup.erosion.exclude_grey_matter(nroi);        
        ERODE=CONN_x.Setup.erosion.erosion_steps(nroi);
        NEIGHB=CONN_x.Setup.erosion.erosion_neighb(nroi);
        if ~CONN_x.Setup.structural_sessionspecific, nsesstemp=1; else nsesstemp=nsess; end
        for nses=1:nsesstemp
            Vmask{nsub}{nroi}{nses}=conn_prepend('e',CONN_x.Setup.rois.files{nsub}{nroi}{nses}{1});
            if FORCERECOMPUTE||(~SKIPRECOMPUTE&&strcmp(lower(REDO),'yes'))||~conn_existfile(Vmask{nsub}{nroi}{nses}),
                
                [nill,nill,ext,fnum]=spm_fileparts(Vmask{nsub}{nroi}{nses});
                switch(ext),
                    case {'.img','.nii','.hdr'},
                        fprintf('%s. File %s\n',CONN_x.Setup.rois.names{nroi}, CONN_x.Setup.rois.files{nsub}{nroi}{nses}{1});
                        V0=spm_vol(CONN_x.Setup.rois.files{nsub}{nroi}{nses}{1}); % mask
                        X0=spm_read_vols(V0);
                        Nstep0=nnz(~isnan(X0)&X0~=0);
                        if nroi>1&&~isempty(THRGM)&&~isnan(THRGM)
                            V2=spm_vol(CONN_x.Setup.rois.files{nsub}{1}{nses}{1}); % grey matter exclusion mask
                            [tx,ty,tz]=ndgrid(1:V0(1).dim(1),1:V0(1).dim(2),1:V0(1).dim(3));
                            txyz=V0(1).mat*[tx(:) ty(:) tz(:) ones(numel(tx),1)]';
                            X2=spm_get_data(V2,pinv(V2(1).mat)*txyz);
                            X0(X2>THRGM)=nan;
                            Nstep0b=nnz(~isnan(X0)&X0~=0);
                            fprintf('   Number of voxels after Grey Matter exclusion mask: %d\n', Nstep0b);
                        end
                        if THRTYPE==1, tTHR=THR;
                            idx1=find(X0(:)>tTHR);
                        else
                            sX0=sort(X0(~isnan(X0)&X0~=0));
                            tTHR=sX0(max(1,min(numel(sX0),round(THR*numel(sX0)))));
                            idx1=find(X0(:)>=tTHR);
                        end
                        Nstep1=numel(idx1);
                        fprintf('   Number of voxels after threhsolding: %d\n',Nstep1);
                        if ~Nstep1, error(sprintf('No suprathreshold voxels in ROI file %s (this typically indicates a problem during normalization/segmentation for this subject; please try re-running normalization/segmentation step after manually re-aligning the structural volumes to better match the default template orientation)',CONN_x.Setup.rois.files{nsub}{nroi}{nses}{1}));
                        end
                        if ERODE>0
                            finished=false;
                            if rem(ERODE,1), tERODE=1;
                            else             tERODE=ERODE;
                            end
                            while ~finished
                                [idxx,idxy,idxz]=ind2sub(size(X0),idx1);
                                idxt=find(idxx>tERODE&idxx<size(X0,1)+1-tERODE&idxy>tERODE&idxy<size(X0,2)+1-tERODE&idxz>tERODE&idxz<size(X0,3)+1-tERODE);
                                [idxdx,idxdy,idxdz]=ndgrid(-tERODE:tERODE);
                                idxd=idxdx+size(X0,1)*idxdy+size(X0,1)*size(X0,2)*idxdz;
                                tidx1=idx1(idxt);
                                valt=zeros(size(tidx1));
                                for n1=1:numel(idxd), valt=valt+(X0(tidx1+idxd(n1))<tTHR); end %for n1=1:length(idxt), valt(n1)=sum(sum(sum(X0(idxx(idxt(n1))+(-tERODE:tERODE),idxy(idxt(n1))+(-tERODE:tERODE),idxz(idxt(n1))+(-tERODE:tERODE))<tTHR,3),2),1); end
                                if rem(ERODE,1), 
                                    svalt=sort(valt);
                                    svalt=svalt(min(numel(svalt),round(ERODE*Nstep1)));
                                    finished=svalt>0;
                                else svalt=NEIGHB;
                                    finished=true;
                                end
                                idxt=idxt(valt<=svalt);
                                idx1=idx1(idxt);
                                tERODE=tERODE+1;
                            end
                            Nstep2=numel(idx1);
%                             [idxx,idxy,idxz]=ind2sub(size(X0),idx1);
%                             idxt=find(idxx>ERODE&idxx<size(X0,1)+1-ERODE&idxy>ERODE&idxy<size(X0,2)+1-ERODE&idxz>ERODE&idxz<size(X0,3)+1-ERODE);
%                             for n1=1:length(idxt), if (sum(sum(sum(X0(idxx(idxt(n1))+(-ERODE:ERODE),idxy(idxt(n1))+(-ERODE:ERODE),idxz(idxt(n1))+(-ERODE:ERODE))<tTHR,3),2),1))>NEIGHB, idxt(n1)=0; end; end
%                             idxt=idxt(idxt>0);
%                             idx1=idx1(idxt);
%                             Nstep2=numel(idx1);
                            fprintf('   Number of voxels after erosion: %d (erosion removed %.1f%% of voxels)\n', Nstep2, 100*(1-Nstep2/Nstep1));
                        else Nstep2=Nstep1;
                        end
                        X1=zeros(size(X0));X1(idx1)=1;
                        V0.fname=regexprep(Vmask{nsub}{nroi}{nses},',\d+$','');
                        spm_write_vol(V0,X1);
                        if ~Nstep2, error(sprintf('No suprathreshold voxels in ROI file %s after erosion (original file contained %d suprathreshold voxels)',CONN_x.Setup.rois.files{nsub}{nroi}{nses}{1},Nstep1));
                        end
                        for qastep2=QASTEP2,
                            roiname=regexprep(CONN_x.Setup.rois.names{nroi},'\s','');
                            if qastep2, roiname=[roiname,'_eroded']; end
                            if ~CONN_x.Setup.structural_sessionspecific, qa_name=['QC_',roiname,'_vol'];
                            else qa_name=['QC_',roiname,'_vol_session',num2str(nses)];
                            end
                            qa_icov=find(strcmp(qa_name,CONN_x.Setup.l2covariates.names(1:end-1)),1);
                            if isempty(qa_icov),
                                qa_icov=numel(CONN_x.Setup.l2covariates.names);
                                CONN_x.Setup.l2covariates.names{qa_icov}=qa_name;
                                CONN_x.Setup.l2covariates.descrip{qa_icov}=['CONN Quality Assurance: # of voxels in ',roiname];
                                CONN_x.Setup.l2covariates.names{qa_icov+1}=' ';
                                for tnsub=1:CONN_x.Setup.nsubjects, CONN_x.Setup.l2covariates.values{tnsub}{qa_icov}=nan; end
                            end
                            if qastep2, CONN_x.Setup.l2covariates.values{nsub}{qa_icov}=Nstep2;
                            else CONN_x.Setup.l2covariates.values{nsub}{qa_icov}=Nstep1;
                            end
                        end
                    otherwise,
                        Vmask{nsub}{nroi}{nses}=CONN_x.Setup.rois.files{nsub}{nroi}{nses}{1};
                end
                
            end
        end
    end
end
