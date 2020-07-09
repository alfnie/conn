function SPM=conn_est_smoothness(SPM)

cwd=pwd;
cd(SPM.swd);
% modified smoothness estimation based on inverse (predictive) second-level model
[Ny,My]=size(SPM.xY.VY);
[nScan nBeta] = size(SPM.xX.X);
if isfield(SPM,'xCon'), X=SPM.xX.X*SPM.xCon(1).c;
else X=SPM.xX.X; 
end
X=X(1:size(X,1)/My,1:size(X,2)/My);
[Nx,Mx]=size(X);
Mr=My;
MAXRES   = spm_get_defaults('stats.maxres');
nSres    = min(Nx*Mr,MAXRES);
DIM=SPM.xVol.DIM;
M=SPM.xVol.M;
XYZ=SPM.xVol.XYZ;
erdf=SPM.xX.erdf;
VM=SPM.VM;
% VResMS = struct('fname',    'ResMS.img',...
%     'dim',      DIM',...
%     'dt',       [spm_type('float64') spm_platform('bigend')],...
%     'mat',      M,...
%     'pinfo',    [1 0 0]',...
%     'descrip',  'spm_spm:Residual sum-of-squares');
% VResMS = spm_create_vol(VResMS);
i_res = round(linspace(1,Nx*Mr,nSres))';        % Indices for residual
VResI(1:nSres) = deal(struct(...
    'fname',    [],...
    'dim',      DIM',...
    'dt',       [spm_type('float64') spm_platform('bigend')],...
    'mat',      M,...
    'pinfo',    [1 0 0]',...
    'descrip',  'spm_spm:StandardisedResiduals'));
for i = 1:nSres
    VResI(i).fname   = sprintf('ResI_%04d.img', i);
    VResI(i).descrip = sprintf('spm_spm:ResI (%04d)', i);
end
VResI = spm_create_vol(VResI);

idxplanes=sparse(1:size(XYZ,2),round(XYZ(3,:)),1,size(XYZ,2),DIM(3));
h=conn_waitbar(0,'re-estimating spatial non-sphericity');
Cxx=X'*X;
iCxx=pinv(Cxx);
for nplane=1:DIM(3)
    
    Q=find(idxplanes(:,nplane));
    Qidx=sub2ind(DIM(1:2),round(XYZ(1,Q)),round(XYZ(2,Q)));
    nvox=numel(Q);
    % Estimates transformed residuals
    %------------------------------------------------------------------
    if nvox
        Y=reshape(spm_get_data(SPM.xY.VY(:),XYZ(:,Q)),[Ny,My,nvox]);
        res=nan([Nx,Mr,nvox]);
        for nvox=1:nvox
            y=Y(:,:,nvox);
            b=iCxx*X'*y;
            e=y-X*b;
            res(:,:,nvox)=e*e(1:Mr,:)';
        end
        res=reshape(res,[Nx*Mr,nvox]);
        ResSS = sum(res.^2,1);                   %-Residual SSQ
        CrResI=res(i_res,:);
        CrResSS=ResSS;
    end    
    %-Write standardised residual images
    %------------------------------------------------------------------
    jj = NaN([DIM(1),DIM(2)]);
    for i = 1:nSres
        if nvox, jj(Qidx) = CrResI(i,:)./sqrt(CrResSS/(erdf/My*Mr)); end
        VResI(i) = spm_write_plane(VResI(i), jj, nplane);
    end
    %-Write ResSS into ResMS (variance) image scaled by tr(RV) above
    %------------------------------------------------------------------
%     if ~isempty(Q), jj(Qidx) = CrResSS/(erdf/My*Mr); end
%     VResMS  = spm_write_plane(VResMS, jj, nplane);
    conn_waitbar(nplane/DIM(3),h);
end
close(h);

if nargout('spm_est_smoothness')==3
    [FWHM,VRpv,R] = spm_est_smoothness(VResI,VM,[Nx*Mr erdf/My*Mr]);
else
    [FWHM,VRpv] = spm_est_smoothness(VResI,VM,[Nx*Mr erdf/My*Mr]);
    R=spm_resels_vol(VM,FWHM)';
end

SPM.xVol.FWHM = FWHM;
SPM.xVol.VRpv = VRpv;
SPM.xVol.R    = R;

j = spm_select('List',SPM.swd,'^ResI_.{4}\..{3}$');
for  k = 1:size(j,1)
    spm_unlink(deblank(j(k,:)));
end
cd(cwd);

end


