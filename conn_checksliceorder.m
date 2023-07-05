function varargout = conn_checksliceorder(fname,ntime,doplot,algo)
% [select, stats] = conn_cheksliceorder(fname)
% checks functional data and returns possible slice-acquisition order(s)
%
% e.g. 
%   conn_checksliceorder rest.nii
%
% note1: only evaluates the following slice orders: 'ascending','descending','interleaved (Siemens)','interleaved (Philips)','interleaved (bottom-up)','interleaved (top-down)','interleaved (middle-top)'
% note2: this procedure has approximately a 75% correct identification rate for an individual functional run.
%        Differences between sequential and interleaved orders are large (and typically noticeable in a single session/subject)
%        while differences between ascending and descending orders are more subtle (and typically require pooling across multiple 
%        sessions/subjects to obtain a reliable estimate)
%

if nargin<2||isempty(ntime), ntime=11; end % number of acquisitions to smooth
if nargin<3||isempty(doplot), doplot=false; end
if nargin<4||isempty(algo), algo=[1 1 1]; end
if ~iscell(fname), fname={fname}; end
if any(conn_server('util_isremotefile',fname)), [varargout{1:nargout}]=conn_server('run',mfilename,conn_server('util_localfile',fname),ntime,doplot,algo); return; end

options={'ascending','descending','interleaved (Siemens)','interleaved (Philips)','interleaved (bottom-up)','interleaved (top-down)','interleaved (middle-top)'};
stats=struct('rows',{fname},'cols',{options},'nslices',[],'rmse',[]);
wmethods=[1 4];
nlevels=4;
alpha=.9;
for nfname=1:numel(fname)
    vol=spm_vol(fname{nfname});
    nslices=vol(1).dim(3);
    nslicehalf=floor(ntime/2);
    globalsignal=zeros(nslices,numel(vol),nlevels);
    inplanemotion=zeros(nslices,numel(vol),nlevels);
    pvol=zeros([vol(1).dim,nlevels]);
    svol=zeros(nslices,numel(vol));
    thr=zeros(nslices,nlevels);
    if 1 % init thr
        nvol=1;
        for nslice=1:nslices
            x=spm_slice_vol(vol(nvol),spm_matrix([0 0 nslice]),vol(nvol).dim(1:2),1);
            x(isnan(x))=0;
            for nlevel=1:nlevels
                if nlevel==1, thr1=mean(x(x>=mean(x(:))/8)); end
                thr2=prctile(x(x>=thr1),100*(nlevel-1)/nlevels);
                thr(nslice,nlevel)=thr2;
            end
        end
    end
    for nvol=1:numel(vol)
        %thrnow=mean(thr,1);
        for nslice=1:nslices
            x=spm_slice_vol(vol(nvol),spm_matrix([0 0 nslice]),vol(nvol).dim(1:2),1);
            x(isnan(x))=0;
            maxx=max(x(:));

            % ~brainmask
            for nlevel=1:nlevels
                if nlevel==1, thr1=mean(x(x>=mean(x(:))/8)); end
                thr2=prctile(x(x>=thr1),100*(nlevel-1)/nlevels);
                thr2=alpha*thr(nslice,nlevel)+(1-alpha)*thr2;
                thr(nslice,nlevel)=thr2;
                thr2=min(maxx,thr2);

                %thr2=min(maxx,thrnow(nlevel));
                p=(x>=thr2);
                %globalsignal2(nslice,nvol)=sum(p(:).*x(:))/sum(p(:));
                p=convn(convn(p,conn_hanning(3),'same'),conn_hanning(3)','same');
                globalsignal(nslice,nvol,nlevel)=sum(sum(p.*x))/max(eps,sum(sum(p)));
                dx=(nvol>1)*(p-pvol(:,:,nslice,nlevel));
                [di,dj]=gradient(pvol(:,:,nslice,nlevel));
                inplanemotion(nslice,nvol,nlevel)=sum(sum( di.*dx + 1i*dj.*dx ));
                pvol(:,:,nslice,nlevel)=p;
                if nlevel==1, svol(nslice,nvol)=sum(sum(p)); end
            end
            % ~gray/csf
%                         %imagesc([x,x.*p]); pause;
%                     case 3
%                         p=x>mean(x(:))/8;
%                         globalsignal(nslice,nvol)=mean(x(p>0));
%                     otherwise,
%                         error
%             end
            %p=tanh(x/(mean(x(:))));
            %x=x(x>mean(x(:))/8);
            %globalsignal(nslice,nvol)=mean(x(:));
        end
    end
    wdata=max(0,mean(svol,2).*mean(globalsignal(:,:,1),2)); wdata=wdata/max(eps,sum(wdata));
    gschange=globalsignal-globalsignal(:,[1,1:size(globalsignal,2)-1],:);
    data=cat(3,... % sources
        gschange,...
        inplanemotion);

%     switch(algo(2))
%         case 1
%             data=globalsignal-globalsignal(:,[1,1:size(globalsignal,2)-1]);
%         case 2
%             data=velocity;
%         case 3
%             data=globalsignal-convn(globalsignal(:,max(1,min(size(globalsignal,2), 1-floor(ntime/2):size(globalsignal,2)+floor(ntime/2)))),conn_hanning(ntime)'/sum(conn_hanning(ntime)),'valid');
%             %globalsignal=globalsignal./convn(globalsignal(:,max(1,min(size(globalsignal,2), 1-floor(ntime/2):size(globalsignal,2)+floor(ntime/2)))),conn_hanning(ntime)'/sum(conn_hanning(ntime)),'valid');
%             %globalsignal=convn(globalsignal,((-floor(ntime/2):floor(ntime/2))==0) - conn_hanning(ntime)'/sum(conn_hanning(ntime)),'valid');
%         case 4
%             temp1=globalsignal;
%             for nvol=1:numel(vol)
%                 temp2=sort(temp1(:,max(1,min(size(temp1,2), nvol-nslicehalf:nvol+nslicehalf ))),2);
%                 data(:,nvol)=globalsignal(:,nvol)-temp2(:,nslicehalf+1); % median filter
%             end
%         case 5
%             data=globalsignal-repmat(mean(globalsignal,2),1,size(globalsignal,2));
%     end
    data=data./repmat(max(eps,std(data,1,2)),[1,size(data,2),1]);
    trmse=[];
    for ntest=1:numel(options),
        switch(ntest)
            case 1, sliceorder=1:nslices;        % ascending
            case 2, sliceorder=nslices:-1:1;     % descending
            case 3, sliceorder=[fliplr(nslices:-2:1) fliplr(nslices-1:-2:1)]; % interleaved (Siemens)
            case 4, sliceorder=cell2mat(arrayfun(@(n)n:round(sqrt(nslices)):nslices,1:round(sqrt(nslices)),'uni',0)); % interleaved (Philips)
            case 5, sliceorder=[1:2:nslices 2:2:nslices]; % interleaved (bottom-up)
            case 6, sliceorder=[nslices:-2:1, nslices-1:-2:1]; % interleaved (top-down)
            case 7, sliceorder=round((nslices-(1:nslices))/2 + (rem((nslices-(1:nslices)),2) * (nslices - 1)/2)) + 1; % interleaved (middle-top)
        end
        x=data(sliceorder,:,:);
        %plot(x,'.-'); title(options{ntest}); pause;
        x=reshape(x,[],size(x,3));
        w=repmat(wdata(sliceorder),size(data,2),1);
        w=(w(1:end-1).*w(2:end));
        w=w/sum(w);
        stats.globalsignal{nfname,ntest}=x;
        %trmse(:,ntest)=sqrt(mean(abs(diff(x,1,1)).^2,1))';
        trmse(:,ntest)=sqrt(w'*(abs(diff(x,1,1)).^2))';
    end
    %trmse=(trmse-repmat(min(trmse,[],2),1,size(trmse,2)))./repmat(max(eps, max(trmse,[],2)-min(trmse,[],2)),1,size(trmse,2)); % scale all sources equally 
    stats.rmse(nfname,:)=kron(wmethods,ones(1,nlevels))*trmse;  % root mean square 
    stats.nslices(nfname)=nslices;

%         switch(algo(3))
%             case 1,
%disp(sqrt(mean(abs(diff(x,1,1)).^2)));
%                 stats.rmse(nfname,ntest)=mean(sqrt(mean(abs(diff(x,1,1)).^2)));  % root mean square
%             case 2,
%                 x=x-convn(x(max(1,min(numel(x), 1-nslicehalf:numel(x)+nslicehalf))), conn_hanning(ntime)/sum(conn_hanning(ntime)),'valid');
%                 stats.rmse(nfname,ntest)=sqrt(mean(abs(x).^2));  
%             case 3,
%                 temp1=x;
%                 for nvol=1:numel(temp1)
%                     temp2=sort(temp1(max(1,min(numel(temp1), nvol-nslicehalf:nvol+nslicehalf ))));
%                     x(nvol)=x(nvol)-temp2(nslicehalf+1); % median filter
%                 end
%                 stats.rmse(nfname,ntest)=sqrt(mean(abs(x).^2));  
%             otherwise,
%                 error
%         end
    if 0,%doplot|~nargout
        [nill,idx]=sort(stats.rmse(nfname,:));
        fprintf('%s (%d slices), RMSE (lower is better) : ',fname{nfname},nslices);
        for n2=reshape(idx,1,[]),
            fprintf('%s=%f ; ',options{n2},stats.rmse(nfname,n2));
        end
        fprintf('\n');
    end
end
[nill,select]=min(stats.rmse,[],2);
varargout={options(select), stats};

if doplot|~nargout
    if numel(fname)>1,
        [uselect,nill,iselect]=unique(select);
        nselect=accumarray(iselect(:),1);
        [nill,idx]=sort(-nselect);
        fprintf('Summary : ');
        for n2=reshape(idx,1,[]),
            fprintf('%s (%d votes) ; ',options{uselect(n2)},nselect(n2));
        end
        fprintf('\n');
    end
    %figure;
    %plot(stats.rmse','o-'); 
    %set(gca,'xtick',1:numel(options),'xticklabel',options);
end