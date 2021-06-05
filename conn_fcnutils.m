function varargout=conn_fcnutils(option,varargin)
% internal function

switch(lower(option))
    case 'demean'
        if any(conn_server('util_isremotevar',varargin{1})), [varargout{1:nargout}]=conn_server('run_keep',mfilename,option,varargin{:}); return; end
        x=varargin{1};
        y=x-repmat(mean(x,1),size(x,1),1);
        varargout={y};
    case 'residual'
        if any(conn_server('util_isremotevar',varargin(1:3))), [varargout{1:nargout}]=conn_server('run_keep',mfilename,option,varargin{:}); return; end
        x=varargin{1};
        y=varargin{2};
        B=varargin{3};
        y=y-x*B;
        varargout={y};
        
    case 'zcorr' %                     [z0,z1,d0, tempA, tempB, dof0, dof1, dof2]=conn_fcnutils('zcorr', X1, xf, CONN_x.Preproc.despiking, CONN_x.Preproc.filter, maxrt); 
        if any(conn_server('util_isremotevar',varargin(1:2))), [varargout{1:nargout}]=conn_server('run',mfilename,option,varargin{:}); return; end
        [X1,xf,Preproc.despiking,Preproc.filter,maxrt]=deal(varargin{1:5});
        x0=X1.sampledata;
        if isfield(X1,'samplexyz')&&numel(X1.samplexyz)==size(x0,2), xyz=cell2mat(X1.samplexyz);
        else xyz=nan(3,size(x0,2));
        end
        %x0=detrend(x0);
        x0orig=x0;
        x0=x0-repmat(mean(x0,1),size(x0,1),1);
        maskx0=~all(abs(x0)<1e-4,1)&~any(isnan(x0),1);
        x0=x0(:,maskx0);
        x0orig=x0orig(:,maskx0);
        xyz=xyz(:,maskx0);
        %if isempty(x0),
        %    conn_disp('Warning! No temporal variation in BOLD signal within sampled grey-matter voxels');
        %end
        x1=x0;
        %fy=mean(abs(fft(x0)).^2,2);
        if isfield(Preproc,'despiking')&&Preproc.despiking==1,
            x1=conn_despike(x1);
            %my=repmat(median(x1,1),[size(x1,1),1]);
            %sy=repmat(4*median(abs(x1-my)),[size(x1,1),1]);
            %x1=my+sy.*tanh((x1-my)./max(eps,sy));
        end
        x1=x1-xf*(pinv(xf'*xf)*(xf'*x1));
        if isfield(Preproc,'despiking')&&Preproc.despiking==2,
            x1=conn_despike(x1);
            %my=repmat(median(x1,1),[size(x1,1),1]);
            %sy=repmat(4*median(abs(x1-my)),[size(x1,1),1]);
            %x1=my+sy.*tanh((x1-my)./max(eps,sy));
        end
        [x1,fy]=conn_filter(maxrt,Preproc.filter,x1);
        fy=mean(abs(fy(1:round(size(fy,1)/2),:)).^2,2);
        %dof=max(0,sum(fy)^2/sum(fy.^2)-size(xf,2)); % change dof displayed to WelchSatterthwaite residual dof approximation
        dof0=size(xf,1)-1;
        dof1=max(0,sum(fy)^2/sum(fy.^2)); % WelchSatterthwaite residual dof approximation
        dof2=[size(xf,1), size(xf,2)];
        z0=corrcoef(x0);z1=corrcoef(x1);d0=shiftdim(sqrt(sum(abs(conn_bsxfun(@minus, xyz,permute(xyz,[1,3,2]))).^2,1)),1);
        maskz=z0~=1&z1~=1;
        z0=z0(maskz);z1=z1(maskz);d0=d0(maskz);
        
        if size(x0,2)==size(xyz,2)&&~all(isnan(xyz(:)))
            [nill,idx]=sort(sum(xyz.^2,1),'descend');
            x0=x0(:,idx); x1=x1(:,idx); xyz=xyz(:,idx);
        end
        temp=[x0 nan(size(x0,1),20) x1]';
        temp=.5+.5*temp/max(abs(temp(:)));
        temp(isnan(temp))=0; 
        tempXYZ=[xyz nan(size(xyz,1),20) xyz]';
        if size(temp,1)>512, 
            temp=temp(round(linspace(1,size(temp,1),512)),:); 
            tempXYZ=tempXYZ(round(linspace(1,size(temp,1),512)),:);
        end

        %tempA=ind2rgb(round(1+(size(CONN_h.screen.colormap,1)/2-1)*temp),CONN_h.screen.colormap);
        %temp=repmat(temp,[1,1,3]);
        %temp=cat(2,temp, nan(size(temp,1),10,3), cat(1,repmat(shiftdim(.75/2*[1,1,1],-1),size(x0,2),10), nan(10,10,3), repmat(shiftdim(1/2*[1,1,0],-1),size(x1,2),10)));
        tempB=[mean(x0orig,2) mean(x1,2)];
        [a0,b0]=hist(z0(:),linspace(-1,1,100));[a1,b1]=hist(z1(:),linspace(-1,1,100));
        if all(isnan(d0)), scatterplotdata={};
        else
            th0=conn_hanning(255); th0=th0/sum(th0); [nill,tidx]=sort(d0(:)); t0=convn(z0(tidx),th0,'valid'); t1=convn(z1(tidx),th0,'valid'); td0=convn(d0(tidx),th0,'valid');
            scatterplotdata={{t0(1:50:end) t1(1:50:end) z0 z1},{td0(1:50:end) td0(1:50:end) d0 d0}};
        end
        z0=struct('mean',mean(z0(:)),'std',std(z0(:)));
        z1=struct('mean',mean(z1(:)),'std',std(z1(:)));
        varargout={scatterplotdata,a0,b0,a1,b1,z0,z1, temp,tempB,tempXYZ, dof0,dof1,dof2};
end
