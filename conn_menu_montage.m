function varargout=conn_menu_montage(h,varargin)
if isempty(h), emphx=.75;
else
    if ishandle(h),
        units=get(h,'units');
        set(h,'units','points');
        emphx=get(h,'position');
        set(h,'units',units);
        emphx=emphx(4)/emphx(3);
    elseif ischar(h)
        switch(h)
            case 'xyz2coords'
                nX=varargin{1};
                xyz=varargin{2};
                [refx,refy]=ndgrid(1:nX(1),1:nX(2));
                ref=[(refx(:)-1)*nX(3) (refy(:)-1)*nX(4)]';
                coords=xyz(1:2,:)+ref(:,max(1,min(nX(6),round(xyz(3,:)))));
                varargout={coords};
            case 'coords2xyz'
                nX=varargin{1};
                coords=varargin{2};
                [refx,refy]=ndgrid(1:nX(1),1:nX(2));
                ref=[(refx(:)-1)*nX(3) (refy(:)-1)*nX(4)]';
                ix=ceil(coords(1,:)/nX(3));
                iy=ceil(coords(2,:)/nX(4));
                xyz=[1+mod(coords(1,:)-1,nX(3)); 1+mod(coords(2,:)-1,nX(4)); ix+(iy-1)*nX(1)]; 
                varargout={xyz};
            case 'figure'
                conn_montage_display(varargin{:});

            case 'plotline'
                % eg sx=[3 3 50 70]; clf; hold on; xyz1=[randi(sx(1)*sx(3));randi(sx(2)*sx(4))]; for n1=1:2e2, xyz2=[randi(sx(1)*sx(3));randi(sx(2)*sx(4))]; [x,y]=conn_menu_montage('plotline',sx,xyz1,xyz2); if 0, hp=plot(x,y,'-','linewidth',4); pause; delete(hp); end; c=xyz2'./sx(1:2)./sx(3:4); c=[c 1-mean(c)]; plot(x,y,'-','color',c); plot(xyz1(1),xyz1(2),'o',xyz2(1),xyz2(2),'o','color',c,'markerfacecolor',c); set(gca,'xtick',(0:sx(1))*sx(3),'ytick',(0:sx(1))*sx(3)); grid on; end; axis equal tight off; set(gcf,'color','k')
                nX=varargin{1};
                coords1=varargin{2};
                coords2=varargin{3};
                if numel(varargin)>=4, Interf=varargin{4};
                else Interf=.15;
                end
                Alpha=8;
                Npts=15;
                DoStraight=false;
                EqualSteps=true;
                Smooth=ceil(Interf*Npts*24);
                xyz1=conn_menu_montage('coords2xyz',nX,coords1);
                xyz2=conn_menu_montage('coords2xyz',nX,coords2);
                i1=round(xyz1(3,:));
                i2=round(xyz2(3,:));
                x=zeros(6*Npts,size(coords1,2));
                y=zeros(6*Npts,size(coords1,2));
                [refx,refy]=ndgrid(1:nX(1),1:nX(2));
                ref=[(refx(:)-1)*nX(3) (refy(:)-1)*nX(4)]';
                if DoStraight
                    ix1=1+mod(i1-1,nX(1));
                    iy1=ceil(i1/nX(1));
                    ix2=1+mod(i2-1,nX(1));
                    iy2=ceil(i2/nX(1));
                    mask=(abs(ix1-ix2)+abs(iy1-iy2)<=0) | (sum(abs(coords1-coords2).^2,1)<=min(nX(3:4)).^2);
                    %mask=sum(abs(coords1-coords2).^2,1)<=max(nX(3:4)).^2/4;
                    if ~isempty(mask)
                        w=linspace(1,0,6*Npts)';
                        x(:,mask)=w*coords1(1,mask)+(1-w)*coords2(1,mask);
                        y(:,mask)=w*coords1(2,mask)+(1-w)*coords2(2,mask);
                    end
                    mask=~mask;
                else mask=true(size(i1));
                end
                if any(mask)
                    pivot{1}=coords1(:,mask); % start
                    pivot{7}=coords2(:,mask); % end
                    D=[0 0 nX(3) nX(3); 0 nX(4) 0 nX(4)];
                    p1=repmat(ref(:,i1(mask)),[1,1,4,4])+repmat(permute(D,[1,3,2,4]),[1,nnz(mask),1,4]);
                    p2=repmat(ref(:,i2(mask)),[1,1,4,4])+repmat(permute(D,[1,4,3,2]),[1,nnz(mask),4,1]);
                    da1=abs(p1-repmat(coords1(:,mask),[1,1,4,4]));
                    da2=abs(p2-repmat(coords2(:,mask),[1,1,4,4]));
                    db=sum(abs(p1-p2),1);
                    [nill,idx]=min(reshape(Alpha*min(da1,[],1)+max(da1,[],1)+Alpha*min(da2,[],1)+max(da2,[],1)+db,nnz(mask),16),[],2);
                    [idx1,idx2]=ind2sub([4,4],idx);
                    pivot{3}=ref(:,i1(mask))+D(:,idx1); % start corner
                    pivot{5}=ref(:,i2(mask))+D(:,idx2); % end corner
                    pivot{2}=pivot{1};t=abs(pivot{1}(1,:)-pivot{3}(1,:))<abs(pivot{1}(2,:)-pivot{3}(2,:));
                    pivot{2}(1,t)=pivot{3}(1,t);
                    pivot{2}(2,~t)=pivot{3}(2,~t); % start edge
                    pivot{6}=pivot{7};t=abs(pivot{7}(1,:)-pivot{5}(1,:))<abs(pivot{7}(2,:)-pivot{5}(2,:));
                    pivot{6}(1,t)=pivot{5}(1,t);
                    pivot{6}(2,~t)=pivot{5}(2,~t); % end edge
                    t=abs(pivot{3}(1,:)-pivot{5}(1,:))<abs(pivot{3}(2,:)-pivot{5}(2,:));
                    pivot{4} = [pivot{5}(1,:); pivot{3}(2,:)];
                    pivot{4}(:,t) = [pivot{3}(1,t); pivot{5}(2,t)]; % mid-point
                    t=all(pivot{5}==pivot{3},1)&min(abs(pivot{2}-pivot{6}),[],1)==0;
                    p0=(pivot{6}(:,t)+pivot{2}(:,t))/2;
                    %p0(1,:)=min(1,max(nX(3).*nX(1), p0(1,:)+randn(1,size(p0,2))*nX(3)/16));
                    %p0(2,:)=min(1,max(nX(4).*nX(2), p0(2,:)+randn(1,size(p0,2))*nX(4)/16));
                    pivot{5}(:,t)=p0;
                    pivot{4}(:,t)=p0;
                    pivot{3}(:,t)=p0;
                    for n=1:numel(pivot)-1
                        w=linspace(1,0,Npts)';
                        x(Npts*(n-1)+(1:Npts),mask)=w*pivot{n}(1,:)+(1-w)*pivot{n+1}(1,:);
                        y(Npts*(n-1)+(1:Npts),mask)=w*pivot{n}(2,:)+(1-w)*pivot{n+1}(2,:);
                    end
                end
                if EqualSteps
                    d=[zeros(1,size(x,2));cumsum(max(1e-10,sqrt(diff(x,1,1).^2+diff(y,1,1).^2)),1)];
                    d=d./repmat(max(eps,max(d,[],1)),size(d,1),1);
                    t=linspace(0,1,size(x,1))';
                    for n=1:size(x,2)
                        xy=interp1(d(:,n),[x(:,n) y(:,n)],t);
                        x(:,n)=xy(:,1); y(:,n)=xy(:,2);
                    end
                end
                if Smooth>0
                    h=conn_hanning(Smooth+1);h=h/sum(h);
                    x(2:end-1,:)=max(1,min(nX(3).*nX(1), x(2:end-1,:)+0*randn(size(x,1)-2,size(x,2))*nX(3)/32));
                    y(2:end-1,:)=max(1,min(nX(4).*nX(2), y(2:end-1,:)+0*randn(size(y,1)-2,size(x,2))*nX(4)/32));
                    x=convn([x(1+zeros(Smooth,1),:);x;x(end+zeros(Smooth,1),:)],h,'valid');
                    y=convn([y(1+zeros(Smooth,1),:);y;y(end+zeros(Smooth,1),:)],h,'valid');
                end
                varargout={x,y};
            otherwise, error('unknown option %s',h);
        end
        return;
    else emphx=h;
    end
end
sX=[size(varargin{end}) 1 1 1];
if sX(4)>1
    [n1,n2]=ndgrid(1:sX(4),1:sX(4)); n1=n1(:); n2=n2(:);
    [nill,idx]=max(min(emphx*sX(2)./n1(:), sX(1)./n2(:))-1e10*(n1(:).*n2(:)<sX(4)));
    n1=n1(idx); n2=n2(idx);
else
    n1=1; n2=1;
end
nX=[n2 n1];
[i2,i1]=ind2sub(nX,1:sX(4));
varargout={};
for n=1:numel(varargin),
    x=nan([sX(1)*n1,sX(2)*n2,sX(3)]);
    for m=1:numel(i1)
        x(sX(1)*(i1(m)-1)+(1:sX(1)),sX(2)*(i2(m)-1)+(1:sX(2)),:)=varargin{n}(:,:,:,min(size(varargin{n},4),m));
    end
    varargout{n}=x; %permute(x,[2,1,3]);
end
varargout{end+1}=[nX sX([2 1 3 4])];

end
