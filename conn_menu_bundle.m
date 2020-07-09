function [Lines,WidthLines] = conn_menu_bundle(Lines,WidthLines,ValueLines,MaxBundleLevel,doplot)
if nargin<2||isempty(WidthLines), WidthLines=.2*ones([size(Lines,1),1,size(Lines,3),size(Lines,4)]); end
if nargin<3||isempty(ValueLines), ValueLines=zeros([1,1,size(Lines,3)]); end
if nargin<4||isempty(MaxBundleLevel), MaxBundleLevel=10; end
if nargin<5||isempty(doplot), doplot=false; end
METHODTEST=false;

if METHODTEST
    widthlines=ones([size(Lines,1),1,size(Lines,3)]);
    wp=abs(linspace(-1,1,size(Lines,1))).^10; 
    wp0=repmat(wp,[size(Lines,3),1]);
    coords=shiftdim(cat(2,Lines(1,:,:,1),Lines(end,:,:,1)),1);
    N=size(coords,2);
    d1=sum(abs(coords).^2,1);
    d=sqrt(max(0,repmat(d1,[N,1])+repmat(d1',[1,N])-2*real(coords'*coords)));
    if doplot, hw=conn_waitbar(0,'Bundling. Please wait',false); end
    BundleLevel=size(Lines,4)-1;
    for bundlerep=1.0*BundleLevel+1:1.0*MaxBundleLevel;
        lines=Lines(:,:,:,1); % Lines: points x dimensions x lines x bundle_levels
        lbundl=bundlerep;
        H=exp(-(d/max(eps,20*lbundl)).^2);
        widthlines=shiftdim(sum(H,1),-1);
        H=H./repmat(max(eps,sum(H,2)),[1,size(H,2)]);
        vertices=permute(lines,[3,1,2]);
        p1=wp0;
        p2=repmat(exp(-(abs(coords(1,:)-coords(end,:))'/20).^2),[1,size(vertices,2)]);
        p=repmat(min(1,p1+p2),[1,1,size(vertices,3)]);
        newvertices=vertices;
        newvertices(:,:)=H*vertices(:,:);
        vertices=vertices.*p + newvertices.*(1-p);
        lines=permute(vertices,[2,3,1]);
        if 1,%~rem(bundlerep,10)==1
            Lines=cat(4,Lines,lines);
            widthlines=.1*(.1+conn_bsxfun(@times,double(1-wp'>.01),40*tanh(max(0,widthlines)/40)));
            WidthLines=cat(4,WidthLines,widthlines);
        end
        if doplot, conn_waitbar((bundlerep-1.0*BundleLevel)/(10*MaxBundleLevel-1.0*BundleLevel),hw); end
    end
else
    lines=Lines(:,:,:,end); % Lines: points x dimensions x lines x bundle_levels
    dirlines=Lines(1,:,:,1)-Lines(end,:,:,1);
    distl=sqrt(max(eps,sum(dirlines.^2,2)));
    dirlines=shiftdim(conn_bsxfun(@rdivide,dirlines,distl),1);
    dirlines=dirlines'*dirlines;
    distl=distl/size(lines,1);
    simlines=sum(abs(conn_bsxfun(@minus,ValueLines,permute(ValueLines,[3,2,1]))).^2,2);
    negdirlines=dirlines<0;
    widthlines=zeros([size(lines,1),1,size(lines,3)]);
    wp=linspace(0,1,size(Lines,1))'; wp=[wp 1-wp];
    wp0=1-max(wp,[],2).^50;
    bundlewidth=10; % 10
    dirweight=100;  % 100  0-1000
    if doplot, hw=conn_waitbar(0,'Bundling. Please wait',false); end
    BundleLevel=size(Lines,4)-1;
    for bundlerep=10*BundleLevel+1:10*MaxBundleLevel;
        newlines1=zeros(size(lines));
        for n1=1:size(lines,3)
            tlines=lines;
            tlines(:,:,negdirlines(:,n1))=tlines(end:-1:1,:,negdirlines(:,n1));
            ltemp=conn_bsxfun(@minus,tlines,lines(:,:,n1));
            dlines1=sum(abs(ltemp).^2,2);
            dlines2=sum(diff(ltemp([1 1:end],:,:),1,1).^2,2);
            dlines=dlines1+dirweight*dlines2;
            dlines=conn_bsxfun(@plus,dlines,1e3*simlines(n1,:,:));
            plines=max(0,exp(-.5*dlines/bundlewidth^2)-0);
            widthlines(:,:,n1)=max(1e-3,sum(plines,3));
            newlines1(:,:,n1)=conn_bsxfun(@rdivide,sum(conn_bsxfun(@times,tlines,plines),3),widthlines(:,:,n1));
        end
        
        newlines=newlines1-lines;
        lines=lines+.1*conn_bsxfun(@times,wp0, convn(newlines,conn_hanning(11)/6,'same'));%-max(-.1,min(.1, newlines));
        
        klines=cumsum(cat(1,zeros(1,1,size(lines,3)),sqrt(sum(diff(lines,1,1).^2,2))),1);
        klines=conn_bsxfun(@rdivide,klines,max(eps,klines(end,:,:)));
        for n1=1:size(lines,3),
            lines(:,:,n1)=interp1(klines(:,:,n1), lines(:,:,n1), linspace(0,1,size(lines,1)));
        end
        if ~rem(bundlerep,10)==1
            Lines=cat(4,Lines,lines);
            widthlines=.1*(.1+conn_bsxfun(@times,double(wp0>.01),40*tanh(max(0,widthlines)/40)));
            WidthLines=cat(4,WidthLines,widthlines);
        end
        if doplot, conn_waitbar((bundlerep-10*BundleLevel)/(10*MaxBundleLevel-10*BundleLevel),hw); end
    end
    if doplot, conn_waitbar('close',hw); end
end

end