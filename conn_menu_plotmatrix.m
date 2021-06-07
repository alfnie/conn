function [hstruct,h]=conn_menu_plotmatrix(x,varargin)

hstruct=[]; h=[];
if ischar(x)
    if numel(varargin)<1||isempty(varargin{1}), hfig=gcf; else hfig=varargin{1}; end
    if numel(varargin)<2||isempty(varargin{2}), nr=20; else nr=varargin{2}; end
    if numel(varargin)<3||isempty(varargin{3}), pos=[.1 .4 .8 .2]; else pos=varargin{3}; end
    if numel(nr)>1, np=nr(2:end); nr=nr(1);
    else np=5; 
    end
    SERIOUSNESS=0; %0,1,2
    if ~ishandle(hfig), return; end
    if numel(np)>2, nv=np(3); np=np(1:2); 
    else nv=inf; 
    end
    if numel(np)==1, np=[1 np]; end
    if strcmp(get(hfig,'type'),'axes'), h=hfig;
    else h=axes('units','norm','position',pos,'parent',hfig);
    end
    if ~SERIOUSNESS, cmap=[.2 .2 .6;.8 .8 .8]; cshape='o';
    else cmap=[.2 .2 .2;.8 .8 .8]; cshape='s';
    end
    axis(h,'equal','off');
    if SERIOUSNESS<2
        for n=1:abs(nr)+0,
            if n>0&&nr>0, pause(.01); drawnow; end
            if isnan(nv), tnv=n; else tnv=nv; end
            conn_menu_plotmatrix(.5*rand(np),'dogray',tnv,'colormap',cmap,'colormapcdata',ceil(2*rand(np)),'shape',cshape,'parent',h,varargin{4:end});
            %hstruct=conn_menu_plotmatrix(1*(n/(abs(nr)+0))^2*ones(np),'colormap',[.2 .2 .2;.8 .8 .8],'colormapcdata',ceil(2*rand(np)));
            %hstruct.vertices(:,3)=n/100;
            %patch(hstruct,'facecolor','flat','edgecolor','none','parent',h,varargin{3:end});
            %axis(h,'equal','tight','off');
        end
    end
    drawnow;
    %delete(h(ishandle(h)));
    return
end

scaled=0;
while nargin>1&&any(strcmp(varargin(1:2:end),'scaled'))
    idx=find(strcmp(varargin(1:2:end),'scaled'),1);
    scaled=varargin{2*idx};
    varargin(2*idx-1+(0:1))=[];
end
scaletoarea=0;
while nargin>1&&any(strcmp(varargin(1:2:end),'scaletoarea'))
    idx=find(strcmp(varargin(1:2:end),'scaletoarea'),1);
    scaletoarea=varargin{2*idx};
    varargin(2*idx-1+(0:1))=[];
end
cmap=[0 0 1;0 0 0;1 0 0];
while nargin>1&&any(strcmp(varargin(1:2:end),'colormap'))
    idx=find(strcmp(varargin(1:2:end),'colormap'),1);
    cmap=varargin{2*idx};
    varargin(2*idx-1+(0:1))=[];
end
colormapcdata=[];
while nargin>1&&any(strcmp(varargin(1:2:end),'colormapcdata'))
    idx=find(strcmp(varargin(1:2:end),'colormapcdata'),1);
    colormapcdata=varargin{2*idx};
    varargin(2*idx-1+(0:1))=[];
end
dogray=inf;
while nargin>1&&any(strcmp(varargin(1:2:end),'dogray'))
    idx=find(strcmp(varargin(1:2:end),'dogray'),1);
    dogray=varargin{2*idx};
    varargin(2*idx-1+(0:1))=[];
end
dotext=false;
while nargin>1&&any(strcmp(varargin(1:2:end),'dotext'))
    idx=find(strcmp(varargin(1:2:end),'dotext'),1);
    dotext=varargin{2*idx};
    varargin(2*idx-1+(0:1))=[];
end
shape='s';
while nargin>1&&any(strcmp(varargin(1:2:end),'shape'))
    idx=find(strcmp(varargin(1:2:end),'shape'),1);
    shape=varargin{2*idx};
    varargin(2*idx-1+(0:1))=[];
end
signed=false;
while nargin>1&&any(strcmp(varargin(1:2:end),'signed'))
    idx=find(strcmp(varargin(1:2:end),'signed'),1);
    signed=varargin{2*idx};
    varargin(2*idx-1+(0:1))=[];
end
scalesigned=.75;
while nargin>1&&any(strcmp(varargin(1:2:end),'scalesigned'))
    idx=find(strcmp(varargin(1:2:end),'scalesigned'),1);
    scalesigned=varargin{2*idx};
    varargin(2*idx-1+(0:1))=[];
end
x0=x;
if scaletoarea, x=sign(x).*sqrt(abs(x)); end
if scaled, x=x/max(eps,max(abs(x(:)))); end

switch(shape)
    case 's'
        sqx=shiftdim([0 -.5 -.5 .5 .5],-1);
        sqy=shiftdim([0 -.5 .5 .5 -.5],-1);
        sqi=shiftdim([1 2 3; 1 3 4; 1 4 5; 1 5 2],-1);
    case 'o'
        n=16;
        sqx=shiftdim([0 .6*cos(2*pi*(0:n-1)/n)],-1);
        sqy=shiftdim([0 .6*sin(2*pi*(0:n-1)/n)],-1);
        sqi=shiftdim([ones(n,1),(2:n+1)',[(3:n+1) 2]'],-1);
    case 'hdo'
        n=100;
        sqx=shiftdim([0 .6*cos(2*pi*(0:n-1)/n)],-1);
        sqy=shiftdim([0 .6*sin(2*pi*(0:n-1)/n)],-1);
        sqi=shiftdim([ones(n,1),(2:n+1)',[(3:n+1) 2]'],-1);
end
sqrx=1:size(x,2);
sqry=(1:size(x,1))';
tempx=conn_bsxfun(@times,max(signed,abs(x)),sqx);
tempy=conn_bsxfun(@times,abs(x),sqy);
if signed, tempx=scalesigned*tempx; tempy=conn_bsxfun(@times,sign(x),max(0,tempy)); end
hstruct=struct(...
    'vertices',[reshape(conn_bsxfun(@plus, sqrx,tempx),[],1) reshape(conn_bsxfun(@plus, sqry,tempy),[],1)],...
    'faces',reshape(conn_bsxfun(@plus,(1:numel(x))',numel(x)*(sqi-1)),[],3),...
    'facevertexcdata', []);
if ~isequal(dotext,false)
    h=text(reshape(repmat(sqrx,size(x,1),1),1,[]), ...
           reshape(repmat(sqry,1,size(x,2)),1,[]), ... 
           ones(1,numel(x0)), ...
           reshape(arrayfun(@num2str,x0,'uni',0),1,[]));
    set(h,'horizontalalignment','center');
    set(h(x0>0),'color','k');
    set(h(x0<0),'color','w');
    if iscell(dotext), set(h,dotext{:}); end
end
if ~isempty(colormapcdata), facevertexcdata=cmap(max(1,min(size(cmap,1),colormapcdata)),:);
else x(isnan(x))=0; facevertexcdata=cmap(2+sign(x),:);
end
facevertexcdata(1:size(facevertexcdata,1)>dogray,:)=repmat(.9+.1*mean(facevertexcdata(1:size(facevertexcdata,1)>dogray,:),2),1,3);
hstruct.facevertexcdata=repmat(facevertexcdata,size(sqi,2),1);
if nargout~=1, try, h=patch(hstruct,'facecolor','flat','edgecolor','none',varargin{:}); end; end
end
