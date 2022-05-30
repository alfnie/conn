function varargout = conn_menu_plothist(x,h,varargin)
% [P,x] = conn_menu_plothist(X,h)
%  X   : column vector of N samples or NxM matrix
%  h   : kernel width
%  P   : kernel density histogram estimate pdf(x|X,h)
%        pdf(x | X,h) = mean( exp(-1/2 * (x-X).^2/h^2) ) / sqrt(2*pi*h^2)

options=struct(...
    'convert','none',...
    'scale',.80,...
    'offset',0,...
    'nsamples',1e2,...
    'plotsamples',true,...
    'plotmeans',true,...
    'plotmedians',false,...
    'colors',[],...
    'edgecolor','none');
if numel(varargin)>0, for n=2:2:numel(varargin), assert(isfield(options,varargin{n-1}),'unrecognized option %s',varargin{n-1}); options.(varargin{n-1})=varargin{n}; end; end
varargout=cell(1,nargout);

if iscell(x), nx=numel(x);
else nx=size(x,2);
end
Pall=[];Xall=[];
for nvar=1:nx
    if iscell(x), tx=x{nvar};
    else tx=x(:,nvar);
    end
    if ~iscell(tx), tx={tx}; end
    tx=cellfun(@(x)x(~isnan(x)),tx,'uni',0);
    if isempty(options.colors), options.colors=repmat(linspace(.25,.75,numel(tx))',1,3); end
    if isequal(options.convert,'logit100'), tx=cellfun(@(x)2*atanh(2*max(.5/numel(x),min(1-.5/numel(x),x/100))-1),tx,'uni',0); end
    minx=min(cellfun(@min,tx));
    maxx=max(cellfun(@max,tx));
    if nargin<2||isempty(h), th=(maxx-minx)/sqrt(max(cellfun(@numel,tx)));
    else th=h;
    end
    X=linspace(minx-3*th,maxx+3*th,options.nsamples);
    P=zeros(numel(tx),options.nsamples);
    for ngroup=1:numel(tx)
        ttx=reshape(tx{ngroup},[],1);
        if numel(ttx)>10*options.nsamples
            P(ngroup,:)=hist(ttx,X);
        elseif numel(ttx)<numel(X)
            p=0;
            for n=1:numel(ttx), p=p+exp(-(X-ttx(n)).^2/th^2/2); end
            P(ngroup,:)=p;
        else
            for n=1:numel(X), P(ngroup,n)=sum(exp(-(X(n)-ttx).^2/th^2/2)); end
        end
    end
    if isequal(options.convert,'logit100'), X=100*(1+tanh(X/2))/2; tx=cellfun(@(x)100*(1+tanh(x/2))/2,tx,'uni',0); end
    P=P./sum(P(:));
    
    if nargout>0
        Pall=cat(3,Pall,P');
        Xall=cat(3,Xall,X');
    else
        maxp=max(sum(P,1));
        pleft=-sum(P,1)/maxp/2;
        for n=1:size(P,1)
            p=P(n,:)/maxp;
            patch(nvar+options.offset+options.scale*[pleft,fliplr(pleft+p)],[X,fliplr(X)],'k','facecolor',options.colors(n,:),'edgecolor',options.edgecolor);
            pleft=pleft+p;
        end
        ttx=[]; for n=1:numel(tx), ttx=[ttx,tx{n}(:)']; end
        if options.plotsamples, hold on; plot(nvar+zeros(numel(ttx),1),sort(ttx(:)),'b.'); hold off; end
        if options.plotmeans, hold on; plot(nvar,mean(ttx(:)),'ko','markerfacecolor','k'); hold off; end
        if options.plotmedians, hold on; plot(nvar,median(ttx(:)),'ko','markerfacecolor','k'); hold off; end
        drawnow
    end
end
if nargout>0
    varargout={Pall,Xall};
else
    hold off;
end
