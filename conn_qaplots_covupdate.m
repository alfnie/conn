function [dataA,labelA,infoA,lineA]=conn_qaplots_covupdate(X,Xnames,Xdescr,nsubs,qafolder,ht,ht_offset,ht_scale)

if nargin<4, nsubs=1:size(X,1); end
if nargin<5, qafolder=[]; end
if nargin<6, ht=[]; end
if nargin<7||isempty(ht_offset), ht_offset=0; end
if nargin<8||isempty(ht_scale), ht_scale=0; end

if size(X,1)==1, X=repmat(X,2,1); end
sx=sort(X,1); IQ=nan(3,size(X,2)); Nsx=sum(~isnan(sx)); 
for n=1:size(X,2), if Nsx(n)>1, IQ(:,n)=interp1(linspace(0,1,Nsx(n))',sx(1:Nsx(n),n),[.25 .5 .75]'); else IQ(:,n)=sx(find(~isnan(sx(:,n)),1),n); end; end
IQR=max(1e-10,IQ(3,:)-IQ(1,:)); IQL=[IQ(1,:)-1.5*IQR; IQ(3,:)+1.5*IQR]; IQ=[IQL(1,:);IQ;IQL(2,:)];
Ka=-IQL(1,:)./max(eps,IQL(2,:)-IQL(1,:));
Kb=1./max(eps,IQL(2,:)-IQL(1,:));
Xdisp=max(-5,min(6, repmat(Ka,size(X,1),1)+repmat(Kb,size(X,1),1).*X )); % scale to same IQR across all measures, all values between 0 and 1
IQdisp=repmat(Ka,size(IQ,1),1)+repmat(Kb,size(IQ,1),1).*IQ;
kt=min([Xdisp(:);IQdisp(:)]); if kt<0, Xdisp=Xdisp-kt; IQdisp=IQdisp-kt; end %Ka=Ka-kt; Xdisp=repmat(Ka,size(X,1),1)+repmat(Kb,size(X,1),1).*X; IQdisp=repmat(Ka,size(IQ,1),1)+repmat(Kb,size(IQ,1),1).*IQ; end
kt=max([Xdisp(:);IQdisp(:)]); if kt>1, Xdisp=Xdisp/kt; IQdisp=IQdisp/kt; end %Ka=Ka/kt; Kb=Kb/kt; Xdisp=repmat(Ka,size(X,1),1)+repmat(Kb,size(X,1),1).*X; IQdisp=repmat(Ka,size(IQ,1),1)+repmat(Kb,size(IQ,1),1).*IQ; end
Npts=200;
hx=linspace(-.1,1.1,Npts)';
py=zeros([Npts,size(Xdisp)]); for ny=1:size(Xdisp,2), py(:,:,ny)=exp(-.5*(repmat(hx,1,size(Xdisp,1))-repmat(Xdisp(:,ny)',Npts,1)).^2/.002); end
ny=sum(~isnan(py),2); %
py(isnan(py))=0;
hy=cumsum(py,2);
hy=permute(cat(2,-hy(:,end,:),2*hy-repmat(hy(:,end,:),[1,size(hy,2),1])),[1,3,2]);
py=permute(py./repmat(max(eps,sum(py,1)),size(py,1),1),[1,3,2]);
Hdisp=.45*hy./repmat(max(eps,max(max(hy,[],1),[],3)),[size(hy,1),1,size(hy,3)]);
dataA={};labelA={};infoA={};lineA={};
for isub=1:numel(nsubs)
    nsub=nsubs(isub);
    Xthis=X(isub,:);
    Xdispthis=Xdisp(isub,:);
    results_patch={[Hdisp(:,:,end);flipud(Hdisp(:,:,1))]+repmat(1:size(Hdisp,2),2*size(Hdisp,1),1), repmat([hx;flipud(hx)],1,size(Hdisp,2))};
    results_line={(1:size(Hdisp,2))+.95*sum(py(:,:,isub).*(Hdisp(:,:,isub+1)+Hdisp(:,:,isub))/2,1),Xdispthis};
    results_info=struct('Values',Xthis,'ValuesDisplay',Xdispthis,'Interquartiles',IQ,'InterquartilesDisplay',IQdisp,'Subjects',{nsubs},'Variables',{Xnames},'Variables_descr',{Xdescr});
    if iscell(nsub), results_label={nsub}; 
    else results_label={{sprintf('Subject %d',nsub)}}; 
    end
    for n1=1:numel(Xthis), results_label{1}{1+n1}=sprintf('%s = %s',Xnames{n1},num2str(Xthis(n1))); end
    results_str={};
    if ~isempty(qafolder)
        filename=fullfile(qafolder,sprintf('QA_COV.subject%03d.mat',nsub));
        conn_savematfile(filename,'results_patch','results_line','results_info','results_label','results_str');
    end
    if nargout>0
        dataA{isub}=results_patch;
        labelA{isub}=results_label;
        infoA{isub}=results_info;
        lineA{isub}=results_line;
    end
    if ~isempty(ht), conn_waitbar(ht_offset+ht_scale*isub/numel(nsubs),ht);
    else fprintf('.');
    end
end
end