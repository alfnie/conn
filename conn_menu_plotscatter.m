function handles=conn_menu_plotscatter(x,y,varargin)
% CONN_MENU_PLOTSCATTER scatter plot with confidence intervals
%    dots       : x/y data
%    black line : linear regression fit
%    inner area : 95% CI of the mean
%    outer area : 95% CI of the data
%
% handles = conn_menu_plotscatter(x,y)
% returns handles to plot objects
%
% e.g. 
% x=randn(10,1);
% y=x+randn(10,1);
% conn_menu_plotscatter(x,y)
%
mask=~isnan(x)&~isnan(y);
x=x(mask); y=y(mask);
N=numel(x);
X=[ones(N,1),x(:)];
B=X\y(:);
yfit=X*B;
mx=mean(x(:));
sx=sum(abs(x(:)-mx).^2);
sy=sum(abs(y(:)-yfit).^2);

mm=[min(x), max(x)];
tx=linspace(mm*[1.1;-.1],mm*[-.1;1.1],100)';
ty=[ones(100,1) tx]*B;
RMSE=sqrt(sy/(N-2));
SE1=RMSE*sqrt(1/N + (tx-mx).^2 ./ sx);
SE2=RMSE*sqrt(1 + 1/N + (tx-mx).^2 ./ sx);
k=spm_invTcdf(.975,N-2); % note: 95% ci's

cla;
hold on;
color=.7-.2*exp(-N/10);
handles.CIdata=patch([tx;flipud(tx)], [ty+k*SE2;flipud(ty-k*SE2)],'w','facecolor',[.925 .925 1],'edgecolor','none','facealpha',.5);
handles.CImean=patch([tx;flipud(tx)], [ty+k*SE1;flipud(ty-k*SE1)],'w','facecolor',[.6 .6 1],'edgecolor','none','facealpha',.5);
handles.data=plot(x,y,'ko','markerfacecolor',color*[1 1 1],'markeredgecolor',color*[1 1 1]);
handles.fit=plot(tx,ty,'b-','color',[.2 .2 1],'linewidth',3);
hold off;
end
