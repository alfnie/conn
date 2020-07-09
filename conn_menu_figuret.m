function [hfig,himg]=conn_menu_figuret(varargin)

keys=find(cellfun(@ischar,varargin));
doremove=[];
match=keys(strcmp(varargin(keys),'color'));
if any(match), mixcolor=varargin{match(end)+1}; doremove=[doremove match(:)' match(:)'+1];
else mixcolor=[1 1 1]; 
end
match=keys(strcmp(varargin(keys),'position'));
if any(match), position=varargin{match(end)+1}; doremove=[doremove match(:)' match(:)'+1];
    match=find(strcmp(varargin(keys(keys<match(end))),'units'));
    if ~isempty(match)&&isempty(strmatch(lower(varargin(match(end)+1)),'normalized')), error('only normalized units supported'); end
else position=[.3 .3 .4 .4]; 
end
match=keys(strcmp(varargin(keys),'mixture'));
if any(match), mixcolorw=varargin{match(end)+1}; doremove=[doremove match(:)' match(:)'+1];
else mixcolorw=.5; 
end
varargin(unique(doremove))=[];
p3=get(0,'screensize');
pos=position.*[p3(3:4) p3(3:4)];
pos(2)=p3(4)-pos(2)-pos(4);
rect = java.awt.Rectangle(pos(1), pos(2), pos(3), pos(4)-0);
robot = java.awt.Robot;
jImage = robot.createScreenCapture(rect);
h = jImage.getHeight;
w = jImage.getWidth;
pixelsData = double(reshape(typecast(jImage.getData.getDataStorage, 'uint8'),4,w,h))/255;
if ~isempty(mixcolor), img = cat(3, mixcolorw*mixcolor(1)+(1-mixcolorw)*reshape(pixelsData(3,:,:),w,h)', mixcolorw*mixcolor(2)+(1-mixcolorw)*reshape(pixelsData(2,:,:),w,h)', mixcolorw*mixcolor(1)+(1-mixcolorw)*reshape(pixelsData(1,:,:),w,h)');
else img = cat(3, reshape(pixelsData(3,:,:),w,h)', reshape(pixelsData(2,:,:),w,h)', reshape(pixelsData(1,:,:),w,h)');
end
img=conn_menu_figuret_filter(img,20); 
hfig = figure(varargin{:},'units','norm','position',position,'menubar','none','resize','off');
hax=axes('units','norm','position',[0 0 1 1]);
himg=image(img,'parent',hax);
set(hax,'visible','off');


%placeholder for location change callback: h = addlistener(hfig,'LocationChanged',@(varargin)disp(randn));
%set(hfig,'windowButtonMotionFcn',@(varargin)disp(randn));
end

function a=conn_menu_figuret_filter(a,fh)
a=cat(1,a,a(end:-1:1,:,:,:));
a=cat(2,a,a(:,end:-1:1,:));
[i1,i2]=ndgrid(1:size(a,1),1:size(a,2)); 
i1=min(i1-1,size(a,1)+1-i1)/size(a,1); 
i2=min(i2-1,size(a,2)+1-i2)/size(a,2);
f=exp(-fh*sqrt(i1.^2+i2.^2)/max([i1(:);i2(:)]));
a=real(ifft2(conn_bsxfun(@times,f,fft2(a))));
a=a(1:end/2,1:end/2,:);
end
