function data = conn_guibackground(option,varargin)
global CONN_gui CONN_h;
if ~nargin, option='setfile'; end
data=[]; 

switch(option)
    case 'clear',
        CONN_gui.background=[];
        
    case {'setfile', 'setfiledefault','setfilecolor'}
        if strcmp(option,'setfiledefault'),%||strcmp(option,'setfilecolor'), 
            if 0
                filename=fullfile(fileparts(which(mfilename)),'conn_guibackground.jpg');
            else
                k=[100 10 5 0.9 1];
                filename=0;
                for k1=[1 2 4 8 12]
                    k(1)=100*k1;
                    k2=k(2);
                    a=rand(3,500,1000);
                    i=randperm(numel(a)/3);
                    i=i(1:end-k2);
                    a(:,i)=0;
                    a=permute(a,[2 3 1]);
                    a=conn_guibackground_filter(a,k(1));
                    a=a/max(a(:));
                    a=max(0,min(1,conn_bsxfun(@times,1-k(5)*rand(1,1,3),max(0,min(1,k(3)*(a-k(4)*mean(a(:))))))));
                    filename=max(filename,a);
                end
            end
        elseif strcmp(option,'setfilecolor'), filename=fullfile(fileparts(which(mfilename)),'conn_guibackground.jpg');
        elseif nargin>1, filename=varargin{1}; 
        else 
            [filename,filepath]=uigetfile({'*',  'All Files (*)'; '*.jpg;*.tif;*.tiff;*.gif;*.bmp;*.png','image files'},'Select image file:');
            if isequal(filename,0), return; end
            filename=fullfile(filepath,filename);
        end
        if ischar(filename), a=imread(filename);
        else a=filename;
        end
        pos=get(CONN_h.screen.hfig,'position');
        if max(a(:))>1, a=double(a)/255; end
        if ~strcmp(option,'setfilecolor')&~strcmp(option,'setfiledefault')
            answer=conn_menu_inputdlg('Background image smoothing? (set to 0 for no smoothing)','',1,{'20'});
            if ~isempty(answer)&&~isempty(str2num(answer{1})), a=conn_guibackground_filter(a,str2num(answer{1})); end
        end
        if size(a,3)>1,
            if nargin>2, answ=varargin{2};
            elseif strcmp(option,'setfilecolor'),  answ='CONN color theme'; 
            else answ=conn_questdlg('Color scheme?','','True color','CONN color theme','True color');
            end
            if isequal(answ,'CONN color theme'), 
                a=a-mean(a(:));
                a=a/max(abs(a(:)))/2;
                a=max(0,min(1, conn_bsxfun(@plus,shiftdim(CONN_gui.backgroundcolor,-1), conn_bsxfun(@times,shiftdim(max(.05,sqrt(CONN_gui.backgroundcolor.*(1-CONN_gui.backgroundcolor))),-1),a)) ));
            end
        elseif strcmp(option,'setfilecolor'),  
            a=max(0,min(1, conn_bsxfun(@plus,shiftdim(CONN_gui.backgroundcolor,-1),conn_bsxfun(@times,a-.5,min(1,.25*shiftdim(max(.01,CONN_gui.backgroundcolor)/mean(max(.01,CONN_gui.backgroundcolor)),-1))))));
        else
            a=max(0,min(1, conn_bsxfun(@plus,shiftdim(CONN_gui.backgroundcolor,-1),conn_bsxfun(@times,a-.5,min(1,.25*shiftdim(max(.01,CONN_gui.backgroundcolor)/mean(max(.01,CONN_gui.backgroundcolor)),-1))))));
        end
        k=pos(3:4)/max(pos(3)/size(a,2),pos(4)/size(a,1));
        CONN_gui.background=uint8(255*a(round(linspace(1,k(2),pos(4))),round(linspace(1,k(1),pos(3))),:));
        data=shiftdim(CONN_gui.backgroundcolor,-1);
        try
            data=double(CONN_gui.background)/255;
        end
        
    case 'cleartrans',
        CONN_gui.background={};
        
    case 'settrans',
        hfig=CONN_h.screen.hfig;
        pos=get(hfig,'position');
        p3=get(0,'screensize');
        %pos(1:2)=pos(1:2)+p3(1:2)-1; % note: fix issue when connecting to external monitor/projector
        pos(2)=p3(4)-pos(2)-pos(4);
        if 1,%~iscell(data)||numel(data)<2||~isequal(data{2},pos)
            rect = java.awt.Rectangle(pos(1), pos(2), pos(3), pos(4));
            robot = java.awt.Robot;
            set(hfig,'visible','off'); 
            if isfield(CONN_h.screen,'hlog')&&ishandle(CONN_h.screen.hlog), set(CONN_h.screen.hlog,'visible','off'); end
            drawnow; pause(.5);
            jImage = robot.createScreenCapture(rect);
            if isfield(CONN_h.screen,'hlog')&&ishandle(CONN_h.screen.hlog), set(CONN_h.screen.hlog,'visible','on'); end
            set(hfig,'visible','on');
            h = jImage.getHeight;
            w = jImage.getWidth;
            pixelsData = reshape(typecast(jImage.getData.getDataStorage, 'uint8'), 4, w, h);
            img = cat(3, reshape(pixelsData(3, :, :), w, h)', reshape(pixelsData(2, :, :), w, h)', reshape(pixelsData(1, :, :), w, h)');
            img = max(0,min(1, bsxfun(@plus,0*shiftdim(CONN_gui.backgroundcolor,-1),1*double(img)/255)));
            %img = max(0,min(1, bsxfun(@plus,.75*shiftdim(CONN_gui.backgroundcolor,-1),.25*double(img)/255)));
            answer=conn_menu_inputdlg('Background image smoothing? (set to 0 for no smoothing)','',1,{'20'});
            if ~isempty(answer)&&~isempty(str2num(answer{1})), img=conn_guibackground_filter(img,str2num(answer{1})); end
            answ=conn_questdlg('Color scheme?','','True color','CONN color theme','True color');
            if isequal(answ,'CONN color theme'),
                img=img-mean(img(:));
                img=img/max(abs(img(:)))/2;
                img=max(0,min(1, conn_bsxfun(@plus,shiftdim(CONN_gui.backgroundcolor,-1), conn_bsxfun(@times,shiftdim(max(.05,sqrt(CONN_gui.backgroundcolor.*(1-CONN_gui.backgroundcolor))),-1),img)) ));
%                 img=mean(double(img),3);
%                 img=max(0,min(1, conn_bsxfun(@plus,shiftdim(CONN_gui.backgroundcolor,-1),conn_bsxfun(@times,img-mean(img(:)),min(1,.5*shiftdim(max(.01,CONN_gui.backgroundcolor)/mean(max(.01,CONN_gui.backgroundcolor)),-1))))));
            end
            CONN_gui.background=uint8(255*img);
        end
        data=shiftdim(CONN_gui.backgroundcolor,-1);
        try
            data=double(CONN_gui.background)/255;
        end
        %CONN_gui.background{2}=pos;
        
    case 'get',
        tpos=varargin{1};
        pos=get(CONN_h.screen.hfig,'position');
        tpos(1:2)=tpos(1:2).*([size(CONN_gui.background,2),size(CONN_gui.background,1)]./pos(3:4));
        tpos(3:4)=tpos(3:4).*([size(CONN_gui.background,2),size(CONN_gui.background,1)]./pos(3:4));
        %p3=get(0,'screensize');
        %tpos(1:2)=tpos(1:2)+p3(1:2)-1; % note: fix issue when connecting to external monitor/projector
        tsiz=varargin{2};
        data=shiftdim(CONN_gui.backgroundcolor,-1);
        try
            data=double(CONN_gui.background(...
                max(1,min(size(CONN_gui.background,1), round(linspace(size(CONN_gui.background,1)+2-tpos(2)-tpos(4),size(CONN_gui.background,1)+1-tpos(2),tsiz(1))) )),...
                max(1,min(size(CONN_gui.background,2), round(linspace(tpos(1),tpos(1)+tpos(3),tsiz(2))))), ...
                :))/255;
        end
end
end

function a=conn_guibackground_filter(a,fh)
a=cat(1,a,a(end:-1:1,:,:,:));
a=cat(2,a,a(:,end:-1:1,:));
[i1,i2]=ndgrid(1:size(a,1),1:size(a,2)); 
i1=min(i1-1,size(a,1)+1-i1)/size(a,1); 
i2=min(i2-1,size(a,2)+1-i2)/size(a,2);
f=exp(-fh*sqrt(i1.^2+i2.^2)/max([i1(:);i2(:)]));
a=abs(a).^2; % minutephysics fix
a=real(ifft2(conn_bsxfun(@times,f,fft2(a))));
a=sqrt(max(0,a));
a=a(1:end/2,1:end/2,:);
end
