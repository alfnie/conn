function [dataM,dataS,dataA,txtA,dataD,dataDmax,dataDidx]=conn_qaplotsexplore_readfiles(files_jpg,files_txt,in, keeprawdata, dogui, skipall);
if nargin<4||isempty(keeprawdata), keeprawdata=true; end
if nargin<5||isempty(dogui), dogui=true; end
if nargin<6||isempty(skipall), skipall=false; end

txtA={};
dataM=0;
dataS=0;
dataD=[];
dataDmax=[];
dataDidx=[];
if any(conn_server('util_isremotefile',files_jpg)), 
    hmsg=conn_msgbox(sprintf('Loading %d images. Please wait...',numel(in)),'',-1);
    [dataM,dataS,dataA,txtA,dataD,dataDmax,dataDidx]=conn_server('run',mfilename,conn_server('util_localfile',files_jpg),conn_server('util_localfile',files_txt), in, false, false, skipall); 
    if ishandle(hmsg), delete(hmsg); end
    dataA=conn_server('util_remotefile',dataA);
    return; 
end
if dogui, ht=conn_waitbar(0,sprintf('Loading %d plots. Please wait...',numel(in)),false); end
if keeprawdata&&~skipall, dataA=[];
else dataA={};
end

for n=1:numel(in),
    if skipall, dataA{n}=files_jpg{in(n)};
    else
        data=conn_fileutils('imread',files_jpg{in(n)});
        if isa(data,'uint8'), data=double(data)/255; end
        if keeprawdata
            if isempty(dataA), dataA=zeros([size(data,1),size(data,2),size(data,3),numel(in)]); dataM=zeros([size(data,1),size(data,2),size(data,3)]); dataS=dataM; end
            if size(data,1)>size(dataA,1), dataA(size(data,1),1,1,1)=0; dataM(size(data,1),1,1)=0; dataS(size(data,1),1,1)=0; end
            if size(data,2)>size(dataA,2), dataA(1,size(data,2),1,1)=0; dataM(1,size(data,2),1)=0; dataS(1,size(data,2),1)=0; end
            if size(data,3)>size(dataA,3), dataA(1,1,size(data,3),1)=0; dataM(1,1,size(data,3))=0; dataS(1,1,size(data,3))=0; end
            if size(dataA,1)>size(data,1), data(size(dataA,1),1,1,1)=0; end
            if size(dataA,2)>size(data,2), data(1,size(dataA,2),1,1)=0; end
            if size(dataA,3)>size(data,3), data(1,1,size(dataA,3),1)=0; end
            %if mean(data(:))<.5, data=1-data; end
            %if size(data,3)==3, data=data.*repmat(shiftdim(figcolor,-1),[size(data,1),size(data,2)]); end
            dataA(:,:,:,n)=data;
        else
            if isempty(dataA), dataM=zeros([size(data,1),size(data,2),size(data,3)]); dataS=dataM; end
            if size(data,1)>size(dataM,1), dataM(size(data,1),1,1)=0; dataS(size(data,1),1,1)=0; end
            if size(data,2)>size(dataM,2), dataM(1,size(data,2),1)=0; dataS(1,size(data,2),1)=0; end
            if size(data,3)>size(dataM,3), dataM(1,1,size(data,3))=0; dataS(1,1,size(data,3))=0; end
            if size(dataM,1)>size(data,1), data(size(dataM,1),1,1,1)=0; end
            if size(dataM,2)>size(data,2), data(1,size(dataM,2),1,1)=0; end
            if size(dataM,3)>size(data,3), data(1,1,size(dataM,3),1)=0; end
            dataA{n}=files_jpg{in(n)};
        end
        dataM=dataM+data;
        dataS=dataS+data.^2;
    end
    descr=''; try, descr = fileread(files_txt{in(n)}); end
    if isempty(descr), txtA{n}={'[empty]'};
    else txtA{n}=regexp(descr,'\n+','split');
    end
    if dogui&&keeprawdata, conn_waitbar(n/numel(in),ht);
    elseif dogui, conn_waitbar(.5*n/numel(in),ht);
    end
end
if ~skipall
    dataM=dataM/numel(in); %mean(dataA,4);
    dataS=sqrt(max(0,dataS/numel(in)-dataM.^2)); %sqrt(mean(dataD.^2,4)); %std(dataA,1,4);
    if 0
        dataD=abs(dataA-repmat(dataM,[1,1,1,size(dataA,4)]));
        %                     flth=1; %conn_hanning(3); flth=flth/sum(flth);
        %                     if numel(flth)==1
        %                         temp=repmat(sum(dataS.^2,3),[1,1,1,size(dataA,4)]);
        %                         dataD=sqrt(sum(dataD.^2,3)./max(eps,.01*mean(temp(:))+temp));
        %                     else
        %                         temp=repmat(convn(convn(sum(dataS.^2,3),flth,'same'),flth','same'),[1,1,1,size(dataA,4)]);
        %                         dataD=sqrt(convn(convn(sum(dataD.^2,3),flth,'same'),flth','same')./max(eps,.01*mean(temp(:))+temp));
        %                     end
        %[dataDmax,dataDidx]=max(dataD,[],4);
    else % skips pre-computation of dataD
        dataD=[];
        dataDmax=-inf(size(dataM));
        dataDidx=nan(size(dataM));
        for n=1:numel(in)
            if keeprawdata
                temp=abs(dataA(:,:,:,n)-dataM); % dataD(:,:,:,n)
            else
                data=conn_fileutils('imread',files_jpg{in(n)});
                if isa(data,'uint8'), data=double(data)/255; end
                if size(dataM,1)>size(data,1), data(size(dataM,1),1,1,1)=0; end
                if size(dataM,2)>size(data,2), data(1,size(dataM,2),1,1)=0; end
                if size(dataM,3)>size(data,3), data(1,1,size(dataM,3),1)=0; end
                temp=abs(data-dataM);
            end
            mask=temp>dataDmax;
            dataDidx(mask)=n;
            dataDmax(mask)=temp(mask);
            if dogui, conn_waitbar(.5+.5*n/numel(in),ht); end
        end
    end
end

if dogui, conn_waitbar('close',ht); end
end
