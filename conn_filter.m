function [y,fy]=conn_filter(rt,filter,x,option,cols)

if nargin<5, cols=[]; end
if nargin<4||isempty(option), option='full'; end
if conn_server('util_isremotevar',x), if nargout>1, [y,fy]=conn_server('run_keep',mfilename,rt,filter,x,option,cols); else y=conn_server('run_keep',mfilename,rt,filter,x,option,cols); end; return; end

USEDCT=true;
if strcmpi(option,'base'), Nx=x; x=eye(Nx); end
if ~isempty(cols), X=x; x=X(:,cols); end
if USEDCT % discrete cosine basis
    Nx=size(x,1);
    fy=fft(cat(1,x,flipud(x)),[],1);
    f=(0:size(fy,1)-1);
    f=min(f,size(fy,1)-f);
    switch(lower(option))
        case {'full','base'}
            idx=find(f<filter(1)*(rt*size(fy,1))|f>=filter(2)*(rt*size(fy,1)));
            %idx=idx(idx>1);
            fy(idx,:)=0;
            k=1; %2*size(fy,1)*(min(.5,filter(2)*rt)-max(0,filter(1)*rt))/max(1,size(fy,1)-numel(idx));
            y=real(ifft(fy,[],1))*k;
            y=y(1:Nx,:);
        case 'partial'
            idx=find(f>=filter(1)*(rt*size(x,1))&f<filter(2)*(rt*size(x,1)));
            %if ~any(idx==1), idx=[1,idx]; end
            y=fy(idx,:);
    end
else % discrete fourier basis
    fy=fft(x,[],1);
    f=(0:size(fy,1)-1);
    f=min(f,size(fy,1)-f);
    switch(lower(option))
        case 'full',
            idx=find(f<filter(1)*(rt*size(fy,1))|f>=filter(2)*(rt*size(fy,1)));
            %idx=idx(idx>1);
            fy(idx,:)=0;
            k=1; %2*size(fy,1)*(min(.5,filter(2)*rt)-max(0,filter(1)*rt))/max(1,size(fy,1)-numel(idx));
            y=real(ifft(fy,[],1))*k;
        case 'partial',
            idx=find(f>=filter(1)*(rt*size(x,1))&f<filter(2)*(rt*size(x,1)));
            %if ~any(idx==1), idx=[1,idx]; end
            y=fy(idx,:);
    end
end
if ~isempty(cols), X(:,cols)=y; y=X; end

