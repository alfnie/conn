function varargout = conn_sliceintensitycorrection(fname,varargin)
% conn_sliceintensitycorrection(fname)
% corrects global differences in intensity across slices (e.g. caused by a scanner transients in a sparse acquisition sequence)
%
% e.g. 
%   conn_sliceintensitycorrection rest.nii
%
% additional options: (use conn_sliceintensitycorrection(filename, option1_name, option1_value, option2_name, option2_value, ...)
%   miniter         : minimum number of iterations [ 100 ]
%   maxiter         : maximum number of iterations [ 1000 ]
%   maxdiff         : break out rule (mean change in optim params below this value) [1e-6]
%   lambda          : regularization factor [ 0.001 ]
%   eps             : iteration step size [ 0.01 ]
%   acc             : iteration acceleration factor [ 0.9 ]
%   maxdk           : maximum scaling (maximum scale = 1+options.maxdk, minimum scale = 1-options.maxdk) [ 0.5 ]
%   domean          : 0: cost function uses individual images (slower); 1: cost function uses average image across timepoints/scans (faster); [ 1 ]
%   doplot          : iterative plots [ false ]
%
% algorithm notes: 
%   accelerated gradient descent with cost function: sum_(x,y,z,t){ | h(z)*s(x,y,z,t) - h(z+1)*s(x,y,z+1,t) | / ( h(z) + h(z+1) ) } + lambda sum_(z){ (h(z)-1)^2 }
%   with s(x,y,z,t) = BOLD signal 4D data
%        h(z) = estimated rescaling factors
%   initial condition: h(z)=1
%

if ~iscell(fname), fname={fname}; end
if any(conn_server('util_isremotefile',fname)), [varargout{1:nargout}]=conn_server('run',mfilename,conn_server('util_localfile',fname),varargin{:}); return; end

options.miniter=1e2; % min iterations
options.maxiter=1e3; % max iterations
options.maxdiff=1e-6;% break out difference in optim params
options.lambda=.001; % regularization factor
options.eps=.01;     % step size
options.acc=.9;      % acceleration factors
options.maxdk=.5;    % maximum scaling (1+-options.maxdk)
options.domean=true;
options.doplot=false;
for n1=1:2:numel(varargin)-1
    assert(isfield(options,lower(varargin{n1})),'unrecognized option %s',varargin{n1});
    options.(lower(varargin{n1}))=varargin{n1+1};
end

fnameout={};
for nfname=1:numel(fname)
    if options.domean
        vol=spm_vol(fname{nfname});
        x=0; mn=inf; mx=-inf; 
        for nt=1:numel(vol), 
            tx=spm_read_vols(vol(nt)); 
            x=x+tx; 
            mn=min(mn,tx);
            mx=max(mx,tx);
        end
        x=x/numel(vol);
    else
        [x,vol]=conn_vol_read(fname{nfname});
    end

    meanx=mean(x(~isnan(x)));
    %dx=max(eps,min(diff(unique(x)))); % possible discretization step
    %x=x/mean(x(~isnan(x)));
    f=zeros(1,size(x,3));

    dftotal=0;
    for nrepeat=1:options.maxiter
        h=1+f;
        %h=exp(f);
        %h=h/mean(h);
        g=x.*repmat(shiftdim(h,-1),[size(x,1),size(x,2),1,size(x,4)]);


        h1=h+h([2:end,end]);
        h2=h+h([1,1:end-1]);
        e1=(g-g(:,:,[2:end,end],:))./repmat(shiftdim(h1,-1),[size(x,1),size(x,2),1,size(x,4)]);
        e2=(g-g(:,:,[1,1:end-1],:))./repmat(shiftdim(h2,-1),[size(x,1),size(x,2),1,size(x,4)]);
        %m1=e1==0; e1(m1)=dx*(rand(nnz(m1),1)-.5);
        %m2=e2==0; e2(m2)=dx*(rand(nnz(m2),1)-.5);
        df= -shiftdim(mean(mean(mean( tanh(100/meanx*e1).*(x - e1) ,1),2),4),1)./max(eps,h1) ...
            -shiftdim(mean(mean(mean( tanh(100/meanx*e2).*(x - e2) ,1),2),4),1)./max(eps,h2);

        madf=max(abs(df));
        if madf>meanx, df=df/madf*meanx; end
        df2=-options.lambda*f+options.eps*df/meanx;
        dftotal=options.acc*dftotal+df2;
        f=f+dftotal;
        %%f=f-mean(f);
        if max(abs(f))>options.maxdk, f=f/max(abs(f))*options.maxdk; end
        if mean(abs(dftotal))<options.maxdiff&nrepeat>options.miniter, break; end

        if options.doplot
            disp(madf/meanx);
            subplot(221);bar(1:numel(f),1+f); set(gca,'ylim',[1-options.maxdk, 1+options.maxdk]); ylabel('re-scaling factor'); xlabel('slice number');
            subplot(223);bar(1:numel(f),df2); ylabel('parameter change'); xlabel('slice number');
            subplot(122);k=ceil(size(g,2)/2); imagesc([flipud(squeeze(x(:,k,:,1))') flipud(squeeze(g(:,k,:,1))')]); axis equal tight;
            drawnow
        elseif ~rem(nrepeat,100) fprintf('.');
        end
    end

    if size(fname{nfname},1)>1, fnameout{nfname}=conn_prepend('i',fname{nfname}(1,:));
    else fnameout{nfname}=conn_prepend('i',fname{nfname});
    end
    fnameout{nfname}=regexprep(fnameout{nfname},',\d+\s*$','');
    if options.domean
        ftype=vol(1).dt(1);
        N=numel(vol);
        mn=mn.*repmat(shiftdim(h,-1),[size(mn,1),size(mn,2),1]);
        mx=mx.*repmat(shiftdim(h,-1),[size(mx,1),size(mx,2),1]);
        pinfo = spm_write_vol_rescale(ftype,min(mn(:)),max(mx(:)));
        a=struct('fname',fnameout{nfname},'mat',vol(1).mat,'dim',vol(1).dim,'n',[1,1],'pinfo',pinfo,'dt',[ftype spm_platform('bigend')]);
        spm_unlink(fnameout{nfname});
        a=repmat(a,1,N); for n=1:N, a(n).n=[n,1]; end
        a=spm_create_vol(a);
        for nt=1:N, 
            tx=spm_read_vols(vol(nt)); 
            tg=tx.*repmat(shiftdim(h,-1),[size(tx,1),size(tx,2),1]);
            a(nt)=spm_write_vol(a(nt),tg); 
        end
    else
        conn_vol_write(conn_prepend('i',fname{nfname}), g, vol);
    end
    try, conn_disp('fprintf','Slice Intensity Correction output file %s (%d iterations, scaling %s)\n',fnameout{nfname}, nrepeat, mat2str(h,4)); end
end

varargout={fnameout,h,x,g};

% conn_sliceintesitycorrection(fname, sliceorder)
% with sliceorder one of the following options:
%                1) leave empty (default) to read this information from a sidecar .json file associated with the input NIFTI file (e.g. a sidecar rest.json file associated with a rest.nii file)
%                2) any of the following keywords: 'ascending','descending','interleaved (Siemens)','interleaved (Philips)','interleaved (bottom-up)','interleaved (top-down)','interleaved (middle-top)'
%                3) a vector of slice indexes from z=1 -first slice in image- to z=? -last slice- in the order they were acquired
%                4) a vector of acquisition times of each slice in milliseconds (e.g. for multiband sequences)
%
if 0 % todo: consider slice acquisition order?
    slicetime=[];
    options={'ascending','descending','interleaved (Siemens)','interleaved (Philips)','interleaved (bottom-up)','interleaved (top-down)','interleaved (middle-top)'};
    vol=spm_vol(fname{nfname});
    nslices=vol(1).dim(3);
    if isempty(slicetime)
        % BIDS json file
        str=conn_jsonread(fname{nfname},'SliceTiming');
        if isempty(str), % SliceTiming information not found in json file
            sliceorder=1:nslices;
            sliceorder=inputdlg(['Slice order? (enter slice indexes from z=1 -first slice in image- to z=',num2str(nslices),' -last slice- in the order they were acquired). Alternatively enter acquisition time of each slice in milliseconds (e.g. for multiband sequences)'],'conn_sliceintensitycorrection',1,{sprintf('%d ',sliceorder)});
            if isempty(sliceorder), return;
            else sliceorder=str2num(regexprep(sliceorder{1},'[a-zA-Z]+',num2str(nslices)));
            end
            if (numel(unique(sliceorder))~=nslices||max(sliceorder)~=nslices||min(sliceorder)~=1) && (numel(sliceorder)~=nslices||any(sliceorder<0)), return; end
        else sliceorder=1000*reshape(str,1,[]); % (ms)
        end
    elseif ischar(slicetime)
        switch(slicetime)
            case 'ascending', sliceorder=1:nslices;        % ascending
            case 'descending', sliceorder=nslices:-1:1;     % descending
            case 'interleaved (Siemens)', sliceorder=[fliplr(nslices:-2:1) fliplr(nslices-1:-2:1)]; % interleaved (Siemens)
            case 'interleaved (Philips)', sliceorder=cell2mat(arrayfun(@(n)n:round(sqrt(nslices)):nslices,1:round(sqrt(nslices)),'uni',0)); % interleaved (Philips)
            case 'interleaved (bottom-up)', sliceorder=[1:2:nslices 2:2:nslices]; % interleaved (bottom-up)
            case 'interleaved (top-down)', sliceorder=[nslices:-2:1, nslices-1:-2:1]; % interleaved (top-down)
            case 'interleaved (middle-top)', sliceorder=round((nslices-(1:nslices))/2 + (rem((nslices-(1:nslices)),2) * (nslices - 1)/2)) + 1; % interleaved (middle-top)
            otherwise, error('unrecognized slice order %s (valid options are %s)',slicetime,sprintf('%s; ',options{:}));
        end
    else sliceorder=slicetime; 
    end
    if (numel(unique(sliceorder))~=nslices||max(sliceorder)~=nslices||min(sliceorder)~=1),  % slice timing (ms)
        fprintf('slice timing information %s\n',mat2str(sliceorder));
        [nill,nill,slicetimerank]=unique(sliceorder);
    elseif 1, % note: convert to slice-timing (ms) syntax for SPM
        fprintf('slice order information %s\n',mat2str(sliceorder));
        [nill,slicetimerank]=sort(sliceorder);
    end
    fprintf('slice ranked time%s\n',mat2str(slicetimerank));
end
end

function pinfo = spm_write_vol_rescale(dt,mn,mx)
    % adapted from spm_write_vol
    pinfo = [1;0;0];
    dt           = dt(1);
    s            = find(dt == [2 4 8 256 512 768]);
    if ~isempty(s)
        dmnmx        = [0 -2^15 -2^31 -2^7 0 0 ; 2^8-1 2^15-1 2^31-1 2^7-1 2^16-1 2^32-1];
        dmnmx        = dmnmx(:,s);
        if isempty(mx), mx = 0; end
        if isempty(mn), mn = 0; end
        if mx ~= mn
            if dmnmx(1) < 0
                pinfo(1) = max(mx/dmnmx(2),mn/dmnmx(1));
            else
                pinfo(1) = mx/dmnmx(2);
            end
            pinfo(2) = 0;
        else
            pinfo(1,1) = mx/dmnmx(2);
            pinfo(2,1) = 0;
        end
    end
end

