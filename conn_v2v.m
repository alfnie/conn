function varargout=conn_v2v(option,varargin)
global CONN_x;
persistent conv_maxstorage conv_kernels conv_data conv_results;
if isempty(conv_maxstorage)
    conv_maxstorage=1;
    conv_kernels=cell(1,conv_maxstorage);
    conv_data=cell(1,conv_maxstorage);
    conv_results=cell(1,conv_maxstorage);
end


switch(lower(option))
    case {'measures','empty','measurenames'}
        measures=struct(...
            'names',{{'group-PCA','group-ICA','group-MVPA','IntrinsicConnectivity','LocalCorrelation','GlobalCorrelation','RadialCorrelation','RadialSimilarity','ALFF','fALFF'}},... % connectome-level Multivariate Pattern Analysis; Integrated Local Correlation; Radial Correlation Contrast; Intrinsic Connectivity Contrast; Radial Similarity Contrast; 
            'alt_names',{{{'group PCA','PCA'},{'group ICA','ICA'},{'group MVPA','MVPA','MCOR'},{'IntrinsicConnectivityContrast','ICC'},{'IntegratedLocalCorrelation','ILC','LCOR'},{'GC','GCOR'},{'RadialCorrelationContrast','RCC'},{'RadialSimilarityContrast','RSC'},{'AmplitudOfLowFrequencyFluctuations'},{'FractionalAmplitudOfLowFrequencyFluctuations','FALFF'}}},...
            'measuretype',{{4,3,2,1,1,5,1,1,6,7}},...
            'localsupport',{{nan,nan,nan,0,25,inf,25,0,nan,nan}},...
            'global',{{nan,nan,nan,1,0,0,0,1,nan,nan}},...
            'deriv',{{nan,nan,nan,0,0,0,1,1,nan,nan}},...
            'norm',{{nan,nan,0,0,0,0,0,0,0,0}},... %{nan,nan,0,1,1,1,1,1,1,1}
            'mask',{{{},{},{},{},{},{},{},{},{},{},{},{}}},...
            'filename',{{'','','','','','','','','',''}},...
            'dimensions_out',{{40,40,4,1,1,1,3,3,1,1}},...
            'dimensions_in',{{64,64,64,64,64,64,64,64,64,64}});
        if strcmpi(option,'measures')
            varargout{1}=measures;
        elseif strcmp(option,'measurenames')
            varargout{1}=measures.names;
        else        
            measures_empty=measures;
            optionsnames=fieldnames(measures);
            for n1=1:numel(optionsnames),
                measures_empty.(optionsnames{n1})={};
            end
            varargout{1}=measures_empty;
        end
        
    case 'list_extended'
        options=varargin{1};
        optionsnames=fieldnames(options);
        measurename={};
        for index=1:numel(options.names)
            for ncomp=1:options.dimensions_out{index}
            %for ncomp=1-(options.dimensions_out{index}>1):options.dimensions_out{index}
                measurename{end+1}='';
                for n1=1:numel(optionsnames),
                    if ischar(options.(optionsnames{n1}){index})&&strcmp(optionsnames{n1},'names'), measurename{end}=[measurename{end},num2str(options.(optionsnames{n1}){index})];
                    elseif ischar(options.(optionsnames{n1}){index}),measurename{end}=[measurename{end},'#<',options.(optionsnames{n1}){index},'#>'];
                    elseif ~iscell(options.(optionsnames{n1}){index}), measurename{end}=[measurename{end},'#<',num2str(options.(optionsnames{n1}){index}),'#>'];
                    end
                end
                measurename{end}=[measurename{end},'#{',num2str(ncomp),'#}'];
            end
        end
        varargout{1}=measurename;

    case 'match_existing'
        measurename=varargin{1};
        idx1=strfind(measurename,'#<');
        if ~isempty(idx1), measurename=measurename(1:idx1(1)-1); end
        idx=find(strncmp([measurename,'#<'],CONN_x.vvAnalyses(CONN_x.vvAnalysis).measures,numel(measurename)+2));
        if ~isempty(idx),
            imeasure=idx;
            isnew=0;
        else
            imeasure=[];
            isnew=1;
        end
        varargout={imeasure,isnew};
        
    case 'match_extended'
        measurename=varargin{1};
        idx1=strfind(measurename,'#{');
        idx2=strfind(measurename,'#}');
        ncomp=str2num(measurename(idx1(1)+2:idx2(1)-1));
        measurename(idx1(1):idx2(1)+1)=[];
        idx=strmatch(measurename,CONN_x.vvAnalyses(CONN_x.vvAnalysis).measurenames,'exact');
        if ~isempty(idx),
            imeasure=idx(1);
            isnew=0;
        else
            imeasure=[];
            isnew=1;
        end
        varargout={imeasure,isnew,ncomp};
        
    case 'fieldtext'
        measurenames=varargin{1};
        if numel(varargin)<2||isempty(varargin{2}), nfield=1; else nfield=varargin{2}; end
        if iscell(measurenames)
            varargout={{}};
            for n1=1:numel(measurenames)
                varargout{1}{n1}=conn_v2v('fieldtext',measurenames{n1},varargin{2:end});
            end
        else
            varargout={''};
            for n=1:nfield-1
                k1=strfind(measurenames,'#<');
                k2=strfind(measurenames,'#>');
                if ~isempty(k1)&&~isempty(k2)&&k1(1)<k2(1), measurenames(1:k2(1)+1)=[]; end
            end
            k1=strfind(measurenames,'#<');
            k2=strfind(measurenames,'#>');
            if ~isempty(k1)&&~isempty(k2)&&k1(1)<k2(1), varargout={measurenames(k1(1)+2:k2(1)-1)}; end
        end
        
    case 'cleartext'
        measurenames=varargin{1};
        if iscell(measurenames)
            varargout={{}};
            for n1=1:numel(measurenames)
                varargout{1}{n1}=conn_v2v('cleartext',measurenames{n1});
            end
        else
            varargout={''};
            k1k2ok=0;
            while ~k1k2ok
                k1=strfind(measurenames,'#<');
                k2=strfind(measurenames,'#>');
                if ~isempty(k1)&&~isempty(k2)&&k1(1)<k2(1), measurenames(k1(1):k2(1)+1)=[];
                else k1k2ok=1; end
            end
            k1=strfind(measurenames,'#{');
            k2=strfind(measurenames,'#}');
            if ~isempty(k1)&&~isempty(k2)&&k1(1)<k2(1), measurenames(k1)='_'; measurenames([k2,k1+1,k2+1])=[]; end
            varargout{1}=measurenames;
        end
        
    case 'pcleartext'
        measurenames=varargin{1};
        if iscell(measurenames)
            varargout={{}};
            for n1=1:numel(measurenames)
                varargout{1}{n1}=conn_v2v('pcleartext',measurenames{n1});
            end
        else
            varargout={''};
            k1k2ok=0;
            while ~k1k2ok
                k1=strfind(measurenames,'#<');
                k2=strfind(measurenames,'#>');
                if ~isempty(k1)&&~isempty(k2)&&k1(1)<k2(1), measurenames(k1)='_'; measurenames([k2,k1+1,k2+1])=[]; 
                %if ~isempty(k1)&&~isempty(k2)&&k1(1)<k2(1), measurenames(k1(1):k2(1)+1)=[];
                else k1k2ok=1; end
            end
            k1=strfind(measurenames,'#{');
            k2=strfind(measurenames,'#}');
            if ~isempty(k1)&&~isempty(k2)&&k1(1)<k2(1), measurenames(k1)='_'; measurenames([k2,k1+1,k2+1])=[]; end
            measurenames=regexprep(measurenames,{'NaN','_+'},{'','_'});
            varargout{1}=measurenames;
        end
        
    case 'match_measures'
        options=varargin{1};
        index=varargin{2};
        optionsnames=fieldnames(options);
        measurename=[];
        for n1=1:numel(optionsnames),
            %assignin('caller',optionsnames{n1},options.(optionsnames{n1}){index}); 
            if ischar(options.(optionsnames{n1}){index})&&strcmp(optionsnames{n1},'names'), measurename=[measurename,options.(optionsnames{n1}){index}];
            elseif ischar(options.(optionsnames{n1}){index}), measurename=[measurename,'#<',options.(optionsnames{n1}){index},'#>'];
            elseif ~iscell(options.(optionsnames{n1}){index}), measurename=[measurename,'#<',num2str(options.(optionsnames{n1}){index}),'#>'];
            end
        end
        if numel(varargin)<3, flag='-'; else flag=varargin{3}; end
        if ~isfield(CONN_x.vvAnalyses(CONN_x.vvAnalysis),'measurenames'),
            CONN_x.vvAnalyses(CONN_x.vvAnalysis).measurenames={};
        end
        idx=strmatch(measurename,CONN_x.vvAnalyses(CONN_x.vvAnalysis).measurenames,'exact');
        if ~isempty(idx),
            imeasure=idx(1);
            isnew=0;
        else
            imeasure=length(CONN_x.vvAnalyses(CONN_x.vvAnalysis).measurenames)+1;
            if flag=='+',
                CONN_x.vvAnalyses(CONN_x.vvAnalysis).measurenames{imeasure}=measurename;
                filemeasurenames=fullfile(CONN_x.folders.firstlevel_vv,CONN_x.vvAnalyses(CONN_x.vvAnalysis).name,'_list_measures.mat');
                filemeasurenames=conn_projectmanager('projectfile',filemeasurenames,CONN_x.pobj,'.mat');
                measurenames=CONN_x.vvAnalyses(CONN_x.vvAnalysis).measurenames;
                save(filemeasurenames,'measurenames');
            end
            isnew=1;
        end
        varargout={imeasure,isnew};
        
    case 'compute_start'
        options=varargin{1};
        index=varargin{2};
        optionsnames=fieldnames(options);
        for n1=1:numel(optionsnames),
            options.(optionsnames{n1})=options.(optionsnames{n1}){index}; 
        end
        options.mat=varargin{3};
        options.issurface=varargin{4};
        if ~options.issurface, options.localsupport=options.localsupport./sqrt(sum(abs(options.mat(1:3,1:3)).^2,1)); end
        if options.measuretype==1
            if options.deriv==1, options.dimensions_out=2+~options.issurface; else options.dimensions_out=1; end
        end
        if options.dimensions_out>1, options.m=repmat({0},[1,options.dimensions_out]);
        else options.m=0;
        end
        varargout{1}=options;
        
    case 'compute_step'
        options=varargin{1};
        x=varargin{2}; % eigenvariate
        d=varargin{3}; % eigenvalue
        dall=varargin{4}; % eigenvalues all
        nx=varargin{5}; % # of voxels

        if all(options.localsupport==0), y=x;
        elseif all(isinf(options.localsupport)), y=sum(x(:))/nx;
        else
            ok=0;
            n=numel(conv_kernels);
            for n1=1:n, % checks in storage to avoid duplication of (time-consuming) convolution operations
                if isequal(conv_kernels{n1},options.localsupport) && isequal(conv_data{n1},x), y=conv_results{n1}; ok=1; break; end
            end
            if ~ok
                y=conn_conv(x,options.localsupport,[],options.issurface);
                nn=1+rem(n,conv_maxstorage);
                conv_kernels{nn}=options.localsupport;
                conv_data{nn}=x;
                conv_results{nn}=y;
            end
        end
%         if options.global==0; k=1;
%         elseif options.global==1, k=d;
%         else k=d.^options.global; 
%         end
        k=d/nx; %numel(x);
        kall=dall/nx; %numel(x);
        idxnull=x==0;
%         if options.measuretype==5
%             alphaKC=.9;
%             options.m=options.m+x.*y*(1/(1-alphaKC*d/max(dall))/k);
        if options.deriv==0
            if options.global==0,
                options.m=options.m+x.*y;
            else
                options.m=options.m+k*abs(y).^2;
            end
        elseif options.deriv==1
            if options.global==0,
                [dx{1:2+~options.issurface}]=conn_deriv(y,1,options.issurface);
                for n=1:2+~options.issurface
                    options.m{n}=options.m{n}+x.*dx{n}/2;
                end
            else
                [dx{1:2+~options.issurface}]=conn_deriv(y,1,options.issurface);
                for n=1:2+~options.issurface
                    options.m{n}=options.m{n}+k*abs(dx{n}/2).^2;
                end
            end
        elseif options.deriv==2
            dx=conn_deriv(y,2);
            if options.global==0,
                options.m=options.m+x.*dx;
            else
                options.m=options.m+k*abs(dx).^2;
            end
        end
        if iscell(options.m)
            for n1=1:numel(options.m)
                options.m{n1}(idxnull)=0;
            end
        else
            options.m(idxnull)=0;
        end
        varargout{1}=options;
        
    case 'compute_end'
        options=varargin{1};
        if options.global>0,
            if iscell(options.m)
                for n1=1:numel(options.m)
                    options.m{n1}=sqrt(options.m{n1});
                end
            else
                options.m=sqrt(options.m);
            end
        end
        if ~isfield(options,'norm')||options.norm
            if ~iscell(options.m)
                idx=find(~isnan(options.m)&options.m~=0);
                [usort,nill,idxsort]=unique(options.m(idx));options.m(idx)=spm_invNcdf(idxsort/(1+numel(usort)));
                %[nill,idxsort]=sort(options.m(idx));idxsort(idxsort)=spm_invNcdf((1:numel(idxsort))/(numel(idxsort)+1));options.m(idx)=idxsort;
                %options.m(idx)=(options.m(idx)-mean(options.m(idx)))/std(options.m(idx));
            else
                for n1=1:numel(options.m)
                    idx=find(~isnan(options.m{n1})&options.m{n1}~=0);
                    [usort,nill,idxsort]=unique(options.m{n1}(idx));options.m{n1}(idx)=spm_invNcdf(idxsort/(1+numel(usort)));
                    %[nill,idxsort]=sort(options.m{n1}(idx));idxsort(idxsort)=spm_invNcdf((1:numel(idxsort))/(numel(idxsort)+1));options.m{n1}(idx)=idxsort;
                    %options.m{n1}(idx)=(options.m{n1}(idx)-mean(options.m{n1}(idx)))/std(options.m{n1}(idx));
                end
            end
        end
        varargout{1}=options.m;
%         if options.deriv==1
%             for n1=1:3,
%                 varargout{1}{n1}=0;
%                 for n2=1:3
%                     varargout{1}{n1}=varargout{1}{n1}+options.mat(n1,n2)*options.m{n2};
%                 end
%             end
%             varargout{1}=sqrt(varargout{1}{1}.^2+varargout{1}{2}.^2+varargout{1}{3}.^2);
%         else
%             varargout{1}=options.m;
%         end
        
    case 'compute_slice'
        options=varargin{1};
        index=varargin{2};
        optionsnames=fieldnames(options);
        for n1=1:numel(optionsnames),
            options.(optionsnames{n1})=options.(optionsnames{n1}){index}; 
        end
        vols=varargin{3};
        options.mat=vols.matdim.mat;
        options.localsupport=options.localsupport./sqrt(sum(abs(options.mat(1:3,1:3)).^2,1));
        d=varargin{4};
        slice=varargin{5};

        %vols=vols(1:min(vols.size.Nt,options.dimensions_in));
        Nt=min(vols.size.Nt,options.dimensions_in);
        d=d(1:min(numel(d),options.dimensions_in));

        ok=0;
        n=numel(conv_kernels);
        for n1=1:n, % checks in storage to avoid duplication of (time-consuming) convolution operations
            if isequal(conv_kernels{n1},vols), x=conv_data{n1}; ok=1; break; end
        end
        if ~ok
            x=zeros([vols.size.Nt,vols.matdim.dim]);
            for n1=1:vols.size.Nt, x(n1,:)=reshape(conn_get_time(vols,n1),1,[]); end
            nn=1+rem(n,conv_maxstorage);
            conv_kernels{nn}=vols;
            conv_data{nn}=x;
            conv_results{nn}=[];
        end
        x=x(1:Nt,:,:,:);
        %k=d.^options.global;
        Nv=prod(vols.matdim.dim);
        k=d/Nv;
        if all(options.localsupport==0), 
            n=double(options.deriv>0);
            y=x(:,:,:,max(1,min(size(x,4),slice+(-n:n))));
        elseif all(isinf(options.localsupport)), 
            %y=bsxfun(@plus,mean(mean(mean(x,2),3),4),zeros(1,size(x,2),size(x,3)));
            y=repmat(mean(mean(mean(x,2),3),4),[1,size(x,2),size(x,3)]);
            if options.deriv>0, y=repmat(y,[1,1,1,3]); end
        else
            [y,h]=conn_conv(x,[0,options.localsupport(1:3)],0);
            n=round((numel(h{4})-1)/2)+(options.deriv>0);
            y=conn_conv(x(:,:,:,max(1,min(size(x,4),slice+(-n:n)))),[0,options.localsupport(1:2),0]);
            if options.deriv==0
                y=sum(conn_bsxfun(@times,y,shiftdim(h{4},-3)),4);
            else
                y=cat(4, sum(conn_bsxfun(@times,y(:,:,:,1:end-2),shiftdim(h{4},-3)),4),...
                    sum(conn_bsxfun(@times,y(:,:,:,2:end-1),shiftdim(h{4},-3)),4),...
                    sum(conn_bsxfun(@times,y(:,:,:,3:end),shiftdim(h{4},-3)),4) );
            end
        end
        x=x(:,:,:,slice);
        idxnull=shiftdim(any(x==0,1),1);
%         if options.measuretype==5
%             v=(y(:,1)./sqrt(d))*Nv;
%             k=sqrt(max(0,Nv-sum(abs(v).^2)));
%             ve=[v;k];
%             de=[d;0];
%             Rpos=(ve*ve'+diag(de))/2;
%             [vpos,nill]=svd(Rpos);
%             vpos=vpos(:,1);
%             if vpos(end)<0, vpos=-vpos; end
%             m=vpos(end)/k+shiftdim(sum(conn_bsxfun(@times,(vpos(1:end-1)-v*vpos(end)/k)./sqrt(d), x), 1),1);
% %             EC = k*Q*w - 1
% %             v = Q?*1
% %             w = v./(1-alpha/2*d)
% %             k = 1/(1-alpha/2*v?*w)
% %             v = (y(:,1)./sqrt(d))*prod(vols.matdim.dim);
% %             alphaKC=10*.99*2/sum((v.^2+d)/2);
% %             w = conn_bsxfun(@rdivide, v, 1-alphaKC/2*d);
% %             k = 1;%/(1-alphaKC/2*(v'*(w-v)+prod(vols.matdim.dim)));
% %             m=shiftdim(sum(conn_bsxfun(@times,k*(w-v)./sqrt(d), x), 1),1);
%             varargout{1}=m;
        if options.deriv==0
            if options.global==0,
                m=shiftdim(sum(x.*y,1),1);
            else
                m=sqrt(shiftdim(sum(conn_bsxfun(@times,k,abs(y).^2),1),1));
            end
            varargout{1}=m;
        elseif options.deriv==1
            if options.global==0,
                dx=y(:,[2:end,end],:,2)-y(:,[1,1:end-1],:,2);
                m{1}=shiftdim(sum(x.*dx,1),1)/2;
                dx=y(:,:,[2:end,end],2)-y(:,:,[1,1:end-1],2);
                m{2}=shiftdim(sum(x.*dx,1),1)/2;
                dx=y(:,:,:,3)-y(:,:,:,1);
                m{3}=shiftdim(sum(x.*dx,1),1)/2;
            else
                dx=y(:,[2:end,end],:,2)-y(:,[1,1:end-1],:,2);
                m{1}=sqrt(shiftdim(sum(conn_bsxfun(@times,k,abs(dx).^2),1),1)/2);
                dx=y(:,:,[2:end,end],2)-y(:,:,[1,1:end-1],2);
                m{2}=sqrt(shiftdim(sum(conn_bsxfun(@times,k,abs(dx).^2),1),1)/2);
                dx=y(:,:,:,3)-y(:,:,:,1);
                m{3}=sqrt(shiftdim(sum(conn_bsxfun(@times,k,abs(dx).^2),1),1)/2);
            end
            for n1=1:3,
                varargout{1}{n1}=0;
                for n2=1:3
                    varargout{1}{n1}=varargout{1}{n1}+options.mat(n1,n2)*m{n2};
                end
            end
            varargout{1}=sqrt(varargout{1}{1}.^2+varargout{1}{2}.^2+varargout{1}{3}.^2);
        elseif options.deriv==2
            dx=y(:,[2:end,end],:,2)+y(:,[1,1:end-1],:,2);
            dx=dx+y(:,:,[2:end,end],2)+y(:,:,[1,1:end-1],2);
            dx=dx+y(:,:,:,3)+y(:,:,:,1);
            if options.global==0,
                m=shiftdim(sum(x.*(y(:,:,:,2)-dx/6),1),1);
            else
                m=sqrt(shiftdim(sum(conn_bsxfun(@times,k,abs(y(:,:,:,2)-dx/6).^2),1),1));
            end                
            varargout{1}=m;
        end
        %if options.measuretype==5, varargout{1}=(1+varargout{1})/2; end
        if iscell(varargout{1})
            for n1=1:numel(varargout{1})
                varargout{1}{n1}(idxnull)=0;
            end
        else
            varargout{1}(idxnull)=0;
        end
%         if ~iscell(varargout{1})
%             idx=find(~isnan(varargout{1})&varargout{1}~=0);
%             varargout{1}(idx)=(varargout{1}(idx)-mean(varargout{1}(idx)))/std(varargout{1}(idx));
%         end
        if ~isfield(options,'norm')||options.norm
            if ~iscell(varargout{1})
                idx=find(~isnan(varargout{1})&varargout{1}~=0);
                [usort,nill,idxsort]=unique(varargout{1}(idx));varargout{1}(idx)=spm_invNcdf(idxsort/(1+numel(usort)));
                %[nill,idxsort]=sort(varargout{1}(idx));idxsort(idxsort)=spm_invNcdf((1:numel(idxsort))/(numel(idxsort)+1));varargout{1}(idx)=idxsort;
                %varargout{1}(idx)=(varargout{1}(idx)-mean(varargout{1}(idx)))/std(varargout{1}(idx));
            else
                for n1=1:numel(varargout{1})
                    idx=find(~isnan(varargout{1}{n1})&varargout{1}{n1}~=0);
                    [usort,nill,idxsort]=unique(varargout{1}{n1}(idx));varargout{1}{n1}(idx)=spm_invNcdf(idxsort/(1+numel(usort)));
                    %[nill,idxsort]=sort(varargout{1}{n1}(idx));idxsort(idxsort)=spm_invNcdf((1:numel(idxsort))/(numel(idxsort)+1));varargout{1}{n1}(idx)=idxsort;
                    %varargout{1}{n1}(idx)=(varargout{1}{n1}(idx)-mean(varargout{1}{n1}(idx)))/std(varargout{1}{n1}(idx));
                end
            end
        end
        
        
    case 'compute_ax'
        x=varargin{1};
        filepath=varargin{2};
        nsubjects=varargin{3};
        nconditions=varargin{4};
        ndims=varargin{5};
        y=zeros(size(x));
        icondition=[];isnewcondition=[];for ncondition=nconditions(:)',[icondition(ncondition),isnewcondition(ncondition)]=conn_conditionnames(CONN_x.Setup.conditions.names{ncondition}); end
        if any(isnewcondition), error(['Some conditions have not been processed yet. Re-run previous step']); end
        for nsub=1:nsubjects,
            for ncond=1:numel(nconditions)
                ncondition=nconditions(ncond);
                fprintf('.')
                %             filename_B1=fullfile(filepath,['vvPC_Subject',num2str(nsub,'%03d'),'_Condition',num2str(ncondition,'%03d'),'.nii']);
                %             Y1=spm_vol(filename_B1);
                filename_B1=fullfile(filepath,['vvPC_Subject',num2str(nsub,'%03d'),'_Condition',num2str(icondition(ncondition),'%03d'),'.mat']);
                Y1=conn_vol(filename_B1);
                if nsub==1&&ncond==1
                    %                 [gridx,gridy,gridz]=ndgrid(1:Y1(1).dim(1),1:Y1(1).dim(2),1:Y1(1).dim(3));xyz=[gridx(:),gridy(:),gridz(:),ones(numel(gridx),1)]';
                    matdim=Y1.matdim.dim;
                else
                    if ~isequal(matdim,Y1.matdim.dim), error('ERROR: unequal volume sizes across subjects when computing group-level measures'); end
                end
                z=0;
                for slice=1:Y1.size.Ns,
                    [y1,idx]=conn_get_slice(Y1,slice);
                    idx=(slice-1)*Y1.matdim.dim(1)*Y1.matdim.dim(2)+idx;
                    z=z+y1(1:min(size(y1,1),ndims),:)*x(idx);
                end
                for slice=1:Y1.size.Ns,
                    [y1,idx]=conn_get_slice(Y1,slice);
                    idx=(slice-1)*Y1.matdim.dim(1)*Y1.matdim.dim(2)+idx;
                    y(idx)=y(idx)+y1(1:min(size(y1,1),ndims),:)'*z;
                end
                figure(2)
                temp=reshape(x,Y1.matdim.dim);subplot(121);imagesc(temp(:,:,30));
                temp=reshape(y,Y1.matdim.dim);subplot(122);imagesc(temp(:,:,30));
                drawnow;
%                 y=zeros(size(x));
%                 [y1,idx]=conn_get_volume(Y1);
%                 y(idx)=y1'*(y1*x(idx));
%                 for ndim=1:ndims
%                     y1=conn_get_time(Y1,ndim);
%                     %y1=spm_get_data(Y1(ndim),xyz);
%                     y=y+y1(:)*(y1(:)'*x);
%                 end
            end
        end
        fprintf('\n')
        varargout{1}=y;
end
