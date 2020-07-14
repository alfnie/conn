function varargout = conn_wilkinson(option, varargin)
% (internal) Wilkinson notation

global CONN_x;

switch(option)
    case 'describe' % 'describe', names [, response, doextend]
        if numel(varargin)<1||isempty(varargin{1}), names={}; else names=varargin{1}; end
        if numel(varargin)<2||isempty(varargin{2}), response='?'; else response=varargin{2}; end
        if numel(varargin)<3||isempty(varargin{3}), cb=[]; else cb=varargin{3}; end
        if numel(varargin)<4||isempty(varargin{4}), cw=[]; else cw=varargin{4}; end
        if numel(varargin)<5||isempty(varargin{5}), other={}; else other=varargin{5}; end
        if numel(varargin)<6||isempty(varargin{6}), doextend=true; else doextend=varargin{6}; end
        
        descrip='';
        str = conn_strjoinstr(names,' + ');
        str = regexprep(str, '\<AllSubjects\>','1');
        str = regexprep(str, '\+?\s+$','');
        %if isempty(regexp(str,'\<1\>')), str=['- 1 + ',str]; end
        if iscell(response), response=conn_strjoinstr(response,' + '); end
        str = ['y ~ ',str];
        if doextend, 
            if isempty(cw)||isequal(cw,1), cw=''; else cw=sprintf('(%s)',regexprep(mat2str(rats(cw)),{'\s+','''','^\s+|\s+$'},{' ','',''})); end
            if isempty(cb)||isequal(cb,1), cb=''; else cb=sprintf('(%s)',regexprep(mat2str(rats(cb)),{'\s+','''','^\s+|\s+$'},{' ','',''})); end
            str={  sprintf('Between-subjects model specification : %s  %s',str, cb),...
                   sprintf('Within-subjects model specification: y ~ %s  %s',response, cw)}; 
            if ~isempty(other), str{end+1}=sprintf('Results stored in %s',regexprep(other,'^.*results.secondlevel.','')); end
        end
        varargout={str};
        
    case 'suggest_within' % 'suggest', names
        if numel(varargin)<1||isempty(varargin{1}), names=CONN_x.Setup.conditions.names(1:end-1); 
        elseif isnumeric(varargin{1}), names=CONN_x.Setup.conditions.names(varargin{1}); 
        else names=varargin{1}; 
        end
        str={}; descrip={}; neffects={}; ceffects={}; priority=[];
        isorig=find(cellfun('length',regexp(names,'x Time\d+$|x Temporal'))==0);
        istime=find(cellfun('length',regexp(names,'x Time\d+$'))>0);
        Tmatch={}; Fmatch={};
        for n1=1:numel(names)
            tnames=names{n1};
            descrip{end+1}=sprintf('functional connectivity during %s',tnames); 
            neffects{end+1}=names(n1);
            ceffects{end+1}=1;
            str{end+1}='connectivity';
            priority(end+1)=1;
            if ~isempty(regexp(tnames,'.?x Temporal Variability$')),
                descrip{end}=regexprep(descrip{end},'^functional connectivity during (.*?)\s*x Temporal Variability$','temporal variability in functional connectivity strength during $1');
                str{end}='variability';
            elseif ~isempty(regexp(tnames,'.?x Temporal Average$')),
                descrip{end}=regexprep(descrip{end},'^functional connectivity during (.*?)\s*x Temporal Average$','temporal average in functional connectivity strength during $1');
                str{end}='connectivity';
            elseif ~isempty(regexp(tnames,'.?x Time\d+$')),
                descrip{end}=regexprep(descrip{end},'^functional connectivity during (.*?)\s*x Time(\d+)$','functional connectivity during $1 sliding-window #$2');
                tmatch=regexp(tnames,'^(.*?)x Time(\d+)$','tokens');
                if numel(tmatch)==1&&numel(tmatch{1})==2, Tmatch=[Tmatch; [tmatch{1},{tnames}]]; end
            elseif ~isempty(regexp(tnames,'\d+$'))
                tmatch=regexp(tnames,'^(.*?)(\d+)$','tokens');
                if numel(tmatch)==1&&numel(tmatch{1})==2, Fmatch=[Fmatch; [tmatch{1},{n1}]]; end
            end
        end
        if ~isempty(Fmatch)
            [ucond,nill,icond]=unique(Fmatch(:,1));
            for n0=1:numel(ucond)
                Fidx=find(icond==n0);
                [nill,idx]=sort(str2double(Fmatch(Fidx,2)));
                if numel(idx)>1
                    Ieffects=cell2mat(Fmatch(Fidx(idx),3)');
                    if numel(idx)<=4
                        for n1=1:numel(idx),
                            for n2=[1:n1-1,n1+1:numel(idx)]
                                ieffects=Ieffects([n1,n2]);
                                descrip{end+1}=sprintf('functional connectivity change between %s and %s',names{ieffects(1)},names{ieffects(2)});
                                neffects{end+1}=names(ieffects);
                                ceffects{end+1}=[-1 1];
                                str{end+1}='connectivity-change';
                                priority(end+1)=100+n1;
                                [nill,iidx]=sort(ieffects); neffects{end}=neffects{end}(iidx); ceffects{end}=ceffects{end}(:,iidx);
                            end
                        end
                    end
                    if numel(idx)>2
                        ieffects=Ieffects;
                        descrip{end+1}=sprintf('functional connectivity changes between %s and %s',conn_strjoinstr(names(ieffects(1:end-1)),', '),names{ieffects(end)});
                        neffects{end+1}=names(ieffects);
                        ceffects{end+1}=diff(eye(numel(idx)));
                        str{end+1}='connectivity-change';
                        priority(end+1)=100+numel(idx)+1;
                    end
                    isorig=setdiff(isorig,Ieffects);
                end
            end
        end
        if numel(isorig)>1
            if numel(isorig)<=4
                for n1=1:numel(isorig),
                    for n2=[1:n1-1,n1+1:numel(isorig)]
                        ieffects=isorig([n1,n2]);
                        descrip{end+1}=sprintf('functional connectivity change between %s and %s',names{ieffects(1)},names{ieffects(2)});
                        neffects{end+1}=names(ieffects);
                        ceffects{end+1}=[-1 1];
                        str{end+1}='connectivity-change';
                        priority(end+1)=200+n1;
                        [nill,iidx]=sort(ieffects); neffects{end}=neffects{end}(iidx); ceffects{end}=ceffects{end}(:,iidx);
                    end
                end
            end
            if numel(isorig)>2
                descrip{end+1}=sprintf('functional connectivity changes between %s and %s',conn_strjoinstr(names(isorig(1:end-1)),', '),names{isorig(end)});
                neffects{end+1}=names(isorig);
                ceffects{end+1}=diff(eye(numel(isorig)));
                str{end+1}='connectivity-change';
                priority(end+1)=200+numel(isorig)+1;
            end
        end
        if ~isempty(Tmatch)
            [ucond,nill,icond]=unique(Tmatch(:,1));
            for n1=1:numel(ucond)
                [nill,idx]=sort(str2double(Tmatch(icond==n1,2)));
                idx(idx)=1:numel(idx);
                if numel(idx)>1
                    descrip{end+1}=sprintf('linear temporal trends in functional connectivity strength during %s',ucond{n1});
                    neffects{end+1}=Tmatch(icond==n1,3)';
                    temp=reshape(-1/2+(idx-1)/(numel(idx)-1),1,[]);
                    ceffects{end+1}=temp/sum(temp.^2);
                    str{end+1}='trend-slope';
                    priority(end+1)=2;
                    
                    descrip{end+1}=sprintf('any temporal changes in functional connectivity strength during %s',ucond{n1});
                    neffects{end+1}=Tmatch(icond==n1,3)';
                    ceffects{end+1}=diff(eye(numel(idx)));
                    str{end+1}='connectivity-change';
                    priority(end+1)=2;                    
                end
            end
        end
        [nill,idx]=sort(priority);
        varargout={struct('str',{str(idx)}, 'descrip',{descrip(idx)}, 'neffects',{neffects(idx)}, 'ceffects',{ceffects(idx)})};
        
    case {'suggest','suggest_between'} % 'suggest', names, values
        if numel(varargin)<1||isempty(varargin{1}), names=CONN_x.Setup.l2covariates.names(1:end-1); names=names(cellfun(@(x)isempty(regexp(x,'^_|^QA_|^QC_')),names)); 
        elseif isnumeric(varargin{1}), names=CONN_x.Setup.l2covariates.names(varargin{1}); 
        else names=varargin{1}; 
        end
        if numel(varargin)<2||isempty(varargin{2}), 
            [ok,idx]=ismember(names, CONN_x.Setup.l2covariates.names(1:end-1));
            names=names(ok);
            idx=idx(ok);
            X=zeros(CONN_x.Setup.nsubjects,length(idx));
            for nsub=1:CONN_x.Setup.nsubjects,
                for ncovariate=1:numel(idx);
                    X(nsub,ncovariate)=CONN_x.Setup.l2covariates.values{nsub}{idx(ncovariate)};
                end
            end
        else X=varargin{2}; 
        end
        Xnan=X;Xnan(isnan(Xnan))=0;
        DODOUBLE=2; % set to 1 to skip duplicate interaction descriptions
        str={}; descrip={}; neffects={}; ceffects={}; ctrl={}; ctrl_label={}; priority=[];
        isinvalid=all(isnan(X)|X==0,1);
        isconstant=all(isnan(X)|X==1,1)&~isinvalid;
        iscateg=all(isnan(X)|ismember(X,[0,1]),1)&~isinvalid&~isconstant;
        isnumer=~isinvalid&~isconstant&~iscateg;
        isinter=false(size(isnumer));
        kinter=sparse(size(X,2),size(X,2)); % [N,N] interaction matrix
        posinter=find(iscateg|isnumer);
        if ~isempty(posinter)
            tidx=find((iscateg|isnumer)&cellfun('length',regexp(names,'\*'))>0);
            for n1=1:numel(tidx),
                [ok,yidx]=ismember(regexp(names{tidx(n1)},'\*','split'),names(posinter));
                if all(ok)
                    isinter(tidx(n1))=true;
                    kinter(posinter(yidx),tidx(n1))=1;
                end
            end
            iscateg(isinter)=false;
            isnumer(isinter)=false;
        end
        iscategornumer=find(iscateg|isnumer|isinter);
        isconstant=find(isconstant);
        iscateg=find(iscateg);
        isnumer=find(isnumer);
        isinter=find(isinter);

        if ~isempty(isconstant)
            tidx=find(strcmp(names(isconstant),'AllSubjects'),1);
            if isempty(tidx), [nill,tidx]=max(sum(~isnan(X(:,isconstant)),1)); end
            tnames=names{isconstant(tidx)};
            %%%%%%%%%%%%
            if isequal(tnames,'AllSubjects'), descrip{end+1}='Does the average connectivity differ from zero?';
            else, descrip{end+1}=sprintf('Does the average connectivity within %s subjects differ from zero?',tnames);
            end
            ieffects=isconstant(tidx);
            neffects{end+1}=names(ieffects);
            ceffects{end+1}=1;
            str{end+1}=conn_wilkinson('describe',neffects{end});
            isorth=mean(abs(Xnan(:,iscategornumer)-Xnan(:,ieffects)*(pinv(Xnan(:,ieffects))*Xnan(:,iscategornumer))).^2,1)>1e-10&all(Xnan(all(Xnan(:,ieffects)==0,2),iscategornumer)==0,1);
            ctrl{end+1}=names(iscategornumer(isorth));
            ctrl_label{end+1}='when estimated at the zero-level of %s';
            priority(end+1)=0;
            for nnumer=1:numel(isnumer)
                %%%%%%%%%%%%
                if isequal(tnames,'AllSubjects'), descrip{end+1}=sprintf('Does the correlation between connectivity and %s differ from zero?',names{isnumer(nnumer)});
                else, descrip{end+1}=sprintf('Does the correlation between connectivity and %s within %s subjects differ from zero?',names{isnumer(nnumer)},tnames);
                end
                ieffects=[isconstant(tidx), isnumer(nnumer)];
                neffects{end+1}=names(ieffects);
                ceffects{end+1}=[0 1];
                str{end+1}=conn_wilkinson('describe',neffects{end});
                isorth=mean(abs(Xnan(:,iscategornumer)-Xnan(:,ieffects)*(pinv(Xnan(:,ieffects))*Xnan(:,iscategornumer))).^2,1)>1e-10&all(Xnan(all(Xnan(:,ieffects)==0,2),iscategornumer)==0,1);
                ctrl{end+1}=names(iscategornumer(isorth));
                ctrl_label{end+1}='after controlling for the influence of %s';
                priority(end+1)=10*nnumer+2;
                [nill,iidx]=sort(ieffects); neffects{end}=neffects{end}(iidx); ceffects{end}=ceffects{end}(:,iidx);
            end
            isconstant=isconstant(tidx);
        end
        if ~isempty(iscateg)
            x=X(:,iscateg);
            x(isnan(x))=0;
            sidx=[];
            for nkeep=1:numel(iscateg)
                ikeep=setdiff(1:nkeep,sidx);
                if ~isempty(ikeep)
                    [ux,nill,ix]=unique(x(:,ikeep),'rows');
                    for n=1:size(ux,1) % individual groups
                        tidx=find(all(repmat(ix==n,1,size(ux,2))==x(:,ikeep),1));
                        if ~isempty(tidx)
                            sidx=[sidx, ikeep(tidx(1))];
                        end
                    end
                end
            end
            % pairs of groups
            c=x(:,sidx)'*x(:,sidx)==0;
            c(1:size(c,1)+1:end)=true;
            maxi1=0; 
            factors={};
            for i1=1:size(c,1)
                if i1>maxi1
                    ok=0;
                    for i2=i1+1:size(c,1)
                        if nnz(c(i1:i2,i1:i2)==0), break;
                        else ok=i2;
                        end
                    end
                    if ok
                        maxi1=ok;
                        tidx=i1:ok;
                        if numel(tidx)==2
                            %%%%%%%%%%%%
                            descrip{end+1}=sprintf('Does the average connectivity differ between %s and %s subjects?',names{iscateg(sidx(tidx(1)))},names{iscateg(sidx(tidx(2)))});
                            ieffects=iscateg(sidx(tidx));
                            neffects{end+1}=names(ieffects);
                            ceffects{end+1}=[-1 1];
                            str{end+1}=conn_wilkinson('describe',neffects{end});
                            isorth=mean(abs(Xnan(:,iscategornumer)-Xnan(:,ieffects)*(pinv(Xnan(:,ieffects))*Xnan(:,iscategornumer))).^2,1)>1e-10&all(Xnan(all(Xnan(:,ieffects)==0,2),iscategornumer)==0,1);
                            ctrl{end+1}=names(iscategornumer(isorth));
                            ctrl_label{end+1}='after controlling for the influence of %s';
                            priority(end+1)=2;
                            [nill,iidx]=sort(ieffects); neffects{end}=neffects{end}(iidx); ceffects{end}=ceffects{end}(:,iidx);
                        else
                            %%%%%%%%%%%%
                            descrip{end+1}=sprintf('Does the average connectivity differ among %s and %s subjects?',conn_strjoinstr(names(iscateg(sidx(tidx(1:end-1)))),', '),names{iscateg(sidx(tidx(end)))});
                            ieffects=iscateg(sidx(tidx));
                            neffects{end+1}=names(ieffects);
                            ceffects{end+1}=diff(eye(numel(tidx)));
                            str{end+1}=conn_wilkinson('describe',neffects{end});
                            isorth=mean(abs(Xnan(:,iscategornumer)-Xnan(:,ieffects)*(pinv(Xnan(:,ieffects))*Xnan(:,iscategornumer))).^2,1)>1e-10&all(Xnan(all(Xnan(:,ieffects)==0,2),iscategornumer)==0,1);
                            ctrl{end+1}=names(iscategornumer(isorth));
                            ctrl_label{end+1}='after controlling for the influence of %s';
                            priority(end+1)=2;
                            [nill,iidx]=sort(ieffects); neffects{end}=neffects{end}(iidx); ceffects{end}=ceffects{end}(:,iidx);
                        end
                        for nnumer=1:numel(isnumer) % group x numer interactions
                            i3=find(sum(kinter,1)==2&kinter(isnumer(nnumer),:));
                            [ok4,i4]=max(kinter(iscateg(sidx(tidx)),i3),[],2);
                            if ~isempty(i3)&all(ok4)
                                for ntidx=1:numel(tidx)
                                    %%%%%%%%%%%%
                                    descrip{end+1}=sprintf('Does the correlation between connectivity and %s within %s subjects differ from zero?',names{isnumer(nnumer)},names{iscateg(sidx(tidx(ntidx)))});
                                    ieffects=[iscateg(sidx(tidx(ntidx))) i3(i4(ntidx))];
                                    neffects{end+1}=names(ieffects);
                                    ceffects{end+1}=[0 1];
                                    str{end+1}=conn_wilkinson('describe',neffects{end});
                                    isorth=mean(abs(Xnan(:,iscategornumer)-Xnan(:,ieffects)*(pinv(Xnan(:,ieffects))*Xnan(:,iscategornumer))).^2,1)>1e-10&all(Xnan(all(Xnan(:,ieffects)==0,2),iscategornumer)==0,1);
                                    ctrl{end+1}=names(iscategornumer(isorth));
                                    ctrl_label{end+1}='after controlling for the influence of %s';
                                    priority(end+1)=10*nnumer+3;
                                    [nill,iidx]=sort(ieffects); neffects{end}=neffects{end}(iidx); ceffects{end}=ceffects{end}(:,iidx);
                                end
                                %%%%%%%%%%%%
                                for ndouble=1:DODOUBLE
                                    if ndouble==1
                                        if numel(tidx)==2, descrip{end+1}=sprintf('Does the correlation between connectivity and %s differ between %s and %s subjects?',names{isnumer(nnumer)},names{iscateg(sidx(tidx(1)))},names{iscateg(sidx(tidx(2)))});
                                        else               descrip{end+1}=sprintf('Does the correlation between connectivity and %s differ among %s and %s subjects?',names{isnumer(nnumer)},conn_strjoinstr(names(iscateg(sidx(tidx(1:end-1)))),', '),names{iscateg(sidx(tidx(end)))});
                                        end
                                    else
                                        if numel(tidx)==2, descrip{end+1}=sprintf('Do the connectivity differences between %s and %s subjects depend on %s?',names{iscateg(sidx(tidx(1)))},names{iscateg(sidx(tidx(2)))},names{isnumer(nnumer)});
                                        else               descrip{end+1}=sprintf('Do the connectivity differences among %s and %s subjects depend on %s?',conn_strjoinstr(names(iscateg(sidx(tidx(1:end-1)))),', '),names{iscateg(sidx(tidx(end)))},names{isnumer(nnumer)});
                                        end
                                    end
                                    ieffects=[reshape(iscateg(sidx(tidx)),1,[]) reshape(i3(i4),1,[])];
                                    neffects{end+1}=names(ieffects);
                                    ceffects{end+1}=[zeros(numel(tidx)-1,numel(tidx)) diff(eye(numel(tidx)))];
                                    str{end+1}=conn_wilkinson('describe',neffects{end});
                                    isorth=mean(abs(Xnan(:,iscategornumer)-Xnan(:,ieffects)*(pinv(Xnan(:,ieffects))*Xnan(:,iscategornumer))).^2,1)>1e-10&all(Xnan(all(Xnan(:,ieffects)==0,2),iscategornumer)==0,1);
                                    ctrl{end+1}=names(iscategornumer(isorth));
                                    ctrl_label{end+1}='after controlling for the influence of %s';
                                    if ndouble==1, priority(end+1)=10*nnumer+3;
                                    else priority(end+1)=4;
                                    end
                                    [nill,iidx]=sort(ieffects); neffects{end}=neffects{end}(iidx); ceffects{end}=ceffects{end}(:,iidx);
                                end
                            end
                        end
                        for nfact=1:numel(factors) % group x factor interactions
                            i3=find(sum(kinter,1)==2);
                            [j1,j2]=ndgrid(factors{nfact},iscateg(sidx(tidx)));
                            [ok4,i4]=max(kinter(j1,i3)&kinter(j2,i3),[],2);
                            if ~isempty(i3)&all(ok4)
                                for ndouble=1:DODOUBLE
                                    %%%%%%%%%%%%
                                    if ndouble==1
                                        if numel(tidx)==2, descrip{end+1}=sprintf('Do the connectivity differences between %s and %s subjects depend on %s-status?',names{iscateg(sidx(tidx(1)))},names{iscateg(sidx(tidx(2)))},conn_strjoinstr(names(factors{nfact}),'/'));
                                        else               descrip{end+1}=sprintf('Do the connectivity differences among %sand %s subjects depend on %s-status?',conn_strjoinstr(names(iscateg(sidx(tidx(1:end-1)))),', '),names{iscateg(sidx(tidx(end)))},conn_strjoinstr(names(factors{nfact}),'/'));
                                        end
                                    else
                                        if numel(factors{nfact})==2, descrip{end+1}=sprintf('Do the connectivity differences between %s and %s subjects depend on %s-status?',names{factors{nfact}(1)},names{factors{nfact}(2)},conn_strjoinstr(names(iscateg(sidx(tidx))),'/'));
                                        else               descrip{end+1}=sprintf('Do the connectivity differences among %s and %s subjects depend on %s-status?',conn_strjoinstr(names(factors{nfact}(1:end-1)),', '),names{factors{nfact}(end)},conn_strjoinstr(names(iscateg(sidx(tidx))),'/'));
                                        end
                                    end
                                    ieffects=reshape(i3(i4),1,[]);
                                    neffects{end+1}=names(ieffects);
                                    ceffects{end+1}=kron(diff(eye(numel(tidx))),diff(eye(numel(factors{nfact}))));
                                    str{end+1}=conn_wilkinson('describe',neffects{end});
                                    isorth=mean(abs(Xnan(:,iscategornumer)-Xnan(:,ieffects)*(pinv(Xnan(:,ieffects))*Xnan(:,iscategornumer))).^2,1)>1e-10&all(Xnan(all(Xnan(:,ieffects)==0,2),iscategornumer)==0,1);
                                    ctrl{end+1}=names(iscategornumer(isorth));
                                    ctrl_label{end+1}='after controlling for the influence of %s';
                                    priority(end+1)=3;
                                    [nill,iidx]=sort(ieffects); neffects{end}=neffects{end}(iidx); ceffects{end}=ceffects{end}(:,iidx);
                                end
                            end
                        end
                        factors{end+1}=iscateg(sidx(tidx));
                        %                 for i2=1:numel(tidx)-1 % pairwise comparisons
                        %                     for i3=i2+1:numel(tidx)
                        %                         str{end+1}=conn_wilkinson('describe',names(iscateg(sidx(tidx([i2,i3])))));
                        %                         descrip{end+1}=sprintf('Does average connectivity differ between %s and %s subjects?',names{iscateg(sidx(tidx(i2)))},names{iscateg(sidx(tidx(i3)))});
                        %                     end
                        %                 end
                    elseif ~isempty(isconstant)
                        descrip{end+1}=sprintf('Does the average connectivity differ between %s and not-%s subjects?',names{iscateg(sidx(i1))},names{iscateg(sidx(i1))});
                        ieffects=[isconstant(1) iscateg(sidx(i1))];
                        neffects{end+1}=names(ieffects);
                        ceffects{end+1}=[0 1];
                        str{end+1}=conn_wilkinson('describe',neffects{end});
                        isorth=mean(abs(Xnan(:,iscategornumer)-Xnan(:,ieffects)*(pinv(Xnan(:,ieffects))*Xnan(:,iscategornumer))).^2,1)>1e-10&all(Xnan(all(Xnan(:,ieffects)==0,2),iscategornumer)==0,1);
                        ctrl{end+1}=names(iscategornumer(isorth));
                        ctrl_label{end+1}='after controlling for the influence of %s';
                        priority(end+1)=2;
                        [nill,iidx]=sort(ieffects); neffects{end}=neffects{end}(iidx); ceffects{end}=ceffects{end}(:,iidx);
                        for nnumer=1:numel(isnumer)
                            i3=find(sum(kinter,1)==2&kinter(isnumer(nnumer),:));
                            [ok4,i4]=max(kinter(iscateg(sidx(i1)),i3),[],2);
                            if ~isempty(i3)&all(ok4)
                                %%%%%%%%%%%%
                                descrip{end+1}=sprintf('Does the correlation between connectivity and %s differ between %s and not-%s subjects?',names{isnumer(nnumer)},names{iscateg(sidx(i1))},names{iscateg(sidx(i1))});
                                ieffects=[isconstant(1) isnumer(nnumer) iscateg(sidx(i1)) i3(i4)];
                                neffects{end+1}=names(ieffects);
                                ceffects{end+1}=[0 0 0 1];
                                str{end+1}=conn_wilkinson('describe',neffects{end});
                                isorth=mean(abs(Xnan(:,iscategornumer)-Xnan(:,ieffects)*(pinv(Xnan(:,ieffects))*Xnan(:,iscategornumer))).^2,1)>1e-10&all(Xnan(all(Xnan(:,ieffects)==0,2),iscategornumer)==0,1);
                                ctrl{end+1}=names(iscategornumer(isorth));
                                ctrl_label{end+1}='after controlling for the influence of %s';
                                priority(end+1)=10*nnumer+2;
                                [nill,iidx]=sort(ieffects); neffects{end}=neffects{end}(iidx); ceffects{end}=ceffects{end}(:,iidx);
                            end
                        end
                    end
                end
            end
            % individual groups
            for n1=1:numel(sidx)
                %%%%%%%%%%%%
                descrip{end+1}=sprintf('Does the average connectivity within %s subjects differ from zero?',names{iscateg(sidx(n1))});
                ieffects=iscateg(sidx(n1));
                neffects{end+1}=names(ieffects);
                ceffects{end+1}=1;
                str{end+1}=conn_wilkinson('describe',neffects{end});
                isorth=mean(abs(Xnan(:,iscategornumer)-Xnan(:,ieffects)*(pinv(Xnan(:,ieffects))*Xnan(:,iscategornumer))).^2,1)>1e-10&all(Xnan(all(Xnan(:,ieffects)==0,2),iscategornumer)==0,1);
                ctrl{end+1}=names(iscategornumer(isorth));
                ctrl_label{end+1}='when estimated at the zero-level of %s';
                priority(end+1)=1;
                [nill,iidx]=sort(ieffects); neffects{end}=neffects{end}(iidx); ceffects{end}=ceffects{end}(:,iidx);
            end
        end
        if ~isempty(isconstant) % numeric*numeric interactions
            i3=find(sum(kinter,1)==2&sum(kinter(isnumer,:),1)==2);
            for nnumer=1:numel(i3)
                i4=find(kinter(:,i3(nnumer)));
                for ndouble=1:DODOUBLE
                    %%%%%%%%%%%%
                    if ndouble==1, descrip{end+1}=sprintf('Does the correlation between connectivity and %s depend on %s?',names{i4(1)},names{i4(2)});
                    else descrip{end+1}=sprintf('Does the correlation between connectivity and %s depend on %s?',names{i4(2)},names{i4(1)});
                    end
                    ieffects=[isconstant(1) i4(1) i4(2) i3(nnumer)];
                    neffects{end+1}=names(ieffects);
                    ceffects{end+1}=[0 0 0 1];
                    str{end+1}=conn_wilkinson('describe',neffects{end});
                    isorth=mean(abs(Xnan(:,iscategornumer)-Xnan(:,ieffects)*(pinv(Xnan(:,ieffects))*Xnan(:,iscategornumer))).^2,1)>1e-10&all(Xnan(all(Xnan(:,ieffects)==0,2),iscategornumer)==0,1);
                    ctrl{end+1}=names(iscategornumer(isorth));
                    ctrl_label{end+1}='after controlling for the influence of %s';
                    if ndouble==1, priority(end+1)=10*find(isnumer==i4(1))+3;
                    else priority(end+1)=10*find(isnumer==i4(2))+3;
                    end
                    [nill,iidx]=sort(ieffects); neffects{end}=neffects{end}(iidx); ceffects{end}=ceffects{end}(:,iidx);
                end
            end
        end
        %[nill,alphrank]=sort(descrip); alphrank(alphrank)=(0:numel(alphrank)-1)/numel(alphrank);
        %[nill,idx]=sort(priority+alphrank);
        [nill,idx]=sort(priority);
        varargout={struct('str',{str(idx)}, 'descrip',{descrip(idx)}, 'neffects',{neffects(idx)}, 'ceffects',{ceffects(idx)}, 'ctrl', {ctrl(idx)}, 'ctrl_label', {ctrl_label(idx)})};
end     
end


function str=conn_strjoinstr(str1,str2)
str=[str1(:)';repmat({str2},1,length(str1))];
str=reshape(str(1:end-1),1,numel(str)-1);
str=[str{:}];
end

