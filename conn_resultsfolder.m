function [newfoldername,foldername,iscryptic]=conn_resultsfolder(ftype,state,varargin)
global CONN_x;

MAXFOLDERLENGTH=64;
switch(ftype)
    case 'subjectsconditions'
        % backwards compatible folder names
        if ismember(state,[1 2])
            foldername={[CONN_x.Analyses(CONN_x.Analysis).name,'.']};
            newfoldername={[CONN_x.Analyses(CONN_x.Analysis).name,filesep]};
            foldername{end+1}=['SUBJECT_EFFECTS_'];
        else
            foldername={'SUBJECT_EFFECTS_'};
            newfoldername={[CONN_x.vvAnalyses(CONN_x.vvAnalysis).name,filesep]};
        end
        doexist=false;
        if ~isempty(varargin)
            [nsubjecteffects,csubjecteffects,nconditions,cconditions]=deal(varargin{1:4});
            [doexist,nill,doexistname]=conn_contrastmanager('check',nsubjecteffects,csubjecteffects,nconditions,cconditions);
            if doexist, 
                altfoldername=newfoldername;
                altfoldername{end+1}=doexistname;
            end
            for n1=1:length(nsubjecteffects),
                if ~any(diff(csubjecteffects)),
                    foldername{end+1}=[CONN_x.Setup.l2covariates.names{nsubjecteffects(n1)},'.'];
                    newfoldername{end+1}=[CONN_x.Setup.l2covariates.names{nsubjecteffects(n1)},'.'];
                else,
                    foldername{end+1}=[CONN_x.Setup.l2covariates.names{nsubjecteffects(n1)},'(',num2str(csubjecteffects(:,n1)','%2.1f'),')','.'];
                    newfoldername{end+1}=[CONN_x.Setup.l2covariates.names{nsubjecteffects(n1)},'(',mat2str(csubjecteffects(:,n1)',3),')','.'];
                end
            end
            newfoldername{end+1}=filesep;
            if isequal(cconditions,'var'),
                foldername{end+1}=['CONDITIONS_var_'];
                newfoldername{end+1}=['var_'];
            else foldername{end+1}=['CONDITIONS_'];
            end
            for n1=1:length(nconditions),
                if isequal(cconditions,'var')||~any(diff(cconditions(:))),
                    foldername{end+1}=[CONN_x.Setup.conditions.names{nconditions(n1)},'.'];
                    newfoldername{end+1}=[CONN_x.Setup.conditions.names{nconditions(n1)},'.'];
                else,
                    foldername{end+1}=[CONN_x.Setup.conditions.names{nconditions(n1)},'(',num2str(cconditions(:,n1)','%2.1f'),')','.'];
                    newfoldername{end+1}=[CONN_x.Setup.conditions.names{nconditions(n1)},'(',mat2str(cconditions(:,n1)',3),')','.'];
                end
            end
        end
        foldername=strcat(foldername{:});
        foldername(foldername==' '|foldername==filesep)='_';
        n=numel(foldername);
        newfoldername=strcat(newfoldername{:});
        
        if doexist,newfoldername={strcat(altfoldername{:}) newfoldername};
        else newfoldername={newfoldername};
        end
        newfoldername=regexprep(newfoldername,'\.(\\)|\.(\/)','$1');
        newfoldername=regexprep(newfoldername,'[^\w\d\s\.\(\)\\\/\-_@&]+|^\.|\.$','');
        newfoldername=regexprep(deblank(newfoldername),'\s+','_');
        if ismember(state,[2,3])&&numel(foldername)>100, foldername=foldername(1:100); end
        for n1=1:numel(newfoldername)
            if isnan(MAXFOLDERLENGTH)&&numel(newfoldername{n1})>250,
                iscryptic(n1)=true;
                newfoldername{n1}=newfoldername{n1}(1:250);
            elseif numel(newfoldername{n1})>MAXFOLDERLENGTH
                iscryptic(n1)=true;
                newfoldername{end+1}=newfoldername{n1};
                n=numel(newfoldername{n1});
                newfoldername{n1}=[newfoldername{n1}(1:max(0,MAXFOLDERLENGTH-16)) char('0'+mod(round(newfoldername{n1}*(sin(6*(1:n)'*(1:16)/n))),10))];
            else
                iscryptic(n1)=false;
            end
        end
        newfoldername=[newfoldername {foldername}];
        iscryptic=[iscryptic false];
        foldername=newfoldername(2:end); % compatibility with older versions
        newfoldername=newfoldername{1};  % current-version
        if doexist, foldername={}; iscryptic=iscryptic(1); end % saved contrasts take precedence always

    case {'sources','sourcesmatched'}
        [sources,nsources,csources]=deal(varargin{1:3});
        REDUCENAMES=(state==2);
        foldername={};
        if state==3
            if strcmp(ftype,'sourcesmatched')
                for n=1:numel(sources)
                    idx=strmatch(sources{n},CONN_x.vvAnalyses(CONN_x.vvAnalysis).measures,'exact');
                    if isempty(idx), idx=strmatch(sources{n},CONN_x.vvAnalyses(CONN_x.vvAnalysis).measures); end
                    if numel(idx)==1, sources{n}=CONN_x.vvAnalyses(CONN_x.vvAnalysis).measures{idx(1)}; end
                end
            end
            fullsources={};for n2=1:length(sources),fullsources{end+1}=conn_v2v('pcleartext',sources{n2}); end
            [nill,idxcharstart]=max(any(diff(double(char(fullsources(nsources))),1,1),1));
        else
            if strcmp(ftype,'sourcesmatched')
                for n=1:numel(sources)
                    idx=strmatch(sources{n},CONN_x.Analyses(CONN_x.Analysis).sources,'exact');
                    if isempty(idx), idx=strmatch(sources{n},CONN_x.Analyses(CONN_x.Analysis).sources); end
                    if numel(idx)==1, sources{n}=CONN_x.Analyses(CONN_x.Analysis).sources{idx(1)}; end
                end
            end
            fullsources=sources;
            [nill,idxcharstart]=max(any(diff(double(char(sources(nsources))),1,1),1));
        end
        if REDUCENAMES
            redsources=regexprep(fullsources,{'^BA\.(\d+) \(([LR])\)\. .*','^\((-?\d+),(-?\d+),(-?\d+)\)$','^SLrois\.|^aal\.|^atlas\.|^networks\.','\s\(([LlRr])\)','([^\(\)]+)\(.+\)\s*$'},{'$1$2','($1 $2 $3)','',' ${lower($1)}','$1'});
            [nill,nill,i2]=unique(redsources);
            validreduced=ismember(i2,find(accumarray(i2(:),1)==1));
            redsources(~validreduced)=fullsources(~validreduced);
        end
        newfoldername=foldername;
        for n2=1:length(nsources),
            if REDUCENAMES, txttmp=redsources{nsources(n2)}; 
            else txttmp=fullsources{nsources(n2)}; 
            end
            if ~REDUCENAMES&&n2>1&&idxcharstart>1
                txttmp=txttmp(idxcharstart:end);
            end            
            if ~any(diff(csources(:))),
                foldername{end+1}=[txttmp,'.'];
                newfoldername{end+1}=[txttmp,'.'];
            else, 
                foldername{end+1}=[txttmp,'(',num2str(csources(:,n2)','%2.1f'),')','.'];
                newfoldername{end+1}=[txttmp,'(',mat2str(csources(:,n2)',3),')','.'];
            end
        end
        foldername=strcat(foldername{:});
        foldername(foldername==' '|foldername==filesep|foldername=='*')='_';
        %newfoldername=char('0'+mod(foldername*sign(sin(6*(1:n)'*(1:16)/n)),10));
        if numel(foldername)>100, foldername=foldername(1:100); end
        newfoldername=strcat(newfoldername{:});
        newfoldername=regexprep(newfoldername,'[^\w\d\s\.\(\)\-_@&]+|^\.|\.$','');
        newfoldername=regexprep(deblank(newfoldername),'\s+','_');
        if isnan(MAXFOLDERLENGTH)&&numel(newfoldername)>250,
            iscryptic=true;
            newfoldername=newfoldername(1:250);
        elseif numel(newfoldername)>MAXFOLDERLENGTH
            iscryptic=true;
            foldername=newfoldername;
            n=numel(newfoldername);
            newfoldername=[newfoldername(1:max(0,MAXFOLDERLENGTH-16)) char('0'+mod(round(newfoldername*(sin(6*(1:n)'*(1:16)/n))),10))];
        else
            iscryptic=false;
        end
        foldername={foldername};
end
end

