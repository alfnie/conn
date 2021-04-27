function [ok,op] = conn_imcalc(filein,fileout,evaluatefunction,designmatrix,cvalues,cnames)

if nargin<4||isempty(designmatrix), designmatrix=''; end % (optional, evaluatefunction=='lin') design matrix
if nargin<5||isempty(cvalues), values=[]; else cvalues=cvalues(:)'; end %(optional, evaluatefunction==@fun) covariate values
if nargin<6||isempty(cnames), cnames=repmat({''},1,size(values,2)); end %(optional, evaluatefunction==@fun) covariate names
if isstruct(evaluatefunction)
    op=evaluatefunction;
elseif isequal(evaluatefunction,'lin'), % lin
    %tstring=CONN_x.Setup.conditions.model{ncondition}{2};
    tstring=designmatrix;
    if ischar(tstring), value=str2num(tstring); else value=tstring; end
    if ischar(tstring)&&isempty(value),
        ok=0;
        for n1=1:3,
            try
                switch(n1)
                    case 1, value=evalin('base',tstring);
                    case 2,
                        %x=cat(2,CONN_x.Setup.l2covariates.values{nsub}{:});
                        %tnames=CONN_x.Setup.l2covariates.names(1:end-1);
                        [nill,idx]=sort(-cellfun('length',cnames));
                        for n1=idx(:)',
                            for n2=fliplr(strfind(tstring,cnames{n1}))
                                tstring=[tstring(1:n2-1) '(' mat2str(cvalues(:,n1)') ')' tstring(n2+numel(cnames{n1}):end)];
                            end
                        end
                        value=evalin('base',tstring);
                    case 3,
                        tstring=regexprep(tstring,'([^\.])(\*)|([^\.])(/)|([^\.])(\^)','$1.$2');
                        value=evalin('base',tstring);
                end
                ok=1;
            end
            if ok, break; end
        end
        if ~ok,
            value=[];
            tstring0=designmatrix;
            tstring0=cellstr(tstring0);tstring0=sprintf('%s;',tstring0{:});
            if isequal(tstring0,tstring), conn_msgbox(['Unable to interpret string ',tstring0(:)'],'',2);
            else conn_msgbox({['Unable to interpret string ',tstring0],['Closest attempt (Matlab string) ',tstring]},'',2);
            end
        end
    end
    op.EvaluateFunction=evaluatefunction;
    op.DesignMatrix=value;
    if size(op.DesignMatrix,1)~=numel(filein), error('mismatch between Design Matrix rows (%d) and number of conditions (%d)',size(op.DesignMatrix,1),numel(filein)); end
    op.CovariateValuesvalues=cvalues;
    op.CovariateNames=cnames;
    %op.CovariateValues=[CONN_x.Setup.l2covariates.values{nsub}{:}];
    %op.CovariateNames=CONN_x.Setup.l2covariates.names(1:end-1);
elseif ~ischar(evaluatefunction) % @fun
    op.EvaluateFunction=evaluatefunction;
    op.CovariateValuesvalues=cvalues;
    op.CovariateNames=cnames;
    %op.CovariateValues=[CONN_x.Setup.l2covariates.values{nsub}{:}];
    %op.CovariateNames=CONN_x.Setup.l2covariates.names(1:end-1);
else % avg,std,...
    op.EvaluateFunction=evaluatefunction;
end

[nill,nill,fext]=fileparts(fileout);
if strcmp(fext,'.mat') % ROI-to-ROI files
    filein=cellstr(filein);
    if ~all(conn_existfile(filein)),ok=false; return; end 
    Vin=cellfun(@load,filein,'uni',0); %load(filename,'Z','regressors','names','names2','xyz','SE','DOF');
    emptycondition=any(cellfun(@(x)all(isnan(x.Z(:))),Vin));
    Vout=Vin{1};
    if emptycondition, Vout.Z(:)=nan; end
    if ~emptycondition,
        y=[];
        for n=1:numel(Vin)
            if ~isequal(Vout.names,Vin{n}.names)||~isequal(Vout.names2,Vin{n}.names2), error('mismatch ROI names in %s and %s. Re-run previous step (Denoising) and try again',filein{1},filein{n}); end
            y=cat(3,y,Vin{n}.Z);
        end
        y=permute(y,[3,1,2]);
        y=y(:,:);
        if ischar(op.EvaluateFunction)
            switch(op.EvaluateFunction)
                case 'avg', t=mean(y,1);
                case 'std', t=std(y,0,1);
                case 'min', t=min(y,[],1);
                case 'max', t=max(y,[],1);
                case 'lin', t=pinv(op.DesignMatrix'*op.DesignMatrix)*op.DesignMatrix'*y; t=t(:,1).';
                otherwise, error('unrecognized function keyword %s (only ''avg'' or ''std'' keywords allowed; use explicit function handle instead)',op.EvaluateFunction);
            end
        elseif nargin(op.EvaluateFunction)>1
            if nargin(op.EvaluateFunction)>2, t=feval(op.EvaluateFunction,y,op.CovariateValues,op.CovariateNames);
            else t=feval(op.EvaluateFunction,y,op.CovariateValues);
            end
        else
            t=feval(op.EvaluateFunction,y);
        end
        Vout.Z=reshape(t,size(Vout.Z));
        ok=true;
    else ok=false;
    end
    Z=Vout.Z; regressors=Vout.regressors; names=Vout.names; names2=Vout.names2; xyz=Vout.xyz; SE=Vout.SE; DOF=Vout.DOF;
    save(fileout,'Z','regressors','names','names2','xyz','SE','DOF');
else % nifti files
    if ~all(conn_existfile(cellstr(filein))),ok=false; return; end 
    Vin=spm_vol(char(filein));
    emptycondition=any(ismember({Vin.descrip},'CONNlabel:MissingData'));
    Vout=struct('mat',Vin(1).mat,'dim',Vin(1).dim,'fname',fileout,'pinfo',[1;0;0],'n',[1,1],'dt',[spm_type('float32') spm_platform('bigend')]);
    if emptycondition, Vout.descrip='CONNlabel:MissingData'; end
    Vout=spm_create_vol(Vout);
    if ~emptycondition,
        [gridx,gridy]=ndgrid(1:Vin(1).dim(1),1:Vin(1).dim(2));
        xyz0=[gridx(:),gridy(:)]';
        for slice=1:Vin(1).dim(3)
            xyz=[xyz0; slice+zeros(1,size(xyz0,2)); ones(1,size(xyz0,2))];
            y=spm_get_data(Vin(:)',xyz);
            if ischar(op.EvaluateFunction)
                switch(op.EvaluateFunction)
                    case 'avg', t=mean(y,1);
                    case 'std', t=std(y,0,1);
                    case 'min', t=min(y,[],1);
                    case 'max', t=max(y,[],1);
                    case 'lin', t=pinv(op.DesignMatrix'*op.DesignMatrix)*op.DesignMatrix'*y; t=t(:,1).';
                    otherwise, error('unrecognized function keyword %s (only ''avg'' or ''std'' keywords allowed; use explicit function handle instead)',op.EvaluateFunction);
                end
            elseif nargin(op.EvaluateFunction)>1
                if nargin(op.EvaluateFunction)>2, t=feval(op.EvaluateFunction,y,op.CovariateValues,op.CovariateNames);
                else t=feval(op.EvaluateFunction,y,op.CovariateValues);
                end
            else
                t=feval(op.EvaluateFunction,y);
            end
            t=reshape(t,Vin(1).dim(1:2));
            Vout=spm_write_plane(Vout,t,slice);
        end
        ok=true;
    else ok=false;
    end
end