function [filename,altout]=conn_get_functional(nsub,nses,nset,flag,alt)
global CONN_x;

if nargin<5||isempty(alt), alt=[]; end
if nargin<4||isempty(flag), flag=false; end
if nargin<3||isempty(nset), nset=0; end
if ischar(nset), nset=conn_datasetlabel(nset,'error'); end

filename=[];altout=[];
if nset>0
    try
        if numel(CONN_x.Setup.secondarydataset)<nset, CONN_x.Setup.secondarydataset(nset)=CONN_x.Setup.secondarydataset(end); end
        if CONN_x.Setup.secondarydataset(nset).functionals_type==4
            if numel(CONN_x.Setup.secondarydataset(nset).functionals_explicit)<nsub||numel(CONN_x.Setup.secondarydataset(nset).functionals_explicit{nsub})<nses, CONN_x.Setup.secondarydataset(nset).functionals_explicit{nsub}{nses}={[],[],[]}; end
            filename=CONN_x.Setup.secondarydataset(nset).functionals_explicit{nsub}{nses}{1};
        else
            if numel(CONN_x.Setup.functional)<nsub||numel(CONN_x.Setup.functional{nsub})<nses, CONN_x.Setup.functional{nsub}{nses}={[],[],[]}; end
            tV=CONN_x.Setup.functional{nsub}{nses}{1};
            if ~isempty(tV), tV=conn_rulebasedfilename(cellstr(tV),CONN_x.Setup.secondarydataset(nset).functionals_type,CONN_x.Setup.secondarydataset(nset).functionals_rule); end
            filename=char(tV);
        end
        if flag, % check if file exists
            existfile=conn_existfile(cellstr(filename)); %existunsmoothed=cellfun(@conn_existfile,VsourceUnsmoothed);
            if ~all(existfile),
                filename=CONN_x.Setup.functional{nsub}{nses}{1};
                conn_disp('fprintf','warning: set-%d data for subject %d session %d not found\n',nset,nsub,nses);
            end
        end
    catch
        if flag, % check if file exists
            filename=CONN_x.Setup.functional{nsub}{nses}{1};
            conn_disp('fprintf','warning: error in dataset definition for set-%d subject %d session %d. Using set-0 functional data instead\n',nset,nsub,nses);
        else
            conn_disp('fprintf','warning: error in dataset definition for set-%d subject %d session %d\n',nset,nsub,nses);
        end
    end
else
    filename=CONN_x.Setup.functional{nsub}{nses}{1};
end
if ~isempty(alt)&&~isempty(filename)
    switch(lower(alt))
        case 'cfile'
            if ~nset, altout=CONN_x.Setup.functional{nsub}{nses};
            elseif CONN_x.Setup.secondarydataset(nset).functionals_type==4, altout=CONN_x.Setup.secondarydataset(nset).functionals_explicit{nsub}{nses};
            else altout=conn_file(filename);
            end
    end
end



end
