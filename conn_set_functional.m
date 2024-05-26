function nV=conn_set_functional(nsub,nses,nset,out)
global CONN_x;

if isempty(nset), nset=0; end
if ischar(nset), nset=conn_datasetlabel(nset,'add'); end

if ~nset, 
    [CONN_x.Setup.functional{nsub}{nses},nV]=conn_file(out); 
    CONN_x.Setup.nscans{nsub}{nses}=nV;
elseif nset<0
    [CONN_x.Setup.structural{nsub}{nses},nV]=conn_file(out); 
else
    CONN_x.Setup.secondarydataset(nset).functionals_type=4;
    [CONN_x.Setup.secondarydataset(nset).functionals_explicit{nsub}{nses},nV]=conn_file(out);
end
end
