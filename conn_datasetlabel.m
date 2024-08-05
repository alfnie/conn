function [nset,isnew] = conn_datasetlabel(label,iffail)
global CONN_x;

if nargin<2||isempty(iffail), iffail='none'; end
nset=find(ismember({CONN_x.Setup.secondarydataset.label},label),1);
if isempty(nset), nset=find(ismember(lower({CONN_x.Setup.secondarydataset.label}),lower(label)),1); end
if isempty(nset), nset=find(ismember(lower({CONN_x.Setup.secondarydataset.label}),regexprep(label,{'^original$','^subjectspace$','^surfacespace$','^mnispace$','^smoothed$'},{'original functional data','subject-space functional data','surface-space functional data','mni-space functional data','smoothed functional data'})),1); end
isnew=isempty(nset);
if isnew
    switch(lower(iffail))
        case 'error', error('unable to match secondary dataset label %s',label);
        case 'add',
            nset=numel(CONN_x.Setup.secondarydataset)+1;
            CONN_x.Setup.secondarydataset(nset)=CONN_x.Setup.secondarydataset(1);
            CONN_x.Setup.secondarydataset(nset).label=label;
        case 'none',
        otherwise, error('unknown option %s',iffail);
    end
end
