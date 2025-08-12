function data = conn_get_l1covariate(nsubject,nl1covariate,nses,varargin)
% internal function

global CONN_x;
if nargin<3, nses=1; end

data=[];
covfilename=CONN_x.Setup.l1covariates.files{nsubject}{nl1covariate}{nses}{1};
switch(covfilename),
    case '[raw values]',
        data=CONN_x.Setup.l1covariates.files{nsubject}{nl1covariate}{nses}{3};
    otherwise,
        data=conn_loadtextfile(covfilename,false,varargin{:});
end
