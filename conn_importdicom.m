function [ok,ERR]=conn_importdicom(filenames,varargin);
global CONN_x;

options=struct('type','functional',...
               'localcopy',false,...
               'nset',0,...
               'bidsname','',...
               'sessions',0,...
               'subjects_id',[],...
               'subjects',[]);
for n1=1:2:nargin-2, if ~isfield(options,lower(varargin{n1})), error('unknown option %s',lower(varargin{n1})); else options.(lower(varargin{n1}))=varargin{n1+1}; end; end
if ~isempty(options.subjects), nsubs=options.subjects; else nsubs=1:CONN_x.Setup.nsubjects; end
if ~isempty(options.bidsname), bidsname=options.bidsname;
elseif ~options.nset,          bidsname='func';
else                           bidsname=sprintf('dataset%dfunc',options.nset);
end

ok=false;
ERR={};
file_name=[];
for isub=1:numel(nsubs),
    if isempty(options.subjects_id), filename=filenames{isub};
    else
        if isempty(file_name), [nill,file_name,nill]=cellfun(@fileparts,filenames,'uni',0); end
        id=options.subjects_id{isub};
        select=cellfun('length',regexp(file_name,['^sub-',id,['^sub-',id,'$|^sub-',id,'_']]))>0;
        filename=filenames(select);
    end
    nsub=nsubs(isub);
    if isempty(options.sessions), nses=[];
    elseif iscell(options.sessions), nses=options.sessions{isub};
    else nses=options.sessions;
    end
    if isequal(nses,0)
        if ~isequal(options.type,'fieldmap'), % match number of sessions to input data
            if numel(CONN_x.Setup.nsessions)==1&&CONN_x.Setup.nsubjects>1, CONN_x.Setup.nsessions=CONN_x.Setup.nsessions+zeros(1,CONN_x.Setup.nsubjects); end
            CONN_x.Setup.nsessions(nsub)=numel(filename);
        end
        nses=[];
    end
    nsess=CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsub));
    if isempty(nses), nsesstemp=1:nsess;
    else nsesstemp=nses;
    end
    switch(options.type)
        case 'functional'
            if numel(filename)~=numel(nsesstemp), ERR{end+1}=sprintf('Subject %d functional import skipped : mismatched number of selected files (%d) and number of target sessions (%d)',nsub,numel(filename),numel(nsesstemp));
            else
                for n2=1:numel(nsesstemp)
                    nses=nsesstemp(n2);
                    if options.localcopy,
                        [nill,nill,nV]=conn_importvol2bids(filename{n2},nsub,nses,bidsname);
                        if ~options.nset, CONN_x.Setup.nscans{nsub}{nses}=nV; end
                    else
                        nV=conn_set_functional(nsub,nses,options.nset,filename{n2});
                    end
                end
            end
        case 'fieldmap'
            for n2=1:numel(nsesstemp)
                nses=nsesstemp(n2);
                if options.localcopy,
                    [nill,nill,nV]=conn_importvol2bids(char(filename),nsub,nses,'fmap','fmap');
                else
                    nV=conn_set_functional(nsub,nses,'fmap',char(filename));
                end
            end
        case 'structural'
            if CONN_x.Setup.structural_sessionspecific
                if numel(filename)~=numel(nsesstemp), ERR{end+1}=sprintf('Subject %d structural import skipped : mismatched number of selected files (%d) and number of target sessions (%d)',nsub,numel(filename),numel(nsesstemp));
                else
                    for n2=1:numel(nsesstemp)
                        nses=nsesstemp(n2);
                        if options.localcopy, conn_importvol2bids(filename{n2},nsub,nses,'anat');
                        else CONN_x.Setup.structural{nsub}{nses}=conn_file(filename{n2});
                        end
                    end
                end
            else
                if numel(filename)~=1, ERR{end+1}=sprintf('Subject %d structural import skipped : selected more than one file (%d)',nsub,numel(filename));
                else
                    for n2=1:numel(nsesstemp)
                        nses=nsesstemp(n2);
                        if options.localcopy, conn_importvol2bids(filename{1},nsub,[1,nses],'anat');
                        else CONN_x.Setup.structural{nsub}{nses}=conn_file(filename{1});
                        end
                    end
                end
            end
            
        otherwise, error('unrecognized type %s (structural|functional)',options.type);
    end
end
ok=isempty(ERR);

end