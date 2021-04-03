function RT = conn_get_rt(nsub,nses,nset, functionals)

global CONN_x;

if nargin<4, functionals={}; end
if nargin<3||isempty(nset), nset=0; end
if nargin<2||isempty(nses), nses=[]; end
if nargin<1||isempty(nsub), nsub=[]; end

if isempty(nsub), nsub=1:CONN_x.Setup.nsubjects; end
if any(isnan(CONN_x.Setup.RT(min(numel(CONN_x.Setup.RT),nsub))))&&conn_projectmanager('inserver'), 
    for insub=nsub(:)'
        for inses=nses(:)'
            filename=conn_get_functional(insub,inses,nset);
            assert(~isempty(filename),'functional data for subject %d session %d not found',insub,inses);
            if size(filename,1)>1, filename=filename(1,:); end
            functionals{insub}{inses}=conn_server('util_localfile',filename);
        end
    end
    RT=conn_server('run',mfilename,nsub,nses,nset,functionals); 
    return; 
end
if numel(nsub)>1
    RT=[];
    for nsub=nsub(:)',
        RT=[RT conn_get_rt(nsub,nses,nset)];
    end
    return
else
    RT = CONN_x.Setup.RT(min(numel(CONN_x.Setup.RT),nsub));
    if isnan(RT)
        if isempty(nses), nses=1:CONN_x.Setup.nsessions(min(numel(CONN_x.Setup.nsessions),nsub)); end
        if numel(nses)>1
            for nses=nses(:)'
                RT=[RT conn_get_rt(nsub,nses,nset)];
            end
            RT=max(RT);
            return
        else
            if isempty(functionals), filename=conn_get_functional(nsub,nses,nset);
            else filename=functionals{nsub}{nses};
            end
            assert(~isempty(filename),'functional data for subject %d session %d not found',nsub,nses);
            if size(filename,1)>1, filename=filename(1,:); end
            [xtemp,failed]=conn_setup_preproc_meanimage(regexprep(filename,'\.gz$|\,\d+$',''),'json');
            if isempty(xtemp)
                try
                    v=nifti(filename);
                    RT=v.timing.tspace;
                catch
                    error('error reading RepetitionTime field; cannot find .json file or timing information associated with data file %s',filename);
                end
            else
                RT=conn_jsonread(xtemp,'RepetitionTime');
            end
            assert(~isempty(RT), 'error reading RepetitionTime field from json file %s',xtemp);
        end
    end
end