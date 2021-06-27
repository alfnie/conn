function varargout=evlab17_run(option,varargin)
% EVLAB17_RUN Ev Lab fMRI preprocessing/analysis procedures
%
% See also evlab17_run_preproc, evlab17_run_model, evlab17_qaplots
%

varargout={[]};
switch(lower(option))
    otherwise
        if ~isempty(which(sprintf('evlab17_run_%s',option))),
            fh=eval(sprintf('@evlab17_run_%s',option));
            if nargout, [varargout{1:nargout}]=feval(fh,varargin{:});
            else feval(fh,varargin{:});
            end
        else
            disp(sprintf('unrecognized option %s or evlab17_run_%s function',option,option));
        end
end
end