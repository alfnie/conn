function varargout=evlab17_module(option,varargin)
%
% Initialize evlab17 software:
%    evlab17_module init
%
% Check evlab17 software internal default values:
%    evlab17_module default fieldname
% with fieldname = 'qa_plots', 'qa_parallel', 'qa_profile', 'mat_format'
%
% other uses for evlab17 internal functions only
%

global CONN_x EVLAB17_info;
persistent filename done_init defaults;

if isempty(defaults), defaults=struct('qa_plots',false,'qa_parallel',0,'qa_profile','Background process (Unix,Mac)','mat_format','-v7.3'); end

switch(lower(option))
    case {'init','initforce'} % misc initialization checks
        if strcmpi(option,'initforce') || isempty(done_init) || ~(isequal(done_init,0)||isequal(done_init,path))
            ATTEMPTSPMCONFLICTRESOLUTION=true;
            verfilename_CONN=fullfile(fileparts(which(mfilename)),'path.CONN.txt');
            if isempty(dir(verfilename_CONN)),
                fprintf('Using implicit CONN paths (path file %s not defined)\n',verfilename_CONN)
                folderOptions_conn={fileparts(fileparts(fileparts(which(mfilename)))), fullfile(fileparts(fileparts(which(mfilename))),'conn'),'/software/conn'}; % folders where to search for CONN
                FORCEVER_CONN=false;
            else
                folderOptions_conn=regexp(fileread(verfilename_CONN),'\n','split');
                FORCEVER_CONN=true; % forces only valid CONN versions to be accepted
            end
            verfilename_SPM=fullfile(fileparts(which(mfilename)),'path.SPM.txt');
            if isempty(dir(verfilename_SPM)),
                fprintf('Using implicit SPM paths (path file %s not defined)\n',verfilename_SPM)
                folderOptions_spm={fullfile(fileparts(fileparts(fileparts(fileparts(which(mfilename))))),'spm12'), fullfile(fileparts(fileparts(which(mfilename))),'spm12'),'/software/spm12'}; % folders where to search for SPM12
                FORCEVER_SPM=false;
            else
                folderOptions_spm=regexp(fileread(verfilename_SPM),'\n','split');
                FORCEVER_SPM=true; % forces only valid SPM versions to be accepted
            end
            verfilename_SPM_SS=fullfile(fileparts(which(mfilename)),'path.SPM_SS.txt');
            if isempty(dir(verfilename_SPM_SS)),
                fprintf('Using implicit SPM_SS paths (path file %s not defined)\n',verfilename_SPM_SS)
                folderOptions_spm_ss={fullfile(fileparts(fileparts(fileparts(fileparts(which(mfilename))))),'spm_ss'), fullfile(fileparts(fileparts(which(mfilename))),'spm_ss'),'/software/spm_ss'}; % folders where to search for SPM_SS
                FORCEVER_SPM_SS=false;
            else
                folderOptions_spm_ss=regexp(fileread(verfilename_SPM_SS),'\n','split');
                FORCEVER_SPM_SS=true; % forces only valid SPM_SS versions to be accepted
            end
            
            % checks versions of SPM/conn/spm_ss
            tfilename=fileparts(which('spm'));
            if ~isempty(tfilename)
                pathall=regexp(path,':','split');
                pathall=pathall(~strncmp(pathall,tfilename,numel(tfilename)));
                conflicts=regexp(sprintf('%s:',pathall{:}),'([^\:]*spm\d+[^\\\/]*)[\\\/]matlabbatch','tokens'); conflicts=unique([conflicts{:}]);
            else conflicts={};
            end
            if isempty(tfilename) || ~isempty(conflicts) || (FORCEVER_SPM&&~ismember(tfilename,folderOptions_spm))
                if ~isempty(tfilename), conflicts=[conflicts {tfilename}]; end
                for nc=1:numel(conflicts)
                    pathall=regexp(path,':','split');
                    fprintf('warning: temporarily removing path entry %s. Please remove this entry from your path to avoid this message in the future\n',conflicts{nc});
                    cellfun(@rmpath,pathall(strncmp(pathall,conflicts{nc},numel(conflicts{nc}))),'uni',0);
                end
                if ~isempty(conflicts)
                    if ATTEMPTSPMCONFLICTRESOLUTION
                        fprintf('warning: conflicting SPM-versions in your path (%s). Please edit your Matlab path and/or edit the list of allowed SPM installations in %s to avoid this warning in the future\n',sprintf('%s ',conflicts{:}),verfilename_SPM);
                        fprintf('IMPORTANT NOTE: switching between different SPM versions within the same Matlab session can lead to unexpected behavior/errors. EVLAB17 will attempt to switch to the correct version now. This requires closing all figures and clearing all variables and class definitions in your Matlab workspace. If this procedure does not finish correctly or results in errors please restart Matlab, change your Matlab paths to avoid these conflicts, and try again\n');
                        delete(findall(get(0,'children')));
                        evalin('base','clear all classes');
                        if ~nargout, evlab17_module(option,varargin{:});
                        else [varargout{1:nargout}]=evlab17_module(option,varargin{:});
                        end
                        return
                    else
                        error('Unable to proceed. Conflicting SPM-versions in your path (%s). Modify your Matlab path and/or edit the file %s to modify the list of allowed SPM installations',sprintf('%s ',conflicts{:}),verfilename_SPM);
                    end
                end
                for n=1:numel(folderOptions_spm),
                    if ~isempty(dir(folderOptions_spm{n})),
                        addpath(folderOptions_spm{n});
                        break
                    end
                end
                if isempty(which('spm')), error('SPM not found in any of the expected locations (%s). Please add the path to SPM toolbox version 12 or above and try again',sprintf('%s ',folderOptions_spm{:})); end
            end
            verfilename_SPM=fullfile(fileparts(which(mfilename)),'ver.SPM.txt');
            if isempty(dir(verfilename_SPM)), if ~FORCEVER_SPM, fprintf('Skipping SPM version/release check (whitelist file %s not defined)\n',verfilename_SPM); end
            else
                ver_SPM=spm('version');
                str=regexp(fileread(verfilename_SPM),'\n','split');
                assert(ismember(lower(ver_SPM),strtrim(lower(str))),'SPM release %s not among your list of allowed SPM versions/releases. Please restart Matlab and update the path to a valid SPM release (%s) before proceeding, or edit the file %s to add "%s" to the list of valid SPM releases, or delete the file %s to skip SPM version-checks in the future',ver_SPM,sprintf('%s ',str{:}),verfilename_SPM,ver_SPM,verfilename_SPM);
            end
            try, spm_get_defaults('mat.format',evlab17_module('default','mat_format')); end
            % checks versions of spm/CONN/spm_ss
            tfilename=fileparts(which('conn'));
            conflicts={};
            if isempty(tfilename) || ~isempty(conflicts) || (FORCEVER_CONN&&~ismember(tfilename,folderOptions_conn))
                if ~isempty(tfilename), conflicts=[conflicts {tfilename}]; end
                for nc=1:numel(conflicts)
                    pathall=regexp(path,':','split');
                    fprintf('warning: temporarily removing path entry %s\n  please remove this entry from your path to avoid this message in the future\n',conflicts{nc});
                    cellfun(@rmpath,pathall(strncmp(pathall,conflicts{nc},numel(conflicts{nc}))),'uni',0);
                end
                for n=1:numel(folderOptions_conn),
                    if ~isempty(dir(folderOptions_conn{n})),
                        addpath(folderOptions_conn{n});
                        break
                    end
                end
                if isempty(which('conn')), error('please add path to CONN toolbox version 18 or above and try again'); end
            end
            verfilename_CONN=fullfile(fileparts(which(mfilename)),'ver.CONN.txt');
            if isempty(dir(verfilename_CONN)), if ~FORCEVER_CONN, fprintf('Skipping CONN version/release check (whitelist file %s not defined)\n',verfilename_CONN); end
            else
                ver_CONN=conn('version');
                str=regexp(fileread(verfilename_CONN),'\n','split');
                assert(ismember(lower(ver_CONN),strtrim(lower(str))),'CONN release %s not among your list of allowed CONN versions/releases. Please restart Matlab and update the path to a valid CONN release (%s) before proceeding, or edit the file %s to add "%s" to the list of valid CONN releases, or delete the file %s to skip CONN version-checks in the future',ver_CONN,sprintf('%s ',str{:}),verfilename_CONN,ver_CONN,verfilename_CONN);
            end
            % checks versions of spm/conn/SPM_SS
            tfilename=fileparts(which('spm_ss'));
            conflicts={};
            if isempty(tfilename) || ~isempty(conflicts) || (FORCEVER_SPM_SS&&~ismember(tfilename,folderOptions_spm_ss))
                if ~isempty(tfilename), conflicts=[conflicts {tfilename}]; end
                for nc=1:numel(conflicts)
                    pathall=regexp(path,':','split');
                    fprintf('warning: temporarily removing path entry %s\n  please remove this entry from your path to avoid this message in the future\n',conflicts{nc});
                    cellfun(@rmpath,pathall(strncmp(pathall,conflicts{nc},numel(conflicts{nc}))),'uni',0);
                end
                for n=1:numel(folderOptions_spm_ss),
                    if ~isempty(dir(folderOptions_spm_ss{n})),
                        addpath(folderOptions_spm_ss{n});
                        break
                    end
                end
                if FORCEVER_SPM_SS&&isempty(which('spm_ss')), error('please add path to SPM_SS toolbox and try again'); end
            end
            verfilename_SPM_SS=fullfile(fileparts(which(mfilename)),'ver.SPM_SS.txt');
            if isempty(dir(verfilename_SPM_SS)), if ~FORCEVER_SPM_SS, fprintf('Skipping SPM_SS version/release check (whitelist file %s not defined)\n',verfilename_SPM_SS); end
            else
                try, 
                    ver_SPM_SS=spm_ss('version');
                    str=regexp(fileread(verfilename_SPM_SS),'\n','split');
                    assert(ismember(lower(ver_SPM_SS),strtrim(lower(str))),'SPM_SS release %s not among your list of allowed SPM_SS versions/releases. Please restart Matlab and update the path to a valid SPM_SS release (%s) before proceeding, or edit the file %s to add "%s" to the list of valid CONN releases, or delete the file %s to skip SPM_SS version-checks in the future',ver_SPM_SS,sprintf('%s ',str{:}),verfilename_SPM_SS,ver_SPM_SS,verfilename_SPM_SS);
                end
            end
            
            if strcmp(pwd,fileparts(which(mfilename))),
                fprintf('warning: adding current directory %s to Matlab path; please include the evlab17 software folder in your Matlab path and cd to a different folder to avoid this warning in the future\n',pwd);
                addpath(pwd);
            end
            conn_setpermissions;
            % checks potential over-shadowing issues
            conn_checkdistributionfiles evlab17 silent;
            conn_checkdistributionfiles spm silent;
            conn_checkdistributionfiles conn silent;
            if ~isempty(which('spm_ss')), conn_checkdistributionfiles spm_ss silent; end
            done_init=path;
        end
        if nargin<=1||~any(strcmpi(varargin(cellfun(@ischar,varargin)),'silent')),
            try, fprintf('EVLAB17 in %s\n',fileparts(which(mfilename))); end
            try, fprintf('SPM version %s in %s\n',spm('version'),fileparts(which('spm')));
            catch, fprintf('SPM version %s in %s\n',spm('ver'),fileparts(which('spm'))); end
            try, fprintf('CONN version CONN%s in %s\n',conn('ver'),fileparts(which('conn'))); end
            try, fprintf('SPM_SS version SPM_SS%s in %s\n',spm_ss('ver'),fileparts(which('spm_ss'))); end
        end
    case 'default'
        if isempty(varargin), varargout={defaults}; 
        elseif isfield(defaults,varargin{1})
            if numel(varargin)>1, defaults.(varargin{1})=varargin{2};
            else varargout={defaults.(varargin{1})};
            end
        else
            if numel(varargin)>1, conn_module('default',varargin{:});
            else varargout={conn_module('default',varargin{1})};
            end
            %error('unrecognized default %s',varargin{1});
        end
    case {'versioncheckreset'} 
        fprintf('EVLAB17 will now continue performing SPM/CONN version number/release/folder checks when run\n');
        done_init=[];
    case {'versioncheckskip'} 
        fprintf('warning: EVLAB17 will now stop performing SPM/CONN version number/release/folder checks when run\n');
        fprintf('Please use "evlab17_module versioncheckreset", or start a new Matlab session, to continue performing SPM/CONN version checks\n');
        done_init=0;
    case 'getinfo'
        for n=1:numel(varargin),
            try, 
                if ~iscell(EVLAB17_info.(varargin{n})), EVLAB17_info.(varargin{n})={EVLAB17_info.(varargin{n})}; end
                varargout{n}=EVLAB17_info.(varargin{n}){end};
            catch, varargout{n}=[];
            end
        end
    case 'setinfo'
        for n=1:2:numel(varargin)-1,
            if isfield(EVLAB17_info,varargin{n})&&~isempty(EVLAB17_info.(varargin{n})), 
                if ~iscell(EVLAB17_info.(varargin{n})), EVLAB17_info.(varargin{n})={EVLAB17_info.(varargin{n})}; end
                EVLAB17_info.(varargin{n}){end+1}=varargin{n+1};
            else EVLAB17_info.(varargin{n})={varargin{n+1}};
            end
        end
    case 'save',
        varargout={false};
        if numel(varargin)>0, filename=conn_prepend('',conn_fullfile(varargin{1}),'.mat'); end
        CONN_x.info=EVLAB17_info;
        save(filename,'CONN_x');
        varargout={true};
    case 'load',
        varargout={false};
        filename=conn_prepend('',conn_fullfile(varargin{1}),'.mat');
        load(filename,'CONN_x');
        conn_updatefolders;
        if isfield(CONN_x,'info'), EVLAB17_info=CONN_x.info;
        else EVLAB17_info=[];
        end
        if isfield(CONN_x,'ispending')&&CONN_x.ispending&&conn_jobmanager('ispending'), evlab17_module('update',filename); end
        varargout={exist('CONN_x','var')};
    case 'update',
        varargout={false};
        if numel(varargin)>0, filename=conn_prepend('',conn_fullfile(varargin{1}),'.mat'); end
        if ~isempty(filename)
            conn('load',filename);
            conn save;
        end
        varargout={exist('CONN_x','var')};
    case 'filename',
        varargout={filename};
    case 'inconnfolders'
        varargout={false};
        try,
            %varargout={isdir(fullfile(conn_prepend('',filename,''),'results','firstlevel'))};
            varargout={~isempty(conn_module('get','filename'))};
        end
    case 'spm'
        evlab17_module init silent;
        if ~nargout, spm(varargin{:});
        else [varargout{1:nargout}]=spm(varargin{:});
        end
    case 'conn'
        evlab17_module init silent;
        if ~nargout, conn(varargin{:});
        else [varargout{1:nargout}]=conn(varargin{:});
        end
    case 'spm_ss'
        evlab17_module init silent;
        if ~nargout, spm_ss(varargin{:});
        else [varargout{1:nargout}]=spm_ss(varargin{:});
        end
    case 'submit'
        evlab17_module init silent;
        if ~nargout, conn('submit',@evlab17,varargin{:}); % e.g. evlab17 submit run_model file.cfg
        else [varargout{1:nargout}]=conn('submit',@evlab17,varargin{:});
        end
    otherwise
        if ~isempty(regexp(lower(option),'^spm_|^conn_|^spm_ss_|^spm_bcc_'))
            evlab17_module init silent;
            if ~isempty(which(option))
                fh=eval(sprintf('@%s',option));
                if ~nargout, feval(fh,varargin{:});
                else [varargout{1:nargout}]=feval(fh,varargin{:});
                end
            else error('unrecognized function %s',option);
            end
        else
            if ~nargout, conn_module(option,varargin{:});
            else [varargout{1:nargout}]=conn_module(option,varargin{:});
            end
        end
end
end
