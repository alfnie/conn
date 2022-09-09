function varargout = evlab17_modelplots(varargin)
% EVLAB17_MODELPLOTS displays first-level analysis results
%   evlab17_modelplots('/myfolder/evlab17_foo.mat'); 
%      displays first-level results associated with dataset evlab17_foo
%   evlab17_modelplots('/myfolder/evlab17_foo.mat',analysis,contrast); 
%      displays selected first-level analysis and contrast (by contrast number or by contrast name)
%   fh = evlab17_modelplots(...,'stats');
%      displays first-level statistics
%

varargout=cell(1,nargout);
evlab17_module init silent;

% interprets info
cwd=pwd;
if numel(varargin)>0&&~isempty(varargin{1}),  
    evlab17_module('load',varargin{1}); 
    if numel(varargin)>1, % select firstlevel analysis
        model_name=regexprep(varargin{2},'^firstlevel_','');
        prependmodelname=true;
        if evlab17_module('inconnfolders'),
            model_folder=fullfile(conn_prepend('',evlab17_module('filename'),''),'results','firstlevel');
            prependmodelname=false;
        else model_folder=fullfile(fileparts(fileparts(evlab17_module('filename'))));
        end
        model_name=fullfile(model_folder,model_name);
        if prependmodelname, model_name=conn_prepend('firstlevel_',model_name); end
        if evlab17_module('get','Setup.nsubjects')>1,
            files=evlab17_module('get','spm');
            subjects=1:evlab17_module('get','Setup.nsubjects');
            for nsubject=subjects(:)'
                files(nsubject)={fullfile(model_name,sprintf('sub-%04d',nsubject),'SPM.mat')};
            end
            evlab17_module('set','spm',files);
        else
            files=evlab17_module('get','spm');
            files(1)={fullfile(model_name,'SPM.mat')};
            evlab17_module('set','spm',files);
        end
    end
end
if numel(varargin)>0&&isequal(varargin{end},'stats'), 
    varargin=varargin(1:end-1);
    files=evlab17_module('get','spm');
    if numel(varargin)>=3, for n=1:numel(files), fh=conn_display(files{n},varargin{3},[],[],1); end
    else for n=1:numel(files), fh=conn_display(files{n},'?',[],[],1); end
    end
    varargout={fh};
else
    conn_displaysubject(varargin{3:end});
end
cd(cwd);


end

