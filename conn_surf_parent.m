function V = conn_surf_parent(V,warnmsg,doerror)
if nargin<2, warnmsg=''; end
if nargin<3||isempty(doerror), doerror=false; end

if isstruct(V), fnames=conn_expandframe(V); 
else
    fnames=cellstr(V);
    V=conn_fileutils('spm_vol',char(fnames));
end

for n=1:numel(fnames),
    if ~isempty(regexp(fnames{n},'\.surf\.(nii|img)(,\d+)?$')),
        fname=regexprep(fnames{n},'\.surf\.(nii|img)(,\d+)?$','.vol.$1$2'); ok=conn_existfile(fname);
        if ~ok, fname=regexprep(fnames{n},'\.surf\.(nii|img)(,\d+)?$','.$1$2'); ok=conn_existfile(fname); end
        if ~ok, [fname_path,fname_name,fname_ext]=fileparts(fname); fname_name=regexprep(fname_name,'^s+',''); fname=fullfile(fname_path,[fname_name,fname_ext]); ok=conn_existfile(fname); end
        %if ~ok, fname=conn_prepend(-1,fname); ok=conn_existfile(fname); end
        %if ~ok, fname=conn_prepend(-1,fname); ok=conn_existfile(fname); end
        if ok,
            if ~isempty(warnmsg), fprintf('%s warning: displaying slices in %s instead of in %s (surface file)\n',warnmsg,fname,fnames{n}); end
            if n==1, V=conn_fileutils('spm_vol',fname);
            else V(n)=conn_fileutils('spm_vol',fname);
            end
        elseif doerror, error('%s error: %s is a non-compatible format (surface file)',warnmsg,fnames{n});
        elseif ~isempty(warnmsg), fprintf('%s warning: %s is a non-compatible format (surface file)\n',warnmsg,fnames{n});
        end
    end
end
end

