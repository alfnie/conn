function out=conn_existfile(filename,acceptdir)
% conn_existfile(filename [,opt])
% check if filename exists
%   opt : 0 (default) filename must be a valid file name (not a folder name)
%         1 filename could be a file or a folder name
%         2 filename must be a valid folder name (not a file name)
% 

if nargin<2||isempty(acceptdir), acceptdir=false; end

out=false;
if iscell(filename)||size(filename,1)>1
    filename=cellfun(@strtrim,cellstr(filename),'uni',0);
    [filepath,filename,fileext,filenum]=cellfun(@spm_fileparts,filename,'uni',0);
    fullfilename=cellfun(@(a,b,c)fullfile(a,[b c]),filepath,filename,fileext,'uni',0);
    [ufullfilename,nill,idx]=unique(fullfilename);
    if acceptdir>1, out=cellfun(@isdir,ufullfilename);                                                       % isdir only
    else
        try
            out=cellfun(@spm_existfile,ufullfilename);                                                              % isfile only
            if acceptdir&&any(~out), out(~out)=cellfun(@isdir,ufullfilename(~out)); end                             % isfile or isdir
        catch
            out=cellfun(@(x)~isempty(dir(x)),ufullfilename);                                                        % isfile or isdir
            if ~acceptdir&&any(out), out(out)=~cellfun(@isdir,ufullfilename(out)); end                              % isfile only
        end
    end
    out=out(idx);
else
    if isempty(deblank(filename)), out=false;
    else
        [filepath,filename,fileext,filenum]=spm_fileparts(deblank(filename));
        if acceptdir>1, out=isdir(fullfile(filepath,[filename,fileext]));
        else
            try
                out=spm_existfile(fullfile(filepath,[filename,fileext]));
                if acceptdir&&~out, out=isdir(fullfile(filepath,[filename,fileext])); end
            catch
                out=~isempty(dir(fullfile(filepath,[filename,fileext])));
                if ~acceptdir&&out, out=~isdir(fullfile(filepath,[filename,fileext])); end
            end
        end
    end
end
