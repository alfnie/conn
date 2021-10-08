function out = conn_jsonread(filename,fieldname,isnumericfield)
if nargin<3||isempty(isnumericfield), isnumericfield=true; end
if nargin<2||isempty(fieldname),fieldname=''; end
if any(conn_server('util_isremotefile',filename)), out=conn_server('run',mfilename,conn_server('util_localfile',filename),fieldname,isnumericfield); return; end

out=[];
if ~isempty(regexp(filename,'\.nii(\.gz)?(\,\d+)?$'))
    [tfilename,failed]=conn_setup_preproc_meanimage(regexprep(filename,'\.gz$|\,\d+$',''),'json');
    if isempty(tfilename), return; end %conn_disp('fprintf','warning: unable to find .json file associated with data file %s',filename); 
    filename=tfilename;
end
if nargin<2||isempty(fieldname),
    try, out=spm_jsonread(filename); end
    if isempty(out), try, out=jsondecode(fileread(filename)); end; end
else % note: faster but only for single-value numeric fields
    str=fileread(filename);
    if isnumericfield
        try, out=str2double(regexp(str,['"',fieldname,'"\s*:\s*([\d\.]+)'],'tokens','once'));
        end
        if isempty(out)
            try, out=str2double(regexp(str,['\.',fieldname,'"\s*:\s*([\d\.]+)'],'tokens','once'));
            end
        end
    else
        try, out=regexp(str,['"',fieldname,'"\s*:\s*"(.*?)"'],'tokens','once');
        end
        if isempty(out)
            try, out=regexp(str,['\.',fieldname,'"\s*:\s*"(.*?)"'],'tokens','once');
            end
        end
    end
    if isempty(out)
        try, 
            data=spm_jsonread(filename);
            out=data.(fieldname);
        end
    end
end
end