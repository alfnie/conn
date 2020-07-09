function filesout=conn_expandframe(filesin)
% expands all frames in 4d nifti file

singlefile=false;
if ~iscell(filesin), 
    if isstruct(filesin) % special-case, passing spm_vol structure
        filesout=cellfun(@(a,b)[a,',',b],{filesin.fname},cellfun(@(x)num2str(x(1)),{filesin.n},'uni',0),'uni',0);
        return
    end
    filesin=cellstr(filesin); % passing single string
    singlefile=true;
end

filesout={};
for n=1:numel(filesin)
    filename=filesin{n};
    if ~isempty(regexp(filename,',\d+$','once'))
        filesout=cat(1,filesout, {filename});
    else
        ni = nifti(filename);
        dm = [ni.dat.dim 1 1 1 1];
        if dm(4)>1
            filesout=cat(1,filesout, arrayfun(@(x)[filename,',',num2str(x)],(1:dm(4))','uni',0) );
        else
            filesout=cat(1,filesout, {filename});
        end
    end
end

if singlefile
    filesout=char(filesout);
end

