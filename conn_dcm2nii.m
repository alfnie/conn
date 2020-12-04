function [filenameout,filetypeout,fileRTout,filedescripout,filename]=conn_dcm2nii(filename,varargin)
% CONN_DCM2NII converts DICOM to NIFTI format
%
% CONN_DCM2NII [rootfilename]-1.dcm
% converts DICOM [rootfilename]-*.dcm files into [rootfilename].nii NIFTI format
%
% CONN_DCM2NII(...,optionname,optionvalue,...) sets CONN_DCM2NII options
%   valid optionname strings are:
%   'folderout' : output folder options
%                   './' to write nii files in same folder as DICOM files
%                   './nii' to write nii files in [DICOMFOLDER]/nii folder
%                   '../nii' to write nii files in [DICOMFOLDER]/../nii folder [default]
%                   use any other string to explicitly define the output folder
%   'overwrite' : 1/0 overwrites target nifti files if they already exist (default 0)
%   'opts'        see spm_dicom_convert help
%   'root_dir'    see spm_dicom_convert help
%   'format'      see spm_dicom_convert help
%   note: options defined using CONN_DCM2NII([],optionname,optionvalue,...) are PERSISTENT for the length of the current Matlab sessions (or until a "clear conn_dcm2nii;" command is issued)
%   note: options defined using CONN_DCM2NII(filename,optionname,optionvalue,...) only appy to the convertion of the file "filename" 
%

persistent saved;

if isempty(saved), 
    saved.folderout='../nii';
    saved.overwrite=false;
    saved.renamefiles=false;
    saved.spm_dicom_convert_opts='all';
    saved.spm_dicom_convert_root_dir='flat';
    saved.spm_dicom_convert_format=spm_get_defaults('images.format');
end
this=saved;

dogui=false;
logfile=[];
if nargin>1,
    for n=1:2:numel(varargin)-1
        switch(varargin{n})
            case 'folderout', this.folderout=varargin{n+1};
            case 'overwrite', this.overwrite=varargin{n+1};
            case 'renamefiles', this.renamefiles=varargin{n+1};
            case 'opts', this.spm_dicom_convert_opts=varargin{n+1};
            case 'root_dir', this.spm_dicom_convert_root_dir=varargin{n+1};
            case 'format', this.spm_dicom_convert_format=varargin{n+1};
            case 'logfile', logfile=varargin{n+1};
            case 'dogui', dogui=varargin{n+1};
            otherwise, error('unrecognized option ',varargin{n});
        end
    end
end
if isempty(filename), saved=this; return; end

conn_disp('fprintf','converting dcm files to nifti format...\n');

if isstruct(filename) % input from conn_dcmdir
else % input from filename
    ischarfilename=ischar(filename);
    if ischarfilename, filename={filename}; end
    filename=cellfun(@strtrim,filename,'uni',0);
    if all(cellfun('length',regexp(filename,'-1\.dcm$'))>0)
        for n0=1:numel(filename)
            allfilename=conn_dir(regexprep(filename{n0},'-1\.dcm$','-*.dcm'),'-R'); % find all -\d.dcm and sort by numbers
            if ~isempty(allfilename), filename{n0}=conn_sortfilenames(allfilename); end
        end
    else
        filename=cellstr(strvcat(filename{:})); % convert to conn_dcmdir output (structure)
        filename=conn_dcmdir(filename);
    end
end

ESSENTIALS=true;
filenameout={};
filetypeout=[];
fileRTout=[];
filedescripout={};
bak_filesout_path=0;
fhlog=[];
if dogui, ht=conn_waitbar(0,'converting dicom files to nifti format. please wait...'); end
if isstruct(filename)&&isfield(filename,'SeriesNumber')
    valid=find(arrayfun(@(n)numel(filename(n).SeriesNumber)==1,1:numel(filename)));
    usn=[filename(valid).SeriesNumber];
    uSN={};
    for n0=1:numel(valid)
        if nnz(usn==filename(valid(n0)).SeriesNumber)==1, uSN{n0}=sprintf('%d',filename(valid(n0)).SeriesNumber);
        else uSN{n0}=sprintf('%d%c',filename(valid(n0)).SeriesNumber,char('a'+sum(usn(1:n0-1)==filename(valid(n0)).SeriesNumber)));
        end
    end
    filename=filename(valid);
end
for n0=1:numel(filename)
    filenameout{n0}=[];
    filetypeout(n0)=0;
    fileRTout(n0)=nan;
    filedescripout{n0}='DICOM series incomplete';
    if isstruct(filename), 
        if isfield(filename,'check')&&(isempty(filename(n0).check)||~filename(n0).check),
            conn_disp('fprintf','skipping conversion of series #%d\n',n0);
            continue
        end
        tfilename=char(filename(n0).files);
    else
        tfilename=char(filename{n0});
    end
    [filesout_path,filesout_name,filesout_ext]=fileparts(deblank(tfilename(1,:)));
    if this.renamefiles
        filesout_name=sprintf('run-%s.nii',uSN{n0}); %filename(n0).SeriesNumber);
    else
        filesout_name=conn_prepend('',regexprep([filesout_name, filesout_ext],'-1\.dcm$',''),'.nii');
    end
    switch(lower(this.folderout))
        case './',
        case './nii',  filesout_path=fullfile(filesout_path,'nii');
            [ok,nill]=mkdir(filesout_path);
        case '../nii', filesout_path=fullfile(fileparts(filesout_path),'nii');
            [ok,nill]=mkdir(filesout_path);
        otherwise,
            if ~isempty(this.folderout)&&this.folderout(1)=='.', filesout_path=fullfile(filesout_path,this.folderout);
            else filesout_path=this.folderout;
            end
            [ok,nill]=mkdir(filesout_path);
    end
    if ~isempty(logfile)&&isequal(bak_filesout_path,0)
        bak_filesout_path=filesout_path;
        fhlog=fopen(logfile,'wt');
    end
    if isempty(logfile)&&~isequal(bak_filesout_path, filesout_path),
        bak_filesout_path=filesout_path;
        try, if ~isempty(fhlog), fclose(fhlog); end; end
        try, fhlog=fopen(fullfile(filesout_path,sprintf('conn_dcm2nii_%s.log',datestr(now,'yyyy_mm_dd_HHMMSSFFF'))),'wt'); end
    end
    filesout_name=conn_prepend('',filesout_name,'.nii');
    filesout_new=fullfile(filesout_path,filesout_name);
    
    if this.overwrite || ~conn_existfile(filesout_new)
        if isstruct(filename), 
            hdrs=filename(n0).headers;
            [info,rt]=conn_dcm2nii_header2json(filename(n0).SeriesInfo); 
        else
            hdrs=spm_dicom_headers(tfilename,false);
            [info,rt]=conn_dcm2nii_header2json(hdrs);
            if ESSENTIALS, for n=2:numel(hdrs), hdrs{n}=spm_dicom_essentials(hdrs{n}); end; end
        end
        fileRTout(n0)=str2double(rt);
        
        if nargin(@spm_dicom_convert)<5, % fix for spm8 and early spm12
            filesout=spm_dicom_convert(hdrs,this.spm_dicom_convert_opts,this.spm_dicom_convert_root_dir,this.spm_dicom_convert_format);
            movefiles=true;
            if isempty(filesout)||isempty(filesout.files), error('Unable to convert DICOM file. Please try installing a more recent SPM version'); end
        else % late spm12 and hopefully beyond
            filesout=spm_dicom_convert(hdrs,this.spm_dicom_convert_opts,this.spm_dicom_convert_root_dir,this.spm_dicom_convert_format,filesout_path);
            movefiles=false;
        end
        if isstruct(filesout)&&isfield(filesout,'files')&&~isempty(filesout.files)
            if movefiles
                for n=1:numel(filesout.files),
                    tfilename=filesout.files{n};
                    [nill,tfilename_name,tfilename_ext]=fileparts(filesout.files{n});
                    newtfilename=fullfile(filesout_path,[tfilename_name,tfilename_ext]);
                    if ~strcmp(newtfilename,tfilename)
                        for ext={'.nii','.mat'}
                            if ispc, [ok,nill]=system(['move "',conn_prepend('',tfilename,ext{1}),'" "',conn_prepend('',newtfilename,ext{1}),'"']);
                            else, [ok,nill]=system(['mv -f ''',conn_prepend('',tfilename,ext{1}),''' ''',conn_prepend('',newtfilename,ext{1}),'''']);
                            end
                        end
                        filesout.files{n}=newtfilename;
                    end
                end
            end
            if numel(filesout.files)==1 % one 3d output file
                if ispc, [ok,nill]=system(['move "',char(filesout.files),'" "',filesout_new,'"']);
                else, [ok,nill]=system(['mv -f ''',char(filesout.files),''' ''',filesout_new,'''']);
                end
                if ispc, [ok,nill]=system(['move "',regexprep(char(filesout.files),'\.nii$','.mat'),'" "',regexprep(filesout_new,'\.nii$','.mat'),'"']);
                else, [ok,nill]=system(['mv -f ''',regexprep(char(filesout.files),'\.nii$','.mat'),''' ''',regexprep(filesout_new,'\.nii$','.mat'),'''']);
                end
                filenameout{n0}=filesout_new;
                filetypeout(n0)=1;
                try,
                    a=spm_vol(filesout_new);
                    str=evalc('disp(hdrs{1})');
                    spm_jsonwrite(conn_prepend('',filesout_new,'.json'),info,struct('indent',' '));
                    fprintf(fhlog,'created 3d nifti file %s\n  Description: %s\n  RT: %s\n  Dimensions: %s\n  Mapping: %s\n%s\n',filesout_new,a.descrip,rt,mat2str(a.dim),mat2str(a.mat),str);
                          conn_disp('fprintf','created 3d nifti file %s\n  Description: %s\n  RT: %s\n  Dimensions: %s\n  Mapping: %s\n',filesout_new,a.descrip,rt,mat2str(a.dim),mat2str(a.mat));
                end
            else
                a=spm_vol(char(filesout.files));
                ok=false;
                try,
                    tdim=cat(1,a.dim);
                    ok=~any(any(diff(tdim,1,1),1),2);
                end
                if ok, %spm_check_orientations(a,false) % one 4d output file
                    spm_file_merge(a,filesout_new);
                    spm_unlink(filesout.files{:});
                    filenameout{n0}=filesout_new;
                    filetypeout(n0)=numel(a);
                    str=evalc('disp(hdrs{1})');
                    spm_jsonwrite(conn_prepend('',filesout_new,'.json'),info,struct('indent',' '));
                    try, fprintf(fhlog,'created 4d nifti file %s (%d volumes)\n  Description: %s\n  RT: %s\n  Dimensions: %s\n  Mapping_first: %s\n  Mapping_last: %s\n%s\n',filesout_new,numel(a),a(1).descrip,rt,mat2str(a(1).dim),mat2str(a(1).mat),mat2str(a(end).mat),str); end
                               conn_disp('fprintf','created 4d nifti file %s (%d volumes)\n  Description: %s\n  RT: %s\n  Dimensions: %s\n  Mapping_first: %s\n  Mapping_last: %s\n',filesout_new,numel(a),a(1).descrip,rt,mat2str(a(1).dim),mat2str(a(1).mat),mat2str(a(end).mat));
                else
                    filenameout{n0}={};
                    for n1=1:numel(filesout.files) % multiple 3d output files
                        tfilesout_new=conn_prepend('',filesout_new,['_',num2str(n1,'%04d'),'.nii']);
                        if ispc, [ok,nill]=system(['move "',char(filesout.files{n1}),'" "',tfilesout_new,'"']);
                        else, [ok,nill]=system(['mv -f ''',char(filesout.files{n1}),''' ''',tfilesout_new,'''']);
                        end
                        if ispc, [ok,nill]=system(['move "',regexprep(char(filesout.files{n1}),'\.nii$','.mat'),'" "',regexprep(tfilesout_new,'\.nii$','.mat'),'"']);
                        else, [ok,nill]=system(['mv -f ''',regexprep(char(filesout.files{n1}),'\.nii$','.mat'),''' ''',regexprep(tfilesout_new,'\.nii$','.mat'),'''']);
                        end
                        filenameout{n0}{n1}=tfilesout_new;
                    end
                    str=evalc('disp(hdrs{1})');
                    spm_jsonwrite(conn_prepend('',filesout_new,'.json'),info,struct('indent',' '));
                    try, fprintf(fhlog,'created multiple 3d nifti files %s-#### (%d volumes)\n  Description: %s\n  RT: %s\n  Dimensions_first: %s\n  Mapping_first: %s\n  Dimensions_last: %s\n  Mapping_last: %s\n%s\n',filesout_new,numel(a),a(1).descrip,rt,mat2str(a(1).dim),mat2str(a(1).mat),mat2str(a(end).dim),mat2str(a(end).mat),str); end
                               conn_disp('fprintf','created multiple 3d nifti files %s-#### (%d volumes)\n  Description: %s\n  RT: %s\n  Dimensions_first: %s\n  Mapping_first: %s\n  Dimensions_last: %s\n  Mapping_last: %s\n',filesout_new,numel(a),a(1).descrip,rt,mat2str(a(1).dim),mat2str(a(1).mat),mat2str(a(end).dim),mat2str(a(end).mat));
                    filenameout{n0}=char(filenameout{n0});
                    filetypeout(n0)=numel(a);
                end
            end
            try, filedescripout{n0}=sprintf('Series #%d : %s (%d%s, %d volumes) %s %s',hdrs{1}.SeriesNumber,hdrs{1}.SeriesDescription,a(1).dim(1),sprintf('x%d',a(1).dim(2:end)),numel(a),hdrs{1}.ImageType,a(1).descrip);
            catch, filedescripout{n0}='DICOM series incomplete';
            end
        else
            filenameout{n0}=[];
            filetypeout(n0)=0;
            fileRTout(n0)=nan;
            filedescripout{n0}='DICOM series incomplete';
        end
    else
        a=spm_vol(filesout_new);
        filenameout{n0}=filesout_new;
        filetypeout(n0)=numel(a);
        fileRTout(n0)=nan;
        filedescripout{n0}='unknown';
    end
    if dogui, conn_waitbar(n0/numel(filename),ht,filedescripout{n0}); end
end
try, fclose(fhlog); end
if ~isstruct(filename)&&ischarfilename,
    while iscell(filenameout)&&numel(filenameout)==1, filenameout=filenameout{1}; end
    try, filenameout=char(filenameout); end;
end
if dogui, conn_waitbar('close',ht); end
conn_disp('fprintf','done\n');
end

function [info,rt]=conn_dcm2nii_header2json(info)
rt='unknown'; 
info0=info;
if iscell(info), info=info{1}; end
if isfield(info,'RepetitionTime'), info.RepetitionTime=round(1e1*info.RepetitionTime)/1e4; rt=mat2str(info.RepetitionTime); end % ms to seconds
if isfield(info,'EchoTime'), info.EchoTime=round(1e1*info.EchoTime)/1e4; end                                                    % ms to seconds
if ~isfield(info,'SliceTiming') % slice-timing information info
    if isfield(info,'MosaicRefAcqTimes')&&~isempty(info.MosaicRefAcqTimes), 
        info.SliceTiming=round(1e1*info.MosaicRefAcqTimes)/1e4;
        info=rmfield(info,'MosaicRefAcqTimes');
    elseif isfield(info,'Private_0019_1029')&&~isempty(info.Private_0019_1029), 
        if ~rem(numel(info.Private_0019_1029),8)&&all(rem(info.Private_0019_1029,1)==0&info.Private_0019_1029>=0&info.Private_0019_1029<=255), info.SliceTiming=typecast(uint8(info.Private_0019_1029),'double')/1e3;
        else info.SliceTiming=round(1e1*info.Private_0019_1029)/1e4;
        end
        info=rmfield(info,'Private_0019_1029');
    end
end
if ~isfield(info,'TotalReadoutTime') % total readout time info
    if isfield(info,'BandwidthPerPixelPhaseEncode')&&numel(info.BandwidthPerPixelPhaseEncode)==1, 
        info.TotalReadoutTime=1/info.BandwidthPerPixelPhaseEncode;
        info=rmfield(info,'BandwidthPerPixelPhaseEncode');
    elseif isfield(info,'Private_0019_1028')&&~isempty(info.Private_0019_1028),
        if numel(info.Private_0019_1028)==8&&all(rem(info.Private_0019_1028,1)==0&info.Private_0019_1028>=0&info.Private_0019_1028<=255), info.TotalReadoutTime=1./typecast(uint8(info.Private_0019_1028),'double');
        else info.TotalReadoutTime=1./info.Private_0019_1028;
        end
        info=rmfield(info,'Private_0019_1028');
    end
end
info.ConversionSoftware=sprintf('CONN%s & %s',conn('ver'),spm('ver'));
info=conn_dcm2nii_expand(info);
end

function info2=conn_dcm2nii_expand(info,fieldname,info2)
if nargin<3, info2=struct; end
if nargin<2, fieldname=''; end
if isempty(info), return; end
if iscell(info)||(isstruct(info)&&numel(info)>1)
    for n=1:numel(info)
        if iscell(info), info2=conn_dcm2nii_expand(info{n},sprintf('%s_%04d',fieldname,n),info2);
        else info2=conn_dcm2nii_expand(info(n),sprintf('%s_%04d',fieldname,n),info2);
        end
    end
    return
end
if isstruct(info)
    names=fieldnames(info);
    for n=1:numel(names)
        if isempty(fieldname), info2=conn_dcm2nii_expand(info.(names{n}),names{n},info2);
        else info2=conn_dcm2nii_expand(info.(names{n}),sprintf('%s_%s',fieldname,names{n}),info2);
        end
    end
    return
end
if ~isempty(fieldname), 
    if ischar(info), info=regexprep(info,'\W+',' '); end
    info2.(fieldname)=info; 
end
end

    

